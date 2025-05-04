#!/usr/bin/env python3
"""
脚本的输入包括：
contigDir 目录下两个子目录：
ends/：存放每个样本原始序列的编号–名称对应文件（按行编号），如 sampleA.ends；
result/：存放每个样本的 contig 文件，格式为 barcodes.<suffix>.contig。
（第二阶段）all_barcodes_raw.fasta：由第一阶段合并产物，可直接用作 QC 流水线输入。
一份引物文件（FASTA 格式），包含 forward 和 reverse 引物序列，用于修剪步骤。
用户可通过命令行参数指定：
覆盖度阈值（如 -cov 5）
片段长度阈值（如 -l 650）
是否开启翻译检查（--cc 关闭）及线粒体密码表编号（-codon 5）。

脚本的输出包括：
合并后的原始条形码 FASTA
stage1 生成：all_barcodes_raw.fasta，所有样本 contig 统一重命名后写入。
质控后的清洗 FASTA
stage2 输出：仅保留经引物修剪、覆盖度/长度/可选翻译检查及 ≥95% 去重后的高质量序列。

日志文件
记录 INFO 与 WARNING 级别的运行详情和统计，文件名带有时间戳。
脚本的主要功能可概括为：
Stage1 合并：将每个样本的 contig 按原始序列编号映射回样本名，并合并到统一 FASTA；
Stage2 QC 流水线：引物修剪 → 覆盖度 & 长度 过滤 →（可选）“两个方向 × 三个阅读框”翻译检查 → 丰度排序 + ≥95% 同一性去重 → 输出最终清洗序列。

日志信息含义示例：
2025-05-03 22:15:55,602 WARNING: Missing order 26 in file 'barcodes.A1.contig' (sample suffix='A1')
表示在 stage1 合并第 26 号 contig 时，ends/A1.ends 中没有对应的第 26 行记录，故跳过该条序列。
2025-05-04 08:57:24,444 WARNING: Translation failed in all frames and RC for sample='A1', seq_id='A1_2-2;7;232-363;k=97_22'
表示在 stage2 的翻译检查中，对于样本 A1 的序列 A1_2-2;7;232-363;k=97_22，其正链和互补链的所有三种阅读框都出现早期终止子，因而被计为“翻译丢弃”并不输出至最终 FASTA。
通过这些输入/输出和日志说明，你可以完整理解脚本如何从原始 contig 目录经过多步 QC，生成高质量的 COI 条形码序列，并通过日志精确追踪任何丢弃原因。
"""
import os
import time
import argparse
import logging
from Bio.Seq import Seq
from natsort import natsorted
import shutil
import sys
import glob

# Utility
def checkDirMkdir(path):
    outdir = os.path.abspath(path)
    os.makedirs(outdir, exist_ok=True)
    return outdir

# Hamming distance
def hamming(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

# Find primer in seq (allow <=1 mismatch)
def find_primer(seq, primer):
    L = len(primer)
    for i in range(len(seq) - L + 1):
        if hamming(seq[i:i+L], primer) <= 1:
            return i
    return None

# Trim primers: extract region between forward and reverse primers
def trim_primers(seq, fwd, rev, logger):
    # try original and reverse complement orientations
    for orientation in (seq, str(Seq(seq).reverse_complement())):
        pos_fwd = find_primer(orientation, fwd)
        if pos_fwd is None:
            continue
        after_fwd = orientation[pos_fwd + len(fwd):]
        # find reverse primer in this region
        pos_rev = find_primer(after_fwd, rev)
        if pos_rev is not None:
            return after_fwd[:pos_rev]
        # try reverse complement of after_fwd
        rc_after = str(Seq(after_fwd).reverse_complement())
        pos_rev = find_primer(rc_after, rev)
        if pos_rev is not None:
            return rc_after[:pos_rev]
    logger.warning("Primers not found for seq: both orientations")
    return None

# Enhanced translation check
def passes_translation(seq, codon_table, logger):
    # try original and reverse complement
    for attempt in range(2):
        s = seq if attempt == 0 else str(Seq(seq).reverse_complement())
        for frame in range(3):
            sub = s[frame:]
            sub = sub[:len(sub) - (len(sub) % 3)]
            prot = Seq(sub).translate(table=codon_table)
            if '*' not in prot:
                return True, s
    # 如果失败，记录是哪个样本、哪个序列出了问题
    # context = getattr(logger, 'context_info', {})
    # sample = context.get('sample', 'unknown_sample')
    # seq_id = context.get('seq_id', 'unknown_seq')
    # logger.warning(
    #     f"Translation failed in all frames and RC for sample='{sample}', sequence_id='{seq_id}'"
    # )
    return False, None

# Deduplication identity
def seq_identity(s1, s2):
    l = min(len(s1), len(s2))
    return sum(a == b for a, b in zip(s1, s2)) / l if l else 0

# Parse primers file
def load_primers(path):
    fwd = rev = None
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                if parts[0].lower() == 'forward':
                    fwd = parts[1]
                elif parts[0].lower() == 'reverse':
                    rev = parts[1]
    if not fwd or not rev:
        raise ValueError("Primer file must contain 'forward' and 'reverse' sequences.")
    return fwd, rev

# Parse FASTA
def parse_fasta(fh):
    header = None
    seq_lines = []
    for line in fh:
        line = line.rstrip()
        if line.startswith('>'):
            if header:
                yield header, ''.join(seq_lines)
            header = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
    if header:
        yield header, ''.join(seq_lines)

# Stage1: merge haplotypes into raw FASTA

def stage1(contigDir, outdir, logger, raw_fa_path):
    with open(raw_fa_path, 'w') as raw_out:
        ends_dir = os.path.join(contigDir, 'ends')
        res_dir = os.path.join(contigDir, 'result')
        for fn in os.listdir(ends_dir):
            if fn.startswith('.'):
                continue
            suffix = fn.split('.')[-1]
            order = {}
            idx = 0
            with open(os.path.join(ends_dir, fn)) as f:
                for line in f:
                    if line.startswith('>'):
                        idx += 1
                        order[idx] = line.strip()[1:]
            cf = os.path.join(res_dir, f"barcodes.{suffix}.contig")
            if not os.path.exists(cf):
                logger.warning(f"Missing {cf}")
                continue
            with open(cf) as fh:
                for head, seq in parse_fasta(fh):
                    parts = head.split()[0].rsplit('_', 1)
                    if len(parts) != 2 or not parts[1].isdigit():
                        logger.warning(f"Unexpected contig header: {head}")
                        continue
                    kmer, idx = parts[0], int(parts[1])
                    orig = order.get(idx)
                    if not orig:
                        # 在警告中同时输出出错的 contig 文件名 和 样本后缀（suffix）
                        logger.warning(
                            f"Missing order {idx} in file 'barcodes.{suffix}.contig' (sample suffix='{suffix}')"
                        )
                        continue
                    raw_out.write(f">{orig};k={kmer}\n{seq}\n")
    return raw_fa_path

# Stage2: polish w/ primer trimming

def stage2(raw_fa, outdir, args, logger):
    polished = os.path.join(outdir, 'all.barcodes.polished.fasta')
    raw_counts, polished_counts, disc = {}, {}, {}
    sample_seqs = {}

    fwd, rev = load_primers(args.primer)
    rev = str(Seq(rev).reverse_complement())

    with open(raw_fa) as fh:
        for head, seq in parse_fasta(fh):
            sample = head.split(';')[0].split('_')[0]
            raw_counts[sample] = raw_counts.get(sample, 0) + 1
            disc.setdefault(sample, {'prim': 0, 'cov': 0, 'len': 0, 'trans': 0, 'dedup': 0})

            # Extract kmer
            kmer = None
            for token in head.split(';'):
                if token.startswith('k='):
                    kmer = token.split('=', 1)[1]
            if not kmer:
                logger.warning(f"kmer missing in header, skipping: {head}")
                disc[sample]['prim'] += 1
                continue

            # Trim primers
            trimmed = trim_primers(seq, fwd, rev, logger)
            if not trimmed:
                disc[sample]['prim'] += 1
                continue

            # Coverage
            cov_item = next((f for f in head.split(';') if '-' in f and all(p.isdigit() for p in f.split('-'))), None)
            if not cov_item:
                logger.warning(f"Coverage missing/malformed in header: {head}")
                disc[sample]['cov'] += 1
                continue
            c1, c2 = map(int, cov_item.split('-'))
            if c1 < args.min_coverage or c2 < args.min_coverage:
                disc[sample]['cov'] += 1
                continue

            if len(trimmed) < args.coi_length:
                disc[sample]['len'] += 1
                continue

            if args.coi_check:
                ok, usable = passes_translation(trimmed, args.codon_table, logger)
                if not ok:
                    # 直接在这里记录 sample 和 seq_id（head）
                    seq_id = head.split()[0]
                    logger.warning(
                        f"Translation failed in all frames and RC for "
                        f"sample='{sample}', seq_id='{seq_id}'"
                    )
                    disc[sample]['trans'] += 1
                    continue
            else:
                usable = trimmed

            abu = (c1 + c2) // 2
            sample_seqs.setdefault(sample, []).append((usable, abu, kmer))

    # Dedup and write
    with open(polished, 'w') as out:
        for sample in natsorted(sample_seqs):
            polished_counts[sample] = 0
            kept = []
            for seq, abu, kmer in sorted(sample_seqs[sample], key=lambda x: x[1], reverse=True):
                if any(seq_identity(seq, k) >= 0.97 for k in kept):
                    disc[sample]['dedup'] += 1
                else:
                    kept.append(seq)
                    tag = 'size' if polished_counts[sample] == 0 else 'abu'
                    out.write(f">{sample};{tag}={abu};k={kmer};l={len(seq)}\n{seq}\n")
                    polished_counts[sample] += 1

    for sample in natsorted(set(raw_counts) | set(polished_counts)):
        logger.info(f"Sample {sample}: raw={raw_counts.get(sample,0)}, polished={polished_counts.get(sample,0)}, discards={disc[sample]}")
    return polished

def collect_results(outpre):
    """
    汇总并移动各阶段结果到统一目录 outpre_results，
    并删除原目录。
    若指定文件未找到，停止运行并提示用户检查。
    如果某目录不存在，则打印警告并跳过该目录。
    """
    # 定义各阶段目录及移动文件，使用通配符匹配 .log
    dirs = [
        (f"{outpre}_assign", ["*.log", "assign_err.fasta", "assign.list", "assign_middle.ssam"]),
        (f"{outpre}_buildend", ["*.log", "buildends.fasta"]),
        (f"{outpre}_chain", ["*.log", "chain.fasta"]),
        # 整个目录移动
        (f"{outpre}_gapfill", None),
        (f"{outpre}_polish", None),
    ]
    # 创建结果目录
    res_dir = f"{outpre}_results"
    os.makedirs(res_dir, exist_ok=True)

    for src, files in dirs:
        if not os.path.exists(src):
            print(f"警告：目录 '{src}' 不存在，已跳过该目录。")
            continue

        if files is None:
            # 整个目录移动
            dst = os.path.join(res_dir, os.path.basename(src))
            try:
                shutil.move(src, dst)
            except Exception as e:
                print(f"移动目录 {src} 时出错：{e}")
                sys.exit(1)
        else:
            # 先检查所有必需文件是否存在
            missing = False
            for pattern in files:
                matches = glob.glob(os.path.join(src, pattern))
                if not matches:
                    print(
                        f"错误：在目录 '{src}' 中未找到匹配 '{pattern}' 的文件。",
                        f"请检查 {src}/{pattern} 是否存在。"
                    )
                    missing = True
            if missing:
                sys.exit(1)

            # 若均存在，执行移动
            for pattern in files:
                for src_path in glob.glob(os.path.join(src, pattern)):
                    try:
                        shutil.move(src_path, res_dir)
                    except Exception as e:
                        print(f"移动文件 {src_path} 时出错：{e}")
                        sys.exit(1)

            # 删除整个目录及其残余文件
            try:
                shutil.rmtree(src)
            except Exception as e:
                print(f"删除目录 {src} 失败：{e}")
                sys.exit(1)

    print(f"所有结果已汇总到目录：{res_dir}")


# Main
if __name__ == '__main__':
    start = time.time()
    p = argparse.ArgumentParser()
    p.add_argument('-d', dest='contigDir', required=True)
    p.add_argument('-outpre', required=True)
    p.add_argument('-primer', dest='primer', required=True)
    p.add_argument('-i', dest='coi_input')
    p.add_argument('-cov', type=int, dest='min_coverage', default=5)
    p.add_argument('-l', type=int, dest='coi_length', default=650)
    p.add_argument('-codon', type=int, dest='codon_table', default=5)
    p.add_argument('-cc', dest='coi_check', action='store_false')
    args = p.parse_args()

    outdir = checkDirMkdir(args.outpre + '_polish')
    logfn = os.path.join(outdir, 'polish_' + time.strftime('%Y%m%d_%H%M%S') + '.log')
    logging.basicConfig(filename=logfn, level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')
    log = logging.getLogger()
    log.info(f"Run start {time.strftime('%Y-%m-%d %H:%M:%S')}")

    raw_fa = args.coi_input if args.coi_input else os.path.join(outdir, 'all_barcodes_raw.fasta')
    if not args.coi_input:
        raw_fa = stage1(args.contigDir, outdir, log, raw_fa)
    polished = stage2(raw_fa, outdir, args, log)

    log.info(f"Run end {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log.info(f"Elapsed {time.time() - start:.2f}s")
    print(f"Done: raw={raw_fa}, polished={polished}, log={logfn}")
    collect_results(args.outpre)