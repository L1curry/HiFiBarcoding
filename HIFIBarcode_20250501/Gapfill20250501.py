'''
间隙 (gap)填充
1) 输入数据:带标签的双端拼接结果作为两端信息输入,中间无标签的拼接结  果统一作为中间信息输入。
2) 组装算法:采用 SOAPBarcode 算法 (图 4) 对条形码序列进行填充,从而  获得高精度的完整的昆虫 COI 基因条形码。
a) 每一对两端数据的前端为起  始点,末端为终点,并用 kmer 构建 de bruijin 图,从起点到终点查找潜在  的连接路径;
b) 采用以下策略保证连接正确性:第一,删除在分叉处 kmer丰度小于平均 kmer 丰度 10%的路径 (蓝色);
第二,如果经过第一步过滤后还有多条出度 (out degree,OD) 存在,则针对不同 ODs 和位于最后分叉前的 kmer 间进行 reads 数目统计,后删除丰度小于平均 reads丰度 10%的ODs(黑色);
第三,删除超出了预先设定的长度的路径 (红色)。
组装过程慢，end sequence 和 middle sequence 序列太多，需要对两者进行精简，end sequence在gapfill这一步按样本分组时进行精简，
而middle sequence应该在chain过程中或chain完成后对其进行精简。
'''
import glob
import os
import time
import subprocess
import argparse
import logging
from datetime import datetime

# 初始化参数解析
parser = argparse.ArgumentParser(
    add_help=True,
    description='Gap filling script with logging and abundance filtering',
    formatter_class=argparse.RawTextHelpFormatter
)
# 参数定义
group = parser.add_argument_group("common arguments")
# 输出前缀
group.add_argument(
    "-outpre", metavar="<STR>", required=True,
    help="Output prefix. All output directories and files will be named using this prefix."
)
# 最大组装长度
group.add_argument(
    "-gmax", metavar="<INT>", type=int, dest="gmax", default=700,
    help="Maximum assembled sequence length (bp). Default: 700."
)
# 最小组装长度
group.add_argument(
    "-gmin", metavar="<INT>", type=int, dest="gmin", default=450,
    help="Minimum assembled sequence length (bp). Default: 450."
)
# 低复杂度 kmer 阈值
group.add_argument(
    "-kl", metavar="<INT>", type=int, dest="lowkmer", default=55,
    help="Low-complexity k-mer abundance threshold. Default: 55."
)
# kmer 区间
group.add_argument(
    "-i", metavar="<INT>", type=int, dest="kmer_interval", default=10,
    help="K-mer size interval. Default: 10."
)
# kmer 长度
group.add_argument(
    "-k", metavar="<INT>", type=int, dest="kmer", default=127,
    help="K-mer length for de Bruijn graph. Default: 127."
)
# 每子文件样本数
group.add_argument(
    "-S", metavar="<INT>", type=int, dest="samp_num", default=1,
    help="Number of samples per split ends file. Default: 1."
)
# 线程数
group.add_argument(
    "-t", metavar="<INT>", type=int, dest="cpu", default=4,
    help="Number of CPU threads to use. Default: 4."
)
# ends 序列文件
group.add_argument(
    "-ends", metavar="<STR>", type=str, required=True, dest="gapfill_ends",
    help="Ends FASTA file path (paired-end merged sequences with For/Rev tags)."
)
# middle 序列文件
group.add_argument(
    "-middle", metavar="<STR>", type=str, required=True, dest="gapfill_mid",
    help="Middle FASTA file path (unlabeled contig sequences)."
)
# barcode 可执行文件路径
group.add_argument(
    "-barcode", metavar="<STR>", required=False, default="barcode",
    help="Path to SOAPBarcode executable. Default 'barcode' (must be in system PATH)."
)
args = parser.parse_args()  # 解析命令行参数

# 创建主要输出目录
# 主输出目录以 outpre_gapfill 命名
gapfill_outdir = os.path.abspath(args.outpre + "_gapfill")
if not os.path.exists(gapfill_outdir): os.mkdir(gapfill_outdir)
# 初始化日志：放在 gapfill_outdir 下
start_time = datetime.now()
log_filename = os.path.join(gapfill_outdir, f"gapfill_{start_time.strftime('%Y%m%d_%H%M%S')}.log")
logging.basicConfig(
    filename=log_filename,
    level=logging.INFO,
    format='%(asctime)s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logging.info("Gapfill run starts")

# 目录创建工具
def checkDirMkdir(fold):
    outdir = os.path.abspath(fold)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    return outdir

# FASTA 解析
def parse_fasta(fh):
    while True:
        header = fh.readline().strip()
        if not header or not header.startswith('>'):
            break
        seq = fh.readline().strip()
        yield header[1:], seq

# 处理 ends 并拆分，包含丰度过滤
def format_ends(ends, subsam, split_ends_dir):
    endfile = []
    seqsf, seqsr, abundance = {}, {}, {}
    all_samples = []
    # 读取序列与丰度
    with open(ends) as fh:
        for head, seq in parse_fasta(fh):
            parts = head.split('_')
            nam = parts[0].replace('For', '').replace('Rev', '')
            if nam not in all_samples: all_samples.append(nam)
            ab = int(parts[2])
            if head.startswith('For'):
                seqsf.setdefault(nam, []).append(seq)
            else:
                seqsr.setdefault(nam, []).append(seq)
            abundance[seq] = ab
    sorted_samples = sorted(all_samples)
    # 记录初始候选对数
    for s in sorted_samples:
        total_pairs = len(seqsf.get(s, [])) * len(seqsr.get(s, []))
        logging.info(f"Sample {s}: initial candidate pairs = {total_pairs}")
    # 拆分逻辑
    sub_len = int(subsam)
    current_len, file_mark = 0, 1
    current_cont = ''
    for s in sorted_samples:
        if s not in seqsf or s not in seqsr:
            logging.warning(f"Sample {s} missing forward or reverse; skipped")
            continue
        current_len += 1
        order_sam = 0
        for f_idx, fseq in enumerate(seqsf[s], 1):
            for r_idx, rseq in enumerate(seqsr[s], 1):
                order_sam += 1
                abf, abr = abundance[fseq], abundance[rseq]
                if abf > abr * 10 or abr > abf * 10:
                    logging.info(f"Sample {s} pair f#{f_idx}-r#{r_idx} dropped: abundance {abf}-{abr}")
                    continue
                header = f">{s}_{f_idx}-{r_idx};{order_sam};{abf}-{abr}\n"
                current_cont += header + fseq + "\n" + rseq + "\n"
        if current_len >= sub_len and current_cont:
            out = f"{split_ends_dir}/ends.{file_mark}"
            with open(out, 'w') as of: of.write(current_cont)
            kept = current_cont.count('>')
            logging.info(f"ends.{file_mark} written, kept pairs = {kept}")
            endfile.append(out)
            current_cont, current_len = '', 0
            file_mark += 1
    # 写入剩余
    if current_cont:
        out = f"{split_ends_dir}/ends.{file_mark}"
        with open(out, 'w') as of: of.write(current_cont)
        kept = current_cont.count('>')
        logging.info(f"ends.{file_mark} written, kept pairs = {kept}")
        endfile.append(out)
    return endfile

# 进一步创建子目录
result_outdir = checkDirMkdir(gapfill_outdir + "/result")
split_ends_dir = checkDirMkdir(gapfill_outdir + "/ends")
shell_outdir = checkDirMkdir(gapfill_outdir + "/shell")
# 写 middle.lis
with open(os.path.join(gapfill_outdir, "middle.lis"), 'w') as mf:
    mf.write(">\n" + "f=" + os.path.abspath(args.gapfill_mid))

# 拆分 ends 文件
if args.samp_num < 1:
    raise ValueError("samp_num must be >=1")
split_ends = format_ends(args.gapfill_ends, args.samp_num, split_ends_dir)
# 记录总输出数
total_out = sum(open(f).read().count('>') for f in split_ends)
logging.info(f"Total output records = {total_out}")

# 顺序执行每个 shell 脚本，确保前一个完成后再启动下一个
for idx, f in enumerate(split_ends, 1):
    shell = f"{shell_outdir}/barcodes.{idx}.sh"
    out_pref = os.path.abspath(f"{result_outdir}/barcodes.{idx}")
    with open(shell, 'w') as sh:
        sh.write(f"#!/bin/sh\n{args.barcode}"
                 f" -e {os.path.abspath(f)}"
                 f" -r {os.path.join(gapfill_outdir, 'middle.lis')}"
                 f" -o {out_pref}"
                 f" -x {args.gmax}"
                 f" -n {args.gmin}"
                 f" -l {args.lowkmer}"
                 f" -v {args.kmer_interval}"
                 f" -k {args.kmer}"
                 f" -t {args.cpu}\n")
    subprocess.run(['chmod', '+x', shell])
    logging.info(f"Running shell script {shell}")
    subprocess.run(['bash', shell], check=True)

# 结束日志
end_time = datetime.now()
elapsed = (end_time - start_time).total_seconds()
logging.info("Gapfill run ends")
logging.info(f"Total elapsed time: {elapsed:.2f}s")