import os
import time
import gzip
import subprocess
import argparse
from Bio.Seq import Seq
import shutil

###############################################################################
#####------------------------- parameters --------------------------------#####
parser = argparse.ArgumentParser(
    add_help=True,
    description="去冗余:将拆分得到的含有标签引物的单端序列按照98%的相似度聚类,减少冗余数据。\n"
                "丰度排序与过滤:根据去冗余的聚类结果将序列按照丰度高低进行排序,并记录相关丰度信息;根据丰度信息保留不同类别序列进行下游分析。\n",
    formatter_class=argparse.RawTextHelpFormatter
)

group = parser.add_argument_group("common arguments")
group.add_argument("-outpre", metavar="<STR>", required=True, help="输出文件名前缀")
group.add_argument("-index", metavar="<INT>", type=int, required=True, help="引物标签长度")
group.add_argument("-vsearch", metavar="<STR>", required=False, default=None, help="vsearch路径，若不指定将自动查找")
group.add_argument("-threads", metavar="<INT>", type=int, default=2, help="vsearch线程数，默认2")
group.add_argument("-cid", metavar="<FLOAT>", type=float, default=0.98, dest="cluster_identity", help="聚类相似度阈值，默认0.98")
group.add_argument("-min", metavar="<INT>", type=int, default=10, dest="min_overlap", help="最小重叠长度，默认10bp")
group.add_argument("-max", metavar="<INT>", type=int, default=150, dest="max_overlap", help="最大重叠长度，默认150bp")
group.add_argument("-tp", metavar="<INT>", type=int, dest="cluster_number_needKeep", help="保留簇数，可手动指定，默认前5")
group.add_argument("-ab", metavar="<INT>", type=int, dest="abundance_threshod", help="丰度阈值，可手动指定，默认保留前5簇")
group.add_argument("-seqs_lim", metavar="<INT>", type=int, default=0, help="ssam文件最大处理序列数，默认0，表示不限制")
group.add_argument("-len", metavar="<INT>", type=int, default=150, dest="standard_length", help="标准读长，默认150bp")
group.add_argument("-ds", dest="drop_short_read", action="store_true", help="丢弃短reads，默认不开启")
group.add_argument("-rc", dest="reads_check", action="store_true", help="对rads进行翻译检查，默认不开启")
group.add_argument("-bom", metavar="<INT>", type=int, dest="buildend_mismatch", default=1, help="允许错配数，默认1个碱基")
group.add_argument("-list", metavar="<FILE>", type=str, required=True, help="输入.ssam文件列表")
args = parser.parse_args()

# 自动查找 vsearch
if args.vsearch is None:
    found = shutil.which("vsearch")
    if found:
        args.vsearch = found
        print(f"[INFO] auto-detected vsearch at {found}")
    else:
        print("[ERROR] vsearch not found; please install or specify -vsearch <path>")
        exit(1)

# 创建输出目录
buildend_outdir = os.path.abspath(args.outpre + "_buildend")
if not os.path.exists(buildend_outdir):
    os.makedirs(buildend_outdir)

# 新日志（仅一次定义）
start_time = time.time()
t = start_time
start_ts = time.strftime("%Y%m%d_%H%M%S", time.localtime(start_time))
runlog = os.path.join(buildend_outdir, f"buildend_{start_ts}.log")

def write_runlog(msg):
    with open(runlog, 'a') as log:
        log.write(msg + "\n")

# 初始化日志
write_runlog(f"START: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}")
write_runlog(f"Parameters: {args}")


# -----------------------functions for common---------------------#
def open_input(file):
    if file.endswith("gz"):
        return gzip.open(file, 'rt')
    else:
        return open(file, 'r')

def checkDirMkdir(fold):
    if not os.path.exists(fold):
        os.mkdir(fold)
    return os.path.abspath(fold)

def print_time(info):
    print(info + " " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

def run_time(info):
    print(f"[INFO]: {info} total run time: {time.time() - t:.2f}s")



# ----------------------functions for buildend------------------#
def complementation(sequence):
    sequence = sequence.upper()
    transtable = str.maketrans('ATCG', 'TAGC')
    return sequence.translate(transtable)

def comp_rev(sequence):
    return complementation(sequence)[::-1]


def parse_ssam(ssam):
    while True:
        tmp = ssam.readline().strip().split()
        name = tmp[0]
        if not name:
            break
        read1, read2 = tmp[1], tmp[3]
        qual1, qual2 = tmp[2], tmp[4]
        if len(read1) == len(read2) and len(qual1) == len(qual2):
            yield name, read1, read2, qual1, qual2


def comp_rev_list(reads_list):
    # make a list of sequences reverse and complement #
    new_reads_list = []
    for read in reads_list:
        read = comp_rev(read)
        new_reads_list.append(read)

    return new_reads_list


def mismatch(str1, str2):
    # ----count matched bases of two sequences----#
    mismatches = 0
    for base in range(len(str1)):
        if str1[base] != str2[base]:
            mismatches += 1

    return mismatches


def translate_dnaseq(seq, codon):
    # ---------translate_dnaseq------------#
    l_dna = len(seq)
    if l_dna % 3 != 0:
        seq = seq[: -(l_dna % 3)]
        # print("your sequence lenght is not tripple" + \
        # "but no worries, I have trimmed well format")
    # coding_dna = Seq(seq, generic_dna)
    # protein = coding_dna.translate(table=codon)
    protein = Seq(seq).translate()
    if "*" in protein:
        return False
    else:
        return True


def connectMetapairReads(record,
                         min_overlap=10,
                         max_overlap=150,
                         overlap_mismatch=1,
                         standard_length=150):
    records = record.strip().split()  # [ForA10_666 ATCG II9I CCGA II-I]
    name = records[0]  # eg ForA10_666
    seq1, seq2 = records[1], comp_rev(records[3])  # 正序列和反向互补序列，使用comp_rev()函数计算反向互补序列，因为反向结果中的reads与正向中的反向互补
    forward_qual, reverse_qual = records[2], records[4][::-1]  # 反向互补序列的质量分数也反转

    # 检查 seq1 和 forward_qual，seq2 和 reverse_qual 的长度是否一致，如果不一致则跳过
    if len(seq1) != len(forward_qual) or len(seq2) != len(reverse_qual):
        print(f"Skipping sequence {name} due to length mismatch between sequence and quality.")
        return None  # 直接跳过该条序列的处理

    overlaps = {}  # 储存不同重叠区域的错配情况
    for s in range(min_overlap, max_overlap + 1):  # 遍历从最小重叠到最大重叠长度的范围，寻找最佳重叠区域[10,150]
        l0 = seq1[-s:]  # 获取正向序列末尾的s个碱基
        l1 = seq2[0:s]  # 获取反向互补序列开头的s个碱基
        tmp_mismatch = mismatch(l0, l1)  # 比较两个重跌区域，计算错配数量
        if tmp_mismatch == 0:  # 如果错配数为0，表示完全匹配
            overlaps[s] = 1  # 完全匹配，赋值1
            # 找到最佳结果，直接跳出循环  #
            break

        elif tmp_mismatch <= overlap_mismatch:  # 如果错配数小于等于允许的错配数，则计算相似度
            tmp_identity = 1 - (tmp_mismatch / standard_length)  # 计算相似度
            overlaps[s] = tmp_identity  # 将相似度记录下来

    # 如果重叠结果为空，跳过该序列
    if not overlaps:
        return None

    # find best overlaping result in all potenial positions
    candidates = sorted(overlaps, key=overlaps.get, reverse=True)
    potenial = candidates[0]  # 选择相似度最高的重叠长度，potential为一个值，10
    s0 = seq1[-potenial:]  # 获取正向序列重叠部分, ATCGATCGAT
    s1 = seq2[0:potenial]  # 获取反向序列重叠部分, ATCGATCGAT
    corrected = ""  # 初始化修正后的序列
    for p in range(len(s0)):  # 遍历重叠区域中的每一个碱基，进行碱基的修正和比较(0,140)
        # site is changed, be careful!#
        tmp_loca0 = len(seq1) - potenial + p  # 标准长情况下：150(不变) - 10(改变) + 0(改变)

        if s0[p] == s1[p]:  # 如果正向碱基和反向碱基相等，直接添加到修正序列中
            corrected += s0[p]  # 添加一致的碱基
        else:  # 如果不相等，根据质量分值选择质量更高的碱基
            if 0 <= tmp_loca0 < len(forward_qual) and 0 <= p < len(reverse_qual):
                for_quality = forward_qual[tmp_loca0]
                rev_quality = reverse_qual[p]
                # print(f"for_quality:{for_quality}~~~~~rev_quality:{rev_quality}")

                if for_quality >= rev_quality:  # 选择质量分值较高的碱基
                    corrected += s0[p]
                else:
                    corrected += s1[p]
            else:
                print(f"[Connect error] Skipping position {p} due to invalid index in quality arrays.")
                print(potenial)
            # 捕获任何异常并记录相关信息
                print(f"Error processing sequence {name} at position {p}")
                print(f"Forward sequence: {s0}")
                print(f"Reverse sequence: {s1}")
                print(f"Forward quality: {forward_qual}")
                print(f"Reverse quality: {reverse_qual}")

# 组合修正后的共识序列
        makeup_consensus = (
                seq1[: standard_length - potenial]  # 正向序列非重叠部分
                + corrected  # 重叠部分修正后的序列
                + seq2[potenial - standard_length:]  # 反向序列非重叠部分
        )
    ## 限制共识序列长度范围在130~270bp之间
    if makeup_consensus and len(makeup_consensus) >= 130 and len(makeup_consensus) <= 270:
        return makeup_consensus  # 返回修正后的共识序列
    else:
        pass

def check_and_open_outhandle(file):
    if os.path.exists(file):
        print("[WARNING]: " + file + " exists! now overwriting")
    else:
        print("[INFO]: " + "open file " + file + "...")
    out = open(file, 'w')
    return out

def write_runlog(msg):
    with open(runlog, 'a') as log: log.write(msg + "\n")

write_runlog(f"vsearch: {args.vsearch}")
write_runlog(f"min_overlap: {args.min_overlap}, max_overlap: {args.max_overlap}")

# -------------------buildend----------------------------#
print_time("[INFO]: Building ends starts:")

# validate overlaps
if args.min_overlap < 10 or args.min_overlap > 120:
    print(f"[ERROR]: min_overlap {args.min_overlap} not in 10-120 range")
    exit(1)
if args.max_overlap < args.min_overlap:
    print("[ERROR]: max_overlap must be >= min_overlap")
    exit(1)

# prepare outputs
buildends_result = os.path.join(buildend_outdir, "buildends.fasta")
fh_out = open(buildends_result, 'w')

# log buildend settings
write_runlog(f"reads_limitation: {args.seqs_lim if args.seqs_lim else 'none'}")
write_runlog(f"input_length: {args.standard_length}")
write_runlog(f"vsearch_path: {args.vsearch}")
write_runlog(f"check_translation: {'yes' if args.reads_check else 'no'}")
write_runlog(f"cluster_identity: {args.cluster_identity}")
write_runlog(f"mismatch_allow: {args.buildend_mismatch}")
write_runlog(f"min_overlap: {args.min_overlap}")
write_runlog(f"max_overlap: {args.max_overlap}")

# --------------main-----------------------#

total_pairs, connected_pairs = 0, 0
with open(args.list) as fh_list:
    lines = [l.strip() for l in fh_list if l.strip()]

for line in lines:
    name = os.path.basename(line).split('.')[0]
    write_runlog(f"// processing {name}")
    if not os.path.exists(line) or os.path.getsize(line) == 0:
        write_runlog(f"! file missing or empty: {line}")
        continue
    success_connected = []
    with open(line) as fh:
        for r in fh:
            total_pairs += 1
            parts = r.split()
            if len(parts) != 5:
                continue
            cons = connectMetapairReads(r, min_overlap=args.min_overlap, max_overlap=args.max_overlap,
                                        overlap_mismatch=args.buildend_mismatch, standard_length=args.standard_length)
            if cons:
                success_connected.append(cons)
                connected_pairs += 1
    if success_connected:
        pid = os.getpid()
        temp_fa = f"temp.fa.{pid}"
        temp_uc = f"temp.uc.{pid}"
        with open(temp_fa,'w') as tf:
            for i,sq in enumerate(success_connected): tf.write(f">{i}\n{sq}\n")
        cmd = f"{args.vsearch} --cluster_smallmem {temp_fa} -usersort --threads {args.threads} --quiet --uc {temp_uc} --id {args.cluster_identity}"
        subprocess.call(cmd, shell=True)
        count = {}
        with open(temp_uc) as uc:
            for L in uc:
                if not L.startswith('H'): continue
                a=L.split()
                count.setdefault(a[9],0)
                count[a[9]]+=1
        sorted_clusters = sorted(count, key=lambda k:(count[k],k), reverse=True)
        if args.cluster_number_needKeep:
            sorted_clusters = sorted_clusters[:args.cluster_number_needKeep]
        elif args.abundance_threshod:
            sorted_clusters = [k for k in sorted_clusters if count[k]>=args.abundance_threshod]
        else:
            sorted_clusters = sorted_clusters[:5]
        for idx,k in enumerate(sorted_clusters,1):
            fh_out.write(f">{name}_{idx}_{count[k]}\n{success_connected[int(k)]}\n")
        write_runlog(f"{name}\tclusters_abundance: {','.join(str(count[k]) for k in sorted_clusters)}")
        os.remove(temp_fa); os.remove(temp_uc)
    else:
        write_runlog(f"[Vsearch WARNING] no connected pairs for {name}")

fh_out.close()

end_time = time.time()
write_runlog(f"END: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}")
write_runlog(f"DURATION: {end_time - start_time:.2f}s")
run_time("Building")