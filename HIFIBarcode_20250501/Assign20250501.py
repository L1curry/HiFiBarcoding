'''
拆分:由于 96 个样品分别对应一对已知的标签序列,因此可以结合测序序  列双末端所包含的标签序列和引物序列的信息,将所有 reads 序列进行拆分。
    可以分为三类:带有 96 种前端 (条形码 5'末端) 标签引物的 reads,带  有 96 种尾端 (条形码 3'末端) 标签引物的 reads 和没有引物序列的 reads。
!!!支持多线程处理，大大缩短运行时间，但不会覆盖生成文件，而是继续追加，因此必需手动删除先前的运行结果后才能再次运行!!!
'''

import os
import time
import gzip
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed # 导入多进程处理模块
###############################################################################
#####------------------------- parameters --------------------------------#####

parser = argparse.ArgumentParser(
    add_help=True,  # 确保帮助信息可用
    description="!!!支持多线程处理，大大缩短运行时间，但不会覆盖生成文件，而是继续追加，因此必需手动删除先前的运行结果后才能再次运行!!!\n"
                "Input:\n"
                "-outpre <STR>: Output file prefix. All output files will be named using this prefix. Example: -outpre output\n"
                "-index <INT>: Length of the primer index (barcode) sequence. Example: -index 6\n"
                "-primer <STR>: Primer list file in the following format: Rev001   AAGCTAAACTTCAGGGTGACCAAAAAATCA For001   AAGCGGTCAACAAATCATAAAGATATTGG. Example: -primer primer_list.txt\n"
                "-tmis <INT>: Number of allowed mismatches in the index (barcode) sequence. Default value is 1. Example: -tmis 1\n"
                "-pmis <INT>: Number of allowed mismatches in the primer sequence. Default value is 1. Example: -pmis 1\n"
                "-fq1 <STR>: Cleaned left-end sequencing file in FASTQ format. Example: -fq1 cleaned_left.fastq.gz\n"
                "-fq2 <STR>: Cleaned right-end sequencing file in FASTQ format. Example: -fq2 cleaned_right.fastq.gz\n"
                "-batch_size <INT>: Number of reads processed per CPU. Default value is 10000. Example: -batch_size 10000\n"
                "-num_threads <INT>: Number of CPU threads to use. Default value is 2. Example: -num_threads 2\n"
                "\n"
                "Output:\n"
                "assign_err.fasta: Records reads that could not be assigned.\n"
                "assign_middle.ssam: Records intermediate processing results.\n"
                "For<sample>.ssam and Rev<sample>.ssam: Assignment result files for each sample.\n"
                "<outpre>_assign/assign.list: Records the file paths of assignment results for each sample.\n"
                "\n"
                "Example:\n"
                "python Assign_threads.py -outpre test -index 6 -primer 96_pairs_primer2.txt -fq1 test1.fq.gz -fq2 test2.fq.gz -batch_size 1000 -num_threads 40",
                formatter_class=argparse.RawTextHelpFormatter  # 使用 RawTextHelpFormatter
)

# 创建一个参数组，用于存放参数
assign_group = parser.add_argument_group("common arguments")

# 添加参数 -outpre
assign_group.add_argument(
    "-outpre",
    metavar="<STR>",
    required=True,
    help="The prefix for output files. All output files will be named with this prefix followed by a specific suffix.\n"
         "For example, if -outpre is set to 'output', the output files might be named as 'output_assign.txt', 'output_stats.txt', etc.",
)

# 添加参数 -index
assign_group.add_argument(
    "-index",
    metavar="<INT>",
    type=int,
    required=True,
    help="The length of the tag sequence in the ends of primers. This value is used to determine the length of the tag sequence when demultiplexing reads.",
)

# 添加参数 -primer
assign_group.add_argument(
    "-primer",
    metavar="<STR>",
    required=True,
    help="The taged-primer list file, which should be in the following format:\n"
         "Rev001   AAGCTAAACTTCAGGGTGACCAAAAAATCA\n"
         "For001   AAGCGGTCAACAAATCATAAAGATATTGG\n"
         "...\n"
         "This format is necessary!",
)

# 添加参数 -tmis
assign_group.add_argument(
    "-tmis",
    metavar="<INT>",
    type=int,
    dest="tag_mismatch",
    default=1,
    help="The allowed number of mismatches in the tag sequence when demultiplexing reads. Default value is 1.",
)

# 添加参数 -pmis
assign_group.add_argument(
    "-pmis",
    metavar="<INT>",
    type=int,
    dest="primer_mismatch",
    default=1,
    help="The allowed number of mismatches in the primer sequence when demultiplexing reads. Default value is 1.",
)

# 添加参数 -fq1
assign_group.add_argument(
    "-fq1",
    metavar="<STR>",
    required=True,
    help="The cleaned fastq1 file, which is the first paired-end sequencing file.",
)

# 添加参数 -fq2
assign_group.add_argument(
    "-fq2",
    metavar="<STR>",
    required=True,
    help="The cleaned fastq2 file, which is the second paired-end sequencing file.",
)

# 添加参数 -batch_size
assign_group.add_argument(
    "-batch_size",
    metavar="<INT>",
    type=int,
    default=10000,
    help="The number of reads handled by one CPU. Default value is 10000.",
)

# 添加参数 -num_threads
assign_group.add_argument(
    "-num_threads",
    metavar="<INT>",
    type=int,
    default=2,
    help="The number of CPU threads used for processing. Default value is 2.",
)

args = parser.parse_args() # 解析命令行的参数并传递给args

# ----------------------functions for assigning-------------------#
def check_and_open_outhandle(file):
    if os.path.exists(file):
        print("[WARNING]: " + file + " exists! so need to sure the is zero kb on the starting before")
    else:
        print("[INFO]: " + "open file " + file + "...")
    out = open(file, 'a')
    return out

def open_input(file):
    if file.endswith("gz"):
        return gzip.open(file, 'rt')
    else:
        return open(file, 'r')

def checkDirMkdir(fold):
    '''
    创建文件储存输出文件
    :param fold: 输出文件夹前缀
    :return:输出文件夹
    '''
    outdir = os.path.abspath(fold)
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
    return outdir

def print_time(info):
    '''
    获脚本开始运行和运行结束的时间
    :param info:是一个字符串，用于说明当前的操作或状态
    :return:打印字符串和当前时间
    '''
    print(info + " " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

def run_time(info):
    '''
    用于计算并打印某段代码的总运行时间
    :param info:字符串，说明当前的操作
    :return:打印info字符串和从t到当前时间的运行总时长，单位是S
    '''
    print("[INFO]: {0} total run time: {1:.2f}".format(info,
                                                       time.time() - t) + "s")
def complementation(sequence):
    '''
    取互补序列
    :param sequence:
    :return:
    '''
    sequence = sequence.upper()
    transtable = str.maketrans('ATCG', 'TAGC')
    sequence = sequence.translate(transtable)
    return sequence


def comp_rev(sequence):
    # 基于上一个函数取反向互补序列 #
    sequence = complementation(sequence)
    return sequence[::-1]


def decode(head_flen, head_rlen, lookup):
    if head_flen in lookup:
        return lookup[head_flen], 'exact'
    if head_rlen in lookup:
        return lookup[head_rlen], 'exact'
    return False, None


def decode_mismatch(f, r, lookup):
    # 首先针对 f 进行 sub/indel 检测
    for seq, sample in lookup.items():
        # substitution 判断
        if len(f) == len(seq):
            tag_mis = primer_mis = 0
            for i, (b1, b2) in enumerate(zip(f, seq)):
                if b1 != b2:
                    if i < args.index:
                        tag_mis += 1
                    else:
                        primer_mis += 1
            if tag_mis <= args.tag_mismatch and primer_mis <= args.primer_mismatch:
                return sample, 'sub'
        # barcode 区 1-bp indel 判断
        if abs(len(f) - len(seq)) == 1:
            longer, shorter = (f, seq) if len(f) > len(seq) else (seq, f)
            for pos in range(args.index):
                if longer[:pos] + longer[pos + 1:] == shorter:
                    return sample, 'indel'
    # 同样流程对 r 再次检测
    for seq, sample in lookup.items():
        # substitution 判断
        if len(r) == len(seq):
            tag_mis = primer_mis = 0
            for i, (b1, b2) in enumerate(zip(r, seq)):
                if b1 != b2:
                    if i < args.index:
                        tag_mis += 1
                    else:
                        primer_mis += 1
            if tag_mis <= args.tag_mismatch and primer_mis <= args.primer_mismatch:
                return sample, 'sub'
        # barcode 区 1-bp indel 判断
        if abs(len(r) - len(seq)) == 1:
            longer, shorter = (r, seq) if len(r) > len(seq) else (seq, r)
            for pos in range(args.index):
                if longer[:pos] + longer[pos + 1:] == shorter:
                    return sample, 'indel'
    return False, None

def dis_index(index_list):
    """
    计算给定条形码列表中每对条形码之间的错配（碱基不同）的个数，返回最小和最大错配数。
    :param barcode_list: 包含多个条形码（字符串形式）的列表
    :return: 最小错配数和最大错配数的元组 (min_dis, max_dis)
    """
    # 初始化一个空列表，用于存储每对index之间的错配数
    dis = []

    # 当 barcode_list 中还有条形码时，继续循环
    while index_list:
        # 从列表中弹出一个index b1
        b1 = index_list.pop()

        # 与剩下的所有index逐一比较
        for b2 in index_list:
            mismatch = 0  # 初始化错配计数器

            # 对两个index的每个碱基逐位进行比较
            for base in range(len(b1)):
                if b1[base] != b2[base]:  # 如果当前位置的碱基不同
                    mismatch += 1  # 错配数加 1

            # 将 b1 和 b2 之间的错配数添加到 dis 列表中
            dis.append(mismatch)

    # 找到 dis 列表中的最小和最大错配数
    min_dis = min(dis)
    max_dis = max(dis)

    # 返回最小错配数和最大错配数的元组
    return (min_dis, max_dis)


# ----------------assign----------------------#
t = time.time()

fq1, fq2 = args.fq1, args.fq2 # 获取两个测序文件

assign_outdir = checkDirMkdir(args.outpre + "_assign")  # 创建输出目录保存分配结果
ErrFile = check_and_open_outhandle(assign_outdir + "/assign_err.fasta")  # 记录出错的reads
Middle_ssam = check_and_open_outhandle(assign_outdir + "/assign_middle.ssam") # 记录COI中间的测序文件

# —— 新增 log 文件 ——
start_time_str = time.strftime("%Y%m%d_%H%M%S", time.localtime(t))
log_path = os.path.join(assign_outdir, f"assign_{start_time_str}.log")
log_fh = open(log_path, "w")
log_fh.write(f"Start Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(t))}\n")

indexlen = args.index  # 读取条形码的索引长度

pris, indp, FH, indexs = {}, {}, {}, [] # {样本名:{For/Rev:引物序列}}；{引物名:引物序列}；{引物序列:引物名}；[For primer index]

with open(args.primer, "r") as p:  # 打开引物文件并解析其中的内容
    primer_lines = 0  # 记录引物文件的行数
    for i in p.readlines():
        primer_lines += 1  # 每读取一行，行数加1
        arr = i.strip().split()  # 去除每行的空格和换行符，并按空格切分
        if len(arr) != 2:
            print("primer set is not well-formated")  # 如果分割结果不是2项，提示引物格式错误
            exit()
        sam, ipr = arr[0], arr[1]  # 引物名, 引物序列
        FH[ipr] = sam  # {ATCG:ForH1, TGAC:RevH1, .....}，
        if sam in indp.keys():  # 检查引物名是否已经存在字典中，避免重复[!!!]
            print(sam + "show twice in primer set")
        else:
            ori = sam[0:3]  # 提取样本名前3个字符，判断是For还是Rev引物
            if ori == "For": # 提取引物名最后3个字符，作为样本名
                num = sam.replace("For","")
                pass
            else:
                num = sam.replace("Rev","")
                pass
            pris[num] = {} # 在字典 pris 中，将键 num（样本名） 对应的值设为一个新的空字典 {H1:{},H2:{}...}
            pris[num][ori] = ipr  # {H1:{Rev:ATCG, For:TAGC}, H2:{Rev:TCGA, For:AGCT}, H3...}
            indp[sam] = ipr  # 与FH字典相反，调换了键值对的位置。{ForH1:ATCG, RevH1:TGAC, .....}

        if ori == "For":
            plenf = len(ipr)  # 记录正向引物的长度
            primerF = ipr[indexlen:]  # 提取正向引物序列，去掉index（一般为前5或6个碱基）
            indexs.append(ipr[0:indexlen])  # 将index保存至列表中

        if ori == "Rev":
            plenr = len(ipr)  # 记录反向引物的长度
            primerR = ipr[indexlen:]  # 提取反向引物索引后的引物序列
# 检查引物文件是否有偶数行（每个样本应有For和Rev引物）
if primer_lines % 2 != 0:
    print("primer lines ({}) is not even number".format(primer_lines))
    print("the primer file need to have each forward and reverse primer")
    exit()

# 分析index之间的距离，确保index错配数合理
(min_dis, max_dis) = dis_index(indexs)
print("[INFO]: min distance among index is {}".format(min_dis))
print("[INFO]: max distance among index is {}".format(max_dis))
print("[INFO]: mismatches allowed in barcodes is {}".format(args.tag_mismatch))

# 检查用户设置的错配数是否合理，不能超过最小index距离减1，由于index序列固定，所以mismatche只能为1，不然都会引起报错
if args.tag_mismatch > (min_dis - 1):
    print("mismatch you set is too large to demultiplex, it must be smaller"
          + " than {},".format(min_dis) + " because min distance among indexs"
          + " is {}".format(min_dis)
          )
    exit()

# 初始化字典，记录分配到的reads数量
count_assigned = {}
for fh in FH.values(): # fh such ForH1
    count_assigned[fh] = 0  # 初始化每个样本的分配计数
count_assigned['middle'] = 0  # 中间区域分配reads数初始化为0

neg_priF = comp_rev(primerF) # 获取正向引物序列的反向互补序列
neg_priR = comp_rev(primerR) # 获取反向引物序列的反向互补序列
assigned_list = assign_outdir + "/assign.list"  # 创建分配列表，记录每个样本的分配文件路径

with open(assigned_list, "w") as ls: # 创建每个样本reads对应的储存文件
    sorted_sample = sorted(pris.keys())  # 按样本编号排序，{H1:{Rev:ATCG, For:TAGC}, H2:{Rev:TCGA, For:AGCT}, H3...}
    for s in sorted_sample:
        ls.write(assign_outdir
                 + "/For" + s + ".ssam" + "\n"
                 + assign_outdir
                 + "/Rev" + s + ".ssam" + "\n")

# 打开所有样本的分配文件，准备写入分配的reads
filehandle = {}
for sam in indp.keys(): #
    filehandle[sam] = check_and_open_outhandle(assign_outdir + "/" + sam + ".ssam")


err, assigned =  0, 0  # 错误计数；序列计数；已分配计数

# 打开过滤后的FASTQ文件进行读取
fh1 = open_input(args.fq1)
fh2 = open_input(args.fq2)
def process_read_pair(name, seq1, seq2, qual1, qual2, plenf, plenr, primerF, primerR, neg_priF, neg_priR, FH):
    """
    处理一对reads，检查引物并分配reads到相应的样本。
    """
    # 初始化正反向匹配结果
    fmatch = rmatch = fmatch_withMis = rmatch_withMis = False
    # 分别获取正反reads的头部序列
    head1f = seq1[:plenf]
    head1r = seq1[:plenr]
    head2f = seq2[:plenf]
    head2r = seq2[:plenr]

    # 去掉引物后的序列
    len_head_cut = max(plenf, plenr)
    tmp1 = seq1[len_head_cut:]
    tmp2 = seq2[len_head_cut:]

    # 检查引物是否存在于序列中，若存在即视为错误
    if any(primer in tmp1 for primer in [primerF, primerR, neg_priF, neg_priR]):
        return "err", (name, seq1, seq2)

    if any(primer in tmp2 for primer in [primerF, primerR, neg_priF, neg_priR]):
        return "err", (name, seq1, seq2)

    # 正反向匹配
    sample_f, kind_f = decode(head1f, head1r, FH)
    if not sample_f:
        sample_f, kind_f = decode_mismatch(head1f, head1r, FH)
    sample_r, kind_r = decode(head2f, head2r, FH)
    if not sample_r:
        sample_r, kind_r = decode_mismatch(head2f, head2r, FH)

    f_has = bool(sample_f)
    r_has = bool(sample_r)

    if f_has and not r_has:
        return "assign", (sample_f, sample_f[:3], "One",
                          name, seq1, qual1, seq2, qual2, kind_f)

    elif (rmatch or rmatch_withMis) and not (fmatch or fmatch_withMis):
        if r_has and not f_has:
            return "assign", (sample_r, sample_r[:3], "Two",
                              name, seq1, qual1, seq2, qual2, kind_r)
    elif (fmatch == False and fmatch_withMis == False and rmatch == False and rmatch_withMis == False):
        return "middle", (name, seq1, qual1, seq2, qual2)
    else:
        return "err", (name, seq1, seq2)


def process_read_batch(read_batch, plenf, plenr, primerF, primerR, neg_priF, neg_priR, FH):
    """
    每个进程处理一批 reads 的逻辑。
    read_batch: List of tuple (name, seq1, seq2, qual1, qual2)
    """
    batch_results = []
    for name, seq1, seq2, qual1, qual2 in read_batch:
        # 假设 process_read_pair 是你定义的用于处理每对 reads 的函数
        result_type, data = process_read_pair(name, seq1, seq2, qual1, qual2, plenf, plenr, primerF, primerR,
                                              neg_priF, neg_priR, FH)
        batch_results.append((result_type, data))
    return batch_results

def process_reads_in_parallel(fh1, fh2, plenf, plenr, primerF, primerR, neg_priF, neg_priR, FH, filehandle, ErrFile,
                              Middle_ssam, batch_size=int(args.batch_size), num_threads=int(args.num_threads)):
    """
    并行处理FASTQ文件中的配对reads，按每1000条为一批次处理。
    """
    # —— 新增统计容器 —— #
    sample_exact = {sam: 0 for sam in FH.values()}
    sample_sub = {sam: 0 for sam in FH.values()}
    sample_indel = {sam: 0 for sam in FH.values()}
    total_sub = total_indel = 0
    assigned = err = middle = 0
    count_assigned = {}
    read_batch = []
    futures = []
    seqnum = 1
    # 创建进程池
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        while True:
            # 逐条读取fastq记录，每次读取四行
            read1 = [fh1.readline().strip() for _ in range(4)]
            read2 = [fh2.readline().strip() for _ in range(4)]

            # 如果读取完毕，退出
            if not read1[0] or not read2[0]:
                break

            name, seq1, seq2, qual1, qual2 = read1[0], read1[1], read2[1], read1[3], read2[3]
            read_batch.append((name, seq1, seq2, qual1, qual2))

            # 如果当前批次已达到指定数量，则提交给进程池处理
            if len(read_batch) >= batch_size:
                futures.append(executor.submit(process_read_batch, read_batch, plenf, plenr, primerF, primerR, neg_priF, neg_priR, FH))
                read_batch = []  # 清空批次，准备下一批

        # 处理最后剩余的reads（不足一个批次的情况）
        if read_batch:
            futures.append(executor.submit(process_read_batch, read_batch, plenf, plenr, primerF, primerR, neg_priF, neg_priR, FH))

        # 处理所有的 futures
        for future in as_completed(futures):
            batch_results = future.result()
            for result_type, data in batch_results:
                if result_type == "err":
                    name, seq1, seq2 = data
                    ErrFile.write(f">{name}_priErr\n{seq1}\n{seq2}\n")
                    seqnum += 1
                    err += 1
                elif result_type == "assign":
                    target, dir, typ, name, seq1, qual1, seq2, qual2, kind = data
                    # —— 统计 exact/sub/indel ——
                    if kind == 'sub':
                        sample_sub[target] += 1;
                        total_sub += 1
                    elif kind == 'indel':
                        sample_indel[target] += 1;
                        total_indel += 1
                    else:
                        sample_exact[target] += 1
                    if typ == "One":
                        if dir == "For":
                            filehandle[target].write(f"{target}_{seqnum} {seq1} {qual1} {seq2} {qual2}\n")
                            seqnum +=1
                            pass
                        else:
                            filehandle[target].write(f"{target}_{seqnum} {seq2} {qual2} {seq1} {qual1}\n")
                            seqnum += 1
                            pass
                        pass
                    elif typ == "Two":
                        if dir == "For":
                            filehandle[target].write(f"{target}_{seqnum} {seq2} {qual2} {seq1} {qual1}\n")
                            seqnum += 1
                            pass
                        else:
                            filehandle[target].write(f"{target}_{seqnum} {seq1} {qual1} {seq2} {qual2}\n")
                            seqnum += 1
                            pass
                        pass
                    else:
                        print("未匹配类型")
                    assigned += 1
                    count_assigned[target] = count_assigned.get(target, 0) + 1

                elif result_type == "middle":
                    middle += 1
                    name, seq1, qual1, seq2, qual2 = data # 这里获取了reads name 但实际并没有使用
                    Middle_ssam.write(f">{seqnum}_middle {seq1} {qual1} {seq2} {qual2}\n")
                    seqnum += 1

    total_reads = err + assigned + middle
    return assigned, err, middle, total_reads, sample_exact, sample_sub, sample_indel, total_sub, total_indel

# 示例调用

if __name__ == '__main__':
    assigned, err, middle, total_reads, sample_exact, sample_sub, sample_indel, total_sub, total_indel = \
                process_reads_in_parallel(fh1, fh2, plenf, plenr, primerF, primerR, neg_priF, neg_priR,
                       FH, filehandle, ErrFile, Middle_ssam)
    summary = f"Assigned: {assigned}, Errors: {err}, Middle: {middle}, Total: {total_reads}"
    print(summary)
    # —— 写入 log 结束信息 ——
    end_t = time.time()
    log_fh.write(f"End   Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_t))}\n")
    log_fh.write(f"Duration: {end_t - t:.2f} s\n")
    log_fh.write("Per-sample counts (exact / sub / indel):\n")
    for s in sorted(sample_exact):
        log_fh.write(f"  {s}: {sample_exact[s]} / {sample_sub[s]} / {sample_indel[s]}\n")
    log_fh.write(f"Total sub mismatches: {total_sub}\n")
    log_fh.write(f"Total indels: {total_indel}\n")
    log_fh.write(f"Errors: {err}\nMiddle: {middle}\n")
    log_fh.close()