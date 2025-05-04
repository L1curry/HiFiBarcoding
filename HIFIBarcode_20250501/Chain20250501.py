'''
无引物reads的双端拼接
1) 拼接:首先,根据 reads 的双端信息把属于同一条片段的两端序列挑选出  来;然后,将每一对 PE reads 按照以下条件拼接起来:
a) 序列之间具有高  于 95%的相似性重叠区;
b) 重叠区域达到一定的重叠长度:前两类有引物  信息的 reads 拼接长度需要在 130-270 bp 间,第三类无引物信息的 reads  拼接后的序列长度限制在 170-270 bp 之间。
2) 深度过滤:统计每种标签对应的拼接序列的各个碱基的深度,将覆盖深度不足5的碱基删除;从而获得可以代表标签对应的昆虫条形码前后两端的一对序列。
支持多线程处理，大大缩短运行时间，但结果文件中的序列比非多线程多了不少，不知道为什么
'''
import os
import time
import gzip
import argparse
from concurrent.futures import ProcessPoolExecutor
import subprocess, shutil
from collections import Counter
import sys

parser = argparse.ArgumentParser(
    add_help=True,  # 禁用默认帮助信息
    description="无引物reads的双端拼接\n"
                "1) 拼接:首先,根据 reads 的双端信息把属于同一条片段的两端序列挑选出来;然后,将每一对PE reads按照以下条件拼接起来:\n"
                "a) 序列之间具有高于95%的相似性重叠区;\n"
                "b) 重叠区域达到一定的重叠长度:前两类有引物信息的reads拼接长度需要在130-270bp间,第三类无引物信息的reads拼接后的序列长度限制在170-270bp之间。\n"
                "2) 深度过滤:统计每种标签对应的拼接序列的各个碱基的深度,将覆盖深度不足5的碱基删除;从而获得可以代表标签对应的昆虫条形码前后两端的一对序列。\n"
                "\n"
                "Input:\n"
                "-outpre parameter (Required): Specifies the prefix for output files.\n"
                "-mf parameter (Required): Specifies the intermediate input file (in SSAM format) containing sequences to be processed.\n"
                "-mi parameter (Optional): Specifies the minimum length of linked fragments, with a default value of 170.\n"
                "-ma parameter (Optional): Specifies the maximum length of linked fragments, with a default value of 270.\n"
                "-bom parameter (Optional): Specifies the allowed number of mismatches in overlap regions, with a default value of 1.\n"
                "-num_threads parameter (Required): Specifies the number of CPU threads to use, with a default value of 2.\n"
                "-batch_size parameter (Required): Specifies the number of reads processed per CPU, with a default value of 10000.\n"
                "\n"
                "Output:\n"
                "chain.fasta: Records reads that consensus sequences, eg: >1(assign_middle.ssam文件中的序列序号)_middle(代表序列中间部分的编号) ATCGATCGATCG···.\n"
                "\n"
                "Example:\n"
                "python Chain_threads.py -outpre test -index 6 -mf test_assign/assign_middle.ssam -num_threads 40 -batch_size 1000"
                ,
    formatter_class=argparse.RawTextHelpFormatter  # 使用 RawTextHelpFormatter
)

# 创建一个参数组，用于存放通用参数
group = parser.add_argument_group("common arguments")

# 添加参数 -outpre
group.add_argument(
    "-outpre",
    metavar="<STR>",
    required=True,
    help="The prefix for output files. All output files will be named with this prefix followed by a specific suffix.\n"
         "For example, if -outpre is set to 'output', the output files might be named as 'output_clustered.fasta', 'output_stats.txt', etc.",
)

# 添加参数 -mf
group.add_argument(
    "-mf",
    metavar="<STR>",
    type=str,
    dest="middle_input",
    required=True,
    help="The input file in ssam format, containing the middle sequences to be processed.",
)

# 添加参数 -mi
group.add_argument(
    "-mi",
    metavar="<INT>",
    type=int,
    dest="min_insertsize",
    default=170,
    help="The minimum length of the connected fragment. Default value is 170.",
)

# 添加参数 -ma
group.add_argument(
    "-ma",
    metavar="<INT>",
    type=int,
    dest="max_insertsize",
    default=270,
    help="The maximum length of the connected fragment. Default value is 270.",
)

# 添加参数 -bom
group.add_argument(
    "-bom",
    metavar="<INT>",
    type=int,
    dest="chain_mismatch",
    default=1,
    help="The number of mismatches allowed in the overlap region during the chain process. Default value is 1.",
)

# 添加参数 -num_threads
group.add_argument(
    "-num_threads",
    metavar="<INT>",
    required=True,
    default=2,
    help="The number of CPU threads to use for processing. Default value is 2.",
)

# 添加参数 -batch_size
group.add_argument(
    "-batch_size",
    metavar="<INT>",
    type=int,
    default=10000,
    help="The number of reads handled by one CPU. Default value is 10000.",
)

# 添加参数 -seqkit，用于指定 seqkit 可执行文件路径
group.add_argument(
    "-seqkit",
    metavar="<STR>",
    dest="seqkit_path",
    help="(Optional) Path to the seqkit executable. If not provided, will try to find 'seqkit' in $PATH.",
)

args = parser.parse_args() # 解析命令行的参数并传递给args

# -----------------------functions for common---------------------#
def open_input(file):
    if file.endswith("gz"):
        return gzip.open(file, 'rt')
    else:
        return open(file, 'r')

def checkDirMkdir(fold):
    outdir = os.path.abspath(fold) # 只会创建_assign文件夹(args.outpre + "_assign")，除非改为接受参数fold
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
    return outdir

def print_time(info):
    print(info + " " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

def mismatch(str1, str2):
    # ----count matched bases of two sequences----#
    mismatches = 0
    for base in range(len(str1)):
        if str1[base] != str2[base]:
            mismatches += 1
    return mismatches

def complementation(sequence):
    # make a sequence complement #
    # replace function of string is too low!
    sequence = sequence.upper()
    transtable = str.maketrans('ATCG', 'TAGC')
    sequence = sequence.translate(transtable)
    return sequence # [BUG20210914]


def comp_rev(sequence):
    # make a sequence complement and reversed #
    sequence = complementation(sequence)
    return sequence[::-1]

def connectMetapairReads(record,
                         min_overlap=10,
                         max_overlap=150,
                         overlap_mismatch=1,
                         standard_length=150):
    results = []
    for item in record:
        name = item[0]  # 序列名称 >1_middle
        seq1, seq2 = item[1], comp_rev(item[3])  # 正序列和反向互补序列，使用comp_rev()函数计算反向互补序列，因为反向结果中的reads与正向中的反向互补
        forward_qual, reverse_qual = item[2], item[4][::-1]  # 反向互补序列的质量分数也反转

    # 检查 seq1 和 forward_qual，seq2 和 reverse_qual 的长度是否一致，如果不一致则跳过
        if len(seq1) != len(forward_qual) or len(seq2) != len(reverse_qual):
            print(f"Skipping sequence {name} due to length mismatch between sequence and quality.")
            continue            # 直接跳过该条序列的处理

        overlaps = {}  # 储存不同重叠区域的错配情况
        for s in range(min_overlap, max_overlap + 1):  # 遍历从最小重叠到最大重叠长度的范围，寻找最佳重叠区域
            l0 = seq1[-s:]  # 获取正向序列末尾的s个碱基
            l1 = seq2[0:s]  # 获取反向互补序列开头的s个碱基
            tmp_mismatch = mismatch(l0, l1)  # 比较两个重跌区域，计算错配数量
            if tmp_mismatch == 0:  # 如果错配数为0，表示完全匹配
                overlaps[s] = 1  # 完全匹配，赋值1
                break

            elif tmp_mismatch <= overlap_mismatch:  # 如果错配数小于等于允许的错配数，则计算相似度
                tmp_identity = 1 - (tmp_mismatch / standard_length)  # 计算相似度
                overlaps[s] = tmp_identity  # 将相似度记录下来
        if not overlaps:# 如果重叠结果为空，跳过该序列
            print(f"Skipping sequence {name} due to no valid overlaps found.")
            continue
        candidates = sorted(
            overlaps, key=overlaps.__getitem__, reverse=True
        )  # 对所有重叠区域的相似度进行排序，选出最好的重叠位置
        if len(candidates) > 0:  # 如果找到了候选重叠区域
            potenial = candidates[0]  # 选择相似度最高的重叠长度，potential为一个值，即最佳重叠范围，如43，代表43个碱基，正反序列重叠最好
            s0 = seq1[-potenial:]  # 获取正向序列的重叠部分
            s1 = seq2[0:potenial]  # 获取反向序列的重叠部分
            corrected = ""  # 初始化修正后的序列
            for p in range(len(s0)):  # 遍历重叠区域中的每一个碱基，进行碱基的修正和比较
                # site is changed, be careful!#
                tmp_loca0 = len(s0) - potenial + p  # 计算正向序列中的位置

                if s0[p] == s1[p]:  # 如果正向碱基和反向碱基相等，直接添加到修正序列中
                    corrected += s0[p]  # 添加一致的碱基
                else:  # 如果不相等，根据质量分值选择质量更高的碱基
                    if 0 <= tmp_loca0 < len(forward_qual) and 0 <= p < len(reverse_qual):
                        for_quality = forward_qual[tmp_loca0]
                        rev_quality = reverse_qual[p]
                        if for_quality >= rev_quality:  # 选择质量分值较高的碱基
                            corrected += s0[p]
                        else:
                            corrected += s1[p]
        # 组合修正后的共识序列
        # 共识序列长度限制在170-270 bp之间
        makeup_consensus = (
                seq1[: standard_length - potenial]  # 正向序列非重叠部分
                + corrected  # 重叠部分修正后的序列
                + seq2[potenial - standard_length:]  # 反向序列非重叠部分
        )
        if makeup_consensus and args.min_insertsize <= len(makeup_consensus) <= args.max_insertsize:
            # 改为收集结果，返回给主进程统一写入
            results.append((name, makeup_consensus))
    return results


# -------------------chain ---------------------#
### 记录运行开始时间 ###
start_time = time.time()
start_str = time.strftime("%Y%m%d_%H%M%S", time.localtime(start_time))
print_time(f"[INFO]: Chain starts at {start_str}")

inputMiddle = args.middle_input
if args.min_insertsize < 150:
    print("[ERROR]: min_insertsize can not be less than 150 bp")
    exit()

chain_outdir = checkDirMkdir(args.outpre + "_chain")
# 日志文件名 chain+开始时间.log
chain_out = os.path.join(chain_outdir, "chain.fasta")
chain_log = os.path.join(chain_outdir, f"chain_{start_str}.log")

cout = open(chain_out, 'w')
clog = open(chain_log, 'w')
# 先写入开始时间
clog.write(f"Start_time\t{start_str}\n")

def chain_main(inputMiddle,num_threads = int(args.num_threads), batch_size = int(args.batch_size), chain_mismatch = int(args.chain_mismatch)):
    total_pairs = 0
    middle_sequence = []
    futures = []
    with open_input(inputMiddle) as fm:
        with ProcessPoolExecutor(max_workers=num_threads) as executor:
            for r in fm:
                total_pairs += 1
                tmp = r.strip().split()
                if len(tmp) != 5:
                    continue
                else:
                    middle_sequence.append(tmp)
                    if len(middle_sequence) >= batch_size:
                        futures.append(executor.submit(connectMetapairReads, middle_sequence, chain_mismatch))
                        middle_sequence = []
                        pass
                    pass
                pass
            if middle_sequence:
                futures.append(executor.submit(connectMetapairReads, middle_sequence, chain_mismatch))
    # 等待所有子进程完成，并收集结果
    for fut in futures:
        for name, seq in fut.result():
            cout.write(name + "\n" + seq + "\n")
    return total_pairs

if __name__ == '__main__':
    total_pairs = chain_main(inputMiddle)
    cout.close()

    # 调用 seqkit：检查可用性
    seqkit_exec = args.seqkit_path or shutil.which("seqkit")
    if not seqkit_exec:
        print("[ERROR]: 'seqkit' not found in PATH. Please install seqkit or specify with -seqkit.")
        exit(1)

    # 排序并去重，输出到同名文件
    sorted_fa = os.path.join(chain_outdir, "chain.sorted.fasta")
    dedup_fa = chain_out  # 直接覆盖 original chain.fasta
    try:
        subprocess.check_call([seqkit_exec, "sort", "--by-length", "--line-width", "0", "-o", sorted_fa, chain_out])
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] seqkit 调用失败：{e}")
        sys.exit(1)
    try:
        subprocess.check_call([seqkit_exec, "-s", "-i", "rmdup", "--line-width", "0", "-o", dedup_fa, sorted_fa])
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] seqkit 调用失败：{e}")
        sys.exit(1)

    os.remove(sorted_fa) # 去掉中间文件

    # 统计长度分布并写入日志（输出 name 和 length）
    out = subprocess.check_output(
        [seqkit_exec, "fx2tab", "-n", "-l", dedup_fa],
        text=True
    )
    # 过滤掉 header（如果有）
    lines = out.splitlines()
    data_lines = [l for l in lines if not l.startswith("name")]
    lengths = [int(line.split("\t")[1]) for line in data_lines]

    freq = Counter(lengths)
    clog.write(f"total_pairs\t{total_pairs}\n")
    clog.write(f"Unique_sequences\t{len(lengths)}\n")
    clog.write("Length\tCount\n")
    for length, count in sorted(freq.items()):
        clog.write(f"{length}\t{count}\n")

    # 记录结束时间和耗时
    end_time = time.time()
    end_str = time.strftime("%Y%m%d_%H%M%S", time.localtime(end_time))
    runtime = end_time - start_time
    clog.write(f"End_time\t{end_str}\n")
    clog.write(f"Runtime_sec\t{runtime:.2f}\n")
    clog.close()

    print_time(f"[INFO]: Chain done at {end_str}")
    print(f"[INFO]: Total run time: {runtime:.2f}s")