# HiFiBarcoding
Methods for Obtaining Full-length DNA Barcodes Using Illumina (250bp pair-150bp) Sequencing
Reference: DOI:10.21769/BioProtoc.1010640
## 文件夹：HIFIBarcode_20250501
96_pairs_primer.txt，index和primer组合文件，对应96孔的各个样本，分配（Assign20250501.py）的时候用。   
primer，primer文件，组装完成后去除tag和primer序列的时候（Polish20250501.py）用。   
barcode，组装软件，用于将两端序列（end.fasta）与中间序列（middle.fasta）组装成完整的COI序列。   
vserach，去重软件，用于聚类各个样本的两端序列（Buildend20250501.py），减少重复序列。   
seqkit，需要提前在系统安装，用于将中间序列进行排序和去重（Chain20250501.py）。   
## Assign20250501.py
### 一、输入（Input）  
脚本通过命令行参数接收以下输入项：

1. `-outpre <STR>`  
   - **说明**：所有输出文件的前缀。  
   - **示例**：`-outpre test` → 输出目录 `test_assign/`，以及文件 `test_assign/assign.list` 等

2. `-index <INT>`  
   - **说明**：条形码（barcode）长度，用于从引物序列中截取前端标签。  
   - **示例**：`-index 6` 表示标签长度为6 bp

3. `-primer <STR>`  
   - **说明**：标签引物列表文件，格式为“样本名 引物序列”（每行两列），正反向各一行。  
   - **示例**：  
     ```text
     RevA1   AAGCTAAACTTCAGGGTGACCAAAAAATCA
     ForA1   AAGCGGTCAACAAATCATAAAGATATTGG
     ......
     ```

4. `-tmis <INT>` (`--tag_mismatch`)  
   - **说明**：标签（barcode）允许的不匹配碱基数，默认 `1`。  
   - **示例**：`-tmis 1`

5. `-pmis <INT>` (`--primer_mismatch`)  
   - **说明**：引物区允许的不匹配碱基数，默认 `1`。  
   - **示例**：`-pmis 1`

6. `-fq1 <STR>`  
   - **说明**：清洗后的左端（forward）FASTQ 文件，可为 `.fastq` 或 `.fastq.gz`。  
   - **示例**：`-fq1 cleaned_left.fastq.gz`

7. `-fq2 <STR>`  
   - **说明**：清洗后的右端（reverse）FASTQ 文件，可为 `.fastq` 或 `.fastq.gz`。  
   - **示例**：`-fq2 cleaned_right.fastq.gz`

8. `-batch_size <INT>`  
   - **说明**：每个子进程一次处理的 reads 数量，默认 `10000`。  
   - **示例**：`-batch_size 10000`

9. `-num_threads <INT>`  
   - **说明**：并行处理使用的线程数（CPU 核心数），默认 `2`。  
   - **示例**：`-num_threads 2`

---

### 二、输出（Output）  
运行结束后，所有结果均保存在以 `-outpre` 参数命名的文件夹下（如 `test_assign/`），主要包括：

1. **`assign_err.fasta`**  
   无法正确匹配引物或标签的 reads 序列，FASTA 格式。

2. **`assign_middle.ssam`**  
   匹配到 COI 中间区域（既非明确分配也非错误）的 reads，SSAM 格式。

3. **`For<样本编号>.ssam`、`Rev<样本编号>.ssam`**  
   按样本分别存放的分配结果，包含 reads 名称、序列和质量值。

4. **`<outpre>_assign/assign.list`**  
   每个样本对应的 `For*.ssam` 与 `Rev*.ssam` 文件路径列表，便于后续批量处理。

5. **`assign_<时间戳>.log`**  
   脚本开始／结束时间、总运行时长、每个样本 exact/sub/indel 三类匹配计数，以及总体错误与中间 reads 数量汇总。

6. **（可选）标准输出汇总**  
   示例：`Assigned: 12345, Errors: 67, Middle: 89, Total: 124...`

---

### 三、主要功能（Functionality）  
1. **多线程并行处理**  
   利用 `concurrent.futures.ProcessPoolExecutor` 实现多进程，将读取的 FASTQ 中的 reads 按 `batch_size` 分批提交至子进程处理，显著提升速度。

2. **条形码＋引物双端 demultiplex**  
   从每条 read 的 5′/3′ 端提取指定长度的标签序列，再去除引物区，判断正向与反向是否匹配：  
   - **exact**：完全一致；  
   - **sub**：允许在标签区／引物区内出现 ≤ `tmis`/`pmis` 个碱基替换；  
   - **indel**：标签区出现 ≤1 碱基的插入或缺失。

3. **错误与中间区域分类**  
   - **assign**：单端或双端匹配成功，将 reads 写入对应样本的 `*.ssam`；  
   - **err**：匹配引物异常或不在预期位置，写入 `assign_err.fasta`；  
   - **middle**：无法分配至任一样本，也非明显错误，写入 `assign_middle.ssam`。

4. **统计与日志**  
   实时统计每类匹配（exact/sub/indel）数目，以及总 reads、assigned/err/middle 数量；  
   在日志文件中记录详细统计和脚本运行时长，便于质量控制与复现。

## Buildend20250501.py
### 一、输入（Input）  
脚本通过命令行参数接收以下输入项：

1. `-outpre <STR>`  
   - **说明**：所有输出文件的前缀（脚本会创建 `<outpre>_buildend/` 目录）。

2. `-index <INT>`  
   - **说明**：引物标签长度，用于解析输入 SSAM 文件中带标签的序列。

3. `-vsearch <STR>`  
   - **说明**：`vsearch` 可执行程序路径；若不指定，则在环境变量中自动查找。

4. `-threads <INT>`  
   - **说明**：`vsearch` 聚类时使用的线程数，默认 `2`。

5. `-cid <FLOAT>` (`--cluster_identity`)  
   - **说明**：`vsearch` 聚类的相似度阈值，默认 `0.98`。

6. `-min <INT>` (`--min_overlap`)  
   - **说明**：配对拼接时最小重叠长度（bp），默认 `10`。

7. `-max <INT>` (`--max_overlap`)  
   - **说明**：配对拼接时最大重叠长度（bp），默认 `150`。

8. `-tp <INT>` (`--cluster_number_needKeep`)  
   - **说明**：每个样本保留的聚类簇数量（手动指定），默认保留前 5 簇。

9. `-ab <INT>` (`--abundance_threshod`)  
   - **说明**：按丰度阈值筛选保留簇（手动指定），默认保留前 5 簇。

10. `-seqs_lim <INT>` (`--seqs_lim`)  
    - **说明**：每个 `.ssam` 文件最大处理序列数，默认 `0`（不限制）。

11. `-len <INT>` (`--standard_length`)  
    - **说明**：标准读长（bp），用于计算相似度和拼接，默认 `150`。

12. `-ds` (`--drop_short_read`)  
    - **说明**：开启后丢弃短于标准长度的 reads，默认不开启。

13. `-rc` (`--reads_check`)  
    - **说明**：开启后对拼接后的序列进行翻译检查（无终止子），默认不开启。

14. `-bom <INT>` (`--buildend_mismatch`)  
    - **说明**：拼接允许的最大错配碱基数，默认 `1`。

15. `-list <FILE>`  
    - **说明**：一个文本文件，列出所有待处理的 `.ssam` 文件路径，一行一个，必需参数。

---

### 二、输出（Output）  
运行结束后，脚本在 `<outpre>_buildend/` 目录下生成：

1. **`buildends.fasta`**  
   - 所有样本拼接并去冗余后、按簇和丰度排序的最终 consensus 序列，FASTA 格式，序列名格式为 `>样本名_簇编号_丰度`。

2. **`buildend_<时间戳>.log`**  
   - 记录脚本开始/结束时间、参数设置、每个样本处理进度、聚类丰度统计及脚本运行时长。

3. **（中间文件）**  
   - 临时生成的 `temp.fa.<pid>` 和 `temp.uc.<pid>` 用于调用 `vsearch cluster_smallmem`，聚类完成后自动删除。

4. **STDOUT 信息**  
   - 自动检测或报错提示 `vsearch` 路径、各步骤开始结束时间和任何警告/错误信息。

---

### 三、主要功能（Functionality）  
1. **自动检测 `vsearch`**  
   - 如果未手动指定，脚本会在系统环境中查找 `vsearch` 可执行文件，并在找不到时报错退出。

2. **读取 `.ssam` 列表**  
   - 从 `-list` 指定的文件中依次读取每个 `.ssam` 文件路径，对每个文件按行解析五列：`<name> <seq1> <qual1> <seq2> <qual2>`。

3. **序列拼接（connectMetapairReads）**  
   - 对每对正反向序列：  
     1. 计算 reverse-complement  
     2. 在 `[min_overlap, max_overlap]` 区间内寻找最佳重叠（≤`buildend_mismatch` 错配）  
     3. 基于质量值选出不一致位点的更高质量碱基  
     4. 拼接成长度限制在 130-270 bp 的 consensus 序列  
     5. （可选）翻译检查保证无内部终止子

4. **去冗余聚类**  
   - 将每个样本成功拼接的 consensus 序列写入临时 FASTA，调用 `vsearch --cluster_smallmem` 按 `cluster_identity` 聚类并生成 `.uc` 文件。

5. **簇丰度统计与筛选**  
   - 解析 `.uc` 文件统计每个簇中序列数（丰度）；  
   - 根据 `-tp` 或 `-ab` 参数决定保留哪些簇（默认保留前 5 簇）。

6. **输出最终序列**  
   - 将每个保留簇的代表序列写入 `buildends.fasta`，序列名中包含样本名、簇编号和簇内序列数。

7. **日志与性能统计**  
   - 全程记录参数、样本处理状态、每步耗时及总运行时长，方便复现和质量控制。

## Chain20250501.py
### 一、输入（Input）  
脚本通过命令行参数接收以下输入项（均由 `argparse` 定义）：

| 参数                 | 变量名             | 类型    | 必需   | 默认     | 说明                                                                                                                                              |
| -------------------- | ------------------ | ------- | ------ | -------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-outpre <STR>`      | `args.outpre`      | String  | 是     | —        | 所有输出文件的前缀，脚本会在当前目录下创建 `<outpre>_chain/` 目录存放结果。|
| `-mf <STR>`          | `args.middle_input`| String  | 是     | —        | 输入的 SSAM 中间文件路径，包含每条记录的五列：`<name> <seq1> <qual1> <seq2> <qual2>`。 |
| `-mi <INT>`          | `args.min_insertsize` | Int     | 否     | 170      | 拼接后 consensus 序列最小长度（bp），脚本要求 ≥150，否则报错退出。 |
| `-ma <INT>`          | `args.max_insertsize` | Int     | 否     | 270      | 拼接后 consensus 序列最大长度（bp）。 |
| `-bom <INT>`         | `args.chain_mismatch` | Int     | 否     | 1        | 重叠区允许最大错配碱基数。 |
| `-num_threads <INT>`  | `args.num_threads` | Int     | 是     | 2        | 并行处理时使用的线程数。 |
| `-batch_size <INT>`  | `args.batch_size`  | Int     | 否     | 10000    | 每个子任务（线程）处理的 reads 对数量。 |
| `-seqkit <STR>`      | `args.seqkit_path` | String  | 否     | —        | 指定外部 `seqkit` 可执行文件路径；若不提供，脚本会在 `$PATH` 中查找。 |

---

### 二、输出（Output）  
脚本运行完成后，所有结果均写入 `<outpre>_chain/` 子目录：

1. **`chain.fasta`**  
   - 拼接并过滤后、去冗余的 consensus 序列（FASTA 格式），每条记录由两行组成：  
     ```  
     >样本名_middle  
     ATCG…  
     ```  

2. **`chain_<YYYYMMDD_HHMMSS>.log`**  
   - 运行日志，包含：  
     - 脚本开始／结束时间  
     - 参数设置  
     - 处理的总 read 对数  
     - 去重后序列数  
     - 序列长度分布统计（Length vs Count）  
     - 总运行时长（秒）

3. **临时中间文件**  
   - `chain.sorted.fasta`：按长度排序的中间 FASTA，用于去重，完成后被删除。  
   - 由 `seqkit` 生成的临时索引文件（`*.fai` 等），自动清理。

4. **Stdout 输出**  
   - 每一步开始／结束时间提示  
   - 错误或警告信息（如 `seqkit` 未找到时报错退出）
---

### 三、主要功能（Functionality）  

1. **参数解析与环境检查**  
   - 利用 `argparse` 定义并读取所有命令行参数；  
   - 检查 `min_insertsize ≥ 150`，否则直接 `exit`。  
   - 自动检测或验证 `seqkit` 可执行路径，找不到则报错退出。

2. **多线程读取与拼接**  
   - `chain_main()`：  
     - 逐行读取 SSAM 文件（支持 `.gz`），将每  `batch_size` 条记录打包为任务提交给 `ProcessPoolExecutor`。  
     - 每条记录包含：`name, seq1, qual1, seq2, qual2`。  
     - 收集所有子进程返回的 `(name, consensus_seq)`，写入 `chain.fasta`。

3. **核心拼接算法**  
   - `connectMetapairReads(record, min_overlap, max_overlap, overlap_mismatch, standard_length)`：  
     1. 计算反向互补序列和质量值反转。  
     2. 在每个可能的重叠长度区间内计算错配数与相似度，选出最佳重叠长度。  
     3. 对重叠区逐位比较：相同直接保留，不同则根据质量值高低选碱基。  
     4. 拼接前后不重叠部分与修正后重叠区，得到 consensus 序列。  
     5. 仅保留长度在 `[min_insertsize, max_insertsize]` 范围内的结果。

4. **排序与去冗余**  
   - 调用 `seqkit sort --by-length` 对 `chain.fasta` 进行长度排序，输出 `chain.sorted.fasta`；  
   - 再调用 `seqkit rmdup` 去除完全相同序列，将结果覆盖回 `chain.fasta`；  
   - 删除中间排序文件。

5. **统计与日志**  
   - 使用 `seqkit fx2tab -n -l` 提取每条序列的长度；  
   - 统计总对数、唯一序列数和各长度出现频次，写入日志文件；  
   - 记录结束时间并计算总耗时（秒）。

## Gapfill20250501.py
### 一、输入（Input）  
脚本通过命令行参数接收以下输入项（由 `argparse` 定义）：

| 参数                 | 变量名             | 类型    | 必需   | 默认   | 说明                                                                                   |
| -------------------- | ------------------ | ------- | ------ | ------ | -------------------------------------------------------------------------------------- |
| `-outpre <STR>`      | `args.outpre`      | String  | 是     | —      | 输出前缀。脚本将创建 `<outpre>_gapfill/` 目录存放所有结果。                              |
| `-gmax <INT>`        | `args.gmax`        | Int     | 否     | 700    | 最大组装序列长度（bp）。                                                               |
| `-gmin <INT>`        | `args.gmin`        | Int     | 否     | 450    | 最小组装序列长度（bp）。                                                               |
| `-kl <INT>`          | `args.lowkmer`     | Int     | 否     | 55     | 低复杂度 k‑mer 丰度阈值。                                                              |
| `-i <INT>`           | `args.kmer_interval` | Int   | 否     | 10     | k‑mer 大小区间步长。                                                                   |
| `-k <INT>`           | `args.kmer`        | Int     | 否     | 127    | 构建 de Bruijn 图的 k‑mer 长度。                                                        |
| `-S <INT>`           | `args.samp_num`    | Int     | 否     | 1      | 每个拆分 ends 子文件包含的样本数。                                                      |
| `-t <INT>`           | `args.cpu`         | Int     | 否     | 4      | 使用的 CPU 线程数。                                                                    |
| `-ends <STR>`        | `args.gapfill_ends`| String  | 是     | —      | ends 序列 FASTA 文件（带 For/Rev 标签的双端拼接序列）。                                  |
| `-middle <STR>`      | `args.gapfill_mid` | String  | 是     | —      | middle 序列 FASTA 文件（无标签的中间拼接结果）。                                        |
| `-barcode <STR>`     | `args.barcode`     | String  | 否     | "barcode" | SOAPBarcode 可执行程序路径，默认需在系统 PATH 中。                                       |

---

### 二、输出（Output）  
运行结束后，脚本在 `<outpre>_gapfill/` 目录下生成：

1. **日志文件**  
   - `gapfill_<YYYYMMDD_HHMMSS>.log`：记录脚本开始／结束时间、参数信息、每步处理统计、警告和总耗时。

2. **拆分的 ends 子文件**  
   - 目录 `ends/` 下的 `ends.1`, `ends.2`, …：按样本分组并经过丰度过滤后的 ends 序列子集。

3. **middle.lis**  
   - 位于主输出目录，内容为 middle 序列文件路径，用于后续填充调用。

4. **Shell 脚本**  
   - 目录 `shell/` 下的 `barcodes.1.sh`, `barcodes.2.sh`, …：每个拆分 ends 对应的 SOAPBarcode 调用脚本。

5. **结果目录**  
   - 目录 `result/` 下的 `barcodes.1*`, `barcodes.2*`, …：SOAPBarcode 输出的填充结果文件。

6. **Stdout 信息**  
   - 终端输出各步骤进度提示和可能的错误/警告。

---

### 三、主要功能（Functionality）  

1. **参数解析与目录准备**  
   - 初始化 `argparse`，读取所有命令行参数。  
   - 在当前工作路径下创建 `<outpre>_gapfill/` 及其子目录 `ends/`、`shell/`、`result/`。  
   - 配置日志系统，将详细信息写入日志文件。

2. **FASTA 解析工具**  
   - 定义 `parse_fasta(fh)` 生成器，按两行一组读取 FASTA 条目，返回 `(header, sequence)`。

3. **ends 序列分组与丰度过滤**  
   - 函数 `format_ends`：  
     - 读取带标签的 ends FASTA，将序列按样本名（去除 “For”/“Rev”）归组，并记录每条序列丰度。  
     - 对每对 forward/reverse 组合，剔除丰度差异超 10× 的对。  
     - 按参数 `-S` 指定的样本数拆分输出子文件，记录保留对数并写入 `ends/`。

4. **middle 列表准备**  
   - 在主目录生成 `middle.lis`，格式为 SOAPBarcode 要求的列表文件。

5. **生成并执行填充脚本**  
   - 遍历所有拆分后的 ends 子文件，为每个文件生成一个可执行的 shell 脚本，调用 SOAPBarcode：  
     ```
     barcode -e <ends.file> -r middle.lis -o result/barcodes.<n> \
             -x gmax -n gmin -l lowkmer -v kmer_interval -k kmer -t cpu
     ```  
   - 赋予执行权限并同步运行，保证前一个完成后再启动下一个。

6. **日志记录与性能统计**  
   - 在各主要步骤记录信息（输入对数、输出对数、警告）。  
   - 脚本结束后计算并记录总耗时，方便性能评估与复现。

## Polish20250501.py
### 一、输入（Input）  
脚本通过命令行参数接收以下输入项（由 `argparse` 定义）：

| 参数               | 变量名             | 类型    | 必需 | 默认    | 说明                                                                                   |
| ------------------ | ------------------ | ------- | ---- | ------- | -------------------------------------------------------------------------------------- |
| `-d <DIR>`         | `args.contigDir`   | String  | 是   | —       | 包含两个子目录 `ends/` 和 `result/` 的 contig 输出目录（stage1输入）。                  |
| `-outpre <STR>`    | `args.outpre`      | String  | 是   | —       | 输出前缀，脚本会创建 `<outpre>_polish/` 目录以及最后汇总到 `<outpre>_results/`。        |
| `-primer <FILE>`   | `args.primer`      | String  | 是   | —       | 引物文件（FASTA 格式），用于 stage2 中的引物修剪。                                      |
| `-i <FILE>`        | `args.coi_input`   | String  | 否   | —       | 如果提供，则跳过 stage1，直接对该 raw FASTA 进行 stage2 polish。                        |
| `-cov <INT>`       | `args.min_coverage`| Int     | 否   | 5       | 最低覆盖度阈值，head 中两端覆盖度均需 ≥ 该值，否则丢弃。                                 |
| `-l <INT>`         | `args.coi_length`  | Int     | 否   | 650     | 最小 COI 片段长度（bp），低于此长度丢弃。                                               |
| `-codon <INT>`     | `args.codon_table` | Int     | 否   | 5       | 翻译检查所用线粒体遗传密码表编号。                                                      |
| `-cc`              | `args.coi_check`   | Flag    | 否   | True    | 是否开启翻译检查；添加 `-cc` 则关闭翻译检查。                                           |

---

### 二、输出（Output）  
1. **Polish 目录及日志**  
   - 在 `<outpre>_polish/` 下生成日志文件 `polish_<YYYYMMDD_HHMMSS>.log`，记录 INFO/WARNING 级别运行详情和每步统计。

2. **Stage1 原始条形码 FASTA**  
   - 文件 `all_barcodes_raw.fasta`（若未提供 `-i`），将 `ends/` ↔ `result/` 中每个 contig 重命名并合并。

3. **Stage2 清洗后 FASTA**  
   - 文件 `all.barcodes.polished.fasta`，保留经引物修剪、覆盖度/长度过滤、可选翻译检查及 ≥97% 去重后的高质量序列。

4. **最终汇总**  
   - 脚本结束后调用 `collect_results`，将 `_assign`、`_buildend`、`_chain`、`_gapfill`、`_polish` 五个阶段目录及其核心产物移动到 `<outpre>_results/`，并清理中间目录。

5. **终端输出**  
   - 打印完成提示，包括 raw FASTA、polished FASTA 和日志文件的路径；以及汇总目录位置或任何缺失警告／错误。

---

### 三、主要功能（Functionality）  

1. **Stage1：合并 Contig → raw FASTA**  
   - 读取 `<contigDir>/ends/*.ends` 中序号→样本名映射；  
   - 解析 `<contigDir>/result/barcodes.<suffix>.contig`，将 contig header 中的序号映射回原始样本名并写入 `all_barcodes_raw.fasta`；  
   - 缺失映射或意外格式均产 Warning 日志。

2. **Stage2：Polish QC 流水线**  
   - **引物修剪**：在正/反链中分别寻找 forward/reverse 引物（允许 ≤1 处错配），截取两引物之间序列；  
   - **覆盖度过滤**：header 中提取 “x-y” 形式的覆盖度，双端均需 ≥ `min_coverage`；  
   - **长度过滤**：保留序列长度 ≥ `coi_length`；  
   - **翻译检查**（可选）：对正链和互补链 3 个阅读框翻译，若任一无终止子即通过，否则丢弃；  
   - **丰度 & 去冗余**：按平均覆盖度排序，保留彼此间序列身份 <97% 的代表序列，其余计入丢弃。

3. **Logging & Reporting**  
   - 针对每个样本统计原始条数、丢弃原因（primers/cov/len/trans/dedup）及最终保留数，并输出 INFO；  
   - Warning 细节记录每次修剪失败、映射缺失或翻译失败的具体样本/序列 ID。

4. **结果汇总**  
   - `collect_results` 按既定阶段顺序检查并移动文件，遇缺失立即报错或报警并退出/跳过；  
   - 最终将五个阶段的产物整合到 `<outpre>_results/`，方便下游分析和存档。  

## 运行命令
### 单命令
Assign20250501.py -fq1 test_1.fq.gz -fq2 test_2.fq.gz -num_threads 40 -index 6 -primer 96_pairs_primer.txt -outpre test     
Buildend20250501.py -outpre test -index 6 -threads 40 -list test_assign/assign.list     
Chain20250501.py -outpre test -num_threads 40 -mf test_assign/assign_middle.ssam     
Gapfill20250501.py -outpre test -ends test_buildend/buildends.fasta -middle test_chain/chain.fasta -t 40     
Polish20250501.py -d test_gapfill  -outpre test -primer primer     
### sh脚本（hifibarcode_py.sh）
hifibarcode_py.sh -fq1 test_1.fq.gz -fq2 test_2.fq.gz -n 40 -i 6 -o test

## 注意事项    
最终结果存储在_polish目录下的all.barcodes.polished.fasta，同时保留了运行过程中的所有log文件，一个脚本会生成2个log文件，其中一个由sh脚本调用命令生成。
![image](https://github.com/user-attachments/assets/6eb8fc60-026b-4cc5-af12-230485c13af0)
该脚本运行较慢，快速组装建议使用perl脚本：https://github.com/liushanlin/metabarcoding    
