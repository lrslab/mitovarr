# mitovarr

[English](README.md) | 简体中文

原生 Rust 线粒体组装与注释软件。Cargo 包名、Rust 库名和命令名统一为 **`mitovarr`**，提供 CaenDR 远程提取、baiting、**BAM-baited de novo assembly**、基于已有 BAM 的 mapping assembly、基因注释、tRNA-Pro 旋转及 NCR 解析。

提取、gzip 读写、图组装、consensus 和注释都在进程内完成，不调用 SPAdes、外部 mapper、samtools 或外部注释程序。Mapping 使用已有 BAM，不负责把原始 FASTQ 全局比对到参考。

[完整流程图](docs/figure_1_workflow/README.md) · [文档目录](docs/README.md) · [实际 CLI 参数](docs/cli.md)

## 构建与安装

要求 Rust 1.90 或以上，使用 edition 2021。TLS 依赖 `ring` 的源码构建需要 C 编译工具链。

```bash
cargo build --release --locked
./target/release/mitovarr --help

cargo install --path . --locked
mitovarr --version
```

Python 和 Poppler 仅用于可选的文档绘图，不是 Rust 程序的构建或运行依赖。

## 输入范围

| 命令 | 输入与含义 |
| --- | --- |
| `fetch-caendr` | 物种和 strain；查询远程已索引线粒体 BAM 区域，恢复可获得的 mates，不下载整个 BAM |
| `assemble-baited` | BAM 提取后的双端 FASTQ；直接组装，验证范围也是提取的 pairs，不进行第二轮 baiting |
| `assemble-denovo` | 原始完整双端 FASTQ；reference-blind 图组装和验证 |
| `assemble-seeded` | 原始 pairs + 显式 seed FASTA；seed 只招募完整 pairs，验证使用全部原始 pairs |
| `bait` | 只进行迭代招募，输出 reads 与收敛状态 |
| `polish` / `assemble-mapping` | 坐标匹配的参考 FASTA 和已排序、已索引 BAM；可加入匹配 FASTQ 做局部结构重建 |
| `annotate` | 已冻结的候选/共识 FASTA 与物种或基因模型；独立注释、旋转和 NCR 定位 |

`assemble-mapping` 是 `polish` 的别名。`cel`、`cbr`、`ctro` 分别指 C. elegans、C. briggsae、C. tropicalis。FASTA/FASTQ 可为普通文本或 gzip；BAM 需要 `.bam.bai` 或 `.bai`，暂不支持 CRAM。输出目录必须不存在或为空。

## CaenDR 提取与 de novo

```bash
mitovarr fetch-caendr --species cbr --strain AF16 --out runs/AF16_reads

mitovarr assemble-baited \
  --r1 runs/AF16_reads/reads.R1.fastq.gz \
  --r2 runs/AF16_reads/reads.R2.fastq.gz --out runs/AF16_denovo
```

提取输出 paired/singleton/unpaired FASTQ、`mt.bam`、索引和报告。HTTP Range 请求只读取索引指向的数据块。输入范围为 **`mt-bam-extracted`**，不是原始全样本 FASTQ；上游比对仍影响哪些 reads 能被提取。无法定位或不存在的 mates、硬剪切丢失序列不能恢复。

包内 [`scripts/fetch_caendr.sh`](scripts/fetch_caendr.sh) 通过 Cargo 构建并启动同一个 Rust 提取器，原样传递参数，可从其他工作目录调用。详见[CaenDR 接口](docs/caendr.md)。

原始 FASTQ 的另外两种入口：

```bash
mitovarr assemble-denovo --r1 original_R1.fastq.gz --r2 original_R2.fastq.gz --out runs/blind

mitovarr assemble-seeded --r1 original_R1.fastq.gz --r2 original_R2.fastq.gz \
  --seed mitochondrial_seed.fasta --out runs/seeded
```

Seed 不进入组装图；招募未完成不能输出最终组装。主图默认 `--k 63`、baseQ 20、k-mer 最少 3 个 pair-record 支持，支持 15-127 的奇数 k；大 k 不一定适合短 reads。仅当低于主 k 时，才用 k=47/31 搜索末端连接。详见[组装方法](docs/de_novo.md)和[双端路径评分](docs/paired_paths.md)。

## Mapping 与局部重建

参考需与 BAM 的 contig 名、长度、坐标一致；多条序列的参考需要显式 `--contig`。

```bash
mitovarr assemble-mapping \
  --ref reference.fasta --bam mt.bam --contig MtDNA \
  --emit-coverage --emit-variants --out runs/mapping
```

普通 SNP 使用质量加权的片段严格多数票，默认 `--min-snp-frac 0.5`，恰好平票保留 N。Indel 默认 `--min-indel-frac 0.7`。重叠 mates 只计一次，repeat-aware 规则始终开启，不用参考碱基填补缺乏支持的位置。`base_votes.tsv` 保留混合等位证据。旧 `--min-alt-frac` 同时设置两个阈值，不能与分别设置的两个参数并用。

加入同一批匹配双端 reads 进行局部 NCR 重建：

```bash
mitovarr assemble-mapping \
  --ref reference.fasta --bam mt.bam --contig MtDNA \
  --local-assembly --r1 reads.R1.fastq.gz --r2 reads.R2.fastq.gz \
  --emit-coverage --emit-variants --out runs/mapping_local
```

局部候选可以包含参考中没有的样本序列。比较整条 contig 时保留原同源拷贝，各候选使用同一批双端片段，替换的两端连接都需有支持。搜索受限或候选竞争未解决时保留不确定性；并非所有 NCR 都能恢复。详见[方法和坐标契约](docs/mapping_local_assembly.md)。

高深度 BAM 可显式增加 `--subsample --max-records 1000000`，默认 seed=42，可用 `--subsample-seed` 修改。按 QNAME 选取完整目标 BAM 片段记录组，不截取坐标排序文件的前一段。局部图和评分使用入选名称对应的完整 FASTQ pairs，并保存为 `subsampled.R1.fastq`、`subsampled.R2.fastq`。缺少完整 FASTQ pair 的 BAM 片段只能贡献 pileup。`report.json.subsampling` 记录数量和范围；coverage/variant 比例描述抽样后的证据。未启用抽样时，超过记录上限仍报错，不静默截断。

## 注释、tRNA-Pro 旋转与 NCR

```bash
mitovarr annotate \
  --fasta runs/AF16_denovo/primary.fasta --species cbr \
  --assembly-ambiguity runs/AF16_denovo/primary_ambiguity.tsv \
  --rotate-to-trna-pro --out runs/AF16_annotation
```

注释先对双链进行核酸同源比对，再检查 CDS 编码特征（遗传密码表 5）与 RNA 边界。每个内置物种面板含 12 CDS、22 tRNA、2 rRNA；检出、部分和完整分别报告。尚未进行 RNA 二级结构评价。

旋转要求输入按环状解释，并有完整且唯一的 tRNA-Pro。自动模式优先读取 `annotation_topology`，否则识别字面 `circular=true`；`--topology circular` 只是注释时的环状假设，不能证明闭环。旋转会同时变换序列、基因和适配的组装歧义坐标。

`--assembly-ambiguity` 接受与 FASTA 坐标一致的 de novo `primary_ambiguity.tsv`、mapping 局部组装 `structural_ambiguity.tsv`，或上一次注释输出的投影表。普通 pileup `ambiguity.tsv` 是参考坐标，不能直接传入。

NCR1-4 由样本中观察到的侧翼基因定义；短/零长度间隔、缺失侧翼和未闭合末端仍记录在 `ncr_loci.tsv`。不是每个物种都有四个长 NCR。注释不会反向选择组装路径或解决重复，跨样本 NCR SV 自动分类尚未实现。详见[注释方法](docs/annotation.md)、[中文完整度规则](docs/annotation.zh-CN.md)和[NCR 定义](docs/ncr.md)。

## 输出与完成状态

| 输出 | 含义 |
| --- | --- |
| `primary.fasta`、`primary_ambiguity.tsv` | de novo 主代表序列和局部不确定性，可能仍为 unresolved |
| `candidates.fasta`、各类 path FASTA | 替代候选，选择证据保留在 JSON 中 |
| `assembly_graph.gfa`、`path_graph.gfa` | 清理后的 unitig 图和实际路径枚举图 |
| `candidate.fasta` | mapping 候选，局部模式包含被接受的样本路径替换 |
| `reference_candidate.fasta`、`local_assembly.gfa` | 局部模式的原始 pileup 骨架和局部候选路径图 |
| `structural_variants.tsv`、`structural_ambiguity.tsv` | 局部结构选择与候选坐标中的歧义 |
| `consensus.fasta` | 仅在该工作流最终规则通过时输出；mapping 完成不等于验证闭环 |
| `report.json`、证据 TSV | 输入范围、阈值、覆盖、投票、选择与不确定性 |
| `annotated.fasta`、GFF3、基因/蛋白 FASTA、NCR 表及 FASTA | 独立注释的序列、特征与位点状态 |

De novo 最终输出要求唯一有支持的主环、匹配碱基覆盖和连接支持；其他无关小分量不自动否决主环。Mapping 要求无 N、无剩余标记歧义，局部模式还要求各窗口解决。完整基因或无 N 不能独立证明 NCR 结构。

退出码 **0** 表示按该命令规则完成；**2** 表示未解决组装/招募或提取没有完整 pairs；**1** 表示输入/运行错误。注释即使缺基因或未旋转也可返回 0。Clap 参数用法错误同样返回 2，需结合终端消息和报告区分。

## 开发

```bash
cargo fmt --check
cargo test --locked
cargo clippy --all-targets --all-features --locked -- -D warnings
```

原包名/命令名为 `mitovar-rust` / `mtasm-rust`，现统一使用 `mitovarr`。为兼容已有报告与基因模型，`mitovar-rust.*` JSON schema 标识保持不变；新 GFF3 的 source 字段改为 `mitovarr`。本地 checkout 目录名不影响包名。
