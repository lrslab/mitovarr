# CaenDR 远程 BAM → 本地 FASTQ

固定入口是包内的 [`scripts/fetch_caendr.sh`](../scripts/fetch_caendr.sh)，提取实现是 [`src/caendr.rs`](../src/caendr.rs)。两者都随 Cargo 源代码包发布。

脚本仅使用 POSIX shell 启动本包的 `mitovarr fetch-caendr`，通过脚本位置定位 `Cargo.toml`。Cargo 首次构建 release 程序，后续复用编译结果；`--locked` 使用随包保存的依赖锁文件。脚本需要 Cargo/Rust 工具链，所有命令参数原样传给 Rust 程序，退出状态也原样传回。

BAM 索引查询、HTTP Range 请求、mate 恢复、gzip FASTQ 和 BAM 输出均由 Rust 实现，不依赖 samtools、SPAdes、Python 或系统 gzip。Rust 程序本身不启动 shell 或其他外部程序。

## 直接运行

从软件包根目录运行，例如 C. tropicalis NIC58：

```bash
./scripts/fetch_caendr.sh \
  --species ctro --strain NIC58 --out runs/NIC58_reads
```

C. elegans N2 和 C. briggsae AF16：

```bash
./scripts/fetch_caendr.sh \
  --species cel --strain N2 --out runs/N2_reads

./scripts/fetch_caendr.sh \
  --species cbr --strain AF16 --out runs/AF16_reads
```

从其他工作目录也可以通过脚本路径启动。`--out` 的相对路径以调用者的工作目录为准；脚本不会改变输出路径的含义。输出目录需新建或为空。需要查看参数时运行：

```bash
./scripts/fetch_caendr.sh --help
```

已有编译好的程序时，可直接运行 `mitovarr fetch-caendr`，不需要这个启动脚本或 Cargo。程序接收的参数与脚本一致。

## 参数

| 参数 | 含义 |
| --- | --- |
| `--species` | `cel`、`ctro`、`cbr`；也接受 `c_elegans`、`c_tropicalis`、`c_briggsae` |
| `--strain` | CaenDR strain 名，如 `N2`、`NIC58`、`AF16` |
| `--out` | 新的或空的输出目录 |
| `--contig` | 可选；覆盖自动识别的线粒体 contig 名 |
| `--timeout-seconds` | 单次 HTTP 请求超时，默认 120 秒 |
| `--max-templates` | 提取过程中保留的模板数量上限，默认 1,000,000 |

例如延长网络超时：

```bash
./scripts/fetch_caendr.sh \
  --species ctro --strain NIC58 --out runs/NIC58_reads \
  --timeout-seconds 300
```

## 提取过程与输出

程序由物种和 strain 组成 CaenDR 的 BAM/BAI 地址，读取索引并查询线粒体 contig，再按已记录的 mate 坐标进行一次恢复。HTTP Range 只读取需要的数据块，不下载完整 BAM，也不进行第二轮迭代 baiting。

仅保留 primary alignments；提取时保留 duplicate/QC-fail 标记的 reads，不设 MAPQ 截断。反向比对的 SEQ/QUAL 恢复到原始 read 方向。无法定位或不存在的 mates、硬剪切丢失的序列不能恢复；缺失 QUAL 写为 Phred 0。

| 输出文件 | 内容 |
| --- | --- |
| `reads.R1.fastq.gz`、`reads.R2.fastq.gz` | 成对、同步的 reads，直接供 `assemble-baited` 使用 |
| `reads.singletons.fastq.gz` | 配对记录中仅恢复到一端的 reads |
| `reads.unpaired.fastq.gz` | 不属于有效 R1/R2 配对槽位的 reads |
| `mt.bam`、`mt.bam.bai` | 提取记录和可用于 mapping consensus 的索引 |
| `report.json`、`recovery.json` | 来源 URL、数量、输出位置和网络请求统计 |

报告不计算数据文件哈希，也不增加完整性审计。正常提取得到 paired reads 时退出码为 `0`；提取结果没有完整 pairs 时为 `2`；运行错误返回非零状态。

这些 FASTQ 的范围是 **`mt-bam-extracted`**，不是全样本原始 FASTQ。此次远程提取就是 baiting 步骤，接下来直接组装：

```bash
cargo run --release --locked --bin mitovarr -- assemble-baited \
  --r1 runs/NIC58_reads/reads.R1.fastq.gz \
  --r2 runs/NIC58_reads/reads.R2.fastq.gz --out runs/NIC58_denovo
```

本脚本没有写死样本或临时二进制路径，原始数据不随软件包发布。
