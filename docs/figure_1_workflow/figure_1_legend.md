# Detailed method schematic legend

## Manuscript legend

**Detailed method schematic | Native mitochondrial assembly and NCR interpretation with mitovarr.**
**a,** BAM-baited de novo assembly uses paired reads recovered from the indexed CaenDR mitochondrial BAM region; validation is restricted to those extracted pairs. Seeded assembly uses iteratively recruited pairs and validates against all original pairs; the seed never enters graph assembly. **b,** Native sample graphs are simplified and compacted, then read threading, smaller-k connections and bounded repeat-copy hypotheses generate candidate paths. Gapped alignments and empirical mate geometry compare alternatives using a common fragment cohort. Representative sequences and alternatives remain available; final de novo output additionally requires resolved structure, coverage and join support. **c,** Existing alignments support fragment-aware pileup consensus. Optional local assembly reconstructs sample-derived sequence between unique flanks and compares whole-contig alternatives while retaining homologous copies. Supported replacements are applied; unresolved loci remain explicit. Mapping completion does not independently establish circularity. **d,** Frozen sequences are annotated using species homology models and coding checks. Optional tRNA-Pro normalization requires a circular-input interpretation and a unique complete Pro model. Gene flanks define positional NCR1-4 intervals, with compatible assembly uncertainty retained in the output coordinates. Gene completeness and annotation topology do not validate repeat structure or physical circularity. Diagrams are schematic; R denotes a repeated segment. Graphs and annotations are reusable outputs, not automatic cross-sample NCR structural-variant calls.

## 中文说明

**详细方法示意图 | mitovarr 的原生线粒体组装与非编码区（NCR）解析流程。**

**a，输入与证据范围。** BAM-baited de novo assembly 从 CaenDR 已索引的线粒体 BAM 区域提取 reads 并恢复可获得的配对端，验证范围限于提取的片段。Seeded 模式以显式 mtDNA seed 迭代招募完整 read pairs，组装使用招募后的样本 reads，验证使用全部原始 pairs；seed 序列不进入组装图。Mapping 路线使用坐标匹配的参考 FASTA 和已排序、已索引 BAM。

**b，样本图 de novo 组装。** 在 Rust 内构建、清理并压缩样本图，通过 read threading、较小 k 的末端连接和有搜索边界的重复拷贝假设生成候选。对同一批片段进行允许错配和 indel 的比对，结合可用的文库方向及片段跨度信息比较路径，每个片段只贡献一次证据。主代表序列与替代路径可以在未解决状态下输出；最终序列还要求结构解析、覆盖和连接证据满足规则。图上成环本身不等于可靠闭环。

**c，参考引导的 mapping assembly。** 从现有 BAM 累积以片段为单位的 pileup，进行质量加权的 SNP/indel 共识调用；普通模式直接输出 pileup 结果。启用匹配双端 FASTQ 的局部组装后，程序检测 soft clips、indel 和歧义窗口，在样本图中搜索独特侧翼间的替代序列，并在保留原同源拷贝的整条 contig 上比较候选及两个连接点。NCR 新序列来自样本路径；仅接受有证据支持的替换，其他情况保留不确定性。Mapping 完成状态不代表独立验证了环化。

**d，组装后的独立注释。** 组装序列与选择状态冻结后，才引入物种同源模型并检查 CDS 编码特征与 RNA 边界。对按环状解释的输入，若存在完整且唯一的 tRNA-Pro，则可统一方向并将其 5′ 端旋转到位置 1。NCR1-4 从样本中识别到的侧翼基因推导，并保留适配后的组装歧义坐标。图中的 12 CDS、22 tRNA、2 rRNA 是内置模型面板的组成，并非每个样本都达到这些完整度。完整基因与注释拓扑不能证明重复结构或物理环化正确。

所有 reads、图结构、候选和基因位置均为原理示意，不代表某个样本的测量结果或真实比例。R 表示重复片段。FASTA、GFA、GFF3、TSV 和 JSON 保留序列、图、注释与证据，供后续 NCR/SV 分析使用；当前并未自动完成跨样本 NCR SV 分类。
