# mitovarr documentation

These guides describe the current code and the `mitovarr` executable.

| Topic | Guide |
| --- | --- |
| Installation, inputs and examples | [English README](../README.md), [中文 README](../README.zh-CN.md) |
| Exact commands and defaults | [CLI help exported from the executable](cli.md) |
| CaenDR BAM extraction and mate recovery | [CaenDR interface](caendr.md) |
| Sample graph, cleanup and final-output rules | [Native de novo assembly](de_novo.md) |
| Gapped alignment, paired scoring and repeat paths | [Paired-path method](paired_paths.md) |
| Pileup, local NCR reconstruction and subsampling | [Mapping local assembly](mapping_local_assembly.md) |
| Homology, CDS/RNA completeness and tRNA-Pro rotation | [Annotation](annotation.md), [中文规则](annotation.zh-CN.md) |
| NCR1-4 flanks, short gaps and assembly uncertainty | [NCR definitions](ncr.md) |
| Complete workflow diagram | [Figure 1: SVG, PDF, PNG and legend](figure_1_workflow/README.md) |

The package, library and binary are named `mitovarr`. Existing
`mitovar-rust.*` JSON schema strings remain format identifiers for compatibility
with saved reports and gene-model panels; they do not name the executable.
GFF3 generated after the rename uses `mitovarr` in the source column.

Assembly candidates, verified final outputs and post-assembly annotation are
different outputs. Use their status and evidence fields together: complete gene
models, an N-free sequence, or an annotation circularity assumption does not
independently establish a resolved NCR structure.
