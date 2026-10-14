# Native NCR annotation

`annotate` derives NCR intervals from the genes detected on each sample contig,
after any tRNA-Pro rotation. It writes both the annotated sequence intervals and
four positional-locus status rows per contig. The implementation runs in Rust,
offline, with the same gene models used by [gene prediction](annotation.md).

```bash
mitovarr annotate --fasta candidate.fasta --species cbr \
  --topology circular --rotate-to-trna-pro --ncr-min-len 20 \
  --out annotation
```

Use `--topology auto` (the default) for de novo candidates. It reads a saved
`annotation_topology=circular|linear` when present, otherwise a literal
`circular=true` requests circular coordinates. Other inputs are linear. The
example's `circular` option is a supplied topology assumption, not a closure
test. Annotation saves that assumption for reuse and never changes assembly
selection or publication status.

## Positional names

The NCR1–4 nomenclature is based on neighboring genes, irrespective of sequence
homology, following [Li et al. 2018, DNA Research](https://pmc.ncbi.nlm.nih.gov/articles/PMC6289777/).
The order below follows gene transcription; the genomic coordinate order reverses
on the minus strand.

| Locus | C. elegans | C. briggsae | C. tropicalis |
| --- | --- | --- | --- |
| NCR1 | tRNA-Ala → tRNA-Pro | tRNA-Ala → tRNA-Pro | tRNA-Ala → tRNA-Pro |
| NCR2 | Gln → Phe position inspected; no baseline long NCR | tRNA-Gln → tRNA-Phe | **tRNA-Arg → tRNA-Gln** |
| NCR3 | ND4 → COX1 | ND4 → COX1 | ND4 → COX1 |
| NCR4 | ND3 → ND5 position inspected; no baseline long NCR | ND3 → ND5 | ND3 → ND5 position inspected; no baseline long NCR |

These species do not all have four qualifying NCRs. A newly observed long gap at
an unexpected species/locus combination receives `NCR2_candidate` or
`NCR4_candidate`; it is not suppressed or called a confirmed insertion.
`expected_in_species_reference` describes the catalog, not a requirement imposed
on a sample.

Short gaps and overlapping gene boundaries remain visible in the locus table. The default
`--ncr-min-len 20` filters emitted sequence intervals; it is configurable and
does not assign absent species labels or require a fixed NCR length.

## Boundaries and sequence handling

The union of observed gene intervals defines intergenic gaps. A partial homolog
covering less than 80% of a model does not mask an NCR when a copy covering at
least 80% is already detected elsewhere. Such fragments remain in the gene
outputs and are listed as `embedded_partial_homologs` when they intersect an NCR.
A partial sole gene copy remains in the mask. This rule does not prove that an
embedded fragment is a pseudogene or that its sequence originated from ND5.

A supported NCR boundary needs one detected flank on each side, consistent
strand/order, at least 80% model coverage, no competing gene identity or limited
model search, and support for the gene endpoint facing the gap. An internal RNA
length variant or a defect at the gene's opposite end does not by itself erase a
supported local boundary. Missing/uncertain facing endpoints produce a
`*_candidate` interval. For example, missing ND4 terminal sequence cannot become
a confidently annotated NCR3 insertion merely because it enlarges the gap.

For a circular hypothesis, a gap crossing position 1 is joined before applying
the minimum length. It has two coordinate segments but one sequence entry. For
a linear contig, head and tail remain independent fragments; a suitably facing
single flank can identify an `NCR*_candidate` terminal fragment. Their lengths
are never summed into a complete NCR length. No detected genes means no inferred
whole-contig NCR. Multiple bounded copies remain separate regions.

NCR sequences are extracted from sample bases, in forward genomic order. Unknown
bases remain unchanged and are counted. `at_fraction` uses only called A/C/G/T
bases and is null for an entirely ambiguous interval. `sequence_origin` is
`not_assessed` unless `--assembly-ambiguity` reports a mapping
`fixed_span_graph_fill` overlapping the interval; it is then
`fixed_span_graph_fill`, meaning the bases were filled from sample-graph walks at
the existing (reference-guided) span and the NCR length was not measured.
Positional annotation does not establish sequence ancestry, regulatory function,
repeat copy number or a read-supported structural variant.

## Files and status fields

| File | Contents |
| --- | --- |
| `annotations.gff3` | Gene features, NCR regions and optional projected assembly evidence, all in `annotated.fasta` coordinates |
| `ncr.gff3` | NCR/gap regions only; supported intervals use `intergenic_region`, uncertain or terminal intervals use `region` |
| `ncr.tsv` | One row per interval: locus, segments, length, genomic left/right genes, boundaries, ambiguity count, AT fraction and embedded partial homologs |
| `ncr_loci.tsv` | Four rows per species contig, including missing, too-short, overlapping and unresolved loci |
| `ncr.fasta` | One sequence per emitted interval; separate entries for separate linear terminal fragments |
| `report.json` | Per-record `ncr` object with catalog, threshold, regions, notes and locus states |

Coordinates are 1-based inclusive. NCR segments follow forward genomic order
through a possible circular origin, unlike gene segments, which follow biological
5′→3′ order. FASTA and GFF IDs agree and are unique across input records.
`features.tsv` continues to contain genes only. GFF includes unclassified gaps
above the threshold as `unclassified_gap`; those are not automatically extra NCRs.

`annotated` and `annotated_with_unknown_bases` indicate supported boundaries, with
unknown sequence explicitly distinguished. `boundary_uncertain` and
`positional_candidate` describe intervals requiring further interpretation.
`below_min_length`, `adjacent_without_gap` and `overlapping_flanks` report the
observed geometry. `unresolved_contig_ends` has no complete NCR length;
`missing_flank_genes` means the needed genes were not detected on this contig.
`multiple_candidates`, `ambiguous_flank_copies` and `flanks_not_adjacent` preserve
unresolved assignments. The locus's `observed_gap_length` is a candidate span,
not an SV size estimate.

With `--assembly-ambiguity`, compatible de novo primary-path or local mapping
structural uncertainty is transformed
alongside the sequence and written to separate `assembly_ambiguity.tsv` and
`assembly_ambiguity.gff3` files. The JSON record's `assembly_evidence.ncr_loci`
keeps each NCR's annotation status and path-evidence status separate. A bounded,
N-free `annotated` NCR can still have `unresolved_path_alternatives`; conversely,
`no_reported_path_ambiguity` does not certify repeat structure. See the
[coordinate handoff](annotation.md#carrying-assembly-uncertainty-into-annotation).

The evidence association also covers gaps below `--ncr-min-len`, touching genes
with zero intergenic bases, and overlapping flanks. A candidate with adjacent
tRNA-Ala and tRNA-Pro can therefore retain unresolved NCR1 evidence even though
it has no NCR1 FASTA/GFF interval. Do not interpret a zero candidate gap as a
read-supported deletion when its junction remains unresolved.

With `--gene-models` instead of `--species`, the command emits generic intervals
but applies no species catalog: `catalog` is null and the locus table has only its
header. It does not infer a species from a panel's accession or gene names.
