use super::{
    ncr::{NcrAnalysis, NcrRegion},
    Feature, Panel, Prediction, Record, Rotation,
};
use crate::io::fasta::write_fasta_records;
use anyhow::Result;
use serde_json::{json, Value};
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

pub(super) struct Annotated {
    pub sequence: Vec<u8>,
    pub prediction: Prediction,
    pub rotation: Rotation,
    pub circular: bool,
    pub ncr: NcrAnalysis,
    pub assembly_evidence: Option<super::uncertainty::Analysis>,
}

fn annotated_header(record: &Record, result: &Annotated) -> String {
    let mut fields: Vec<_> = record
        .header
        .split_whitespace()
        .enumerate()
        .filter(|(i, s)| {
            *i == 0
                || !(s.starts_with("annotation_topology=") || s.starts_with("annotation_origin="))
        })
        .map(|(_, s)| s.to_owned())
        .collect();
    fields.push(format!(
        "annotation_topology={}",
        if result.circular {
            "circular"
        } else {
            "linear"
        }
    ));
    if result.rotation.applied {
        fields.push("annotation_origin=tRNA-Pro".into());
    }
    fields.join(" ")
}
pub(super) fn counts(features: &[Feature], full: bool) -> Value {
    let count = |kind: &str| {
        features
            .iter()
            .filter(|f| f.kind == kind && (!full || f.full_length_supported))
            .count()
    };
    json!({"CDS": count("CDS"), "rRNA": count("rRNA"), "tRNA": count("tRNA")})
}
fn escape(s: &str) -> String {
    let mut result = String::new();
    for b in s.bytes() {
        if b.is_ascii_alphanumeric() || b".:^*$@!+_?-|".contains(&b) {
            result.push(b as char);
        } else {
            result.push_str(&format!("%{b:02X}"));
        }
    }
    result
}
fn tsv(s: &str) -> String {
    s.replace(['\t', '\n', '\r'], " ")
}
fn gff_feature(writer: &mut impl Write, contig: &str, f: &Feature) -> Result<()> {
    let mut consumed = 0;
    for (part, segment) in f.segments.iter().enumerate() {
        let phase = f
            .cds
            .as_ref()
            .map(|c| {
                // GFF phase is the bases to skip at the biological beginning of this
                // part. It continues through a split codon across the circular origin.
                ((c.phase + 3 - consumed % 3) % 3).to_string()
            })
            .unwrap_or_else(|| ".".into());
        writeln!(writer, "{}\tmitovarr\t{}\t{}\t{}\t{}\t{}\t{}\tID={};Name={};product={};part={}/{};status={};full_length_supported={};model_coverage={:.6};identity={:.6}{}",
            escape(contig), f.kind, segment.start, segment.end, f.alignment_score, f.strand, phase,
            f.id, escape(&f.name), escape(&f.product), part+1, f.segments.len(), f.status,
            f.full_length_supported, f.query_coverage, f.identity,
            if f.kind == "CDS" { ";transl_table=5" } else { ";rna_structure=not_assessed" })?;
        consumed += segment.end - segment.start + 1;
    }
    Ok(())
}

fn gff_ncr(writer: &mut impl Write, contig: &str, r: &NcrRegion) -> Result<()> {
    for (i, s) in r.segments.iter().enumerate() {
        let kind = if r.boundary_status == "supported" {
            "intergenic_region"
        } else {
            "region"
        };
        writeln!(writer,"{}\tmitovarr\t{}\t{}\t{}\t.\t.\t.\tID={};Name={};part={}/{};boundary_status={};ambiguous_bases={};left_gene={};right_gene={};sequence_origin=not_assessed",
            escape(contig),kind,s.start,s.end,r.id,r.name,i+1,r.segments.len(),r.boundary_status,r.ambiguous_bases,
            r.left_genes.iter().map(|x|escape(x)).collect::<Vec<_>>().join(","),r.right_genes.iter().map(|x|escape(x)).collect::<Vec<_>>().join(","))?;
    }
    Ok(())
}

fn gff_assembly(
    writer: &mut impl Write,
    contig: &str,
    r: &super::uncertainty::Region,
) -> Result<()> {
    for (i, s) in r.segments.iter().enumerate() {
        writeln!(writer,"{}\tmitovarr\tregion\t{}\t{}\t.\t.\t.\tID={};source_evidence={};reason={};alternatives={};part={}/{};ncr_loci={}",
            escape(contig),s.start,s.end,r.id,escape(&r.evidence),escape(&r.reason),r.alternatives,i+1,r.segments.len(),
            r.overlapping_ncrs.iter().map(|x|escape(x)).collect::<Vec<_>>().join(","))?;
    }
    Ok(())
}

fn segment_text(segments: &[super::Segment]) -> String {
    segments
        .iter()
        .map(|s| format!("{}..{}", s.start, s.end))
        .collect::<Vec<_>>()
        .join(",")
}
pub(super) fn write(
    out: &Path,
    records: &[Record],
    results: &[Annotated],
    panel: &Panel,
) -> Result<()> {
    let mut gff = BufWriter::new(File::create(out.join("annotations.gff3"))?);
    writeln!(gff,"##gff-version 3\n# model_accession={}\n# coordinates=annotated.fasta; method=homology_assisted; RNA_structure=not_assessed", escape(&panel.accession))?;
    let mut table = BufWriter::new(File::create(out.join("features.tsv"))?);
    let mut ncr_gff = BufWriter::new(File::create(out.join("ncr.gff3"))?);
    writeln!(ncr_gff,"##gff-version 3\n# coordinates=annotated.fasta; labels=positional; sequence_origin=not_assessed")?;
    let mut ncr_table = BufWriter::new(File::create(out.join("ncr.tsv"))?);
    writeln!(ncr_table,"contig\tid\tname\tlocus\tsegments_1based\tlength\tleft_genes\tright_genes\tboundary_status\tterminal_fragment\tambiguous_bases\tat_fraction_called_bases\tsequence_origin\tembedded_partial_homologs")?;
    let mut loci_table = BufWriter::new(File::create(out.join("ncr_loci.tsv"))?);
    writeln!(loci_table,"contig\tlocus\texpected_in_species_reference\tflank_genes\tstatus\tobserved_gap_length\toverlap_bases\tregion_ids\tmissing_flank_genes")?;
    writeln!(table,"contig\tid\tgene\ttype\tstrand\tsegments_1based\tlength\tstatus\tfull_length_supported\tquery_start\tquery_end\tmodel_length\tquery_coverage\tidentity\tgap_bases\tambiguous_bases\tboundary_5prime\tboundary_3prime\tRNA_structure\tCDS_stop\tCDS_internal_stops\tCDS_non_triplet_gaps\tcompeting_models")?;
    let mut genes = Vec::new();
    let mut proteins = Vec::new();
    let mut annotated = Vec::new();
    let mut original = Vec::new();
    let mut ncr_sequences = Vec::new();
    let mut assembly_files = if results.iter().any(|r| r.assembly_evidence.is_some()) {
        let mut table = BufWriter::new(File::create(out.join("assembly_ambiguity.tsv"))?);
        let mut gff = BufWriter::new(File::create(out.join("assembly_ambiguity.gff3"))?);
        writeln!(table,"sequence\tstart\tend\tevidence\treason\talternatives\tid\tinput_segments_1based\tgenes\tncr_loci")?;
        writeln!(gff,"##gff-version 3\n# coordinates=annotated.fasta; evidence=projected_assembly_intervals; structural_variants=not_inferred")?;
        Some((table, gff))
    } else {
        None
    };
    for (record, result) in records.iter().zip(results) {
        writeln!(
            gff,
            "##sequence-region {} 1 {}",
            escape(&record.name),
            record.sequence.len()
        )?;
        original.push((record.header.clone(), record.sequence.clone()));
        let header = annotated_header(record, result);
        annotated.push((header, result.sequence.clone()));
        if let (Some(evidence), Some((table, evidence_gff))) =
            (&result.assembly_evidence, &mut assembly_files)
        {
            writeln!(
                evidence_gff,
                "##sequence-region {} 1 {}",
                escape(&record.name),
                record.sequence.len()
            )?;
            for r in &evidence.regions {
                gff_assembly(&mut gff, &record.name, r)?;
                gff_assembly(evidence_gff, &record.name, r)?;
                for s in &r.segments {
                    writeln!(
                        table,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        tsv(&record.name),
                        s.start,
                        s.end,
                        tsv(&r.evidence),
                        tsv(&r.reason),
                        r.alternatives,
                        r.id,
                        segment_text(&r.input_segments),
                        r.overlapping_genes
                            .iter()
                            .map(|s| tsv(s))
                            .collect::<Vec<_>>()
                            .join(","),
                        r.overlapping_ncrs.join(",")
                    )?;
                }
            }
        }
        writeln!(
            ncr_gff,
            "##sequence-region {} 1 {}",
            escape(&record.name),
            record.sequence.len()
        )?;
        for r in &result.ncr.regions {
            gff_ncr(&mut gff, &record.name, r)?;
            gff_ncr(&mut ncr_gff, &record.name, r)?;
            let segments = r
                .segments
                .iter()
                .map(|s| format!("{}..{}", s.start, s.end))
                .collect::<Vec<_>>()
                .join(",");
            writeln!(
                ncr_table,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                tsv(&record.name),
                r.id,
                r.name,
                r.locus.as_deref().unwrap_or(""),
                segments,
                r.length,
                r.left_genes.join(","),
                r.right_genes.join(","),
                r.boundary_status,
                r.terminal_fragment,
                r.ambiguous_bases,
                r.at_fraction.map(|x| format!("{x:.6}")).unwrap_or_default(),
                r.sequence_origin,
                r.embedded_partial_homologs.join(",")
            )?;
            ncr_sequences.push((
                format!(
                    "{} name={} contig={} segments={} boundary_status={} ambiguous_bases={}",
                    r.id,
                    r.name,
                    escape(&record.name),
                    segments,
                    r.boundary_status,
                    r.ambiguous_bases
                ),
                r.sequence.clone(),
            ));
        }
        for locus in &result.ncr.loci {
            writeln!(
                loci_table,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                tsv(&record.name),
                locus.name,
                locus.expected_in_species_reference,
                locus.flank_genes.join(","),
                locus.status,
                locus
                    .observed_gap_length
                    .map(|x| x.to_string())
                    .unwrap_or_default(),
                locus.overlap_bases,
                locus.region_ids.join(","),
                locus.missing_flank_genes.join(",")
            )?;
        }
        for f in &result.prediction.features {
            gff_feature(&mut gff, &record.name, f)?;
            let parts = f
                .segments
                .iter()
                .map(|s| format!("{}..{}", s.start, s.end))
                .collect::<Vec<_>>()
                .join(",");
            let stops = f
                .cds
                .as_ref()
                .map(|c| {
                    c.internal_stop_codons
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .unwrap_or_default();
            writeln!(table,"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                tsv(&record.name), f.id, tsv(&f.name), f.kind, f.strand, parts, f.length, f.status, f.full_length_supported,
                f.query_start, f.query_end, f.model_length, f.query_coverage, f.identity, f.gap_bases, f.ambiguous_bases,
                f.five_prime_boundary, f.three_prime_boundary, f.rna_structure.as_deref().unwrap_or("not_applicable"),
                f.cds.as_ref().map(|c|c.stop.as_str()).unwrap_or("not_applicable"), stops,
                f.cds.as_ref().map(|c|c.non_triplet_alignment_gaps.to_string()).unwrap_or_default(), f.competing_models.join(","))?;
            let header = format!(
                "{} gene={} contig={} strand={} status={}",
                f.id,
                escape(&f.name),
                escape(&record.name),
                f.strand,
                f.status
            );
            genes.push((header.clone(), f.sequence.clone()));
            if let Some(cds) = &f.cds {
                proteins.push((header, cds.protein.as_bytes().to_vec()));
            }
        }
    }
    gff.flush()?;
    table.flush()?;
    ncr_gff.flush()?;
    ncr_table.flush()?;
    loci_table.flush()?;
    if let Some((mut table, mut gff)) = assembly_files {
        table.flush()?;
        gff.flush()?;
    }
    write_fasta_records(&out.join("input.fasta"), &original)?;
    write_fasta_records(&out.join("annotated.fasta"), &annotated)?;
    write_fasta_records(&out.join("genes.fna"), &genes)?;
    write_fasta_records(&out.join("proteins.faa"), &proteins)?;
    write_fasta_records(&out.join("ncr.fasta"), &ncr_sequences)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotation::{predict, GeneModel, Panel};
    #[test]
    fn gff_cds_phase_continues_across_origin() {
        let q = b"ATGGCTGCCGATCGTGGTACGCTCGCATTTGGGGCTGCTGCTTAA";
        let panel = Panel {
            schema: "mitovar-rust.gene-models.v1".into(),
            accession: "test".into(),
            genetic_code: 5,
            genes: vec![GeneModel {
                name: "x".into(),
                kind: "CDS".into(),
                sequence: String::from_utf8(q.to_vec()).unwrap(),
                product: String::new(),
            }],
        };
        let mut target = [q.as_slice(), b"CGCGCGCGCGCGCGCGCGCG"].concat();
        target.rotate_left(10);
        for seq in [
            target.clone(),
            crate::annotation::reverse_complement(&target),
        ] {
            let p = predict(&seq, &panel, true);
            let f = &p.features[0];
            assert_eq!(f.segments.len(), 2);
            let mut bytes = Vec::new();
            gff_feature(&mut bytes, "x", f).unwrap();
            let text = String::from_utf8(bytes).unwrap();
            let phases: Vec<_> = text
                .lines()
                .map(|l| l.split('\t').nth(7).unwrap())
                .collect();
            assert_eq!(phases, vec!["0", "2"]);
        }
        let mut partial_panel = Panel::bundled(crate::cli::Species::Cel).unwrap();
        partial_panel.genes.retain(|g| g.name == "ND4L");
        let coding_query = partial_panel.genes[0].sequence.as_bytes();
        let p = predict(&coding_query[1..], &partial_panel, false);
        let f = &p.features[0];
        assert!(!f.full_length_supported);
        let mut bytes = Vec::new();
        gff_feature(&mut bytes, "x", f).unwrap();
        let text = String::from_utf8(bytes).unwrap();
        assert_eq!(text.lines().next().unwrap().split('\t').nth(7), Some("2"));
    }
}
