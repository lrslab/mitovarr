use clap::{Args, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    name = "mitovarr",
    version,
    about = "Native Rust mtDNA de novo and mapping assembly, annotation and tRNA-Pro rotation"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Fetch mitochondrial reads and available mates from an indexed CaenDR BAM.
    FetchCaendr(FetchCaendrArgs),
    /// Iteratively recruit original complete read pairs with seed kmers.
    Bait(BaitArgs),
    /// Compatibility entry for direct assembly of original pairs without recruitment.
    #[command(hide = true)]
    AssembleDenovo(AssemblyArgs),
    /// BAM-baited de novo assembly from extracted pairs; no second recruitment.
    AssembleBaited(AssemblyArgs),
    /// Seeded de novo assembly; seed recruits reads only, validation uses all originals.
    AssembleSeeded(SeededArgs),
    /// Mapping assembly from a reference and indexed BAM, with optional local reconstruction.
    #[command(alias = "polish")]
    AssembleMapping(MappingArgs),
    /// Predict and assess genes on frozen candidates independently of NCR1 completeness.
    Annotate(AnnotationArgs),
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum Topology {
    /// Use annotation_topology when present, otherwise a literal circular=true header.
    Auto,
    Linear,
    /// Treat each input as a circular hypothesis; this does not validate closure.
    Circular,
}

#[derive(Debug, Args)]
pub struct AnnotationArgs {
    /// Candidate/consensus FASTA (all records, plain or gzip); never fed back to assembly.
    #[arg(long)]
    pub fasta: PathBuf,
    /// Bundled, pinned CDS/rRNA/tRNA homology models.
    #[arg(
        long,
        value_enum,
        required_unless_present = "gene_models",
        conflicts_with = "gene_models"
    )]
    pub species: Option<Species>,
    /// Custom mitovar-rust.gene-models.v1 JSON panel (translation table 5).
    #[arg(long)]
    pub gene_models: Option<PathBuf>,
    #[arg(long)]
    pub out: PathBuf,
    #[arg(long, value_enum, default_value = "auto")]
    pub topology: Topology,
    /// Orient a circular input so its unique, complete tRNA-Pro starts at base 1 on +.
    #[arg(long)]
    pub rotate_to_trna_pro: bool,
    /// Project assembly ambiguity intervals from matching FASTA coordinates.
    #[arg(long)]
    pub assembly_ambiguity: Option<PathBuf>,
    /// Minimum intergenic length for NCR intervals; shorter named loci remain in ncr_loci.tsv.
    #[arg(long, default_value_t = 20, value_parser = clap::value_parser!(u32).range(1..))]
    pub ncr_min_len: u32,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum Species {
    #[value(alias = "c_elegans")]
    Cel,
    #[value(alias = "c_tropicalis")]
    Ctro,
    #[value(alias = "c_briggsae")]
    Cbr,
}

impl Species {
    pub fn name(self) -> &'static str {
        match self {
            Self::Cel => "c_elegans",
            Self::Ctro => "c_tropicalis",
            Self::Cbr => "c_briggsae",
        }
    }
}

#[derive(Debug, Args)]
pub struct FetchCaendrArgs {
    #[arg(long, value_enum)]
    pub species: Species,
    #[arg(long)]
    pub strain: String,
    #[arg(long)]
    pub out: PathBuf,
    /// Override automatic mitochondrial contig-name detection.
    #[arg(long)]
    pub contig: Option<String>,
    /// Timeout for each HTTP request.
    #[arg(long, default_value_t = 120)]
    pub timeout_seconds: u64,
    /// Bound the number of mitochondrial templates retained for mate recovery.
    #[arg(long, default_value_t = 1_000_000)]
    pub max_templates: usize,
}

#[derive(Debug, Clone, Args)]
pub struct Reads {
    #[arg(long)]
    pub r1: PathBuf,
    #[arg(long)]
    pub r2: PathBuf,
}

#[derive(Debug, Clone, Args)]
pub struct BaitSettings {
    #[arg(long, default_value_t = 21, value_parser = parse_k)]
    pub bait_k: usize,
    #[arg(long, default_value_t = 30)]
    pub max_rounds: usize,
    #[arg(long, default_value_t = 2_000_000)]
    pub max_pairs: usize,
}

#[derive(Debug, Clone, Args)]
pub struct Limits {
    /// Hard limit on distinct kmers in each recruitment or assembly set.
    #[arg(long, default_value_t = 2_000_000)]
    pub max_kmers: usize,
    /// Optional hard base-quality cutoff; 0 disables it. Qualities still inform read scores.
    #[arg(long, default_value_t = 0, value_parser = clap::value_parser!(u8).range(0..=93))]
    pub min_baseq: u8,
}

#[derive(Debug, Args)]
pub struct BaitArgs {
    #[command(flatten)]
    pub reads: Reads,
    #[arg(long)]
    pub seed: PathBuf,
    #[arg(long)]
    pub out: PathBuf,
    #[command(flatten)]
    pub bait: BaitSettings,
    #[command(flatten)]
    pub limits: Limits,
}

#[derive(Debug, Clone, Args)]
pub struct GraphSettings {
    /// Also propose a coverage-fitted circular layout from sample graph paths; never certifies a final assembly.
    #[arg(long)]
    pub coverage_layout: bool,
    /// Add cross-k layout proposals and read-context-guided walks; hypotheses remain unverified.
    #[arg(long, requires = "coverage_layout")]
    pub multi_k_layout: bool,
    /// Main odd k in 15..127; smaller k=47/31 graphs also propose open-end connections.
    #[arg(long, default_value_t = 63, value_parser = parse_graph_k)]
    pub k: usize,
    /// Minimum pair-record support for a retained kmer (mates counted once).
    #[arg(long, default_value_t = 3)]
    pub min_kmer_count: u32,
    /// Maximum weak/strong support ratio for tips, weak links and repeat routes; 0 disables all cleanup.
    #[arg(long, default_value_t = 0.2)]
    pub max_error_fraction: f64,
    /// Maximum weak/strong mean coverage ratio for bounded local bubbles; does not relax repeat routes.
    #[arg(long, default_value_t = 0.49)]
    pub max_bubble_fraction: f64,
    /// Minimum matching-base fragment depth after gapped validation.
    #[arg(long, default_value_t = 10)]
    pub min_depth: u32,
    /// Minimum named fragments supporting each repeat bridge and the output junction.
    #[arg(long, default_value_t = 3)]
    pub min_junction_fragments: u32,
    /// Hard limit on unique names retained by final read validation.
    #[arg(long, default_value_t = 2_000_000)]
    pub max_validation_pairs: usize,
}

#[derive(Debug, Args)]
pub struct AssemblyArgs {
    #[command(flatten)]
    pub reads: Reads,
    #[arg(long)]
    pub out: PathBuf,
    #[command(flatten)]
    pub graph: GraphSettings,
    #[command(flatten)]
    pub limits: Limits,
}

#[derive(Debug, Args)]
pub struct SeededArgs {
    #[command(flatten)]
    pub assembly: AssemblyArgs,
    #[arg(long)]
    pub seed: PathBuf,
    #[command(flatten)]
    pub bait: BaitSettings,
}

#[derive(Debug, Clone, Args)]
pub struct MappingArgs {
    /// Build a sample graph from paired FASTQ and output a fast consensus with unresolved regions.
    #[arg(long, requires_all = ["r1", "r2"])]
    pub local_assembly: bool,
    /// Enumerate local structures and compute complete paired-fragment scores (slower; opt-in).
    #[arg(long, requires = "local_assembly")]
    pub deep_validation: bool,
    /// Retain a usable consensus even if structure remains unresolved; enabled by default in fast local mode.
    #[arg(long, alias = "legacy-bridge", requires = "local_assembly")]
    pub best_effort: bool,
    /// Infer separate coverage-guided local layouts; bounded candidates never certify final structure.
    #[arg(long, requires = "deep_validation")]
    pub coverage_layout: bool,
    /// Add cross-k local layout proposals with long flanks and shared read scoring.
    #[arg(long, requires = "coverage_layout")]
    pub multi_k_layout: bool,
    /// Complete paired reads for local assembly (e.g. the matching fetch-caendr FASTQ).
    #[arg(long, requires_all = ["r2", "local_assembly"])]
    pub r1: Option<PathBuf>,
    #[arg(long, requires_all = ["r1", "local_assembly"])]
    pub r2: Option<PathBuf>,
    /// Sample graph k for local assembly; independent of the de novo main k.
    #[arg(long, default_value_t = 31, value_parser = parse_k)]
    pub local_k: usize,
    /// Maximum sequence length between local anchors, including the anchors.
    #[arg(long, default_value_t = 2000)]
    pub local_max_span: usize,
    /// Maximum local paths; reaching a search limit leaves the window unresolved.
    #[arg(long, default_value_t = 128)]
    pub local_max_paths: usize,
    /// Worker threads for deep validation; results keep input order. Capped by available CPUs.
    #[arg(long, default_value_t = 4, value_parser = clap::value_parser!(u16).range(1..=64))]
    pub scoring_threads: u16,
    #[arg(long = "ref")]
    pub reference: PathBuf,
    #[arg(long)]
    pub bam: PathBuf,
    #[arg(long)]
    pub contig: Option<String>,
    #[arg(long)]
    pub out: PathBuf,
    #[arg(long, default_value_t = 20, value_parser = clap::value_parser!(u8).range(0..=60))]
    pub min_mapq: u8,
    /// Optional hard base-quality cutoff; 0 disables it. Quality-weighted voting and scoring remain.
    #[arg(long, default_value_t = 0, value_parser = clap::value_parser!(u8).range(0..=93))]
    pub min_baseq: u8,
    #[arg(long, default_value_t = 10)]
    pub min_depth: u32,
    /// Minimum weighted SNP fraction; a strict majority (> 0.5) is always required.
    #[arg(long, default_value_t = 0.5)]
    pub min_snp_frac: f64,
    /// Minimum fragment support fraction for an insertion or deletion.
    #[arg(long, default_value_t = 0.7)]
    pub min_indel_frac: f64,
    /// Legacy override for both SNP and indel fractions; prefer the separate options.
    #[arg(long, conflicts_with_all = ["min_snp_frac", "min_indel_frac"])]
    pub min_alt_frac: Option<f64>,
    #[arg(long, default_value_t = 100)]
    pub window: usize,
    /// Bound records retained for deterministic fragment-aware pileup.
    #[arg(long, default_value_t = 500_000)]
    pub max_records: usize,
    /// Deterministically subsample complete named fragments to the record cap.
    #[arg(long)]
    pub subsample: bool,
    /// Seed for fragment priority sampling (not a data integrity fingerprint).
    #[arg(long, default_value_t = 42, requires = "subsample")]
    pub subsample_seed: u64,
    #[arg(long)]
    pub emit_variants: bool,
    #[arg(long)]
    pub emit_coverage: bool,
}

fn parse_k(value: &str) -> Result<usize, String> {
    let k = value.parse::<usize>().map_err(|_| "k must be an integer")?;
    if (15..=31).contains(&k) {
        Ok(k)
    } else {
        Err("k must be in 15..=31".into())
    }
}

fn parse_graph_k(value: &str) -> Result<usize, String> {
    let k = value.parse::<usize>().map_err(|_| "k must be an integer")?;
    if (15..=127).contains(&k) && k % 2 == 1 {
        Ok(k)
    } else {
        Err("assembly k must be odd in 15..=127".into())
    }
}
