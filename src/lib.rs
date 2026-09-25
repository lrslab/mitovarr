use anyhow::{ensure, Result};
use serde_json::{json, Value};
use std::{
    fs,
    io::{BufWriter, Write},
    path::Path,
};

pub mod annotation;
pub mod bait;
pub mod caendr;
pub mod cli;
pub mod consensus;
pub mod denovo;
pub mod io;
pub mod kmer;
mod mapping;
pub mod pileup;
pub mod polish;
pub mod repeat;
pub mod report;
mod subsample;

pub fn prepare_output(out: &Path) -> Result<()> {
    ensure!(
        !out.exists() || (out.is_dir() && fs::read_dir(out)?.next().is_none()),
        "output directory must be new or empty: {}",
        out.display()
    );
    fs::create_dir_all(out)?;
    Ok(())
}

pub fn write_json(path: &Path, value: &Value) -> Result<()> {
    let temporary = path.with_extension("json.tmp");
    let mut writer = BufWriter::new(fs::File::create(&temporary)?);
    serde_json::to_writer_pretty(&mut writer, value)?;
    writeln!(writer)?;
    writer.flush()?;
    fs::rename(temporary, path)?;
    Ok(())
}

pub fn run(cli: cli::Cli) -> Result<bool> {
    use cli::Commands;
    let out = match &cli.command {
        Commands::FetchCaendr(a) => &a.out,
        Commands::Bait(a) => &a.out,
        Commands::AssembleDenovo(a) => &a.out,
        Commands::AssembleBaited(a) => &a.out,
        Commands::AssembleSeeded(a) => &a.assembly.out,
        Commands::Polish(a) => &a.out,
        Commands::Annotate(a) => &a.out,
    }
    .clone();
    prepare_output(&out)?;
    let result = (|| match cli.command {
        Commands::FetchCaendr(args) => caendr::fetch(&args),
        Commands::Bait(args) => {
            let recruited = bait::recruit(&args.reads, &args.seed, &out, &args.bait, &args.limits)?;
            write_json(&out.join("report.json"), &recruited.report)?;
            Ok(recruited.complete)
        }
        Commands::AssembleDenovo(args) => denovo::assemble(
            &args.reads,
            &args.reads,
            &out,
            &args.graph,
            &args.limits,
            None,
            false,
        ),
        Commands::AssembleBaited(args) => denovo::assemble(
            &args.reads,
            &args.reads,
            &out,
            &args.graph,
            &args.limits,
            None,
            true,
        ),
        Commands::AssembleSeeded(args) => {
            let a = &args.assembly;
            let recruited = bait::recruit(
                &a.reads,
                &args.seed,
                &out.join("bait"),
                &args.bait,
                &a.limits,
            )?;
            if !recruited.complete {
                write_json(
                    &out.join("report.json"),
                    &json!({"schema":"mitovar-rust.assembly.v1","status":"unresolved","reason":"recruitment_not_usable","final_published":false,"recruitment":recruited.report}),
                )?;
                return Ok(false);
            }
            denovo::assemble(
                &recruited.reads,
                &a.reads,
                &out,
                &a.graph,
                &a.limits,
                Some(recruited.report),
                false,
            )
        }
        Commands::Polish(args) => polish::polish(&args),
        Commands::Annotate(args) => annotation::run(&args),
    })();
    if let Err(error) = &result {
        // This directory was proven empty before execution; never touch prior runs.
        let _ = fs::remove_file(out.join("consensus.fasta"));
        let _ = write_json(
            &out.join("report.json"),
            &json!({"schema":"mitovar-rust.failure.v1","status":"failed","final_published":false,"error":format!("{error:#}")}),
        );
    }
    result
}
