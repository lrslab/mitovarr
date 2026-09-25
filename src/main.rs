use clap::Parser;

fn main() -> anyhow::Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env().unwrap_or_else(|_| "info".into()),
        )
        .with_writer(std::io::stderr)
        .init();
    let complete = mitovarr::run(mitovarr::cli::Cli::parse())?;
    if !complete {
        std::process::exit(2);
    }
    Ok(())
}
