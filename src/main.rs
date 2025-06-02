use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Array2;
use rayon::prelude::*;
use std::path::PathBuf;

mod ms_data;
mod similarity;

use ms_data::{MSRun, Spectrum};
use similarity::{cosine_similarity, modified_cosine_similarity, SimilarityMethod};

#[derive(Parser)]
#[command(name = "ms-similarity-tool")]
#[command(about = "A tool for untargeted mass spectrometry data similarity analysis")]
struct Cli {
    /// Input mzML or mzXML file
    #[arg(short, long)]
    input: PathBuf,

    /// Output file for similarity matrix (CSV format)
    #[arg(short, long)]
    output: PathBuf,

    /// Similarity method to use
    #[arg(short, long, default_value = "cosine")]
    method: SimilarityMethod,

    /// Minimum intensity threshold for peaks
    #[arg(long, default_value = "1000.0")]
    min_intensity: f64,

    /// Mass tolerance for modified cosine similarity (in Da)
    #[arg(long, default_value = "0.02")]
    mass_tolerance: f64,

    /// Maximum number of spectra to process (for testing)
    #[arg(long)]
    max_spectra: Option<usize>,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

fn main() -> Result<()> {
    let args = Cli::parse();

    if args.verbose {
        println!("Loading MS data from: {:?}", args.input);
    }

    // Load MS data
    let ms_run = MSRun::from_file(&args.input)
        .context("Failed to load MS data")?;

    if args.verbose {
        println!("Loaded {} spectra", ms_run.spectra.len());
    }

    // Filter spectra by intensity and limit if requested
    let filtered_spectra: Vec<_> = ms_run.spectra
        .into_iter()
        .filter(|spec| spec.max_intensity() >= args.min_intensity)
        .take(args.max_spectra.unwrap_or(usize::MAX))
        .collect();

    if args.verbose {
        println!("Processing {} spectra after filtering", filtered_spectra.len());
    }

    if filtered_spectra.is_empty() {
        anyhow::bail!("No spectra found meeting the criteria");
    }

    // Calculate similarity matrix
    let similarity_matrix = calculate_similarity_matrix(
        &filtered_spectra,
        args.method,
        args.mass_tolerance,
        args.verbose,
    )?;

    // Save results
    save_similarity_matrix(&similarity_matrix, &filtered_spectra, &args.output)
        .context("Failed to save similarity matrix")?;

    if args.verbose {
        println!("Similarity matrix saved to: {:?}", args.output);
        println!("Matrix dimensions: {}x{}", similarity_matrix.nrows(), similarity_matrix.ncols());
    }

    Ok(())
}

fn calculate_similarity_matrix(
    spectra: &[Spectrum],
    method: SimilarityMethod,
    mass_tolerance: f64,
    verbose: bool,
) -> Result<Array2<f64>> {
    let n = spectra.len();
    let mut matrix = Array2::zeros((n, n));
    
    let progress_bar = if verbose {
        let pb = ProgressBar::new((n * (n - 1) / 2) as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
                .unwrap()
                .progress_chars("#>-"),
        );
        Some(pb)
    } else {
        None
    };

    // Calculate upper triangle of similarity matrix in parallel
    let similarities: Vec<_> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            (i + 1..n).into_par_iter().map(move |j| {
                let similarity = match method {
                    SimilarityMethod::Cosine => cosine_similarity(&spectra[i], &spectra[j]),
                    SimilarityMethod::ModifiedCosine => {
                        modified_cosine_similarity(&spectra[i], &spectra[j], mass_tolerance)
                    }
                };
                (i, j, similarity)
            })
        })
        .collect();

    // Fill the matrix
    for (i, j, similarity) in similarities {
        matrix[[i, j]] = similarity;
        matrix[[j, i]] = similarity; // Symmetric matrix
        
        if let Some(ref pb) = progress_bar {
            pb.inc(1);
        }
    }

    // Set diagonal to 1.0 (self-similarity)
    for i in 0..n {
        matrix[[i, i]] = 1.0;
    }

    if let Some(pb) = progress_bar {
        pb.finish_with_message("Similarity calculation complete");
    }

    Ok(matrix)
}

fn save_similarity_matrix(
    matrix: &Array2<f64>,
    spectra: &[Spectrum],
    output_path: &PathBuf,
) -> Result<()> {
    let mut writer = csv::Writer::from_path(output_path)?;

    // Write header with retention times
    let mut header = vec!["RT".to_string()];
    for spec in spectra {
        header.push(format!("{:.3}", spec.retention_time));
    }
    writer.write_record(&header)?;

    // Write each row
    for (i, spec) in spectra.iter().enumerate() {
        let mut row = vec![format!("{:.3}", spec.retention_time)];
        for j in 0..matrix.ncols() {
            row.push(format!("{:.6}", matrix[[i, j]]));
        }
        writer.write_record(&row)?;
    }

    writer.flush()?;
    Ok(())
}