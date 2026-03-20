// use std::collections::HashMap;
use hashbrown::HashMap;
use ndarray::Array2;
use std::process;
use std::fmt; // Import fmt for Display implementation
use std::path::Path;
use std::io::{self, Write}; // Import io and Write traits for flushing output
use std::time::Instant;
use std::error::Error;
use std::fs::File;

mod ms_io;
use ms_io::{import_mzml, OutputFormat, write_similarity_matrix};
mod similarity;
use similarity::{
    compute_pairwise_similarity_matrix_sparse,
    prune_background_bins_sparse,
    prune_background_similarity_regions,
    compute_sparse_vec_map,
};

mod visual;
use visual::{plot_similarity_heatmap, ImageFormat, ColorTheme, ThemeName, ChromatogramType};

mod processing;
use processing::{ProcessingConfig, NormMethod, ThresholdMethod};

// Custom error types for better error handling
#[derive(Debug)]
enum CliError {
    InvalidInput(String),
    FileError(String),
    ProcessingError(String),
    IoError(io::Error),
}

impl fmt::Display for CliError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CliError::InvalidInput(msg) => write!(f, "Invalid input: {}", msg),
            CliError::FileError(msg) => write!(f, "File error: {}", msg),
            CliError::ProcessingError(msg) => write!(f, "Processing error: {}", msg),
            CliError::IoError(err) => write!(f, "IO error: {}", err),
        }
    }
}

impl Error for CliError {}

impl From<io::Error> for CliError {
    fn from(err: io::Error) -> Self {
        CliError::IoError(err)
    }
}

type Result<T> = std::result::Result<T, CliError>;

fn main() {
    let art = "  ____   _    ____ ____         ____ _     ___ \n |  _ \\ / \\  / ___/ ___|       / ___| |   |_ _|\n | |_) / _ \\ \\___ \\___ \\ _____| |   | |    | | \n |  __/ ___ \\ ___) |__) |_____| |___| |___ | | \n |_| /_/   \\_\\____/____/       \\____|_____|___|\n";
    println!("{art}");

    //First, give a welcome and description of what this command line interface does.
    println!("**Welcome to PASS-CLI (Pairwise Analyzer for Spectral Similarity)**
----------------------------------------------------------------------------------------------
||  -A Rust-powered CLI for mass spectrometry analysis.                                     
||  -Parsing `.mzML` files for spectral data                                                
||  -Computing pairwise similarity with cosine and modified cosine metrics                  
||  -Building exportable similarity matrices for downstream visualization and analysis      
||  -Supports file format: `.mzML`
----------------------------------------------------------------------------------------------
||  Please follow the directions to process your data.
----------------------------------------------------------------------------------------------"
    );

    // Use a proper error handling approach with early returns
    if let Err(e) = run_cli() {
        eprintln!("Error: {}", e);
        eprintln!("Program terminated with errors.");
        process::exit(1);
    }
}

fn run_cli() -> Result<()> {
    let input_path = prompt_input_path()?;
    let (output_path, output_format) = prompt_output_path()?;
    let ms1_similarity_metrics = prompt_similarity_metrics(1.0)?;
    let ms2_similarity_metrics = prompt_similarity_metrics(2.0)?;
    let ms1_threshold = prompt_threshold_method(true)?;
    let ms2_threshold = prompt_threshold_method(false)?;
    let mass_tolerance = prompt_mass_tolerance()?;

    // Visuals config
    let heatmap_enabled = prompt_generate_heatmap()?;
    
    // Background pruning config
    let pre_similarity_prune_fraction = prompt_pre_similarity_prune_fraction()?;
    let (enable_pruning, similarity_threshold, min_background_fraction) = prompt_enable_background_pruning()?;

    // If heatmap is enabled, get all config at once including chromatogram options
    let (image_path, image_format, theme, chromatogram_type, enable_smoothing) = if heatmap_enabled {
        let path = prompt_output_image_path()?;
        let format = prompt_image_format()?;
        let theme = prompt_color_theme()?;
        let chrom_type = prompt_chromatogram_type()?;
        let smoothing = prompt_enable_smoothing()?;
        (path, format, theme, chrom_type, smoothing)
    } else {
        ("".to_string(), ImageFormat::Png, ThemeName::Classic.get_theme(), ChromatogramType::BPI, false)
    };

    // Print a final confirmation
    print_configuration(&input_path,
        &output_path,
        &output_format,
        &ms1_similarity_metrics,
        &ms2_similarity_metrics,
        ms1_threshold,
        ms2_threshold,
        mass_tolerance,
        heatmap_enabled,
        Some(&image_path),
        Some(&image_format),
        Some(&theme),
        Some(chromatogram_type),
        Some(enable_smoothing)
    );

    // Write the configuration to a file
    let config = Config {
        input_path: input_path.clone(),
        output_path: output_path.clone(),
        output_format: output_format.clone(),
        ms1_similarity_metrics: ms1_similarity_metrics.clone(),
        ms2_similarity_metrics: ms2_similarity_metrics.clone(),
        ms1_threshold,
        ms2_threshold,
        mass_tolerance,
        heatmap_enabled,
        image_path: Some(image_path.clone()),
        image_format: Some(image_format.clone()),
        theme: Some(theme.clone()),
    };
    write_config_to_file(&config)?;

    if !confirm_processing()? {
        println!("**Processing cancelled by user.**");
        return Ok(());
    }

    let start = Instant::now();
    
    // Process the data with proper error handling
    process_spectral_data(
        &input_path, &output_path, &image_path, output_format,
        &ms1_similarity_metrics, &ms2_similarity_metrics,
        ms1_threshold, ms2_threshold,
        mass_tolerance,
        heatmap_enabled, Some(&image_format), Some(&theme),
        chromatogram_type, enable_smoothing,
        pre_similarity_prune_fraction,
        enable_pruning, similarity_threshold, min_background_fraction
    )?;

    let duration = start.elapsed();
    println!("||    Success! Processing completed in {duration:.2?}");
    Ok(())
}

fn process_spectral_data(
    input_path: &str,
    output_path: &str,
    image_path: &str,
    output_format: OutputFormat,
    ms1_similarity_metrics: &[&str],
    ms2_similarity_metrics: &[&str],
    ms1_threshold: ThresholdMethod,
    ms2_threshold: ThresholdMethod,
    mass_tolerance: f32,
    heatmap_enabled: bool,
    image_format: Option<&ImageFormat>,
    theme: Option<&ColorTheme>,
    chromatogram_type: ChromatogramType,
    enable_smoothing: bool,
    pre_similarity_prune_fraction: f64,
    enable_pruning: bool,
    similarity_threshold: f64,
    min_background_fraction: f64
) -> Result<()> {
    // Importing the data file with better error handling
    println!("||    Loading spectral data from: {}", input_path);
    let spec_map_result = import_mzml(input_path)
        .map_err(|e| CliError::FileError(format!("Failed to import mzML file '{}': {:?}", input_path, e)))?;
    let (spec_map, spec_metadata) = spec_map_result;

    if spec_map.is_empty() {
        return Err(CliError::ProcessingError("No spectra found in the input file".to_string()));
    }
    println!("||    Parsed {} spectra successfully.", spec_map.len());

    // Partition into MS1/MS2 in a single pass (avoids cloning the full map twice)
    let (ms1_spec_map, ms2_spec_map): (hashbrown::HashMap<_, _>, hashbrown::HashMap<_, _>) =
        spec_map.into_iter().partition(|(scan_id, _)| {
            spec_metadata.get(scan_id).map(|m| m.ms_level == 1).unwrap_or(false)
        });
    println!("||    Num. MS1 spectra: {}, Num. MS2 spectra: {}", ms1_spec_map.len(), ms2_spec_map.len());

    // Compute max m/z before preprocessing (used for binning config and later)
    let max_ms1_mz = ms1_spec_map.values()
        .map(|p| p.iter().map(|p| p.mz).fold(0., f32::max))
        .fold(0., f32::max);
    let max_ms2_mz = ms2_spec_map.values()
        .map(|p| p.iter().map(|p| p.mz).fold(0., f32::max))
        .fold(0., f32::max);

    // ------------------ Processing: noise → (optional) background → normalization ------------------

    // Choose sensible defaults (can move to CLI/config later)
    let proc_cfg = ProcessingConfig {
        bin_width: mass_tolerance,               // reuse your CLI param
        max_mz: max_ms1_mz.max(1.0),             // guard
        threshold: ms1_threshold,
        snr_min: 5.0,                            // MAD-based S/N ≥ 5
        norm: NormMethod::PQN,                   // robust global scaling
        blank_ratio_max: 0.20,
        sample_over_blank_min: 3.0,
    };

    // If you have blank scans, list their scan IDs here; otherwise pass None.
    let blank_scan_ids: Option<Vec<String>> = None; // or Some(vec!["scan=1".into(), "scan=5".into(), ...]);

    let proc_result = processing::preprocess_after_import(
        &ms1_spec_map,
        &spec_metadata,
        &proc_cfg,
        blank_scan_ids.as_ref().map(|v| v.as_slice())
    );

    // Replace MS1 with cleaned & normalized MS1
    let ms1_spec_map = proc_result.cleaned_ms1;

    // Apply the *same* normalization factors to MS2 so MS1/MS2 stay on a consistent scale per scan.
    let ms2_spec_map = processing::apply_scale_factors(ms2_spec_map, &proc_result.scale_factors);
    // Apply adaptive threshold to MS2 (SNR + intensity gate, no normalization).
    let ms2_spec_map = processing::apply_threshold_to_spec_map(ms2_spec_map, ms2_threshold, 5.0);

    // Quick log
    println!("||    Processing summary:");
    println!("||        - scans: {}", proc_result.stats.n_scans);
    println!("||        - noise-removed peaks: {}", proc_result.stats.n_removed_peaks_noise);
    println!("||        - background bins: {}", proc_result.stats.n_background_bins);
    println!("||        - norm method: {:?}, median factor: {:.3}", proc_result.stats.norm_method, proc_result.stats.norm_factors_median);
    println!("------------------------------------------------------------------------------");

    if max_ms1_mz == 0. && max_ms2_mz == 0. {
        return Err(CliError::ProcessingError("**No MS1 or MS2 data available in the input file**".to_string()));
    }

    let max_mz = f32::max(max_ms1_mz, max_ms2_mz);
    let vector_length = (max_mz / mass_tolerance).ceil() as usize;
    
    println!("||    Processing parameters:");
    println!("||        - m/z range: 0.0 to {:.2}", max_mz);
    println!("||        - Bin width: {:.4}", mass_tolerance);
    println!("||        - Vector size: {}", vector_length);
    println!("------------------------------------------------------------------------------");

    // Compute sparse binary vectors with progress indication
    println!("||    Computing sparse binary vectors...");
    let mut ms1_bits_map = compute_sparse_vec_map(&ms1_spec_map, mass_tolerance, max_ms1_mz, 0.0)
        .map_err(|e| CliError::ProcessingError(format!("**Failed to compute MS1 sparse vectors: {:?}**", e)))?;
    let mut ms2_bits_map = compute_sparse_vec_map(&ms2_spec_map, mass_tolerance, max_ms2_mz, 0.0)
        .map_err(|e| CliError::ProcessingError(format!("**Failed to compute MS2 sparse vectors: {:?}**", e)))?;
    println!("||    ✔ Computed sparse binary vectors successfully.");
    println!("------------------------------------------------------------------------------");

    // Filter background signals
    println!("||    Pruning background signals (threshold: {:.0}%)...", pre_similarity_prune_fraction * 100.0);
    prune_background_bins_sparse(&mut ms1_bits_map, pre_similarity_prune_fraction)
        .map_err(|e| CliError::ProcessingError(format!("**Failed to prune MS1 background bins: {:?}**", e)))?;
    prune_background_bins_sparse(&mut ms2_bits_map, pre_similarity_prune_fraction)
        .map_err(|e| CliError::ProcessingError(format!("**Failed to prune MS2 background bins: {:?}**", e)))?;
    println!("||    ✔ Pruned background signals successfully.");
    println!("------------------------------------------------------------------------------");

    // Compute pairwise cosine similarities
    let mut ms1_results: HashMap<String, (Vec<String>, Array2<f32>)> = HashMap::new();
    let mut ms2_results: HashMap<String, (Vec<String>, Array2<f32>)> = HashMap::new();

    // Process MS1 metrics
    for metric in ms1_similarity_metrics {
        println!("||    Computing MS1 pairwise similarity matrix using {}...", metric);
        let (ms1_scans, ms1_mat) = compute_pairwise_similarity_matrix_sparse(
            &ms1_bits_map, &spec_metadata, metric.to_string(), 1, mass_tolerance
        ).map_err(|e| CliError::ProcessingError(format!("**Failed to compute MS1 {} matrix: {:?}**", metric, e)))?;
        ms1_results.insert(metric.to_string(), (ms1_scans, ms1_mat));
    }
    println!("||    ✔ Computed MS1 pairwise similarity matrices successfully.");
    println!("------------------------------------------------------------------------------");

    // Process MS2 metrics
    for metric in ms2_similarity_metrics {
        if ms2_bits_map.is_empty() {
            println!("**Skipping MS2 {} - no MS2 spectra available**", metric);
            continue;
        }
        println!("||    Computing MS2 pairwise similarity matrix using {}...", metric);
        let (ms2_scans, ms2_mat) = compute_pairwise_similarity_matrix_sparse(
            &ms2_bits_map, &spec_metadata, metric.to_string(), 2, mass_tolerance
        ).map_err(|e| CliError::ProcessingError(format!("**Failed to compute MS2 {} matrix: {:?}**", metric, e)))?;
        ms2_results.insert(metric.to_string(), (ms2_scans, ms2_mat));
    }
    println!("||    ✔ Computed MS2 pairwise similarity matrices successfully.");
    println!("------------------------------------------------------------------------------");

    // Apply post-similarity background pruning if enabled
    if enable_pruning {
        println!("||    Applying background pruning to similarity matrices...");
        
        // Prune MS1 results
        let mut total_ms1_pruned = 0;
        for (metric, (_, ms1_mat)) in ms1_results.iter_mut() {
            match prune_background_similarity_regions(ms1_mat, similarity_threshold, min_background_fraction) {
                Ok(pruned_count) => {
                    total_ms1_pruned += pruned_count;
                    if pruned_count > 0 {
                        println!("||    Zeroed {} background spectra in MS1 {} matrix", pruned_count, metric);
                    }
                },
                Err(e) => println!("||    Warning: Failed to prune MS1 {} matrix: {:?}", metric, e),
            }
        }
        
        // Prune MS2 results
        let mut total_ms2_pruned = 0;
        for (metric, (_, ms2_mat)) in ms2_results.iter_mut() {
            match prune_background_similarity_regions(ms2_mat, similarity_threshold, min_background_fraction) {
                Ok(pruned_count) => {
                    total_ms2_pruned += pruned_count;
                    if pruned_count > 0 {
                        println!("||    Zeroed {} background spectra in MS2 {} matrix", pruned_count, metric);
                    }
                },
                Err(e) => println!("||    Warning: Failed to prune MS2 {} matrix: {:?}", metric, e),
            }
        }
        
        if total_ms1_pruned > 0 || total_ms2_pruned > 0 {
            println!("||    ✔ Background pruning complete. Total zeroed: {} MS1, {} MS2 spectra", total_ms1_pruned, total_ms2_pruned);
        } else {
            println!("||    No background regions detected with current thresholds");
        }
        println!("------------------------------------------------------------------------------");
    }

    // Export results
    println!("||    Exporting similarity matrices...");
    
    let out_base = strip_extension(output_path);
    let out_ext = output_format.to_string();

    for (metric, (ms1_scans, ms1_mat)) in &ms1_results {
        let ms1_output_path = format!("{}_ms1_{}.{}", out_base, metric, out_ext);
        write_similarity_matrix(ms1_scans, ms1_mat, &ms1_output_path, output_format)
            .map_err(|e| CliError::FileError(format!("**Failed to write MS1 {} matrix: {}**", metric, e)))?;
        println!("||    Exported MS1 {} matrix to: {}", metric, ms1_output_path);
    }

    for (metric, (ms2_scans, ms2_mat)) in &ms2_results {
        let ms2_output_path = format!("{}_ms2_{}.{}", out_base, metric, out_ext);
        write_similarity_matrix(ms2_scans, ms2_mat, &ms2_output_path, output_format)
            .map_err(|e| CliError::FileError(format!("**Failed to write MS2 {} matrix: {}**", metric, e)))?;
        println!("||    Exported MS2 {} matrix to: {}", metric, ms2_output_path);
    }

    if heatmap_enabled {
        println!("||    Creating heatmaps...");

        // Calculate global scan range from all data (MS1 and MS2) for consistent scaling
        let mut all_scan_ids: Vec<String> = Vec::new();
        for (_, (scans, _)) in &ms1_results {
            all_scan_ids.extend(scans.clone());
        }
        for (_, (scans, _)) in &ms2_results {
            all_scan_ids.extend(scans.clone());
        }
        
        let global_scan_range = if !all_scan_ids.is_empty() {
            let mut min_scan = usize::MAX;
            let mut max_scan = 0;
            for scan_id in &all_scan_ids {
                if let Ok(scan_num) = scan_id.parse::<usize>() {
                    min_scan = min_scan.min(scan_num);
                    max_scan = max_scan.max(scan_num);
                }
            }
            Some((min_scan, max_scan))
        } else {
            None
        };
        
        // println!("||    Global scan range for consistent heatmap scaling: {:?}", global_scan_range);

        let img_base = strip_extension(image_path);
        let img_fmt = image_format.unwrap_or(&ImageFormat::Png);

        for (metric, (ms1_scans, ms1_mat)) in &ms1_results {
            let ext = get_extension(img_fmt);
            let file_path = format!("{}_ms1_{}.{}", img_base, metric, ext);

            // Extract MS1 chromatogram data based on selected type
            let ms1_chromatogram: Vec<f32> = ms1_scans
                .iter()
                .map(|scan_id| {
                    spec_metadata.get(scan_id)
                        .map(|meta| match chromatogram_type {
                            ChromatogramType::BPI => meta.base_peak_intensity,
                            ChromatogramType::TIC => meta.total_ion_current,
                        })
                        .unwrap_or(0.0)
                })
                .collect();

            plot_similarity_heatmap(
                ms1_scans,
                ms1_mat,
                &ms1_chromatogram,
                chromatogram_type,
                enable_smoothing,
                &file_path,
                img_fmt,
                theme.unwrap_or(&ThemeName::Classic.get_theme()),
                global_scan_range,
                None // MS1 uses its own scan IDs directly
            )
            .map_err(|e| CliError::FileError(format!("**Failed to create MS1 {} heatmap: {}**", metric, e)))?;
            println!("||    Created MS1 {} heatmap: {}", metric, file_path);
        }

        for (metric, (ms2_scans, ms2_mat)) in &ms2_results {
            let ext = get_extension(img_fmt);
            let file_path = format!("{}_ms2_{}.{}", img_base, metric, ext);

            // Use the SAME MS1 chromatogram data and scan IDs directly 
            // This ensures identical chromatogram representation across all heatmaps
            let (ms1_scan_ids, ms1_chromatogram): (Vec<String>, Vec<f32>) = ms1_results
                .values()
                .next() // Take any MS1 result (they all have the same scan IDs)
                .map(|(ms1_scans, _)| {
                    let scan_ids = ms1_scans.clone();
                    let chromatogram = ms1_scans
                        .iter()
                        .map(|scan_id| {
                            spec_metadata.get(scan_id)
                                .map(|meta| match chromatogram_type {
                                    ChromatogramType::BPI => meta.base_peak_intensity,
                                    ChromatogramType::TIC => meta.total_ion_current,
                                })
                                .unwrap_or(0.0)
                        })
                        .collect();
                    (scan_ids, chromatogram)
                })
                .unwrap_or_default();

            plot_similarity_heatmap(
                ms2_scans,
                ms2_mat,
                &ms1_chromatogram, // Use MS1 chromatogram data
                chromatogram_type,
                enable_smoothing,
                &file_path,
                img_fmt,
                theme.unwrap_or(&ThemeName::Classic.get_theme()),
                global_scan_range,
                Some(ms1_scan_ids) // Pass MS1 scan IDs for proper chromatogram mapping
            )
            .map_err(|e| CliError::FileError(format!("**Failed to create MS2 {} heatmap: {}**", metric, e)))?;
            println!("||    Created MS2 {} heatmap: {}", metric, file_path);
        }
    }

    Ok(())
}

pub fn print_configuration(
    input_path: &str,
    output_path: &str,
    output_format: &OutputFormat,
    ms1_metrics: &[&str],
    ms2_metrics: &[&str],
    ms1_threshold: ThresholdMethod,
    ms2_threshold: ThresholdMethod,
    mass_tolerance: f32,
    heatmap_enabled: bool,
    image_path: Option<&str>,
    image_format: Option<&ImageFormat>,
    theme: Option<&ColorTheme>,
    chromatogram_type: Option<ChromatogramType>,
    enable_smoothing: Option<bool>,
) {
    println!("----------------------------------------------------------------------------------------------");
    println!(" CONFIGURATION SUMMARY");
    println!("----------------------------------------------------------------------------------------------");
    println!("||    Input file: {}", input_path);
    println!("||    Output file: {}", output_path);
    println!("||    Output format: {}", output_format);
    println!("||    MS1 Similarity metrics: {:?}", ms1_metrics);
    println!("||    MS2 Similarity metrics: {:?}", ms2_metrics);
    println!("||    MS1 intensity threshold: {:?}", ms1_threshold);
    println!("||    MS2 intensity threshold: {:?}", ms2_threshold);
    println!("||    Mass tolerance: {:.4}", mass_tolerance);

    if heatmap_enabled {
        println!("||");
        println!("||    Heatmap generation: ENABLED");
        if let Some(path) = image_path {
            println!("||    Heatmap output path: {}", path);
        }
        if let Some(format) = image_format {
            println!("||    Heatmap image format: {:?}", format);
        }
        if let Some(chrom_type) = chromatogram_type {
            println!("||    Chromatogram type: {}", chrom_type.as_str());
        }
        if let Some(smoothing) = enable_smoothing {
            println!("||    Chromatogram smoothing: {}", if smoothing { "ENABLED" } else { "DISABLED" });
        }
        if let Some(t) = theme {
            println!("||    Color theme:");
            println!("||       - Background: rgb({}, {}, {})", t.background.r, t.background.g, t.background.b);
            println!("||       - Low similarity: rgb({}, {}, {})", t.low.r, t.low.g, t.low.b);
            println!("||       - High similarity: rgb({}, {}, {})", t.high.r, t.high.g, t.high.b);
        }
    } else {
        println!("||    Heatmap generation: DISABLED");
    }

    println!("----------------------------------------------------------------------------------------------");
}

fn confirm_processing() -> Result<bool> {
    print!("||    Continue with processing? [Y/n]: ");
    io::stdout().flush()?;
    
    let mut input = String::new();
    io::stdin().read_line(&mut input)?;
    
    let response = input.trim().to_lowercase();
    Ok(response.is_empty() || response == "y" || response == "yes")
}

fn detect_input_file_type(file_path: &str) -> Option<&str> {
    match file_path.rsplitn(2, '.').next() {
        Some("mzML") | Some("mzml") => Some("mzML"),
        Some("mzXML") | Some("mzxml") => Some("mzXML"),
        _ => None,
    }
}

fn detect_output_file_type(file_path: &str) -> Option<OutputFormat> {
    match file_path.rsplitn(2, '.').next() {
        Some("csv") => Some(OutputFormat::Csv),
        Some("tsv") => Some(OutputFormat::Tsv),
        Some("json") => Some(OutputFormat::Json),
        _ => Some(OutputFormat::Csv), // Default to CSV if no extension
    }
}

/// Implement Display so the format can be printed easily.
impl fmt::Display for OutputFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            OutputFormat::Csv => "csv",
            OutputFormat::Tsv => "tsv",
            OutputFormat::Json => "json",
        };
        write!(f, "{}", s)
    }
}

// Ordered metric lists — index + 1 is the number the user types.
const MS1_METRICS: &[&str] = &[
    "cosine", "weighted-cosine", "spectral-entropy", "hellinger", "bray-curtis",
];
const MS2_METRICS: &[&str] = &[
    "cosine", "weighted-cosine", "spectral-entropy", "hellinger", "bray-curtis",
    "modified-cosine", "modified-weighted-cosine", "modified-spectral-entropy",
    "modified-hellinger", "modified-bray-curtis",
];

fn detect_similarity_metrics(input: &str, ms_level: f32) -> Option<Vec<&'static str>> {
    if input.is_empty() {
        return Some(vec!["cosine"]);
    }

    let numbered = if ms_level == 1.0 { MS1_METRICS } else { MS2_METRICS };

    let mut metrics = Vec::new();
    for token in input.split(',').map(|s| s.trim()) {
        // Accept a bare number (e.g. "3") as an index into the numbered list.
        if let Ok(n) = token.parse::<usize>() {
            if n >= 1 && n <= numbered.len() {
                metrics.push(numbered[n - 1]);
            } else {
                return None; // out-of-range number
            }
            continue;
        }
        // Accept the metric name (or common aliases).
        match token.to_lowercase().as_str() {
            "cosine"                                              => metrics.push("cosine"),
            "weighted-cosine"  | "weighted_cosine"  | "weighted"  => metrics.push("weighted-cosine"),
            "spectral-entropy" | "spectral_entropy" | "entropy"   => metrics.push("spectral-entropy"),
            "hellinger"                                           => metrics.push("hellinger"),
            "bray-curtis"      | "bray_curtis"      | "bray"      => metrics.push("bray-curtis"),
            "modified-cosine"           | "modified_cosine"           | "mod-cosine"    => metrics.push("modified-cosine"),
            "modified-weighted-cosine"  | "modified_weighted_cosine"  | "mod-weighted"  => metrics.push("modified-weighted-cosine"),
            "modified-spectral-entropy" | "modified_spectral_entropy" | "mod-entropy"   => metrics.push("modified-spectral-entropy"),
            "modified-hellinger"        | "modified_hellinger"        | "mod-hellinger" => metrics.push("modified-hellinger"),
            "modified-bray-curtis"      | "modified_bray_curtis"      | "mod-bray"      => metrics.push("modified-bray-curtis"),
            _ => return None,
        }
    }

    if metrics.is_empty() { None } else { Some(metrics) }
}

fn read_input_line() -> Result<String> {
    let mut input = String::new();
    io::stdin().read_line(&mut input)?;
    Ok(input.trim().to_string())
}

fn prompt_input_path() -> Result<String> {
    const MAX_ATTEMPTS: u32 = 3;
    let mut attempts = 0;
    
    loop {
        attempts += 1;
        println!("----------------------------------------------------------------------------------------------");
        print!("||    Please enter the full path to your input file: ");
        io::stdout().flush()?;

        let input = read_input_line()?;
        
        // Handle empty input
        if input.is_empty() {
            clear_current_line();
            println!("||      Path cannot be empty. Please try again.");
            if attempts >= MAX_ATTEMPTS {
                return Err(CliError::InvalidInput("**Too many invalid attempts**".to_string()));
            }
            continue;
        }

        // Trim quotation marks
        let path = input.trim_matches('"').trim_matches('\'');

        // Check if file exists
        if !Path::new(path).exists() {
            clear_current_line();
            println!("||     File does not exist: {}", path);
            if attempts >= MAX_ATTEMPTS {
                return Err(CliError::FileError(format!("**File not found after {} attempts: {}**", MAX_ATTEMPTS, path)));
            }
            continue;
        }

        // Check if it's a file (not a directory)
        if !Path::new(path).is_file() {
            clear_current_line();
            println!("||     Path is not a file: {}", path);
            if attempts >= MAX_ATTEMPTS {
                return Err(CliError::FileError("**Path is not a valid file**".to_string()));
            }
            continue;
        }

        // Try to detect file type
        if let Some(ft) = detect_input_file_type(path) {
            clear_current_line();
            println!("||     File selected: {}", path);
            println!("||     File type detected: {}", ft);
            return Ok(path.to_string());
        } else {
            clear_current_line();
            println!("||     Unsupported file type. Supported formats: .mzML, .mzXML");
            if attempts >= MAX_ATTEMPTS {
                return Err(CliError::InvalidInput("**Unsupported file format after multiple attempts**".to_string()));
            }
        }
    }
}

fn prompt_output_path() -> Result<(String, OutputFormat)> {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        println!("||    Please enter the desired path to your output file.");
        println!("||    Supported output formats: csv, tsv, json.");
        println!("||    Leave blank to default to 'output.csv'.");
        print!("||    Path: ");
        io::stdout().flush()?;

        let input = read_input_line()?;
        
        let path = if input.is_empty() { 
            "output.csv".to_string() 
        } else { 
            input.clone()
        };

        // For the default case (output.csv), skip directory validation since it goes to current directory
        if path == "output.csv" {
            clear_current_line();
            println!("||     Using default output file: {}", path);
            println!("||     Output format: csv");
            return Ok((path, OutputFormat::Csv));
        }

        // Trim quotation marks
        let path = input.trim_matches('"').trim_matches('\'');

        // Check if the output directory exists and is writable for non-default paths
        if let Some(parent) = Path::new(&path).parent() {
            // Only check if parent is not the current directory
            if parent != Path::new("") && parent != Path::new(".") && !parent.exists() {
                clear_current_line();
                println!("||     Output directory does not exist: {:?}", parent);
                continue;
            }
        }

        // Validate format
        if let Some(ft) = detect_output_file_type(&path) {
            clear_current_line();
            println!("||     Output file: {}", path);
            println!("||     Output format: {}", ft);
            return Ok((path.to_string(), ft));
        } else {
            clear_current_line();
            println!("||    Error: Could not determine output format. Using CSV as default.");
            return Ok((format!("{}.csv", path), OutputFormat::Csv));
        }
    }
}

fn prompt_similarity_metrics(ms_level: f32) -> Result<Vec<&'static str>> {
    const DIVIDER: &str = "----------------------------------------------------------------------------------------------";

    loop {
        println!("{}", DIVIDER);
        println!("||    MS{} similarity metrics — enter numbers or names, comma-separated (e.g. 1,3 or cosine,hellinger):", ms_level as u8);
        println!("||      1) cosine           — standard dot-product cosine (default)");
        println!("||      2) weighted-cosine  — mz² × I^0.5 weighted cosine (MassBank standard)");
        println!("||      3) spectral-entropy — information-theoretic; robust to noise/missing peaks");
        println!("||      4) hellinger        — Bhattacharyya coeff on sqrt-prob vectors; ESI-tuned");
        println!("||      5) bray-curtis      — compositional dissimilarity; tolerant to outlier peaks");
        if ms_level == 2.0 {
            println!("||    MS2-only — precursor-mass shift alignment applied before scoring:");
            println!("||      6) modified-cosine");
            println!("||      7) modified-weighted-cosine");
            println!("||      8) modified-spectral-entropy");
            println!("||      9) modified-hellinger");
            println!("||     10) modified-bray-curtis");
        }
        print!("||    Metric(s) [1 / cosine]: ");
        io::stdout().flush()?;

        let input = read_input_line()?;

        match detect_similarity_metrics(&input, ms_level) {
            Some(metrics) if ms_level == 1.0 && metrics.iter().any(|m| m.starts_with("modified-")) => {
                clear_current_line();
                println!("||    modified-* metrics are only available for MS2 spectra.");
            }
            Some(metrics) => {
                clear_current_line();
                println!("||     MS{} similarity metric(s): {:?}", ms_level as u8, metrics);
                return Ok(metrics);
            }
            None => {
                clear_current_line();
                println!("||     Unrecognised input. Enter numbers (e.g. 1,3) or names (e.g. cosine,hellinger).");
            }
        }
    }
}

fn clear_current_line() {
    print!("\r\x1B[K");
}

fn prompt_threshold_method(is_ms1: bool) -> Result<ThresholdMethod> {
    let level = if is_ms1 { "MS1" } else { "MS2" };
    loop {
        println!("----------------------------------------------------------------------------------------------");
        println!("||    {} intensity threshold — enter a number or name:", level);
        println!("||      1) absolute (abs)        — fixed intensity floor (default)");
        println!("||      2) percent  (pbp)         — keep peaks >= X% of scan's base peak");
        println!("||      3) topn     (top-n)        — keep N most intense peaks per scan");
        print!("||    Choice [1 / absolute]: ");
        io::stdout().flush()?;

        let input = read_input_line()?;
        let raw = if input.is_empty() { "1".to_string() } else { input.to_lowercase() };

        // Normalise both numbers and name aliases to a canonical key.
        let choice = match raw.trim() {
            "1" | "absolute" | "abs"                                   => "absolute",
            "2" | "percent"  | "pbp" | "percentbasepeak"
                | "percent-base-peak" | "percent_base_peak"            => "percent",
            "3" | "topn"     | "top-n" | "top_n" | "top"              => "topn",
            _ => "unknown",
        };

        match choice {
            "absolute" => {
                let default = if is_ms1 { 1.0_f32 } else { 0.0_f32 };
                print!("||    Absolute floor (default: {:.1}): ", default);
                io::stdout().flush()?;
                let val_input = read_input_line()?;
                let val_str = if val_input.is_empty() { default.to_string() } else { val_input };
                match val_str.parse::<f32>() {
                    Ok(v) if v >= 0.0 => {
                        clear_current_line();
                        println!("||     {} threshold: Absolute({:.2})", level, v);
                        return Ok(ThresholdMethod::Absolute(v));
                    }
                    _ => { clear_current_line(); println!("||     Invalid value — must be >= 0."); }
                }
            }
            "percent" => {
                print!("||    Fraction of base peak (0.0–1.0, default: 0.01): ");
                io::stdout().flush()?;
                let val_input = read_input_line()?;
                let val_str = if val_input.is_empty() { "0.01".to_string() } else { val_input };
                match val_str.parse::<f32>() {
                    Ok(v) if v > 0.0 && v <= 1.0 => {
                        clear_current_line();
                        println!("||     {} threshold: PercentBasePeak({:.4})", level, v);
                        return Ok(ThresholdMethod::PercentBasePeak(v));
                    }
                    _ => { clear_current_line(); println!("||     Invalid fraction — must be in (0.0, 1.0]."); }
                }
            }
            "topn" => {
                print!("||    Top N peaks to keep per scan (default: 20): ");
                io::stdout().flush()?;
                let val_input = read_input_line()?;
                let val_str = if val_input.is_empty() { "20".to_string() } else { val_input };
                match val_str.parse::<usize>() {
                    Ok(n) if n >= 1 => {
                        clear_current_line();
                        println!("||     {} threshold: TopN({})", level, n);
                        return Ok(ThresholdMethod::TopN(n));
                    }
                    _ => { clear_current_line(); println!("||     N must be >= 1."); }
                }
            }
            _ => {
                clear_current_line();
                println!("||     Unrecognised input. Enter 1/absolute/abs, 2/percent/pbp, or 3/topn/top-n.");
            }
        }
    }
}

fn prompt_mass_tolerance() -> Result<f32> {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        print!("||    Mass tolerance (default: 0.01): ");
        io::stdout().flush()?;

        let input = read_input_line()?;
        let input_str = if input.is_empty() { "0.01" } else { &input };

        match input_str.parse::<f32>() {
            Ok(value) if value > 0.0 => {
                clear_current_line();
                println!("||     Mass tolerance: {:.4}", value);
                return Ok(value);
            }
            Ok(_) => {
                clear_current_line();
                println!("||     Mass tolerance must be positive");
            }
            Err(_) => {
                clear_current_line();
                println!("||     Invalid number format");
            }
        }
    }
}

fn prompt_generate_heatmap() -> Result<bool> {
    println!("----------------------------------------------------------------------------------------------");
    print!("||    Generate similarity matrix heatmap? [Y/n]: ");
    io::stdout().flush()?;

    let response = read_input_line()?.to_lowercase();
    Ok(response.is_empty() || response == "y" || response == "yes")
}

fn prompt_image_format() -> Result<ImageFormat> {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        println!("||    Output image format: png, svg, or jpeg");
        print!("||    Format (default: png): ");
        io::stdout().flush()?;

        let input = read_input_line()?;
        let format_str = if input.is_empty() { "png" } else { input.as_str() };

        if let Some(fmt) = ImageFormat::from_ext(format_str) {
            println!("||     Output format set to: {}", format_str.to_lowercase());
            return Ok(fmt);
        } else {
            println!("||     Invalid format. Try: png, svg, jpeg.");
        }
    }
}

fn prompt_color_theme() -> Result<ColorTheme> {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        println!("||    Choose a color theme:");
        for name in ThemeName::list() {
            println!("||      - {name}");
        }

        print!("||    Theme (default: classic): ");
        io::stdout().flush()?;

        let input = read_input_line()?;
        let name = if input.is_empty() { "classic" } else { &input };

        if let Some(theme) = ThemeName::from_str(name) {
            println!("||     Selected theme: {}", name.to_lowercase());
            return Ok(theme.get_theme());
        } else {
            println!("||     Invalid theme. Please try again.");
        }
    }
}

fn prompt_output_image_path() -> Result<String> {
    println!("----------------------------------------------------------------------------------------------");
    println!("||    Path to save the heatmap image (default: similarity_heatmap.png): ");
    print!("||    Path: ");
    io::stdout().flush()?;

    let input = read_input_line()?;
    Ok(if input.is_empty() {
        "similarity_heatmap.png".to_string()
    } else {
        input
    })
}

fn prompt_chromatogram_type() -> Result<ChromatogramType> {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        println!("||    Select chromatogram data type to display on heatmap axes:");
        println!("||    Available options: {}", ChromatogramType::list().join(", "));
        println!("||    (BPI = Base Peak Intensity, TIC = Total Ion Current)");
        print!("||    Chromatogram type (default: BPI): ");
        io::stdout().flush()?;

        let input = read_input_line()?;
        let choice = if input.is_empty() { "BPI" } else { &input };

        if let Some(chrom_type) = ChromatogramType::from_str(choice) {
            println!("||     Selected chromatogram type: {}", chrom_type.as_str());
            return Ok(chrom_type);
        } else {
            println!("||     Invalid selection. Please choose from: {}", ChromatogramType::list().join(", "));
        }
    }
}

fn prompt_pre_similarity_prune_fraction() -> Result<f64> {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        println!("||    Pre-similarity background pruning: remove m/z bins present in too many spectra.");
        print!("||    Fraction threshold (0.0-1.0, default: 0.3): ");
        io::stdout().flush()?;

        let input = read_input_line()?;
        let val_str = if input.is_empty() { "0.3" } else { &input };

        match val_str.parse::<f64>() {
            Ok(val) if val >= 0.0 && val <= 1.0 => {
                clear_current_line();
                println!("||     Pre-similarity prune fraction: {:.2}", val);
                return Ok(val);
            }
            _ => {
                clear_current_line();
                println!("||     Invalid value. Please enter a number between 0.0 and 1.0.");
            }
        }
    }
}

fn prompt_enable_background_pruning() -> Result<(bool, f64, f64)> {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        println!("||    Enable post-similarity background pruning?");
        println!("||    This removes background noise regions that have high similarity to many other spectra.");
        print!("||    Enable background pruning? (y/n, default: n): ");
        io::stdout().flush()?;

        let input = read_input_line()?;
        let choice = if input.is_empty() { "n" } else { &input };

        match choice.trim().to_lowercase().as_str() {
            "y" | "yes" | "true" | "1" => {
                println!("||     Background pruning enabled");
                
                // Get similarity threshold
                let similarity_threshold = loop {
                    print!("||    Similarity threshold for background detection (0.0-1.0, default: 0.7): ");
                    io::stdout().flush()?;
                    let input = read_input_line()?;
                    let threshold_str = if input.is_empty() { "0.7" } else { &input };
                    
                    match threshold_str.parse::<f64>() {
                        Ok(val) if val >= 0.0 && val <= 1.0 => {
                            println!("||     Similarity threshold set to: {}", val);
                            break val;
                        }
                        _ => println!("||     Invalid threshold. Please enter a value between 0.0 and 1.0."),
                    }
                };
                
                // Get minimum background fraction  
                let min_background_fraction = loop {
                    print!("||    Minimum fraction of high similarities for background (0.0-1.0, default: 0.8): ");
                    io::stdout().flush()?;
                    let input = read_input_line()?;
                    let fraction_str = if input.is_empty() { "0.8" } else { &input };
                    
                    match fraction_str.parse::<f64>() {
                        Ok(val) if val >= 0.0 && val <= 1.0 => {
                            println!("||     Background fraction set to: {}", val);
                            break val;
                        }
                        _ => println!("||     Invalid fraction. Please enter a value between 0.0 and 1.0."),
                    }
                };
                
                return Ok((true, similarity_threshold, min_background_fraction));
            }
            "n" | "no" | "false" | "0" => {
                println!("||     Background pruning disabled");
                return Ok((false, 0.7, 0.8)); // Default values, not used
            }
            _ => {
                println!("||     Invalid input. Please enter 'y' for yes or 'n' for no.");
            }
        }
    }
}

fn prompt_enable_smoothing() -> Result<bool> {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        println!("||    Enable chromatogram smoothing to reduce noise?");
        println!("||    Smoothing uses a 3-point moving average to create cleaner chromatograms.");
        print!("||    Enable smoothing? (y/n, default: y): ");
        io::stdout().flush()?;

        let input = read_input_line()?;
        let choice = if input.is_empty() { "y" } else { &input };

        match choice.trim().to_lowercase().as_str() {
            "y" | "yes" | "true" | "1" => {
                println!("||     Smoothing enabled");
                return Ok(true);
            }
            "n" | "no" | "false" | "0" => {
                println!("||     Smoothing disabled");
                return Ok(false);
            }
            _ => {
                println!("||     Invalid input. Please enter 'y' for yes or 'n' for no.");
            }
        }
    }
}

/// Strip a known file extension from a path, returning the base.
fn strip_extension(path: &str) -> &str {
    for ext in &[".csv", ".tsv", ".json", ".png", ".svg", ".jpeg", ".jpg"] {
        if let Some(base) = path.strip_suffix(ext) {
            return base;
        }
    }
    path
}

fn get_extension(format: &ImageFormat) -> &'static str {
    match format {
        ImageFormat::Png => "png",
        ImageFormat::Svg => "svg",
        ImageFormat::Jpeg => "jpeg",
    }
}

#[derive(serde::Serialize)]
struct Config {
    input_path: String,
    output_path: String,
    output_format: OutputFormat,
    ms1_similarity_metrics: Vec<&'static str>,
    ms2_similarity_metrics: Vec<&'static str>,
    ms1_threshold: ThresholdMethod,
    ms2_threshold: ThresholdMethod,
    mass_tolerance: f32,
    heatmap_enabled: bool,
    image_path: Option<String>,
    image_format: Option<ImageFormat>,
    theme: Option<ColorTheme>,
}

fn write_config_to_file(config: &Config) -> Result<()> {
    let mut file = File::create("config.toml")?;

    let toml_str = toml::ser::to_string(config)
        .map_err(|e| CliError::FileError(format!("Failed to serialize config to TOML: {}", e)))?;
    use std::io::Write;
    file.write_all(toml_str.as_bytes())?;

    Ok(())
}