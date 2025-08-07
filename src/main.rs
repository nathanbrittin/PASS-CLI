// use std::collections::HashMap;
use hashbrown::HashMap;
use ndarray::Array2;
use std::process;
use std::fmt; // Import fmt for Display implementation
use std::path::Path;
use std::io::{self, Write}; // Import io and Write traits for flushing output
use std::time::Instant;
use std::error::Error;

mod ms_io;
pub use ms_io::{
    import_mzml,
    OutputFormat, 
    write_similarity_matrix, 
    filter_by_ms_level
};
mod similarity;
pub use similarity::{
    compute_pairwise_similarity_matrix_sparse, 
    prune_background_columns, 
    prune_background_bins_sparse, 
    prune_unused_bins, 
    spectrum_to_dense_vec, 
    compute_sparse_vec_map, 
    compute_dense_vec_map, 
    cosine_similarity, 
    compute_pairwise_similarity_matrix_ndarray
};

mod visual;
use visual::{plot_similarity_heatmap, ImageFormat, ColorTheme, ThemeName};

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
||  -Supports file formats: `.mzML` (default) and `.mzXML` (future)                                          
----------------------------------------------------------------------------------------------
||  Please follow the directions to process your data.                                      
||  Alternatively, provide a config.yaml file to automate processing.                       
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
    let ms1_minimum_intensity = prompt_min_intensity(true)?;
    let ms2_minimum_intensity = prompt_min_intensity(false)?;
    let noise_threshold = prompt_noise_threshold()?;
    let mass_tolerance = prompt_mass_tolerance()?;

    // Visuals config
    let heatmap_enabled = prompt_generate_heatmap()?;

    // If heatmap is enabled, get all config at once
    let (image_path, image_format, theme) = if heatmap_enabled {
        let path = prompt_output_image_path()?;
        let format = prompt_image_format()?;
        let theme = prompt_color_theme()?;
        (path, format, theme)
    } else {
        ("".to_string(), ImageFormat::Png, ThemeName::Classic.get_theme())
    };

    // Print a final confirmation
    print_configuration(&input_path, 
        &output_path, 
        &output_format, 
        &ms1_similarity_metrics, 
        &ms2_similarity_metrics,
        ms1_minimum_intensity,
        ms2_minimum_intensity,
        noise_threshold,
        mass_tolerance,
        heatmap_enabled,
        Some(&image_path),
        Some(&image_format),
        Some(&theme)
    );

    if !confirm_processing()? {
        println!("**Processing cancelled by user.**");
        return Ok(());
    }

    let start = Instant::now();
    
    // Process the data with proper error handling
    process_spectral_data(
        &input_path, &output_path, output_format,
        &ms1_similarity_metrics, &ms2_similarity_metrics,
        ms1_minimum_intensity, ms2_minimum_intensity,
        noise_threshold, mass_tolerance,
        heatmap_enabled, Some(&image_format), Some(&theme),
    )?;

    let duration = start.elapsed();
    println!("||    Success! Processing completed in {duration:.2?}");
    Ok(())
}

fn process_spectral_data(
    input_path: &str,
    output_path: &str,
    output_format: OutputFormat,
    ms1_similarity_metrics: &[&str],
    ms2_similarity_metrics: &[&str],
    ms1_minimum_intensity: f32,
    ms2_minimum_intensity: f32,
    _noise_threshold: u32,
    mass_tolerance: f32,
    heatmap_enabled: bool,
    image_format: Option<&ImageFormat>,
    theme: Option<&ColorTheme>,
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

    // Filter to only MS1 and MS2
    let ms1_spec_map = filter_by_ms_level(spec_map.clone(), spec_metadata.clone(), 1);
    let ms2_spec_map = filter_by_ms_level(spec_map.clone(), spec_metadata.clone(), 2);
    println!("||    Num. MS1 spectra: {}, Num. MS2 spectra: {}", ms1_spec_map.len(), ms2_spec_map.len());


    // Determine a max_mz for binning (e.g., highest m/z across all spectra)
    let max_ms1_mz = ms1_spec_map.values()
        .map(|p| p.iter().map(|p| p.mz).fold(0., f32::max))
        .fold(0., f32::max);

    let max_ms2_mz = ms2_spec_map.values()
        .map(|p| p.iter().map(|p| p.mz).fold(0., f32::max))
        .fold(0., f32::max);

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
    let mut ms1_bits_map = compute_sparse_vec_map(&ms1_spec_map, mass_tolerance, max_ms1_mz, ms1_minimum_intensity)
        .map_err(|e| CliError::ProcessingError(format!("**Failed to compute MS1 sparse vectors: {:?}**", e)))?;
    let mut ms2_bits_map = compute_sparse_vec_map(&ms2_spec_map, mass_tolerance, max_ms2_mz, ms2_minimum_intensity)
        .map_err(|e| CliError::ProcessingError(format!("**Failed to compute MS2 sparse vectors: {:?}**", e)))?;
    println!("||    ✔ Computed sparse binary vectors successfully.");
    println!("------------------------------------------------------------------------------");

    // Filter background signals
    println!("||    Pruning background signals...");
    prune_background_bins_sparse(&mut ms1_bits_map, 0.5)
        .map_err(|e| CliError::ProcessingError(format!("**Failed to prune MS1 background bins: {:?}**", e)))?;
    prune_background_bins_sparse(&mut ms2_bits_map, 0.5)
        .map_err(|e| CliError::ProcessingError(format!("**Failed to prune MS2 background bins: {:?}**", e)))?;
    println!("||    ✔ Pruned background signals successfully.");
    println!("------------------------------------------------------------------------------");

    // Compute pairwise cosine similarities
    let mut ms1_results: HashMap<String, (Vec<String>, Array2<f32>)> = HashMap::new();
    let mut ms2_results: HashMap<String, (Vec<String>, Array2<f32>)> = HashMap::new();

    // Process MS1 metrics
    for metric in ms1_similarity_metrics {
        println!("||    Computing MS1 pairwise similarity matrix using {}...", metric);
        let (ms1_scans, ms1_mat) = match *metric {
            "cosine" => compute_pairwise_similarity_matrix_sparse(&ms1_bits_map, &spec_metadata, metric.to_string(), 1, mass_tolerance)
                .map_err(|e| CliError::ProcessingError(format!("**Failed to compute MS1 similarity matrix: {:?}**", e)))?,
            _ => return Err(CliError::ProcessingError(format!("**Unsupported MS1 similarity metric: {}**", metric))),
        };
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
        let (ms2_scans, ms2_mat) = match *metric {
            "cosine" => compute_pairwise_similarity_matrix_sparse(&ms2_bits_map, &spec_metadata, metric.to_string(), 2, mass_tolerance)
                .map_err(|e| CliError::ProcessingError(format!("**Failed to compute MS2 similarity matrix: {:?}**", e)))?,
            "modified-cosine" => compute_pairwise_similarity_matrix_sparse(&ms2_bits_map, &spec_metadata, metric.to_string(), 2, mass_tolerance)
                .map_err(|e| CliError::ProcessingError(format!("**Failed to compute MS2 similarity matrix: {:?}**", e)))?,
            _ => return Err(CliError::ProcessingError(format!("**Unsupported MS2 similarity metric: {}**", metric))),
        };
        ms2_results.insert(metric.to_string(), (ms2_scans, ms2_mat));
    }
    println!("||    ✔ Computed MS2 pairwise similarity matrices successfully.");
    println!("------------------------------------------------------------------------------");

    // Export results
    println!("||    Exporting similarity matrices...");
    
    for (metric, (ms1_scans, ms1_mat)) in &ms1_results {
        let ms1_output_path = output_path.replace(".csv", &format!("_ms1_{}.csv", metric));
        write_similarity_matrix(ms1_scans, ms1_mat, &ms1_output_path, output_format)
            .map_err(|e| CliError::FileError(format!("**Failed to write MS1 {} matrix: {}**", metric, e)))?;
        println!("||    Exported MS1 {} matrix to: {}", metric, ms1_output_path);
    }

    for (metric, (ms2_scans, ms2_mat)) in &ms2_results {
        let ms2_output_path = output_path.replace(".csv", &format!("_ms2_{}.csv", metric));
        write_similarity_matrix(ms2_scans, ms2_mat, &ms2_output_path, output_format)
            .map_err(|e| CliError::FileError(format!("**Failed to write MS2 {} matrix: {}**", metric, e)))?;
        println!("||    Exported MS2 {} matrix to: {}", metric, ms2_output_path);
    }

    if heatmap_enabled {
        println!("||    Creating heatmaps...");

        for (metric, (ms1_scans, ms1_mat)) in &ms1_results {
            let ext = get_extension(image_format.unwrap_or(&ImageFormat::Png));
            let file_path = output_path.replace(".csv", &format!("_ms1_{}_heatmap.{}", metric, ext));

            plot_similarity_heatmap(
                ms1_scans,
                ms1_mat,
                &file_path,
                image_format.unwrap_or(&ImageFormat::Png),
                theme.unwrap_or(&ThemeName::Classic.get_theme())
            )
            .map_err(|e| CliError::FileError(format!("**Failed to create MS1 {} heatmap: {}**", metric, e)))?;
            println!("||    Created MS1 {} heatmap: {}", metric, file_path);
        }

        for (metric, (ms2_scans, ms2_mat)) in &ms2_results {
            let ext = get_extension(image_format.unwrap_or(&ImageFormat::Png));
            let file_path = output_path.replace(".csv", &format!("_ms2_{}_heatmap.{}", metric, ext));

            plot_similarity_heatmap(
                ms2_scans,
                ms2_mat,
                &file_path,
                image_format.unwrap_or(&ImageFormat::Png),
                theme.unwrap_or(&ThemeName::Classic.get_theme())
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
    ms1_min_intensity: f32,
    ms2_min_intensity: f32,
    noise_threshold: u32,
    mass_tolerance: f32,
    heatmap_enabled: bool,
    image_path: Option<&str>,
    image_format: Option<&ImageFormat>,
    theme: Option<&ColorTheme>,
) {
    println!("----------------------------------------------------------------------------------------------");
    println!(" CONFIGURATION SUMMARY");
    println!("----------------------------------------------------------------------------------------------");
    println!("||    Input file: {}", input_path);
    println!("||    Output file: {}", output_path);
    println!("||    Output format: {}", output_format);
    println!("||    MS1 Similarity metrics: {:?}", ms1_metrics);
    println!("||    MS2 Similarity metrics: {:?}", ms2_metrics);
    println!("||    MS1 Minimum intensity threshold: {:.2}", ms1_min_intensity);
    println!("||    MS2 Minimum intensity threshold: {:.2}", ms2_min_intensity);
    println!("||    Noise threshold: {}", noise_threshold);
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
        if let Some(t) = theme {
            println!("||    Color theme:");
            println!("||       - Background: rgb({}, {}, {})", t.background.0, t.background.1, t.background.2);
            println!("||       - Low similarity: rgb({}, {}, {})", t.low.0, t.low.1, t.low.2);
            println!("||       - High similarity: rgb({}, {}, {})", t.high.0, t.high.1, t.high.2);
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

fn detect_similarity_metrics(input: &str) -> Option<Vec<&'static str>> {
    if input.is_empty() {
        return Some(vec!["cosine"]);  // Default case
    }
    
    let mut metrics = Vec::new();
    let inputs: Vec<&str> = input.split(',').map(|s| s.trim()).collect();
    
    for metric in inputs {
        match metric.to_lowercase().as_str() {
            "cosine" => metrics.push("cosine"),
            "modified-cosine" | "modified_cosine" => metrics.push("modified-cosine"),
            _ => return None,  // Invalid metric found
        }
    }
    
    if metrics.is_empty() {
        None
    } else {
        Some(metrics)
    }
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
            input
        };

        // For the default case (output.csv), skip directory validation since it goes to current directory
        if path == "output.csv" {
            clear_current_line();
            println!("||     Using default output file: {}", path);
            println!("||     Output format: csv");
            return Ok((path, OutputFormat::Csv));
        }

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
            return Ok((path, ft));
        } else {
            clear_current_line();
            println!("||    Error: Could not determine output format. Using CSV as default.");
            return Ok((format!("{}.csv", path), OutputFormat::Csv));
        }
    }
}

fn prompt_similarity_metrics(ms_level: f32) -> Result<Vec<&'static str>> {
    const DIVIDER: &str = "----------------------------------------------------------------------------------------------";
    const PROMPT_MS1: &str = "||    MS1 similarity metric(s) [cosine(default)]: ";
    const PROMPT_MS2: &str = "||    MS2 similarity metric(s) [cosine(default), modified-cosine]: ";
    
    loop {
        println!("{}", DIVIDER);
        let prompt = if ms_level == 1.0 { PROMPT_MS1 } else { PROMPT_MS2 };
        print!("{}", prompt);
        io::stdout().flush()?;

        let input = read_input_line()?;
        
        match detect_similarity_metrics(&input) {
            // If the ms_level is 1, do not allow modified-cosine
            Some(metrics) if ms_level == 1.0 && metrics.contains(&"modified-cosine") => {
                clear_current_line();
                println!("||    modified-cosine is not supported for MS1 spectra");
            }
            Some(metrics) => {
                clear_current_line();
                if ms_level == 1.0 {
                    println!("||     MS1 similarity metric(s): {:?}", metrics);
                } else {
                    println!("||     MS2 similarity metric(s): {:?}", metrics);
                }
                return Ok(metrics);
            }
            None => {
                clear_current_line();
                println!("||     Invalid metrics. Supported: cosine{}",
                        if ms_level == 2.0 { ", modified-cosine" } else { "" });
            }
        }
    }
}

fn clear_current_line() {
    print!("\r\x1B[K");
}

fn prompt_min_intensity(is_ms1: bool) -> Result<f32> {
    let default_value = if is_ms1 { 1.0 } else { 0.0 };
    
    loop {
        println!("----------------------------------------------------------------------------------------------");
        if is_ms1 {
            print!("||    MS1 minimum intensity threshold (default: {:.1}): ", default_value);
        } else {
            print!("||    MS2 minimum intensity threshold (default: {:.1}): ", default_value);
        }
        io::stdout().flush()?;

        let input = read_input_line()?;
        let input_str = if input.is_empty() {
            default_value.to_string()
        } else {
            input
        };

        match input_str.parse::<f32>() {
            Ok(value) if value >= 0.0 => {
                clear_current_line();
                println!("||     Minimum intensity threshold: {:.2}", value);
                return Ok(value);
            }
            Ok(_) => {
                clear_current_line();
                println!("||     Intensity must be non-negative");
            }
            Err(_) => {
                clear_current_line();
                println!("||     Invalid number format");
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

fn prompt_noise_threshold() -> Result<u32> {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        print!("||    Noise threshold (default: 100): ");
        io::stdout().flush()?;

        let input = read_input_line()?;
        let input_str = if input.is_empty() { "100" } else { &input };

        match input_str.parse::<u32>() {
            Ok(value) => {
                clear_current_line();
                println!("||     Noise threshold: {}", value);
                return Ok(value);
            }
            Err(_) => {
                clear_current_line();
                println!("||     Invalid number format. Please enter a positive integer");
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

fn get_extension(format: &ImageFormat) -> &'static str {
    match format {
        ImageFormat::Png => "png",
        ImageFormat::Svg => "svg",
        ImageFormat::Jpeg => "jpeg",
    }
}