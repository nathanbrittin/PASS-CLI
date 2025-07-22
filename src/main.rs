// use std::collections::HashMap;
use std::process;
use std::fmt; // Import fmt for Display implementation
// use std::path::Path;
use std::io::{self, Write}; // Import io and Write traits for flushing output
use std::time::Instant;
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

    let input_path = prompt_input_path();

    let (output_path, output_format) = prompt_output_path();

    let similarity_metric = prompt_similarity_metric();

    let ms1_minimum_intensity = prompt_min_intensity(true);

    let ms2_minimum_intensity = prompt_min_intensity(false);

    let noise_threshold = prompt_noise_threshold();

    // let max_peaks = prompt_max_peaks();

    let mass_tolerance = prompt_mass_tolerance();

    // let verbose = prompt_verbose();

    // Print a final confirmation
    println!("----------------------------------------------------------------------------------------------");
    println!("||    Input file: {input_path}");
    println!("||    Output file: {output_path}");
    println!("||    Similarity metric: {similarity_metric}");
    println!("||    MS1 Minimum intensity threshold: {ms1_minimum_intensity:.2}");
    println!("||    MS2 Minimum intensity threshold: {ms2_minimum_intensity:.2}");
    println!("||    Noise threshold: {noise_threshold:.2}");
    println!("||    Mass tolerance: {mass_tolerance:.2}");
    // println!("||    Maximum number of peaks: {}", max_peaks);
    // println!("||    Verbose flag: {}", verbose);

    // Wait for the user to press enter
    println!("----------------------------------------------------------------------------------------------");
    println!("||    Press enter to continue...");
    let mut input = String::new();
    io::stdin()
        .read_line(&mut input)
        .expect("Failed to read line");

    let start = Instant::now();
    // Importing the data file
    let (spec_map, spec_metadata) = import_mzml(&input_path).expect("Failed to import mzml file. Check formatting.");
    println!("Parsed {} spectra successfully.", spec_map.len());

    // Filter to only MS1 and MS2
    let ms1_spec_map = filter_by_ms_level(spec_map.clone(), spec_metadata.clone(), 1);
    let ms2_spec_map = filter_by_ms_level(spec_map.clone(), spec_metadata.clone(), 2);
    println!("Number of MS1 spectra: {}", ms1_spec_map.len());
    println!("Number of MS2 spectra: {}", ms2_spec_map.len());

    // Print the first MS1 spectrum
    println!("First MS1 Spectrum:");
    println!("{:?}", ms1_spec_map.values().next().unwrap()[0]);

    // Print the first MS2 spectrum
    println!("First MS2 Spectrum:");
    println!("{:?}", ms2_spec_map.values().next().unwrap()[0]);

    // Determine a max_mz for binning (e.g., highest m/z across all spectra)
    let max_ms1_mz = ms1_spec_map.values().map(|p| p.iter().map(|p| p.mz).fold(0., f32::max)).fold(0., f32::max);
    // Get vector length for binning
    let vector_length = (max_ms1_mz / mass_tolerance).ceil() as usize;
    println!("Min m/z: 0.0, Max m/z: {max_ms1_mz}, Bin width: {mass_tolerance:.2}, Vector Size: {vector_length}");

    // Compute dense binary vectors
    println!("Computing MS1 spare binary vectors...");
    let mut ms1_bits_map = compute_sparse_vec_map(&ms1_spec_map, mass_tolerance, max_ms1_mz, ms1_minimum_intensity);
    println!("Computing MS2 spare binary vectors...");
    let mut ms2_bits_map = compute_sparse_vec_map(&ms2_spec_map, mass_tolerance, max_ms1_mz, ms2_minimum_intensity);
    println!("Computed dense binary vectors successfully.");

    // Filter background signals
    println!("Pruning background signals...");
    prune_background_bins_sparse(&mut ms1_bits_map, 0.5);
    prune_background_bins_sparse(&mut ms2_bits_map, 0.5);
    println!("Pruned background signals successfully.");

    // // Print vector size
    // println!("Vector size: {}", ms1_bits_map.values().next().unwrap().len());

    // // Prune bits_map
    // println!("Pruning unused bins...");
    // let ms1_bits_map = prune_unused_bins(&ms1_bits_map);
    // let ms2_bits_map = prune_unused_bins(&ms2_bits_map);
    // println!("Pruned vector size: {}", ms1_bits_map.values().next().unwrap().len());  

    // Compute pairwise cosine similarities
    println!("Computing MS1 pairwise similarity matrix...");
    // let (ms1_scans, ms1_mat) = compute_pairwise_similarity_matrix_ndarray(&ms1_bits_map);
    let (ms1_scans, ms1_mat) = compute_pairwise_similarity_matrix_sparse(&ms1_bits_map);
    println!("Computing MS2 pairwise similarity matrix...");
    // let (ms2_scans, ms2_mat) = compute_pairwise_similarity_matrix_ndarray(&ms2_bits_map);
    let (ms2_scans, ms2_mat) = compute_pairwise_similarity_matrix_sparse(&ms2_bits_map);
    println!("Computed pairwise similarity matrix successfully.");

    println!("Exporting similarity matricies...");
    let ms1_output_path = output_path.replace(".csv", "_ms1.csv");
    let ms2_output_path = output_path.replace(".csv", "_ms2.csv");
    let _ = write_similarity_matrix(&ms1_scans, &ms1_mat, &ms1_output_path, output_format);
    let _ = write_similarity_matrix(&ms2_scans, &ms2_mat, &ms2_output_path, output_format);
    println!("Exported MS1 similarity matrix successfully.");

    println!("Success! File written to: {output_path}");
    let duration = start.elapsed();
    println!("Time elapsed: {duration:.2?}");
    process::exit(0);

}

fn detect_input_file_type(file_path: &str) -> Option<&str> {
    match file_path.split('.').next_back() {
        Some("mzML") | Some("mzml") => Some("mzML"),
        Some("mzXML") | Some("mzxml") => Some("mzXML"),
        _ => None,
    }
}

fn detect_output_file_type(file_path: &str) -> Option<OutputFormat> {
    match file_path.split('.').next_back() {
        Some("") => Some(OutputFormat::Csv),
        Some("csv") => Some(OutputFormat::Csv),
        Some("tsv") => Some(OutputFormat::Tsv),
        Some("json") => Some(OutputFormat::Json),
        _ => None,
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
        write!(f, "{s}")
    }
}

fn detect_similarity_metric(similarity_metric: &str) -> Option<&str> {
    // If blank, default to cosine
    match similarity_metric {
        "" => Some("cosine"),
        "cosine" => Some("cosine"),
        "modified-cosine" => Some("modified-cosine"),
        _ => None,
    }
}

fn prompt_input_path() -> String {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        print!("||    Please enter the full path to your input file: ");
        io::stdout().flush().unwrap();

        // Read into a temporary buffer
        let mut buffer = String::new();
        io::stdin()
            .read_line(&mut buffer)
            .expect("Failed to read path.");

        // Trim whitespace/newlines
        let trimmed = buffer.trim();

        // Trim quotation marks
        let trimmed = trimmed.trim_matches('"');
        let trimmed = trimmed.trim_matches('\'');

        // Try to detect file type
        if let Some(ft) = detect_input_file_type(trimmed) {
            // Clear the current prompt line and confirm selection
            print!("\r\x1B[K");
            println!("||    File selected: {trimmed}");
            println!("||    File type detected: {ft}");
            return trimmed.to_string();
        } else {
            // Clear the prompt and show an error, then loop again
            print!("\r\x1B[K");
            println!("||    !!! File type not detected. Please verify the file path. !!!");
        }
    };
}

fn prompt_output_path() -> (String, OutputFormat) {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        println!("||    Please enter the desired path to your output file.");
        println!("||    Supported output formats are: csv, tsv, json.");
        println!("||    Leave blank to default to 'output.csv'.");
        print!("||    Path: ");
        io::stdout().flush().unwrap();

        // Read into a fresh buffer
        let mut buffer = String::new();
        io::stdin()
            .read_line(&mut buffer)
            .expect("Failed to read path.");

        // Trim and apply default
        let trimmed = {
            let t = buffer.trim();
            if t.is_empty() { "output.csv" } else { t }
        };

        // Validate format
        if let Some(ft) = detect_output_file_type(trimmed) {
            // Clear the prompt line and confirm
            print!("\r\x1B[K");
            println!("||    File selected: {trimmed}");
            println!("||    File type detected: {ft}");
            return (trimmed.to_string(), ft);
        } else {
            // Clear and show error, then loop again
            print!("\r\x1B[K");
            println!("||    !!! Unsupported format. Please verify the file path. !!!");
        }
    }
}

fn prompt_similarity_metric() -> String {
    loop {
        // Print the divider and context lines
        println!("----------------------------------------------------------------------------------------------");
        // Use `print!` for the input prompt, then flush so it appears immediately
        print!("||    Please enter the desired similarity metric [cosine(default) or modified-cosine]: ");
        let mut similarity_metric = String::new();
        io::stdout().flush().unwrap();

        // Clear buffer and read the user’s response
        similarity_metric.clear();
        io::stdin()
            .read_line(&mut similarity_metric)
            .expect("Failed to read similarity metric.");

        // Trim and detect
        let metric_input = similarity_metric.trim();
        let metric = detect_similarity_metric(metric_input);

        if let Some(m) = metric {
            // Erase the prompt line, then print the selected metric
            print!("\r\x1B[K");
            println!("||    Similarity metric selected: {m}");
            return m.to_string();
        } else {
            // Erase the prompt line, then print an error
            print!("\r\x1B[K");
            println!("||    !!! Similarity metric not detected. Please verify your input. !!!");
            // loop back to re-prompt
        }
    }
}

fn prompt_min_intensity(is_ms1: bool) -> f32 {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        // Prompt for minimum intensity threshold
        if is_ms1 {
            print!("||    Please enter the desired MS1 minimum intensity threshold (default: 1.0): ");
        } else {
            print!("||    Please enter the desired MS2 minimum intensity threshold (default: 0.0): ");
        }
        let mut min_int_input = String::new();
        io::stdout().flush().unwrap();

        // Read into our buffer
        min_int_input.clear();
        io::stdin()
            .read_line(&mut min_int_input)
            .expect("Failed to read minimum intensity threshold.");

        // Trim & apply default
        let trimmed_min_int = min_int_input.trim();
        let input_min_int_str = if trimmed_min_int.is_empty() {
            "1.0"
        } else {
            trimmed_min_int
        };

        // Parse
        match input_min_int_str.parse::<f32>() {
            Ok(min_intensity) if min_intensity > 0.0 => {
                // Erase the prompt+input line
                print!("\r\x1B[K");
                // Print the selected value
                println!("||    Minimum intensity threshold selected: {min_intensity:.2}");
                return min_intensity;
            }
            _ => {
                // Erase the prompt+input line
                print!("\r\x1B[K");
                // Print error
                println!("||    !!! Invalid intensity. Please enter a positive number. !!!");
                // loop back
            }
        }
    }
}

fn prompt_mass_tolerance() -> f32 {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        // Prompt for mass tolerance
        print!("||    Please enter the desired mass tolerance (default: 0.01): ");
        let mut mass_tol_input = String::new();
        io::stdout().flush().unwrap();

        // Read into our buffer
        mass_tol_input.clear();
        io::stdin()
            .read_line(&mut mass_tol_input)
            .expect("Failed to read mass tolerance.");

        // Trim & apply default
        let trimmed_mass_tol = mass_tol_input.trim();
        let input_mass_tol_str = if trimmed_mass_tol.is_empty() {
            "0.01"
        } else {
            trimmed_mass_tol
        };

        // Parse
        match input_mass_tol_str.parse::<f32>() {
            Ok(mass_tolerance) if mass_tolerance > 0.0 => {
                // Erase the prompt+input line
                print!("\r\x1B[K");
                // Print the selected value
                println!("||    Mass tolerance selected: {mass_tolerance:.2}");
                return mass_tolerance;
            }
            _ => {
                // Erase the prompt+input line
                print!("\r\x1B[K");
                // Print error
                println!("||    !!! Invalid mass tolerance. Please enter a positive number. !!!");
                // loop back
            }
        }
    }
}

fn prompt_noise_threshold() -> u32 {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        // Prompt for the intensity of the noise threshold
        print!("||    Please enter the desired intensity of the noise threshold (default: 100): ");
        let mut noise_thresh_input = String::new();
        io::stdout().flush().unwrap();

        // Read into our buffer
        noise_thresh_input.clear();
        io::stdin()
            .read_line(&mut noise_thresh_input)
            .expect("Failed to read noise threshold.");

        // Trim & apply default
        let trimmed_noise_thresh = noise_thresh_input.trim();
        let input_noise_thresh_str = if trimmed_noise_thresh.is_empty() {
            "100"
        } else {
            trimmed_noise_thresh
        };

        // Parse
        match input_noise_thresh_str.parse::<u32>() {
            Ok(noise_threshold) if noise_threshold > 0 => {
                // Erase the prompt+input line
                print!("\r\x1B[K");
                // Print the selected value
                println!("||    Noise threshold selected: {}", noise_threshold);
                return noise_threshold;
            }
            _ => {
                // Erase the prompt+input line
                print!("\r\x1B[K");
                // Print error
                println!("||    !!! Invalid number. Please enter a positive integer. !!!");
                // loop back
            }
        }
    }
}

// fn prompt_verbose() -> bool {
//     loop {
//         println!("----------------------------------------------------------------------------------------------");
//         // Prompt for verbose flag
//         print!("||    Please enter the desired verbose flag (default: true): ");
//         let mut verbose_input = String::new();
//         io::stdout().flush().unwrap();

//         // Read into our buffer
//         verbose_input.clear();
//         io::stdin()
//             .read_line(&mut verbose_input)
//             .expect("Failed to read verbose flag.");

//         // Trim & apply default
//         let trimmed_verbose = verbose_input.trim();
//         let input_verbose_str = if trimmed_verbose.is_empty() {
//             "true"
//         } else {
//             trimmed_verbose
//         };

//         // Parse
//         match input_verbose_str.parse::<bool>() {
//             Ok(verbose) => {
//                 // Erase the prompt+input line
//                 print!("\r\x1B[K");
//                 // Print the selected value
//                 println!("||    Verbose flag selected: {}", verbose);
//                 return verbose;
//             }
//             _ => {
//                 // Erase the prompt+input line
//                 print!("\r\x1B[K");
//                 // Print error
//                 println!("||    !!! Invalid flag. Please enter true or false. !!!");
//                 // loop back
//             }
//         }
//     }
// }