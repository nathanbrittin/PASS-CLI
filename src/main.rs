use std::process;
use std::path::Path;
use std::io::{self, Write}; // Import io and Write traits for flushing output
mod ms_io;
use ms_io::import_mzml;

fn main() {

    let art = " ____                  _             __  __             \n/ ___| _ __   ___  ___| |_ _ __ __ _|  \\/  | __ _ _ __  \n\\___ \\| '_ \\ / _ \\/ __| __| '__/ _` | |\\/| |/ _` | '_ \\ \n ___) | |_) |  __/ (__| |_| | | (_| | |  | | (_| | |_) |\n|____/| .__/ \\___|\\___|\\__|_|  \\__,_|_|  |_|\\__,_| .__/ \n      |_|                                        |_|    \n";
    println!("{}", art);

    //First, give a welcome and description of what this command line interface does.
    println!("**Welcome to SpectraMap**
----------------------------------------------------------------------------------------------
||  -A Rust-powered CLI for mass spectrometry analysis.                                     
||  -Parsing `.mzML` files for spectral data                                                
||  -Computing pairwise similarity with cosine and modified cosine metrics                  
||  -Building exportable similarity matrices for downstream visualization and analysis      
||  -Supports file formats: `.mzML` and `.mzXML`                                            
----------------------------------------------------------------------------------------------
||  Please follow the directions to process your data.                                      
||  Alternatively, provide a config.yaml file to automate processing.                       
----------------------------------------------------------------------------------------------"
    );

    let input_path = prompt_input_path();

    let output_path = prompt_output_path();

    let similarity_metric = prompt_similarity_metric();

    let minimum_intensity = prompt_min_intensity();

    let max_peaks = prompt_max_peaks();

    let mass_tolerance = prompt_mass_tolerance();

    let verbose = prompt_verbose();

    // Print a final confirmation
    println!("----------------------------------------------------------------------------------------------");
    println!("||    Input file: {}", input_path);
    println!("||    Output file: {}", output_path);
    println!("||    Similarity metric: {}", similarity_metric);
    println!("||    Minimum intensity threshold: {:.2}", minimum_intensity);
    println!("||    Mass tolerance: {:.2}", mass_tolerance);
    println!("||    Maximum number of peaks: {}", max_peaks);
    println!("||    Verbose flag: {}", verbose);

    // Wait for the user to press enter
    println!("----------------------------------------------------------------------------------------------");
    println!("||    Press enter to continue...");
    let mut input = String::new();
    io::stdin()
        .read_line(&mut input)
        .expect("Failed to read line");

    // Importing the data file
    let test_file: &str = &input_path;
    // let test_file: &str = "tests/data/1min.mzML"; // Placeholder for testing
    // let test_file: &str = "tiny.mzML";
    // // Verify the file exists
    // if !verify_file_exists(test_file) {
    //     println!("Error: The file '{}' does not exist.", test_file);
    //     process::exit(1);
    // } else {
    //     println!("The file '{}' exists.", test_file);
    // }

    match import_mzml(test_file) {
        Ok(map) => {
            println!("Parsed {} spectra successfully.", map.len());
        }
        Err(errs) => {
            eprintln!("Failed to import {} with {} errors:", test_file, errs.len());
            for e in errs {
                eprintln!("  - {}", e);
            }
            process::exit(1);
        }
    }

    println!("Success so far!");
    process::exit(0);

}

fn detect_input_file_type(file_path: &str) -> Option<&str> {
    match file_path.split('.').last() {
        Some("mzML") | Some("mzml") => Some("mzML"),
        Some("mzXML") | Some("mzxml") => Some("mzXML"),
        _ => None,
    }
}

fn detect_output_file_type(file_path: &str) -> Option<&str> {
    match file_path.split('.').last() {
        Some("") => Some("csv"),
        Some("csv") => Some("csv"),
        Some("tsv") => Some("tsv"),
        Some("json") => Some("json"),
        _ => None,
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

fn verify_file_exists(file_path: &str) -> bool {
    Path::new(file_path).exists()
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

        // Try to detect file type
        if let Some(ft) = detect_input_file_type(trimmed) {
            // Clear the current prompt line and confirm selection
            print!("\r\x1B[K");
            println!("||    File selected: {}", trimmed);
            println!("||    File type detected: {}", ft);
            return trimmed.to_string();
        } else {
            // Clear the prompt and show an error, then loop again
            print!("\r\x1B[K");
            println!("||    !!! File type not detected. Please verify the file path. !!!");
        }
    };
}

fn prompt_output_path() -> String {
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
            println!("||    File selected: {}", trimmed);
            println!("||    File type detected: {}", ft);
            return trimmed.to_string();
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
            println!("||    Similarity metric selected: {}", m);
            return m.to_string();
        } else {
            // Erase the prompt line, then print an error
            print!("\r\x1B[K");
            println!("||    !!! Similarity metric not detected. Please verify your input. !!!");
            // loop back to re-prompt
        }
    }
}

fn prompt_min_intensity() -> f64 {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        // Prompt for minimum intensity threshold
        print!("||    Please enter the desired minimum intensity threshold (default: 1000.0): ");
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
            "1000.0"
        } else {
            trimmed_min_int
        };

        // Parse
        match input_min_int_str.parse::<f64>() {
            Ok(min_intensity) if min_intensity > 0.0 => {
                // Erase the prompt+input line
                print!("\r\x1B[K");
                // Print the selected value
                println!("||    Minimum intensity threshold selected: {:.2}", min_intensity);
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

fn prompt_mass_tolerance() -> f64 {
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
        match input_mass_tol_str.parse::<f64>() {
            Ok(mass_tolerance) if mass_tolerance > 0.0 => {
                // Erase the prompt+input line
                print!("\r\x1B[K");
                // Print the selected value
                println!("||    Mass tolerance selected: {:.2}", mass_tolerance);
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

fn prompt_max_peaks() -> u32 {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        // Prompt for maximum number of peaks
        print!("||    Please enter the desired maximum number of peaks (default: 1000): ");
        let mut max_peaks_input = String::new();
        io::stdout().flush().unwrap();

        // Read into our buffer
        max_peaks_input.clear();
        io::stdin()
            .read_line(&mut max_peaks_input)
            .expect("Failed to read maximum number of peaks.");

        // Trim & apply default
        let trimmed_max_peaks = max_peaks_input.trim();
        let input_max_peaks_str = if trimmed_max_peaks.is_empty() {
            "1000"
        } else {
            trimmed_max_peaks
        };

        // Parse
        match input_max_peaks_str.parse::<u32>() {
            Ok(max_peaks) if max_peaks > 0 => {
                // Erase the prompt+input line
                print!("\r\x1B[K");
                // Print the selected value
                println!("||    Maximum number of peaks selected: {}", max_peaks);
                return max_peaks;
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

fn prompt_verbose() -> bool {
    loop {
        println!("----------------------------------------------------------------------------------------------");
        // Prompt for verbose flag
        print!("||    Please enter the desired verbose flag (default: true): ");
        let mut verbose_input = String::new();
        io::stdout().flush().unwrap();

        // Read into our buffer
        verbose_input.clear();
        io::stdin()
            .read_line(&mut verbose_input)
            .expect("Failed to read verbose flag.");

        // Trim & apply default
        let trimmed_verbose = verbose_input.trim();
        let input_verbose_str = if trimmed_verbose.is_empty() {
            "true"
        } else {
            trimmed_verbose
        };

        // Parse
        match input_verbose_str.parse::<bool>() {
            Ok(verbose) => {
                // Erase the prompt+input line
                print!("\r\x1B[K");
                // Print the selected value
                println!("||    Verbose flag selected: {}", verbose);
                return verbose;
            }
            _ => {
                // Erase the prompt+input line
                print!("\r\x1B[K");
                // Print error
                println!("||    !!! Invalid flag. Please enter true or false. !!!");
                // loop back
            }
        }
    }
}