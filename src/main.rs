use std::process;
use std::io::{self, Write}; // Import io and Write traits for flushing output

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

    let mut input_path = String::new();
    // loop {
    //     println!("----------------------------------------------------------------------------------------------");
    //     println!("||    Please enter the full path to your input file:");
    //     std::io::stdin().read_line(&mut input_path).expect("Failed to read path, please verify that it has proper formatting.");

    //     let input_path = input_path.trim();
    //     let file_type = detect_input_file_type(input_path);

    //     if file_type.is_some() {
    //         println!("||    File selected: {}", input_path);
    //         println!("||    File type detected: {}", file_type.unwrap());
    //         break;
    //     } else {
    //         println!("
    //         *************File type not detected. Please verify that the file path is correct.*************
    //         ");
    //     }
    
        loop {
        // Print the prompt without a newline
        println!("----------------------------------------------------------------------------------------------");
        print!("||    Please enter the full path to your input file: ");
        io::stdout().flush().unwrap();

        // Read the line (this echoes the user’s typing)
        input_path.clear();                   // clear the buffer each time!
        io::stdin()
            .read_line(&mut input_path)
            .expect("Failed to read path.");

        let input_path = input_path.trim();
        let file_type = detect_input_file_type(input_path);

        if let Some(ft) = file_type {
            // Move cursor back to start of line, clear it, then print our formatted text
            // \r   → carriage return (go to col 0)
            // \x1B[K → ANSI “erase to end of line”
            print!("\r\x1B[K");             
            println!("||    File selected: {}", input_path);
            println!("||    File type detected: {}", ft);
            break;
        } else {
            // similarly overwrite the prompt line with an error message
            print!("\r\x1B[K");
            println!("||    !!! File type not detected. Please verify the file path. !!!");
            // then loop back and re‑prompt
        }
    }

    let mut output_path = String::new();
    loop {
        // Print multi-line prompt...
        println!("----------------------------------------------------------------------------------------------");
        println!("||    Please enter the desired path to your output file.");
        println!("||    Supported output formats are: csv, tsv, json.");
        println!("||    You can also leave it blank, which defaults to '/current_dir/output.csv'.");

        print!("||    Path: ");
        io::stdout().flush().unwrap();

        // Clear the buffer and read
        output_path.clear();
        io::stdin()
            .read_line(&mut output_path)
            .expect("Failed to read path, please verify formatting.");

        // Trim & handle default
        let mut trimmed = output_path.trim();
        if trimmed.is_empty() {
            trimmed = "output.csv";
        }
        let output_file_type = detect_output_file_type(trimmed);

        if let Some(ft) = output_file_type {
            print!("\r\x1B[K");
            println!("||    File selected: {}", trimmed);
            println!("||    File type detected: {}", ft);
            break;
        } else {
            // Erase the input line, then print your error
            print!("\r\x1B[K");
            println!("||    !!! File type not detected. Please verify the file path. !!!");
            // loop around and re-prompt
        }
    }

    //Ask the user for the similarity metric.
    let mut similarity_metric = String::new();
    loop {
        // Print the divider and context lines
        println!("----------------------------------------------------------------------------------------------");
        // Use `print!` for the input prompt, then flush so it appears immediately
        print!("||    Please enter the desired similarity metric [cosine(default) or modified-cosine]: ");
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
            break;
        } else {
            // Erase the prompt line, then print an error
            print!("\r\x1B[K");
            println!("||    !!! Similarity metric not detected. Please verify your input. !!!");
            // loop back to re-prompt
        }
    }

    //Ask the user for the minimum intensity threshold.
    let mut min_int_input = String::new();
    loop {
        // 1) Divider
        println!("----------------------------------------------------------------------------------------------");
        // 2) Prompt on one line
        print!("||    Please enter the desired minimum intensity threshold (default: 1000.0): ");
        io::stdout().flush().unwrap();

        // 3) Read into our buffer
        min_int_input.clear();
        io::stdin()
            .read_line(&mut min_int_input)
            .expect("Failed to read minimum intensity threshold.");

        // 4) Trim & apply default
        let trimmed_min_int = min_int_input.trim();
        let input_min_int_str = if trimmed_min_int.is_empty() {
            "1000.0"
        } else {
            trimmed_min_int
        };

        // 5) Parse
        match input_min_int_str.parse::<f64>() {
            Ok(min_intensity) if min_intensity > 0.0 => {
                // Erase the prompt+input line
                print!("\r\x1B[K");
                // Print the selected value
                println!("||    Minimum intensity threshold selected: {:.2}", min_intensity);
                break;
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