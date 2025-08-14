# PASS-CLI

<p align="center">
  <a href="https://github.com/<your‑username>/<your‑repo>">
    <img src="assets/pass-logo.png" alt="PASS logo" width="400" />
  </a>
</p>

<!-- [![Build Status](https://img.shields.io/github/actions/workflow/status/YourUser/{PASS}-CLI/ci.yml)](https://github.com/YourUser/PASS-CLI/actions) [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE) [![Tests](https://img.shields.io/github/actions/workflow/status/YourUser/PASS-CLI/tests.yml?label=tests)](https://github.com/YourUser/PASS-CLI/actions) -->

**PASS-CLI** (**P**airwise **A**nalyzer for **S**pectral **S**imilarity) is a high-performance, cross-platform command-line tool designed for comprehensive spectral similarity analysis in untargeted mass spectrometry workflows. Built with Rust for maximum speed and reliability, PASS-CLI efficiently processes MS/MS spectra to compute pairwise similarity matrices, enabling researchers to identify spectral relationships, detect molecular families, and explore chemical space within their datasets.

## What PASS-CLI Does

- **Analyzes spectral relationships**: Computes similarity scores between all MS/MS spectra pairs in your dataset using proven cosine similarity algorithms
- **Handles large datasets**: Leverages parallel processing to efficiently analyze thousands of spectra with optimal memory usage  
- **Provides flexible output**: Exports similarity matrices in multiple formats (CSV, TSV, JSON) for downstream analysis in your preferred tools
- **Visualizes results**: Generates publication-ready heatmaps (PNG, SVG, JPEG) to immediately visualize spectral clustering patterns
- **Streamlines workflows**: Interactive prompts guide users through analysis parameters without complex command-line syntax

Perfect for metabolomics researchers, natural product chemists, and mass spectrometry practitioners who need fast, reliable spectral comparison tools for data exploration and hypothesis generation.

---

## Table of Contents

1. [Features](#features)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Releases & Precompiled Binaries](releases-&-precompiled-binaries)
5. [Usage](#usage)
6. [Quick Start Example](#quick-start-example)
7. [Example Data & Expected Output](#example-data--expected-output)
8. [Output](#output)
9. [Similarity Methods](#similarity-methods)
10. [Performance](#performance)
11. [Roadmap](#roadmap)
12. [Tests](#tests)
13. [Contact & Support](#contact--support)
14. [Contributing](#contributing)
15. [License](#license)

---

## Features

### **Multi-format Input Support**
* **Native mzML & mzXML parsing**: Direct support for standard mass spectrometry formats without external dependencies
* **Robust file handling**: Automatic format detection and validation with detailed error reporting
* **Large file optimization**: Efficient memory management for processing datasets with thousands of spectra

### **Advanced Similarity Algorithms**
* **Cosine similarity**: Industry-standard vector cosine similarity with configurable m/z tolerance (default: 0.01 Da)
* **Modified cosine similarity**: Advanced algorithm accounting for neutral losses and precursor mass shifts for enhanced metabolite family detection

### **High-Performance Computing**
* **Parallel processing**: Multi-threaded computation using [Rayon](https://crates.io/crates/rayon) for optimal CPU utilization
* **Memory efficient**: Smart algorithms that scale gracefully with dataset size (O(N²) complexity optimization)
* **Progress tracking**: Real-time progress indicators and performance metrics during analysis

### **Flexible Data Processing**
* **Intelligent filtering**: Remove low-intensity peaks and noise with customizable thresholds
* **Spectrum selection**: Limit analysis to specific spectra subsets for rapid prototyping and testing
* **Quality control**: Built-in validation to ensure data integrity throughout the analysis pipeline

### **Comprehensive Output Options**
* **Multiple matrix formats**: Export similarity matrices as `CSV`, `TSV`, or `JSON` for seamless integration with R, Python, Excel, and other analysis tools
* **Publication-ready visualizations**: Generate high-quality heatmaps in multiple formats (`PNG`, `SVG`, `JPEG`)
* **Customizable themes**: Choose from 13 built-in color themes or create custom visualizations

### **User-Friendly Interface**
* **Interactive prompts**: Guided workflow eliminates the need to memorize complex command-line arguments
* **Smart defaults**: Sensible parameter defaults based on mass spectrometry best practices
* **Cross-platform compatibility**: Runs natively on Windows, macOS, and Linux with identical functionality

## Prerequisites for Developers

* Rust (1.60 or later) and `cargo`
* Operating systems: Windows, macOS, Linux

## Installation

Clone the repository and compile:

```bash
# Clone the project
git clone https://github.com/nathanbrittin/PASS-CLI.git
cd PASS-CLI
```

Compile either directly using cargo or using the build-all.sh bash script for automatically building binaries for all OS.

1. Build in cargo release mode
```
# Build in cargo release mode
cargo build --release
```

The above command creates an optimized executable in the `target/release/` directory.

2. Optionally execute build-all.sh in Unix or WSL command line.

```
.\build-all.sh
```

The resulting executables are:

* `target/release/pass-cli-<version>-linux-x64.tar.gz` (Unix/Linux/macOS)
* `target\release\pass-cli-<version>-windows-x64.zip` (Windows)

## Releases & Precompiled Binaries

Precompiled executables are available on the [GitHub Releases page](https://github.com/nathanbrittin/PASS-CLI/releases) for major platforms:

- **Linux/macOS**: `pass-cli-<version>-linux-x64.tar.gz`
- **Windows**: `pass-cli-<version>-windows-x64.zip`

### Installation from a Release

1. Download the appropriate `.tar.gz` (Linux/macOS) or `.zip` (Windows) file from the [latest release](https://github.com/nathanbrittin/PASS-CLI/releases).
2. Extract the archive and run the executable:

**Linux/macOS**
```sh
tar -xzf pass-cli-vX.Y.Z-linux.tar.gz
./pass-cli
```

**Windows**
1. Double-click "pass-cli.exe"
or
2. Open Command Prompt or PowerShell, navigate to the directory, and run the executable. Example:
```sh
.\pass-cli.exe
```

## Usage

PASS-CLI features an interactive, prompt-driven workflow to simplify operation for users who may be unfamiliar with complex command-line flags or syntax. On launch, the tool will guide you through each configuration step—there’s no need to remember options or refer back to documentation mid-run.

To get started, just execute the program:

```bash
# Unix
./target/release/pass-cli

# Windows
.\\target\\release\\pass-cli.exe
#or move the file to your directory
pass-cli.exe
```

You will be prompted to enter:

1. **Input file path** (mzML or mzXML)
2. **Output file path** (CSV, TSV, or JSON)
3. **Similarity metric** (`cosine` or `modified-cosine`)
4. **Minimum peak intensity** (e.g. `1000.0`)
5. **Noise Threshold** (eg. '1000.0')
6. **Mass tolerance** in Da (for modified cosine, e.g. `0.02`)

After entry, PASS‑CLI computes and saves the matrix.

## Quick Start Example

### Example Data & Expected Output

A few example mzML are provided under `tests\data`. These files have all been verified to work.

Assuming you use the small test file at `pass-cli\tests\data\FeatureFinderMetaboIdent_1_input.mzML`:

```bash
# Run the tool
./target/release/pass-cli
# or
pass-cli.exe

# Enter when prompted:
# Input file path: tests\data\FeatureFinderMetaboIdent_1_input.mzML
# Output file path: tests\data\FeatureFinderMetaboIdent_1_output.csv
# MS1 Similarity metric: cosine
# MS2 Similarity metric(s): cosine, modified-cosine
# MS1 Minimum peak intensity: 1.0
# MS2 Minimum peak intensity: 0.0
# Noise Threshold: 100.0
# Mass tolerance: 0.01
# Generate similarity matrix heatmap? [Y/n]: Y
# Path to save the heatmap image: tests\data\heatmap.png
# Output image format: png
# Choose a color theme: classic
# Continue with processing? [Y/n]: Y
```

Inspect `examples/test_similarity.csv`—it should look like:

```csv
  0,      2,        3
0,1.0000,0.0967,0.2710
2,0.0967,1.0000,0.1434
3,0.2710,0.1434,1.0000
```

### Output

The output file contains an N×N similarity matrix, where N is the number of spectra:

* The first row/column holds the scan identifiers or the scan index.
* Each cell `[i,j]` is the similarity score (0.0 – 1.0).

## Similarity Methods

PASS-CLI implements spectral similarity algorithms optimized for mass spectrometry data analysis. Both methods support configurable parameters to accommodate different experimental setups and research objectives.

### Cosine Similarity

**Standard vector cosine similarity** - the gold standard for spectral comparison in metabolomics.

**How it works:**
- Converts MS/MS spectra into high-dimensional vectors based on m/z and intensity values
- Computes the cosine of the angle between spectrum vectors in multi-dimensional space
- Produces similarity scores from 0.0 (completely different) to 1.0 (identical)

**Key features:**
- **Configurable m/z tolerance**: Default 0.01 Da, adjustable for different mass accuracy requirements
- **Intensity normalization**: Automatic peak intensity normalization for fair comparison
- **Peak matching**: Intelligent algorithm to align peaks across spectra within tolerance windows

**Best used for:**
- General spectral comparison and clustering
- Quality control and replicate analysis  
- Initial data exploration and pattern discovery
- Comparing spectra from the same ionization mode and similar fragmentation conditions

### Modified Cosine Similarity

**Advanced algorithm** designed specifically for metabolomics applications where neutral losses and mass shifts are biologically relevant.

**How it works:**
- Extends standard cosine similarity by considering shifted peak patterns
- Accounts for neutral losses (H₂O, NH₃, CO₂, etc.) commonly observed in metabolite fragmentation
- Compares not only direct peak matches but also shifted versions representing potential neutral losses
- Automatically detects and weights biologically relevant mass differences

**Key features:**
- **Precursor mass compensation**: Accounts for differences in precursor masses between related compounds
- **Configurable mass tolerance**: Fine-tune sensitivity for different analytical platforms

**Best used for:**
- **Metabolite family detection**: Identifying related compounds with similar core structures
- **Analog discovery**: Finding structural analogs and derivatives in complex mixtures
- **Natural product analysis**: Detecting compound series with systematic mass differences
- **Untargeted metabolomics**: Exploring chemical space for novel compound relationships

**Technical parameters:**
- **Mass tolerance**: Configurable tolerance for peak matching and neutral loss detection
- **Shift range**: Automatic detection of relevant mass shifts based on common neutral losses
- **Minimum intensity threshold**: Filters low-intensity peaks that may introduce noise

### Algorithm Performance Comparison

| Feature | Cosine Similarity | Modified Cosine Similarity |
|---------|-------------------|----------------------------|
| **Speed** | Very Fast | Fast |
| **Memory Usage** | Low | Moderate |
| **Best for Replicates** | Excellent | Excellent |
| **Best for Analogs** | Good | Excellent |
| **Parameter Sensitivity** | Low | Moderate |
| **Biological Relevance** | High | Very High |

### Choosing the Right Method

- **Start with Cosine Similarity** for initial data exploration, quality control, and when comparing highly similar spectra
- **Use Modified Cosine Similarity** for comprehensive metabolite annotation, analog discovery, and when exploring chemical diversity in complex samples
- **Consider your research goals**: Use both methods and compare results to gain complementary insights into your data


## Visualization & Heatmap Options

PASS-CLI automatically generates publication-friendly heatmaps to visualize spectral similarity patterns in your data. The interactive workflow guides you through visualization options, making it easy to create compelling figures for presentations and publications.

### Heatmap Features

- **Multiple output formats**: PNG, SVG, JPEG for different use cases
- **High resolution**: Optimized for both screen display and print publication
- **Automatic scaling**: Smart color mapping based on your similarity score distribution
- **Professional styling**: Clean, readable layouts with proper axis labels and legends

### Available Color Themes

PASS-CLI includes 13 carefully designed color themes optimized for scientific visualization. Some themes are designed to be colorblind-friendly and print-ready.

<table>
<tr>
<td align="center">
<br><strong>MS1 Cosine Score Heatmap</strong>
</td>
<td align="center">
<br><strong>MS2 Cosine Score Heatmap</strong>
</td>
<td align="center">
<br><strong>MS2 Modified Cosine Score Heatmap</strong>
</td>
<tr>
<td align="center">
<img src="examples\classic_output_ms1_cosine_heatmap.png" alt="Classic Theme" width="200"/>
<br><strong>-------------------</strong>
</td>
<td align="center">
<img src="examples\classic_output_ms2_cosine_heatmap.png" alt="Classic Theme" width="200"/>
<br><strong>Classic</strong>
</td>
<td align="center">
<img src="examples\classic_output_ms2_modified-cosine_heatmap.png" alt="Classic Theme" width="200"/>
<br><strong>-------------------</strong>
</td>
</tr>
<tr>
<td align="center">
<img src="examples\darkblue_output_ms1_cosine_heatmap.png" alt="Darkblue Theme" width="200"/>
<br><strong>-------------------</strong>
</td>
<td align="center">
<img src="examples\darkblue_output_ms2_cosine_heatmap.png" alt="Darkblue Theme" width="200"/>
<br><strong>Darkblue</strong>
</td>
<td align="center">
<img src="examples\darkblue_output_ms2_modified-cosine_heatmap.png" alt="Darkblue Theme" width="200"/>
<br><strong>-------------------</strong>
</td>
</tr>
<tr>
<td align="center">
<img src="examples\jet_output_ms1_cosine_heatmap.png" alt="Jet Theme" width="200"/>
<br><strong>-------------------</strong>
</td>
<td align="center">
<img src="examples\jet_output_ms2_cosine_heatmap.png" alt="Jet Theme" width="200"/>
<br><strong>Jet</strong>
</td>
<td align="center">
<img src="examples\jet_output_ms2_modified-cosine_heatmap.png" alt="Jet Theme" width="200"/>
<br><strong>-------------------</strong>
</td>
</tr>
<tr>
<td align="center">
<img src="examples\viridis_output_ms1_cosine_heatmap.png" alt="Viridis Theme" width="200"/>
<br><strong>-------------------</strong>
</td>
<td align="center">
<img src="examples\viridis_output_ms2_cosine_heatmap.png" alt="Viridis Theme" width="200"/>
<br><strong>Viridis</strong>
</td>
<td align="center">
<img src="examples\viridis_output_ms2_modified-cosine_heatmap.png" alt="Viridis Theme" width="200"/>
<br><strong>-------------------</strong>
</td>
</tr>
</table>

### Theme Recommendations

- **classic**: Traditional blue gradient, excellent for presentations
- **viridis**: Perceptually uniform, colorblind-safe, ideal for scientific publications
- **plasma/inferno/magma**: High-contrast themes perfect for highlighting strong similarities
- **blues/reds/greens/purples**: Monochromatic themes for specific aesthetic requirements
- **grays**: Grayscale option for print publications or accessibility needs

## Roadmap

* [x] Fix error handling of files not containing MS2 data
* [x] Implement the modified cosine score for effective comaprison of MS2 spectra
* [ ] Preprocessing filters (e.g., baseline subtraction, feature extraction)
* [ ] Chromatogram overlay with similarity mapping
* [ ] Export to network formats (e.g., GEXF, GraphML)
* [x] Output visualization scripts for heatmaps and network graphs
* [ ] Add colorbar option to heatmap generation
* [ ] Utilize ML/AI based comparative tools such as MS2DeepScore/Spec2Vec/DeepMASS/DreaMS

## Tests

Run unit and integration tests:

```bash
cargo test -- --nocapture
```

Integration tests live in `tests/`, using sample data in `tests/data`. CI runs these on each push.

## Contact & Support

* **GitHub Issues:** [https://github.com/nathanbrittin/PASS-CLI/issues](https://github.com/nathanbrittin/PASS-CLI/issues)
* **Email:** Nathan Brittin [nathanbrittin@gmail.com](mailto:nathanbrittin@gmail.com)

## Contributing

1. Fork the repo
2. Create a feature branch (`git checkout -b feature/YourFeature`)
3. Commit your changes (`git commit -m "Add feature X"`)
4. Push to your branch (`git push origin feature/YourFeature`)
5. Open a Pull Request

Please follow Rust formatting conventions (`cargo fmt`) and include tests where applicable.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE.md) for details.
