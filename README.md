# PASS-CLI

<p align="center">
  <a href="https://github.com/nathanbrittin/PASS-CLI">
    <img src="assets/pass-logo.png" alt="PASS logo" width="400" />
  </a>
</p>

[![CI](https://img.shields.io/github/actions/workflow/status/nathanbrittin/PASS-CLI/rust.yml?label=CI)](https://github.com/nathanbrittin/PASS-CLI/actions) [![Release](https://img.shields.io/github/v/release/nathanbrittin/PASS-CLI)](https://github.com/nathanbrittin/PASS-CLI/releases) [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)

**PASS-CLI** (**P**airwise **A**nalyzer for **S**pectral **S**imilarity) is a high-performance, cross-platform command-line tool for spectral similarity analysis in untargeted mass spectrometry workflows. Built with Rust for maximum speed and reliability, PASS-CLI computes pairwise similarity matrices across all MS1 and MS2 spectra in a dataset, helping researchers identify spectral relationships, detect molecular families, and explore chemical space.

## What PASS-CLI Does

- **Analyzes spectral relationships**: Computes pairwise similarity scores across all spectra using multiple algorithms simultaneously
- **Handles large datasets**: Multi-threaded parallel computation scales to thousands of spectra
- **Adaptive noise filtering**: Three intensity thresholding modes remove noise without requiring instrument-specific knowledge
- **Flexible output**: Exports similarity matrices as CSV, TSV, or JSON for downstream analysis
- **Visualizes results**: Generates publication-ready heatmaps with chromatogram overlays in PNG, SVG, or JPEG

---

## Table of Contents

1. [Features](#features)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Releases & Precompiled Binaries](#releases--precompiled-binaries)
5. [Usage](#usage)
6. [Quick Start Example](#quick-start-example)
7. [Output](#output)
8. [Similarity Methods](#similarity-methods)
9. [Visualization & Heatmap Options](#visualization--heatmap-options)
10. [Performance](#performance)
11. [Roadmap](#roadmap)
12. [Tests](#tests)
13. [Contact & Support](#contact--support)
14. [Contributing](#contributing)
15. [License](#license)

---

## Features

### Multi-format Input Support
- Native mzML and mzXML parsing with no external dependencies
- Automatic format detection and validation with detailed error reporting
- Efficient memory management for large datasets with thousands of spectra

### Ten Similarity Algorithms
- Five standard metrics for both MS1 and MS2 spectra
- Five modified (precursor-mass shift) variants for MS2 analog detection
- Select any combination of metrics in a single run

### Adaptive Intensity Thresholding
- **Absolute**: Fixed floor value (backward-compatible, good when the noise floor is known)
- **Percent Base Peak (PBP)**: Keeps peaks above a fraction of each scan's own maximum — adapts automatically to signal level
- **Top-N**: Retains only the N most intense peaks per scan — robust to dynamic range differences

### High-Performance Computing
- Multi-threaded pairwise computation via [Rayon](https://crates.io/crates/rayon)
- Sparse vector representation for memory-efficient similarity calculation
- MAD-based SNR noise filtering before similarity computation

### Comprehensive Output Options
- Similarity matrices in CSV, TSV, or JSON
- Publication-ready heatmaps in PNG, SVG, or JPEG
- 13 built-in color themes including colorblind-safe options
- Chromatogram overlay (BPI or TIC) on heatmaps

### User-Friendly Interface
- Interactive prompts guide through all parameters — no flags to memorize
- Accept metric names or numbers interchangeably (`cosine`, `1`, `entropy`, etc.)
- Sensible defaults based on mass spectrometry best practices
- Runs identically on Windows, macOS, and Linux

---

## Prerequisites

To build from source you need Rust 1.60 or later and `cargo`. Download from [rustup.rs](https://rustup.rs).

Precompiled binaries require no dependencies — download and run.

---

## Installation

### From Source

```bash
git clone https://github.com/nathanbrittin/PASS-CLI.git
cd PASS-CLI
cargo build --release
```

The optimized executable is placed at `target/release/pass-cli` (Linux/macOS) or `target/release/pass-cli.exe` (Windows).

---

## Releases & Precompiled Binaries

Every tagged release on GitHub includes precompiled binaries for four platforms, built automatically by GitHub Actions:

| Platform | Archive |
|----------|---------|
| Windows x64 | `pass-cli-vX.Y.Z-windows-x64.zip` |
| Linux x64 | `pass-cli-vX.Y.Z-linux-x64.tar.gz` |
| macOS Intel (x86_64) | `pass-cli-vX.Y.Z-macos-x64.tar.gz` |
| macOS Apple Silicon (aarch64) | `pass-cli-vX.Y.Z-macos-arm64.tar.gz` |

Download from the [latest release](https://github.com/nathanbrittin/PASS-CLI/releases), extract, and run.

**Linux/macOS:**
```sh
tar -xzf pass-cli-vX.Y.Z-linux-x64.tar.gz
./pass-cli
```

**Windows** — double-click `pass-cli.exe`, or from a terminal:
```sh
.\pass-cli.exe
```

---

## Usage

PASS-CLI uses an interactive prompt-driven workflow. Launch the executable and follow the prompts:

```bash
# Linux/macOS
./pass-cli

# Windows
.\pass-cli.exe
```

The prompts walk you through:

1. **Input file path** — mzML or mzXML
2. **Output file path** — CSV, TSV, or JSON
3. **MS level** — MS1, MS2, or both
4. **Similarity metrics** — enter names or numbers (can select multiple for MS2)
5. **MS1 intensity threshold** — choose Absolute / PercentBasePeak / TopN and enter the value
6. **MS2 intensity threshold** — same three options
7. **Mass tolerance** in Da (default 0.01)
8. **Heatmap options** — format, color theme, output path

---

## Quick Start Example

Using the small test file at `tests/data/FeatureFinderMetaboIdent_1_input.mzML`:

```
Input file path: tests/data/FeatureFinderMetaboIdent_1_input.mzML
Output file path: output.csv
MS level to analyze: both
MS1 similarity metric: cosine
MS2 similarity metric(s): cosine, modified-cosine
MS1 intensity threshold method: 2  (PercentBasePeak)
MS1 PBP fraction (e.g. 0.01): 0.01
MS2 intensity threshold method: 2  (PercentBasePeak)
MS2 PBP fraction: 0.01
Mass tolerance (Da): 0.01
Generate heatmap? [Y/n]: Y
Heatmap output path: heatmap.png
Output format: png
Color theme: classic
Continue? [Y/n]: Y
```

---

## Output

Each selected metric produces its own similarity matrix file. For example, selecting `cosine` and `modified-cosine` for MS2 produces `output_ms2_cosine.csv` and `output_ms2_modified-cosine.csv`.

Each matrix is N×N where N is the number of spectra:

- Row and column headers are scan identifiers
- Each cell `[i,j]` holds the similarity score (0.0–1.0, except Bray-Curtis which is a dissimilarity)

A `config.toml` is also written alongside each run, recording all parameters used.

---

## Similarity Methods

PASS-CLI implements ten metrics. The five standard metrics run on both MS1 and MS2 spectra. The five modified variants are MS2-only — they shift one spectrum's m/z bins by the precursor mass difference before computing the score, allowing analogs with the same fragmentation pattern but different precursor masses to score highly.

### Standard Metrics (MS1 + MS2)

| # | Name | Key | Description |
|---|------|-----|-------------|
| 1 | Cosine | `cosine` | Classic dot product of L2-normalized intensity vectors. Fast and widely used. |
| 2 | Weighted Cosine | `weighted-cosine` | Weights each peak by m/z² × √intensity before comparison, emphasizing high-mass, high-intensity ions. |
| 3 | Spectral Entropy | `entropy` | Uses Shannon entropy of the intensity distribution. More tolerant to missing peaks and noise than cosine. Score of 1.0 means identical. |
| 4 | Hellinger | `hellinger` | Bhattacharyya coefficient on intensity probability distributions. Naturally bounded [0,1], robust to outlier intensities. |
| 5 | Bray-Curtis | `bray-curtis` | Returns **dissimilarity** (0 = identical, 1 = no overlap). Useful for ecological-style diversity analysis of spectral features. |

### Modified Variants (MS2 only)

Each modified metric applies a precursor mass shift to one spectrum before computing the corresponding standard score:

| # | Name | Key |
|---|------|-----|
| 6 | Modified Cosine | `modified-cosine` |
| 7 | Modified Weighted Cosine | `modified-weighted-cosine` |
| 8 | Modified Spectral Entropy | `modified-entropy` |
| 9 | Modified Hellinger | `modified-hellinger` |
| 10 | Modified Bray-Curtis | `modified-bray-curtis` |

### Choosing a Metric

- **General exploration / QC**: `cosine` — fast, interpretable, well-established
- **Noise tolerance / missing peaks**: `entropy` — less sensitive to low-intensity contaminants
- **Analog and metabolite family detection**: `modified-cosine` or `modified-entropy`
- **Emphasizing defining high-mass ions**: `weighted-cosine`
- **Diversity analysis**: `bray-curtis` (note: dissimilarity, not similarity)
- **Multiple perspectives**: select several metrics in one run and compare the resulting matrices

---

## Visualization & Heatmap Options

PASS-CLI generates publication-ready heatmaps with a chromatogram overlay for each output matrix. The heatmap workflow is guided interactively.

### Heatmap Features

- Output formats: PNG, SVG, JPEG
- Chromatogram overlay: Base Peak Intensity (BPI) or Total Ion Current (TIC) with smoothing
- Automatic color scaling based on the score distribution
- Clean axis labels with scan identifiers

### Available Color Themes

13 built-in themes: `classic`, `darkblue`, `jet`, `viridis`, `plasma`, `inferno`, `magma`, `blues`, `reds`, `greens`, `purples`, `grays`, `rainbow`.

<table>
<tr>
<td align="center"><br><strong>MS1 Cosine</strong></td>
<td align="center"><br><strong>MS2 Cosine</strong></td>
<td align="center"><br><strong>MS2 Modified Cosine</strong></td>
</tr>
<tr>
<td align="center">
<img src="examples/classic_output_ms1_cosine_heatmap.png" alt="Classic MS1" width="200"/>
<br><strong>classic</strong>
</td>
<td align="center">
<img src="examples/classic_output_ms2_cosine_heatmap.png" alt="Classic MS2 Cosine" width="200"/>
<br><strong>classic</strong>
</td>
<td align="center">
<img src="examples/classic_output_ms2_modified-cosine_heatmap.png" alt="Classic MS2 Modified" width="200"/>
<br><strong>classic</strong>
</td>
</tr>
<tr>
<td align="center">
<img src="examples/darkblue_output_ms1_cosine_heatmap.png" alt="Darkblue MS1" width="200"/>
<br><strong>darkblue</strong>
</td>
<td align="center">
<img src="examples/darkblue_output_ms2_cosine_heatmap.png" alt="Darkblue MS2 Cosine" width="200"/>
<br><strong>darkblue</strong>
</td>
<td align="center">
<img src="examples/darkblue_output_ms2_modified-cosine_heatmap.png" alt="Darkblue MS2 Modified" width="200"/>
<br><strong>darkblue</strong>
</td>
</tr>
<tr>
<td align="center">
<img src="examples/jet_output_ms1_cosine_heatmap.png" alt="Jet MS1" width="200"/>
<br><strong>jet</strong>
</td>
<td align="center">
<img src="examples/jet_output_ms2_cosine_heatmap.png" alt="Jet MS2 Cosine" width="200"/>
<br><strong>jet</strong>
</td>
<td align="center">
<img src="examples/jet_output_ms2_modified-cosine_heatmap.png" alt="Jet MS2 Modified" width="200"/>
<br><strong>jet</strong>
</td>
</tr>
<tr>
<td align="center">
<img src="examples/viridis_output_ms1_cosine_heatmap.png" alt="Viridis MS1" width="200"/>
<br><strong>viridis</strong>
</td>
<td align="center">
<img src="examples/viridis_output_ms2_cosine_heatmap.png" alt="Viridis MS2 Cosine" width="200"/>
<br><strong>viridis</strong>
</td>
<td align="center">
<img src="examples/viridis_output_ms2_modified-cosine_heatmap.png" alt="Viridis MS2 Modified" width="200"/>
<br><strong>viridis</strong>
</td>
</tr>
</table>

**Theme recommendations:**
- `classic`: Traditional blue gradient, good for presentations
- `viridis`: Perceptually uniform, colorblind-safe — preferred for publications
- `plasma` / `inferno` / `magma`: High contrast, good for highlighting strong similarities
- `grays`: Grayscale for print or accessibility

---

## Performance

PASS-CLI is optimized for speed:

- **Sparse vector computation** — only shared m/z bins are evaluated, skipping zero-intensity regions
- **Rayon parallelism** — all N×N pairs computed in parallel across CPU cores
- **MAD noise filtering** — reduces vector dimensionality before similarity computation
- **Release profile** — `opt-level = 3`, `lto = true`, `codegen-units = 1`

---

## Roadmap

- [x] Modified cosine similarity for MS2 analog detection
- [x] Chromatogram overlay on heatmaps
- [x] Spectral entropy, weighted cosine, Hellinger, and Bray-Curtis metrics
- [x] Modified variants of all new metrics
- [x] Adaptive intensity thresholding (Absolute, PercentBasePeak, TopN)
- [x] Metric and threshold selection by name or number
- [x] Cross-platform CI and automated GitHub Release builds
- [ ] Retention time filter (focus on spectra within a specific RT window)
- [ ] Baseline subtraction and additional preprocessing filters
- [ ] Export to network formats (GEXF, GraphML) for molecular networking
- [ ] Colorbar option for heatmap generation
- [ ] Config file support to skip interactive prompts for repeated analyses

---

## Tests

Run all tests with output:

```bash
cargo test -- --nocapture
```

Integration tests live in `tests/` and use real mzML files from `tests/data/`. The CI workflow runs tests on Ubuntu, macOS, and Windows on every push.

---

## Contact & Support

- **GitHub Issues**: [https://github.com/nathanbrittin/PASS-CLI/issues](https://github.com/nathanbrittin/PASS-CLI/issues)
- **Email**: Nathan Brittin — [nathanbrittin@gmail.com](mailto:nathanbrittin@gmail.com)

---

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Commit your changes (`git commit -m "Add feature X"`)
4. Push to your branch (`git push origin feature/your-feature`)
5. Open a Pull Request

Please run `cargo fmt` and `cargo test` before submitting.

---

## License

This project is licensed under the MIT License. See [LICENSE.md](LICENSE.md) for details.
