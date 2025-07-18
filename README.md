# SpectraMap-CLI

[![Build Status](https://img.shields.io/github/actions/workflow/status/YourUser/SpectraMap-CLI/ci.yml)](https://github.com/YourUser/SpectraMap-CLI/actions) [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE) [![Tests](https://img.shields.io/github/actions/workflow/status/YourUser/SpectraMap-CLI/tests.yml?label=tests)](https://github.com/YourUser/SpectraMap-CLI/actions)

SpectraMap-CLI is a fast, cross‑platform Rust command‑line tool for untargeted mass spectrometry data spectral similarity analysis. It computes pairwise similarity scores between all MS/MS spectra in a run and exports the resulting similarity matrix in CSV, TSV, or JSON formats.

---

## Table of Contents

1. [Features](#features)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Quick Start Example](#quick-start-example)
6. [Example Data & Expected Output](#example-data--expected-output)
7. [Output](#output)
8. [Similarity Methods](#similarity-methods)
9. [Performance](#performance)
10. [Roadmap](#roadmap)
11. [Tests](#tests)
12. [Contact & Support](#contact--support)
13. [Contributing](#contributing)
14. [License](#license)

---

## Features

* **Multi-format support**: Read `mzML` and `mzXML` files out of the box.
* **Similarity metrics**:

  * Standard cosine similarity
  * Modified cosine similarity (accounts for neutral losses and precursor mass shifts)
* **Parallel processing**: Built on [Rayon](https://crates.io/crates/rayon) for efficient multi‑threaded computation.
* **Flexible filtering**: Discard low‑intensity peaks and limit the number of spectra for rapid prototyping.
* **Output formats**: Export a full similarity matrix as `CSV`, `TSV`, or `JSON`.

## Prerequisites

* Rust (1.60 or later) and `cargo`
* Operating systems: Windows, macOS, Linux

## Installation

Clone the repository and compile:

```bash
# Clone the project
git clone https://github.com/YourUser/SpectraMap-CLI.git
cd SpectraMap-CLI

# Build in release mode
cargo build --release
```

The resulting executable is:

* `target/release/spectramap-cli` (Unix)
* `target\release\spectramap-cli.exe` (Windows)

## Usage

SpectraMap-CLI features an interactive, prompt-driven workflow to simplify operation for users who may be unfamiliar with complex command-line flags or syntax. On launch, the tool will guide you through each configuration step—there’s no need to remember options or refer back to documentation mid-run.

To get started, just execute the program:

```bash
# Unix
./target/release/spectramap-cli

# Windows
.\\target\\release\\spectramap-cli.exe
```

You will be prompted to enter:

1. **Input file path** (mzML or mzXML)
2. **Output file path** (CSV, TSV, or JSON)
3. **Similarity metric** (`cosine` or `modified-cosine`)
4. **Minimum peak intensity** (e.g. `1000.0`)
5. **Mass tolerance** in Da (for modified cosine, e.g. `0.02`)
6. **Maximum number of spectra** (or leave blank to process all)
7. **Output format** (`csv`, `tsv`, or `json`)
8. **Verbose mode** (`y`/`n`)

After entry, SpectraMap‑CLI computes and saves the matrix.

## Quick Start Example

Assuming you have a small test file at `examples/test_run.mzML`:

```bash
# Run the tool
./target/release/spectramap-cli

# Enter when prompted:
# Input file path: examples/test_run.mzML
# Output file path: examples/test_similarity.csv
# Similarity metric: cosine
# Minimum peak intensity: 1000.0
# Mass tolerance: 0.02
# Maximum number of spectra: [ENTER]
# Output format: csv
# Verbose mode: n
```

Inspect `examples/test_similarity.csv`—it should look like:

```csv
RT,45.1,82.3,120.5
45.1,1.000000,0.812345,0.234567
82.3,0.812345,1.000000,0.456789
120.5,0.234567,0.456789,1.000000
```

## Example Data & Expected Output

A tiny sample mzML is provided under `examples/`. The expected similarity matrices (`.csv`, `.tsv`, and `.json`) live alongside it, so you can verify correct installation without your own data.

## Output

The output file contains an N×N similarity matrix, where N is the number of spectra:

* The first row/column holds spectrum retention times or scan identifiers.
* Each cell `[i,j]` is the similarity score (0.0 – 1.0).

## Similarity Methods

### Cosine Similarity

* Standard vector cosine similarity with default peak m/z tolerance of `0.01` Da.
* Good for overall spectral shape comparisons.

### Modified Cosine Similarity

* Compensates for mass shifts and neutral losses between MS/MS spectra.
* Suitable for detecting related fragmentation patterns across different precursor masses.
* User‑configurable tolerance (`--mass-tolerance`).

## Performance

* Utilizes Rayon to parallelize O(N²) comparisons.
* Memory usage scales quadratically with the number of spectra.
* For large runs (> 10,000 spectra), consider downsampling or limiting via prompts.
* `Verbose` mode prints per‑step timings.

## Roadmap

* [ ] Preprocessing filters (e.g., baseline subtraction, decharging)
* [ ] Chromatogram overlay with similarity mapping
* [ ] Export to network formats (e.g., GEXF, GraphML)
* [ ] Output visualization scripts for heatmaps and network graphs

## Tests

Run unit and integration tests:

```bash
cargo test -- --nocapture
```

Integration tests live in `tests/`, using sample data in `examples/`. CI runs these on each push.

## Contact & Support

* **GitHub Issues:** [https://github.com/YourUser/SpectraMap-CLI/issues](https://github.com/YourUser/SpectraMap-CLI/issues)
* **Email:** Nathan Brittin [nathan.brittin@wisc.edu](mailto:nathan.brittin@wisc.edu)

## Contributing

1. Fork the repo
2. Create a feature branch (`git checkout -b feature/YourFeature`)
3. Commit your changes (`git commit -m "Add feature X"`)
4. Push to your branch (`git push origin feature/YourFeature`)
5. Open a Pull Request

Please follow Rust formatting conventions (`cargo fmt`) and include tests where applicable.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
