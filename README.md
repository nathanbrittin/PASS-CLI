# MS Similarity Tool

A cross-platform command-line tool for untargeted mass spectrometry data similarity analysis. This tool performs similarity comparisons between all collected spectra in a single MS run and outputs the results as a CSV similarity matrix.

## Features

- **File Format Support**: mzML and mzXML files
- **Similarity Methods**: 
  - Cosine similarity
  - Modified cosine similarity (accounts for mass shifts and neutral losses)
- **Performance**: Multi-threaded processing using Rayon for fast computation
- **Cross-platform**: Works on Windows, macOS, and Linux
- **Flexible Filtering**: Filter spectra by intensity threshold and limit processing count

## Installation

### From Source

1. Install Rust (if not already installed):
   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   ```

2. Clone and build:
   ```bash
   git clone <repository-url>
   cd ms-similarity-tool
   cargo build --release
   ```

3. The executable will be available at `target/release/ms-similarity-tool`

## Usage

### Basic Usage

```bash
ms-similarity-tool -i input.mzML -o similarity_matrix.csv
```

### Advanced Options

```bash
ms-similarity-tool \
  --input data.mzML \
  --output results.csv \
  --method modified-cosine \
  --min-intensity 1000.0 \
  --mass-tolerance 0.02 \
  --max-spectra 1000 \
  --verbose
```

### Command Line Options

- `-i, --input <FILE>`: Input mzML or mzXML file (required)
- `-o, --output <FILE>`: Output CSV file for similarity matrix (required)
- `-m, --method <METHOD>`: Similarity method [`cosine`, `modified-cosine`] (default: `cosine`)
- `--min-intensity <FLOAT>`: Minimum intensity threshold for peaks (default: 1000.0)
- `--mass-tolerance <FLOAT>`: Mass tolerance for modified cosine similarity in Da (default: 0.02)
- `--max-spectra <INT>`: Maximum number of spectra to process (for testing)
- `-v, --verbose`: Enable verbose output
- `-h, --help`: Print help information

## Output Format

The tool outputs a CSV file containing the similarity matrix where:
- First row contains retention times as column headers
- First column contains retention times for each spectrum
- Each cell contains the similarity score between the corresponding spectra
- Values range from 0.0 (no similarity) to 1.0 (identical spectra)

Example output structure:
```csv
RT,1.234,2.456,3.789
1.234,1.000000,0.750000,0.200000
2.456,0.750000,1.000000,0.450000
3.789,0.200000,0.450000,1.000000
```

## Similarity Methods

### Cosine Similarity
- Standard cosine similarity between spectral vectors
- Matches peaks within 0.01 Da tolerance
- Good for general spectral comparison

### Modified Cosine Similarity
- Enhanced version that accounts for precursor mass differences
- Useful for comparing spectra with mass shifts or neutral losses
- Particularly effective for MS2 fragmentation spectra
- Uses user-defined mass tolerance (default: 0.02 Da)

## Performance Considerations

- The tool uses parallel processing for similarity calculations
- Memory usage scales with O(n²) where n is the number of spectra
- For large datasets (>10,000 spectra), consider using `--max-spectra` for testing
- Progress bar shows estimation of completion time when `--verbose` is enabled

## Future Development

This CLI tool serves as the foundation for a full interactive visualization tool that will include:
- Interactive heatmap visualization
- Chromatogram display with similarity mapping
- Real-time spectral comparison
- Advanced preprocessing options
- Extracted ion chromatogram support

## Dependencies

- `mzdata`: Mass spectrometry data parsing
- `ndarray`: N-dimensional array operations
- `rayon`: Data parallelism
- `clap`: Command line argument parsing
- `csv`: CSV file handling
- `indicatif`: Progress bars

## License

[Add your license information here]

## Contributing

[Add contribution guidelines here]