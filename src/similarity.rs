use crate::ms_data::Spectrum;
use clap::ValueEnum;

#[derive(Debug, Clone, ValueEnum)]
pub enum SimilarityMethod {
    #[value(name = "cosine")]
    Cosine,
    #[value(name = "modified-cosine")]
    ModifiedCosine,
}

/// Calculate cosine similarity between two spectra
pub fn cosine_similarity(spec1: &Spectrum, spec2: &Spectrum) -> f64 {
    if spec1.peaks.is_empty() || spec2.peaks.is_empty() {
        return 0.0;
    }

    let mut dot_product = 0.0;
    let mut norm1 = 0.0;
    let mut norm2 = 0.0;

    // Simple approach: find matching peaks within tolerance
    let mass_tolerance = 0.01; // 0.01 Da tolerance

    for peak1 in &spec1.peaks {
        for peak2 in &spec2.peaks {
            if (peak1.mz - peak2.mz).abs() <= mass_tolerance {
                dot_product += peak1.intensity * peak2.intensity;
                break; // Only match each peak once
            }
        }
    }

    // Calculate norms
    for peak in &spec1.peaks {
        norm1 += peak.intensity * peak.intensity;
    }
    for peak in &spec2.peaks {
        norm2 += peak.intensity * peak.intensity;
    }

    norm1 = norm1.sqrt();
    norm2 = norm2.sqrt();

    if norm1 == 0.0 || norm2 == 0.0 {
        0.0
    } else {
        dot_product / (norm1 * norm2)
    }
}

/// Calculate modified cosine similarity between two spectra
/// This version accounts for mass shifts and neutral losses
pub fn modified_cosine_similarity(spec1: &Spectrum, spec2: &Spectrum, mass_tolerance: f64) -> f64 {
    if spec1.peaks.is_empty() || spec2.peaks.is_empty() {
        return 0.0;
    }

    // For MS2 spectra, we can use precursor information
    let precursor_diff = match (spec1.precursor_mz, spec2.precursor_mz) {
        (Some(p1), Some(p2)) => p1 - p2,
        _ => 0.0,
    };

    let mut dot_product = 0.0;
    let mut norm1 = 0.0;
    let mut norm2 = 0.0;

    // Modified cosine: try to match peaks with and without precursor mass shift
    for peak1 in &spec1.peaks {
        let mut best_match_intensity = 0.0;
        
        for peak2 in &spec2.peaks {
            // Try direct matching
            if (peak1.mz - peak2.mz).abs() <= mass_tolerance {
                best_match_intensity = best_match_intensity.max(peak2.intensity);
            }
            
            // Try matching with precursor mass difference
            if precursor_diff != 0.0 {
                let shifted_mz = peak2.mz + precursor_diff;
                if (peak1.mz - shifted_mz).abs() <= mass_tolerance {
                    best_match_intensity = best_match_intensity.max(peak2.intensity);
                }
            }
        }
        
        if best_match_intensity > 0.0 {
            dot_product += peak1.intensity * best_match_intensity;
        }
    }

    // Calculate norms
    for peak in &spec1.peaks {
        norm1 += peak.intensity * peak.intensity;
    }
    for peak in &spec2.peaks {
        norm2 += peak.intensity * peak.intensity;
    }

    norm1 = norm1.sqrt();
    norm2 = norm2.sqrt();

    if norm1 == 0.0 || norm2 == 0.0 {
        0.0
    } else {
        dot_product / (norm1 * norm2)
    }
}

/// Calculate spectral entropy (for future use in entropy-based similarity)
pub fn spectral_entropy(spectrum: &Spectrum) -> f64 {
    if spectrum.peaks.is_empty() {
        return 0.0;
    }

    let total_intensity: f64 = spectrum.peaks.iter().map(|p| p.intensity).sum();
    if total_intensity == 0.0 {
        return 0.0;
    }

    let mut entropy = 0.0;
    for peak in &spectrum.peaks {
        let probability = peak.intensity / total_intensity;
        if probability > 0.0 {
            entropy -= probability * probability.log2();
        }
    }

    entropy
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ms_data::Peak;

    fn create_test_spectrum(peaks: Vec<(f64, f64)>) -> Spectrum {
        Spectrum {
            retention_time: 1.0,
            precursor_mz: None,
            peaks: peaks.into_iter().map(|(mz, intensity)| Peak { mz, intensity }).collect(),
            ms_level: 1,
        }
    }

    #[test]
    fn test_cosine_similarity_identical() {
        let spec = create_test_spectrum(vec![(100.0, 1000.0), (200.0, 2000.0), (300.0, 1500.0)]);
        let similarity = cosine_similarity(&spec, &spec);
        assert!((similarity - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cosine_similarity_orthogonal() {
        let spec1 = create_test_spectrum(vec![(100.0, 1000.0), (200.0, 2000.0)]);
        let spec2 = create_test_spectrum(vec![(300.0, 1000.0), (400.0, 2000.0)]);
        let similarity = cosine_similarity(&spec1, &spec2);
        assert!(similarity.abs() < 1e-10);
    }

    #[test]
    fn test_cosine_similarity_partial_overlap() {
        let spec1 = create_test_spectrum(vec![(100.0, 1000.0), (200.0, 2000.0), (300.0, 1500.0)]);
        let spec2 = create_test_spectrum(vec![(100.0, 1000.0), (400.0, 2000.0)]);
        let similarity = cosine_similarity(&spec1, &spec2);
        assert!(similarity > 0.0 && similarity < 1.0);
    }
}