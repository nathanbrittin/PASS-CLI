use std::collections::HashMap;
// use std::error::Error;

use crate::ms_io::{import_mzml, Peak};

/// Compute a binary bit vector representation of a spectrum's peaks.
pub fn spectrum_to_dense_vec(peaks: &[Peak], bin_width: f32, max_mz: f32) -> Vec<f32> {
    let nbins = (max_mz / bin_width).ceil() as usize;
    let mut bins = vec![0.0; nbins];
    for p in peaks {
        if (0.0..=max_mz).contains(&p.mz) {
            let idx = (p.mz / bin_width) as usize;
            if idx < nbins {
                bins[idx] = p.intensity as f32;
            }
        }
    }
    bins
}

/// Compute Cosine similarity between two bit vectors.
/// Cosine similarity = dot(v1, v2) / (||v1|| * ||v2||)
pub fn cosine_similarity(v1: &[f32], v2: &[f32]) -> f32 {
    let n = v1.len().min(v2.len());
    let mut dot = 0.0;
    let mut norm1 = 0.0;
    let mut norm2 = 0.0;
    for i in 0..n {
        dot += v1[i] * v2[i];
        norm1 += v1[i] * v1[i];
        norm2 += v2[i] * v2[i];
    }
    if norm1 == 0.0 || norm2 == 0.0 {
        0.0
    } else {
        dot / (norm1.sqrt() * norm2.sqrt())
    }
}

/// Given a map of scan_id -> peaks, compute pairwise similarity matrix using cosine.
/// Returns a map of (scan_i, scan_j) -> similarity.
pub fn compute_pairwise_similarity(
    spec_map: &HashMap<String, Vec<Peak>>,
    bin_width: f32,
    max_mz: f32,
) -> HashMap<(String, String), f32> {
    let mut sims = HashMap::new();
    // Precompute bit vectors
    let mut bits_map: HashMap<String, Vec<f32>> = HashMap::new();
    for (scan, peaks) in spec_map {
        let bv = spectrum_to_dense_vec(peaks, bin_width, max_mz);
        bits_map.insert(scan.clone(), bv);
    }
    let scans: Vec<_> = bits_map.keys().cloned().collect();
    for i in 0..scans.len() {
        for j in (i + 1)..scans.len() {
            let a = &scans[i];
            let b = &scans[j];
            let v1 = &bits_map[a];
            let v2 = &bits_map[b];
            let sim = cosine_similarity(v1, v2);
            sims.insert((a.clone(), b.clone()), sim);
            sims.insert((b.clone(), a.clone()), sim);
        }
        // self-similarity = 1.0
        sims.insert((scans[i].clone(), scans[i].clone()), 1.0);
    }
    sims
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;
    use std::error::Error;

    #[test]
    fn test_bitvec_and_cosine() -> Result<(), Vec<String>> {
        // Load sample data
        let map = import_mzml("tests/data/1min.mzML")?;
        let peaks = map.get("39").expect("Scan 39 missing");
        // Determine max_mz from peaks
        let max_mz = peaks.iter().map(|p| p.mz).fold(0./0., f32::max);
        let bin_width = 0.01;
        let bv = spectrum_to_dense_vec(peaks, bin_width, max_mz);
        // Ensure the bit vector length is correct
        assert_eq!(bv.len(), (max_mz / bin_width).ceil() as usize);
        // Self similarity must be 1
        let sim = cosine_similarity(&bv, &bv);
        assert!((sim - 1.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn test_pairwise_all_files_cosine() -> Result<(), Box<dyn Error>> {
        let data_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests").join("data");
        let mut full_map: HashMap<String, Vec<Peak>> = HashMap::new();

        for entry_res in std::fs::read_dir(&data_dir)? {
            let entry = entry_res?;
            let p = entry.path();
            if p.extension().and_then(|e| e.to_str()) == Some("mzML") {
                let s = import_mzml(p.to_str().unwrap())
                    .map_err(|errs| errs.join(", "))?;
                for (scan, peaks) in s {
                    let key = format!("{}:{}", p.file_name().unwrap().to_string_lossy(), scan);
                    full_map.insert(key, peaks);
                }
            }
        }

        let sims = compute_pairwise_similarity(&full_map, 1.0, 2000.0);
        
        let self_count = sims.keys().filter(|(a, b)| a == b).count();
        assert_eq!(full_map.len(), self_count,
            "Expected {} self-similarity entries, found {}", full_map.len(), self_count);
        Ok(())
    }
}
