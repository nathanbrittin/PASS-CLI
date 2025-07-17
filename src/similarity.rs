use std::collections::HashMap;
// use std::error::Error;
use ndarray::{Array2};
use crate::ms_io::Peak;
use rayon::prelude::*;
use std::time::Instant;

/// Converts a spectrum from a sparse representation of peaks to a dense representation as a vector.
///
/// The sparse representation is a list of peaks, where each peak is a pair of m/z and intensity.
/// The dense representation is a vector of floats, where each index represents a range of m/z values.
/// The m/z ranges are determined by the `bin_width` argument, and the index of each peak is determined
/// by dividing its m/z by the `bin_width` and taking the floor.
///
/// The value at each index is the intensity of the peak at that index, or 0.0 if there is no peak.
///
/// The vector is of length `nbins`, where `nbins` is the smallest integer greater than or equal to
/// `max_mz / bin_width`. The vector is initialized with all elements set to 0.0.
///
/// ### Arguments
///
/// * `peaks` - A slice of `Peak` values representing the sparse spectrum.
/// * `bin_width` - A float representing the width of each m/z range.
/// * `max_mz` - A float representing the maximum m/z of the spectrum.
///
/// ### Returns
///
/// * A `Vec<f32>` representing the dense spectrum.
pub fn spectrum_to_dense_vec(peaks: &[Peak], bin_width: f32, max_mz: f64, minimum_intensity: f64) -> Vec<f32> {
    // Compute number of bins as ceil(max_mz / bin_width)
    let nbins = (max_mz / bin_width as f64).ceil() as usize;
    // Create a Vec<f32> initialized to 0.0
    let mut bins = vec![0.0_f32; nbins];
    for p in peaks {
        // Only consider peaks with intensity >= minimum_intensity
        if p.intensity <= minimum_intensity as f64 {
            continue;
        }
        // Only consider peaks in [0.0, max_mz]
        if (0.0_f64..=max_mz as f64).contains(&p.mz) {
            // Determine which bin this m/z belongs in
            let idx = ((p.mz as f64) / (bin_width as f64)).floor() as usize;
            if idx < nbins {
                bins[idx] = p.intensity as f32;
            }
        }
    }
    bins
}

/// Computes the cosine similarity between two vectors.
/// 
/// The cosine similarity is a measure of similarity between two non-zero vectors of an inner product space.
/// It is calculated as the dot product of the vectors divided by the product of their magnitudes.
/// 
/// ### Arguments
/// 
/// * `v1` - A slice of f32 values representing the first vector.
/// * `v2` - A slice of f32 values representing the second vector.
/// 
/// ### Returns
/// 
/// * A f32 value representing the cosine similarity between the two vectors. The value ranges from 0.0 to 1.0.
/// 
/// ### Notes
/// 
/// * If the vectors are not of the same length, the function returns 0.0.
/// * If either vector has a magnitude of zero, the function returns 0.0.
pub fn cosine_similarity(v1: &[f32], v2: &[f32]) -> f32 {
    // Check if vectors are the same length
    if v1.len() != v2.len() {
        return 0.0;
    }
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

/// Computes the dot product of two vectors.
///
/// The dot product is the sum of the products of the corresponding entries of the two sequences of numbers.
///
/// ### Arguments
///
/// * `v1` - A slice of f32 values representing the first vector.
/// * `v2` - A slice of f32 values representing the second vector.
///
/// ### Returns
///
/// * A f32 value representing the dot product of the two vectors.
///
/// ### Notes
///
/// * The vectors must be of the same length.
pub fn dot(v1: &[f32], v2: &[f32]) -> f32 {
    // assuming v1.len() == v2.len()
    v1.iter().zip(v2).map(|(a, b)| a * b).sum()
}

/// Normalizes a vector in-place.
///
/// The vector is divided by its magnitude, leaving it unchanged if the magnitude is zero.
///
/// ### Arguments
///
/// * `v` - The vector to normalize.
fn normalize(v: &mut [f32]) {
    let norm = v.iter().map(|x| x * x).sum::<f32>().sqrt();
    if norm > 0.0 {
        for x in v {
            *x /= norm;
        }
    }
}

/// Computes a pairwise similarity matrix for spectra based on their binary bit vectors.
///
/// This function takes a map from scan IDs to their binary bit vector 
/// representations, computes the cosine similarity between each pair of spectra, 
/// and returns both the ordered list of scan IDs and the resulting similarity
/// matrix. The matrix is symmetric, with diagonal elements set to 1.0 (indicating 
/// identical spectra).
///
/// ### Arguments
///
/// * `bits_map` - A reference to a HashMap where keys are scan IDs and values are 
///                their corresponding binary bit vector representations.
///
/// ### Returns
///
/// A tuple containing:
/// - A `Vec<String>` of scan IDs in the order they appear in the matrix.
/// - An `Array2<f32>` representing the pairwise similarity matrix, where each 
///   cell contains the cosine similarity score between the corresponding spectra.
pub fn compute_pairwise_similarity_matrix_ndarray(
    bits_map: &HashMap<String, Vec<f32>>,
) -> (Vec<String>, Array2<f32>) {
    let start = Instant::now();
    // 1) Order and sort scan IDs numerically
    let mut scans: Vec<String> = bits_map.keys().cloned().collect();
    scans.sort_by(|a, b| a.parse::<f64>().unwrap().partial_cmp(&b.parse::<f64>().unwrap()).unwrap());
    let elapsed = start.elapsed();
    println!("Sorted scan IDs in {:.2?}", elapsed);
    let start = Instant::now();
    // 2) Cache references to the spectra slices inside an Arc for thread-safe sharing
    use std::sync::Arc;
    let spectra: Arc<Vec<&[f32]>> = Arc::new(scans.iter().map(|id| bits_map[id].as_slice()).collect());

    let n = scans.len();
    let mut mat = Array2::<f32>::zeros((n, n));
    let elapsed = start.elapsed();
    println!("Cached scan IDs in {:.2?}", elapsed);
    let start = Instant::now();
    // 3) Parallel compute upper-triangle similarities
    let entries: Vec<(usize, usize, f32)> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            // Clone Arc for this thread
            let spectra = Arc::clone(&spectra);
            (i + 1..n).into_par_iter().map(move |j| {
                // spectra now is safely shared
                let sim = dot(spectra[i], spectra[j]);
                (i, j, sim)
            })
        })
        .collect();
    let elapsed = start.elapsed();
    println!("Computed cosine similarities in {:.2?}", elapsed);
    let start = Instant::now();
    // 4) Fill matrix: diagonal and symmetrize
    for i in 0..n {
        mat[(i, i)] = 1.0;
    }
    for (i, j, sim) in entries {
        mat[(i, j)] = sim;
        mat[(j, i)] = sim;
    }
    let elapsed = start.elapsed();
    println!("Filled matrix in {:.2?}", elapsed);

    (scans, mat)
}

/// Compute a map from scan IDs to binary bit vector representations of their peaks
/// at the given bin width and maximum m/z.
///
/// The binary bit vector representation is a fixed-size vector of length
/// `ceil(max_mz / bin_width)` where each element is the intensity of the peak that falls in the
/// corresponding m/z range, and 0.0 otherwise. The `spectrum_to_dense_vec`
/// function is used to compute each binary bit vector. The vector is then normalized and added
/// to the map.
///
/// This function is used to preprocess the data before computing the similarity
/// matrix. The result can be passed to `compute_pairwise_similarity_matrix_ndarray`.
/// 
/// ### Arguments
///
/// * `spec_map` - A reference to a HashMap where keys are scan IDs and values are 
///                their corresponding spectrum data.
/// * `bin_width` - The width of each bin in m/z units.
/// * `max_mz` - The maximum m/z value to consider.
///
/// ### Returns
///
/// A `HashMap<String, Vec<f32>>` where keys are scan IDs and values are their 
/// corresponding binary bit vector representations.
pub fn compute_dense_vec_map(
    spec_map: &HashMap<String, Vec<Peak>>,
    bin_width: f32,
    max_mz: f64,
    minimum_intensity: f64
) -> HashMap<String, Vec<f32>> {
    let mut bits_map = HashMap::with_capacity(spec_map.len());
    for (scan, peaks) in spec_map {
        let mut bv = spectrum_to_dense_vec(peaks, bin_width, max_mz, minimum_intensity);
        normalize(&mut bv);
        bits_map.insert(scan.clone(), bv);
    }
    bits_map
}

/// Remove any vector positions (bins) that are zero across all spectra.
///
/// # Arguments
///
/// * `bits_map` - map from scan ID to its dense vector
///
/// # Returns
///
/// A new `HashMap` where each vector has been pruned to only the indices
/// where at least one spectrum had a non-zero intensity.
pub fn prune_unused_bins(
    bits_map: &HashMap<String, Vec<f32>>,
) -> HashMap<String, Vec<f32>> {
    if bits_map.is_empty() {
        return HashMap::new();
    }
    let n = bits_map.values().next().unwrap().len();
    // Determine which indices are used by any vector
    let mut used = vec![false; n];
    for vec in bits_map.values() {
        for (i, &val) in vec.iter().enumerate() {
            if val != 0.0 {
                used[i] = true;
            }
        }
    }
    // Collect active positions
    let active_indices: Vec<usize> = used.iter()
        .enumerate()
        .filter_map(|(i, &u)| if u { Some(i) } else { None })
        .collect();

    // Build pruned map
    let mut pruned = HashMap::with_capacity(bits_map.len());
    for (scan, vec) in bits_map {
        let pruned_vec: Vec<f32> = active_indices
            .iter()
            .map(|&i| vec[i])
            .collect();
        pruned.insert(scan.clone(), pruned_vec);
    }
    pruned
}

/// Set columns to zero if they have a hit in at least `min_frac * n_spec` spectra.
///
/// This is useful for removing background noise from the similarity matrix.
///
/// # Arguments
///
/// * `bits_map` - Mutable map from scan ID to its dense vector
/// * `min_frac` - Minimum fraction of spectra that must contain a hit in the column
///                for it to be considered a background column. e.g. 0.5 for 50%
pub fn prune_background_columns(
    bits_map: &mut HashMap<String, Vec<f32>>,
    min_frac: f64 // e.g. 0.5 for 50%
) {
    let n_spec = bits_map.len();
    if n_spec == 0 { return; }

    // assume all vectors have the same length
    let len = bits_map.values().next().unwrap().len();
    let threshold = (n_spec as f64 * min_frac).ceil() as usize;

    // 1) Count spectra with a “hit” in each column
    let mut counts = vec![0usize; len];
    for vec in bits_map.values() {
        for (j, &v) in vec.iter().enumerate() {
            if v != 0.0 {
                counts[j] += 1;
            }
        }
    }

    // 2) Mark background columns
    let bg_cols: Vec<usize> = counts.iter()
        .enumerate()
        .filter_map(|(j, &c)| if c >= threshold { Some(j) } else { None })
        .collect();

    // 3a) Option A: zero out background columns
    for vec in bits_map.values_mut() {
        for &j in &bg_cols {
            vec[j] = 0.0;
        }
    }

    // 3b) Option B: *drop* background columns entirely
    // (uncomment if you want to shrink vectors)
    /*
    for vec in bits_map.values_mut() {
        // build a new Vec without bg columns
        let mut pruned = Vec::with_capacity(len - bg_cols.len());
        let mut bg_iter = bg_cols.iter().peekable();
        for (j, &v) in vec.iter().enumerate() {
            if Some(&j) == bg_iter.peek() {
                bg_iter.next();
                continue;
            }
            pruned.push(v);
        }
        *vec = pruned;
    }
    */
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::ms_io::import_mzml;

    /// Tests the generation of a map of binary bit vectors from a map of spectra,
    /// and then computes the cosine similarity matrix from the binary bit vectors.
    /// The similarity matrix is then printed to the console.
    /// 
    /// The test uses the sample mzML file `tests\\data\\FeatureFinderCentroided_1_input.mzML`,
    /// and prints out the number of spectra in the file, the maximum m/z, and the
    /// bin width used. The resulting similarity matrix is then printed to the
    /// console.
    #[test]
    fn test_bitvec_and_cosine() -> Result<(), Vec<String>> {
        // Load sample data
        println!("Importing test mzML file...");
        let (map, _) = import_mzml("tests\\data\\FeatureFinderCentroided_1_input.mzML")?;
        println!("Number of Spectra in test file: {}", map.len());
        // Determine the max m/z
        let max_mz = map.values().map(|p| p.iter().map(|p| p.mz).fold(0., f64::max)).fold(0., f64::max);
        println!("Max m/z: {}", max_mz);
        // Determine bin width
        let bin_width = 0.01;
        println!("Bin width: {}", bin_width);
        // Generate map of bin Vectors
        let bits_map = compute_dense_vec_map(&map, bin_width, max_mz, 0.0);
        // Compute cosine similarity matrix
        let sims = compute_pairwise_similarity_matrix_ndarray(&bits_map);
        // Print confirmation if sims is not empty
        if sims.0.len() > 0 {
            println!("Cosine Similarity Matrix Successfully Generated!");
        }
        Ok(())
    }
}
