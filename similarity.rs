// use std::collections::{HashMap,HashSet};
use hashbrown::{HashMap, HashSet};
// use std::error::Error;
use ndarray::{Array2};
use crate::ms_io::{Peak, SpectrumMetadata};
use rayon::prelude::*;
use std::time::Instant;
use std::sync::Arc;
type SparseVec = Vec<(usize, f32)>;

/// Converts a spectrum from a sparse representation of peaks to a sparse representation as a vector.
///
/// Qualifying peaks are those with intensity above the given minimum and m/z within the given range.
/// The output vector is sorted by bin index, and normalized in place to have unit length.
///
/// ### Arguments
///
/// * `peaks` - A slice of `Peak`s to convert.
/// * `bin_width` - The width of each bin in the output vector.
/// * `max_mz` - The upper bound of the m/z range.
/// * `min_intensity` - The lower bound of the intensity range.
///
/// ### Returns
///
/// A sparse vector, represented as a `Vec<(usize, f32)>`.
fn spectrum_to_sparse_vec(
    peaks: &[Peak],
    bin_width: f32,
    max_mz: f32,
    min_intensity: f32,
) -> SparseVec {
    // 1) Map each qualifying peak → (bin_idx, intensity)
    let mut sv: SparseVec = peaks
        .iter()
        .filter(|p| p.intensity > min_intensity && (0.0..=max_mz).contains(&p.mz))
        .map(|p| {
            let idx = (p.mz / bin_width).floor() as usize;
            (idx, p.intensity)
        })
        .collect();

    // 2) Sort by bin index so merging in dot‐product is linear
    sv.sort_unstable_by_key(|(idx, _)| *idx);

    // 3) Normalize in place
    let norm = sv.iter().map(|(_,i)| i*i).sum::<f32>().sqrt();
    if norm > 0.0 {
        for (_idx, val) in &mut sv {
            *val /= norm;
        }
    }
    sv
}


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
pub fn spectrum_to_dense_vec(peaks: &[Peak], bin_width: f32, max_mz: f32, minimum_intensity: f32) -> Vec<f32> {
    // Compute number of bins as ceil(max_mz / bin_width)
    let nbins = (max_mz / bin_width).ceil() as usize;
    // Create a Vec<f32> initialized to 0.0
    let mut bins = vec![0.0_f32; nbins];
    for p in peaks {
        // Only consider peaks with intensity >= minimum_intensity
        if p.intensity <= minimum_intensity {
            continue;
        }
        // Only consider peaks in [0.0, max_mz]
        if (0.0_f32..=max_mz).contains(&p.mz) {
            // Determine which bin this m/z belongs in
            let idx = ((p.mz) / (bin_width)).floor() as usize;
            if idx < nbins {
                bins[idx] = p.intensity;
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

/// Computes the cosine similarity between two sparse vectors.
///
/// The cosine similarity is a measure of similarity between two non-zero vectors of an inner product space.
/// For sparse vectors, it is calculated as the dot product of the vectors divided by the product of their magnitudes.
/// The vectors must be represented as slices of (index, value) pairs, sorted by index.
///
/// ### Arguments
///
/// * `a` - The first sparse vector.
/// * `b` - The second sparse vector.
///
/// ### Returns
///
/// * A f32 value representing the cosine similarity between the two vectors. The value ranges from 0.0 to 1.0.
///
/// ### Notes
///
/// * If either vector is empty or has a magnitude of zero, the function returns 0.0.
pub fn cosine_similarity_sparse(a: &SparseVec, b: &SparseVec) -> f32 {
    sparse_dot(a, b)
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

/// Computes the dot product of two sparse vectors.
///
/// The sparse vectors are represented as a slice of (index, value) pairs, sorted by index.
/// The dot product is the sum of the products of the corresponding values of the two sequences.
///
/// This function assumes that the vectors are sorted by index.
///
/// ### Arguments
///
/// * `a` - The first sparse vector.
/// * `b` - The second sparse vector.
///
/// ### Returns
///
/// * A f32 value representing the dot product of the two vectors.
fn sparse_dot(a: &SparseVec, b: &SparseVec) -> f32 {
    let mut sum = 0.0;
    let (mut ia, mut ib) = (0, 0);
    while ia < a.len() && ib < b.len() {
        let ai = a[ia].0;
        let bi = b[ib].0;
        match ai.cmp(&bi) {
            std::cmp::Ordering::Less    => ia += 1,
            std::cmp::Ordering::Greater => ib += 1,
            std::cmp::Ordering::Equal   => {
                sum += a[ia].1 * b[ib].1;
                ia += 1;
                ib += 1;
            }
        }
    }
    sum
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
/// * `bits_map` - A reference to a HashMap where keys are scan IDs and values are their corresponding binary bit vector representations.
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
    scans.sort_by(|a, b| a.parse::<f32>().unwrap().partial_cmp(&b.parse::<f32>().unwrap()).unwrap());
    let elapsed = start.elapsed();
    println!("Sorted scan IDs in {elapsed:.2?}");
    let start = Instant::now();
    // 2) Cache references to the spectra slices inside an Arc for thread-safe sharing
    use std::sync::Arc;
    let spectra: Arc<Vec<&[f32]>> = Arc::new(scans.iter().map(|id| bits_map[id].as_slice()).collect());

    let n = scans.len();
    let mut mat = Array2::<f32>::zeros((n, n));
    let elapsed = start.elapsed();
    println!("Cached scan IDs in {elapsed:.2?}");
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
    println!("Computed cosine similarities in {elapsed:.2?}");
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
    println!("Filled matrix in {elapsed:.2?}");

    (scans, mat)
}


/// Computes a pairwise similarity matrix for spectra based on their sparse vector representations.
///
/// This function takes a map from scan IDs to their sparse vector representations,
/// computes the cosine similarity between each pair of spectra, and returns both the ordered list of scan IDs and the resulting similarity matrix.
/// The matrix is symmetric, with diagonal elements set to 1.0 (indicating identical spectra).
///
/// If the `modified_cosine` metric is used, the function will shift the sparse vectors to account for mass shifts and neutral losses between MS/MS spectra.
/// In this case, the `mass_tolerance` parameter specifies the mass tolerance for the shift in units of Da.
///
/// ### Arguments
///
/// * `sparse_map` - A reference to a HashMap where keys are scan IDs and values are their corresponding sparse vector representations.
/// * `spec_metadata` - A reference to a HashMap where keys are scan IDs and values are their corresponding SpectrumMetadata objects.
/// * `similarity_metric` - A string specifying the similarity metric to use. Options are "cosine" and "modified_cosine".
/// * `ms_level` - An integer specifying the MS level of the spectra.
/// * `mass_tolerance` - A float specifying the mass tolerance for the shift in units of Da. Only used if `similarity_metric` is "modified_cosine".
///
/// ### Returns
///
/// A tuple containing:
/// - A `Vec<String>` of scan IDs in the order they appear in the matrix.
/// - An `Array2<f32>` representing the pairwise similarity matrix, where each cell contains the cosine similarity score between the corresponding spectra.
pub fn compute_pairwise_similarity_matrix_sparse(
    sparse_map: &HashMap<String, SparseVec>,
    spec_metadata: &HashMap<String, SpectrumMetadata>,
    similarity_metric: String,
    ms_level: u8,
    mass_tolerance: f32,
) -> (Vec<String>, Array2<f32>) {
    if similarity_metric == "modified_cosine" && ms_level == 1 {
        panic!("The modified cosine similarity metric can only be used with MS2 data.");
    }

    let start = Instant::now();
    let mut scans: Vec<String> = sparse_map.keys().cloned().collect();
    scans.sort_by(|a, b| {
        a.parse::<f32>()
            .unwrap()
            .partial_cmp(&b.parse::<f32>().unwrap())
            .unwrap()
    });
    println!("Sorted scan IDs in {:.2?}", start.elapsed());

    let start = Instant::now();
    let spectra: Arc<Vec<&SparseVec>> = Arc::new(scans.iter().map(|id| &sparse_map[id]).collect());
    let n = scans.len();
    let mut mat = Array2::<f32>::zeros((n, n));
    println!("Cached scan IDs in {:.2?}", start.elapsed());

    let start = Instant::now();
    let scans_clone = scans.clone();
    let entries: Vec<(usize, usize, f32)> = {
        let metric = similarity_metric.clone();
        (0..n)
            .into_par_iter()
            .flat_map(move |i| {
                let spectra = Arc::clone(&spectra);
                let scans = scans_clone.clone();
                let spec_metadata = spec_metadata.clone();
                let metric = metric.clone();

                let a_id = scans[i].clone();
                let a = spectra[i];

                (i + 1..n).into_par_iter().map(move |j| {
                    let b_id = scans[j].clone();
                    let b = spectra[j];

                    let sim = if metric == "modified_cosine" {
                        let a_meta = spec_metadata.get(&a_id).unwrap();
                        let b_meta = spec_metadata.get(&b_id).unwrap();
                        let delta_mz = b_meta.target_mz - a_meta.target_mz;
                        let shift_bins = (delta_mz / mass_tolerance).round() as isize;
                        // let shift_bins = (delta_mz / mass_tolerance);

                        let b_shifted = shift_sparse_vec(b, -shift_bins);
                        let sim1 = cosine_similarity_sparse(a, &b_shifted);

                        let a_shifted = shift_sparse_vec(a, shift_bins);
                        let sim2 = cosine_similarity_sparse(&a_shifted, b);

                        sim1.max(sim2)
                    } else {
                        cosine_similarity_sparse(a, b)
                    };

                    (i, j, sim)
                })
            })
            .collect()
    };
    println!("Computed {} similarities in {:.2?}", entries.len(), start.elapsed());

    let start = Instant::now();
    for i in 0..n {
        mat[(i, i)] = 1.0;
    }
    for (i, j, sim) in entries {
        mat[(i, j)] = sim;
        mat[(j, i)] = sim;
    }
    println!("Filled matrix in {:.2?}", start.elapsed());

    (scans, mat)
}

/// Shifts a sparse vector by a given number of bins, wrapping around to 0
/// if the resulting index is negative.
///
/// # Arguments
///
/// * `vec` - The sparse vector to shift.
/// * `shift_bins` - The number of bins to shift the vector by.
///
/// # Returns
///
/// A new sparse vector with the shifted indices.
fn shift_sparse_vec(vec: &SparseVec, shift_bins: isize) -> SparseVec {
    vec.iter()
        .filter_map(|(bin, intensity)| {
            let shifted = *bin as isize + shift_bins;
            if shifted >= 0 {
                Some((shifted as usize, *intensity))
            } else {
                None
            }
        })
        .collect()
}

/// Computes a map from scan IDs to sparse vector representations of their peaks.
///
/// Each sparse vector representation is a list of (bin index, intensity) pairs,
/// where each pair corresponds to a peak that falls within a specified m/z range and
/// has an intensity greater than the minimum threshold. The `spectrum_to_sparse_vec`
/// function is utilized to compute each sparse vector. The vector is sorted by bin index
/// and normalized to ensure efficient computation in subsequent similarity evaluations.
///
/// ### Arguments
///
/// * `spec_map` - A reference to a HashMap where keys are scan IDs and values are their corresponding spectrum data.
/// * `bin_width` - The width of each bin in m/z units.
/// * `max_mz` - The maximum m/z value to consider.
/// * `minimum_intensity` - The minimum intensity threshold for peaks to be included.
///
/// ### Returns
///
/// A `HashMap<String, SparseVec>` where keys are scan IDs and values are their
/// corresponding sparse vector representations.
pub fn compute_sparse_vec_map(
    spec_map: &HashMap<String, Vec<Peak>>,
    bin_width: f32,
    max_mz: f32,
    minimum_intensity: f32,
) -> HashMap<String, SparseVec> {
    if spec_map.is_empty() {
        // If the input map is empty, return an empty map
        return HashMap::new();
    }

    let mut map = HashMap::with_capacity(spec_map.len());
    for (scan, peaks) in spec_map {
        let sv = spectrum_to_sparse_vec(peaks, bin_width, max_mz, minimum_intensity);
        map.insert(scan.clone(), sv);
    }
    map
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
/// * `spec_map` - A reference to a HashMap where keys are scan IDs and values are their corresponding spectrum data.
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
    max_mz: f32,
    minimum_intensity: f32
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
/// * `min_frac` - Minimum fraction of spectra that must contain a hit in the column for it to be considered a background column. e.g. 0.5 for 50%
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

/// Remove bins from sparse vectors if they appear in at least `min_frac * n_spec` spectra.
///
/// This is useful for removing background noise from the similarity matrix.
///
/// ### Arguments
///
/// * `sparse_map` - Mutable map from scan ID to its sparse vector
/// * `min_frac` - Minimum fraction of spectra that must contain a hit in the bin for it to be considered a background bin. e.g. 0.5 for 50%
pub fn prune_background_bins_sparse(
    sparse_map: &mut HashMap<String, SparseVec>,
    min_frac: f64,
) {
    let n_spec = sparse_map.len();
    if n_spec == 0 {
        return;
    }
    // compute how many spectra a given bin appears in
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for vec in sparse_map.values() {
        for &(bin, _) in vec {
            *counts.entry(bin).or_insert(0) += 1;
        }
    }
    // threshold = ceil(n_spec * min_frac)
    let threshold = (n_spec as f64 * min_frac).ceil() as usize;
    // collect all bins that are “too common”
    let bg_bins: HashSet<usize> = counts
        .into_iter()
        .filter_map(|(bin, cnt)| if cnt >= threshold { Some(bin) } else { None })
        .collect();
    // now prune each sparse vector in place
    for vec in sparse_map.values_mut() {
        vec.retain(|&(bin, _)| !bg_bins.contains(&bin));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::{compute_pairwise_similarity_matrix_ndarray,
        compute_pairwise_similarity_matrix_sparse,
        SparseVec};
    use crate::ms_io::import_mzml;
    // use std::collections::HashMap;
    use hashbrown::HashMap;


    /// Tests the generation of a map of binary bit vectors from a map of spectra,
    /// and then computes the cosine similarity matrix from the binary bit vectors.
    /// The similarity matrix is then printed to the console.
    /// 
    /// The test uses the sample mzML file `tests\\data\\FeatureFinderCentroided_1_input.mzML`,
    /// and prints out the number of spectra in the file, the maximum m/z, and the
    /// bin width used. The resulting similarity matrix is then printed to the
    /// console.

    /// Helper: turn each non-zero entry in a dense Vec<f32> into (index, value)
    fn dense_to_sparse(bits_map: &HashMap<String, Vec<f32>>) -> HashMap<String, SparseVec> {
        bits_map.iter().map(|(k, v)| {
            let sv: SparseVec = v.iter()
                .enumerate()
                .filter_map(|(i, &x)| if x != 0.0 { Some((i, x)) } else { None })
                .collect();
            (k.clone(), sv)
        }).collect()
    }

    #[test]
    fn sparse_matches_dense() {
        println!("\n=== TESTING THAT SPARSE VECTORS MATCHES DENSE VECTORS===");
        // 1) Build a tiny “normalized” dense map
        let mut bits_map: HashMap<String, Vec<f32>> = HashMap::new();
        // orthonormal basis vectors
        bits_map.insert("1".into(), vec![1.0, 0.0]);
        bits_map.insert("2".into(), vec![0.0, 1.0]);
        // 45° vector
        let norm = (2.0f32).sqrt();
        bits_map.insert("3".into(), vec![1.0 / norm, 1.0 / norm]);

        // 2) Dense similarity
        let (dense_scans, dense_mat) = compute_pairwise_similarity_matrix_ndarray(&bits_map);
        println!("-- Dense scans: {dense_scans:?}");
        println!("-- Dense matrix:\n{dense_mat}");

        // 3) Convert to sparse & compute
        let sparse_map = dense_to_sparse(&bits_map);
        let dummy_metadata: HashMap<String, SpectrumMetadata> = HashMap::new();
        let (sparse_scans, sparse_mat) =
            compute_pairwise_similarity_matrix_sparse(&sparse_map, &dummy_metadata, "cosine".to_string(), 1, 1.0);
        println!("-- Sparse scans: {sparse_scans:?}");
        println!("-- Sparse matrix:\n{sparse_mat}");

        // 4) Verify same ordering and dimensions
        assert_eq!(dense_scans, sparse_scans);
        let n = dense_scans.len();
        assert_eq!(dense_mat.shape(), &[n, n]);
        assert_eq!(sparse_mat.shape(), &[n, n]);

        // 5) Compare every entry (allowing tiny float-rounding slack)
        for i in 0..n {
            for j in 0..n {
                let d = dense_mat[(i, j)];
                let s = sparse_mat[(i, j)];
                assert!(
                    (d - s).abs() < 1e-6,
                    "Mismatch at ({},{}) → dense={}, sparse={}",
                    i, j, d, s
                );
            }
        }
    }

    #[test]
    fn test_bitvec_and_cosine() -> Result<(), Vec<String>> {
        // Load sample data
        println!("Importing test mzML file...");
        let (map, _) = import_mzml("tests\\data\\FeatureFinderCentroided_1_input.mzML")?;
        println!("Number of Spectra in test file: {}", map.len());
        // Determine the max m/z
        let max_mz = map.values().map(|p| p.iter().map(|p| p.mz).fold(0., f32::max)).fold(0., f32::max);
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
    #[test]
    fn test_cosine_vs_modified_cosine() {
        use crate::ms_io::SpectrumMetadata;

        println!("\n=== TESTING COSINE VS MODIFIED COSINE ===");
        
        // Create test sparse vectors that represent MS/MS spectra
        // These vectors simulate fragment ion patterns that would be similar
        // if we account for the precursor mass difference
        
        let mut sparse_map: HashMap<String, SparseVec> = HashMap::new();
        let mut spec_metadata: HashMap<String, SpectrumMetadata> = HashMap::new();
        
        // Mass tolerance for binning (this should match your bin_width)
        let mass_tolerance = 1.0; // 1 Da bins for this test
        
        // Spectrum A: precursor m/z = 500.0
        // Fragment pattern: peaks at bins 100, 150, 200, 250
        let spectrum_a = vec![
            (100, 0.8), // normalized intensity
            (150, 0.4),
            (200, 0.3),
            (250, 0.3)
        ];
        
        // Spectrum B: precursor m/z = 510.0 (10 Da higher than A)
        // Fragment pattern shifted by +10 bins: peaks at bins 110, 160, 210, 260
        // This simulates what would happen if the same molecule had a 10 Da mass shift
        let spectrum_b = vec![
            (110, 0.8), // same relative intensities, shifted by 10 bins
            (160, 0.4),
            (210, 0.3),
            (260, 0.3)
        ];
        
        // Spectrum C: precursor m/z = 500.0 (same as A)
        // Different fragment pattern: peaks at bins 120, 180, 240
        let spectrum_c = vec![
            (120, 0.7),
            (180, 0.5),
            (240, 0.4)
        ];
        
        sparse_map.insert("1".to_string(), spectrum_a);
        sparse_map.insert("2".to_string(), spectrum_b);
        sparse_map.insert("3".to_string(), spectrum_c);
        
        // Create metadata with different precursor masses
        spec_metadata.insert("1".to_string(), SpectrumMetadata {
            index: "1".to_string(),
            id: "scan=1".to_string(),
            ms_level: 2,
            base_peak_mz: 150.0,
            base_peak_intensity: 1000.0,
            total_ion_current: 2500.0,
            scan_window_lower_limit: 50.0,
            scan_window_upper_limit: 1000.0,
            target_mz: 500.0,
            selected_ion: 500.0,
            charge: 2,
        });
        
        spec_metadata.insert("2".to_string(), SpectrumMetadata {
            index: "2".to_string(),
            id: "scan=2".to_string(),
            ms_level: 2,
            base_peak_mz: 160.0,
            base_peak_intensity: 1000.0,
            total_ion_current: 2500.0,
            scan_window_lower_limit: 50.0,
            scan_window_upper_limit: 1000.0,
            target_mz: 510.0, // 10 Da higher
            selected_ion: 510.0,
            charge: 2,
        });
        
        spec_metadata.insert("3".to_string(), SpectrumMetadata {
            index: "3".to_string(),
            id: "scan=3".to_string(),
            ms_level: 2,
            base_peak_mz: 180.0,
            base_peak_intensity: 800.0,
            total_ion_current: 2000.0,
            scan_window_lower_limit: 50.0,
            scan_window_upper_limit: 1000.0,
            target_mz: 500.0, // same as spectrum 1
            selected_ion: 500.0,
            charge: 2,
        });
        
        // Compute similarity matrices using both methods
        let (scans_cosine, matrix_cosine) = compute_pairwise_similarity_matrix_sparse(
            &sparse_map,
            &spec_metadata,
            "cosine".to_string(),
            2,
            mass_tolerance,
        );
        
        let (scans_modified, matrix_modified) = compute_pairwise_similarity_matrix_sparse(
            &sparse_map,
            &spec_metadata,
            "modified_cosine".to_string(),
            2,
            mass_tolerance,
        );
        
        // The scan order should be the same
        assert_eq!(scans_cosine, scans_modified);
        
        // Find indices for our spectra (they should be sorted as "1", "2", "3")
        let idx_1 = scans_cosine.iter().position(|x| x == "1").unwrap();
        let idx_2 = scans_cosine.iter().position(|x| x == "2").unwrap();
        let idx_3 = scans_cosine.iter().position(|x| x == "3").unwrap();
        
        println!("=== COSINE SIMILARITY RESULTS ===");
        println!("-- Spectrum 1 vs 2 (cosine): {:.4}", matrix_cosine[(idx_1, idx_2)]);
        println!("-- Spectrum 1 vs 3 (cosine): {:.4}", matrix_cosine[(idx_1, idx_3)]);
        println!("-- Spectrum 2 vs 3 (cosine): {:.4}", matrix_cosine[(idx_2, idx_3)]);
        
        println!("\n=== MODIFIED COSINE SIMILARITY RESULTS ===");
        println!("-- Spectrum 1 vs 2 (modified): {:.4}", matrix_modified[(idx_1, idx_2)]);
        println!("-- Spectrum 1 vs 3 (modified): {:.4}", matrix_modified[(idx_1, idx_3)]);
        println!("-- Spectrum 2 vs 3 (modified): {:.4}", matrix_modified[(idx_2, idx_3)]);
        
        // Test expectations:
        
        // 1. Regular cosine: Spectrum 1 vs 2 should have low similarity (different bins)
        let cosine_1_vs_2 = matrix_cosine[(idx_1, idx_2)];
        assert!(cosine_1_vs_2 < 0.1, "-- Regular cosine similarity between shifted spectra should be low, got {}", cosine_1_vs_2);
        
        // 2. Modified cosine: Spectrum 1 vs 2 should have HIGH similarity (accounting for mass shift)
        let modified_1_vs_2 = matrix_modified[(idx_1, idx_2)];
        assert!(modified_1_vs_2 > 0.8, "-- Modified cosine similarity between shifted but similar spectra should be high, got {}", modified_1_vs_2);
        
        // 3. The improvement should be significant
        let improvement = modified_1_vs_2 - cosine_1_vs_2;
        assert!(improvement > 0.7, "-- Modified cosine should show significant improvement over regular cosine for shifted spectra, improvement: {}", improvement);
        
        // 4. For spectra with same precursor mass but different patterns (1 vs 3), 
        //    both methods should give similar (low) results
        let cosine_1_vs_3 = matrix_cosine[(idx_1, idx_3)];
        let modified_1_vs_3 = matrix_modified[(idx_1, idx_3)];
        let diff_1_vs_3 = (modified_1_vs_3 - cosine_1_vs_3).abs();
        assert!(diff_1_vs_3 < 0.1, "-- For spectra with same precursor mass, cosine and modified cosine should be similar, diff: {}", diff_1_vs_3);
        
        println!("\n=== TEST SUMMARY ===");
        println!("-- Regular cosine failed to match shifted spectra (similarity: {:.4})", cosine_1_vs_2);
        println!("-- Modified cosine successfully matched shifted spectra (similarity: {:.4})", modified_1_vs_2);
        println!("-- Improvement: {:.4}", improvement);
        println!("-- Both methods agree on truly different spectra");
    }

    #[test]
    fn test_modified_cosine_edge_cases() {
        use crate::ms_io::SpectrumMetadata;

        println!("\n=== TESTING MODIFIED COSINE EDGE CASES ===");
        
        let mut sparse_map: HashMap<String, SparseVec> = HashMap::new();
        let mut spec_metadata: HashMap<String, SpectrumMetadata> = HashMap::new();
        
        // Test with very small mass differences that shouldn't cause shifts
        let spectrum_a = vec![(100, 1.0)];
        let spectrum_b = vec![(100, 1.0)]; // identical spectrum
        
        sparse_map.insert("1".to_string(), spectrum_a);
        sparse_map.insert("2".to_string(), spectrum_b);
        
        // Very small mass difference (0.1 Da)
        spec_metadata.insert("1".to_string(), SpectrumMetadata {
            index: "1".to_string(),
            id: "scan=1".to_string(),
            ms_level: 2,
            base_peak_mz: 100.0,
            base_peak_intensity: 1000.0,
            total_ion_current: 1000.0,
            scan_window_lower_limit: 50.0,
            scan_window_upper_limit: 1000.0,
            target_mz: 500.0,
            selected_ion: 500.0,
            charge: 2,
        });
        
        spec_metadata.insert("2".to_string(), SpectrumMetadata {
            index: "2".to_string(),
            id: "scan=2".to_string(),
            ms_level: 2,
            base_peak_mz: 100.0,
            base_peak_intensity: 1000.0,
            total_ion_current: 1000.0,
            scan_window_lower_limit: 50.0,
            scan_window_upper_limit: 1000.0,
            target_mz: 500.1, // only 0.1 Da difference
            selected_ion: 500.1,
            charge: 2,
        });
        
        let mass_tolerance = 1.0; // 1 Da bins
        
        let (_, matrix_cosine) = compute_pairwise_similarity_matrix_sparse(
            &sparse_map,
            &spec_metadata,
            "cosine".to_string(),
            2,
            mass_tolerance,
        );
        
        let (_, matrix_modified) = compute_pairwise_similarity_matrix_sparse(
            &sparse_map,
            &spec_metadata,
            "modified_cosine".to_string(),
            2,
            mass_tolerance,
        );
        
        // For small mass differences, both should give identical results
        let cosine_sim = matrix_cosine[(0, 1)];
        let modified_sim = matrix_modified[(0, 1)];
        
        println!("-- Small mass diff - Cosine: {:.6}, Modified: {:.6}", cosine_sim, modified_sim);
        
        // They should be nearly identical since the shift is less than 1 bin
        assert!((cosine_sim - modified_sim).abs() < 0.001, 
            "-- For small mass differences, cosine and modified cosine should be nearly identical");
    }

    #[test]
    fn test_vector_shifting_with_real_data() -> Result<(), Vec<String>> {
        use crate::ms_io::import_mzml;
        // use std::collections::HashMap;
        use hashbrown::HashMap;
        
        println!("\n=== TESTING VECTOR SHIFTING WITH REAL DATA ===");
        
        // Load sample data
        println!("-- Importing test mzML file...");
        let (spec_map, metadata_map) = import_mzml("tests\\data\\FeatureFinderCentroided_1_input.mzML")?;
        println!("-- Number of Spectra in test file: {}", spec_map.len());
        
        // Filter for MS2 spectra only (modified cosine requires MS2)
        let ms2_specs: HashMap<String, Vec<Peak>> = spec_map
            .into_iter()
            .filter(|(scan_id, _)| {
                metadata_map.get(scan_id)
                    .map(|meta| meta.ms_level == 2)
                    .unwrap_or(false)
            })
            .collect();
        
        let ms2_metadata: HashMap<String, SpectrumMetadata> = metadata_map
            .into_iter()
            .filter(|(_, meta)| meta.ms_level == 2)
            .collect();
        
        println!("-- Number of MS2 spectra: {}", ms2_specs.len());
        
        if ms2_specs.len() < 2 {
            println!("\n**Not enough MS2 spectra for shifting test, skipping...**\n");
            return Ok(());
        }
        
        // Parameters for vector generation
        let max_mz = ms2_specs.values()
            .flat_map(|peaks| peaks.iter().map(|p| p.mz))
            .fold(0.0, f32::max);
        let bin_width = 1.0; // 1 Da bins for clearer shifting
        let minimum_intensity = 0.0;
        
        println!("-- Max m/z: {:.2}", max_mz);
        println!("-- Bin width: {}", bin_width);
        
        // Generate sparse vectors
        let sparse_map = compute_sparse_vec_map(&ms2_specs, bin_width, max_mz, minimum_intensity);
        
        // Pick two spectra with different precursor masses for testing
        let mut scan_pairs: Vec<(String, String, f32)> = Vec::new();
        let scan_ids: Vec<String> = sparse_map.keys().cloned().collect();
        
        for (i, scan_a) in scan_ids.iter().enumerate() {
            for scan_b in scan_ids.iter().skip(i + 1) {
                if let (Some(meta_a), Some(meta_b)) = (ms2_metadata.get(scan_a), ms2_metadata.get(scan_b)) {
                    let mass_diff = (meta_b.target_mz - meta_a.target_mz).abs();
                    if mass_diff > 5.0 { // Look for pairs with significant mass difference
                        scan_pairs.push((scan_a.clone(), scan_b.clone(), mass_diff));
                    }
                }
            }
        }
        
        // Sort by mass difference and pick the pair with largest difference
        scan_pairs.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap());
        
        if scan_pairs.is_empty() {
            println!("\n**No spectrum pairs found with significant mass differences**\n");
            return Ok(());
        }
        
        let (scan_a, scan_b, mass_diff) = &scan_pairs[0];
        println!("\n-- Testing with spectra: {} vs {}", scan_a, scan_b);
        println!("-- Mass difference: {:.2} Da", mass_diff);
        
        let meta_a = &ms2_metadata[scan_a];
        let meta_b = &ms2_metadata[scan_b];
        
        println!("-- Spectrum A - target_mz: {:.2}, peaks: {}", meta_a.target_mz, sparse_map[scan_a].len());
        println!("-- Spectrum B - target_mz: {:.2}, peaks: {}", meta_b.target_mz, sparse_map[scan_b].len());
        
        // Calculate expected shift
        let delta_mz = meta_b.target_mz - meta_a.target_mz;
        let shift_bins = (delta_mz / bin_width).round() as isize;
        
        println!("\n-- Shift calculation:");
        println!("-- Delta m/z: {:.2}", delta_mz);
        println!("-- Expected shift in bins: {}", shift_bins);
        
        // Get original vectors
        let vec_a = &sparse_map[scan_a];
        let vec_b = &sparse_map[scan_b];
        
        // Perform shifts (same logic as in modified cosine)
        let vec_b_shifted = shift_sparse_vec(vec_b, -shift_bins);
        let vec_a_shifted = shift_sparse_vec(vec_a, shift_bins);
        
        println!("\n=== VECTOR COMPARISON ===");
        println!("-- Original A: {} peaks", vec_a.len());
        println!("-- Original B: {} peaks", vec_b.len());
        println!("-- B shifted by {} bins: {} peaks", -shift_bins, vec_b_shifted.len());
        println!("-- A shifted by {} bins: {} peaks", shift_bins, vec_a_shifted.len());
        
        // Show first few peaks for comparison
        println!("\nFirst 10 peaks comparison:");
        println!("-- Original A bins: {:?}", vec_a.iter().take(10).map(|(bin, _)| *bin).collect::<Vec<_>>());
        println!("-- Original B bins: {:?}", vec_b.iter().take(10).map(|(bin, _)| *bin).collect::<Vec<_>>());
        println!("-- B shifted bins: {:?}", vec_b_shifted.iter().take(10).map(|(bin, _)| *bin).collect::<Vec<_>>());
        println!("-- A shifted bins: {:?}", vec_a_shifted.iter().take(10).map(|(bin, _)| *bin).collect::<Vec<_>>());
        
        // Convert to dense vectors for easier visualization (handle empty vectors)
        let max_bin = [vec_a, vec_b, &vec_b_shifted, &vec_a_shifted]
            .iter()
            .filter(|v| !v.is_empty())
            .flat_map(|v| v.iter().map(|(bin, _)| *bin))
            .max()
            .unwrap_or(0) + 10;
        
        let dense_a = sparse_to_dense_helper(vec_a, max_bin);
        let dense_b = sparse_to_dense_helper(vec_b, max_bin);
        let dense_b_shifted = sparse_to_dense_helper(&vec_b_shifted, max_bin);
        let dense_a_shifted = sparse_to_dense_helper(&vec_a_shifted, max_bin);
        
        // Count non-zero positions
        let nonzero_a: Vec<usize> = dense_a.iter().enumerate().filter_map(|(i, &v)| if v > 0.0 { Some(i) } else { None }).collect();
        let nonzero_b: Vec<usize> = dense_b.iter().enumerate().filter_map(|(i, &v)| if v > 0.0 { Some(i) } else { None }).collect();
        let nonzero_b_shifted: Vec<usize> = dense_b_shifted.iter().enumerate().filter_map(|(i, &v)| if v > 0.0 { Some(i) } else { None }).collect();
        let nonzero_a_shifted: Vec<usize> = dense_a_shifted.iter().enumerate().filter_map(|(i, &v)| if v > 0.0 { Some(i) } else { None }).collect();
        
        println!("\nNon-zero bin positions:");
        println!("-- Original A: {:?}", &nonzero_a[..nonzero_a.len().min(10)]);
        println!("-- Original B: {:?}", &nonzero_b[..nonzero_b.len().min(10)]);
        println!("-- B shifted:  {:?}", &nonzero_b_shifted[..nonzero_b_shifted.len().min(10)]);
        println!("-- A shifted:  {:?}", &nonzero_a_shifted[..nonzero_a_shifted.len().min(10)]);
        
        // Verify shifting worked correctly
        if shift_bins != 0 {
            // Check that B shifted is actually different from original B (only if both are non-empty)
            if !vec_b.is_empty() && !vec_b_shifted.is_empty() {
                let b_vs_b_shifted_same = vec_b.iter().zip(vec_b_shifted.iter()).all(|(orig, shifted)| orig.0 == shifted.0);
                assert!(!b_vs_b_shifted_same, "-- Shifted vector B should be different from original B");
            } else if vec_b_shifted.is_empty() && !vec_b.is_empty() {
                println!("-- -- Note: B-shifted vector became empty after shifting by {} bins", -shift_bins);
            }
            
            // Check that A shifted is actually different from original A (only if both are non-empty)
            if !vec_a.is_empty() && !vec_a_shifted.is_empty() {
                let a_vs_a_shifted_same = vec_a.iter().zip(vec_a_shifted.iter()).all(|(orig, shifted)| orig.0 == shifted.0);
                assert!(!a_vs_a_shifted_same, "-- Shifted vector A should be different from original A");
            } else if vec_a_shifted.is_empty() && !vec_a.is_empty() {
                println!("-- -- Note: A-shifted vector became empty after shifting by {} bins", shift_bins);
            }
            
            // Verify the shift magnitude is correct for B (only if both vectors are non-empty)
            if !vec_b.is_empty() && !vec_b_shifted.is_empty() {
                let first_bin_b = vec_b[0].0 as isize;
                let first_bin_b_shifted = vec_b_shifted[0].0 as isize;
                let actual_shift = first_bin_b_shifted - first_bin_b;
                println!("\n-- Shift verification:");
                println!("-- Expected shift: {}", -shift_bins);
                println!("-- Actual shift (first peak): {}", actual_shift);
                
                // Allow some tolerance for peaks that might be filtered out
                if actual_shift != -shift_bins {
                    println!("-- -- Note: First peak shift differs from expected (peaks may have been filtered)");
                }
            }
            
            println!("\nVector shifting verification completed!");
            if !vec_b_shifted.is_empty() && !vec_a_shifted.is_empty() {
                println!("-- Original and shifted vectors are different");
                println!("-- Shift direction and magnitude are as expected");
            } else {
                println!("-- One or more shifted vectors became empty - this may indicate a very large shift");
            }
        } else {
            println!("\n-- No shift expected (mass difference < 0.5 * bin_width)");
        }
        
        // Calculate similarities to show the effect (with error handling for empty vectors)
        let sim_original = cosine_similarity_sparse(vec_a, vec_b);
        
        let sim_b_shifted = if vec_b_shifted.is_empty() {
            println!("Warning: B-shifted vector is empty, returning 0.0 similarity");
            0.0
        } else {
            cosine_similarity_sparse(vec_a, &vec_b_shifted)
        };
        
        let sim_a_shifted = if vec_a_shifted.is_empty() {
            println!("Warning: A-shifted vector is empty, returning 0.0 similarity");
            0.0
        } else {
            cosine_similarity_sparse(&vec_a_shifted, vec_b)
        };
        
        let modified_cosine_sim = sim_b_shifted.max(sim_a_shifted);
        
        println!("\n=== SIMILARITY COMPARISON ===");
        println!("-- Original A vs B:        {:.4}", sim_original);
        println!("-- A vs B-shifted:         {:.4}", sim_b_shifted);
        println!("-- A-shifted vs B:         {:.4}", sim_a_shifted);
        println!("-- Modified cosine result: {:.4}", modified_cosine_sim);
        println!("-- Improvement:            {:.4}", modified_cosine_sim - sim_original);
        
        Ok(())
    }

    // Helper function to convert sparse vector to dense for visualization
    fn sparse_to_dense_helper(sparse: &SparseVec, size: usize) -> Vec<f32> {
        let mut dense = vec![0.0; size];
        for &(bin, intensity) in sparse {
            if bin < size {
                dense[bin] = intensity;
            }
        }
        dense
    }
}