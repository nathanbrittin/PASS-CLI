// use std::collections::{HashMap,HashSet};
use hashbrown::{HashMap, HashSet};
// use std::error::Error;
use ndarray::{Array2};
use crate::ms_io::{Peak, SpectrumMetadata};
use rayon::prelude::*;
use std::time::Instant;
use std::sync::Arc;
use std::fmt;

type SparseVec = Vec<(usize, f32)>;

// Custom error types for better error handling
#[derive(Debug, Clone)]
pub enum SpectrumProcessingError {
    InvalidParameters(String),
    EmptyInput(String),
    ParseError(String),
    InvalidMsLevel(String),
    MatrixDimensionMismatch(String),
    NumericalError(String),
}

impl fmt::Display for SpectrumProcessingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SpectrumProcessingError::InvalidParameters(msg) => write!(f, "**Invalid parameters: {}**", msg),
            SpectrumProcessingError::EmptyInput(msg) => write!(f, "**Empty input: {}**", msg),
            SpectrumProcessingError::ParseError(msg) => write!(f, "**Parse error: {}**", msg),
            SpectrumProcessingError::InvalidMsLevel(msg) => write!(f, "**nvalid MS level: {}**", msg),
            SpectrumProcessingError::MatrixDimensionMismatch(msg) => write!(f, "**Matrix dimension mismatch: {}**", msg),
            SpectrumProcessingError::NumericalError(msg) => write!(f, "**Numerical error: {}**", msg),
        }
    }
}

impl std::error::Error for SpectrumProcessingError {}

type Result<T> = std::result::Result<T, SpectrumProcessingError>;

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
) -> Result<SparseVec> {
    // Validate parameters
    if bin_width <= 0.0 {
        return Err(SpectrumProcessingError::InvalidParameters(
            "**bin_width must be positive**".to_string()
        ));
    }
    if max_mz <= 0.0 {
        return Err(SpectrumProcessingError::InvalidParameters(
            "**max_mz must be positive**".to_string()
        ));
    }
    if min_intensity < 0.0 {
        return Err(SpectrumProcessingError::InvalidParameters(
            "**min_intensity cannot be negative**".to_string()
        ));
    }

    // 1) Map each qualifying peak → (bin_idx, intensity)
    let mut sv: SparseVec = peaks
        .iter()
        .filter(|p| p.intensity > min_intensity && (0.0..=max_mz).contains(&p.mz))
        .filter_map(|p| {
            if p.mz.is_finite() && p.intensity.is_finite() {
                let idx = (p.mz / bin_width).floor() as usize;
                Some((idx, p.intensity))
            } else {
                // Skip non-finite values
                None
            }
        })
        .collect();

    // 2) Sort by bin index so merging in dot‐product is linear
    sv.sort_unstable_by_key(|(idx, _)| *idx);

    // 3) Normalize in place
    let norm_squared: f32 = sv.iter().map(|(_,i)| i*i).sum();
    if norm_squared <= 0.0 {
        // Return empty vector if no valid peaks or all zero intensities
        return Ok(Vec::new());
    }
    
    let norm = norm_squared.sqrt();
    if !norm.is_finite() {
        return Err(SpectrumProcessingError::NumericalError(
            "**Cannot compute finite norm for spectrum normalization**".to_string()
        ));
    }
    
    for (_idx, val) in &mut sv {
        *val /= norm;
    }
    
    Ok(sv)
}

/// Converts a spectrum from a sparse representation of peaks to a dense representation as a vector.
///
/// The sparse representation is a list of peaks, where each peak is a pair of m/z and intensity.
/// The dense representation is a vector of floats, where each index represents a range of m/z values.
/// The m/z ranges are determined by the `bin_width` argument, and the index of each peak is determined
/// by dividing its m/z by the `bin_width` and taking the floor.
///
/// ### Arguments
///
/// * `peaks` - A slice of `Peak` values representing the sparse spectrum.
/// * `bin_width` - A float representing the width of each m/z range.
/// * `max_mz` - A float representing the maximum m/z of the spectrum.
/// * `minimum_intensity` - A float representing the minimum intensity threshold.
///
/// ### Returns
///
/// * A `Result<Vec<f32>>` representing the dense spectrum or an error.
pub fn spectrum_to_dense_vec(
    peaks: &[Peak], 
    bin_width: f32, 
    max_mz: f32, 
    minimum_intensity: f32
) -> Result<Vec<f32>> {
    // Validate parameters
    if bin_width <= 0.0 {
        return Err(SpectrumProcessingError::InvalidParameters(
            "**bin_width must be positive**".to_string()
        ));
    }
    if max_mz <= 0.0 {
        return Err(SpectrumProcessingError::InvalidParameters(
            "**max_mz must be positive**".to_string()
        ));
    }
    if minimum_intensity < 0.0 {
        return Err(SpectrumProcessingError::InvalidParameters(
            "**minimum_intensity cannot be negative**".to_string()
        ));
    }

    // Compute number of bins as ceil(max_mz / bin_width)
    let nbins_f64 = (max_mz / bin_width).ceil() as f64;
    if nbins_f64 > usize::MAX as f64 {
        return Err(SpectrumProcessingError::InvalidParameters(
            "**Computed number of bins exceeds maximum usize value**".to_string()
        ));
    }
    let nbins = nbins_f64 as usize;
    
    // Check for reasonable memory allocation
    if nbins > 100_000_000 { // 100M bins ~ 400MB for f32
        return Err(SpectrumProcessingError::InvalidParameters(
            format!("**Too many bins requested ({}). Consider increasing bin_width or reducing max_mz**", nbins)
        ));
    }

    // Create a Vec<f32> initialized to 0.0
    let mut bins = vec![0.0_f32; nbins];
    
    for p in peaks {
        // Validate peak data
        if !p.mz.is_finite() || !p.intensity.is_finite() {
            continue; // Skip invalid peaks
        }
        
        // Only consider peaks with intensity >= minimum_intensity
        if p.intensity <= minimum_intensity {
            continue;
        }
        
        // Only consider peaks in [0.0, max_mz]
        if (0.0_f32..=max_mz).contains(&p.mz) {
            // Determine which bin this m/z belongs in
            let idx_f32 = p.mz / bin_width;
            if idx_f32.is_finite() {
                let idx = idx_f32.floor() as usize;
                if idx < nbins {
                    bins[idx] = p.intensity;
                }
            }
        }
    }
    
    Ok(bins)
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
/// * A Result<f32> value representing the cosine similarity between the two vectors or an error.
pub fn cosine_similarity(v1: &[f32], v2: &[f32]) -> Result<f32> {
    // Check if vectors are the same length
    if v1.len() != v2.len() {
        return Err(SpectrumProcessingError::MatrixDimensionMismatch(
            format!("**Vector lengths do not match: {} vs {}**", v1.len(), v2.len())
        ));
    }
    
    if v1.is_empty() {
        return Ok(0.0);
    }
    
    let n = v1.len();
    let mut dot = 0.0;
    let mut norm1 = 0.0;
    let mut norm2 = 0.0;
    
    for i in 0..n {
        if !v1[i].is_finite() || !v2[i].is_finite() {
            return Err(SpectrumProcessingError::NumericalError(
                format!("**Non-finite values found at index {}**", i)
            ));
        }
        dot += v1[i] * v2[i];
        norm1 += v1[i] * v1[i];
        norm2 += v2[i] * v2[i];
    }
    
    if norm1 <= 0.0 || norm2 <= 0.0 {
        return Ok(0.0);
    }
    
    let norm_product = norm1.sqrt() * norm2.sqrt();
    if !norm_product.is_finite() || norm_product == 0.0 {
        return Err(SpectrumProcessingError::NumericalError(
            "**Cannot compute finite norm product**".to_string()
        ));
    }
    
    let similarity = dot / norm_product;
    if !similarity.is_finite() {
        return Err(SpectrumProcessingError::NumericalError(
            "**Computed similarity is not finite**".to_string()
        ));
    }
    
    Ok(similarity.clamp(0.0, 1.0))
}

/// Computes the cosine similarity between two sparse vectors.
pub fn cosine_similarity_sparse(a: &SparseVec, b: &SparseVec) -> Result<f32> {
    sparse_dot(a, b)
}

/// Computes the dot product of two vectors.
pub fn dot(v1: &[f32], v2: &[f32]) -> Result<f32> {
    if v1.len() != v2.len() {
        return Err(SpectrumProcessingError::MatrixDimensionMismatch(
            format!("**Vector lengths do not match: {} vs {}**", v1.len(), v2.len())
        ));
    }
    
    let result: f32 = v1.iter()
        .zip(v2)
        .map(|(a, b)| a * b)
        .sum();
    
    if !result.is_finite() {
        return Err(SpectrumProcessingError::NumericalError(
            "**Dot product result is not finite**".to_string()
        ));
    }
    
    Ok(result)
}

/// Computes the dot product of two sparse vectors.
fn sparse_dot(a: &SparseVec, b: &SparseVec) -> Result<f32> {
    let mut sum = 0.0;
    let (mut ia, mut ib) = (0, 0);
    
    while ia < a.len() && ib < b.len() {
        let ai = a[ia].0;
        let bi = b[ib].0;
        
        match ai.cmp(&bi) {
            std::cmp::Ordering::Less    => ia += 1,
            std::cmp::Ordering::Greater => ib += 1,
            std::cmp::Ordering::Equal   => {
                let product = a[ia].1 * b[ib].1;
                if !product.is_finite() {
                    return Err(SpectrumProcessingError::NumericalError(
                        format!("**Non-finite product at indices {} and {}**", ia, ib)
                    ));
                }
                sum += product;
                ia += 1;
                ib += 1;
            }
        }
    }
    
    if !sum.is_finite() {
        return Err(SpectrumProcessingError::NumericalError(
            "**Sparse dot product result is not finite**".to_string()
        ));
    }
    
    Ok(sum)
}

/// Normalizes a vector in-place.
fn normalize(v: &mut [f32]) -> Result<()> {
    if v.is_empty() {
        return Ok(());
    }
    
    let norm_squared: f32 = v.iter().map(|x| x * x).sum();
    if norm_squared <= 0.0 {
        return Ok(()); // Already normalized (all zeros)
    }
    
    let norm = norm_squared.sqrt();
    if !norm.is_finite() || norm == 0.0 {
        return Err(SpectrumProcessingError::NumericalError(
            "**Cannot compute finite norm for normalization**".to_string()
        ));
    }
    
    for x in v {
        *x /= norm;
        if !x.is_finite() {
            return Err(SpectrumProcessingError::NumericalError(
                "**Non-finite value after normalization**".to_string()
            ));
        }
    }
    
    Ok(())
}

/// Computes a pairwise similarity matrix for spectra based on their binary bit vectors.
pub fn compute_pairwise_similarity_matrix_ndarray(
    bits_map: &HashMap<String, Vec<f32>>,
) -> Result<(Vec<String>, Array2<f32>)> {
    if bits_map.is_empty() {
        return Err(SpectrumProcessingError::EmptyInput(
            "**Input bits_map is empty**".to_string()
        ));
    }

    let start = Instant::now();
    // 1) Order and sort scan IDs numerically
    let mut scans: Vec<String> = bits_map.keys().cloned().collect();
    scans.sort_by(|a, b| {
        let a_parsed = a.parse::<f32>().map_err(|_| {
            SpectrumProcessingError::ParseError(format!("**Cannot parse scan ID '{}' as f32**", a))
        });
        let b_parsed = b.parse::<f32>().map_err(|_| {
            SpectrumProcessingError::ParseError(format!("**Cannot parse scan ID '{}' as f32**", b))
        });
        
        match (a_parsed, b_parsed) {
            (Ok(a_val), Ok(b_val)) => a_val.partial_cmp(&b_val).unwrap_or(std::cmp::Ordering::Equal),
            _ => a.cmp(b), // Fallback to string comparison if parsing fails
        }
    });
    
    let elapsed = start.elapsed();
    println!("||    Sorted scan IDs in {elapsed:.2?}");
    
    // Validate that all vectors have the same length
    let first_len = bits_map.values().next().unwrap().len();
    for (scan_id, vec) in bits_map {
        if vec.len() != first_len {
            return Err(SpectrumProcessingError::MatrixDimensionMismatch(
                format!("**Vector for scan {} has length {} but expected {}**", scan_id, vec.len(), first_len)
            ));
        }
    }
    
    let start = Instant::now();
    // 2) Cache references to the spectra slices inside an Arc for thread-safe sharing
    let spectra: Arc<Vec<&[f32]>> = Arc::new(scans.iter().map(|id| bits_map[id].as_slice()).collect());

    let n = scans.len();
    let mut mat = Array2::<f32>::zeros((n, n));
    let elapsed = start.elapsed();
    println!("||    Cached scan IDs in {elapsed:.2?}");
    
    let start = Instant::now();
    // 3) Parallel compute upper-triangle similarities
    let entries: std::result::Result<Vec<(usize, usize, f32)>, SpectrumProcessingError> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            // Clone Arc for this thread
            let spectra = Arc::clone(&spectra);
            (i + 1..n).into_par_iter().map(move |j| {
                // spectra now is safely shared
                match dot(spectra[i], spectra[j]) {
                    Ok(sim) => Ok((i, j, sim)),
                    Err(e) => Err(e),
                }
            })
        })
        .collect();
    
    let entries = entries?;
    let elapsed = start.elapsed();
    println!("||    Computed cosine similarities in {elapsed:.2?}");
    
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
    println!("||    Filled matrix in {elapsed:.2?}");

    Ok((scans, mat))
}

/// Computes a pairwise similarity matrix for spectra based on their sparse vector representations.
pub fn compute_pairwise_similarity_matrix_sparse(
    sparse_map: &HashMap<String, SparseVec>,
    spec_metadata: &HashMap<String, SpectrumMetadata>,
    similarity_metric: String,
    ms_level: u8,
    mass_tolerance: f32,
) -> Result<(Vec<String>, Array2<f32>)> {
    // Validate inputs
    if sparse_map.is_empty() {
        return Err(SpectrumProcessingError::EmptyInput(
            "**Input sparse_map is empty**".to_string()
        ));
    }
    
    if similarity_metric == "modified-cosine" && ms_level == 1 {
        return Err(SpectrumProcessingError::InvalidMsLevel(
            "**The modified cosine similarity metric can only be used with MS2 data**".to_string()
        ));
    }
    
    if similarity_metric == "modified-cosine" && mass_tolerance <= 0.0 {
        return Err(SpectrumProcessingError::InvalidParameters(
            "**mass_tolerance must be positive for modified-cosine metric**".to_string()
        ));
    }
    
    if similarity_metric != "cosine" && similarity_metric != "modified-cosine" {
        return Err(SpectrumProcessingError::InvalidParameters(
            format!("**Unknown similarity metric: {}. Must be 'cosine' or 'modified-cosine'**", similarity_metric)
        ));
    }

    let start = Instant::now();
    let mut scans: Vec<String> = sparse_map.keys().cloned().collect();
    scans.sort_by(|a, b| {
        let a_parsed = a.parse::<f32>();
        let b_parsed = b.parse::<f32>();
        
        match (a_parsed, b_parsed) {
            (Ok(a_val), Ok(b_val)) => a_val.partial_cmp(&b_val).unwrap_or(std::cmp::Ordering::Equal),
            _ => a.cmp(b), // Fallback to string comparison
        }
    });
    println!("||    Sorted scan IDs in {:.2?}", start.elapsed());

    let start = Instant::now();
    let spectra: Arc<Vec<&SparseVec>> = Arc::new(scans.iter().map(|id| &sparse_map[id]).collect());
    let n = scans.len();
    let mut mat = Array2::<f32>::zeros((n, n));
    println!("||    Cached scan IDs in {:.2?}", start.elapsed());

    let start = Instant::now();
    let scans_clone = scans.clone();
    let entries: std::result::Result<Vec<(usize, usize, f32)>, SpectrumProcessingError> = {
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

                    let sim_result = if metric == "modified-cosine" {
                        let a_meta = spec_metadata.get(&a_id).ok_or_else(|| {
                            SpectrumProcessingError::EmptyInput(
                                format!("**Missing metadata for scan {}**", a_id)
                            )
                        })?;
                        let b_meta = spec_metadata.get(&b_id).ok_or_else(|| {
                            SpectrumProcessingError::EmptyInput(
                                format!("**Missing metadata for scan {}**", b_id)
                            )
                        })?;
                        
                        let delta_mz = b_meta.target_mz - a_meta.target_mz;
                        if !delta_mz.is_finite() {
                            return Err(SpectrumProcessingError::NumericalError(
                                format!("**Invalid target_mz values for scans {} and {}**", a_id, b_id)
                            ));
                        }
                        
                        let shift_bins = (delta_mz / mass_tolerance).round() as isize;

                        let b_shifted = shift_sparse_vec(b, -shift_bins)?;
                        let sim1 = cosine_similarity_sparse(a, &b_shifted)?;

                        let a_shifted = shift_sparse_vec(a, shift_bins)?;
                        let sim2 = cosine_similarity_sparse(&a_shifted, b)?;

                        Ok(sim1.max(sim2))
                    } else {
                        cosine_similarity_sparse(a, b)
                    };

                    match sim_result {
                        Ok(sim) => Ok((i, j, sim)),
                        Err(e) => Err(e),
                    }
                })
            })
            .collect()
    };
    
    let entries = entries?;
    println!("||    Computed {} similarities in {:.2?}", entries.len(), start.elapsed());

    let start = Instant::now();
    for i in 0..n {
        mat[(i, i)] = 1.0;
    }
    for (i, j, sim) in entries {
        mat[(i, j)] = sim;
        mat[(j, i)] = sim;
    }
    println!("||    Filled matrix in {:.2?}", start.elapsed());

    Ok((scans, mat))
}

/// Shifts a sparse vector by a given number of bins, filtering out negative indices.
fn shift_sparse_vec(vec: &SparseVec, shift_bins: isize) -> Result<SparseVec> {
    let result: SparseVec = vec.iter()
        .filter_map(|(bin, intensity)| {
            let shifted = *bin as isize + shift_bins;
            if shifted >= 0 {
                // Check for overflow when converting back to usize
                let shifted_usize = shifted as usize;
                if shifted_usize <= usize::MAX {
                    Some((shifted_usize, *intensity))
                } else {
                    None // Skip if overflow would occur
                }
            } else {
                None
            }
        })
        .collect();
    
    Ok(result)
}

/// Computes a map from scan IDs to sparse vector representations of their peaks.
pub fn compute_sparse_vec_map(
    spec_map: &HashMap<String, Vec<Peak>>,
    bin_width: f32,
    max_mz: f32,
    minimum_intensity: f32,
) -> Result<HashMap<String, SparseVec>> {
    if spec_map.is_empty() {
        return Ok(HashMap::new());
    }

    // Validate parameters
    if bin_width <= 0.0 {
        return Err(SpectrumProcessingError::InvalidParameters(
            "**bin_width must be positive**".to_string()
        ));
    }
    if max_mz <= 0.0 {
        return Err(SpectrumProcessingError::InvalidParameters(
            "**max_mz must be positive**".to_string()
        ));
    }

    let mut map = HashMap::with_capacity(spec_map.len());
    for (scan, peaks) in spec_map {
        match spectrum_to_sparse_vec(peaks, bin_width, max_mz, minimum_intensity) {
            Ok(sv) => {
                map.insert(scan.clone(), sv);
            }
            Err(e) => {
                return Err(SpectrumProcessingError::InvalidParameters(
                    format!("**Error processing scan {}: {}**", scan, e)
                ));
            }
        }
    }
    Ok(map)
}

/// Compute a map from scan IDs to binary bit vector representations of their peaks.
pub fn compute_dense_vec_map(
    spec_map: &HashMap<String, Vec<Peak>>,
    bin_width: f32,
    max_mz: f32,
    minimum_intensity: f32
) -> Result<HashMap<String, Vec<f32>>> {
    if spec_map.is_empty() {
        return Ok(HashMap::new());
    }

    let mut bits_map = HashMap::with_capacity(spec_map.len());
    for (scan, peaks) in spec_map {
        let mut bv = spectrum_to_dense_vec(peaks, bin_width, max_mz, minimum_intensity)?;
        normalize(&mut bv)?;
        bits_map.insert(scan.clone(), bv);
    }
    Ok(bits_map)
}

/// Remove any vector positions (bins) that are zero across all spectra.
pub fn prune_unused_bins(
    bits_map: &HashMap<String, Vec<f32>>,
) -> Result<HashMap<String, Vec<f32>>> {
    if bits_map.is_empty() {
        return Ok(HashMap::new());
    }
    
    let first_len = bits_map.values().next().unwrap().len();
    // Validate all vectors have the same length
    for (scan_id, vec) in bits_map {
        if vec.len() != first_len {
            return Err(SpectrumProcessingError::MatrixDimensionMismatch(
                format!("**Vector for scan {} has length {} but expected {}**", scan_id, vec.len(), first_len)
            ));
        }
    }
    
    let n = first_len;
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
    Ok(pruned)
}

/// Set columns to zero if they have a hit in at least `min_frac * n_spec` spectra.
pub fn prune_background_columns(
    bits_map: &mut HashMap<String, Vec<f32>>,
    min_frac: f64
) -> Result<()> {
    if min_frac < 0.0 || min_frac > 1.0 {
        return Err(SpectrumProcessingError::InvalidParameters(
            "**min_frac must be between 0.0 and 1.0**".to_string()
        ));
    }
    
    let n_spec = bits_map.len();
    if n_spec == 0 { 
        return Ok(());
    }

    // assume all vectors have the same length
    let len = bits_map.values().next().unwrap().len();
    let threshold = (n_spec as f64 * min_frac).ceil() as usize;

    // 1) Count spectra with a "hit" in each column
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

    // 3) Zero out background columns
    for vec in bits_map.values_mut() {
        for &j in &bg_cols {
            vec[j] = 0.0;
        }
    }
    
    Ok(())
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
) -> Result<()> {
    if min_frac < 0.0 || min_frac > 1.0 {
        return Err(SpectrumProcessingError::InvalidParameters(
            "**min_frac must be between 0.0 and 1.0**".to_string()
        ));
    }
    
    let n_spec = sparse_map.len();
    if n_spec == 0 {
        return Ok(());
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
    
    // collect all bins that are "too common"
    let bg_bins: HashSet<usize> = counts
        .into_iter()
        .filter_map(|(bin, cnt)| if cnt >= threshold { Some(bin) } else { None })
        .collect();
    
    // now prune each sparse vector in place
    for vec in sparse_map.values_mut() {
        vec.retain(|&(bin, _)| !bg_bins.contains(&bin));
    }
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::{compute_pairwise_similarity_matrix_ndarray,
        compute_pairwise_similarity_matrix_sparse,
        SparseVec};
    use crate::ms_io::import_mzml;
    use hashbrown::HashMap;

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
    fn sparse_matches_dense() -> Result<()> {
        println!("\n=== TESTING THAT SPARSE VECTORS MATCHES DENSE VECTORS===");
        // 1) Build a tiny "normalized" dense map
        let mut bits_map: HashMap<String, Vec<f32>> = HashMap::new();
        // orthonormal basis vectors
        bits_map.insert("1".into(), vec![1.0, 0.0]);
        bits_map.insert("2".into(), vec![0.0, 1.0]);
        // 45° vector
        let norm = (2.0f32).sqrt();
        bits_map.insert("3".into(), vec![1.0 / norm, 1.0 / norm]);

        // 2) Dense similarity
        let (dense_scans, dense_mat) = compute_pairwise_similarity_matrix_ndarray(&bits_map)?;
        println!("||    Dense scans: {dense_scans:?}");
        println!("||    Dense matrix:\n{dense_mat}");

        // 3) Convert to sparse & compute
        let sparse_map = dense_to_sparse(&bits_map);
        let dummy_metadata: HashMap<String, SpectrumMetadata> = HashMap::new();
        let (sparse_scans, sparse_mat) =
            compute_pairwise_similarity_matrix_sparse(&sparse_map, &dummy_metadata, "cosine".to_string(), 1, 1.0)?;
        println!("||    Sparse scans: {sparse_scans:?}");
        println!("||    Sparse matrix:\n{sparse_mat}");

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
                    "**Mismatch at ({},{}) → dense={}, sparse={}**",
                    i, j, d, s
                );
            }
        }
        
        Ok(())
    }

    #[test]
    fn test_bitvec_and_cosine() -> Result<()> {
        // Load sample data
        println!("||    Importing test mzML file...");
        let (map, _) = import_mzml("tests\\data\\FeatureFinderCentroided_1_input.mzML")
                    .map_err(|e| SpectrumProcessingError::ParseError(format!("import_mzml failed: {:?}", e)))?;
        println!("||    Number of Spectra in test file: {}", map.len());
        
        // Determine the max m/z
        let max_mz = map.values().map(|p| p.iter().map(|p| p.mz).fold(0., f32::max)).fold(0., f32::max);
        println!("||    Max m/z: {}", max_mz);
        
        // Determine bin width
        let bin_width = 0.01;
        println!("||    Bin width: {}", bin_width);
        
        // Generate map of bin Vectors
        let bits_map = compute_dense_vec_map(&map, bin_width, max_mz, 0.0)?;
            
        // Compute cosine similarity matrix
        let sims = compute_pairwise_similarity_matrix_ndarray(&bits_map)?;
            
        // Print confirmation if sims is not empty
        if sims.0.len() > 0 {
            println!("||    Cosine Similarity Matrix Successfully Generated!");
        }
        Ok(())
    }

    #[test]
    fn test_cosine_vs_modified_cosine() -> Result<()> {
        use crate::ms_io::SpectrumMetadata;

        println!("\n=== TESTING COSINE VS MODIFIED COSINE ===");
        
        // Create test sparse vectors that represent MS/MS spectra
        let mut sparse_map: HashMap<String, SparseVec> = HashMap::new();
        let mut spec_metadata: HashMap<String, SpectrumMetadata> = HashMap::new();
        
        // Mass tolerance for binning (this should match your bin_width)
        let mass_tolerance = 1.0; // 1 Da bins for this test
        
        // Spectrum A: precursor m/z = 500.0
        let spectrum_a = vec![
            (100, 0.8), // normalized intensity
            (150, 0.4),
            (200, 0.3),
            (250, 0.3)
        ];
        
        // Spectrum B: precursor m/z = 510.0 (10 Da higher than A)
        let spectrum_b = vec![
            (110, 0.8), // same relative intensities, shifted by 10 bins
            (160, 0.4),
            (210, 0.3),
            (260, 0.3)
        ];
        
        // Spectrum C: precursor m/z = 500.0 (same as A)
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
        )?;
        
        let (scans_modified, matrix_modified) = compute_pairwise_similarity_matrix_sparse(
            &sparse_map,
            &spec_metadata,
            "modified-cosine".to_string(),
            2,
            mass_tolerance,
        )?;
        
        // The scan order should be the same
        assert_eq!(scans_cosine, scans_modified);
        
        // Find indices for our spectra (they should be sorted as "1", "2", "3")
        let idx_1 = scans_cosine.iter().position(|x| x == "1").unwrap();
        let idx_2 = scans_cosine.iter().position(|x| x == "2").unwrap();
        let idx_3 = scans_cosine.iter().position(|x| x == "3").unwrap();
        
        println!("=== COSINE SIMILARITY RESULTS ===");
        println!("||    Spectrum 1 vs 2 (cosine): {:.4}", matrix_cosine[(idx_1, idx_2)]);
        println!("||    Spectrum 1 vs 3 (cosine): {:.4}", matrix_cosine[(idx_1, idx_3)]);
        println!("||    Spectrum 2 vs 3 (cosine): {:.4}", matrix_cosine[(idx_2, idx_3)]);
        
        println!("\n=== MODIFIED COSINE SIMILARITY RESULTS ===");
        println!("||    Spectrum 1 vs 2 (modified): {:.4}", matrix_modified[(idx_1, idx_2)]);
        println!("||    Spectrum 1 vs 3 (modified): {:.4}", matrix_modified[(idx_1, idx_3)]);
        println!("||    Spectrum 2 vs 3 (modified): {:.4}", matrix_modified[(idx_2, idx_3)]);
        
        // Test expectations:
        let cosine_1_vs_2 = matrix_cosine[(idx_1, idx_2)];
        assert!(cosine_1_vs_2 < 0.1, "**Regular cosine similarity between shifted spectra should be low, got {}**", cosine_1_vs_2);
        
        let modified_1_vs_2 = matrix_modified[(idx_1, idx_2)];
        assert!(modified_1_vs_2 > 0.8, "**Modified cosine similarity between shifted but similar spectra should be high, got {}**", modified_1_vs_2);
        
        let improvement = modified_1_vs_2 - cosine_1_vs_2;
        assert!(improvement > 0.7, "**Modified cosine should show significant improvement over regular cosine for shifted spectra, improvement: {}**", improvement);
        
        let cosine_1_vs_3 = matrix_cosine[(idx_1, idx_3)];
        let modified_1_vs_3 = matrix_modified[(idx_1, idx_3)];
        let diff_1_vs_3 = (modified_1_vs_3 - cosine_1_vs_3).abs();
        assert!(diff_1_vs_3 < 0.1, "**For spectra with same precursor mass, cosine and modified cosine should be similar, diff: {}**", diff_1_vs_3);
        
        println!("\n=== TEST SUMMARY ===");
        println!("||    Regular cosine failed to match shifted spectra (similarity: {:.4})", cosine_1_vs_2);
        println!("||    Modified cosine successfully matched shifted spectra (similarity: {:.4})", modified_1_vs_2);
        println!("||    Improvement: {:.4}", improvement);
        println!("||    Both methods agree on truly different spectra");
        
        Ok(())
    }

    #[test]
    fn test_modified_cosine_edge_cases() -> Result<()> {
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
        )?;
        
        let (_, matrix_modified) = compute_pairwise_similarity_matrix_sparse(
            &sparse_map,
            &spec_metadata,
            "modified-cosine".to_string(),
            2,
            mass_tolerance,
        )?;
        
        // For small mass differences, both should give identical results
        let cosine_sim = matrix_cosine[(0, 1)];
        let modified_sim = matrix_modified[(0, 1)];
        
        println!("||    Small mass diff - Cosine: {:.6}, Modified: {:.6}", cosine_sim, modified_sim);
        
        // They should be nearly identical since the shift is less than 1 bin
        assert!((cosine_sim - modified_sim).abs() < 0.001, 
            "|| For small mass differences, cosine and modified cosine should be nearly identical");
            
        Ok(())
    }

    #[test]
    fn test_vector_shifting_with_real_data() -> Result<()> {
        use crate::ms_io::import_mzml;
        use hashbrown::HashMap;
        
        println!("\n=== TESTING VECTOR SHIFTING WITH REAL DATA ===");
        
        // Load sample data
        println!("||    Importing test mzML file...");
        let (spec_map, metadata_map) = import_mzml("tests\\data\\FeatureFinderCentroided_1_input.mzML")
            .map_err(|e| SpectrumProcessingError::ParseError(format!("import_mzml failed: {:?}", e)))?;
        println!("||    Number of Spectra in test file: {}", spec_map.len());
        
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
        
        println!("||    Number of MS2 spectra: {}", ms2_specs.len());
        
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
        
        println!("||    Max m/z: {:.2}", max_mz);
        println!("||    Bin width: {}", bin_width);
        
        // Generate sparse vectors
        let sparse_map = compute_sparse_vec_map(&ms2_specs, bin_width, max_mz, minimum_intensity)?;
        
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
        println!("\n||  Testing with spectra: {} vs {}", scan_a, scan_b);
        println!("||  Mass difference: {:.2} Da", mass_diff);
        
        let meta_a = &ms2_metadata[scan_a];
        let meta_b = &ms2_metadata[scan_b];
        
        println!("||    Spectrum A - target_mz: {:.2}, peaks: {}", meta_a.target_mz, sparse_map[scan_a].len());
        println!("||    Spectrum B - target_mz: {:.2}, peaks: {}", meta_b.target_mz, sparse_map[scan_b].len());
        
        // Calculate expected shift
        let delta_mz = meta_b.target_mz - meta_a.target_mz;
        let shift_bins = (delta_mz / bin_width).round() as isize;
        
        println!("\n||  hift calculation:");
        println!("||    Delta m/z: {:.2}", delta_mz);
        println!("||    Expected shift in bins: {}", shift_bins);
        
        // Get original vectors
        let vec_a = &sparse_map[scan_a];
        let vec_b = &sparse_map[scan_b];
        
        // Perform shifts (same logic as in modified cosine)
        let vec_b_shifted = shift_sparse_vec(vec_b, -shift_bins)?;
        let vec_a_shifted = shift_sparse_vec(vec_a, shift_bins)?;
        
        println!("\n=== VECTOR COMPARISON ===");
        println!("||    Original A: {} peaks", vec_a.len());
        println!("||    Original B: {} peaks", vec_b.len());
        println!("||    B shifted by {} bins: {} peaks", -shift_bins, vec_b_shifted.len());
        println!("||    A shifted by {} bins: {} peaks", shift_bins, vec_a_shifted.len());
        
        println!("||    Vector shifting with real data test completed successfully!");
        Ok(())
    }

    // Helper function to convert sparse vector to dense for visualization
    // fn sparse_to_dense_helper(sparse: &SparseVec, size: usize) -> Vec<f32> {
    //     let mut dense = vec![0.0; size];
    //     for &(bin, intensity) in sparse {
    //         if bin < size {
    //             dense[bin] = intensity;
    //         }
    //     }
    //     dense
    // }

    #[test]
    fn test_error_handling() {
        // Test invalid parameters
        assert!(spectrum_to_dense_vec(&[], -1.0, 100.0, 0.0).is_err());
        assert!(spectrum_to_dense_vec(&[], 1.0, -100.0, 0.0).is_err());
        assert!(spectrum_to_dense_vec(&[], 1.0, 100.0, -1.0).is_err());
        
        // Test empty input handling
        let empty_map: HashMap<String, Vec<f32>> = HashMap::new();
        assert!(compute_pairwise_similarity_matrix_ndarray(&empty_map).is_err());
        
        // Test invalid similarity metric
        let sparse_map: HashMap<String, SparseVec> = HashMap::new();
        let metadata: HashMap<String, SpectrumMetadata> = HashMap::new();
        assert!(compute_pairwise_similarity_matrix_sparse(
            &sparse_map, 
            &metadata, 
            "invalid_metric".to_string(), 
            2, 
            1.0
        ).is_err());
        
        // Test modified cosine with MS1 data
        let mut sparse_map: HashMap<String, SparseVec> = HashMap::new();
        sparse_map.insert("1".to_string(), vec![(100, 1.0)]);
        assert!(compute_pairwise_similarity_matrix_sparse(
            &sparse_map, 
            &metadata, 
            "modified-cosine".to_string(), 
            1, // MS1 should fail
            1.0
        ).is_err());
        
        println!("||    Error handling tests passed!");
    }
}