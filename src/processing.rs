use hashbrown::{HashMap, HashSet};
use crate::ms_io::{Peak, SpectrumMetadata};
use crate::similarity::spectrum_to_dense_vec;

/// How to normalize scans
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormMethod { #[allow(dead_code)] TIC, PQN }

/// How to threshold peak intensities before similarity computation
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize)]
pub enum ThresholdMethod {
    Absolute(f32),
    PercentBasePeak(f32),  // fraction in (0.0, 1.0]
    TopN(usize),
}

impl ThresholdMethod {
    /// Absolute floor for internal dense-vector helpers (PQN/background detection).
    /// For adaptive methods, returns 0.0 because peaks are already gated.
    pub fn dense_floor(&self) -> f32 {
        match self {
            ThresholdMethod::Absolute(v) => *v,
            _ => 0.0,
        }
    }
}

/// Config for the post-import processing
#[derive(Debug, Clone)]
pub struct ProcessingConfig {
    pub bin_width: f32,           // use your mass_tolerance
    pub max_mz: f32,              // computed from MS1
    pub threshold: ThresholdMethod,  // replaces minimum_intensity
    pub snr_min: f32,             // e.g., 5.0 (MAD-based S/N)
    pub norm: NormMethod,         // TIC or PQN
    // Optional blank filtering (bin-level background removal)
    pub blank_ratio_max: f32,         // e.g., 0.20
    pub sample_over_blank_min: f32,   // e.g., 3.0
}

/// Stats returned to print in logs
#[derive(Debug, Clone)]
pub struct ProcessingStats {
    pub n_scans: usize,
    #[allow(dead_code)]
    pub n_kept_peaks: usize,
    pub n_removed_peaks_noise: usize,
    #[allow(dead_code)]
    pub n_removed_peaks_background: usize,
    pub norm_factors_median: f32,
    pub norm_method: NormMethod,
    pub n_background_bins: usize,
}

/// Output bundle
#[derive(Debug, Clone)]
pub struct ProcessingResult {
    pub cleaned_ms1: HashMap<String, Vec<Peak>>,
    pub scale_factors: HashMap<String, f32>,        // apply to MS2 as well
    #[allow(dead_code)]
    pub background_bins: HashSet<usize>,            // which m/z bins were blank-like
    pub stats: ProcessingStats,
}

/// Main entrypoint:
/// 1) MAD/SNR noise filter per scan
/// 2) Optional background bin detection from blanks, remove those peaks
/// 3) Compute normalization factors (TIC or PQN) on MS1
/// NOTE: normalization factors are returned so you can apply to MS2 too.
pub fn preprocess_after_import(
    ms1_spec_map: &HashMap<String, Vec<Peak>>,
    metadata: &HashMap<String, SpectrumMetadata>,
    cfg: &ProcessingConfig,
    blank_scan_ids: Option<&[String]>,   // pass None if not using blanks
) -> ProcessingResult {
    // --- 0. Guard ---
    if ms1_spec_map.is_empty() {
        return ProcessingResult {
            cleaned_ms1: HashMap::new(),
            scale_factors: HashMap::new(),
            background_bins: HashSet::new(),
            stats: ProcessingStats {
                n_scans: 0, n_kept_peaks: 0, n_removed_peaks_noise: 0,
                n_removed_peaks_background: 0, norm_factors_median: 1.0,
                norm_method: cfg.norm, n_background_bins: 0,
            },
        };
    }

    // --- 1. Per-scan noise filtering using MAD-based SNR ---
    let mut cleaned: HashMap<String, Vec<Peak>> = HashMap::with_capacity(ms1_spec_map.len());
    let mut removed_noise = 0usize;
    let mut kept = 0usize;

    for (scan, peaks) in ms1_spec_map {
        let (_thr, noise) = mad_threshold(peaks);
        let sn_floor = cfg.snr_min.max(0.0) * noise;

        // Phase A: MAD/SNR floor
        let snr_passed: Vec<Peak> = peaks.iter()
            .cloned()
            .filter(|p| p.intensity.is_finite() && p.intensity >= sn_floor)
            .collect();

        // Phase B: adaptive intensity gate
        let filtered: Vec<Peak> = match cfg.threshold {
            ThresholdMethod::Absolute(min_int) => snr_passed.into_iter()
                .filter(|p| p.intensity >= min_int)
                .collect(),

            ThresholdMethod::PercentBasePeak(fraction) => {
                let base_peak = snr_passed.iter()
                    .map(|p| p.intensity)
                    .fold(0.0_f32, f32::max);
                if base_peak <= 0.0 {
                    Vec::new()
                } else {
                    let pbp_floor = base_peak * fraction;
                    snr_passed.into_iter()
                        .filter(|p| p.intensity >= pbp_floor)
                        .collect()
                }
            }

            ThresholdMethod::TopN(n) => {
                if n == 0 {
                    Vec::new()
                } else {
                    let mut ranked = snr_passed;
                    ranked.sort_unstable_by(|a, b| {
                        b.intensity.partial_cmp(&a.intensity)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    });
                    ranked.truncate(n);
                    ranked
                }
            }
        };

        removed_noise += peaks.len().saturating_sub(filtered.len());
        kept += filtered.len();
        cleaned.insert(scan.clone(), filtered);
    }

    // --- 2. Optional blank-based background bin removal (bin logic uses unnormalized dense vectors) ---
    let mut background_bins: HashSet<usize> = HashSet::new();
    let mut removed_bg_count = 0usize;
    if let Some(blank_ids) = blank_scan_ids {
        if !blank_ids.is_empty() {
            background_bins = detect_background_bins(
                &cleaned,
                blank_ids,
                cfg.bin_width,
                cfg.max_mz,
                cfg.threshold.dense_floor(),
                cfg.blank_ratio_max,
                cfg.sample_over_blank_min
            );
            // remove peaks that fall into background bins; count removals before modifying
            for (_scan, peaks) in cleaned.iter_mut() {
                peaks.retain(|p| {
                    let bin = (p.mz / cfg.bin_width).floor() as isize;
                    let keep = bin >= 0 && !background_bins.contains(&(bin as usize));
                    if !keep { removed_bg_count += 1; }
                    keep
                });
            }
            kept = cleaned.values().map(|v| v.len()).sum();
        }
    }

    // --- 3. Compute normalization factors on MS1 ---
    let scale_factors = match cfg.norm {
        NormMethod::TIC => tic_norm_factors(&cleaned, metadata),
        NormMethod::PQN => pqn_norm_factors(&cleaned, cfg.bin_width, cfg.max_mz, cfg.threshold.dense_floor()),
    };

    // Apply factors to MS1 now (we’ll apply to MS2 in main.rs using the same factors)
    apply_scale_factors_in_place(&mut cleaned, &scale_factors);

    let median_factor = median(scale_factors.values().copied().collect());
    let n_background_bins = background_bins.len();

    ProcessingResult {
        cleaned_ms1: cleaned,
        scale_factors,
        background_bins,
        stats: ProcessingStats {
            n_scans: ms1_spec_map.len(),
            n_kept_peaks: kept,
            n_removed_peaks_noise: removed_noise,
            n_removed_peaks_background: removed_bg_count,
            norm_factors_median: median_factor,
            norm_method: cfg.norm,
            n_background_bins,
        }
    }
}

// ---------- helpers ----------

fn median(mut v: Vec<f32>) -> f32 {
    if v.is_empty() { return 0.0; }
    v.sort_by(|a,b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = v.len();
    if n % 2 == 1 { v[n/2] } else { 0.5*(v[n/2-1] + v[n/2]) }
}

/// Compute median & MAD, return (threshold, MAD_as_noise).
/// Here noise is MAD*1.4826 (≈σ if Gaussian). We’ll use snr_min * noise as floor.
fn mad_threshold(peaks: &[Peak]) -> (f32, f32) {
    if peaks.is_empty() { return (0.0, 0.0); }
    let mut vals: Vec<f32> = peaks.iter().map(|p| p.intensity.max(0.0)).collect();
    let med = median(vals.clone());
    for x in vals.iter_mut() { *x = ( *x - med ).abs(); }
    let mad = median(vals) * 1.4826_f32;
    let thr = med + 3.0*mad; // not used directly for SNR floor, but useful reference
    (thr, mad.max(1e-12))
}

/// TIC-based normalization factors: factor = (median_TIC / TIC_i)
fn tic_norm_factors(
    ms1: &HashMap<String, Vec<Peak>>,
    _meta: &HashMap<String, SpectrumMetadata>
) -> HashMap<String, f32> {
    let mut tics: Vec<(String, f32)> = Vec::with_capacity(ms1.len());
    for (scan, peaks) in ms1 {
        let tic: f32 = peaks.iter().map(|p| p.intensity.max(0.0)).sum();
        tics.push((scan.clone(), tic.max(1e-12)));
    }
    let median_tic = median(tics.iter().map(|(_,tic)| *tic).collect()).max(1e-12);
    tics.into_iter().map(|(scan, tic)| (scan, median_tic / tic)).collect()
}

/// PQN factors on *unnormalized* binned vectors (no L2 normalization).
/// factor = 1 / median_i( v_i / ref_i ), ref = median spectrum across scans.
fn pqn_norm_factors(
    ms1: &HashMap<String, Vec<Peak>>,
    bin_width: f32, max_mz: f32, min_intensity: f32
) -> HashMap<String, f32> {
    // Build dense bins for each scan
    let mut dense: HashMap<String, Vec<f32>> = HashMap::with_capacity(ms1.len());
    for (scan, peaks) in ms1 {
        if let Ok(vec) = spectrum_to_dense_vec(peaks, bin_width, max_mz, min_intensity) {
            dense.insert(scan.clone(), vec);
        }
    }
    if dense.is_empty() { return HashMap::new(); }

    let n_bins = dense.values().next().map(|v| v.len()).unwrap_or(0);
    // reference spectrum = median across scans per bin.
    // Accumulate bin values in a single row-major pass (cache-friendly) rather than
    // accessing column b of every scan vector in a nested loop.
    let mut bin_vals: Vec<Vec<f32>> = vec![Vec::new(); n_bins];
    for v in dense.values() {
        for (b, &x) in v.iter().enumerate() {
            if x.is_finite() && x > 0.0 {
                bin_vals[b].push(x);
            }
        }
    }
    let mut ref_spec = vec![0.0_f32; n_bins];
    for (b, vals) in bin_vals.into_iter().enumerate() {
        ref_spec[b] = if vals.is_empty() { 0.0 } else { median(vals) };
    }

    // compute factor for each scan
    let mut factors = HashMap::with_capacity(dense.len());
    for (scan, vec) in dense {
        let mut ratios: Vec<f32> = Vec::new();
        for (x, r) in vec.iter().zip(ref_spec.iter()) {
            if *x>0.0 && *r>0.0 { ratios.push(x / r); }
        }
        let q = if ratios.is_empty() { 1.0 } else { median(ratios) };
        let factor = if q > 0.0 { 1.0 / q } else { 1.0 };
        factors.insert(scan, factor);
    }
    factors
}

/// Identify "background" bins that are abundant in blanks relative to samples,
/// using median blank/samples vectors.
fn detect_background_bins(
    cleaned_ms1: &HashMap<String, Vec<Peak>>,
    blank_ids: &[String],
    bin_width: f32, max_mz: f32, min_intensity: f32,
    blank_ratio_max: f32, sample_over_blank_min: f32
) -> HashSet<usize> {
    // Build O(1) lookup set for blank scan IDs
    let blank_set: HashSet<&String> = blank_ids.iter().collect();

    // split dense maps
    let mut blank_vecs: Vec<Vec<f32>> = Vec::new();
    let mut sample_vecs: Vec<Vec<f32>> = Vec::new();
    for (scan, peaks) in cleaned_ms1 {
        if let Ok(v) = spectrum_to_dense_vec(peaks, bin_width, max_mz, min_intensity) {
            if blank_set.contains(scan) {
                blank_vecs.push(v);
            } else {
                sample_vecs.push(v);
            }
        }
    }
    let n_bins = sample_vecs.get(0).map(|v| v.len()).or_else(|| blank_vecs.get(0).map(|v| v.len())).unwrap_or(0);
    let median_of = |vecs: &Vec<Vec<f32>>, idx: usize| -> f32 {
        let vals: Vec<f32> = vecs.iter().map(|v| v[idx]).filter(|x| x.is_finite()).collect();
        if vals.is_empty() { 0.0 } else { median(vals) }
    };

    let mut bad = HashSet::new();
    for b in 0..n_bins {
        let mb = median_of(&blank_vecs, b);
        let ms = median_of(&sample_vecs, b);
        // (i) blank is ≥ blank_ratio_max of sample, OR
        // (ii) sample/blank < sample_over_blank_min
        let cond1 = ms > 0.0 && (mb / (ms + 1e-12)) >= blank_ratio_max;
        let cond2 = mb > 0.0 && (ms / (mb + 1e-12)) < sample_over_blank_min;
        if cond1 || cond2 { bad.insert(b); }
    }
    bad
}

/// Multiply all intensities by per-scan factors (in place).
pub fn apply_scale_factors_in_place(
    spec_map: &mut HashMap<String, Vec<Peak>>,
    factors: &HashMap<String, f32>
) {
    for (scan, vec) in spec_map.iter_mut() {
        if let Some(f) = factors.get(scan) {
            for p in vec.iter_mut() {
                p.intensity = p.intensity * *f;
            }
        }
    }
}

/// Apply scale factors to another map (e.g., MS2) and return it.
pub fn apply_scale_factors(
    mut spec_map: HashMap<String, Vec<Peak>>,
    factors: &HashMap<String, f32>
) -> HashMap<String, Vec<Peak>> {
    apply_scale_factors_in_place(&mut spec_map, factors);
    spec_map
}

/// Apply SNR + adaptive threshold to a spec map without normalization (for MS2).
pub fn apply_threshold_to_spec_map(
    spec_map: HashMap<String, Vec<Peak>>,
    method: ThresholdMethod,
    snr_min: f32,
) -> HashMap<String, Vec<Peak>> {
    spec_map.into_iter().map(|(scan, peaks)| {
        let (_thr, noise) = mad_threshold(&peaks);
        let sn_floor = snr_min.max(0.0) * noise;

        let snr_passed: Vec<Peak> = peaks.into_iter()
            .filter(|p| p.intensity.is_finite() && p.intensity >= sn_floor)
            .collect();

        let filtered: Vec<Peak> = match method {
            ThresholdMethod::Absolute(min_int) => snr_passed.into_iter()
                .filter(|p| p.intensity >= min_int).collect(),
            ThresholdMethod::PercentBasePeak(fraction) => {
                let base_peak = snr_passed.iter().map(|p| p.intensity).fold(0.0_f32, f32::max);
                if base_peak <= 0.0 { Vec::new() } else {
                    let floor = base_peak * fraction;
                    snr_passed.into_iter().filter(|p| p.intensity >= floor).collect()
                }
            }
            ThresholdMethod::TopN(n) => {
                if n == 0 { return (scan, Vec::new()); }
                let mut ranked = snr_passed;
                ranked.sort_unstable_by(|a, b| b.intensity.partial_cmp(&a.intensity)
                    .unwrap_or(std::cmp::Ordering::Equal));
                ranked.truncate(n);
                ranked
            }
        };
        (scan, filtered)
    }).collect()
}
