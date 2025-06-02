use anyhow::{Context, Result};
use mzdata::io::prelude::*;
use mzdata::spectrum::MultiLayerSpectrum;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Peak {
    pub mz: f64,
    pub intensity: f64,
}

#[derive(Debug, Clone)]
pub struct Spectrum {
    pub retention_time: f64,
    pub precursor_mz: Option<f64>,
    pub peaks: Vec<Peak>,
    pub ms_level: u8,
}

impl Spectrum {
    pub fn max_intensity(&self) -> f64 {
        self.peaks
            .iter()
            .map(|p| p.intensity)
            .fold(0.0, f64::max)
    }

    pub fn normalize(&mut self) {
        let max_intensity = self.max_intensity();
        if max_intensity > 0.0 {
            for peak in &mut self.peaks {
                peak.intensity /= max_intensity;
            }
        }
    }

    pub fn filter_peaks(&mut self, min_intensity: f64) {
        self.peaks.retain(|p| p.intensity >= min_intensity);
    }

    pub fn sort_peaks_by_mz(&mut self) {
        self.peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());
    }
}

#[derive(Debug)]
pub struct MSRun {
    pub spectra: Vec<Spectrum>,
    pub base_peak_chromatogram: Vec<(f64, f64)>, // (retention_time, intensity)
}

impl MSRun {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let extension = path
            .extension()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_lowercase();

        match extension.as_str() {
            "mzml" => Self::from_mzml(path),
            "mzxml" => Self::from_mzxml(path),
            _ => Err(anyhow::anyhow!(
                "Unsupported file format. Expected .mzML or .mzXML"
            )),
        }
    }

    fn from_mzml<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut reader = MZReader::open_path(path)
            .context("Failed to open mzML file")?;

        let mut spectra = Vec::new();
        let mut base_peak_chromatogram = Vec::new();

        for spectrum_result in reader.iter() {
            let spectrum = spectrum_result
                .context("Failed to read spectrum")?;

            let retention_time = spectrum
                .start_time()
                .unwrap_or(0.0);

            let ms_level = spectrum.ms_level().unwrap_or(1);

            // Extract peaks
            let mut peaks = Vec::new();
            let arrays = spectrum.arrays.as_ref().unwrap();
            let mz_array = &arrays.mzs;
            let intensity_array = &arrays.intensities;

            for (mz, intensity) in mz_array.iter().zip(intensity_array.iter()) {
                if *intensity > 0.0 {
                    peaks.push(Peak {
                        mz: *mz,
                        intensity: *intensity,
                    });
                }
            }

            if peaks.is_empty() {
                continue;
            }

            // Get precursor m/z for MS2+ spectra
            let precursor_mz = if ms_level > 1 {
                spectrum.precursor_ion()
                    .and_then(|p| p.mz())
            } else {
                None
            };

            // Calculate base peak intensity for chromatogram
            let base_peak_intensity = peaks
                .iter()
                .map(|p| p.intensity)
                .fold(0.0, f64::max);

            base_peak_chromatogram.push((retention_time, base_peak_intensity));

            let mut spec = Spectrum {
                retention_time,
                precursor_mz,
                peaks,
                ms_level,
            };

            spec.sort_peaks_by_mz();
            spectra.push(spec);
        }

        // Sort spectra by retention time
        spectra.sort_by(|a, b| a.retention_time.partial_cmp(&b.retention_time).unwrap());
        
        // Sort chromatogram by retention time
        base_peak_chromatogram.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        Ok(MSRun {
            spectra,
            base_peak_chromatogram,
        })
    }

    fn from_mzxml<P: AsRef<Path>>(path: P) -> Result<Self> {
        // mzdata crate handles mzXML through the same interface
        Self::from_mzml(path)
    }

    pub fn get_ms1_spectra(&self) -> Vec<&Spectrum> {
        self.spectra
            .iter()
            .filter(|s| s.ms_level == 1)
            .collect()
    }

    pub fn get_ms2_spectra(&self) -> Vec<&Spectrum> {
        self.spectra
            .iter()
            .filter(|s| s.ms_level == 2)
            .collect()
    }

    pub fn get_spectra_in_rt_range(&self, start_rt: f64, end_rt: f64) -> Vec<&Spectrum> {
        self.spectra
            .iter()
            .filter(|s| s.retention_time >= start_rt && s.retention_time <= end_rt)
            .collect()
    }
}