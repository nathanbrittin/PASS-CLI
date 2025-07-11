use std::fs;
use std::collections::HashMap;
use std::error::Error;
use std::io::{Cursor, Read};

use roxmltree::Document;
use base64::engine::general_purpose::STANDARD;
use base64::Engine;
use byteorder::{BigEndian, ReadBytesExt};
use flate2::read::ZlibDecoder;

/// A single m/z-intensity pair.
#[derive(Debug, PartialEq)]
pub struct Peak {
    pub mz: f32,
    pub intensity: f32,
}

/// All data for one scan.
#[derive(Debug, PartialEq)]
pub struct SpectrumData {
    pub scan_id: String,
    pub peaks: Vec<Peak>,
}

/// Parses an mzML XML file, returning a map from scan ID to its peaks.
/// Returns Err if any spectrum has mismatched m/z vs intensity lengths.
pub fn import_mzml(
    file_path: &str,
) -> Result<HashMap<String, Vec<Peak>>, Vec<String>> {
    // Read the file
    let xml = fs::read_to_string(file_path).map_err(|e| vec![e.to_string()])?;
    // Parse
    let doc = Document::parse(&xml).map_err(|e| vec![e.to_string()])?;

    // Iterate over each spectrum
    let mut result: HashMap<String, Vec<Peak>> = HashMap::new();
    let mut errors: Vec<String> = Vec::new();

    for spec in doc.descendants().filter(|n| n.tag_name().name() == "spectrum") {
        let scan = spec
            .attribute("scanNumber")
            .or_else(|| spec.attribute("id"))
            .unwrap_or("unknown")
            .to_string();

        let mut mzs = Vec::new();
        let mut ints = Vec::new();
        
        // Iterate over each binary data array
        for bda in spec.children().filter(|n| n.tag_name().name() == "binaryDataArray") {
            let mut is_mz = false;
            let mut is_int = false;
            let mut compressed = false;
            
            // Iterate over each cvParam to check for the type
            for cv in bda.children().filter(|n| n.tag_name().name() == "cvParam") {
                match cv.attribute("accession") {
                    Some("MS:1000514") => is_mz = true,
                    Some("MS:1000515") => is_int = true,
                    Some("MS:1000574") => compressed = true,
                    Some("MS:1000576") => (),
                    _ => (),
                }
            }

            // Extract the data
            if let Some(blob) = bda
                .children()
                .find(|n| n.tag_name().name() == "binary")
                .and_then(|bin| bin.text())
            {
                // Decode
                let raw = STANDARD.decode(blob.trim()).unwrap_or_default();
                let data = if compressed {
                    let mut decoder = ZlibDecoder::new(Cursor::new(raw));
                    let mut buf = Vec::new();
                    decoder.read_to_end(&mut buf).unwrap_or(0);
                    buf
                } else {
                    raw
                };

                let mut rdr = Cursor::new(data);
                let mut vals = Vec::new();
                while let Ok(v) = rdr.read_f32::<BigEndian>() {
                    vals.push(v);
                }

                // Check if we have m/z or intensity data
                if is_mz {
                    mzs = vals;
                } else if is_int {
                    ints = vals;
                }
            }
        }

        // Check lengths match or record error
        if mzs.len() != ints.len() {
            errors.push(format!(
                "Spectrum {} has mismatched lengths: m/z={} intensity={} bytes",
                scan,
                mzs.len(),
                ints.len()
            ));
        } else {
            // zip into Peak structs and insert into map
            let peaks: Vec<Peak> = mzs
                .into_iter()
                .zip(ints.into_iter())
                .map(|(mz, intensity)| Peak { mz, intensity })
                .collect();
            result.insert(scan, peaks);
        }
    }

    // Return or error
    if errors.is_empty() {
        Ok(result)
    } else {
        Err(errors)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_import_on_sample() {
        let data = import_mzml("tests/data/1min.mzML");
        assert!(data.is_ok(), "Expected successful import");
        let map = data.unwrap();
        // spectrum 39 should have 11 peaks
        let peaks = map.get("39").expect("Scan 39 missing");
        assert_eq!(peaks.len(), 70);
        // Check first peak roughly
        assert!((peaks[0].mz - 400.8245).abs() < 1e-3);
        assert!((peaks[0].intensity - 1036.0).abs() < 1e-3);
    }
}
