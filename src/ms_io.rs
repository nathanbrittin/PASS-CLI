use std::fs;
use std::collections::HashMap;
// use std::error::Error;
use std::io::{Cursor, Read};

use roxmltree::Document;
use base64::engine::general_purpose::STANDARD;
use base64::Engine;
use byteorder::{BigEndian, ReadBytesExt};
use flate2::read::ZlibDecoder;

/// A single m/z-intensity pair.
#[derive(Debug, PartialEq, Clone)]
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
        println!("------------------------------------------------------------------------------");
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
                    println!("Found m/z data for scan {}: {}: {:?}", scan, mzs.len(), mzs);
                } else if is_int {
                    ints = vals;
                    println!("Found intensity data for scan {}: {}: {:?}", scan, ints.len(), ints);
                }
            }
        }

        // Check lengths match or record error
        if mzs.len() != ints.len() {
            errors.push(format!(
                "Spectrum {} has mismatched lengths: m/z={} intensity={}",
                scan,
                mzs.len(),
                ints.len()
            ));
        } else {
            // zip into Peak structs and insert into map
            let mut peaks: Vec<Peak> = mzs
                .into_iter()
                .zip(ints.into_iter())
                .map(|(mz, intensity)| Peak { mz, intensity })
                .collect();

            // Sort peaks
            peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());

            result.insert(scan, peaks);
        }
    }

    // Sort the HashMap by scan ID
    result = result.into_iter().map(|(k, v)| (k.clone(), v.clone())).collect();

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
    fn test_import_on_samples() {
        // Iterate through all .mzML files and check that each spectrum has the correct number of spectra
        let mut file_num_spectra_map = HashMap::new();
        file_num_spectra_map.insert(r"tests\data\1min.mzML".to_string(), 39);
        file_num_spectra_map.insert(r"tests\data\2min.mzML".to_string(), 38);
        file_num_spectra_map.insert(r"tests\data\tiny.msdata.mzML0.99.9.mzML".to_string(), 2);
        file_num_spectra_map.insert(r"tests\data\tiny.msdata.mzML0.99.10.mzML".to_string(), 2);
        file_num_spectra_map.insert(r"tests\data\tiny.pwiz.1.1.mzML".to_string(), 4);
        file_num_spectra_map.insert(r"tests\data\tiny.pwiz.mzML0.99.9.mzML".to_string(), 2);
        file_num_spectra_map.insert(r"tests\data\tiny.pwiz.mzML0.99.10.mzML".to_string(), 2);
        file_num_spectra_map.insert(r"tests\data\tiny1.mzML0.99.0.mzML".to_string(), 2);
        file_num_spectra_map.insert(r"tests\data\tiny1.mzML0.99.1.mzML".to_string(), 2);
        // file_num_spectra_map.insert(r"tests\data\tiny2_SRM.mzML0.99.0.mzML".to_string(), 2);
        // file_num_spectra_map.insert(r"tests\data\tiny2_SRM.mzML0.99.1.mzML".to_string(), 2);
        file_num_spectra_map.insert(r"tests\data\tiny4_LTQ-FT.mzML0.99.0.mzML".to_string(), 2);
        file_num_spectra_map.insert(r"tests\data\tiny4_LTQ-FT.mzML0.99.1.mzML".to_string(), 2);
        for (file, expected_count) in file_num_spectra_map {
            let data = import_mzml(&file);
            assert!(data.is_ok(), "Failed to import {}", file);
            let map = data.unwrap();
            assert_eq!(map.len(), expected_count, "File {} has incorrect number of spectra", file);
        }
        println!("All tests passed!");
    }
}
