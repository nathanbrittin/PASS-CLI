use std::fs::File;
use std::fs;
use std::collections::HashMap;
use std::error::Error;
use std::path::Path;
use std::io::{Cursor, Read};
use core::clone::Clone;
use serde_json;
use roxmltree::Document;
use base64::engine::general_purpose;
use base64::Engine;
use byteorder::{LittleEndian, ReadBytesExt};
use flate2::read::ZlibDecoder;
use csv::WriterBuilder;
use ndarray::{Array2, Axis};

/// A single m/z-intensity pair.
#[derive(Debug, PartialEq, Clone)]
pub struct Peak {
    pub mz: f64,
    pub intensity: f64,
}

#[derive(Clone)]
pub struct SpectrumMetadata {
    pub index: String,
    pub id: String,
    pub ms_level: u8,
    pub base_peak_mz: f64,
    pub base_peak_intensity: f64,
    pub total_ion_current: f64,
    pub scan_window_lower_limit: f64,
    pub scan_window_upper_limit: f64,
    pub target_mz: f64,
    pub selected_ion: f64,
    pub charge: u8
}

/// Parses an mzML XML file, returning a map from scan ID to its peaks.
/// Returns Err if any spectrum has mismatched m/z vs intensity lengths.
pub fn import_mzml(file_path: &str) -> Result<(HashMap<String, Vec<Peak>>, HashMap<String, SpectrumMetadata>), Vec<String>> {
    // Read the file
    let xml = fs::read_to_string(file_path).map_err(|e| vec![e.to_string()])?;
    // Parse
    let doc = Document::parse(&xml).map_err(|e| vec![e.to_string()])?;

    // Iterate over each spectrum
    let mut result: HashMap<String, Vec<Peak>> = HashMap::new();
    let mut result_metadata: HashMap<String, SpectrumMetadata> = HashMap::new();
    let mut errors: Vec<String> = Vec::new();

    for spec in doc.descendants().filter(|n| n.tag_name().name() == "spectrum") {
        let scan = spec
            .attribute("index")
            .or_else(|| spec.attribute("id"))
            .unwrap_or("unknown")
            .to_string();

        let id = spec.attribute("id").unwrap_or("unknown").to_string();
        let index = spec.attribute("index").unwrap_or("unknown").to_string();

        let mut base_peak_mz = 0.0;
        let mut base_peak_intensity = 0.0;
        let mut total_ion_current = 0.0;
        let mut scan_window_lower_limit = 0.0;
        let mut scan_window_upper_limit = 0.0;
        let mut target_mz = 0.0;
        let mut selected_ion = 0.0;
        let mut charge = 0;

        // Get metadata from cvParams
        let mut ms_level = 1;
        for cv in spec.descendants().filter(|n| n.tag_name().name() == "cvParam") {
            // MS:1000511 = ms level, MS:1000504 = base peak m/z, MS:1000505 = base peak intensity, MS:1000285 = total ion current
            // MS:1000501 = scan window lower limit, MS:1000500 = scan window upper limit
            // MS:1000827 = target m/z, MS:1000744 = selected ion, MS:1000041 = charge
            match cv.attribute("accession") {
                Some("MS:1000511") => ms_level = cv.attribute("value").unwrap().parse().unwrap(),
                Some("MS:1000504") => base_peak_mz = cv.attribute("value").unwrap().parse().unwrap(),
                Some("MS:1000505") => base_peak_intensity = cv.attribute("value").unwrap().parse().unwrap(),
                Some("MS:1000285") => total_ion_current = cv.attribute("value").unwrap().parse().unwrap(),
                Some("MS:1000501") => scan_window_lower_limit = cv.attribute("value").unwrap().parse().unwrap(),
                Some("MS:1000500") => scan_window_upper_limit = cv.attribute("value").unwrap().parse().unwrap(),
                Some("MS:1000827") => target_mz = cv.attribute("value").unwrap().parse().unwrap(),
                Some("MS:1000744") => selected_ion = cv.attribute("value").unwrap().parse().unwrap(),
                Some("MS:1000041") => charge = cv.attribute("value").unwrap().parse().unwrap(),
                _ => (),
            }
        }
        
        let spec_metadata = SpectrumMetadata {
            index,
            id,
            ms_level,
            base_peak_mz,
            base_peak_intensity,
            total_ion_current,
            scan_window_lower_limit,
            scan_window_upper_limit,
            target_mz,
            selected_ion,
            charge
        };

        result_metadata.insert(scan.clone(), spec_metadata);

        let mut mzs = Vec::new();
        let mut ints = Vec::new();

        // Iterate over each binary data array
        for bda in spec.descendants().filter(|n| n.tag_name().name() == "binaryDataArray") {
            let mut is_mz = false;
            let mut is_int = false;
            let mut compressed = false;
            let mut is_64bit   = false;
            
            // Iterate over each cvParam to check for the type
            for cv in bda.descendants().filter(|n| n.tag_name().name() == "cvParam") {
                match cv.attribute("accession") {
                    Some("MS:1000514") => is_mz = true,
                    Some("MS:1000515") => is_int = true,
                    Some("MS:1000521") => (),    // 32-bit
                    Some("MS:1000523") => is_64bit = true,  // 64-bit
                    Some("MS:1000574") => compressed = true,
                    Some("MS:1000576") => (),    // no compression
                    _ => (),
                }
            }

            // Extract the data
            if let Some(blob) = bda
                .descendants()
                .find(|n| n.tag_name().name() == "binary")
                .and_then(|bin| bin.text())
            {
                if is_mz {
                    mzs = decode_mz_array(blob, is_64bit, compressed);
                } else if is_int {
                    ints = decode_int_array(blob, is_64bit, compressed);
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

    // Return result or errors
    if errors.is_empty() {
        Ok((result, result_metadata))
    } else {
        Err(errors)
    }
}

fn decode_mz_array(base64_seq: &str, is_64bit: bool, is_zlib: bool) -> Vec<f64> {
    // 1) Base64 decode
    let decoded = general_purpose::STANDARD
        .decode(base64_seq)
        .expect("base64 decode failed");

    // 2) If compressed, zlib-decompress; if that fails, fall back to the raw bytes
    let data: Vec<u8> = if is_zlib {
        let mut decompressor = ZlibDecoder::new(&decoded[..]);
        let mut buf = Vec::new();
        match decompressor.read_to_end(&mut buf) {
            Ok(_) => buf,
            Err(_) => decoded,
        }
    } else {
        decoded
    };

    // 3) Read as 32 or 64-bit little-endian floats
    let mut cursor = Cursor::new(&data);
    let mut floats = Vec::new();
    if is_64bit {
        while let Ok(val) = cursor.read_f64::<LittleEndian>() {
            floats.push(val);
        }
    } else {
        while let Ok(val) = cursor.read_f32::<LittleEndian>() {
            floats.push(val as f64);
        }
    }
    floats
}

fn decode_int_array(base64_seq: &str, is_64bit: bool, is_zlib: bool) -> Vec<f64> {
    // Base64 decode
    let decoded = general_purpose::STANDARD
        .decode(base64_seq)
        .expect("base64 decode failed");

   // If compressed, zlib-decompress; if that fails, fall back to the raw bytes
    let data: Vec<u8> = if is_zlib {
        let mut decompressor = ZlibDecoder::new(&decoded[..]);
        let mut buf = Vec::new();
        match decompressor.read_to_end(&mut buf) {
            Ok(_) => buf,
            Err(_) => decoded,
        }
    } else {
        decoded
    };

    // Read as 32 or 64-bit little-endian integers
    let mut cursor = Cursor::new(&data);
    let mut floats = Vec::new();
    if is_64bit {
        while let Ok(val) = cursor.read_f64::<LittleEndian>() {
            floats.push(val);
        }
    } else {
        while let Ok(val) = cursor.read_f32::<LittleEndian>() {
            floats.push(val as f64);
        }
    }
    floats
}

/// Supported output formats for the similarity matrix.
#[derive(Clone, Copy)]
pub enum OutputFormat {
    Csv,
    Tsv,
    Json,
}

/// Write the similarity matrix to a file in CSV, TSV, or JSON.
///
/// # Arguments
///
/// * `scans` - Slice of scan IDs in row/column order.
/// * `mat` - The n×n similarity matrix.
/// * `output_path` - File path to write to.
/// * `format` - Desired output format.
pub fn write_similarity_matrix<P: AsRef<Path>>(
    scans: &[String],
    mat: &Array2<f32>,
    output_path: P,
    format: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    match format {
        OutputFormat::Csv | OutputFormat::Tsv => {
           let delimiter = if let OutputFormat::Tsv = format { b'\t' } else { b',' };
            let mut wtr = WriterBuilder::new()
                .delimiter(delimiter)
                .from_path(output_path)?;

            // Write header: empty corner + scan IDs
            let mut header = Vec::with_capacity(scans.len() + 1);
            header.push(String::new());
            header.extend_from_slice(scans);
            wtr.write_record(&header)?;

            // Write each row: scan ID + row values
            for (i, scan_id) in scans.iter().enumerate() {
                let mut record = Vec::with_capacity(scans.len() + 1);
                record.push(scan_id.clone());
                // row as f32 strings
                for val in mat.index_axis(Axis(1), i) {
                    // Make the value only 4 digits after decimal point
                    let val = format!("{:.4}", val);
                    record.push(val.to_string());
                }
                wtr.write_record(&record)?;
            }

            wtr.flush()?;
        }
        OutputFormat::Json => {
            // Flatten matrix rows as owned Vec<f32>
            let rows: Vec<Vec<f32>> = mat
                .axis_iter(Axis(0))
                .map(|row| row.to_vec())
                .collect();
            let container = serde_json::json!({
                "scans": scans,
                "matrix": rows
            });
            let file = File::create(output_path)?;
            serde_json::to_writer_pretty(file, &container)?;
        }
    }
    Ok(())
}


    /// Filter a map of scan IDs to their corresponding peak data by MS level.
    ///
    /// # Arguments
    ///
    /// * `map` - The map from scan IDs to their peak data.
    /// * `map_metadata` - The map from scan IDs to their corresponding metadata.
    /// * `ms_level` - The desired MS level to filter by.
    ///
    /// # Returns
    ///
    /// A new `HashMap` where each key-value pair is a filtered version of the
    /// input map. Only scans with the specified MS level are included in the
    /// output.
pub fn filter_by_ms_level(map: HashMap<String, Vec<Peak>>, map_metadata: HashMap<String, SpectrumMetadata>, ms_level: u8) -> HashMap<String, Vec<Peak>> {
    let mut filtered_map = HashMap::new();
    for (scan_id, peaks) in map {
        if let Some(metadata) = map_metadata.get(&scan_id) {
            if metadata.ms_level == ms_level {
                filtered_map.insert(scan_id, peaks);
            }
        }
    }
    filtered_map
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::path::Path;

    #[test]
    fn test_import_on_samples() {
        // Map of pyOpenMS sample file names to expected spectrum counts
        let mut file_num_spectra_map: HashMap<&str, usize> = HashMap::new();
        let data_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests").join("data");
        println!("Data directory: {}", data_dir.display());
        file_num_spectra_map.insert("190509_Ova_native_25ngul_R.mzML", 194);
        file_num_spectra_map.insert("BSA1.mzML", 1684);
        file_num_spectra_map.insert("FeatureFinderCentroided_1_input.mzML", 112);
        file_num_spectra_map.insert("FeatureFinderMetaboIdent_1_input.mzML", 380);
        file_num_spectra_map.insert("Metabolomics_1_aligned.mzML", 1052);
        file_num_spectra_map.insert("Metabolomics_1.mzML", 1052);
        file_num_spectra_map.insert("Metabolomics_2_aligned.mzML", 963);
        file_num_spectra_map.insert("Metabolomics_2.mzML", 963);
        file_num_spectra_map.insert("peakpicker_tutorial_1_baseline_filtered.mzML", 1);
        file_num_spectra_map.insert("PeakPickerHiRes_input.mzML", 5);
        file_num_spectra_map.insert("PrecursorPurity_input.mzML", 7);
        file_num_spectra_map.insert("small.mzML", 236);

        for (fname, &expected_count) in &file_num_spectra_map {
            let path = data_dir.join(fname);
            let path_str = path.to_str().expect("Invalid path");
            println!("----------------------------------------------------------------------------------");
            println!("Testing file: {}", path_str);

            // Verify that the file exists and is readable
            assert!(Path::new(path_str).exists(), "File {} does not exist", fname);

            // Use the public API import_mzml_to_map
            let result = import_mzml(path_str);
            assert!(result.is_ok(), "Failed to import {}", fname);
            let (map, _) = result.unwrap();

            // Ensure the number of spectra matches the expected count
            let num_spectra = map.len();
            assert_eq!(num_spectra, expected_count, "Incorrect number of spectra in {}", fname);

            println!("File {} has {} spectra", fname, map.len());

            // Print a summary of the first spectra (scan number, mzs, ints)
            let (scan, peak) = map.iter().next().unwrap();
            // Grab the first mz and int
            let mz = peak[0].mz;
            let int = peak[0].intensity;
            println!(
                "First spectrum: scan {}, ex. m/z: {}, int: {}",
                scan, mz, int
            );
        }
    }

    #[test]
    fn test_decoding() {
        let mut binary_map = HashMap::new();
        // Add binary sequences to map
        binary_map.insert("example_1", "AAAAYP4ZhEAAAAAgPyKEQAAAAKAIKoRAAAAAQDkuhEAAAADATzKEQAAAAOAINoRAAAAAwNw5hEAAAAAAKEKEQAAAAKBWSoRAAAAAYBdahEAAAABgKWKEQAAAACAAaoRAAAAAYBZyhEAAAACADXaEQAAAAMAKeoRAAAAAQEKChEAAAAAAQoaEQAAAAAAdioRAAAAAoN+NhEAAAABgAJKEQA==");
        binary_map.insert("example_2", "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");
        binary_map.insert("example_3", "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");
        binary_map.insert("example_4", "AAAAwNwVaUA=");
        binary_map.insert("example_5", "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");

        // Iterate and decode
        let mut decoded_map = HashMap::new();
        for (name, b64) in binary_map {
            let mzs = decode_mz_array(b64, true, true);
            decoded_map.insert(name, mzs.clone());
            // println!("{} → {} points, {:?}", name, mzs.clone().len(), &mzs);
        }
        // Confirm with the known values
        let mut confirmation_map = HashMap::new();
            
        confirmation_map.insert("example_1", vec![643.2492065429688, 644.2808227539063, 645.2542114257813, 645.7779541015625, 646.2889404296875, 646.7543334960938, 647.2327880859375, 648.26953125, 649.2922973632813, 651.2614135742188, 652.2702026367188, 653.2500610351563,654.2609252929688, 654.756591796875, 655.2552490234375, 656.2823486328125,656.7822265625, 657.26416015625, 657.7341918945313, 658.2501831054688,]);
        confirmation_map.insert("example_2", vec![1781.9692735887363, 1782.0406065523073, 1782.11194379921, 1782.183285329787, 1805.08903473707, 1805.161760452934, 1805.234490564, 1805.307225070622, 1805.3799639731549, 1805.452707271952, 1805.5254549673682, 1805.5982070597577, 1805.670963549475, 1805.743724436874, 1805.8164904696157, 1805.8892601534424, 1805.9620342360142, 1806.0348127176858, 1820.6794132501532, 1820.7530831976412, 1820.8267576165833, 1820.9004365073406, 1820.974119870276, 1821.0478077057508, 1821.121500014127, 1821.1951967957666, 1821.268898051032, 1821.3426037802853, 1821.4163141003544, 1821.4900287786695, 1821.5637479320587, 1821.6374715608843, 1823.1128845372227, 1823.1867022339559, 1823.2605244141039, 1823.33435107803, 1823.4081822260976, 1823.4820178586692, 1823.5558579761087, 1823.6297025787787, 1823.7035516670428, 1823.777405710447, 1823.8512637709891, 1823.9251263182152, 1823.998993352489, 1832.8957210594679, 1832.9701337200231, 1833.0445509122399, 1833.118972636486, 1833.1933988931296, 1833.2678296825386, 1833.3422650050813, 1833.4167048611257, 1833.49114925104, 1833.5655981751925, 1833.6400528866673, 1833.7145108804004, 1833.7889734094765, 1833.8634404742638, 1858.9915183883586, 1859.067525885866, 1859.143538044996, 1859.219554866129, 1859.2955763496475, 1859.371602495932, 1859.4476333053638, 1859.5236687783247, 1859.5997089151956, 1859.6757562059292, 1859.7518056717652, 1859.827859802656, 1859.903918598983, 1910.8222153489364, 1910.9014237032532, 1910.9806369827652, 1911.0598551878818, 1911.1390783190104, 1911.2183063765603, 1911.2975393609395, 1911.3767772725566, 1911.4560201118202, 1911.5352678791385, 1911.614520635173, 1911.6937782598275, 1911.7730408137634, 1911.8523082973888, 1912.0901403304929, 1912.1694275369664, 1912.2487196751745, 1912.3280167455255, 1912.4073187484291, 1912.4866256842943, 1912.56593755353, 1912.6452543565454, 1912.7245760937503, 1912.8039027655532, 1912.8832352418553, 1912.9625717840831, 1913.0419132621375, 1913.1212596764276, 1930.2156324651946, 1930.2960497493061, 1930.3764720590975, 1930.4568993949877, 1930.5373317573956, 1930.6177691467396, 1930.6982115634394, 1930.7786590079136, 1930.859111480581, 1930.9395537438447, 1931.0200162741562, 1931.1004838339184, 1931.18095642355, 1938.525071643414, 1938.606008780176, 1938.6869509859825, 1938.7678982612574, 1938.8488506064239, 1938.9298080219053, 1939.0107705081248, 1939.091738065506, 1939.1727106944725, 1939.253688811444, 1939.334671584852, 1939.4156594311162, 1939.4966523506603, 1948.9263717405013, 1949.0079611707718, 1949.0895557246315, 1949.171155402509, 1949.2527602048335, 1949.3343701320343, 1949.4159851845402, 1949.4976053627806, 1949.579230667185, 1949.660861098182, 1949.742497185222, 1949.824137870693, 1949.905783684045, 1949.987434625708, 1962.2109763263074, 1962.2934014040559, 1962.3758316754818, 1962.4582671410217, 1962.5407078011121, 1962.6231536561893, 1962.7056047066897, 1962.7880609530498, 1962.8705223957063, 1962.952989074917, 1963.0354609114765, 1963.117937945642, 1963.2004201778511, 1963.28290760854, 1963.3654002381456, 1963.4478980671054, 1963.5304010958557, 1963.6129093248337, 1963.6954227544766, 1963.7779413852213, 1963.8604652175047, 1963.9429942910622, 1964.025528527735, 1964.108067967259, 1964.190612610071, 1964.2731624566081, 1964.3557175073083, 1964.4382777626083, 1964.5208432229463, 1964.6034138887594, 1964.6859897604857, 1964.7685708385623, 1964.851157123427, 1964.933748654924, 1965.0163453546786, 1965.0989472625345, 1965.1815543789298, 1965.2641667043022, 1965.34678423909, 1965.4294069837315, 1965.5120349386639, 1965.5946681043258, 1965.6773064811555, 1965.7599500695908, 1965.8425988700699, 1965.9252529606078, 1966.0079121864903, 1966.0905766257315, 1966.1732462787702, 1967.0002296645025, 1967.082956698262, 1967.165688951084, 1967.2484264234067, 1967.3311691156698, 1967.413917028312, 1967.4966701617727, 1967.5794285164914, 1967.6621920929067, 1967.7449608914585, 1967.8277349125858, 1967.9105141962364, 1967.993298663833, 1968.0760883553235, 1968.1588832711475, 1968.2416834117444, 1968.324488777554, 1968.407299369016, 1968.49011518657, 1968.5729362306556, 1968.6557625017128, 1968.7385940001816, 1968.8214307265016, 1968.90427373147, 1968.9871209148127, 1969.0699733273266, 1969.152830969452, 1994.3327088295175, 1994.4171661662447, 1994.501628868089, 1994.5860969355049, 1994.670570368947, 1994.7550491688696, 1994.8395333357275, 1994.924022869975, 1995.008517772067, 1995.0930180424582, 1995.1775253551536, 1995.2620363635074, 1995.346552741525, 1995.4310744896611, 2036.2814090318054, 2036.3685450785615, 2036.455686718514, 2036.542833952142, 2036.6299867799244, 2036.7171452023397, 2036.8043092198668, 2036.8914788329848, 2036.9786540421724, 2037.0658362404106, 2037.1530226431748, 2037.240214643446, 2037.3274122417033, 2071.14031502956, 2071.229698168522, 2071.319087093825, 2071.4084818059673, 2071.4978823054485, 2071.587288592769, 2071.6767006684277, 2071.766118532925, 2071.8555421867595, 2071.9449726273915, 2072.0344078614016, 2072.123848886249, 2072.213295702434, 2096.1223973304327, 2096.2134025686123, 2096.3044137335687, 2096.3954308258153, 2096.486453845868, 2096.57748279424, 2096.6685176714477, 2096.7595584780047, 2096.850605214426, 2096.941658807901, 2097.0327174055974, 2097.1237819347034, 2097.214852395735, 2119.3370393172654, 2119.429560586319, 2119.522087914159, 2119.614621301312, 2119.7071607483094, 2119.799706255679, 2119.8922578239512, 2119.984815453654, 2120.077379145317, 2120.169949092516, 2120.262524909688, 2120.35510679041, 2120.4476947352096, 2124.0633602149474, 2124.156191157774, 2124.2490281864407, 2124.341871301479, 2124.434720503421, 2124.5275757927984, 2124.620437170144, 2124.71330463599, 2124.806178190868, 2124.899057835312, 2124.9919435698516, 2125.084835436247, 2125.1777333525793, 2125.270637360606, 2125.3635474608604, 2125.4564636538744, 2125.549385940181, 2125.6423143203137, 2125.735248794804, 2125.8281893641865, 2125.9211360289924, 2126.014088789756, 2126.107047684438, 2126.2000126387156, 2126.29298369055, 2126.3859608404746, 2126.478944089022, 2126.571933436727, 2126.6649288841218, 2126.7579304317396, 2126.8509380801156, 2126.968945674774, 2127.0619655262653, 2127.1549914801144, 2127.2480235368553, 2129.01679287479, 2129.1099471015573, 2129.20310744244, 2129.2962738979722, 2129.38944646869, 2129.4826251551285, 2129.5758099578225, 2129.669000877308, 2129.7621979141195, 2129.855401068793, 2129.948610341864, 2130.0418257338674, 2130.138651119114, 2130.2318787505906, 2130.3251125026068, 2130.4183523756983, 2130.5115983704013, 2130.604850487252, 2130.6981087267864, 2130.7913730895393, 2130.884643576048, 2130.977920186848, 2131.0843019987024, 2131.177590859694, 2131.2708858465867, 2131.3641869599155, 2131.4574942002173, 2131.5508075680286, 2131.644127063886, 2131.737452688326, 2131.8307844418855, 2131.924122325101, 2132.0174663385096, 2132.069119925127, 2132.1624762005317, 2132.2558386077403, 2132.3492071472883, 2132.4425818197146, 2132.535962625556, 2132.629349565349, 2132.7227426396307, 2132.8161418489394, 2132.9095482723415, 2133.002959753315, 2133.0963773709273, 2133.189801125716, 2159.116137536818, 2159.211275910654, 2159.3064205728588, 2159.4015715239866, 2159.4967287645914, 2159.591892295228, 2159.6870621164503, 2159.782238228813, 2159.8774206328712, 2159.972609481523, 2160.067804470634, 2160.1630057531042, 2160.258213329488, 2162.831199019935, 2162.9265830513773, 2163.0219733928484, 2163.117370044904, 2163.212773008101, 2163.308182282996, 2163.4035978701454, 2163.4990197701068, 2163.594447983437, 2163.6898825106923, 2163.7853235152447, 2163.880770672022, 2163.976224144397, 2164.0716839329248, 2166.938417366632, 2167.0340732326435, 2167.129735432672, 2167.2254039672766, 2167.3210788370166, 2167.416760042452, 2167.5124475841417, 2167.608141462645, 2167.703841678522, 2167.7995483982677, 2167.89526129057, 2167.9909805219254, 2168.086706092892, 2171.057346297642, 2171.1532750321417, 2171.2492101247617, 2171.345151576064, 2171.4410993866095, 2171.537053556962, 2171.6330140876817, 2171.728980979332, 2171.8249542324743, 2171.9209338476717, 2172.016920284199, 2172.112912625193, 2172.208911329929, 2172.3049163989695, 2182.6143501895035, 2182.7110459278847, 2182.807748092259, 2182.9044566831963, 2183.0011717012667, 2183.097893147039, 2183.1946210210826, 2183.2913553239678, 2183.3880960562637, 2183.484843991902, 2183.58159758473, 2183.678357608678, 2183.7751240643165, 2201.984142062365, 2202.0821278735093, 2202.1801202252127, 2202.2781191180575, 2202.3761245526257, 2202.4741365295004, 2202.5721550492626, 2202.6701801124955, 2202.768211719781, 2202.8662401989645, 2202.964284896104, 2203.062336139044, 2203.1603939283673, 2668.312432062636, 2668.4431390949876, 2668.5738557316286, 2668.704581973501, 2668.8353178215443, 2668.9660632767013, 2669.0968183399127, 2669.22758301212, 2669.3583572942644, 2669.4891416878745, 2669.619935192719, 2669.750738310326, 2669.8815510416384, 2685.3852166997826, 2685.517180221396, 2685.6491534705997, 2685.781136448349, 2685.9131291556023, 2686.0451315933133, 2686.177143762441, 2686.30916566394, 2686.441197298768, 2686.573238667882, 2686.7052897949466, 2686.8373506355024, 2686.969421213215, 2687.1015015290413, 2687.2335915839394, 2687.365691378867, 2687.4978009470783, 2687.6299202249365, 2687.7620492456967, 2687.894188010317, 2688.026336519756, 2688.1584947749716, 2688.2906627769216, 2688.422840526566, 2688.5550280248613, 2688.6872252727676, 2688.8194322712434, 2688.9516490212477, 2689.0838755237396, 2689.2161117796777, 3636.459320951778, 3636.6672742744454, 3636.8752454356636, 3637.083234437472]);
        confirmation_map.insert("example_3", vec![300.06577092599525, 300.08932376441896, 300.2018173360483, 300.28960754478305, 300.2981327863004, 300.3346816032124, 301.1415329911221, 301.21654231553987, 301.9879869591199, 302.0453317036539, 302.14468923704726, 302.16034243700244, 303.08239388782465, 303.1121322675646, 303.2321950715736, 303.28969877168703, 303.98475623793445, 304.0242798872948, 304.0607795757647, 304.08140951060466, 304.09729659736263, 304.1755344266212, 304.24871526835136, 304.28498516972763, 305.12736148798683, 305.171087617278, 305.24793123134737, 305.2516968114553, 305.268882796776, 305.96745383771906, 306.04031753309573, 306.0765388865728, 306.08786505738226, 306.17431214358896, 306.1918933381472, 306.26416201074625, 306.3005546190294, 307.07998314707777, 307.0841429771381, 307.1329099048425, 307.1902745021192, 307.2273219916964, 307.26405430648833, 308.0162544397057, 308.0561261656696, 308.09197484651986, 308.27995890254005, 309.2273140178521, 309.27936780888496, 309.2827527515315, 309.7022124175072, 310.0342192839173, 310.2007830928396, 310.203990196222, 311.1449553557436, 311.1506720431911, 311.16035999678786, 311.29507493900087, 311.3316653113528, 312.0296582973308, 312.04025304963454, 312.06572597276363, 312.25362217838915, 313.0329055338333, 313.1439146819056, 313.2744789310042, 313.31059592635074, 314.04510660697883, 314.1470893565481, 314.3134894860743, 315.04826981713313, 315.2325175663815, 315.2530807334047, 316.06098578038313, 316.1761233278128, 316.19698890464065, 317.11530021832164, 317.2004194559206, 317.2476934896519, 317.30506128208935, 318.0764063688618, 318.087463876401, 318.0964379641782, 318.200481800701, 319.2043485382809, 319.26319435766965, 320.0920049759792, 320.1026601793058, 321.09635366340393, 321.15232648172525, 321.2060623281531, 321.2437072705161, 321.84244333920657, 322.0183516578768, 322.1076333254574, 322.274243909091, 324.0202615225021, 324.1119586674887, 324.15565944258606, 325.3106090265444, 325.34719185408034, 326.045108725874, 326.08121980390115, 326.1269636107289, 326.17097260120545, 326.26912317267943, 327.0084510790779, 327.02391701004717, 327.07845355967436, 327.1114752548257, 327.3259782444802, 328.06078791623133, 328.06675129943244, 328.08156017589044, 329.00542480277466, 329.0263721624682, 330.05095011023195, 330.0764616358369, 330.15944224621904, 330.1697177069163, 330.1758692701768, 330.1914254480724, 330.3369239907588, 331.0025044629897, 332.0553941207956, 332.09156454544666, 332.1200783409621, 332.19168907939434, 332.25883387209785, 332.2952392130187, 332.5379967850014, 332.8729002011879, 332.93646986430826, 333.2424106077791, 333.30065612339604, 334.0458427555852, 334.0714286237368, 334.10742826236583, 334.3319265983347, 335.0470080567035, 335.22182321052367, 335.2797884534358, 336.0466331702311, 336.08655794677696, 336.1292292093568, 336.2255256974554, 337.12587839151337, 337.25882298186855, 339.36238132351457, 340.10614445804407, 341.094166297876, 341.26896194888053, 341.3416885983188, 342.07631651013537, 342.1211175408017, 343.1544219317576, 343.2844950600877, 344.0912941324163, 344.1050801425514, 344.2281701911388, 344.2879362789383, 345.1083305003673, 345.15272376240796, 345.2024857352375, 346.1078242775908, 346.2060109770754, 346.29534180547716, 347.0364888440506, 347.20090887348186, 347.2174395871094, 347.29793456394015, 348.12322767311275, 349.23783128284134, 350.1268442957349, 350.32633783311803, 352.00909059734516, 352.14282814779176, 352.23815402007483, 353.268918151187, 353.3778559148341, 354.2852251149236, 354.7056368006603, 355.0702776750471, 355.2076338256516, 355.2884448903345, 356.06983769639373, 356.2798953864984, 356.93610197107, 357.04961389865434, 357.06639496164536, 357.28641721942665, 357.3001666423217, 358.0665668927001, 358.07018994287205, 358.3034633216913, 360.15026028776396, 360.22312615489903, 361.10832222788997, 361.2332353899544, 363.04469090169675, 363.17226146539406, 363.25314897652044, 363.293488236833, 364.9255296061541, 366.12158825744064, 367.22292326464674, 368.1377912661735, 369.12542612272046, 369.1825811779155, 370.1286058239166, 370.1532620052534, 371.3157188257477, 372.1008639841235, 372.31915796909584, 373.0807020771242, 373.0982676706859, 373.3232062358784, 374.06634243562496, 374.0976116620109, 374.94718566267255, 375.09511915687176, 376.1544137973484, 376.2174610827652, 377.1275406837615, 377.2689993760731, 379.2877619125361, 381.0554345507501, 381.4090507197432, 381.6866405867304, 383.14505458660346, 385.1599325212386, 386.119111117379, 386.15182398578787, 386.2211940839725, 386.2898686817156, 388.127988134391, 388.3425264630013, 389.112163300556, 389.3456648234618, 390.0611045096754, 390.10720027804854, 390.1335213395728, 390.17008630370225, 391.25061135630233, 391.2557723265717, 391.27552927758086, 391.2840663410693, 391.3150333874384, 392.103495061887, 392.2435280478371, 392.287520454517, 392.9579723021781, 393.0992000725768, 393.2909602946618, 393.9943296755873, 396.04983116786855, 396.31446217841676, 396.92750823245007, 397.2945734252651, 398.29774898614596, 399.250522645929, 399.34667253893963, 400.25354614689365, 404.15943650222243, 404.18602043403246, 405.054940919055, 405.2998967964835, 406.05715441156383, 407.0569017584098, 407.20371394030985, 408.2075042682958, 411.171946148813, 413.26609402688484, 413.3625218846327, 414.2696157941232, 418.2014983899513, 419.3155065263834, 420.3189979137288, 421.0115255463586, 421.15719662618096, 421.3223139791545, 422.826521905188, 423.1593474760457, 423.18976280204765, 424.23927321509973, 425.81481842876644, 426.14916069395673, 426.4797600369745, 426.76253190394584, 427.3783096028602, 428.381327441505, 429.24068501700293, 430.27997215630626, 430.8884295694197, 432.8851290133681, 433.331285407009, 434.3345999702576, 435.0204134069146, 435.2356111654327, 438.74412862360646, 441.3213409262396, 444.4046751574796, 445.2300130261354, 445.7320569463724, 446.1421756736548, 446.23236856600806, 446.24892203629963, 447.34658786343755, 448.11559608142807, 448.12230273317437, 448.34986182342107, 449.09517002483847, 449.11370880415643, 452.11245640379235, 454.12802958607324, 455.764998840321, 457.7235689938244, 458.22540232179335, 458.3107413042411, 459.1718511458954, 461.36246443349955, 462.14636331462066, 463.1308036322965, 463.14578481718934, 463.2098753204567, 464.1251123825406, 464.1433722808555, 465.12578962078896, 465.141791033475, 466.13993562776295, 466.8521002683997, 467.1016219589245, 469.0891826504268, 476.19826863027095, 478.17824009803434, 479.34773689469483, 482.9534548047927, 483.5389779757922, 496.0084248607909, 499.1920214170151, 499.51958317622575, 499.85400892690245, 500.1879422917905, 503.1075234084895, 504.1061342865272, 504.2279889783473, 505.08651922010984, 505.1051332636861, 505.11636709610036, 505.22460275080823, 505.61850347582924, 506.08529533244996, 507.083057229636, 507.22275999147257, 509.84904488013103, 510.18333648958554, 510.5176609620802, 511.1838526053851, 514.1893910400889, 514.5235041857638, 514.8553115303058, 518.8549620675617, 518.9901790758638, 519.1879984834034, 519.5219142491613, 519.8542316195463, 520.189934566082, 521.120007280286, 522.1129383375669, 522.1346587340684, 523.1119093175676, 523.1318437403029, 525.8956342135886, 529.889660653381, 530.1782348813133, 530.5124066457134, 530.8459554062916, 531.1770837297595, 531.5116303102745, 531.8451328156815, 533.191213571393, 536.1652069802592, 537.165690554377, 538.147084414892, 538.1623801905953, 538.1725654571201, 539.1614924959924, 539.182314961809, 539.1997143218365, 539.5337683223456, 539.867622103337, 540.1580973919104, 540.167816030319, 540.1996714752104, 541.1204547651669, 541.1574866558914, 541.911328425844, 542.4049688043835, 543.9058347911914, 544.9065176004398, 545.9042329026179, 547.7647561976321, 547.9003886425804, 548.2673403223193, 548.766127243682, 548.900997380026, 550.2171423671628, 550.5067609675161, 561.9166621009176, 562.9181872051932, 564.2379233958118, 565.053162433659, 565.2415693075473, 565.9113589690614, 570.4363882361984, 577.1259253345877, 579.1055011042679, 580.1047280690273, 581.1021503403105, 582.0798640492886, 591.2228912398559, 591.7248366944465, 592.2189163489313, 593.1577753351431, 594.1578819133756, 595.1543891585054, 596.1322113538263, 596.154231735106, 597.132134079257, 597.1501644490644, 598.1288943548636, 610.1842554452198, 611.1846140627105, 612.1809998090316, 612.1927850234154, 613.1803506057414, 614.1768154933488, 618.7209184834993, 639.0722219701454, 653.1247074715898, 654.1240223733496, 655.1209769667864, 656.0995735366292, 667.1768176471533, 668.1766571875644, 669.1729029060152, 669.1856651466671, 670.1510680665256, 670.1732488101055, 671.150685718646, 671.1704786198087, 672.1479316045459, 684.2033126065634, 685.2031009213292, 686.1995668033246, 686.4631604382547, 687.199385874575, 688.1964449210457, 718.7597078518974, 727.1438185642227, 741.1959406159153, 742.1959084381465, 743.1924600500697, 744.1701032202724, 744.1922455445751, 745.1880201392116, 747.7633460111631, 748.2647403962501, 748.7665570522424, 757.7994719456809, 758.2221851231714, 759.2218739253365, 760.2188716561647, 761.2194731735834, 777.7801693789845, 790.0187210547588, 794.7648690699343, 795.2668737297647]);
        confirmation_map.insert("example_4", vec![200.68319702148438]);
        confirmation_map.insert("example_5", vec![200.00523376464844, 200.0917510986328, 200.12815856933594, 200.13958740234375, 200.18386840820313, 200.20091247558594, 200.97235107421875, 201.00831604003906, 201.10983276367188, 201.11219787597656, 201.1234130859375, 201.14866638183594, 201.20440673828125, 201.90542602539063, 201.93319702148438, 201.97247314453125, 201.9757537841797, 202.1800079345703, 202.21661376953125, 202.97109985351563, 202.98727416992188, 203.1390380859375, 203.9384765625, 204.13829040527344, 204.92803955078125, 204.96681213378906, 205.045166015625, 205.08599853515625, 205.0896759033203, 205.09693908691406, 205.12245178222656, 205.15476989746094, 205.1585693359375, 205.1855926513672, 205.1950225830078, 205.97036743164063, 205.9740447998047, 206.08116149902344, 206.15402221679688, 206.9879913330078, 207.1017303466797, 207.11282348632813, 207.15899658203125, 207.1743927001953, 207.96990966796875, 207.9852294921875, 208.09681701660156, 208.1282196044922, 208.1310577392578, 208.13319396972656, 208.17762756347656, 208.9437255859375, 208.9550018310547, 208.98484802246094, 208.98878479003906, 209.0810089111328, 209.11720275878906, 209.12828063964844, 209.13662719726563, 209.15365600585938, 209.9820098876953, 210.02394104003906, 210.11245727539063, 210.12789916992188, 210.94818115234375, 211.060302734375, 211.10780334472656, 211.1329803466797, 211.1441650390625, 211.16928100585938, 211.670166015625, 211.9274444580078, 211.9403839111328, 212.09156799316406, 212.10299682617188, 212.1432647705078, 212.96420288085938, 213.12315368652344, 213.1486358642578, 213.18505859375, 213.9226531982422, 214.21665954589844, 215.08139038085938, 215.10264587402344, 215.1280059814453, 215.13890075683594, 215.16419982910156, 215.96412658691406, 216.01148986816406, 216.10182189941406, 216.12307739257813, 216.96334838867188, 217.06588745117188, 217.077880859375, 217.10702514648438, 217.1546630859375, 217.1954803466797, 217.95950317382813, 218.08116149902344, 218.1103973388672, 218.13868713378906, 218.95956420898438, 218.98263549804688, 219.0445556640625, 219.04745483398438, 219.0765838623047, 219.1379852294922, 219.1743927001953, 219.21060180664063, 219.95635986328125, 220.13331604003906, 220.17767333984375, 221.016845703125, 221.080810546875, 221.11721801757813, 221.1536865234375, 221.18980407714844, 221.6698455810547, 222.11251831054688, 222.1488494873047, 222.52590942382813, 222.96095275878906, 223.06361389160156, 223.09390258789063, 223.11585998535156, 223.15225219726563, 223.16941833496094, 223.20542907714844, 223.9195556640625, 223.98838806152344, 224.0918426513672, 224.12814331054688, 224.92723083496094, 224.95889282226563, 225.09083557128906, 225.10980224609375, 225.13145446777344, 225.13449096679688, 225.14859008789063, 225.1849822998047, 225.19615173339844, 225.995849609375, 226.0975341796875, 226.13206481933594, 226.13499450683594, 226.18003845214844, 226.90992736816406, 226.95428466796875, 227.0732879638672, 227.12509155273438, 227.1389617919922, 227.1643524169922, 227.1754608154297, 227.20042419433594, 227.9542236328125, 227.95762634277344, 227.9827423095703, 228.0866241455078, 228.12315368652344, 228.23220825195313, 229.14341735839844, 229.4440460205078, 231.1128692626953, 231.93325805664063, 232.90487670898438, 232.9354705810547, 232.9801025390625, 233.09205627441406, 233.1284637451172, 233.15377807617188, 233.96177673339844, 234.07623291015625, 234.13185119628906, 234.1488800048828, 234.18540954589844, 234.94326782226563, 234.94650268554688, 234.95167541503906, 234.96234130859375, 235.1081085205078, 235.13308715820313, 235.16932678222656, 235.20553588867188, 235.90455627441406, 235.9512176513672, 235.96995544433594, 235.98019409179688, 236.1129608154297, 236.1281280517578, 236.1727752685547, 236.9386444091797, 237.12335205078125, 237.18504333496094, 237.22137451171875, 237.90834045410156, 237.9383087158203, 238.0647735595703, 238.10757446289063, 238.14376831054688, 238.22518920898438, 238.9354248046875, 239.08895874023438, 239.14889526367188, 239.16433715820313, 239.20054626464844, 239.2368621826172, 239.90721130371094, 240.12294006347656, 240.97157287597656, 241.0682373046875, 241.15499877929688, 241.17979431152344, 241.21624755859375, 242.10218811035156, 243.13368225097656, 243.17047119140625, 243.17453002929688, 243.21788024902344, 243.92218017578125, 243.96644592285156, 244.173828125, 244.92047119140625, 244.93594360351563, 244.965087890625, 246.13392639160156, 246.91763305664063, 246.9354705810547, 246.94403076171875, 247.1075897216797, 247.16954040527344, 247.94329833984375, 248.90049743652344, 249.12339782714844, 249.18495178222656, 249.9439239501953, 250.12689208984375, 250.1439208984375, 250.90130615234375, 250.95574951171875, 250.9712677001953, 251.0001678466797, 251.04685974121094, 251.0703125, 251.16429138183594, 251.18515014648438, 251.20065307617188, 251.91436767578125, 252.00311279296875, 252.159423828125, 252.94378662109375, 252.95399475097656, 252.96267700195313, 252.9972686767578, 253.11831665039063, 253.16287231445313, 253.17971801757813, 253.21624755859375, 254.1833038330078, 254.94927978515625, 254.99432373046875, 255.11036682128906, 255.13381958007813, 255.19552612304688, 255.23184204101563, 255.9490509033203, 256.2635192871094, 256.8817443847656, 256.9084777832031, 256.9537658691406, 257.2272033691406, 257.2479248046875, 257.2669372558594, 257.96771240234375, 258.9671936035156, 259.16943359375, 259.1905517578125, 259.89453125, 259.9051818847656, 259.9644775390625, 260.1939392089844, 261.13079833984375, 262.12860107421875, 262.14410400390625, 262.946533203125, 263.16448974609375, 263.2368469238281, 263.8995361328125, 264.1594543457031, 264.872802734375, 265.1183776855469, 265.14361572265625, 265.18023681640625, 265.25250244140625, 265.8281555175781, 265.9031066894531, 265.9305114746094, 266.788818359375, 266.86761474609375, 266.8733825683594, 266.9110107421875, 267.17449951171875, 267.1954650878906, 267.505615234375, 267.9520263671875, 268.17767333984375, 268.9516296386719, 268.9659118652344, 269.2112731933594, 269.2475891113281, 269.8804626464844, 269.9250183105469, 271.1329345703125, 271.1695861816406, 271.2267761230469, 271.8468933105469, 272.2301940917969, 272.8810119628906, 272.9598693847656, 273.1485595703125, 273.8936462402344, 273.9019775390625, 274.1514587402344, 274.9013671875, 274.909912109375, 274.91876220703125, 274.9442138671875, 275.20074462890625, 275.888916015625, 275.94354248046875, 275.9783630371094, 276.896728515625, 277.1799011230469, 277.8969421386719, 277.92633056640625, 278.8970642089844, 279.0933532714844, 279.1015930175781, 279.1592102050781, 279.2319641113281, 279.9093933105469, 280.09661865234375, 280.1051025390625, 280.13665771484375, 280.23541259765625, 280.86798095703125, 280.89215087890625, 280.9836730957031, 281.1400146484375, 281.1749572753906, 281.2474365234375, 281.41839599609375, 281.8756103515625, 281.9101867675781, 282.02301025390625, 282.1700134277344, 282.2216491699219, 282.2505798339844, 282.88336181640625, 283.17529296875, 283.2253723144531, 283.2632141113281, 284.2947998046875, 284.8340148925781, 284.89398193359375, 284.90289306640625, 284.9284973144531, 285.2790222167969, 285.2984924316406, 286.8307189941406, 286.8543395996094, 286.93231201171875, 287.88897705078125, 288.9005126953125, 288.94647216796875, 289.6796875, 290.91253662109375, 290.9542236328125, 291.40277099609375, 292.8670959472656, 292.9407958984375, 293.1091613769531, 293.2471618652344, 293.8549499511719, 294.8629455566406, 295.2269287109375, 295.863037109375, 295.87066650390625, 295.8966369628906, 295.90252685546875, 295.9366760253906, 296.13153076171875, 296.8599853515625, 296.91790771484375, 297.01922607421875, 297.08251953125, 297.1349182128906, 297.2424011230469, 297.8753662109375, 298.0225830078125, 298.2459411621094, 298.8774108886719, 299.0163269042969, 299.0616455078125, 299.2012634277344, 299.258056640625, 299.8890380859375, 299.9081115722656, 300.0199279785156, 300.2893371582031, 300.9384765625, 301.013427734375, 301.179931640625, 301.888427734375, 301.8965759277344, 302.9046936035156, 302.939208984375, 303.1954040527344, 303.8838195800781, 303.8918762207031, 303.92108154296875, 304.8918762207031, 304.90472412109375, 304.9317321777344, 304.98046875, 305.2476501464844, 305.9683837890625, 306.2511901855469, 308.8690490722656, 308.8794250488281, 308.9001770019531, 309.0435791015625, 309.2034912109375, 309.8762512207031, 311.29339599609375, 311.9450988769531, 312.12664794921875, 313.2738952636719, 313.8474426269531, 313.882080078125, 314.8501892089844, 315.8798828125, 315.8857727050781, 316.8768310546875, 317.1964111328125, 317.8767395019531, 319.8626708984375, 321.8948974609375, 323.8785400390625, 324.91259765625, 324.9211120605469, 325.0434265136719, 325.1436462402344, 325.9207763671875, 326.8686218261719, 326.9935302734375, 327.0082092285156, 327.0782775878906, 327.2529296875, 327.9969177246094, 328.0116882324219, 328.3212890625, 328.9906311035156, 329.00531005859375, 329.267578125, 329.8915710449219, 329.9214782714844, 329.9936828613281, 331.0023498535156, 331.28460693359375, 331.8850402832031, 331.89251708984375, 331.9161376953125, 332.8868103027344, 332.9103698730469, 332.9991149902344, 334.87664794921875, 335.9451904296875, 337.8932189941406, 341.30535888671875, 341.84259033203125, 343.8456726074219, 343.88043212890625, 344.8714904785156, 344.9076843261719, 345.84295654296875, 345.879638671875, 346.83306884765625, 346.89898681640625, 349.9142761230469, 352.3396911621094, 352.9161376953125, 352.9216613769531, 353.343017578125, 353.8750305175781, 353.95574951171875, 354.33673095703125, 354.8741760253906, 354.9322509765625, 355.0701904296875, 355.2845458984375, 356.8590087890625, 356.87396240234375, 357.0677490234375, 359.3158874511719, 359.8534851074219, 359.8930358886719, 360.870361328125, 360.91131591796875, 361.86395263671875, 362.8712463378906, 363.78985595703125, 364.9100036621094, 365.90118408203125, 369.3004150390625, 369.35223388671875, 371.10137939453125, 371.31591796875, 371.7725524902344, 371.8408203125, 372.3192138671875, 372.881103515625, 373.18389892578125, 373.3219299316406, 374.8504333496094, 374.90386962890625, 375.8482666015625, 378.37353515625, 379.87445068359375, 379.94793701171875, 380.88238525390625, 381.8006896972656, 381.86993408203125, 382.89239501953125, 383.3313903808594, 383.811279296875, 385.3468933105469, 389.8586120605469, 391.2847900390625, 391.8271179199219, 392.28765869140625, 392.8355407714844, 392.8544921875, 392.861572265625, 392.9057312011719, 393.29779052734375, 393.8634033203125, 394.3013610839844, 394.346923828125, 395.91864013671875, 397.91094970703125, 398.9122009277344, 399.3473205566406, 400.3507995605469, 400.9036865234375, 401.822021484375, 401.85638427734375, 402.8907165527344, 402.89874267578125, 403.8983154296875, 405.01031494140625, 408.7881774902344, 408.85479736328125, 412.85430908203125, 413.8075866699219, 413.8333740234375, 414.8410949707031, 414.8921203613281, 416.3737487792969, 419.3160095214844, 424.8317565917969, 424.8725891113281, 426.7982177734375, 427.3783874511719, 428.38189697265625, 429.089599609375, 429.8114013671875, 430.8095703125, 430.8646240234375, 430.88775634765625, 431.85223388671875, 432.8852233886719, 433.20635986328125, 434.74676513671875, 434.8361511230469, 434.8548278808594, 434.8624572753906, 435.2961120605469, 435.822509765625, 439.8412170410156, 441.3945007324219, 443.8013610839844, 445.1207275390625, 445.9025573730469, 446.12005615234375, 446.9018859863281, 447.1170959472656, 447.34710693359375, 447.89886474609375, 448.8992614746094, 450.8581237792969, 452.3221130371094, 452.8730163574219, 453.88128662109375, 454.8696594238281, 455.8603210449219, 462.14739990234375, 462.8592529296875, 464.81072998046875, 465.7989501953125, 467.8119201660156, 467.83990478515625, 468.79461669921875, 469.8352355957031, 470.7239074707031, 471.8377380371094, 472.77825927734375, 472.84051513671875, 473.82037353515625, 478.7914733886719, 478.877685546875, 480.78692626953125, 481.7822265625, 482.7300720214844, 482.8171081542969, 482.8523864746094, 484.7293395996094, 485.81304931640625, 490.7427062988281, 493.7796936035156, 493.7958068847656, 494.8507385253906, 495.84893798828125, 499.8193664550781, 500.7660217285156, 501.815673828125, 503.1086730957031, 503.92724609375, 507.832763671875, 508.7808837890625, 513.8093872070313, 513.9700317382813, 514.8727416992188, 519.139404296875, 520.139404296875, 521.1190185546875, 521.1361083984375, 521.7451171875, 522.7479858398438, 523.7822265625, 524.758544921875, 524.8025512695313, 526.756591796875, 532.74560546875, 534.7874145507813, 536.1658935546875, 536.7882080078125, 537.1550903320313, 537.166015625, 538.1632690429688, 541.8035278320313, 542.7215576171875, 548.7364501953125, 552.7523193359375, 553.4191284179688, 559.8115844726563, 560.7835693359375, 566.3176879882813, 569.73046875, 585.7570190429688, 588.797607421875, 593.1588745117188, 593.6732788085938, 594.15771484375, 601.7276611328125, 607.7799072265625, 610.1847534179688, 611.183837890625, 616.8135986328125, 622.7173461914063, 623.6702270507813, 634.7681274414063, 639.3335571289063, 648.707763671875, 657.654296875, 663.5307006835938, 678.6911010742188, 683.6181030273438, 684.2041625976563, 693.0982055664063, 704.960205078125, 715.6765747070313, 727.599609375, 728.6154174804688, 747.0206909179688, 758.2227172851563, 759.2227783203125, 759.5880126953125, 796.6764526367188]);
        
        // Compare the vectors between the two scans
        for (name, expected) in &confirmation_map {
            let got = &decoded_map[name];
            // 1) same length
            assert_eq!(
                got.len(), expected.len(),
                "{}: decoded length {} != expected length {}",
                name, got.len(), expected.len()
            );

            // 2) same values (within a tiny epsilon)
            for (i, (a, b)) in got.iter().zip(expected.iter()).enumerate() {
                let diff = (a - b).abs();
                assert!(
                    diff < 1e-8,
                    "{}[{}] = {} but expected {}; diff {} > eps",
                    name, i, a, b, diff
                );
            }
        }
        println!("Decoding Unit Test Passed Successfully!");
    }
}

