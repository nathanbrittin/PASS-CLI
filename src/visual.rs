use plotters::prelude::*;
use ndarray::Array2;
use plotters::style::RGBColor;
use std::collections::HashMap;

/// Chromatogram data type selection
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChromatogramType {
    BPI,  // Base Peak Intensity
    TIC,  // Total Ion Current
}

impl ChromatogramType {
    pub fn from_str(input: &str) -> Option<Self> {
        match input.trim().to_lowercase().as_str() {
            "bpi" | "base peak intensity" => Some(ChromatogramType::BPI),
            "tic" | "total ion current" => Some(ChromatogramType::TIC),
            _ => None,
        }
    }

    pub fn list() -> Vec<&'static str> {
        vec!["BPI", "TIC"]
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            ChromatogramType::BPI => "BPI",
            ChromatogramType::TIC => "TIC",
        }
    }
}

/// Interpolate chromatogram data to create smoother lines - enhanced combination approach
fn interpolate_chromatogram(data: &[f32]) -> Vec<f32> {
    let mut interpolated = data.to_vec();
    
    // Step 1: Fill gaps between data points using linear interpolation
    let mut last_valid_idx = None;
    let mut last_valid_value = 0.0;
    
    for i in 0..interpolated.len() {
        if interpolated[i] > 0.0 {
            // Found a data point
            if let Some(prev_idx) = last_valid_idx {
                // Interpolate between previous and current point
                let steps = i - prev_idx;
                if steps > 1 {
                    let step_size = (interpolated[i] - last_valid_value) / steps as f32;
                    for j in 1..steps {
                        interpolated[prev_idx + j] = last_valid_value + (step_size * j as f32);
                    }
                }
            }
            last_valid_idx = Some(i);
            last_valid_value = interpolated[i];
        }
    }
    
    // Step 2: Create a denser dataset by upsampling (adding more points)
    let upsample_factor = 3; // Triple the density
    let mut upsampled = Vec::with_capacity(interpolated.len() * upsample_factor);
    
    for i in 0..interpolated.len() - 1 {
        upsampled.push(interpolated[i]);
        
        // Add intermediate points between each pair of consecutive points
        for j in 1..upsample_factor {
            let t = j as f32 / upsample_factor as f32;
            let interpolated_value = interpolated[i] * (1.0 - t) + interpolated[i + 1] * t;
            upsampled.push(interpolated_value);
        }
    }
    // Add the final point
    if let Some(&last) = interpolated.last() {
        upsampled.push(last);
    }
    
    // Step 3: Apply aggressive smoothing to the upsampled data
    let smoothed = smooth_chromatogram(&upsampled, 10); // 10-point moving average for aggressive smoothing

    // Step 4: Downsample back to original resolution but keep the smoothness
    let mut final_result = Vec::with_capacity(interpolated.len());
    for i in 0..interpolated.len() {
        let upsampled_idx = i * upsample_factor;
        if upsampled_idx < smoothed.len() {
            final_result.push(smoothed[upsampled_idx]);
        } else {
            final_result.push(interpolated[i]);
        }
    }
    
    final_result
}

/// Smooth chromatogram data using moving average smoothing
fn smooth_chromatogram(data: &[f32], window_size: usize) -> Vec<f32> {
    if data.is_empty() || window_size == 0 {
        return data.to_vec();
    }
    
    let half_window = window_size / 2;
    let mut smoothed = Vec::with_capacity(data.len());
    
    for i in 0..data.len() {
        let start = i.saturating_sub(half_window);
        let end = (i + half_window + 1).min(data.len());
        let sum: f32 = data[start..end].iter().sum();
        let avg = sum / (end - start) as f32;
        smoothed.push(avg);
    }
    
    smoothed
}

/// Normalize data to 0-1 range
fn normalize_data(data: &[f32]) -> Vec<f32> {
    if data.is_empty() {
        return Vec::new();
    }
    
    let max_val = data.iter().cloned().fold(0.0f32, f32::max);
    let min_val = data.iter().cloned().fold(f32::INFINITY, f32::min);
    
    if max_val == min_val {
        return vec![0.0; data.len()];
    }
    
    data.iter()
        .map(|&x| (x - min_val) / (max_val - min_val))
        .collect()
}

/// Parse scan IDs to find the min and max numeric scan IDs
fn get_scan_range(scan_ids: &[String]) -> (usize, usize) {
    let mut min_scan = usize::MAX;
    let mut max_scan = 0;
    
    for scan_id in scan_ids {
        if let Ok(scan_num) = scan_id.parse::<usize>() {
            min_scan = min_scan.min(scan_num);
            max_scan = max_scan.max(scan_num);
        }
    }
    
    (min_scan, max_scan)
}

/// Create a mapping from scan ID to its position in the full range
fn create_scan_position_map(scan_ids: &[String], min_scan: usize) -> HashMap<String, usize> {
    let mut position_map = HashMap::new();
    
    for scan_id in scan_ids {
        if let Ok(scan_num) = scan_id.parse::<usize>() {
            position_map.insert(scan_id.clone(), scan_num - min_scan);
        }
    }
    
    position_map
}

/// Supported output formats
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub enum ImageFormat {
    Png,
    Svg,
    Jpeg,
}

impl ImageFormat {
    pub fn from_ext(ext: &str) -> Option<Self> {
        match ext.to_lowercase().as_str() {
            "png" => Some(ImageFormat::Png),
            "svg" => Some(ImageFormat::Svg),
            "jpg" | "jpeg" => Some(ImageFormat::Jpeg),
            _ => None,
        }
    }
}

/// A serializable color type that can be converted to/from RGBColor
#[derive(Clone, Copy, Debug, serde::Serialize, serde::Deserialize)]
pub struct SerializableColor {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

impl SerializableColor {
    pub fn new(r: u8, g: u8, b: u8) -> Self {
        Self { r, g, b }
    }

    pub fn to_rgb_color(&self) -> RGBColor {
        RGBColor(self.r, self.g, self.b)
    }

    pub fn from_rgb_color(color: RGBColor) -> Self {
        Self {
            r: color.0,
            g: color.1,
            b: color.2,
        }
    }
}

#[derive(Clone, Copy, Debug, serde::Serialize, serde::Deserialize)]
pub struct ColorTheme {
    pub background: SerializableColor,
    pub low: SerializableColor,
    pub high: SerializableColor,
}

#[derive(serde::Serialize, serde::Deserialize, Debug, Clone, PartialEq, Eq, Hash)]
pub enum ThemeName {
    Classic,
    DarkBlue,
    Inferno,
    Viridis,
    Plasma,
    GrayScale,
    WhiteRed,
    WhiteBlue,
    WhiteGreen,
    WhiteOrange,
    WhitePurple,
    BlackYellow,
    Jet,
}

impl ThemeName {
    pub fn from_str(input: &str) -> Option<Self> {
        match input.trim().to_lowercase().as_str() {
            "classic" => Some(ThemeName::Classic),
            "darkblue" => Some(ThemeName::DarkBlue),
            "inferno" => Some(ThemeName::Inferno),
            "viridis" => Some(ThemeName::Viridis),
            "plasma" => Some(ThemeName::Plasma),
            "grayscale" => Some(ThemeName::GrayScale),
            "whitered" | "white-red" => Some(ThemeName::WhiteRed),
            "whiteblue" | "white-blue" => Some(ThemeName::WhiteBlue),
            "whitegreen" | "white-green" => Some(ThemeName::WhiteGreen),
            "whiteorange" | "white-orange" => Some(ThemeName::WhiteOrange),
            "whitepurple" | "white-purple" => Some(ThemeName::WhitePurple),
            "blackyellow" | "black-yellow" => Some(ThemeName::BlackYellow),
            "jet" => Some(ThemeName::Jet),
            _ => None,
        }
    }

    pub fn list() -> Vec<&'static str> {
        vec![
            "classic", "darkblue", "inferno", "viridis", "plasma", "grayscale",
            "whitered", "whiteblue", "whitegreen", "whiteorange", "whitepurple",
            "blackyellow", "jet"
        ]
    }

    pub fn get_theme(&self) -> ColorTheme {
        match self {
            ThemeName::Classic => ColorTheme {
                background: SerializableColor::new(255, 255, 255),
                low: SerializableColor::new(255, 255, 255),
                high: SerializableColor::new(0, 0, 255),
            },
            ThemeName::DarkBlue => ColorTheme {
                background: SerializableColor::new(20, 20, 30),
                low: SerializableColor::new(30, 30, 50),
                high: SerializableColor::new(0, 200, 255),
            },
            ThemeName::Inferno => ColorTheme {
                background: SerializableColor::new(10, 10, 10),
                low: SerializableColor::new(40, 0, 20),
                high: SerializableColor::new(255, 200, 50),
            },
            ThemeName::Viridis => ColorTheme {
                background: SerializableColor::new(20, 20, 20),
                low: SerializableColor::new(68, 1, 84),
                high: SerializableColor::new(253, 231, 37),
            },
            ThemeName::Plasma => ColorTheme {
                background: SerializableColor::new(30, 5, 40),
                low: SerializableColor::new(13, 8, 135),
                high: SerializableColor::new(240, 249, 33),
            },
            ThemeName::GrayScale => ColorTheme {
                background: SerializableColor::new(255, 255, 255),
                low: SerializableColor::new(255, 255, 255),
                high: SerializableColor::new(50, 50, 50),
            },
            // New white background themes
            ThemeName::WhiteRed => ColorTheme {
                background: SerializableColor::new(255, 255, 255),
                low: SerializableColor::new(255, 255, 255),
                high: SerializableColor::new(220, 20, 20),
            },
            ThemeName::WhiteBlue => ColorTheme {
                background: SerializableColor::new(255, 255, 255),
                low: SerializableColor::new(255, 255, 255),
                high: SerializableColor::new(20, 20, 220),
            },
            ThemeName::WhiteGreen => ColorTheme {
                background: SerializableColor::new(255, 255, 255),
                low: SerializableColor::new(255, 255, 255),
                high: SerializableColor::new(20, 180, 20),
            },
            ThemeName::WhiteOrange => ColorTheme {
                background: SerializableColor::new(255, 255, 255),
                low: SerializableColor::new(255, 255, 255),
                high: SerializableColor::new(255, 140, 20),
            },
            ThemeName::WhitePurple => ColorTheme {
                background: SerializableColor::new(255, 255, 255),
                low: SerializableColor::new(255, 255, 255),
                high: SerializableColor::new(160, 20, 200),
            },
            // Black background theme
            ThemeName::BlackYellow => ColorTheme {
                background: SerializableColor::new(0, 0, 0),
                low: SerializableColor::new(40, 40, 0),
                high: SerializableColor::new(255, 255, 20),
            },
            // Jet colormap (blue -> cyan -> yellow -> red)
            ThemeName::Jet => ColorTheme {
                background: SerializableColor::new(255, 255, 255),
                low: SerializableColor::new(0, 0, 140),
                high: SerializableColor::new(140, 0, 0),
            },
        }
    }
}

/// Generic function that works with any backend - renders heatmap with full scan range
fn render_heatmap_full_range<B>(
    backend: B,
    scan_ids: &[String],
    similarity_matrix: &Array2<f32>,
    chromatogram_data: &[f32],
    _chromatogram_type: ChromatogramType,
    enable_smoothing: bool,
    theme: &ColorTheme,
    is_jet_theme: bool,
    global_scan_range: Option<(usize, usize)>, // Optional global range override
    chromatogram_scan_ids: Option<&[String]>, // Optional scan IDs that correspond to chromatogram_data
) -> Result<(), Box<dyn std::error::Error>> 
where
    B: DrawingBackend,
    B::ErrorType: 'static + Send + Sync,
{
    let root = backend.into_drawing_area();
    root.fill(&theme.background.to_rgb_color())?;

    let n = similarity_matrix.nrows();
    assert_eq!(n, similarity_matrix.ncols(), "Matrix must be square");
    assert_eq!(n, scan_ids.len(), "Mismatch between matrix and labels");
    
    // Note: chromatogram_data may be longer than scan_ids when using MS1 chromatogram for MS2 heatmap
    // This is intentional for consistent chromatogram display across heatmaps

    // Get the scan range - use global range if provided, otherwise calculate from scan_ids
    let (min_scan, max_scan) = if let Some((global_min, global_max)) = global_scan_range {
        (global_min, global_max)
    } else {
        get_scan_range(scan_ids)
    };
    let full_range = max_scan - min_scan + 1;
    
    // Create position mapping for available scans
    let position_map = create_scan_position_map(scan_ids, min_scan);
    
    // println!("||    Heatmap will show full scan range: {} to {} ({} positions)", min_scan, max_scan, full_range);
    // println!("||    Available scans: {} (showing {} data points)", n, n);

    // Prepare chromatogram data: apply user-controlled interpolated smoothing if enabled, then normalize
    let processed_data = if enable_smoothing {
        interpolate_chromatogram(chromatogram_data) // User-controlled interpolated smoothing
    } else {
        chromatogram_data.to_vec()
    };
    let normalized_data = normalize_data(&processed_data);

    // Create full-range chromatogram data (with zeros for missing scans)
    let mut full_chromatogram = vec![0.0; full_range];
    
    // Handle two cases:
    // 1. Normal case: chromatogram_data matches scan_ids (same length)
    // 2. Cross-reference case: chromatogram_data is from MS1, scan_ids are from MS2
    if normalized_data.len() == scan_ids.len() {
        // Normal case: direct mapping
        for (i, scan_id) in scan_ids.iter().enumerate() {
            if let Some(&pos) = position_map.get(scan_id) {
                full_chromatogram[pos] = normalized_data[i];
            }
        }
    } else {
        // Cross-reference case: chromatogram_data is from MS1, map using provided scan IDs
        if let Some(chrom_scan_ids) = chromatogram_scan_ids {
            // Use the provided scan IDs that correspond to chromatogram_data
            for (i, chrom_scan_id) in chrom_scan_ids.iter().enumerate() {
                if let Ok(scan_num) = chrom_scan_id.parse::<usize>() {
                    if scan_num >= min_scan && scan_num <= max_scan && i < normalized_data.len() {
                        let pos = scan_num - min_scan;
                        full_chromatogram[pos] = normalized_data[i];
                    }
                }
            }
        } else {
            // Fallback: Create a map from scan number to chromatogram intensity
            let chromatogram_map: HashMap<usize, f32> = (0..normalized_data.len())
                .filter_map(|i| {
                    let scan_num = min_scan + i; // Reconstruct scan number from index
                    if i < normalized_data.len() {
                        Some((scan_num, normalized_data[i]))
                    } else {
                        None
                    }
                })
                .collect();
            
            // Fill the full_chromatogram using the map
            for scan_num in min_scan..=max_scan {
                let pos = scan_num - min_scan;
                if let Some(&intensity) = chromatogram_map.get(&scan_num) {
                    full_chromatogram[pos] = intensity;
                }
            }
        }
    }

    // Generate caption based on chromatogram type and smoothing
    // let caption = format!("Similarity Matrix with {} {} (Scans {}-{})", 
    //     chromatogram_type.as_str(),
    //     if enable_smoothing { "(Smoothed)" } else { "" },
    //     min_scan,
    //     max_scan
    // );

    // Create the main heatmap using expanded coordinates to accommodate chromatograms
    let chromatogram_scale = 250.0; // Scale factor for chromatogram
    let chromatogram_margin = 250; // Fixed reasonable margin for chromatogram
    let chart_width = full_range + chromatogram_margin;
    let chart_height = full_range + chromatogram_margin;
    
    let mut chart = ChartBuilder::on(&root)
        .margin(60)
        // .caption(&caption, ("sans-serif", 24))
        .build_cartesian_2d(0..chart_width, 0..chart_height)?;

    chart.configure_mesh()
        .disable_mesh()
        .x_labels(8)
        .y_labels(8)
        .x_label_formatter(&|x| {
            if *x < full_range {
                format!("{}", min_scan + x)
            } else {
                "".to_string() // Don't show labels in chromatogram area
            }
        })
        .y_label_formatter(&|y| {
            if *y < full_range {
                format!("{}", min_scan + y)
            } else {
                "".to_string() // Don't show labels in chromatogram area
            }
        })
        .draw()?;

    // Draw the main heatmap - use rectangles with transparency for traditional heatmap appearance
    // Set rectangle size to 1% of heatmap length with minimum of 1
    let rect_size = ((full_range as f64 * 0.01).round() as usize).max(1);
    let alpha = 0.7; // Transparency level (0.0 = fully transparent, 1.0 = fully opaque)
    
    // println!("||    Using rectangle size: {} (1% of range {} = {})", rect_size, full_range, full_range as f64 * 0.01);
    
    for (i, scan_id_i) in scan_ids.iter().enumerate() {
        if let Some(&pos_i) = position_map.get(scan_id_i) {
            for (j, scan_id_j) in scan_ids.iter().enumerate() {
                if let Some(&pos_j) = position_map.get(scan_id_j) {
                    let value = similarity_matrix[(i, j)];
                    
                    // Skip very low similarity values to improve performance and visibility
                    if value < 0.1 {
                        continue;
                    }

                    // Skip the first 20 scans
                    if pos_i < 20 || pos_j < 20 {
                        continue;
                    }

                    let base_color = if is_jet_theme {
                        interpolate_jet_colormap(value as f64)
                    } else {
                        interpolate_rgb(value as f64, &theme.low.to_rgb_color(), &theme.high.to_rgb_color())
                    };
                    
                    // Convert to RGBA with alpha transparency
                    let transparent_color = RGBAColor(
                        base_color.0, 
                        base_color.1, 
                        base_color.2, 
                        alpha
                    );
                    
                    let y_coord = full_range.saturating_sub(pos_i);
                    let x_right = (pos_j + rect_size).min(full_range);
                    let y_bottom = if pos_i + rect_size < full_range {
                        full_range - pos_i - rect_size
                    } else {
                        0
                    };
                    
                    // Use rectangles with transparency for traditional heatmap look
                    chart.draw_series(std::iter::once(Rectangle::new(
                        [(pos_j, y_bottom), (x_right, y_coord)],
                        transparent_color.filled(),
                    )))?;
                }
            }
        }
    }

    // Add chromatogram as overlay lines at the edges of the heatmap
    // Plot chromatograms OUTWARD from the heatmap for better visibility
    let chromatogram_offset = 10; // Small offset to separate chromatogram from heatmap edge
    
    // Apply interpolation for smoother chromatogram lines
    let interpolated_chromatogram = interpolate_chromatogram(&full_chromatogram);
    
    // Simple line color based on background brightness
    let bg_avg = (theme.background.r as u16 + theme.background.g as u16 + theme.background.b as u16) / 3;
    let line_color = if bg_avg < 128 {
        RGBColor(255, 255, 255) // White for dark backgrounds
    } else {
        RGBColor(0, 0, 0) // Black for light backgrounds
    };
    
    // Top edge: horizontal chromatogram (plotted above the heatmap)
    let top_chromatogram_data: Vec<(usize, usize)> = interpolated_chromatogram
        .iter()
        .enumerate()
        .map(|(i, &val)| (i, full_range + chromatogram_offset + (val * chromatogram_scale) as usize))
        .collect();
    
    if !top_chromatogram_data.is_empty() {
        chart.draw_series(LineSeries::new(top_chromatogram_data, line_color.stroke_width(1)))?;
    }

    // Right edge: vertical chromatogram (plotted to the right of the heatmap, rotated 90°)
    let right_chromatogram_data: Vec<(usize, usize)> = interpolated_chromatogram
        .iter()
        .enumerate()
        .map(|(i, &val)| (full_range + chromatogram_offset + (val * chromatogram_scale) as usize, full_range - i - 1))
        .collect();
    
    if !right_chromatogram_data.is_empty() {
        chart.draw_series(LineSeries::new(right_chromatogram_data, line_color.stroke_width(1)))?;
    }

    root.present()?;
    
    Ok(())
}

/// Main plotting function - handles lifetime issues properly
pub fn plot_similarity_heatmap(
    scan_ids: &[String],
    similarity_matrix: &Array2<f32>,
    chromatogram_data: &[f32],
    chromatogram_type: ChromatogramType,
    enable_smoothing: bool,
    output_file: &str,
    image_format: &ImageFormat,
    theme: &ColorTheme,
    global_scan_range: Option<(usize, usize)>, // Optional global range for consistent scaling
    chromatogram_scan_ids: Option<Vec<String>>, // Optional scan IDs that correspond to chromatogram_data
) -> Result<(), Box<dyn std::error::Error>> {
    let image_size = (1024, 1024);
    
    // Create owned string to satisfy lifetime requirements
    let output_path = output_file.to_owned();
    
    // Check if this is the Jet theme (special handling needed)
    let is_jet_theme = theme.low.r == 0 && theme.low.g == 0 && theme.low.b == 140 &&
                      theme.high.r == 140 && theme.high.g == 0 && theme.high.b == 0;

    match image_format {
        ImageFormat::Png | ImageFormat::Jpeg => {
            // Create backend with owned string
            let backend = BitMapBackend::new(output_path.as_str(), image_size);
            render_heatmap_full_range(backend, scan_ids, similarity_matrix, chromatogram_data, chromatogram_type, enable_smoothing, theme, is_jet_theme, global_scan_range, chromatogram_scan_ids.as_deref())?;
        }
        ImageFormat::Svg => {
            // Create backend with owned string  
            let backend = SVGBackend::new(output_path.as_str(), image_size);
            render_heatmap_full_range(backend, scan_ids, similarity_matrix, chromatogram_data, chromatogram_type, enable_smoothing, theme, is_jet_theme, global_scan_range, chromatogram_scan_ids.as_deref())?;
        }
    }

    // println!("||    Heatmap saved to {output_file}");
    Ok(())
}

fn interpolate_rgb(score: f64, low: &RGBColor, high: &RGBColor) -> RGBColor {
    let clamped = score.clamp(0.0, 1.0);
    let r = (low.0 as f64 * (1.0 - clamped) + high.0 as f64 * clamped).round() as u8;
    let g = (low.1 as f64 * (1.0 - clamped) + high.1 as f64 * clamped).round() as u8;
    let b = (low.2 as f64 * (1.0 - clamped) + high.2 as f64 * clamped).round() as u8;
    RGBColor(r, g, b)
}

fn interpolate_jet_colormap(score: f64) -> RGBColor {
    let clamped = score.clamp(0.0, 1.0);
    
    // Jet colormap: Blue -> Cyan -> Yellow -> Red
    if clamped < 0.25 {
        // Blue to Cyan
        let t = clamped * 4.0;
        RGBColor(0, (255.0 * t) as u8, 255)
    } else if clamped < 0.5 {
        // Cyan to Green  
        let t = (clamped - 0.25) * 4.0;
        RGBColor(0, 255, (255.0 * (1.0 - t)) as u8)
    } else if clamped < 0.75 {
        // Green to Yellow
        let t = (clamped - 0.5) * 4.0;
        RGBColor((255.0 * t) as u8, 255, 0)
    } else {
        // Yellow to Red
        let t = (clamped - 0.75) * 4.0;
        RGBColor(255, (255.0 * (1.0 - t)) as u8, 0)
    }
}

#[allow(dead_code)]
fn draw_color_legend<D>(
    area: &DrawingArea<D, plotters::coord::Shift>,
    low: &SerializableColor,
    high: &SerializableColor,
) -> Result<(), Box<dyn std::error::Error>> 
where
    D: DrawingBackend,
    D::ErrorType: 'static + Send + Sync,
{
    let bar_width = 300;
    
    let mut legend = ChartBuilder::on(area)
        .margin(1)
        .build_cartesian_2d(0..bar_width, 0..15)?;

    legend.configure_mesh().disable_mesh().draw()?;

    // Draw a very thin horizontal color gradient bar (15px total height)
    legend.draw_series((0..bar_width).map(|x| {
        let v = x as f64 / (bar_width - 1) as f64;
        let color = interpolate_rgb(v, &low.to_rgb_color(), &high.to_rgb_color());
        Rectangle::new([(x, 2), (x + 1, 8)], color.filled())
    }))?;

    // Add minimal labels
    legend.draw_series(std::iter::once(Text::new("0.0", (5, 9), ("sans-serif", 9))))?;
    legend.draw_series(std::iter::once(Text::new("1.0", (bar_width - 20, 9), ("sans-serif", 9))))?;
    legend.draw_series(std::iter::once(Text::new("Similarity", (bar_width / 2 - 20, 9), ("sans-serif", 9))))?;

    Ok(())
}