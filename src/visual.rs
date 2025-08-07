use plotters::prelude::*;
use ndarray::Array2;
use plotters::style::RGBColor;

/// Supported output formats
#[derive(Debug, Clone)]
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

pub struct ColorTheme {
    pub background: RGBColor,
    pub low: RGBColor,
    pub high: RGBColor,
}

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
                background: RGBColor(255, 255, 255),
                low: RGBColor(255, 255, 255),
                high: RGBColor(0, 0, 255),
            },
            ThemeName::DarkBlue => ColorTheme {
                background: RGBColor(20, 20, 30),
                low: RGBColor(30, 30, 50),
                high: RGBColor(0, 200, 255),
            },
            ThemeName::Inferno => ColorTheme {
                background: RGBColor(10, 10, 10),
                low: RGBColor(40, 0, 20),
                high: RGBColor(255, 200, 50),
            },
            ThemeName::Viridis => ColorTheme {
                background: RGBColor(20, 20, 20),
                low: RGBColor(68, 1, 84),
                high: RGBColor(253, 231, 37),
            },
            ThemeName::Plasma => ColorTheme {
                background: RGBColor(30, 5, 40),
                low: RGBColor(13, 8, 135),
                high: RGBColor(240, 249, 33),
            },
            ThemeName::GrayScale => ColorTheme {
                background: RGBColor(255, 255, 255),
                low: RGBColor(255, 255, 255),
                high: RGBColor(50, 50, 50),
            },
            // New white background themes
            ThemeName::WhiteRed => ColorTheme {
                background: RGBColor(255, 255, 255),
                low: RGBColor(255, 255, 255),
                high: RGBColor(220, 20, 20),
            },
            ThemeName::WhiteBlue => ColorTheme {
                background: RGBColor(255, 255, 255),
                low: RGBColor(255, 255, 255),
                high: RGBColor(20, 20, 220),
            },
            ThemeName::WhiteGreen => ColorTheme {
                background: RGBColor(255, 255, 255),
                low: RGBColor(255, 255, 255),
                high: RGBColor(20, 180, 20),
            },
            ThemeName::WhiteOrange => ColorTheme {
                background: RGBColor(255, 255, 255),
                low: RGBColor(255, 255, 255),
                high: RGBColor(255, 140, 20),
            },
            ThemeName::WhitePurple => ColorTheme {
                background: RGBColor(255, 255, 255),
                low: RGBColor(255, 255, 255),
                high: RGBColor(160, 20, 200),
            },
            // Black background theme
            ThemeName::BlackYellow => ColorTheme {
                background: RGBColor(0, 0, 0),
                low: RGBColor(40, 40, 0),
                high: RGBColor(255, 255, 20),
            },
            // Jet colormap (blue -> cyan -> yellow -> red)
            ThemeName::Jet => ColorTheme {
                background: RGBColor(255, 255, 255),
                low: RGBColor(0, 0, 140),
                high: RGBColor(140, 0, 0),
            },
        }
    }
}

/// Generic function that works with any backend
fn render_heatmap<B>(
    backend: B,
    scan_ids: &[String],
    similarity_matrix: &Array2<f32>,
    theme: &ColorTheme,
    is_jet_theme: bool,
) -> Result<(), Box<dyn std::error::Error>> 
where
    B: DrawingBackend,
    B::ErrorType: 'static + Send + Sync,
{
    let root = backend.into_drawing_area();
    root.fill(&theme.background)?;

    let n = similarity_matrix.nrows();
    assert_eq!(n, similarity_matrix.ncols(), "Matrix must be square");
    assert_eq!(n, scan_ids.len(), "Mismatch between matrix and labels");

    // Temporarily remove legend - use full area for heatmap
    let mut chart = ChartBuilder::on(&root)
        .margin(50)
        .caption("Similarity Matrix Heatmap", ("sans-serif", 30))
        .build_cartesian_2d(0..n, 0..n)?;

    chart.configure_mesh()
        .disable_mesh()
        .x_labels(10)
        .y_labels(10)
        .x_label_formatter(&|x| scan_ids.get(*x).unwrap_or(&"?".to_string()).to_string())
        .y_label_formatter(&|y| scan_ids.get(*y).unwrap_or(&"?".to_string()).to_string())
        .draw()?;

    for row in 0..n {
        for col in 0..n {
            let value = similarity_matrix[(row, col)];
            let color = if is_jet_theme {
                interpolate_jet_colormap(value as f64)
            } else {
                interpolate_rgb(value as f64, &theme.low, &theme.high)
            };
            chart.draw_series(std::iter::once(Rectangle::new(
                [(col, n - row - 1), (col + 1, n - row)],
                color.filled(),
            )))?;
        }
    }

    // Commented out the legend for now
    // draw_color_legend(&legend_area, &theme.low, &theme.high)?;
    root.present()?;
    
    Ok(())
}

/// Main plotting function - handles lifetime issues properly
pub fn plot_similarity_heatmap(
    scan_ids: &[String],
    similarity_matrix: &Array2<f32>,
    output_file: &str,
    image_format: &ImageFormat,
    theme: &ColorTheme,
) -> Result<(), Box<dyn std::error::Error>> {
    let image_size = (1024, 1024);
    
    // Create owned string to satisfy lifetime requirements
    let output_path = output_file.to_owned();
    
    // Check if this is the Jet theme (special handling needed)
    let is_jet_theme = theme.low.0 == 0 && theme.low.1 == 0 && theme.low.2 == 140 &&
                      theme.high.0 == 140 && theme.high.1 == 0 && theme.high.2 == 0;

    match image_format {
        ImageFormat::Png | ImageFormat::Jpeg => {
            // Create backend with owned string
            let backend = BitMapBackend::new(output_path.as_str(), image_size);
            render_heatmap(backend, scan_ids, similarity_matrix, theme, is_jet_theme)?;
        }
        ImageFormat::Svg => {
            // Create backend with owned string  
            let backend = SVGBackend::new(output_path.as_str(), image_size);
            render_heatmap(backend, scan_ids, similarity_matrix, theme, is_jet_theme)?;
        }
    }

    println!("||    Heatmap saved to {output_file}");
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

fn draw_color_legend<D>(
    area: &DrawingArea<D, plotters::coord::Shift>,
    low: &RGBColor,
    high: &RGBColor,
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
        let color = interpolate_rgb(v, low, high);
        Rectangle::new([(x, 2), (x + 1, 8)], color.filled())
    }))?;

    // Add minimal labels
    legend.draw_series(std::iter::once(Text::new("0.0", (5, 9), ("sans-serif", 9))))?;
    legend.draw_series(std::iter::once(Text::new("1.0", (bar_width - 20, 9), ("sans-serif", 9))))?;
    legend.draw_series(std::iter::once(Text::new("Similarity", (bar_width / 2 - 20, 9), ("sans-serif", 9))))?;

    Ok(())
}