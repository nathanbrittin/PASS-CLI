use plotters::prelude::*;
use ndarray::Array2;
use plotters::style::RGBColor;

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
    root.fill(&theme.background.to_rgb_color())?;

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
                interpolate_rgb(value as f64, &theme.low.to_rgb_color(), &theme.high.to_rgb_color())
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
    let is_jet_theme = theme.low.r == 0 && theme.low.g == 0 && theme.low.b == 140 &&
                      theme.high.r == 140 && theme.high.g == 0 && theme.high.b == 0;

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