//! RUISDAP adapters for the shared native Rustplotlib browser components.
use wasm_bindgen::JsValue;
#[wasm_bindgen::prelude::wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace=window, js_name=ruisdapField, catch)]
    fn native_field(
        id: &str,
        values: &[f32],
        nx: usize,
        ny: usize,
        spec: JsValue,
    ) -> Result<JsValue, JsValue>;
    #[wasm_bindgen(js_namespace=window, js_name=ruisdapCurve, catch)]
    fn native_curve(id: &str, spec: JsValue) -> Result<(), JsValue>;
}
pub fn shared_field(
    id: &str,
    values: &[f32],
    nx: usize,
    ny: usize,
    spec: &serde_json::Value,
) -> Result<(), JsValue> {
    native_field(id, values, nx, ny, js_sys::JSON::parse(&spec.to_string())?)?;
    Ok(())
}
pub fn shared_curve(id: &str, spec: &str) {
    if let Err(e) = js_sys::JSON::parse(spec).and_then(|v| native_curve(id, v)) {
        web_sys::console::error_1(&e);
    }
}

pub fn rectilinear(id: &str, values: &[f32], x: &[f64], y: &[f64]) -> Result<(), JsValue> {
    fn edges(a: &[f64], unit: f64) -> Vec<f64> {
        if a.len() == 1 {
            return vec![(a[0] - 0.5) / unit, (a[0] + 0.5) / unit];
        }
        std::iter::once((a[0] - (a[1] - a[0]) / 2.) / unit)
            .chain(a.windows(2).map(|w| (w[0] + w[1]) / 2. / unit))
            .chain(std::iter::once(
                (a[a.len() - 1] + (a[a.len() - 1] - a[a.len() - 2]) / 2.) / unit,
            ))
            .collect()
    }
    let xe = edges(x, 1e6);
    let ye = edges(y, 1000.);
    let limit = crate::percentile(values, 0.995) as f64;
    shared_field(
        id,
        values,
        x.len(),
        y.len(),
        &serde_json::json!({"bounds":[xe[0],xe[xe.len()-1],ye[0],ye[ye.len()-1]],"x_edges":xe,"y_edges":ye,"limits":[-limit,limit],"mode":2,"palette":1,"time":false,"external_navigation":true,"xlabel":"Frequency offset / MHz","ylabel":"Slant range / km"}),
    )
}

// These constructors keep the library's wasm-bindgen exports linked even though
// application drawing is called through the shared JavaScript component.
#[wasm_bindgen::prelude::wasm_bindgen]
pub fn native_field_view(
    canvas: web_sys::HtmlCanvasElement,
) -> Result<rustplotlib_web::FieldView, JsValue> {
    rustplotlib_web::FieldView::new(canvas)
}
#[wasm_bindgen::prelude::wasm_bindgen]
pub fn native_series_view(
    canvas: web_sys::HtmlCanvasElement,
) -> Result<rustplotlib_web::SeriesView, JsValue> {
    rustplotlib_web::SeriesView::new(canvas)
}
#[wasm_bindgen::prelude::wasm_bindgen]
pub fn native_plot_viewport() -> rustplotlib_web::PlotViewport {
    rustplotlib_web::PlotViewport::new()
}
