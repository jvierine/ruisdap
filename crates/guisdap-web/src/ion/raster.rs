use crate::{doc, el, plotting};
use wasm_bindgen::{JsCast, JsValue};
use web_sys::{CanvasRenderingContext2d, HtmlCanvasElement};
pub fn draw(
    id: &str,
    values: &[f32],
    flags: &[u8],
    width: usize,
    height: usize,
    spec: &serde_json::Value,
) -> Result<(), JsValue> {
    plotting::shared_field(id, values, width, height, spec)?;
    let overlay_id = format!("{id}-quality");
    if doc().get_element_by_id(&overlay_id).is_none() {
        let c = doc().create_element("canvas")?;
        c.set_id(&overlay_id);
        c.set_attribute("style", "pointer-events:none;background:transparent")?;
        el(id).parent_element().unwrap().append_child(&c)?;
    }
    let c: HtmlCanvasElement = el(&overlay_id).dyn_into()?;
    c.set_width(width as u32);
    c.set_height(height as u32);
    let ctx: CanvasRenderingContext2d = c.get_context("2d")?.unwrap().dyn_into()?;
    ctx.clear_rect(0., 0., width as f64, height as f64);
    ctx.set_fill_style_str("#cccccc38");
    for iy in 0..height {
        for ix in 0..width {
            if flags[iy * width + ix] != 0 && (ix + height - 1 - iy) % 16 == 0 {
                ctx.fill_rect(ix as f64, (height - 1 - iy) as f64, 1., 1.);
            }
        }
    }
    Ok(())
}
