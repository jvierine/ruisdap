use crate::{doc, error};
use wasm_bindgen::{JsCast, JsValue};
pub fn setup() -> Result<(), JsValue> {
    let location = web_sys::window().unwrap().location();
    let path = location.pathname()?;
    let processing = path.ends_with("live.html");
    if let Some(link) = doc().query_selector(if processing {
        "[data-nav=processing]"
    } else {
        "[data-nav=browser]"
    })? {
        link.set_attribute("aria-current", "page")?;
    }
    let query = location.search()?;
    let experiment = query
        .trim_start_matches('?')
        .split('&')
        .find(|s| s.starts_with("experiment="))
        .unwrap_or("");
    let links = doc().query_selector_all("[data-dataset-view]")?;
    for i in 0..links.length() {
        let link = links
            .item(i)
            .unwrap()
            .dyn_into::<web_sys::Element>()
            .map_err(|_| error("navigation element"))?;
        let route = if link.get_attribute("data-dataset-view").as_deref() == Some("ion") {
            "ion.html"
        } else {
            "plasma.html"
        };
        link.set_attribute(
            "href",
            &format!(
                "{route}{}",
                if experiment.is_empty() {
                    String::new()
                } else {
                    format!("?{experiment}")
                }
            ),
        )?;
    }
    Ok(())
}
