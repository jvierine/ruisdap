use super::*;

pub(super) fn spectrum(v: &Viewer) -> String {
    let g = &selected(v).gates[v.gate];

    let frequency: Vec<_> = (0..=800).map(|i| -40000. + i as f64 * 100.).collect();
    let lags: Vec<_> = g.lag_us.iter().map(|t| t * 1e-6).collect();
    // BEATA's exported lag grid is 10 us. Reject incompatible products rather
    // than silently deriving a different normalization from a sparse subset.
    let estimate = guisdap_ion::dft::AcfDft::new(&lags, &frequency, 10e-6).and_then(|op| {
        let values = op.transform(&g.observed)?;
        let errors = op.diagonal_standard_error(&g.sigma)?;
        Ok((op, values, errors))
    });
    let Ok((op, values, errors)) = estimate else {
        return empty_curve("Measured ACF spectrum unavailable");
    };
    let fitted = g
        .result
        .as_ref()
        .and_then(|f| op.transform(&f.predicted).ok());
    let low = values
        .iter()
        .zip(&errors)
        .map(|(a, b)| a - b)
        .chain(fitted.iter().flatten().copied())
        .fold(0., f64::min);
    let high = values
        .iter()
        .zip(&errors)
        .map(|(a, b)| a + b)
        .chain(fitted.iter().flatten().copied())
        .fold(1., f64::max);
    let pad = (high - low) * 0.06;
    let scale = 10f64.powf(high.abs().max(low.abs()).max(1.).log10().floor());
    let series:Vec<_>=[(Some(&values),"#647b89"),(fitted.as_ref(),"#236b8e")].into_iter().filter_map(|(data,color)|data.map(|data|serde_json::json!({"xy":frequency.iter().zip(data).flat_map(|(f,v)|[*f/1000.,*v/scale]).collect::<Vec<_>>(),"color":color}))).collect();
    let band: Vec<_> = (0..values.len())
        .flat_map(|i| {
            [
                frequency[i] / 1000.,
                (values[i] - errors[i]) / scale,
                (values[i] + errors[i]) / scale,
            ]
        })
        .collect();
    serde_json::json!({"title":"Measured ACF spectrum","bounds":[-40.,40.,(low-pad)/scale,(high+pad)/scale],"xlabel":"Frequency offset / kHz","ylabel":format!("ACF DFT / {scale:.0e} m⁻³ Hz⁻¹"),"series":series,"bands":[{"values":band,"color":"#647b894d"}],"caption":format!("{} measured lags · {} missing · {}",op.unique_lags,op.missing_lags,if op.has_zero_lag{"zero lag measured"}else{"baseline undetermined"})}).to_string()
}
