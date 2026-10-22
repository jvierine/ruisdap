//! Conditional-candidate stress test: zero-mean Gaussian ion-line scatter,
//! without any coherent target. Not an end-to-end CFAR false-alarm estimate.
use anyhow::{ensure, Result};
use guisdap_ion::spectrum::{IonSpecies, PlasmaState, SpectrumModel};
use guisdap_millstone::{
    debris::{cancel, template, Detection, Settings},
    C,
};
use rustfft::FftPlanner;
fn noise(n: usize, seed: u64) -> Vec<C> {
    let mut state = seed;
    let mut u = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        ((state >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    };
    (0..n)
        .map(|_| C::from_polar((-2. * u().ln()).sqrt(), std::f64::consts::TAU * u()))
        .collect()
}
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().collect();
    ensure!(
        (a.len() == 3 || (a.len() == 4 && a[3] == "--spectrum-only"))
            && !std::path::Path::new(&a[1]).exists(),
        "new report.h5 correlation_lag_s [--spectrum-only] required"
    );
    let n = 4096;
    let sr = 1e6;
    let model = SpectrumModel {
        radar_frequency_hz: 440.2e6,
        ions: vec![IonSpecies {
            mass_amu: 16.,
            fraction: 1.,
        }],
        electron_collision_s: 0.,
    };
    let state = PlasmaState {
        electron_density_m3: 3e11,
        ion_temperature_k: 1000.,
        electron_temperature_k: 1800.,
        velocity_m_s: 0.,
        ion_collision_s: 0.,
    };
    let f: Vec<_> = (0..n)
        .map(|i| {
            if i < n / 2 {
                i as f64 * sr / n as f64
            } else {
                (i as f64 - n as f64) * sr / n as f64
            }
        })
        .collect();
    let mut psd = model.evaluate(state, &f)?;
    for (p, f) in psd.iter_mut().zip(&f) {
        if f.abs() > 20000. {
            *p = 0.;
        }
    }
    let norm = (2. * psd.iter().sum::<f64>()).sqrt();
    let fft = FftPlanner::<f64>::new().plan_fft_inverse(n);
    let settings = Settings {
        correlation_lag_s: a[2].parse()?,
        refine_delay_samples: 0.,
        refine_doppler_hz: 0.,
        ..Settings::default()
    };
    let report = hdf5::File::create(&a[1])?;
    report.new_attr_builder().with_data("Forced candidate; no coherent target; Gaussian ion-line PSD restricted to +/-20 kHz; no CFAR selection, range-extended plasma or real-archive false-alarm-rate claim").create("scope")?;
    report
        .new_dataset_builder()
        .with_data(&serde_json::to_vec(&(model, state))?)
        .create("physics_json")?;
    report
        .new_dataset_builder()
        .with_data(&f)
        .create("frequency_hz")?;
    report
        .new_dataset_builder()
        .with_data(&psd)
        .create("ion_line_psd")?;
    report
        .new_attr::<f64>()
        .create("sample_rate_hz")?
        .write_scalar(&sr)?;
    if a.get(3).is_some_and(|flag| flag == "--spectrum-only") {
        report
            .new_attr::<bool>()
            .create("spectrum_only")?
            .write_scalar(&true)?;
        return Ok(());
    }
    for family in ["long", "coded"] {
        let mut tx = vec![C::default(); n];
        for i in 64..544 {
            tx[i] = C::new(
                if family == "long" || ((i - 64) / 30) % 2 == 0 {
                    1.
                } else {
                    -1.
                },
                0.,
            );
        }
        let detection = Detection {
            delay_samples: 600.,
            doppler_hz: 0.,
        };
        let t = template(&tx, detection.delay_samples, 0., sr, n, 16);
        for (power_index, power) in [0., 0.03, 0.1, 0.3, 1.].iter().enumerate() {
            let mut rows = Vec::new();
            let mut applied = 0;
            let mut confirmed = 0;
            let trials = 100;
            for seed in 1..=trials {
                let mut colored = noise(n, seed + 10000);
                for (z, p) in colored.iter_mut().zip(&psd) {
                    *z *= p.sqrt() / norm;
                }
                fft.process(&mut colored);
                let bg = noise(n, seed);
                let y: Vec<_> = bg
                    .iter()
                    .zip(&colored)
                    .zip(&t)
                    .map(|((w, z), t)| *w / 2f64.sqrt() + *t * *z * (*power as f64).sqrt())
                    .collect();
                let c = cancel(&y, &tx, &[detection.clone()], sr, [3072, 4000], &settings)?;
                applied += usize::from(c.diagnostic.applied && c.diagnostic.accepted);
                confirmed += usize::from(c.diagnostic.confirmed);
                rows.push(c.diagnostic);
            }
            let g = report.create_group(&format!("{family}/power_{power_index}"))?;
            g.new_attr::<f64>()
                .create("diffuse_to_receiver_power")?
                .write_scalar(power)?;
            g.new_attr::<u64>()
                .create("trials")?
                .write_scalar(&trials)?;
            g.new_attr::<u64>()
                .create("accepted_subtractions")?
                .write_scalar(&(applied as u64))?;
            g.new_attr::<u64>()
                .create("confirmed_subtractions")?
                .write_scalar(&(confirmed as u64))?;
            g.new_dataset_builder()
                .with_data(&serde_json::to_vec(&rows)?)
                .create("diagnostics_json")?;
            println!("family={family} diffuse_power={power} trials={trials} false_applied={applied} false_confirmed={confirmed}");
        }
    }
    Ok(())
}
