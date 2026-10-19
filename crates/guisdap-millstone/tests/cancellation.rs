use guisdap_millstone::{
    debris::{cancel, template, Detection, Settings},
    C,
};
fn noise(n: usize, seed: u64) -> Vec<C> {
    let mut state = seed;
    let mut u = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        ((state >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    };
    (0..n)
        .map(|_| {
            let r = (-2. * u().ln()).sqrt();
            let t = std::f64::consts::TAU * u();
            C::from_polar(r, t)
        })
        .collect()
}
fn tx() -> Vec<C> {
    let mut tx = vec![C::default(); 2048];
    let code = [1., -1., 1., 1., -1., -1., 1., -1.];
    for (i, z) in tx.iter_mut().enumerate().take(320).skip(64) {
        *z = C::new(code[(i - 64) / 32], 0.);
    }
    tx
}
#[test]
fn injection_recovers_complex_amplitude_and_preserves_noise() {
    let tx = tx();
    let d = Detection {
        delay_samples: 500.25,
        doppler_hz: 1250.,
    };
    let t = template(&tx, d.delay_samples, d.doppler_hz, 1e6, tx.len(), 16);
    let bg = noise(tx.len(), 1);
    let amp = C::new(200., -100.);
    let y: Vec<_> = bg.iter().zip(&t).map(|(n, t)| n + amp * t).collect();
    let mut s = Settings::default();
    s.refine_delay_samples = 0.;
    s.refine_doppler_hz = 0.;
    let c = cancel(&y, &tx, &[d], 1e6, [1600, 2000], &s).unwrap();
    assert!(c.diagnostic.accepted, "{:?}", c.diagnostic);
    assert!((c.diagnostic.amplitudes[0][0] - amp.re).abs() < 0.4);
    assert!((c.diagnostic.amplitudes[0][1] - amp.im).abs() < 0.4);
    assert!(c.diagnostic.heldout_rms_over_noise > 0.8 && c.diagnostic.heldout_rms_over_noise < 1.2);
    assert_eq!(c.cleaned[1500], y[1500]);
}
#[test]
fn clipped_echo_is_never_accepted() {
    let tx = tx();
    let d = Detection {
        delay_samples: 500.,
        doppler_hz: 1000.,
    };
    let t = template(&tx, d.delay_samples, d.doppler_hz, 1e6, tx.len(), 16);
    let y: Vec<_> = noise(tx.len(), 2)
        .iter()
        .zip(&t)
        .map(|(n, t)| {
            let a = n + 50000. * t;
            C::new(a.re.clamp(-32768., 32767.), a.im.clamp(-32768., 32767.))
        })
        .collect();
    let c = cancel(&y, &tx, &[d], 1e6, [1600, 2000], &Settings::default()).unwrap();
    assert!(!c.diagnostic.accepted);
    assert_eq!(c.diagnostic.reason, "clipped_adc");
    assert!(c.usable.iter().any(|x| !*x));
    assert!(c.usable[1500]);
}
#[test]
fn wrong_doppler_or_bit_shape_fails_heldout_test() {
    let tx = tx();
    let d = Detection {
        delay_samples: 500.,
        doppler_hz: 1000.,
    };
    let t = template(&tx, d.delay_samples, 2500., 1e6, tx.len(), 16);
    let y: Vec<_> = noise(tx.len(), 3)
        .iter()
        .zip(&t)
        .map(|(n, t)| n + 100. * t)
        .collect();
    let c = cancel(&y, &tx, &[d], 1e6, [1600, 2000], &Settings::default()).unwrap();
    assert!(!c.diagnostic.accepted);
    assert!(c.diagnostic.clip_fraction == 0.);
}
#[test]
fn simultaneous_overlapping_targets_are_fitted_jointly() {
    let tx = tx();
    let ds = [
        Detection {
            delay_samples: 500.,
            doppler_hz: 0.,
        },
        Detection {
            delay_samples: 600.,
            doppler_hz: 3000.,
        },
    ];
    let mut y = noise(tx.len(), 9);
    for (d, a) in ds.iter().zip([C::new(100., 50.), C::new(50., -80.)]) {
        for (z, t) in y.iter_mut().zip(template(
            &tx,
            d.delay_samples,
            d.doppler_hz,
            1e6,
            tx.len(),
            16,
        )) {
            *z += a * t
        }
    }
    let mut s = Settings::default();
    s.refine_delay_samples = 0.;
    s.refine_doppler_hz = 0.;
    let c = cancel(&y, &tx, &ds, 1e6, [1600, 2000], &s).unwrap();
    assert!(c.diagnostic.accepted, "{:?}", c.diagnostic);
    assert_eq!(c.diagnostic.fitted_complex_parameters, 2);
}

#[test]
fn weak_subtraction_is_distinct_from_confirmed_removal_and_extremely_weak_is_skipped() {
    let tx = tx();
    let d = Detection {
        delay_samples: 500.,
        doppler_hz: 0.,
    };
    let t = template(&tx, d.delay_samples, d.doppler_hz, 1e6, tx.len(), 16);
    let mut s = Settings::default();
    s.refine_delay_samples = 0.;
    s.refine_doppler_hz = 0.;
    // About 9 dB matched-filter amplitude SNR, despite per-sample power far below noise.
    let y: Vec<_> = noise(tx.len(), 123)
        .iter()
        .zip(&t)
        .map(|(n, t)| n + C::new(0.36, 0.12) * t)
        .collect();
    let c = cancel(&y, &tx, &[d.clone()], 1e6, [1600, 2000], &s).unwrap();
    assert!(c.diagnostic.applied, "{:?}", c.diagnostic);
    assert!(!c.diagnostic.confirmed);
    assert!(c.diagnostic.accepted);
    assert_eq!(c.diagnostic.reason, "weak_applied");
    let y = noise(tx.len(), 123);
    let c = cancel(&y, &tx, &[d], 1e6, [1600, 2000], &s).unwrap();
    assert!(!c.diagnostic.applied);
    assert!(!c.diagnostic.confirmed);
    assert_eq!(c.diagnostic.reason, "skipped_low_snr");
    assert_eq!(c.cleaned, y);
    assert!(c.usable.iter().all(|v| *v));
}
#[test]
fn echo_over_background_rejects_pulse_without_fabricated_calibration() {
    let tx = tx();
    let d = Detection {
        delay_samples: 1500.,
        doppler_hz: 1000.,
    };
    let t = template(&tx, d.delay_samples, d.doppler_hz, 1e6, tx.len(), 16);
    let y: Vec<_> = noise(tx.len(), 71)
        .iter()
        .zip(&t)
        .map(|(n, t)| n + 200. * t)
        .collect();
    let c = cancel(&y, &tx, &[d], 1e6, [1600, 2000], &Settings::default()).unwrap();
    assert_eq!(c.cleaned, y);
    assert!(c.usable.iter().all(|&v| !v));
    assert!(!c.diagnostic.accepted);
    assert!(!c.diagnostic.applied);
    assert_eq!(c.diagnostic.reason, "background_gate_contaminated");
    let qa = serde_json::to_value(c.diagnostic).unwrap();
    assert!(qa["noise_rms"].is_null());
    assert!(qa["heldout_rms_over_noise"].is_null());
}

#[test]
fn dead_recorded_background_rejects_only_that_pulse() {
    let tx = tx();
    let y = vec![C::new(0.212, 0.221); tx.len()];
    let c = cancel(&y, &tx, &[], 1e6, [1600, 2000], &Settings::default()).unwrap();
    assert_eq!(c.cleaned, y);
    assert!(c.usable.iter().all(|v| !*v));
    assert!(!c.diagnostic.accepted);
    assert!(!c.diagnostic.applied);
    assert_eq!(c.diagnostic.reason, "invalid_background_noise");
    assert!(serde_json::to_value(&c.diagnostic).unwrap()["noise_rms"].is_null());
    let good = noise(tx.len(), 44);
    let next = cancel(&good, &tx, &[], 1e6, [1600, 2000], &Settings::default()).unwrap();
    assert_eq!(next.cleaned, good);
    assert!(next.usable.iter().all(|v| *v));
}
#[test]
fn legacy_settings_keep_identical_serialized_source_configuration() {
    let value = serde_json::to_value(Settings::default()).unwrap();
    assert!(value.get("correlation_lag_s").is_none());
}
