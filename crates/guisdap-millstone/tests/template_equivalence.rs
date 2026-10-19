use guisdap_millstone::{debris::template, C};
fn reference(tx: &[C], delay: f64, doppler: f64, sr: f64, n: usize, half: usize) -> Vec<C> {
    fn sinc(x: f64) -> f64 {
        if x.abs() < 1e-12 {
            1.
        } else {
            (std::f64::consts::PI * x).sin() / (std::f64::consts::PI * x)
        }
    }
    (0..n)
        .map(|i| {
            let t = i as f64 - delay;
            let centre = t.floor() as isize;
            let mut z = C::default();
            for j in centre - half as isize + 1..=centre + half as isize {
                if j >= 0 && (j as usize) < tx.len() {
                    let u = t - j as f64;
                    if u.abs() < half as f64 {
                        z += tx[j as usize] * (sinc(u) * sinc(u / half as f64));
                    }
                }
            }
            z * C::from_polar(1., 2. * std::f64::consts::PI * doppler * i as f64 / sr)
        })
        .collect()
}
#[test]
fn optimized_templates_preserve_every_output_bit() {
    for sparse in [false, true] {
        let mut tx = vec![C::default(); 1024];
        if sparse {
            for (i, z) in tx.iter_mut().enumerate().take(320).skip(64) {
                *z = C::new(
                    if i % 3 == 0 { -1. } else { 1. },
                    if i % 7 == 0 { 0.25 } else { 0. },
                );
            }
        }
        for delay in [-40.25, 0., 0.5, 500.25, 900.75, 1200.] {
            for doppler in [-35000., 0., 1250.] {
                let a = reference(&tx, delay, doppler, 1e6, 1024, 16);
                let b = template(&tx, delay, doppler, 1e6, 1024, 16);
                for (i, (a, b)) in a.iter().zip(&b).enumerate() {
                    assert_eq!(
                        (a.re.to_bits(), a.im.to_bits()),
                        (b.re.to_bits(), b.im.to_bits()),
                        "sample {i}, delay {delay}, Doppler {doppler}"
                    );
                }
            }
        }
    }
}
