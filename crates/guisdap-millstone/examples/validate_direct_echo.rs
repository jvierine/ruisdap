//! Bounded original-archive cancellation check: one catalogue echo, no RF output.
use anyhow::{ensure, Result};
use guisdap_millstone::{
    archive::{Archive, Manifest},
    debris, C,
};
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().collect();
    let m: Manifest = serde_json::from_slice(&std::fs::read(&a[1])?)?;
    let k: u64 = a[2].parse()?;
    let mut archive = Archive::open(&m)?;
    let (observations, _, _) = archive.observations(m.start_unix_s)?;
    for p in observations {
        if let Some(i) = p.samples.iter().position(|&s| s == k) {
            let dc = C::new(-0.212, -0.221);
            let echo: Vec<_> = p.echo[i].iter().map(|z| *z - dc).collect();
            let mut tx = p.transmit[i].clone();
            for z in &mut tx[p.config.lpi.tx_gate[0]..p.config.lpi.tx_gate[1]] {
                *z -= dc;
            }
            let result = debris::cancel(
                &echo,
                &tx,
                &p.detections[i],
                p.config.lpi.sample_rate_hz,
                p.config.lpi.noise_gate,
                &p.config.cancellation,
            )?;
            ensure!(
                result.diagnostic.accepted && result.diagnostic.confirmed,
                "known real echo failed"
            );
            println!("{}", serde_json::to_string_pretty(&result.diagnostic)?);
            return Ok(());
        }
    }
    anyhow::bail!("reference echo not selected")
}
