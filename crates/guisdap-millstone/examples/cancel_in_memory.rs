//! Diagnostic IPC helper: voltages enter and leave anonymous pipes only.
//! Each frame is u64 JSON-byte length, JSON header, then little-endian complex128
//! receive and transmit arrays. Responses contain diagnostic JSON, cleaned
//! complex128 and a byte mask. This is not a stored RF format or production API.
use anyhow::{ensure, Result};
use guisdap_millstone::{
    debris::{self, Detection, Settings},
    C,
};
use serde::Deserialize;
use std::io::{Read, Write};

#[derive(Deserialize)]
struct Request {
    n: usize,
    sample_rate_hz: f64,
    noise_gate: [usize; 2],
    detections: Vec<Detection>,
    #[serde(default)]
    settings: Settings,
}
fn read_complex(input: &mut impl Read, n: usize) -> Result<Vec<C>> {
    let mut raw = vec![0u8; n * 16];
    input.read_exact(&mut raw)?;
    Ok(raw
        .chunks_exact(16)
        .map(|v| {
            C::new(
                f64::from_le_bytes(v[..8].try_into().unwrap()),
                f64::from_le_bytes(v[8..].try_into().unwrap()),
            )
        })
        .collect())
}
fn main() -> Result<()> {
    let mut input = std::io::stdin().lock();
    let mut output = std::io::BufWriter::new(std::io::stdout().lock());
    loop {
        let mut size = [0u8; 8];
        // Clean EOF between frames is the only successful shutdown.
        match input.read(&mut size[..1])? {
            0 => return Ok(()),
            _ => input.read_exact(&mut size[1..])?,
        };
        let size = u64::from_le_bytes(size) as usize;
        ensure!(
            size > 0 && size <= 1_000_000,
            "Invalid diagnostic frame size"
        );
        let mut json = vec![0u8; size];
        input.read_exact(&mut json)?;
        let r: Request = serde_json::from_slice(&json)?;
        ensure!(r.n > 0 && r.n <= 100_000, "Invalid diagnostic sample count");
        let y = read_complex(&mut input, r.n)?;
        let tx = read_complex(&mut input, r.n)?;
        let c = match debris::cancel(
            &y,
            &tx,
            &r.detections,
            r.sample_rate_hz,
            r.noise_gate,
            &r.settings,
        ) {
            Ok(c) => c,
            Err(e) => {
                let json = serde_json::to_vec(&serde_json::json!({"error":e.to_string()}))?;
                output.write_all(&(json.len() as u64).to_le_bytes())?;
                output.write_all(&json)?;
                output.flush()?;
                continue;
            }
        };
        let json = serde_json::to_vec(&c.diagnostic)?;
        output.write_all(&(json.len() as u64).to_le_bytes())?;
        output.write_all(&json)?;
        for z in c.cleaned {
            output.write_all(&z.re.to_le_bytes())?;
            output.write_all(&z.im.to_le_bytes())?;
        }
        for value in c.usable {
            output.write_all(&[u8::from(value)])?;
        }
        output.flush()?;
    }
}
