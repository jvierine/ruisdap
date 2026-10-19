//! Direct-reader validation boundary. Writes only the requested test output,
//! not an intermediate archive; production workers use Reader in memory.
use anyhow::{ensure, Result};
use guisdap_digital_rf::Reader;
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().collect();
    ensure!(a.len() > 5, "output.h5 channel start length root...");
    let mut r = Reader::open(&a[5..])?;
    let start = a[3].parse()?;
    let length = a[4].parse()?;
    let out = r.read_vector(start, length, &a[2], None)?;
    let shape = out.shape.clone();
    let z = out.into_complex64();
    let h = hdf5::File::create(&a[1])?;
    let data: Vec<_> = z.iter().flat_map(|z| [z.re, z.im]).collect();
    let d = h
        .new_dataset::<f64>()
        .shape([z.len(), 2])
        .create("vector")?;
    d.write_raw(&data)?;
    h.new_dataset_builder().with_data(&shape).create("shape")?;
    h.new_attr::<u64>()
        .create("cache_bytes")?
        .write_scalar(&(r.cache_stats().last_chunk_cache_bytes as u64))?;
    Ok(())
}
