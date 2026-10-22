//! Direct local Digital RF 2.x reads from multiple top-level directories.
//! Absolute sample indices stay u64; file-time arithmetic uses checked u128.
//! Missing files and indexed gaps are returned explicitly, never zero-filled.
//! No RF conversion, pack files, or persistent index is created.
use anyhow::{bail, ensure, Context, Result};
use chrono::{DateTime, Utc};
use hdf5::{Dataset, File, H5Type};
use num_complex::Complex64;
use std::{
    collections::{BTreeMap, VecDeque},
    fs,
    ops::Range,
    path::{Path, PathBuf},
};

#[repr(C)]
#[derive(Clone, Copy, Debug)]
struct ComplexValue {
    r: f64,
    i: f64,
}
unsafe impl H5Type for ComplexValue {
    fn type_descriptor() -> hdf5::types::TypeDescriptor {
        use hdf5::types::{CompoundField, CompoundType, TypeDescriptor};
        TypeDescriptor::Compound(CompoundType {
            fields: vec![
                CompoundField::typed::<f64>("r", 0, 0),
                CompoundField::typed::<f64>("i", 8, 1),
            ],
            size: 16,
        })
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Properties {
    pub sample_rate_numerator: u64,
    pub sample_rate_denominator: u64,
    pub num_subchannels: usize,
    pub is_complex: bool,
    pub epoch: String,
    pub format: Format,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Real32,
    Real64,
    Complex32,
    Complex64,
}
#[derive(Debug)]
pub enum Samples {
    Real32(Vec<f32>),
    Real64(Vec<f64>),
    Complex32(Vec<num_complex::Complex32>),
    Complex64(Vec<Complex64>),
}
#[derive(Debug)]
pub struct Vector {
    pub samples: Samples,
    pub shape: Vec<usize>,
}
impl Vector {
    pub fn into_complex64(self) -> Vec<Complex64> {
        match self.samples {
            Samples::Real32(v) => v
                .into_iter()
                .map(|x| Complex64::new(f64::from(x), 0.))
                .collect(),
            Samples::Real64(v) => v.into_iter().map(|x| Complex64::new(x, 0.)).collect(),
            Samples::Complex32(v) => v
                .into_iter()
                .map(|z| Complex64::new(f64::from(z.re), f64::from(z.im)))
                .collect(),
            Samples::Complex64(v) => v,
        }
    }
}
#[derive(Clone, Debug)]
pub struct CacheOptions {
    pub minimum_bytes: usize,
    pub max_open_files: usize,
}
impl Default for CacheOptions {
    fn default() -> Self {
        Self {
            minimum_bytes: 64 * 1024 * 1024,
            max_open_files: 8,
        }
    }
}
#[derive(Clone, Copy, Debug, Default)]
pub struct CacheStats {
    pub file_opens: u64,
    pub file_hits: u64,
    pub last_chunk_cache_bytes: usize,
}
#[derive(Clone)]
struct Location {
    path: PathBuf,
    file_ms: u64,
    subdir_s: u64,
}
struct Channel {
    properties: Properties,
    locations: Vec<Location>,
}
struct Cached {
    path: PathBuf,
    _file: File,
    data: Dataset,
    index: Vec<(u64, u64)>,
    rows: u64,
}
#[derive(Clone, Debug)]
pub struct SourceSpan {
    pub samples: Range<u64>,
    pub path: PathBuf,
}
#[derive(Clone, Debug)]
pub struct SampleBlock {
    pub start: u64,
    pub data: Vec<Complex64>,
}
#[derive(Debug)]
pub struct MissingSamples {
    pub gaps: Vec<Range<u64>>,
}
impl std::fmt::Display for MissingSamples {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "DigitalRF sample gaps: {:?}", self.gaps)
    }
}
impl std::error::Error for MissingSamples {}
#[derive(Clone, Debug)]
pub struct Read {
    pub requested: Range<u64>,
    pub blocks: Vec<SampleBlock>,
    pub gaps: Vec<Range<u64>>,
    pub sources: Vec<SourceSpan>,
}
impl Read {
    pub fn is_complete(&self) -> bool {
        self.gaps.is_empty()
    }
    pub fn into_contiguous(self) -> Result<Vec<Complex64>> {
        if !self.gaps.is_empty() {
            return Err(MissingSamples { gaps: self.gaps }.into());
        }
        if self.requested.is_empty() {
            return Ok(Vec::new());
        }
        ensure!(
            self.blocks.len() == 1 && self.blocks[0].start == self.requested.start,
            "invalid contiguous read"
        );
        Ok(self.blocks.into_iter().next().unwrap().data)
    }
}
/// Read-only DigitalRF reader. Supply recording roots containing channels,
/// recording roots containing `rf_data/`, or individual channel directories.
/// Duplicate roots are deduplicated. Overlapping samples must agree exactly;
/// inconsistent sample rates, epochs, or subchannel layouts are errors.
pub struct Reader {
    channels: BTreeMap<String, Channel>,
    cache: VecDeque<Cached>,
    cache_files: usize,
    minimum_bytes: usize,
    stats: CacheStats,
}
fn attr_u64(h: &File, n: &str) -> Result<u64> {
    let v = h.attr(n)?.read_raw::<u64>()?;
    ensure!(v.len() == 1, "non-scalar {n}");
    Ok(v[0])
}
fn attr_text(h: &File, n: &str) -> Result<String> {
    let a = h.attr(n)?;
    match a.dtype()?.to_descriptor()? {
        hdf5::types::TypeDescriptor::FixedAscii(_) => Ok(a
            .read_raw::<hdf5::types::FixedAscii<256>>()?[0]
            .as_str()
            .to_owned()),
        hdf5::types::TypeDescriptor::VarLenUnicode => Ok(a
            .read_raw::<hdf5::types::VarLenUnicode>()?[0]
            .as_str()
            .to_owned()),
        hdf5::types::TypeDescriptor::VarLenAscii => Ok(a.read_raw::<hdf5::types::VarLenAscii>()?
            [0]
        .as_str()
        .to_owned()),
        _ => bail!("unsupported string attribute {n}"),
    }
}
fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b > 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}
impl Reader {
    pub fn open(roots: impl IntoIterator<Item = impl AsRef<Path>>) -> Result<Self> {
        let mut channels: BTreeMap<String, Channel> = BTreeMap::new();
        let mut seen = std::collections::BTreeSet::new();
        for root in roots {
            let root = root.as_ref();
            ensure!(
                root.is_dir(),
                "DigitalRF root does not exist: {}",
                root.display()
            );
            let root = if root.join("rf_data").is_dir() {
                root.join("rf_data")
            } else {
                root.to_path_buf()
            };
            let candidates = if root.join("drf_properties.h5").is_file() {
                vec![root]
            } else {
                fs::read_dir(root)?
                    .filter_map(|e| e.ok().map(|x| x.path()))
                    .filter(|p| p.join("drf_properties.h5").is_file())
                    .collect()
            };
            for path in candidates {
                let path = path.canonicalize()?;
                if !seen.insert(path.clone()) {
                    continue;
                }
                let h = File::open(path.join("drf_properties.h5"))?;
                let version = attr_text(&h, "digital_rf_version")?;
                ensure!(
                    version.split('.').next() == Some("2"),
                    "unsupported DigitalRF version {version}"
                );
                let n = attr_u64(&h, "sample_rate_numerator")?;
                let d = attr_u64(&h, "sample_rate_denominator")?;
                ensure!(n > 0 && d > 0, "invalid sample rate");
                let divisor = gcd(n, d);
                let class = attr_u64(&h, "H5Tget_class")?;
                let bytes = attr_u64(&h, "H5Tget_size")?;
                ensure!(
                    matches!(class, 0 | 1) && matches!(bytes, 1 | 2 | 4 | 8),
                    "unsupported DigitalRF numeric type"
                );
                let complex = attr_u64(&h, "is_complex")? != 0;
                let small = if class == 0 { bytes <= 2 } else { bytes <= 4 };
                let format = match (complex, small) {
                    (false, true) => Format::Real32,
                    (false, false) => Format::Real64,
                    (true, true) => Format::Complex32,
                    (true, false) => Format::Complex64,
                };
                let properties = Properties {
                    format,
                    sample_rate_numerator: n / divisor,
                    sample_rate_denominator: d / divisor,
                    num_subchannels: usize::try_from(attr_u64(&h, "num_subchannels")?)?,
                    is_complex: attr_u64(&h, "is_complex")? != 0,
                    epoch: attr_text(&h, "epoch")?,
                };
                ensure!(
                    properties.num_subchannels > 0 && properties.epoch == "1970-01-01T00:00:00Z",
                    "unsupported subchannel count or epoch"
                );
                let location = Location {
                    path: path.clone(),
                    file_ms: attr_u64(&h, "file_cadence_millisecs")?,
                    subdir_s: attr_u64(&h, "subdir_cadence_secs")?,
                };
                ensure!(
                    location.file_ms > 0 && location.subdir_s > 0,
                    "invalid file cadence"
                );
                let name = path
                    .file_name()
                    .unwrap()
                    .to_str()
                    .context("non-UTF8 channel")?
                    .to_owned();
                let channel = channels.entry(name.clone()).or_insert_with(|| Channel {
                    properties: properties.clone(),
                    locations: Vec::new(),
                });
                ensure!(
                    channel.properties == properties,
                    "incompatible DigitalRF properties for channel {name}"
                );
                channel.locations.push(location);
            }
        }
        ensure!(!channels.is_empty(), "no DigitalRF channels found");
        Ok(Self {
            channels,
            cache: VecDeque::new(),
            cache_files: 8,
            minimum_bytes: 64 * 1024 * 1024,
            stats: CacheStats::default(),
        })
    }
    pub fn set_cache_options(&mut self, options: CacheOptions) -> Result<()> {
        ensure!(
            options.minimum_bytes > 0 && options.max_open_files > 0,
            "nonpositive RF cache budget"
        );
        self.clear_cache();
        self.cache_files = options.max_open_files;
        self.minimum_bytes = options.minimum_bytes;
        Ok(())
    }
    pub fn cache_stats(&self) -> CacheStats {
        self.stats
    }
    /// DigitalRF read_vector(start_sample, vector_length, channel, sub_channel).
    /// Contiguous reads are strict: a missing sample is an error for this vector,
    /// which the campaign handles by auditing/skipping only the affected pulse.
    /// All subchannels are returned row-major when None; single dimensions squeeze.
    /// Numeric promotion follows NumPy promote_types(f4/c8, stored component).
    pub fn read_vector(
        &mut self,
        start: u64,
        length: usize,
        channel: &str,
        subchannel: Option<usize>,
    ) -> Result<Vector> {
        ensure!(
            length > 0,
            "Number of samples requested must be greater than 0"
        );
        let end = start
            .checked_add(u64::try_from(length)?)
            .context("vector sample overflow")?;
        let p = self.properties(channel)?.clone();
        let channels: Vec<_> =
            subchannel.map_or_else(|| (0..p.num_subchannels).collect(), |s| vec![s]);
        ensure!(
            channels.iter().all(|&c| c < p.num_subchannels),
            "invalid subchannel"
        );
        // Match the upstream squeeze/length check, including the single-row,
        // multiple-subchannel case rejected by its read_vector_raw implementation.
        ensure!(
            length != 1 || channels.len() == 1,
            "Requested 1 samples, but got {} after squeeze",
            channels.len()
        );
        let columns: Vec<_> = channels
            .iter()
            .map(|&c| self.read(channel, start..end, c)?.into_contiguous())
            .collect::<Result<_>>()?;
        let size = length
            .checked_mul(columns.len())
            .context("vector allocation overflow")?;
        let mut values = Vec::with_capacity(size);
        for row in 0..length {
            for column in &columns {
                values.push(column[row]);
            }
        }
        let shape = if channels.len() == 1 {
            vec![length]
        } else {
            vec![length, channels.len()]
        };
        let samples = match p.format {
            Format::Real32 => Samples::Real32(values.iter().map(|z| z.re as f32).collect()),
            Format::Real64 => Samples::Real64(values.iter().map(|z| z.re).collect()),
            Format::Complex32 => Samples::Complex32(
                values
                    .iter()
                    .map(|z| num_complex::Complex32::new(z.re as f32, z.im as f32))
                    .collect(),
            ),
            Format::Complex64 => Samples::Complex64(values),
        };
        Ok(Vector { samples, shape })
    }
    pub fn channels(&self) -> Vec<&str> {
        self.channels.keys().map(String::as_str).collect()
    }
    pub fn properties(&self, channel: &str) -> Result<&Properties> {
        Ok(&self
            .channels
            .get(channel)
            .context("unknown DigitalRF channel")?
            .properties)
    }
    pub fn clear_cache(&mut self) {
        self.cache.clear();
    }
    fn file(&mut self, path: &Path, subchannels: usize) -> Result<Option<&Cached>> {
        if let Some(i) = self.cache.iter().position(|c| c.path == path) {
            let c = self.cache.remove(i).unwrap();
            self.cache.push_back(c);
            self.stats.file_hits += 1;
            return Ok(self.cache.back());
        }
        if !path.is_file() {
            return Ok(None);
        }
        let h = match File::open(path) {
            Ok(h) => h,
            Err(_) if !path.exists() => return Ok(None),
            Err(e) => return Err(e).with_context(|| format!("opening {}", path.display())),
        };
        let probe = h.dataset("rf_data")?;
        let logical = probe
            .shape()
            .iter()
            .try_fold(probe.dtype()?.size(), |a, &b| {
                a.checked_mul(b).context("RF dataset size overflow")
            })?;
        let bytes = self.minimum_bytes.max(logical);
        let chunks = probe
            .chunk()
            .map(|chunk| chunk.iter().product::<usize>())
            .unwrap_or(1)
            .max(1);
        let chunk_count = probe.size().div_ceil(chunks);
        let slots = 1_000_003usize.max(chunk_count.saturating_mul(10).saturating_add(1));
        drop(probe);
        drop(h);
        let h = File::with_options()
            .with_fapl(|p| p.chunk_cache(slots, bytes, 0.75))
            .open(path)?;
        self.stats.file_opens += 1;
        self.stats.last_chunk_cache_bytes = bytes;
        let data = h.dataset("rf_data")?;
        let shape = data.shape();
        ensure!(
            shape.len() == 2 && shape[1] == subchannels,
            "invalid rf_data shape in {}",
            path.display()
        );
        let idx = h.dataset("rf_data_index")?;
        ensure!(
            idx.shape().len() == 2 && idx.shape()[1] == 2,
            "invalid rf_data_index"
        );
        let index: Vec<_> = idx
            .read_raw::<u64>()?
            .chunks_exact(2)
            .map(|x| (x[0], x[1]))
            .collect();
        let rows = shape[0] as u64;
        ensure!(
            !index.is_empty() && index[0].1 == 0,
            "missing first RF block"
        );
        for (i, &(sample, offset)) in index.iter().enumerate() {
            let end = index.get(i + 1).map_or(rows, |x| x.1);
            ensure!(offset < end && end <= rows, "invalid block offsets");
            let stop = sample
                .checked_add(end - offset)
                .context("RF sample index overflow")?;
            if let Some(next) = index.get(i + 1) {
                ensure!(stop <= next.0, "overlapping rf_data_index blocks");
            }
        }
        self.cache.push_back(Cached {
            path: path.to_owned(),
            _file: h,
            data,
            index,
            rows,
        });
        while self.cache.len() > self.cache_files {
            self.cache.pop_front();
        }
        Ok(self.cache.back())
    }
    /// Half-open global sample range. File cadence guides lookup; rf_data_index
    /// is authoritative for block boundaries and gaps inside each HDF5 file.
    pub fn read(&mut self, channel: &str, samples: Range<u64>, subchannel: usize) -> Result<Read> {
        ensure!(samples.start <= samples.end, "reversed sample interval");
        let c = self
            .channels
            .get(channel)
            .context("unknown DigitalRF channel")?;
        let p = c.properties.clone();
        let locations = c.locations.clone();
        ensure!(subchannel < p.num_subchannels, "invalid subchannel");
        let mut pieces: Vec<SampleBlock> = Vec::new();
        let mut sources = Vec::new();
        if !samples.is_empty() {
            let millis = |s: u64| -> Result<u64> {
                Ok(u64::try_from(
                    u128::from(s) * u128::from(p.sample_rate_denominator) * 1000
                        / u128::from(p.sample_rate_numerator),
                )?)
            };
            for l in locations {
                let mut ms = millis(samples.start)? / l.file_ms * l.file_ms;
                let last = millis(samples.end - 1)? / l.file_ms * l.file_ms;
                loop {
                    let second = ms / 1000;
                    let directory = second / l.subdir_s * l.subdir_s;
                    let date = DateTime::<Utc>::from_timestamp(i64::try_from(directory)?, 0)
                        .context("timestamp outside calendar range")?;
                    let path = l
                        .path
                        .join(date.format("%Y-%m-%dT%H-%M-%S").to_string())
                        .join(format!("rf@{}.{:03}.h5", ms / 1000, ms % 1000));
                    if let Some(f) = self.file(&path, p.num_subchannels)? {
                        for (i, &(start, offset)) in f.index.iter().enumerate() {
                            let stop =
                                start + (f.index.get(i + 1).map_or(f.rows, |x| x.1) - offset);
                            let lo = start.max(samples.start);
                            let hi = stop.min(samples.end);
                            if lo >= hi {
                                continue;
                            }
                            let a = usize::try_from(offset + lo - start)?;
                            let b = usize::try_from(offset + hi - start)?;
                            let data = if p.is_complex {
                                f.data
                                    .read_slice_1d::<ComplexValue, _>((a..b, subchannel))?
                                    .iter()
                                    .map(|z| Complex64::new(z.r, z.i))
                                    .collect()
                            } else {
                                f.data
                                    .read_slice_1d::<f64, _>((a..b, subchannel))?
                                    .iter()
                                    .map(|x| Complex64::new(*x, 0.))
                                    .collect()
                            };
                            pieces.push(SampleBlock { start: lo, data });
                            sources.push(SourceSpan {
                                samples: lo..hi,
                                path: path.clone(),
                            });
                        }
                    }
                    if ms == last {
                        break;
                    }
                    ms = ms
                        .checked_add(l.file_ms)
                        .context("file timestamp overflow")?;
                    ensure!(ms <= last, "invalid file cadence traversal");
                }
            }
        }
        pieces.sort_by_key(|b| b.start);
        let mut blocks: Vec<SampleBlock> = Vec::new();
        for piece in pieces {
            if let Some(last) = blocks.last_mut() {
                let end = last.start + last.data.len() as u64;
                if piece.start <= end {
                    let overlap = (end - piece.start).min(piece.data.len() as u64) as usize;
                    let position = usize::try_from(piece.start - last.start)?;
                    ensure!(
                        last.data[position..position + overlap]
                            .iter()
                            .zip(&piece.data[..overlap])
                            .all(|(a, b)| a.re.to_bits() == b.re.to_bits()
                                && a.im.to_bits() == b.im.to_bits()),
                        "conflicting DigitalRF overlap at sample {}",
                        piece.start
                    );
                    last.data.extend_from_slice(&piece.data[overlap..]);
                    continue;
                }
            }
            blocks.push(piece);
        }
        let mut gaps = Vec::new();
        let mut next = samples.start;
        for b in &blocks {
            if b.start > next {
                gaps.push(next..b.start);
            }
            next = b.start + b.data.len() as u64;
        }
        if next < samples.end {
            gaps.push(next..samples.end);
        }
        Ok(Read {
            requested: samples,
            blocks,
            gaps,
            sources,
        })
    }
    /// Inclusive first/last stored sample across all directories. Only file
    /// indices are inspected, never RF voltage datasets. Gaps remain possible.
    pub fn bounds(&mut self, channel: &str) -> Result<Option<(u64, u64)>> {
        let c = self
            .channels
            .get(channel)
            .context("unknown DigitalRF channel")?;
        let locations = c.locations.clone();
        let ns = c.properties.num_subchannels;
        let mut first = None;
        let mut last = None;
        for l in locations {
            let mut files = Vec::new();
            for d in fs::read_dir(&l.path)? {
                let d = d?.path();
                if !d.is_dir() {
                    continue;
                }
                for f in fs::read_dir(d)? {
                    let p = f?.path();
                    let n = p.file_name().unwrap().to_string_lossy();
                    if let Some(t) = n.strip_prefix("rf@").and_then(|s| s.strip_suffix(".h5")) {
                        if let Some((s, ms)) = t.split_once('.') {
                            if let (Ok(s), Ok(ms)) = (s.parse::<u64>(), ms.parse::<u64>()) {
                                files.push((u128::from(s) * 1000 + u128::from(ms), p));
                            }
                        }
                    }
                }
            }
            files.sort_by_key(|f| f.0);
            for (_, path) in &files {
                if let Some(f) = self.file(path, ns)? {
                    let a = f.index[0].0;
                    first = Some(first.map_or(a, |v: u64| v.min(a)));
                    break;
                }
            }
            for (_, path) in files.iter().rev() {
                if let Some(f) = self.file(path, ns)? {
                    let &(s, o) = f.index.last().unwrap();
                    let b = s + f.rows - o - 1;
                    last = Some(last.map_or(b, |v: u64| v.max(b)));
                    break;
                }
            }
        }
        Ok(first.zip(last))
    }
}
