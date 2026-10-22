use guisdap_digital_rf::{CacheOptions, Reader, Samples};
use hdf5::{
    types::{CompoundField, CompoundType, FixedAscii, TypeDescriptor},
    File, H5Type,
};
use std::path::{Path, PathBuf};
#[repr(C)]
#[derive(Clone, Copy)]
struct IQ {
    r: i16,
    i: i16,
}
unsafe impl H5Type for IQ {
    fn type_descriptor() -> TypeDescriptor {
        TypeDescriptor::Compound(CompoundType {
            fields: vec![
                CompoundField::typed::<i16>("r", 0, 0),
                CompoundField::typed::<i16>("i", 2, 1),
            ],
            size: 4,
        })
    }
}
fn properties(root: &Path, n: u64, d: u64, channels: usize) -> PathBuf {
    let p = root.join("ch");
    std::fs::create_dir_all(&p).unwrap();
    let h = File::create(p.join("drf_properties.h5")).unwrap();
    for (k, v) in [
        ("sample_rate_numerator", n),
        ("sample_rate_denominator", d),
        ("num_subchannels", channels as u64),
        ("is_complex", 1),
        ("file_cadence_millisecs", 1000),
        ("subdir_cadence_secs", 3600),
        ("H5Tget_class", 0),
        ("H5Tget_size", 2),
    ] {
        h.new_attr::<u64>()
            .create(k)
            .unwrap()
            .write_scalar(&v)
            .unwrap();
    }
    h.new_attr::<FixedAscii<5>>()
        .create("digital_rf_version")
        .unwrap()
        .write_scalar(&FixedAscii::<5>::from_ascii("2.6.0").unwrap())
        .unwrap();
    h.new_attr::<FixedAscii<20>>()
        .create("epoch")
        .unwrap()
        .write_scalar(&FixedAscii::<20>::from_ascii("1970-01-01T00:00:00Z").unwrap())
        .unwrap();
    p
}
fn file(channel: &Path, ms: u64, index: &[[u64; 2]], logical: &[i16], ns: usize) {
    let dir = channel.join(
        chrono::DateTime::from_timestamp((ms / 1000 / 3600 * 3600) as i64, 0)
            .unwrap()
            .format("%Y-%m-%dT%H-%M-%S")
            .to_string(),
    );
    std::fs::create_dir_all(&dir).unwrap();
    let h = File::create(dir.join(format!("rf@{}.{:03}.h5", ms / 1000, ms % 1000))).unwrap();
    let data: Vec<_> = logical
        .iter()
        .flat_map(|&v| {
            (0..ns).map(move |c| IQ {
                r: v + c as i16 * 1000,
                i: -v,
            })
        })
        .collect();
    h.new_dataset::<IQ>()
        .shape([logical.len(), ns])
        .chunk([logical.len().min(100).max(1), ns])
        .deflate(1)
        .create("rf_data")
        .unwrap()
        .write_raw(&data)
        .unwrap();
    h.new_dataset::<u64>()
        .shape([index.len(), 2])
        .create("rf_data_index")
        .unwrap()
        .write_raw(&index.iter().flatten().copied().collect::<Vec<_>>())
        .unwrap();
}
#[test]
fn crosses_files_roots_and_hour_directory_without_packing() {
    let t = tempfile::tempdir().unwrap();
    let a = t.path().join("a");
    let b = t.path().join("b");
    let ac = properties(&a, 10, 1, 1);
    let bc = properties(&b, 10, 1, 1);
    file(&ac, 3599000, &[[35990, 0]], &(0..10).collect::<Vec<_>>(), 1);
    file(
        &bc,
        3600000,
        &[[36000, 0]],
        &(10..20).collect::<Vec<_>>(),
        1,
    );
    let mut r = Reader::open([&a, &b]).unwrap();
    assert_eq!(r.channels(), vec!["ch"]);
    assert_eq!(r.bounds("ch").unwrap(), Some((35990, 36009)));
    let v = r.read_vector(35995, 10, "ch", None).unwrap();
    assert_eq!(v.shape, vec![10]);
    assert!(matches!(v.samples, Samples::Complex32(_)));
    let z = v.into_complex64();
    assert_eq!(
        z.iter().map(|v| v.re as i16).collect::<Vec<_>>(),
        (5..15).collect::<Vec<_>>()
    );
}
#[test]
fn indexed_gaps_remain_gaps_and_subchannels_match() {
    let t = tempfile::tempdir().unwrap();
    let a = t.path().join("a");
    let b = t.path().join("b");
    let ac = properties(&a, 10, 1, 2);
    let bc = properties(&b, 10, 1, 2);
    file(&ac, 0, &[[0, 0], [7, 4]], &[0, 1, 2, 3, 7, 8, 9], 2);
    file(&bc, 1000, &[[10, 0]], &[10, 11, 12, 13, 14], 2);
    let mut r = Reader::open([a, b]).unwrap();
    let out = r.read("ch", 2..15, 1).unwrap();
    assert_eq!(out.gaps, vec![4..7]);
    assert_eq!(out.blocks.len(), 2);
    assert_eq!(out.blocks[1].start, 7);
    assert_eq!(out.blocks[1].data.len(), 8);
    assert!(r.read_vector(2, 13, "ch", None).is_err());
    let v = r.read_vector(7, 8, "ch", None).unwrap();
    assert_eq!(v.shape, vec![8, 2]);
    let z = v.into_complex64();
    assert_eq!(z[0].re, 7.);
    assert_eq!(z[1].re, 1007.);
    assert_eq!(z[15].re, 1014.);
    assert!(r.read_vector(7, 0, "ch", None).is_err());
    assert!(r.read_vector(7, 1, "ch", None).is_err());
    assert_eq!(r.read_vector(7, 1, "ch", Some(1)).unwrap().shape, vec![1]);
    assert!(r.read_vector(7, 2, "ch", Some(2)).is_err());
}
#[test]
fn missing_file_does_not_prevent_later_reads() {
    let t = tempfile::tempdir().unwrap();
    let c = properties(t.path(), 10, 1, 1);
    file(&c, 0, &[[0, 0]], &(0..10).collect::<Vec<_>>(), 1);
    file(&c, 2000, &[[20, 0]], &(20..30).collect::<Vec<_>>(), 1);
    let mut r = Reader::open([t.path()]).unwrap();
    assert_eq!(r.read("ch", 5..25, 0).unwrap().gaps, vec![10..20]);
    assert!(r.read_vector(5, 20, "ch", Some(0)).is_err());
    assert_eq!(
        r.read_vector(20, 5, "ch", None).unwrap().into_complex64()[0].re,
        20.
    );
}
#[test]
fn compressed_sparse_reads_reuse_full_file_sized_chunk_cache() {
    let t = tempfile::tempdir().unwrap();
    let c = properties(t.path(), 10000, 1, 1);
    file(&c, 0, &[[0, 0]], &(0..10000).collect::<Vec<_>>(), 1);
    let mut r = Reader::open([t.path()]).unwrap();
    r.set_cache_options(CacheOptions {
        minimum_bytes: 128,
        max_open_files: 2,
    })
    .unwrap();
    for start in [0, 9000, 20, 9990, 300, 8000] {
        assert_eq!(
            r.read_vector(start, 5, "ch", None)
                .unwrap()
                .into_complex64()[0]
                .re,
            start as f64
        );
    }
    let s = r.cache_stats();
    assert_eq!(s.file_opens, 1);
    assert_eq!(s.file_hits, 5);
    assert!(s.last_chunk_cache_bytes >= 40000);
    assert!(r
        .set_cache_options(CacheOptions {
            minimum_bytes: 0,
            max_open_files: 1
        })
        .is_err());
}
#[test]
fn absolute_sample_indices_do_not_round_through_f64() {
    let t = tempfile::tempdir().unwrap();
    let c = properties(t.path(), 1000000000, 1, 1);
    let start = 1712599210000000007;
    assert!(start > 1 << 53);
    file(&c, 1712599210000, &[[start, 0]], &[4, 5, 6], 1);
    let mut r = Reader::open([t.path()]).unwrap();
    assert_eq!(r.bounds("ch").unwrap(), Some((start, start + 2)));
    assert_eq!(
        r.read_vector(start + 1, 1, "ch", None)
            .unwrap()
            .into_complex64()[0]
            .re,
        5.
    );
    assert!(r.read_vector(u64::MAX, 2, "ch", None).is_err());
}
#[test]
fn rational_sample_rates_and_duplicate_directory_arguments() {
    let t = tempfile::tempdir().unwrap();
    let c = properties(t.path(), 3, 2, 1);
    file(&c, 4000, &[[7, 0]], &[7], 1);
    file(&c, 5000, &[[8, 0]], &[8], 1);
    let mut r = Reader::open([t.path(), t.path()]).unwrap();
    assert_eq!(
        r.read_vector(7, 2, "ch", None)
            .unwrap()
            .into_complex64()
            .iter()
            .map(|z| z.re)
            .collect::<Vec<_>>(),
        vec![7., 8.]
    );
}
#[test]
fn incompatible_channel_properties_are_rejected() {
    let t = tempfile::tempdir().unwrap();
    let a = t.path().join("a");
    let b = t.path().join("b");
    properties(&a, 10, 1, 1);
    properties(&b, 20, 1, 1);
    assert!(Reader::open([a, b]).is_err());
}
#[test]
fn conflicting_overlap_is_rejected() {
    let t = tempfile::tempdir().unwrap();
    let a = t.path().join("a");
    let b = t.path().join("b");
    let ac = properties(&a, 10, 1, 1);
    let bc = properties(&b, 10, 1, 1);
    file(&ac, 0, &[[0, 0]], &[0, 1, 2], 1);
    file(&bc, 0, &[[0, 0]], &[0, 9, 2], 1);
    let mut r = Reader::open([a, b]).unwrap();
    assert!(r.read_vector(0, 3, "ch", None).is_err());
}
