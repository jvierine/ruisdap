# Rust measured-waveform radar processing

Rust float32/complex64 range–Doppler–Doppler-rate matched filtering, following
[jvierine/hard_target](https://github.com/jvierine/hard_target), reference commit
`1b1d2b0a1f4b3b3ee99fd2f2e9ab0539f839d4bf` (MIT).

CPU: fused AVX complex products with runtime fallback, Rayon range searches and
planned RustFFT transforms. GPU: NVIDIA driver PTX kernels, bounded range tiles
and single-precision batched cuFFT. A pointwise kernel handles undecimated data;
warp reductions handle decimated products. CUDA is optional. Ion LPI uses the
existing complex128 weighted solver for conditioning; filtering and measured
products stay complex64, accumulators use float64, products are float32 HDF5.

## One staged recording, four products

The Arecibo dispatcher stages each one-second raw file once. Workers run meteor
CFAR, space-object CFAR, measured-transmit ion-line LPI and plasma-line matched
spectra against that file, then collect its products before deleting the staged
copy. Ion LPI and plasma processing run concurrently. Raw files remain untouched.
The work laptop j only transports data, manages queues and publishes products.
Scientific processing excludes j, avaruus and revontuli.

`arecibo_scan input.h5 output cpu|gpu meteor|space` writes CFAR detections only:
pulse UTC, exact source sample, slant range, IQ Doppler and SNR. Acceleration
search and CFAR settings are metadata. No hard-target matched-filter maps are
retained. The GPU downloads compact detections instead of full transformed
maps. Spectral and range CA-CFAR tests both must pass; the nominal per-cell PFA
assumes independent exponential training samples. Pulse-window correlations and
the total trial count are not included in that probability. These are detections,
not verified object classifications. Measured TX correlation and pulse support
set suppression scales for nearby maxima. Fixed-delay searches do not correct
range migration.

`arecibo_products input.h5 output cpu|gpu ion|plasma` produces:

* `ion/ion-lpi.h5`: signed real and imaginary LPI ACFs, variance, range/lag axes,
  filtered IQ, receiver mask, FIR taps and measured transmit templates. It uses
  actual measured TX waveforms, finite receive support and empirical pulse
  variance. No SVD or missing-lag interpolation. ADC rails are excluded. Moderate
  CFAR peaks must form a persistent track before they mask ion-line data; weak
  unclassified peaks do not censor the power estimate. The WLS solver uses a
  finite range grid; unconstrained columns remain NaN.
* `plasma/plasma-line.h5`: mean measured-waveform matched periodograms and pulse
  standard errors. Frequency is relative to the captured transmitter. The full
  25 MHz recording bandwidth is retained before display-bin averaging, and
  partially masked pulse energy is accounted for. These requested plasma-line
  products are distinct from the discarded hard-target maps.

Ion LPI is experimental: its provisional wideband-white receiver-noise estimate
comes from first differences and is propagated through the FIR. Calibration,
range boundaries, unambiguous-range contamination and covariance assumptions
need validation before plasma-parameter inference. LPI amplitude is in receiver
ADU², not calibrated electron density. Plasma spectra include receiver background
and the measured pulse ambiguity; display-grid spacing is not physical resolution.

The campaign recipe is restricted to measured 440 µs TX support, 10 ms IPP and
25 MHz DigitalRF v1 recordings. Unknown support or layout fails. Source sample
indices stay uint64; UTC seconds are float64. Positive Doppler means positive
recorded-IQ frequency shift. Slant range is never labelled altitude.

## Products and browsing

`bundle.json` records each product's status, duration, source recording and exact
executable SHA-256. Executables are pinned per bundle so deployments cannot change
an active analysis. The source HDF5 products have descriptive dataset names and
UTF-8 JSON provenance. Binary float32 tiles are display exports only.

The shared RUISDAP Rust/WASM/WebGL/rustplotlib plots browse individual chunks as
they arrive. The root index carries only the latest ten minutes; earlier chunks
are indexed by acquisition hour. CFAR views use detection scatter plots and a
four-column table, not full range–Doppler maps. Counts and links are also merged
into the main Rust data catalog.

## Validation

`cargo test -p ruisdap-hard-target --release` checks signed complex recovery,
CFAR local backgrounds, missing-equation WLS and conservative contamination masks.
`validate_cfar` recovers an injected delay/Doppler on CPU and GPU; `validate_power`
compares the pointwise GPU path and unequal frequency-bin averaging against CPU.
Real one-second recordings have produced all four products. Passing these checks
does not establish absolute calibration or validate the nominal global false-alarm
rate on this experiment.

Plasma-line diagnostic background suppression (`spectral-median-v1`): each range
spectrum has a reflected running frequency median subtracted (approximately
0.5 MHz full window, exact width in the manifest). HDF5 retains the original
matched PSD and standard error, background, signed residual PSD, and fractional
excess `(PSD-background)/background`. Negative residuals are retained. This is
not range–Doppler deconvolution or a detection-significance estimate. Broad
physical features can be suppressed; background-estimation uncertainty is not
included in the original standard errors. `plasma_highpass DIRECTORY` adds these
products to an existing plasma chunk without reopening raw IQ.

## MPI and clusters

See [CLUSTER.md](CLUSTER.md) for optional MPI delay-row distribution, bounded
per-rank CPU/CUDA resources, rank-local output, Slurm launching and validation.
