# Millstone Hill

Develop and test the GUI on avaruus, in `/home/j/src/guisdap26`.
Use `http://avaruus.uit.no/guisdap/`. Exported measurements live in
`/srv/guisdap/data`; source recordings remain on Revontuli.

The data browser opens a dataset once. Its plasma-line and ion-line views keep
one selected time and the loaded Rustplotlib renderers. All instruments use
`web/scientific-plots.js`: Rustplotlib fields, axes, colour bars, curves and
default Pagella typography. Millstone's plasma view shows the spectrum collapsed
over 185–560 km, the selected native frequency–range image, raw Ne map, density
profile and reviewed calibration points. The ion view shows Ne, Ti, Te, velocity and
the signed raw 10 microsecond lag, as maps and selected profiles. Range is slant
range. Existing density products remain explicitly uncalibrated proxies.

The April 8 quicklook contains 145 fits, 2,900 raw profiles and 573 plasma-line
integrations for ZENITH; MISA has 95, 1,907 and 380. Short gaps extend neighbouring
measured cells to their midpoint; gaps of one hour or more remain blank.
The first view covers 15–21 UTC; Full day and the normal Rustplotlib pan/zoom
controls expose the rest. Colour limits start from the IQR. Plasma contrast
preserves signed source estimates; sliders can include negative values.
Changing a view, colour limit or viewport does not fetch measurements again.
Selecting a new plasma integration or native ACF fetches that product on demand.
Scientific labels use the default Rustplotlib Palatino/Pagella family.

## Products

`scripts/isr_product.py` imports processed Millstone and legacy EISCAT HDF5 into
`ruisdap-isr-v1`. `/records/<id>` stores explicit instrument, UTC interval,
frequency convention, independent coordinates, fits and spectra. `/sources`
retains source trees and hashes. `/calibration` stores reviewed results.
`validate`, `export` and `calibrate` are separate commands. Imports and calibrated
versions refuse to overwrite existing products. The existing raw EISCAT MAT
reader remains available; unsupported layouts fail explicitly.

The Millstone density proxy is the source fitted zero-lag amplitude multiplied
by `(1 + Te/Ti) * range_km² / transmitter_W`. Raw Ne uses the real 10 microsecond
lag instead, without the temperature-ratio factor or receiver-gain division.
The source fitted amplitude already has receiver gain removed. A fitted-proxy
calibration must therefore not be applied directly to the raw proxy.
Neither is silently labelled
as calibrated electron density. `dTe/Ti` is ratio uncertainty, not `dTe`.

`scripts/export-millstone-day.py ROOT OUTPUT --date YYYY-MM-DD` makes the linked
quicklook from existing fit, LPI and plasma-line products. Its HDF5 stores all
displayed arrays, units and source metadata. Binary float32 arrays and JSON are
browser exports. `web/millstone-config.json` registers antenna/date exports;
`millstone-datasets.json` adds the observation days to the main catalogue.

`scripts/index-millstone-viewer.py SOURCE QUICKLOOK OUTPUT --date YYYY-MM-DD`
registers native products in the catalogue service's private
`millstone-sources.json`. Native ACFs are receiver-gain corrected; their spectra
use the same ACF Fourier transform as EISCAT. These individual integrations are
not the legacy fit's time/range average. The source does not retain model ACFs
or fitted ion composition, so an exact fitted spectrum is unavailable.

The plasma-line view estimates the absolute Millstone magic constant C from
included calibration points automatically after loading or changing the saved review, without launching reprocessing. New points are included by default; individual exclusions are retained. The source proxy is
`p = A0_K (1 + Te/Ti) range_km² / transmitter_W`. The definition is
`C p = Ne (1 + Te/Ti) / [(1+q)(1+q+Te/Ti)]`, with `q = k² lambda_D²`.
Thus C has units `W m^-3 K^-1 km^-2` for this chain; receiver gain has already
been removed from A0. In the small-Debye limit, `Ne = C p`.

Transmit power is **already divided out**, using each source fit's recorded
`P_tx` in watts, not the transmitter's nominal rating. Calibration validates
that this power is finite and positive, copies it to each checked point, and
retains it per fit in the calibrated HDF5 product. It must not divide the fitted
proxy by power a second time. Doubling signal and transmit power at fixed
plasma must leave the inferred C unchanged.

The legacy fitter uses an arithmetic mean of source `P_tx`, whereas it averages
ACFs with inverse-variance weights. This is not an exact power correction when
power and weights vary together. Historical products retain that definition;
correcting the averaging requires the source ACFs and a new fit, rather than an
extra scalar division in the calibration GUI. In a new fit, power belongs in
the measurement model (including its variance and range/lag response) before
combining measurements.

`POST millstone/estimate` verifies the exact saved notebook revision and derives
C from trusted source fits and plasma frequencies. It returns per-point C,
median C within each fit, the median across independent fits, and scaled MAD
between fits. That scatter is not total uncertainty; frequency, temperature and
systematic uncertainties are excluded. One fit cannot establish between-fit
scatter. Estimation creates no products and changes no active selector. Applying
C uses the same calculation in the authenticated density-recalculation endpoint.
Calibrated HDF5 products retain C and its units.

Calibration points save to the existing revisioned calibration service on
avaruus. Clicking an individual spectral feature adds an included point. Untick Use to
exclude it, or remove it from a linked plot. Numbered, outlined reticles keep
the echo centre visible. Click a point’s UTC time in the table to locate its
individual spectrum. The C plot labels the chain-specific units explicitly. The authenticated
`millstone/calibrate` endpoint recomputes gains from the saved revision and
trusted source fits. It publishes an immutable `ruisdap-isr-v1` HDF5 product and
new browser export, retaining the original source. Density uses the finite-Debye
ion-power correction at fixed source temperatures. This is a density
recalculation, not a new ACF fit. `catalog/millstone-config.json` selects the
active product; deployment preserves this server-owned selector.

Use `conda run -n base python` for the Python commands. On a compute host without
Conda, use a dedicated environment with NumPy/h5py rather than system Python.

## Native processing

`guisdap-millstone` reads an explicit `ruisdap-pulses-v1` HDF5 pulse pack. It
contains receive/transmit complex voltages, absolute sample keys, real paired
pulses, detection coordinates and explicit mode settings. The small DigitalRF
I/O bridges are `pack-millstone-pulses.py` for cancellation inspection and
`pack-millstone-train.py` for a bounded paired alternating-code control.

The Rust path includes measured-response fractional delay/Doppler cancellation,
noncircular filtering, sum decimation, lag products, clutter differences,
robust cuts, weighted lag inversion and native range–Doppler diagnostics.
Cancellation precedes lag products. Failed subtraction masks the echo support,
with guards expanded through the filter and both lag partners. Weak subtraction
and confirmed removal are distinct. Extremely weak amplitude estimates are left
unchanged. Receiver rail hits reject removal; analog compression is not inferred
solely from a large residual.

Build with `cargo build --release -p guisdap-millstone --features mpi`. The MPI
manifest is a small JSON object with `inputs`, `output`, and optional
`live_directory`. Each rank owns whole integrations and writes separate atomic
HDF5 outputs. Run `mpirun -n N target/release/guisdap-millstone process manifest.json`.
A non-MPI build refuses a multi-rank launcher. `cancel INPUT OUTPUT` inspects one
pack without pretending to have ground-clutter partners.

The live endpoint is `/guisdap/api/catalog/live`, integrated in `catalog-store`.
`RUISDAP_LIVE_ROOTS` lists shared status directories, separated by colons;
`RUISDAP_LEGACY_STATUS` lists existing status/queue files. Probing renews a
12-second lease. Small checkpoints continue independently, but workers generate
detailed arrays only while the lease is valid. The GUI pauses probes while
hidden or explicitly paused. Give the two host users access only to the shared
status directory when their numeric user IDs differ.

## Evidence and limits

- All workspace tests pass, including five voltage-cancellation cases.
- A 74-pulse real control agrees with independent NumPy/SciPy lag inversion to
  1.23e-13 relative to peak, including variances, receiver gain and retention.
- One- and two-rank Revontuli runs match exactly for 28 numerical datasets.
- Separate real 14.8 and 35 dB echoes pass removal checks. The 50 and 78.7 dB
  cases fail residual checks and are masked. No ADC rail hits occur in these
  four echoes; the bright failure mechanism remains unresolved.
- The six plain reference tables are in `jvierine/isr_spec_validation` and use
  Kudeki–Milla Part I. The native theory also agrees with an independent
  implementation over 162 collisionless, unmagnetized Maxwellian cases.

The full-day driver covers the declared 16 × 30 microsecond ZENITH mode.
Support for all other pulse modes and exact parity with all legacy spectral
averaging/rejection rules remain unfinished. Cancellation's
signal-loss and covariance effects still need end-to-end injections. The
recalculation queue freezes reviewed EISCAT calibration revisions; a validated
Millstone production refit recipe is not yet registered. These are not supplied
by a successful display or a short MPI pilot.


## Native run on 19 October 2026

The explicit driver `scripts/prepare-millstone-native-run.py` selects observed
monostatic coded pulses, validates recorded power, retains real same-code
clutter partners, and imports the existing range–Doppler/CFAR detections. It
measures the TX-to-receiver delay instead of assuming eleven samples. The
current recipe is the 16 × 30 microsecond coded mode, 1 MHz sampling, a
30-sample range basis and sixteen measured lags from 10 to 160 microseconds.
Other modes still require their own declared timing recipe.

On avaruus, 128 distinct one-second integrations (9,472 paired pulses) were run
with 128 MPI ranks bound to 128 physical cores. The slowest rank took 58.5 s,
excluding preparation. All requested ACF cells were finite. Of 1,096 seeded
pulse echoes, 924 passed the removal checks, 37 received weak subtraction,
and 135 failed (90 held-out excess, 44 coherent residual, one edge residual).
Failed supports were masked from lag inversion. These are algorithm decisions,
not independent ground truth. No target window hit the recorded ADC rails.

Pulse normalization is explicit in `lpi.code_energy_samples`, one value per
pulse at unit peak transmit power. The unit-energy ambiguity row is multiplied
by that energy and `range_cell_km = c * decimation / (2 * sample_rate)`.
The recovered ACF has units K/km; the corresponding calibration has units
W m^-3 K^-1 km^-1. For a fixed-code old product,
`C_native = C_legacy * code_energy_samples * range_cell_km`. Per-pulse scaling
supports mixed code lengths. The inverse is checked against independent
peak-power forward sums for 30 and 100 microsecond codes, mixed codes, and
10/30/100-sample range bases. This does not validate antenna efficiencies or
all receiver effects, and historical products keep their original definition.

The native `fit MANIFEST` command also supports MPI. Its manifest supplies
inputs, output, integration_s, magic_constant_by_antenna, ions, range_km,
initial plasma state and fit settings. It removes each source's recorded
transmit power before inverse-variance averaging, includes the finite receive
FIR in the common ISR spectrum model, and writes canonical HDF5 fit arrays,
model ACFs, full parameter covariance and standard errors. Calibration,
composition, and shared-lag covariance are not included in those error bars.
The first comparison produced four 180-second ZENITH profiles over 200–540 km,
with explicit O+ composition. MISA was not assigned the ZENITH calibration.
The products are experimental; the historical GUI dataset has not yet been
replaced by these new fit products.

BEATA native power is exported by
`cargo run --release -p guisdap-worker --example export_beata_power -- INPUT_FRAMES NEW_OUTPUT`.
The HDF5 retains 40,166 native profiles and their original altitude grids.
Five resolutions of 256-column browser tiles contain arithmetic means before
logarithmic colour mapping. The browser fetches only the visible resolution,
keeps at most twelve tiles, and exposes approximately 0.357 s profiles when
zoomed in. Slice timing is assigned uniformly within recorded five-second
bounds; pulse-level timestamps are not available in those products.

## Viewing historical and native analyses

The shared Millstone viewer now lists server-registered analyses in its Analysis
selector. Historical opens the full UTC day; the native trial opens its available
coverage. Switching preserves selected time and slant range. Both use the same
Rustplotlib parameter maps, profiles and instrument diagnostics.

`export-millstone-native-view.py` exports canonical native HDF5 fits into that
viewer. Native Ne and its errors are already absolute and are not multiplied by
the historical calibration again. Matching measured/model ACFs and component
errors come from the exact fit integration; the displayed spectrum applies the
same measured-lag DFT to data and model. The native negative-lag convention is
conjugated for the shared positive-lag DFT. HDF5 downloads retain full covariance
and processing metadata. Historical raw Ne and plasma spectra remain explicitly
separate context. Calibration review editing continues in Historical.

The four-profile trial is registered as `native-20261019`. A full-day ZENITH
30-us-code run was launched separately with `process-millstone-full-day.py`:
preparation scans all 24 hours with one reader per hour; reduction and fitting use 128
MPI ranks on physical cores, with experimental CFAR-seeded subtraction enabled.
It publishes `native-full-20261019` only after all eligible inputs reduce and fit
successfully. This extends time coverage for the declared recipe, not support
for other pulse modes or a MISA calibration. Failed processing does not publish
a partially complete analysis, and the service has no automatic restart.

Explicit `--resume-products` verifies the frozen input SHA256, common-format
schema and native executable identity before retaining completed products. The
remaining inputs are distributed by descending CFAR-detection work estimate.
Resume audits are HDF5; no completed waveform or ACF product is overwritten.

The original 8 April sampled trial is registered as `native-full-20261019`:
5,540 selected integrations produced ACFs, and 295 time buckets produced native
fits, spanning 00:16:40–23:59:56 UTC. The viewer was checked against canonical
HDF5 arrays and every input is represented once. Across 22,125 fitted range
gates, 11,863 converged, 1,905 were low-information, 4,399 did not converge and
3,958 reached parameter boundaries. These flags are retained; completion is not
a claim that every gate is a usable fit or that cancellation bias is validated.


### 2026-10-19 coverage audit and corrected eclipse processing

The original `native-full-20261019` export is an incomplete sampled trial,
not an all-pulse day reduction. Its input builder required an exact per-second
power-meter file, while the archive records power at five-second cadence. It
therefore retained only 5,540 one-second coded packs (406,172 pulses; 5,385 s
wall exposure), omitted the long-pulse waveform, used 16 lags, decoded only
126–535 km and fitted only 200–540 km. The catalogue now labels it incomplete;
its immutable products remain available for inspection.

The corrected driver is `scripts/process-millstone-eclipse.py`, launched in
`ruisdap-millstone-eclipse-all-20261019.service`. It processes April 8, April 7
and April 9 sequentially, with 128 physical MPI cores and experimental
CFAR-seeded complex-voltage satellite subtraction. The server archive covers
April 7 about 10 UTC through April 9 about 11 UTC; the boundary days are partial.
A frozen executable and full source pulse-pack hashes identify the numerical
run. Input preparation interpolates the recorded power-meter series, retains
real cross-boundary clutter partners without counting them as science pulses,
and records an exact unique science-pulse audit in HDF5.

A checked 18:00:10–18:00:20 UTC pilot contains all 748 coded and 374 long science
pulses in the metadata window, plus three separately identified clutter partners.
Both native reductions contain 240 range cells (0–1,074.756 km) and 45 measured
lags (10–450 microseconds). At 200–500 km the median real short-lag amplitude
ratio to the response-normalized historical ACF was 0.959 in this pilot. That is
an input/normalization check, not validation of full-day fitted plasma parameters
or of satellite-removal bias. Negative and unidentifiable estimates remain
stored with their quality flags.

Campaign fits use the actual historical integration windows, the original
range-averaging radii 2/3/5 at 0–300/300–700/700–1500 km, below-300-km paired
clutter ACFs, the Holt 32/16-amu molecular-ion prescription, three initial guesses
and historical Ti, Te/Ti and velocity bounds. Above 700 km the current Rust
implementation profiles O+/H+ on a composition grid. This is explicitly distinct
from the old continuous composition fit. The physical finite-filter,
finite-Debye forward model and propagated diagonal measurement errors are also
retained rather than copying historical free-amplitude and error approximations.
The legacy slice is end-exclusive: its actual interior averages contain 4/6/10 cells, despite the stored nominal widths of 5/7/11 cells. The matching path uses the actual slices. The archive's fourth-moment, shared-lag and range covariance remain incomplete.

Native analyses now support independent calibration notebooks, automatic C
estimation using the measured complex ACF and the exact fitted receiver/model
response, and authenticated application through a new 128-core ACF refit.
The estimator evaluates the Bohm–Gross plasma-line density at conditional
source Ti, Te, velocity and composition and solves the weighted real amplitude
projection for C. It verifies source HDF5 hashes, rejects unmatched/nonconverged
points explicitly, groups by native fit, and reports between-fit scaled MAD.
Applying C preserves the original analysis, writes canonical HDF5 fits and a
calibration audit, and adds a separate catalogue choice. It does not merely
rescale fitted Ne. Native C has units W m^-3 K^-1 km^-1; initial campaign C is
transferred from the frozen reviewed April 8 same-antenna calibration. Other
days' context is exported from their own archive products and labels that transfer.

Instrument traces share the common ISR plot, bridge gaps shorter than an hour,
and use blue azimuth, purple elevation, green Tsys and orange transmit power.
The shared view retains the same inner time-axis extent as the parameter maps.

Native calibration refits run as tracked server jobs; the UI polls only while the calibration view is visible, resumes an existing job after reload, and reports interrupted jobs after server restart. The campaign progress heartbeat exits with its explicitly launched systemd job and never schedules a restart.

## Bounded real-data benchmark (19 October 2026)

Memo 39 records a ten-second, 1,122-science-pulse ZENITH workload: 748 coded
and 374 long pulses, 45 measured lags and 240 requested range cells. On one
pinned EPYC 9754 physical core, three-repeat median LPI times were 21.403 s
(coded) and 10.503 s (long). An independent NumPy/SciPy replay of the same
operator took 53.077 and 25.438 s. The ACF peak-relative discrepancies were
5.99e-14 and 2.67e-11; variances, clutter ACFs, retained fractions and receiver
gains were also checked. This is not a direct historical-controller benchmark.

The reproducible harnesses are `examples/benchmark_lpi.rs`,
`examples/benchmark_cancellation.rs`, `scripts/benchmark-millstone-lpi-reference.py`
and `scripts/benchmark-millstone-pipeline.py`. `scripts/collect-millstone-benchmark.py`
freezes timing arrays, configuration, source hashes and numerical checks into
HDF5. The bounded snapshot uses explicit avaruus pilot paths; use a fresh
output version when repeating it.

Two complete process-command repeats took 33.847–33.938 s with one MPI rank
and 23.380–23.448 s with two ranks. There are only two unequal jobs: do not
extrapolate this measurement to 128 ranks. The pilot has zero CFAR seeds.
Separate five-repeat cancellation measurements of the three Memo 38 echoes
were 0.793, 0.789 and 0.774 s/echo, retaining the accept/reject/reject decisions.
These exclude catalogue construction and input I/O. Fresh archive preparation
was 47.25 s including three hours of metadata discovery; this is startup cost,
not an hourly-worker steady-state cost. A single short-pilot plasma fit took
134.47 s; the low convergence count is not representative of scientific
200–300 s integrations.

## Production campaign: direct DigitalRF reads

Use `scripts/process-millstone-eclipse-direct.py`, passing one or more `--root`
recording directories. The native `process-archive` command reads the original
DigitalRF files by absolute sample index; no campaign-wide preparation, pulse
packs, combined RF directory, converted RF or persistent RF cache is used.
The older packing controllers are historical experiments and must not be used
for production campaign processing.

`guisdap-digital-rf::Reader` accepts multiple top-level, recording or channel
directories. It exposes `read_vector(start, length, channel, subchannel)` with
upstream contiguous-vector length, subchannel squeezing and safe floating-point
promotion semantics. File cadence is used only to find candidates;
`rf_data_index` determines the actual stored blocks and gaps. Sample indices
remain u64 and cadence arithmetic uses u128, including rational sample rates.
Missing files and indexed gaps do not become zeros. A vector containing a gap
returns a typed `MissingSamples` error; the campaign audits that pulse, skips it
and continues with subsequent valid data. Corrupt or unsupported layouts remain
explicit errors. Local filesystem/NFS DigitalRF 2.x with rational-rate properties
is supported; network URL transport and pre-rational-rate property layouts are
not implemented.

HDF5 files remain open in a bounded LRU. Each raw chunk cache is at least the
uncompressed logical `rf_data` size and at least 64 MiB by default, with hash
slots sized to avoid chunk-cache thrashing. Campaign workers retain up to 32
files, keeping both receiver and transmit files warm across a ten-second
observation. This is an in-memory read/decompression cache, not a disk archive.
Direct DigitalMetadata reads preserve the last antenna state until the next
actual event, interpolate nearby recorded power, and use catalogue delay/Doppler
seeds. The CFAR cancellation and LPI kernels are shared with the validated tests.

MPI workers dynamically claim ten-second windows, so busy satellite passes do
not leave work permanently assigned to one rank. Only reduced ACFs, variances,
receiver calibration, cancellation diagnostics, mask intervals and coverage
audits are written. Atomic per-window HDF5 checkpoints verify output hashes
before resuming. Fit checkpoints bind both the exact fit manifest and executable
hash, and verify their reduced source inputs. The original RF is never modified.

Validation compares the native reader to upstream `digital_rf` 2.6.14 across
compressed synthetic multi-directory/hour boundaries, missing files, real
Millstone file boundaries and the actual 28-second receive gap. The real
1,122-science-pulse pilot gives exact ACF, clutter-ACF, variance, retained-row and
receiver-gain equality with the independently validated operator. A known 35 dB
real satellite echo remains accepted with a held-out RMS/background ratio 0.957.
The reader tests also cover indices above 2^53, rational sample rates, incompatible
properties and sparse compressed reads with cache reuse.

### Independent pulse coverage and invalid receiver backgrounds

`audit-millstone-source-coverage.py` independently reconstructs eligible pulses
using upstream DigitalRF/DigitalMetadata readers, original antenna/power/code
metadata and continuous RF block indices. It compares exact candidate and
paired-science sample keys with every ten-second reduction audit. Full-day
publication requires a verified complete-day HDF5 reference report. The auditor
can resume verified chunks, checking the hashes of previous reduction audits;
it reads no RF voltage samples and writes no RF representation. Synthetic tests
exercise a real indexed RF gap, loss of its same-code partner, and rejection of
a modified audit checkpoint.

Flat recorded receiver backgrounds are explicit per-pulse quality failures.
They produce an all-false usable mask, preserve the input voltage, and have
undefined noise/validation ratios recorded as null. They do not supply receiver
calibration statistics or abort unrelated windows. At the April 8 receive-gap
boundary, 76 such pulses were checked against the original upstream RF reader;
their background samples are actually constant in the recording. Ordinary coded
and long-pulse control ACFs, variances and receiver gains remain bitwise equal.

### Subtraction review and cancellation timing

The current campaign subtraction remains experimental. Two complete matched
integrations, using identical source pulse sets and fitting kernels, show median
temperature changes below 0.3% over 250–500 km. In the separate 525–575 km check,
median ion-temperature changes are 5–8%, and velocity changes reach 31 m/s.
Parameter-boundary fits remain with subtraction both enabled and disabled.
Neither branch is a ground-truth plasma measurement.

`check_diffuse_null` synthesizes zero-mean Gaussian ion-line fluctuations from
the native spectrum model with no coherent target. With a *forced* candidate,
the legacy canceller accepted subtractions in 26/100 long-pulse and 22/100 coded
trials at diffuse/receiver power ratio 0.1. This is a conditional-candidate
stress test, not an end-to-end CFAR or real-archive false-alarm rate.
Alternating sixteen-sample training/holdout blocks do not ensure independence
for colored plasma scatter. An opt-in Bartlett HAC amplitude-variance prototype
(`correlation_lag_s`, in seconds) reduces those counts to 8/100 for each family
at a 250 µs bandwidth, but still falsely confirms some higher-power null cases.
It remains disabled by default and in the campaign. Complete-day publication is
held pending the detector/subtraction review. Weak amplitudes remain labeled
weak in the viewer even when a different echo in the same pulse is strong.

Cancellation interpolation now skips only identically zero transmit support
and reuses the separable delay interpolation and Doppler phase grids. Regression
tests compare every template output bit with the original numerical definition,
including fractional shifts, noncircular edges and signed zeros.
`benchmark_direct_cancellation` measures only cancellation, excluding archive I/O
and LPI, and writes diagnostics, timing and output/mask hashes without voltages.
Three repeated synthetic and real coded/long single-pulse cases on one allocated
core showed approximately 30× kernel speedup with exact output and diagnostic
equality. This is not a measured 30× speedup of the whole campaign pipeline.

### Incremental fits, legacy timing and debris diagnostics

`stream-millstone-direct.py` fits completed integrations while archive reduction
continues. The integration algorithm is taken from `fit_lpi.fit_lpifiles`:
anchor at the first available ACF product and include source starts strictly
below anchor + 300 seconds. Antenna switching commonly limits ZENITH exposure
to about 200 seconds. Historical `t1` is the **start** of the last ten-second
source product, so its common 190-second recorded span represents about 200
seconds of observations. Native fits record the actual science-pulse extent.
The last integration is flushed at a day boundary; the old Python loop omitted
its trailing group. A missing reduction audit blocks the contiguous watermark;
an explicitly audited missing-RF window does not. This prevents publishing an
incomplete integration as if it had no further observations.

The live scheduler uses eight lower-priority fit workers on physical cores
already allocated to the 128-rank reduction. Final draining uses up to 128
workers. Fits have immutable, per-integration plans and checkpoint hashes;
per-integration leases prevent duplicate fitting between the live and final
drivers. Public exports use immutable revision directories and switch the
manifest atomically after all assets are present. The viewer checks for new
revisions every 15 seconds while visible, preserves selection and viewport,
and labels the analysis live until the audited complete day is published.
No new RF reader or RF staging step is introduced by this scheduler.

The native exporter writes a compact, UTC-relative four-column browser array
of debris detections: time offset in seconds, slant range in km, status and
signed Doppler in Hz. Original pulse keys and full diagnostics remain in reduced
HDF5. Only science pulses are exported, excluding repeated clutter partners.
The delay is relative to the measured leakage template with the measured
TX/RX channel timing offset removed. The viewer places the detection panel
below raw electron density, with identical inner horizontal plot bounds and
time navigation. It distinguishes confirmed subtractions, weak subtractions,
failed-validation masks and retained low-SNR candidates. These are CFAR-seeded
candidates and subtraction quality states, not independent orbit confirmations.
At broad zoom, one actual event per screen bin and category is drawn; zooming
reveals the native event coordinates. No observational values are averaged or
removed from the authoritative HDF5 products.

The first live-fit comparison exposed a high-density initialization failure:
on the same 301 km ACF in the first April 8 integration, a fixed 1e11 m^-3
starting density entered a hot-ion local minimum with reduced chi-square about
50.5. A start near 1e12 m^-3 found a much lower residual near 5.0. Historical-mode
fits now include density starts profiled from the measured complex amplitude
through the identical finite receiver/model response. Two amplitude updates
account approximately for density-dependent Debye response before the full
nonlinear solve. The checked real-data regression yields reduced chi-square
5.08 and Ti about 1303 K; its remaining temperature-ratio boundary flag is
retained. This is evidence of improved initialization, not proof that remaining
large residuals, lag covariance or composition uncertainty are resolved.
The fitting-kernel revision is frozen separately from reduction kernels, so
valid ACF reductions can be reused without repeating RF processing.

Background-gate contamination now rejects the entire affected pulse and omits
it from receiver noise/injection gain estimation, rather than terminating the
window. Rejected pulse keys, masks and undefined diagnostics are explicit in
HDF5; undefined ratios are JSON null. An ordinary coded/long pilot remained
bitwise identical in ACF, variance, clutter ACF and receiver gain after this
guard. A previously failing archive window completed with two fully rejected
pulses. Every reduction kernel is identified in the product/coverage ledger.
Work counters include their manifest start/end interval, so bounded repairs
cannot reset the full-day MPI queue. Earlier duplicate claims are rejected by
output guards; checkpoint reuse verifies immutable product hashes.

### 2026-10-20 coverage and residual checks

The direct April 8 reduction passed the independent full-day metadata and RF
block-index audit: 3,099,908 science pulses from 3,103,238 candidates. The
reduction separately records 82 science pulses rejected in full for unusable
background calibration. All reduced products are represented exactly once
across 145 completed fit integrations. This establishes source coverage, not
scientific validity of every fit. Complete-day publication remains held while
the subtraction review is open; the controller proceeds with April 7 and 9.

The external CFAR detector was also run unchanged on synthetic diffuse-only
Gaussian ion-line returns. Thin-layer coded returns can produce candidates
without coherent targets; distributing the same nominal return over a broad
range profile greatly reduces this effect. These model-dependent tests do not
measure the campaign's empirical false-alarm rate. Direct checks of four
ten-second source windows found marginal excess RF power around 0.13–0.23 of
receiver background near 300 km and about 0.01–0.04 near 550 km in selected
no-seed pulses. Absence of catalogue seeds does not establish absence of debris.

`check-millstone-residual-cfar.py` reconstructs recorded subtraction directly
from the original RF, refined templates and recorded complex amplitudes in
memory. It requires agreement with the native support RMS to relative 1e-9,
then runs the unchanged detector on original and residual voltages. It writes
only scalar QA, candidate coordinates and provenance in HDF5. In eight targeted
windows selected for high counts of accepted or masked 525–575 km candidates,
39 of 511 accepted subtractions retained a CFAR-level residual in that range.
None of those residuals differed by more than 4 kHz from every fitted Doppler.
Among 329 failed-validation/masked pulses, 133 retained residuals, including 24
with that Doppler separation. None of 55 no-seed controls had a layer residual.
These are targeted samples, not representative campaign fractions.

The detector's range-only peak suppression can omit components within one
pulse-length range neighbourhood. The native held-out check only tests the
selected templates, so a small projection onto them does not prove an empty
residual ambiguity surface. Neither finding justifies adding arbitrary
residual candidates in production: such additions can also remove diffuse ISR.
The `check_residual_candidates` example provides a bounded joint-fit experiment
on the original source without changing production seeds or writing voltages.
In that experiment, adding the observed residual components cleared the
CFAR-level layer residual in 37 of the 39 targeted cases; all 39 joint fits
passed the native checks and two still had a detectable residual. Additional
amplitude SNR ranged from 11.9 to 24.0 dB under the existing white-noise
uncertainty assumption. This demonstrates a possible fit improvement, not
plasma preservation or a confirmed second object. Diffuse-null and injected
target controls remain required before adopting this in campaign processing.

The subsequent CFAR/native-cancellation control run contains 2,400 synthetic
trials across coded/long pulses, thin/broad Gaussian diffuse layers, six diffuse
power levels, and single/two-target injections. The native canceller receives
voltages through anonymous pipes; neither RF voltages nor pulse packs are
written. At nominal diffuse/receiver power 0.2, the coded thin-layer null
produced nine accepted subtractions in fifty trials. With one genuine injected
target, its paired real 10 microsecond ACF error was -18.7% of the diffuse
reference before residual augmentation and -33.1% afterwards. The corresponding
broad-layer coded single-target result was -4.6% with no augmentation. These
are direct-voltage lag moments, not complete LPI inversions or measured
campaign density errors. They reject automatic residual augmentation as a
validated solution and retain the completed-day scientific hold.

For fixed templates the amplitude stage is a linear projection, even if a
coherent target really exists. If `Q = I - P`, then the remaining diffuse
covariance is `Q R Q^H`; using the original unprojected plasma response can
therefore bias the recovered parameters. Interleaved within-pulse training
and validation samples do not provide independent plasma realizations. A
small residual and a successful native validation flag do not remove this
effect. Geometry refinement and candidate selection introduce further
data-dependent effects beyond this fixed-template identity.

The saved-amplitude phase diagnostic found normalized inter-pulse coherence
about 0.93 in one targeted F-region coded track and 0.94–0.99 in two targeted
550 km coded tracks. Each exceeded all 64 random-phase controls under the
same frequency/chirp search. That is evidence of coherent structure, not an
object identification, a resolved Doppler alias or a campaign false-alarm
rate. Long-pulse endpoints lacked the required sixteen neighbouring pulses.
It would be incorrect to identify every F-region candidate as diffuse plasma
solely from its location.

An independent-pulse amplitude control supplies the known injected geometry
and aligned phase, estimates each pulse's amplitude from the other pulses,
and records paired ACF changes. Independent groups, rather than pulses that
share a training pool, are the bootstrap units. This isolates the projection
effect but does not supply a deployable tracker: phase/geometry prediction,
false-candidate control and full LPI/fit validation on actual observations are
still required. Production subtraction settings remain unchanged.

The oracle test completed 24 cases, each with eight independent groups of
sixteen pulses. With fixed known geometry, same-pulse amplitude estimation had
mean 10 microsecond ACF errors from -3.8% to -36.7% of the diffuse reference.
Leave-one-pulse-out prediction changed these point estimates to +0.2% through
+3.7%; the bootstrap interval extrema across cases were -0.7% and +5.8%.
This demonstrates removal of the negative projection effect under the stated
oracle assumptions, with a remaining positive estimation-noise contribution.
It does not validate phase tracking, candidate confirmation or full LPI fits
on real observations.

### Receiver calibration exclusions

April 7's first direct reduction stopped at 21:20:40 UTC when the long-pulse
injection measurement was below the receiver background. The calibration
guard now distinguishes this measured failure from a numerical or unsupported
layout error. It keeps the coded product and excludes only the uncalibratable
long-pulse family, recording the original pulse keys, configuration, complete
cancellation diagnostics and measured powers in the window audit. It does not
invent or interpolate a gain. Unrelated windows continue.

In that window, the filtered noise and injection powers were 216.0171 and
145.4731 ADC-squared respectively, giving gain -0.0601911 ADC-squared/K.
An independent upstream DigitalRF read and NumPy FIR/calibration replay
reproduced the failure for 374 excluded long-pulse science records. The coded
product's SHA256 stayed unchanged. Positive-gain arithmetic is unchanged, and
an invalid injection-temperature configuration still raises an error.

The source-coverage verifier compares the independently expected pulse set to
the union of reduced science pulses and explicitly calibration-rejected
pulses. It independently reconstructs every claimed calibration exclusion
from original RF, including checks that valid calibration pulses were not
silently omitted. Day products separately record readable, reduced and
calibration-rejected pulse counts. A quality exclusion therefore remains
visible without masquerading as a missing source recording or valid density.

April 9's independent audit caught a second exclusion-accounting case at
10:30:50 UTC. Missing RF left 74 paired coded science pulses and 37 long
science pulses, but only eight coded and five long pulses had non-flat
receiver-background gates. A complete family requires sixteen usable
calibration pulses. The reducer already skipped those families; its audit
now explicitly records all 111 excluded science keys, the retained
calibration keys, configuration and pulse diagnostics. The upstream reader
independently replays every claimed flat background and verifies that no
valid calibration pulse was omitted. The existing reduced products and
fits are unchanged. These exclusions are separate from the day's 374
nonpositive-gain exclusions.

### Independent-pulse subtraction development

The experimental Rust track forecaster removes the queried pulse before
association, geometry regression and complex phase/amplitude estimation.
It associates nearby CFAR-derived events from other pulses, fits quadratic
delay and within-pulse Doppler, and predicts a linear complex envelope with
a frequency/chirp phase model. The phase search uses sparse FFT sums of
amplitude metadata; absent pulse epochs carry zero weight, with no invented
RF samples. A local Nelder–Mead refinement follows the FFT grid search.
The phase frequency is aliased on the measured pulse clock; it is not
substituted for the independently fitted within-pulse Doppler.

At least sixteen other pulses are required. Tests use neighbourhoods of
0.35 seconds for coded pulses and 0.5 seconds for long pulses, at most
48 training events, and training explained fraction at least 0.9. Unsupported
timing, insufficient support and multiple associated components in a
training pulse produce explicit unsupported results. Held-pulse metadata
mutation tests leave the forecast unchanged. On 64 fixed real-data queries,
Rust and an independent Python implementation agree on support/training
membership and forecast tolerances; the largest relative amplitude
difference is 5.8e-8. This numerical agreement is not scientific validation.

A non-oracle detector control passes synthetic voltages through the unchanged
CFAR detector and native single-pulse amplitude estimator, then gives only
the resulting event metadata to the native forecaster. Four independent
groups of 32 pulses per case cover coded/long pulses, thin/broad diffuse
layers, and absent/present moving targets. All 512 injected-target pulses
receive forecasts; none of the 512 target-free pulses receive a subtraction.
The limited null sample does not establish a campaign false-alarm rate.
These controls save paired lag moments, not raw RF or full LPI products.

Two further controls inject a declared moving target at slant ranges 300
and 550 km into the same quiet original-RF window, 2024-04-08 18:00:10 UTC.
All 750 coded and 375 long pulses receive independent forecasts. These
candidate-controlled tests use coarse seeds; they are not CFAR detection
efficiency tests. They run all 45 native lags through the actual LPI,
including the clutter-pair branch, recorded noise-injection calibration,
physical pulse-energy/range normalization and propagated diagonal variance.
Original ACFs, clutter ACFs and variances match the canonical campaign
products exactly. Calibration gains are identical across all four branches.
No raw or cleaned voltages are written.

For example, over the 30 km neighbourhood of the 300 km injection, the
coded real-ACF change at the measured 20 microsecond lag is -15.44% with
same-pulse subtraction and +0.66% with independent prediction, relative to
the mean absolute original real ACF in that neighbourhood. At 10
microseconds those figures are -7.66% and -0.13%. Keeping lag coordinates
explicit matters: the first native lag is 10 microseconds, not zero.

The frozen campaign fitter was then run on the original, uncorrected
injection, same-pulse subtraction and independent-prediction ACF branches.
These are ten-second controls with the same legacy range averaging and
calibration, not full campaign integrations. Median fitted changes from
the identical original RF within 30 km of each injected range were:

| Injected slant range | Method | Density change | Ti change / K | Te change / K | vi change / m/s |
| --- | --- | ---: | ---: | ---: | ---: |
| 300 km | Same-pulse | +6.94% | -255.6 | +166.1 | +16.31 |
| 300 km | Independent | -0.82% | +52.5 | +10.0 | -0.88 |
| 550 km | Same-pulse | -1.81% | -17.1 | -74.6 | +15.58 |
| 550 km | Independent | -2.57% | +24.3 | +35.1 | -4.44 |

These results show improvement in several parameters, with remaining errors
and no uniform improvement in every quantity. They use one quiet measured
window, not independent campaign-wide trials. Prediction variance is
conditional on the fitted phase/geometry and does not supply a complete
uncertainty model. Held-RF pilot checks also retain some detectable residuals.
The forecaster remains outside production preparation; the scientific
publication hold stays in place pending full-integration real-event,
false-candidate and coverage tests in a fresh immutable product version.

### Requested conservative same-pulse removal

Neighbouring-pulse echo prediction has been rejected for campaign processing.
The processing path fits and subtracts echoes from the pulse being processed.
The diagnostic experiments above do not define the production method.

New archive runs require the recorded CFAR ratio to reach twice the recorded
CFAR threshold, a 3.01 dB increase in power threshold. They only admit targets
at **400 km altitude or above**. The manifest explicitly records
`debris_cfar_threshold_multiplier` (default 2) and
`debris_min_altitude_km` (default 400). Original candidate statistics and
selection decisions are retained in product metadata.

This reader supports monostatic vertical ZENITH. For that geometry, altitude
is slant range plus the campaign station height of 0.146 km. The cutoff
therefore corresponds to 399.854 km slant range. It is applied after the
measured transmitter/receiver timing correction. Delay refinement is
constrained by the same boundary; targets below it are neither subtracted
nor masked as debris. This is a target-altitude selection, not an ion-fit
range restriction. It does not claim that a finite waveform has no response
at neighbouring plasma gates.

Each output directory has a cancellation-policy ledger. Changed policies or
legacy checkpoints require a fresh product version rather than silently
reusing products generated with different thresholds. The direct DigitalRF
reader, pulse selection, calibration, lag inversion and fitting settings
remain unchanged. No neighbouring-pulse option is accepted by the manifest.

## isr_analysis analyses (Henrik, 26 October 2026)

`scripts/export-isr-analysis-view.py` writes one day of an
[isr_analysis](https://github.com/ionodev/isr_analysis) run (the FYS-3730
satellite and RFI mitigation work) in the quicklook schema, as an extra
analysis beside Historical and the native runs. It reads that run's
`fit_lpi.py` products (`lpi_60/<channel>/pp-*.h5`) and its 10-second
inversions (`lpi-*.h5`). Ne is the fitted proxy times the run's own density
constant (7.15e8 W m^-3 K^-1 km^-2 for its lpi_60 products, from the
digisonde, isr_analysis memo 43), with the viewer's finite-Debye correction;
raw Ne is the run's own 10 us lag with Te/Ti = 1. The plasma-line panel
references the day's existing export. `scripts/register-millstone-analysis.py`
adds or removes one analysis in `catalog/millstone-config.json` under the
config's lock file.

First exports (7 and 8 April, ZENITH, 11-23 UTC), staged on revontuli in
`/mnt/data/henrik/viewer_exports_2026-10/`: `isr-nosat-v1` (variance-weighted
LPI with a background column, no satellite handling), `isr-detected-v1`
(satellite echoes as extra unknowns on the detected pulses), `isr-tracked-v1`
(the same with pass tracking, isr_analysis memo 45). Publishing them needs
write access to `data/` and `catalog/`.

`web/compare.html` compares two analyses of one day before and after
mitigation: each parameter map is drawn from both, on the same axes and
colour limits, and a vertical bar (drag it, or use the slider) shows the
"before" analysis left of it and the "after" analysis right of it; pan and
zoom are shared. Where the range gates differ, both are shown on the finer
grid, each cell taking the value of the gate it lies in. At most four
parameters are shown at once: each map uses two WebGL contexts, and browsers
drop the oldest beyond about sixteen. Native (Rust) analyses are compared
with their absolute Ne as stored; their raw-Ne row is left out, since a native
analysis shows its source's raw Ne rather than its own.
Example: `compare.html?date=2024-04-08&before=isr-nosat-v1&after=isr-tracked-v1`.

Quality flags and the eclipse (27 October 2026). With `--isr-analysis
<checkout>` (and `--echo-rows`, a .npy of the catalogue's echoes, or
`--echo-catalogue`) the export writes each fit's quality flags as
isr_analysis's `fit_quality.py` defines them (`fits[].quality_flags`, one
integer per gate; FAILED 1, BOUND 2, NOISY 4, SATELLITE 8 for information;
the definition in the manifest's `quality_flags`). The exports with flags
are the `-v2` folders (parameters identical to `-v1`), and the private copy's
menu points at them (ids `isr-*-v2`). `compare.html` has a checkbox "Hide
flagged gates" (off by default; `&flags=hide` in the URL): each side hides
its own flagged gates, FAILED, BOUND and NOISY for isr_analysis, any
`fit_status` other than "converged" for native analyses; Historical has no
flags. On 8 April every map carries dashed lines at the eclipse's first
contact, maximum and last contact at Millstone Hill (ground: 18:15:50,
19:29:30 with 94 % obscured, 20:38:50 UTC; isr_analysis `eclipse_timing.py`).
Example: `compare.html?date=2024-04-08&before=isr-nosat-v2&after=isr-v2rfi-v2&flags=hide`.
