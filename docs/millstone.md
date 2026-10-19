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

The whole-day native antenna/power/gap controller and exact parity with all
legacy spectral averaging/rejection rules are not finished. Cancellation's
signal-loss and covariance effects still need end-to-end injections. The
recalculation queue freezes reviewed EISCAT calibration revisions; a validated
Millstone production refit recipe is not yet registered. These are not supplied
by a successful display or a short MPI pilot.
