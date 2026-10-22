# Compute and static hosting

Develop and test the GUI directly on avaruus so changes are visible at the
real site. Do not use a local GUI preview as the development target.

Authoritative development checkout: `avaruus.uit.no:/home/j/src/guisdap26`.
Compute products: `/home/j/data/guisdap/2022-11-02-v4/frames/*.h5`.
Static export: `/home/j/src/guisdap26/web/data/2022-11-02-v4/`.

The event copy currently used is under
`/nfs/revontuli/data/oliver/02112022/beata_cp1_2.2u_AA@uhf/20221102_17`.
The archive root supplied for future discovery is `/urdr/data/beata` on revontuli;
on avaruus the same urdr mount is `/nfs/urdr/data/beata`. That archive does not
contain this event at the exact `beata_cp1_2.2u_AA@uhf/20221102_17` path on the
inspection date. Do not substitute another experiment version silently.

## Avaruus

Apache serves `http://avaruus.uit.no/guisdap/` directly from `/srv/guisdap/`.
It is enabled at boot. GUI assets and exported products are deployed there;
the catalogue and calibration APIs are proxied to local user services.
`scripts/apache-guisdap.conf` defines the alias,
disables directory listing, and sets the WebAssembly MIME type. The legacy
user service on port 8088 remains useful for development previews.

```sh
# One-time administrator setup after installing apache2:
sudo install -d -o j -g j -m 755 /srv/guisdap
sudo install -m 644 scripts/apache-guisdap.conf /etc/apache2/conf-available/guisdap.conf
sudo a2enmod headers
sudo a2enconf guisdap
sudo apache2ctl configtest
sudo systemctl reload apache2
# If UFW is active with default-deny incoming:
sudo ufw allow 80/tcp comment "GUISDAP Apache HTTP"
# Build and publish from the avaruus checkout:
scripts/build-web.sh
scripts/deploy-web.sh
```

The data browser is the root page. Experiments open `dataset.html`, with paired
ion-line and plasma-line views. Millstone exports are registered by antenna and
date in `millstone-config.json`; see [Millstone processing](millstone.md).

The processing monitor shares the catalogue API, at `/guisdap/api/catalog/live`.
It needs no additional Apache route. Configure `RUISDAP_LIVE_ROOTS` on the
catalogue service to include writable status directories on each shared compute
resource. The monitor's 12-second lease is readable by the worker; detailed arrays
are generated only during a valid lease. Progress checkpoints remain independent.

Deployment preserves `data/`, `catalog/` and existing deployment configuration.
Do not overwrite those directories with a development export.

## Thin host, e.g. juha.no

Copy the `web/` assets, including `light.css`, `scientific-plots.js`,
`plot-data.js`, `embedded.js`, `vendor/` and generated `pkg/`, into the GUI
document directory, e.g. `/var/www/guisdap/`. Preserve its deployment selectors.
Export processed data into `/mnt/shovel/guisdap/YYYY-MM-DD/` and serve that root
through a read-only static alias. Set `config.json` to:

```json
{"data_root":"/guisdap/data/2022-11-02-v4/"}
```

For Apache, adapt `scripts/apache-guisdap.conf` and add a more specific `/guisdap/data/` Alias to `/mnt/shovel/guisdap/` before the GUI alias, with a matching read-only Directory block. An nginx alternative is `scripts/nginx-thin.conf`.
No compute connection is required. If GUI and data use different origins, the
product server must allow the GUI origin through CORS. Avoid public directory
listing. HDF5 files remain the internal canonical products; the static export is
little-endian float32 arrays and a JSON manifest.

Export into a staging day directory and switch it into place after completion.
Within an export the JSON index is written last. Existing date exports should not
be rewritten while clients are reading them; publish a new versioned directory
when reprocessing, then change `config.json` atomically.

The ion-fit inspector is `ion.html`, linked from the main viewer. Its static
product URL is configured independently in `ion-config.json`. Copy the JSON
products and this configuration to the thin host along with the same WASM assets.
No raw input, experiment-definition MAT file, HDF5 reader, or backend is required
by the browser. Internal ion-fit HDF5 products remain on the compute/data host.
