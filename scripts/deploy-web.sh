#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
destination="${1:-/srv/guisdap}"
# Catalogue, products, and deployment-specific selectors belong to the server.
rsync -a --chmod=D755,F644 --exclude data --exclude catalog \
  --exclude config.json --exclude ion-config.json --exclude calibration-config.json \
  web/ "$destination/"
# Initialise a new deployment without overwriting an existing selector.
for config in config.json ion-config.json calibration-config.json; do
  if [[ ! -f "$destination/$config" && -f "web/$config" ]]; then
    install -m 644 "web/$config" "$destination/$config"
  fi
done
