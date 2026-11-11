#!/usr/bin/env bash
# Download the public released N-glycan dataset used for the real-data case study and
# verify each file against the MD5 checksums published on Zenodo.
#
# Dataset: R. Smith, "Released N-glycans from IgG and Fetuin on Agilent QTOF and Orbitrap",
# Zenodo (2021), doi:10.5281/zenodo.4705546, CC BY 4.0. Data prepared and donated by the
# Biologics MS group at Charles River Labs.
#
# On Sapelo2, run this on the transfer node (xfer.gacrc.uga.edu), not a login node:
#   bash cluster/download_glycan_data.sh /scratch/$USER/glycan_data
set -euo pipefail
DEST="${1:-data/raw}"
mkdir -p "$DEST"
cd "$DEST"
BASE="https://zenodo.org/records/4705546/files"
while read -r md5 name; do
  if [[ -f "$name" ]] && echo "$md5  $name" | md5sum -c --quiet - 2>/dev/null; then
    echo "OK (already present): $name"; continue
  fi
  echo "Downloading $name"
  curl -L --fail --retry 3 -o "$name" "$BASE/$name?download=1"
  echo "$md5  $name" | md5sum -c -
done <<'LIST'
dcbdf46866f2d17fde18ff6250594d1c MASOP00220_AB_CRL_07JAN2021_Fetuin_001.d.zip
72a269831600c6818e3ec2b54fad7875 MASOP00011_AB_0332_16DEC2020_hIgG_NoDE_40ul_001.d.zip
LIST
for z in *.d.zip; do [[ -d "${z%.zip}" ]] || unzip -q "$z"; done
echo "Done. Convert the .d folders to mzML next (see cluster/README.md)."
