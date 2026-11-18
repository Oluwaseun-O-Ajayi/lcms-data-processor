"""List and download files from a PRIDE Archive project, with size and checksum checks.

Run on the GACRC transfer node (xfer.gacrc.uga.edu). Works with Python 3.6 or later (the
transfer node has 3.6.8). Downloads use curl, which resumes a partial file after a dropped
connection. For long runs, start it in the background so it survives logging out:
    nohup python3 -u cluster/pride_files.py ... > download.log 2>&1 & The workspace that wrote this script
cannot reach EBI, so the API response is checked at run time: the raw JSON is always saved,
and missing fields are reported instead of guessed.

Step 1, list (writes <dest>/<accession>_files.csv and the raw JSON):
    python cluster/pride_files.py PXD078897 --dest /scratch/$USER/glycan_data/pride

Step 2, download the files you chose (exact names, or a regular expression):
    python cluster/pride_files.py PXD078897 --dest /scratch/$USER/glycan_data/pride \
        --get "Glycan_standard_01.raw"
    python cluster/pride_files.py PXD078897 --dest /scratch/$USER/glycan_data/pride \
        --pattern "(?i)glyc.*\\.raw$" --max 1
"""

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional

API = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}/files?pageSize={size}&page={page}"
PAGE_SIZE = 100


def fetch_json(url: str, retries: int = 3, wait: float = 3.0):
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as exc:  # network errors and bad JSON both retry
            if attempt == retries:
                raise SystemExit(f"Failed to fetch {url}: {exc}")
            time.sleep(wait * attempt)


def unwrap(page):
    """Return the list of file records from one API page, whatever the wrapper."""
    if isinstance(page, list):
        return page
    if isinstance(page, dict):
        for key in ("_embedded", "content", "files", "list"):
            if key in page:
                inner = page[key]
                if isinstance(inner, dict):
                    for v in inner.values():
                        if isinstance(v, list):
                            return v
                if isinstance(inner, list):
                    return inner
    return []


def https_location(rec) -> Optional[str]:
    """Pick a download URL; PRIDE FTP paths are also served over HTTPS."""
    locs = rec.get("publicFileLocations") or []
    for loc in locs:
        value = loc.get("value", "") if isinstance(loc, dict) else str(loc)
        if value.startswith("ftp://ftp.pride.ebi.ac.uk"):
            return "https://" + value[len("ftp://"):]
        if value.startswith("https://"):
            return value
    return None


def list_files(acc: str, dest: Path) -> List[Dict]:
    records, page = [], 0
    while True:
        data = fetch_json(API.format(acc=acc, size=PAGE_SIZE, page=page))
        (dest / f"{acc}_files_page{page}.json").write_text(json.dumps(data, indent=1))
        batch = unwrap(data)
        if not batch:
            break
        records.extend(batch)
        if len(batch) < PAGE_SIZE:
            break
        page += 1

    rows, missing = [], set()
    for rec in records:
        row = {
            "file_name": rec.get("fileName"),
            "category": (rec.get("fileCategory") or {}).get("value")
            if isinstance(rec.get("fileCategory"), dict) else rec.get("fileCategory"),
            "size_bytes": rec.get("fileSizeBytes"),
            "checksum": rec.get("checksum"),
            "url": https_location(rec),
        }
        missing.update(k for k, v in row.items() if v in (None, ""))
        rows.append(row)

    out = dest / f"{acc}_files.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["file_name"], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{acc}: {len(rows)} files -> {out}")
    if missing:
        print(f"WARNING: some records lack: {', '.join(sorted(missing))}. "
              f"Check the saved JSON before downloading.")
    return rows


def verify(path: Path, row: Dict) -> bool:
    ok = True
    size = row.get("size_bytes")
    if size not in (None, ""):
        if path.stat().st_size != int(size):
            print(f"  SIZE MISMATCH: {path.stat().st_size} vs {size}")
            ok = False
    checksum = (row.get("checksum") or "").lower()
    algo = {32: "md5", 40: "sha1", 64: "sha256"}.get(len(checksum))
    if algo:
        h = hashlib.new(algo)
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() != checksum:
            print(f"  {algo.upper()} MISMATCH")
            ok = False
        else:
            print(f"  {algo.upper()} OK")
    else:
        print("  no checksum published; size check only")
    return ok


def download(rows: List[Dict], dest: Path) -> int:
    failures = 0
    for row in rows:
        name, url = row["file_name"], row["url"]
        target = dest / name
        if not url:
            print(f"SKIP {name}: no download URL in the API record")
            failures += 1
            continue
        if target.exists() and verify(target, row):
            print(f"OK (already present): {name}")
            continue
        print(f"Downloading {name} ({int(row['size_bytes'] or 0) / 1e6:.1f} MB)", flush=True)
        tmp = target.with_suffix(target.suffix + ".part")
        # curl resumes a partial file (-C -); retry up to 20 times after a dropped connection
        for attempt in range(1, 21):
            rc = subprocess.call(["curl", "-L", "--fail", "-sS", "-C", "-", "-o", str(tmp), url])
            if rc == 0:
                break
            print(f"  curl exit {rc}; retry {attempt}/20 in 30 s", flush=True)
            time.sleep(30)
        else:
            print(f"  FAILED after 20 attempts: {name} (partial file kept: {tmp.name})", flush=True)
            failures += 1
            continue
        tmp.rename(target)
        if not verify(target, row):
            failures += 1
    return failures


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("accession")
    p.add_argument("--dest", type=Path, required=True)
    p.add_argument("--get", nargs="+", default=[], help="exact file names to download")
    p.add_argument("--pattern", help="regular expression on file names to download")
    p.add_argument("--max", type=int, default=None, help="download at most this many matches")
    a = p.parse_args()

    a.dest.mkdir(parents=True, exist_ok=True)
    rows = list_files(a.accession, a.dest)
    if not (a.get or a.pattern):
        return

    chosen = [r for r in rows if r["file_name"] in a.get]
    if a.pattern:
        rx = re.compile(a.pattern)
        chosen += [r for r in rows if r["file_name"] and rx.search(r["file_name"]) and r not in chosen]
    unknown = set(a.get) - {r["file_name"] for r in rows}
    if unknown:
        print(f"Not in project: {', '.join(sorted(unknown))}")
    if a.max is not None:
        chosen = chosen[: a.max]
    if not chosen:
        raise SystemExit("Nothing matched; check the names in the CSV.")
    sys.exit(1 if download(chosen, a.dest) else 0)


if __name__ == "__main__":
    main()
