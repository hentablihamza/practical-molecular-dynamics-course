"""Re-download the PDB entries in data/structures/ from the RCSB PDB and record the date in data/provenance.json.

    python scripts/fetch_structures.py
"""
import datetime as dt
import json
import urllib.request
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"
ENTRIES = {"5AWL": "chignolin CLN025, X-ray", "1UAO": "chignolin, NMR ensemble", "1M17": "EGFR kinase + erlotinib, X-ray"}

meta_file = DATA / "provenance.json"
meta = json.loads(meta_file.read_text()) if meta_file.exists() else {}
for pdb_id, desc in ENTRIES.items():
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    target = DATA / "structures" / f"{pdb_id}.pdb"
    urllib.request.urlretrieve(url, target)
    meta[f"{pdb_id}.pdb"] = {"description": desc, "url": url, "retrieved": dt.date.today().isoformat()}
    print("downloaded", pdb_id, target.stat().st_size, "bytes")
meta_file.write_text(json.dumps(meta, indent=2) + "\n")
