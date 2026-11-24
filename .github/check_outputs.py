"""CI helper: fail if any notebook has no outputs or contains an error output."""
import glob
import sys

import nbformat

bad = []
for f in sorted(glob.glob("notebooks/*.ipynb")):
    nb = nbformat.read(f, 4)
    outputs = [o for c in nb.cells if c.cell_type == "code" for o in c.get("outputs", [])]
    if not outputs:
        bad.append(f"{f}: no outputs")
    if any(o.get("output_type") == "error" for o in outputs):
        bad.append(f"{f}: error output")
if bad:
    print("\n".join(bad))
    sys.exit(1)
print("All notebooks executed with outputs and no errors.")
