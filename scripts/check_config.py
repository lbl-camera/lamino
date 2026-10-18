#!/usr/bin/env python3
"""Check that the files referenced by [[input]] entries in a <config>.toml exist.

Checks the keys filename, angles and alignment of every [[input]] table.
filename and angles are required; alignment is optional.

Usage: check_config.py config.toml [...]
Exit status is 1 if any file is missing.
"""

import sys
import tomllib
from pathlib import Path

REQUIRED = ("filename", "angles")
OPTIONAL = ("alignment",)


def check(config):
    with open(config, "rb") as f:
        cfg = tomllib.load(f)

    nmissing = 0
    for i, entry in enumerate(cfg.get("input", [])):
        for key in REQUIRED + OPTIONAL:
            label = f"input[{i}].{key}"
            if key not in entry:
                if key in REQUIRED:
                    print(f"[MISSING] {label}: key not set")
                    nmissing += 1
                continue
            path = Path(entry[key]).expanduser()
            ok = path.exists()
            print(f"[{'OK' if ok else 'MISSING'}] {label}: {path}")
            nmissing += not ok
    return nmissing


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    total = sum(check(c) for c in sys.argv[1:])
    print(f"\n{total} missing" if total else "\nall present")
    sys.exit(1 if total else 0)


if __name__ == "__main__":
    main()
