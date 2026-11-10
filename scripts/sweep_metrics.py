#!/usr/bin/env python3
"""Collect wall-width and data-fit metrics of the emma/665 lambda sweep.

For every <sweep>/tv<a>_l2<b>_n<N>/ that holds recon{0,1,2}.tiff:
  data fit    forward-project the recon (build/release/forward, same alignment as the
              data), compare with the measured projections inside the projected support
              (compare_forward.py) -> mean Pearson r and nrmse over all angles and both
              projection sets; cached in <run>/misfit.json
  wall width  wall_width.py with one fixed reference (--ref) so that all runs are fitted at
              the same wall points; the common-points delta is the one to compare.

Also the extra stems given with --extra (e.g. the existing noreg / constmag recons) are
analysed with wall_width.py (not forward projected unless --extra-misfit).
Writes <sweep>/metrics.json: a list of dicts (run, lambda_tv, lambda_l2, n, delta_fit_common,
iqr_common, delta_g, noise_deg, r, nrmse).

Usage: sweep_metrics.py sweep_dir --ref stem --mask-dir DIR [--half-length 10]
                        [--forward-template forward/Emma665_tv_constmag_513.toml]
"""

import argparse
import json
import re
import subprocess
import sys
import tomllib
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, str(Path(__file__).parent))
import compare_forward as cf  # noqa: E402

RUN = re.compile(r"tv([\d.]+)_l2([\d.]+)_n(\d+)$")


def misfit(run, template, mask_dir, fwd_root, gpu):
    cache = run / "misfit.json"
    if cache.exists():
        return json.loads(cache.read_text())
    fdir = fwd_root / run.name
    fdir.mkdir(parents=True, exist_ok=True)
    cfg = fdir / "forward.toml"
    text = template.read_text()
    text = re.sub(r'(?m)^data_path *=.*', f'data_path = "{run.resolve()}"', text)
    for i in range(3):
        text = re.sub(rf'(?m)^component_{i + 1} *=.*', f'component_{i + 1} = "recon{i}.tiff"', text)
    text = re.sub(r'(?m)^basedir *=.*', f'basedir = "{fdir}"', text)
    cfg.write_text(text)
    env = {"CUDA_VISIBLE_DEVICES": str(gpu), "PATH": "/usr/bin:/bin"}
    subprocess.run(["build/release/forward", str(cfg)], check=True, env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    recon, fwd = tomllib.loads((run / "config.toml").read_text()), tomllib.loads(text)
    rs, es = [], []
    for inp, pr in zip(recon["input"], fwd["projections"]):
        data = tifffile.imread(inp["filename"]).astype(np.float32)
        model = tifffile.imread(fdir / pr["filename"]).astype(np.float32)
        mask = cf.load_mask(mask_dir, pr["filename"], model.shape)
        oy, ox = cf.crop_offsets(data, model, mask)
        data = data[:, oy : oy + model.shape[1], ox : ox + model.shape[2]]
        rows = cf.compare(data, model, mask)
        rs += [x[1] for x in rows]
        es += [x[3] for x in rows]
    out = {"r": float(np.mean(rs)), "nrmse": float(np.mean(es))}
    cache.write_text(json.dumps(out))
    for f in fdir.glob("gamma*.tiff"):  # the projections are large; the numbers are cached
        f.unlink()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("sweep")
    ap.add_argument("--ref", required=True, help="stem whose walls are fitted in every recon")
    ap.add_argument("--mask-dir", required=True, help="projected support masks")
    ap.add_argument("--mask", default="recons/emma/665/mask_tv.tiff")
    ap.add_argument("--half-length", type=int, default=10)
    ap.add_argument("--thickness", type=int, default=5)
    ap.add_argument("--forward-template", default="forward/Emma665_tv_constmag_513.toml")
    ap.add_argument("--fwd-root", default="/homes/dinesh/tmp/claude/sweep_fwd")
    ap.add_argument("--gpu", type=int, default=0)
    ap.add_argument("--extra", nargs="*", default=[], help="additional stems (walls only)")
    ap.add_argument("--no-misfit", action="store_true")
    args = ap.parse_args()

    sweep = Path(args.sweep)
    runs = sorted(r for r in sweep.iterdir() if RUN.match(r.name) and (r / "recon2.tiff").exists())
    if not runs:
        sys.exit(f"no finished runs in {sweep}")
    stems = [str(r / "recon") for r in runs] + args.extra
    wj = sweep / "walls.json"
    cmd = [sys.executable, str(Path(__file__).parent / "wall_width.py"), *stems,
           "--ref", args.ref, "--mask", args.mask, "--thickness", str(args.thickness),
           "--half-length", str(args.half_length), "--json", str(wj)]
    # the reference must be among the analysed stems for the pairing to include it
    if args.ref not in stems:
        cmd.insert(3, args.ref)
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    walls = {w["name"]: w for w in json.loads(wj.read_text())}

    rows = []
    for r in runs:
        m = RUN.match(r.name)
        row = dict(run=r.name, lambda_tv=float(m[1]), lambda_l2=float(m[2]), n=int(m[3]))
        w = walls[str(r / "recon")]
        row.update({k: w.get(k) for k in ("delta_fit", "delta_fit_common", "iqr_common",
                                          "delta_g", "noise_deg", "n_ok", "n_common")})
        if not args.no_misfit:
            row.update(misfit(r, Path(args.forward_template), args.mask_dir,
                              Path(args.fwd_root), args.gpu))
        rows.append(row)
        print(row, flush=True)
    for s in args.extra + [args.ref]:
        if s in walls:
            rows.append(dict(run=s, **{k: walls[s].get(k) for k in
                                       ("delta_fit", "delta_fit_common", "iqr_common",
                                        "delta_g", "noise_deg", "n_ok", "n_common")}))
    (sweep / "metrics.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
