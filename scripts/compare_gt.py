"""Compare vector reconstructions (recon{0,1,2}.tiff, shape (z,y,x)) with an OVF ground truth.

Usage: python compare_gt.py gt.ovf name=path/to/basename [name=path/to/basename ...]
"""
import sys
import numpy as np
import tifffile


def read_ovf(path):
    hdr, n = {}, 0
    with open(path) as f:
        for n, line in enumerate(f):
            if "Begin: Data Text" in line:
                break
            if line.startswith("#") and ":" in line:
                k, v = line[1:].split(":", 1)
                hdr[k.strip()] = v.strip()
    nx, ny, nz = (int(hdr[k]) for k in ("xnodes", "ynodes", "znodes"))
    d = np.loadtxt(path, comments="#", skiprows=n + 1, dtype=np.float32)
    return d.reshape(nz, ny, nx, 3)  # x fastest -> (z, y, x, comp)


def metrics(r, g):
    mask = np.linalg.norm(g, axis=-1) > 0.5
    rv, gv = r[mask], g[mask]
    out = {}
    out["corr"] = [np.corrcoef(rv[:, c], gv[:, c])[0, 1] for c in range(3)]
    out["rmse"] = np.sqrt(np.mean((rv - gv) ** 2))
    s = np.sum(rv * gv) / np.sum(rv * rv)  # best global scale of recon
    out["scale"] = s
    out["rmse_scaled"] = np.sqrt(np.mean((s * rv - gv) ** 2))
    out["rel_err"] = np.linalg.norm(s * rv - gv) / np.linalg.norm(gv)
    cos = np.sum(rv * gv, -1) / (np.linalg.norm(rv, axis=-1) * np.linalg.norm(gv, axis=-1) + 1e-12)
    ang = np.degrees(np.arccos(np.clip(cos, -1, 1)))
    out["ang_mean"], out["ang_med"] = ang.mean(), np.median(ang)
    out["mag_mean"] = np.linalg.norm(rv, axis=-1).mean()
    return out


gt = read_ovf(sys.argv[1])
print(f"GT shape {gt.shape}, mean |m| = {np.linalg.norm(gt, axis=-1).mean():.3f}")
for arg in sys.argv[2:]:
    name, base = arg.split("=")
    rec = np.stack([tifffile.imread(f"{base}{i}.tiff") for i in range(3)], -1)
    if rec.shape != gt.shape:
        print(name, "shape mismatch", rec.shape); continue
    m = metrics(rec, gt)
    print(f"\n{name}: corr(mx,my,mz) = " + ", ".join(f"{c:.3f}" for c in m["corr"]))
    print(f"  RMSE {m['rmse']:.4f} | best scale {m['scale']:.3f} -> RMSE {m['rmse_scaled']:.4f}, rel err {m['rel_err']:.3f}")
    print(f"  angular err mean {m['ang_mean']:.1f} deg, median {m['ang_med']:.1f} deg; mean |m| {m['mag_mean']:.3f}")
