#!/usr/bin/env python3
"""Plots and tables of the emma/665 lambda_tv x lambda_l2 sweep (from sweep_metrics.py).

Reads <sweep>/metrics.json and writes into outdir:
  sweep_heatmaps.png   delta_g, delta_fit, number of fitted points, data nrmse over
                       (lambda_tv, lambda_l2); lambda = 0 is the first row/column
  sweep_lines.png      delta_g against lambda_tv (one line per lambda_l2) and against
                       lambda_l2 (one line per lambda_tv), reference recons as dashed lines
  sweep_tradeoff.png   delta_g against data nrmse (colour lambda_tv, marker lambda_l2)
  sweep_table.md       all runs, and the n = 500 vs n = 2000 comparison

Usage: plot_lambda_sweep.py sweep_dir outdir
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def grid(rows, key, tvs, l2s):
    g = np.full((len(tvs), len(l2s)), np.nan)
    for r in rows:
        v = r.get(key)
        if v is not None:
            g[tvs.index(r["lambda_tv"]), l2s.index(r["lambda_l2"])] = v
    return g


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("sweep")
    ap.add_argument("outdir")
    ap.add_argument("--n", type=int, default=500, help="iteration count of the grid")
    args = ap.parse_args()

    allrows = json.loads((Path(args.sweep) / "metrics.json").read_text())
    runs = [r for r in allrows if "lambda_tv" in r]
    refs = [r for r in allrows if "lambda_tv" not in r]
    rows = [r for r in runs if r["n"] == args.n]
    tvs = sorted({r["lambda_tv"] for r in rows})
    l2s = sorted({r["lambda_l2"] for r in rows})
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    # heatmaps
    keys = [("delta_g", "delta_g (no wall selection) [px]"), ("delta_fit", "delta (tanh fit, accepted walls) [px]"),
            ("n_ok", "accepted fit points"), ("nrmse", "data nrmse")]
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    for ax, (k, title) in zip(axs.ravel(), keys):
        g = grid(rows, k, tvs, l2s)
        im = ax.imshow(g, origin="lower", cmap="viridis")
        for i in range(len(tvs)):
            for j in range(len(l2s)):
                if np.isfinite(g[i, j]):
                    ax.text(j, i, f"{g[i, j]:.2f}", ha="center", va="center", color="w", fontsize=8)
        ax.set_xticks(range(len(l2s)), [f"{v:g}" for v in l2s])
        ax.set_yticks(range(len(tvs)), [f"{v:g}" for v in tvs])
        ax.set_xlabel("lambda_l2")
        ax.set_ylabel("lambda_tv")
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(out / "sweep_heatmaps.png", dpi=130)

    # lines
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    g, e = grid(rows, "delta_g", tvs, l2s), np.zeros((len(tvs), len(l2s)))
    for j, l2 in enumerate(l2s):
        axs[0].errorbar(range(len(tvs)), g[:, j], yerr=e[:, j] / 2, marker="o", ms=3, capsize=2,
                        label=f"lambda_l2 = {l2:g}")
    axs[0].set_xticks(range(len(tvs)), [f"{v:g}" for v in tvs])
    axs[0].set_xlabel("lambda_tv")
    for i, tv in enumerate(tvs):
        axs[1].errorbar(range(len(l2s)), g[i], yerr=e[i] / 2, marker="o", ms=3, capsize=2,
                        label=f"lambda_tv = {tv:g}")
    axs[1].set_xticks(range(len(l2s)), [f"{v:g}" for v in l2s])
    axs[1].set_xlabel("lambda_l2")
    for ax in axs:
        for r, ls in zip(refs, ("--", ":", "-.")):
            if r.get("delta_g") is not None:
                ax.axhline(r["delta_g"], color="gray", ls=ls, lw=1,
                           label=Path(r["run"]).name.replace("recon_", ""))
        ax.set_ylabel("delta_g [px]")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "sweep_lines.png", dpi=130)

    # trade-off
    fig, ax = plt.subplots(figsize=(6.5, 5))
    marks = "osD^vP*X"
    for j, l2 in enumerate(l2s):
        pts = [r for r in rows if r["lambda_l2"] == l2 and "nrmse" in r]
        sc = ax.scatter([r["nrmse"] for r in pts], [r["delta_g"] for r in pts],
                        c=[np.log10(r["lambda_tv"] + 0.005) for r in pts], cmap="plasma",
                        marker=marks[j % len(marks)], s=40, label=f"lambda_l2 = {l2:g}")
    fig.colorbar(sc, ax=ax, label="log10(lambda_tv + 0.005)")
    ax.set_xlabel("data nrmse (mean over angles)")
    ax.set_ylabel("delta_g [px]")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "sweep_tradeoff.png", dpi=130)

    # table
    def f(v, p=3):
        return "" if v is None else f"{v:.{p}f}"

    lines = ["| lambda_tv | lambda_l2 | n | delta_fit [px] | n fit | delta_g [px] | noise | nrmse | r |",
             "|---|---|---|---|---|---|---|---|---|"]
    for r in sorted(runs, key=lambda r: (r["n"], r["lambda_tv"], r["lambda_l2"])):
        lines.append(f"| {r['lambda_tv']:g} | {r['lambda_l2']:g} | {r['n']} | "
                     f"{f(r.get('delta_fit'), 2)} | {r.get('n_ok', '')} | "
                     f"{f(r.get('delta_g'), 2)} | {f(r.get('noise_deg'), 1)} | "
                     f"{f(r.get('nrmse'))} | {f(r.get('r'))} |")
    lines += ["", "References (walls only):", ""]
    for r in refs:
        lines.append(f"- {r['run']}: delta_fit {f(r.get('delta_fit'), 2)} px, "
                     f"delta_g {f(r.get('delta_g'), 2)}, noise {f(r.get('noise_deg'), 1)}")
    conv = [r for r in runs if r["n"] != args.n]
    if conv:
        lines += ["", "Iteration count (n = %d vs the other n, same lambdas):" % args.n, "",
                  "| lambda_tv | lambda_l2 | n | delta_fit | delta_g | nrmse |", "|---|---|---|---|---|---|"]
        for r in conv:
            for q in [r] + [b for b in rows if (b["lambda_tv"], b["lambda_l2"]) ==
                            (r["lambda_tv"], r["lambda_l2"])]:
                lines.append(f"| {q['lambda_tv']:g} | {q['lambda_l2']:g} | {q['n']} | "
                             f"{f(q.get('delta_fit'), 2)} | {f(q.get('delta_g'), 2)} | "
                             f"{f(q.get('nrmse'))} |")
    (out / "sweep_table.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {out}/sweep_heatmaps.png, sweep_lines.png, sweep_tradeoff.png, sweep_table.md")


if __name__ == "__main__":
    main()
