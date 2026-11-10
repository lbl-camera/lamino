#!/usr/bin/env python3
"""Does TV help or hurt the recovery of Bloch / Neel walls? (thin Mumax lambda_tv sweep)

Reads the tables of wall_chirality.py --json for the central slices and the surface slices
(entries: the thin truth, then recons/Mumax/thin_tv/{clean,noisy}/<run>/recon with run =
noreg, c<tv> (const_mag + TV) or t<tv> (TV only)) and plots, against lambda_tv:

  f_B          Bloch fraction <I_t^2> / (<I_t^2> + <I_n^2>) (truth: dashed line)
  I_t / truth  integrated in-plane moment along the walls
  I_n / truth  integrated in-plane moment across the walls
  r            mean Pearson r of m_x, m_y with the truth (central slice)

Blue = const_mag + TV, orange = TV only (lambda_tv = 0 of that family: const_mag only /
noreg); solid = clean projections, dashed hollow = emma-matched noise. Writes
out_prefix.png and out_prefix.md (all numbers).

Usage: plot_thin_tv.py central.json surface.json out_prefix
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"  # validated pair (validate_palette.js, light)
TVS = ["0", "0.02", "0.05", "0.1", "0.2", "0.5"]


def parse(rows):
    """-> truth row, {(family, level): [value per TVS or None]} for each key in the row."""
    truth, runs = None, {}
    for r in rows:
        parts = Path(r["name"]).parts
        if parts[-1] != "recon":
            truth = r
            continue
        level, run = parts[-3], parts[-2]
        runs[(level, run)] = r
    return truth, runs


def series(runs, level, family, key, truth_val=None):
    ys = []
    for tv in TVS:
        name = "noreg" if (family == "t" and tv == "0") else f"{family}{tv}"
        r = runs.get((level, name))
        v = None if r is None else r[key]
        ys.append(v if v is None or truth_val is None else v / truth_val)
    return ys


def main():
    central, surface, prefix = sys.argv[1:4]
    sets = [("central slices", *parse(json.load(open(central)))),
            ("surface slices", *parse(json.load(open(surface))))]
    panels = [("f_b", "Bloch fraction f_B", None),
              ("rms_it", "I_t / truth (along the wall)", "ratio"),
              ("rms_in", "I_n / truth (across the wall)", "ratio"),
              ("r_xy", "in-plane r with truth", None)]

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                         "figure.facecolor": SURFACE, "axes.facecolor": SURFACE})
    fig, axs = plt.subplots(2, 4, figsize=(14, 6.4), sharex=True)
    lines = []
    for i, (title, truth, runs) in enumerate(sets):
        for r in runs.values():  # mean in-plane correlation
            r["r_xy"] = 0.5 * (r["r_x"] + r["r_y"])
        for j, (key, label, mode) in enumerate(panels):
            ax = axs[i, j]
            tv = truth[key] if (mode == "ratio") else None
            for family, colour in (("c", BLUE), ("t", ORANGE)):
                for level, ls, filled in (("clean", "-", True), ("noisy", "--", False)):
                    ys = series(runs, level, family, key, tv)
                    xs = [x for x, y in enumerate(ys) if y is not None]
                    ln, = ax.plot(xs, [ys[x] for x in xs], ls, color=colour, lw=2, marker="o",
                                  ms=7, mfc=colour if filled else SURFACE, mec=colour, mew=1.6,
                                  label=f"{'const_mag + TV' if family == 'c' else 'TV only'}, "
                                        f"{level}")
                    if i == 0 and j == 0:
                        lines.append(ln)
            if key == "f_b":
                ax.axhline(truth[key], color=INK2, lw=1.2, ls=":")
            elif mode == "ratio":
                ax.axhline(1.0, color=INK2, lw=1.2, ls=":")
            elif key == "r_xy":
                ax.axhline(1.0, color=INK2, lw=1.2, ls=":")
            ax.set_title(f"{label}" + (f"  ({title})" if j == 0 else ""), fontsize=10,
                         loc="left", color=INK)
            ax.grid(axis="y", color=GRID, lw=0.8)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if i == 1:
                ax.set_xticks(range(len(TVS)), TVS)
                ax.set_xlabel("lambda_tv")
    fig.legend(handles=lines, loc="upper center", ncol=4, frameon=False, fontsize=10)
    fig.text(0.995, 0.005, "dotted: truth (thin Mumax)", ha="right", va="bottom", color=INK2,
             fontsize=8)
    fig.tight_layout(rect=(0, 0.02, 1, 0.93))
    fig.savefig(f"{prefix}.png", dpi=130)

    out = []
    for title, truth, runs in sets:
        out += [f"### {title}", "",
                f"truth: f_B {truth['f_b']:.2f}, I_t {truth['rms_it']:.2f}, I_n {truth['rms_in']:.2f}",
                "", "| level | run | f_B | I_t | I_n | floor | r_x | r_y |", "|---|---|---|---|---|---|---|---|"]
        for (level, run), r in sorted(runs.items(), key=lambda kv: (kv[0][0], kv[0][1][0], float(kv[0][1][1:]) if kv[0][1] != 'noreg' else -1.0)):
            out.append(f"| {level} | {run} | {r['f_b']:.2f} | {r['rms_it']:.2f} | {r['rms_in']:.2f} | "
                       f"{r['floor']:.2f} | {r['r_x']:.2f} | {r['r_y']:.2f} |")
        out.append("")
    Path(f"{prefix}.md").write_text("\n".join(out) + "\n")
    print(f"wrote {prefix}.png, {prefix}.md")


if __name__ == "__main__":
    main()
