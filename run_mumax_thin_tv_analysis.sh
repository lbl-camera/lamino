#! /bin/bash
# Analysis of run_mumax_thin_tv.sh: integrated in-plane wall test of every recon against the thin
# truth (central slices and surface slices), then the figure and table.
# Usage (inside tmux, after the sweep):  ./run_mumax_thin_tv_analysis.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
export MPLCONFIGDIR=${MPLCONFIGDIR:-$HOME/tmp/claude/mpl}
py="$HOME/.venvs/lamino/bin/python"
base=recons/Mumax/thin_tv
truth=/homes/dinesh/data/mumax/sims_thin9/m_thin9.ovf
mkdir -p $base/plots

stems=()
for level in clean noisy; do
    for d in $base/$level/*/; do
        [ -f "$d/recon2.tiff" ] && stems+=("${d%/}/recon")
    done
done
echo "${#stems[@]} reconstructions"

$py scripts/wall_chirality.py --ovf $truth "${stems[@]}" --mask-none --erode 6 \
    --json $base/plots/central.json | tee $base/plots/central.txt
$py scripts/wall_chirality.py --ovf $truth "${stems[@]}" --mask-none --erode 6 \
    --slices 0 1 7 8 --json $base/plots/surface.json | tee $base/plots/surface.txt
$py scripts/plot_thin_tv.py $base/plots/central.json $base/plots/surface.json $base/plots/thin_tv
echo ANALYSIS DONE
