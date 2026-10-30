#! /bin/bash
# Analysis of the z-masked emma/665 sweep (run_emma665_lambda_sweep.sh with MASK_T=9, 7, 11, 13):
# wall width + data fit (sweep_metrics.py), integrated in-plane wall test (wall_chirality.py)
# and comparison figures against the unmasked sweep, which has the same lambdas.
# Usage (inside tmux, after the sweeps):  ./run_emma665_masked_analysis.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
export MPLCONFIGDIR=${MPLCONFIGDIR:-$HOME/tmp/claude/mpl}
py="$HOME/.venvs/lamino/bin/python"
E=recons/emma/665
ref=$E/constmag/recon_constmag_tv
maskdir=$HOME/tmp/claude/projmask   # project_mask.py output of mask_tv.tiff

unm=()
for tv in 0.02 0.05 0.1 0.2 0.5; do for l2 in 0 0.1; do unm+=($E/sweep/tv${tv}_l2${l2}_n500/recon); done; done

echo "== walls + data fit, z-masked grid (t = 9)"
$py scripts/sweep_metrics.py $E/sweep_mask_t9 --ref $ref --mask-dir "$maskdir" --gpu 7 \
    --extra $E/noreg/recon_noreg $E/constmag/recon_constmag "${unm[@]}" > /dev/null
$py scripts/plot_lambda_sweep.py $E/sweep_mask_t9 $E/sweep_mask_t9/plots

for t in 7 11 13; do
    echo "== walls + data fit, thickness check t = $t"
    $py scripts/sweep_metrics.py $E/sweep_mask_t$t --ref $ref --mask-dir "$maskdir" --gpu 7 \
        --extra $E/sweep_mask_t9/tv0.05_l20.1_n500/recon $E/sweep/tv0.05_l20.1_n500/recon > /dev/null
done

echo "== integrated in-plane wall test (first entry gives the walls)"
$py scripts/wall_chirality.py $E/sweep_mask_t9/tv0.05_l20.1_n500/recon \
    $E/sweep/tv0.05_l20.1_n500/recon $E/sweep_mask_t9/tv0.1_l20.1_n500/recon \
    $E/sweep/tv0.1_l20.1_n500/recon $E/sweep_mask_t9/tv0.2_l20.1_n500/recon \
    $E/noreg/recon_noreg $E/constmag/recon_constmag --mask $E/mask_tv.tiff --erode 6

echo "== figures"
cmp=($E/noreg/recon_noreg $E/sweep/tv0.05_l20.1_n500/recon $E/sweep_mask_t9/tv0.05_l20.1_n500/recon
     $E/sweep_mask_t9/tv0.1_l20.1_n500/recon $E/sweep_mask_t9/tv0.2_l20.1_n500/recon)
mkdir -p $E/sweep_mask_t9/plots
$py scripts/plot_sweep_yz.py $E/sweep_mask_t9/plots/yz.png "${cmp[@]}" --mask $E/mask_tv.tiff
$py scripts/plot_sweep_slices.py $E/sweep_mask_t9/plots/slices "${cmp[@]}" --mask $E/mask_tv.tiff
echo ANALYSIS DONE
