#! /bin/bash
# Sweep of lambda_tv x lambda_l2 of the combined (const_mag + TV + L2) regularizer on emma/665 (GPU).
# Configs are generated from recons/emma/665/config_tv_constmag.toml (only lambda_tv, lambda_l2,
# max_iters and the output path differ) into recons/emma/665/sweep/<run>/config.toml; one queue
# per GPU, results in recons/emma/665/sweep/<run>/ with the log in recon.log.
# Usage (inside tmux):  ./run_emma665_lambda_sweep.sh [grid|conv]    (default: grid)
#   grid  lambda_tv x lambda_l2 at n = 500 (35 runs)
#   conv  corners and the current 0.05/0.1 at n = 2000 (convergence check)
# Env: GPUS="0 1 2 3 4 5 6 7"  TV_LIST, L2_LIST, ITERS override the grid.
#      MASK_T=N  limit the magnitude support to N central slices (mask_dims = [N, 399, 399] with
#      the 399 x 399 box of recons/emma/665/mask_tv.tiff, fiducial holes cut out); runs go to
#      recons/emma/665/sweep_mask_t<N>/. The emma sample is 8-9 px thick; recon_dims[0] is 21.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"

cmake --build build/release

mode="${1:-grid}"
template="recons/emma/665/config_tv_constmag.toml"
basedir="recons/emma/665/sweep"
mask_sed=()
if [ -n "${MASK_T:-}" ]; then
    basedir="recons/emma/665/sweep_mask_t${MASK_T}"
    maskfile="recons/emma/665/mask_tv_399.tiff"
    if [ ! -f "$maskfile" ]; then
        "$HOME/.venvs/lamino/bin/python" -c "
import tifffile
m = tifffile.imread('recons/emma/665/mask_tv.tiff')
tifffile.imwrite('$maskfile', m[57:456, 57:456])"
    fi
    mask_sed=(-e "/^recon_dims/a mask = \"$maskfile\"" -e "/^recon_dims/a mask_dims = [${MASK_T}, 399, 399]")
fi
gpus=(${GPUS:-0 1 2 3 4 5 6 7})

case "$mode" in
grid)
    tvs="${TV_LIST:-0 0.01 0.02 0.05 0.1 0.2 0.5}"
    l2s="${L2_LIST:-0 0.01 0.03 0.1 0.3}"
    iters="${ITERS:-500}"
    ;;
conv)
    tvs="${TV_LIST:-0 0.5 0.05}"
    l2s="${L2_LIST:-0 0.3 0.1}"
    iters="${ITERS:-2000}"
    ;;
*)
    echo "unknown mode $mode" >&2
    exit 1
    ;;
esac

mkdir -p "$basedir"
runs=()
for tv in $tvs; do
    for l2 in $l2s; do
        # conv: only the corners (0,0), (0.5,0.3), (0,0.3), (0.5,0) and (0.05,0.1)
        if [ "$mode" = conv ]; then
            case "$tv/$l2" in
            0/0 | 0.5/0.3 | 0/0.3 | 0.5/0 | 0.05/0.1) ;;
            *) continue ;;
            esac
        fi
        run="tv${tv}_l2${l2}_n${iters}"
        dir="$basedir/$run"
        mkdir -p "$dir"
        cfg="$dir/config.toml"
        sed -e "s#^filename = \"recons/.*#filename = \"$dir/recon.tiff\"#" \
            -e "s#^max_iters = .*#max_iters = ${iters}#" \
            -e "s#^lambda_tv = .*#lambda_tv = ${tv}#" \
            -e "s#^lambda_l2 = .*#lambda_l2 = ${l2}#" "${mask_sed[@]}" "$template" > "$cfg"
        if ! python3 scripts/check_config.py "$cfg" >/dev/null 2>&1; then
            echo "sanity check failed: $cfg" >&2
            exit 1
        fi
        runs+=("$run")
    done
done

echo "all ${#runs[@]} config files exist, and look well-formed"

# lambda = 0 may be rejected by the solver parse; fail early rather than after the sweep
if [ "${DRYRUN:-0}" = 1 ]; then
    printf '%s\n' "${runs[@]}"
    exit 0
fi

ngpu=${#gpus[@]}
for g in $(seq 0 $((ngpu - 1))); do
    (
        for i in $(seq "$g" "$ngpu" $((${#runs[@]} - 1))); do
            run="${runs[$i]}"
            dir="$basedir/$run"
            echo "branch: $(git branch --show-current)" >> "$dir/recon.log"
            echo "reconstructing $dir/config.toml on GPU ${gpus[$g]} ($(date))" >> "$dir/recon.log"
            rc=0
            CUDA_VISIBLE_DEVICES="${gpus[$g]}" ./build/release/recon "$dir/config.toml" \
                >> "$dir/recon.log" 2>&1 || rc=$?
            echo "done $run gpu=${gpus[$g]} exit=$rc" >> "$basedir/chains.log"
        done
    ) &
done
wait
echo "sweep finished: $(cat "$basedir/chains.log" | wc -l) runs logged in $basedir/chains.log"
