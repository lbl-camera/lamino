#! /bin/bash
# lambda_tv sweep on the thin Mumax case (9 voxels of 8 nm, see run_mumax_thin.sh), clean and with
# white noise matched to emma/665 (scripts/add_noise.py), to see whether TV helps or hurts the
# recovery of Bloch / Neel walls (scripts/wall_chirality.py against the thin truth).
#   c<tv>   const_mag + TV(tv) + L2 0, magnitude support = the 9 true slices (tv = 0: const_mag only)
#   t<tv>   TV only (no const_mag), tv > 0
#   noreg   no regularizer, 60 iterations
# Output: recons/Mumax/thin_tv/{clean,noisy}/<name>/recon{0,1,2}.tiff, chains.log in recons/Mumax/thin_tv.
# Env: GPUS="0 1 ... 7"  NOISE_SCALE (multiplies the matched noise, default 1)  ITERS (default 500)
# Usage (inside tmux, after run_mumax_thin.sh has made the clean projections):  ./run_mumax_thin_tv.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"

cmake --build build/release

py="$HOME/.venvs/lamino/bin/python"
clean="/homes/dinesh/data/mumax/projs_thin9"
noisy="/homes/dinesh/data/mumax/projs_thin9_noisy"
base="recons/Mumax/thin_tv"
iters="${ITERS:-500}"
gpus=(${GPUS:-0 1 2 3 4 5 6 7})
tvs="0.02 0.05 0.1 0.2 0.5"
emma="/homes/dinesh/data/emma/665"

mkdir -p "$noisy" "$base"
for g in gamma0 gamma45; do
    seed=0; [ "$g" = gamma45 ] && seed=1
    if [ ! -f "$noisy/$g.tiff" ]; then
        $py scripts/add_noise.py "$clean/$g.tiff" "$noisy/$g.tiff" --seed $seed --scale "${NOISE_SCALE:-1}" \
            --like "$emma/Gamma0/ualigned_XMCD_665_Gamma0.tiff" \
                   "$emma/GammaNeg45/unaligned_XMCD_665_Gammaneg45.tiff" | tee -a "$noisy/noise.log"
    fi
done
cp -n "$clean/angles.txt" "$noisy/angles.txt"

inputs() { # projection directory
    cat <<TOML
[[input]]
filename = "$1/gamma0.tiff"
angles = "$1/angles.txt"
gamma = 0
beta = 0

[[input]]
filename = "$1/gamma45.tiff"
angles = "$1/angles.txt"
gamma = 45
beta = 0

TOML
}

mk() { # level(clean|noisy) name max_iters recon_dims mask_dims regularizer-block
    local level=$1 name=$2 n=$3 dims=$4 mdims=$5 reg=$6
    local dir="$base/$level/$name" pdir="$clean"
    [ "$level" = noisy ] && pdir="$noisy"
    mkdir -p "$dir"
    {
        inputs "$pdir"
        echo "[output]"
        echo "filename = \"$dir/recon.tiff\""
        echo
        echo "[recon_params]"
        echo "max_iters = $n"
        echo "tol = 1e-5"
        echo "xtol = 1e-5"
        echo "recon_dims = $dims"
        [ -n "$mdims" ] && echo "mask_dims = $mdims"
        [ -n "$reg" ] && { echo; echo "$reg"; }
    } > "$dir/config.toml"
    python3 scripts/check_config.py "$dir/config.toml" >/dev/null
    jobs_list+=("$level/$name")
}

reg() { # const_mag(true|false) lambda_tv
    cat <<TOML
[recon_params.regularizer]
method = "combined"

[recon_params.regularizer.combined]
const_mag = $1
lambda_tv = $2
lambda_l2 = 0.0
mu = 10.0
inner_iters = 5
rho = 10.0
rho_growth = 1.01
rho_max = 200.0
estimate_Ms = true
TOML
}

jobs_list=()
for level in clean noisy; do
    mk $level noreg 60 "[21, 511, 511]" "" ""
    for tv in 0 $tvs; do
        mk $level "c$tv" "$iters" "[21, 449, 449]" "[9, 319, 319]" "$(reg true $tv)"
    done
    for tv in $tvs; do
        mk $level "t$tv" "$iters" "[21, 449, 449]" "" "$(reg false $tv)"
    done
done
echo "all ${#jobs_list[@]} config files exist, and look well-formed"
[ "${DRYRUN:-0}" = 1 ] && { printf '%s\n' "${jobs_list[@]}"; exit 0; }

ngpu=${#gpus[@]}
for g in $(seq 0 $((ngpu - 1))); do
    (
        for i in $(seq "$g" "$ngpu" $((${#jobs_list[@]} - 1))); do
            job="${jobs_list[$i]}"
            dir="$base/$job"
            echo "branch: $(git branch --show-current)" >> "$dir/recon.log"
            echo "reconstructing $dir/config.toml on GPU ${gpus[$g]} ($(date))" >> "$dir/recon.log"
            rc=0
            CUDA_VISIBLE_DEVICES="${gpus[$g]}" ./build/release/recon "$dir/config.toml" \
                >> "$dir/recon.log" 2>&1 || rc=$?
            echo "done $job gpu=${gpus[$g]} exit=$rc" >> "$base/chains.log"
        done
    ) &
done
wait
echo "finished: $(wc -l < "$base/chains.log") runs logged in $base/chains.log"
