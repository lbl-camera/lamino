#! /bin/bash
# Thin Mumax case (GPU): the 40-slice simulation (2 nm cells in z) is block-averaged into 9
# voxels of 8 nm (72 nm, like emma/665 which is 8-9 px thick), forward projected (gamma 0 and
# 45, angles -70..70) and reconstructed with 21 slices as in the emma sweep:
#   noreg          no regularizer, 60 iterations
#   ctv_zfull      const_mag + TV + L2 (0.1, 0.1), magnitude support = all 21 slices (as the
#                  emma configs: the 2D mask is applied to every slice)
#   ctv_zmatched   same, magnitude support restricted to the 9 true slices (mask_dims)
# Data: /homes/dinesh/data/mumax/sims_thin9/m_thin9.ovf, projections in .../projs_thin9;
# configs and recons in recons/Mumax/thin/. Truth for comparison: the thin OVF itself.
# Usage (inside tmux):  ./run_mumax_thin.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"

cmake --build build/release

simdir="/homes/dinesh/data/mumax/sims_thin9"
projdir="/homes/dinesh/data/mumax/projs_thin9"
base="recons/Mumax/thin"
ovf="$simdir/m_thin9.ovf"
mkdir -p "$simdir" "$projdir" "$base"

if [ ! -f "$ovf" ]; then
    "$HOME/.venvs/lamino/bin/python" scripts/make_thin_mumax.py /homes/dinesh/data/mumax/sims/m000000.ovf "$ovf"
fi

fwd="$base/forward_thin9.toml"
cat > "$fwd" <<TOML
[input]
data_path = "$simdir"
component_1 = "m_thin9.ovf"
component_2 = ""
component_3 = ""

[angles]
begin = -70
end = 71
num_projs = 141

[output]
basedir = "$projdir"
dims = [449, 449]

[[projections]]
gamma = 0
filename = "gamma0.tiff"

[[projections]]
gamma = 45
filename = "gamma45.tiff"
TOML

if [ ! -f "$projdir/gamma45.tiff" ]; then
    mkdir -p "$projdir"
    echo "forward projecting $fwd  ($(date))" | tee -a "$projdir/forward.log"
    ./build/release/forward "$fwd" 2>&1 | tee -a "$projdir/forward.log"
fi

inputs="[[input]]
filename = \"$projdir/gamma0.tiff\"
angles = \"$projdir/angles.txt\"
gamma = 0
beta = 0

[[input]]
filename = \"$projdir/gamma45.tiff\"
angles = \"$projdir/angles.txt\"
gamma = 45
beta = 0
"

regtv='[recon_params.regularizer]
method = "combined"

[recon_params.regularizer.combined]
const_mag = true
lambda_tv = 0.1
lambda_l2 = 0.1
mu = 10.0
inner_iters = 5
rho = 10.0
rho_growth = 1.01
rho_max = 200.0
estimate_Ms = true'

mk() { # name, stem, max_iters, recon_dims, mask_dims ("" for none), regularizer block
    local name=$1 stem=$2 iters=$3 dims=$4 mdims=$5 reg=$6
    mkdir -p "$base/$name"
    {
        echo "$inputs"
        echo "[output]"
        echo "filename = \"$base/$name/$stem.tiff\""
        echo
        echo "[recon_params]"
        echo "max_iters = $iters"
        echo "tol = 1e-5"
        echo "xtol = 1e-5"
        echo "recon_dims = $dims"
        [ -n "$mdims" ] && echo "mask_dims = $mdims"
        [ -n "$reg" ] && { echo; echo "$reg"; }
    } > "$base/$name/config.toml"
    python3 scripts/check_config.py "$base/$name/config.toml" >/dev/null
}

mk noreg recon 60 "[21, 511, 511]" "" ""
mk ctv_zfull recon 500 "[21, 449, 449]" "[21, 319, 319]" "$regtv"
mk ctv_zmatched recon 500 "[21, 449, 449]" "[9, 319, 319]" "$regtv"

echo "all config files exist, and look well-formed"

gpu=0
for name in noreg ctv_zfull ctv_zmatched; do
    (
        outdir="$base/$name"
        echo "branch: $(git branch --show-current)" | tee -a "$outdir/recon.log"
        echo "reconstructing $outdir/config.toml on GPU $gpu ($(date))" | tee -a "$outdir/recon.log"
        CUDA_VISIBLE_DEVICES=$gpu ./build/release/recon "$outdir/config.toml" >> "$outdir/recon.log" 2>&1
        echo "done $name exit=$?" >> "$base/chains.log"
    ) &
    gpu=$((gpu + 1))
done
wait
echo "finished: $(cat "$base/chains.log")"
