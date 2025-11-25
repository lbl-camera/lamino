#! /bin/bash

# set paths
SCRIPTS_DIR=$(dirname "$(realpath "$0")")

if [ ! -f "proj.tif" ]; then
    sh $SCRIPTS_DIR/forward.sh
fi

cat << EOF > bwd.json 
{
  "filename": "proj.tif",
  "angles": "proj_angles.txt",
  "thickness": 21,
  "gamma": 0.0,
  "output": "adj.tif"
}
EOF
# build 
cmake --build debug
./debug/tests/backproj bwd.json
tiffview adj.tif
