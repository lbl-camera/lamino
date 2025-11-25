#! /bin/bash

# setup paths
SCRIPT_DIR="$(dirname "$(realpath "$0")")"

echo $SCRIPT_DIR
# create phantom
if [ ! -f "shepp.tif" ]; then
    ./debug/tests/shepp3d
fi

# create angels 
if [ ! -f "proj_angles.txt" ]; then
    . $SCRIPT_DIR/angles.sh
    write_angles_file -70 70
fi

# write input file
cat << EOF > fwd.json 
{
  "filename": "shepp.tif",
  "angles": "proj_angles.txt",
  "gamma": 0.0,
  "output": "proj.tif"
}
EOF

# build and run
cmake --build debug
./debug/tests/forward fwd.json
tiffview proj.tif
