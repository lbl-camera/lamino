#! /bin/bash

write_angles_file() {
    t0=$1
    t1=$2
    outfile="proj_angles.txt"
    awk -v a="$t0" -v b="$t1" '
    BEGIN {
    for (x = a; x <= b;  x += 1.0)
        print x
    }' > "$outfile"
}
