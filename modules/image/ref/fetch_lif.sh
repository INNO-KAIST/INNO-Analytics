#!/bin/sh
# 실측 Leica LIF 시험 파일 (GitHub 공개 시험 자료: nimne/readlif, bioio-devs/bioio-lif(git LFS)) — 시험에만 씀
# 사용: sh fetch_lif.sh <폴더> → python3 -I gen_lif.py <폴더> <bioformats_package.jar> <BFDump 폴더>
#  bioformats_package.jar: PyPI python-bioformats 4.1.0 휠 안 (Bio-Formats 7.3.1), BFDump: javac -cp <jar> -d <폴더> bf/BFDump.java
set -e; D=${1:-real/lif}; mkdir -p "$D"; cd "$D"
g(){ curl -sSfL -o "$2" "$1" && echo "$2"; }
R=https://raw.githubusercontent.com/nimne/readlif/master/tests/artifacts
for f in misc/LeicaLASX_wavelength-sweep_example.lif misc/new_lasx.lif xyzt-example/xyzt-example.lif xz-example/xz-example.lif; do g "$R/$f" "readlif_$(basename $f)"; done
B=https://media.githubusercontent.com/media/bioio-devs/bioio-lif/main/bioio_lif/tests/resources
for f in merged-tiles s_1_t_1_c_2_z_1 s_1_t_4_c_2_z_1 single_tile_scan tiled; do g "$B/$f.lif" "bioio_$f.lif"; done
