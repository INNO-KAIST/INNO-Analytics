#!/bin/sh
# 실측 AFM 시험 파일 내려받기 (GitHub 공개 시험 자료 — 앱에 넣지 않고 시험에만 씀)
#  pycroscopy/SciFiDatasets, AFM-SPM/AFMReader, AFM-SPM/TopoStats
# 사용: sh fetch_afm.sh <폴더>  → python3 -I gen_afm.py ../test/io/afm <폴더>
set -e; D=${1:-real/afm}; mkdir -p "$D"; cd "$D"
g(){ curl -sSfL -o "$2" "https://raw.githubusercontent.com/$1" && echo "$2"; }
B=pycroscopy/SciFiDatasets/main/data/microscopy/spm
g $B/GwiddionReader_bipm.gwy gwy_bipm.gwy; g $B/GwiddionReader_chip.gwy gwy_chip.gwy; g $B/GwiddionReader_chip.gsf gsf_chip.gsf; g $B/GwiddionReader_samples.gwy gwy_samples.gwy; g $B/afm/PTO_110_Virgin0001.gwy gwy_PTO.gwy
g $B/afm/BrukerReader_180524_MAPI-HPA_KPFM_02.0_00006_modified.spm bruker_kpfm.spm; g $B/afm/BrukerReader_Image.001 bruker_image.001; g $B/afm/BrukerReader_ForceCurve_Sapphire_TAP525.001 bruker_force.001
g AFM-SPM/AFMReader/main/tests/resources/sample_0.gwy afmr_sample0.gwy; g AFM-SPM/AFMReader/main/tests/resources/sample_0.spm afmr_sample0.spm
for f in old_bruker.002 old_bruker.004 plasmids.spm; do g AFM-SPM/TopoStats/main/tests/resources/$f ts_$f; done
