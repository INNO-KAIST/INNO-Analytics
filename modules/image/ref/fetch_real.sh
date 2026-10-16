#!/bin/sh
# 실측 nd2·czi 시험 파일 내려받기 (GitHub 의 공개 시험 자료 — 앱에 넣지 않고 시험에만 씀)
#  nd2: bioio-devs/bioio-nd2 (원래 tlambert03/nd2 시험 세트, git LFS), soft-matter/pims_nd2
#  czi: sebi06/czitools (Zeiss), AllenCellModeling/aicspylibczi, ZEISS/pylibczirw (git LFS)
# 사용: sh fetch_real.sh <폴더>  (Argo_sphere_water.nd2 는 약 1 GB)
set -e; D=${1:-real}; mkdir -p "$D"; cd "$D"
g(){ curl -sSfL -o "$2" "$1" && echo "$2"; }
N=https://media.githubusercontent.com/media/bioio-devs/bioio-nd2/main/bioio_nd2/tests/resources
for f in Argo_sphere_water ND2_aryeh_but3_cont200-1 ND2_dims_c2y32x32 ND2_dims_p1z5t3c2y32x32 ND2_dims_p2z5t3-2c4y32x32 ND2_dims_p4z5t3c2y32x32 ND2_dims_rgb ND2_dims_rgb_t3p2c2z3x64y64 ND2_dims_t3c2y32x32 ND2_jonas_header_test2 ND2_maxime_BF007; do g "$N/$f.nd2" "$f.nd2"; done
g https://raw.githubusercontent.com/soft-matter/pims_nd2/master/pims_nd2/cluster.nd2 pims_cluster.nd2
C=https://raw.githubusercontent.com/sebi06/czitools/main/data
for f in 100x100 1_tileregion Airyscan Al2O3_SE_020_sp CellDivision_T3_Z5_CH2_X240_Y170 CellDivision_T3_Z5_CH2_X240_Y170_comp1 CellDivision_T3_Z5_CH2_X240_Y170_comp2 CellDivision_T3_Z6_CH1_X300_Y200_DCV_ZSTD DAPI_GFP FOV7_HV110_P0500510000 LLS7_small LSM880_AiryScan S3_1Pos_2Mosaic_T2Z3_CH2_sm Tumor_HE_RGB WellD6-7_S2 newCZI_compressed nuc_small_new_red "w96_B2_S=3_B3_S=2_B4_S=1" S2_3x3_CH2 overlapp; do g "$C/$f.czi" "$f.czi"; done
A=https://raw.githubusercontent.com/AllenCellModeling/aicspylibczi/main/c_tests/resources
for f in RGB-multichannel tiles_time_z mosaic_test; do g "$A/$f.czi" "aics_$f.czi"; done
Z=https://media.githubusercontent.com/media/ZEISS/pylibczirw/main/test_data
for f in c1_bgr48 c1_bgr96float c2_gray8_gray16 c1_gray8_s2_overlapping_bounding_boxes c2_gray8_t3_z5_s2; do g "$Z/$f.czi" "zeiss_$f.czi"; done
