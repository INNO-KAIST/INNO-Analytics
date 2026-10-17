#!/bin/sh
# Drop Shape 실측 시험 그림 내려받기 — OpenDrop 예제 (jdber1/opendrop, GPL-3.0)
#  저작권·라이선스 때문에 저장소·앱에 넣지 않고 시험에만 쓴다 (test/drop_real 은 배포 제외 목록).
#  기준값(OpenDrop 문서 화면 docs/usage/images/ift_results.png · conan_results.png):
#   water_in_air001: 바늘 0.7176 mm, ρ 1000/0 kg/m³, g 9.8 → IFT 73.31 mN/m, V 11.15 mm³, A 24.05 mm², Wo 0.6609, Bo 0.2157
#   drop_on_surface: 접선 왼 120.0° · 오른 124.2°
#  OpenDrop: Huang, Gao, Berry et al., JOSS 6 (2021) 2604, doi:10.21105/joss.02604 · Berry et al., J. Colloid Interface Sci. 454 (2015) 226
# 사용: sh ref/fetch_drop.sh            → test/drop_real/*.png (md5 확인)
set -e; H=$(cd "$(dirname "$0")" && pwd); D="$H/../test/drop_real"; mkdir -p "$D"; cd "$D"
U=https://raw.githubusercontent.com/jdber1/opendrop/61b227b70a92cf2f6e97d9ae20ba7a14fdd41bfd/example_images
while read m f; do curl -sSfL -o "$f" "$U/$f"; echo "$m  $f" | md5sum -c - ; done <<'L'
ebd294697b26825aaebf84b70d3d587c drop_on_surface.png
fdf05b95012e8b67e51bedd1ac51b049 drop_on_surface_with_needle.png
38862d7170c81d4546eb00c611453a7a water_in_air001.png
f90b4d21e301d8f53e5e8a6b2fc14bdc water_in_air002.png
69b11fc2d7657fdfcfd6cec54a47d05d water_in_air003.png
da8a17c4954ba17c471377490c56aa7d water_in_air004.png
2009efad2e4dd795b2ec1dfc0f9f84eb water_in_air005.png
L
