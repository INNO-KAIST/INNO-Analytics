#!/usr/bin/env bash
# Run every XPS test suite one after another and print one summary table (name · pass · fail · skip · time).
#   test/run_all.sh [--quick] [--seeds A-B] [--html "INNO Analytics.html" [--fig DIR]] [--only REGEX] [--logs DIR]
#   --quick        QUICK=1 for the suites that support it (test_xps_casa, test_guide, test_multi, test_arxps, test_proc)
#   --seeds A-B    run the seed-dependent Monte-Carlo suites (test_multi, test_arxps) once per SEED_OFFSET A..B (rows name@sN)
#   --html FILE    also run the browser scripts test/ui_xps*.py against FILE (figures to --fig DIR, default <logs>/fig)
#   --only REGEX   run only the suites whose name matches REGEX (bash =~)
# Suites run strictly sequentially (test_shapes2 has timing checks). Exit status 1 if any suite fails, 0 otherwise.
# A real-data suite whose private data are absent counts as SKIP (exit 0, "SKIPPED (data absent)" line), not as a pass.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; MOD="$(dirname "$HERE")"
QUICK_FLAG=""; SEEDS=""; HTML=""; FIG=""; ONLY=""; LOGS=""
while [ $# -gt 0 ]; do
  case "$1" in
    --quick) QUICK_FLAG=1 ;;
    --seeds) SEEDS="$2"; shift ;;
    --html) HTML="$2"; shift ;;
    --fig) FIG="$2"; shift ;;
    --only) ONLY="$2"; shift ;;
    --logs) LOGS="$2"; shift ;;
    -h|--help) sed -n '2,10p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac; shift
done
[ -n "$LOGS" ] || LOGS="$(mktemp -d "${TMPDIR:-/tmp}/xps_run_all.XXXXXX")"
mkdir -p "$LOGS"; [ -n "$FIG" ] || FIG="$LOGS/fig"
cd "$MOD" || exit 2

ROWS=(); ANYFAIL=0
# parse a log → "pass fail skip" (pass/fail '-' when the suite prints no counts)
parse_counts() {
  local log="$1" line p="-" f="-" s=0
  line="$(grep -aE '[0-9]+ passed, [0-9]+ failed|[0-9]+ 통과 · [0-9]+ 실패|^[0-9]+/[0-9]+ passed' "$log" | tail -1)"
  if [[ "$line" =~ ([0-9]+)\ passed,\ ([0-9]+)\ failed(,\ ([0-9]+)\ skipped)? ]]; then p="${BASH_REMATCH[1]}"; f="${BASH_REMATCH[2]}"; s="${BASH_REMATCH[4]:-0}"
  elif [[ "$line" =~ ([0-9]+)\ 통과\ ·\ ([0-9]+)\ 실패(\ ·\ ([0-9]+)\ 건너뜀)? ]]; then p="${BASH_REMATCH[1]}"; f="${BASH_REMATCH[2]}"; s="${BASH_REMATCH[4]:-0}"
  elif [[ "$line" =~ ^([0-9]+)/([0-9]+)\ passed ]]; then p="${BASH_REMATCH[1]}"; f=$(( BASH_REMATCH[2] - BASH_REMATCH[1] ))
  else  # fall back to per-check lines ("  ok   …" / "  FAIL …")
    local no nf; no=$(grep -acE '^ +ok( |$)' "$log"); nf=$(grep -acE '^ +FAIL( |$)' "$log")
    if [ "$no" -gt 0 ] || [ "$nf" -gt 0 ]; then p="$no"; f="$nf"; fi
  fi
  echo "$p $f $s"
}
run_suite() {  # name, command...
  local name="$1"; shift
  if [ -n "$ONLY" ] && ! [[ "$name" =~ $ONLY ]]; then return; fi
  local log="$LOGS/${name//[^A-Za-z0-9_@.-]/_}.log" t0 t1 rc st
  printf '%-22s … ' "$name" >&2
  t0=$(date +%s.%N); "$@" >"$log" 2>&1; rc=$?; t1=$(date +%s.%N)
  local dt; dt=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.1f", b-a}')
  read -r p f s <<<"$(parse_counts "$log")"
  if grep -aq '^SKIPPED (data absent)' "$log" && [ "$p" = "-" ] && [ "$rc" -eq 0 ]; then st=SKIP; s=1
  elif [ "$rc" -ne 0 ] || { [ "$f" != "-" ] && [ "$f" -gt 0 ]; }; then st=FAIL; ANYFAIL=1
  else st=PASS; fi
  [ "$st" = FAIL ] && [ "$f" = "-" ] && f="rc=$rc"
  echo "$st (${dt} s)" >&2
  ROWS+=("$(printf '%-22s %-5s %6s %6s %5s %8s' "$name" "$st" "$p" "$f" "$s" "$dt")")
}
Q() { if [ -n "$QUICK_FLAG" ]; then QUICK=1 "$@"; else "$@"; fi; }

run_suite test_xps        node test/test_xps.js
run_suite test_xps_casa   Q node test/test_xps_casa.js
run_suite test_casa_v3    node test/test_casa_v3.js
run_suite test_shapes2    node test/test_shapes2.js
for f in test/test_real_*.js; do [ -e "$f" ] || continue; n="$(basename "$f" .js)"; run_suite "$n" node "$f"; done
run_suite test_guide      Q node test/test_guide.js
if [ -n "$SEEDS" ]; then
  a="${SEEDS%-*}"; b="${SEEDS#*-}"
  for ((o=a; o<=b; o++)); do
    run_suite "test_multi@s$o"  env SEED_OFFSET=$o ${QUICK_FLAG:+QUICK=1} node test/test_multi.js
    run_suite "test_arxps@s$o"  env SEED_OFFSET=$o ${QUICK_FLAG:+QUICK=1} node test/test_arxps.js
  done
else
  run_suite test_multi    Q node test/test_multi.js
  run_suite test_arxps    Q node test/test_arxps.js
fi
run_suite test_proc       Q node test/test_proc.js
run_suite test_image      node test/test_image.js
for f in test/test_robust_*.js; do [ -e "$f" ] || continue; n="$(basename "$f" .js)"; run_suite "$n" node "$f"; done
if [ -n "$HTML" ]; then
  for f in test/ui_xps*.py; do [ -e "$f" ] || continue; n="$(basename "$f" .py)"; mkdir -p "$FIG/$n"; run_suite "$n" python3 "$f" "$HTML" "$FIG/$n"; done
fi

echo
printf '%-22s %-5s %6s %6s %5s %8s\n' suite status pass fail skip "time s"
printf '%s\n' "------------------------------------------------------------"
NP=0; NF=0; NS=0
for r in "${ROWS[@]}"; do echo "$r"; set -- $r; case "$2" in PASS) NP=$((NP+1));; FAIL) NF=$((NF+1));; SKIP) NS=$((NS+1));; esac; done
printf '%s\n' "------------------------------------------------------------"
echo "suites: ${#ROWS[@]} — ${NP} passed, ${NF} failed, ${NS} skipped (data absent).  logs: $LOGS"
exit $ANYFAIL
