#!/usr/bin/env bash
# Check the per-variant scan of nearby bases that differ from the reference (getWindowsDiffRef in
# ParsingBam.cpp) on reads with mismatches, deletions and insertions, with the same alignment written
# with M and with =/X CIGAR operations. Links tests/diffref_cigar_walk_test.cpp against the object
# files of a configured build.
# Usage: tests/test_diffref_cigar_walk.sh   (after `make`)
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
build="$(mktemp -d)"
trap 'rm -rf "$build"' EXIT

objs=()
for o in Haplotag.o ParsingBam.o Util.o HaplotagProcess.o PhasingProcess.o Phasing.o PhasingGraph.o \
         MethylXgbModel.o MethylXgbFeatureExtraction.o ModCall.o ModCallParsingBam.o ModCallProcess.o; do
    [ -f "$root/$o" ] || { echo "FAIL: $o not found; run make first"; exit 1; }
    objs+=("$root/$o")
done
# the libraries the longphase-to link uses (HTSLIB_LIB from config.mk, JEMLIB from the Makefile)
libs="$(cd "$root" && make -s --no-print-directory -f Makefile -f - print-link-libs <<'MK'
print-link-libs: ; @echo $(HTSLIB_LIB) $(JEMLIB)
MK
)"

cd "$root"
"${CXX:-g++}" -std=c++11 -Wall -fopenmp -I. -Ihtslib -o "$build/diffref_cigar_walk_test" \
    tests/diffref_cigar_walk_test.cpp "${objs[@]}" $libs
"$build/diffref_cigar_walk_test"
