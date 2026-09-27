#!/usr/bin/env bash
# Build and run the unit tests. Needs only a C++11 compiler (no htslib).
# Usage: tests/run_unit_tests.sh   (from any directory)
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
build="$(mktemp -d)"
trap 'rm -rf "$build"' EXIT

"${CXX:-g++}" -std=c++11 -Wall -I"$root" \
    -o "$build/test_methyl_xgb_feature_extraction" \
    "$root/tests/unit/test_methyl_xgb_feature_extraction.cpp" \
    "$root/MethylXgbFeatureExtraction.cpp"
"$build/test_methyl_xgb_feature_extraction"
