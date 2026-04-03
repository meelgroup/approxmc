#!/usr/bin/env bash

set -euo pipefail

rm -rf .cmake
rm -rf lib*
rm -rf Test*
rm -rf tests*
rm -rf include
rm -rf tests
rm -rf approxmc*
rm -rf apx-src
rm -rf CM*
rm -rf cmake*
rm -rf deps
rm -rf _deps
CXX=clang++ cmake -DENABLE_TESTING=ON -DSANITIZE=ON \
    -Dcadical_DIR=../../cadical/build \
    -Dcadiback_DIR=../../cadiback \
    -Dcryptominisat5_DIR=../../cryptominisat/build \
    -Dsbva_DIR=../../sbva/build \
    -Dtreedecomp_DIR=../../treedecomp/build \
    -Darjun_DIR=../../arjun/build \
    ..
make -j$(nproc)
make test
