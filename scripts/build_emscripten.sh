#!/bin/bash
set -e
rm -rf lib*
rm -rf Test*
rm -rf tests*
rm -rf include
rm -rf tests
rm -rf scalmc*
rm -rf CM*
rm -rf cmake*
rm -rf approxmc*
rm -rf apx-src
rm -rf deps
rm -rf _deps
emcmake cmake -DCMAKE_INSTALL_PREFIX=$EMINSTALL \
    -DBUILD_SHARED_LIBS=OFF \
    -Dcadical_DIR="${SOLVERS_DIR}/cadical/build" \
    -Dcadiback_DIR="${SOLVERS_DIR}/cadiback/build" \
    -Dcryptominisat5_DIR="${SOLVERS_DIR}/cryptominisat/build" \
    -Dsbva_DIR="${SOLVERS_DIR}/sbva/build" \
    -Dtreedecomp_DIR="${SOLVERS_DIR}/treedecomp/build" \
    -Darjun_DIR="${SOLVERS_DIR}/arjun/build" \
    ..
emmake make -j$(nproc)
emmake make install
cp approxmc.wasm ../html
cp $EMINSTALL/bin/approxmc.js ../html
