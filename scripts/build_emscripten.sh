#!/bin/bash
set -e
rm -rf lib* Test* tests* include tests scalmc* CM* cmake* approxmc* apx-src
emcmake cmake -DCMAKE_INSTALL_PREFIX=$EMINSTALL ..
emmake make -j$(nproc)
emmake make install
cp approxmc.wasm ../html
cp $EMINSTALL/bin/approxmc.js ../html
