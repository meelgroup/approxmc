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
emcmake cmake -DCMAKE_INSTALL_PREFIX=$EMINSTALL ..
emmake make -j$(nproc)
emmake make install
cp approxmc.wasm ../html
cp $EMINSTALL/bin/approxmc.js ../html
