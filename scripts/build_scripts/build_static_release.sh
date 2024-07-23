#!/bin/bash

set -e

rm -rf lib* Test* tests* include tests scalmc* CM* cmake*
cmake -DCMAKE_BUILD_TYPE=Release -DSTATICCOMPILE=ON ..
make -j6
make test
