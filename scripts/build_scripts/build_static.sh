#!/bin/bash

set -e

rm -rf lib* Test* tests* include tests scalmc* CM* cmake*
cmake -DSTATICCOMPILE=ON ..
make -j6
make test
