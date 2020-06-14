#!/bin/bash

set -e

rm -rf lib* Test* tests* include tests scalmc* CM* cmake*
CXX=clang++ cmake -DENABLE_TESTING=ON -DSANITIZE=ON ..
make -j6
make test
