#!/bin/bash

set -e

rm -rf lib* Test* tests* include tests scalmc* CM* cmake* approxmc* apx-src
cmake -DENABLE_TESTING=ON ..
make -j16
make test
