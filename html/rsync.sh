#!/bin/bash
set -euxo pipefail

rsync -vaP approxmc.js approxmc.wasm  index.html msoos.org:/var/www/approxmc/
