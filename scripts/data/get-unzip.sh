#!/usr/bin/env bash

set -euo pipefail
set -x

if [[ $# -ne 1 ]]; then
  echo "usage: $0 <id>" >&2
  exit 1
fi

id="$1"

if ! [[ "$id" =~ ^[0-9]+$ ]]; then
  echo "error: '$id' is not a number" >&2
  exit 1
fi

dirs=(out-approxmc-"${id}"-*/)
if [[ -d "${dirs[0]}" ]]; then
  echo "out-approxmc-${id}-* already present, skipping rsync"
else
  rsync -vaP "trillium:/scratch/msoos/outfiles/*${id}*" .
fi

find out-approxmc-"${id}"-*/ -type f -name '*.xz' -print0 \
  | xargs -0 -r -P "$(nproc)" -n 1 unxz

./get_data_approxmc.py --files "out-approxmc-${id}-*/*"
