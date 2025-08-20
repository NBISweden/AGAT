#!/usr/bin/env bash
set -euo pipefail
# Use local/ for this ONE command, then exit (perfect for isolated agent runs)
eval "$(perl -Mlocal::lib=--deactivate-all,local)"
exec "$@"
