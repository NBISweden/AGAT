#!/usr/bin/env bash
# Execute any command with the repo-local Perl env active.
set -euo pipefail
cd "$(git rev-parse --show-toplevel 2>/dev/null || echo .)"
eval "$(perl -Mlocal::lib=local)"
export PATH="$PWD/local/bin:$PATH"
exec "$@"
