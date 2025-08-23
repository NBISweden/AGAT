#!/usr/bin/env bash
# Minimal, fast bootstrap for AGAT hacking (Debian/Ubuntu).
# Installs system libs (if available), sets up local::lib, bootstraps cpm,
# installs deps WITHOUT running CPAN tests, then builds the Makefile.
# By default, it does NOT run the full test suite (see RUN_TESTS env).

set -euo pipefail

# 0) Repo root
cd "$(git rev-parse --show-toplevel 2>/dev/null || echo .)"

# 1) System libs (skip if you lack privileges — most modules are pure-Perl)
if command -v apt-get >/dev/null 2>&1; then
  if command -v sudo >/dev/null 2>&1; then SUDO=sudo; else SUDO=""; fi
  $SUDO apt-get update -y
  $SUDO apt-get install -y --no-install-recommends \
    build-essential curl ca-certificates libdb-dev libexpat1-dev libxml2-dev cpanminus
fi

# 2) Local install root (no global pollution)
mkdir -p local/bin .agents
eval "$(perl -Mlocal::lib=local)"
export PATH="$PWD/local/bin:$PATH"

# 3) Fast installer: cpm (single-file bootstrap; no cpanmetadb needed)
if ! command -v cpm >/dev/null 2>&1; then
  curl -fsSL https://raw.githubusercontent.com/skaji/cpm/main/cpm -o .agents/cpm
  chmod +x .agents/cpm
  ln -sf "$PWD/.agents/cpm" local/bin/cpm
fi

# 4) Configure-time prerequisite (avoids Makefile.PL abort)
# If you also list it under `on 'configure'` in cpanfile, this is instant.
cpm install -L local --no-test File::ShareDir::Install

# 5) Project dependencies (runtime + tests as needed)
# For maximal speed we skip CPAN's own test suites.
# Set WITH_DEVELOP=1 to also pull dev tools (Perl::Critic, Devel::Cover, etc.).
DEVEL_FLAG=()
[[ "${WITH_DEVELOP:-0}" = "1" ]] && DEVEL_FLAG+=(--with-develop)
cpm install -L local --no-test "${DEVEL_FLAG[@]}" \
  --workers "$(nproc 2>/dev/null || echo 2)" \
  --show-build-log-on-failure

# 6) Build (Makefile + compile steps as needed)
perl Makefile.PL
make -j"$(nproc 2>/dev/null || echo 2)"

# 7) Optional smoke tests (fast). Full suite is opt-in.
# if [[ "${RUN_TESTS:-0}" = "1" ]]; then
#   # Run a quick pass; keep it short by default
#   if [[ -d t/smoke ]]; then
#     prove -lr -j"$(nproc 2>/dev/null || echo 2)" t/smoke
#   else
#     prove -lr -j"$(nproc 2>/dev/null || echo 2)" t
#   fi
# fi

echo "Bootstrap complete. Use .agents/with-perl-local.sh to run tools in this env."
