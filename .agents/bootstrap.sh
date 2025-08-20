#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------
# System deps (Linux headers & toolchain)
# --------------------------------------------
# Installs only when a package manager is available and we have permissions.
# Skip by setting: AGENTS_SKIP_SYSTEM_DEPS=1
install_system_deps() {
  # Packages you asked for:
  #   libdb-dev make gcc libexpat1-dev libxml2-dev
  # Equivalent names for other distros are mapped below.

  # If explicitly skipped, return.
  if [[ "${AGENTS_SKIP_SYSTEM_DEPS:-0}" == "1" ]]; then
    echo "[bootstrap] Skipping system deps per AGENTS_SKIP_SYSTEM_DEPS=1"
    return 0
  fi

  echo "[bootstrap] Attempting to install system deps (if permitted)..."

  # Detect root/sudo capability
  _as_root() {
    if [[ "$(id -u)" -eq 0 ]]; then
      "$@"
    elif command -v sudo >/dev/null 2>&1; then
      sudo "$@"
    else
      return 126  # "permission denied" sentinel
    fi
  }

  if command -v apt-get >/dev/null 2>&1; then
    export DEBIAN_FRONTEND=noninteractive
    if _as_root apt-get update -y; then
      _as_root apt-get install -y --no-install-recommends \
        libdb-dev make gcc libexpat1-dev libxml2-dev \
        pkg-config ca-certificates curl cpanminus
      echo "[bootstrap] System deps installed via apt-get."
      return 0
    fi
  elif command -v apk >/dev/null 2>&1; then
    # Alpine equivalents
    #   libdb-dev     -> db-dev
    #   libexpat1-dev -> expat-dev
    #   libxml2-dev   -> libxml2-dev
    #   make,gcc      -> build-base includes both
    if _as_root apk add --no-cache db-dev expat-dev libxml2-dev build-base \
        pkgconfig ca-certificates curl perl-app-cpanminus; then
      echo "[bootstrap] System deps installed via apk."
      return 0
    fi
  elif command -v dnf >/dev/null 2>&1; then
    # Fedora/RHEL (dnf)
    if _as_root dnf install -y libdb-devel expat-devel libxml2-devel gcc make \
        pkgconf-pkg-config ca-certificates curl perl-App-cpanminus; then
      echo "[bootstrap] System deps installed via dnf."
      return 0
    fi
  elif command -v yum >/dev/null 2>&1; then
    # Older RHEL/CentOS (yum)
    if _as_root yum install -y libdb-devel expat-devel libxml2-devel gcc make \
        pkgconfig ca-certificates curl perl-App-cpanminus; then
      echo "[bootstrap] System deps installed via yum."
      return 0
    fi
  elif command -v zypper >/dev/null 2>&1; then
    # openSUSE
    if _as_root zypper -n install libdb-devel libexpat-devel libxml2-devel gcc make \
        pkg-config ca-certificates curl perl-App-cpanminus; then
      echo "[bootstrap] System deps installed via zypper."
      return 0
    fi
  fi

  echo "[bootstrap] Could not install system deps (no perms or unknown OS)."
  echo "[bootstrap] Continuing without system packages; builds that need headers may fail."
  return 0
}

install_system_deps || true

# --------------------------------------------
# Perl deps (contained in ./local, no sudo)
# --------------------------------------------
mkdir -p local/bin

# Ensure cpanm is available (downloaded locally if not provided by system)
if ! command -v cpanm >/dev/null 2>&1; then
  echo "[bootstrap] Downloading cpanm..."
  if command -v curl >/dev/null 2>&1; then
    curl -fsSL https://raw.githubusercontent.com/miyagawa/cpanminus/master/cpanm -o local/bin/cpanm \
      || curl -fsSL https://cpanmin.us -o local/bin/cpanm
  else
    wget -O local/bin/cpanm https://raw.githubusercontent.com/miyagawa/cpanminus/master/cpanm \
      || wget -O local/bin/cpanm https://cpanmin.us
  fi
  chmod +x local/bin/cpanm
fi

# Activate local::lib for THIS shell only (ephemeral; perfect for isolated agents)
eval "$(perl -Mlocal::lib=--deactivate-all,local)"

# Common cpanm options: use a known mirror and avoid cpanmetadb
CPANM=(cpanm --mirror https://cpan.metacpan.org --mirror-only --notest --local-lib-contained=local)

# Install build-time tools required to parse Makefile.PL
"${CPANM[@]}" File::ShareDir::Install

# 1) Install your dist deps from Makefile.PL/META (authoritative for CI/CPAN)
"${CPANM[@]}" --installdeps --with-recommends .

# 2) Install agent-only extras (tests/dev tools) from .agents/cpanfile
if [[ -f ".agents/cpanfile" ]]; then
  "${CPANM[@]}" --installdeps --cpanfile .agents/cpanfile .
fi

# 3) Build the distribution to generate shared files and fixtures
perl Makefile.PL
make
make install

echo "[bootstrap] Done. Use '.agents/with-perl-local.sh <cmd>' to run with this env."

