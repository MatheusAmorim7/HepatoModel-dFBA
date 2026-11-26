#!/bin/bash
# Installs the system packages and Python modules the project needs (Ubuntu/Debian, WSL, Google Colab, Docker).
# - installs only what is missing; skips apt entirely when everything is already there
# - uses root/sudo only if it will not hang (root, passwordless sudo, or an interactive terminal); otherwise tells you what to install
set -e
APT_PKGS="g++ make git coinor-libclp-dev coinor-libcoinutils-dev libsbml5-dev liblapack-dev libxml2-dev libbz2-dev zlib1g-dev python3 python3-pip"
MISSING=""
for p in $APT_PKGS; do dpkg -s "$p" >/dev/null 2>&1 || MISSING="$MISSING $p"; done
if [ -z "$MISSING" ]; then
  echo "system packages: all installed"
else
  SUDO=""; CAN=0
  if [ "$(id -u)" = "0" ]; then CAN=1
  elif command -v sudo >/dev/null 2>&1 && { sudo -n true 2>/dev/null || [ -t 0 ]; }; then SUDO="sudo"; CAN=1; fi
  if [ "$CAN" = "1" ]; then
    $SUDO apt-get update -y
    DEBIAN_FRONTEND=noninteractive $SUDO apt-get install -y --no-install-recommends $MISSING
  else
    echo "Missing system packages and no root/sudo available here. Ask your administrator to install:"; echo "  $MISSING"
    echo "(or use the Docker image of this repository)"; exit 1
  fi
fi
PYBIN=${PY:-python3}
if $PYBIN -c "import numpy, scipy, pandas, cobra, matplotlib" 2>/dev/null; then echo "python modules: all installed"
else
  $PYBIN -m pip install --quiet numpy scipy pandas cobra matplotlib 2>/dev/null || $PYBIN -m pip install --quiet --break-system-packages numpy scipy pandas cobra matplotlib
fi
echo "dependencies ok"
