#!/usr/bin/env bash
set -e
export DEBIAN_FRONTEND=noninteractive

# --- 1) Repos y clave de CRAN para Ubuntu 24.04 (noble) ---
apt-get update
apt-get install -y --no-install-recommends ca-certificates gpg curl
install -d -m 0755 /etc/apt/keyrings
curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
  | gpg --dearmor -o /etc/apt/keyrings/cran.gpg
echo "deb [signed-by=/etc/apt/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" \
  > /etc/apt/sources.list.d/cran-r.list

# --- 2) Instalar R (mínimo) ---
apt-get update
apt-get install -y --no-install-recommends r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev

# --- 3) Paquetes R básicos para check ---
R -q -e "options(repos='https://cloud.r-project.org'); install.packages(c('rcmdcheck','devtools','testthat'))"

# --- 4) (Opcional, rápido) Instalar sólo Depends/Imports del paquete ---
# si el repo es un paquete R con DESCRIPTION en la raíz:
R -q -e "if (file.exists('DESCRIPTION')) remotes::install_deps(dependencies=c('Depends','Imports'))"

# --- 5) Limpieza para acelerar arranques posteriores ---
apt-get clean
rm -rf /var/lib/apt/lists/*
