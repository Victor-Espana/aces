#!/usr/bin/env bash
set -e
export DEBIAN_FRONTEND=noninteractive

# 1) System update + herramientas básicas
apt-get update
apt-get install -y --no-install-recommends \
build-essential gfortran \
libblas-dev liblapack-dev \
libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
pandoc qpdf ca-certificates

# 2) Instalar R (desde los repos del sistema)
apt-get install -y --no-install-recommends r-base r-base-dev

# 3) Paquetes R básicos para desarrollo/chequeo
R -q -e "options(repos='https://cloud.r-project.org'); \
         install.packages(c('remotes','devtools','rcmdcheck','testthat','roxygen2','BiocManager'))"

# 4) Instalar las dependencias del paquete según DESCRIPTION
R -q -e "remotes::install_deps(dependencies = TRUE)"
