#!/bin/bash

echo "rstudio ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers

apt-get update && \
apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    nano \
    libpng-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libjpeg-dev \
    libxml2-dev \
    libcairo2-dev \
    libgit2-dev \
    default-libmysqlclient-dev \
    libpq-dev \
    libsasl2-dev \
    libsqlite3-dev \
    libssh2-1-dev \
    libxtst6 \
    libtiff5-dev \
    unixodbc-dev \
    zlib1g-dev \
    libfontconfig1-dev \
    libudunits2-dev \
    libgdal-dev \
    libmagick++-dev
    
# Rscript -e "install.packages(c('remotes', 'tidyverse', 'rmarkdown', 'BiocManager'))"
# consider devtools
    
## UNTESTED - dplyr database backends
# install2.r --error --skipmissing --skipinstalled -n "$NCPUS" \
#     arrow \
#     dbplyr \
#     DBI \
#     dtplyr \
#     duckdb \
#     nycflights13 \
#     Lahman \
#     RMariaDB \
#     RPostgres \
#     RSQLite \
#     fst

# Clean up
apt-get autoremove -y && \
apt-get autoclean -y && \
rm -rf /var/lib/apt/lists/*