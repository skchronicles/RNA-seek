# Base image for differential gene 
# expression analysis, uses rocker
# r-ver:v4.4.3 as base image which
# uses Ubuntu/24.04 LTS (noble)
FROM rocker/r-ver:4.4.3

############### ABOUT ###############
# Dependencies of differential gene expression analysis:
# R>=4.4.0           apt-get, installs: 4.3.3
#  CRAN
#    • pander           installed via CRAN
#    • plotly           installed via CRAN
#    • ggplot2          installed via CRAN
#    • DT               installed via CRAN
#    • knitr            installed via CRAN
#    • RColorBrewer     installed via CRAN
#    • tidyverse        installed via CRAN
#    • dplyr            installed via CRAN
#    • rmarkdown        installed via CRAN
#    • stringr          installed via CRAN
#    • corrplot         installed via CRAN
#  Biconductor
#    • edgeR            installed via Bioconductor
#    • limma            installed via Bioconductor
#    • EnhancedVolcano  installed via Bioconductor
LABEL maintainer=kuhnsa@nih.gov

############### INIT ################
# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with the host's
# filesystem, i.e. /opt and /data
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
        && echo $TZ > /etc/timezone

############### SETUP ################
# This section installs system packages
# required for your project. If you need
# extra system packages add them here.
# Get system library dependencies of
# R packages by running within R:
# pak::pkg_sysreqs(c("devtools","packageX","packageY","packageZ"), sysreqs_platform = "ubuntu-24.04")
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        bzip2 \
        ca-certificates \
        curl \
        gawk \
        git \
        gzip \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libgit2-dev \
        libharfbuzz-dev \
        libicu-dev \
        libjpeg-dev \
        libpng-dev \
        libssl-dev \
        libtiff-dev \
        libwebp-dev \
        libx11-dev \
        libxml2-dev \
        locales \
        make \
        pandoc \
        # installs python/3.12.3
        python3 \
        python3-numpy \
        python3-pandas \
        python3-pip \
        python3-scipy \
        unzip \
        wget \
        zlib1g-dev \    
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Make python3 the default interpreter
# and install Python Packages 
RUN ln -sf /usr/bin/python3 /usr/bin/python

############### INSTALL ################
# Install required CRAN packages,
# use install.packages() for CRAN packages
RUN Rscript -e "use_ncpus <- max(parallel::detectCores()-2, 2); \
        cran_packages = c( \
                'BiocManager', \
                'corrplot', \
                'devtools', \
                'dplyr', \
                'DT', \
                'ggplot2', \
                'knitr', \
                'pander', \
                'plotly', \
                'RColorBrewer', \
                'rmarkdown', \
                'stringr', \
                'tidyverse' \
        ); \
        install.packages(cran_packages, repos = 'https://cloud.r-project.org/', Ncpus = use_ncpus)" 

# Install required Bioconductor packages,
RUN Rscript -e "use_ncpus <- max(parallel::detectCores()-2, 2); \
        bioc_packages = c( \
                'edgeR', \
                'EnhancedVolcano', \
                'limma' \
        ); \
        BiocManager::install(bioc_packages, Ncpus = use_ncpus)"

# Install GNU which, increases compatability
# with which distros included with Debian/Rocky
# linux installations, this version of which
# includes extra options that the ubuntu version 
# does not have. The gnu verison is backwards 
# compatiable with the ubunutu version.
RUN wget --no-check-certificate -O /opt2/which.tar.gz https://ftp.gnu.org/gnu/which/which-2.21.tar.gz \
        && tar -xvf /opt2/which.tar.gz -C /opt2/ \
        && rm /opt2/which.tar.gz \
        && cd /opt2/which-2.21/ \
        && ./configure --prefix=/opt2/which-2.21 \
        && make \
        && make install
WORKDIR /opt2

################ POST #################
# Add Dockerfile and export important 
# environment variables
ADD Dockerfile /opt2/rna-seek_dge_v0-1-0.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2/which-2.21/bin:${PATH}"
WORKDIR /data2
