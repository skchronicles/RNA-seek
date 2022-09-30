# Dockerfile for MultiQC/1.12
# https://github.com/skchronicles/RNA-seek/issues/26
# Using Ubuntu Focal (20.04 LTS) as base image
FROM ubuntu:20.04
LABEL maintainer=kuhnsa@nih.gov

# Overview of requirements
#  • python/3.8
#    - multiqc/1.12   (pypi)
#    - argparse       (pypi)
#    - pandas         (pypi)
#    - scipy          (pypi)
#    - matplotlib     (pypi)
#    - XlsxWriter     (pypi)

# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with host filesyetem 
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

# This section installs system packages 
# required for your project. If you need 
# extra system packages add them here.
# Installs python/3.8.10
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        locales \
        python3 \
        python3-pip \
        wget \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Make python3 the default interpreter
# and install python packages
RUN ln -sf /usr/bin/python3.8 /usr/bin/python \
    && ln -sf /usr/bin/python3.8 /usr/bin/python3
RUN pip3 install --upgrade pip \
    && pip3 install multiqc==1.12 \
    && pip3 install argparse \
    && pip3 install pandas \
    && pip3 install scipy \
    && pip3 install matplotlib \
    && pip3 install XlsxWriter

# Add Dockerfile and export env variables
ADD Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2 
ENV PATH="/opt2:$PATH"
WORKDIR /data2