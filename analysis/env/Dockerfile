# Accept R_VERSION and DIGEST as build arguments
ARG R_VERSION=4.3.1
ARG RVER_DIGEST
ARG CUSTOM_BASE  
# ARG CUSTOM_BASE=josephrich98/base_image:rocker_rver_4.3.1  # uncomment (and manually set version) OR provide as argument with --build-args to use my base instead

# Use custom image if specified, otherwise use rocker; and if using rocker, use DIGEST if specified, otherwise use R_VERSION
FROM ${CUSTOM_BASE:-rocker/rstudio:${RSTUDIO_DIGEST:-$R_VERSION}}

WORKDIR /usr/src/app

COPY ./scripts/install_tidyverse_rstudio.sh /tmp/

RUN chmod +x /tmp/install_tidyverse_rstudio.sh && \
	/tmp/install_tidyverse_rstudio.sh && \
	rm /tmp/install_tidyverse_rstudio.sh


### Define miniconda build arguments ###
ARG miniconda_python_version_for_r_with_conda=py310
ARG miniconda_version_for_r_with_conda=23.5.2

COPY ./scripts/install_miniconda_rstudio.sh /tmp/

RUN chmod +x /tmp/install_miniconda_rstudio.sh && \
    /tmp/install_miniconda_rstudio.sh $miniconda_python_version_for_r_with_conda $miniconda_version_for_r_with_conda && \
    rm /tmp/install_miniconda_rstudio.sh

WORKDIR /
