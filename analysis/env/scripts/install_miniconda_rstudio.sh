#!/bin/bash

miniconda_python_version_for_r_with_conda=$1
miniconda_version_for_r_with_conda=$2

if [ "$(uname -m)" = "x86_64" ]; then \
    wget https://repo.anaconda.com/miniconda/Miniconda3-${miniconda_python_version_for_r_with_conda}_${miniconda_version_for_r_with_conda}-0-Linux-x86_64.sh -O /home/miniconda.sh; \
    elif [ "$(uname -m)" = "aarch64" ]; then \
    wget https://repo.anaconda.com/miniconda/Miniconda3-${miniconda_python_version_for_r_with_conda}_${miniconda_version_for_r_with_conda}-0-Linux-aarch64.sh -O /home/miniconda.sh; \
    fi

chmod +x /home/miniconda.sh && \
    /bin/bash /home/miniconda.sh -b -p /home/rstudio/miniconda && \
    rm -rf /home/miniconda.sh && \
    /home/rstudio/miniconda/bin/conda clean -tipsy

echo 'export PATH="/home/rstudio/miniconda/bin:${PATH}"' >> /root/.bashrc
echo "source /home/rstudio/miniconda/bin/activate" >> /root/.bashrc

echo 'export PATH="/home/rstudio/miniconda/bin:${PATH}"' >> /home/rstudio/.bashrc
echo "source /home/rstudio/miniconda/bin/activate" >> /home/rstudio/.bashrc

ENV PATH=/home/miniconda/bin:$PATH



# Add Conda channels (add in reverse priority)
/home/miniconda/bin/conda config --add channels r
/home/miniconda/bin/conda config --add channels bioconda
/home/miniconda/bin/conda config --add channels conda-forge
/home/miniconda/bin/conda config --add channels defaults
