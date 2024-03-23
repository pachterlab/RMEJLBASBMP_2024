This is the repository which contains all code for the preprint The impact of package selection and versioning on single-cell RNA-seq analysis by Joseph Rich, Lambda Moses, et al.

To run the docker image (based on rocker/rstudio, with installed conda, system packages, R packages, and prebuilt conda environments):
docker run -d -p 8787:8787 -e PASSWORD=yourpassword -v /path/to/project:/workspace josephrich98/rmejlbasbmp_2024:1.0_seu5

The dockerfile can be found in ./analysis/env

There are two main folders
- analysis (analysis of Seurat vs. Scanpy, read/cell downsampling, and version control)
- matrix_generation (count matrix generation from fastq files)

Each folder has its own README.md file with more information.
