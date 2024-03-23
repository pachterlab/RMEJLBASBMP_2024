# Seurat vs Scanpy project

This is a project to compare package selection and version control for single-cell RNA-seq analysis.

To run the docker image (based on rocker/rstudio, with installed conda, system packages, and python/R packages):
docker run -d -p 8787:8787 -e PASSWORD=yourpassword -v /path/to/project:/workspace josephrich98/scrnaseq_packages_and_versioning:1.0

The dockerfile can be found in ./analysis/env

There are two main folders
- matrix_generation
- analysis

Each folder has its own README.md file with more information.
