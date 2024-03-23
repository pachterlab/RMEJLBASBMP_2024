import nbformat
import os
import subprocess
import glob

new_cell_2_code = """if (!requireNamespace("reticulate", quietly = TRUE)) remotes::install_version("reticulate", version = "1.34.0", upgrade = "never")

using_colab <- reticulate::py_run_string("
try:
    import google.colab
    using_colab = True
except ImportError:
    using_colab = False
using_colab
")$using_colab

if (using_colab) {
    system("git clone https://github.com/josephrich98/scrnaseq_packages_and_versioning.git", intern = FALSE)
}"""
    
new_cell_4_code = """seurat_version_for_download <- gsub("_", ".", seurat_version)
scanpy_version_for_download <- gsub("_", ".", scanpy_version)

if (using_colab) {
    py_command <- sprintf("import subprocess; subprocess.run(['pip', 'install', 'scanpy==%s', 'python-igraph==0.10.8', 'leidenalg==0.10.1', 'anndata==0.10.2', 'hdf5plugin==4.2.0', 'kb-python==0.27.3', 'umap-learn==0.5.2', 'louvain==0.8.1', 'git+https://github.com/has2k1/scikit-misc.git@269f61e'])", scanpy_version_for_download)
    
    reticulate::py_run_string(py_command)
}

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

if (!requireNamespace("tidyverse", quietly = TRUE)) remotes::install_version("tidyverse", version = "2.0.0", upgrade = "never")
if (!requireNamespace("rmarkdown", quietly = TRUE)) remotes::install_version("rmarkdown", version = "2.25", upgrade = "never")

# if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
if (!require(pak, quietly = TRUE)) {
    install.packages("pak", repos = sprintf(
        "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
        .Platform$pkgType,
        R.Version()$os,
        R.Version()$arch
    ))
}
if (!requireNamespace("igraph", quietly = TRUE)) pak::pak("igraph/rigraph")
if (!requireNamespace("Seurat", quietly = TRUE)) remotes::install_version("Seurat", version = seurat_version_for_download, upgrade = "never")
if (!requireNamespace("Matrix", quietly = TRUE)) remotes::install_version("Matrix", version = "1.6.4", upgrade = "never")
if (!requireNamespace("patchwork", quietly = TRUE)) remotes::install_version("patchwork", version = "1.1.3", upgrade = "never")
if (!requireNamespace("eulerr", quietly = TRUE)) remotes::install_version("eulerr", version = "7.0.0", upgrade = "never")
if (!requireNamespace("scattermore", quietly = TRUE)) remotes::install_version("scattermore", version = "1.2", upgrade = "never")
if (!requireNamespace("assertthat", quietly = TRUE)) remotes::install_version("assertthat", version = "0.2.1", upgrade = "never")
if (!requireNamespace("pheatmap", quietly = TRUE)) remotes::install_version("pheatmap", version = "1.0.12", upgrade = "never")
if (!requireNamespace("ggforce", quietly = TRUE)) remotes::install_version("ggforce", version = "0.4.1", upgrade = "never")
if (!requireNamespace("ggplotify", quietly = TRUE)) remotes::install_version("ggplotify", version = "0.1.2", upgrade = "never")
if (!requireNamespace("mclust", quietly = TRUE)) remotes::install_version("mclust", version = "6.0.1", upgrade = "never")
if (!requireNamespace("ggalluvial", quietly = TRUE)) remotes::install_version("ggalluvial", version = "0.12.5", upgrade = "never")
if (!requireNamespace("UpSetR", quietly = TRUE)) remotes::install_version("UpSetR", version = "1.4.0", upgrade = "never")
if (!requireNamespace("ggpointdensity", quietly = TRUE)) remotes::install_version("ggpointdensity", version = "0.1.0", upgrade = "never")
if (!requireNamespace("dbscan", quietly = TRUE)) remotes::install_version("dbscan", version = "1.1.12", upgrade = "never")
if (!requireNamespace("presto", quietly = TRUE)) remotes::install_github("immunogenomics/presto@31dc97f", upgrade = "never")


if (!requireNamespace("BiocManager", quietly = TRUE)) remotes::install_version("BiocManager", version = "1.30.22", upgrade = "never")
bioconductor_version <- "3.18"

if (!requireNamespace("BUSpaRse", quietly = TRUE)) BiocManager::install("BUSpaRse", version = bioconductor_version, update = FALSE)
if (!requireNamespace("DropletUtils", quietly = TRUE)) BiocManager::install("DropletUtils", version = bioconductor_version, update = FALSE)
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt", version = bioconductor_version, update = FALSE)
if (!requireNamespace("bluster", quietly = TRUE)) BiocManager::install("bluster", version = bioconductor_version, update = FALSE)"""


def rmd_to_ipynb(source_dir, dest_dir):
    rmd_files = glob.glob(os.path.join(source_dir, "*.Rmd"))
    
    # Loop through each .Rmd file and convert it to a .ipynb file using notedown
    for file_path in rmd_files:
        filename = os.path.basename(file_path).replace(".Rmd", "")
        output_file_path = os.path.join(dest_dir, f"{filename}.ipynb")
        # Construct the notedown command
        command = f"notedown {file_path} --to notebook --output {output_file_path}"
        # Execute the command
        subprocess.run(command, shell=True)


def process_notebook(input_filename, output_filename):
    # Load the notebook
    with open(input_filename, 'r', encoding='utf-8') as f:
        notebook = nbformat.read(f, as_version=4)

    # Initialize a list to hold the modified and new cells
    modified_cells = []

    for i, cell in enumerate(notebook.cells):
        if cell.cell_type == 'code':
            # Check if the cell is an R cell (starts with %%R)
            if cell.source.startswith('%%R'):
                # Remove the %%R magic
                cell.source = cell.source.replace('%%R\n', '', 1)
                modified_cells.append(cell)
            else:
                # Prepare the cell source for R by escaping single quotes and encapsulating with single quotes
                escaped_source = cell.source.replace("'", "\\'")
                cell.source = f"py_run_string('{escaped_source}')"
                modified_cells.append(cell)
        else:
            modified_cells.append(cell)
        
        # Insert new cells after processing existing cells at specified positions
        if i == 1:  # After the original second cell
            new_cell_2 = nbformat.v4.new_code_cell(source=new_cell_2_code)
            modified_cells.append(new_cell_2)
        elif i == 3:  # After the original fourth cell
            new_cell_4 = nbformat.v4.new_code_cell(source=new_cell_4_code)
            modified_cells.append(new_cell_4)

    # Replace the notebook cells with the modified list
    notebook.cells = modified_cells

    # Write the modified notebook
    with open(output_filename, 'w', encoding='utf-8') as f:
        nbformat.write(notebook, f)
        
def process_all_notebooks(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".ipynb"):
            input_filename = os.path.join(directory, filename)
            output_filename = os.path.join(directory, filename)
            process_notebook(input_filename, output_filename)
            print(f"Processed {input_filename} into {output_filename}")

# Specify the directory containing the .ipynb files
if __name__ == "__main__":
    rmd_dir = "/workspace/analysis/rmd"
    ipynb_dir = "/workspace/analysis/ipynb"
    rmd_to_ipynb(rmd_dir, ipynb_dir)
    process_all_notebooks(ipynb_dir)
    print("Remember to make any manual adjustments")


