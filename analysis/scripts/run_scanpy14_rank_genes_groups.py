import os
import sys
import shutil
import scanpy as sc
import hdf5plugin

adata_input_path = sys.argv[1]
adata_output_path = sys.argv[2]
scanpy_clustering_algorithm = sys.argv[3] if len(sys.argv) > 3 else "leiden"
scanpy_correction_method = sys.argv[4] if len(sys.argv) > 4 else "benjamini-hochberg"

adata = sc.read_h5ad(adata_input_path)
sc.tl.rank_genes_groups(adata, scanpy_clustering_algorithm, use_raw=True, method='wilcoxon', corr_method=scanpy_correction_method)
adata.write_h5ad(adata_output_path, compression = hdf5plugin.FILTERS["zstd"])
