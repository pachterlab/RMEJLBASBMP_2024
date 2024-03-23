import os
import sys
import shutil
import scanpy as sc
import hdf5plugin

adata_input_path = sys.argv[1]
adata_output_path = sys.argv[2]
umap_seed = sys.argv[3]

adata = sc.read_h5ad(adata_input_path)
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
sc.tl.umap(adata, random_state = int(umap_seed))
adata.write_h5ad(adata_output_path, compression = hdf5plugin.FILTERS["zstd"])


# scan_n_neighbors = sys.argv[3] if len(sys.argv) > 3 else 15
# scan_num_pcs = sys.argv[4] if len(sys.argv) > 4 else None
# sc.pp.neighbors(adata, n_neighbors=int(scan_n_neighbors), n_pcs=int(scan_num_pcs))

# dir_name, file_name = os.path.split(adata_input_path)
# new_file_name = "adata_umap.h5ad"
# output_path = os.path.join(dir_name, new_file_name)
