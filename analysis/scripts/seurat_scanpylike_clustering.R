input_path <- "/workspace/analysis/output/SC3_v3_NextGem_SI_PBMC_10K/seuratv4.3.0_vs_scanpyv1.9.5/methods_scanpy_like/input_scanpy/kb0_28_0_raw_generated_sa/cell_fraction_1_0/read_fraction_1_0/data_files/seu.rds"
output_path <- "/workspace/analysis/output/SC3_v3_NextGem_SI_PBMC_10K/seuratv4.3.0_vs_scanpyv1.9.5/methods_scanpy_like/input_scanpy/kb0_28_0_raw_generated_sa/cell_fraction_1_0/read_fraction_1_0/data_files/seu.rds"


seu <- readRDS(input_path)
seu <- Seurat::FindClusters(seu, verbose = FALSE, algorithm = 4, resolution = 1)
saveRDS(seu, file = output_path)