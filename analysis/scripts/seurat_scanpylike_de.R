seu_input_path <- "/workspace/analysis/output/SC3_v3_NextGem_SI_PBMC_10K/seuratv4.3.0_vs_scanpyv1.9.5/methods_scanpy_like/input_scanpy/kb0_28_0_raw_generated_sa/cell_fraction_1_0/read_fraction_1_0/data_files/seu.rds"
markers_output_path <- "/workspace/analysis/output/SC3_v3_NextGem_SI_PBMC_10K/seuratv4.3.0_vs_scanpyv1.9.5/methods_scanpy_like/input_scanpy/kb0_28_0_raw_generated_sa/cell_fraction_1_0/read_fraction_1_0/data_files/markers.rds"

seu <- readRDS(seu_input_path)
markers <- Seurat::FindAllMarkers(seu, logfc.threshold = 0, min.pct = 0, return.thresh = 1.0001)
markers$p_val_adj <- p.adjust(markers$p_val, method = "BH")
saveRDS(markers, file = markers_output_path)
