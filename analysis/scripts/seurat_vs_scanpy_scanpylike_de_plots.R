group1_name <- "Seurat"
group2_name <- "Scanpy"
output_base_path <- "/workspace/analysis/output/SC3_v3_NextGem_SI_PBMC_10K/seuratv4.3.0_vs_scanpyv1.9.5/methods_scanpy_like/input_scanpy/kb0_28_0_raw_generated_sa/cell_fraction_1_0/read_fraction_1_0"   # "/workspace/analysis/output/SC3_v3_NextGem_SI_PBMC_10K/scanpyv1.9.5/input_scan1/kb0_28_0_raw_generated_sa/cell_fraction_1_0/read_fraction_scan1_1_0_vs_scan2_0_04_seed100"
# rds_path_markers_full <- glue::glue("{output_base_path}/data_files/markers_Full.rds")  # "/workspace/analysis/output/SC3_v3_NextGem_SI_PBMC_10K/seuratv4.3.0_vs_scanpyv1.9.5/methods_scanpy_like/input_scanpy/kb0_28_0_raw_generated_sa/cell_fraction_1_0/read_fraction_1_0/data_files/markers2.rds"
# rds_path_markers_downsampled <- glue::glue("{output_base_path}/data_files/markers_Downsampled_reads.rds")


rds_path_markers_seu <- glue::glue("{output_base_path}/data_files/markers.rds")  # "/workspace/analysis/output/SC3_v3_NextGem_SI_PBMC_10K/seuratv4.3.0_vs_scanpyv1.9.5/methods_scanpy_like/input_scanpy/kb0_28_0_raw_generated_sa/cell_fraction_1_0/read_fraction_1_0/data_files/markers2.rds"
rds_path_markers_scan <- glue::glue("{output_base_path}/data_files/results_scan.rds")

markers <- readRDS(rds_path_markers_seu)
result <- readRDS(rds_path_markers_scan)

library(tidyverse)
library(Seurat)
library(ggforce)
library(ggplotify)
library(ggalluvial)
library(glue)
library(eulerr)
theme_set(theme_bw())

source("/workspace/analysis/scripts/plotting_and_stats.R")

file_paths <- list()
file_paths$logFC_scatterplot_file_path <- glue::glue("/{output_base_path}/plots/logFC_scatterplot.tiff")
file_paths$wilcoxon_scatterplot_file_path <- glue::glue("/{output_base_path}/plots/wilcoxon_scatterplot.tiff")
file_paths$logFC_scatterplot_file_path_with_legend <- glue::glue("/{output_base_path}/plots/logFC_scatterplot_with_legend.tiff")
file_paths$upset_markers_genes_only <- glue::glue("{output_base_path}/plots/upset_marker_genes_only.tiff")
file_paths$upset_markers <- glue::glue("{output_base_path}/plots/upset_markers.tiff")
file_paths$de_stats_file <- glue::glue("{output_base_path}/stats/de_stats.txt")



seu_filtered_markers <- markers %>%
    filter(p_val_adj < 0.05)

# vectorized_seu_unfiltered_markers <- unique(markers$gene)
vectorized_seu_filtered_markers <- unique(seu_filtered_markers$gene)

scan_filtered_markers <- result %>%
    filter(p_value_adj < 0.05)

vectorized_scan_filtered_markers <- unique(scan_filtered_markers$gene)

markers_euler_genes_only <- make_euler_seurat_vs_scanpy(vectorized_seu_filtered_markers, vectorized_scan_filtered_markers, comparison = "Marker Gene", save_plot = FALSE, save_stats = file_paths$de_stats_file)
markers_euler_genes_only

upset_marker_gene_only <- make_upset_seurat_vs_scanpy(vectorized_seu_filtered_markers, vectorized_scan_filtered_markers, comparison = "Marker Gene", save = file_paths$upset_markers_genes_only)

seu_markers_df <- markers %>% select(gene = gene, cluster = cluster)
scan_markers_df <- result %>% select(gene = gene, cluster = cluster)

vectorized_seu_markers <- paste(seu_markers_df$gene, seu_markers_df$cluster, sep = "-")
vectorized_scan_markers <- paste(scan_markers_df$gene, scan_markers_df$cluster, sep = "-")

markers_euler <- make_euler_seurat_vs_scanpy(vectorized_seu_markers, vectorized_scan_markers, comparison = "Marker", save_plot = FALSE, save_stats = file_paths$de_stats_file)
markers_euler

upset_markers_all <- make_upset_seurat_vs_scanpy(vectorized_seu_markers, vectorized_scan_markers, comparison = "Marker", save = file_paths$upset_markers)


markers2 <- markers |>
    inner_join(result, by = c("cluster", "gene"))

markers2 <- markers2 |>
    dplyr::rename(
        p_val_r = p_val, logFC_r = avg_log2FC, p_val_adj_r = p_val_adj,
        p_val_py = p_value, p_val_adj_py = p_value_adj,
        logFC_py = log_fc
    )

markers2 <- markers2 |>
    mutate(cluster = factor(cluster, levels = as.character(seq_len(length(unique(cluster))) - 1)))

markers2 <- markers2 |>
    group_by(cluster) |>
    mutate(rank_r = seq_along(gene))

markers2$FC_r <- 2^markers2$logFC_r
markers2$FC_py <- 2^markers2$logFC_py

markers2 <- calculate_de_stats(markers2, save = file_paths$de_stats_file)

markers2$p_val_adj_r[markers2$p_val_adj_r == 0] <- .Machine$double.xmin
markers2$p_val_adj_py[markers2$p_val_adj_py == 0] <- .Machine$double.xmin

seurat_vs_scanpy_logFC_scatterplot <- plot_scatterplot_de_logfc(markers2, ccc = markers2$CCC[1], save = file_paths$logFC_scatterplot_file_path, outliers_excluded = FALSE)
seurat_vs_scanpy_pvaladj_scatterplot <- plot_scatterplot_de_wilcoxon(markers2, save = file_paths$wilcoxon_scatterplot_file_path, outliers_excluded = FALSE)

seurat_vs_scanpy_logFC_scatterplot_with_legend <- plot_scatterplot_de_logfc(markers2, ccc = markers2$CCC[1], save = file_paths$logFC_scatterplot_file_path_with_legend, outliers_excluded = FALSE, show_legend = TRUE)

subset_markers2 <- markers2[, c("gene", "cluster", "logFC_py", "logFC_r", "p_val_adj_r", "p_val_adj_py", "logFC_difference_magnitude", "logFC_difference_signed", "pvaladj_difference_magnitude", "pvaladj_difference_signed")]

saveRDS(subset_markers2, file = glue("{output_base_path}/data_files/markers2.rds"))
