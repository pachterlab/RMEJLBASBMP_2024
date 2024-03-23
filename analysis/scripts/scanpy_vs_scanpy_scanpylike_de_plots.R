group1_name <- "Full"
group2_name <- "Downsampled_reads"
output_base_path <- "/workspace/analysis/output/SC3_v3_NextGem_SI_PBMC_10K/scanpyv1.9.5/input_scan1/kb0_28_0_raw_generated_sa/cell_fraction_1_0/read_fraction_scan1_1_0_vs_scan2_0_04_seed100"  # "/workspace/analysis/output/SC3_v3_NextGem_SI_PBMC_10K/scanpyv1.9.5/input_scan1/kb0_28_0_raw_generated_sa/cell_fraction_1_0/read_fraction_scan1_1_0_vs_scan2_0_04_seed100"
rds_path_markers_full <- glue::glue("{output_base_path}/data_files/markers_Full.rds")
rds_path_markers_downsampled <- glue::glue("{output_base_path}/data_files/markers_Downsampled_reads.rds")

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
file_paths$logFC_scatterplot_file_path <- glue::glue("/workspace/analysis/logFC_scatterplot.tiff")
file_paths$wilcoxon_scatterplot_file_path <- glue::glue("/workspace/analysis/wilcoxon_scatterplot.tiff")
file_paths$logFC_scatterplot_file_path_with_legend <- glue::glue("/workspace/analysis/logFC_scatterplot_with_legend.tiff")
file_paths$upset_markers_genes_only <- glue::glue("{output_base_path}/plots/upset_marker_genes_only.tiff")
file_paths$upset_markers <- glue::glue("{output_base_path}/plots/upset_markers.tiff")
file_paths$de_stats_file <- glue::glue("{output_base_path}/stats/de_stats.txt")


scan1_name <- group1_name
scan2_name <- group2_name

result1 <- readRDS(rds_path_markers_full)
result2 <- readRDS(rds_path_markers_downsampled)

scan1_filtered_markers <- result1 %>% filter(p_val_adj < 0.05)
scan2_filtered_markers <- result2 %>% filter(p_val_adj < 0.05)

vectorized_scan1_filtered_markers <- unique(scan1_filtered_markers$gene)
vectorized_scan2_filtered_markers <- unique(scan2_filtered_markers$gene)


markers_euler_genes_only <- make_euler_scanpy(vectorized_scan1_filtered_markers, vectorized_scan2_filtered_markers, comparison = "Marker Gene", group_names = scanpy_group_names, save_plot = file_paths$euler_after_qc_marker_genes_only, save_stats = file_paths$de_stats_file)
markers_euler_genes_only

upset_marker_gene_only <- make_upset_scanpy(vectorized_scan1_filtered_markers, vectorized_scan2_filtered_markers, comparison = "Marker Gene", group_names = scanpy_group_names, save = file_paths$upset_markers_genes_only)

result1_markers_df <- result1 %>% select(gene = gene, cluster = cluster)
result2_markers_df <- result2 %>% select(gene = gene, cluster = cluster)

vectorized_scan1_markers <- paste(result1_markers_df$gene, result1_markers_df$cluster, sep = "-")
vectorized_scan2_markers <- paste(result2_markers_df$gene, result2_markers_df$cluster, sep = "-")

markers_euler <- make_euler_scanpy(vectorized_scan1_markers, vectorized_scan2_markers, comparison = "Marker", group_names = scanpy_group_names, save_plot = file_paths$euler_after_qc_marker_file_path, save_stats = file_paths$de_stats_file)
markers_euler

upset_markers_all <- make_upset_scanpy(vectorized_scan1_markers, vectorized_scan2_markers, comparison = "Marker", group_names = scanpy_group_names, save = file_paths$upset_markers)


markers2 <- result1 |>
    inner_join(result2, by = c("cluster", "gene"), suffix = c(glue(".{scan1_name}"), glue(".{scan2_name}")))

markers2 <- markers2 |>
    mutate(cluster = factor(cluster, levels = as.character(seq_len(length(unique(cluster))) - 1)))

markers2 <- markers2 |>
    group_by(cluster) |>
    mutate(rank_r = seq_along(gene))


markers2[[glue("FC.{scan1_name}")]] <- 2^markers2[[glue("avg_log2FC.{scan1_name}")]]
markers2[[glue("FC.{scan2_name}")]] <- 2^markers2[[glue("avg_log2FC.{scan2_name}")]]

markers2 <- calculate_de_stats(markers2, group1_name = scan1_name, group2_name = scan2_name, save = FALSE)

# markers2 <- readRDS(rds_path)

if ((group1_name == "Seurat" && group2_name == "Scanpy") || (group2_name == "Seurat" && group1_name == "Scanpy")) {
    markers2$p_val_adj_r[markers2$p_val_adj_r == 0] <- .Machine$double.xmin
    markers2$p_val_adj_py[markers2$p_val_adj_py == 0] <- .Machine$double.xmin
} else {
    markers2[[glue("p_val_adj.{group1_name}")]][markers2[[glue("p_val_adj.{group1_name}")]] == 0] <- .Machine$double.xmin
    markers2[[glue("p_val_adj.{group2_name}")]][markers2[[glue("p_val_adj.{group2_name}")]] == 0] <- .Machine$double.xmin
}

logFC_scatterplot <- plot_scatterplot_de_logfc(markers2, group1_name = group1_name, group2_name = group2_name, ccc = markers2$CCC[1], save = file_paths$logFC_scatterplot_file_path, outliers_excluded = FALSE)
pvaladj_scatterplot <- plot_scatterplot_de_wilcoxon(markers2, group1_name = group1_name, group2_name = group2_name, save = file_paths$wilcoxon_scatterplot_file_path, outliers_excluded = FALSE)

logFC_scatterplot_with_legend <- plot_scatterplot_de_logfc(markers2, group1_name = group1_name, group2_name = group2_name, ccc = markers2$CCC[1], save = file_paths$logFC_scatterplot_file_path_with_legend, outliers_excluded = FALSE, show_legend = TRUE)
