# dpi <- 300
# dpi <- 500
axis_text_size <- 1.7
axis_numbering_size <- 1.4


ditto_colors <- c(
    "#D55E00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#E69F00", "#CC79A7", "#666666", "#AD7700", "#1C91D4", "#007756", "#D5C711", "#005685",
    "#A04700", "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71", "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C", "#FFCB57", "#9AD2F2",
    "#2CFFC6", "#F6EF8E", "#38B7FF", "#FF9B4D", "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45", "#AA9F0D", "#00446B", "#803800", "#8D3666",
    "#3D3D3D"
)

if (!exists("group1_color")) {
    group1_color <- "#D55E00"
}

if (!exists("group2_color")) {
    group2_color <- "#56B4E9"
}
seurat_v_scanpy_baseline_color <- "gray60"

seurat_group_names_default <- c(Seurat1 = "Seurat1", Seurat2 = "Seurat2")

pal <- c(Seurat = group1_color, Scanpy = group2_color, Seurat1 = group1_color, Seurat2 = group2_color, Scanpy1 = group1_color, Scanpy2 = group2_color)

output_base_path_default <- getwd()

file_paths_default <- list(
    euler_stats_before_QC_file = glue("{output_base_path_default}/stats/euler_stats.txt"),
    euler_stats_after_QC_file = glue("{output_base_path_default}/stats/euler_stats.txt"),
    pca_knn_clustering_umap_file = glue::glue("{output_base_path_default}/stats/pca_knn_clustering_umap_stats.txt"),
    de_stats_file = glue("{output_base_path_default}/stats/de_stats.txt"),
    knee_plot = glue("{output_base_path_default}/plots/knee_plot.tiff"),
    umi_scatterplot = glue("{output_base_path_default}/plots/umi_scatterplot.tiff"),
    violin_file_path = glue("{output_base_path_default}/plots/violin_plot.tiff"),
    upset_cells = glue::glue("{output_base_path_default}/plots/upset_cells.tiff"),
    upset_genes = glue::glue("{output_base_path_default}/plots/upset_genes.tiff"),
    upset_hvgs = glue::glue("{output_base_path_default}/plots/upset_hvgs.tiff"),
    upset_markers_genes_only = glue::glue("{output_base_path_default}/plots/upset_marker_genes_only.tiff"),
    upset_markers = glue::glue("{output_base_path_default}/plots/upset_markers.tiff"),
    euler_before_qc_cell_file_path = glue("{output_base_path_default}/plots/euler_cells_beforeQC.tiff"),
    euler_before_qc_gene_file_path = glue("{output_base_path_default}/plots/euler_genes_beforeQC.tiff"),
    euler_after_qc_cell_file_path = glue("{output_base_path_default}/plots/euler_cells_afterQC.tiff"),
    euler_after_qc_gene_file_path = glue("{output_base_path_default}/plots/euler_genes_afterQC.tiff"),
    euler_after_qc_hvg_file_path = glue("{output_base_path_default}/plots/euler_hvgs_afterQC.tiff"),
    euler_after_qc_marker_file_path = glue("{output_base_path_default}/plots/euler_markers.tiff"),
    euler_after_qc_marker_manual_bonferroni_file_path = glue("{output_base_path_default}/plots/euler_markers_manual_bonferroni.tiff"),
    pca_elbow_filepath_combined = glue("{output_base_path_default}/plots/pca_elbow_combined.tiff"),
    pca_12_filepath = glue("{output_base_path_default}/plots/pca_scatterplot_12.tiff"),
    pca_34_filepath = glue("{output_base_path_default}/plots/pca_scatterplot_34.tiff"),
    pca_loading_diffs = glue("{output_base_path_default}/plots/pc_loading_diffs.tiff"),
    pca_eigs_diff = glue("{output_base_path_default}/plots/pc_eig_diff.tiff"),
    pca_cluster_filepath_seu = glue("{output_base_path_default}/plots/pca_scatterplot_clusters_seu.tiff"),
    pca_cluster_filepath_scan = glue("{output_base_path_default}/plots/pca_scatterplot_clusters_scan.tiff"),
    combined_pc_variance_loadings_plot = glue::glue("{output_base_path_default}/plots/combined_pc_variance_loadings_plot.tiff"),
    jaccards = glue("{output_base_path_default}/plots/jaccards.tiff"),
    knn_scatterplot = glue("{output_base_path_default}/plots/knn_scatterplot.tiff"),
    jaccard_degree_scatterplot = glue::glue("{output_base_path_default}/plots/jaccard_degree_scatterplot.tiff"),
    pheatmap = glue("{output_base_path_default}/plots/cluster_pheatmap.tiff"),
    alluvial = glue("{output_base_path_default}/plots/cluster_alluvial.tiff"),
    umap_seu = glue("{output_base_path_default}/plots/umap_seu.tiff"),
    umap_scan = glue("{output_base_path_default}/plots/umap_scan.tiff"),
    umap_centroid_distances = glue("{output_base_path_default}/plots/umap_centroid)distances.tiff"),
    umap_jaccard_knn_density = glue::glue("{output_base_path}/plots/umap_jaccard_knn_density.tiff"),
    logFC_histogram_magnitude_file_path = glue("{output_base_path_default}/plots/logFC_histogram_magnitude.tiff"),
    logFC_histogram_signed_file_path = glue("{output_base_path_default}/plots/logFC_histogram_signed.tiff"),
    wilcoxon_histogram_magnitude_file_path = glue("{output_base_path_default}/plots/wilcoxon_histogram_magnitude.tiff"),
    wilcoxon_histogram_signed_file_path = glue("{output_base_path_default}/plots/wilcoxon_histogram_signed.tiff"),
    logFC_scatterplot_file_path = glue("{output_base_path_default}/plots/logFC_scatterplot.tiff"),
    wilcoxon_scatterplot_file_path = glue("{output_base_path_default}/plots/wilcoxon_scatterplot.tiff"),
    logFC_scatterplot_outliers_removed_file_path = glue("{output_base_path_default}/plots/logFC_scatterplot_no_outliers.tiff"),
    wilcoxon_scatterplot_outliers_removed_file_path = glue("{output_base_path_default}/plots/wilcoxon_scatterplot_no_outliers.tiff"),
    logFC_boxplot_magnitude_file_path = glue("{output_base_path_default}/plots/logFC_boxplot_magnitude.tiff"),
    logFC_boxplot_signed_file_path = glue("{output_base_path_default}/plots/logFC_boxplot_signed.tiff"),
    wilcoxon_boxplot_magnitude_file_path = glue("{output_base_path_default}/plots/wilcoxon_boxplot_magnitude.tiff"),
    wilcoxon_boxplot_signed_file_path = glue("{output_base_path_default}/plots/wilcoxon_boxplot_signed.tiff")
)


make_save_path <- function(filepath = NULL, default_filepath = NULL) {
    if (filepath == TRUE) {
        filepath <- default_filepath
    } else if (is.character(filepath) && (substr(filepath, 1, 1) != "/")) {
        if (!is.character(filepath)) {
            output_base_path_default <- getwd()
        }
        filepath <- file.path(output_base_path_default, filepath)
    }

    if (!dir.exists(dirname(filepath))) {
        dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
    }

    return(filepath)
}


make_knee_plot <- function(bc_rank, save = FALSE) {
    options(repr.plot.width = 9, repr.plot.height = 6)

    knee_plt <- tibble(
        rank = bc_rank[["rank"]],
        total = bc_rank[["total"]]
    ) %>%
        distinct() %>%
        dplyr::filter(total > 0)
    annot <- tibble(
        inflection = metadata(bc_rank)[["inflection"]],
        rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]])
    )
    p <- ggplot(knee_plt, aes(rank, total)) +
        geom_line() +
        geom_hline(aes(yintercept = inflection), data = annot, linetype = 2) +
        geom_vline(aes(xintercept = rank_cutoff), data = annot, linetype = 2) +
        scale_x_log10() +
        scale_y_log10() +
        annotation_logticks() +
        labs(x = "Rank", y = "Total UMIs")

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$knee_plot)
        ggsave(filepath, plot = p, dpi = dpi)
    }

    return(p)
}


get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}


make_umi_scatterplot <- function(res_mat1, res_mat2, UMI_cutoff1 = NULL, UMI_cutoff2 = NULL, point_density = FALSE, res_mat1_name = "kb", res_mat2_name = "cellranger", color_points = FALSE, save = FALSE) {
    common_cells <- intersect(colnames(res_mat1), colnames(res_mat2))

    umi_counts_res_mat <- colSums(res_mat1[, common_cells])
    umi_counts_res_mat2 <- colSums(res_mat2[, common_cells])

    data_for_plot <- data.frame(
        Cell = common_cells,
        UMI_Counts_res_mat1 = umi_counts_res_mat,
        UMI_Counts_res_mat2 = umi_counts_res_mat2
    )

    minor_breaks <- rep(1:9, 21) * (10^rep(-10:10, each = 9))

    if (!color_points) {
        p <- ggplot(data_for_plot, aes(x = UMI_Counts_res_mat1, y = UMI_Counts_res_mat2))
    } else {
        res_mat1_filtered <- res_mat1[, tot_counts > UMI_cutoff]
        res_mat1_filtered <- res_mat1_filtered[Matrix::rowSums(res_mat1_filtered) > 0, ]

        res_mat2_filtered <- res_mat2[, tot_counts > UMI_cutoff]
        res_mat2_filtered <- res_mat2_filtered[Matrix::rowSums(res_mat2_filtered) > 0, ]

        res_mat1_common <- res_mat1_filtered[, colnames(res_mat1_filtered) %in% common_cells]
        res_mat2_common <- res_mat2_filtered[, colnames(res_mat2_filtered) %in% common_cells]

        data_for_plot$in_res_mat1_filtered <- data_for_plot$Cell %in% colnames(res_mat1_common)
        data_for_plot$in_res_mat2_filtered <- data_for_plot$Cell %in% colnames(res_mat2_common)

        data_for_plot$color_group <- with(data_for_plot, factor(paste(in_res_mat1_filtered, in_res_mat2_filtered)))

        p <- ggplot(data_for_plot, aes(x = UMI_Counts_res_mat1, y = UMI_Counts_res_mat2, color = color_group)) +
            scale_color_manual(values = c("TRUE TRUE" = "orange", "TRUE FALSE" = "green", "FALSE TRUE" = "blue", "FALSE FALSE" = "brown"))
    }

    if (point_density) {
        data_for_plot$density <- get_density(data_for_plot$UMI_Counts_res_mat1, data_for_plot$UMI_Counts_res_mat2, n = 100)

        p <- p +
            geom_point(aes(color = data_for_plot$density), alpha = 0.2, size = 0.5) +
            # ggpointdensity::geom_pointdensity(size = 0.5, alpha = 1, adjust = 0.2) +
            scico::scale_color_scico(palette = "grayC", direction = -1, end = 0.8)
    } else {
        p <- p +
            geom_point(alpha = 0.2, size = 0.5)
    }

    p <- p +
        scale_x_log10(
            minor_breaks = minor_breaks,
            labels = function(x) format(x, scientific = TRUE)
        ) +
        scale_y_log10(
            minor_breaks = minor_breaks,
            labels = function(x) format(x, scientific = TRUE)
        ) +
        annotation_logticks() +
        geom_abline(slope = 1, intercept = 0, show.legend = FALSE, linewidth = 0.3, color = "gray30", linetype = 2) +
        theme(
            text = element_text(family = "Arial"),
            legend.position = "none",
            plot.margin = margin(l = 5, r = 15, t = 4),
            axis.text = element_text(size = rel(axis_numbering_size)), # Increase axis tick labels size
            axis.title = element_text(size = rel(axis_text_size)) # Increase axis titles size
        ) +
        labs(
            x = glue("UMI Counts in {res_mat1_name}"),
            y = glue("UMI Counts in {res_mat2_name}")
        )

    if (!is.null(UMI_cutoff1)) {
        p <- p +
            geom_vline(aes(xintercept = UMI_cutoff1), color = "gray30", linetype = 2) +
            annotate("text", x = UMI_cutoff1 * 0.9, y = max(data_for_plot$UMI_Counts_res_mat2), label = glue("{res_mat1_name} UMI cutoff"), hjust = 1, vjust = 0.88, color = "gray30", size = 5)
    }

    if (!is.null(UMI_cutoff2)) {
        p <- p +
            geom_hline(aes(yintercept = UMI_cutoff2), linetype = 2) +
            annotate("text", x = 0, y = UMI_cutoff2 * 1.13, label = glue("{res_mat2_name} UMI cutoff"), hjust = -0.05, vjust = 0, color = "gray30", size = 5)
    }

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$umi_scatterplot)
        ggsave(filepath, plot = p, dpi = dpi, width = 2100, height = 2100, units = "px")
    }

    return(p)
}


make_violin_plot <- function(seu, show_points = FALSE, color = NULL, save = FALSE) {
    pt.size <- ifelse(show_points == TRUE, 0.01, 0)
    p1 <- VlnPlot(seu, features = "nFeature_RNA", group.by = "orig.ident", pt.size = pt.size, cols = color) +
        coord_cartesian(ylim = c(0, 8250)) +
        scale_y_continuous(breaks = seq(0, 8500, by = 2000)) +
        theme(text = element_text(family = "Arial"), legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

    p2 <- VlnPlot(seu, features = "nCount_RNA", group.by = "orig.ident", pt.size = pt.size, cols = color) +
        coord_cartesian(ylim = c(0, 64000)) +
        scale_y_continuous(breaks = seq(0, 60000, by = 20000)) +
        theme(text = element_text(family = "Arial"), legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

    p3 <- VlnPlot(seu, features = "pct_mt", group.by = "orig.ident", pt.size = pt.size, cols = color) +
        coord_cartesian(ylim = c(0, 80)) +
        scale_y_continuous(breaks = seq(0, 80, by = 10)) +
        theme(text = element_text(family = "Arial"), legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

    # Arrange the plots into one plot with 3 columns
    combined_plot <- (p1 | p2 | p3) + plot_layout(ncol = 3)

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$violin_file_path)
        ggsave(filepath, plot = combined_plot, dpi = dpi)
    }

    return(combined_plot)
}


upset_plot_general <- function(data, group1_name, group2_name, comparison, before_filtering = FALSE, as_ggplot = FALSE, save = FALSE, default_plotpath) {
    all_elements <- unique(c(data[[group1_name]], data[[group2_name]]))

    df <- setNames(
        data.frame(
            element = all_elements,
            all_elements %in% data[[group1_name]],
            all_elements %in% data[[group2_name]]
        ),
        c("element", group1_name, group2_name)
    )


    df <- df %>%
        mutate(!!sym(group1_name) := as.integer(!!sym(group1_name)),
               !!sym(group2_name) := as.integer(!!sym(group2_name)))

    if (before_filtering) {
        y_label <- glue("{comparison} Intersection (before UMI filtering)")
    } else {
        y_label <- glue("{comparison} Intersection")
    }
    
    if (comparison == "HVG") {
        set_size_scale_max_value <- 1.55  # 1.35
        intersection_title_size <-  3.3
    } else if (comparison == "Marker") {
        set_size_scale_max_value <- 1.95  # 1.65
        intersection_title_size <-  4
    } else {
        set_size_scale_max_value <- 1.65  # 1.4
        intersection_title_size <-  4
    }

    p <- UpSetR::upset(df[, -1], set_size.show = TRUE, set_size.scale_max = set_size_scale_max_value * max(length(data[[group1_name]]), length(data[[group2_name]])), set_size.numbers_size = 10.5, mainbar.y.label = y_label, sets.x.label = glue("{comparison}s"), text.scale = c(intersection_title_size, 3, 2.25, 1.3, 3.7, 4), empty.intersections = TRUE)  # text.scale = c(2.4, 1.4, 1.85, 1, 1.85, 1.95)

    if (as_ggplot) {
        p <- as.ggplot(p)
        if (save == TRUE || is.character(save)) {
            filepath <- make_save_path(filepath = save, default_filepath = default_plotpath)
            ggsave(filepath, plot = p, dpi = dpi, bg = "white", width = 2300, height = 2100, units = "px")
        }
    } else {
        if (save == TRUE || is.character(save)) {
            filepath <- make_save_path(filepath = save, default_filepath = default_plotpath)
            tiff(filepath, width = 2300, height = 2100, res = dpi, bg = "white", units = "px")
            print(p)
            dev.off()
        }
    }

    return(p)
}


make_upset_seurat_vs_scanpy <- function(seu, adata, comparison, as_ggplot = FALSE, save = FALSE, default_plotpath) {
    if (comparison == "Gene") {
        data <- list(
            Seurat = rownames(seu@assays$RNA$counts),
            Scanpy = unlist(adata$var_names$tolist())
        )
        default_plotpath <- file_paths_default$upset_cells
    } else if (comparison == "Cell") {
        data <- list(
            Seurat = colnames(seu@assays$RNA$counts),
            Scanpy = unlist(adata$obs_names$tolist())
        )
        default_plotpath <- file_paths_default$upset_gene
    } else if (comparison == "HVG") {
        py_run_string("scanpy_highly_variable_genes = adata.var.index[adata.var['highly_variable']]")
        scanpy_highly_variable_genes <- unlist(py$scanpy_highly_variable_genes$tolist())
        data <- list(
            Seurat = VariableFeatures(seu),
            Scanpy = scanpy_highly_variable_genes
        )
        default_plotpath <- file_paths_default$upset_hvg
    } else if (comparison == "Marker" || comparison == "Marker Gene") {
        data <- list(
            Seurat = seu,
            Scanpy = adata
        )
        default_plotpath <- default_plotpath <- file_paths_default$upset_markers
    }

    p <- upset_plot_general(data, group1_name = "Seurat", group2_name = "Scanpy", comparison = comparison, before_filtering = FALSE, as_ggplot = as_ggplot, save = save, default_plotpath = default_plotpath)

    print(p)

    return(p)
}


make_euler_seurat_vs_scanpy <- function(seu, adata, comparison, before_QC = FALSE, save_plot = FALSE, save_stats = FALSE) {
    if (comparison == "Gene") {
        data <- list(
            Seurat = rownames(seu@assays$RNA$counts),
            Scanpy = unlist(adata$var_names$tolist())
        )
        default_plotpath <- ifelse(before_QC, file_paths_default$euler_before_qc_gene_file_path, file_paths_default$euler_after_qc_gene_file_path)
    } else if (comparison == "Cell") {
        data <- list(
            Seurat = colnames(seu@assays$RNA$counts),
            Scanpy = unlist(adata$obs_names$tolist())
        )
        default_plotpath <- ifelse(before_QC, file_paths_default$euler_before_qc_cell_file_path, file_paths_default$euler_after_qc_cell_file_path)
    } else if (comparison == "HVG") {
        py_run_string("scanpy_highly_variable_genes = adata.var.index[adata.var['highly_variable']]")
        scanpy_highly_variable_genes <- unlist(py$scanpy_highly_variable_genes$tolist())
        data <- list(
            Seurat = VariableFeatures(seu),
            Scanpy = scanpy_highly_variable_genes
        )
        default_plotpath <- file_paths_default$euler_after_qc_hvg_file_path
    } else if (comparison == "Marker" || comparison == "Marker Gene") {
        data <- list(
            Seurat = seu,
            Scanpy = adata
        )
        default_plotpath <- file_paths_default$euler_after_qc_marker_file_path
    }

    euler_data <- euler(data)
    seu_unique <- euler_data$original.values["Seurat"]
    scan_unique <- euler_data$original.values["Scanpy"]
    total_overlap <- euler_data$original.values["Seurat&Scanpy"]


    if (save_stats == TRUE || is.character(save_stats)) {
        default_statpath <- ifelse(before_QC, file_paths_default$euler_stats_before_QC_file, file_paths_default$euler_stats_after_QC_file)
        filepath <- make_save_path(filepath = save_stats, default_filepath = default_statpath)
        sink(filepath, append = TRUE, split = TRUE)
    }

    print(glue("Total {comparison}s in Seurat: {seu_unique+total_overlap}"))
    print(glue("Total {comparison}s in Scanpy: {scan_unique+total_overlap}"))
    print(glue("Total {comparison}s overlapping: {total_overlap}"))
    print(glue("Fraction of {comparison}s in Seurat overlapping in Scanpy: {total_overlap/(total_overlap + scan_unique)}"))
    print(glue("Fraction of {comparison}s in Scanpy overlapping in Seurat: {total_overlap/(total_overlap + seu_unique)}"))
    print(glue("{comparison}s Jaccard: {total_overlap/(total_overlap + seu_unique + scan_unique)}"))

    if (save_stats == TRUE || is.character(save_stats)) {
        sink()
    }

    euler_plot <- plot(euler_data, fills = pal, quantities = list(cex = 2), labels = list(cex = 2))

    euler_plot_ggplot_compatible <- as.ggplot(euler_plot)

    if (before_QC) {
        title <- glue("{comparison}s (before QC)")
    } else {
        title <- glue("{comparison}s")
    }

    euler_plot_ggplot_compatible <- euler_plot_ggplot_compatible +
        ggtitle(glue("{comparison}s")) +
        theme(plot.title = element_text(size = 40, hjust = 0.5))

    if (save_plot == TRUE || is.character(save_plot)) {
        filepath <- make_save_path(filepath = save_plot, default_filepath = default_plotpath)
        ggsave(filepath, plot = euler_plot_ggplot_compatible, dpi = dpi, bg = "white", width = 2300, height = 2100, units = "px")
    }

    return(euler_plot_ggplot_compatible)
}


plot_differences_histogram_seurat_vs_scanpy <- function(df, column, title, median_or_variance = NULL, x_label = "Differences", save = FALSE) {
    column_sym <- sym(column)

    p <- ggplot(df, aes(x = !!column_sym)) +
        geom_histogram(aes(y = after_stat(density)), binwidth = 0.1, fill = seurat_v_scanpy_baseline_color, color = seurat_v_scanpy_baseline_color, alpha = 0.7) +
        ggtitle(title) +
        xlab(x_label) +
        ylab("Density") +
        theme(plot.title = element_text(hjust = 0.5))

    if (column == "FC_difference_magnitude") {
        p <- p + xlim(0, 10)
    } else if (column == "FC_difference_signed") {
        p <- p + xlim(-10, 2.5)
    }

    if (grepl("FC", column)) {
        if (column == "logFC_difference_magnitude") {
            p <- p + xlim(0, 10)
            inset_xlim <- list(min = 25, max = 30)
            inset_bin_width <- 0.1
            inset_plot_coordinates <- list(xmin = 7.5, xmax = 10.4, ymin = 0.3, ymax = 0.6)
            x_scale_breaks <- waiver()
        } else if (column == "logFC_difference_signed") {
            p <- p + xlim(-10, 10)
            inset_xlim <- list(min = 25, max = 30)
            inset_bin_width <- 0.1
            inset_plot_coordinates <- list(xmin = 5, xmax = 10.8, ymin = 0.1, ymax = 0.3)
            x_scale_breaks <- waiver()
        } else if (column == "FC_difference_signed") {
            p <- p + xlim(-5, 5)
            inset_xlim <- list(min = -65, max = -5)
            inset_bin_width <- 10
            inset_plot_coordinates <- list(xmin = -5, xmax = -2, ymin = 0.2, ymax = 0.5)
            x_scale_breaks <- c(-65, -45, -25, -5)
        }

        df_inset <- df %>% filter((!!column_sym >= inset_xlim$min) & (!!column_sym <= inset_xlim$max))
        p_inset <- ggplot(df_inset, aes(x = !!column_sym)) +
            geom_histogram(aes(y = after_stat(density)), binwidth = inset_bin_width, fill = seurat_v_scanpy_baseline_color, color = seurat_v_scanpy_baseline_color, alpha = 0.7) +
            xlim(inset_xlim$min, inset_xlim$max) +
            xlab(x_label) + # Set x-axis label based on variable x_label
            scale_x_continuous(breaks = x_scale_breaks) +
            theme(axis.title.x = element_text(size = 8))

        # Create a grob from the inset plot
        inset_grob <- ggplotGrob(p_inset)

        p <- p +
            annotation_custom(grob = inset_grob, xmin = inset_plot_coordinates$xmin, xmax = inset_plot_coordinates$xmax, ymin = inset_plot_coordinates$ymin, ymax = inset_plot_coordinates$ymax)
    }

    if (column == "logFC_difference_magnitude") {
        default_plotpath <- file_paths_default$logFC_histogram_magnitude_file_path
    } else if (column == "logFC_difference_signed") {
        default_plotpath <- file_paths_default$logFC_histogram_signed_file_path
    } else if (column == "pvaladj_difference_magnitude") {
        default_plotpath <- file_paths_default$wilcoxon_histogram_magnitude_file_path
    } else if (column == "pvaladj_difference_signed") {
        default_plotpath <- file_paths_default$wilcoxon_histogram_signed_file_path
    } else {
        default_plotpath <- NULL
    }

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = default_plotpath)
        ggsave(filepath, plot = p, dpi = dpi)
    }

    return(p)
}


plot_differences_boxplot_seurat_vs_scanpy <- function(df, column, title, x_label = "Differences", save = FALSE) {
    column_sym <- sym(column)

    p <- ggplot(markers2, aes(x = "", y = !!column_sym)) +
        geom_boxplot() +
        xlab("") +
        ylab(x_label) +
        ggtitle(title)

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = default_plotpath)
        ggsave(filepath, plot = p, dpi = dpi)
    }

    if (column == "logFC_difference_magnitude") {
        default_plotpath <- file_paths_default$logFC_boxplot_magnitude_file_path
    } else if (column == "logFC_difference_signed") {
        default_plotpath <- file_paths_default$logFC_boxplot_signed_file_path
    } else if (column == "pvaladj_difference_magnitude") {
        default_plotpath <- file_paths_default$wilcoxon_boxplot_magnitude_file_path
    } else if (column == "pvaladj_difference_signed") {
        default_plotpath <- file_paths_default$wilcoxon_boxplot_signed_file_path
    } else {
        default_plotpath <- NULL
    }

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = default_plotpath)
        ggsave(filepath, plot = p, dpi = dpi)
    }

    return(p)
}


scatterplot_naming <- function(group1_name, group2_name) {
    if (group1_name == "Seurat") {
        logFC_group1 <- "logFC_r"
        p_val_adj_group1 <- "p_val_adj_r"
    } else if (group1_name == "Scanpy") {
        logFC_group1 <- "logFC_py"
        p_val_adj_group1 <- "p_val_adj_py"
    } else {
        logFC_group1 <- glue("avg_log2FC.{group1_name}")
        p_val_adj_group1 <- glue("p_val_adj.{group1_name}")
    }

    if (group2_name == "Seurat") {
        logFC_group2 <- "logFC_r"
        p_val_adj_group2 <- "p_val_adj_r"
    } else if (group2_name == "Scanpy") {
        logFC_group2 <- "logFC_py"
        p_val_adj_group2 <- "p_val_adj_py"
    } else {
        logFC_group2 <- glue("avg_log2FC.{group2_name}")
        p_val_adj_group2 <- glue("p_val_adj.{group2_name}")
    }

    return(list(logFC_group1 = logFC_group1, logFC_group2 = logFC_group2, p_val_adj_group1 = p_val_adj_group1, p_val_adj_group2 = p_val_adj_group2))
}


plot_scatterplot_de_wilcoxon <- function(markers2, metric, outliers_excluded = FALSE, show_legend = FALSE, group1_name = "Seurat", group2_name = "Scanpy", spearman = NULL, save = FALSE) {
    scatterplot_names <- scatterplot_naming(group1_name, group2_name)

    p_val_adj_group1 <- scatterplot_names$p_val_adj_group1
    p_val_adj_group2 <- scatterplot_names$p_val_adj_group2

    markers2$log_p_val_adj_1 <- -log(markers2[[p_val_adj_group1]], base = 10)
    markers2$log_p_val_adj_2 <- -log(markers2[[p_val_adj_group2]], base = 10)

    markers2$log_p_val_adj_1 <- ifelse(is.infinite(markers2$log_p_val_adj_1), 330, markers2$log_p_val_adj_1)
    markers2$log_p_val_adj_2 <- ifelse(is.infinite(markers2$log_p_val_adj_2), 330, markers2$log_p_val_adj_2)

    max_log_p <- max(
        max(markers2$log_p_val_adj_2, na.rm = TRUE),
        max(markers2$log_p_val_adj_1, na.rm = TRUE)
    )

    if (outliers_excluded) {
        plot_title <- glue("{group1_name} vs. {group2_name} Adjusted p-value (excluding outliers)")
    } else {
        plot_title <- glue("{group1_name} vs. {group2_name} Adjusted p-value")
    }

    p <- ggplot(markers2, aes(log_p_val_adj_2, log_p_val_adj_1)) +
        ggpointdensity::geom_pointdensity(size = 0.3, alpha = 1, adjust = 0.2, show.legend = FALSE) +
        scico::scale_color_scico(palette = "grayC", direction = -1, end = 0.8) +
        coord_equal() +
        coord_fixed(ratio = 1, xlim = c(0, max_log_p), ylim = c(0, max_log_p)) +
        labs(
            x = bquote(-log[10](.(group2_name) ~ "Adjusted p-value")),
            y = bquote(-log[10](.(group1_name) ~ "Adjusted p-value")),
            title = plot_title
        ) +
        theme(
            text = element_text(family = "Arial"), 
            plot.title = element_blank(),          # Increase plot title size
            axis.text = element_text(size = rel(axis_numbering_size)), # Increase axis tick labels size
            axis.title = element_text(size = rel(axis_text_size))         # Increase axis text size
        )
    
    # show spearman correlation
    if (!is.null(spearman)) {
        p <- p +
            annotate(
                "text",
                x = Inf,  # Places the annotation on the far right
                y = -Inf, # Places the annotation on the bottom
                label = glue::glue("Spearman Correlation = {sprintf('%.3f', spearman)}"),
                hjust = 1.04, # Align text to the right
                vjust = -0.28,  # Align text to the bottom
                size = 3.5
            )
        
        ymin_inset <- -5
        ymax_inset <- 145
        
    } else {
        ymin_inset <- -10
        ymax_inset <- 140
    }

    df_inset <- markers2 %>% filter((log_p_val_adj_1 <= 10) & (log_p_val_adj_2 <= 10))

    p_inset <- ggplot(df_inset, aes(log_p_val_adj_2, log_p_val_adj_1)) +
        ggpointdensity::geom_pointdensity(size = 0.3, alpha = 1, adjust = 0.2, show.legend = FALSE) +
        scico::scale_color_scico(palette = "grayC", direction = -1, end = 0.8) +
        xlim(0, 10) +
        ylim(0, 10) +
        geom_hline(yintercept = -log(0.05, base = 10), linetype = "dashed", color = "black") + # Horizontal line at y = 0.05
        geom_vline(xintercept = -log(0.05, base = 10), linetype = "dashed", color = "black") +
        theme(
            text = element_text(family = "Arial"), 
            axis.title.x = element_blank(), # Remove x-axis label
            axis.title.y = element_blank(),
            panel.background = element_rect(fill = "transparent", colour = NA), # Transparent background
            plot.background = element_rect(fill = "transparent", colour = NA), # Transparent plot background
            axis.text = element_text(size = 10)
        ) # Adjust text size for axis labels

    # Create a grob from the inset plot
    inset_grob <- ggplotGrob(p_inset)

    p <- p +
        annotation_custom(grob = inset_grob, xmin = 170, xmax = 320, ymin = ymin_inset, ymax = ymax_inset) # 170, 320   # 320, 470   # ymin_inset and ymax_inset adjusted based on whether or not spearman shows

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$wilcoxon_scatterplot_file_path)
        ggsave(filepath, plot = p, dpi = dpi, width = 2100, height = 2100, units = "px")
    }

    return(p)
}


plot_scatterplot_de_logfc <- function(markers2, outliers_excluded = FALSE, show_legend = FALSE, group1_name = "Seurat", group2_name = "Scanpy", ccc = NULL, save = FALSE) {
    scatterplot_names <- scatterplot_naming(group1_name, group2_name)
    logFC_group1 <- scatterplot_names$logFC_group1
    logFC_group2 <- scatterplot_names$logFC_group2
    
    data_pca <- data.frame(x = markers2[[logFC_group2]], 
                           y = markers2[[logFC_group1]])
    
    pca_result <- prcomp(data_pca, center = TRUE, scale. = FALSE)
    
    data_mean <- colMeans(data_pca)
    
    pc1_direction <- pca_result$rotation[,1]
    
    point1 <- data_mean - pc1_direction
    point2 <- data_mean + pc1_direction
    
    m <- (point2[2] - point1[2]) / (point2[1] - point1[1])
    
    b <- data_mean[2] - m * data_mean[1]

    if (!outliers_excluded) {
        markers2 |>
            filter(abs(.data[[logFC_group1]]) < 15) |>
            filter(abs(.data[[logFC_group2]]) < 15) -> markers2_filtered

        data_pca_filtered <- data.frame(x = markers2_filtered[[logFC_group2]], 
                               y = markers2_filtered[[logFC_group1]])
        
        pca_result_filtered <- prcomp(data_pca_filtered, center = TRUE, scale. = FALSE)
        
        data_mean_filtered <- colMeans(data_pca_filtered)
        
        pc1_direction_filtered <- pca_result_filtered$rotation[,1]
        
        point1_filtered <- data_mean_filtered - pc1_direction_filtered
        point2_filtered <- data_mean_filtered + pc1_direction_filtered
        
        m_filtered <- (point2_filtered[2] - point1_filtered[2]) / (point2_filtered[1] - point1_filtered[1])
        
        b_filtered <- data_mean_filtered[2] - m_filtered * data_mean_filtered[1]
    }

    bottom_margin <- 0
    top_margin <- 5
    
    # if (show_legend) {
    #     bottom_margin <- 0
    #     top_margin <- 0
    # } else {
    #     bottom_margin <- 0  # -40
    #     top_margin <- 5
    # }
    
    max_value <- max(max(abs(markers2[[logFC_group1]])), max(abs(markers2[[logFC_group2]])))
    
    p <- ggplot(markers2, aes(!!sym(logFC_group2), !!sym(logFC_group1))) +
        labs(x = bquote(log[2](.(group2_name) ~ "Fold Change")), y = bquote(log[2](.(group1_name) ~ "Fold Change"))) +
        coord_cartesian(xlim = c(-max_value, max_value), ylim = c(-max_value, max_value)) +
        theme(
            text = element_text(family = "Arial"), 
            plot.title = element_blank(),          # Increase plot title size
            plot.margin = margin(l = 5, r = 11.5, b = bottom_margin, t = top_margin),
            axis.text = element_text(size = rel(axis_numbering_size)),
            axis.title = element_text(size = rel(axis_text_size))
        )
    
    if (max_value > 20) {
        rounded_max_value <- round(max_value / 5) * 5
        p <- p + 
            scale_x_continuous(
                limits = c(-rounded_max_value, rounded_max_value), # Set limits based on max_value
                breaks = seq(-rounded_max_value, rounded_max_value, by = 5), # Major breaks every 5 units
                minor_breaks = setdiff(seq(-rounded_max_value, rounded_max_value, by = 2.5), seq(-rounded_max_value, rounded_max_value, by = 5)) # Minor breaks every 2.5, excluding major breaks
            ) +
            scale_y_continuous(
                limits = c(-rounded_max_value, rounded_max_value), # Set limits based on max_value
                breaks = seq(-rounded_max_value, rounded_max_value, by = 5), # Major breaks every 5 units
                minor_breaks = setdiff(seq(-rounded_max_value, rounded_max_value, by = 2.5), seq(-rounded_max_value, rounded_max_value, by = 5)) # Minor breaks every 2.5, excluding major breaks
            )
    }

    if (!is.null(ccc)) {
        p <- p +
            annotate(
                "text",
                x = Inf,  # Places the annotation on the far right
                y = -Inf, # Places the annotation on the bottom
                label = glue::glue("CCC = {sprintf('%.2f', ccc)}"),
                hjust = 1.04, # Align text to the right
                vjust = -0.28,  # Align text to the bottom
                size = 5.2
            )
    }
    
    if (!outliers_excluded && ((m / m_filtered < 0.95 || m / m_filtered > 1.05) || (b / b_filtered < 0.95 || b / b_filtered > 1.05))) {
        baseline_pca_model_filtered <- glue("PCA fit (excluding outliers): y={sprintf('%.2f', m_filtered)}x+{sprintf('%.2f', b_filtered)}")
    }
    
    baseline_pca_model <- glue("PCA fit: y={sprintf('%.2f', m)}x+{sprintf('%.2f', b)}")
    
    if (show_legend) {
        p <- p +
            geom_abline(aes(slope = 1, intercept = 0, linetype = "y=x", color = "y=x"), show.legend = FALSE, linewidth = 0.5) +
            theme(text = element_text(family = "Arial"), legend.text = element_text(size = 14))
        if (!outliers_excluded && ((m / m_filtered < 0.95 || m / m_filtered > 1.05) || (b / b_filtered < 0.95 || b / b_filtered > 1.05))) {
            p <- p +
                geom_abline(linewidth = 0.5, show.legend = TRUE, aes(slope = m_filtered, intercept = b_filtered, linetype = baseline_pca_model_filtered, color = baseline_pca_model_filtered))
        } else {
            p <- p +
                geom_abline(linewidth = 0.5, show.legend = TRUE, aes(slope = m, intercept = b, linetype = baseline_pca_model, color = baseline_pca_model))
        }
        p <- p +
            scale_color_manual(name = "", values = c("black", "black")) +
            scale_linetype_manual(name = "", values = c(1, 2))
    } else {
        p <- p +
            geom_abline(slope = 1, intercept = 0, linetype = 2, color = "black", show.legend = FALSE, linewidth = 0.5)
        if (!outliers_excluded && ((m / m_filtered < 0.95 || m / m_filtered > 1.05) || (b / b_filtered < 0.95 || b / b_filtered > 1.05))) {
            p <- p +
                geom_abline(slope = m_filtered, intercept = b_filtered, linewidth = 0.5, color = "black", linetype = 1)
        } else {
            p <- p +
                geom_abline(slope = m, intercept = b, linewidth = 0.5, linetype = 1, show.legend = FALSE, color = "black")
        }
    }
    
    if (show_legend) {
        p <- p +
            geom_point(alpha = 0.2, size = 0.3, color = "gray60") +
            theme(text = element_text(family = "Arial"), legend.position = "bottom", legend.direction = "vertical", legend.background = element_blank(), legend.box.spacing = grid::unit(-0.8, "lines"), ) +
            guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))
    } else {
        p <- p +
            ggpointdensity::geom_pointdensity(size = 0.3, alpha = 1, adjust = 0.2, show.legend = FALSE) +
            scico::scale_color_scico(palette = "grayC", direction = -1, end = 0.8)
    }
    
    # p <- p +
        # coord_equal()
        # theme(plot.title = element_text(hjust = 0.48)) # Center the title

    if (!outliers_excluded) {
        # p <- p + labs(title = glue("{group1_name} vs. {group2_name} logFC"))
        default_plotpath <- file_paths_default$logFC_scatterplot_file_path
    } else {
        # p <- p + labs(title = glue("{group1_name} vs. {group2_name} logFC (excluding outliers)"))
        default_plotpath <- file_paths_default$logFC_scatterplot_outliers_removed_file_path
    }


    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = default_plotpath)
        ggsave(filepath, plot = p, dpi = dpi, width = 2100, height = 2100, units = "px")
    }

    return(p)
}


make_bar_plot <- function(df, metric, save = FALSE, filename = NULL) {
    p <- ggplot(df, aes(x = !!sym("Categories"), y = !!sym(metric))) +
        geom_bar(stat = "identity", aes(fill = !!sym("Categories"))) +
        scale_fill_manual(values = ditto_colors) +
        labs(y = metric, x = "Category", title = metric) +
        theme_minimal(base_family = "Arial") +
        theme(text = element_text(family = "Arial"), legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

    print(p)

    if (save) {
        ggsave(filename, plot = p, dpi = dpi)
    }

    return(p)
}


plot_var_explained <- function(eigs_df, npcs = 20, group_names = waiver(), save = FALSE) {
    # eigs_df: each column is eigenvalues from one package, and
    # there must be PC column indicating which PC
    pcs_ve <- eigs_df |>
        pivot_longer(cols = -PC, names_to = "package", values_to = "value")
    required_rows <- pcs_ve %>%
        group_by(package) %>%
        slice_head(n = npcs) %>%
        ungroup()

    p <- ggplot(required_rows, aes(PC, value, color = package, shape = package)) +
        geom_point(size = 3) +
        scale_color_manual(values = pal, labels = group_names) +
        scale_shape_manual(values = c(16, 17), labels = group_names) +
        scale_x_continuous(
            minor_breaks = setdiff(seq(0, 50, by = 1), seq(0, 50, by = 5)),
            breaks = seq(0, 50, by = 5)
        ) +
        # breaks = scales::breaks_pretty(Q = c(1,5,2,4,3))) +
        labs(
            y = "Variance Explained",
            # title = "PC Similarity",
            color = "Package", # Change legend title for color
            shape = "Package"
        ) + # Change legend title for shape
        theme(
            text = element_text(family = "Arial"), 
            legend.text = element_text(size = rel(1.5)),
            legend.title = element_blank(),
            legend.position = c(0.97, 0.93), # Adjust coordinates for top right position
            legend.justification = c("right", "top"), # Anchor point of the legend
            legend.box.just = "right",
            legend.margin = margin(-10, -10, -10, -10),
            plot.title = element_text(size = rel(1), hjust = 0.5),
            axis.title.x = element_text(size = rel(1)),
            axis.text.x = element_text(size = rel(1)),
            axis.title.y = element_text(size = rel(1.4)), # X axis label size
            axis.text.y = element_text(size = rel(1.15)) # Y axis tick number size
        ) +
        guides(
            color = guide_legend(override.aes = list(size = 3)), # Increase legend symbols size for color
            shape = guide_legend(override.aes = list(size = 3)) # Increase legend symbols size for shape
        )

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$pca_elbow_filepath_combined)
        ggsave(filepath, plot = p, dpi = dpi)
    }

    return(p)
}


plot_pca_compare <- function(embeddings1, embeddings2,
                             pcs = 1:2, group1_name = "Seurat", group2_name = "Scanpy", group_labels = waiver(), legend_position = "outside", save = FALSE) {
    # See if needs to be flipped
    embeddings2 <- flip_pcs(embeddings1, embeddings2, pcs)
    df1 <- make_pca_emb_df(embeddings1, group1_name, pcs)
    df2 <- make_pca_emb_df(embeddings2, group2_name, pcs)

    df <- bind_rows(df1, df2)
    dfp <- df |>
        pivot_wider(names_from = package, values_from = c(x, y))
    namesx <- paste("x", c(group1_name, group2_name), sep = "_")
    namesy <- paste("y", c(group1_name, group2_name), sep = "_")
    p <- ggplot(df) +
        geom_point(aes(x, y, color = package, shape = package),
            size = 1.8, alpha = 0.7
        ) +
        geom_segment(
            data = dfp, aes(
                x = .data[[namesx[1]]], y = .data[[namesy[1]]],
                xend = .data[[namesx[2]]], yend = .data[[namesy[2]]]
            ),
            alpha = 0.25, linewidth = 0.1
        ) +
        scale_color_manual(values = pal, labels = group_labels) +
        scale_shape_manual(values = c(16, 17), labels = group_labels) +
        labs( # title = glue("PCA embeddings, {group1_name} vs. {group2_name}"),
            x = paste0("PC", pcs[1]),
            y = paste0("PC", pcs[2]), color = "Package", shape = "Package"
        ) +
        guides(shape = guide_legend(title = ""), color = guide_legend(title = "", override.aes = list(alpha = 1))) +
        theme(
            text = element_text(family = "Arial"), 
            axis.text = element_text(size = rel(axis_numbering_size)), # Increase axis tick labels size
            axis.title = element_text(size = rel(axis_text_size)), # Increase axis titles size
            legend.text = element_text(size = rel(1.5)),
            legend.title = element_text(size = rel(1.5))
        )
    
    if (is.character(legend_position) && legend_position != "outside") {
        p <- p +
            theme(
                text = element_text(family = "Arial"), 
                legend.margin = margin(-10, -10, -10, -10),  # Optional: Adjust margin to move closer to or further from the edges
                legend.box.margin = margin(0, 0, 0, 0),  # Optional: Adjust box margin
                legend.background = element_blank(),  # Optional: Remove background
                legend.key = element_blank()
            )
        
        if (legend_position == "TL") {
            p <- p +
                theme(
                    text = element_text(family = "Arial"), 
                    legend.position = c(0.025, 1.02),
                    legend.justification = c(0, 1),
                    legend.box.just = "left"
                )
        } else if (legend_position == "BR") {
            p <- p +
                theme(
                    text = element_text(family = "Arial"), 
                    legend.position = c(0.96, 0.14),
                    legend.justification = c(1, 1),
                    legend.box.just = "right"
                )
        } else if (legend_position == "BL") {
            p <- p +
                theme(
                    text = element_text(family = "Arial"), 
                    legend.position = c(0.025, 0.14),
                    legend.justification = c(0, 1),
                    legend.box.just = "left"
                )
        } else {  # put it in the top-right
            p <- p +
                theme(
                    text = element_text(family = "Arial"), 
                    legend.position = c(0.96, 1.02),  # Top-right corner
                    legend.justification = c(1, 1),  # Anchor the legend at the top-right
                    legend.box.just = "right"  # Justify the legend box at the right
                )
        }
    }

    if (identical(pcs, 1:2)) {
        default_filepath <- file_paths_default$pca_12_filepath
    } else if (identical(pcs, 3:4)) {
        default_filepath <- file_paths_default$pca_34_filepath
    } else {
        default_filepath <- NULL
    }

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$knee_plot)
        ggsave(filepath, plot = p, width = 2100, height = 2100, dpi = dpi, units = "px")
    }

    return(p)
}


plot_loading_diffs <- function(df, mean_loadings_diff = NULL, save = FALSE) {
    df$differences[is.na(df$differences)] <- .Machine$double.eps
    df$differences[(df$differences == 0)] <- .Machine$double.eps

    min_y <- min(sqrt(.Machine$double.eps), min(df$differences))

    p <- ggplot(df, aes(PC, differences, color = type)) +
        geom_point(color = "black") +
        # geom_col(fill = seurat_v_scanpy_baseline_color, color = "black") +
        scale_x_continuous(
            minor_breaks = setdiff(seq(0, 50, by = 1), seq(0, 50, by = 5)),
            breaks = scales::breaks_width(5)
        ) +
        labs(
            y = "Sine of PCA Eigenvectors θ"
        ) +
        theme(
            text = element_text(family = "Arial"), 
            plot.title = element_text(size = rel(1.5), hjust = 0.5),
            axis.title.x = element_text(size = rel(axis_text_size)), # X axis label size
            axis.title.y = element_text(size = rel(1.4)), # X axis label size
            axis.text.x = element_text(size = rel(axis_numbering_size)), # Y axis tick number size
            axis.text.y = element_text(size = rel(1.15)) # Y axis tick number size
        ) +
        scale_color_manual(values = ditto_colors) +
        guides(color = FALSE) +
        scale_y_log10(limits = c(min_y, 1)) +
        annotation_logticks(sides = "l") +
        geom_hline(yintercept = sqrt(.Machine$double.eps), linetype = 2, color = "gray30") +
        annotate("text", x = 37, y = sqrt(.Machine$double.eps) * 1.4, label = "Double Precision Limit", hjust = 0, vjust = 0, color = "gray30", size = 4)

    if (!is.null(mean_loadings_diff)) {
        mean_loadings_diff <- format(mean_loadings_diff, scientific = TRUE, digits = 2)
        p <- p + annotate("text", x = 20, y = 3.5e-09, label = glue("PC1-3 mean value: {mean_loadings_diff}"), hjust = 0, vjust = 0, color = "black", size = 4)
    }

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$pca_loading_diffs)
        ggsave(filepath, plot = p, dpi = dpi)
    }

    return(p)
}


plot_eigs_diffs <- function(df, save = FALSE) {
    p <- ggplot(df, aes(PC, value, color = type)) +
        geom_point(color = "black") +
        # geom_col(fill = seurat_v_scanpy_baseline_color, color = "black") +
        annotation_logticks(sides = "l") +
        scale_color_manual(values = ditto_colors) +
        geom_hline(
            yintercept = sqrt(.Machine$double.eps), linetype = 2,
            color = "gray30"
        ) +
        labs(
            # title = "Proportion of variance explained",
            y = "Absolute differences",
        ) +
        theme(text = element_text(family = "Arial"), plot.title = element_text(size = rel(1.5), hjust = 0.5)) +
        guides(color = FALSE)

    if (!all(df$value == 0)) {
        p <- p + scale_y_log10()
    }

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$pca_eigs_diff)
        ggsave(filepath, plot = p, dpi = dpi)
    }

    return(p)
}


make_jaccard_plot <- function(jaccards, median_jaccard = NULL, save = FALSE) {
    if (all(jaccards$Jaccard == 1)) {
        # All values are 1, proceed with plotting

        # Create a simple data frame for plotting
        plot_data <- data.frame(Jaccard = 1, Fraction = 1)

        # Plot
        jaccard_plot <- ggplot(plot_data, aes(x = Jaccard, y = Fraction)) +
            geom_point(size = 4, color = "black") + # Adjust point size as needed
            xlim(0.9, 1.1) + # Adjust limits to focus on the point
            ylim(0.9, 1.1) +
            xlab("Jaccard") +
            ylab("Fraction of Data Points") +
            theme_minimal(base_family = "Arial")
    } else {
        jaccard_plot <- ggplot(jaccards, aes(Jaccard)) +
            geom_density(fill = NA, color = seurat_v_scanpy_baseline_color) +
            labs(y = "Density") +
            scale_color_manual(values = ditto_colors)
    }

    jaccard_plot <- jaccard_plot +
        theme(
            text = element_text(family = "Arial"), 
            axis.text = element_text(size = rel(axis_numbering_size)), # Increase axis tick labels size
            axis.title = element_text(size = rel(axis_text_size)), # Increase axis titles size
            plot.margin = margin(5.5, 17, 5.5, 5.5)
        )

    if (!is.null(median_jaccard)) {
        max_density <- max(density(jaccards$Jaccard)$y)
        median_jaccard <- sprintf("%.2f", median_jaccard)
        jaccard_plot <- jaccard_plot + annotate("text", x = 0.57, y = max_density, label = glue("Median jaccard: {median_jaccard}"), hjust = 0, vjust = 0.23, color = "black", size = 6.8)
    }

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$pca_eigs_diff)
        ggsave(filepath, plot = jaccard_plot, dpi = dpi)
    }

    return(jaccard_plot)
}


make_knn_scatterplot <- function(nei_pairs, save = FALSE) {
    knn_scatterplot <- ggplot(nei_pairs, aes(value1, value2)) +
        geom_point(size = 0.5, alpha = 0.3, color = seurat_v_scanpy_baseline_color) +
        geom_abline(slope = 1, intercept = 0, color = "gray30") +
        coord_equal() +
        labs(x = "Degree (Group 1)", y = "Degree (Group 2)") +
        theme(
            text = element_text(family = "Arial"), 
            axis.text = element_text(size = rel(axis_numbering_size)), # Increase axis tick labels size
            axis.title = element_text(size = rel(axis_text_size))
        )

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$knn_scatterplot)
        ggsave(filepath, plot = knn_scatterplot, dpi = dpi)
    }

    return(knn_scatterplot)
}


make_combined_pc_variance_loadings_plot <- function(combined_pc_variance, loading_diffs, save = FALSE) {
    combined_pc_variance_mod <- combined_pc_variance +
        theme(
            axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(), # Remove x-axis labels
            axis.ticks.x = element_blank()
        )

    combined_plot <- combined_pc_variance_mod / loading_diffs

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$combined_pc_variance_loadings_plot)
        ggsave(filepath, plot = combined_plot, dpi = dpi)
    }

    return(combined_plot)
}


make_snn_jaccard_degree_scatterplot <- function(jaccards, neighbor_space = "knn", save = FALSE) {
    xmax <- ceiling(max(abs(jaccards$logged_degree_ratio)))
    xmin <- -xmax
    
    if (!(all(jaccards$Jaccard == 1) && all(jaccards$degree_ratio == 1))) {
        point_size = 0.3
    } else {
        point_size = 1.5
    }

    p <- ggplot(jaccards, aes(x = logged_degree_ratio, y = Jaccard)) +
        ggpointdensity::geom_pointdensity(size = point_size, alpha = 1, adjust = 0.2) +
        scico::scale_color_scico(palette = "grayC", direction = -1, end = 0.8) +
        geom_vline(xintercept = 0, color = "black", linetype = "solid") +
        geom_hline(yintercept = 0, color = "black", linetype = "solid") +
        scale_y_continuous(
            minor_breaks = seq(0, 1, by = 0.1),
            breaks = seq(0, 1, by = 0.2)
        ) +
        scale_x_continuous(breaks = seq(xmin, xmax, by = 1)) +
        coord_cartesian(
            xlim = c(xmin, xmax),
            ylim = c(0, 1)
        ) +
        theme_minimal(base_family = "Arial") +
        stat_function(fun = function(x) 2^x, color = "grey30", linetype = 2, xlim = c(xmin, 0)) +
        stat_function(fun = function(x) 2^(-x), color = "grey30", linetype = 2, xlim = c(0, xmax)) +
        theme(
            text = element_text(family = "Arial"), 
            legend.position = "none",
            panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
            panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
            axis.text = element_text(size = rel(axis_numbering_size)), # Increase axis tick labels size
            axis.title = element_text(size = rel(axis_text_size)),
        )
    # geom_density_2d()
    
    if (!(all(jaccards$Jaccard == 1) && all(jaccards$degree_ratio == 1))) {
        p <- p +
            annotate("text", x = xmin, y = 0.13, label = "y == 2^{x}", parse = TRUE, color = "grey30", size = 4.5) +
            annotate("text", x = xmax, y = 0.13, label = "y == 2^-{x}", parse = TRUE, color = "grey30", size = 4.5)
    }

    if (neighbor_space == "umap") {
        p <- p +
            labs(x = bquote(log[2]("Ratio of UMAP SNN Graph Degrees")), y = "UMAP SNN Graph Neighborhood Jaccard Index")
    } else {
        p <- p +
            labs(x = bquote(log[2]("Ratio of SNN Graph Degrees")), y = "SNN Graph Neighborhood Jaccard Index")
    }

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$jaccard_degree_scatterplot)
        ggsave(filepath, plot = p, dpi = dpi, bg = "white")
    }

    return(p)
}


make_umap_jaccard_plot <- function(jaccards_df, facet = NULL, save = FALSE) {
    if (all(jaccards_df$JaccardIndex == 1)) {
        # Create a simple data frame for plotting
        plot_data <- data.frame(JaccardIndex = 1, Density = 1)

        # Plot
        umap_jaccard_plot <- ggplot(plot_data, aes(x = JaccardIndex, y = Density)) +
            geom_point(size = 4, color = "black") + # Adjust point size as needed
            xlim(0.9, 1.1) + # Adjust limits to focus on the point
            ylim(0.9, 1.1) +
            xlab("UMAP KNN Jaccard") +
            ylab("Density") +
            theme_minimal(base_family = "Arial")
    } else {
        umap_jaccard_plot <- ggplot(jaccards_df, aes(x = JaccardIndex)) +
            geom_density(fill = NA, color = "black") +
            labs(y = "Density", x = "UMAP KNN Jaccard") +
            scale_x_continuous(
                breaks = seq(0, 1, by = 0.1), # Major breakpoints
                minor_breaks = seq(0, 1, by = 0.05) # Minor breakpoints
            ) +
            scale_color_manual(values = ditto_colors) +
            coord_cartesian(xlim = c(0, 1)) +
            theme(
                text = element_text(family = "Arial"), 
                axis.text = element_text(size = rel(axis_numbering_size)), # Increase axis tick labels size
                axis.title = element_text(size = rel(axis_text_size)), # Increase axis titles size
                plot.margin = margin(5.5, 17, 5.5, 5.5)
            )

        if (!is.null(facet)) {
            umap_jaccard_plot <- umap_jaccard_plot +
                facet_wrap(~ .data[[facet]]) +
                theme(text = element_text(family = "Arial"), axis.text = element_text(size = 5.5))
        }
    }

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$umap_jaccard_knn_density)
        ggsave(filepath, plot = umap_jaccard_plot, dpi = dpi, bg = "white", width = 2100, height = 2100, units = "px")
    }

    return(umap_jaccard_plot)
}


plot_heatmap <- function(jacc_seu_sc3, ari_value = NULL, show_axis_titles = FALSE, show_trees = TRUE, save = FALSE) {
    if (show_trees) {
        treeheight_row <- 30
        treeheight_col <- 30
    } else {
        treeheight_row <- 0
        treeheight_col <- 0
    }

    cluster_pheatmap <- as.ggplot(pheatmap::pheatmap(jacc_seu_sc3, color = scales::viridis_pal(end = max(jacc_seu_sc3))(255), treeheight_row = treeheight_row, treeheight_col = treeheight_col))
    if (show_axis_titles || !is.null(ari_value)) {
        cluster_pheatmap <- cluster_pheatmap + coord_cartesian(clip = "off") +
            theme(plot.margin = margin(15, 10, 23, 25, "pt"))
        if (show_axis_titles) {
            cluster_pheatmap <- cluster_pheatmap +
                annotate("text", x = 0.35, y = 0, label = "Seurat", hjust = 0, vjust = 0.8, color = "black", size = 10) +
                annotate("text", x = 0, y = 0.5, label = "Scanpy", hjust = 0.6, vjust = -0.1, color = "black", size = 10, angle = 90)
        }
        if (!is.null(ari_value)) {
            ari_value <- sprintf("%.2f", ari_value)
            cluster_pheatmap <- cluster_pheatmap + annotate("text", x = 0.81, y = 1, label = glue("ARI: {ari_value}"), hjust = 0.6, vjust = 0.2, color = "black", size = 8)
        }
    }

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$pheatmap)
        ggsave(filepath, plot = cluster_pheatmap, dpi = dpi, bg = "white")
    }
    return(cluster_pheatmap)
}


find_group2_colors <- function(clus_df_gather, group1_name, group2_name) {
    clus_df_filtered <- clus_df_gather[, c(group1_name, group2_name, "value")]
    
    clus_df_filtered[[group1_name]] <- paste0("G1_", clus_df_filtered[[group1_name]])
    clus_df_filtered[[group2_name]] <- paste0("G2_", clus_df_filtered[[group2_name]])
    
    g <- igraph::graph_from_data_frame(d = clus_df_filtered, directed = FALSE)
    igraph::V(g)$type <- ifelse(igraph::V(g)$name %in% clus_df_filtered[[group1_name]], TRUE, FALSE)
    
    matching <- igraph::max_bipartite_match(g, weights = clus_df_filtered$value)
    
    keys <- names(matching$matching)
    number_group1_clusters <- length(sub("^G1_", "", keys[grep("^G1_", keys)]))
    number_group2_clusters <- length(sub("^G2_", "", keys[grep("^G2_", keys)]))
    
    # Extract and filter the matching pairs
    non_na_pairs <- matching$matching[!is.na(matching$matching)]
    
    # Filter out pairs where S is matched to C (excluding C matched to S or NA)
    g1_to_g2_pairs <- non_na_pairs[grep("^G1_", names(non_na_pairs))]
    
    
    # Extract numeric indices from the filtered pairs
    g1_indices <- as.numeric(sub("G1_", "", names(g1_to_g2_pairs)))
    g2_indices <- as.numeric(sub("G2_", "", g1_to_g2_pairs))
    
    # Initialize the new colors vector
    group2_colors <- vector("character", number_group2_clusters)
    
    # Assign colors based on the matching
    for (i in seq_along(g1_indices)) {
        # Check if the C index is within the bounds of ditto_colors
        if (g2_indices[i] <= length(ditto_colors)) {
            group2_colors[g2_indices[i]] <- ditto_colors[g1_indices[i]]
        }
    }
    
    remaining_colors <- ditto_colors[number_group1_clusters+1:length(ditto_colors)]
    group2_colors[(group2_colors == "")] <- remaining_colors[1:sum(group2_colors == "")]
    
    return (group2_colors)
}

plot_alluvial <- function(clus_df_gather, group1_name = "Seurat", group2_name = "Scanpy", group1_name_mapping = "Seurat", group2_name_mapping = "Scanpy", color_boxes = TRUE, color_bands = FALSE, alluvial_alpha = 0.5, match_colors = TRUE, save = FALSE, include_labels_in_boxes = FALSE, include_axis_titles = FALSE, show_group_2_box_labels_in_ascending = FALSE) {
    num_levels_group1 <- length(levels(clus_df_gather[[group1_name]]))
    num_levels_group2 <- length(levels(clus_df_gather[[group2_name]]))
    
    if (show_group_2_box_labels_in_ascending) {
        group2_name_mapping <- group2_name
    }

    # Extract colors for each factor, assuming ditto_colors is long enough
    colors_group1 <- ditto_colors[1:num_levels_group1]
    
    if (match_colors) {
        colors_group2 <- find_group2_colors(clus_df_gather, group1_name, group2_name)
    } else {
        colors_group2 <- ditto_colors[1:num_levels_group2]
    }

    colors_group1_reverse <- rev(colors_group1)
    colors_group2_reverse <- rev(colors_group2)

    # Combine the colors
    combined_colors <- c(colors_group1, colors_group2)
    combined_colors_reverse <- c(colors_group1_reverse, colors_group2_reverse)

    # uncomment to attempt mapping
    # p <- ggplot(data = clus_df_gather, aes(axis1 = !!sym(group1_name_mapping), axis2 = !!sym(group2_name_mapping), y = value))  # commented out as of Jan 2025
    p <- ggplot(data = clus_df_gather, aes(axis1 = !!sym(group1_name), axis2 = !!sym(group2_name), y = value))  # uncommented as of Jan 2025

    if (color_bands) {
        if (num_levels_group2 > num_levels_group1) {
            p <- p +
                geom_alluvium(aes(fill = !!sym(group2_name)), alpha = alluvial_alpha) +
                scale_fill_manual(values = colors_group2) +
                labs(fill = NULL)
        } else {
            p <- p +
                geom_alluvium(aes(fill = !!sym(group1_name)), alpha = alluvial_alpha) +
                scale_fill_manual(values = colors_group1) +
                labs(fill = NULL)
        }
    } else {
        p <- p + geom_alluvium()
    }

    if (color_boxes) {
        p <- p + geom_stratum(fill = combined_colors_reverse)
    } else {
        p <- p + geom_stratum()
    }
    
    if (include_labels_in_boxes) {
        p <- p + 
            geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, color = "black")
    } 
    
    if (include_axis_titles) {
        p <- p +
            annotate("text", x = 1, y = max(clus_df_gather$value) + 510, label = group1_name, size = 5, hjust = 0.5) +
            annotate("text", x = 2, y = max(clus_df_gather$value) + 510, label = group2_name, size = 5, hjust = 0.5)
    } 

    p <- p +
        # geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
        theme_void() +
        annotate("text", x = 1.023, y = 1, label = num_levels_group1, hjust = 1, vjust = 1.35, size = 5) + # Adjust x, y for Seurat
        annotate("text", x = 1.978, y = 1, label = num_levels_group2, hjust = 0, vjust = 1.35, size = 5) + # Adjust x, y for Scanpy
        theme(
            text = element_text(family = "Arial"), 
            legend.text = element_text(size = rel(axis_text_size))
        )

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$alluvial)
        ggsave(filepath, plot = p, dpi = dpi, bg = "white")
    }

    return(p)
}


make_umap_plot <- function(dataframe, package, colors = ditto_colors, title = "UMAP", overall_min_dim1 = 0, overall_max_dim1 = 20, overall_min_dim2 = 0, overall_max_dim2 = 20, show_legend = FALSE) {
    centroids <- dataframe %>%
        group_by(cluster) %>%
        summarise(
            UMAP_1 = mean(UMAP_1),
            UMAP_2 = mean(UMAP_2)
        )
    
    p <- ggplot(dataframe, aes(x = UMAP_1, y = UMAP_2, color = as.factor(cluster))) +
        geom_point(size = 0.1) +
        xlim(overall_min_dim1, overall_max_dim1) +
        ylim(overall_min_dim2, overall_max_dim2) +
        ggrepel::geom_text_repel(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = cluster), alpha = 0.6, 
                         vjust = "center", hjust = "center", color = "black", size = 4.5, force = 2.5e-4) +
        # ggrepel::geom_label_repel(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = cluster, fill = cluster), alpha = 0.6, 
                   # vjust = "center", hjust = "center", color = "black", size = 4) +
        # scale_fill_manual(values = ditto_colors, name = "Cluster") +
        labs(
            x = "", # Remove x-axis title
            y = "", # Remove y-axis title
            title = title
        ) +
        theme(
            # text = element_text(family = "Arial"), 
            axis.text.x = element_blank(), # Turn off x-axis numbers
            axis.text.y = element_blank(), # Turn off y-axis numbers
            axis.ticks = element_blank(), # Optionally, turn off axis ticks as well
            axis.text = element_text(size = rel(axis_numbering_size)), # Increase axis tick labels size
            plot.title = element_text(size = rel(2), hjust = 0.5) # Center the title
        ) +
        scale_color_manual(values = colors, name = "Cluster") +
        guides(color = guide_legend(override.aes = list(size = 3)))

    if (!show_legend) {
        p <- p +
            theme(text = element_text(family = "Arial"), 
                  legend.position = "none")
    }

    return(p)
}


plot_umap <- function(group1_umap_info, group1_clusters, group2_umap_info, group2_clusters, colors_group1 = ditto_colors, colors_group2 = ditto_colors, group1 = "Seurat", group2 = "Scanpy", group1_title = "", group2_title = "", show_legend = FALSE, save = FALSE) {
    if ((nrow(group1_umap_info) == 0) || (length(group2_umap_info) == 0)) {
        return (list(NA, NA))
    }
    
    umap_df_group1 <- as.data.frame(group1_umap_info)
    names(umap_df_group1) <- c("UMAP_1", "UMAP_2")
    umap_df_group1$cluster <- group1_clusters
    
    if ("0" %in% levels(umap_df_group1$cluster)) {
        new_dataframe <- as.numeric(as.character(umap_df_group1$cluster)) + 1
        new_levels <- as.numeric(levels(umap_df_group1$cluster)) + 1
        umap_df_group1$cluster <- factor(new_dataframe, levels = new_levels, labels = new_levels)
    }
    
    group1_max_cluster_id <- max(as.numeric(as.character(group1_clusters[!is.na(group1_clusters)])))
    group1_min_cluster_id <- min(as.numeric(as.character(group1_clusters[!is.na(group1_clusters)])))
    umap_df_group1$cluster <- factor(umap_df_group1$cluster, levels = as.character(group1_min_cluster_id:group1_max_cluster_id))

    
    
    
    umap_df_group2 <- as.data.frame(group2_umap_info)
    names(umap_df_group2) <- c("UMAP_1", "UMAP_2")
    umap_df_group2$cluster <- group2_clusters
    
    if ("0" %in% levels(umap_df_group2$cluster)) {
        new_dataframe <- as.numeric(as.character(umap_df_group2$cluster)) + 1
        new_levels <- as.numeric(levels(umap_df_group2$cluster)) + 1
        umap_df_group2$cluster <- factor(new_dataframe, levels = new_levels, labels = new_levels)
    }
    
    group2_max_cluster_id <- max(as.numeric(as.character(group2_clusters[!is.na(group2_clusters)])))
    group2_min_cluster_id <- min(as.numeric(as.character(group2_clusters[!is.na(group2_clusters)])))
    umap_df_group2$cluster <- factor(umap_df_group2$cluster, levels = as.character(group2_min_cluster_id:group2_max_cluster_id))


    group1_min_dim1 <- min(group1_umap_info[, 1])
    group1_min_dim2 <- min(group1_umap_info[, 2])
    group1_max_dim1 <- max(group1_umap_info[, 1])
    group1_max_dim2 <- max(group1_umap_info[, 2])

    group2_min_dim1 <- min(group2_umap_info[, 1])
    group2_min_dim2 <- min(group2_umap_info[, 2])
    group2_max_dim1 <- max(group2_umap_info[, 1])
    group2_max_dim2 <- max(group2_umap_info[, 2])

    overall_min_dim1 <- min(group1_min_dim1, group2_min_dim1)
    overall_min_dim2 <- min(group1_min_dim2, group2_min_dim2)
    overall_max_dim1 <- max(group1_max_dim1, group2_max_dim1)
    overall_max_dim2 <- max(group1_max_dim2, group2_max_dim2)

    if (group1_title == "") {
        group1_title <- tools::toTitleCase(group1)
    }

    if (group2_title == "") {
        group2_title <- tools::toTitleCase(group2)
    }

    p1 <- make_umap_plot(umap_df_group1, package = group1, colors = colors_group1, title = group1_title, overall_min_dim1, overall_max_dim1, overall_min_dim2, overall_max_dim2, show_legend = show_legend)
    p2 <- make_umap_plot(umap_df_group2, package = group2, colors = colors_group2, title = group2_title, overall_min_dim1, overall_max_dim1, overall_min_dim2, overall_max_dim2, show_legend = show_legend)


    if (length(save) > 1) {
        save_1 <- save[1]
        save_2 <- save[2]
    } else {
        save_1 <- save
        save_2 <- save
    }

    if (save_1 == TRUE || (is.character(save_1)) && !is.na(save_1)) {
        filepath_1 <- make_save_path(filepath = save_1, default_filepath = file_paths_default$umap_seu)
        ggsave(filepath_1, plot = p1, dpi = dpi, bg = "white", width = 2100, height = 2100, units = "px")
    }

    if (save_2 == TRUE || (is.character(save_2)) && !is.na(save_2)) {
        filepath_2 <- make_save_path(filepath = save_2, default_filepath = file_paths_default$umap_scan)
        ggsave(filepath_2, plot = p2, dpi = dpi, bg = "white", width = 2100, height = 2100, units = "px")
    }

    return(list(p1, p2))
}


calculate_individual_de_stats <- function(markers2, column_name, column_equation) {
    markers2[[column_name]] <- column_equation

    # Calculating the mean magnitude of the difference
    mean_stat <- mean(markers2[[column_name]], na.rm = TRUE)
    print(glue("Mean magnitude of {column_name}: {mean_stat}"))

    median_stat <- median(markers2[[column_name]], na.rm = TRUE)
    print(glue("Median magnitude of {column_name}: {median_stat}"))

    variance_stat <- var(markers2[[column_name]], na.rm = TRUE)
    print(glue("Variance of magnitude of {column_name}: {variance_stat}"))

    print("------------------------------------------------------------------------------")

    return(markers2)
}


calculate_de_stats <- function(markers2, group1_name = "Seurat", group2_name = "Scanpy", save = FALSE) {
    de_names <- scatterplot_naming(group1_name, group2_name)
    logFC_group1 <- de_names$logFC_group1
    logFC_group2 <- de_names$logFC_group2
    p_val_adj_group1 <- de_names$p_val_adj_group1
    p_val_adj_group2 <- de_names$p_val_adj_group2

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$de_stats_file)
        sink(filepath, split = TRUE, append = TRUE)
    }
    
    filtered_markers2 <- markers2[abs(markers2[[logFC_group1]]) <= 20 & abs(markers2[[logFC_group2]]) <= 20, ]
    
    logFC_ccc_results <- CCC(x = filtered_markers2[[logFC_group1]], y = filtered_markers2[[logFC_group2]])
    logFC_ccc <- logFC_ccc_results$rho.c$est
    print(glue("logFC CCC: {logFC_ccc}"))
    
    markers2 <- markers2 %>% mutate(CCC = logFC_ccc)

    logFC_difference_magnitude_equation <- abs(markers2[[logFC_group1]] - markers2[[logFC_group2]])
    markers2 <- calculate_individual_de_stats(markers2, column_name = "logFC_difference_magnitude", column_equation = logFC_difference_magnitude_equation)

    logFC_difference_signed_equation <- (markers2[[logFC_group1]] - markers2[[logFC_group2]])
    markers2 <- calculate_individual_de_stats(markers2, column_name = "logFC_difference_signed", column_equation = logFC_difference_signed_equation)


    # filtered_markers2 <- subset(markers2, abs(logFC_r) < 20 & abs(logFC_py) < 20)
    # filtered_logFC_difference_magnitude_equation <- abs(filtered_markers2$logFC_r - filtered_markers2$logFC_py)
    # filtered_logFC_difference_signed_equation <- (filtered_markers2$logFC_r - filtered_markers2$logFC_py)
    #
    # filtered_markers2 <- calculate_individual_stat(filtered_markers2, column_name="filtered_logFC_difference_magnitude", column_equation=filtered_logFC_difference_magnitude_equation)
    # filtered_markers2 <- calculate_individual_stat(filtered_markers2, column_name="filtered_logFC_difference_signed", column_equation=filtered_logFC_difference_signed_equation)

    # p_val_adj
    pvaladj_difference_magnitude_equation <- abs(markers2[[p_val_adj_group1]] - markers2[[p_val_adj_group2]])
    markers2 <- calculate_individual_de_stats(markers2, column_name = "pvaladj_difference_magnitude", column_equation = pvaladj_difference_magnitude_equation)

    pvaladj_difference_signed_equation <- (markers2[[p_val_adj_group1]] - markers2[[p_val_adj_group2]])
    markers2 <- calculate_individual_de_stats(markers2, column_name = "pvaladj_difference_signed", column_equation = pvaladj_difference_signed_equation)

    p_both_small <- nrow(markers2 %>% filter(!!sym(p_val_adj_group1) <= 0.05, !!sym(p_val_adj_group2) <= 0.05))

    p1_small <- nrow(markers2 %>% filter(!!sym(p_val_adj_group1) <= 0.05, !!sym(p_val_adj_group2) > 0.05))

    p2_small <- nrow(markers2 %>% filter(!!sym(p_val_adj_group1) > 0.05, !!sym(p_val_adj_group2) <= 0.05))

    p_both_big <- nrow(markers2 %>% filter(!!sym(p_val_adj_group1) > 0.05, !!sym(p_val_adj_group2) > 0.05))


    fraction_flipped <- (p1_small + p2_small) / nrow(markers2)

    print(glue("Adjusted p value, fraction <=0.05 in both groups: {p_both_small/nrow(markers2)}"))
    print(glue("Adjusted p value, fraction <=0.05 in {group1_name} but >0.05 in {group2_name}: {p1_small/nrow(markers2)}"))
    print(glue("Adjusted p value, fraction <=0.05 in {group2_name} but >0.05 in {group1_name}: {p2_small/nrow(markers2)}"))
    print(glue("Adjusted p value, fraction >0.05 in both groups: {p_both_big/nrow(markers2)}"))
    print(glue("Adjusted p value, fraction that flipped across 0.05 threshold: {fraction_flipped}"))
    
    # Add new columns with ranks for p_val_adj_py and p_val_adj_r (for Spearman)
    markers2_ranked <- markers2 %>%
        mutate(
            rank_p_val_adj_group1 = rank(!!sym(p_val_adj_group1), ties.method = "min"),
            rank_p_val_adj_group2 = rank(!!sym(p_val_adj_group2), ties.method = "min"),
        )
    
    # Compute Spearman correlation between the two rank columns
    spearman_correlation <- cor(markers2_ranked$rank_p_val_adj_group1, markers2_ranked$rank_p_val_adj_group2, method = "spearman")
    print(glue("Adjusted p value, Spearman correlation: {spearman_correlation}"))
    
    markers2 <- markers2 %>% mutate(spearman = spearman_correlation)

    # # FC (exponentiated)
    # FC_difference_magnitude_equation <- abs(markers2$FC_r - markers2$FC_py)
    # markers2 <- calculate_individual_de_stats(markers2, column_name="FC_difference_magnitude", column_equation=FC_difference_magnitude_equation)
    #
    # FC_difference_signed_equation <- (markers2$FC_r - markers2$FC_py)
    # markers2 <- calculate_individual_de_stats(markers2, column_name="FC_difference_signed", column_equation=FC_difference_signed_equation)

    if (save == TRUE || is.character(save)) {
        sink()
    }

    return(markers2)
}


make_upset_scanpy <- function(adata1, adata2, comparison = comparison, group_names = seurat_group_names_default, before_filtering = FALSE, as_ggplot = FALSE, save = FALSE) {
    if (comparison == "Gene") {
        data <- list(
            Scanpy1 = unlist(adata1$var_names$tolist()),
            Scanpy2 = unlist(adata2$var_names$tolist())
        )
        default_plotpath <- file_paths_default$upset_genes
    } else if (comparison == "Cell") {
        data <- list(
            Scanpy1 <- unlist(adata1$obs_names$tolist()),
            Scanpy2 <- unlist(adata2$obs_names$tolist())
        )
        default_plotpath <- file_paths_default$upset_cells
    } else if (comparison == "HVG") {
        py_run_string("scanpy_highly_variable_genes1 = adata1.var.index[adata1.var['highly_variable']]")
        py_run_string("scanpy_highly_variable_genes2 = adata2.var.index[adata2.var['highly_variable']]")
        scanpy_highly_variable_genes1 <- unlist(py$scanpy_highly_variable_genes1$tolist())
        scanpy_highly_variable_genes2 <- unlist(py$scanpy_highly_variable_genes2$tolist())
        data <- list(
            Scanpy1 = scanpy_highly_variable_genes1,
            Scanpy2 = scanpy_highly_variable_genes2
        )
        default_plotpath <- file_paths_default$upset_hvg
    } else if (comparison == "Marker" || comparison == "Marker Gene") {
        data <- list(
            Scanpy1 = adata1,
            Scanpy2 = adata2
        )
        default_plotpath <- default_plotpath <- file_paths_default$upset_markers
    }

    data <- setNames(data, unlist(group_names))

    p <- upset_plot_general(data, group1_name = group_names$Scanpy1, group2_name = group_names$Scanpy2, comparison = comparison, before_filtering = before_filtering, as_ggplot = as_ggplot, save = save, default_plotpath = default_plotpath)

    print(p)

    return(p)
}


make_upset_seurat <- function(group1, group2, comparison = comparison, group_names = seurat_group_names_default, before_filtering = FALSE, as_ggplot = FALSE, save = FALSE) {
    if (comparison == "Gene") {
        data <- list(
            Seurat1 <- rownames(group1),
            Seurat2 <- rownames(group2)
        )
        default_plotpath <- file_paths_default$upset_genes
    } else if (comparison == "Cell") {
        data <- list(
            Seurat1 <- colnames(group1),
            Seurat2 <- colnames(group2)
        )
        default_plotpath <- file_paths_default$upset_cells
    } else if (comparison == "HVG") {
        data <- list(
            Seurat1 = VariableFeatures(group1),
            Seurat2 = VariableFeatures(group2)
        )
        default_plotpath <- file_paths_default$upset_hvg
    } else if (comparison == "Marker" || comparison == "Marker Gene") {
        data <- list(
            Seurat1 = group1,
            Seurat2 = group2
        )
        default_plotpath <- default_plotpath <- file_paths_default$upset_markers
    }

    data <- setNames(data, unlist(group_names))

    p <- upset_plot_general(data, group1_name = group_names$Seurat1, group2_name = group_names$Seurat2, comparison = comparison, before_filtering = before_filtering, as_ggplot = as_ggplot, save = save, default_plotpath = default_plotpath)

    print(p)

    return(p)
}


make_euler_scanpy <- function(adata1, adata2, comparison, before_QC = FALSE, group_names = seurat_group_names_default, save_plot = FALSE, save_stats = FALSE) {
    if (comparison == "Gene") {
        data <- list(
            Scanpy1 = unlist(adata1$var_names$tolist()),
            Scanpy2 = unlist(adata2$var_names$tolist())
        )
        default_plotpath <- file_paths_default$upset_genes
    } else if (comparison == "Cell") {
        data <- list(
            Scanpy1 <- unlist(adata1$obs_names$tolist()),
            Scanpy2 <- unlist(adata2$obs_names$tolist())
        )
        default_plotpath <- file_paths_default$upset_cells
    } else if (comparison == "HVG") {
        py_run_string("scanpy_highly_variable_genes1 = adata1.var.index[adata1.var['highly_variable']]")
        py_run_string("scanpy_highly_variable_genes2 = adata2.var.index[adata2.var['highly_variable']]")
        scanpy_highly_variable_genes1 <- unlist(py$scanpy_highly_variable_genes1$tolist())
        scanpy_highly_variable_genes2 <- unlist(py$scanpy_highly_variable_genes2$tolist())
        data <- list(
            Scanpy1 = scanpy_highly_variable_genes1,
            Scanpy2 = scanpy_highly_variable_genes2
        )
        default_plotpath <- file_paths_default$upset_hvg
    } else if (comparison == "Marker" || comparison == "Marker Gene") {
        data <- list(
            Scanpy1 = adata1,
            Scanpy2 = adata2
        )
        default_plotpath <- default_plotpath <- file_paths_default$upset_markers
    }

    data <- setNames(data, unlist(group_names))
    scanpy1_name <- unlist(group_names)["Scanpy1"]
    scanpy2_name <- unlist(group_names)["Scanpy2"]
    overlap_name <- glue("{scanpy1_name}&{scanpy2_name}")

    euler_data <- euler(data)
    scan1_unique <- euler_data$original.values[scanpy1_name]
    scan2_unique <- euler_data$original.values[scanpy2_name]
    total_overlap <- euler_data$original.values[overlap_name]


    if (save_stats == TRUE || is.character(save_stats)) {
        default_statpath <- ifelse(before_QC, file_paths_default$euler_stats_before_QC_file, file_paths_default$euler_stats_after_QC_file)
        filepath <- make_save_path(filepath = save_stats, default_filepath = default_statpath)
        sink(filepath, append = TRUE, split = TRUE)
    }

    print(glue("Total {comparison}s in {scanpy1_name}: {scan1_unique+total_overlap}"))
    print(glue("Total {comparison}s in {scanpy2_name}: {scan2_unique+total_overlap}"))
    print(glue("Total {comparison}s overlapping: {total_overlap}"))
    print(glue("Fraction of {comparison}s in {scanpy1_name} overlapping in {scanpy2_name}: {total_overlap/(total_overlap + scan2_unique)}"))
    print(glue("Fraction of {comparison}s in {scanpy2_name} overlapping in {scanpy1_name}: {total_overlap/(total_overlap + scan1_unique)}"))
    print(glue("{comparison}s Jaccard: {total_overlap/(total_overlap + scan1_unique + scan2_unique)}"))

    if (save_stats == TRUE || is.character(save_stats)) {
        sink()
    }

    euler_plot <- plot(euler_data, fills = pal, quantities = list(cex = 2), labels = list(cex = 2))

    euler_plot_ggplot_compatible <- as.ggplot(euler_plot)

    if (before_QC) {
        title <- glue("{comparison}s (before QC)")
    } else {
        title <- glue("{comparison}s")
    }

    euler_plot_ggplot_compatible <- euler_plot_ggplot_compatible +
        ggtitle(glue("{comparison}s")) +
        theme(plot.title = element_text(size = 40, hjust = 0.5))

    if (save_plot == TRUE || is.character(save_plot)) {
        filepath <- make_save_path(filepath = save_plot, default_filepath = default_plotpath)
        ggsave(filepath, plot = euler_plot_ggplot_compatible, dpi = dpi, bg = "white", width = 2300, height = 2100, units = "px")
    }

    return(euler_plot_ggplot_compatible)
}


make_euler_seurat <- function(seu1, seu2, comparison, before_QC = FALSE, group_names = seurat_group_names_default, save_plot = FALSE, save_stats = FALSE) {
    if (comparison == "Gene") {
        data <- list(
            Seurat1 = rownames(seu1@assays$RNA$counts),
            Seurat2 = rownames(seu2@assays$RNA$counts)
        )
        default_plotpath <- ifelse(before_QC, file_paths_default$euler_before_qc_cell_file_path, file_paths_default$euler_after_qc_cell_file_path)
    } else if (comparison == "Cell") {
        data <- list(
            Seurat1 = colnames(seu1@assays$RNA$counts),
            Seurat2 = colnames(seu2@assays$RNA$counts)
        )
        default_plotpath <- ifelse(before_QC, file_paths_default$euler_before_qc_gene_file_path, file_paths_default$euler_after_qc_gene_file_path)
    } else if (comparison == "HVG") {
        data <- list(
            Seurat1 = VariableFeatures(seu1),
            Seurat2 = VariableFeatures(seu2)
        )
        default_plotpath <- file_paths_default$euler_after_qc_hvg_file_path
    } else if (comparison == "Marker" || comparison == "Marker Gene") {
        data <- list(
            Seurat1 = seu1,
            Seurat2 = seu2
        )
        default_plotpath <- file_paths_default$euler_after_qc_marker_file_path
    }

    data <- setNames(data, unlist(group_names))
    seurat1_name <- unlist(group_names)["Seurat1"]
    seurat2_name <- unlist(group_names)["Seurat2"]
    overlap_name <- glue("{seurat1_name}&{seurat2_name}")

    euler_data <- euler(data)
    seu1_unique <- euler_data$original.values[seurat1_name]
    seu2_unique <- euler_data$original.values[seurat2_name]
    total_overlap <- euler_data$original.values[overlap_name]


    if (save_stats == TRUE || is.character(save_stats)) {
        default_statpath <- ifelse(before_QC, file_paths_default$euler_stats_before_QC_file, file_paths_default$euler_stats_after_QC_file)
        filepath <- make_save_path(filepath = save_stats, default_filepath = default_statpath)
        sink(filepath, append = TRUE, split = TRUE)
    }

    print(glue("Total {comparison}s in {seurat1_name}: {seu1_unique+total_overlap}"))
    print(glue("Total {comparison}s in {seurat2_name}: {seu2_unique+total_overlap}"))
    print(glue("Total {comparison}s overlapping: {total_overlap}"))
    print(glue("Fraction of {comparison}s in {seurat1_name} overlapping in {seurat2_name}: {total_overlap/(total_overlap + seu2_unique)}"))
    print(glue("Fraction of {comparison}s in {seurat2_name} overlapping in {seurat1_name}: {total_overlap/(total_overlap + seu1_unique)}"))
    print(glue("{comparison}s Jaccard: {total_overlap/(total_overlap + seu1_unique + seu2_unique)}"))

    if (save_stats == TRUE || is.character(save_stats)) {
        sink()
    }

    euler_plot <- plot(euler_data, fills = pal, quantities = list(cex = 2), labels = list(cex = 2))

    euler_plot_ggplot_compatible <- as.ggplot(euler_plot)

    if (before_QC) {
        title <- glue("{comparison}s (before QC)")
    } else {
        title <- glue("{comparison}s")
    }

    euler_plot_ggplot_compatible <- euler_plot_ggplot_compatible +
        ggtitle(glue("{comparison}s")) +
        theme(plot.title = element_text(size = 40, hjust = 0.5))

    if (save_plot == TRUE || is.character(save_plot)) {
        filepath <- make_save_path(filepath = save_plot, default_filepath = default_plotpath)
        ggsave(filepath, plot = euler_plot_ggplot_compatible, dpi = dpi, bg = "white", width = 2300, height = 2100, units = "px")
    }

    return(euler_plot_ggplot_compatible)
}


make_violin_nfeatures_seu <- function(seu1, seu2, group1_name = "Group 1", group2_name = "Group 2", save = FALSE) {
    max_features <- max(seu1$nFeature_RNA, seu2$nFeature_RNA)
    
    vln_seu1 <- VlnPlot(seu1, features = "nFeature_RNA", pt.size = 0, cols = group1_color) +
        coord_cartesian(ylim = c(0, max_features), clip = "off") +
        scale_y_continuous(breaks = seq(0, max_features, by = 2000)) +
        annotate("text", x = 1, y = -Inf, label = group1_name, vjust = 1.7, hjust = 0.5, size = 9) +
        theme(
            text = element_text(family = "Arial"), 
            legend.position = "none",
            plot.margin = margin(r = 0, b = 35, t = 10, unit = "pt"),
            plot.title = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = rel(axis_numbering_size)), # Increase axis tick labels size
            axis.title.x = element_blank(),
            axis.title.y = element_blank()
        )

    vln_seu2 <- VlnPlot(seu2, features = "nFeature_RNA", pt.size = 0, cols = group2_color) +
        coord_cartesian(ylim = c(0, max_features), clip = "off") +
        scale_y_continuous(breaks = seq(0, max_features, by = 2000)) +
        annotate("text", x = 1, y = -Inf, label = group2_name, vjust = 1.7, hjust = 0.5, size = 9) +
        theme(
            text = element_text(family = "Arial"), 
            legend.position = "none",
            plot.margin = margin(l = 0, b = 35, r = 30, unit = "pt"),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            plot.title = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank()
        )

    # Combine the plots
    combined_plot <- vln_seu1 + vln_seu2 + plot_layout(ncol = 2, widths = c(0.3, 0.3))
    
    combined_plot <- combined_plot + plot_annotation(title = "Number of genes per cell", 
                                                     theme = theme(
                                                         text = element_text(family = "Arial"), plot.title = element_text(hjust = 0.8, vjust = 0, size = 30)))

    if (save == TRUE || is.character(save)) {
        filepath <- make_save_path(filepath = save, default_filepath = file_paths_default$violin_counts_comparison)
        ggsave(filepath, plot = combined_plot, dpi = 300, width = 2100, height = 2100, units = "px")
    }

    return(combined_plot)
}



# From DescTools v0.99.53
CCC <- function(x, y, ci = "z-transform", conf.level = 0.95, na.rm = FALSE){
    
    dat <- data.frame(x, y)
    
    if(na.rm) dat <- na.omit(dat)
    #   id <- complete.cases(dat)
    #   nmissing <- sum(!complete.cases(dat))
    #   dat <- dat[id,]
    
    
    N. <- 1 - ((1 - conf.level) / 2)
    zv <- qnorm(N., mean = 0, sd = 1)
    lower <- "lwr.ci"
    upper <- "upr.ci"
    
    k <- length(dat$y)
    yb <- mean(dat$y)
    sy2 <- var(dat$y) * (k - 1) / k
    sd1 <- sd(dat$y)
    
    xb <- mean(dat$x)
    sx2 <- var(dat$x) * (k - 1) / k
    sd2 <- sd(dat$x)
    
    r <- cor(dat$x, dat$y)
    sl <- r * sd1 / sd2
    
    sxy <- r * sqrt(sx2 * sy2)
    p <- 2 * sxy / (sx2 + sy2 + (yb - xb)^2)
    
    delta <- (dat$x - dat$y)
    rmean <- apply(dat, MARGIN = 1, FUN = mean)
    blalt <- data.frame(mean = rmean, delta)
    
    # Scale shift:
    v <- sd1 / sd2
    # Location shift relative to the scale:
    u <- (yb - xb) / ((sx2 * sy2)^0.25)
    # Variable C.b is a bias correction factor that measures how far the best-fit line deviates from a line at 45 degrees (a measure of accuracy). No deviation from the 45 degree line occurs when C.b = 1. See Lin (1989 page 258).
    # C.b <- (((v + 1) / (v + u^2)) / 2)^-1
    
    # The following taken from the Stata code for function "concord" (changed 290408):
    C.b <- p / r
    
    # Variance, test, and CI for asymptotic normal approximation (per Lin (March 2000) Biometrics 56:325-5):
    sep = sqrt(((1 - ((r)^2)) * (p)^2 * (1 - ((p)^2)) / (r)^2 + (2 * (p)^3 * (1 - p) * (u)^2 / r) - 0.5 * (p)^4 * (u)^4 / (r)^2 ) / (k - 2))
    ll = p - zv * sep
    ul = p + zv * sep
    
    # Statistic, variance, test, and CI for inverse hyperbolic tangent transform to improve asymptotic normality:
    t <- log((1 + p) / (1 - p)) / 2
    set = sep / (1 - ((p)^2))
    llt = t - zv * set
    ult = t + zv * set
    llt = (exp(2 * llt) - 1) / (exp(2 * llt) + 1)
    ult = (exp(2 * ult) - 1) / (exp(2 * ult) + 1)
    
    if(ci == "asymptotic"){
        rho.c <- as.data.frame(cbind(p, ll, ul))
        names(rho.c) <- c("est", lower, upper)
        rval <- list(rho.c = rho.c, s.shift = v, l.shift = u, C.b = C.b, blalt = blalt ) # , nmissing = nmissing)
    }
    
    else if(ci == "z-transform"){
        rho.c <- as.data.frame(cbind(p, llt, ult))
        names(rho.c) <- c("est", lower, upper)
        rval <- list(rho.c = rho.c, s.shift = v, l.shift = u, C.b = C.b, blalt = blalt) #, nmissing = nmissing)
    }
    return(rval)
}



create_volcano_plots <- function(data, output_base_path, package, dpi = 300) {
    dir.create(output_base_path, recursive = TRUE, showWarnings = FALSE)
    
    if (package == "Seurat") {
        logfc_col <- "logFC_r"
        pval_col <- "p_val_adj_r"
    } else if (package == "Scanpy") {
        logfc_col <- "logFC_py"
        pval_col <- "p_val_adj_py"
    } else{
        stop("Invalid package specified. Please use 'Seurat' or 'Scanpy'.")
    }
    
    # Dynamically access the specified columns
    logFC <- sym(logfc_col)
    pval <- sym(pval_col)
    
    # Get unique clusters
    clusters <- unique(data$cluster)
    
    for (clust in clusters) {
        # Subset data for the current cluster
        subset_data <- subset(data, cluster == clust)
        
        # Cap y-values at 5 for display purposes
        subset_data <- subset_data %>% mutate(neg_log10_pval = pmin(-log10(!!pval), 5))
        
        # Create a new color column
        subset_data$color <- "gray"  # Default color
        subset_data$color[subset_data[[logfc_col]] >= 1 & subset_data[[pval_col]] < 0.05] <- "red"
        subset_data$color[subset_data[[logfc_col]] <= -1 & subset_data[[pval_col]] < 0.05] <- "blue"
        
        # Create volcano plot
        p <- ggplot(subset_data, aes(x = !!logFC, y = neg_log10_pval)) +
            geom_point(aes(color = color), alpha = 0.6) +
            scale_color_identity() +  # Use color directly from the data
            geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Dashed lines for logFC thresholds
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Dashed line for p-value threshold
            scale_x_continuous(
                limits = c(-10, 10),
                breaks = seq(-10, 10, 2),  # Major ticks every 2
                minor_breaks = seq(-10, 10, 1)  # Minor ticks every 1
            ) +
            scale_y_continuous(
                limits = c(0, 5),
                breaks = seq(0, 5, 1),
                labels = c(0, 1, 2, 3, 4, "≤5")  # Custom label for the top value
            ) +
            labs(
                title = glue("{package} Cluster {clust}"),
                x = expression(Log[2] ~ "Fold Change"),
                y = expression(-Log[10] ~ "Adjusted P-Value")
            ) +
            theme_minimal(base_size = 14) +
            theme(
                plot.background = element_rect(fill = "white", color = NA),  # White background
                plot.title = element_text(hjust = 0.5),  # Center title
                panel.grid.minor = element_line(color = "lightgray", linetype = "dotted")  # Minor grid lines
            )
        
        # Save plot as TIFF
        output_file <- file.path(output_base_path, paste0("volcano_", package, "_cluster_", clust, ".tiff"))
        ggsave(filename = output_file, plot = p, dpi = dpi, width = 6, height = 4, units = "in")
    }
    
    message("Volcano plots saved to: ", output_base_path)
    
    # Create a stand-alone legend
    create_legend <- function() {
        dummy_data <- data.frame(
            x = c(1, -1, 0),
            y = c(1, 1, 1),
            significance = factor(
                c("Up", 
                  "Down", 
                  "Normal"),
                levels = c("Up", 
                           "Down", 
                           "Normal"),
            )
        )
        
        legend_plot <- ggplot(dummy_data, aes(x = x, y = y, color = significance)) +
            geom_point(size = 4) +
            scale_color_manual(
                values = c(
                    "Up" = "red",
                    "Down" = "blue",
                    "Normal" = "gray"
                )
            ) +
            labs(color = "Significance") +
            theme_void() +  # Remove all axes and backgrounds
            theme(
                legend.position = "right",
                legend.title = element_text(size = 14, face = "bold"),
                legend.text = element_text(size = 12)
            )
        
        # Extract and save the legend
        legend <- cowplot::get_legend(legend_plot)
        legend_file <- file.path(output_base_path, "volcano_legend.tiff")
        ggsave(legend_file, plot = gridExtra::grid.arrange(legend), dpi = dpi, width = 3, height = 2, units = "in")
        message("Legend saved to: ", legend_file)
    }
    
    # Call the function to create and save the legend
    create_legend()
}


run_gget_enrichr <- function(genes_list, databases, filename_base, output_base_path_biological_output, ensembl=FALSE, custom_background_list=NULL, save_kegg=FALSE) {
    if (class(genes_list) != "list") {
        genes_list <- as.list(genes_list)
    }
    
    if (is.null(custom_background_list)) {
        background <- TRUE
        ensembl_bkg <- FALSE
    } else {
        background <- FALSE
        ensembl_bkg <- ensembl
        if (class(custom_background_list) != "list") {
            custom_background_list <- as.list(unique(custom_background_list))
        }
    }
    
    for (database in databases) {
        if (save_kegg == TRUE) {
            if (database == "pathway" || grepl("^KEGG", database)) {
                kegg_out_path <- file.path(output_base_path_biological_output, paste0("kegg_", database, ".png"))
            } else {
                kegg_out_path <- NULL
            }
        } else {
            kegg_out_path <- NULL
        }
        
        gget$enrichr(genes_list, database=database, species='human', background=background, background_list=custom_background_list, ensembl=ensembl, ensembl_bkg=ensembl_bkg, plot=TRUE, save=TRUE, kegg_out=kegg_out_path)
        # Define paths
        png_destination_path <- file.path(output_base_path_biological_output, paste0(filename_base, "_", database, ".png"))
        csv_destination_path <- file.path(output_base_path_biological_output, paste0(filename_base, "_", database, ".csv"))
        dir.create(output_base_path_biological_output, recursive = TRUE, showWarnings = FALSE)
        
        # Define source files
        png_source_path <- "./gget_enrichr_results.png"
        csv_source_path <- "./gget_enrichr_results.csv"
        
        # Move the .png file if it exists
        if (file.exists(png_source_path)) {
            file.rename(png_source_path, png_destination_path)
            print(glue("{png_source_path} moved to {png_destination_path}"))
        } else {
            print(sprintf("File not found: %s", png_source_path))
        }
        
        # Move the .csv file if it exists
        if (file.exists(csv_source_path)) {
            file.rename(csv_source_path, csv_destination_path)
            print(glue("{csv_source_path} moved to {csv_destination_path}"))
        } else {
            print(sprintf("File not found: %s", csv_source_path))
        }
    }
}


make_volcano_go_plots <- function(markers2, databases = NULL, ensembl=FALSE, custom_background_list=NULL, save_kegg=FALSE) {
    dir.create(output_base_path, recursive = TRUE, showWarnings = FALSE)
    
    gene_name_column <- if (ensembl) markers2$gene else markers2$gene_symbol
    
    # Get unique clusters
    clusters <- unique(markers2$cluster)
    
    # Loop over each cluster
    for (clust in clusters) {
        # Subset data for the current cluster
        subset_data <- markers2[markers2$cluster == clust, ]
        
        # Get gene sets for Seurat and Scanpy thresholds
        genes_seurat <- unlist(list(unique(gene_name_column[abs(subset_data$logFC_r) >= 1 & subset_data$p_val_adj_r < 0.05])))
        genes_seurat <- unique(genes_seurat)
        genes_seurat <- genes_seurat[!is.null(genes_seurat) & !is.na(genes_seurat)]
        
        genes_scanpy <- unlist(list(unique(gene_name_column[abs(subset_data$logFC_py) >= 1 & subset_data$p_val_adj_py < 0.05])))
        genes_scanpy <- unique(genes_scanpy)
        genes_scanpy <- genes_scanpy[!is.null(genes_scanpy) & !is.na(genes_scanpy)]
        
        enrichr_filename_base = glue("enrichr_seurat_cluster{clust}")
        run_gget_enrichr(genes_seurat, databases=databases, filename_base=enrichr_filename_base, output_base_path_biological_output=output_base_path_biological_output, custom_background_list=custom_background_list, save_kegg=save_kegg, ensembl=ensembl)  # custom_background_list=background_list
        
        enrichr_filename_base = glue("enrichr_scanpy_cluster{clust}")
        run_gget_enrichr(genes_scanpy, databases=databases, filename_base=enrichr_filename_base, output_base_path_biological_output=output_base_path_biological_output, custom_background_list=custom_background_list, save_kegg=save_kegg, ensembl=ensembl)  # custom_background_list=background_list
        
        print(glue("Done with cluster {clust}"))
    }
}

# Function to calculate Jaccard index
calculate_jaccard <- function(set1, set2) {
    intersection_size <- length(intersect(set1, set2))
    union_size <- length(union(set1, set2))
    if (union_size == 0) return(0)  # Avoid division by zero
    return(intersection_size / union_size)
}

make_volcano_jaccards_and_upsets <- function(markers2, make_upset_plots = TRUE) {
    dir.create(output_base_path, recursive = TRUE, showWarnings = FALSE)
    
    # Initialize variables to store global sets
    global_genes_seurat <- c()
    global_genes_scanpy <- c()
    
    # Get unique clusters
    clusters <- unique(markers2$cluster)
    
    # Loop over each cluster
    jaccard_indices <- list()  # Store per-cluster Jaccard indices
    for (clust in clusters) {
        # Subset data for the current cluster
        subset_data <- markers2[markers2$cluster == clust, ]
        
        # Get gene sets for Seurat and Scanpy thresholds
        genes_seurat <- subset_data$gene[abs(subset_data$logFC_r) >= 1 & subset_data$p_val_adj_r < 0.05]
        genes_scanpy <- subset_data$gene[abs(subset_data$logFC_py) >= 1 & subset_data$p_val_adj_py < 0.05]
        
        # Add to global sets
        global_genes_seurat <- union(global_genes_seurat, genes_seurat)
        global_genes_scanpy <- union(global_genes_scanpy, genes_scanpy)
        
        # Calculate Jaccard index for the current cluster
        jaccard_index <- calculate_jaccard(genes_seurat, genes_scanpy)
        jaccard_indices[[as.character(clust)]] <- jaccard_index
        
        data <- list(Seurat = genes_seurat, Scanpy = genes_scanpy)
        
        if (make_upset_plots) {
            p <- upset_plot_general(data, group1_name = "Seurat", group2_name = "Scanpy", comparison = "Gene", as_ggplot = FALSE, save = glue("{output_base_path_biological_output}/upset_volcano_clust{clust}.tiff"))  #!!!! uncomment   
        }
    }
    
    # Calculate the global Jaccard index across all clusters
    global_jaccard_index <- calculate_jaccard(global_genes_seurat, global_genes_scanpy)
    
    # Prepare the output file path
    output_file <- glue::glue("{output_base_path_biological_output}/jaccard_indices.txt")
    
    # Write results to the file
    output_lines <- c(
        "Jaccard Index for each cluster:",
        paste(names(jaccard_indices), jaccard_indices, sep = ": "),
        "",
        glue::glue("Global Jaccard Index across all clusters: {global_jaccard_index}")
    )
    writeLines(output_lines, con = output_file)
}