download_box <- function(file_link, filename, data_path_root, file_ext = '.tar.gz') {
    if (!grepl("shared", file_link)) {
        if (file_ext == ".tar.gz") {
            file_ext_box_link <- ".gz"
        } else {
            file_ext_box_link <- file_ext
        }
        file_id <- sub("https://caltech.box.com/s/(.*)", "\\1", file_link)
        file_link <- paste0("https://caltech.box.com/shared/static/", file_id, file_ext_box_link)
    }
    
    output_path <- paste(data_path_root, filename, sep = "/")
    
    if (!grepl(".{0,10}\\.[^\\.]+$", filename)) {
        output_path_with_extension <- paste0(output_path, file_ext)
    } else {
        output_path_with_extension <- output_path
    }
    
    if (!dir.exists(output_path)) {
        dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
    }
    
    download.file(file_link, destfile = output_path_with_extension, mode = "wb")
    
    if (file_ext == '.tar.gz') {
        untar(tarfile = output_path_with_extension, exdir = dirname(output_path))
    }
    
    return (output_path)
}

read_count_output_modified <- function(dir, name, unspliced = FALSE, batch = FALSE, tcc = FALSE) {
    dir <- normalizePath(dir, mustWork = TRUE)
    if (unspliced) {
        m <- readMM(paste0(dir, "/", name, ".total.mtx"))
    } else {
        if (file.exists(paste0(dir, "/", name, ".cell.mtx"))) {
            # If the first file exists, read from it
            m <- readMM(paste0(dir, "/", name, ".cell.mtx"))
        } else {
            m <- readMM(paste0(dir, "/", name, ".mtx"))
        }
    }
    m <- Matrix::t(m)
    m <- as(m, "dgCMatrix")
    # The matrix read has cells in rows
    ge <- if (tcc) ".ec.txt" else ".genes.txt"
    con_genes <- file(paste0(dir, "/", name, ge))
    if (batch) {
        con_bcs <- file(paste0(dir, "/", name, ".barcodes.combined.txt"))
    } else {
        con_bcs <- file(paste0(dir, "/", name, ".barcodes.txt"))
    }
    genes <- readLines(con_genes)
    barcodes <- readLines(con_bcs)
    colnames(m) <- barcodes
    rownames(m) <- genes
    close(con_genes)
    close(con_bcs)
    return(m)
}


get_pc_diffs <- function(mat1, mat2, npcs = 20) {
    diffs <- numeric(npcs)
    for (i in seq_len(npcs)) {
        pc1 <- mat1[, i]
        pc2 <- mat2[, i]
        # diffs[i] <- abs(abs(sum(pc1 * pc2) / sqrt(sum(pc1^2)*sum(pc2^2))) - 1)

        # Cosine similarity
        cos_sim <- abs(sum(pc1 * pc2) / (sqrt(sum(pc1^2) * sum(pc2^2))))

        # Sine similarity
        sine_sim <- sqrt(1 - cos_sim^2)

        diffs[i] <- abs(sine_sim)
    }
    diffs
}


make_pc_diffs_df <- function(loadings_list, npcs) {
    # Generate all pairwise combinations of the loadings_list names
    names_combinations <- combn(names(loadings_list), 2, simplify = FALSE)

    # Create pairwise comparison data frames
    dfs <- lapply(names_combinations, function(x) {
        loadings1 <- loadings_list[[x[1]]]
        loadings2 <- loadings_list[[x[2]]]
        genes_use <- intersect(rownames(loadings1), rownames(loadings2))
        diffs <- get_pc_diffs(loadings1[genes_use, ], loadings2[genes_use, ], npcs = npcs)
        tibble(
            differences = diffs,
            PC = seq_len(npcs),
            type = paste(x[1], x[2], sep = " vs. ")
        )
    })

    # Combine all data frames into one
    bind_rows(dfs)
}


make_pairwise_df <- function(data_df) {
    # Get all pairwise combinations of the columns
    column_names <- colnames(data_df)
    combs <- t(combn(column_names, 2))

    # Create pairwise comparison data frames
    dfs <- apply(combs, 1, function(x) {
        tibble(
            value1 = data_df[[x[1]]],
            value2 = data_df[[x[2]]],
            package1 = x[1],
            package2 = x[2]
        )
    })

    # Combine all data frames into one
    bind_rows(dfs)
}


# Remove the lower left panel in pairwise comparison plots with 2x2 panels
# Or equal number of rows and columns
# Not that flexible, just make sure that that panel is empty by using fct_relevel
# Returns grob
rm_bl_panel <- function(p, ncol = 2) {
    # https://stackoverflow.com/a/49525552/8916916
    grob <- ggplotGrob(p)
    # Remove lower triangle, basically i > j
    inds_rm <- expand_grid(
        i = seq_len(ncol),
        j = seq_len(ncol)
    ) |>
        filter(i > j)
    panels_rm <- paste("panel", inds_rm$i, inds_rm$j, sep = "-")
    inds <- which(grob$layout$name %in% panels_rm)
    for (i in inds) {
        grob$grobs[[i]] <- nullGrob()
    }
    # Move the bottom axes
    indsj <- unique(inds_rm$j)
    for (i in indsj) {
        axis_b_rm <- paste("axis-b", i, sep = "-")
        k <- which(grob$layout$name == axis_b_rm)
        grob$layout[k, c("t", "b")] <- grob$layout[k, c("t", "b")] - 2 * (ncol - i)
    }
    indsi <- unique(inds_rm$i)
    for (i in indsi) {
        axis_l_rm <- paste("axis-l", i, sep = "-")
        k <- which(grob$layout$name == axis_l_rm)
        grob$layout[k, c("l", "r")] <- grob$layout[k, c("l", "r")] + 2 * (i - 1)
    }
    grob
}

calculate_angle <- function(v1, v2) {
    dot_product <- sum(v1 * v2)
    
    magnitude_vector1 <- sqrt(sum(v1^2))
    magnitude_vector2 <- sqrt(sum(v2^2))
    
    cos_theta <- dot_product / (magnitude_vector1 * magnitude_vector2)
    
    if (isTRUE(all.equal(cos_theta, 1))) {
        angle_radians <- 0
    } else if (isTRUE(all.equal(cos_theta, -1))) {
        angle_radians <- pi
    } else {
        angle_radians <- acos(cos_theta)
    }
    
    return (angle_radians)
}

flip_pcs <- function(v1, v2, pcs) {
    for (i in pcs) {
        # v <- c(
        #     max(v1[, i] - v2[, i]),
        #     max(v1[, i] + v2[, i])
        # )
        # if (which.min(v) == 2L) { # do flip
        #     v2[, i] <- -v2[, i]
        # }
        
        theta1_2 <- calculate_angle(v1[, i], v2[, i])
        theta1_neg2 <- calculate_angle(v1[, i], -v2[, i])
        
        if (theta1_neg2 < theta1_2) {
            v2[, i] <- -v2[, i]
        }
    }
    v2
}


make_pca_emb_df <- function(embeddings, package, pcs) {
    df <- as_tibble(embeddings[, pcs]) |>
        mutate(
            package = package,
            ID = seq_len(nrow(embeddings))
        )
    names(df)[1:2] <- c("x", "y")
    df
}


mat2list <- function(m, get_wts = FALSE) {
    apply(m, 1, function(x) {
        inds <- which(x > 0)
        if (get_wts) {
            out <- x[inds]
        } else {
            out <- inds
        }
        out
    })
}


find_jaccards <- function(ll) {
    # ll is a list of two sets of neighborhood lists, must be named.

    # Extract the two lists
    nb1 <- ll[[1]]
    nb2 <- ll[[2]]

    # Calculate the Jaccard index
    jaccard_values <- vapply(seq_along(nb1), function(i) {
        l1 <- unclass(nb1[[i]])
        l2 <- unclass(nb2[[i]])
        length(intersect(l1, l2)) / length(union(l1, l2))
    }, FUN.VALUE = numeric(1))

    # Create a data frame for the output
    out <- data.frame(
        type = paste(names(ll)[1], "vs.", names(ll)[2]),
        Jaccard = jaccard_values
    )

    out
}


calculate_knn_jaccards <- function(knn1_ids, knn2_ids) {
    jaccards_all_cells <- c()
    for (i in 1:nrow(knn1_ids)) {
        row1 <- as.integer(knn1_ids[i, ])
        row2 <- as.integer(knn2_ids[i, ])
        intersection <- length(intersect(row1, row2))
        union <- length(union(row1, row2))
        jaccard <- intersection / union
        jaccards_all_cells <- c(jaccards_all_cells, jaccard) # Append jaccard value
    }
    return(jaccards_all_cells)
}


find_wts_corrs <- function(ll_wts, ll_inds, combs, method = "pearson") {
    out <- apply(combs, 2, function(x) {
        nb1 <- ll_inds[[x[1]]]
        nb2 <- ll_inds[[x[2]]]
        wts1 <- ll_wts[[x[1]]]
        wts2 <- ll_wts[[x[2]]]
        # they should have the same length, shouldn't be singletons
        vapply(seq_along(nb1), function(i) {
            l1 <- unclass(nb1[[i]])
            l2 <- unclass(nb2[[i]])
            inds <- intersect(l1, l2)
            inds1 <- which(l1 %in% inds)
            inds2 <- which(l2 %in% inds)
            cor(wts1[[i]][inds1], wts2[[i]][inds2], method = method)
        }, FUN.VALUE = numeric(1))
    }, simplify = FALSE)
    out_names <- paste(combs[1, ], combs[2, ], sep = " vs. ")
    names(out) <- out_names
    out <- as.data.frame(out, optional = TRUE) |>
        pivot_longer(everything(), names_to = "type", values_to = "cor")
    out
}


get_alluvial_df <- function(df) {
    # Each column of df is clustering results from one package
    df <- df |>
        mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
        group_by_all() |>
        dplyr::count(name = "value")
    gather_set_data(df, 1:2)
}


# Code from https://gitlab.svi.edu.au/biocellgen-public/mage_2020_marker-gene-benchmarking/-/blob/master/code/run-scanpy.R
get_py_de_results <- function(adata_name) {
    name_lookup <- list(
        # FIXME: Not all are strictly t-values...
        names = "gene",
        scores = "scaled_statistic",
        pvals = "p_value",
        pvals_adj = "p_value_adj",
        logfoldchanges = "log_fc",
        pts_rec_array = "pts",
        pts_rest_rec_array = "pts_rest"
    )
    
    if (py_eval(glue("'pts_rec_array' not in {adata_name}.uns['rank_genes_groups']"))) {
        name_lookup$pts_rec_array <- NULL
    }
    
    if (py_eval(glue("'pts_rest_rec_array' not in {adata_name}.uns['rank_genes_groups']"))) {
        name_lookup$pts_rest_rec_array <- NULL
    }
    
    out <- list()
    for (ele in names(name_lookup)) {
        # We have to do this on the 'python side' to prevent array type 20 errors.
        py_run_string(
            glue("{ele} = pd.DataFrame({adata_name}.uns['rank_genes_groups']['{ele}'])")
        )
        out <- c(
            out,
            list(tidyr::pivot_longer(
                py[[ele]],
                cols = everything(),
                names_to = "cluster",
                values_to = name_lookup[[ele]]
            ))
        )
    }

    cluster <- tibble::tibble(cluster = out[[1]]$cluster)
    out <- lapply(out, function(x) dplyr::select(x, -cluster))

    result <- dplyr::bind_cols(cluster, !!!out) |>
        group_by(cluster) |>
        mutate(rank_py = seq_along(gene))

    # In the past the output from Scanpy grouped the results for each cluster
    # together with genes ranked. At some point however this changed for an
    # unknown reason, so we rectify this here. We DO NOT change the ranking of
    # scanpy.
    result <- dplyr::arrange(result, cluster)
    result
}


reorder_clusters_descending <- function(clusters) {
    # Count the size of each cluster
    cluster_sizes <- table(clusters)

    # Sort the clusters by size and get the ordered names
    ordered_cluster_names <- names(sort(cluster_sizes, decreasing = TRUE))

    # Create a mapping from old to new cluster numbers
    cluster_mapping <- setNames(seq_along(ordered_cluster_names), ordered_cluster_names)

    # Apply the mapping to renumber the clusters
    renumbered_clusters <- cluster_mapping[as.character(clusters)]

    # Convert the factor to numeric
    renumbered_clusters_factor <- factor(as.numeric(as.character(renumbered_clusters)))

    return(renumbered_clusters_factor)
}


increment_if_zeros <- function(clus_df_gather, column) {
    clus_df_gather <- clus_df_gather %>% mutate(group_numeric = as.numeric(as.character(.data[[column]])))

    if (any(clus_df_gather$group_numeric == 0, na.rm = TRUE)) {
        clus_df_gather$group_numeric <- clus_df_gather$group_numeric + 1
        clus_df_gather <- clus_df_gather %>% mutate(!!column := factor(group_numeric))
    }

    clus_df_gather <- clus_df_gather %>% select(-group_numeric)

    return(clus_df_gather)
}


sort_clusters_by_agreement <- function(clus_df_gather, stable_column = "Seurat", reordered_column = "Scanpy") {
    for (n in 1:2) {
        reordered_column_original_clusters_name <- paste0(reordered_column, "_original_clusters")
        
        clus_df_gather <- increment_if_zeros(clus_df_gather, stable_column)
        clus_df_gather <- increment_if_zeros(clus_df_gather, reordered_column)
        clus_df_gather <- increment_if_zeros(clus_df_gather, "y")
        
        # Initialize variables
        half_rows <- nrow(clus_df_gather) / 2
        
        subset_data <- clus_df_gather %>%
            ungroup() %>%
            dplyr::slice((half_rows + 1):nrow(clus_df_gather))
        
        subset_data <- subset_data %>% mutate(
            !!reordered_column_original_clusters_name := as.numeric(as.character(.data[[reordered_column]])),
            !!reordered_column := -as.numeric(as.character(.data[[reordered_column]])),
            y := -as.numeric(as.character(y)),
            best_cluster_agreement := .data[[reordered_column]]
        )
        
        # Loop over each unique cluster number in Group 2
        for (cluster_number in sort(unique(subset_data[[reordered_column_original_clusters_name]]))) {
            # Subset the data for the current cluster number
            subset_data2 <- subset_data[subset_data[[reordered_column_original_clusters_name]] == cluster_number, ]
            
            # Find the row with the largest overlap (value)
            largest_overlap_row <- subset_data2[which.max(subset_data2$value), ]
            
            # Check if the corresponding Group 1 number is available
            new_cluster_number <- as.numeric(as.character(largest_overlap_row[[stable_column]]))
            
            best_cluster_number <- new_cluster_number
            
            subset_data$best_cluster_agreement[subset_data[[reordered_column_original_clusters_name]] == cluster_number] <- best_cluster_number
            
            any(subset_data[[reordered_column]] == best_cluster_number)
            
            while (any(subset_data[[reordered_column]] == new_cluster_number)) {
                new_cluster_number <- new_cluster_number + 1
                if (!any(subset_data$best_cluster_agreement[subset_data[[reordered_column]] == new_cluster_number] <= best_cluster_number)) {
                    subset_data[[reordered_column]] <- ifelse(subset_data[[reordered_column]] >= new_cluster_number, subset_data[[reordered_column]] + 1, subset_data[[reordered_column]])
                    subset_data$y <- ifelse(subset_data$y >= new_cluster_number, subset_data$y + 1, subset_data$y)
                    break
                }
            }
            
            # Assign the new cluster number
            subset_data[[reordered_column]][subset_data[[reordered_column_original_clusters_name]] == cluster_number] <- new_cluster_number
            
            for (i in (new_cluster_number + 1):(max(as.numeric(as.character(clus_df_gather$y))))) {
                if (!any(subset_data[[reordered_column]] == i)) {
                    subset_data[[reordered_column]] <- ifelse(subset_data[[reordered_column]] > i, subset_data[[reordered_column]] - 1, subset_data[[reordered_column]])
                    subset_data$y <- ifelse(subset_data$y > i, subset_data$y - 1, subset_data$y)
                    break
                }
            }
        }
        
        mapping <- setNames(
            seq_along(sort(unique(subset_data[[reordered_column]]))),
            sort(unique(subset_data[[reordered_column]]))
        )
        
        subset_data <- subset_data %>% mutate(!!reordered_column := mapping[as.character(.data[[reordered_column]])])
        
        clus_df_gather[[reordered_column]][(1:half_rows)] <- subset_data[[reordered_column]]
        clus_df_gather[[reordered_column]][((half_rows + 1):nrow(clus_df_gather))] <- subset_data[[reordered_column]]
        clus_df_gather$y[((half_rows + 1):nrow(clus_df_gather))] <- subset_data[[reordered_column]]
        
        sorted_levels <- sort(as.numeric(levels(clus_df_gather[[reordered_column]])))
        sorted_levels <- as.character(sorted_levels)
        clus_df_gather[[reordered_column]] <- factor(clus_df_gather[[reordered_column]], levels = sorted_levels)
        
        sorted_levels <- sort(as.numeric(levels(clus_df_gather$y)))
        sorted_levels <- as.character(sorted_levels)
        clus_df_gather$y <- factor(clus_df_gather$y, levels = sorted_levels)
    }
        
    return(clus_df_gather)
}
