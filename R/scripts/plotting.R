#' @title gene set specificity - functions to use for tool 1
#' @author Mads Porse Pedersen

#' @export
plot_bar <- function(output,
                     value_to_plot,
                     n_tissues = 10,
                     font_size = 8){
  if(value_to_plot == "es.value") {
    ggplot2::ggplot(data = head(output, n_tissues), mapping = ggplot2::aes(x = reorder(tissue.cell, -es.value), y = es.value)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::xlab("Cell type") +
      ggplot2::ylab(value_to_plot) +
      ggplot2::ggtitle(paste0("Bar plot of cell type specific expression (top ", n_tissues, ") - ", value_to_plot)) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
                     axis.text.x = ggplot2::element_text(angle = 270, vjust = 0.5, hjust = 0, size = font_size))
  } else if(value_to_plot == "z.value") {
    ggplot2::ggplot(data = head(output, n_tissues), mapping = ggplot2::aes(x = reorder(tissue.cell, -z.value), y = z.value)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::xlab("Cell type") +
      ggplot2::ylab(value_to_plot) +
      ggplot2::ggtitle(paste0("Bar plot of cell type specific expression (top ", n_tissues, ") - ", value_to_plot)) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
                     axis.text.x = ggplot2::element_text(angle = 270, vjust = 0.5, hjust = 0, size = font_size))
  } else if(value_to_plot == "p.value") {
    ggplot2::ggplot(data = head(output, n_tissues), mapping = ggplot2::aes(x = reorder(tissue.cell, -(-log10(p.value))), y = -log10(p.value))) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::xlab("Cell type") +
      ggplot2::ylab(paste0("-log10(", value_to_plot, ")")) +
      ggplot2::ggtitle(paste0("Bar plot of cell type specific expression (top ", n_tissues, ") - ", value_to_plot)) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
                     axis.text.x = ggplot2::element_text(angle = 270, vjust = 0.5, hjust = 0, size = font_size))
  } else if(value_to_plot == "p.value.emp") {
    ggplot2::ggplot(data = head(output, n_tissues), mapping = ggplot2::aes(x = reorder(tissue.cell, -(-log10(p.value.emp))), y = -log10(p.value.emp))) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::xlab("Cell type") +
      ggplot2::ylab(paste0("-log10(", value_to_plot, ")")) +
      ggplot2::ggtitle(paste0("Bar plot of cell type specific expression (top ", n_tissues, ") - ", value_to_plot)) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
                     axis.text.x = ggplot2::element_text(angle = 270, vjust = 0.5, hjust = 0, size = font_size))
  }
  ggplot2::ggsave(paste0("output_barplot_", value_to_plot,".png"))
}

#' @export
plot_box <- function(output, subset, cellex_data, n_tissues = 10, param = "ESmu", font_size = 8){ # background
  top_cells <- head(output, n_tissues)$tissue.cell
  subset_new <- cellex_data %>%
    dplyr::filter(gene %in% subset) %>%
    dplyr::select(gene, top_cells) %>%
    tidyr::gather(cell_type, value, -gene)
  ggplot2::ggplot(data = subset_new, mapping = ggplot2::aes(x = reorder(cell_type, -value), y = value)) +
    ggplot2::geom_boxplot() +
    ggplot2::ggtitle(paste0("Boxplot of gene set ", param," values")) +
    ggplot2::xlab("cell type") +
    ggplot2::ylab(param) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
                   axis.text.x = ggplot2::element_text(angle = 270, vjust = 0.5, hjust = 0, size = font_size))
  # background_new <- cellex_data %>%
  #   dplyr::filter(gene %in% background) %>%
  #   dplyr::select(gene, top_cells) %>%
  #   tidyr::gather(cell_type, value, -gene)
  # p2 <- ggplot2::ggplot(data = background_new, mapping = ggplot2::aes(x = reorder(cell_type, -value), y = value)) +
  #   ggplot2::geom_boxplot() +
  #   ggplot2::ggtitle(paste0("Boxplot of background genes ", param," values")) +
  #   ggplot2::xlab("cell type") +
  #   ggplot2::ylab(param) +
  #   ggplot2::theme_bw() +
  #   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
  #                  axis.text.x = ggplot2::element_text(angle = 270, vjust = 0.5, hjust = 0, size = font_size))
  # gridExtra::grid.arrange(p1, p2, nrow = 1)
  ggplot2::ggsave("output_boxplot.png")
}


# plot_heatmap <- function(output, cellex_data, n_tissues = 10){
#   top_tissues <- head(output, n_tissues)$tissue.cell
#   top_tissues <- data %>% dplyr::select(top_tissues)
#   top_tissues <- as.matrix(top_tissues)
#   rownames(top_tissues) <- data$gene
#
#   pheatmap::pheatmap(top_tissues, main = paste0("Heatmap of CELLEX mu values for top ", n_tissues, " tissues"),
#                      cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, scale = "none")
# }


# plot_ppi_graph <- function(ppi_data_graph, gene_expr){
#   unique_proteins <- unique(c(ppi_data_graph$unique_A, ppi_data_graph$unique_B))
#   gene_expr_filter <- gene_expr %>% dplyr::filter(uniprotswissprot %in% unique_proteins)
#   g <- igraph::graph_from_data_frame(ppi_data_graph, vertices = gene_expr_filter, directed = FALSE)
#   #g <- igraph::simplify(g)
#   expr_vec <- gene_expr_filter$expr_val
#   expr_vec <- (expr_vec) / max(expr_vec)
#
#   color_fun <- grDevices::colorRampPalette(c("red", "blue"))
#   color_vec <- color_fun(1000000)
#   expr_vec <- round(expr_vec*1000000)
#   expr_vec <- sapply(expr_vec, FUN = function(x){
#     if(x != 1000000){
#       return(x+1)
#     } else {
#       return(x)
#     }
#   })
#   expr_vec <- unlist(expr_vec, use.names = FALSE)
#   igraph::V(g)$color <- color_vec[expr_vec]
#
#   igraph::V(g)$label.cex = 0.7
#   #set.seed(1)
#   igraph::plot.igraph(g, label.dist = 1, vertex.size = 5,
#                       vertex.label = NA, main = "iGraph plot")
# }
