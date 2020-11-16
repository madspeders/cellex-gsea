#' @title gene set specificity
#' @author Mads Porse Pedersen
#' @usage Takes CELLEX data and an input gene set (in ensembl, uniprot, or hgnc format),
#' then extracts the ensembl gene names and performs a user-specified test to determine
#' gene specific expression for various cell types.
#' @return Output dataframe with p-values, z-values, and/or enrichment scores (es-values)
#' for the various cell types.

#' @importFrom magrittr %>%



#' @export
get_alias <- function(input,
                      input_type = c("ensembl", "uniprot", "gene"),
                      ensembl_dataset = "hsapiens_gene_ensembl",
                      ensembl_verbose = TRUE,
                      return_only_ensembl = FALSE
){

  loop <- TRUE

  while(loop){
    ### Regarding biomaRt ###
    try(
      if(input_type[1] == "ensembl"){
        ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = ensembl_dataset, verbose = ensembl_verbose)
        alias_df <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "uniprotswissprot"),
                                   filters = "ensembl_gene_id", values = input, mart = ensembl)
        alias_df <- alias_df %>%
          dplyr::arrange(desc(hgnc_symbol), desc(uniprotswissprot)) %>%
          dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
        missing <- input[!(input %in% alias_df$ensembl_gene_id)]
        alias_df <- alias_df %>% dplyr::add_row(data.frame(ensembl_gene_id = missing,
                                                           hgnc_symbol = rep("", length(missing)),
                                                           uniprotswissprot = rep("", length(missing))))
      } else if(input_type[1] == "gene"){
        ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = ensembl_dataset, verbose = ensembl_verbose)
        alias_df <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "uniprotswissprot"),
                                   filters = "hgnc_symbol", values = input, mart = ensembl)
        alias_df <- alias_df %>%
          dplyr::arrange(desc(hgnc_symbol), desc(uniprotswissprot)) %>%
          dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
        missing <- input[!(input %in% alias_df$hgnc_symbol)]
        alias_df <- alias_df %>% dplyr::add_row(data.frame(ensembl_gene_id = length(missing)),
                                                hgnc_symbol = rep("", missing,
                                                                  uniprotswissprot = rep("", length(missing))))
      } else if(input_type[1] == "uniprot"){
        ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = ensembl_dataset, verbose = ensembl_verbose)
        alias_df <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "uniprotswissprot"),
                                   filters = "uniprotswissprot", values = input, mart = ensembl)
        alias_df <- alias_df %>%
          dplyr::arrange(desc(hgnc_symbol), desc(uniprotswissprot)) %>%
          dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
        missing <- input[!(input %in% alias_df$uniprotswissprot)]
        alias_df <- alias_df %>% dplyr::add_row(data.frame(ensembl_gene_id = length(missing)),
                                                hgnc_symbol = rep("", length(missing)),
                                                uniprotswissprot = rep("", missing))
      }
    )


    if(exists("alias_df")){
      loop <- FALSE
    }

  }

  alias_df <- alias_df %>% dplyr::filter(ensembl_gene_id != "")

  if(return_only_ensembl){
    return(alias_df$ensembl_gene_id)
  } else {
    return(alias_df)
  }
}



#' @export
first_order_network <- function(input,
                                ppi_network = "InBio_Map_core_2016_09_12/core_v3.psimitab",
                                col_names = c("unique_A", "unique_B"),
                                ensembl_dataset = "hsapiens_gene_ensembl",
                                ensembl_verbose = TRUE,
                                return_only_ensembl = FALSE
){
  # Read in ppi data, generate graph and detect first order proteins.
  ppi_data <- readr::read_tsv(file = ppi_network, col_names = col_names)

  for(pro in input){
    if(!exists("ppi_list")){
      ppi_list <- ppi_data %>% dplyr::filter(str_detect(unique_A, pro))
    } else {
      ppi_temp <- ppi_data %>% dplyr::filter(str_detect(unique_A, pro))
      ppi_list <- dplyr::bind_rows(ppi_list, ppi_temp)
      rm(ppi_temp)
    }
  }
  ppi_data <- ppi_list
  rm(ppi_list)

  ppi_data$unique_A <- sub("uniprotkb:", "", ppi_data$unique_A)
  ppi_data$unique_B <- sub("uniprotkb:", "", ppi_data$unique_B)
  first_order_network <- unique(c(unique(ppi_data$unique_A), unique(ppi_data$unique_B)))



  ### Regarding biomaRt ###
  loop <- TRUE

  while(loop){
    try(
      ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = ensembl_dataset, verbose = ensembl_verbose)
    )

    if(exists("ensembl")){
      loop <- FALSE
    }
  }

  alias_df <- biomaRt::getBM(attributes = c('ensembl_gene_id', "uniprotswissprot", 'hgnc_symbol'),
                             filters = "uniprotswissprot", values = first_order_network, mart = ensembl)
  alias_df <- alias_df %>%
    dplyr::arrange(desc(hgnc_symbol), desc(uniprotswissprot), ensembl_gene_id) %>%
    dplyr::distinct(uniprotswissprot, .keep_all = TRUE) %>%
    dplyr::filter(ensembl_gene_id != "", uniprotswissprot != "", hgnc_symbol != "")
  missing <- first_order_network[!(first_order_network %in% alias_df$uniprotswissprot)]
  alias_df <- alias_df %>% dplyr::add_row(data.frame(ensembl_gene_id = rep("", length(missing))),
                                          hgnc_symbol = rep("", length(missing)),
                                          uniprotswissprot = missing)

  alias_df <- alias_df %>% dplyr::filter(ensembl_gene_id != "")

  if(return_only_ensembl){
    return(alias_df$ensembl_gene_id)
  } else {
    return(alias_df)
  }
}



#' @export
fastRndWalk <- function(gSetIdx, geneRanking, j, Ra) {
  m <- length(geneRanking)
  k <- length(gSetIdx)
  idxs <- sort.int(fastmatch::fmatch(gSetIdx, geneRanking))

  stepCDFinGeneSet2 <-
    sum(Ra[geneRanking[idxs], j] * (m - idxs + 1)) /
    sum((Ra[geneRanking[idxs], j]))

  stepCDFoutGeneSet2 <- (m * (m + 1) / 2 - sum(m - idxs + 1)) / (m - k)
  walkStat <- stepCDFinGeneSet2 - stepCDFoutGeneSet2

  return(walkStat)
}



#' @export
ssgsea <- function(data, geneSet, alpha=0.25, n_cores = 1){
  #data <- tibble::column_to_rownames(.data = data, var = "gene")
  p <- nrow(data)
  n <- ncol(data)

  R <- apply(data, 2, function(x, p) as.integer(rank(x, ties.method = "min")), p) #"average", "min", "max"
  Ra <- abs(R)^alpha

  gSetIdx <- which(rownames(data) %in% geneSet)

  es <- matrix(NA,
               nrow=1,
               ncol=ncol(data))

  es <- parallel::mclapply(as.list(1:n), mc.cores = n_cores, FUN = function(j) {
    geneRanking <- sort.list(R[,j], decreasing=TRUE)
    es_sample <- lapply(list(gSetIdx), fastRndWalk, geneRanking, j, Ra)

    unlist(es_sample)
  })

  es <- do.call("cbind", es)
  es <- matrix(es, nrow=1)

  # Normalization of values.
  es <- apply(es, 2, function(x, es) x / max(es), es)
  es <- matrix(es, nrow=1)

  # Output.
  es <- data.frame(tissue.cell = colnames(data), es.value = es[1,])

  return(es)
}



#' @export
statistical_test <- function(data,
                             subset,
                             background,
                             stat_test = c("W", "KS", "T"),
                             p_val = TRUE,
                             p_val_adjust = TRUE,
                             emp_p_val = FALSE,
                             es_val = FALSE,
                             plot = FALSE,
                             n_rep = 1000,
                             n_cores = 1,
                             n_background = 0){

  # Set seed.
  set.seed(42)
  data <- tibble::column_to_rownames(.data = data, var = "gene")
  names <- colnames(data)

  if(emp_p_val & n_background == 0){
    to_sample <- replicate(n_rep, {
      sample(x = rownames(data), size = length(subset))
    })

    if(stat_test[1] == "W"){
      stat_df_list <- parallel::mclapply(X = 1:n_rep, mc.cores = n_cores, FUN = function(i){
        stat_vec <- numeric(length(names))
        for(j in 1:length(stat_vec)){
          stat_vec[j] <- wilcox.test(x = data[to_sample[,i],j],
                                     y = data[!(rownames(data) %in% to_sample[,i]),j])$statistic
        }
        return(stat_vec)
      })

      if(n_rep > 1) {
        stat_df <- Reduce(rbind, stat_df_list)
        colnames(stat_df) <- names
        stat_df <- tibble::tibble(as.data.frame(stat_df))
      } else {
        stat_df <- t(Reduce(rbind, stat_df_list))
        colnames(stat_df) <- names
        stat_df <- tibble::tibble(as.data.frame(stat_df))
      }


    } else if(stat_test[1] == "T"){
      stat_df_list <- parallel::mclapply(X = 1:n_rep, mc.cores = n_cores, FUN = function(i){
        stat_vec <- numeric(length(names))
        for(j in 1:length(stat_vec)){
          stat_vec[j] <- t.test(x = data[(to_sample[,i]),j],
                                y = data[!(rownames(data) %in% to_sample[,i]),j])$statistic
        }
        return(stat_vec)
      })

      if(n_rep > 1) {
        stat_df <- Reduce(rbind, stat_df_list)
        colnames(stat_df) <- names
        stat_df <- tibble::tibble(as.data.frame(stat_df))
      } else {
        stat_df <- t(Reduce(rbind, stat_df_list))
        colnames(stat_df) <- names
        stat_df <- tibble::tibble(as.data.frame(stat_df))
      }

    } else if(stat_test[1] == "KS"){
      stat_df_list <- parallel::mclapply(X = 1:n_rep, mc.cores = n_cores, FUN = function(i){
        stat_vec <- numeric(length(names))
        for(j in 1:length(stat_vec)){
          stat_vec[j] <- ks.test(x = data[(to_sample[,i]),j],
                                 y = data[!(rownames(data) %in% to_sample[,i]),j])$statistic
        }
        return(stat_vec)
      })

      if(n_rep > 1) {
        stat_df <- Reduce(rbind, stat_df_list)
        colnames(stat_df) <- names
        stat_df <- tibble::tibble(as.data.frame(stat_df))
      } else {
        stat_df <- t(Reduce(rbind, stat_df_list))
        colnames(stat_df) <- names
        stat_df <- tibble::tibble(as.data.frame(stat_df))
      }

    }


  } else if(emp_p_val & n_background != 0) {
    # Use only a smaller selection of the background genes instead of all of them.
    # For instance, n_background = 1000.
    # Reduces accuracy of output the smaller the background gene set, though.

    to_sample <- replicate(n_rep, {
      sample(x = rownames(data), size = length(subset))
    })

    to_sample_back <- apply(X = to_sample, MARGIN = 2, FUN = function(x){
      gene_names <- rownames(data)[!(rownames(data) %in% x)]
      gene_names <- sample(x = gene_names, size = n_background)
      return(gene_names)
    })


    if(stat_test[1] == "W"){
      stat_df_list <- parallel::mclapply(X = 1:n_rep, mc.cores = n_cores, FUN = function(i){
        stat_vec <- numeric(length(names))
        for(j in 1:length(stat_vec)){
          stat_vec[j] <- wilcox.test(x = data[to_sample[,i],j],
                                     y = data[to_sample_back[,i],j])$statistic
        }
        return(stat_vec)
      })

      if(n_rep > 1) {
        stat_df <- Reduce(rbind, stat_df_list)
        colnames(stat_df) <- names
        stat_df <- tibble::tibble(as.data.frame(stat_df))
      } else {
        stat_df <- t(Reduce(rbind, stat_df_list))
        colnames(stat_df) <- names
        stat_df <- tibble::tibble(as.data.frame(stat_df))
      }

    } else if(stat_test[1] == "T"){
      stat_df_list <- parallel::mclapply(X = 1:n_rep, mc.cores = n_cores, FUN = function(i){
        stat_vec <- numeric(length(names))
        for(j in 1:length(stat_vec)){
          stat_vec[j] <- t.test(x = data[(to_sample[,i]),j],
                                y = data[to_sample_back[,i],j])$statistic
        }
        return(stat_vec)
      })

      if(n_rep > 1) {
        stat_df <- Reduce(rbind, stat_df_list)
        colnames(stat_df) <- names
        stat_df <- tibble::tibble(as.data.frame(stat_df))
      } else {
        stat_df <- t(Reduce(rbind, stat_df_list))
        colnames(stat_df) <- names
        stat_df <- tibble::tibble(as.data.frame(stat_df))
      }

    } else if(stat_test[1] == "KS"){
      stat_df_list <- parallel::mclapply(X = 1:n_rep, mc.cores = n_cores, FUN = function(i){
        stat_vec <- numeric(length(names))
        for(j in 1:length(stat_vec)){
          stat_vec[j] <- ks.test(x = data[(to_sample[,i]),j],
                                 y = data[to_sample_back[,i],j])$statistic
        }
        return(stat_vec)
      })

      if(n_rep > 1) {
        stat_df <- Reduce(rbind, stat_df_list)
        colnames(stat_df) <- names
        stat_df <- tibble::tibble(as.data.frame(stat_df))
      } else {
        stat_df <- t(Reduce(rbind, stat_df_list))
        colnames(stat_df) <- names
        stat_df <- tibble::tibble(as.data.frame(stat_df))
      }

    }

  }


  if(p_val | emp_p_val) {
    # Generate analytical test statistics.
    output <- parallel::mclapply(names, mc.cores = n_cores, FUN = function(col){
      if(stat_test[1] == "W"){
        test <- wilcox.test(x = data[subset, col],
                            y = data[background, col])
        df.test <- broom::tidy(test)
      } else if(stat_test[1] == "T"){
        test <- t.test(x = data[subset, col],
                       y = data[background, col])
        df.test <- broom::tidy(test)
      } else if(stat_test[1] == "KS"){
        test <- ks.test(x = data[subset, col],
                        y = data[background, col])
        df.test <- broom::tidy(test)
      }

      df.test <- df.test %>% tibble::add_column(tissue.cell = col, .before = "statistic")
      return(df.test)

    })

    output <- tibble::tibble(Reduce(rbind, output))

    if(p_val){
      output <- output %>% dplyr::select(tissue.cell, statistic, p.value)
    } else {
      output <- output %>% dplyr::select(tissue.cell, statistic)
    }

  }



  if(emp_p_val){
    p_z_val_out <- parallel::mclapply(colnames(stat_df), mc.cores = n_cores, FUN = function(col){
      stat_an <- output %>% dplyr::filter(tissue.cell == col) %>% dplyr::select(statistic)
      stat_an <- stat_an[[1]]
      stat_null <- stat_df[,col]
      stat_null <- stat_null[[1]]
      #p_value <- (sum(stat_null > stat_an) + 1) / (n_rep) # one-sided
      p_value <- ( sum(stat_null < -abs(stat_an)) + sum(stat_null > abs(stat_an)) ) / (n_rep) # two-sided
      z_value <- (stat_an - mean(stat_null)) / sd(stat_null)
      return(c(p_value, unname(z_value)))
    })

    p_z_val_out <- tibble::tibble(as.data.frame(Reduce(rbind, p_z_val_out)))
    names(p_z_val_out) <- c("p.value", "Z.value")
    p_val_out <- p_z_val_out$p.value
    z_val_out <- p_z_val_out$Z.value

    p_val_out[p_val_out > 1] <- 1.0
    p_val_out[p_val_out == 0] <- 1 / n_rep
    names(p_val_out) <- output$tissue.cell
    names(z_val_out) <- output$tissue.cell

  }


  # Create output data frame and sort by adjusted P-values (multiple testing).
  if(p_val_adjust & p_val){
    output$p.value <- p.adjust(output$p.value, method = "fdr")
  }

  if(emp_p_val){
    output <- output %>% tibble::add_column(p.value.emp = p_val_out)
    output <- output %>% tibble::add_column(z.value = z_val_out, .after = "p.value.emp")
  }

  # Calculate ES values (from ssGSEA method) and append to output.
  if(es_val & !(p_val | emp_p_val)) {
    output <- tibble::tibble(ssgsea(data, geneSet = subset, n_cores = n_cores))
  } else if(es_val) {
    es <- ssgsea(data, geneSet = subset, n_cores = n_cores)
    output <- output %>% tibble::add_column(es.value = es$es.value)
  }

  # Sort the output values.
  if(es_val & p_val & emp_p_val) {
    output <- output %>% dplyr::arrange(desc(es.value), p.value, desc(z.value))
  } else if(es_val & p_val) {
    output <- output %>% dplyr::arrange(desc(es.value), p.value)
  } else if(es_val & emp_p_val) {
    output <- output %>% dplyr::arrange(desc(es.value), desc(z.value))
  } else if(emp_p_val & p_val) {
    output <- output %>% dplyr::arrange(p.value, desc(z.value))
  } else if(emp_p_val) {
    output <- output %>% dplyr::arrange(p.value.emp, desc(z.value))
  } else if(p_val) {
    output <- output %>% dplyr::arrange(p.value)
  } else if(es_val) {
    output <- output %>% dplyr::arrange(desc(es.value))
  }


  if(plot & emp_p_val){
    # Plot the calculated statistic(s) as a histogram and put a line
    # where the analytical value is - only for liver cells (hepatocytes).
    # NB: "stat_df" is necessary for this!
    tissue <- output$tissue.cell[1]
    intercept <- output %>% dplyr::filter(tissue.cell == tissue) %>% dplyr::select(statistic)
    intercept <- unname(intercept[[1]])

    ggplot2::ggplot(data = stat_df, aes_string(x = tissue)) +
      ggplot2::geom_histogram() +
      ggplot2::geom_vline(xintercept = intercept, color = "red", linetype = "longdash", size = 1) +
      ggplot2::xlab(paste0(stat_test[1], " statistic")) +
      ggplot2::ggtitle(paste0("Histogram of statistic values (", tissue, ")")) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = element_text(size = 12, hjust = 0.5))
    ggplot2::ggsave("stat_histogram.png")
  }


  if(p_val | emp_p_val) {
    # Remove the "statistic" column.
    output <- output %>% dplyr::select(-statistic)
    # Return the output dataframe.
    return(output)
  } else if(es_val) {
    # Return the output dataframe.
    return(output)
  } else {
    return(NULL)
  }
}



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





#' Perform Enrichment Analysis
#'
#' This function takes as input a gene set (in ENSEMBL, UNIPROT, or HNGC format) and performs
#' enrichment analysis / gene specificity analysis and outputs a dataframe with p-values and/or
#' scores indicating how specifically expressed the gene set is in different cell types. The
#' p-values are calculated by comparing the CELLEX values for the input gene set to the CELLEX values
#' for the background gene set.
#'
#' @param input_set Gene set or protein set, as a character vector.
#' @param input_type What format the input_set is in. One of "ensembl" (ENSEMBL format), "uniprot" (UNIPROT format), or "gene" (HGNC format). Defaults to "ensembl" if no user specified value is selected.
#' @param cellex_data The CELLEX dataset to use in the package - must be in .csv format, and can also be gzip compressed as well. The path to a user-provided CELLEX dataset can be put here. Otherwise, two CELLEX datasets are included in the package, which can also be used. Set "cellex_data = 1" (default) for tabula_muris cellex dataset, set "cellex_data = 2" for gtex_v8 cellex dataset.
#' @param first_order If TRUE, obtain the first order network for the input genes / proteins, and use for the analysis. Defaults to FALSE.
#' @param all_genes_as_background If TRUE, use all genes in the CELLEX dataset as background for analysis. If FALSE, use all genes in the CELLEX dataset which are not in the input_set, as background for analysis.
#' @param statistic Which test statistic to use. One of "W" (Wilcoxon test), "KS" (Kolmogorov–Smirnov test), or "T" (Student's t-test). Defaults to "W".
#' @param p_value If TRUE, calculates and outputs the p-value (for the selected test statistic) for the expression specificity analysis. Defaults to TRUE.
#' @param p_value_adjust If TRUE, adjusts the calculated p-values for multiple testing, using the FDR approach. Defaults to TRUE.
#' @param emp_p_value If TRUE, calculates empirical p-values by randomly sampling a number of genes from the CELLEX dataset equivalent to the number of genes in the input_set, which is repeated a number of times (see parameter "reps") to obtain a null distribution. The null distribution is then compared to the test statistic for the genes in the input_set. Related parameters: "reps", "num_cores", "num_background_genes", "statistic_plot".
#' @param es_value If TRUE, calculates enrichment scores for the input_set for each cell type / tissue type in the CELLEX. Uses ssGSEA approach implemented in the R package GSVA.
#' @param reps Default is 1000. The number of random samplings to perform when computing the null distribution that is used to calculate the empirical p-value (see parameter "emp_p_value").
#' @param num_cores Default is 1. The number of cores/threads to use when computing the null distribution for the empirical p-value (see parameter "emp_p_value").
#' @param num_background_genes Default is 0 (which means to use all background genes). The number of background genes to use when computing the null distribution for the empirical p-value (see parameter "emp_p_value").
#' @param statistic_plot Default is FALSE. Whether to plot (TRUE) or not plot (FALSE) the null distribution and the analytical test statistic for the most significant cell type.
#' @param save_output Default is TRUE. Whether to save (TRUE) or not save (FALSE) the output of the tool. If TRUE, the tool saves the output dataframe as a .csv file.
#'
#' @export
cellex_analysis <- function(input_set, # input gene set or protein set.
                            input_type = c("ensembl", "uniprot", "gene"),
                            cellex_data = c(1, 2),
                            p_value = TRUE,
                            p_value_adjust = TRUE,
                            emp_p_value = FALSE,
                            es_value = FALSE,
                            first_order = FALSE,
                            all_genes_as_background = FALSE,
                            statistic = c("W", "KS", "T"),
                            reps = 1000,
                            num_cores = 1,
                            num_background_genes = 0,
                            statistic_plot = FALSE,
                            save_output = TRUE
                            ) {

  set.seed(42)

  # Obtain aliases for the input gene set.
  #source("gene_set.R")
  message("Using BioMart to obtain ENSEMBL, UNIPROT and HNGC names...")
  gene_set <- get_alias(input = input_set, input_type = input_type)
  message("Done!")

  # Get PPI data matching the protein names from the alias dataframe.
  # Also get expression values for each ES?.

  if(first_order){
    #source("ppi.R")
    message("\nObtaining first order network...")
    proteins <- unique(gene_set$uniprotswissprot)
    proteins <- proteins[!(proteins == "")]
    gene_set <- first_order_network(input = proteins)
    message("Done!")
  }

  # Load cellex data.
  message("\nLoading CELLEX data...")
  if(cellex_data[1] == 1) {
    cellex_data <- system.file("extdata", "tabula_muris.esmu.csv.gz", package = "cellex.analysis")
    cellex_data <- readr::read_csv(gzfile(cellex_data))
  } else if(cellex_data[1] == 2) {
    cellex_data <- system.file("extdata", "gtex_v8_filtered_normLog_all.esmu.csv.gz", package = "cellex.analysis")
    cellex_data <- readr::read_csv(gzfile(cellex_data))
  } else {
    if((substr(cellex_data, start = nchar(cellex_data)-2, stop = nchar(cellex_data)) == ".gz") | (substr(cellex_data, start = nchar(cellex_data)-2, stop = nchar(cellex_data)) == ".GZ")){
      cellex_data <- readr::read_csv(gzfile(cellex_data))
    } else {
      cellex_data <- readr::read_csv(cellex_data)
    }
  }

  # Choose subset genes and background genes.
  subset_genes <- cellex_data$gene[cellex_data$gene %in% gene_set$ensembl_gene_id]
  if(all_genes_as_background){
    background_genes <- cellex_data$gene
  } else {
    background_genes <- cellex_data$gene[!(cellex_data$gene %in% gene_set$ensembl_gene_id)]
  }

  gene_set <- gene_set %>% dplyr::filter(ensembl_gene_id %in% subset_genes)
  message("Done!")

  # Perform a statistical test for each tissue in the cellex data.
  # Possibility to perform permutation tests as well.
  # Load statistical tests functions.
  #source("statistical_tests.R")

  # Perform statistical tests.
  message("\nRunning analysis...")
  output <- statistical_test(data = cellex_data,
                             subset = subset_genes,
                             background = background_genes,
                             stat_test = statistic[1],
                             p_val = p_value,
                             p_val_adjust = p_value_adjust,
                             emp_p_val = emp_p_value,
                             es_val = es_value,
                             plot = statistic_plot,
                             n_rep = reps,
                             n_cores = num_cores,
                             n_background = num_background_genes)
  message("Done!")

  if(is.null(output)){
    message("\nERROR!\nNo value or score was selected for computation, so no output was generated.\nExiting...")
  } else {
    # Whether to save the output to a file.
    if(save_output){
      message("\nSaving output...")
      #saveRDS(output, file = "output.rds")
      readr::write_csv(output, file = "output.csv")
    }
    message("Done!")

    message("\nSaving plots...")
    #source("plotting.R")

    for(value in c("es.value", "z.value", "p.value", "p.value.emp")){
      if(value %in% colnames(output)){
        plot_bar(output, value)
      }
    }
    plot_box(output, subset_genes, cellex_data)

    message("Done!")

    message("\nOutput:\n")
    output

    return(output)
  }
}
