#' @title CELLEX-GSEA
#' @author Mads Porse Pedersen
#' @usage Takes CELLEX data and an input gene set (in ensembl, uniprot, or hgnc format),
#' then extracts the ensembl gene names and performs a user-specified statistical test
#' to determine gene specific expression for various cell types.
#' @return Outputs dataframe with either empirical p-values and z-values, p-values
#' (non-empirical) or enrichment score (es) values for the various cell types.

#' @importFrom magrittr %>%


#' @export
get_alias <- function(input,
                      input_type = c("ensembl", "uniprot", "hgnc"),
                      ensembl_dataset = "hsapiens_gene_ensembl",
                      ensembl_verbose = TRUE,
                      return_only_ensembl = FALSE,
                      download_biomart_data = TRUE
){

  if(download_biomart_data){
    loop <- 5

    while(loop != 0){
      ### Regarding biomaRt ###
      try(
        ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = ensembl_dataset, verbose = ensembl_verbose)
      )

      if(exists("ensembl")){
        loop <- 0
      } else {
        loop <- loop - 1
      }
    }
  }

  if(!exists("ensembl")){
    ensembl <- readRDS(system.file("biomart_data", "biomart_2020.rds", package = "cellex.gsea"))
  }

  if(input_type[1] == "ensembl"){
    alias_df <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "uniprotswissprot"),
                               filters = "ensembl_gene_id", values = input, mart = ensembl)
    alias_df <- alias_df %>%
      dplyr::arrange(desc(hgnc_symbol), desc(uniprotswissprot)) %>%
      dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
    missing <- input[!(input %in% alias_df$ensembl_gene_id)]
    alias_df <- alias_df %>% dplyr::add_row(data.frame(ensembl_gene_id = missing,
                                                       hgnc_symbol = rep("", length(missing)),
                                                       uniprotswissprot = rep("", length(missing))))
  } else if(input_type[1] == "hgnc"){
    alias_df <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "uniprotswissprot"),
                               filters = "hgnc_symbol", values = input, mart = ensembl)
    alias_df <- alias_df %>%
      dplyr::arrange(hgnc_symbol, desc(uniprotswissprot)) %>%
      dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
    missing <- input[!(input %in% alias_df$hgnc_symbol)]
    alias_df <- alias_df %>% dplyr::add_row(data.frame(ensembl_gene_id = rep("", length(missing)),
                                            hgnc_symbol = missing,
                                            uniprotswissprot = rep("", length(missing))))
  } else if(input_type[1] == "uniprot"){
    alias_df <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "uniprotswissprot"),
                               filters = "uniprotswissprot", values = input, mart = ensembl)
    alias_df <- alias_df %>%
      dplyr::arrange(desc(hgnc_symbol), desc(uniprotswissprot)) %>%
      dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
    missing <- input[!(input %in% alias_df$uniprotswissprot)]
    alias_df <- alias_df %>% dplyr::add_row(data.frame(ensembl_gene_id = rep("", length(missing)),
                                            hgnc_symbol = rep("", length(missing)),
                                            uniprotswissprot = missing))
  }

  alias_df <- alias_df %>% dplyr::filter(ensembl_gene_id != "")

  if(return_only_ensembl){
    return(alias_df$ensembl_gene_id)
  } else {
    return(alias_df)
  }
}



first_order_network <- function(input,
                                ppi_network = "InBio_Map_core_2016_09_12/core_v3.psimitab",
                                col_names = c("unique_A", "unique_B"),
                                ensembl_dataset = "hsapiens_gene_ensembl",
                                ensembl_verbose = TRUE,
                                return_only_ensembl = FALSE,
                                download_biomart_data = FALSE
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
  if(download_biomart_data){
    loop <- 20

    while(loop != 0){
      ### Regarding biomaRt ###
      try(
        ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = ensembl_dataset, verbose = ensembl_verbose)
      )

      if(exists("ensembl")){
        loop <- 0
      } else {
        loop <- loop - 1
      }
    }
  } else {
    ensembl <- readRDS(system.file("biomart_data", "biomart_2020.rds", package = "cellex.gsea"))
  }

  if(!exists("ensembl")){
    ensembl <- readRDS(system.file("biomart_data", "biomart_2020.rds", package = "cellex.gsea"))
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
  rm(ensembl)

  if(return_only_ensembl){
    return(alias_df$ensembl_gene_id)
  } else {
    return(alias_df)
  }
}



fastRndWalk <- function(gSetIdx, geneRanking, j, Ra) {
  # Function used for the ssgsea approach of calculating ES values.
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



ssgsea <- function(data, geneSet, alpha=0.25, n_cores = 1, es_val_norm = FALSE){
  # Function for calculating the ES values, as used in the ssgsea approach of the GSVA package.
  p <- nrow(data)
  n <- ncol(data)

  R <- apply(data, 2, function(x, p) as.integer(rank(x, ties.method = "min")), p) #"average", "min", "max"
  Ra <- abs(R)^alpha

  gSetIdx <- which(rownames(data) %in% geneSet)

  es <- matrix(NA, nrow=1, ncol=ncol(data))

  es <- parallel::mclapply(as.list(1:n), mc.cores = n_cores, FUN = function(j) {
    geneRanking <- sort.list(R[,j], decreasing=TRUE)
    es_sample <- lapply(list(gSetIdx), fastRndWalk, geneRanking, j, Ra)

    unlist(es_sample)
  })

  es <- do.call("cbind", es)
  es <- matrix(es, nrow=1)

  # Normalization of values (to values in the range 0 to 1).
  if(es_val_norm){
    es <- apply(es, 2, function(x, es) (x - min(es)) / (max(es) - min(es)), es)
    es <- matrix(es, nrow=1)
  }

  # Output.
  es <- data.frame(tissue.cell = colnames(data), es.value = es[1,])

  return(es)
}



statistical_test <- function(data,
                             subset,
                             background,
                             stat_test = c("KS", "T", "W", "ES"),
                             emp_p_val = TRUE,
                             p_val_adjust = FALSE,
                             n_reps = 1000,
                             n_cores = 1,
                             n_background = 0,
                             get_null_dist = TRUE){
  # Function to perform statistical tests and calculate p-values.

  # Set seed.
  set.seed(42)

  # Load data and gene set, and perform analysis.
  data <- tibble::column_to_rownames(.data = data, var = "gene")
  names <- colnames(data)

  if(emp_p_val & n_background == 0){

    if(get_null_dist & length(subset) <= 1005) {
      if(length(subset) <= 1000) {
        numbers <- c(1:99, seq(from = 100, to = 1000, by = 5))
        if(length(subset) %in% numbers){
          number <- length(subset)
        } else {
          if(length(subset) %% 5 == 1) {
            number <- length(subset) - 1
          } else if(length(subset) %% 5 == 2) {
            number <- length(subset) - 2
          } else if(length(subset) %% 5 == 3) {
            number <- length(subset) + 2
          } else if(length(subset) %% 5 == 4) {
            number <- length(subset) + 1
          }
        }
        rm(numbers)

      } else {
        number <- 1000
      }

      dir.create("./tmp", showWarnings = FALSE)
      token <- readRDS(system.file("token", "token.rds", package = "cellex.gsea"))

      if(stat_test[1] == "KS") {
        rdrop2::drop_download(path = paste0("/null_dists/KS_less/", number, ".rds"),
                              local_path = paste0("./tmp/", number, ".rds"),
                              overwrite = TRUE,
                              dtoken = token)
      } else {
        rdrop2::drop_download(path = paste0("/null_dists/", stat_test[1], "/", number, ".rds"),
                              local_path = paste0("./tmp/", number, ".rds"),
                              overwrite = TRUE,
                              dtoken = token)
      }

      stat_df <- readRDS(paste0("./tmp/", number, ".rds"))
      file.remove(paste0("./tmp/", number, ".rds"))
      unlink("./tmp", recursive = TRUE)
      rm(number)
      rm(token)

    } else {
      to_sample <- replicate(n_reps, {
        sample(x = rownames(data), size = length(subset))
        })
      if(length(to_sample) == n_reps){
        to_sample <- matrix(to_sample, ncol = n_reps)
        }

      if(stat_test[1] == "W"){
        # New approach to calculate null test statistics.
        stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
          stat_vec <- matrixTests::col_wilcoxon_twosample(x = data[to_sample[,i],],
                                                          y = data[!(rownames(data) %in% to_sample[,i]),],
                                                          alternative = "greater")$statistic
          return(stat_vec)
        })

      } else if(stat_test[1] == "T"){
        # New approach to calculate null test statistics.
        stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
          stat_vec <- matrixTests::col_t_welch(x = data[to_sample[,i],],
                                               y = data[!(rownames(data) %in% to_sample[,i]),],
                                               alternative = "greater")$statistic
          return(stat_vec)
        })

      } else if(stat_test[1] == "KS"){
        stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
          stat_vec <- numeric(length(names))
          for(j in 1:length(stat_vec)){
            stat_vec[j] <- ks.test(x = data[(to_sample[,i]),j],
                                   y = data[!(rownames(data) %in% to_sample[,i]),j],
                                   alternative = "less")$statistic
          }
          return(stat_vec)
        })

      } else if(stat_test[1] == "ES") {
        stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
          es_vec <- ssgsea(data, geneSet = to_sample[,i], n_cores = n_cores, es_val_norm = FALSE)$es.value
          return(es_vec)
        })

      }

      if(n_reps > 1) {
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

    to_sample <- replicate(n_reps, {
      sample(x = rownames(data), size = length(subset))
    })
    if(length(to_sample) == n_reps){
      to_sample <- matrix(to_sample, ncol = n_reps)
    }

    to_sample_back <- apply(X = to_sample, MARGIN = 2, FUN = function(x){
      gene_names <- rownames(data)[!(rownames(data) %in% x)]
      gene_names <- sample(x = gene_names, size = n_background)
      return(gene_names)
    })

    if(stat_test[1] == "W"){
      # New approach to calculate null test statistics.
      stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
        stat_vec <- matrixTests::col_wilcoxon_twosample(x = data[to_sample[,i],],
                                                        y = data[to_sample_back[,i],])$statistic
        return(stat_vec)
      })

    } else if(stat_test[1] == "T"){
      # New approach to calculate null test statistics.
      stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
        stat_vec <- matrixTests::col_t_welch(x = data[to_sample[,i],],
                                             y = data[to_sample_back[,i],])$statistic
        return(stat_vec)
      })

    } else if(stat_test[1] == "KS"){
      stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
        stat_vec <- numeric(length(names))
        for(j in 1:length(stat_vec)){
          stat_vec[j] <- ks.test(x = data[to_sample[,i],j],
                                 y = data[to_sample_back[,i],j])$statistic
        }
        return(stat_vec)
      })

    } else if(stat_test[1] == "ES") {
      stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
        es_vec <- ssgsea(data[to_sample_back[,i],], geneSet = to_sample[,i], n_cores = 1, es_val_norm = FALSE)$es.value
        return(es_vec)
      })

    }

    if(n_reps > 1) {
      stat_df <- Reduce(rbind, stat_df_list)
      colnames(stat_df) <- names
      stat_df <- tibble::tibble(as.data.frame(stat_df))
    } else {
      stat_df <- t(Reduce(rbind, stat_df_list))
      colnames(stat_df) <- names
      stat_df <- tibble::tibble(as.data.frame(stat_df))
    }

  }



  # Generate analytical test statistics.
  if(stat_test[1] == "W"){
    output <- matrixTests::col_wilcoxon_twosample(x = data[subset,],
                                                  y = data[background,],
                                                  alternative = "greater")
    output <- output %>%
      tibble::rownames_to_column(var = "tissue.cell") %>%
      tibble::as_tibble() %>%
      dplyr::rename(p.value = pvalue)
  } else if(stat_test[1] == "T") {
    output <- matrixTests::col_t_welch(x = data[subset,],
                                       y = data[background,],
                                       alternative = "greater")
    output <- output %>%
      tibble::rownames_to_column(var = "tissue.cell") %>%
      tibble::as_tibble() %>%
      dplyr::rename(p.value = pvalue)
  } else if(stat_test[1] == "KS") {
    output <- parallel::mclapply(names, mc.cores = n_cores, FUN = function(col){
      test <- ks.test(x = data[subset, col],
                      y = data[background, col],
                      alternative = "less")
      df.test <- broom::tidy(test)
      df.test <- df.test %>% tibble::add_column(tissue.cell = col, .before = "statistic")
      return(df.test)
    })
    output <- tibble::tibble(Reduce(rbind, output))
  } else if(stat_test[1] == "ES") {
    if(emp_p_val) {
      output <- tibble::as_tibble(ssgsea(data, geneSet = subset, n_cores = n_cores, es_val_norm = FALSE))
    } else {
      output <- tibble::as_tibble(ssgsea(data, geneSet = subset, n_cores = n_cores, es_val_norm = TRUE))
    }
  }

  # Which values to extract.
  if(emp_p_val & stat_test[1] != "ES"){
    output <- output %>% dplyr::select(tissue.cell, statistic)
  } else if(stat_test[1] != "ES") {
    output <- output %>% dplyr::select(tissue.cell, p.value)
  } else if(stat_test[1] == "ES") {
    output <- output %>% dplyr::select(tissue.cell, es.value)
  }

  # Calculate empirical p-values and z-values from the previously computed statistical values.
  if(emp_p_val){

    p_z_val_out <- parallel::mclapply(colnames(stat_df), mc.cores = n_cores, FUN = function(col){

      if(stat_test[1] == "ES") {
        stat_an <- output %>% dplyr::filter(tissue.cell == col) %>% dplyr::select(es.value)
      } else {
        stat_an <- output %>% dplyr::filter(tissue.cell == col) %>% dplyr::select(statistic)
      }

      stat_an <- stat_an[[1]]
      stat_null <- stat_df[,col]
      stat_null <- stat_null[[1]]

      # Calculate empirival p-values, and z-values.
      p_value <- ( sum(stat_null > stat_an) +1 ) / (n_reps+1) # one-sided
      z_value <- (stat_an - mean(stat_null)) / sd(stat_null)

      return(c(p_value, unname(z_value)))
    })

    p_z_val_out <- tibble::tibble(as.data.frame(Reduce(rbind, p_z_val_out)))
    names(p_z_val_out) <- c("p.value", "Z.value")
    p_val_out <- p_z_val_out$p.value
    z_val_out <- p_z_val_out$Z.value

    # p_val_out[p_val_out > 1] <- 1.0
    # p_val_out[p_val_out == 0] <- 1 / n_reps
    names(p_val_out) <- output$tissue.cell
    names(z_val_out) <- output$tissue.cell

  }


  # Create output data frame and sort by adjusted P-values (multiple testing).
  if(emp_p_val){
    output <- output %>% tibble::add_column(p.value = p_val_out)
    output <- output %>% tibble::add_column(z.value = z_val_out, .after = "p.value")
  }

  # Whether to adjust the calculated p-values.
  if(p_val_adjust){
    output$p.value <- p.adjust(output$p.value, method = "fdr")
  }

  # Sort the output values.
  if(emp_p_val) {
    output <- output %>% dplyr::arrange(p.value, desc(z.value))
  } else if(stat_test == "ES" & !emp_p_val) {
    output <- output %>% dplyr::arrange(desc(es.value))
  } else {
    output <- output %>% dplyr::arrange(p.value)
  }

  if(emp_p_val) {
    output_list <- list(output = output, stat_df = stat_df)
  } else {
    output_list <- list(output = output, stat_df = NULL)
  }

  return(output_list)
}



#' @export
plot_bar <- function(output,
                     value_to_plot,
                     n_tissues = 10,
                     font_size = 8,
                     save_plots = TRUE){
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
  }
  if(save_plots) {
    dir.create("./plots", showWarnings = FALSE)
    ggplot2::ggsave(paste0("plots/output_barplot_", value_to_plot,".png"))
  }
}



#' @export
plot_box <- function(output,
                     subset,
                     cellex_data,
                     n_tissues = 10,
                     param = "ESmu",
                     font_size = 8,
                     save_plots = TRUE){ # background
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
  if(save_plots) {
    dir.create("./plots", showWarnings = FALSE)
    ggplot2::ggsave("plots/output_boxplot.png")
  }
}

#' @export
plot_hist <- function(output,
                      stat_df,
                      statistic = c("KS", "T", "W", "ES"),
                      save_plots = TRUE) { # "output" parameter needs to have a column with "statistic" values!
  tissues <- output$tissue.cell[1:12]
  if(statistic[1] != "ES") {
    intercept <- output %>% dplyr::filter(tissue.cell %in% tissues) %>% dplyr::select(statistic)
  } else {
    intercept <- output %>% dplyr::filter(tissue.cell %in% tissues) %>% dplyr::select(es.value)
  }
  intercept <- unname(intercept[[1]])
  intercepts <- c()
  for(number in intercept) {
    intercepts <- c(intercepts, rep(x = number, times = nrow(stat_df)))
  }
  stat_df <- stat_df %>% dplyr::select(tissues)
  stat_df <- reshape2::melt(stat_df, )
  colnames(stat_df) <- c("tissue.cell", "statistic")
  stat_df <- stat_df %>% tibble::add_column(intercept = intercepts)


  ggplot2::ggplot(data = stat_df, ggplot2::aes(x = statistic)) +
    ggplot2::geom_histogram() +
    ggplot2::xlab(paste0(statistic[1], " statistic")) +
    ggplot2::ggtitle(paste0("Histograms of statistic values")) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, hjust = 0.5)) +
    ggplot2::facet_wrap(~tissue.cell, scales = "free") + # maybe use 'scales = "fixed" ' instead?
    ggplot2::geom_vline(data = dplyr::filter(stat_df, tissue.cell == tissues[1]), ggplot2::aes(xintercept = unique(intercept)), color = "red", linetype = "longdash", size = 1) +
    ggplot2::geom_vline(data = dplyr::filter(stat_df, tissue.cell == tissues[2]), ggplot2::aes(xintercept = unique(intercept)), color = "red", linetype = "longdash", size = 1) +
    ggplot2::geom_vline(data = dplyr::filter(stat_df, tissue.cell == tissues[3]), ggplot2::aes(xintercept = unique(intercept)), color = "red", linetype = "longdash", size = 1) +
    ggplot2::geom_vline(data = dplyr::filter(stat_df, tissue.cell == tissues[4]), ggplot2::aes(xintercept = unique(intercept)), color = "red", linetype = "longdash", size = 1) +
    ggplot2::geom_vline(data = dplyr::filter(stat_df, tissue.cell == tissues[5]), ggplot2::aes(xintercept = unique(intercept)), color = "red", linetype = "longdash", size = 1) +
    ggplot2::geom_vline(data = dplyr::filter(stat_df, tissue.cell == tissues[6]), ggplot2::aes(xintercept = unique(intercept)), color = "red", linetype = "longdash", size = 1) +
    ggplot2::geom_vline(data = dplyr::filter(stat_df, tissue.cell == tissues[7]), ggplot2::aes(xintercept = unique(intercept)), color = "red", linetype = "longdash", size = 1) +
    ggplot2::geom_vline(data = dplyr::filter(stat_df, tissue.cell == tissues[8]), ggplot2::aes(xintercept = unique(intercept)), color = "red", linetype = "longdash", size = 1) +
    ggplot2::geom_vline(data = dplyr::filter(stat_df, tissue.cell == tissues[9]), ggplot2::aes(xintercept = unique(intercept)), color = "red", linetype = "longdash", size = 1) +
    ggplot2::geom_vline(data = dplyr::filter(stat_df, tissue.cell == tissues[10]), ggplot2::aes(xintercept = unique(intercept)), color = "red", linetype = "longdash", size = 1)
    # ggplot2::geom_vline(data = dplyr::filter(stat_df, tissue.cell == tissues[11]), ggplot2::aes(xintercept = unique(intercept)), color = "red", linetype = "longdash", size = 1) +
    # ggplot2::geom_vline(data = dplyr::filter(stat_df, tissue.cell == tissues[12]), ggplot2::aes(xintercept = unique(intercept)), color = "red", linetype = "longdash", size = 1)

  if(save_plots) {
    dir.create("./plots", showWarnings = FALSE)
    ggplot2::ggsave("plots/stat_histogram.png", units = "cm", width = 40, height = 22.5)
  }
}



#' export
filter_results <- function(input, p_threshold = c(0.001, 0.005, 0.01, 0.05)){

  output <- input %>% dplyr::filter(p.value <= p_threshold)

  if("z.value" %in% colnames(input)) {
    output <- output %>% dplyr::arrange(desc(z.value))
  } else {
    output <- output %>% dplyr::arrange(desc(p.value))
  }

  return(output)
}



#' Perform Enrichment Analysis
#'
#' This function takes as input a gene set (in ENSEMBL, UNIPROT, or HGNC format) and performs
#' enrichment analysis / gene specificity analysis and outputs a dataframe with p-values and/or
#' scores indicating how specifically expressed the gene set is in different cell types. The
#' p-values are calculated by comparing the CELLEX values for the input gene set to the CELLEX values
#' for the background gene set.
#'
#' @param input_set Gene set or protein set, as a character vector.
#' @param input_type What format the input_set is in. One of "ensembl" (ENSEMBL format), "uniprot" (UNIPROT format), or "gene" (HGNC format). Defaults to "ensembl" if no user specified value is selected.
#' @param cellex_data The CELLEX dataset to use in the package - must be in .csv format, and can also be gzip compressed as well. The path to a user-provided CELLEX dataset can be put here. Otherwise, three CELLEX datasets are included in the package, which can also be used. Set "cellex_data = 1" (default) for tabula_muris cellex dataset, set "cellex_data = 2" for gtex_v8 cellex dataset, or set "cellex_data = 3 for human cell landscape (HCL) cellex dataset.
#' @param statistic Which test statistic to use. One of "KS" (Kolmogorov-Smirnov test), "W" (Wilcoxon test), "T" (Student's t-test), and "ES" (the ssGSEA approach). Defaults to "KS".
#' @param emp_p_value If TRUE, calculates empirical p-values by randomly sampling a number of genes from the CELLEX dataset equivalent to the number of genes in the input_set, which is repeated a number of times (see parameter "reps") to obtain a null distribution. The null distribution is then compared to the test statistic for the genes in the input_set. Related parameters: "reps", "num_cores", "num_background_genes", "statistic_plot".
#' @param p_threshold The p-value threshold used for filtering the results in the outputted table. Is only used if the "p_threshold" parameter is set to TRUE.
#' @param p_value_adjust If TRUE, adjusts the calculated p-values for multiple testing, using the FDR approach. Defaults to FALSE.
#' @param n_cores Default is 1. The number of cores/threads to use when computing the null distribution for the empirical p-value (see parameter "emp_p_value").
#' @param n_reps Default is 1000. The number of random samplings to perform when computing the null distribution that is used to calculate the empirical p-value (see parameter "emp_p_value").
#' @param n_background Default is 0 (which means to use all background genes). The number of background genes to use when computing the null distribution for the empirical p-value (see parameter "emp_p_value").
#' @param all_genes_as_background If TRUE, use all genes in the CELLEX dataset as background for analysis. If FALSE, use all genes in the CELLEX dataset which are not in the input_set, as background for analysis.
#' @param first_order If TRUE, obtain the first order network for the input genes / proteins, and use for the analysis. Defaults to FALSE.
#' @param generate_plots If TRUE, generates plots when running the tool. Defaults to TRUE.
#' @param save_plots If TRUE, saves the plots the tool created to the working directory, to a folder named "plots". Defaults to TRUE.
#' @param save_output If TRUE, saves the generated table to the working directory, to a folder named "outputs". The table is saved as both a .csv and .rds file. Defaults to TRUE.
#' @param download_biomart If TRUE, uses the biomaRt package to connect to an online database of gene names, for gene names conversion - NB: Can increase computation time. If FALSE, uses an offline "snapshot" of the database, which comes with the tool. Defaults to FALSE.
#' @param get_null_dist Used in combination with the "emp_p_value" parameter. If "get_null_dist" is TRUE, downloads the appropriate null distribution of test statistics for the gene set size used when running the tool. If FALSE, a new null distribution is calculated locally. Defaults to TRUE.
#' @param del_stat_vals Parameter reserved for use in online tool. DO NOT CHANGE.
#' @param output_stat_df Parameter reserved for use in online tool. DO NOT CHANGE.


#' @export
gsea_analysis <- function(input_set, # input gene set or protein set.
                          input_type = c("ensembl", "uniprot", "hgnc"),
                          cellex_data = c(1, 2, 3),
                          statistic = c("KS", "T", "W", "ES"),
                          emp_p_val = TRUE,
                          filter_output = FALSE,
                          p_threshold = 0.001,
                          p_val_adjust = FALSE,
                          n_cores = 1,
                          n_reps = 1000,
                          n_background = 0,
                          all_genes_as_background = FALSE,
                          first_order = FALSE,
                          generate_plots = TRUE,
                          save_plots = TRUE,
                          save_output = TRUE,
                          download_biomart = FALSE,
                          get_null_dist = TRUE,
                          del_stat_vals = TRUE,
                          output_stat_df = FALSE
                          ) {

  set.seed(42)

  # Obtain aliases for the input gene set.
  message("Using BioMart to obtain ENSEMBL, UNIPROT and HGNC names...")
  gene_set <- get_alias(input = input_set,
                        input_type = input_type,
                        download_biomart_data = download_biomart)
  message("Done!")

  # Get PPI data matching the protein names from the alias dataframe.
  # Also get expression values for each ES?.

  if(first_order){
    #source("ppi.R")
    message("\nObtaining first order network...")
    proteins <- unique(gene_set$uniprotswissprot)
    proteins <- proteins[!(proteins == "")]
    gene_set <- first_order_network(input = proteins,
                                    download_biomart_data = download_biomart)
    message("Done!")
  }

  # Load cellex data.
  message("\nLoading CELLEX data...")
  if(cellex_data[1] == 1) {
    cellex_data <- readr::read_csv(system.file("cellex_data", "tabula_muris.gz", package = "cellex.gsea"))
  } else if(cellex_data[1] == 2) {
    cellex_data <- readr::read_csv(system.file("cellex_data", "gtex_v8.gz", package = "cellex.gsea"))
  } else if(cellex_data[1] == 3) {
    cellex_data <- readr::read_csv(system.file("cellex_data", "hcl.gz", package = "cellex.gsea"))
  } else {
    cellex_data <- readr::read_csv(cellex_data)
  }

  # Choose subset genes and background genes.
  subset_genes <- cellex_data$gene[cellex_data$gene %in% gene_set$ensembl_gene_id]

  if(all_genes_as_background){
    background_genes <- cellex_data$gene
  } else {
    background_genes <- cellex_data$gene[!(cellex_data$gene %in% subset_genes)]
  }

  gene_set <- gene_set %>% dplyr::filter(ensembl_gene_id %in% subset_genes)
  message("Done!")

  # Perform a statistical test for each tissue in the cellex data.
  # Possibility to perform permutation tests as well.
  # Load statistical tests functions.
  #source("statistical_tests.R")

  # Perform statistical tests.
  message("\nRunning analysis...")
  output_list <- statistical_test(data = cellex_data,
                                 subset = subset_genes,
                                 background = background_genes,
                                 stat_test = statistic[1],
                                 p_val_adjust = p_val_adjust,
                                 emp_p_val = emp_p_val,
                                 n_reps = n_reps,
                                 n_cores = n_cores,
                                 n_background = n_background,
                                 get_null_dist = get_null_dist)

  output <- output_list$output
  stat_df <- output_list$stat_df
  rm(output_list)
  message("Done!")

  if(is.null(output)) {
    message("\nERROR!\nNo output was generated.\nExiting...")
  } else {

    if(generate_plots & emp_p_val) {
      # Histograms.
      message("\nGenerating histograms of statistics values...")
      plot_hist(output = output, stat_df = stat_df, statistic = statistic[1], save_plots = save_plots)
      message("Done!")
    }


    if(emp_p_val) {
      if(statistic[1] != "ES" & del_stat_vals) {
        # Remove the "statistic" column.
        output <- output %>% dplyr::select(-statistic)
      } else if(statistic[1] == "ES" & del_stat_vals) {
        # Remove the "es.value" column.
        output <- output %>% dplyr::select(-es.value)
      }
    }


    # Whether to filter output before returning it.
    if(filter_output) {
      output <- filter_results(input = output, p_threshold = p_threshold)
    }

    # Whether to save the output to a file.
    if(save_output) {
      message("\nSaving output...")
      dir.create("./outputs", showWarnings = FALSE)
      saveRDS(output, file = "outputs/output.rds")
      readr::write_csv(output, file = "outputs/output.csv")
    }
    message("Done!")

    # Whether to save plots.
    if(generate_plots) {
      message("\nGenerating plots...")

      for(value in c("es.value", "z.value", "p.value")) {
        if(value %in% colnames(output)){
          plot_bar(output, value, save_plots = save_plots)
        }
      }
      plot_box(output, subset_genes, cellex_data, save_plots = save_plots)

      message("Done!\n")
    }



    if(output_stat_df) {
      output_list <- list(output = output, stat_df = stat_df)
      return(output_list)
    } else {
      #output_list <- list(output = output, stat_df = NULL)
      #return(output_list)
      return(output)
    }
  }
}
