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
                      return_only_ensembl = FALSE,
                      download_biomart_data = TRUE
){

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
    ensembl <- readRDS(system.file("biomart_data", "biomart_2020.rds", package = "cellex.analysis"))
  }

  if(!exists("ensembl")){
    ensembl <- readRDS(system.file("biomart_data", "biomart_2020.rds", package = "cellex.analysis"))
  }

  if(input_type[1] == "ensembl"){
    alias_df <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "uniprotswissprot"),
                               filters = "ensembl_gene_id", values = input, mart = ensembl)
    alias_df <- alias_df %>%
      dplyr::arrange(desc(hgnc_symbol), desc(uniprotswissprot)) %>%
      dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
    # alias_df <- alias_df %>%
    #   dplyr::arrange(hgnc_symbol, ensembl_gene_id) %>%
    #   dplyr::distinct(hgnc_symbol, .keep_all = TRUE)
    missing <- input[!(input %in% alias_df$ensembl_gene_id)]
    alias_df <- alias_df %>% dplyr::add_row(data.frame(ensembl_gene_id = missing,
                                                       hgnc_symbol = rep("", length(missing)),
                                                       uniprotswissprot = rep("", length(missing))))
  } else if(input_type[1] == "gene"){
    alias_df <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "uniprotswissprot"),
                               filters = "hgnc_symbol", values = input, mart = ensembl)
    alias_df <- alias_df %>%
      dplyr::arrange(hgnc_symbol, desc(uniprotswissprot)) %>%
      dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
    # alias_df <- alias_df %>%
    #   dplyr::arrange(hgnc_symbol, ensembl_gene_id) %>%
    #   dplyr::distinct(hgnc_symbol, .keep_all = TRUE)
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
    # alias_df <- alias_df %>%
    #   dplyr::arrange(hgnc_symbol, ensembl_gene_id) %>%
    #   dplyr::distinct(hgnc_symbol, .keep_all = TRUE)
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



#' @export
first_order_network <- function(input,
                                ppi_network = "InBio_Map_core_2016_09_12/core_v3.psimitab",
                                col_names = c("unique_A", "unique_B"),
                                ensembl_dataset = "hsapiens_gene_ensembl",
                                ensembl_verbose = TRUE,
                                return_only_ensembl = FALSE,
                                download_biomart_data = TRUE
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
    ensembl <- readRDS(system.file("biomart_data", "biomart_2020.rds", package = "cellex.analysis"))
  }

  if(!exists("ensembl")){
    ensembl <- readRDS(system.file("biomart_data", "biomart_2020.rds", package = "cellex.analysis"))
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
  #es <- apply(es, 2, function(x, es) x / max(es), es) # approach 1.
  es <- apply(es, 2, function(x, es) (x - min(es)) / (max(es) - min(es)), es) # approach 2.
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
                             emp_p_val = FALSE,
                             p_val_adjust = FALSE,
                             es_val = FALSE,
                             plot = FALSE,
                             z_val_norm = FALSE,
                             n_reps = 1000,
                             n_cores = 1,
                             n_background = 0,
                             get_null_dist = NULL,
                             save_null_dist = FALSE){

  # Set seed.
  set.seed(42)
  data <- tibble::column_to_rownames(.data = data, var = "gene")
  names <- colnames(data)

  # Create folder to store null distributions in.
  dir.create(paste0(getwd(), "/null_dists"), showWarnings = FALSE)

  if(emp_p_val & n_background == 0){

    if(paste0("null_dist_", stat_test[1], "_", length(subset), ".rds") %in% dir(paste0(getwd(), "/null_dists"))){
      stat_df <- readRDS(paste0(getwd(), "/null_dists/", paste0("null_dist_", stat_test[1], "_", length(subset), ".rds")))
    } else if(is.null(get_null_dist)){
      to_sample <- replicate(n_reps, {
        sample(x = rownames(data), size = length(subset))
      })

      if(stat_test[1] == "W"){
        # Old approach to calculate null test statistics.
        # stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
        #   stat_vec <- numeric(length(names))
        #   for(j in 1:length(stat_vec)){
        #     stat_vec[j] <- wilcox.test(x = data[to_sample[,i],j],
        #                                y = data[!(rownames(data) %in% to_sample[,i]),j])$statistic
        #   }
        #   return(stat_vec)
        # })

        # New approach to calculate null test statistics.
        stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
          stat_vec <- matrixTests::col_wilcoxon_twosample(x = data[to_sample[,i],],
                                                          y = data[!(rownames(data) %in% to_sample[,i]),])$statistic
          return(stat_vec)
        })

      } else if(stat_test[1] == "T"){
        # Old approach to calculate null test statistics.
        # stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
        #   stat_vec <- numeric(length(names))
        #   for(j in 1:length(stat_vec)){
        #     stat_vec[j] <- t.test(x = data[(to_sample[,i]),j],
        #                           y = data[!(rownames(data) %in% to_sample[,i]),j])$statistic
        #   }
        #   return(stat_vec)
        # })

        # New approach to calculate null test statistics.
        stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
          stat_vec <- matrixTests::col_t_welch(x = data[to_sample[,i],],
                                               y = data[!(rownames(data) %in% to_sample[,i]),])$statistic
          return(stat_vec)
        })

      } else if(stat_test[1] == "KS"){
        stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
          stat_vec <- numeric(length(names))
          for(j in 1:length(stat_vec)){
            stat_vec[j] <- ks.test(x = data[(to_sample[,i]),j],
                                   y = data[!(rownames(data) %in% to_sample[,i]),j])$statistic
          }
          return(stat_vec)
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

      if(save_null_dist){
        saveRDS(stat_df, paste0(getwd(), "/null_dists/null_dist_", stat_test[1], "_", length(subset), ".rds"))
      }

    } else {
      stat_df <- get_null_dist
    }


  } else if(emp_p_val & n_background != 0) {
    # Use only a smaller selection of the background genes instead of all of them.
    # For instance, n_background = 1000.
    # Reduces accuracy of output the smaller the background gene set, though.

    to_sample <- replicate(n_reps, {
      sample(x = rownames(data), size = length(subset))
    })

    to_sample_back <- apply(X = to_sample, MARGIN = 2, FUN = function(x){
      gene_names <- rownames(data)[!(rownames(data) %in% x)]
      gene_names <- sample(x = gene_names, size = n_background)
      return(gene_names)
    })

    if(stat_test[1] == "W"){
      # Old approach to calculate null test statistics.
      # stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
      #   stat_vec <- numeric(length(names))
      #   for(j in 1:length(stat_vec)){
      #     stat_vec[j] <- wilcox.test(x = data[to_sample[,i],j],
      #                                y = data[to_sample_back[,i],j])$statistic
      #   }
      #   return(stat_vec)
      # })

      # New approach to calculate null test statistics.
      stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
        stat_vec <- matrixTests::col_wilcoxon_twosample(x = data[to_sample[,i],],
                                                        y = data[to_sample_back[,i],])$statistic
        return(stat_vec)
      })

    } else if(stat_test[1] == "T"){
      # Old approach to calculate null test statistics.
      # stat_df_list <- parallel::mclapply(X = 1:n_reps, mc.cores = n_cores, FUN = function(i){
      #   stat_vec <- numeric(length(names))
      #   for(j in 1:length(stat_vec)){
      #     stat_vec[j] <- t.test(x = data[to_sample[,i],j],
      #                           y = data[to_sample_back[,i],j])$statistic
      #   }
      #   return(stat_vec)
      # })

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
                                                  y = data[background,])
    output <- output %>%
      tibble::rownames_to_column(var = "tissue.cell") %>%
      tibble::as_tibble() %>%
      dplyr::rename(p.value = pvalue)
  } else if(stat_test[1] == "KS") {
    output <- matrixTests::col_t_welch(x = data[subset,],
                                       y = data[background,])
    output <- output %>%
      tibble::rownames_to_column(var = "tissue.cell") %>%
      tibble::as_tibble() %>%
      dplyr::rename(p.value = pvalue)
  } else if(stat_test[1] == "T") {
    output <- parallel::mclapply(names, mc.cores = n_cores, FUN = function(col){
      test <- ks.test(x = data[subset, col],
                      y = data[background, col])
      df.test <- broom::tidy(test)
      df.test <- df.test %>% tibble::add_column(tissue.cell = col, .before = "statistic")
      return(df.test)
    })
    output <- tibble::tibble(Reduce(rbind, output))
  }

  # Which values to extract.
  if(emp_p_val){
    output <- output %>% dplyr::select(tissue.cell, statistic)
  } else {
    output <- output %>% dplyr::select(tissue.cell, p.value)
  }

  # Calculate empirical p-values and z-values from the previously computed statistical values.
  if(emp_p_val){
    p_z_val_out <- parallel::mclapply(colnames(stat_df), mc.cores = n_cores, FUN = function(col){
      stat_an <- output %>% dplyr::filter(tissue.cell == col) %>% dplyr::select(statistic)
      stat_an <- stat_an[[1]]
      stat_null <- stat_df[,col]
      stat_null <- stat_null[[1]]
      #p_value <- (sum(stat_null > stat_an) + 1) / (n_reps) # one-sided
      p_value <- ( sum(stat_null < -abs(stat_an)) + sum(stat_null > abs(stat_an)) ) / (n_reps) # two-sided
      z_value <- (stat_an - mean(stat_null)) / sd(stat_null)
      return(c(p_value, unname(z_value)))
    })

    p_z_val_out <- tibble::tibble(as.data.frame(Reduce(rbind, p_z_val_out)))
    names(p_z_val_out) <- c("p.value", "Z.value")
    p_val_out <- p_z_val_out$p.value
    z_val_out <- p_z_val_out$Z.value
    if(z_val_norm){
      #z_val_out <- z_val_out / max(z_val_out) # approach 1
      z_val_out <- (z_val_out - min(z_val_out)) / (max(z_val_out) - min(z_val_out)) # approach 2
    }

    p_val_out[p_val_out > 1] <- 1.0
    p_val_out[p_val_out == 0] <- 1 / n_reps
    names(p_val_out) <- output$tissue.cell
    names(z_val_out) <- output$tissue.cell

  }


  # Create output data frame and sort by adjusted P-values (multiple testing).
  if(emp_p_val){
    output <- output %>% tibble::add_column(p.value = p_val_out)
    output <- output %>% tibble::add_column(z.value = z_val_out, .after = "p.value")
  }

  # Whether to adjust the calculated p-values.
  if(!is.null(p_val_adjust)){
    output$p.value <- p.adjust(output$p.value, method = "fdr")
  }

  # Calculate ES values (from ssGSEA method) and append to output.
  if(es_val) {
    es <- ssgsea(data, geneSet = subset, n_cores = n_cores)
    output <- output %>% tibble::add_column(es.value = es$es.value)
  }

  # Sort the output values.
  if(es_val & emp_p_val) {
    output <- output %>% dplyr::arrange(desc(es.value), desc(z.value))
  } else if(es_val & !emp_p_val) {
    output <- output %>% dplyr::arrange(desc(es.value), p.value)
  } else {
    output <- output %>% dplyr::arrange(p.value)
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


  if(emp_p_val) {
    # Remove the "statistic" column.
    output <- output %>% dplyr::select(-statistic)
    # Return the output dataframe.
    return(output)
  } else {
    # Return the output dataframe.
    return(output)
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
#' This function takes as input a gene set (in ENSEMBL, UNIPROT, or HGNC format) and performs
#' enrichment analysis / gene specificity analysis and outputs a dataframe with p-values and/or
#' scores indicating how specifically expressed the gene set is in different cell types. The
#' p-values are calculated by comparing the CELLEX values for the input gene set to the CELLEX values
#' for the background gene set.
#'
#' @param input_set Gene set or protein set, as a character vector.
#' @param input_type What format the input_set is in. One of "ensembl" (ENSEMBL format), "uniprot" (UNIPROT format), or "gene" (HGNC format). Defaults to "ensembl" if no user specified value is selected.
#' @param cellex_data The CELLEX dataset to use in the package - must be in .csv format, and can also be gzip compressed as well. The path to a user-provided CELLEX dataset can be put here. Otherwise, three CELLEX datasets are included in the package, which can also be used. Set "cellex_data = 1" (default) for tabula_muris cellex dataset, set "cellex_data = 2" for gtex_v8 cellex dataset, or set "cellex_data = 3 for human cell landscape (HCL) cellex dataset.
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
                            cellex_data = c(1, 2, 3),
                            emp_p_val = FALSE,
                            p_val_adjust = FALSE,
                            es_val = TRUE,
                            z_val_norm = FALSE,
                            first_order = FALSE,
                            all_genes_as_background = FALSE,
                            statistic = c("W", "KS", "T"),
                            n_reps = 1000,
                            n_cores = 1,
                            n_background = 0,
                            plot = FALSE,
                            save_output = TRUE,
                            download_biomart = TRUE,
                            get_null_dist = NULL,
                            save_null_dist = FALSE
                            ) {

  set.seed(42)

  # Obtain aliases for the input gene set.
  #source("gene_set.R")
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
    cellex_data <- readr::read_csv(system.file("cellex_data", "tabula_muris.gz", package = "cellex.analysis"))
  } else if(cellex_data[1] == 2) {
    cellex_data <- readr::read_csv(system.file("cellex_data", "gtex_v8.gz", package = "cellex.analysis"))
  } else if(cellex_data[1] == 3) {
    cellex_data <- readr::read_csv(system.file("cellex_data", "hcl.gz", package = "cellex.analysis"))
  } else {
    cellex_data <- readr::read_csv(cellex_data)
  }

  # Criteria for which genes to use:
  # idx <- sapply(1:nrow(cellex_data), FUN = function(row_num){
  #   row <- sort(unlist(cellex_data[row_num,-1], use.names = FALSE), decreasing = TRUE)
  #   return( any(row > 0.5) & any((max(row) / mean(row)) > 10) & (row[1] / row[2]) > 1.1 & (sum(test != 0) / sum(test == 0)) < 1/3 )
  # })
  # good_genes <- unlist(cellex_data[idx,1], use.names = FALSE)

  # Choose subset genes and background genes.
  subset_genes <- cellex_data$gene[cellex_data$gene %in% gene_set$ensembl_gene_id]
  # subset_genes <- intersect(subset_genes, good_genes) # MAYBE
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
  output <- statistical_test(data = cellex_data,
                             subset = subset_genes,
                             background = background_genes,
                             stat_test = statistic[1],
                             p_val_adjust = p_val_adjust,
                             emp_p_val = emp_p_val,
                             es_val = es_val,
                             z_val_norm = z_val_norm,
                             plot = plot,
                             n_reps = n_reps,
                             n_cores = n_cores,
                             n_background = n_background,
                             get_null_dist = get_null_dist,
                             save_null_dist = save_null_dist)
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

    for(value in c("es.value", "z.value", "p.value")){
      if(value %in% colnames(output)){
        plot_bar(output, value)
      }
    }
    plot_box(output, subset_genes, cellex_data)

    message("Done!")

    return(output)
  }
}
