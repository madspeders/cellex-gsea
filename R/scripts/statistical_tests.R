#' @title gene set specificity - functions to use for tool 1
#' @author Mads Porse Pedersen

# For loop to perform the Wilcoxon test or KS test for each tissue in the
# Tabula Muris ES? data. Possibility to perfom permutation tests as well.

# Function for ssGSEA (borrowed from the GSVA implementation).
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

# ssGSEA approach.
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
