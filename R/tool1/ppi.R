#' @title gene set specificity - functions to use for tool 1
#' @author Mads Porse Pedersen

# Function to load in the full PPI network data from Intomics, and given an input (vector) of
# proteins it extracts all proteins all these proteins are connected to. The function also connects
# to the ENSEMBL database and extracts gene IDs for the proteins in the network.
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


# get_ppi_data_subset <- function(ppi_data,
#                                 data_subset,
#                                 swissprot. = swissprot,
#                                 val = c("diff_max_min", "median", "mean", "max", "min")
# ){
#   CELLEX_genes <- data_subset$gene
#   swissprot_subset <- swissprot %>% dplyr::filter(ensembl_gene_id %in% CELLEX_genes)
#   #swissprot_subset <- swissprot_subset[!duplicated(swissprot_subset$ensembl_gene_id),]
#   expr_vec <- c()
#   for(gene_name in swissprot_subset$ensembl_gene_id){
#     expr <- unlist(data_subset %>% dplyr::filter(gene == gene_name) %>% dplyr::select(-gene))
#     if(val[1] == "diff_max_min") {
#       expr <- max(expr[expr != 0]) - min(expr[expr != 0])
#     } else if(val[1] == "median") {
#       expr <- median(expr[expr != 0])
#     } else if(val[1] == "mean") {
#       expr <- mean(expr[expr != 0])
#     } else if(val[1] == "max") {
#       expr <- max(expr[expr != 0])
#     } else if(val[1] == "min") {
#       expr <- min(expr[expr != 0])
#     }
#     expr_vec <- c(expr_vec, expr)
#   }
#
#   swissprot_subset <- dplyr::bind_cols(swissprot_subset, expr_val = expr_vec) # use either "diff_max_min", "mean", "max or "min"
#   gene_expr_df <- swissprot_subset %>% dplyr::select(uniprotswissprot, expr_val)
#   ppi_data_subset <- ppi_data %>% dplyr::filter(unique_B %in% swissprot_subset$uniprotswissprot)
#   ppi_data_subset <- ppi_data_subset %>% dplyr::filter(unique_A %in% swissprot_subset$uniprotswissprot)
#   ppi_data_subset <- ppi_data_subset[,1:2]
#
#   ppi_expr_list <- list("ppi_data_subset" = ppi_data_subset, "gene_expr_df" = gene_expr_df)
#
#   return(ppi_expr_list)
# }

# Function to create PPI graph object using the "igraph" package.

# create_ppi_graph <- function(ppi_data,
#                              data_subset,
#                              swissprot. = swissprot){
#
#   proteins_total <- unique(unlist(ppi_data_subset))
#   CELLEX_genes <- swissprot %>% dplyr::filter(uniprotswissprot %in% proteins_total)
#   CELLEX_genes <- CELLEX_genes$ensembl_gene_id
#   ppi_data_graph <- igraph::graph_from_data_frame(ppi_data_subset)
#   return(ppi_data_graph)
# }
