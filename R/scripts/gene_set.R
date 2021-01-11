#' @title gene set specificity - functions to use for tool 1
#' @author Mads Porse Pedersen

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
# load_gene_set <- function(input,
#                           test_genes = FALSE,
#                           test_proteins = FALSE,
#                           tissue_group = c("Enriched", "Elevated", "Only")){
#
#   if(test_genes) {
#     # Select specific gene (if selected above).
#     gene_set <- input
#     return(gene_set)
#   } else if(test_proteins) {
#     alias_df <- get_alias(input, input_type = "uniprot")
#     gene_set <- alias_df$ensembl_gene_id
#     return(gene_set)
#   } else {
#     # Choose tissue to test.
#     if(typeof(tissue_group) == "character" & length(tissue_group) > 1){
#       group <- tissue_group[1]
#       tissue <- input
#     } else if(typeof(tissue_group) == "character" & length(tissue_group) == 1) {
#       group <- tissue_group
#       tissue <- input
#     }
#     gene_set <- readr::read_tsv(gzfile(paste0("tissue_data/Protein_atlas/", group, "/", tissue, ".tsv.gz")))
#     gene_set <- gene_set$Ensembl
#     return(gene_set)
#   }
# }
