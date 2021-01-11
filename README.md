# cellex.analysis
R package to calculate cell specific expression of genes, using an input gene set and a CELLEX dataset.  

## Installation instructions
Use the function `install_github` from the package `devtools` to install this R package:
```R
devtools::install_github("Mporse/cellex-gsea")
```

## Short tutorial
There are several functions in the package, however you should only use one for the analysis (all the other functions are helper functions to make the tool work correctly):
```R
gsea_analysis()
```

Follow these to run an analysis of a input gene set using the function:  
1. Load the package using:  
```R
library(cellex.gsea)
```

2. Save your gene set (in ENSEMBL or HGNC format) or protein set (in UNIPROT format) as a variable. A gene set that is enriched in the liver (obtained from the Human Protein Atlas, [proteinatlas.org](https://www.proteinatlas.org/)) is already included in the package, and can be loaded with the following command:
```R
gene_set <- readRDS(system.file("gene_set_data", "liver.rds", package = "cellex.gsea"))
```

3. Run the tool #1 analysis. If the input is in ENSEMBL format, you only need to supply the input gene set, as such:
```R
gsea_analysis(gene_set)
```

That's it!  

## Important function arguments  

### Using other gene set input formats  
The input gene set is not always in ENSEMBL format. Therefore, the argument `input_type` in the function needs to be changed:  
* **If the input set is in UNIPROT format:**  
```R
gsea_analysis(gene_set, input_type = "uniprot")
```

* **If the input set is in HGNC format:**
```R
gsea_analysis(gene_set, input_type = "gene")
```

### Using other CELLEX datasets  
The used CELLEX dataset can be provided by the user (so setting the argument `cellex_data` equal to the path to the dataset), or one of two datasets included in the package. Set the argument `cellex_data` equal to **1** to use the Tabula Muris CELLEX dataset, or set equal to **2** to use the GTEX_v8 CELLEX dataset.  

### Specifying the preferred output  
Several output values can be calculated using the tool:  
* Set argument `p_value` equal to `TRUE` (default) to calculate p-values for each cell type, otherwise set equal to `FALSE`.  
* Set argument `p_value_adjust` equal to `TRUE` (default) to correct the p-values using the FDR approach, otherwise set equal to `FALSE`.  
* Set argument `emp_p_value` equal to `TRUE` to calculate empirical p-values for each cell type (this can take some time), otherwise set equal to `FALSE` (default).  
* Set argument `es_value` equal to `TRUE` to calculate enrichment scores (using the ssGSEA approach from the GSVA package), otherwise set equal to `FALSE` (default).  
  
### Specifying the test statistics  
By default, the p-values are calculated using the Wilcoxon test. This can be changed by changing the value of the argument `statistic`:
* `statistic = "W"` to use the Wilcoxon test.  
* `statistic = "KS"` to use the Kolmogorov–Smirnov test.  
* `statistic = "T"` to use the Student's t-test.  
* `statistic = "ES"` to use the ssGSEA algorithm from the GSVA package to calculate enrichment scores (ES).

### Other arguments  
To see which other arguments are available and which values they accept, check the R documentation for the function:  
```R
?gsea_analysis
```
