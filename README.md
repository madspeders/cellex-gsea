# CELLEX-GSEA
R package to calculate cell specific expression of genes, using an input gene set and a CELLEX dataset.  

## Installation instructions
Use the function `install_github` from the package `devtools` to install this R package (for the code to function properly, please update any other outdated packages in the process when/if prompted to):
```R
devtools::install_github("https://github.com/Mporse/cellex-gsea")
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

3. Run the enrichment analysis. If the input is in ENSEMBL format, you only need to supply the input gene set, as such:
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
gsea_analysis(gene_set, input_type = "hgnc")
```

### Using other CELLEX datasets  
The used CELLEX dataset can be provided by the user (so setting the argument `cellex_data` equal to the path to the dataset), or one of three datasets included in the package. Set the argument `cellex_data` equal to **1** to use the Tabula Muris CELLEX dataset, set equal to **2** to use the GTEX_v8 CELLEX dataset, or set sequal to **3** to use the HCL CELLEX dataset.  

### Specifying the preferred output  
It can be specified whether or not to calulate empirical p-values:  
* Set argument `emp_p_value` equal to `TRUE` to calculate empirical p-values for each cell type (this can take some time), otherwise set equal to `FALSE` (default).  

**NB:** If calculating empirical P-values, I recommend using the Tabula Muris CELLEX dataset. Pre-computed null distributions have been made available online only for this dataset, and the package will download the proper null distribution (according to the size of the gene set). If calculating empirical P-values without using the pre-computed null distributions, make sure to increase the number of utilized cores through the `n_cores` argument (only possible on Unix systems), otherwise the computation can take a long time.  
  
### Specifying the test statistics  
By default, the p-values are calculated using the Wilcoxon test. This can be changed by changing the value of the argument `statistic`:
* `statistic = "KS"` to use the Kolmogorov–Smirnov test.  
* `statistic = "T"` to use the Welch's T-test.  
* `statistic = "W"` to use the Mann-Whitney U (Wilcoxon rank-sum) test.  
* `statistic = "ES"` to use the ssGSEA approach from the GSVA package to calculate enrichment scores (ES).

### Other arguments  
To see which other arguments are available and which values they accept, check the R documentation for the function:  
```R
?gsea_analysis
```

## Web application  
The tool has also been made available as a Shiny web application. It is being hosted on _shinyapps.io_ using the free tier, which means the app is limited to 25 active hours per month. Thus, the use of the web application is quite limited and I would recommend using the R package version instead whenever possible, since it yields the same results.  
* Shiny web application: https://porse.shinyapps.io/cellex-gsea/

## Credits  
The tool makes use of ESµ values from the CELLEX tool:  
* GitHub: https://github.com/perslab/CELLEX
* Journal article (_P. N. Timshel, J. J. Thompson, and T. H. Pers_): https://elifesciences.org/articles/55851

The tool also makes use of code from the GSVA package for performing the ssGSEA approach:  
* GitHub: https://github.com/rcastelo/GSVA
* Bioconductor: https://www.bioconductor.org/packages/release/bioc/html/GSVA.html
* Journal article (_S. Hänzelmann, R. Castelo, and J. Guinney_): https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7
