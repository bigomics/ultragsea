# Getting Started with ultragsea

## Introduction

ultragsea is a novel, ultrafast and memory optimized gene set enrichment
scoring algorithm not unlike fGSEA. ultragsea demonstrates highly
correlated gene set scores and p-values and typically is 100 times
faster than fGSEA.

## Example

### Preparing gensets

For this vignette, we retrieve the Hallmark gene sets using the
`msigdbr` package. Of course you can use your own gene sets or download
from somewhere else.

``` r

library(ultragsea)
library(msigdbr)
gs <- msigdbr::msigdbr(collection = "H")
#> Downloading gene sets (first use only, may take a few minutes)...
gmt <- tapply(gs$gene_symbol,gs$gs_name,list)
length(gmt)
#> [1] 50
```

### Performing the statistical test

Just for this example we use some random values for `fc` but in general
this is a fold change from some statistical comparison of your
experiment.

``` r

G <- gmt2mat(gmt)
fc <- rnorm(nrow(G))
names(fc) <- rownames(G)
res <- ultragsea(G, fc, format='as.gsea', method='cor')
res <- ultragsea(G, fc, format='as.gsea', method='ztest')
head(res)
#>                               pathway       pval      padj log2err          ES
#>                                <char>      <num>     <num>   <num>       <num>
#> 1:           HALLMARK_HEME_METABOLISM 0.96527162 0.9652716      NA  0.01423022
#> 2:     HALLMARK_INFLAMMATORY_RESPONSE 0.68056175 0.8725151      NA  0.13455730
#> 3:         HALLMARK_KRAS_SIGNALING_UP 0.63095925 0.8725151      NA  0.15700453
#> 4: HALLMARK_OXIDATIVE_PHOSPHORYLATION 0.67483574 0.8725151      NA -0.13711441
#> 5:              HALLMARK_ADIPOGENESIS 0.13486035 0.6686212      NA  0.48868756
#> 6:       HALLMARK_ALLOGRAFT_REJECTION 0.07052317 0.6686212      NA -0.59109455
#>           NES  size                                 leadingEdge
#>         <num> <int>                                      <list>
#> 1:  0.0384407   200    CA2,BNIP3L,CDC27,H1-0,BCAM,BLVRA,...[33]
#> 2:  0.3634855   200    IL4R,CCL2,HBEGF,INHBA,ITGB3,EREG,...[36]
#> 3:  0.4241232   200      CA2,HBEGF,INHBA,KLF4,EREG,AKT2,...[36]
#> 4: -0.3703931   200 ACAA1,SDHC,ECHS1,GPI,GRPEL1,ATP6V1F,...[27]
#> 5:  1.3201130   200  ALDOA,ECH1,GPX4,CRAT,DHCR7,GADD45A,...[41]
#> 6: -1.5967495   200    CCND2,TIMP1,LIF,IFNGR2,CCL7,LCP2,...[41]
```

For the fastest speed, especially when computing enrichment for multiple
foldchanges, it is better to call the low-level functions for geneset
correlation or geneset z-test directly. They do not return the nicely
formatted output tables as
[`ultragsea()`](https://bigomics.github.io/ultragsea/reference/ultragsea.md)
but they are much faster as they are optimized for matrix inputs `F`.

``` r

F <- cbind(fc,fc,fc,fc,fc)
res1 <- gset.cortest(G, F, compute.p=TRUE)
res2 <- gset.ztest(G, F)
```

### Plotting

You can proceed by plotting the gene set using e.g. the plotting
functions of the fgsea package.

``` r

#library(fgsea)
#k <- res$pathway[which.min(res$pval)]
#fgsea::plotEnrichment(gmt[[k]],fc,ticksSize = 0.2)
```

## Session info

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] msigdbr_26.1.1   ultragsea_0.1.12 BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.7-6        jsonlite_2.0.0      dplyr_1.2.1        
#>  [4] compiler_4.6.2      BiocManager_1.30.27 tidyselect_1.2.1   
#>  [7] slam_0.1-56         parallel_4.6.2      assertthat_0.2.1   
#> [10] jquerylib_0.1.4     systemfonts_1.3.3   textshaping_1.0.5  
#> [13] yaml_2.3.12         fastmap_1.2.0       lattice_0.23-1     
#> [16] R6_2.6.1            generics_0.1.4      qlcMatrix_0.9.9    
#> [19] curl_8.1.0          knitr_1.52          tibble_3.3.1       
#> [22] bookdown_0.48       desc_1.4.3          bslib_0.12.0       
#> [25] pillar_1.11.1       rlang_1.3.0         cachem_1.1.0       
#> [28] xfun_0.62           fs_2.1.0            sass_0.4.10        
#> [31] otel_0.3.0          cli_3.6.6           pkgdown_2.2.1      
#> [34] withr_3.0.3         magrittr_2.0.5      digest_0.6.39      
#> [37] grid_4.6.2          docopt_0.7.2        lifecycle_1.0.5    
#> [40] vctrs_0.7.3         sparsesvd_0.2-3     data.table_1.18.6.1
#> [43] evaluate_1.0.5      glue_1.8.1          ragg_1.5.2         
#> [46] rmarkdown_2.32      matrixStats_1.5.0   tools_4.6.2        
#> [49] pkgconfig_2.0.3     htmltools_0.5.9
```
