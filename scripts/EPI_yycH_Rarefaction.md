Epidome yycH rarefaction
================
Anna Ingham, SSI Copenhagen
February 2020

``` r
library(phyloseq)
require(ggplot2)
```

    ## Loading required package: ggplot2

``` r
require(scales)
```

    ## Loading required package: scales

``` r
require(reshape2)
```

    ## Loading required package: reshape2

``` r
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-5

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ranacapa)
library(openxlsx)
library(tibble)
theme_set(theme_bw())
```

Load dada2 output from workspace
EPI\_yycH\_Bayes\_tax\_190920\_workspace\_mb35.RData and metadata
Epidome\_metadata.xlsx

#### Generate phyloseq object

``` r
sam_df <- sam_df %>% remove_rownames() %>% tibble::column_to_rownames('Sample')

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(tax), sample_data(sam_df))
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 122 taxa and 52 samples ]
    ## sample_data() Sample Data:       [ 52 samples by 2 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 122 taxa by 1 taxonomic ranks ]

``` r
ps <- prune_taxa(taxa_sums(ps) != 0, ps)
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 122 taxa and 52 samples ]
    ## sample_data() Sample Data:       [ 52 samples by 2 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 122 taxa by 1 taxonomic ranks ]

#### Generate rarefaction curves

``` r
set.seed(123)
p1 <- ggrare(ps, step = 500, se = FALSE, label = "Sample" ) 
```

    ## rarefying sample Extraction_control_1
    ## rarefying sample Extraction_control_2
    ## rarefying sample P01_nose_1
    ## rarefying sample P01_nose_2
    ## rarefying sample P01_skin_1
    ## rarefying sample P01_skin_2
    ## rarefying sample P02_nose_1
    ## rarefying sample P02_nose_2
    ## rarefying sample P02_skin_1
    ## rarefying sample P02_skin_2
    ## rarefying sample P03_nose_1
    ## rarefying sample P03_nose_2
    ## rarefying sample P03_skin_1
    ## rarefying sample P03_skin_2
    ## rarefying sample P04_nose_1
    ## rarefying sample P04_nose_2
    ## rarefying sample P04_skin_1
    ## rarefying sample P04_skin_2
    ## rarefying sample P05_nose_1
    ## rarefying sample P05_nose_2
    ## rarefying sample P05_skin_1
    ## rarefying sample P05_skin_2
    ## rarefying sample P06_nose_1
    ## rarefying sample P06_nose_2
    ## rarefying sample P06_skin_1
    ## rarefying sample P06_skin_2
    ## rarefying sample P07_nose_1
    ## rarefying sample P07_nose_2
    ## rarefying sample P07_skin_1
    ## rarefying sample P07_skin_2
    ## rarefying sample P08_nose_1
    ## rarefying sample P08_nose_2
    ## rarefying sample P08_skin_1
    ## rarefying sample P08_skin_2
    ## rarefying sample P09_nose_1
    ## rarefying sample P09_nose_2
    ## rarefying sample P09_skin_1
    ## rarefying sample P09_skin_2
    ## rarefying sample P10_nose_1
    ## rarefying sample P10_nose_2
    ## rarefying sample P10_skin_1
    ## rarefying sample P10_skin_2
    ## rarefying sample P11_nose_1
    ## rarefying sample P11_nose_2
    ## rarefying sample P11_skin_1
    ## rarefying sample P11_skin_2
    ## rarefying sample even-mock3-1_S258_L001
    ## rarefying sample even-mock3-2_S282_L001
    ## rarefying sample even-mock3-3_S199_L001
    ## rarefying sample staggered-mock3-1_S270_L001
    ## rarefying sample staggered-mock3-2_S211_L001
    ## rarefying sample staggered-mock3-3_S223_L001

![](EPI_yycH_Rarefaction_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
p1 + facet_wrap(.~site, ncol = 1)  + theme(panel.background = element_blank(), axis.title.x = element_text(size =14, face = "bold"), axis.title.y = element_text(size =14, face = "bold"), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16), strip.text.x = element_text(angle = 0, face = "bold", size = 12), strip.background =element_rect(fill="white")) + xlab("Post-QC library size per sample") + ylab("# of observed ASVs") + ggtitle("yycH gene - rarefaction curves") + geom_vline(xintercept=5000, color = "red", size = 0.8) + geom_vline(xintercept=10000, color = "orange", size = 0.8) + geom_vline(xintercept=2500, color = "darkred", size = 0.8) + geom_text(aes(x=2050, label="2500", y=20), colour="darkred", angle=90) + geom_text(aes(x=4550, label="5000", y=20), colour="red", angle=90) + geom_text(aes(x=9550, label="10000", y=20), colour="orange", angle=90) 
```

![](EPI_yycH_Rarefaction_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
p1 + scale_x_continuous(limits = c(0,25000), breaks = seq(0, 25000, 1000)) + geom_vline(xintercept=5000, color = "red", size = 2) + geom_vline(xintercept=10000, color = "orange", size = 2) + geom_vline(xintercept=2500, color = "darkred", size = 2) + xlab("Post-QC library size per sample") + ylab("# of observed ASVs") + ggtitle("yycH gene - rarefaction curves") 
```

    ## Warning: Removed 23 rows containing missing values (geom_text).

    ## Warning: Removed 992 rows containing missing values (geom_path).

![](EPI_yycH_Rarefaction_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

#### Zoom

``` r
p1 + scale_x_continuous(limits = c(0,12000), breaks = seq(0, 12000, 1000)) + geom_vline(xintercept=5000, color = "red", size = 2) + geom_vline(xintercept=10000, color = "orange", size = 2) + geom_vline(xintercept=2500, color = "darkred", size = 2) + xlab("Post-QC library size per sample") + ylab("# of observed ASVs") + ggtitle("yycH gene - rarefaction curves") 
```

    ## Warning: Removed 38 rows containing missing values (geom_text).

    ## Warning: Removed 1779 rows containing missing values (geom_path).

![](EPI_yycH_Rarefaction_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
