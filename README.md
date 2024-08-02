# wlab.block

Here I provide an R package to quickly implement the images in the previous article in our laboratory.

## Installation

You can easily install wlab.block like so:

``` r
# Open an R session and install the wlab.block R package
if(!require(devtools)) install.packages("devtools")
devtools::install_github("cheding/wlab.block")
```

## Example pipeline

This is a basic example pipeline which shows you how to process DiMSum report, and finally generate heatmaps:

``` r
## basic example code
library(wlab.block)
```

### Step1

Complete the amino acid sequence based on the wild-type sequence.

``` r
aa_complement(required_file,wt_aa)
```

### Step2

Normalize data between different blocks and select single mutations.

``` r
nor_fit<-nor_fitness(block1="path/to/block1",block2="path/to/block2")
nor_fit_pos<-pos_id(nor_fit,wt_aa)
nor_fit_single<-nor_fitness_single_mut(input=nor_fit_pos)
```

### Step3

Draw fitness heatmaps based on normalized data.

``` r
fitness_heatmap(input=nor_fit_single,wt_aa,title="fitness")
ggplot2::ggsave("fitness_heatmap.pdf", device = cairo_pdf,height = 4,width=20)
```
## Other pipeline

[Real backgrounds](docs/real_background.md)

