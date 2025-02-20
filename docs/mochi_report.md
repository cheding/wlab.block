# MoCHI Reports

MOCHI is a neural network model that we usually use to quantify energy changes in deep mutational scanning data. This pipeline will show you how to process the mochi results and generate some pictures.

## Example pipeline

This is a basic example pipeline which shows you how to use mochi results:
``` r
## basic example code
library(wlab.block)
setwd("path/to/mochi_result/weignts")
```

### generate binding/folding energy change heatmaps
Drawing a heatmap using Mochi is simpler, just give the file location and change the title.
``` r
ddG_heatmap(
    input="./weights_Folding.txt",
    wt_aa=wt_aa,
    title="folding free energy change"
)
ggplot2::ggsave("folding_heatmap.pdf", device = cairo_pdf,height = 4,width=20)
```
### observed fitness vs. predicted fitness
Mochi uses thermodynamic equations to predict fitness through binding and folding energy change, so comparing the observed fitness with the predicted fitness can demonstrate the credibility of the model.
Normalize the model prediction results, and the fitness in the output only corresponds to the fineness of "assay name".
``` r
# normalize prediction data
pre_nor<-nor_predict(
    pre_file="../predictions/predicted_phenotypes_all.txt",
    assay="assay_name",
    abundance1="path/to/abundance1.RData",
    abundance2="path/to/abundance2.RData",
    abundance3="path/to/abundance3.RData",
    abundance4...
)
# compare observed and predicted fitness
plot2D_versus(
    input=pre_nor,
    assay="assay_name",
    block=1
)
```
At the same time, there are also some more intuitive graphics to reflect the reliability of the model's predictions.
``` r
# abundance fitness vs. folding energy change
plot2D_fold(
    input=pre_nor,
    model="path/to/model/you want to test",
    block=1
)
ggplot2::ggsave("./plot2D_fold.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")

# binding fitness vs. binding and folding energy change
plot3D_bind(
    input=pre_nor,
    model="path/to/model/you want to test",
    assay="assay_name",
    block=1,
    output_file="./fig.pdf"
)
```
Don't forget to normalize again when using different assay_name.
