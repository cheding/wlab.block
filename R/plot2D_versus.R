#' A Function to plot2D folding predicted fitness against observed fitness
#'
#' @description
#' This function allows you to plot folding predicted fitness against observed fitness.
#'
#' @param input normalized predicted phenotypes data table
#' @param assay assay name
#' @param block target block number
#'
#' @return Nothing
#' @export
#' @import data.table
plot2D_versus<-function(
    input,
    assay,
    block
){
  pre_nor<-input

  if(assay=="Abundance"){
    nb=0
  }else{
    nb<-which(colnames(pre_nor)==paste0("Binding1_",assay))-which(colnames(pre_nor)=="Abundance1")
  }

  phenotypen<-block+nb

  lm_mochi<-lm(pre_nor_fitness~ob_nor_fitness,pre_nor[phenotype==phenotypen,])
  ggplot2::ggplot()+
    ggplot2::stat_binhex(data=pre_nor[phenotype==phenotypen,],ggplot2::aes(x=ob_nor_fitness,y=pre_nor_fitness),
                         bins = 50,linewidth=0,color="black") +
    ggplot2::scale_fill_gradient(low="white",high="black",trans="log10",guide = ggplot2::guide_colorbar(barwidth = 0.5,barheight = 1.5)) +
    ggplot2::geom_hline(yintercept=0)+
    ggplot2::geom_vline(xintercept=0)+
    ggplot2::geom_abline(intercept = 0,slope=1,linetype="dashed")+
    ggplot2::annotate("text",x=-0.8,y=0.3,
                      label = paste0("R\u00B2 = ",round(summary(lm_mochi)$r.squared,2)),
                      size=7*0.35 )+
    ggplot2::theme_classic()+
    ggplot2::xlab("Observed fitness")+
    ggplot2::ylab("Predicted fitness")+
    ggplot2::theme(text = ggplot2::element_text(size=7),
                   axis.text = ggplot2::element_text(size=7),
                   legend.text = ggplot2::element_text(size=7),
                   legend.key.size = ggplot2::unit(1, "cm"),
                   plot.title = ggplot2::element_text(size=7))+
    ggplot2::coord_fixed()
}
