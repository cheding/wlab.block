#' A Function to plot2D folding fitness against folding ddG
#'
#' @description
#' This function allows you to plot folding fitness against folding ddG.
#'
#' @param input normalized predicted phenotypes data table
#' @param model MoCHI linears_weights model
#' @param block target block number
#'
#' @return plot about folding fitness against folding ddG
#' @export
#' @import data.table
plot2D_fold<-function(
    input,
    model,
    block
){
  mochi<-fread(model)

  pre_nor<-input
  pre_nor_input<-pre_nor[phenotype==block,]

  folding_range<-pre_nor_input[Nham_aa>0,range(fold_1_additive_trait0*RT,na.rm=T)]
  folding_energy_grid<-seq(folding_range[1],
                           folding_range[2],
                           (folding_range[2]-folding_range[1])/500)

  folding_fraction<-as.numeric(1/(1+exp(folding_energy_grid/RT)))
  fitness_folding<- folding_fraction * mochi[fold==1,kernel] + mochi[fold==1,bias]

  pred_fitness_dt <- data.table(
    f_dg_pred = folding_energy_grid,
    observed_fitness = fitness_folding)

  ggplot2::ggplot()+
    ggplot2::stat_binhex(data=pre_nor_input[Nham_aa>0,],ggplot2::aes(x=fold_1_additive_trait0*RT,y=fitness),bins=50,size=0,color="black")+
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::geom_line(data = pred_fitness_dt,ggplot2::aes(x=f_dg_pred,y=observed_fitness), color = "red")+
    ggplot2::geom_vline(xintercept=0,size=0.1,color="black")+
    ggplot2::geom_hline(yintercept=0,size=0.1,color="black")+
    ggpubr::theme_classic2()+
    ggplot2::ylab("Abundance(observed)")+
    ggplot2::xlab("Folding ddG (inferred)")+
    ggplot2::theme(text = ggplot2::element_text(size=7),
                   legend.position="right",
                   legend.text = ggplot2::element_text(size=7),
                   axis.text.x = ggplot2::element_text(size =7, vjust=.5, hjust=.5),
                   axis.text.y = ggplot2::element_text(size=7, vjust = .5,hjust = .5,margin=ggplot2::margin(0,-0.5,0,0,"mm")),
                   legend.key.height= ggplot2::unit(3.1, 'mm'),
                   legend.key.width = ggplot2::unit(3.1, 'mm'),
                   legend.key.size = ggplot2::unit(1,"mm"),
                   plot.margin=ggplot2::margin(0,0,0,0))
}
