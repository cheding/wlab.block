#' A Function to plot3D binding fitness against folding ddG and binding ddG
#'
#' @description
#' This function allows you to plot binding fitness against folding ddG and binding ddG.
#'
#' @param input normalized predicted phenotypes data table
#' @param model MoCHI linears_weights model
#' @param assay assay name
#' @param block target block number
#' @param output_file plot output path
#'
#' @return Nothing
#' @export
#' @import data.table
plot3D_bind<-function(
    input,
    model,
    assay,
    block,
    output_file="./",
    RT = 0.001987*(273+30)
){
  mochi<-fread(model)

  pre_nor<-input
  bind_ddg1<-pre_nor[eval(parse(text=paste0("Binding",block,"_",assay)))==1,.(nt_seq,Nham_aa,aa_seq,fold_1_additive_trait1,fitness)]
  fold_ddg1<-pre_nor[eval(parse(text=paste0("Abundance",block)))==1,.(nt_seq,Nham_aa,aa_seq,fold_1_additive_trait0)]
  ddg_merge1<-merge(fold_ddg1,bind_ddg1,by=c("nt_seq","aa_seq","Nham_aa"))

  binding_range1<-ddg_merge1[Nham_aa>0,range(fold_1_additive_trait1)]

  binding_energy_grid1<-seq(binding_range1[1],
                            binding_range1[2],
                            (binding_range1[2]-binding_range1[1])/15)
  folding_range1<-ddg_merge1[Nham_aa>0,range(fold_1_additive_trait0)]
  folding_energy_grid1<-seq(folding_range1[1],
                            folding_range1[2],
                            (folding_range1[2]-folding_range1[1])/15)

  energy_grid_dt<-as.data.table(expand.grid(folding_energy_grid = folding_energy_grid1,
                                            binding_energy_grid = binding_energy_grid1))


  fraction_bound_fun<- function(
    folding_energy,
    binding_energy
  ){
    return(1/(1+exp(binding_energy/RT)*(1+exp(folding_energy/RT))))
  }
  fraction_bound<-fraction_bound_fun(folding_energy=energy_grid_dt[,folding_energy_grid],
                                     binding_energy = energy_grid_dt[,binding_energy_grid])

  fitness_binding<- fraction_bound * mochi[fold==1,kernel] + mochi[fold==1,bias]


  # fitness_binding
  pred_fitness_dt_binding <- data.table(
    f_dg_pred = energy_grid_dt[,folding_energy_grid],
    b_dg_pred = energy_grid_dt[,binding_energy_grid],
    observed_fitness = fitness_binding
  )

  Cairo::CairoPDF(file = output_file)
  plot3D::persp3D(x=binding_energy_grid1*RT,
                  y=folding_energy_grid1*RT,
                  z=matrix(data=pred_fitness_dt_binding[,observed_fitness],
                           nrow=length(binding_energy_grid1),ncol=length(folding_energy_grid1)),
                  r=2, shade=0.4, axes=TRUE,scale=TRUE, box=TRUE, nticks=5, ticktype="detailed",
                  colvar=F, col="white", alpha = 0, border="#F4270C", lwd=0.2,
                  cex.lab=1,cex.main=1,cex.axis=1,
                  xlab="ddG Folding",
                  ylab="ddG Binding",
                  zlab="Binding(observed)"
  )
  plot3D::scatter3D(x=ddg_merge1[,fold_1_additive_trait1*RT],
                    y=ddg_merge1[,fold_1_additive_trait0*RT],
                    z=ddg_merge1[,fitness],
                    add = T, col = "black", alpha = 0.2, cex = 0.2
  )
  dev.off()

}
