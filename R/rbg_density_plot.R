#' A function to draw density plots of false vresus true mutations
#'
#' @param input_dir path to "filter_variants.RData"
#' @param output_dir output path (default=input_dir)
#'
#' @return Nothing
#' @export
#' @import data.table
rbg_density_plot<-function(
    input_dir,
    output_dir=NULL
){
  if(is.null(output_dir)){
    output_dir=input_dir
  }
  #Load real backgrounds
  load(file.path(input_dir, "filter_variants.RData"))
  name<-as.data.frame(strsplit(names(mc_doubles_list), '_'))
  block<-unique(as.character(name[nrow(name),]))
  num<-length(mc_doubles_list)/length(block)
  #Real doubles based on real backgrounds
  plot_dt <- rbindlist(mc_doubles_list)
  plot_dt[, block_plot := paste0("block = ", substr(phenotype, nchar(phenotype), nchar(phenotype)))]
  plot_dt[, phenotypes_plot := substr(phenotype, 0, nchar(phenotype)-2)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(x = log10(mean_count+1), color = real_bg)) +
    ggplot2::geom_density() +
    ggplot2::xlab("log10(mean input count + 1)") +
    ggplot2::ylab("Density") +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(block_plot~phenotypes_plot)
  ggplot2::ggsave(file.path(output_dir, "real_backgrounds_mean_input_count.pdf"), d, width = num+2, height = length(block)+1, useDingbats=FALSE)

  plot_class<-plot_dt[,.(sum_count=sum(mean_count)),by=.(phenotype,real_bg)]
  p <- ggplot2::ggplot(plot_class,ggplot2::aes(x = factor(phenotype),y=sum_count, fill =factor(real_bg))) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::xlab("experiment") +
    ggplot2::ylab("Percentage of reads with mutations") +
    ggplot2::labs(fill=NULL)+
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(output_dir, "real_backgrounds_mean_input_percentage.pdf"), p, width = num+2.5, height = 4, useDingbats=FALSE)

  #Real doubles based on NNK
  plot_dt <- rbindlist(mc_nnk_list)
  plot_dt[, block_plot := paste0("block = ", substr(phenotype, nchar(phenotype), nchar(phenotype)))]
  plot_dt[, phenotypes_plot := substr(phenotype, 0, nchar(phenotype)-2)]
  d <- ggplot2::ggplot(plot_dt[Nham_aa==2],ggplot2::aes(x = log10(mean_count+1), color = real_bg)) +
    ggplot2::geom_density() +
    ggplot2::xlab("log10(mean input count + 1)") +
    ggplot2::ylab("Density") +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(block_plot~phenotypes_plot)
  ggplot2::ggsave(file.path(output_dir, "real_nnk_doubles_mean_input_count.pdf"), d, width = num+2, height = length(block)+1, useDingbats=FALSE)

  plot_class<-plot_dt[Nham_aa==2,.(sum_count=sum(mean_count)),by=.(phenotype,real_bg)]
  p <- ggplot2::ggplot(plot_class,ggplot2::aes(x = factor(phenotype),y=sum_count, fill =factor(real_bg))) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::xlab("experiment") +
    ggplot2::ylab("Percentage of reads with mutations") +
    ggplot2::labs(fill=NULL)+
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(output_dir, "real_nnk_doubles_mean_input_percentage.pdf"), p, width = num+2.5, height = 4, useDingbats=FALSE)

  #Real singles based on NNK
  plot_dt <- rbindlist(mc_nnk_list)
  plot_dt[, block_plot := paste0("block = ", substr(phenotype, nchar(phenotype), nchar(phenotype)))]
  plot_dt[, phenotypes_plot := substr(phenotype, 0, nchar(phenotype)-2)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(x = log10(mean_count+1), color = real_bg)) +
    ggplot2::geom_density() +
    ggplot2::xlab("log10(mean input count + 1)") +
    ggplot2::ylab("Density") +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(block_plot~phenotypes_plot)
  ggplot2::ggsave(file.path(output_dir, "real_nnk_singles_mean_input_count.pdf"), d, width = num+2, height = length(block)+1, useDingbats=FALSE)

  plot_class<-plot_dt[Nham_aa==1,.(sum_count=sum(mean_count)),by=.(phenotype,real_bg)]
  p <- ggplot2::ggplot(plot_class[Nham_aa==1],ggplot2::aes(x = factor(phenotype),y=sum_count, fill =factor(real_bg))) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::xlab("experiment") +
    ggplot2::ylab("Percentage of reads with mutations") +
    ggplot2::labs(fill=NULL)+
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(output_dir, "real_nnk_singles_mean_input_percentage.pdf"), p, width = num+2.5, height = 4, useDingbats=FALSE)
}
