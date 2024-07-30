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
  ggplot2::ggsave(file.path(output_dir, "real_backgrounds_mean_input_count.pdf"), d, width = 6, height = 4, useDingbats=FALSE)

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
  ggplot2::ggsave(file.path(output_dir, "real_nnk_doubles_mean_input_count.pdf"), d, width = 6, height = 4, useDingbats=FALSE)

  #Real singles based on NNK
  plot_dt <- rbindlist(mc_nnk_list)
  plot_dt[, block_plot := paste0("block = ", substr(phenotype, nchar(phenotype), nchar(phenotype)))]
  plot_dt[, phenotypes_plot := substr(phenotype, 0, nchar(phenotype)-2)]
  d <- ggplot2::ggplot(plot_dt[Nham_aa==1],ggplot2::aes(x = log10(mean_count+1), color = real_bg)) +
    ggplot2::geom_density() +
    ggplot2::xlab("log10(mean input count + 1)") +
    ggplot2::ylab("Density") +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(block_plot~phenotypes_plot)
  ggplot2::ggsave(file.path(output_dir, "real_nnk_singles_mean_input_count.pdf"), d, width = 6, height = 4, useDingbats=FALSE)
}
