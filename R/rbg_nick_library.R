#' A Function to obtain mutations distribution from nick library
#'
#' @param input_dir_list list contain path to "variant_data_merge.tsv"
#' @param output_dir output path (default="./")
#' @param block block information of input_dir_list (default=c('1','2','3'))
#'
#' @return Nothing
#' @export
#' @import data.table
rbg_nick_library<-function(
    input_dir_list,
    output_dir="./",
    block=c('1','2','3')
){
  input_files <- file.path(unlist(input_dir_list), sapply(input_dir_list, list.files, pattern = "*Q20_variant_data_merge.tsv"))
  dir.create(output_dir)
  num<-length(input_dir_list)/length(block)-1
  #Load count tables
  count_list <- lapply(input_files, 'fread')
  names(count_list) <- names(input_dir_list)

  ###########################
  ### Singles frequency in doubles
  ###########################
  #Singles list
  min_mean_count <- 20
  singles_list = list()
  for(i in block){
    print(paste0('block', i))
    block_i = grep(paste0('_', i, '$'), names(count_list))
    #Get table of singles in doubles for all phenotypes
    singles_dt <- data.table()
    for(j in block_i){
      print(j)
      doubles <- rbg_get_doubles(count_list[[j]])
      bg_tab <- table(unlist(doubles[mean_count>min_mean_count,.(Mut1_codon, Mut2_codon)]))
      temp_dt <- data.table(
        id = names(bg_tab),
        count = as.integer(bg_tab))
      names(temp_dt)[2] <- names(count_list)[j]
      if(nrow(singles_dt)==0){
        singles_dt <- temp_dt
      }else{
        singles_dt <- merge(singles_dt, temp_dt, by = 'id', all = T)
      }
    }
    #Save
    singles_list[[i]] <- singles_dt
  }

  #Save
  save(singles_list, file = file.path(output_dir, "singles_list.RData"))

  ###########################
  ### Real backgrounds
  ###########################

  #Load singles list
  load(file.path(output_dir, "singles_list.RData"))

  #Number of backgrounds based on thresholds
  rbg_list <- list()
  for(i in block){
    print(paste0('block', i))
    singles_dt <- singles_list[[i]]
    for(j in 0:num){
      for(min_obs in 11:300){
        rbg_list <- c(rbg_list,
                      list(data.frame(
                        count = sum(apply(as.data.frame(singles_dt[,.SD,,.SDcols = !grepl('id', names(singles_dt))]>min_obs), 1, sum, na.rm = T)>j),
                        min_obs = min_obs,
                        min_phenotypes = j+1,
                        block = i)))
      }
    }
  }
  rgb_dt <- rbindlist(rbg_list)

  #Save
  save(rgb_dt, file = file.path(output_dir, "rgb_dt.RData"))

  ###########################
  ### Plot #real backgrounds vs. minimum observations in doubles
  ###########################

  #Load real backgrounds
  load(file.path(output_dir, "rgb_dt.RData"))

  plot_dt <- rgb_dt
  plot_dt[, block_plot := paste0("block = ", block)]
  plot_dt[, min_phenotypes_plot := paste0("min_phen. = ", min_phenotypes)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(x = min_obs, y = count)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = 200, linetype = 2, col = 'red') +
    ggplot2::geom_hline(yintercept = 10, linetype = 3, col = 'grey') +
    ggplot2::geom_hline(yintercept = 20, linetype = 2, col = 'grey') +
    ggplot2::xlab("Minimum observed doubles") +
    ggplot2::ylab("High confidence backgrounds") +
    ggplot2::coord_cartesian(ylim = c(1, 1000)) +
    ggplot2::scale_y_continuous(trans='log10') +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(block_plot~min_phenotypes_plot)
  ggplot2::ggsave(file.path(output_dir, "real_backgrounds_lineplots.pdf"), d, width = num+3, height = num+1, useDingbats=FALSE)
}
