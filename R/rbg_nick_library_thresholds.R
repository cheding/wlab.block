#' A function to determine mutation background in nick library based on frequency
#'
#' @param input_files_thre libraries used for filtering
#' @param input_dir rbg_nick_library output directory
#' @param output_dir output path (default=input_dir)
#' @param block block information of input_files_thre (default=c('1','2','3'))
#' @param min_obs minimum permitted mutation reads
#' @param min_phenotypes minimum permitted frequency in different library
#'
#' @return Nothing
#' @export
#' @import data.table
rbg_nick_library_thresholds<-function(
    input_files_thre,
    input_dir,
    output_dir=NULL,
    block=c('1','2','3'),
    min_obs=200,
    min_phenotypes=NULL
){
  input_files<- file.path(unlist(input_files_thre), sapply(input_files_thre, list.files, pattern = "*Q20_variant_data_merge.tsv"))
  #Load singles list
  load(file.path(input_dir, "singles_list.RData"))
  load(file.path(input_dir, "rgb_dt.RData"))
  if(is.null(output_dir)){
    output_dir=input_dir
  }
  num<-ncol(singles_list[[1]])-1
  if(is.null(min_phenotypes)){
    min_phenotypes=num
  }
  #Number of backgrounds based on thresholds
  rbg_list <- list()
  for(i in block){
    singles_dt <- singles_list[[i]]
    rbg_list[[i]] <- singles_dt[apply(as.data.frame(singles_dt[,.SD,,.SDcols = !grepl('id', names(singles_dt))]>min_obs), 1, sum, na.rm = T)>=min_phenotypes,'id']
  }

  ###########################
  ### Filter variants
  ###########################

  #Load count tables
  count_list <- lapply(input_files, 'fread')
  names(count_list) <- names(input_files_thre)

  #Loop over all count tables
  mc_doubles_list <- list()
  mc_nnk_list <- list()
  for(i in 1:length(count_list)){
    print(i)
    block <- rev(unlist(strsplit(names(count_list)[i], '_')))[1]
    final_backgrounds <- rbg_list[[block]][,id]
    #WT
    wt <- count_list[[i]][WT==T]
    #Annotate singles (including synonymous mutations)
    singles <- rbg_get_singles(copy(count_list[[i]]))
    #Mean counts of singles with/without G/T in 3rd position
    mc_nnk_list <- c(mc_nnk_list,
                     list(singles[,.(mean_count = mean_count, real_bg = (!Mut1_codon %in% final_backgrounds & substr(Mut1_codon, nchar(Mut1_codon), nchar(Mut1_codon)) %in% c('g', 't')), phenotype = names(count_list)[i], Nham_aa = Nham_aa)]))
    #Filter singles for those matching NNK (or real backgrounds)
    singles <- singles[((Mut1_codon %in% final_backgrounds) | (!Mut1_codon %in% final_backgrounds & substr(Mut1_codon, nchar(Mut1_codon), nchar(Mut1_codon)) %in% c('g', 't')))]
    #Annotate doubles
    doubles <- rbg_get_doubles(copy(count_list[[i]]))
    #Mean counts of doubles with/without real background
    mc_doubles_list[[names(count_list)[i]]] <- doubles[,.(mean_count = mean_count, real_bg = (Mut1_codon %in% final_backgrounds | Mut2_codon %in% final_backgrounds), phenotype = names(count_list)[i])]
    #Filter doubles for those with real backgrounds
    doubles <- doubles[Mut1_codon %in% final_backgrounds | Mut2_codon %in% final_backgrounds]
    #Mean counts of doubles with/without G/T in 3rd position
    mc_nnk_list <- c(mc_nnk_list, list(doubles[,.(mean_count = mean_count, real_bg = ((Mut1_codon %in% final_backgrounds & substr(Mut2_codon, nchar(Mut2_codon), nchar(Mut2_codon)) %in% c('g', 't')) | (Mut2_codon %in% final_backgrounds & substr(Mut1_codon, nchar(Mut1_codon), nchar(Mut1_codon)) %in% c('g', 't'))), phenotype = names(count_list)[i], Nham_aa = Nham_aa)]))
    #Filter doubles for those matching NNK
    doubles <- doubles[((Mut1_codon %in% final_backgrounds & substr(Mut2_codon, nchar(Mut2_codon), nchar(Mut2_codon)) %in% c('g', 't')) | (Mut2_codon %in% final_backgrounds & substr(Mut1_codon, nchar(Mut1_codon), nchar(Mut1_codon)) %in% c('g', 't')))]
    #Merge with remaining variants
    all_variants <- rbind(
      wt,
      singles[,.SD,,.SDcols = names(count_list[[i]])],
      doubles[,.SD,,.SDcols = names(count_list[[i]])])
    #Remove unnecessary columns
    all_variants <- all_variants[,.SD,,.SDcols = names(all_variants)[grepl('nt_seq|input|output', names(all_variants))]]
    #Rename columns
    names(all_variants) <- gsub(".*input", "input", names(all_variants))
    names(all_variants) <- gsub(".*output", "output", names(all_variants))
    names(all_variants) <- c("nt_seq", sapply(strsplit(names(all_variants)[2:length(names(all_variants))], '_'), '[', 1))
    #Save count table
    write.table(all_variants, file = file.path(output_dir, paste0("variantCounts_", names(count_list)[i], '.tsv')), sep = "\t", row.names = F, quote = F)
  }

  #Save
  save(mc_doubles_list, mc_nnk_list, file = file.path(output_dir, "filter_variants.RData"))
}
