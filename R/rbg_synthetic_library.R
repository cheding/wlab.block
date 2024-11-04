#' A Function to remove false mutations based on mutation backgrounds in synthetic library
#'
#' @param input_dir_list list contain path to "variant_data_merge.tsv"
#' @param output_dir output path (default="./")
#' @param position list contain mutation background position
#' @param codon code type (default="NNK")
#'
#' @return Nothing
#' @export
#' @import data.table
rbg_synthetic_library<-function(
    input_dir_list,
    output_dir="./",
    position,
    codon="NNK"
){
  get_codon <- function(codon){
    code1<-substr(codon,1,1)
    code2<-substr(codon,2,2)
    code3<-substr(codon,3,3)
    codon<-expand.grid(unlist(nuc_codes[code1]),unlist(nuc_codes[code2]),unlist(nuc_codes[code3]))
    codon<- apply(codon, 1, function(x) paste(x, collapse = ""))
    return(codon)
  }

  get_backgrounds <- function(position, codon){
    backgrounds<-expand.grid(position, '_', codon)
    backgrounds<- apply(backgrounds, 1, function(x) paste(x, collapse = ""))
    return(backgrounds)
  }

  #IUPAC nucleic acid codon
  nuc_codes <- list(
    "A" = "a",
    "C" = "c",
    "G" = "g",
    "T" = "t",
    "R" = c("a", "g"),
    "Y" = c("c", "t"),
    "S" = c("c", "g"),
    "W" = c("a", "t"),
    "K" = c("g", "t"),
    "M" = c("a", "c"),
    "B" = c("c", "g", "t"),
    "D" = c("a", "g", "t"),
    "H" = c("a", "c", "t"),
    "V" = c("a", "c", "g"),
    "N" = c("a", "c", "g", "t")
  )

  input_files <- file.path(unlist(input_dir_list), sapply(input_dir_list, list.files, pattern = "*Q20_variant_data_merge.tsv"))
  dir.create(output_dir)
  #Load count tables
  count_list <- lapply(input_files, 'fread')
  names(count_list) <- names(input_dir_list)

  rbg_list <- list()
  #Codon range
  codon <- get_codon(codon)
  #Real backgrounds based on library design
  for (i in 1:length(count_list)) {
    print(names(count_list[i]))
    pos <- unlist(position[names(count_list[i])])
    backgrounds <- get_backgrounds(pos,codon)
    rbg_list[names(count_list[i])] <- as.data.frame(backgrounds)
  }

  #Save
  save(rbg_list, file = file.path(output_dir, "rbg_list.RData"))

  ###########################
  ### Filter variants
  ###########################

  mc_doubles_list <- list()
  mc_nnk_list <- list()
  for(i in 1:length(count_list)){
    print(names(count_list[i]))
    final_backgrounds <- unlist(rbg_list[names(count_list[i])])
    #WT
    wt <- count_list[[i]][WT==T]
    #Annotate singles (including synonymous mutations)
    singles <- rbg_get_singles(copy(count_list[[i]]))
    #Mean counts of singles with/without G/T in 3rd position
    mc_nnk_list <- c(mc_nnk_list,
                     list(singles[,.(mean_count = mean_count, real_bg = ((Mut1_codon %in% final_backgrounds) | ((!Mut1_codon %in% final_backgrounds) & (substr(Mut1_codon, nchar(Mut1_codon), nchar(Mut1_codon)) %in% c('g', 't')))), phenotype = names(count_list)[i], Nham_aa = Nham_aa)]))
    #Filter singles for those matching NNK (or real backgrounds)
    singles <- singles[((Mut1_codon %in% final_backgrounds) | ((!Mut1_codon %in% final_backgrounds) & (substr(Mut1_codon, nchar(Mut1_codon), nchar(Mut1_codon)) %in% c('g', 't'))))]
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
    # names(all_variants) <- gsub(".*input", "input", names(all_variants))
    # names(all_variants) <- gsub(".*output", "output", names(all_variants))
    names(all_variants) <- c("nt_seq", sapply(strsplit(names(all_variants)[2:length(names(all_variants))], '_'), '[', 1))
    #Save count table
    write.table(all_variants, file = file.path(output_dir, paste0("variantCounts_", names(count_list)[i], '.tsv')), sep = "\t", row.names = F, quote = F)
  }

  #Save
  save(mc_doubles_list, mc_nnk_list, file = file.path(output_dir, "filter_variants.RData"))
}
