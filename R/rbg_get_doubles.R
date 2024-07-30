#' A Function to find double mutations from DiMSum "variant_data_merge.tsv"
#'
#' @param input_dt data table of "variant_data_merge.tsv"
#'
#' @return data table with double mutations
#' @export
#' @import data.table
rbg_get_doubles <- function(input_dt){
  #Mean count
  input_dt[, mean_count := rowMeans(.SD),,.SDcols = grepl('input', names(input_dt))]
  #WT seqs
  wt_ntseq <- input_dt[WT==T,nt_seq]
  wt_ntseq_split <- strsplit(input_dt[WT==T,nt_seq],"")[[1]]
  wt_AAseq_split <- strsplit(input_dt[WT==T,aa_seq],"")[[1]]
  #Double AA mutants
  doubles <- input_dt[Nham_aa==2]
  #Add position, mutant AA, WT AA
  doubles[,Pos1 := unique(ceiling(which(strsplit(nt_seq,"")[[1]] !=wt_ntseq_split)/3))[1],nt_seq]
  doubles[,Pos2 := unique(ceiling(which(strsplit(nt_seq,"")[[1]] !=wt_ntseq_split)/3))[2],nt_seq]
  doubles[,Mut1 := strsplit(aa_seq,"")[[1]][Pos1],aa_seq]
  doubles[,Mut2 := strsplit(aa_seq,"")[[1]][Pos2],aa_seq]
  doubles[,WT_AA1 := wt_AAseq_split[Pos1],aa_seq]
  doubles[,WT_AA2 := wt_AAseq_split[Pos2],aa_seq]
  #Add mutant codon, WT codon
  doubles[,Mut1_codon := paste0(Pos1, '_', substr(nt_seq,(Pos1-1)*3+1, (Pos1-1)*3+3)),nt_seq]
  doubles[,Mut2_codon := paste0(Pos2, '_', substr(nt_seq,(Pos2-1)*3+1, (Pos2-1)*3+3)),nt_seq]
  doubles[,WT_AA1_codon := paste0(Pos1, '_', substr(wt_ntseq,(Pos1-1)*3+1, (Pos1-1)*3+3)),nt_seq]
  doubles[,WT_AA2_codon := paste0(Pos2, '_', substr(wt_ntseq,(Pos2-1)*3+1, (Pos2-1)*3+3)),nt_seq]
  return(doubles)
}
