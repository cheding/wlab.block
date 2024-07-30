#' A Function to find single mutations from DiMSum "variant_data_merge.tsv"
#'
#' @param input_dt data table of "variant_data_merge.tsv"
#'
#' @return data table with single mutations
#' @export
#' @import data.table
rbg_get_singles <- function(
    input_dt
){
  #Mean count
  input_dt[, mean_count := rowMeans(.SD),,.SDcols = grepl('input', names(input_dt))]
  #WT seqs
  wt_ntseq <- input_dt[WT==T,nt_seq]
  wt_ntseq_split <- strsplit(input_dt[WT==T,nt_seq],"")[[1]]
  wt_AAseq_split <- strsplit(input_dt[WT==T,aa_seq],"")[[1]]
  #Single AA mutants
  singles <- input_dt[Nmut_codons==1]
  #Add position, mutant AA, WT AA
  singles[,Pos1 := ceiling(which(strsplit(nt_seq,"")[[1]] !=wt_ntseq_split)[1]/3),nt_seq]
  singles[,Mut1 := strsplit(aa_seq,"")[[1]][Pos1],aa_seq]
  singles[,WT_AA1 := wt_AAseq_split[Pos1],aa_seq]
  #Add mutant codon, WT codon
  singles[,Mut1_codon := paste0(Pos1, '_', substr(nt_seq,(Pos1-1)*3+1, (Pos1-1)*3+3)),nt_seq]
  singles[,WT_AA1_codon := paste0(Pos1, '_', substr(wt_ntseq,(Pos1-1)*3+1, (Pos1-1)*3+3)),nt_seq]
  return(singles)
}
