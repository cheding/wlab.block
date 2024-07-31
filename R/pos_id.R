#' A Function to find mutations from aa sequence
#'
#' @description
#' This function allows you to find mutation from dimsum output data.
#'
#' @param input merged dimsum output data table
#' @param wt_aa wt amino acid sequence
#'
#' @return data.table with mutation positions
#' @export
#' @import data.table
pos_id<-function(
    input,
    wt_aa
){
  output<-copy(input)
  num<-nchar(wt_aa)
  #obtain mutant locations
  output[,AA_Pos1 :=  which(unlist(strsplit(aa_seq, ""))!=unlist(strsplit(wt_aa, ""))[1:nchar(aa_seq)])[1],aa_seq]
  output[,AA_Pos2 :=  which(unlist(strsplit(aa_seq, ""))!=unlist(strsplit(wt_aa, ""))[1:nchar(aa_seq)])[2],aa_seq]
  #substr(s, first, last)
  for(i in 1:num){output[AA_Pos1==i,wtcodon1 := substr(wt_aa,i,i)]}
  for(i in 1:num){output[AA_Pos2==i,wtcodon2 := substr(wt_aa,i,i)]}
  output[,codon1 := substr(aa_seq,AA_Pos1,AA_Pos1)]
  output[,codon2 := substr(aa_seq,AA_Pos2,AA_Pos2)]
  #real mutant locations
  output[,AA_Pos1 :=  AA_Pos1+1]
  output[,AA_Pos2 :=  AA_Pos2+1]

  output[,mt1 := paste0(wtcodon1,AA_Pos1,codon1)]
  output[,mt2 := paste0(wtcodon2,AA_Pos2,codon2)]
  return(output)
}
