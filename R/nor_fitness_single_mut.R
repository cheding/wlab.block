#' A Function to get single aa mutants fitness data after normalization
#'
#' @description
#' This function allows you to get single aa mutants fitness data,
#' which are required for downstream  heatmap.
#'
#' @param input fitness data table after normalization
#' @return data.table with single aa mutants
#' @export
#' @import data.table
nor_fitness_single_mut<-function(
    input
){
  #select single aa mutants data
  input_single<-input[Nham_aa==1,]
  #overlap fitness normalization
  input_single[,nor_fitness_nooverlap:=sum(.SD[[1]]/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = c("nor_fitness","nor_fitness_sigma"),by="aa_seq"]
  input_single[,nor_fitness_nooverlap_sigma:=sqrt(1/sum(1/.SD[[1]]^2, na.rm = T)),.SDcols = c("nor_fitness_sigma"),by="aa_seq"]
  #overlap growth rate normalization
  input_single[,nor_gr_nooverlap:=sum(.SD[[1]]/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = c("nor_gr","nor_gr_sigma"),by="aa_seq"]
  input_single[,nor_gr_nooverlap_sigma:=sqrt(1/sum(1/.SD[[1]]^2, na.rm = T)),.SDcols = c("nor_gr_sigma"),by="aa_seq"]
  #remove duplicated mutant
  input_single<-input_single[!duplicated(aa_seq),]
  return(input_single)
}
