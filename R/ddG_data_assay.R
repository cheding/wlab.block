#' A Function to obtain statistical values of interest from ddG data
#'
#'#' @description
#' This function allows you to get statistical values of interest  from ddG data.
#'
#' @param input path to MoCHI output ".txt" file
#' @param wt_aa wt amino acid sequence
#'
#' @return data table with summarized statistical values
#' @export
#' @import data.table
ddG_data_assay<-function(
    input,
    wt_aa
){
  ddG<-fread(input)
  num<-nchar(wt_aa)+1
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))

  #get real mutant position
  ddG[,Pos_real:=Pos_ref+1]
  ddG[id!="WT",wt_codon:=substr(id,1,1)]
  ddG[id!="WT",mt_codon:=substr(id,nchar(id),nchar(id))]
  ddG[,mt:=paste0(wt_codon,Pos_real,mt_codon)]

  #create assistant data table
  heatmap_tool<-data.table(wt_codon = rep(unlist(strsplit(wt_aa,"")),each=20),
                           Pos_real = rep(2:num,each=20),
                           mt_codon = unlist(aa_list))
  ddG<-merge(ddG,heatmap_tool,by=c("Pos_real","wt_codon","mt_codon"),all=T)

  #get statistical values of interest
  codon<-ddG[Pos_real>1,unique(.SD[[1]]),.SDcols = c("wt_codon"),by="Pos_real"]
  setnames(codon,"V1","codon")
  mean<-ddG[Pos_real>1,sum(.SD[[1]]/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = c("mean_kcal/mol","std_kcal/mol"),by="Pos_real"]
  setnames(mean,"V1","mean")
  abs_mean<-ddG[Pos_real>1,sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = c("mean_kcal/mol","std_kcal/mol"),by="Pos_real"]
  setnames(abs_mean,"V1","abs_mean")
  sigma<-ddG[Pos_real>1,sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),.SDcols = c("mean_kcal/mol","std_kcal/mol"),by="Pos_real"]
  setnames(sigma,"V1","sigma")
  max<-ddG[Pos_real>1&!is.na(`mean_kcal/mol`),max(.SD[[1]]),.SDcols = c("mean_kcal/mol"),by="Pos_real"]
  setnames(max,"V1","max")
  min<-ddG[Pos_real>1&!is.na(`mean_kcal/mol`),min(.SD[[1]]),.SDcols = c("mean_kcal/mol"),by="Pos_real"]
  setnames(min,"V1","min")
  count<-ddG[Pos_real>1&!is.na(`mean_kcal/mol`),nrow(.SD),.SDcols = c("mean_kcal/mol"),by="Pos_real"]
  setnames(count,"V1","count")

  #summarize and return data
  output<-merge(codon,mean,by="Pos_real",all=T)
  output<-merge(output,abs_mean,by="Pos_real",all=T)
  output<-merge(output,sigma,by="Pos_real",all=T)
  output<-merge(output,max,by="Pos_real",all=T)
  output<-merge(output,min,by="Pos_real",all=T)
  output<-merge(output,count,by="Pos_real",all=T)
  return(output)
}
