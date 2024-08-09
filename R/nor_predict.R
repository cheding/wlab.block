#' A Function to normalize predicted fitness from MoCHI
#'
#'
#' @param pre_file "predicted_phenotypes_all.txt" from MoCHI
#' @param assay assay name
#' @param ... paths to "replicate.RData"
#'
#' @return data.table with predicted fitness
#' @export
#' @import data.table
nor_predict<-function(
    pre_file,
    assay="Abundance",
    ...
){
  required_file <- list(...)

  fit_lm<-data.frame(matrix(nrow = 2,ncol = 0))
  for(i in 1:length(required_file)){
    load(required_file[i][[1]])
    #fitness
    #normalize growth rate
    all_variants$fitness_over_sigmasquared<-all_variants$fitness/(all_variants$sigma)**2
    all_variants$one_over_fitness_sigmasquared<-1/(all_variants$sigma)**2
    #average stop mutants fitness
    dead_fitness<-all_variants[all_variants$STOP==T,]
    stop1_fitness<-sum(dead_fitness$fitness_over_sigmasquared, na.rm = TRUE)/sum(dead_fitness$one_over_fitness_sigmasquared, na.rm = TRUE)
    #average WT_aa fitness
    wt_fitness<-all_variants[all_variants$WT==T,]
    wt1_fitness<-sum(wt_fitness$fitness_over_sigmasquared, na.rm = TRUE)/sum(wt_fitness$one_over_fitness_sigmasquared, na.rm = TRUE)

    fit_lm<-dplyr::bind_cols(fit_lm,c(stop=stop1_fitness,wt=wt1_fitness))
    colnames(fit_lm)[i]<-names(required_file)[i]
  }


  pre_nor<-fread(pre_file)

  n1<-which(colnames(pre_nor)=="fold_1")
  n2<-which(colnames(pre_nor)=="Fold")
  if(assay=="Abundance"){
    nb=0
  }else{
    nb<-which(colnames(pre_nor)==paste0("Binding1_",assay))-which(colnames(pre_nor)=="Abundance1")
  }

  extract_prediction<-function(row){return(row[n1+as.numeric(row[n2])-1])}
  pre_nor$predicted_fitness<-apply(pre_nor,MARGIN = 1,FUN = extract_prediction)
  #create a new column for predicted fitness. In a matrix, 1 in MARGIN means rows
  pre_nor$predicted_fitness<-as.numeric(pre_nor$predicted_fitness)

  # if(assay=="Abundance"){
  #   extract_additive_trait0<-function(row){return(row[n2+as.numeric(row[n2])*2-1])}
  #   extract_additive_trait1<-function(row){return(row[n2+as.numeric(row[n2])*2])}
  #
  #   pre_nor$additive_trait0<-apply(pre_nor,MARGIN = 1,FUN = extract_additive_trait0)
  #   #create a new column for predicted fitness. In a matrix, 1 in MARGIN means rows
  #   pre_nor$additive_trait0<-as.numeric(pre_nor$additive_trait0)
  #   pre_nor$additive_trait1<-apply(pre_nor,MARGIN = 1,FUN = extract_additive_trait1)
  #   #create a new column for predicted fitness. In a matrix, 1 in MARGIN means rows
  #   pre_nor$additive_trait1<-as.numeric(pre_nor$additive_trait1)
  #   pre_nor[,additive_trait:=additive_trait0+additive_trait1]
  # }

  #normalize phenotype1
  pre_nor[phenotype==1+nb,pre_nor_mean_fitness:=mean]
  pre_nor[phenotype==1+nb,pre_nor_fitness_sigma:=std]
  pre_nor[phenotype==1+nb,ob_nor_fitness:=fitness]
  pre_nor[phenotype==1+nb,ob_nor_fitness_sigma:=sigma]
  pre_nor[phenotype==1+nb,pre_nor_fitness:=predicted_fitness]

  #normalize other blocks
  for(i in 2:length(required_file)){

    #fitness lm parameter
    formula2 <- as.formula(paste0(colnames(fit_lm)[1],"~", colnames(fit_lm)[i]))
    c1c2<-lm(formula = formula2,data = fit_lm)
    c1c2<-summary(c1c2)
    d2<-c1c2$coefficients[[2]]
    e2<-c1c2$coefficients[[1]]
    #normalize block[i]
    pre_nor[phenotype==i+nb,pre_nor_mean_fitness:=mean*d2+e2]
    pre_nor[phenotype==i+nb,pre_nor_fitness_sigma:=std*d2]
    pre_nor[phenotype==i+nb,ob_nor_fitness:=fitness*d2+e2]
    pre_nor[phenotype==i+nb,ob_nor_fitness_sigma:=sigma*d2]
    pre_nor[phenotype==i+nb,pre_nor_fitness:=predicted_fitness*d2+e2]

  }
  return(pre_nor)
}
