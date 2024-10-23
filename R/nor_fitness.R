#' A Function to normalize fitness from different blocks
#'
#' @description
#' This function allows you to normalize fitness from
#' different blocks dimsum output data.
#'
#' @param ... paths to "replicate.RData"
#'
#' @return data.table with normalized fitness and sigma
#' @export
#' @import data.table
nor_fitness <- function(...) {
  required_file <- list(...)

  all_data<-list()
  gr_lm<-data.frame(matrix(nrow = 2,ncol = 0))
  fit_lm<-data.frame(matrix(nrow = 2,ncol = 0))
  for(i in 1:length(required_file)){
    load(required_file[i][[1]])
    ######
    # growth rate
    ######
    #normalize growth rate
    all_variants$gr_over_sigmasquared<-all_variants$growthrate/(all_variants$growthrate_sigma)**2
    all_variants$one_over_sigmasquared<-1/(all_variants$growthrate_sigma)**2
    #average stop mutants growth rate
    dead_gr<-all_variants[STOP==T,]
    stop1<-sum(dead_gr$gr_over_sigmasquared, na.rm = TRUE)/sum(dead_gr$one_over_sigmasquared, na.rm = TRUE)
    #average WT_aa growth rate
    wt_gr<-all_variants[WT==T,]
    wt1<-sum(wt_gr$gr_over_sigmasquared, na.rm = TRUE)/sum(wt_gr$one_over_sigmasquared, na.rm = TRUE)

    gr_lm<-dplyr::bind_cols(gr_lm,c(stop=stop1,wt=wt1))
    colnames(gr_lm)[i]<-names(required_file)[i]
    ######
    # fitness
    ######
    #normalize fitness
    all_variants$fitness_over_sigmasquared<-all_variants$fitness/(all_variants$sigma)**2
    all_variants$one_over_fitness_sigmasquared<-1/(all_variants$sigma)**2
    #average stop mutants fitness
    dead_fitness<-all_variants[STOP==T,]
    stop1_fitness<-sum(dead_fitness$fitness_over_sigmasquared, na.rm = TRUE)/sum(dead_fitness$one_over_fitness_sigmasquared, na.rm = TRUE)
    #average WT_aa fitness
    wt_fitness<-all_variants[WT==T,]
    wt1_fitness<-sum(wt_fitness$fitness_over_sigmasquared, na.rm = TRUE)/sum(wt_fitness$one_over_fitness_sigmasquared, na.rm = TRUE)

    fit_lm<-dplyr::bind_cols(fit_lm,c(stop=stop1_fitness,wt=wt1_fitness))
    colnames(fit_lm)[i]<-names(required_file)[i]
    #combine
    all_data[[names(required_file[i])]]<-all_variants
  }
  #normalize block1
  all_data[[1]][,nor_gr:=growthrate]
  all_data[[1]][,nor_gr_sigma:=growthrate_sigma]
  all_data[[1]][,nor_fitness:=fitness]
  all_data[[1]][,nor_fitness_sigma:=sigma]
  #normalize other blocks
  if(length(required_file)>1){
    for(i in 2:length(required_file)){
      #growth rate lm parameter
      formula1 <- as.formula(paste0(colnames(gr_lm)[1],"~", colnames(gr_lm)[i]))
      b1b2<-lm(formula = formula1,data = gr_lm)
      b1b2<-summary(b1b2)
      a2<-b1b2$coefficients[[2]]
      b2<-b1b2$coefficients[[1]]
      #fitness lm parameter
      formula2 <- as.formula(paste0(colnames(fit_lm)[1],"~", colnames(fit_lm)[i]))
      c1c2<-lm(formula = formula2,data = fit_lm)
      c1c2<-summary(c1c2)
      d2<-c1c2$coefficients[[2]]
      e2<-c1c2$coefficients[[1]]
      #normalize block[i]
      all_data[[i]][,nor_gr:=growthrate*a2+b2]
      all_data[[i]][,nor_gr_sigma:=growthrate_sigma*a2]
      all_data[[i]][,nor_fitness:=fitness*d2+e2]
      all_data[[i]][,nor_fitness_sigma:=sigma*d2]
    }
  }
  data_after_nor<-dplyr::bind_rows(all_data,.id='block')
  return(data_after_nor)
}


