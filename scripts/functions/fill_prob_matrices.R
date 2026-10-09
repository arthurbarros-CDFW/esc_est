#' @title fill_prob_matrices
#'
#' @description builds matrices of survival and capture probabilities. Called by
#' total_escapement().
#' 
#' @param ch matrix of capture histories.
#' 
#' @param beta list of four model parameters.
#' 
#' @param cap_X matrix of capture related covariates.
#' 
#' @param surv_X matrix of survival related covariates.
#' 
#' @param ints confidence intervals, defaults to 95%.
#' 
#' @return list with matrices of estimate capture and survival probabilities for
#' each carcass at each observation point.
#' 

fill_prob_matrices<-function(ch,beta,cap_X,surv_X,ints){
  nan=nrow(ch)
  ns=ncol(ch)
  p_hat<-s_hat<-matrix( 0, nan, ns ) #create empty matrices
  for(i in 1:nan){
    for(j in 1:ns){ #fill each cell with pro_capsur() function
      p_hat[i,j]<-pro_capsur(i,j,ch,beta,
                             cap_X,
                             surv_X,ints)$p.hat
      s_hat[i,j]<-pro_capsur(i,j,ch,beta,
                             cap_X,
                             surv_X,ints)$s.hat
    }
  }
  p_hat[,1]<-NA
  
  est_list<-list('p_hat'=p_hat,'s_hat'=s_hat)
  return(est_list)
}