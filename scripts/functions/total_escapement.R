#' @title total_escapement
#'
#' @description Produces estimate of total escapement for survey period.
#' 
#' @param ch matrix of capture histories.
#' 
#' @param beta list of four model parameters drawn from cpp_optim.
#' 
#' @param cap_X matrix of capture related covariates.
#' 
#' @param surv_X matrix of survival related covariates.
#' 
#' @param ints confidence intervals, defaults to 95%.
#' 
#' @param subsampling_weeks matrix of weeks in which subsampling occured.
#' 
#' @param observed_per_week matrix of number of carcasses observed each week. Derived from 
#' number tagged each week * sub.sampling.
#' 
#' @param tagged_per_week matrix of number of carcasses tagged each week. Derived
#' as colSums(ch>=1,na.rm=T).
#' 
#' @return ans list with p_hat and s_hat matrices, N_hat for each sampling period,
#' Bstar, and estimate of total escapement for survey period.
#' 

total_escapement<-function(ch,beta,cap_X,
                           surv_X,ints,
                           subsampling_weeks,
                           observed_per_week,
                           tagged_per_week){
  
  p_hat<-fill_prob_matrices(ch,beta,
                            cap_X,
                            surv_X,ints)$p_hat
  s_hat<-fill_prob_matrices(ch,beta,
                            cap_X ,
                            surv_X,ints)$s_hat
  
  N_hat=Horvitz_Thompson(p_hat,ch,subsampling_weeks,
                         observed_per_week,
                         tagged_per_week)
  
  Bstar=B_star(ch,s_hat,p_hat,
               subsampling_weeks,observed_per_week,tagged_per_week)
  
  escapement<-N_hat[[2]]*(log(mean(s_hat[,1]))/(mean(s_hat[,1])-1)) + 
    sum(Bstar,na.rm=T)
  
  #very rarely we get values of s_hat that are 1, which results in
  #Bstar of NA, and escapement of NA, so if that is the case,
  #replace super population adjustment with the mean of N_hat
  if(is.na(escapement)){
    escapement=mean(unlist(N_hat))
  }
  
  ans<-list("p_hat"=p_hat,"s_hat"=s_hat,"N_hat"=N_hat,
            "Bstar"=Bstar,"escapement"=escapement)
  
  return(ans)
}