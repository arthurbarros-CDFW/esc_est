#' @title B_star
#'
#' @description Calculates R_j, the number of carcasses released with tags during period j, 
#' and n_j the number of carcasses captured during j. Utilizes the Horvitz_Thompson() 
#' function to estimate N_j. Uses B1 and s_hat to calculate Bstar_j, the number
#' of births between j and j+1.
#' 
#' @param ch matrix of capture histories.
#' 
#' @param s_hat matrix of survival probability estimates for each carcass i at
#'  sampling period j.
#' 
#' @param p_hat matrix of capture probability estimates for each carcass i at
#'  sampling period j.
#' 
#' @param subsampling_weeks matrix of weeks in which subsampling occurred.
#' 
#' @param observed_per_week matrix of number of carcasses observed each week. Derived from 
#' number tagged each week * sub.sampling.
#' 
#' @param tagged_per_week matrix of number of carcasses tagged each week. Derived
#' as colSums(ch>=1,na.rm=T).
#' 
#' @return Bstar Estimate of Bstar, or number of births between j and j+1.
#' 
B_star<-function(ch,s_hat,p_hat,
                 subsampling_weeks,observed_per_week,tagged_per_week){
  R=n=list()
  nan=nrow(ch)
  ns=ncol(ch)
  for(j in 1:ns){
    R[j]<-as.numeric(length(which(ch[j]==1))) #carcasses released with tags
    
    if(!is.null(observed_per_week) && j<=length(observed_per_week)){
      n[j] <- observed_per_week[j] #incorporate subsampling if applicable
    } else {
      n[j]<-as.numeric(length(which(ch[j]==1))+length(which(ch[j]==2))) #total captured carcasses
    }
  }
  
  N_hat<-Horvitz_Thompson(p_hat,ch,subsampling_weeks,
                          observed_per_week,
                          tagged_per_week)
  #next B1, or total number of births for each period
  B1<-list()
  for(j in 2:(ns-2)){
    B1[j]<-N_hat[[j+1]]-mean(s_hat[,j])*(N_hat[[j]]-(n[[j]]-R[[j]]))
  }
  
  #next Bstar, number of births adjusted for those entering the system between
  #j and j+1, but not surviving to j+1
  Bstar<-NULL
  for(j in 2:(ns-2)){
    Bstar[j]<-as.numeric(B1[[j]]*(log(mean(s_hat[,j]))/(mean(s_hat[,j])-1)))
  }
  Bstar<-Bstar[-(1)]
  return(Bstar)
}