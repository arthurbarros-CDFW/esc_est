##########################################
#B_star
##########################################
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