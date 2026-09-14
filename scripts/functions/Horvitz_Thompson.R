##########################################
#Horvitz_Thompson
##########################################
Horvitz_Thompson<-function(p_hat,ch,
                           subsampling_weeks,observed_per_week,
                           tagged_per_week){
  nan=nrow(ch)
  ns=ncol(ch)
  N_hat<-list()
  n_mat<-matrix(NA,nan,ns)
  for(j in 1:ns){
    for(i in 1:nan){
      n_mat[[i,j]]<-if(ch[[i,j]]>=1){
        1/p_hat[[i,j]]
      } else {0}
    }
    N_tagged<- sum(n_mat[, j], na.rm = TRUE)
    
    #if this is a subsampled week, adjust for skipped carcasses
    if(subsampling_weeks[j]>1){
      if(observed_per_week[j]>tagged_per_week[j]){
        skipped_count<-observed_per_week[j]-tagged_per_week[j]
        skipped_contribution<-matrix()
        for(s in 1:skipped_count){
          #for each skipped carcass, sample from p_hat and use to estimate the skipped contribution
          p_imputed<-sample(p_hat[,j],1)#sample by column j in case we eventuality figure out time covariates
          skipped_contribution[s]<-1/p_imputed
        }
      }
      N_hat[[j]]=N_tagged+sum(skipped_contribution)
    }else{
      N_hat[[j]]=N_tagged
    }
  }
  return(N_hat)
}