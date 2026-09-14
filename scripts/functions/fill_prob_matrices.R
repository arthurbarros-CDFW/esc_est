##########################################
#fill_prob_matrices
##########################################
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