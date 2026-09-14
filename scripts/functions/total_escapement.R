##########################################
#Total escapement
##########################################
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
  
  ans<-list("p_hat"=p_hat,"s_hat"=s_hat,"N_hat"=N_hat,
            "Bstar"=Bstar,"escapement"=escapement)
  
  return(ans)
}