##########################################
#pro_capsur
##########################################
pro_capsur<-function(i,j,ch, beta,cap_X,surv_X,ints){
  nan=nrow(ch)
  ns=ncol(ch)
  
  p=length(beta)
  
  #purpose: evaluate probability of capture and survival for each animal i
  cap_beta<-beta[1:(p/2)]
  surv_beta<-beta[((p/2)+1):p]
  
  zp<-exp(cap_beta[1]*1+cap_beta[2]*cap_X[i,j])
  zs<-exp(surv_beta[1]*1+surv_beta[2]*surv_X[i,j])
  
  p.hat<-zp/(1+zp)
  s.hat<-zs/(1+zs)
  
  s.hat<-s.hat**ints[j]
  
  est_list<-list('p.hat'=p.hat,'s.hat'=s.hat)
  return(est_list)
}