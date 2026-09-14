CJS_bootstrap<-function(iterations,ch,cap_X,surv_X,ints,
                        subsampling_weeks,observed_per_week,
                        tagged_per_week,
                        progress_callback=NULL){
  #Part 1: set initial values
  results<-data.frame()
  initial_beta=numeric(4)
  
  #Part 2: run the for loop for each iteration
  for(r in 1:iterations){
    iter_start<-Sys.time()
    
    #Part 2.1: index and sample the capture histories
    index=1:dim(ch)[1]
    samp<- sample(index, replace = T)
    
    #Part 2.2: create new capture histories and
    #cap_X and surv_X matrices based on sampled indices
    ch_iteration = ch[samp,]
    cap_X_iteration=cap_X[samp,]
    surv_X_iteration=surv_X[samp,]
    
    #Part 2.3: use cpp_optim to find the optimizal beta
    #parameters for the given iteration data
    {gc()
      optim_iter<-cpp_optim(beta=initial_beta,
                            ch=as.matrix(ch_iteration),
                            cap_X=as.matrix(cap_X_iteration),
                            surv_X=as.matrix(surv_X_iteration),
                            ints=ints)
    }
    iter_beta<-optim_iter$par
    iter_lik<-optim_iter$value
    
    #Part 2.4: estimate total escapement for iteration
    est_esc_iter<-total_escapement(ch_iteration,
                                   iter_beta,
                                   cap_X_iteration,
                                   surv_X_iteration,
                                   ints,
                                   subsampling_weeks,observed_per_week,
                                   tagged_per_week)$escapement
    
    iter_end<-Sys.time()
    iter_time<-iter_end-iter_start
    
    #Part 2.5: store relevant iteration information in results
    d<-data.frame("iteration"=r,
                  "log-likelihood"=iter_lik,
                  "escapement"=est_esc_iter,
                  "time"=iter_time)
    results<-results%>%rbind(d)
    print(paste("bootstrap iteration: ",r)) #print iteration number for progress tracking
    if(!is.null(progress_callback)) {
      progress_callback(r)
    }
  }
  return(results)
}