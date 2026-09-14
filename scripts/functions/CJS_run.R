##########################################
#CJS_run
##########################################
CJS_run<-function(){
  #ask for data if needed, check if already exists
  refreshData <- TRUE
  if (exists("prepped_data")){
    new.data <- select.list(c("Keep the 'CURRENT' data file sets",
                              "Load 'NEW' data file sets"),
                            multiple = FALSE,
                            graphics = FALSE,
                            title = "Which Data Files to Use?")
    refreshData <- substr(new.data,1,4)=="Load"
  } 
  
  
  #run CJS_data_prep()
  if(refreshData==TRUE){
    prepped_data<-CJS_data_prep()
    assign("prepped_data", prepped_data, pos=.GlobalEnv) # save for next round, possibly
  }
  
  #select one or more model
  models = c(	"constant capture and survival rates",
              "constant capture rate and survival related the sex",
              "constant capture rate and survival related the length",
              "capture related to sex and constant survival rate",
              "capture related to length and constant survival rate",
              "capture related to sex and survival related to length",
              "capture related to length and survival related to sex",
              "capture related to sex and survival related to sex",
              "capture related to length and survival related to length")
  model_fit<-CJS_model_select(prepped_data$covars_ask,
                              prepped_data$sex_matrix,
                              prepped_data$lengths_matrix,
                              prepped_data$ch)
  
  
  n_boot<-get_bootstrap_iterations()
  
  c_hat = 1
  initial_beta=numeric(4)
  model_results=data.frame()
  
  for(i in model_fit$models_ran){
    #est_escapement = NaN when i=6 (cap~sex, surv~length)
    #this seems to be because the fill_prob_matrices is producing all 1s for length covariates?
    #this is produced in the pro_capsur function
    #seems to be because of exp(cap_beta[1]*1+cap_beta[2]*cap_X[i,j])
    #the beta parameters are way high when using length
    #this only seems to happen for this model variant for this data?
    #note it doesn't happen with Trent's application
    #skip for now (3/26/2025)
    print(paste("beginning model:",i," (",
                models[i],")",sep=""))
    model_num=i
    starttime<-Sys.time()
    {gc()
      optim_results<-cpp_optim(beta=initial_beta,
                               ch=as.matrix(prepped_data$ch),
                               cap_X=as.matrix(model_fit$cap_X[[i]]),
                               surv_X=as.matrix(model_fit$surv_X[[i]]),
                               ints=as.matrix(prepped_data$intervals))
    }
    endtime<-Sys.time()
    optim_speed<-endtime-starttime 
    
    ans<-total_escapement(ch=prepped_data$ch,
                          beta=optim_results$par,
                          cap_X=model_fit$cap_X[[i]],
                          surv_X=model_fit$surv_X[[i]],
                          ints=as.matrix(prepped_data$intervals),
                          subsampling_weeks=prepped_data$subsampling_weeks,
                          observed_per_week=prepped_data$observed_per_week,
                          tagged_per_week=prepped_data$tagged_per_week)
    est_escapement<-ans$escapement
    
    p_hat<-ans$p_hat
    s_hat<-ans$s_hat
    ch=prepped_data$ch
    beta=optim_results$par
    cap_X=model_fit$cap_X[[i]]
    surv_X=model_fit$surv_X[[i]]
    ints=as.matrix(prepped_data$intervals)
    
    ######################################################
    #calculate model fit statistics (aic, qaic, qaicc)
    ######################################################
    loglik<-optim_results$value
    ic<-prepped_data$ch
    nan=nrow(ic)
    ns=ncol(ic)
    
    fit_results <- cjs_fit_simple(ch, beta, cap_X, surv_X, ints)
    c_hat<-fit_results
    #df<-fit_results$idfgt
    n_params<-length(beta)
    
    AIC=2*loglik+2*n_params
    AICc=AIC+((2*n_params)*(n_params+1))/(nan-n_params-1)
    QAIC=((2*loglik)/c_hat)+(2*n_params)
    QAICC=QAIC+(2*n_params*(n_params+1))/(nan-n_params-1)
    
    
    if(n_boot!=0){
      boot_start<-Sys.time()
      boot_results<-CJS_bootstrap(iterations=n_boot,
                                  ch=ic,
                                  cap_X=model_fit$cap_X[[i]],
                                  surv_X=model_fit$surv_X[[i]],
                                  ints=prepped_data$intervals,
                                  subsampling_weeks=prepped_data$subsampling_weeks,
                                  observed_per_week=prepped_data$observed_per_week,
                                  tagged_per_week=prepped_data$tagged_per_week)
      boot_end<-Sys.time()
      boot_speed=boot_end-boot_start
      #confidence intervals
      conf.level = 95
      alpha = 1 - conf.level/100
      lower = alpha/2
      upper = 1 - alpha/2
      mid=.5
      ci<-boot_results%>%
        summarise(lower_ci=ceiling(quantile(escapement,probs = c(lower),na.rm=T)),
                  mid_ci=ceiling(quantile(escapement,probs = c(mid),na.rm=T)),
                  upper_ci=ceiling(quantile(escapement,probs = c(upper),na.rm=T)))
      p<-ggplot(boot_results,aes(x=escapement))+
        geom_histogram(color = "#000000", fill = "#0099F8")+
        geom_segment(data=ci,aes(x=lower_ci,xend=lower_ci,y=0,yend=Inf),
                     linewidth=1,linetype='dashed')+
        geom_segment(data=ci,aes(x=upper_ci,xend=upper_ci,y=0,yend=Inf),
                     linewidth=1,linetype='dashed')+
        geom_segment(data=ci,aes(x=est_escapement,
                                 xend=est_escapement,y=0,yend=Inf),
                     linewidth=1,linetype='dashed',color='red')+
        #scale_x_continuous(breaks = seq(0,10000,500)) +
        theme_classic()
      ggsave(paste("outputs/plot_model-",i,"_iter-",n_boot,"_",Sys.Date(),".png",
                   sep=""),p,scale=4)
    }
    
    d<-data.frame(est_escapement,
                  optim_speed,
                  cap_beta1=optim_results$par[1],
                  cap_beta2=optim_results$par[2],
                  surv_beta1=optim_results$par[3],
                  surv_beta2=optim_results$par[4],
                  loglik,
                  model=models[i],
                  AIC,
                  AICc,
                  QAIC,
                  QAICC,
                  c_hat
    )
    if(n_boot!=0){
      d<-data.frame(d,lower_ci=ci$lower_ci,
                    upper_ci=ci$upper_ci,
                    n_boot,
                    boot_speed)
    }
    model_results<-model_results%>%
      rbind(d)
    
  }
  
  model_results <- model_results %>%
    mutate(
      delta_AIC = AIC - min(AIC),
      delta_AICc = AICc - min(AICc),
      AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC)),
      AICc_weight = exp(-0.5 * delta_AICc) / sum(exp(-0.5 * delta_AICc))
    ) %>%
    arrange(AIC)  # Sort by best model first
  
  write.csv(model_results,paste("outputs/CJS_outputs-",
                                Sys.Date(),
                                ".csv",
                                sep=""),
            row.names = F)
  return(model_results)
}