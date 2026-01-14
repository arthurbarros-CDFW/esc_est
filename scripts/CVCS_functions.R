library(knitr)
library(tinytex)
library(tidyverse)
library(ggplot2)
library(Rcpp)
sourceCpp('scripts/CJS_functions.cpp')

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

##########################################
#CJS_data_prep
##########################################
CJS_data_prep<-function(){
  suppressWarnings(rm("ch","chops","covars"))#remove objects
  #Part 1: select data prompts
  #ask for CH
  ch<-read.csv(choose.files(caption="Select your CAPTURE HISTORY data (.csv file)"), head=T, as.is=T)
  
  #ask for chops
  chops.ask<-select.list(c(	"YES", 
                 "NO"), 
              graphics = FALSE, 
              multiple = FALSE, title = "Do you have CHOPS on 1st Capture data NOT represented as '2' in capture histories?")
  if(chops.ask == "YES"){
    chops<-read.csv(choose.files(caption="Select your CHOPS data (.csv file)"), head=T, as.is=T)
  }
  #ask for covars
  covars.ask<-select.list(c(	"YES", 
                            "NO"), 
                         graphics = FALSE, 
                         multiple = FALSE, title = "Do you have COVARIATE data?")
  if(covars.ask == "YES"){
    covars<-read.csv(choose.files(caption="Select your COVARIATE data (.csv file)"), head=T, as.is=T)
    #standardize field names
    names(covars) <- tolower(names(covars))
    
    if( nrow(covars) < nrow(ch) ){
      dropcovars.ask<-select.list(c(	"YES", 
                                "NO"), 
                             graphics = FALSE, 
                             multiple = FALSE, title = "The number of covariate records does not match the number of capture history records. Would you like to drop capture history records missing covariate records?")
      if(dropcovars.ask == "YES"){
        missing_in_covars <- ch$disctag[!ch$disctag %in% covars$disctag]
        ch<-ch%>%
          filter(disctag!=missing_in_covars)
      } else {
        stop("Number of covariate records must be the same as the number of capture history records")
      }
    } else if(nrow(covars)>nrow(ch)){
      stop("Number of covariate records must be the same as the number of capture history records")
    }
    
    #check covars for na values
    if(any(is.na(covars$length))){
      missing_lengths<-covars%>%filter(is.null(length)|is.na(length))
      message(paste("Found", sum(is.na(covars$length)), 
                    "NA values in length. These rows will be removed."))
      # Remove rows with NA length
      covars <- covars%>%filter(!disctag %in% missing_lengths$disctag)
      
      # Also remove corresponding rows from ch
      ch <- ch%>%filter(!disctag %in% missing_lengths$disctag)
    }
    if(any(is.na(covars$sex))){
      missing_sex<-covars%>%filter(is.null(sex)|is.na(sex))
      message(paste("Found", sum(is.na(covars$sex)), 
                    "NA values in sex. These rows will be removed."))
      # Remove rows with NA length
      covars <- covars%>%filter(!disctag %in% missing_sex$disctag)
      
      # Also remove corresponding rows from ch
      ch <- ch%>%filter(!disctag %in% missing_sex$disctag)
    }
  } else {
    sex <- NULL
    length <- NULL
  }
  
  #ask for unequal timing
  #see Mrawin.f90 (https://github.com/tmcd82070/MRA/blob/master/src/Mrawin.f90) 
  #line 1270 incorporated in prosur
  ints.equal <- select.list(c("YES",
                              "NO"),
                            multiple = FALSE,
                            graphics = FALSE, 
                            title = "Are all intervals between occasions equal?")
  if(ints.equal == "NO"){
    ints = read.csv(choose.files(caption="Select the file containing lengths of un-equal time intervals (.csv file)"),
                    head=T, as.is=T)
    ints <- ints[,grep("[0-9]",names(ints))]
    ints <- unlist(ints)  # must be a vector, not a data frame
    # at this point ch may have other columns in it
    if( length(ints) != (ncol(ch)-2) ){
      stop("Number of time intervals must be 1 less than number of sampling occasions")
    }
    if( any(ints == 0) ){
      stop("Zero time intervals are not allowed.")
    }
  } else {
    ints = rep(1,ncol(ch)-2)
  }
  
  #ask about subsampling
  #this is driven by LAR having to subsample every other fish in high return years
  sampling.unequal <- select.list(c("YES",
                              "NO"),
                            multiple = FALSE,
                            graphics = FALSE, 
                            title = "Where there any periods when subsampling was performed?")
  if(sampling.unequal == "YES"){
    sub.sampling = read.csv(choose.files(caption="Select the file containing records of subsampling (.csv file)"),
                    head=T, as.is=T)
    sub.sampling <- sub.sampling[,grep("[0-9]",names(sub.sampling))]
    sub.sampling <- unlist(sub.sampling)  # must be a vector, not a data frame
    # at this point ch may have other columns in it
    if( length(sub.sampling) != (ncol(ch)-1) ){
      stop("Number of records must equal number of sampling periods.")
    }
    if( any(ints == 0) ){
      stop("Zero sampling periods is not allowed.")
    }
  } else {
    sub.sampling = rep(1,ncol(ch)-1)
  }
  
  #Part 2: prepare ch and covars data
  ch=ch[-1] #remove disctag vector from capture histories
  null_matrix<-matrix(1,nrow=nrow(ch),ncol=ncol(ch))
  
  if(covars.ask == "YES"){
    covars$sex<-as.numeric(ifelse(covars$sex%in%c('F',"f"),1,0)) #change sex to numeric value 
    covars$length<-as.numeric(covars$length)
  }
  
  #Part 3: deal with subsampling periods
  tagged_per_week<-colSums(ch>=1,na.rm=T)
  observed_per_week <- tagged_per_week*sub.sampling
  
  #Part 4: prep chops data
  if(chops.ask == "YES"){
    chops<-chops[-1]
    clean_chops<-matrix(ncol=ncol(chops))
    for(i in 1:ncol(chops)){
      d<-chops[i]
      r<-rep(0,ncol(chops))
      r[i]=2
      r<-matrix(rep((r),d),ncol=ncol(chops),byrow=T)
      clean_chops<-clean_chops%>%rbind(r)
    }
    clean_chops<-clean_chops[-1, ,drop=FALSE]
    colnames(clean_chops)<-colnames(ch)
  }
  
  #Part 5: generate covariate data for chops
  #here we'll do things differently than the escapeMR app
  #instead of taking the mean values for sex and length and assigning them to chops
  #we will randomly sample existing covars and assign to chop data
  #I like this better because it incorporates variability into the chops
  #consider including this in the bootstrapping somehow?
  if(covars.ask == "YES" & chops.ask == "YES"){
    n<-nrow(clean_chops)
    chops_covars<-covars[sample(nrow(covars),n,replace=T),] #sample for number of chops
    covars<-covars%>%rbind(chops_covars)
  }
  
  if(chops.ask == "YES"){
    #add in chops to ch
    ch<-ch%>%rbind(clean_chops)
  }
   if(covars.ask == "YES"){
     lengths_matrix<-matrix(covars$length,
                            nrow=nrow(ch),ncol=ncol(ch))
     sex_matrix<-matrix(as.numeric(covars$sex),
                            nrow=nrow(ch),ncol=ncol(ch))
   } else {
     lengths_matrix<-null_matrix
     sex_matrix<-null_matrix
   }
  
  return(list('ch'=ch,'lengths_matrix'=lengths_matrix,
              'sex_matrix'=sex_matrix,"intervals"=ints,
              'covars_ask'=covars.ask,
              'tagged_per_week'=tagged_per_week,
              'observed_per_week'=observed_per_week,
              'sampling.unequal'=sampling.unequal,
              'subsampling_weeks'=sub.sampling))
}

##########################################
#CJS_model_select
##########################################
CJS_model_select<-function(covars.ask,sex_matrix,lengths_matrix,ch){
  models = c(	"constant capture and survival rates",
              "constant capture rate and survival related the sex",
              "constant capture rate and survival related the length",
              "capture related to sex and constant survival rate",
              "capture related to length and constant survival rate",
              "capture related to sex and survival related to length",
              "capture related to length and survival related to sex",
              "capture related to sex and survival related to sex",
              "capture related to length and survival related to length")
  model.rank = c(1:9)
  if(covars.ask == "YES"){
    models.to.fit <-select.list(models, 
                                multiple = TRUE, 
                                graphics = TRUE,
                                title = "Select one or more models to fit")
    models.to.fit <- which( models %in% models.to.fit )
  } else {
    models.to.fit <-c(1)
  }
  model_covars = list(cap_X=NULL, surv_X=NULL)
  for(m in models.to.fit){
    if(m==1){
      model_covars$cap_X[[m]] = matrix(1,nrow=nrow(ch),ncol=ncol(ch))
      model_covars$surv_X[[m]]= matrix(1,nrow=nrow(ch),ncol=ncol(ch))
    }
    if(m==2){
      model_covars$cap_X[[m]] = matrix(1,nrow=nrow(ch),ncol=ncol(ch))
      model_covars$surv_X[[m]]= sex_matrix
    }
    if(m==3){
      model_covars$cap_X[[m]] = matrix(1,nrow=nrow(ch),ncol=ncol(ch))
      model_covars$surv_X[[m]]= lengths_matrix
    }
    if(m==4){
      model_covars$cap_X[[m]] = sex_matrix
      model_covars$surv_X[[m]]= matrix(1,nrow=nrow(ch),ncol=ncol(ch))
    }
    if(m==5){
      model_covars$cap_X[[m]] = lengths_matrix
      model_covars$surv_X[[m]]= matrix(1,nrow=nrow(ch),ncol=ncol(ch))
    }
    if(m==6){
      model_covars$cap_X[[m]] = sex_matrix
      model_covars$surv_X[[m]]= lengths_matrix
    }
    if(m==7){
      model_covars$cap_X[[m]] = lengths_matrix
      model_covars$surv_X[[m]]= sex_matrix
    }
    if(m==8){
      model_covars$cap_X[[m]] = sex_matrix
      model_covars$surv_X[[m]]= sex_matrix
    }
    if(m==9){
      model_covars$cap_X[[m]] = lengths_matrix
      model_covars$surv_X[[m]]= lengths_matrix
    }
  }
  return(list('models_ran'=models.to.fit,
              'ch'=ch,
              'cap_X'=model_covars$cap_X,
              'surv_X'=model_covars$surv_X))
}

CJS_loglik_wrapper <- function(beta, cap_X, surv_X, ch) {
  -CJS_loglik(beta, cap_X, surv_X, ch)  # Return the negative log-likelihood
}

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

##########################################
#Total escapement
##########################################
total_escapement<-function(ch,beta,cap_X,surv_X,ints,
                           subsampling_weeks,observed_per_week,
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

##########################################
#CJS_bootstrap
##########################################
get_bootstrap_iterations <- function() {
  #ask if user wants to bootstrap
  perform_boot <- tolower(readline("Perform bootstrapping for confidence intervals? (y/n): "))
  
  #validate yes/no response
  while(!perform_boot %in% c("y", "n", "yes", "no")) {
    message("Error: Please answer 'y' or 'n'.")
    perform_boot <- tolower(readline("Perform bootstrapping? (y/n): "))
  }
  
  if (perform_boot %in% c("n", "no")) {
    return(0)  #return 0 if no bootstrapping
  }
  
  #if yes, get bootstrap iterations
  while(TRUE) {
    n_boot <- readline("Enter bootstrap iterations (integer 25-1000): ")
    
    #check if numeric
    if (is.na(suppressWarnings(as.numeric(n_boot)))) {
      message("Error: Input must be a number.")
      next
    }
    
    n_boot <- as.integer(n_boot)
    
    #check if integer
    if (is.na(n_boot) || n_boot != as.numeric(n_boot)) {
      message("Error: Input must be an integer.")
      next
    }
    
    # Check range
    if (n_boot < 25 || n_boot > 1000) {
      message("Error: Input must be between 25 and 1000.")
      next
    }
    
    return(n_boot)
  }
}

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

##########################################
#model selection for shiny app
##########################################
#wrapper function for 
CJS_model_select_app <- function(covars_used, sex_matrix, lengths_matrix, ch, selected_models) {
  models = c(
    "constant capture and survival rates",
    "constant capture rate and survival related to sex",
    "constant capture rate and survival related to length",
    "capture related to sex and constant survival rate",
    "capture related to length and constant survival rate",
    "capture related to sex and survival related to length",
    "capture related to length and survival related to sex",
    "capture related to sex and survival related to sex",
    "capture related to length and survival related to length"
  )
  
  if (covars_used) {
    if (is.null(selected_models)) {
      # Default to first model if none selected but covariates are available
      models.to.fit <- 1
    } else {
      models.to.fit <- which(models %in% selected_models)
    }
  } else {
    # If no covariates, only allow the first model
    models.to.fit <- 1
  }
  
  model_covars = list(cap_X = NULL, surv_X = NULL)
  
  for (m in models.to.fit) {
    if (m == 1) {
      model_covars$cap_X[[m]] = matrix(1, nrow = nrow(ch), ncol = ncol(ch))
      model_covars$surv_X[[m]] = matrix(1, nrow = nrow(ch), ncol = ncol(ch))
    }
    if (m == 2) {
      model_covars$cap_X[[m]] = matrix(1, nrow = nrow(ch), ncol = ncol(ch))
      model_covars$surv_X[[m]] = sex_matrix
    }
    if (m == 3) {
      model_covars$cap_X[[m]] = matrix(1, nrow = nrow(ch), ncol = ncol(ch))
      model_covars$surv_X[[m]] = lengths_matrix
    }
    if (m == 4) {
      model_covars$cap_X[[m]] = sex_matrix
      model_covars$surv_X[[m]] = matrix(1, nrow = nrow(ch), ncol = ncol(ch))
    }
    if (m == 5) {
      model_covars$cap_X[[m]] = lengths_matrix
      model_covars$surv_X[[m]] = matrix(1, nrow = nrow(ch), ncol = ncol(ch))
    }
    if (m == 6) {
      model_covars$cap_X[[m]] = sex_matrix
      model_covars$surv_X[[m]] = lengths_matrix
    }
    if (m == 7) {
      model_covars$cap_X[[m]] = lengths_matrix
      model_covars$surv_X[[m]] = sex_matrix
    }
    if (m == 8) {
      model_covars$cap_X[[m]] = sex_matrix
      model_covars$surv_X[[m]] = sex_matrix
    }
    if (m == 9) {
      model_covars$cap_X[[m]] = lengths_matrix
      model_covars$surv_X[[m]] = lengths_matrix
    }
  }
  
  return(list(
    'models_ran' = models.to.fit,
    'ch' = ch,
    'cap_X' = model_covars$cap_X,
    'surv_X' = model_covars$surv_X
  ))
}

##########################################
# cjs_fit for c_hat model fitting
##########################################
#purpose: Calculate goodness-of-fit for CJS models using Test 2 and Test 3
#returns: c-hat (vif), chi-square statistic, and degrees of freedom
#inputs:
#ng=1 number of "groups"
#ic=capture history (ch)
#ig=dimensions of groups
#min_expected= minimum expected count for chi-square


cjs_fit <- function(ch, beta, cap_X, surv_X, ints, min_expected = 1) {
  
  # Calculate fitted probabilities from your model
  prob_matrices <- fill_prob_matrices(ch, beta, cap_X, surv_X, ints)
  p_hat <- prob_matrices$p_hat
  s_hat <- prob_matrices$s_hat
  
  nan <- nrow(ch)
  ns <- ncol(ch)
  
  # ------------------------------------------------------------------
  # TEST 2: Compare observed vs expected m-array
  # ------------------------------------------------------------------
  
  # Create OBSERVED m-array
  m_obs <- matrix(0, nrow = ns, ncol = ns)
  releases <- numeric(ns)
  
  for (i in 1:nan) {
    for (j in 1:(ns-1)) {
      if (ch[i, j] >= 1) {
        releases[j] <- releases[j] + 1
        found_next <- FALSE
        for (k in (j+1):ns) {
          if (ch[i, k] >= 1) {
            m_obs[j, k] <- m_obs[j, k] + 1
            found_next <- TRUE
            break
          }
        }
      }
    }
  }
  
  # Create EXPECTED m-array from model predictions
  m_exp <- matrix(0, nrow = ns, ncol = ns)
  
  for (j in 1:(ns-1)) {
    if (releases[j] > 0) {
      for (i in 1:nan) {
        if (ch[i, j] >= 1) {
          # Probability of surviving from j to k and being captured at k
          # but not captured between j and k
          for (k in (j+1):ns) {
            surv_prob <- 1
            for (t in j:(k-1)) {
              surv_prob <- surv_prob * s_hat[i, t]
            }
            cap_prob <- p_hat[i, k]
            not_cap_prob <- 1
            for (t in (j+1):(k-1)) {
              not_cap_prob <- not_cap_prob * (1 - p_hat[i, t])
            }
            m_exp[j, k] <- m_exp[j, k] + surv_prob * cap_prob * not_cap_prob
          }
        }
      }
    }
  }
  
  # Calculate Test 2 chi-square
  chi_sq_2 <- 0
  df_2 <- 0
  
  if (ns >= 4) {
    for (j in 2:(ns-2)) {
      # Compare newly marked vs previously marked
      prev_marked_obs <- sum(m_obs[1:(j-1), (j+1):ns])
      new_marked_obs <- sum(m_obs[j, (j+1):ns])
      prev_marked_exp <- sum(m_exp[1:(j-1), (j+1):ns])
      new_marked_exp <- sum(m_exp[j, (j+1):ns])
      
      # Skip if expected counts too small
      if (prev_marked_exp < min_expected || new_marked_exp < min_expected) {
        next
      }
      
      # Chi-square contribution (simplified - actual would be more detailed)
      chi_sq_2 <- chi_sq_2 + (prev_marked_obs - prev_marked_exp)^2 / prev_marked_exp
      chi_sq_2 <- chi_sq_2 + (new_marked_obs - new_marked_exp)^2 / new_marked_exp
      df_2 <- df_2 + 1
    }
  }
  
  # ------------------------------------------------------------------
  # TEST 3: Compare observed vs expected capture histories
  # ------------------------------------------------------------------
  
  chi_sq_3 <- 0
  df_3 <- 0
  
  for (j in 2:(ns-1)) {
    # Create 2x2 table: captured at j vs captured at j+1
    obs_table <- matrix(0, nrow = 2, ncol = 2)
    exp_table <- matrix(0, nrow = 2, ncol = 2)
    
    for (i in 1:nan) {
      # Only animals alive at j
      alive_at_j <- FALSE
      for (t in j:ns) {
        if (ch[i, t] >= 1) {
          alive_at_j <- TRUE
          break
        }
      }
      
      if (alive_at_j) {
        captured_j <- as.numeric(ch[i, j] >= 1)
        captured_j1 <- as.numeric(ch[i, j+1] >= 1)
        
        obs_table[captured_j + 1, captured_j1 + 1] <- 
          obs_table[captured_j + 1, captured_j1 + 1] + 1
        
        # Expected probability
        if (captured_j == 1) {
          # Animal captured at j
          exp_prob <- s_hat[i, j] * p_hat[i, j+1]
        } else {
          # Animal not captured at j but could be captured at j+1
          exp_prob <- s_hat[i, j] * p_hat[i, j+1] * (1 - p_hat[i, j])
        }
        
        # Simplified expected counts
        exp_table[captured_j + 1, captured_j1 + 1] <- 
          exp_table[captured_j + 1, captured_j1 + 1] + exp_prob
      }
    }
    
    # Calculate chi-square if expected counts sufficient
    if (sum(exp_table) > min_expected * 4) {
      for (r in 1:2) {
        for (c in 1:2) {
          if (exp_table[r, c] > min_expected) {
            chi_sq_3 <- chi_sq_3 + (obs_table[r, c] - exp_table[r, c])^2 / exp_table[r, c]
          }
        }
      }
      df_3 <- df_3 + 1
    }
  }
  
  # ------------------------------------------------------------------
  # Calculate c-hat
  # ------------------------------------------------------------------
  
  total_chi_sq <- chi_sq_2 + chi_sq_3
  total_df <- df_2 + df_3
  
  if (total_df > 0) {
    c_hat <- total_chi_sq / total_df
    c_hat <- max(1.0, c_hat)  # c-hat should be at least 1
  } else {
    c_hat <- 1.0
  }
  
  return(list(
    vif = c_hat,
    chigt = total_chi_sq,
    idfgt = total_df,
    diagnostic = list(
      test2_chi_sq = chi_sq_2,
      test2_df = df_2,
      test3_chi_sq = chi_sq_3,
      test3_df = df_3
    )
  ))
}

# ------------------------------------------------------------------
# Helper function for TEST 2
# ------------------------------------------------------------------
calculate_test2 <- function(ns, m, releases, min_expected, use_chat_rot, chat_rot) {
  
  chi_sq <- 0
  df <- 0
  
  # For each release occasion j from 2 to ns-2
  for (j in 2:(ns-2)) {
    
    # Create contingency table for release occasion j
    # Rows: Previously marked (releases before j) vs Newly marked (releases at j)
    # Columns: Recaptured at occasions j+1 to ns
    
    # Previously marked: sum of releases before j that were recaptured
    prev_marked <- numeric(ns - j)
    new_marked <- numeric(ns - j)
    
    for (k in (j+1):ns) {
      # Previously marked: animals released before j and recaptured at k
      prev_marked[k-j] <- sum(m[1:(j-1), k])
      # Newly marked: animals released at j and recaptured at k
      new_marked[k-j] <- m[j, k]
    }
    
    # Row totals
    prev_total <- sum(prev_marked)
    new_total <- sum(new_marked)
    col_totals <- prev_marked + new_marked
    
    # Column totals (total recaptures at each occasion)
    total_animals <- prev_total + new_total
    
    # Skip if insufficient data
    if (total_animals < min_expected) {
      next
    }
    
    # Calculate expected values and chi-square
    for (k in 1:length(prev_marked)) {
      # Expected counts under independence
      exp_prev <- prev_total * col_totals[k] / total_animals
      exp_new <- new_total * col_totals[k] / total_animals
      
      # Skip if expected count too small
      if (exp_prev < min_expected || exp_new < min_expected) {
        next
      }
      
      # Chi-square contribution
      chi_sq <- chi_sq + ((prev_marked[k] - exp_prev)^2 / exp_prev)
      chi_sq <- chi_sq + ((new_marked[k] - exp_new)^2 / exp_new)
      df <- df + 1
    }
    
    # Adjust degrees of freedom (for each column we lose 1 df)
    df <- df - 1  # Because column totals sum to row totals
  }
  
  return(list(chi_sq = chi_sq, df = max(0, df)))
}

# ------------------------------------------------------------------
# Helper function for TEST 3
# ------------------------------------------------------------------
calculate_test3 <- function(nan, ns, ic, min_expected, use_chat_rot, chat_rot) {
  
  chi_sq <- 0
  df <- 0
  
  # For each occasion j from 2 to ns-1
  for (j in 2:(ns-1)) {
    
    # Create 2x2 contingency table for occasion j
    # Rows: Captured at j (yes/no)
    # Columns: Captured at j+1 (yes/no)
    # This tests if capture probability depends on previous capture
    
    contingency <- matrix(0, nrow = 2, ncol = 2)
    
    for (i in 1:nan) {
      # Only consider animals alive at j (captured at or after j)
      if (any(ic[i, j:ns] >= 1)) {
        captured_j <- as.numeric(ic[i, j] >= 1)
        captured_j1 <- as.numeric(ic[i, j+1] >= 1)
        
        contingency[captured_j + 1, captured_j1 + 1] <- 
          contingency[captured_j + 1, captured_j1 + 1] + 1
      }
    }
    
    # Skip if any cell has insufficient counts
    if (any(contingency < min_expected)) {
      next
    }
    
    # Calculate expected values
    row_totals <- rowSums(contingency)
    col_totals <- colSums(contingency)
    total <- sum(contingency)
    
    # Skip if totals are too small
    if (total < min_expected * 4) {  # All 4 cells need minimum
      next
    }
    
    # Calculate chi-square
    for (r in 1:2) {
      for (c in 1:2) {
        expected <- row_totals[r] * col_totals[c] / total
        if (expected >= min_expected) {
          chi_sq <- chi_sq + ((contingency[r, c] - expected)^2 / expected)
        }
      }
    }
    
    df <- df + 1  # Each 2x2 table contributes 1 degree of freedom
  }
  
  return(list(chi_sq = chi_sq, df = max(0, df)))
}

# ------------------------------------------------------------------
# Simplified wrapper function for model fitting
# ------------------------------------------------------------------
cjs_fit_simple <- function(ch, beta, cap_X, surv_X, ints) {
  # Simple wrapper that matches how you're calling it in CJS_run()
  result <- calculate_simple_c_hat(
    ch = as.matrix(ch),
    beta = beta,
    cap_X = as.matrix(cap_X),
    surv_X = as.matrix(surv_X),
    ints = as.matrix(ints)
  )
  return(result)
}

calculate_simple_c_hat <- function(ch, beta, cap_X, surv_X, ints) {
  nan <- nrow(ch)
  ns <- ncol(ch)
  
  # Get predicted probabilities
  probs <- fill_prob_matrices(ch, beta, cap_X, surv_X, ints)
  p_hat <- probs$p_hat
  
  # Calculate Pearson residuals
  residuals <- matrix(0, nan, ns)
  for (i in 1:nan) {
    for (j in 1:ns) {
      if (!is.na(p_hat[i, j]) && p_hat[i, j] > 0 && p_hat[i, j] < 1) {
        observed <- as.numeric(ch[i, j] >= 1)
        expected <- p_hat[i, j]
        residuals[i, j] <- (observed - expected) / sqrt(expected * (1 - expected))
      }
    }
  }
  
  # Simple c-hat estimate
  c_hat <- sum(residuals^2, na.rm = TRUE) / (nan * ns - length(beta))
  return(max(1.0, c_hat))
}
