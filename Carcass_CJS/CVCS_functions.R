library(knitr)
library(tinytex)
library(tidyverse)
library(ggplot2)
library(Rcpp)
sourceCpp('CJS_functions.cpp')

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
  
  s.hat<-s.hat*ints[j]
  
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
  } else {
    sex <- NULL
    length <- NULL
  }
  
  #ask for unequal timing
  #this still needs to be incorporated into the CJS_run side?
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
    if( length(ints) != (length(grep("survey", names(ch)))-1) ){
      stop("Number of time intervals must be 1 less than number of occasions")
    }
    if( any(ints == 0) ){
      stop("Zero time intervals are not allowed.")
    }
  } else {
    ints = NULL
  }
  
  #Part 2: prepare ch and covars data
  ch=ch[-1] #remove disctag vector from capture histories
  null_matrix<-matrix(1,nrow=nrow(ch),ncol=ncol(ch))
  
  if(covars.ask == "YES"){
    covars$sex<-as.numeric(ifelse(covars$sex%in%c('F',"f"),1,0)) #change sex to numeric value 
  }
  
  #Part 3: prep chops data
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
    clean_chops<-clean_chops[-1,]
    colnames(clean_chops)<-colnames(ch)
  }
  
  #Part 4: generate covariate data for chops
  #here we'll do things differently than the escapeMR app
  #instead of taking the mean values for sex and length and assigning them to chops
  #we will randomly sample existing covars and assign to chop data
  #I like this better because it incorporates variability into the chops
  #consider including this in the bootstrapping somehow?
  if(covars.ask == "YES" & chops.ask == "YES"){
    n<-nrow(clean_chops)
    chops_covars<-as.matrix(covars[sample(nrow(covars),n,replace=T),]) #sample for number of chops
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
              'covars_ask'=covars.ask))
}

##########################################
#Model_fit
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
B_star<-function(ch,s_hat,p_hat){
  R=n=list()
  nan=nrow(ch)
  ns=ncol(ch)
  for(j in 1:ns){
    d<-ch #select just the capture matrix data
    R[j]<-as.numeric(length(which(d[j]==1)))
    n[j]<-as.numeric(length(which(d[j]==1))+length(which(d[j]==2)))
  }
  
  N_hat<-Horvitz_Thompson(p_hat,ch)
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
total_escapement<-function(ch,beta,cap_X,surv_X,ints){
  
  p_hat<-fill_prob_matrices(ch,beta,
                            cap_X,
                            surv_X,ints)$p_hat
  s_hat<-fill_prob_matrices(ch,beta,
                            cap_X ,
                            surv_X,ints)$s_hat
  
  N_hat=Horvitz_Thompson(p_hat,ch)
  
  Bstar=B_star(ch,s_hat,p_hat)
  
  escapement<-N_hat[[2]]*(log(mean(s_hat[,1]))/(mean(s_hat[,1])-1)) + 
    sum(Bstar,na.rm=T)
  
  return(ceiling(escapement))
}

##########################################
#Horvitz_Thompson
##########################################
Horvitz_Thompson<-function(p_hat,ch){
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
    N_hat[[j]]=sum(n_mat[,j])
  }
  return(N_hat)
}

##########################################
#cjs_fit
##########################################
#inputs:
#nan=number of animals or records
#ns=number of surveys
#ic=the capture history (nan*ns)
#ng= 1= number of parameters, based on group value from F.cjs.estim?
#however, in the esc_model.R script they don't assign group
#so in Trents work, group defaults to rep(1,nan) and ng=1
#ig=dimensions of groups=rep(1, nan)

#outputs:
#chigt=total chi-squared
#vif=c_hat, measure of overdispersion=chigt/idfgt
#idfgt=degrees of freedom for test
#ic=ch

cjs_fit <- function(nan, ns, ic, ng, ig) {
  # Constants
  chat_rot <- 5  #c_hat rule of thumb
  
  #initialize outputs
  vif <- 1.0
  chigt <- 0.0
  idfgt <- 0
  
  #test 2 "group by group"
  #this is artifact from MRA, probably not necessary but keeping
  #test 2 and 3 are covered in Amstrup et al handbook 2005
  
  for (l in 1:ng) {
    #initialize matrix array
    nang <- 0
    relese <- integer(ns)
    m <- matrix(0, nrow = ns, ncol = ns)
    
    #compute matrix array
    #this is essentially the chi-squared contingency table
    #amstrup 2005 page 236
    for (i in 1:nan) {
      if (ig[i] == l) {
        nang <- nang + 1
        for (j in 1:(ns-1)) {
          if (ic[i, j] >= 1) {
            relese[j] <- relese[j] + 1
            for (k in (j+1):ns) {
              if (ic[i, k] >= 1) {
                m[j, k] <- m[j, k] + 1
                break
              }
            }
          }
        }
      }
    }
    
    if (nang > 0) {
      test2_result <- cjs_test2(ns, m, chat_rot) #run test2 function
      chigt <- chigt + test2_result$chitot
      idfgt <- idfgt + test2_result$idftot
    }
  }
  
  #test 3 group by group
  for (l in 1:ng) {
    #calculate resight and subsequent components for releases 2 to NS-1
    for (j in 2:(ns-1)) {
      #initialize components for contingency table
      n11 <- n12 <- n21 <- n22 <- 0
      m11 <- m12 <- m21 <- m22 <- 0
      
      for (i in 1:nan) {
        if (ig[i] == l && ic[i, j] >= 1) {
          #calculate IB (how many times was animal seen before j)
          ib <- sum(ic[i, 1:(j-1)])
          
          #calculate IA (how many times was animal seen after j)
          ia <- sum(ic[i, (j+1):ns])
          
          #counts for R component (resighting)
          if (ib > 0 && ia > 0) {
            n11 <- n11 + 1
          } else if (ib > 0 && ia == 0) {
            n12 <- n12 + 1
          } else if (ib == 0 && ia > 0) {
            n21 <- n21 + 1
          } else {
            n22 <- n22 + 1
          }
          
          #counts for S component (subsequent components?) (only if j < ns-1)
          if (j < ns-1) {
            if (ib > 0 && ic[i, j+1] >= 1) {
              m11 <- m11 + 1
            } else if (ib > 0 && ia > 0) {
              m12 <- m12 + 1
            } else if (ib == 0 && ic[i, j+1] >= 1) {
              m21 <- m21 + 1
            } else if (ib == 0 && ia > 0) {
              m22 <- m22 + 1
            }
          }
        }
      }
      
      #chi-squared for resight component
      iuser <- 1
      r1 <- n11 + n12
      r2 <- n21 + n22
      c1 <- n11 + n21
      c2 <- n12 + n22
      
      if (r1 < chat_rot || r2 < chat_rot || c1 < chat_rot || c2 < chat_rot) {
        iuser <- 0
      }
      
      tot <- r1 + r2
      if (r1 == 0 || r2 == 0 || c1 == 0 || c2 == 0) {
        compr <- 0.0
        idfr <- 0
      } else {
        e11 <- r1 * c1 / tot
        e12 <- r1 * c2 / tot
        e21 <- r2 * c1 / tot
        e22 <- r2 * c2 / tot
        compr <- (n11 - e11)^2 / e11 + (n12 - e12)^2 / e12 + 
          (n21 - e21)^2 / e21 + (n22 - e22)^2 / e22
        idfr <- 1
      }
      
      # Chi-squared for S component
      iuses <- 1
      r1 <- m11 + m12
      r2 <- m21 + m22
      c1 <- m11 + m21
      c2 <- m12 + m22
      
      if (r1 < chat_rot || r2 < chat_rot || c1 < chat_rot || c2 < chat_rot) {
        iuses <- 0
      }
      
      tot <- r1 + r2
      if (r1 == 0 || r2 == 0 || c1 == 0 || c2 == 0) {
        comps <- 0.0
        idfs <- 0
      } else {
        e11 <- r1 * c1 / tot
        e12 <- r1 * c2 / tot
        e21 <- r2 * c1 / tot
        e22 <- r2 * c2 / tot
        comps <- (m11 - e11)^2 / e11 + (m12 - e12)^2 / e12 + 
          (m21 - e21)^2 / e21 + (m22 - e22)^2 / e22
        idfs <- 1
      }
      
      chigt <- chigt + compr * iuser + comps * iuses
      idfgt <- idfgt + idfr * iuser + idfs * iuses
    }
  }
  
  #calculate variance inflation factor (remember vif=c_hat)
  if (idfgt > 0) {
    vif <- chigt / idfgt
    vif <- max(vif, 1.0)
  } else {
    vif <- 1.0
  }
  
  return(list(vif = vif, chigt = chigt, idfgt = idfgt))
}

#helper function for TEST2
cjs_test2 <- function(ns, m, chat_rot) {
  # Initialize outputs
  iuse <- integer(ns)
  idf <- integer(ns)
  chisq <- numeric(ns)
  chitot <- 0.0
  idftot <- 0
  
  #exit test2 if not possible
  if (ns < 4) {
    return(list(chisq = chisq, idf = idf, chitot = chitot, idftot = idftot))
  }
  
  #compute components one by one
  for (l in 2:(ns-2)) {
    iuse[l] <- 1
    
    #find frequencies for contingency table
    n <- matrix(0, nrow = 2, ncol = ns)
    for (j in (l+1):ns) {
      n[1, j] <- sum(m[1:(l-1), j])
      n[2, j] <- m[l, j]
    }
    
    #find row and column totals
    r <- numeric(2)
    c <- numeric(ns)
    
    for (i in 1:2) {
      r[i] <- sum(n[i, (l+1):ns])
    }
    
    for (j in (l+1):ns) {
      c[j] <- sum(n[, j])
    }
    
    tot <- sum(r)
    
    #check minimum counts
    if (r[1] < chat_rot || r[2] < chat_rot) {
      iuse[l] <- 0
    }
    
    for (j in (l+1):ns) {
      if (c[j] < chat_rot) {
        iuse[l] <- 0
      }
    }
    
    #calculate chi-squared
    if (r[1] <= 0 || r[2] <= 0) {
      chisq[l] <- 0.0
      idf[l] <- 0
    } else {
      chisq[l] <- 0.0
      idf[l] <- ns - l - 1
      
      for (j in (l+1):ns) {
        if (c[j] <= 0) {
          idf[l] <- idf[l] - 1
        } else {
          for (i in 1:2) {
            exp <- r[i] * c[j] / tot
            chisq[l] <- chisq[l] + (n[i, j] - exp)^2 / exp
          }
        }
      }
      
      if (idf[l] <= 0) {
        idf[l] <- 0
        chisq[l] <- 0.0
        iuse[l] <- 0
      }
    }
    
    chitot <- chitot + chisq[l] * iuse[l]
    idftot <- idftot + idf[l] * iuse[l]
  }
  
  return(list(chisq = chisq, idf = idf, chitot = chitot, idftot = idftot))
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
                                   ints)
    
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
                               ints=rep(1,ncol(prepped_data$ch)))
    }
    endtime<-Sys.time()
    optim_speed<-endtime-starttime 
    
    est_escapement<-total_escapement(prepped_data$ch,
                                     optim_results$par,
                                     model_fit$cap_X[[i]],
                                     model_fit$surv_X[[i]],
                                     ints=rep(1,ncol(prepped_data$ch)))
    
    ######################################################
    #calculate model fit statistics (aic, qaic, qaicc)
    ######################################################
    loglik<-optim_results$value
    ic<-prepped_data$ch
    nan=nrow(ic)
    ns=ncol(ic)
    
    fit_results<-cjs_fit(nan,ns,ic,ng=1,ig=rep(1, nan))
    df<-fit_results$idfgt
    c_hat<-fit_results$vif
    #AIC needs to be calculated for each model separatly
      #but QAIC needs to be done outside of this after all AIC is calc
      #reference Amstrup et al Handbook
    AIC=2*loglik+2*df
    
    #AICC is corrected for small samples (<40 per parameter)
    AICC=AIC+((2*df)*(df+1))/(nan-df-1)
    
    #QAIC
      #c_hat=variance inflation factor (vif)
      #c_hat measures the overdispersion, or excess variance, of the data
      #due to model assumption violations
      #this changes the analysis from a true "maximum likelihood"
      #to a "quasi-likelihood" thus "quasi-AIC" or QAIC
      #we do this with tow methods called "Test 2" and "Test 3"
      #page 236 of Amstrup et al Handbook

    QAIC=(2*loglik)/c_hat+(2*df)
    QAICC=QAIC+(2*df*(df+1))/(nan-df-1)
    
    if(n_boot!=0){
      boot_start<-Sys.time()
      boot_results<-CJS_bootstrap(n_boot,
                                  ic,
                                  model_fit$cap_X[[i]],
                                  model_fit$surv_X[[i]])
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
        scale_x_continuous(breaks = seq(0,10000,500)) +
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
  return(model_results)
  write.csv(model_results,paste("outputs/CJS_outputs-",
                                Sys.Date(),
                                ".csv",
                                sep=""),
            row.names = F)
}

##########################################
#model selection for shiny ap
##########################################
#wrapper function for 
CJS_model_select_app <- function(covars_used, sex_matrix, lengths_matrix, ch, selected_models) {
  models = c(
    "constant capture and survival rates",
    "constant capture rate and survival related the sex",
    "constant capture rate and survival related the length",
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

