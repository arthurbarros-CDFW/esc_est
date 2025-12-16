library(tidyverse)
library(escapeMR)
library(Rcpp)
sourceCpp('scripts/CJS_functions.cpp')
source("scripts/CVCS_functions.R")

#Methods comparisons

#automated data loading for testing
path<-"data/testing data/"
files <- list.files(path, full.names = TRUE)
test_data<-data.frame()

for(i in 1: length(files)){
  file_name<-basename(files[i])
  river<-gsub("_.*","",file_name)
  year<-gsub(".*_(\\d{4})_.*", "\\1",file_name)
  data_type<-gsub(".*?_.*?_(.*?)\\.csv", "\\1",file_name)
  d<-data.frame(file_name,river,year,data_type)
  test_data<-test_data%>%rbind(d)
}

test_data$id<-paste(test_data$river,test_data$year)

prepped_list<-list()

id_list<-unique(test_data$id)
#the end goal here is a list that has each year and river combo
#with an output of ch (including chops), lengths_matrix, and sex_matrix
for(i in 1:length(id_list)){
  target_id<-id_list[i]
  target_data<-test_data%>%
    filter(id==target_id)
  
  target_ch<-target_data%>%filter(data_type=="CH")%>%select(file_name)
  target_chops<-target_data%>%filter(data_type=="Chops")%>%select(file_name)
  target_covars<-target_data%>%filter(data_type=="Covars")%>%select(file_name)
  
  ch<-read.csv(file = paste("data/testing data/",target_ch,sep=""), head=T, as.is=T)
  chops<-read.csv(file = paste("data/testing data/",target_chops,sep=""), head=T, as.is=T)
  covars<-read.csv(file = paste("data/testing data/",target_covars,sep=""), head=T, as.is=T)
  
  #Part 2: prepare ch and covars data
  ch=ch[-1] #remove disctag vector from capture histories
  covars$sex<-as.numeric(ifelse(covars$sex%in%c('F',"f"),1,0)) #change sex to numeric value

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
  
  n<-nrow(clean_chops)
  chops_covars<-as.matrix(covars[sample(nrow(covars),n,replace=T),]) #sample for number of chops
  covars<-covars%>%rbind(chops_covars)
  
  #add in chops to ch
  ch<-ch%>%rbind(clean_chops)
  
  lengths_matrix<-matrix(covars$length,
                         nrow=nrow(ch),ncol=ncol(ch))
  sex_matrix<-matrix(as.numeric(covars$sex),
                     nrow=nrow(ch),ncol=ncol(ch))
  
  prepped_list[[target_id]]<-list(
    ch=ch,
    lengths_matrix=lengths_matrix,
    sex_matrix=sex_matrix
  )
  
}

saveRDS(prepped_list,"outputs/testing/prepped_data.rds")

#Run a loop that gets the following for both R optim methods and mra methods:

#model structure
#N_est
#loglik
#AIC
#QAICC
#cap_beta1 and 2
#surv_beta1 and 2
#optim_speed
#n_boot
#total_boot_time
#lower_ci
#upper_ci
#number of optimization fails

n_boot=500
models = c(	"constant capture and survival rates",
            "constant capture rate and survival related the sex",
            "constant capture rate and survival related the length",
            "capture related to sex and constant survival rate",
            "capture related to length and constant survival rate",
            "capture related to sex and survival related to length",
            "capture related to length and survival related to sex",
            "capture related to sex and survival related to sex",
            "capture related to length and survival related to length")

Roptim_results=data.frame()
mra_total_results<-data.frame()
boot_list=list()

for(r in 4){
  data_id<-id_list[r]
  data<-prepped_list[[data_id]]
  ch=data$ch
  lengths_matrix=data$lengths_matrix
  sex_matrix=data$sex_matrix
  
  #structure data for all model forms
  model_covars = list(cap_X=NULL, surv_X=NULL)
  data_id_results=data.frame()
  for(m in 1:length(models)){
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
  
  c_hat = 1
  initial_beta=numeric(4)
  model_results<-data.frame()
  mra_results<-data.frame()
  for(i in 1:length(models)){
    print(paste(data_id," beginning model:",i," (",
                models[i],")",sep=""))
    ######################
    ###R optim
    starttime<-Sys.time()
    {gc()
      optim_results<-cpp_optim(beta=initial_beta,
                               ch=as.matrix(ch),
                               cap_X=as.matrix(model_covars$cap_X[[i]]),
                               surv_X=as.matrix(model_covars$surv_X[[i]]),
                               ints=rep(1,ncol(ch)))
    }
    endtime<-Sys.time()
    optim_speed<-endtime-starttime 
    est_escapement<-total_escapement(ch,
                                     optim_results$par,
                                     model_covars$cap_X[[i]],
                                     model_covars$surv_X[[i]],
                                     ints=rep(1,ncol(ch)))
    loglik<-optim_results$value
    ic<-ch
    nan=nrow(ic)
    ns=ncol(ic)
    
    fit_results<-cjs_fit(nan,ns,ic,ng=1,ig=rep(1, nan))
    df<-fit_results$idfgt
    c_hat<-fit_results$vif
    AIC=2*loglik+2*df
    AICC=AIC+((2*df)*(df+1))/(nan-df-1)
    
    QAIC=(2*loglik)/c_hat+(2*df)
    QAICC=QAIC+(2*df*(df+1))/(nan-df-1)
    
    boot_time<-0
    if(n_boot!=0){
      boot_start<-Sys.time()
      boot_results<-CJS_bootstrap(n_boot,
                                  ch,
                                  model_covars$cap_X[[i]],
                                  model_covars$surv_X[[i]],
                                  ints=rep(1,ncol(ch)))
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
      boot_time=boot_time+boot_speed
    }
    total_boot_time=
    d<-data.frame(
      data_set=data_id,
      method="R optim",
      model_num=i,
      model=models[i],
      N_est=est_escapement,
      loglik=loglik,
      optim_speed=optim_speed,
      AIC=AIC,
      QAICC=QAICC,
      cap_beta1=optim_results$par[1],
      cap_beta2=optim_results$par[2],
      surv_beta1=optim_results$par[3],
      surv_beta2=optim_results$par[4],
      n_boot=n_boot,
      boot_time=boot_time,
      lower_ci=ci$lower_ci,
      mid_ci=ci$mid_ci,
      upper_ci=ci$upper_ci,
      optim_fails=sum(is.na(boot_results$escapement)==TRUE)
    )
    model_results=model_results%>%rbind(d)
    data_id_results=data_id_results%>%rbind(d)
    
    saveRDS(model_results,paste("outputs/testing/Roptim_",data_id,"_model_",i,'.rds',sep=""))
    
    ######################
    ###mra
    iteration_sample<-list()
    #Run the for loop for each iteration
    for(r in 1:length(n_boot)){
      iter_start<-Sys.time()
      
      #Part 2.1: index and sample the capture histories
      index=1:dim(ch)[1]
      samp<- sample(index, replace = T)
      iteration_sample[[r]]<-samp
      #samp<-iteration_sample[[373]]
      
      #Part 2.2: create new capture histories and
      #cap_X and surv_X matrices based on sampled indices
      ch_iteration = ch[samp,]
      cap_X_iteration=model_covars$cap_X[[i]][samp,]
      surv_X_iteration=model_covars$surv_X[[i]][samp,]
      
      #Part 2.3: use cpp_optim to find the optimizal beta
      #parameters for the given iteration data
      source("scripts/escapeMR_testing.R")
      iter_beta<-ans$parameters
      iter_lik<-ans$loglik
      iter_lik
      
      #Part 2.4: estimate total escapement for iteration
      est_esc_iter<-total_escapement(ch_iteration,
                                     iter_beta,
                                     cap_X_iteration,
                                     surv_X_iteration,
                                     ints=rep(1,ncol(ch)))
      
      iter_end<-Sys.time()
      iter_time<-iter_end-iter_start
      
      #Part 2.5: store relevant iteration information in results
      d<-data.frame(data_set=id_list[r],
                    method="R optim",
                    model=models[i],
                    model_num=i,
                    "iteration"=r,
                    "log_likelihood"=iter_lik,
                    "beta1"=ans$parameters[1],
                    "beta2"=ans$parameters[1],
                    "beta3"=ans$parameters[1],
                    "beta4"=ans$parameters[1],
                    "escapement"=est_esc_iter,
                    "time"=iter_time)
      mra_results<-mra_results%>%rbind(d)
      print(r)
      saveRDS(mra_results,paste("outputs/testing/mra_",data_id,"_model_",i,'.rds',sep=""))
    }
    saveRDS(model_results,paste("outputs/testing/Roptim_",data_id,'.rds',sep=""))
    saveRDS(mra_results,paste("outputs/testing/mra_",data_id,'.rds',sep=""))
  }
  Roptim_results<-Roptim_results%>%
    rbind(model_results)
  mra_total_results<-mra_total_results%>%
    rbind(mra_results)
  
}

saveRDS(Roptim_results,"outputs/testing/Roptim_results.rds")
saveRDS(mra_total_results,"outputs/testing/mra_results.rds")
