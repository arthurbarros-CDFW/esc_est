#testing analysis
library(tidyverse)


#what do we want for each dataset?
#1) number of capture history records
#2) number of survey periods
#3) number of chops
#4) optim duration for each model and approach (mra/optim)
#5) estimation of escapement for each model and approach
# ) coefficient estimates for each model and approach
#6) upper and lower CI for escapement estimate for each model and approach
#7) AIC? model selection criteria
#8) number of failures to converge
#9) total bootstrapping time

optim_files<-c("Roptim_AR 2021",
               "Roptim_Butte 2020",
               "Roptim_Butte 2022",
               "Roptim_FR 2021",
               "Roptim_SRfall 2021",
               "Roptim_SRwinter 2021")

mra_files<-c("mra_AR 2021",
             "mra_Butte 2020",
             "mra_Butte 2022",
             "mra_FR 2021",
             "mra_SRfall 2021",
             "mra_SRwinter 2021")

#load mra files
mra_results<-optim_results<-data.frame()
for(i in 1:length(mra_files)){
  optim_name<-optim_files[i]
  mra_name<-mra_files[i]
  
  optim_d<-readRDS(file = paste("outputs/testing/",
                                      optim_name,".rds",sep=""))
  mra_d<-readRDS(file = paste("outputs/testing/",
                                    mra_name,".rds",sep=""))
  optim_d<-optim_d%>%
    mutate(data_set=optim_name)%>%
    rename(loglik_fails=optim_fails)
  mra_d<-mra_d%>%
    mutate(data_set=mra_name,
           method="mra")
  optim_results<-optim_results%>%rbind(optim_d)
  mra_results<-mra_results%>%rbind(mra_d)
}

#summarize mra results
#confidence intervals
conf.level = 95
alpha = 1 - conf.level/100
lower = alpha/2
upper = 1 - alpha/2

mra_summary<-mra_results%>%
  group_by(data_set,method,model,model_num)%>%
  summarise(N_est=ceiling(mean(escapement)),
            lower_ci=ceiling(quantile(escapement,probs = c(lower),na.rm=T)),
            upper_ci=ceiling(quantile(escapement,probs = c(upper),na.rm=T)),
            loglik_fails=sum(is.nan(log_likelihood)),
            n_boot=n(),
            boot_time=sum(time))

combined_results<-optim_results%>%
  select(data_set,method,model,model_num,N_est,boot_time,
         lower_ci,upper_ci,loglik_fails,n_boot)%>%
  rbind(mra_summary)

combined_results<-combined_results%>%
  mutate(data_set=sub(".*_", "", data_set))

model_list<-unique(combined_results$model)
data_list<-unique(combined_results$data_set)
#make some plots
for(i in 1:length(data_list)){
  d<-combined_results%>%
    filter(data_set==data_list[i])
  
  dodge <- position_dodge(width=.5) 
  p<-ggplot(data=d,aes(x = factor(model_num),
                       y = N_est,
                       color = method,
                       group = method))+
    geom_point(aes(x=factor(model_num),y=N_est),
               position = dodge)+
    geom_errorbar(aes(x=factor(model_num),
                      ymax=upper_ci,
                      ymin=lower_ci),
                  position = dodge)+
    geom_text(aes(x=factor(model_num),
                  y=upper_ci,
                  label=loglik_fails),
              position = dodge,
              vjust=-.5,
              show.legend = FALSE)+
    theme_bw()+
    labs(x="model number",
         y="estimated escapement",
         title=data_list[i],
         subtitle = paste("boot strap iterations:",
                          unique(d$n_boot),sep=" "),
         color="optimization method")
  ggsave(paste("outputs/testing_plots/",
               data_list[i],".png",sep=""),
         scale=2)
}
