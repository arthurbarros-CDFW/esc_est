#testing analysis
library(tidyverse)

mra_test<-readRDS("outputs/testing/mra_LAR 2023_model_9.rds")
optim_test<-readRDS("outputs/testing/Roptim_FR 2021.rds")

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
               "Roptim_Butte 2021",
               "Roptim_FR 2021",
               "Roptim_SRfall 2021",
               "Roptim_SRwinter 2021")

mra_files<-c("mra_AR 2021",
             "mra_Butte 2020",
             "mra_Butte 2021",
             "mra_FR 2021",
             "mra_SRfall 2021",
             "mra_SRwinter 2021")

#load mra files
mra_results<-data.frame()
for(i in 1:length(mra_files)){
  
}