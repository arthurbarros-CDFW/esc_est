#quick esc_est testing
#created: 12/23/2025
#last update: 12/23/2025
rm( list = ls()) #clear env
library(tidyverse)
library(escapeMR)
library(Rcpp)

source("scripts/CVCS_functions.R")

results<-CJS_run()

#CJSscript()
