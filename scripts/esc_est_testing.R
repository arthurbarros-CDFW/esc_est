#quick esc_est testing
#created: 12/23/2025
#last update: 09/03/2026
rm( list = ls()) #clear env
library(tidyverse)
library(escapeMR)
library(Rcpp)

#Load functions
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)
sourceCpp('scripts/CJS_functions.cpp')

results<-CJS_run()

#CJSscript()
ch<-read.csv('data/demo_CH.csv')
