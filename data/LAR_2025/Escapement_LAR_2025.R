#This script was prepared for Jennifer O'Brien by Vanessa Kollmar
#If you have any questions please reach out to Vanessa Kollmar (CDFW R4) via TEAMS, email or phone

# Installation of Packages ------------------------------------------------
##Installation only needs to happen the first time you use this program
##After you intall these packages put a pound sign (hashtag) in front of "install" to deactivate this line of code.
#install.packages(c("shiny", "mra", "DT", "withr", "RODBC", "devtools", "readxl", "tidyverse", "lubridate", "dplyr", "janitor", "ggplot2", "cowplot"))


# Library and Work Directory ----------------------------------------------
##The library lines of code will turn on the packages you need to run for this analysis.
##You will need to run this code every time you use the program
library(readxl)
library(tidyverse)
library(lubridate)
library(dplyr)
library(janitor)
library(stringr)
##Setting your working directory is just telling the program in which folder you have stored your escapement data
##and where you would like to export your files.make sure to use all forward slashes "/". IF you copy and paste
##your path you will have to replace all backslashes "\" with forward slashes.

setwd("U:/241-FISHERIES/R2 Low Elevation Fisheries/American River & NFH/LAR Carcass Survey/2025/Reports/Annual/Estimate")


# Importing and Formatting Your Data --------------------------------------
##The word that the arrow points is the name you are giving your data set. Remember that R does not like spaces!
##Don't use spaces when naming data sets, column names or row names. Also make sure the file you import isn't already
##named with spaces. Erase spaces or use an underscore "_" instead. 
##Important! Before your import your data, open the excel file and replace the word "Not chopped" with "Notchopped" in the 
##recovery sheet. If you don't make this change the recovery analysis wont work.
Tags<-read_excel("EstimateQueries2026.01.05.xls", sheet= "QQrySurveyIndividualFish")
Rec<-read_excel("EstimateQueries2026.01.05.xls", sheet= "QQrySurveyRecoveredTags") 
Chops<-read_excel("EstimateQueries2026.01.05.xls", sheet= "QQrySurveyChops")

##This section will format your data to produce your covariate file. "Cov23" is the name I chose since
##this is covariate data for 2023. In subsequent years you can change the 23 to 24, 25, etc to represent the 
##escapement year you are working on. All the "23" have to be changed where you see then below. 
##I assigned dummy tags to the fish you sex and measured but never tagged. This is necessary if you want to include then
##in the analysis. The dummy tags can be differentiated by the letter "A" at the end. 
Cov25<-Tags%>%
  mutate(Month= month(Date), SurveyYear = if_else(month(Date) %in% c(1:2), year(Date) - 1, year(Date)))%>% #this line is just establishing your survey year with an extra column at the end of the file
  transmute(disctag= as.numeric(DiscTagApplied), Week, Sex, Length= FLmm) #this line is renaming your columns so that they match the column names escapeMR needs for the analysis. It is all only keeping the 3 columns for the covariate run

# Find number of NAs in the disctag column that needs to be replaced with a dummy value
LengthNAs <- length(which(is.na(Cov25$disctag)))
# Beginning at the max value used in the dsctag column, count by one to the number of NA values that needs dummy values
DummyIndex <- (max(Cov25$disctag, na.rm = T) + 1):(max(Cov25$disctag, na.rm = T) + LengthNAs)
# Fill in those NAs with the dummy values, these should be exactly the same length
Cov25$disctag[is.na(Cov25$disctag)] <- DummyIndex

##This next line of code will export your covariate file.   
write_excel_csv(select(Cov25, -Week), path="Cov25LAR.csv")

##This section will format your data to produce the chops table. 

week<-""  #Will use later to add the Week header back into the data frame

Chops25<-Chops%>%
  mutate(Month= month(Date), SurveyYear = if_else(month(Date) %in% c(1:2), year(Date) - 1, year(Date))) %>% #the only reason I have this line of code here is if you ever had a multiyear database and wanted to filter it by survey year, you could. 
  select(Week, ChopCount)%>%
  na.omit()%>%
  group_by(Week)%>%
  mutate(Sum_of_Chops_per_Wk= sum(ChopCount))%>%
  select(Week, Sum_of_Chops_per_Wk)%>%
  unique()%>%
  arrange(Week)%>%
  pivot_wider(names_from = Week, values_from = Sum_of_Chops_per_Wk )

Chops25<- cbind(week,Chops25)

##This next line of code will export your chop file.
write_excel_csv(Chops25, path="Chop25LAR.csv")

##This section will format your data to produce your capture history
##This line of code will show you if a single tag code was used more than once during the season. You will need to 
##correct this error in the access database before you move forward. If this code returns nothing, then you have no 
##repeat tag code errors.
TagQC<-Tags%>%
  transmute(disctag= DiscTagApplied, Week)%>%
  na.omit()%>%
  filter(duplicated(disctag))
  
##This line of code will assign a number 1 to the tag number on the week it was tagged. Remember that 1's in the
##the analysis represent a capture event. 
T25<-Tags%>% 
  mutate(Month= month(Date), SurveyYear = if_else(month(Date) %in% c(1:2), year(Date) - 1, year(Date))) %>% #Again this line of code can be ignored.
  transmute(disctag= as.numeric(DiscTagApplied), Week)%>%
  mutate(Chop= as.numeric("1"))%>%
  na.omit()

##These few lines of code will name the chops, you took covariate data for, with the same dummy tag numbering system used before. 
##Because these fish were chopped and not used in the mark-recapture part of the study they will be assigned a "2" which represents a chop in the analysis

# Find all instances of disctag that has a dummy index and assign a Chop = 2 to them
D25 <- Cov25 %>%
  filter(disctag %in% DummyIndex) %>%
  transmute(disctag, Week, Chop = 2)

##This line of code will assign the appropriate numbers (1 or 2) to recapture and chop events 
R25<-Rec%>%  
  mutate(Month= month(Date), SurveyYear = if_else(month(Date) %in% c(1:2), year(Date) - 1, year(Date))) %>% #Again you can ignore this
  transmute(Week, disctag= as.numeric(TagRecovered), Chop= Disposition)%>%
  mutate(Chop= ifelse(Chop == "Notchopped", 1, 2 ))

##This line of code will tell you which tag codes were recaptured more than once in the same week. 
##This has to be corrected in your database or else you cannot proceed.
CapQC<-bind_rows(T25,D25,R25)%>%     
  group_by(Week)%>%
  filter(duplicated(disctag))%>%
  select(disctag, Week)

# This will give you a warning of non-unique "Chop" values if you have duplicated disctag values in your Access DB
Cap25<-bind_rows(T25,D25,R25)%>%     
  mutate(Week= factor(Week, levels=1:10))%>%  #the levels will need to be changed based on the number of survey weeks 
  complete(Week,disctag, fill = list(Chop=0))%>%
  pivot_wider(names_from = Week, values_from = Chop )

#This file will not export if you have unresolved database errors from the TagQC and CapQC analysis
write_excel_csv(Cap25, path="CapHis25LAR.csv")

######DONE WITH THIS CODE IF YOU ARE USING THE SHINY APP (MUCH EASIER)
######ONLY USE THIS NEXT SECTION IF YOU CAN'T USE THE SHINY APP
# escapeMR ----------------------------------------------------------------

library(escapeMR)
escapeMR("script")
