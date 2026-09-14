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
      #remove rows with NA length
      covars <- covars%>%filter(!disctag %in% missing_lengths$disctag)
      
      #also remove corresponding rows from ch
      ch <- ch%>%filter(!disctag %in% missing_lengths$disctag)
    }
    if(any(is.na(covars$sex))){
      missing_sex<-covars%>%filter(is.null(sex)|is.na(sex))
      message(paste("Found", sum(is.na(covars$sex)), 
                    "NA values in sex. These rows will be removed."))
      #remove rows with NA length
      covars <- covars%>%filter(!disctag %in% missing_sex$disctag)
      
      #also remove corresponding rows from ch
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