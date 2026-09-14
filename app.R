# app.R
# CDFW Escapement Estimator
# for carcass mark-recapture surveys

library(shiny)
library(tidyverse)
library(ggplot2)
library(shinythemes)
library(Rcpp)
library(DT)

#load functions
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)
sourceCpp('scripts/CJS_functions.cpp')

#define the ui
ui <- fluidPage(
  theme = shinytheme("cosmo"),
  # title
  titlePanel("CDFW CVCS Escapement Model"),
  
  sidebarLayout(
    sidebarPanel(
      
      #################
      #file inputs
      #################
      #cap history input
      fileInput("ch_input", "Upload Capture History (CSV)*",
                accept = c(".csv")),
      
      #chops input
      checkboxInput("use_chop", "Include chop file", value = FALSE),
      conditionalPanel(
        condition = "input.use_chop == true",
        fileInput("chop_input", "Upload Chop File (CSV)",
                  accept = c(".csv"))
      ),
      
      #covariate input
      checkboxInput("use_cov", "Include covariates file", value = FALSE),
      conditionalPanel(
        condition = "input.use_cov == true",
        fileInput("cov_input", "Upload Covariates File (CSV)",
                  accept = c(".csv"))
      ),
      
      #intervals input
      checkboxInput("use_ints", "Include unequal intervals", value = FALSE),
      conditionalPanel(
        condition = "input.use_ints == true",
        fileInput("ints_input", "Upload Unequal Intervals File (CSV)",
                  accept = c(".csv"))
      ),
      
      #subsampling input
      checkboxInput("use_subsampling", "Was any subsampling conducted?", value = FALSE),
      conditionalPanel(
        condition = "input.use_subsampling == true",
        fileInput("subsampling_input", "Upload Subsampling File (CSV)",
                  accept = c(".csv"))
      ),
      
      #help text
      helpText("* Required field"),
      helpText("Note: All files should be in CSV format with headers."),
      
      #################
      #model selection
      #################
      conditionalPanel(
        condition = "input.use_cov == true",
        checkboxGroupInput("use_model", "Select one or more model",
                           c(	"constant capture and survival rates",
                              "constant capture rate and survival related the sex",
                              "constant capture rate and survival related the length",
                              "capture related to sex and constant survival rate",
                              "capture related to length and constant survival rate",
                              "capture related to sex and survival related to length",
                              "capture related to length and survival related to sex",
                              "capture related to sex and survival related to sex",
                              "capture related to length and survival related to length"))
      ),
      
      #button to run analysis
      actionButton("run_models", "run CJS models", 
                   class = "btn-primary")
      
      
    ),
    
    #################
    #main panel
    #################
    mainPanel(

      #tabs for outputs
      tabsetPanel(
        tabPanel("Model Results",
                 uiOutput("model_results"),
                 #enter number of bootstrap reps
                 h4(strong("Boostrap for confidence intervals?")),
                 checkboxInput("use_boots", "Y/N", value = FALSE),
                 conditionalPanel(
                   condition = "input.use_boots == true",
                   numericInput("boot_input", "Enter a number of bootstrap replications to perform",
                                10,min=10,max=1000),
                   helpText("Enter numeric value between 10 - 1000")
                 ),
                 actionButton("run_selected_models",
                              "Estimate Escapement",
                              class = "btn-primary")
                 ),
        tabPanel("Escapement Plot", 
                 plotOutput("p_esc",height = "400px"),
                 downloadButton("download_esc_plot", "Download Plot")
                 ),
        tabPanel("Data Preview", 
                 dataTableOutput("ch_preview"),
                 conditionalPanel(
                   condition = "input.use_chop == true",
                   dataTableOutput("chop_preview")
                 ),
                 conditionalPanel(
                   condition = "input.use_cov == true",
                   dataTableOutput("cov_preview")
                 ))
      )
    )
    
    #################
  )
)

#################
#server side
#################
server <- function(input, output, session) {
  
  #################
  #reactive value to store uploaded data
  #################
  uploaded_data <- reactiveValues(
    ch = NULL,
    chop = NULL,
    cov = NULL,
    ints = NULL
  )
  
  #observe file uploads
  observe({
    req(input$ch_input)
    uploaded_data$ch <- read.csv(input$ch_input$datapath)
  })
  
  observe({
    if(input$use_chop) {
      req(input$chop_input)
      uploaded_data$chop <- read.csv(input$chop_input$datapath)
    } else {
      uploaded_data$chop <- NULL
    }
  })
  
  observe({
    if(input$use_ints) {
      req(input$ints_input)
      uploaded_data$ints_input <- read.csv(input$ints_input$datapath)
    } else {
      uploaded_data$ints_input <- NULL
    }
  })
  
  observe({
    if(input$use_subsampling) {
      req(input$subsampling_input)
      uploaded_data$subsampling_input <- read.csv(input$subsampling_input$datapath)
    } else {
      uploaded_data$subsampling_input <- NULL
    }
  })
  
  observe({
    if(input$use_cov) {
      req(input$cov_input)
      uploaded_data$cov <- read.csv(input$cov_input$datapath)
    } else {
      uploaded_data$cov <- NULL
    }
  })
  
  #########################
  #data preparation
  #########################
  prepare_data<-reactive({
    req(uploaded_data$ch)
    
    #initial variables
    ch=uploaded_data$ch
    chops=uploaded_data$chop
    covars=uploaded_data$cov
    
    #handle intervals
    if(input$use_ints && !is.null(uploaded_data$ints_input)) {
      ints <- uploaded_data$ints_input
      
      #ensure intervals have correct lengths
      if(length(ints) != (ncol(ch)-1)) {
        showNotification("Number of intervals must be one less than number of occasions, defaulting to equal intervals", 
                         type = "error")
        ints <- rep(1, ncol(ch)-1)
      }
    } else {
      ints <- rep(1, ncol(ch)-1)  # Default to equal intervals
    }
    
    #handle subsampling
    if(input$use_subsampling && !is.null(uploaded_data$subsampling_input)){
      sub.sampling<-uploaded_data$subsampling_input
      sub.sampling <- sub.sampling[,grep("[0-9]",names(sub.sampling))]
      sub.sampling <- unlist(sub.sampling)  # must be a vector, not a data frame
      # at this point ch may have other columns in it
      if( length(sub.sampling) != (ncol(ch)-1) ){
        stop("Number of records must equal number of sampling periods.")
      }
    } else {
      sub.sampling = rep(1,ncol(ch)-1)
    }
    
    #prepare ch and covar data
    ch=ch[-1]
    null_matrix<-matrix(1,nrow=nrow(ch),ncol=ncol(ch))
    
    #process covariates
    if(input$use_cov && !is.null(covars)){
      covars$sex <- as.numeric(ifelse(covars$sex %in% c('F', "f"), 1, 0))
    }
    
    #prepare tagged and observed_per_week for total_escapement()
    tagged_per_week<-colSums(ch>=1,na.rm=T)
    observed_per_week <- tagged_per_week*sub.sampling
    
    #prepare chops data
    if(input$use_chop && !is.null(chops)){
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
    
    #generate chops covariate data
    #here we'll do things differently than the escapeMR app
    #instead of taking the mean values for sex and length and assigning them to chops
    #we will randomly sample existing covars and assign to chop data
    #I like this better because it incorporates variability into the chops
    if(input$use_cov && input$use_chop && !is.null(covars) && !is.null(chops)){
      n <- nrow(clean_chops)
      chops_covars <- as.matrix(covars[sample(nrow(covars), n, replace = TRUE),])
      covars <- covars %>% rbind(chops_covars)
    }
    
    if(input$use_chop && !is.null(chops)){
      # add in chops to ch
      ch <- ch %>% rbind(clean_chops)
    }
    
    if(input$use_cov && !is.null(covars)){
      lengths_matrix <- matrix(covars$length,
                               nrow = nrow(ch), ncol = ncol(ch))
      sex_matrix <- matrix(as.numeric(covars$sex),
                           nrow = nrow(ch), ncol = ncol(ch))
    } else {
      lengths_matrix <- null_matrix
      sex_matrix <- null_matrix
    }
    
    #prepare model covariates based on selected models
    model_data <- CJS_model_select_app(
      covars_used = input$use_cov,
      sex_matrix = sex_matrix,
      lengths_matrix = lengths_matrix,
      ch = ch,
      selected_models = input$use_model
    )
    
    #return formatted data
    list(
      'ch' = ch,
      "covars"=covars,
      'lengths_matrix' = lengths_matrix,
      'sex_matrix' = sex_matrix,
      'intervals' = ints,
      'covars_used' = input$use_cov,
      'chops_used' = input$use_chop,
      'model_data' = model_data,
      'tagged_per_week'=tagged_per_week,
      'observed_per_week'=observed_per_week,
      'subsampling_weeks'=sub.sampling
    )
  })
  
  #########################
  #data previews
  #########################
  output$ch_preview <- renderDT({
    (uploaded_data$ch)
  })
  
  output$chop_preview <- renderDT({
    if(input$use_chop && !is.null(uploaded_data$chop)) {
      (uploaded_data$chop)
    }
  })
  
  output$cov_preview <- renderDT({
    if(input$use_cov && !is.null(uploaded_data$cov)) {
      (uploaded_data$cov)
    }
  })
  
  #########################
  #run CJS model when button pressed
  #########################
  observeEvent(input$run_models,{
    req(prepare_data())
    prepped_data <- prepare_data()
    
    model_list<-c(	"constant capture and survival rates",
                   "constant capture rate and survival related the sex",
                   "constant capture rate and survival related the length",
                   "capture related to sex and constant survival rate",
                   "capture related to length and constant survival rate",
                   "capture related to sex and survival related to length",
                   "capture related to length and survival related to sex",
                   "capture related to sex and survival related to sex",
                   "capture related to length and survival related to length")
    
    #create a progress object
    progress <- Progress$new(session, min=0, max=1)
    on.exit(progress$close()) #progress bar closes when done
    
    #ensure intervals are properly formatted
    intervals <- if(!is.null(prepped_data$intervals)) {
      as.numeric(prepped_data$intervals)
    } else {
      rep(1, ncol(prepped_data$ch)-1)
    }
    
    #access model data
    model_data <- prepped_data$model_data
    
    #pull subsampling data
    subsampling_weeks<-prepped_data$subsampling_weeks
    observed_per_week<-prepped_data$observed_per_week
    tagged_per_week<-prepped_data$tagged_per_week
    
    #display 
    output$model_output <- renderPrint({
      cat("Selected models:", paste(model_data$models_ran, collapse = ", "), "\n")
      cat("Number of capture histories:", nrow(model_data$ch), "\n")
    })
    
    #set initial parameters
    c_hat = 1
    initial_beta=numeric(4)
    model_results=data.frame()
    total_models <- length(model_data$models_ran)
    
    for(i in model_data$models_ran){
      current_model <- model_data$models_ran[i]
      #update progress
      progress$set(
        message = paste("Processing Model", model_list[i]),
        detail = paste("Model:", i),
        value = (i-1)/total_models
      )
      
      model_num=i
      starttime<-Sys.time()
      {gc()
        optim_results <- cpp_optim(
          beta = initial_beta,
          ch = as.matrix(model_data$ch),
          cap_X = as.matrix(model_data$cap_X[[i]]),
          surv_X = as.matrix(model_data$surv_X[[i]]),
          ints = intervals)
      }
      endtime<-Sys.time()
      optim_speed<-endtime-starttime 
      
      est_escapement<-total_escapement(model_data$ch,
                                       optim_results$par,
                                       model_data$cap_X[[i]],
                                       model_data$surv_X[[i]],
                                       ints=intervals,
                                       subsampling_weeks=subsampling_weeks,
                                       observed_per_week=observed_per_week,
                                       tagged_per_week=tagged_per_week)
      
      ######################################################
      #calculate model fit statistics (aic, qaic, qaicc)
      ######################################################
      loglik<-optim_results$value
      ic<-model_data$ch
      nan=nrow(ic)
      ns=ncol(ic)
      
      fit_results <- cjs_fit_simple(ch=prepped_data$ch,
                                    beta=optim_results$par,
                                    cap_X=model_data$cap_X[[i]],
                                    surv_X=model_data$surv_X[[i]],
                                    ints=intervals)

      c_hat<-fit_results
      
      n_params<-length(beta)
      
      AIC=2*loglik+2*n_params
      AICc=AIC+((2*n_params)*(n_params+1))/(nan-n_params-1)
      QAIC=((2*loglik)/c_hat)+(2*n_params)
      QAICC=QAIC+(2*n_params*(n_params+1))/(nan-n_params-1)
      
      if(input$use_boots == TRUE && input$boot_input!=0){
        
        #update progress for bootstrap
        progress$set(
          message = paste("Processing Model", model_list[i]),
          detail = paste("Bootstrapping Model:", i),
          value = (i-0.5)/total_models
        )
        
        boot_start<-Sys.time()
        
        #add bootstrap progress updates
        withProgress(message = 'Running bootstrap...', value = 0, {
          boot_results<-CJS_bootstrap(input$boot_input,
                                      ic,
                                      model_data$cap_X[[i]],
                                      model_data$surv_X[[i]],
                                      ints=intervals,
                                      subsampling_weeks=subsampling_weeks,
                                      observed_per_week=observed_per_week,
                                      tagged_per_week=tagged_per_week,
                                      progress_callback = function(iter) {
                                        incProgress(1/input$boot_input, 
                                                    detail = paste("Bootstrap iteration", iter))
                                      })
        })
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
        
      }
      
      d<-data.frame("escapement"=ceiling(est_escapement$escapement),
                    "optim_speed"=round(optim_speed,2),
                    cap_beta1=signif(optim_results$par[1],3),
                    cap_beta2=signif(optim_results$par[2],3),
                    surv_beta1=signif(optim_results$par[3],3),
                    surv_beta2=signif(optim_results$par[4],3),
                    loglik=signif(loglik,3),
                    model=model_data$models_ran[i],
                    "AIC"=signif(AIC,3),
                    "AICc"=signif(AICc,3),
                    "QAIC"=signif(QAIC,3),
                    "QAICC"=signif(QAICC,3),
                    "c_hat"=signif(c_hat,3)
      )
      if(input$boot_input!=0 && input$use_boots==T){
        d<-data.frame(d,lower_ci=ci$lower_ci,
                      upper_ci=ci$upper_ci,
                      input$boot_input,
                      boot_speed)
      }
      model_results<-model_results%>%
        rbind(d)
      
      #update progress
      progress$set(
        value = i/total_models,
        detail = paste("Completed Model:", current_model)
      )
    }
    
    #final progress update
    progress$set(
      value = 1,
      detail = "All models processed!"
    )
    
    #model results output
    output$model_results <- renderDT({
      (model_results)
    })

  })
  
  observeEvent(input$run_selected_models,{
    if(input$use_boots == TRUE && input$boot_input!=0){
      
      #update progress for bootstrap
      progress$set(
        message = paste("Processing Model", model_list[i]),
        detail = paste("Bootstrapping Model:", i),
        value = (i-0.5)/total_models
      )
      
      boot_start<-Sys.time()
      
      #add bootstrap progress updates
      withProgress(message = 'Running bootstrap...', value = 0, {
        boot_results<-CJS_bootstrap(input$boot_input,
                                    ic,
                                    model_data$cap_X[[i]],
                                    model_data$surv_X[[i]],
                                    ints=intervals,
                                    subsampling_weeks=subsampling_weeks,
                                    observed_per_week=observed_per_week,
                                    tagged_per_week=tagged_per_week,
                                    progress_callback = function(iter) {
                                      incProgress(1/input$boot_input, 
                                                  detail = paste("Bootstrap iteration", iter))
                                    })
      })
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
      
    }
  })
  
}

#run the app
shinyApp(ui = ui, server = server)
