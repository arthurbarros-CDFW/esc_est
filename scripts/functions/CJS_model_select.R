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