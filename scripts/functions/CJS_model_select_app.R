#' @title CJS_model_select_app
#'
#' @description function that sets model parameters using ch and covariate matrices.
#' 
#' @param covars_used checks uploaded data to see if covariates.
#' 
#' @param sex_matrix matrix of sex covariate data.
#' 
#' @param lengths_matrix matrix of lengths covariate data.
#' 
#' @param ch matrix of capture histories.
#' 
#' @return list with capture histories and covariate matrices for each model possible
#' for available covariates.

CJS_model_select_app <- function(covars_used,
                                 sex_matrix,
                                 lengths_matrix, 
                                 ch) {
  models = c(
    "constant capture and survival rates",
    "constant capture rate and survival related to sex",
    "constant capture rate and survival related to length",
    "capture related to sex and constant survival rate",
    "capture related to length and constant survival rate",
    "capture related to sex and survival related to length",
    "capture related to length and survival related to sex",
    "capture related to sex and survival related to sex",
    "capture related to length and survival related to length"
  )
  
  model_covars = list(cap_X = NULL, surv_X = NULL)
  
  for (m in 1:length(models)) {
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
    'ch' = ch,
    'cap_X' = model_covars$cap_X,
    'surv_X' = model_covars$surv_X
  ))
}