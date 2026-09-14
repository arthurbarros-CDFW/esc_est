##########################################
#get_bootstrap_iterations
##########################################
get_bootstrap_iterations <- function() {
  #ask if user wants to bootstrap
  perform_boot <- tolower(readline("Perform bootstrapping for confidence intervals? (y/n): "))
  
  #validate yes/no response
  while(!perform_boot %in% c("y", "n", "yes", "no")) {
    message("Error: Please answer 'y' or 'n'.")
    perform_boot <- tolower(readline("Perform bootstrapping? (y/n): "))
  }
  
  if (perform_boot %in% c("n", "no")) {
    return(0)  #return 0 if no bootstrapping
  }
  
  #if yes, get bootstrap iterations
  while(TRUE) {
    n_boot <- readline("Enter bootstrap iterations (integer 25-1000): ")
    
    #check if numeric
    if (is.na(suppressWarnings(as.numeric(n_boot)))) {
      message("Error: Input must be a number.")
      next
    }
    
    n_boot <- as.integer(n_boot)
    
    #check if integer
    if (is.na(n_boot) || n_boot != as.numeric(n_boot)) {
      message("Error: Input must be an integer.")
      next
    }
    
    # Check range
    if (n_boot < 25 || n_boot > 1000) {
      message("Error: Input must be between 25 and 1000.")
      next
    }
    
    return(n_boot)
  }
}