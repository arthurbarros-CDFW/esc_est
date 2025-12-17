DRAFT. Code to replicate and build upon the CVCS escapement estimating methods utilizing the Cormack-Jolly-Seber modeling methods established by Trent McDonald in the escapeMR package (now unlisted on CRAN). Still in a working draft form.
Updates since October 2024:
1) Added functionality to select multiple various covariate model structures for comparison.
2) Incorporated a 'CJS_run()' function to streamline the data loading, model selection, bootstrapping, and output of model results.
3) Ran multiple tests using carcass survey data from various years and sources to compare outputs, and model failure, between these methods and the escapeMR/MRA approach.
4) Currently working on incorporating model fit tests (AIC, QAICC), which is a work in progress.
