#This script is taken from the escapeMR package and repurposed here to be used for methods comparisons so that 
#our bootstrap iteration data can be run multiple times and provide parameter outputs and loglik estimates.

sex=as.numeric(surv_X_iteration[,1])
length=as.numeric(cap_X_iteration[,1])
intervals=NULL
if(is.null(intervals)){
  intervals <- rep(1,ncol(ch_iteration)-1)
}

#function calls
cap.model = as.formula("~ivar(length,ns)")
surv.model =as.formula("~ivar(sex,ns)")

survival=surv.model
capture=cap.model
histories=as.matrix(ch_iteration)

cap.init=NULL
nhat.v.meth=1
sur.init=NULL
c.hat=-1.0
link="logit"
group=NULL
chat=-1
ns = dim(ch_iteration)[2]
nan = dim(ch_iteration)[1]
Chat_user = chat

control=mra.control()
covars <- F.cr.model.matrix(capture, survival, nan, ns )  
nx <- covars$n.cap.covars  #nx and ny include intercept.  This is total number of parameters
ny <- covars$n.sur.covars

if( length(union( unique(histories), c(0,1,2))) > 3 ) stop("Capture histories must consist of 0's, 1's, and 2's only.")

if( missing(cap.init) ){
  cap.init <- rep(0,nx)
} else if(length(cap.init) < (nx) ){
  cap.init <- c(cap.init, rep(0, nx-length(cap.init)))
} 
if( missing(sur.init) ){
  sur.init <- rep(0,ny)
} else if(length(sur.init) < (ny) ){
  sur.init <- c(sur.init, rep(0, ny-length(sur.init)))
} 
#   Set up the tolerance vector, if not specified, or if not long enough
if( length(control$tol) < (nx+ny) ){
  control$tol <- rep(control$tol, trunc((nx+ny) / length(control$tol))+1)[1:(nx+ny)]
} else if( length(control$tol > (nx+ny)) ){
  control$tol <- control$tol[1:(nx+ny)]
}
if( missing( group )){
  group <- rep(1, nan)
  ng <- 1
} else {
  ng <- length(unique(group))
}
vif <- c.hat

#   Do the estimation, but first allocate room for answers
loglik <- chisq.vif <- df.vif <- 0
parameters <- se.param <- rep(0, nx + ny )
covariance <- matrix( 0, nx+ny, nx+ny )
p.hat <- se.p.hat <- s.hat <- se.s.hat <- matrix( 0, nan, ns )
exit.code <- cov.code <- 0
n.hat <- se.n.hat <- rep(0, ns)
# on entry to .Fortran, maxfn is maximum number of function evals.  On return, maxfn is actual number of evaluations.
maxfn <- control$maxfn  

#   Re-code the link specification to integers
if( link=="logit" ){
  link.code <- 1
} else if( link == "sine" ){
  link.code <- 2
} else if( link == "hazard" ){
  link.code <- 3
} else {
  stop("Unknown link function specified.")
}

if(control$trace) cat( "Calling MRA DLL to maximize likelihood.  Please wait...\n")
ans <- .Fortran( "cjsmod", 
                 nan         = as.integer(nan), 
                 ns          = as.integer(ns), 
                 nx          = as.integer(nx), 
                 ny          = as.integer(ny), 
                 ng          = as.integer(ng), 
                 histories   = as.integer(histories), 
                 group       = as.integer(group), 
                 algorithm   = as.integer(control$algorithm), 
                 cov.meth    = as.integer(control$cov.meth), 
                 link        = as.integer(link.code),
                 nhat.v.meth = as.integer(nhat.v.meth), 
                 capX        = as.double(covars$capX), 
                 survX       = as.double(covars$survX), 
                 cap.init    = as.double(cap.init), 
                 sur.init    = as.double(sur.init), 
                 maxfn       = as.integer(maxfn),
                 beta.tol.vec= as.double(control$tol), 
                 loglik      = as.double(loglik), 
                 vif         = as.double(vif), 
                 chisq.vif   = as.double(chisq.vif), 
                 df.vif      = as.double(df.vif), 
                 parameters  = as.double(parameters),
                 se.param    = as.double(se.param), 
                 covariance  = as.double(covariance), 
                 p.hat       = as.double(p.hat), 
                 se.p.hat    = as.double(se.p.hat), 
                 s.hat       = as.double(s.hat), 
                 se.s.hat    = as.double(se.s.hat), 
                 n.hat       = as.double(n.hat), 
                 se.n.hat    = as.double(se.n.hat), 
                 exit.code   = as.integer(exit.code), 
                 cov.code    = as.integer(cov.code), 
                 intervals   = as.double(intervals), 
                 PACKAGE="mra" 
) 


