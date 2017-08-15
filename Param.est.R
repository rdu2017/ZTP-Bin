
mle.lambda <- function(u,lambda0){
  
  # Maximum likelihood estimator (MLE) of lambda under a zero-truncated Poisson distribution
  #
  # Args:
  #   u: a vector; an element of it is the number of reads fragmented from a CDS associated to the COG of study
  #   lambda0: an initial value for solving the ML equation
  #
  # Return:
  #   a integer, the obtained MLE of lambda

  f<-function(u,lambda) # the log-likelihood function
  {
    lambda-mean(u)*(1-exp(-lambda))  
  }
  
  fprime<-function(u,lambda) # first derivative of the log-likelihood function
  {
    1-mean(u)*exp(-lambda)
  }
  
  lambda=lambda0 # set the initial value to start the Newton-Raphson iteration to find MLE
  repeat
  { 
    lambdaK=lambda-f(u,lambda)/fprime(u,lambda) 
    diff=abs(lambdaK-lambda)

    cat(paste('MLE is obtained, the absolute difference is: ', 
              as.character(diff),
              sep=''))
    cat('\n')
    
    if (diff<1e-6) break else lambda<-lambdaK # the iteration stops when the absolute difference between two consecutive solulions is less than 1e-6
  }
  
  lambdaK # return the obtained MLE of lambda
}


