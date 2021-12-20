#' @title varreduction
#' @description Implement a function to generate samples from standard normal distribution,and then use antithetic variables.In the end we will get the percent of variance reduction
#' @param m times of repetition
#' @param x upper limit of Monte Carlo Integration
#' @param R the sample size
#' @return the reduction of the variance through using the antithetic variables in Monte Carlo Integration
#' @examples
#' \dontrun{
#' m<-1000
#' R<-10000
#' x<-2
#' var.redu(m,x,R)
#' }
#' @export
var.redu<-function(m,x,R){
  zt1 <- zt(x, anti = FALSE)
  zt2 <- zt(x)
  zt1 <- zt2 <- numeric(m)
  for (i in 1:m) {
    zt1[i] <- zt(x, R=1000, anti = FALSE)
    zt2[i] <- zt(x, R=1000)
  }
  print((var(zt1) - var(zt2))/var(zt1))
}
zt <- function(x,R=1000, antithetic = TRUE) {
  u <- runif(R/2)
  if (!antithetic) 
    v <- runif(R/2) 
  else
    v <- 1 - u
  u <- c(u, v)
  cdf <- numeric(length(x))
  for (i in 1:length(x)) {
    g <- x[i]*exp(-(u*x[i]^2/2))
    cdf[i] <- mean(g)/sqrt(2*pi)+1/2
  }
  cdf
}
#' @title rejpercent
#' @description use the Metropolis-Hastings sampler to generate a sample from a Rayleigh distribution,and then evaluate the rejection percent of the candidate point
#' @param m the length of chain
#' @param sigma the parameter of Rayleigh distribution
#' @return the rejection percent of the final chain
#' @examples
#' \dontrun{
#' m<-1000
#' sigma<-5
#' rej.perc(m,sigma)
#' }
#' @export
rej.perc<-function(m,sigma){
  x<-numeric(m)
  x[1]<-rchisq(1,df=10)
  k<-0
  u<-runif(m)
  for(i in 2:m){
    xt<-x[i-1]
    y<-rchisq(1,df=xt)
    gyj1<-rayleigh(y,sigma)*dchisq(xt,df=y)
    gyj2<-rayleigh(xt,sigma)*dchisq(y,df=xt)
    if(u[i]<=gyj1/gyj2)
      x[i]<-y
    else{
      x[i]<-xt
      k<-k+1
    }
  }
  print(k/m)
}
rayleigh<-function(x,sigma){
  if(any(x<0))
    return(0)
  stopifnot(sigma>0)
  return((x/sigma^2)*exp(-x^2/(2*sigma^2)))
}
