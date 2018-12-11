gaussPeak <- function(x, sigConv, miConv) {
  a <- 1.0 / (sigConv * ((2 * pi) ^ 0.5))
  a * exp(-0.5 * ((x - miConv) ^ 2) / (sigConv ^ 2))
}

# generate gauss distr
generateGauss <- function(x, params=list(), a=1, sigma=1, mu=0) {
  if(!is.na(params)) {
    a=params$a
    sigma=params$sigma
    mu=params$mu
  }
  a*exp(-1/2*(x-mu)^2/sigma^2)
}

# use standart implementation
fitGauss <- function(x, y, mu0, sigma0, a0) {
  # converting the data to valid std form
  dataSpectrum <- data.frame(x=x, y=y)
  res <- nls( y ~ a*exp(-0.5*(x-mu)^2/(sigma^2)), start=c(a=a0, mu=mu0,sigma=sigma0) , data = dataSpectrum, control=list(maxiter=100))
  coeffs <- as.vector(coef(res))
  #TODO: check!
  return(list(a=coeffs[1], mu=coeffs[2], sigma=coeffs[3]))
}