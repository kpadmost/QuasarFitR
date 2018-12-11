fitElementGaussian <- function(spectrumPart, spectrumFull, element) {
  gParams <- fitGauss(
    spectrumPart$lambda,
    spectrumPart$flux,
    sigma0=element$c,
    mu0=element$b,
    a0=element$a
  )
  # extract an interesting params

  elementGauss <- generateGauss(spectrumFull$lambda, params=gParams)
  list(elementGauss=elementGauss, fitGaussParams=gParams)
}
