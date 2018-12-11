library(ggplot2)
library(tidyr)

drawSpectrum <- function(spectrumSpecs, continuumSpecs=NULL, feTemplate=NULL) {
  plt <- ggplot(spectrumSpecs$lambda, spectrumSpecs$flux)
  if(!is.null(continuumSpecs))
    plt <- plt + ggplot(spectrumSpecs$lambda, spectrumSpecs$flux)
  if(!is.null(feTemplate))
    plt <- plt + ggplot(feTemplate$lambda, feTemplate$flux)
}

drawFit <- function(spectrum, continuum, template, fitElements) {
  df <- data.frame(lambda = spectrum$wavelength, 
                   flux = spectrum$flux, 
                   c = continuum,
                   fe=template,
                   feC = continuum + template,
                   fFi = spectrum$flux - continuum - template
                   )
  df %>% gather(key, value, flux, c, feC) %>% ggplot(aes(x=lambda, y=value, color=key))
  grid.arrange(pl0, pl1, plc, plf)
}