expandFeTemplate <- function(spectrum, template, feFitParams) {

  # sigmaConv

  fwhmConv <- (feFitParams$fwhmn ^ 2 - feFitParams$fwhmt ^ 2) ^ 0.5

  sigmaConv <- fwhmConv / (2 * ((2 * log(2)) ^ 0.5)) / C * 1e3
  # miConv is set on a middle of a vector
  lambdaL <- log10(template$lambda)
  miConv <- lambdaL[length(lambdaL) / 2]

  expandedL <- gaussPeak(lambdaL, sigmaConv, miConv)
  expandedFlux <- filter(template$flux, expandedL, method = 'convolution', circular = T)
  interP <- plotsInterpolation(spectrum$lambda, rep(0, length(spectrum$lambda)), template$lambda, expandedFlux, actionFunc = function(x, y) x + y, spectrum$flux)
  interP[is.na(interP)] <- 0.0

  list(lambda=spectrum$lambda, flux=interP)
}

scaleFeTemplate <- function(spectrum, templateExpanded) {
  # find scaling factor
  a <- leastq(spectrum$flux, templateExpanded$flux)$a
  scalingParams <- nleastq(templateExpanded$flux, spectrum$flux, a)
  list(a=scalingParams['a'])
}

filterTemplate <- function(spectrum, template, ferrumEmissionWindows, method = 'fwin') {
  # filter from pure fe emission
  filteredSpectrumFlux <- filterFromWindows(spectrum$lambda, spectrum$flux, ferrumEmissionWindows)
  filteredSpectrumError <- filterFromWindows(spectrum$lambda, spectrum$error, ferrumEmissionWindows)
  filteredTemplate <- filterFromWindows(template$lambda, template$flux, ferrumEmissionWindows)

  filteredWavelength <- filteredSpectrumFlux$x
  filteredFlux <- filteredSpectrumFlux$y
  filteredFluxError <- filteredSpectrumError$y

  filteredTemplateWavelength <- filteredTemplate$x
  filteredTemplateFlux <- filteredTemplate$y

  if (method == 'fwin') {
    paragon <- filteredTemplate$y
    filteredWavelength <- filterFromValues(paragon, filteredWavelength)
    filteredFlux <- filterFromValues(paragon, filteredSpectrumFlux$y)
    filteredFluxError <- filterFromValues(paragon, filteredFluxError)
    filteredTemplateWavelength <- filterFromValues(paragon, filteredTemplateWavelength)
    filteredTemplateFlux <- filterFromValues(paragon, filteredTemplateFlux)
  }
  list(
    spectrum = list(
      lambda = filteredWavelength,
      flux = filteredFlux,
      error = filteredFluxError
    ),
    template = list(
      lambda = filteredTemplateWavelength,
      flux = filteredTemplateFlux
    )
  )
}
