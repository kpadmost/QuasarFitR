library(QuasarRfit)

continuumFitting <- function(wavelength, flux, error, lambdaAmplitude, continuumWindows) {
    filtered <- filtering(wavelength, flux, error, continuumWindows)
    # change to log scale
    # use lm, check how works substitute with least square interface
    logf <- log10(filtered$flux)
    logl <- log10(filtered$lambda)
    continuumReglin <- reglin(logl, logf)
    #get reglin at lapt
    if(lambdaAmplitude > 0) {
      logla <- logl - log10(lambdaAmplitude)
      ampReglin <- reglin(logla, logf)
      ampReglin <- reglinFix(ampReglin, mode='amp')
    }
    # TODO: check coeffs whether is correct
    continuumReglin <- reglinFix(continuumReglin)

    cfunFiltered <- continuumReglin$b * (filtered$lambda ^ continuumReglin$a)
    cfun <- continuumReglin$b * (wavelength ^ continuumReglin$a)
    dcfun <- continuumReglin$sib * (wavelength ^ continuumReglin$a) + continuumReglin$sib * continuumReglin$a * log(wavelength)

    continuumReglin <- reglinFix(continuumReglin)
    result <- list(
      cfun = cfun,
      dcfun = dcfun,
      fitParams = list(
        cReglin = continuumReglin,
        chisq = chisqsDistribution(filtered$flux, cfunFiltered, filtered$error, len =length(filtered$flux), reduced=T)
      )
    )
    if(lambdaAmplitude > 0)
      result$fitParams$reglin <- ampReglin
    result
}


ferrumTemplateFitting <- function(spectrum, continuum, feFitParams) {
    template <- loadDefaultFeTemplate()
    ferrumEmissionWindows <- loadDefaultFeWindows()
   # interpolation + convolute template
    expandedTemplate <- expandFeTemplate(spectrum, template, feFitParams)
    # filtering from ferrumEmissionWindows ^ zeros ^ pure ferrum emission
    filtered <- filterTemplate(spectrum, expandedTemplate, ferrumEmissionWindows, feFitParams$fitType)
    scale <- scaleFeTemplate(filtered$spectrum, filtered$template)

    feTemplateFitted <- expandedTemplate$flux * scale$a

    #  fitting template params
    chisqReducedFull <- chisqsDistribution(spectrum$flux, feTemplateFitted, spectrum$error, len=length(spectrum$flux), reduced=T)
    chisqReducedFiltered <- chisqsDistribution(filtered$spectrum$flux, filtered$template$flux * scale$a, filtered$spectrum$error, len=length(spectrum$flux), reduced=T)
    ewFull <- equiwalentWidth(list(lambda = spectrum$lambda, flux=feTemplateFitted), continuum)
    # fitting template in range params
    fitRange <- feFitParams$feFitRange
    spectrumFRange <- slice(a = fitRange[1], b = fitRange[2], spectrum$lambda, spectrum$flux)
    indices <- spectrumFRange$indices

    wavelengthsFRange <- spectrumFRange$x
    templateFRange <- feTemplateFitted[indices]
    continuumFRange <- continuum[indices]
    fluxErrorFRange <- spectrum$flux[indices]
    chisqsReducetFRange <- chisqsDistribution(spectrumFRange$y,
                                              templateFRange,
                                              fluxErrorFRange,
                                              len = length(templateFRange),
                                              reduced = T)
    ewFRange <- equiwalentWidth(list(lambda = wavelengthsFRange, flux=templateFRange), continuumFRange)
    list(
      feTemplate = feTemplateFitted,
      fitParams = list(
        scaleRate = scale$a,
        reducedChisqsFull = chisqReducedFull,
        reducedChisqsWindows = chisqReducedFiltered,
        ewsFull = ewFull,
        reducedChisqsFeRange = chisqsReducetFRange
      )
    )
}









elementFitting <- function(element, spectrum, continuum, fluxError) {
  # consider indices as update
  spectrumPartsInfo <- getSpectrumElementPart(c(element$left, element$right), spectrum)
  indices <- spectrumPartsInfo$indices
  if(!spectrumPartsInfo$fit) {

    return(list(fit=F))
  }
  # fit gaussian
  fitGaussParams <- fitElementGaussian(spectrumPartsInfo$spectrum, spectrumFull = spectrum, element)
  gaussParams <- fitGaussParams$fitGaussParams
  gaussFlux <- fitGaussParams$elementGauss

  fwhm <- aTokms(gaussParams$mu, sigmaToFWHM(gaussParams$sigma))
  equWidth <- equiwalentWidth(list(lambda=spectrum$lambda[indices], flux=gaussFlux[indices]), continuum[indices])
  chisqsDist <- chisqsDistribution(spectrum$flux, gaussFlux, fluxError)
  list(
    fwhm = fwhm,
    ew = equWidth,
    chisqs = chisqsDist,
    fit=T
  )
}



rFitting <- function(filename) {
    # load params
    quasarParameters <- readQuasarInfo(filename)
    # get redshift
    z <- quasarParameters$specs$z
    # exclude skylines
    quasarParameters$values <- skylinesElimination(quasarParameters$values)
    wavelength <- quasarParameters$values$lambda
    # eliminate  redshift
    wavelength <- wavelength / (z + 1)

    # srednia kroczaca
    flux <- movingAverage(quasarParameters$values$flux, 10)
    fluxError <- quasarParameters$values$noise
    #plot(flux)
    # iterations with fitting cont & remove iron emissions
    lambAmp <- 3000
    contWindows <- loadDefaultContinuumWindows()
    # TODO : add result quality check
    feFitted <- list(feTemplate = rep(0, length(wavelength)))
    feFitParams <- loadDefaultFeFitParams()
    for(i in seq(3)) {
        # fit in cont windows
        continuumParams <- continuumFitting(wavelength = wavelength, flux=flux - feFitted$feTemplate, error=fluxError, lambdaAmplitude = lambAmp, continuumWindows = contWindows)
        # remove flux
        #fluxCless <- flux - continuumInfo$flux
        # remove iron emissions
        feFitted <- ferrumTemplateFitting(list(lambda = wavelength, flux=flux - continuumParams$cfun, error=fluxError), continuumParams$cfun, feFitParams)
    }
    globalFitParams <- append(list(continuum = continuumParams$fitParams), list(fe=feFitted$fitParams))
    fluxFitted <- flux - continuumParams$cfun - feFitted$feTemplate
    # for each possible emission element
    elements <- loadDefaultSpectralLines()

    elementsFitParams <- list()
    for(row in 1:nrow(elements)) {
      currentElement <- elements[row,]
      fitParams <- elementFitting(
        currentElement,
        spectrum = list(lambda = wavelength, flux=fluxFitted),
        continuum = continuumParams$cfun,
        fluxError = fluxError
      )
      if(fitParams$fit) {
        fitParams$fit <- NULL
        elementsFitParams <- append(elementsFitParams, list(list(element=currentElement$name, params=fitParams)))
      } else
        print(paste('element', currentElement$name, 'does not fit into object!'))
    }
    for(e in elementsFitParams) {
      print(e$element)
    }
    globalFitParams <- append(globalFitParams, list(elements = elementsFitParams))
    #drawFit(list(wavelength=wavelength, flux=flux), continuumParams$cfun, feFitted$feTemplate, elementsFitParams)
    saveResults(globalFitParams, sub('.fit', '', filename))
}

main <- function() {
    # read file
    rFitting('spSpec-51792-0354-147.fit')
}
main()
