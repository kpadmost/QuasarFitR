library(QuasarRfit)





filtering <- function(wavelength, flux, error, continuumWindows) {
    # filter wavelengths with continuum windows

    filteredFluxW <- filterFromWindows(wavelength, flux, continuumWindows)
    # filter error vec
    filteredFluxErrorW <- filterFromWindows(wavelength, error, continuumWindows)

    # clear from nonlogar wavel
    paragon <- filteredFluxW$y
    filteredWavelength <- filterFromValues(paragon, filteredFluxW$x)
    filteredFlux <- filterFromValues(paragon, filteredFluxW$y)
    filteredFluxError <- filterFromValues(paragon, filteredFluxErrorW$y)
   list(
     lambda=filteredWavelength,
     flux=filteredFlux,
     error=filteredFluxError
   )
}

#fitContinuum(wavelength, flux, continuumFlux)

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




templ <- function(template, expandedF) {
    plot(template$lambda, template$flux)
    plot(template$lambda, expandedF)
}



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

#TODO: add error matrix fix
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
