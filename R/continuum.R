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
