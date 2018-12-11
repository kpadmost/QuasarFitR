
skylinesElimination <- function(params) {
  EMISSION_LINE <- 0x40000000
  bitmask <- params$bitmask
  eliminateElements <- which(bitmask != 0 & bitmask != EMISSION_LINE)
  params$lambda <- params$lambda[-eliminateElements]
  params$noise <- params$noise[-eliminateElements]
  params$flux <- params$flux[-eliminateElements]
  return(params)
}

filterFromValues <- function(paragon, unfiltered, value= 0.0, comparator = function(a, b) a > b) {
  filtered <- comparator(paragon, value)
  unfiltered[filtered]
}

filterFromWindows <- function(x, y, windows) {
  xResult <- c()
  yResult <- c()
  for(window in windows) {
    filtered <- slice(window[1], window[2], x, y)
    xResult <- c(xResult, filtered$x)
    yResult <- c(yResult, filtered$y)
  }
  list(x=xResult, y=yResult)
}

# function which returns a slice of limited array x, y(x), where values of function are within [a, b]
slice <- function(a, b, x, y) {
  indices <- x >= a & x <= b
  
  result <- list(x=x[indices], y=y[indices], indices=indices)
}

# TODO: filter with wavelenght window
getSpectrumElementPart <- function(wavelengthWindow, spectrum) {
  if( wavelengthWindow[1] < spectrum$lambda[1] ||
      wavelengthWindow[2] > spectrum$lambda[length(spectrum$lambda)]) {
    return(list(fit=F))
  }
  indexes <- spectrum$lambda >= wavelengthWindow[1] & spectrum$lambda <= wavelengthWindow[2]
  
  wavelengthPart <- spectrum$lambda[indexes]
  fluxPart <- spectrum$flux[indexes]
  list(spectrum=list(lambda=wavelengthPart,flux=fluxPart), fit=T, indices=indexes)
}