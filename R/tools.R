# Title     : Tools for fitting loading
# Objective : TODO
# Created by: konstantin
# Created on: 03.04.18

#' @import pracma
#' @import forecast
RAD <- 15.0
C <- 299792458.0


boundaries <- function(a, b, l) l[l >= a & l <= b]





plotsInterpolation <- function(x, y, s, t, actionFunc, yR) {
    yInterpolated <- approx(s, t, xout = x)
    yInterpolated$y[is.na(yInterpolated)] <- 0.0
    yA <- actionFunc(y, yInterpolated$y)
    yA
}

# elements fitting group

# convert sigma to FWHM
sigmaToFWHM <- function(sigma) 2 * sigma * (2 * log(2)) ^ 0.5

# convert a to km/s
aTokms <- function(a0, da) {
  interm <- (1 + da/a0) ^ 2
  # speed of light, check
  C * (interm - 1) / (interm + 1) * 1e-3
}

# equiwalent width calculations
equiwalentWidth <- function(spectrum, continuumFlux) {

  interpol <- plotsInterpolation(spectrum$lambda, spectrum$flux, spectrum$lambda, rep(0, length(spectrum$lambda)), function(x, y) x - y)
  equWidth <- plotsInterpolation(spectrum$lambda, interpol, spectrum$lambda, continuumFlux, function(x, y) x / y)
  trapz(x = spectrum$lambda, y = equWidth)
}



# leas squares method group
# sucks, TODO  test with library
movingAverage <- function(v, pt=3) {
    len <- length(v)
    result <- rep(0, len)
    for(i in seq(len)) {
        if(i - pt <= 0)
            slice = v[1:(i + pt)]
        else if(i + pt > len)
            slice = v[(i - pt):len]
        else
             slice = v[(i - pt):(i + pt)]
        result[i] <- mean(slice)
    }

    return(result)
}

# movingAverage <- function(v, n = 3) {
#
#   vFiltered <- ma(v, order = n, centre = T)
#   startInds <- 1:(n / 2)
#   endInds <- (length(v) - (n / 2) + 1):length(v)
#   vFiltered[startInds] <- v[startInds]
#   vFiltered[endInds] <- v[endInds]
#   vFiltered
# }

hoursConv <- function(degrees) {
    temp <- degrees / RAD
    return(temp + temp %% 1)
}

# convert degrees to HMS
degreesToSI <- function(degrees) {
    valConv <- function(func, val) return(as.integer(func(val)))
    hours <- as.integer(degrees)
    minutes <- valConv(func=function(x) (x %% 1) * 60, degrees)
    seconds <- degrees - (degrees ) %% 1 * 60 %% 1 * 60
    return(list(H=hours, M=minutes, S=seconds))
}
