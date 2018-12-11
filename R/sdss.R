# Title     : TODO
# Objective : TODO
# Created by: konstantin
# Created on: 21.01.18

#read data from sdss server


UNIT <- 1e-17

extractFITSHeaderValue <- function(header, key) {
    valueIndex <- which(header == key)
    return(header[valueIndex + 1])
}


# restore lambda vector
# FIX error in restoring
# FITS keeps data in log form: l = 10 ^ (a + b* i)
restoreLambda <- function(wavelength, a, b) {
    len <- 0:(wavelength - 1)
    return(10.0 ^ (a + b * len))
}




getHeaderInfo <- function(fitFile) {
    header <- fitFile$hdr
    a <- as.numeric(extractFITSHeaderValue(header, 'COEFF0'))
    b <- as.numeric(extractFITSHeaderValue(header, 'COEFF1'))
    mjd <- as.numeric(extractFITSHeaderValue(header, "MJD"))
    fiber <- as.numeric(extractFITSHeaderValue(header, "FIBERID"))
    plate <- as.numeric(extractFITSHeaderValue(header, "PLATEID"))
    z <- as.numeric(extractFITSHeaderValue(header, "Z"))
    ra <- as.numeric(extractFITSHeaderValue(header, "RAOBJ"))
    dec <- as.numeric(extractFITSHeaderValue(header, "DECOBJ"))
    type <- as.numeric(extractFITSHeaderValue(header, "OBJTYPE"))
    list(mjd=mjd, a=a, b=b, fiber=fiber, plate=plate, z=z, ra=ra, dec=dec, type=type)
}

getImageInfo <- function(fitFile) {
    noise <- fitFile$imDat[, 3] * UNIT
    flux <- fitFile$imDat[, 1] * UNIT
    bitmask <- fitFile$imDat[, 4]
    list(flux=flux, noise=noise, bitmask=bitmask)
}

#' @export
readQuasarInfo <- function(filepath) {
    # since we need first header, is the best choise
    fitFile <- FITSio::readFITS(filepath, hdu = 1)
    headerInfo <- getHeaderInfo(fitFile)
    imageInfo <- getImageInfo(fitFile)

    # restoring lambda from params
    lambda <- restoreLambda(length(imageInfo$flux), as.double(headerInfo$a), as.double(headerInfo$b))
    imageInfo$lambda <- lambda
    list(specs=headerInfo, values=imageInfo)
}





# function to download files from server

# function to collect data from ???


# function to convert quasar specs to filename
qsoParamsToFilename <- function(mjd, plate, fiber) {
    return(sprintf(
        "sp-Spec-%d-%.5d-%.4d.fits",
        as.integer(mjd),
        as.integer(plate),
        as.integer(fiber)
    ))
}


#substitute for getspecbyname
getQuasarSpecsAdress <- function(mjd, plate, fiber, dataRelease='default') {
    sdssServerCatalogue <- config['data_release'][paste('DR_', dataRelease)]
    fitsFilename <- qsoParamsToFilename(mjd, plate, fiber)
    address <- file.path(config['servers']['main_server'], sdssServerCatalogue, plate, '1d', fitsFilename)
    return(address)
}

# separated from getspecbyname
downloadFitsFile <- function(address, folder) {
    fileName <- basename(address)
    dir.create(folder, reqursively=T)
    fileDestination <- file.path(folder, fileName)
    download.file(address, fileDestination)
    return(fileDestination)
}

getQuasarFile <- function(mjd, plate, fiber, dataRelease='default', folder=getwd()) {
    address <- getQuasarSpecsAdress(mjd, plate, fiber, dataRelease)
    return(downloadFitsFile(address, folder))
}



getFITSInfo <- function(filename) {
    header <- parseHdr(filename)
    z <- header['Z']
    zError <- header['Z_ERR']
    rightAssention <- header['RAOBJ']
    declination <- header['DECOBJ']
    mjd <- header['MJD']
	plate <- header['PLATEID']
	fiber <- header['FIBERID']
	type <- header['OBJTYPE']
    objectName <- qsoParamsToFilename(mjd, plate, fiber)
    return(list(name=objectName, z=z, zError=zError, ra=rightAssention, dec=declination, type=type))
}

getRedshift <- function(filename) {
    header <- getFITSInfo(filename)
    return(header$z)
}

