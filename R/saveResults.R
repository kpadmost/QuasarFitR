#' @importFrom RJSONIO toJSON
saveResults <- function(fitResults, objectName) {
  write(toJSON(fitResults), paste(objectName, '.json', sep=''))
}
