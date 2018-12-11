library(stats)
# simplest version for this
leastq <- function(x, y) {
  model <- lm(y ~ x)
  list(a=model$coefficients[2], b=model$coefficients[1])
}


nleastq <- function(x, y, a0) {
  df <- data.frame(x = x, y = y)
  model <- nls(y ~ a * x, data = df, start=list(a=a0))
  coefficients(model)
}




reglin <- function(x, y) {
  model <- lm(y ~ x)
  smodel <- summary(model)
  coeffis <- smodel$coefficients
  a <- coeffis[2,1]
  aSigma <- coeffis[2, 2]
  b <- coeffis[1,1]
  bSigma <- coeffis[1, 2]
  ySigma <- sigma(model)
  list(a = a, b = b, sia = aSigma, sib=bSigma, siy = ySigma)
}

#mode - cont |
reglinFix <- function(reglin, mode='amp') {
  reglin$b <- 10.0 ^ reglin$b
  if(mode == 'amp') {
    return(reglin)
  }
  reglin$a <- reglin$a
  reglin$sib <- reglin$b * log(10) * reglin$sib ^ 0.5
  reglin$sia <- reglin$sia ^ 0.5
  reglin$siy <- reglin$siy  ^ 0.5
  reglin
}


chisqsDistribution <- function(x, y, sigma, len=1, reduced = F) {
  chisqs <- sum(((x - y) / sigma) ^ 2)
  if(reduced)
    chisqs <- chisqs / len
  chisqs
}



