library(ggplot2)

# pl <- list(y, color)
drawPlots <- function(x, pl1, pl2, c1='blue', c2='red') {
    df1 <- data.frame(x = x, y = pl1)
    df2 <- data.frame(x = x, y = pl2)
    pl <- ggplot(df1, aes(x, y)) + geom_line(color = c1) + geom_line(data = df2, aes(x, y), color = c2)
    show(pl)
}

plotl <- function(x, y) {
    plot(x, y, type = 'n')
    lines(x, y, type='l')
}