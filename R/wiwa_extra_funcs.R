


# here I have hacked floating.pie from plotrix so it doesn't change
# colors when you have zeroes in your vector.
my.floating.pie <- function (xpos, ypos, x, edges = 200, radius = 1, col = NULL, 
    startpos = 0, shadow = FALSE, shadow.col = c("#ffffff", "#cccccc"), 
    ...) 
{
    if (!is.numeric(x)) 
        stop("floating.pie: x values must be numeric.")
    validx <- which(!is.na(x) & x > 0)
    x <- c(0, cumsum(x[validx])/sum(x[validx]))
    dx <- diff(x)
    nx <- length(dx)
    if (is.null(col)) 
        col <- rainbow(nx)
    else if (length(col) < nx) 
        col <- rep(col, nx)
    else 
    		col <- col[validx]
    		
    xylim <- par("usr")
    plotdim <- par("pin")
    yradius <- radius * (xylim[4] - xylim[3])/(xylim[2] - xylim[1]) * 
        plotdim[1]/plotdim[2]
    bc <- 2 * pi * (x[1:nx] + dx/2) + startpos
    if (shadow) {
        xc <- c(cos(seq(0, 2 * pi, length = edges)) * radius + 
            xpos)
        yc <- c(sin(seq(0, 2 * pi, length = edges)) * yradius + 
            ypos)
        polygon.shadow(xc, yc, col = shadow.col)
    }
    for (i in 1:nx) {
        n <- max(2, floor(edges * dx[i]))
        t2p <- 2 * pi * seq(x[i], x[i + 1], length = n) + startpos
        xc <- c(cos(t2p) * radius + xpos, xpos)
        yc <- c(sin(t2p) * yradius + ypos, ypos)
        polygon(xc, yc, col = col[i], ...)
        t2p <- 2 * pi * mean(x[i + 0:1]) + startpos
        xc <- cos(t2p) * radius
        yc <- sin(t2p) * radius
    }
    invisible(bc)
}
