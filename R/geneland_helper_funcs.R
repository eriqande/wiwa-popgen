## Hack their posterior mode function to change the color of the dots,
## and also to add a map to it.
PostMode2 <- function (coordinates, path.mcmc, plotit = TRUE, format = "pdf", 
    new.dev = TRUE, printit = FALSE, file, main.title = "", dot.pch=21, dot.cex=1.6, dot.col="blue") 
{
    coordinates <- as.matrix(coordinates)
    fileparam <- paste(path.mcmc, "parameters.txt", sep = "")
    param <- as.matrix(read.table(fileparam))
    delta.coord <- as.numeric(param[param[, 1] == "delta.coord", 
        3])
    npopmax <- as.numeric(param[param[, 1] == "npopmax", 3])
    param.postprocess <- as.matrix(read.table(paste(path.mcmc, 
        "postprocess.parameters.txt", sep = "")))
    nxdom <- as.numeric(param.postprocess[1, 3])
    nydom <- as.numeric(param.postprocess[2, 3])
    print(nxdom)
    print(nydom)
    s <- coordinates
    filedom <- paste(path.mcmc, "proba.pop.membership.txt", sep = "/")
    data <- as.matrix(read.table(filedom))
    dom.post <- data[, -(1:2)]
    coord.grid <- data[, (1:2)]
    s[, 1] <- s[, 1] - min(s[, 1])
    s[, 2] <- s[, 2] - min(s[, 2])
    map.dom <- t(apply(dom.post, 1, order))[, npopmax]
    if (plotit) {
        if (new.dev) 
            dev.new()
        frame <- max(max(coordinates[, 1]) - min(coordinates[, 
            1]), max(coordinates[, 2]) - min(coordinates[, 2]))/40
        image(seq(min(coordinates[, 1] - delta.coord/2), max(coordinates[, 
            1] + delta.coord/2), length = nxdom), seq(min(coordinates[, 
            2] - delta.coord/2), max(coordinates[, 2] + delta.coord/2), 
            length = nydom), matrix(map.dom, nrow = nxdom, ncol = nydom, 
            byrow = TRUE), xlab = "", ylab = "", main = "", cex = 1.5, 
            cex.lab = 1.5, col = terrain.colors(npopmax), xlim = c(min(coordinates[, 
                1] - delta.coord/2 - frame), max(coordinates[, 
                1] + delta.coord/2 + frame)), ylim = c(min(coordinates[, 
                2] - delta.coord/2 - frame), max(coordinates[, 
                2] + delta.coord/2 + frame)), asp = 1)
        points(coordinates, pch = dot.pch, cex = dot.cex, col=dot.col)
        title(sub = "Estimated cluster membership")
        title(main = main.title, pch = 16)
    }
    if (printit) {
        if (format == "ps") {
            postscript(file)
        }
        if (format == "pdf") {
            pdf(file)
        }
        frame <- max(max(coordinates[, 1]) - min(coordinates[, 
            1]), max(coordinates[, 2]) - min(coordinates[, 2]))/40
        image(seq(min(coordinates[, 1] - delta.coord/2), max(coordinates[, 
            1] + delta.coord/2), length = nxdom), seq(min(coordinates[, 
            2] - delta.coord/2), max(coordinates[, 2] + delta.coord/2), 
            length = nydom), matrix(map.dom, nrow = nxdom, ncol = nydom, 
            byrow = TRUE), xlab = "", ylab = "", main = "", cex = 1.5, 
            cex.lab = 1.5, col = terrain.colors(npopmax), xlim = c(min(coordinates[, 
                1] - delta.coord/2 - frame), max(coordinates[, 
                1] + delta.coord/2 + frame)), ylim = c(min(coordinates[, 
                2] - delta.coord/2 - frame), max(coordinates[, 
                2] + delta.coord/2 + frame)), asp = 1)
        points(coordinates,  pch = dot.pch, cex = dot.cex, col=dot.col)
        title(main = main.title, sub = "Estimated cluster membership")
        
        # now ECA adds some map stuff in there.
        library(maps)
        map("world", add=T)
        
        
        dev.off()
    }
}

