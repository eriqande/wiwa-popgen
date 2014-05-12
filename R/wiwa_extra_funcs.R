





# Here is a function to write the data frames to gsi_sim format.  Supply a df with just the genetic
# data in it with the individual IDs as the rownames.  If you supply a list of pop-names 
# (preferably as a factor with levels ordered as you want the pops to appear) to the 
# pop.ize.them variable, then it splits them into pops, ordered as in the factor
df2gsi.sim <- function(d, pop.ize.them=NULL, outfile="gsi_sim_file.txt") {
	write(dim(d)/c(1,2), outfile)
	write(cbind(names(d)[seq(1, ncol(d), 2)]), outfile, append=T)
	if(is.null(pop.ize.them)) {
		write("POP  Mixture", outfile, append=T)
		write.table(d, outfile, append=T, quote=F, row.names=T, col.names=F)
	}
	else {
		idx.list <- split(1:nrow(d), pop.ize.them)
		catchit <- lapply(names(idx.list), function(x) {
				write(paste("POP", x), outfile, append=T)
				write.table(d[idx.list[[x]],] , outfile, append=T, quote=F, row.names=T, col.names=F)
			}
		)
	}	
}




# read the output of gsi_sim self-assignment run and extract the results into a list of data frames,
# one for posteriors and the other for logls.  This assumes that the individual IDs in the gsi_sim
# output are what they would be for the pipeline, so you can get the population name by stripping
# the numbers off them.
gsi_simSelfAss2DF <- function(file="self-ass-output.txt") {
	x <- readLines(file)
	x1<-x[grep("UNSORTED_SELF_ASS_LIKE_GC_CSV:", x)]  # get just the lines we want
	x2<-read.table(textConnection(x1), sep=";", stringsAsFactors=F)
	numpops <- (ncol(x2)-4)/4  # remove one column for each of ID, NumMissingLoci, NonMissingLocusNames, and MissingLocusNames
	popnames <- as.character(x2[1,seq(from=2, by=3, length.out=numpops)])  # names of all the pops in the baseline, in that order
	IDs <- gsub("UNSORTED_SELF_ASS_LIKE_GC_CSV:/", "", x2$V1)
	FromPops <- gsub("[0-9]*","", IDs) 
	FromPops <- factor(FromPops, levels=popnames)  # make them a factor that preserves the order they went into gsi sim with
	Posteriors <- x2[,seq(from=3, by=3, length.out=numpops)]/100.0
	LogLs <- x2[, seq(from=2+3*numpops, length.out=numpops)]
	NumLoci <- x2[, ncol(x2)-2]

	# send result back as a list of data frames, either Posteriors or LogLs
	lapply(list(Post=Posteriors, LogLs=LogLs), function(z) {
		names(z) <- popnames
		ll <- cbind(PopulationOfOrigin=FromPops, NumberOfLoci=NumLoci, z)
		rownames(ll) <- IDs
		ll
		}
	)
}


# given a data frame x that has columns represententing scores (like posterior
# probs) for each population in the pops.str (which should be a character
# vector of population names) which are themselves grouped into reporting
# units according to rep.units (preferably a factor), this returns a new data frame
# with all the non-pop columns first and then the pop columns aggregated by summing
# them up.
gsi_aggByScoresByRepUnit <- function(x, pops.str, rep.units) {
	other.cols <- x[ !(names(x) %in% pops.str) ]
	grps <- split(pops.str, rep.units)
	sapply(grps, function(z) rowSums(x[z]))
	cbind(other.cols, sapply(grps, function(z) rowSums(x[z])))
}


# pass this a data frame d that has the posterior probabilities in columns
# with indexes (or names) in variable columns, and this will return a vector
# of the names of the column in which the max posterior occurs.  And it also 
# returns the maximum posterior associated with that assignment
gsi_simMaxColumnAndPost <- function(d, columns) {
	tmp <- as.matrix(d[,columns])
	nam <- colnames(tmp)
	maxcolumns <- nam[max.col(tmp)]  # hold this temporarily
	maxcolumns <- factor(maxcolumns, levels=nam)  # preserve the order they were in d
	data.frame(MaxColumn=maxcolumns, MaxPosterior=apply(tmp, 1, max))
}


# given a factor fp giving the "FromPopulation" and a factor mc given the "MaxColumn"
# unit assigned to, and a numeric vector mp of the maximum posterior assignment scores,
# and a numeric vector of cutoffs for the max posterior, this returns a list indexed by
# cutoffs of 
#   - The number of individuals assigned
#   - The assignment table of the assigned individuals
gsi_simAssTableVariousCutoffs <- function(fp, mc, mp, cutoffs=seq(0,.95, by=.05)) {
	names(cutoffs) <- paste("Cutoff",cutoffs, sep="_")
	lapply(cutoffs, function(x) {
			f <- fp[mp>=x]
			m <- mc[mp>=x]
			AssTable=table(f,m)
			names(dimnames(AssTable)) <- c("From", "To")
			list(NumAssigned=sum(mp>=x), AssTable=AssTable)
		}
	)
}


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
