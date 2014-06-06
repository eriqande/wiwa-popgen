## This is a collection of code to deal with the geneland output.
## Note that the geneland output is currently on a jump drive.


load("outputs/WIWA-main-carryover-variables.Rda")


#GLOUT<-"/Volumes/ERIC/WIWA_SEP30_Geneland_Runs/GL-Run2--Sept30-version"
#GLOUT<-"/Users/eriq/Documents/work/assist/kristenruegg/WIWA_popgen/WIWA_popgen_analysis/geneland_outputs/GL-Run-Oct21"

GLOUT<-"/Volumes/eriq/Documents/UnsyncedSimulationOutputs/WIWA-Geneland/GL-Run-Oct21"

OUTP <- "/tmp"

library(Geneland)

source("./R/geneland_helper_funcs.R")

x<-1
post.file <- paste(GLOUT, "/GL-run-sub-",x,"/", "proba.pop.membership.txt", sep="")

# get the sampling coordinates
#coord <- read.table("../data/Sept30_version/Latlong.txt")
coord <- gl.coord


# now read in the posteriors and clip it etc:
pfile <- read.table(post.file, header=T)

# make spatial po


library(maps)
catchit <- lapply(3:10, function(x) {
  PostMode2(
    coordinates=coord, 
    path.mcmc=paste(GLOUT, "/GL-run-sub-",x,"/", sep=""), 
    file=paste(OUTP, "PostModeMap-", x, ".pdf", sep=""),
    dot.cex=.8
  )
  map("world", add=T)
  dev.copy2pdf(file=paste(OUTP, "PostModeMap-", x, ".pdf", sep="/"))
}
)
