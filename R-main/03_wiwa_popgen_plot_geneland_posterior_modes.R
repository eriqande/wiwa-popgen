## This is a collection of code to deal with the geneland output.
## Note that the geneland output is currently on a jump drive.



#### Load packages, get stored data, source necessary R-scripts, define directories ####
stopifnot(
  library(Geneland, logical.return = TRUE),
  library(maps, logical.return = TRUE)
)
load("./outputs/WIWA-main-carryover-variables.Rda")

if(!all(c("gl.geno", "gl.coord") %in% got.these)) {
  stop("Failed to successfully load variables gl.geno and/or gl.coord from file WIWA-main-carryover-variables.Rda")
}

source("./R/geneland_helper_funcs.R")


#GLOUT<-"/Volumes/ERIC/WIWA_SEP30_Geneland_Runs/GL-Run2--Sept30-version"
#GLOUT<-"/Users/eriq/Documents/work/assist/kristenruegg/WIWA_popgen/WIWA_popgen_analysis/geneland_outputs/GL-Run-Oct21"

GLOUT<-"./outputs"

OUTP <- "./outputs/geneland-posterior-mode-plots"
dir.create(OUTP)


#### NOT SURE WHY THIS IS HERE. I MIGHT BE ABLE TO TOSS THIS CHUNK ####
x<-1
post.file <- paste(GLOUT, "/GeneLandRun-",x,"/", "proba.pop.membership.txt", sep="")

# get the sampling coordinates
coord <- gl.coord


# now read in the posteriors and clip it etc:
pfile <- read.table(post.file, header=T)


#### PRODUCE ALL THE POSTERIOR MODE PLOTS ####
catchit <- lapply(3:10, function(x) {
  PostMode2(
    coordinates=coord, 
    path.mcmc=paste(GLOUT, "/GeneLandRun-",x,"/", sep=""), 
    file=paste(OUTP, "PostModeMap-", x, ".pdf", sep=""),
    dot.cex=.8
  )
  map("world", add=T)
  dev.copy2pdf(file=paste(OUTP, "PostModeMap-", x, ".pdf", sep="/"))
}
)
