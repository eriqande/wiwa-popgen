# this should be run on a machine with multiple cores to make it go faster

#### Load libraries and get data ####
stopifnot(
  library(parallel, logical.return = TRUE),
  library(Geneland, logical.return = TRUE)
)


# get the data.  they are in gl.coord and gl.geno after you load the following
got.these <- load("outputs/WIWA-main-carryover-variables.Rda")

if(!all(c("gl.geno", "gl.coord") %in% got.these)) {
	stop("Failed to successfully load variables gl.geno and/or gl.coord from file WIWA-main-carryover-variables.Rda")
}


#### DEFINE A FUNCTION TO DO THE GENELAND RUNS ####
# this function creates the directory subdir-suffix, changes to and
# then runs geneland there.  This is a wrapper to let you run it
# using multicore.  You have to define 
runGeneLand <- function(subdir="GL-run", suffix="0", geno, coord) {
	
	curd <- getwd()
	dn <- paste(subdir, suffix, sep="-")
	dir.create(dn)
	setwd(dn)

	# then run the MCMC function
	# these appear to be largely the values in kristen parameter file
	MCMC(coordinates=coord, 
		geno.dip.codom=geno,
		varnpop=TRUE, 
		npopmax=10, 
		spatial=TRUE, 
		freq.model="Uncorrelated", 
		nit=2.2e6, 
		thinning=100, 
		path.mcmc="./",
		filter.null.alleles = FALSE
		)
	
	
	
	# then post-process (# note they have a bug in the manual, you don't put geno.dip.codom=geno in there!)
	PostProcessChain(
		coordinates=coord, 
		path.mcmc="./", 
		nxdom=250, 
		nydom=250, 
		burnin=0.5e4  # this is the number of thinned samples. Toss out the first quarter (rouhgly) of them
		)
		
		setwd(curd)
}


#### DO THE GENELAND RUNS  ####
# and here we spawn ourselves 10 processes running geneland
mclapply(1:10, function(x) runGeneLand(
		subdir="./outputs/GeneLandRun", 
		suffix=x, 
		geno=gl.geno,
		coord=gl.coord),
	mc.cores = 10
)

