# This is the main file that goes through everything from the genotypes input file
# note that it assumes that the working directory is the one in 
# which the directories "R", "R-main", and "data" reside.


##### LOAD PACKAGES, INSTALL GPIPER, SOURCE SOME FILES  ####
stopifnot(
  library(plyr, logical.return = TRUE),
  library(lubridate, logical.return = TRUE),
  library(maps, logical.return = TRUE),
  library(mapdata, logical.return = TRUE),
  library(devtools, logical.return = TRUE)
)

# get and install gpiper from github.  I put in a reference for it in case 
# we change gpiper in ways that break things here.
if("gpiper" %in% rownames(installed.packages()) == FALSE) {
  install_github(repo = "gpiper", username="eriqande", ref="36be155a125b7cf1975a6e5c1f30bd05fb312fd7")
}
stopifnot(library(gpiper, logical.return = TRUE))

# source some files with functions and other variables
source("./R/wiwa_extra_funcs.R")  # for my floating pies
source("./R/wiwa_colors.R")

##### PREPARE THE WORK AREA.  BASICALLY CREATE AN OUTPUT DIRECTORY  ####
# note that this directory will not be under version control and will be ignored
# by git.
dir.create("outputs", showWarnings = F)

##### READ IN ALL THE DATA AND REMOVE SOME INCOMPLETE CASES AND DUPLICATED NAMES  ####
# this is the "WIWA-All columns"
# read it like this and then remove the duplicated mMXMX06 because there was a naming glitch on the plate
WA <- read.csv("./data/wiwa-snp-genotypes.csv",  stringsAsFactors=F, na.strings=c(""))

# finde those that don't have lat-longs
no_lat_longs_logical <- (is.na(WA$Latitude) | is.na(WA$Longitude))
# and print them out here:
WA$Short_Name[no_lat_longs_logical]

# finally, drop them from the data set
WA <- WA[!no_lat_longs_logical, ]


# find ones with duplicated short names (this seems to be only one...the last one)
duplicated_short_names_logical <- duplicated(WA$Short_Name)
# print them out here
WA$Short_Name[duplicated_short_names_logical]

# and then drop them
WA <- WA[!duplicated_short_names_logical, ]

# now I make the short_names the rownames and then I toss the first column which is the sample names
rownames(WA) <- WA$Short_Name
WA <- WA[,-1]

# we add a column on the end that is the letters of their names:
WA$Pop <- gsub("[0-9]*", "", rownames(WA) )

#### NOW DEFINE THE COLUMNS THAT HAVE THE LOCI AND DISCARD ONES WITH TOO MUCH MISSING DATA  ####
loc.idx <- 1:192  # these are just the indices of the loci in the data frame

# now, we are going to drop the individuals that have more than 6 missing loci:
too_much_missing_logical <- (rowSums(WA[,loc.idx]==0)/2)>6
WA.MD <- WA[!too_much_missing_logical, ]  # WA.MD ==> Wiwa All Missers Dropped

# let us make a note of those individuals that got tossed (had more than 6 missing loci)
TheseGotDropped <- WA[too_much_missing_logical, ] 

# actually we are going to add a column in WA called Was.Included
# and put Yes if it was included and No otherwise
WA$Was.Included <- ifelse(too_much_missing_logical, "No", "Yes")

# now write that out for Kristen
write.table(cbind(ShortName=rownames(WA), Was.Included=WA$Was.Included), file="outputs/WIWA_WhichWasUsed.txt", quote=F, row.names=F, sep="\t")

##### SET ORDER OF POPS, SPLIT INTO BREEDING AND NON-BREEDING, AND FIND LAT-LONG CENTROIDS OF GROUPS ####
# here are the breeding "Pops" in the order that Kristen wants them in.
# and we also define some reporting units here
breed.ord <- c("wAKDE", "wAKYA", "wAKUG", "wAKJU", "wABCA", "wBCMH", "wWADA", 
  "wORHA", "wORMB", "wCAEU", "wCAHM", "wCABS", "wCASL", "wCATE", "wCACL", "wCAHU",
  "wMTHM", "wOREL", "wCOPP", "wCOGM", "eQCCM", "eONHI", "eNBFR"
  )
	
if(TRUE) { 
rep.names <- c("AK.EastBC.AB", "AK.EastBC.AB", "AK.EastBC.AB", "AK.EastBC.AB", "AK.EastBC.AB", "AK.EastBC.AB", "Wa.To.NorCalCoast", 
	"Wa.To.NorCalCoast", "Wa.To.NorCalCoast", "Wa.To.NorCalCoast", "CentCalCoast", "CentCalCoast", "CentCalCoast", "CalSierra", "CalSierra", "CalSierra", 
	"Basin.Rockies", "Basin.Rockies", "Basin.Rockies", "Basin.Rockies", "Eastern", "Eastern", "Eastern")	
}
if(FALSE) {  # these were a different collection tried at one point
	rep.names <- c("Alaska", "Alaska", "Alaska", "Alaska", "EastBC.Alberta", "EastBC.Alberta", "Wa.To.NorCalCoast", 
		"Wa.To.NorCalCoast", "Wa.To.NorCalCoast", "Wa.To.NorCalCoast", "CentCalCoast", "CentCalCoast", "CentCalCoast", "CalSierra", "CalSierra", "CalSierra", 
		"Basin.Rockies", "Basin.Rockies", "Basin.Rockies", "Basin.Rockies", "Eastern", "Eastern", "Eastern")		
}
rep.units <- factor(rep.names, levels=unique(rep.names))



# let's pick out the breeders, and then maybe sort them as desired:
WA.B <- WA.MD[ WA.MD$Pop %in% breed.ord, ] 
WA.B$Pop <- factor(WA.B$Pop, levels=breed.ord)  # make a factor out of them in the specified order
WA.B <- WA.B[order(WA.B$Pop),]  # order them by factor

# now print out a table with the number of breeders by unique lat/long
write.table(
	count(WA.B, c("Pop", "Latitude", "Longitude", "Area_General", "Area_Specific", "State_Province", "Country")), 
	file="outputs/breeding-location-counts.txt", 
	row=F, quote=F, sep="\t"
)

# now, within each "short-name" group/population, we are going to want to get the average Lat-Long
# So I can print a map that has those letters on it, so Kristen knows where to place them:
tmp <- aggregate(cbind(Longitude, Latitude)~Pop, data=WA.B, FUN=mean)
Pop.Centers <- cbind(Letter=LETTERS[1:nrow(tmp)], tmp)
names(Pop.Centers)[3:4] <- c("lon", "lat")


# and then get all the non-breeders too:
WA.WM <- WA.MD[ !(WA.MD$Pop %in% breed.ord), ]   # WA.WM ==> Wiwa-All.Wintering-Migrants

#### SELF-ASSIGNMENT TESTS WITH GSI_SIM  #####

# make the gsi_sim baseline file:
gPdf2gsi.sim(WA.B[,loc.idx], WA.B$Pop, "outputs/baseline.txt")


# do a self-assignment gsi_sim run with just the baseline, and after that read the text file
# output back into a data frame and manipulate it
gsi_Run_gsi_sim(arg.string=c("-b", "outputs/baseline.txt", "--self-assign"), stdout.to = "outputs/self-ass-output.txt")
SA.df <- gsi.simSelfAss2DF("outputs/self-ass-output.txt")$Post
SA.to.Rep <- gsi_aggScoresByRepUnit(SA.df, breed.ord, rep.units)
SA.to.Rep <- cbind(SA.to.Rep, gsi_simMaxColumnAndPost(SA.to.Rep, 3:ncol(SA.to.Rep)) )
SA.rep.cuts <- gsi_simAssTableVariousCutoffs(SA.to.Rep$PopulationOfOrigin, SA.to.Rep$MaxColumn, SA.to.Rep$MaxPosterior)
# I sent the above as text to Kristen.

#### NOW ASSIGN INDIVIDUALS FROM THE WINTERING AND MIGRATORY SAMPLES VIA GSI_SIM  ####

# make the gsi-sim mixture and reporting units files
gPdf2gsi.sim(WA.WM[,loc.idx], outfile="outputs/mixture.txt")
gsi_WriteReportingUnitsFile(breed.ord, rep.units, "outputs/repunits.txt")

# change into the outputs directory so we don't bomb the cwd with all the gsi_sim output:
basedir <- getwd()
setwd("outputs")

# now run gsi-sim. Note that I could add more reps for the z-scores. I just want this to be fast right now
system(paste(gsi_simBinaryPath(), " -b baseline.txt -t mixture.txt -r repunits.txt --mix-logl-sims 10 0"), ignore.stdout=T)

# then change back to the old working directory (get out out outputs)
setwd(basedir)

# now we capture all the wintering and migrant gsi-sim results in a single data frame (dropping the genetic data for simplicity)
WM.gr <- cbind(WA.WM[, -loc.idx], read.table("outputs/rep_unit_pofz_full_em_mle.txt", header=T), read.table("outputs/em_mixture_logl_summary.txt", header=T))
WM.gr$AssignedTo <- factor(WM.gr$AssignedTo, levels=breed.ord)  # put these factors in the right order


# and I give that data frame to kristen so she can peruse it.
write.table(WM.gr, "outputs/GsiResults1.txt", quote=F, sep="\t", row.names=T, col.names=NA)

# now, find the highest posterior max repu in each case:
WM.gr$MaxRepu <- factor(levels(rep.units)[max.col(WM.gr[,levels(rep.units)])], levels=levels(rep.units))
# and add the max posterior as well:
WM.gr$PostMax <- apply(WM.gr[,levels(rep.units)], 1, max)

#### PLOT THE DISTRIBUTION OF POSTERIOR PROBS FOR THE ASSIGNMENTS ####
# check out the distribution of posterior prabilities to different groups
pdf(width = 6, height = 6, file = "outputs/posterior_prob_boxplots.pdf")
par(mar=c(5, 4, 4, 2)+c(9,.1,.1,.1))
plot(WM.gr$MaxRepu, WM.gr$PostMax, las=2)
dev.off()


#### SPLIT RESULTS INTO WW (WINTERING BIRDS) AND MM (MIGRATING BIRDS) and write out to files for later #####
# OK, this is cool. Now we want to break these into wintering and migrating:
# The migrants either are from Pop=="WIWA" or have Pop starting with "m"
migrant.idx <- WM.gr$Pop=="WIWA" | grepl("^m", WM.gr$Pop)
WW <- WM.gr[ !migrant.idx, ]  ## WW =  Wintering Wintering!!
MM <- WM.gr[  migrant.idx, ]  ## MM = Migrating Migrating !!

write.csv(MM, file="outputs/assignment-results-migrants.csv")
write.csv(WW, file="outputs/assignment-results-wintering.csv")


#### ANALYZE THE WINTERING BIRDS GSI RESULTS (WW) ####
# here are some fixes that I have to make.  Perhaps Kristen can update the data base.
# one of the baja birds does not have an entry for "Area_General"
WW$Area_General[WW$Latitude==22.883 & is.na(WW$Area_General)] <- "San Jose del Cabo"


# let us just have a look at where these all are from:
count(WW, c("Pop", "Latitude", "Longitude", "Area_General", "Area_Specific", "State_Province", "Country"))


# OK, now we need to aggregate some of the wintering birds into groups tied to a cluster
# of nearby lat-longs, which is obvious from the above.  Here is what we do:
# 1. write out a table of all the unique lat-long and locations
write.table(count(WW, c("Pop", "Latitude", "Longitude", "Area_General", "Area_Specific", "State_Province", "Country")), 
            quote=F, sep="\t", file="outputs/winter_locations.txt")
# 2. Then eric and kristen by hand lumped those unique combinations into different groups. That
# file is read back here
winter.lumps <- read.csv("./data/eric_lumps_winter_locations.csv", stringsAsFactors=F)
# 3. Then we put another column into WW called GeoGroup that tells which Geographic group each wintering bird is in
wilmps <- winter.lumps$Wintering.Group
names(wilmps) <- winter.lumps$Area_General
WW$GeoGroup <- wilmps[WW$Area_General]


# now, let us write those counts out for the table.  # this really needs to be more automated in the
# future!!  I should see what Kristen wants and then generate it in R.
LumpCounts <- do.call(rbind, lapply(split(WW, WW$GeoGroup), function(x) table(x$MaxRepu)))
write.csv(LumpCounts, file="outputs/LumpCounts.csv")

# and we also might want the assignments to every specific lat-long.  There must be a better way to do this:
ugly <- split(WW, paste(WW$Latitude, WW$Longitude, sep="---"))
LatLongCounts <- do.call(rbind, lapply(ugly[sapply(ugly, nrow)>0], function(x) table(x$MaxRepu)))

# and then we create a variable that we will write out for making the pretty map in script 04
WW.list <- split(WW, WW$GeoGroup)
tmp <- lapply(WW.list, 
              function(x) {
                lat <- mean(x$Latitude)
                long <- mean(x$Longitude)
                Locats <- paste(unique(x$Area_General), collapse="/") 
                tab <- t(data.frame((table(x$MaxRepu)[])))
                cnt <- nrow(x)
                cbind(tab, n=cnt, data.frame(lat=lat, long=long, Locats=Locats))
              }
)
WW.assign.df <- do.call(rbind, tmp)


#### PLOT THE WINTERING BIRDS RESULTS IN VARIOUS FORMATS  ####
WW.ass.list <- lapply(split(WW$MaxRepu, WW$GeoGroup), sort)  # assignments of birds in each GeoGroup, sorted by reporting unit

# Here, we are going to represent each bird as a dot of a certain color in a grid like so:
# We actually ended up not using it...
dot.squares <- function(x=WW.ass.list)  {
	w <- sapply(x, function(z) ceiling(sqrt(length(z))))  # how wide each square will be (in number of dots)
	
	xs <- c(0,cumsum(w))  # x starting values for each square
	names(xs) <- names(w)  # shift names down
	xs <- xs[-length(xs)]  # drop the last element which is irrelevant
	
	ys <- cumsum(w)  # starting y values for each square
	
	# now, the position of each dot in each square can be determined:
	dp <- lapply(names(w), function(z) {
		xx<-xs[z]
		yy<-ys[z]
		ww<-w[z]
		cbind(X=xx+(0:(length(x[[z]])-1) %% ww) , Y=yy-(0:(length(x[[z]])-1) %/% ww) )
		}
	)
	dp.mat <- do.call(rbind, dp)
	d <- data.frame(X=dp.mat[,"X"], Y=dp.mat[,"Y"], Color=pie.cols[unlist(sapply(x, as.numeric))], stringsAsFactors=F)
	plot(d$X, d$Y, col=d$Color, pch=19)
}

pdf(width = 9, height = 6, file = "outputs/dot-squares.pdf")
dot.squares()
dev.off()

# Here we make pie charts that Kristen ended up putting onto the map in photoshop
# or illustrator.  It allowed us to have them vectorized while sitting atop a 
# raster-based map.  To get them right, you would have to open up a graphics window
# to the desired size. 
pdf(width = 9, height = 6, file = "outputs/wintering-pies-with-sample-sizes.pdf")
WW.ass.tab <- lapply(WW.ass.list, table)
par(mar=c(0,0,0,0))
# make a blank space
plot(c(-1,2.7*length(WW.ass.tab)), c(-1.2, 2.2), type="n")
xx<-0
lapply(names(WW.ass.tab), function(z) {
	xx<<-xx+2.5; 
	x <- WW.ass.tab[[z]]
	my.floating.pie(xx, 0, x, col=pie.cols)
	text(xx, 1.0, label=z, srt=90, adj=0)
	text(xx-.2, 0, label=sum(x), adj=1, cex=.8)
	}
)
dev.off()

#### NOW, ANALYZE THE MIGRANT BIRDS RESULTS  ####
# now, we analyze the migrants:
MM$date <- dmy(MM$Collection_Date)
MM$Year <- year(MM$date)

# For general interest, count up assignments in different weeks, in different years
wk.tab <- table(MM$MaxRepu, week(MM$date), year(MM$date), MM$Area_Specific)  

# make a summary table for kristen:
write.table(count(MM, c("Pop", "Latitude", "Longitude", "Area_General", "Area_Specific", "State_Province", "Country")),
	quote=F, sep="\t", file="outputs/migrant_locations.txt")

write.table(count(MM, c("Pop", "Latitude", "Longitude", "Area_General", "Area_Specific", "State_Province", "Country", "Year")),
	quote=F, sep="\t", file="outputs/migrant_locations_by_year_too.txt")


# count up number of birds going through each Area_Specific across all weeks and years:
AreaByRepuTab <- table(MM$Area_Specific, MM$MaxRepu)
AreaByRepuTabProportions <- apply(AreaByRepuTab, 1, function(x) x/sum(x))
#dev.copy2pdf(file="outputs/area-specific-naked-arrows.pdf")


# now, we want to put the Cibola birds in there in separate years (2008 and 2009)
# but all the other ones we just lump all the different years
MMSplit <- split(MM, MM$Area_Specific)  # split into a bunch of data frames by location
CiboByYears <- split(MMSplit[["Cibola"]], year(MMSplit[["Cibola"]]$date)) # this is splitting Cibola by years
names(CiboByYears) <- paste("Cibola", names(CiboByYears), sep=".") # here we name them with their years


# now we make a list of only the data frames that we want.  We drop the eastern and really small samples
MMSplit <- c(CiboByYears, MMSplit[c(3,6,7,9)])

# confirm that there are only 4 reporting units represented amongst these:
lapply(MMSplit, function(x) levels(droplevels(x$MaxRepu)))

# table up the time series of these for each week in the fall migration:
WeekTabs <- lapply(MMSplit, function(x) {
		mrep <- factor(x$MaxRepu, levels=c("AK.EastBC.AB", "Wa.To.NorCalCoast", "CentCalCoast", "CalSierra"))
		table(mrep, factor(week(x$date), levels=11:21))
	}
)

#### PLOT THE MIGRANT BIRD RESULTS ####
pdf(width = 8, height = 9, file = "outputs/week-tables-six-spots.pdf")
# now plot those
# make a separate plot for each location:

# drop some and make them north south
WeekTabsLight<-WeekTabs #[c(4,6,2,1)]
counter<<-0
par(mar=c(.4,3.2,1.4,.1), las=1, mfrow=c(length(WeekTabsLight),1), oma=c(5,4,0,0))
lapply(names(WeekTabsLight), function(x) {
	counter<<-counter + 1
	y<-WeekTabs[[x]]; 
	if(counter<length(WeekTabsLight))  nameys <- rep("", ncol(y))
	else nameys <- colnames(y)
	barplot(y, beside=T, col=pie.cols[1:4], , xlab="", names.arg=nameys);
	text(0, 0.75*max(y), x, pos=4) 
	#legend("top", legend=rownames(y), fill=pie.cols); 
	#z <- gsub(",-", "-", gsub("  *", "-", x));
	#dev.copy2pdf(file=paste("mig-barplot-",z,".pdf", sep="")) 
	}
)
mtext("Week of the year", side=1, line=2.7, outer=T)
mtext("Birds encountered per week", side=2, line=4, las=0, adj=-2)
dev.off()


# OK, NOW KRISTEN WANTS A FIGURE WHERE (a) is two panels (Cibola 2008 and then 2009) and
# (b) is the Oneill Forebay. I will make those separately
pdf(width = 10, height = 6, file = "outputs/cibola-bar-plot.pdf")
CB <- WeekTabs[1:2]  # this is the Cib08 and 09
par(mar=c(.4,1.9,1.4,.1), las=1, mfrow=c(length(CB),1), oma=c(5,4,0,0))
lapply(names(CB), function(x) {
	y<-WeekTabs[[x]]; 
	barplot(y, beside=T, col=pie.cols[1:4], xlab="", names.arg=rep("",11), space=c(0,0));
	axis(side=1, at=seq(1, length.out=ncol(y), by=4), labels=F)
	}
)
# at the bottom here we put the annotation:
week.centers <- mdy("1-1-2009") + weeks(as.numeric(colnames(CB[[1]]))) + days(3)
week.cent.str <- paste(month(week.centers, label=T, abbr=T), mday(week.centers), sep="-")
axis(side=1, at=2+seq(1, length.out=11, by=4), tick=F, labels=week.cent.str)
mtext(text="Date of Week Midpoint", side=1, line=2.6, cex=1.3)
mtext(text="Number of Birds", side=2, cex=1.3, las=3, adj=.45, line=.3, outer=T)
dev.off()



# Now, the Oneill Forebay stuff:
pdf(width = 8, height = 4, file = "outputs/oneill-bar-plot.pdf")
Oneill <- MMSplit$"O'neill Forbay Wildlife Area, Merced County, CA"
OneTab <- table(factor(Oneill$MaxRepu, levels=c("AK.EastBC.AB", "Wa.To.NorCalCoast", "CentCalCoast", "CalSierra")), factor(week(Oneill$date), 13:40))
par(mar=c(.4,1.9,1.4,.1), las=1, oma	=c(5,3,0,0))
barplot(OneTab, beside=T, col=pie.cols[1:4], xlab="", names.arg=rep("",ncol(OneTab)), space=c(0,0));
axis(side=1, at=seq(1, length.out=ncol(OneTab), by=4), labels=F)

# at the bottom here we put the annotation:
week.centers <- mdy("1-1-2009") + weeks(as.numeric(colnames(OneTab))) + days(3)
week.cent.str <- paste(month(week.centers, label=T, abbr=T), mday(week.centers), sep="-")
axis(side=1, at=2+seq(1, length.out=ncol(OneTab), by=4), tick=F, labels=week.cent.str, las=3)
mtext(text="Date of Week Midpoint", side=1, line=4.3, cex=1.3)
mtext(text="Number of Birds", side=2, cex=1.3, las=3, adj=.38, line=.3, outer=T)
dev.off()



# migrant pie figure
# quartz(width=9.5, height=4.6)
pdf(width = 9.5, height = 4.6, file = "outputs/migrant-two-year-and-onfb-pies.pdf")
par(mfrow=c(1,1))
par(mar=c(0,0,0,0))
yspot <- .57
radi <- .4
plot(0:12,seq(-1,1,length.out=13), type="n", axes=F, xlab="", ylab="", asp=1)  # make an empty space
for(i in colnames(WeekTabs$Cibola.2008)) {
	x <- as.numeric(i) - 10
	vec08 <- WeekTabs$Cibola.2008[, i]
	vec09 <- WeekTabs$Cibola.2009[, i]
	onfb <-  WeekTabs$"O'neill Forbay Wildlife Area, Merced County, CA"[, i]
	if(sum(vec08)>0) {
		my.floating.pie(xpos=x, ypos=yspot, x=vec08, edges=200, radius=radi, col=pie.cols, startpos=0, shadow=FALSE)
		text(x+radi*.19, yspot+radi*0.05, sum(vec08), pos=2)
	}
	if(sum(vec09)>0) {
		my.floating.pie(xpos=x, ypos=-yspot, x=vec09, edges=200, radius=radi, col=pie.cols, startpos=0, shadow=FALSE)
		text(x+radi*.19, -yspot+radi*0.05, sum(vec09), pos=2)
	}
	if(sum(onfb)>0) {
		my.floating.pie(xpos=x, ypos=-yspot-(2*yspot), x=onfb, edges=200, radius=radi, col=pie.cols, startpos=0, shadow=FALSE)
		text(x+radi*.19, -yspot-(2*yspot)+radi*0.05, sum(onfb), pos=2)
	}
}

# line below was for week numbers along the top, but we don't need them, so it is commented out
#text(as.numeric(colnames(wk.tab[,,"2008"]))-10, yspot+radi*1.3, labels=colnames(wk.tab[,,"2008"]), pos=3, cex=1.3)
text(x=c(1,1)-radi*1.2, y=c(1,-1)*yspot, labels=c("2008", "2009"), pos=2, cex=1.45)  # years
# now put the week-midpoint dates on:
week.centers <- mdy("1-1-2009") + weeks(11:21) + days(3)
week.cent.str <- paste(month(week.centers, label=T, abbr=T), mday(week.centers), sep="-")
text(x=as.numeric(colnames(WeekTabs$Cibola.2008))-10, y=-yspot-2*yspot-radi*1.73, labels=week.cent.str, pos=1, cex=1.1) 
dev.off()


#### MAKE A TABLE OF MIGRANTS IN ARIZONA  ####
# now we are going to put the same information in a big ol' table:
az.tab <- rbind(cbind(year=2008, t(wk.tab[,,"2008", "Cibola"])), cbind(year=2009, t(wk.tab[,,"2009", "Cibola"])))
az.tab <- cbind(week.num=as.numeric(rownames(az.tab)), az.tab)
rownames(az.tab)<-NULL
az.tab <- data.frame(az.tab)
wc <- mdy(paste("1-1",az.tab$year, sep="-")) + weeks(az.tab$week.num) + days(3)
az.tab <- cbind(week.centers=paste(month(wc, label=T, abbr=T), mday(wc), sep="-"), az.tab)
write.table(az.tab, file="outputs/arizona-mirants-pie-ingredients.txt", row=F, quote=F, sep="\t")

#### PREPARE A BUNCH OF STUFF FOR THE GENELAND ANALYSIS  ####
# now, let's prep stuff for geneland.
gl.coord <- WA.B[, c("Longitude", "Latitude")]  # coordinates for geneland
gl.geno <- WA.B[, loc.idx]  # genetic data for geneland
gl.geno[gl.geno==0] <- NA 	# then the 0s into bona-fide missing data


# here we put two fake birds in the geneland data set at the upper right and 
# lower left corners of things (so it will infer across the whole range.). We
# let these birds be missing at all loci except one which is diagnostic for 
# east vs west, and hence essentially known for those two "fake" birds.
# To make the birds that are all missing data except one locus to put on the extreme ends of
# the map, I want to find a locus that is fixed in east v west.  A good locus for that is East_West_01.
# it is allele 3 in eastern pops and 1 in western pops
fakies<-data.frame(matrix(NA, nrow=2, ncol=length(loc.idx)))
names(fakies) <- names(gl.geno)
rownames(fakies) <- c("wFake", "eFake")
fakies[1, c("East_West_01.1", "East_West_01.2")] <- 1
fakies[2, c("East_West_01.1", "East_West_01.2")] <- 3

gl.geno <- rbind(gl.geno, fakies)

# and the positions we want to put them at are:
# The west one should be at (-170, 72) and I think
# that the fake east one would be fine at (-50, 47).
gl.coord["wFake", ] <- c(-170, 72)
gl.coord["eFake", ] <- c(-50, 47)

#### AND THEN LET US PRINT SOME DATA OUT IN PIPELINE FORMAT TO RUN STRUCTURE ON IT ####
slg <- WA.B[, loc.idx]
names(slg)<-gsub("\\.[12]$", "", names(slg))
write.table(slg, file="outputs/wibreed-pipe-genos.txt", quote=F, col.names=NA, sep="\t")
write(breed.ord, "outputs/wibreed-pipe-pops.txt")
write(names(slg)[c(T,F)], "outputs/wibreed-pipe-locs.txt")

#### FINALLY, SAVE SOME VARIABLES IN AN RDA FILE FOR LATER USE ####
save(Pop.Centers, rep.units, WA.B, gl.geno, gl.coord, WM.gr, WW.assign.df,file="outputs/WIWA-main-carryover-variables.Rda")





