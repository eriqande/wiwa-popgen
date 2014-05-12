# This is the main file that goes through everything from our standardized 
# data set forward.
library(plyr)
library(lubridate)
library(maps)
library(mapdata)
#library(gpiper)
library(devtools)

setwd("/Users/eriq/Documents/work/assist/kristenruegg/WIWA_popgen/WIWA_popgen_analysis/")

source("./R/wiwa_extra_funcs.R")  # for my floating pies

load_all("/Users/eriq/Documents/git-repos/gpiper")  # this loads the gpiper library (which is still under development)



# source some files with functions and other variables
source("./R/wiwa_colors.R")


# this is the "WIWA-All columns"
#WA <- read.csv("./data/WIWA_SNP_and_metadata_Oct.csv", row=1, stringsAsFactors=F)
# read it like this and then remove the duplicated mMXMX06 because there was a naming glitch on the plate
WA <- read.csv("./data/WIWA_SNP_Jan2014.mer",  stringsAsFactors=F, na.strings=c(""))

# drop those that don't have lat-longs
WA <- WA[!(is.na(WA$Latitude) | is.na(WA$Longitude)),]

# toss ones with duplicated short names (this seems to be only one...the last one)
WA <- WA[!duplicated(WA$Short_Name),]

# now I make the short_names the rownames and then I toss the first column which is the sample names
rownames(WA) <- WA$Short_Name
WA <- WA[,-1]

# we add a column on the end that is the letters of their names:
WA$Pop <- gsub("[0-9]*", "", rownames(WA) )

loc.idx <- 1:192  # these are just the indices of the loci in the data frame


# now, we are going to drop the individuals that have more than 6 missing loci:
WA.MD <- WA[!(rowSums(WA[,loc.idx]==0)/2)>6, ]  # WA.MD ==> Wiwa All Missers Dropped

# let us make a note of those individuals that got tossed (had more than 6 missing loci)
TheseGotDropped <- WA[(rowSums(WA[,loc.idx]==0)/2)>6, ] 

# actually we are going to add a column in WA called Was.Included
# and put Yes if it was included and No otherwise
WA$Was.Included <- ifelse( (rowSums(WA[,loc.idx]==0)/2)>6, "No", "Yes")

# now write that out for Kristen
write.table(cbind(ShortName=rownames(WA), Was.Included=WA$Was.Included), file="WIWA_WhichWasUsed.txt", quote=F, row.names=F, sep="\t")


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
	file="breeding-location-counts.txt", 
	row=F, quote=F, sep="\t"
)

# now, within each "short-name" group/population, we are going to want to get the average Lat-Long
# So I can print a map that has those letters on it, so Kristen knows where to place them:
tmp <- aggregate(cbind(Longitude, Latitude)~Pop, data=WA.B, FUN=mean)
Pop.Centers <- cbind(Letter=LETTERS[1:nrow(tmp)], tmp)
names(Pop.Centers)[3:4] <- c("lon", "lat")


# and then get all the non-breeders too:
WA.WM <- WA.MD[ !(WA.MD$Pop %in% breed.ord), ]   # WA.WM ==> Wiwa-All.Wintering-Migrants




# make the gsi_sim files:
gPdf2gsi.sim(WA.B[,loc.idx], WA.B$Pop, "baseline.txt")
gPdf2gsi.sim(WA.WM[,loc.idx], outfile="mixture.txt")



# now, write out a reporting units file too!
gsi_WriteReportingUnitsFile(breed.ord, rep.units, "repunits.txt")


#### SELF-ASSIGNMENT tests with GSI_SIM
# do a self-assignment gsi_sim run with just the baseline, and after that read the text file
# output back into a data frame and manipulate it
gsi_Run_gsi_sim(arg.string=c("-b", "baseline.txt", "--self-assign"), stdout.to = "self-ass-output.txt")
SA.df <- gsi.simSelfAss2DF("self-ass-output.txt")$Post
SA.to.Rep <- gsi_aggScoresByRepUnit(SA.df, breed.ord, rep.units)
SA.to.Rep <- cbind(SA.to.Rep, gsi_simMaxColumnAndPost(SA.to.Rep, 3:ncol(SA.to.Rep)) )
SA.rep.cuts <- gsi_simAssTableVariousCutoffs(SA.to.Rep$PopulationOfOrigin, SA.to.Rep$MaxColumn, SA.to.Rep$MaxPosterior)
# I sent the above as text to Kristen.

# now run gsi-sim. Note that I could add more reps for the z-scores. I just want this to be fast right now
system(paste(gsi_simBinaryPath(), " -b baseline.txt -t mixture.txt -r repunits.txt --mix-logl-sims 10 0"), ignore.stdout=T)


# now we capture all the migrant gsi-sim results in a single data frame (dropping the genetic data for simplicity)
WM.gr <- cbind(WA.WM[, -loc.idx], read.table("rep_unit_pofz_full_em_mle.txt", header=T), read.table("em_mixture_logl_summary.txt", header=T))
WM.gr$AssignedTo <- factor(WM.gr$AssignedTo, levels=breed.ord)  # put these factors in the right order


# and I give that data frame to kristen so she can peruse it.
write.table(WM.gr, "GsiResults1.txt", quote=F, sep="\t", row.names=T, col.names=NA)

# now, find the highest posterior max repu in each case:
WM.gr$MaxRepu <- factor(levels(rep.units)[max.col(WM.gr[,levels(rep.units)])], levels=levels(rep.units))
# and add the max posterior as well:
WM.gr$PostMax <- apply(WM.gr[,levels(rep.units)], 1, max)




# check out the distribution of posterior prabilities to different groups
par(mar=c(5, 4, 4, 2)+c(9,.1,.1,.1))
plot(WM.gr$MaxRepu, WM.gr$PostMax, las=2)


# OK, this is cool. Now we want to break these into wintering and migrating:
# The migrants either are from Pop=="WIWA" or have Pop starting with "m"
migrant.idx <- WM.gr$Pop=="WIWA" | grepl("^m", WM.gr$Pop)
WW <- WM.gr[ !migrant.idx, ]  ## WW =  Wintering Wintering!!
MM <- WM.gr[  migrant.idx, ]  ## MM = Migrating Migrating !!


# here are some fixes that I have to make.  Perhaps Kristen can update the data base.
# one of the baja birds does not have an entry for "Area_General"
WW$Area_General[WW$Latitude==22.883 & is.na(WW$Area_General)] <- "San Jose del Cabo"


# let us just have a look at where these all are from:
count(WW, c("Pop", "Latitude", "Longitude", "Area_General", "Area_Specific", "State_Province", "Country"))



# OK, now we need to aggregate some of the wintering birds into groups tied to a cluster
# of nearby lat-longs, which is obvious from the above.  Here is what we do:
# 1. write out a table of all the unique lat-long and locations
write.table(count(WW, c("Pop", "Latitude", "Longitude", "Area_General", "Area_Specific", "State_Province", "Country")), 
	quote=F, sep="\t", file="winter_locations.txt")
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
# and we also want the assignments to every specific lat-long.  There must be a better way to do this:
ugly <- split(WW, paste(WW$Latitude, WW$Longitude, sep="---"))
LatLongCounts <- do.call(rbind, lapply(ugly[sapply(ugly, nrow)>0], function(x) table(x$MaxRepu)))
write.csv(LumpCounts, file="LumpCounts.csv")



# now, we are going to represent each bird as a dot of a certain color in a grid like so:
WW.ass.list <- lapply(split(WW$MaxRepu, WW$GeoGroup), sort)  # assignments of birds in each GeoGroup, sorted by reporting unit
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


dot.squares()


###################### HERE WE MAKE SOME PIES:
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
dev.copy2pdf(file="wintering-pies-with-sample-sizes.pdf")

############################################################################################################



################ DOING THE MIGRANT BIRDS  #####################
# now, we analyze the migrants:
MM$date <- dmy(MM$Collection_Date)
MM$Year <- year(MM$date)

# For general interest, count up assignments in different weeks, in different years
wk.tab <- table(MM$MaxRepu, week(MM$date), year(MM$date), MM$Area_Specific)  

# make a summary table for kristen:
write.table(count(MM, c("Pop", "Latitude", "Longitude", "Area_General", "Area_Specific", "State_Province", "Country")),
	quote=F, sep="\t", file="migrant_locations.txt")

write.table(count(MM, c("Pop", "Latitude", "Longitude", "Area_General", "Area_Specific", "State_Province", "Country", "Year")),
	quote=F, sep="\t", file="migrant_locations_by_year_too.txt")


# count up number of birds going through each Area_Specific across all weeks and years:
AreaByRepuTab <- table(MM$Area_Specific, MM$MaxRepu)
AreaByRepuTabProportions <- apply(AreaByRepuTab, 1, function(x) x/sum(x))
dev.copy2pdf(file="area-specific-naked-arrows.pdf")


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


if(FALSE) {
		################# working up a stacked area thing here:
		library(ggplot2)
		Cib <- MMSplit$Cibola.2009
		Cib$Week <- week(Cib$date)
		Cib$MaxRepu <- factor(Cib$MaxRepu, levels=levels(MM$MaxRepu)[c(4,3,2,1)]) # drop missing levels and order by rarity

		CibTab.df <- count(Cib, c("MaxRepu", "Week"))

		# now make the plot
		p <- ggplot(CibTab.df, aes( Week, freq)) + scale_color_manual(values = pie.cols[c(4,3,2,1)]) + scale_fill_manual(values = pie.cols[c(4,3,2,1)])
		p + geom_area(aes(colour = MaxRepu, fill= MaxRepu), position = 'stack') 


		# now try doing both years together as facets and doing it by day of the year.  This is pretty cool.
		Cib <- MM[MM$Area_Specific=="Cibola",]
		Cib$Week <- week(Cib$date)
		Cib$DayOfYear <- yday(Cib$date)
		Cib$Year <- year(Cib$date)
		Cib$MaxRepu <- factor(Cib$MaxRepu, levels=levels(MM$MaxRepu)[c(4,3,2,1)]) # drop missing levels and order by rarity

		CibTab.df <- count(Cib, c("MaxRepu", "DayOfYear", "Year"))

		# now make the plot
		p <- ggplot(CibTab.df, aes( DayOfYear, freq)) + 
					scale_color_manual(values = pie.cols[c(4,3,2,1)]) + 
					scale_fill_manual(values = pie.cols[c(4,3,2,1)]) + 
					geom_area(aes(colour = MaxRepu, fill= MaxRepu), position = 'stack') +
					facet_grid(Year ~ .)

		p
}



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

dev.copy2pdf(file="week-tables-six-spots.pdf")


# OK, NOW KRISTEN WANTS A FIGURE WHERE (a) is two panels (Cibola 2008 and then 2009) and
# (b) is the Oneill Forebay. I will make those separately
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
dev.copy2pdf(file="cibola-bar-plot.pdf")



# Now, the Oneill Forebay stuff:
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
dev.copy2pdf(file="oneill-bar-plot.pdf")



# migrant pie figure
quartz(width=9.5, height=4.6)
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

dev.copy2pdf(file="migrant-two-year-and-onfb-pies.pdf")


# now we are going to put the same information in a big ol' table:
az.tab <- rbind(cbind(year=2008, t(wk.tab[,,"2008"])), cbind(year=2009, t(wk.tab[,,"2009"])))
az.tab <- cbind(week.num=as.numeric(rownames(az.tab)), az.tab)
rownames(az.tab)<-NULL
az.tab <- data.frame(az.tab)
wc <- mdy(paste("1-1",az.tab$year, sep="-")) + weeks(az.tab$week.num) + days(3)
az.tab <- cbind(week.centers=paste(month(wc, label=T, abbr=T), mday(wc), sep="-"), az.tab)
write.table(az.tab, file="arizona-mirants-pie-ingredients.txt", row=F, quote=F, sep="\t")

###### NOW, WE PREPARE A BUNCH OF STUFF FOR THE GENELAND ANALYSIS
# now, let's prep stuff for geneland.
gl.coord <- WA.B[, c("Long", "Lat")]  # coordinates for geneland
gl.geno <- WA.B[, loc.idx]  # genetic data for geneland
gl.geno[gl.geno==0] <- NA 	# then the 0s into bona-fide missing data

# now, to make some birds that are all missing data except one locus to put on the extreme ends of
# the map, I want to find a locus that is fixed in east v west.  A good locus for that is East_West_01.
# it is allele 3 in eastern pops and 1 in western pops
fakies<-data.frame(matrix(NA, nrow=2, ncol=length(loc.idx)))
names(fakies) <- names(gl.geno)
rownames(fakies) <- c("wFake", "eFake")
fakies[1, c("East_West_01", "East_West_01.1")] <- 1
fakies[2, c("East_West_01", "East_West_01.1")] <- 3

gl.geno <- rbind(gl.geno, fakies)

# and the positions we want to put them at are:
# The west one should be at (-170, 72) and I think
# that the fake east one would be fine at (-50, 47).
gl.coord["wFake", ] <- c(-170, 72)
gl.coord["eFake", ] <- c(-50, 47)


##### AND THEN LET US PRINT SOME DATA OUT IN PIPELINE FORMAT TO RUN STRUCTURE ON IT
slg <- WA.B[, loc.idx]
names(slg)<-gsub("\\.1$", "", names(slg))
write.table(slg, file="wibreed-pipe-genos.txt", quote=F, col.names=NA, sep="\t")
write(breed.ord, "wibreed-pipe-pops.txt")
write(names(slg)[c(T,F)], "wibreed-pipe-locs.txt")




save(Pop.Centers, rep.units, WW.assign.df, WA.B, gl.geno, gl.coord, AssLatLong, WM.gr, file="WIWA-main-carryover-variables.Rda")