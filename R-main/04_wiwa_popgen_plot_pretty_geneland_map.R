# set to true if you want to reproject the base-map and/or regenerated the 
# rasters of the polygons giving the breeding and the wintering range of 
# wiwa's.  (If you do that, you will have to get the basemap and the 
# shapefiles for the bird ranges.) 
REGENERATE_BASE_MAP <- FALSE
REGENERATE_POLY_RASTS <- FALSE


#### LOAD LIBRARIES, DEFINE SMALL HELPER FUNCTIONS, GET DATA, SOURCE NECESSARY FILES  ####
stopifnot(
  library(raster, logical.return = TRUE),
  library(rgdal, logical.return = TRUE),
  library(geosphere, logical.return = TRUE),
  library(RCurl, logical.return = TRUE)
)


# here is a function to download stuff via RCurl
bdown=function(url, file){
  f = CFILE(file, mode="wb")
  a = curlPerform(url = url, writedata = f@ref, noprogress=FALSE)
  close(f)
  return(a)
}


# get some of the results from wiwa_analysis_main.R
load("outputs/WIWA-main-carryover-variables.Rda")

source("./R/wiwa_colors.R")
source("./R/wiwa_extra_funcs.R")


# define here the proba.pop.membership file you wish to use:
ProbAFile <- "./intermediates/proba.pop.membership.txt"
# and tell it how to permute the colors, which is specific to this particular run 
# of geneland:
col.perm <- c(2, 6, 1, 5, 3, 4)

# now, permute the colors:
wiwa.cols <- wiwa.cols[col.perm]




# For a projection
# I ended up using Lamberts Conformal Conic as follows to allow us to snip a square part out of the map with
# illustrator.
wiwa.crs <- CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-110 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")




if(REGENERATE_BASE_MAP==TRUE) {
  
  ## NOTE: IF YOU WANT TO REGENERATE THE BASEMAP YOU WILL HAVE TO DOWNLOAD THE BASE MAP
  ## FROM NATURAL EARTH DATA.  
	b<-brick("~/Documents/NaturalEarthData/HYP_LR_SR_W_DR/HYP_LR_SR_W_DR.tif")
		
	b1 <- crop(b, extent(-180, -5, -8, 90))  # this is how we crop it to give us a huge area so we can cut out
																					# a rectangular part of it after projecting.  There must be a more 
																					# efficient way to do that...

	# if it were desired to downsample the resolution to make it easier to handle the files, then it
	# would be possible with the following line (and comment out b2 <- b1 in that case)
	b2 <- aggregate(b1, 8)
	#b2 <- b1  # or don't downsample it
	
	
	# For historical interest, here was how I looked at the format for a few different projections
	# like the Albers Equal Area....
	#EPSG <- make_EPSG()  # load the data frame of various used projections
	#EPSG[grep("Albers", EPSG$note), ]  # find those that mention Albers
	#CRS("+init=epsg:3309")  # this appears to be a california centric one, which might be good
	
	# project this huge raster file to a very large swath	
	b3.oversize <-  projectRaster(b2, crs=wiwa.crs)  # this takes a minute or two, even with the lo-res map...
																									 # and it takes even longer in full resolution

	# now, we crop our oversized b3 to be rectangular according to the extent of wiwa.SPDF.proj
	# which includes the breeding and wintering ranges
	b3 <- crop(b3.oversize, extent(wiwa.SPDF.proj))
	
	
	# this puts a 1.5 Gb file down to disk
	#writeRaster(b3.oversize, "/tmp/b3_oversize")  
	# I ended up putting that as well as the cropped version (in ncdf format as well as native raster) 
	# on barshis machine in:\\
	#/Users/eriq/Documents/UnsyncedSimulationOutputs/WIWA-MapStuff/b3_cropped.grd
	#/Users/eriq/Documents/UnsyncedSimulationOutputs/WIWA-MapStuff/b3_cropped.gri
	#/Users/eriq/Documents/UnsyncedSimulationOutputs/WIWA-MapStuff/b3_cropped.nc
	#/Users/eriq/Documents/UnsyncedSimulationOutputs/WIWA-MapStuff/b3_oversize.grd
	#/Users/eriq/Documents/UnsyncedSimulationOutputs/WIWA-MapStuff/b3_oversize.gri
	
	# then I copied the ncdf file into:
	#/Users/eriq/Documents/work/assist/kristenruegg/WIWA_popgen/WIWA_popgen_analysis/basemap/b3_cropped.nc


}  else {
  dir.create("mapstuff")  # make a directory to put the stuff
  
  message("\n\nStarting download of basemap: b3_cropped.nc\n\n")
  bdown(url = "https://dl.dropboxusercontent.com/u/19274778/mapstuff/b3_cropped.nc", 
        file = "mapstuff/b3_cropped.nc")
	b3 <- brick("mapstuff/b3_cropped.nc")
}



# now, another task that takes a long time, so we will not re-do it if not necessary
if(REGENERATE_POLY_RASTS==TRUE) {
  ## NOTE, IF YOU WANT TO MAKE THESE RASTERS OF THE POLYGONS, YOU WILL NEED TO 
  ## GET THOSE SHAPEFILES FROM THE BIRD RANGE FOLKS.  SEE THE MANUSCRIPT AND TERMS-BLI.md IN
  ## THE GITHUB REPOSITORY FOR A CITATION.
  
  # here are some maneuvers with the range map shapefiles.
  # Here get the distribution maps as shapefiles.  We read them in as spatialPolygonsDataFrames
  # We will use these to clip the Geneland results so as to just include the breeding range.
  wiwa.SPDF <- readOGR(dsn="./shapefiles/Wilsonia_pusilla_9151_NS", layer="Wilsonia_pusilla_9151_NS")
  
  # once we have done that, we need to project it
  wiwa.SPDF.proj <- spTransform(wiwa.SPDF,  wiwa.crs) 
  
  # then whittle it down to just the breeding range and also get the wintering range:
  wiwa.breed.SPDF <- wiwa.SPDF.proj[wiwa.SPDF.proj$SEASONAL==2,]
  wiwa.winter.SPDF <- wiwa.SPDF.proj[wiwa.SPDF.proj$SEASONAL==3,]
  

	# now, let us see if we can rasterize that
	wibreed.rast <- rasterize(wiwa.breed.SPDF, b3)
	wiwinter.rast <- rasterize(wiwa.winter.SPDF, b3)
	
	# and write those into the shapefiles directory:
	writeRaster(wibreed.rast, "shapefiles/wibreed_rast", "CDF")
	writeRaster(wiwinter.rast, "shapefiles/wiwinter_rast", "CDF")
	
} else {
	# if not regenerating, then just download them and read them in
  dir.create("mapstuff")  # make a directory to put the stuff
    
  message("\n\nStarting download of wibreed_rast.nc\n\n")
  bdown(url = "https://dl.dropboxusercontent.com/u/19274778/mapstuff/wibreed_rast.nc", 
                file = "mapstuff/wibreed_rast.nc")
  
  message("\n\nStarting download of wiwinter_rast.nc\n\n")
  bdown(url = "https://dl.dropboxusercontent.com/u/19274778/mapstuff/wiwinter_rast.nc", 
        file = "mapstuff/wiwinter_rast.nc")
  
	wibreed.rast <- raster("mapstuff/wibreed_rast.nc")
	wiwinter.rast <- raster("mapstuff/wiwinter_rast.nc")
}

# get the sampling locations properly projected:
# these are the weighted average sampling locations for each place
Pop.Centers.sp <- SpatialPointsDataFrame(Pop.Centers[,c("lon", "lat")], Pop.Centers, coords.nrs=1:2, proj4string=CRS("+proj=longlat"))
# these are all the actual unique lat longs:
BreedSampLocations <- SpatialPoints(unique(cbind(WA.B$Long, WA.B$Lat)), proj4string=CRS("+proj=longlat"))
# these are the weighted average lat longs for the different wintering groups
WintGroupLocs <- SpatialPointsDataFrame(WW.assign.df[,c("long", "lat")], WW.assign.df,  proj4string=CRS("+proj=longlat"))

# now project them all:
PC.proj <- spTransform(Pop.Centers.sp, wiwa.crs)
BSL.proj <- spTransform(BreedSampLocations, wiwa.crs)
WGL.proj <- spTransform(WintGroupLocs, wiwa.crs)

# OK! That is quite pretty so far.  Now I need to read in the geneland output
# and put it into a raster brick (one layer for each cluster---retaining all the posterior
# values).  
# then we have to disaggregate it so that is closer in size to the wibreed raster.
# Then mask it according to its intersection with wibreed.rast.  And then plot them 
# with the transparency varying according to some function of the posterior probability. 



  # the resampling step can take quite long, maybe 5 minutes, (much less than some of the 
  # earlier maneuvers.  But I don't think it will make sense to store these at this point
	# first get the posteriors:
	pp <- read.table(ProbAFile, header=T)
	
	# then try rasterizing.  Get a whole list of them and turn them into a rasterBrick
	ppb <- brick(lapply(3:8, function(x)  rasterFromXYZ(pp[,c(1,2,x)], crs=CRS("+proj=longlat"))))
	
	
	# then we have to project it to wiwa.crs:
	ppb1 <- projectRaster(ppb, crs=wiwa.crs)
	
	# and we might have to resample it to be like b3
	ppb2 <- resample(x=ppb1, y=b3)
	
	# and then mask it according to the breeding habitat
	ppb3 <- mask(ppb2, wibreed.rast)


# here are some crazy hacks I have to do to get the plots to come out
na3 <- b3
na3[] <- NA  # make a raster like the earth map but all NAs
wib.na <- wibreed.rast
wib.na[] <- NA # make another dummy that I can plot with no real result


# now, we plot na3 (and nothing shows up), the we
# plot nothing three times.  For some inexplicable reason, this resets the plot
# area.  So thing get plotted in a different place.
#png(file="try.png", width=4098, height=3047)  # pdf looks terrible.  This looks great and is only 3.7 Mb or so
jpeg(file="Oct21-6-pops-basemap-way-bigger.jpg", width=floor(4098*2.51), height=floor(3047*2.47), quality=100, res=600)
#tiff(file="Oct21-6-pops-basemap.tiff", width=floor(4098*1.02), height=floor(3047*1.05), res=600, compression="none")

plotRGB(na3) # plot the whole thing.
for(i in 1:3) {plot(wib.na, add=T, legend=F)} 

# now that we have tricked it into resetting the plot area (is that a bug in raster?) 
# we can go ahead and plot our real data by adding to it.
plotRGB(b3, maxpixels=ncell(b3), add=T)
plot(wibreed.rast, maxpixels=ncell(wibreed.rast), add=T, col=rgb(255, 255, 255, 80, maxColorValue=255), legend=F)
plot(wiwinter.rast, maxpixels=ncell(wiwinter.rast), add=T, col=rgb(255, 255, 255, 80, maxColorValue=255), legend=F)


# in order to properly add these it seems to have to plot them as rasterLayers not as parts of a rasterBrick
for(y in 1:6) {
	plot(ppb3[[y]], add=T, maxpixels=1.4*ncell(ppb3[[y]]), col=wiwa.cols[[y]], legend=F, useRaster=F, asp=1)
	print(par("plt"))
}

# here we can put the letter labels on there, if need be
text(PC.proj, as.character(PC.proj$Letter), cex=.02)

# breeding sampling location points
#plot(BSL.proj, add=T, pch=19, col="black", cex=1.0)
#plot(BSL.proj, add=T, pch=1, col="white", lwd=6, cex=1.0)

# wintering group average locations.  I just print out their number...
text(WGL.proj, as.character(1:nrow(WGL.proj)), cex=.1)


# here is some code to make some lines
if(FALSE) {
#	AssLatLon.spdf <- spTransform(SpatialPointsDataFrame(AssLatLong[,c("lon", "lat")], AssLatLon, coords.nrs=1:2, proj4string=CRS("+proj=longlat")), wiwa.crs)
  AssLatLon.spdf <- SpatialPointsDataFrame(AssLatLong[,c("lon", "lat")], AssLatLong, coords.nrs=1:2, proj4string=CRS("+proj=longlat"))

	WM.gr.no.miss <- WM.gr[!is.na(WM.gr$Long), ]
	AssLatLon.spdf <- AssLatLon.spdf[!is.na(WM.gr$Long), ]  # toss the ones with no origin information, too
#	WM.gr.spdf <- spTransform(SpatialPointsDataFrame(WM.gr.no.miss[,c("Long", "Lat")], WM.gr.no.miss, coords.nrs=c(7,6), proj4string=CRS("+proj=longlat")), wiwa.crs)
	WM.gr.spdf <- SpatialPointsDataFrame(WM.gr.no.miss[,c("Long", "Lat")], WM.gr.no.miss, coords.nrs=c(7,6), proj4string=CRS("+proj=longlat")) 

	# here we compute the great circle points (100 for each bird)
	gCircs <- gcIntermediate(WM.gr.spdf, AssLatLon.spdf, n=100, addStartEnd=T)
	gCircs.tf <- lapply(gCircs, function(x) spTransform(SpatialPoints(x, proj4string=CRS("+proj=longlat")), wiwa.crs))
	lapply(gCircs.tf, function(x) lines(x$lat, x$lon, lwd=.1))
	# this makes it look like the projections are way off...
}




#plot(coords.proj, add=T, col="white", cex=2.45)  # I have to scale these so they show up!! Not sure how to do that automatically with png
#plot(coords.proj, add=T, col="red", cex=2.17)
dev.off()


