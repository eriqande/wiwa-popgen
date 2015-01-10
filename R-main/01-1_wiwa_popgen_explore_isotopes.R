
library(mapdata)
library(maptools)
library(ggplot2)
library(dplyr)



#### Get data generated in 01
# get the data.  they are in gl.coord and gl.geno after you load the following
got.these <- load("outputs/WIWA-main-carryover-variables.Rda")

if(!all(c("gl.geno", "gl.coord") %in% got.these)) {
  stop("Failed to successfully load variables gl.geno and/or gl.coord from file WIWA-main-carryover-variables.Rda")
}


#### Get the isotope data   ####
# get all the isotope data (wintering and breeding)
isotopes_all <- read.table("data/wiwa-isotopes.txt", header = T, row.names = 1)

# pick out values just for the breeders in gl.geno
isotopes <- cbind(isotopes_all[rownames(gl.geno), ])
rownames(isotopes) <- rownames(gl.geno)


# Then make a data frame that has their locations on it
iso <- cbind(gl.coord, isotopes)




#### Make a simple ggplot
ggplot(iso, aes(x = Longitude, y = Latitude, color = isotopes)) + 
  geom_jitter(position = position_jitter(5,5), size = 3) + 
  scale_color_gradientn(colours = rev(rainbow(7))) + coord_fixed(1.3)

# figure out what the lat-long bounds should be:
range(iso$Longitude)
range(iso$Latitude)


# now prepare a base map
if (!rgeosStatus()) gpclibPermit()
gshhs.f.b <- "/Users/eriq/Maps/gshhg-bin-2.3.3/gshhs_l.b"
noram <- getRgshhsMap(gshhs.f.b, xlim = c(-175, -50), ylim = c(30, 74)) %>%
  fortify()


# then plot that basemap in gray in the background and put the
# points over it
isopoints <- ggplot() + geom_polygon(data = noram, aes(x=long, y = lat, group = group), fill = "grey91") + 
  coord_fixed(1.5, xlim = c(-175,-50), ylim = c(30,74)) + 
  geom_jitter(data = iso, aes(x=Longitude, y = Latitude, color = isotopes), position = position_jitter(1,1), size = 1) + 
  scale_color_gradientn(colours = rev(rainbow(7))) + 
  theme_bw()


ggsave("wiwa-iso-points.pdf", isopoints, width = 11, height = 8)
