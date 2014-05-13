

# here are some colors for the pie charts
pie.cols <- c(
	rgb(152, 78, 163, maxColorValue=255), # purple
	rgb(77, 175, 74, maxColorValue=255), # green
	rgb(255, 255, 51, maxColorValue=255),  # yellow
	rgb(247, 129, 191, maxColorValue=255), # pink
	rgb(255, 127, 0, maxColorValue=255),  # orange
	rgb(228, 26, 28, maxColorValue=255) # red
)


# set up some colors with fading transparencies.  These are for the raster maps. 
wiwa.cols <- list(
	rgb(228, 26, 28, c(rep(0,100), 0:200), maxColorValue=255), # red
	rgb(247, 129, 191, c(rep(0,100), 0:200), maxColorValue=255), # pink
	rgb(77, 175, 74, c(rep(0,100), 0:200), maxColorValue=255), # green
	rgb(255, 255, 51, c(rep(0,100), 0:200), maxColorValue=255),  # yellow
	rgb(152, 78, 163, c(rep(0,100), 0:200), maxColorValue=255), # purple
	rgb(255, 127, 0, c(rep(0,100), 0:200), maxColorValue=255)  # orange
	)

 