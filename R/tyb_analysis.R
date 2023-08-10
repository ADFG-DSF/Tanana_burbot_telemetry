## Analysis script for Tanana Burbot Telemetry
## Matt Tyers, fall 2023

# load packages
library(riverdist)   # for spatial river network analysis (note: used v0.16.0)
library(tidyverse)   # for streamlined data manipulation

# load data
load(file="data/tyb.Rdata")  # rivernetwork created from shapefile
telemdata <- read.csv("data/TananaBurbotTelem.csv")
head(telemdata)

## convert long/lat to Alaska Albers (shapefile projection)
# projection string
akalbers <- "+proj=aea +lat_1=55 +lat_2=65 
    +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs 
    +ellps=GRS80"

# conversion
telem_albers <- telemdata %>% 
  select(c("longitude", "latitude")) %>%    # extracting coord matrix
  sf::sf_project(pts=., to=akalbers)        # re-projecting to Alaska Albers
colnames(telem_albers) <- c("AlbersX","AlbersY")
# head(telem_albers)



## quick plot to verify that the points & rivers are in the right place
plot(tyb)
points(telem_albers, pch=16,
       col=1+as.numeric(as.factor(telemdata$use)))
legend("topleft", pch=16, col=c(3,2), title="Use", legend=c("Y","N"))
# kinda looks like a grayling drinking from a silly straw



#### MT note: since we aren't using the fish locations far up/down the Yukon,
#### we could remove the corresponding segments in the river network.
#### This would be the place to do it.

tyb_trim <- trimriver(trim=c(79, 81, 82, 80, 83, 86, 87, 85), rivers=tyb)
# tyb_trim <- trimriver(trim=c(79, 81, 82, 80), rivers=tyb)  # just lower yukon
plot(tyb_trim)
points(telem_albers[telemdata$use=="Y",])

## redefining the river mouth, since we chopped it off (might not be needed)
# showends(80,tyb_trim)
tyb_trim$mouth$mouth.seg <- 80
tyb_trim$mouth$mouth.vert <- 1316


## finishing the data stuff ASSUMING WE USE tyb_trim
# converting to river coordinates
telem_segvert <- xy2segvert(x=telem_albers[,1], y=telem_albers[,2], rivers=tyb_trim)

# bundling albers & segvert with original data.frame
telemdata <- cbind(telemdata, as.data.frame(telem_albers), telem_segvert)
