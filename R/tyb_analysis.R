#### Analysis script for Tanana Burbot Telemetry
#### Matt Tyers, fall 2023


#### Data & setup


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

# coordinate conversion
telem_albers <- telemdata %>% 
  select(c("longitude", "latitude")) %>%    # extracting coord matrix
  sf::sf_project(pts=., to=akalbers) %>%   # re-projecting to Alaska Albers
  as.data.frame                             # make it a data.frame
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
telemdata <- cbind(telemdata, telem_albers, telem_segvert) %>%
  filter(use=="Y" & total_length!="TEST TAG")   #### important!!!
telemdata$total_length <- as.numeric(telemdata$total_length)
head(telemdata)

# grabbing the distance from lowest point on the river, as a single metric of
# migratory position on the river
telemdata$mouthdist <- with(telemdata, mouthdist(seg=seg, vert=vert, rivers=tyb_trim))
telemdata$mouthdist <- telemdata$mouthdist/1000  # re-expressing in km
### idea: what if instead of measuring from an arbitrary junction, we measure
### from the Tanana-Yukon confluence, with negative numbers meaning below this?



# data: looking for cases in which multiple records exist for unique combinations
# of individual x flight
thetab <- table(telemdata$unique_id_num, telemdata$flight_num)
sum(thetab>1)

# returning the offending rows...
badrows <- list()
k <- 1
for(i in sort(unique(telemdata$unique_id_num))) {
  for(j in sort(unique(telemdata$flight_num))) {
    if(sum(telemdata$unique_id_num==i & telemdata$flight_num==j)>1) {
      badrows[[k]] <- telemdata[telemdata$unique_id_num==i & telemdata$flight_num==j, ]
      k <- k+1
    }
  }
}

# print the duplicate entries
badrows

# how far apart are each pair of points?  spatially, which to keep might not matter
sapply(badrows, function(x) dist(select(x, c("AlbersX", "AlbersY"))))
# max 1.6km, most <1km

# much difference in power?
sapply(badrows, function(x) min(x$power)/max(x$power))
# min ratio of 75%


# Transforming from long-format data.frame to "3d" list
# This might make it easier to grab longitudinal slices of a given variable
df23d <- function(df, x1, x2, stopifmultiple=T, tiebreaker=NULL, tiemax=T) {  # make this operation a function so i can use it later!
  # df: input data.frame
  # x1: variable 1 name (as character) - will correspond to dim 1 of list elements
  # x2: variable 2 name (as character) - will correspond to dim 2 of list elements
  # stopifmultiple: stop if multiple records are encountered for a given combination
  # tiebreaker: optional variable name to use as a tiebreaker if multiple records
  # tiemax: whether to keep max value, as opposed to min
  
  var1 <- df[,names(df)==x1]
  var2 <- df[,names(df)==x2]
  unique1 <- sort(unique(var1))
  unique2 <- sort(unique(var2))
  
  if(!is.null(tiebreaker)) tievar <- df[,names(df)==tiebreaker]
  
  outlist <- list() 
  for(i in 1:ncol(df)) {
    outlist[[i]] <- matrix(nrow=length(unique1), ncol=length(unique2))
    for(j in 1:length(unique1)) {   ### I should probably rework this with a reshape function!! this is slow
      for(k in 1:length(unique2)) {
        theval <- df[var1==unique1[j] & var2==unique2[k], i]
        if(is.null(theval)) {
          outlist[[i]][j,k] <- NA
        }
        if(length(theval) > 1) {
          if(stopifmultiple) stop(paste("Multiple records exist for x1 =", unique1[j], "and x2 =", unique2[k]))
          else {
            if(i==1) { # don't need to be warned multiple times
              warning(paste("Multiple records exist for x1 =", unique1[j], "and x2 =", unique2[k]))
            }
          }
        } 
        if(length(theval) == 1 | (!stopifmultiple & length(theval) > 1)) {
          if(!is.null(tiebreaker)) {
            whichone <- ifelse(tiemax, which.max(tievar[var1==unique1[j] & var2==unique2[k]]),
                                       which.min(tievar[var1==unique1[j] & var2==unique2[k]]))
          } else {
            whichone <- 1
          }
          outlist[[i]][j,k] <- theval[whichone]  
        }
      }
    }
    rownames(outlist[[i]]) <- unique1
    colnames(outlist[[i]]) <- unique2
  }
  names(outlist) <- colnames(df)
  return(outlist) 
}
telem3d <- df23d(df=telemdata, x1="unique_id_num", x2="flight_num", 
                 stopifmultiple = F, tiebreaker="power")  # wow this is slow





#### Analysis!

### length composition

# now we have a way of getting records for each individual fish
# check if we can just grab the first flight number
colMeans(!is.na(telem3d$total_length))  # all good: all rows are not NA

# grab unique lengths
lengths <- telem3d$total_length[,1]
tagging_locations <- telem3d$tagging_location[,1]

# histogram form
hist(lengths, xlab="Total Length (mm)", main="")   # insert title in main=

# any difference by tagging location?
boxplot(lengths~tagging_locations, ylab="Total Length (mm)", xlab="Tagging Location")
lm(lengths~tagging_locations) %>% anova
# yes!! fish tagged lower down were bigger

aov(lengths~tagging_locations) %>% TukeyHSD
# family-wise significant differences between: 
# - Middle & Lower
# - Upper & Lower


## to make:
## - table output for all fish: n, mean, sd, se
## - table output for each tagging loc: n, mean, sd, se
n <- c(length(lengths), table(tagging_locations))
mn <- c(mean(lengths), tapply(lengths, tagging_locations, mean))
sds <- c(sd(lengths), tapply(lengths, tagging_locations, sd))
ses <- sds/sqrt(n)
round_digits <- 1   #### can change number of digits for table here
ci95 <- paste0("(",round(mn + qnorm(0.025)*ses, round_digits), " - ",
               round(mn + qnorm(0.975)*ses, round_digits), ")")
length_tbl_mns <- data.frame(n, round(mn, round_digits), 
                             round(sds, round_digits), 
                             round(ses, round_digits), ci95)
rownames(length_tbl_mns) <- c("All Sections", names(n)[2:4])
colnames(length_tbl_mns) <- c("n", "Mean", "SD", "SE", "95% CI for Mean")
length_tbl_mns                # print to console
knitr::kable(length_tbl_mns)  # print to Markdown doc

## - table output for binned lengths: n, prop, se
## - table output for binned lengths by tagging loc: n, prop, se
# length_bins <- seq(600, 1050, by=50)  #### can change length bins here
length_bins <- c(600, 700, 800, 900, 1050)
lengths_binned <- cut(lengths, length_bins)
lengthtab <- table(lengths_binned, tagging_locations)
par(mfrow=c(1,1))
lengthtab  # print to console
mosaicplot(t(lengthtab), las=1, main="",
           ylab="Total Length (mm)", xlab="Tagging Location",
           col=grey.colors(nrow(lengthtab), rev=T))
chisq.test(lengthtab)
# no surprise, a Chi-squared test gives strong evidence of differences in length
# categories by tagging location

length_tbl_props <- matrix(nrow=4, ncol=3*nrow(lengthtab))
n_props <- cbind(table(lengths_binned), lengthtab)
p_props <- t(t(n_props)/colSums(n_props))
se_props <- sqrt(p_props*(1-p_props)/(n_props-1))
round_digits <- 2  #### can change number of digits for table here
for(i in 1:4) {
  length_tbl_props[, (i-1)*3 + 1] <- n_props[i,]
  length_tbl_props[, (i-1)*3 + 2] <- round(p_props[i,], round_digits)
  length_tbl_props[, (i-1)*3 + 3] <- round(se_props[i,], round_digits)
  # colnames(length_tbl_props)[(i-1)*3 + 1] <- paste("n", rownames(n_props)[i])
  # colnames(length_tbl_props)[(i-1)*3 + 2] <- paste("p_hat", rownames(n_props)[i])
  # colnames(length_tbl_props)[(i-1)*3 + 3] <- paste("se(p_hat)", rownames(n_props)[i])
}
rownames(length_tbl_props) <- rownames(length_tbl_mns)
colnames(length_tbl_props) <- apply(expand.grid(c("n","p_hat","se(p_hat)"), rownames(n_props)), 1, paste, collapse=" ")

length_tbl_props  # print to console
knitr::kable(length_tbl_props)  # print to Markdown - THIS WILL NEED TO BE REFORMATTED TO FIT

# alternative way to print this
length_tbl_props_alt <- rbind(rep(c("n","p_hat","se(p_hat)"), nrow(n_props)),
                              length_tbl_props)
newcolnames <- rep("", ncol(length_tbl_props)) 
newcolnames[seq(1,length(newcolnames), by=3)] <- rownames(n_props)
colnames(length_tbl_props_alt) <- newcolnames
length_tbl_props_alt
knitr::kable(length_tbl_props_alt)  # print to Markdown



# looking at relationship between length and river location, by survey.
# generally fish lower down are bigger, but the smallest fish are somewhere upper middle?
par(mfrow=c(4,4))
for(i in 1:ncol(telem3d$mouthdist)) {
  plot(lengths ~ telem3d$mouthdist[,i], main=colnames(telem3d$mouthdist)[i],
       xlim=range(telem3d$mouthdist, na.rm=T), 
       xlab="River Position (km)", ylab="Total Length (mm)")
}
# did they disperse?
par(mfrow=c(4,4))
for(i in 1:ncol(telem3d$mouthdist)) {
  plot(lengths ~ telem3d$mouthdist[,i], main=colnames(telem3d$mouthdist)[i],
       xlim=range(telem3d$mouthdist, na.rm=T), 
       xlab="River Position (km)", ylab="Total Length (mm)", 
       col=as.numeric(as.factor(telem3d$tagging_location[,1]))+1, pch=16)
} # not as much as I expected
### THINK ABOUT HOW TO VISUALIZE THIS BETTER, PERHAPS SPATIALLY


# proportion



## run timing??



## estimate proportions within each river section, per survey



## identify probable spawning areas in the mainstem Tanana River during late January



## describe seasonal distributions and migrations (this will be the big one)

# maybe generate list by survey date? or wide format
