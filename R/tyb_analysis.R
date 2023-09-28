#### Preliminary analysis script for Tanana Burbot Telemetry
#### Matt Tyers, fall 2023


#### Data & setup


# load packages 
# note: all packages are available on CRAN and can be installed using install.packages()
library(riverdist)   # for spatial river network analysis (note: used v0.16.0)
library(tidyverse)   # for streamlined data manipulation
library(knitr)       # for table output to Markdown
library(jagsUI)      # for running Bayesian survival analysis 
library(jagshelper)  # for plotting & model diagnostics
# Note: JAGS can be downloaded here: https://sourceforge.net/projects/mcmc-jags/


# load data (from folder available on Github)
load(file="data/tyb.Rdata")  # rivernetwork created from shapefile
# telemdata <- read.csv("data/TananaBurbotTelem.csv")
telemdata <- read.csv("data/TananaBurbotTelem2.csv")
head(telemdata)  # looking at the first few rows


### some data manipulation & conversion...

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
  filter(use=="Y" & total_length!="TEST TAG" & !(flight_num %in% c(2.5, 3.5)))   #### important!!!
# misc data reformatting
telemdata$total_length <- as.numeric(telemdata$total_length)
telemdata$date <- as.Date(telemdata$date, format="%m/%d/%Y")
head(telemdata)

# grabbing the distance from lowest point on the river, as a single metric of
# migratory position on the river
telemdata$mouthdist <- with(telemdata, mouthdist(seg=seg, vert=vert, rivers=tyb_trim))
telemdata$mouthdist <- telemdata$mouthdist/1000  # re-expressing in km

### idea: what if instead of measuring from an arbitrary junction, we measure
### from the Tanana-Yukon confluence, with negative numbers meaning below this?

# telemdata$upriver <- NA
# for(i in 1:nrow(telemdata)) {
#   telemdata$upriver[i] <- upstream(startseg=52, startvert=1119,
#                                    endseg=telemdata$seg[i],
#                                    endvert=telemdata$vert[i],
#                                    rivers=tyb_trim)/1000 # re-expressing in km
# } # actually this is a slow way of doing it.
telemdata$upriver <- telemdata$mouthdist - tyb_trim$lengths[80]/1000
# all.equal(telemdata$upriver, telemdata$upriver2)  # check to make sure the above worked



## We will end up creating a restructured data object, described later on.
## First, we need to see if there are cases in which multiple records exist for 
## unique combinations of individual x flight

# quick tabluate
thetab <- table(telemdata$unique_id_num, telemdata$flight_num)
sum(thetab>1)  # how many times there are non-unique entries

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


## Big question: will we lose information if we keep only one entry where there
## are multiple records for individual x flight?

# how far apart are each pair of points?  spatially, which to keep might not matter
sapply(badrows, function(x) dist(select(x, c("AlbersX", "AlbersY"))))
# max 1.6km, most <1km  - NOT VERY FAR APART

# is there much difference in receiver power?
sapply(badrows, function(x) min(x$power)/max(x$power))
# min ratio (min/max) of 75%  - MIGHT MAKE SENSE TO USE THIS AS TIEBREAKER METRIC


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

## first, add a vector of mean flight dates associated with each number
## (easier for interpretation)
mnDates <- as.Date(round(tapply(telemdata$date, telemdata$flight_num, mean)), origin="1970-01-01")
mnDates_df <- data.frame(flight_num=as.numeric(names(mnDates)), mnDate=mnDates)
telemdata <- left_join(telemdata, mnDates_df)  # join these to original data

## finally, make a 3d list!!
telem3d <- df23d(df=telemdata, x1="unique_id_num", x2="flight_num", 
                 stopifmultiple = F, tiebreaker="power")  # wow this is slow






#### Analysis!

### length composition

# now we have a way of getting records for each individual fish

# check if we can just grab the first flight number (marking event) for length
colMeans(!is.na(telem3d$total_length))  # all good: all rows are not NA

# making sure total_length doesn't change for a given individual
all(apply(telem3d$total_length, 1, function(x) length(unique(x[!is.na(x)])))==1)


# grab unique lengths
lengths <- telem3d$total_length[,1]
tagging_locations <- telem3d$tagging_location[,1]

# histogram form
hist(lengths, xlab="Total Length (mm)", main="")   # insert title in main=

# one-way ANOVA: is there any difference by tagging location?
boxplot(lengths~tagging_locations, ylab="Total Length (mm)", xlab="Tagging Location")
lm(lengths~tagging_locations) %>% anova
## yes!! evidence of difference in mean lengths by tagging area
## fish tagged lower down were bigger

# Multiple comparison (Tukey's HSD)
aov(lengths~tagging_locations) %>% TukeyHSD
# family-wise significant differences between: 
# - Middle & Lower
# - Upper & Lower


## MAKING ASSOCIATED TABLES

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

# alternative way to print this (also needs to be fixed in Excel/Word)
length_tbl_props_alt <- rbind(rep(c("n","p_hat","se(p_hat)"), nrow(n_props)),
                              length_tbl_props)
newcolnames <- rep("", ncol(length_tbl_props)) 
newcolnames[seq(1,length(newcolnames), by=3)] <- rownames(n_props)
colnames(length_tbl_props_alt) <- newcolnames
length_tbl_props_alt       # print to console
knitr::kable(length_tbl_props_alt)  # print to Markdown


### SOME EXPLORATORY PLOTS - I think this was handled better later on

# looking at relationship between length and river location, by survey.
# generally fish lower down are bigger, but the smallest fish are somewhere upper middle?
par(mfrow=c(4,4))
for(i in 2:ncol(telem3d$upriver)) {  # starting at 2 to not include tagging date
  plot(lengths[telem3d$current_state[,i]!="dead"] ~ telem3d$upriver[,i][telem3d$current_state[,i]!="dead"], main=mnDates[i],
       xlim=range(telem3d$upriver, na.rm=T), 
       xlab="River Position (km)", ylab="Total Length (mm)")
}
# did they disperse?
par(mfrow=c(4,4))
for(i in 2:ncol(telem3d$upriver)) {  # starting at 2 to not include tagging date
  plot(lengths[telem3d$current_state[,i]!="dead"] ~ telem3d$upriver[,i][telem3d$current_state[,i]!="dead"], main=mnDates[i],
       xlim=range(telem3d$upriver, na.rm=T), 
       xlab="River Position (km)", ylab="Total Length (mm)", 
       col=as.numeric(as.factor(telem3d$tagging_location[,1]))[telem3d$current_state[,i]!="dead"]+1, pch=16)
} 
### THINK ABOUT HOW TO VISUALIZE THIS BETTER, PERHAPS SPATIALLY



## run timing??




# ## identify probable spawning areas in the mainstem Tanana River during late January
# ###### THIS IS LONG RUNNING AND IS CURRENTLY COMMENTED OUT

# Making a couple new data objects
telemdata1 <- telemdata %>% filter(flight_num>=1 & current_state!="dead")
telem_alive <- telemdata %>% filter(current_state!="dead")

# 
# # empirical k-functions 
# # review this methodology!  it is very slow and doesn't seem to gain anything
# par(mfrow=c(4,4))
# with(telemdata1, kfunc(seg=seg, vert=vert, survey=mnDate,
#       rivers=tyb_trim, envreps=100, maxdist=200000))
# 
# # kernel density (network)
# # darn, this doesn't seem to gain anything either
# # network_dens <- with(telemdata1, makeriverdensity(seg=seg, vert=vert, survey=mnDate,
# #                                        rivers=tyb_trim, bw=20*1000))
# network_dens <- with(telemdata1, makeriverdensity(seg=seg, vert=vert, survey=mnDate,
#                                        rivers=tyb_trim, bw=20*1000,
#                                        kernel="rect"))
# # plotriverdensitypoints(network_dens)  # maybe more points than needed
# par(mfrow=c(4,4))
# plot(network_dens, scalebyN=F)
# par(mfrow=c(1,1))
# plot(network_dens, points=F, scalebyN = F)  # maybe I should try animating this


## kernel density (network) version 2.0
## First objective: show kernel density for surveys identified as spawning season
## for all individuals (this will not be the same as survey by itself). Identify
## the survey preceding spawning time for all individuals, and contrast the two.
## Contrast pre/spawning for all three spawning seasons, and contrast all seasons together.

# using data.frame telem_alive, define a new column for spawning season (where spawning_use is Y)
telem_alive$spawn <- with(telem_alive, ifelse(!is.na(spawning_use) & spawning_use=="Y",
                          paste(spawning_season, "spawn", sep="_"), NA))

# then go by individual, and back-fill prespawn, defined as last observation before spawn
indivs <- sort(unique(telem_alive$unique_id)) # vector of unique individuals
telem_alive$spawn_diff <- telem_alive$prespawn_diff <- NA # initializing columns
for(i_indiv in seq_along(indivs)[-300]) {
  df_indiv <- subset(telem_alive, unique_id==indivs[i_indiv])
  df_indiv$prespawn <- NA        # this will be prespawn season
  df_indiv$prespawn_diff <- df_indiv$spawn_diff <- NA   # this will be number of days between prespawn and spawn
  spawnseasons <- unique(df_indiv$spawn[!is.na(df_indiv$spawn)])
  df_indiv <- df_indiv[order(df_indiv$date),]  # to make sure these are in order by date
  for(i_season in seq_along(spawnseasons)) {  # backfilling prespawn...
    ## backfilling to last-known location...
    # df_indiv$prespawn[which(df_indiv$spawn==spawnseasons[i_season])-1] <- 
    #   paste(substr(spawnseasons[i_season], 1, 9), "prespawn", sep="_")
    # df_indiv$spawn_diff[which(df_indiv$spawn==spawnseasons[i_season])] <-
    #   diff(df_indiv$date[which(df_indiv$spawn==spawnseasons[i_season])-1:0])
    # df_indiv$prespawn_diff[which(df_indiv$spawn==spawnseasons[i_season])-1] <-
    #   diff(df_indiv$date[which(df_indiv$spawn==spawnseasons[i_season])-1:0])
    # needs to be before spawning survey AND closest to diff of 128
    
    ## backfilling to most-comparable date to 2018-2019 spawning season
    ## that is, closest to 128-day difference
    diffs <- df_indiv$date[which(df_indiv$spawn==spawnseasons[i_season])] - df_indiv$date
    df_indiv$prespawn[which.min(abs((diffs-128)[diffs>0]))] <- 
      paste(substr(spawnseasons[i_season], 1, 9), "prespawn", sep="_")
    df_indiv$prespawn_diff[which.min(abs((diffs-128)[diffs>0]))] <- 
      diffs[which.min(abs((diffs-128)[diffs>0]))]
    df_indiv$spawn_diff[which(df_indiv$spawn==spawnseasons[i_season])] <- 
      diffs[which.min(abs((diffs-128)[diffs>0]))]
  }
  # i_indiv of 302 is a good test case, i_season of 1
  # print these to test: 
  # df_indiv$orderdate <- order(df_indiv$date)
  # print(df_indiv %>% select(c("spawning_use", "spawning_season", "spawn", "prespawn",
  #                             "spawn_diff", "prespawn_diff", "date", "orderdate")))
  
  telem_alive$spawn[as.numeric(rownames(df_indiv))] <-
    ifelse(!is.na(df_indiv$spawn), df_indiv$spawn, df_indiv$prespawn)
  telem_alive$spawn_diff[as.numeric(rownames(df_indiv))] <- df_indiv$spawn_diff
  telem_alive$prespawn_diff[as.numeric(rownames(df_indiv))] <- df_indiv$prespawn_diff
}

# creating a subset data.frame from just prespawn/spawn surveys
# can also try a max prespawn/spawn_diff number of days?
telem_spawn <- subset(telem_alive, !is.na(spawn))

# creating (network) kernel density thing!
network_dens_spawn <- with(telem_spawn, makeriverdensity(seg=seg, vert=vert, survey=spawn,
                                                        rivers=tyb_trim, bw=50*1000,
                                                        kernel="rect"))
par(mfrow=c(3,2))
plot(network_dens_spawn, scalebyN = F, maxlwd = 15)

# hacking the structure of a riverdensity object to create an object of contrasts
network_dens_contrasts <- network_dens_spawn  # copying the existing structure
densities_contrasts <- list(list(), list(), list())                   # initializing an object that will become the new $densities
for(i_line in seq_along(network_dens_spawn$densverts)) {
  densities_contrasts[[1]][[i_line]] <- network_dens_spawn$densities[[2]][[i_line]] - network_dens_spawn$densities[[1]][[i_line]]
  densities_contrasts[[2]][[i_line]] <- network_dens_spawn$densities[[4]][[i_line]] - network_dens_spawn$densities[[3]][[i_line]]
  densities_contrasts[[3]][[i_line]] <- network_dens_spawn$densities[[6]][[i_line]] - network_dens_spawn$densities[[5]][[i_line]]
}
network_dens_contrasts$densities <- densities_contrasts
network_dens_contrasts$survey <- as.factor(c("2018_2019 contrast", "2019_2020 contrast", "2020_2021 contrast"))

# plotting!! I hope
# plot(network_dens_contrasts, bycol=T, bylwd=F, ramp="stoplight", lwd = 3,
#      main=c("2018_2019 contrast", "2019_2020 contrast", "2020_2021 contrast"),
#      showN=F, points=F,pwr=1)

## modifying riverdist:::plot.riverdensity to produce a diverging color scheme to 
## show contrast.  this is an extreme hack!!
plotriverdensity1 <- function(x,whichplots=NULL,points=TRUE,bycol=TRUE,bylwd=TRUE,maxlwd=10,pwr=0.7,scalebyN=TRUE,ramp="grey",lwd=1,linecol="black",denscol="black",alpha=1,dark=1,showN=TRUE,main=NULL,xlab="",ylab="",add=FALSE,scalebar=TRUE,...) {
  if(!inherits(x, "riverdensity")) stop("Argument x must be an object returned from makeriverdensity().")
  if(dark>1 | dark<0) dark <-1
  if(alpha>1 | alpha<0) alpha <-1
  densities <- x$densities
  endptverts <- x$endptverts
  densverts <- x$densverts
  seg <- x$pointsegs
  vert <- x$pointverts
  survey <- x$survey
  if(!is.factor(survey)) survey <- as.factor(survey)
  rivers <- x$rivers
  if(length(main)==1) main <- rep(main,length(unique(survey)))
  # if(is.null(main) & length(unique(survey))>1) main <- sort(unique(as.character(survey)))
  if(is.null(main) & length(unique(survey))>1) main <- levels(survey)
  
  lines <- rivers$lines
  length <- length(lines)
  xmin <- min(lines[[1]][,1])
  xmax <- max(lines[[1]][,1])
  ymin <- min(lines[[1]][,2])
  ymax <- max(lines[[1]][,2])
  if(length>1) {
    for(j in 2:length) {
      if(min(lines[[j]][,1])<xmin) xmin <- min(lines[[j]][,1])
      if(max(lines[[j]][,1])>xmax) xmax <- max(lines[[j]][,1])
      if(min(lines[[j]][,2])<ymin) ymin <- min(lines[[j]][,2])
      if(max(lines[[j]][,2])>ymax) ymax <- max(lines[[j]][,2])
    }
  }
  
  nsize <- NA
  isurvey <- 1
  # for(surveyi in sort(unique(survey))) {
  for(surveyi in levels(survey)) {
    nsize[isurvey] <- length(seg[survey==surveyi])
    isurvey <- isurvey+1
  }
  
  iisurvey <- 1
  if(is.null(whichplots)) whichplots <- 1:length(unique(survey))
  # whichplotsurvey <- (sort(unique(survey)))[whichplots]
  whichplotsurvey <- (levels(survey))[whichplots]
  
  if(!scalebyN){ 
    isurvey <- 1
    # for(surveyi in sort(unique(survey))) {
    for(surveyi in levels(survey)) {
      for(segi in 1:length(rivers$lines)) {
        densities[[isurvey]][[segi]] <- densities[[isurvey]][[segi]]*max(nsize[whichplots])/nsize[isurvey]
      }
      isurvey <- isurvey+1
    }
  }
  
  for(surveyi in whichplotsurvey) {
    isurvey <- whichplots[iisurvey]
    if(showN) mainforplot <- paste0(main[isurvey],"  (n=",length(seg[survey==surveyi]),")")
    if(!showN) mainforplot <- main[isurvey]
    if(!add) plot(c(xmin,xmax),c(ymin,ymax),col="white",cex.axis=.6,asp=1,xlab=xlab,ylab=ylab,main=mainforplot,...=...)
    
    for(segi in 1:length(rivers$lines)) {
      quants <- (densities[[isurvey]][[segi]]/max(unlist(densities)))^pwr
      if(bycol) {
        # if(ramp=="grey" | ramp=="gray") {
        #   cols <- grey((1-quants)*.8)
        #   denscol <- 1
        #   linecol <- grey(.8)
        # }
        # if(ramp=="red") {
        #   cols <- rgb(1,(1-quants)*.8,(1-quants)*.8)
        #   denscol <- 2
        #   linecol <- rgb(1,.8,.8)
        # }
        # if(ramp=="green") {
        #   cols <- rgb((1-quants)*.8,1,(1-quants)*.8)
        #   denscol <- 3
        #   linecol <- rgb(.8,1,.8)
        # }
        # if(ramp=="blue") {
        #   cols <- rgb((1-quants)*.8,(1-quants)*.8,1)
        #   denscol <- 4
        #   linecol <- rgb(.8,.8,1)
        # }
        # if(ramp=="heat") {
        #   cols <- heat.colors(1000)[ceiling(900*(1-quants))+1]
        #   denscol <- heat.colors(1000)[1]
        #   linecol <- heat.colors(1000)[901]
        # }
        # if(ramp=="stoplight") {
        #   cols <- rainbow(1000)[ceiling(300*(1-quants))+1] 
        #   denscol <- rainbow(1000)[1] 
        #   linecol <- rainbow(1000)[301]
        # }
        # if(ramp=="rainbow") {
        #   cols <- rainbow(1000)[ceiling(700*(1-quants))+1]
        #   denscol <- rainbow(1000)[1] 
        #   linecol <- rainbow(1000)[701]
        # }
        
        nnn <- max(abs(unlist(densities)))
        nnn <- max(abs(unlist(densities[[isurvey]])))
        allthecols <- c(colorRampPalette(c("blue","grey90"))(nnn), colorRampPalette(c("grey90","red"))(nnn)[-1])
        cols <- allthecols[densities[[isurvey]][[segi]]+nnn]
        linecol="grey90"
      }
      if(bylwd) {
        lwds <- maxlwd*quants
      }
      if(!bycol) cols <- rep(denscol,(length(endptverts[[segi]])-1))
      if(!bylwd) lwds <- rep(lwd,(length(endptverts[[segi]])-1))
      
      if(alpha<1) cols <- adjustcolor(cols,alpha.f=alpha)
      if(alpha<1) linecol <- adjustcolor(linecol,alpha.f=alpha)
      if(dark<1) cols <- adjustcolor(cols,red.f=dark,green.f=dark,blue.f=dark)
      if(dark<1) linecol <- adjustcolor(linecol,red.f=dark,green.f=dark,blue.f=dark)
      if(dark<1) denscol <- adjustcolor(denscol,red.f=dark,green.f=dark,blue.f=dark)
      
      for(vertsi in 1:(length(endptverts[[segi]])-1)) {
        if(dim(matrix(rivers$lines[[segi]][(endptverts[[segi]][vertsi]):(endptverts[[segi]][vertsi+1]),],ncol=2))[1] > 1) {
          lines(rivers$lines[[segi]][(endptverts[[segi]][vertsi]):(endptverts[[segi]][vertsi+1]),],lwd=1,col=linecol,lend=1)
          if(densities[[isurvey]][[segi]][vertsi] != 0) {
            # lines(rivers$lines[[segi]][(endptverts[[segi]][vertsi]):(endptverts[[segi]][vertsi+1]),],lwd=lwds[vertsi],col=cols[vertsi],lend=1)
            lines(rivers$lines[[segi]][(endptverts[[segi]][vertsi]):(endptverts[[segi]][vertsi+1]),],
                  lwd=lwd*(abs(densities[[isurvey]][[segi]])[vertsi]/nnn)^pwr,col=cols[vertsi],lend=1)
          }
        }
      }
    }
    if(points) riverpoints(seg=seg[survey==surveyi],vert=vert[survey==surveyi],rivers=rivers,pch=21,bg=0,col=denscol)
    if(scalebar) riverdist:::scalebar(rivers)
    iisurvey <- iisurvey+1
  }
}
par(mfrow=c(2,2))
# par(mfrow=c(1,1))
plotriverdensity1(network_dens_contrasts, bycol=T, bylwd=F, lwd = 10,
     main=c("2018_2019 contrast", "2019_2020 contrast", "2020_2021 contrast"),
     showN=F, points=F,pwr=.7)



## kernel density (linear)
## THIS IS A SIMPLIFICATION - instead of on the full river network,
## density is expressed just in terms of upriver position.  The hope was that
## areas of high/low density would be easier to spot.
## It should be noted that this is calculated in terms of network distance 
## from/below the Tanana/Yukon confluence as opposed to position on the Tanana 
## mainstem - this is relevant in the case of fish on tributaries

bandwidth <- 50  # kernel bandwidth, in km (trial and error)

# generating a list of kernel density objects
lin_denses <- with(telemdata1, tapply(upriver, flight_num, density, bw=bandwidth))

# extracting max y-value for plotting
maxy <- max(sapply(lin_denses, function(x) max(x$y)))

# plotting
par(mfrow=c(4,4))
for(i in 1:length(lin_denses)) {
  plot(lin_denses[[i]], xlim=range(telemdata1$upriver), ylim=c(0, maxy), main=mnDates[i])
  xx <- with(telemdata1, upriver[flight_num==names(lin_denses)[i]])
  points(x=xx, y=abs(jitter(rep(0, length(xx)), 0.005)), col=adjustcolor(1, alpha.f=.3))
}
par(mfrow=c(1,1))
for(i in 1:length(lin_denses)) {
  plot(lin_denses[[i]], xlim=range(telemdata1$upriver), ylim=c(0, maxy), main=mnDates[i])
  xx <- with(telemdata1, upriver[flight_num==names(lin_denses)[i]])
  points(x=xx, y=abs(jitter(rep(0, length(xx)), 0.005)), col=adjustcolor(1, alpha.f=.3))
} # try animating this one too


## overlaying all densities
for(i in 1) {
  plot(lin_denses[[i]], xlim=range(telemdata1$upriver), ylim=c(0, maxy), col=adjustcolor(1, alpha.f=.3), main="")
  xx <- with(telemdata1, upriver[flight_num==names(lin_denses)[i]])
  points(x=xx, y=abs(jitter(rep(0, length(xx)), 0.005)), col=adjustcolor(1, alpha.f=.3))
}
for(i in 2:length(lin_denses)) {
  lines(lin_denses[[i]], xlim=range(telemdata1$upriver), ylim=c(0, maxy), col=adjustcolor(1, alpha.f=.3))
  xx <- with(telemdata1, upriver[flight_num==names(lin_denses)[i]])
  points(x=xx, y=abs(jitter(rep(0, length(xx)), 0.005)), col=adjustcolor(1, alpha.f=.3))
}

## Then showing spawning seasons differently!!
spawnseasons <- paste(2018:2020, 2019:2021, "spawn", sep="_")
for(i_season in seq_along(spawnseasons)) {
  thedensity <- density(telem_spawn$upriver[telem_spawn$spawn==spawnseasons[i_season]], bw=bandwidth)
  lines(thedensity, lwd=3)
  text(x=thedensity$x[which.max(thedensity$y)], 
       y=thedensity$y[which.max(thedensity$y)], 
       labels=spawnseasons[i_season], pos=c(2,4,4)[i_season], offset=2.5)
}

## Comparing pre-spawn to spawn one more time
par(mfrow=c(3,1))
preseasons <- paste(2018:2020, 2019:2021, "prespawn", sep="_")
for(i_season in seq_along(spawnseasons)) {
  thedensity <- density(telem_spawn$upriver[telem_spawn$spawn==spawnseasons[i_season]], bw=bandwidth)
  plot(thedensity, xlim=range(telemdata1$upriver), ylim=c(0, maxy), main=spawnseasons[i_season])
  thedensity <- density(telem_spawn$upriver[telem_spawn$spawn==preseasons[i_season]], bw=bandwidth)
  lines(thedensity, lty=2)
}

## any utility from KS tests?
for(i_season in seq_along(spawnseasons)) {
  print(ks.test(telem_spawn$upriver[telem_spawn$spawn==spawnseasons[i_season]],
          telem_spawn$upriver[telem_spawn$spawn==preseasons[i_season]]))
}


## Data checks.  I would like to see how much time elapsed between prespawn and spawn
## for each season, and how consistent this is among years.  Also, when (calendar-wise)
## was identified as prespawn for each year.

# boxplot of time difference for each spawning season
par(mfrow=c(1,2))
with(telem_spawn, boxplot(spawn_diff~spawn, ylab="Difference (days) between prespawn and spawn"))

# median difference
with(subset(telem_spawn, !is.na(spawn_diff)), tapply(spawn_diff, spawn, median, na.rm=T))

# sd of difference
with(subset(telem_spawn, !is.na(spawn_diff)), tapply(spawn_diff, spawn, sd, na.rm=T))


# boxplot of the dates themselves
with(telem_spawn, boxplot(date~spawn, ylab="Dates", col=c(2,4)))

# median dates
with(telem_spawn, as.Date(tapply(date, spawn, median)), origin="1970-01-01")

# sd of dates
with(telem_spawn, tapply(date, spawn, sd))

## Take-home message: pre-spawn dates were not consistent among years.  
## 2018-19 median was late September, as opposed to January in other years.
## Thus migration patterns probably won't be consistent between years 

## Update: after redefining pre-spawn dates in the second two seasons to be more 
## consistent with the first, the inferences don't really change appreciably.



## SHOWING THIS ALL IN ONE PLOT
## two versions of a stacked dotplot (flipping X/Y coordinates)

## Individuals observed in sequential flights are linked by solid lines
## Gaps in observations are indicated by dashed lines

upm1 <- telem3d$upriver  # making the name smaller for less typing
upm1[telem3d$current_state=="dead"] <- NA   ### TAKING OUT DEAD FISH ###
daym <- telem3d$date
parmar <- par("mar")  # storing original margin settings
par(mar=c(6.1, 4.1, 4.1, 2.1))  # tweaking margin for this plot
par(mfrow=c(1,1))
plot(NA, xlim=range(daym, na.rm=T), ylim=range(upm1, na.rm=T),
     ylab="Upriver position (km)", xaxt='n', xlab="")
for(i in 1:nrow(upm1)) {
  points(daym[i,], upm1[i,], col=adjustcolor(1, alpha.f=.5))
  lines(daym[i,], upm1[i,], col=adjustcolor(1, alpha.f=.5))
  lines(daym[i,][!is.na(upm1[i,])], upm1[i,!is.na(upm1[i,])], lty=3, col=adjustcolor(1,alpha.f=.3))
}
axis(side=1, at=as.numeric(mnDates[-1]), mnDates[-1], las=2)
# note: axis ticks are not shown for tagging dates

par(mar=c(5.1, 6.1, 4.1, 2.1))  # tweaking margin for this plot
plot(NA, ylim=c(max(daym, na.rm=T), min(daym, na.rm=T)), xlim=range(upm1, na.rm=T),
     xlab="Upriver position (km)", yaxt='n', ylab="")
for(i in 1:nrow(upm1)) {
  points(upm1[i,], daym[i,], col=adjustcolor(1, alpha.f=.5))
  lines(upm1[i,], daym[i,], col=adjustcolor(1, alpha.f=.5))
  lines(upm1[i,!is.na(upm1[i,])], daym[i,][!is.na(upm1[i,])], lty=3, col=adjustcolor(1,alpha.f=.3))
}
axis(side=2, at=as.numeric(mnDates[-1]), mnDates[-1], las=2)
# note: axis ticks are not shown for tagging dates
par(mar=parmar)  # resetting margins to default state




### estimate proportions within each river section, per survey
# note: combined Yukon with Lower!!

telemdata1$section[telemdata1$section=="Yukon"] <- "Lower"
sectiontable_alive <- with(telemdata1, table(flight_num, section))
mosaicplot(sectiontable_alive, col=grey.colors(3, rev=T), main="", xlab="Flight Number", ylab="")
chisq.test(sectiontable_alive)  # weak evidence of difference in proportions

## generating a table
n_section <- rowSums(sectiontable_alive)
p_section <- sectiontable_alive/n_section
se_section <- sqrt(p_section*(1-p_section)/n_section)

proptable <- matrix(NA,nrow=nrow(p_section), ncol=3*ncol(p_section))
proptable[,c(1,4,7)] <- sectiontable_alive
proptable[,c(2,5,8)] <- round(p_section, 3)
proptable[,c(3,6,9)] <- round(se_section, 3)

rownames(proptable) <- rownames(sectiontable_alive)
colnames(proptable) <- apply(expand.grid(c("n","p_hat","SE(p_hat)"),colnames(sectiontable_alive)), 1, paste, collapse=" ")

proptable         # print to console
kable(proptable)  # print to Markdown



### estimate proportions within each river section, per survey 
### FOR EACH TAGGED STOCK - (Lower/Middle/Upper)

sectiontables_bystock <- with(telemdata1, table(flight_num, section, tagging_location))
sectiontables_bystock  # print to console
apply(sectiontables_bystock, 3, kable)  # curious if this works
for(i in 1:3) print(kable(sectiontables_bystock[,,i]))  # or if this works

# some simple mosaic plots
par(mfrow=c(2,2))
for(i in 1:3) {
  mosaicplot(sectiontables_bystock[,,i], col=T, 
             ylab="Section observed", xlab="Flight number", 
             main=paste("Tagging location:", c("Lower","Middle","Upper")[i]))
}
# par(mfrow=c(3,3))
# for(i in 1:13) {
#   mosaicplot(t(sectiontables_bystock[i,,]), col=T, 
#              ylab="Section observed", xlab="Tagging Location",
#              main=paste("Flight", i))
# }  # I don't think this one is meaningful

# constructing a simplfied(ish) table for this
n_arr <- array(dim=dim(sectiontables_bystock))
for(j in 1:3) {
  for(k in 1:3) {
    n_arr[,j,k] <- rowSums(sectiontables_bystock[,,k])
  }
}
p_arr <- sectiontables_bystock/n_arr
se_arr <- sqrt(p_arr*(1-p_arr)/(n_arr-1))
paste0(sectiontables_bystock,", ",p_arr, ", (",se_arr,")")
toprint <- matrix(nrow=dim(p_arr)[1], ncol=prod(dim(p_arr)[2:3]))
for(j in 1:3) {
  for(k in 1:3) {
    toprint[,(k-1)*3 + j] <- paste0(round(p_arr[,j,k],3), " (",
                                    round(se_arr[,j,k],3), ")")
  }
  # toprint[,2*j-1] <- c("asdf",sectiontables_bystock[,j])
}
toprint <- rbind(c("Lower","Middle","Upper"), toprint)
rownames(toprint) <- c("Section observed", 1:dim(p_arr)[1])
colnames(toprint) <- c("Tagged - Lower","","","Tagged - Middle","","","Tagged - Upper","","")
toprint         # print to console
kable(toprint)  # print to Markdown




## describe seasonal distributions and migrations (this will be the big one)
## INVESTIGATING CALCULATING A SUMMARY METRIC OF TRAVEL FOR EACH INDIVIDUAL
## (such that there is one value for each individual)

# - calculating minimum homerange for each individual
hr <- homerange(unique=telem_alive$unique_id_num,
                survey=telem_alive$flight_num,
                seg=telem_alive$seg,
                vert=telem_alive$vert,
                rivers=tyb_trim)
hr_table <- hr$ranges  # extracting table of homeranges
hr_table$range <- hr_table$range/1000  # converting to km
hist(hr_table$range, main="", xlab="Minimum homerange (km)")

# looking at the locations of the top mover (seems reasonable)
hr_table[which.max(hr_table$range),]
par(mfrow=c(1,1))
plot(tyb_trim)
riverpoints(seg=telem_alive$seg[telem_alive$unique_id_num==115],
            vert=telem_alive$vert[telem_alive$unique_id_num==115],
            rivers=tyb_trim, pch=16)

#### NOTE FOUND A BUG in plot.homerange: tries to plot for indiv 28



## what if we did total observed distance INSTEAD of homerange?
## this is defined as the absolute distance between each possible pairing of
## sequential observations (or non-sequential if individual was not observed)

# to calculate, for each individual:
# - subset
# - sort by date 
# - calculate non-missing distances sequentially
indiv <- sort(unique(telem_alive$unique_id_num))
nobs <- cumuldist <- ndays <- rep(NA, length(indiv))
for(i in 1:length(indiv)) {
  # print(i)
  d1 <- telem_alive[telem_alive$unique_id_num==indiv[i],]
  d2 <- d1[order(d1$date),]
  nobs[i] <- nrow(d2)  # number of observations per individual
  ndays[i] <- d2$date %>% range %>% diff %>% as.numeric  # total number of elapsed days
  cumuldist[i] <- 0
  if(nrow(d2)>1) {
    for(irow in 2:nrow(d2)) {
      cumuldist[i] <- cumuldist[i] + riverdistance(startseg=d2$seg[irow-1],
                                                   startvert=d2$vert[irow-1],
                                                   endseg=d2$seg[irow],
                                                   endvert=d2$vert[irow],
                                                   rivers=tyb_trim)/1000 # make it km
    }
  }
}
dtab <- data.frame(nobs, ndays, cumuldist)  # bundle summary metrics
rownames(dtab) <- indiv

## exploratory plots to look at behavior of possible summary metrics
# cumulative distance
hist(dtab$cumuldist)  
plot(nobs, cumuldist)
boxplot(cumuldist~nobs)

# cumulative distance per possible pair of observations
hist(cumuldist/(nobs-1))  
plot(nobs, cumuldist/(nobs-1))
boxplot(cumuldist/(nobs-1) ~ nobs)  # seems consistent enough to use as metric

# cumulative distance per day
boxplot(cumuldist/ndays ~ nobs)  # less consistent, one big outlier
hist(cumuldist/ndays, breaks=10)


## looking at the top mover again
dtab[which.max(cumuldist),]
par(mfrow=c(1,1))
plot(tyb_trim, empty=T, linecol="grey")
riverpoints(seg=telemdata$seg[telemdata$unique_id_num==115],
            vert=telemdata$vert[telemdata$unique_id_num==115],
            rivers=tyb_trim, pch=16)
text(x=telemdata$AlbersX[telemdata$unique_id_num==115],
     y=telemdata$AlbersY[telemdata$unique_id_num==115],
     labels=telemdata$flight_num[telemdata$unique_id_num==115],
     pos=1, cex=1.5)


## ok, bundle all summary metrics
dtab$dist_per_obs <- dtab$cumuldist/(dtab$nobs-1)   # distance per pair of observations
dtab$dist_per_day <- dtab$cumuldist/ndays
# this is a little awkward, have to include individuals with zero homerange
dtab$homerange <- 0
for(i in 1:nrow(hr_table)) dtab$homerange[rownames(dtab)==hr_table$ID[i]] <- hr_table$range[i]




### INVESTIGATING RELATIONSHIPS BETWEEN DISTANCE STUFF AND OTHER VARIABLES
### This is extremely incomplete, and is still mostly an exploratory analysis.
### At the moment, there are four possible metrics for how much an individual moved:
### - homerange
### - cumulative distance
### - cumulative distance per observation
### - cumulative distance per day
### And there are multiple individual-level variables:
### - River section (observed)
### - Tagging location
### - tagging habitat
### - life history
### - Total length (also considered as categorical)
### - Avg upriver position (also considered as categorical)

# defining per-individual summary metrics
themode <- function(x) unname(names(sort(table(x), decreasing=T))[1])
sectionmode <- unname(apply(telem3d$section, 1, themode))  # river section most often observed
sectionmode[sectionmode=="Yukon"] <- "Lower"       # combining Yukon with Lower
avgupriver <- rowMeans(telem3d$upriver, na.rm=T)   # avg upriver position (km)

## lots of exploratory plots
for(imetric in c(6,3:5)) { # this was just a logical order of columns
  par(mfrow=c(2,4))
  boxplot(dtab[,imetric] ~ sectionmode, xlab="River section", ylab=names(dtab)[imetric])
  boxplot(dtab[,imetric] ~ telem3d$tagging_location[,1], xlab="Tagging location", ylab=names(dtab)[imetric])
  boxplot(dtab[,imetric] ~ telem3d$tagging_habitat[,1], xlab="Tagging habitat", ylab=names(dtab)[imetric])
  boxplot(dtab[,imetric] ~ telem3d$life_history[,1], xlab="Life history", ylab=names(dtab)[imetric])
  plot(dtab[,imetric] ~ telem3d$total_length[,1], xlab="Total length (mm)", ylab=names(dtab)[imetric])
  plot(dtab[,imetric] ~ telem3d$total_length[,1], xlab="Total length (mm) - log scale", ylab=paste(names(dtab)[imetric], "- log scale"), log="xy")
  plot(dtab[,imetric] ~ avgupriver, xlab="Avg upriver position (km)", ylab=names(dtab)[imetric])
  plot(dtab[,imetric] ~ avgupriver, xlab="Avg upriver position (km)", ylab=paste(names(dtab)[imetric], "- log scale"), log="y")
}
# additional relationships
par(mfrow=c(2,2))
boxplot(telem3d$total_length[,1] ~ sectionmode, xlab="River section", ylab="Total length (mm)")
boxplot(telem3d$total_length[,1] ~ telem3d$life_history[,1], xlab="Life history", ylab="Total length (mm)")
boxplot(telem3d$total_length[,1] ~ telem3d$tagging_habitat[,1], xlab="Tagging habitat", ylab="Total length (mm)")
mosaicplot(table(sectionmode, telem3d$life_history[,1]), xlab="River section", ylab="Life history", main="",col=T)

## considering total length and upriver position as categorical
lengthcut <- cut(telem3d$total_length[,1], breaks=c(600, 700, 800, 900, 1050))
uprivercut <- cut(avgupriver, breaks=c(0, 300, 500, 800, 1000))
par(mfrow=c(2,4))
par(mar=c(7.1, 4.1, 4.1, 2.1))
for(imetric in c(6,3:5)) { # this was just a logical order of columns
  boxplot(dtab[,imetric] ~ lengthcut, xlab="", main="by Total length (mm)", ylab=names(dtab)[imetric], las=2)
}
for(imetric in c(6,3:5)) { # this was just a logical order of columns
  boxplot(dtab[,imetric] ~ uprivercut, xlab="", main="by Upriver position (km)", ylab=names(dtab)[imetric], las=2)
}
par(mar=parmar)




### bringing back the time series plot (line 412ish), but also showing 
### by-individual stuff from this section 
make_a_ts <- function(colvec=rep(1,309), lwdvec=rep(1,309), alpha=c(.5,.3), ...) {
  par(mar=c(5.1, 6.1, 4.1, 2.1))  # tweaking margin for this plot
  plot(NA, ylim=c(max(daym, na.rm=T), min(daym, na.rm=T)), xlim=range(upm1, na.rm=T),
       xlab="Upriver position (km)", yaxt='n', ylab="", ...=...)
  for(i in 1:nrow(upm1)) {
    points(upm1[i,], daym[i,], col=adjustcolor(colvec[i], alpha.f=alpha[1]))
    lines(upm1[i,], daym[i,], col=adjustcolor(colvec[i], alpha.f=alpha[1]), lwd=lwdvec[i])
    lines(upm1[i,!is.na(upm1[i,])], daym[i,][!is.na(upm1[i,])], lty=3, 
          col=adjustcolor(colvec[i], alpha.f=alpha[2]), lwd=lwdvec[i])
  }
  axis(side=2, at=as.numeric(mnDates[-1]), mnDates[-1], las=2)
  # note: axis ticks are not shown for tagging dates
  par(mar=parmar)  # resetting margins to default state
}

par(mfrow=c(1,1))
# make_a_ts()
make_a_ts(colvec=1+as.numeric(as.factor(sectionmode)), alpha=c(.8,.4), 
          main="River section most often occupied")
make_a_ts(colvec=1+as.numeric(as.factor(telem3d$tagging_location[,1])), alpha=c(.8,.4), 
          main="Section tagged")
make_a_ts(colvec=1+as.numeric(as.factor(telem3d$tagging_habitat[,1])), alpha=c(.8,.4), 
          main="Tagging habitat")
make_a_ts(colvec=1+as.numeric(as.factor(telem3d$life_history[,1])), alpha=c(.8,.4), 
          main="Life history")
make_a_ts(lwdvec=as.numeric(as.factor(lengthcut)),
          main="Length")
make_a_ts(lwdvec=as.numeric(cut(dtab$homerange,3)),
          main="Homerange")
make_a_ts(lwdvec=as.numeric(cut(dtab$cumuldist,4)),
          main="Cumulative distance")
make_a_ts(lwdvec=as.numeric(cut(dtab$dist_per_obs,3)),
          main="Cumulative distance per survey")



### Interactive plot to highlight individual!!  Not for the report, but might
### be helpful for storytelling

### IMPORTANT!! - If you close the app but it keeps running and you are prompted 
### to restart R, click No, then click on the console and hit <Esc> 

library(shiny)
server <- shinyServer(function(input, output) {
  output$thePlot <- renderPlot(
    make_a_ts(lwdvec = 1+4*((1:309)==input$indiv), 
              colvec = 2-((1:309)==input$indiv), 
              alpha=c(.8,.6), main=telem3d$unique_id[input$indiv,1]))
  
})
ui <- shinyUI(fluidPage(
  titlePanel("Add Title Here"),
  sidebarLayout(
    sidebarPanel(numericInput("indiv","Individual:",min=1,max=309,value=1)),
    
    #beginning of main section
    mainPanel(plotOutput("thePlot", height="800px", width="700px"))
  )
))
shinyApp(ui = ui, server = server)



### investigating a hierarchical clustering algorithm to group similar fish
### Are discrete movement/behavior patterns detected?
hc <- hclust(dist(upm1))  # could take first column out
plot(hc)
plot(rev(hc$height))   # plotting the height of each split
plot(-diff(rev(hc$height)), log="y", type="b", xlim=c(0,50))  # difference in heights (log scale)

grp <- cutree(hc, k=14)
plot(hc)
rect.hclust(hc, k=14, border=1:14)
make_a_ts(colvec = jagshelper::rcolors(100)[grp], lwdvec = 2)


### The same thing as an interactive Shiny app, in which the user gets to choose
### the number of clusters.  Panels indicating "- NEW" in the title are the groups
### that were split at this level.  Interpreting these results, one can ask
### "is this a behaviorally meaningful split?"  I think the split at k=14 is 
### actually quite meaningful, fwiw.

### I would be curious if the groups identified have any relationship with 
### known variables (life history, size, etc)

library(shiny)
server <- shinyServer(function(input, output) {
  output$thePlot <- renderPlot(
    {
      grp <- cutree(hc, k=input$nclust)
      thetab <- table(cutree(hc, k=input$nclust-1), grp)
      newvec <- rep("", ncol(thetab))
      for(j in 1:nrow(thetab)) {
        if(sum(thetab[j,] > 0) > 1) newvec[thetab[j,] > 0] <- "- NEW"
      }
      par(mfrow=c(4,5))
      plot(hc)
      rect.hclust(hc, k=input$nclust, border=rainbow(input$nclust))
      for(i in 1:input$nclust) {
        make_a_ts(lwdvec = 1+2*(grp==i), 
                  colvec = 2-(grp==i), 
                  alpha=c(.8,.6),
                  main=paste("n =", sum(grp==i), newvec[i]))
        }
    })
  
})
ui <- shinyUI(fluidPage(
  titlePanel("Add Title Here"),
  sidebarLayout(
    sidebarPanel(sliderInput("nclust","Number of Clusters:",min=2,max=19,value=14)),
    
    #beginning of main section
    mainPanel(plotOutput("thePlot", height="1000px", width="1200px"))
  )
))
shinyApp(ui = ui, server = server)


# plot by itself for the k=14 clustering result
par(mfrow=c(4,4))
plot(hc)
rect.hclust(hc, k=14, border=rainbow(14))
for(i in 1:14) {
  make_a_ts(lwdvec = 1+2*(grp==i), 
            colvec = 2-(grp==i), 
            alpha=c(.8,.6),
            main=paste("n =", sum(grp==i)))
}

par(mfrow=c(1,1))  # making a plot of group 2 by itself to address a discussion question
for(i in 2) {
  make_a_ts(lwdvec = 1+2*(grp==i), 
            colvec = 2-(grp==i), 
            alpha=c(.8,.6),
            main=paste("n =", sum(grp==i)))
}



### Trying a suite of multiple regression models, with each distance metric as
### the response variable and logical combinations of by-individual variables
### as predictors.
### From the exploratory plots, it's apparent that there are stories to tell -
### my hope is that MLR might isolate what the cleanest stories are.

# First pulling out potential explanatory variables and grouping them by set
# note: "lake" was recoded as "tributary" since sample sizes were small for lake

tagging_loc <- telem3d$tagging_location[,1]
# avgupriver
# uprivercut
# sectionmode

tagging_hab <- telem3d$tagging_habitat[,1]
tagging_hab[tagging_hab=="lake"] <- "tributary"

life_hist <- telem3d$life_history[,1]
life_hist[life_hist=="lake"] <- "tributary"

total_length <- telem3d$total_length[,1]
# lengthcut


### Statistically, this might not be the best practice; programmatically, it's
### so cool I can hardly stand it.
### The algorithm below fits multiple regression models for all possible 
### combinations of each set of explanatory variables:
### - v1 are possible variables for general location
### - v2 are possible variables for life history
### - v3 are possible variables for length.

v1 <- c("tagging_loc","avgupriver","uprivercut","sectionmode")
v2 <- c("tagging_hab","life_hist")
v3 <- c("total_length","lengthcut")
o1 <- c("+","*")

# # including all possible interactions
# alltheformulas <- c(expand.grid(v1,o1,v2,o1,v3) %>%
#   apply(1, paste, collapse=" "),
#   expand.grid(v1,o1,v2) %>%
#     apply(1, paste, collapse=" "),
# expand.grid(v1,o1,v3) %>%
#   apply(1, paste, collapse=" "),
# expand.grid(v2,o1,v3) %>%
#   apply(1, paste, collapse=" "),
# v1, v2, v3) %>%
#   paste("y ~", .) 

# # only including main effects
# alltheformulas <- c(expand.grid(v1,v2,v3) %>%
#   apply(1, paste, collapse=" + "),
#   expand.grid(v1,v2) %>%
#     apply(1, paste, collapse=" + "), 
#   expand.grid(v1,v3) %>%
#     apply(1, paste, collapse=" + "), 
#   expand.grid(v2,v3) %>%
#     apply(1, paste, collapse=" + "),
#   v1, v2, v3) %>%
#   paste("y ~", .) 

### let's see if I can generalize this - SUCCESS!!
### Input will be a list of variable names
### NOTE: I would like to develop this further, but should probably move forward
### with the project!
makealltheformulas <- function(varlist, interactions=F) {
  whichvars <- list()
  for(ii in seq_along(varlist)) whichvars[[ii]] <- combn(seq_along(varlist), ii)
  
  # incremental element 1
  thing1 <- sapply(whichvars, function(y) apply(y, 2, function(x) expand.grid(varlist[x])))  # trying this
  # str(thing1)
  
  if(interactions) o1 <- c("*", "+")
  if(!interactions) o1 <- "+"
  
  thing2 <- list()  # thing2 will match list structure of thing1
  for(i1 in seq_along(thing1)) {
    thing2[[i1]] <- list()
    for(i2 in seq_along(thing1[[i1]])) {
      # thing2[[i1]][[i2]] <- thing1[[i1]][[i2]] # this is just to make sure indexing worked
      thedf <- thing1[[i1]][[i2]]
      if(ncol(thedf)>1) {
        omat <- expand.grid(replicate(ncol(thedf)-1, o1, simplify=F))  # matrix defined by operators (* or +)
        linker <- expand.grid(1:nrow(thedf), 1:nrow(omat))             # links variables and operators
        splicer <- matrix(nrow=nrow(linker), ncol=ncol(thedf)+ncol(omat))  # splices variables and operators
        for(i in 1:ncol(thedf)) {
          splicer[,2*i-1] <- as.character(thedf[,i][linker[,1]])
        }
        for(i in 1:ncol(omat)) {
          splicer[,2*i] <- as.character(omat[,i][linker[,2]])
        }
        thing2[[i1]][[i2]] <- paste("y ~", apply(splicer, 1, paste, collapse=" "))
      } else {
        thing2[[i1]][[i2]] <- paste("y ~", as.character(thedf[,1]))
      }
    }
  }
  thing3 <- unlist(thing2)
  return(thing3)
}
alltheformulas <- makealltheformulas(varlist=list(v1, v2, v3), interactions=T)



### fitting all candidate response variables to the suite of all predictors...

y <- log(dtab$homerange+1)

allthelms <- sapply(alltheformulas, function(x) lm(as.formula(x)))
alltheAICs <- sapply(allthelms, AIC)
par(mfrow=c(2,3))
plot(alltheAICs)
abline(h=min(alltheAICs)+2, lty=2)
alltheAICs[alltheAICs < min(alltheAICs)+2]
lm1 <- allthelms[[which.min(alltheAICs)]]
AIC(lm1)
plot(lm1)
anova(lm1)
summary(lm1)


y <- log(dtab$cumuldist+1)

allthelms <- sapply(alltheformulas, function(x) lm(as.formula(x)))
alltheAICs <- sapply(allthelms, AIC)
par(mfrow=c(2,3))
plot(alltheAICs)
abline(h=min(alltheAICs)+2, lty=2)
alltheAICs[alltheAICs < min(alltheAICs)+2]
lm1 <- allthelms[[which.min(alltheAICs)]]
AIC(lm1)
plot(lm1)
anova(lm1)
summary(lm1)


y <- log(dtab$dist_per_obs+1)  ### this is by far the cleanest in terms of diagnostics

allthelms <- sapply(alltheformulas, function(x) lm(as.formula(x)))
alltheAICs <- sapply(allthelms, AIC)
par(mfrow=c(2,3))
plot(alltheAICs)
abline(h=min(alltheAICs)+2, lty=2)
data.frame(alltheAICs, 1:length(alltheAICs))[alltheAICs < min(alltheAICs)+2,]
# lm1 <- allthelms[[which.min(alltheAICs)]]
lm1 <- allthelms[[96]]
AIC(lm1)
plot(lm1)
anova(lm1)
summary(lm1)


### Fitting an equivalent Bayesian model, for the purpose of visualizing effect sizes.
### This model is parameterized differently from the model above in two ways:
### - Instead of considering main effects & interaction terms for the one-way
###   interaction between life_hist and sectionmode, a single predictor variable
###   was constructed by concatenating the two.
### - Perhaps more importantly, this model is parameterized in terms of a global
###   mean or baseline, plus effects for [length] and [life_hist & sectionmode].
###   The global mean should not be interpreted as the mean across all fish, but
###   rather across all category means.
### Note: the previous version of this model included an interaction between 
### avg upriver position and life history, and model code should still be updated.

# specify model, which is written to a temporary file
dist_jags <- tempfile()
cat('model {
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0 + b_length[lengthcut[i]] + b_uplife[upriver_lifehist[i]]
  }

  for(j_length in 1:(nlength-1)) {
    b_length[j_length] ~ dnorm(0, 0.001)
  }
  b_length[nlength] <- -sum(b_length[1:(nlength-1)])

  for(j_uplife in 1:(nuplife-1)) {
    b_uplife[j_uplife] ~ dnorm(0, 0.001)
  }
  b_uplife[nuplife] <- -sum(b_uplife[1:(nuplife-1)])

  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)
  b0 ~ dnorm(0, 0.001)
}', file=dist_jags)


# bundle data to pass into JAGS
dist_data <- list(lengthcut=as.numeric(as.factor(lengthcut[!is.na(y) & !is.na(lengthcut)])),
                  # upriver_lifehist=as.numeric(as.factor(paste(uprivercut, life_hist)[!is.na(y) & !is.na(lengthcut)])),
                  # upriver_lifehist=as.numeric(as.factor(paste(life_hist, uprivercut)[!is.na(y) & !is.na(lengthcut)])),
                  upriver_lifehist=as.numeric(as.factor(paste(life_hist, sectionmode)[!is.na(y) & !is.na(lengthcut)])),
                  y=y[!is.na(y) & !is.na(lengthcut)])
# dist_data$lengthcut <- as.numeric(as.factor(dist_data$lengthcut_raw))
dist_data$nlength <- length(unique(dist_data$lengthcut))
# dist_data$upriver_lifehist <- as.numeric(as.factor(dist_data$upriver_lifehist_raw))
dist_data$nuplife <- length(unique(dist_data$upriver_lifehist))
dist_data$n <- length(dist_data$y)

# JAGS controls
niter <- 10000
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

{  ## running JAGS
  tstart <- Sys.time()
  print(tstart)
  dist_jags_out <- jagsUI::jags(model.file=dist_jags, data=dist_data,
                                parameters.to.save=c("b0","b_length","b_uplife","sig"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}
## a few JAGS diagnostic plots
# nbyname(dist_jags_out)
# plotRhats(dist_jags_out)
# traceworstRhat(dist_jags_out)

## pulling out data.frames associated with posteriors of length & life history/position effects
b_length <- jags_df(dist_jags_out, p="b_length")
b_uplife <- jags_df(dist_jags_out, p="b_uplife")


## Plotting effect sizes!!

par(mar=c(10.1, 4.1, 2.1, 2.1))  ## tweaking margins so plots below will work

par(mfrow=c(2,2))
caterpillar(b_length, xax=rep("",4), main="Effect sizes for length category")
axis(1, 1:4, label=levels(as.factor(lengthcut)), las=2)
abline(h=0, lty=3)
caterpillar(b_uplife, xax=rep("",ncol(b_uplife)), main="Effect sizes for section mode & life history")
# axis(1, 1:ncol(b_uplife), label=levels(as.factor(paste(life_hist,uprivercut))), las=2)
axis(1, 1:ncol(b_uplife), label=levels(as.factor(paste(life_hist,sectionmode))), las=2)
abline(h=0, lty=3)
caterpillar(exp(b_length), xax=rep("",4), 
            main="Exponentiated effect sizes for length category",
            ylab="Multiplicative change from baseline")
axis(1, 1:4, label=levels(as.factor(lengthcut)), las=2)
abline(h=1, lty=3)
caterpillar(exp(b_uplife), xax=rep("",ncol(b_uplife)), 
            main="Exponentiated effect sizes for section mode & life history",
            ylab="Multiplicative change from baseline")
# axis(1, 1:ncol(b_uplife), label=levels(as.factor(paste(life_hist,uprivercut))), las=2)
axis(1, 1:ncol(b_uplife), label=levels(as.factor(paste(life_hist,sectionmode))), las=2)
abline(h=1, lty=3)

## Looking at how much of the effect size posteriors are above/below zero
## This can be taken as evidence of which effects are different from zero 
colMeans(b_length>0)
colMeans(b_uplife>0)

## can interpret these as p-values(ish)
sapply(colMeans(b_length>0), function(x) 2*min(x,1-x))
sapply(colMeans(b_uplife>0), function(x) 2*min(x,1-x))


## Boxplots of raw data, broken out by factors identified by the model

# boxplot(y ~ lengthcut, main="log(dist per obs) by length category", 
#         las=2, xlab="")
# boxplot(y ~ paste(life_hist,sectionmode), main="log(dist per obs) by section mode & life history", 
#         las=2, xlab="")
boxplot(dtab$dist_per_obs+1 ~ lengthcut, main="Dist per obs, by length category", 
        las=2, xlab="", log="y")
boxplot(dtab$dist_per_obs+1 ~ paste(life_hist, sectionmode), main="Dist per obs, by section mode & life history", 
        las=2, xlab="", log="y")
boxplot(dtab$dist_per_obs ~ lengthcut, main="Dist per obs, by length category", 
        las=2, xlab="")
boxplot(dtab$dist_per_obs ~ paste(life_hist, sectionmode), main="Dist per obs, by section mode & life history", 
        las=2, xlab="")

par(mar=parmar)  # re-setting margins


### making time-series plots to show up/downriver movement
### broken out by life history / river section
par(mfrow=c(2,3))
for(i in 1:6) {
make_a_ts(lwdvec=1+2*(as.numeric(as.factor(paste(life_hist, sectionmode)))==i),
          colvec=2-1*(as.numeric(as.factor(paste(life_hist, sectionmode)))==i),
          main=names(table(paste(life_hist, sectionmode)))[i])
}



### STARTING TO LOOK AT DISTANCES BETWEEN SEQUENTIAL SURVEYS
# Just looking at sequential pairs - this does not account for individuals
# that were missed in a given flight
# Second and third capture events were included as part of flights that happened 
# at approximately the same time: flights 3 and 5

# data stuff to prepare for calculation
tag_breaks <- as.Date(c("2018-09-01","2019-01-01","2019-09-01","2020-01-01"))
# table(telem_alive$date[telem_alive$flight_num==0.1])
# table(telem_alive$date, telem_alive$flight_num)
telem_alive <- telemdata %>% filter(current_state!="dead")
telem_alive$flight_num[telem_alive$flight_num==0.1 & telem_alive$date>tag_breaks[2] & telem_alive$date<tag_breaks[3]] <- 3  # for plot
telem_alive$flight_num[telem_alive$flight_num==0.1 & telem_alive$date>tag_breaks[3] & telem_alive$date<tag_breaks[4]] <- 5  # for plot

# - table of distances between sequential surveys
sequentialdists <- riverdistanceseq(unique=telem_alive$unique_id_num,
                                    survey=telem_alive$flight_num,
                                    seg=telem_alive$seg,
                                    vert=telem_alive$vert,
                                    rivers=tyb_trim)/1000 # make it km

# - table of upstream distances between sequential surveys
upstreamdists <- upstreamseq(unique=telem_alive$unique_id_num,
                                    survey=telem_alive$flight_num,
                                    seg=telem_alive$seg,
                                    vert=telem_alive$vert,
                                    rivers=tyb_trim)/1000 # make it km

## plotting distances and upstream distances
par(mfrow=c(1,2))
plotseq(sequentialdists, ylab="River distance (km)")
plotseq(upstreamdists, ylab="Upstream rive distance (km)")

## making a table of n/mean/SD/SE
mean_tbl <- function(x) {
  n <- apply(!is.na(x), 2, sum)
  mn <- round(apply(x, 2, mean, na.rm=T), 2)
  SD <- round(apply(x, 2, sd, na.rm=T), 2)
  SE <- round(SD/sqrt(n), 2)
  return(cbind(n,mn,SD,SE))
}
mean_tbl(sequentialdists)
kable(mean_tbl(sequentialdists))
mean_tbl(upstreamdists)
kable(mean_tbl(upstreamdists))



### overlap of spatial extent for each stock (defined as tagged section)
### for each survey date
### THIS IS DEFINED AS INTERSECTION/UNION FOR MINIMUM SPATIAL EXTENT

## calculating and plotting for each flight

overlaps <- list()
flights <- sort(unique(telem_alive$flight_num))
for(i in 2:length(flights)) {
  hri <- with(subset(telem_alive, flight_num==flights[i]),
              homerange(unique=tagging_location, seg=seg, vert=vert, rivers=tyb_trim))
  par(mfrow=c(2,2))
  plot(NA, xlim=0:1, ylim=0:1, xlab="", ylab="",main="")
  text(.5, .7, paste("flight",flights[i]), cex=2)
  text(.5, .4, mnDates[i], cex=2)
  plot(hri)
  overlaps[[i-1]] <- homerangeoverlap(hri)$prop_both
}
overlaps_tbl <- data.frame(LowerMiddle=round(sapply(overlaps, function(x) x[1,2]),3),
                           MiddleUpper=round(sapply(overlaps, function(x) x[2,3]),3),
                           LowerUpper=round(sapply(overlaps, function(x) x[1,3]),3))
rownames(overlaps_tbl) <- flights[-1]

# tables
overlaps_tbl          # print to console
kable(overlaps_tbl)   # print to Markdown

# and a plot to go with the table
par(mfrow=c(1,1))
plot(overlaps_tbl$LowerMiddle, xaxt="n", ylab="Percent Overlap", xlab="Flight number", 
     pch=18, col=2, ylim=c(0,.6))
lines(overlaps_tbl$LowerMiddle, col=2)
points(overlaps_tbl$MiddleUpper, pch=15, col=3)
lines(overlaps_tbl$MiddleUpper, col=3)
points(overlaps_tbl$LowerUpper, pch=16, col=4)
lines(overlaps_tbl$LowerUpper, col=4)
axis(side=1, labels=flights[-1], at=seq_along(flights[-1]))
legend("topleft", pch=c(18,15,16), lty=c(1,1,1), col=2:4, 
       legend=c("Lower - Middle","Middle - Upper", "Lower - Upper"))



# - percent overlap between each pair of surveys (homerange) - might not do this





#### survival analysis thing!!
#### I developed a nifty Bayesian model to estimate survival probabilities
#### (1-mortality probability) for each time period between flights.
#### The cool thing about this model is it treats survival (1=alive, 0=dead) as
#### an unknown parameter when it is not observed for an individual in a given
#### flight. Mortality then propagates sequentially (once an individual is dead,
#### it stays dead).

#### Before defining and running the model, there was a lot of data stuff.
#### The main input to the model was a survival matrix (ones and zeroes) in
#### which each row was an individual, and each column was a flight.  Additional
#### ones and zeroes were filled in, given known information: if an individual
#### was observed to be dead in a given flight, it was considered to be dead
#### at every sequential time period.  Conversely, if an individual was observed
#### to be alive, it was considered alive at every time period prior, following 
#### the date when it was captured.  All unknown values of survival were coded as NA.

#### Note: for this analysis, I wanted to be a bit more explicit with the dates,
#### so I inserted columns (dates) associated with capture events 2 and 3 rather
#### than considering them to be concurrent with flights.


## DATA STRATEGY:
# - insert columns corresponding to dates of tag deployment, insert alive then
# - define a vector of flight date when first present
# -    second tag deployment is basically concurrent with 3rd flight: say 4th flight
# -    third tag deployment slightly before 5th flight: say 6th flight
tagdate <- as.Date(as.character(telem3d$tagging_date[,1]), format="%m/%d/%Y")
tagdatecut <- as.numeric(cut(tagdate, as.Date(c("2018-09-01","2019-01-01","2019-09-01","2020-01-01"))))
# table(tagdate, tagdatecut)  # check that this cutting strategy worked

# define survival table
# - each row is an individual
# - each column is a flight event (not counting tagging)
# - 0=alive, 1=dead, NA=not found - THIS COMES FROM RAW DATA SPREADSHEET (convert later)
survtable <- telem3d$current_state_num[,-1]            # take out tagging event
survtable <- cbind(ifelse(tagdatecut==1, 0, NA), survtable)  # insert survival corresponding to first tagging event
survtable <- cbind(survtable[,1:3],
                   ifelse(tagdatecut==2, 0, NA),
                   survtable[,4:ncol(survtable)])
survtable <- cbind(survtable[,1:6],
                   ifelse(tagdatecut==3, 0, NA),
                   survtable[,7:ncol(survtable)])

# this will return flight first fully present
# that is, the first flight AFTER entry 
firstpresent <- c(2, 5, 8)[tagdatecut]  

# take out indivs that were never found
# firstpresent <- firstpresent[rowMeans(is.na(survtable))<1]
# survtable <- survtable[rowMeans(is.na(survtable))<1,]  

# define a vector as the first event each individual was encountered dead
# once an individual is dead, it stays dead: does not contribute anything to prob model
# also, back-filling "alive" when an individual is found in later survey
firstdead <- firstalive <- rep(NA, nrow(survtable))
for(i in 1:nrow(survtable)) {
  founddead <- F
  j <- 1
  while(!founddead) {
    if(!is.na(survtable[i,j])) {
      if(survtable[i,j]==1) {
        founddead<-T
        firstdead[i] <- j
        survtable[i,j:ncol(survtable)] <- 1  # once it's dead, it stays dead
      }
      if(survtable[i,j]==0) {
        if(is.na(firstalive[i])) firstalive[i] <- j
        survtable[i, firstalive[i]:j] <- 0   # backfilling alive to point of entry
      }
    }
    j <- j+1
    if(j>ncol(survtable)) founddead <- T
  }
}
firstdead[is.na(firstdead)] <- ncol(survtable)  # necessary for Bayesian model to function

survtable <- 1-survtable ## make 0=dead, 1=alive !!!!!!


# Bayesian model
# - Of interest is parameter vector p[], which represents the probability of
#   mortality for each sequential pair of surveys.
# - If an individual is not observed (has value NA) in a given survey, the unobserved
#   state will be treated as an unknown parameter

# skeleton("surv")
# specify model
surv_jags <- tempfile()
cat('model {
  for(i in 1:n) {  
    for(j in firstpresent[i]:firstdead[i]) {          # for each survey
      survtable[i,j] ~ dbin(p[j-1], survtable[i,j-1])   # for each event present   
    }
  }
  for(j in 1:np) {
    p[j] ~ dbeta(.5,.5)
  }
}', file=surv_jags)

# bundle data to pass into JAGS
surv_data <- list(survtable=survtable,  
                  firstdead=firstdead,
                  firstpresent=firstpresent,
                  n=nrow(survtable),
                  np=ncol(survtable)-1)

# JAGS controls
niter <- 100000
ncores <- min(10, parallel::detectCores()-1)  # number of cores to use

{
  tstart <- Sys.time()
  print(tstart)
  surv_jags_out <- jagsUI::jags(model.file=surv_jags, data=surv_data,
                                parameters.to.save=c("p","survtable"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

# checking output & assessing convergence (all looks good)
nbyname(surv_jags_out)     # how many parameters are returned
par(mfrow=c(2,2))
plotRhats(surv_jags_out)   # plotting Rhat for each parameter: near 1 is good
traceworstRhat(surv_jags_out)  # trace plots for worst parameters


## plot 50 & 95% credible intervals for each survival parameter p!
par(mfrow=c(1,1))
parmar <- par("mar")  # storing original margin settings
par(mar=c(6.1, 4.1, 4.1, 2.1))
caterpillar(surv_jags_out, p="p", xax=rep("", surv_data$np))  # 50 & 95% credible intervals for each p!
mnDates1 <- as.character(mnDates)
plotdates <- c("Capture 1", mnDates1[2:3],"Capture 2", mnDates1[4:5],"Capture 3",mnDates1[6:length(mnDates)])
axis(side=1, at=0:surv_data$np, plotdates, las=2)


## plotting survival & survival probability for all individuals
## 0=dead, 1=alive, somewhere in between represents PROBABILITY of survival!
## note: points/lines are jittered slightly for better visibility
par(mfrow=c(1,1))
par(mar=c(6.1, 4.1, 4.1, 2.1))
plot(NA, xlim=c(1,ncol(survtable)), ylim=0:1,
     xlab="", ylab="Survival Probability", xaxt="n")
for(i in 1:nrow(survtable)) {
  jmeans <- surv_jags_out$mean$survtable[i,] + runif(ncol(surv_jags_out$mean$survtable), -0.01, 0.01)
  jdates <- 1:ncol(surv_jags_out$mean$survtable) + runif(ncol(surv_jags_out$mean$survtable), -0.05, 0.05)
  lines(jdates, jmeans, col=adjustcolor(1, alpha.f = .4))
  points(jdates, jmeans, col=adjustcolor(1, alpha.f = .4))
}
axis(side=1, at=1:length(plotdates), plotdates, las=2)
par(mar=parmar)


## making a table of survival parameters & CI's
pmed <- surv_jags_out$q50$p
pse <- surv_jags_out$sd$p
plo <- surv_jags_out$q2.5$p
phi <- surv_jags_out$q97.5$p
ndigits <- 3
ci95 <- paste0("(",round(plo, ndigits)," - ", round(phi, ndigits),")")
ptable <- data.frame(p_est=round(pmed, ndigits),
                     SE=round(pse, ndigits),
                     CI95=ci95)
rownames(ptable) <- plotdates[-1]
ptable   # print to console
knitr::kable(ptable)  # print to markdown



### THIS SECTION WILL PROBABLY LOOK A LITTLE WEIRD
### and is mostly included as a sanity check to validate the model and make sure
### it behaved as intended.
### The intent was to calculate "raw" binomial point estimates of survival probability
### at each time step, defined as the proportion observed to survive between sequential
### observations j-1 and j, out of all individuals alive at j-1 and observed (dead
### OR alive) at both j-1 and j.
### I also wanted to calculate an "adjusted" point estimate of survival probability
### for each time step, essentially weighted averages of the numbers observed and
### the numbers unobserved and the individual-level probabilities of survival.

### The hope was that there would be total agreement between the parameter estimates
### from the Bayesian model and the adjusted point estimates that also used output
### from the Bayesian model (this would mean the model was internally consistent),
### and to see if there were differences in inference between the Bayesian model
### and the "raw" binomial estimates.

### The raw binomial estimates were consistently higher than the Bayesian estimates!
### I haven't figured out why yet.  However, the product of all survival parameters
### from the Bayesian model agrees more with the total survival of all individuals
### from the beginning to the end of the study, so it's likely that the Bayesian
### estimates are better.  Thankfully the "adjusted" point estimates were consistent
### with the Bayesian parameter estimates from the same model!


## Tabulating all possible cases for each individual at each time step:
# - known alive at j-1, known alive at j
# - known alive at j-1, known dead at j
# - known alive at j-1, unknown at j
# - known dead at j-1, known dead at j
# - unknown at j-1, known dead at j
# - unknown at j-1, unknown at j
transition_arr <- array(dim=c(nrow(survtable), ncol(survtable)-1, 6))
for(j in 2:ncol(survtable)) {
  transition_arr[,j-1,1] <- survtable[,j-1]==1 & survtable[,j]==1 & !is.na(survtable[,j-1]) & !is.na(survtable[,j])
  transition_arr[,j-1,2] <- survtable[,j-1]==1 & survtable[,j]==0 & !is.na(survtable[,j-1]) & !is.na(survtable[,j])
  transition_arr[,j-1,3] <- survtable[,j-1]==1 & !is.na(survtable[,j-1]) & is.na(survtable[,j])
  transition_arr[,j-1,4] <- survtable[,j-1]==0 & survtable[,j]==0 & !is.na(survtable[,j-1]) & !is.na(survtable[,j])
  transition_arr[,j-1,5] <- survtable[,j]==0 & is.na(survtable[,j-1]) & !is.na(survtable[,j])
  transition_arr[,j-1,6] <- is.na(survtable[,j-1]) & is.na(survtable[,j])
}
for(i in 1:length(firstpresent)) {  # accounting for entry to the system
  if(firstpresent[i] > 2) transition_arr[i, 1:(firstpresent[i]-2), ] <- NA
}

# summary table should have:
# - number known alive at j-1, known alive at j
# - number known alive at j-1, known dead at j
# - number known alive at j-1, unknown at j (avg prob alive at j)
# - number known dead at j-1, known dead at j
# - number unknown at j-1, known dead at j (avg prob alive at j-1)
# - number unknown at j-1, unknown at j (avg prob alive at j-1, avg prob alive at j)
transition_mat <- apply(transition_arr, 2:3, sum, na.rm=T)  # raw sums
survprob <- surv_jags_out$mean$survtable  # grabbing the survival matrix from jags output
case3_avgprob_j <- case5_avgprob_j1 <- case6_avgprob_j1 <- case6_avgprob_j <- rep(NA, nrow(transition_mat))
for(j in 2:(nrow(transition_mat)+1)) {
  case3_avgprob_j[j-1] <- mean(survprob[transition_arr[,j-1,3], j], na.rm=T)
  case5_avgprob_j1[j-1] <- mean(survprob[transition_arr[,j-1,5], j-1], na.rm=T)
  case6_avgprob_j1[j-1] <- mean(survprob[transition_arr[,j-1,6], j-1], na.rm=T)
  case6_avgprob_j[j-1] <- mean(survprob[transition_arr[,j-1,6], j], na.rm=T)
}
case3_avgprob_j[is.nan(case3_avgprob_j)] <- 0
case5_avgprob_j1[is.nan(case5_avgprob_j1)] <- 0
case6_avgprob_j1[is.nan(case6_avgprob_j1)] <- 0
case6_avgprob_j[is.nan(case6_avgprob_j)] <- 0
transition_tbl <- data.frame(case1=transition_mat[,1],
                             case2=transition_mat[,2],
                             case3=paste0(transition_mat[,3], " (", round(100*case3_avgprob_j),"%)"),
                             case4=transition_mat[,4],
                             case5=paste0(transition_mat[,5], " (", round(100*case5_avgprob_j1),"%)"),
                             case6=paste0(transition_mat[,6], " (", round(100*case6_avgprob_j1),"%, ", round(100*case6_avgprob_j),"%)"),
                             p_raw=round(transition_mat[,1]/rowSums(transition_mat[,1:2]), 3),
                             p_adj=round((transition_mat[,1] + (transition_mat[,3]*case3_avgprob_j) + (transition_mat[,6]*case6_avgprob_j))/
                               (rowSums(transition_mat[,1:3]) + (transition_mat[,5]*case5_avgprob_j1) + (transition_mat[,6]*case6_avgprob_j1)), 3),
                             p_bayes=round(pmed, 3)
)
rownames(transition_tbl) <- plotdates[-1]
transition_tbl  # print to console
kable(transition_tbl)  # print to Markdown


## how do raw and adjusted probability estimates compare to Bayesian model?
# REPEATING OLD PLOT OF SURVIVAL PARAMETER p 
par(mfrow=c(1,1))
par(mar=c(6.1, 4.1, 4.1, 2.1))
caterpillar(surv_jags_out, p="p", xax=rep("", surv_data$np))  # 50 & 95% credible intervals for each p!
mnDates1 <- as.character(mnDates)
plotdates <- c("Capture 1", mnDates1[2:3],"Capture 2", mnDates1[4:5],"Capture 3",mnDates1[6:length(mnDates)])
axis(side=1, at=0:surv_data$np, plotdates, las=2)

## overlay point estimates - raw and adjusted
points(transition_tbl$p_raw)
points(transition_tbl$p_adj, pch=16)
legend("bottomright", legend=c("p raw", "p adjusted"), pch=c(1,16))


# # should look at how many days elapsed between events??
# mnTagdates <- tapply(tagdate, tagdatecut, mean) %>% 
#   round %>% 
#   unname %>%
#   as.Date(origin=as.Date("1970-01-01"))
# 
# allMnDates <- sort(c(mnDates[-1], mnTagdates)) %>% unname
# diffs <- diff(allMnDates) %>% as.numeric
# 
# plot((1-pmed)~diffs)



# #### AT THIS POINT I TRIED EXPANDING THE MODEL TO ESTIMATE SURVIVAL/MORTALITY
# #### SEPARATELY BY SECTION.  THE MODEL SEEMED TO WORK, BUT I WASN'T SURE IF
# #### THE RESULTS WERE MEANINGFUL OR INTERPRETABLE.
# 
# ### does survival/mortality differ by section?? harvest certainly does
# 
# # build section matrix from telem3d$section to correspond with survtable
# # - insert columns for capture events (take entries from original column 1)
# # - for NA, impute highest-frequency section for each row (individual)
# # - combine Lower and Yukon
# 
# # recoding
# sectiontable <- telem3d$section
# sectiontable[sectiontable=="Lower"] <- 1
# sectiontable[sectiontable=="Yukon"] <- 1
# sectiontable[sectiontable=="Middle"] <- 2
# sectiontable[sectiontable=="Upper"] <- 3
# sectiontable <- matrix(as.numeric(sectiontable), nrow=nrow(sectiontable), ncol=ncol(sectiontable))
# 
# # adding columns for capture events
# sectiontagged <- sectiontable[,1]
# sectiontable[,1] <- NA
# sectiontable <- cbind(sectiontable[,1:3], NA, sectiontable[,4:5], NA, sectiontable[,6:ncol(sectiontable)])
# sectiontable[tagdatecut==1, 1] <- sectiontagged[tagdatecut==1]
# sectiontable[tagdatecut==2, 4] <- sectiontagged[tagdatecut==2]
# sectiontable[tagdatecut==3, 7] <- sectiontagged[tagdatecut==3]
# 
# # imputing missing values
# themode <- function(x) as.numeric(names(sort(table(x), decreasing=T))[1])
# for(i in 1:nrow(sectiontable)) {
#   sectiontable[i, is.na(sectiontable[i,])] <- themode(sectiontable[i,])
# }
# 
# ## the next Bayesian model!!
# survsection_jags <- tempfile()
# cat('model {
#   for(i in 1:n) {  
#     for(j in firstpresent[i]:firstdead[i]) {          # for each survey
#       survtable[i,j] ~ dbin(p[sectiontable[i,j], j-1], survtable[i,j-1])   # for each event present   
#     }
#   }
#   for(i in 1:3) {
#     for(j in 1:np) {
#       p[i,j] ~ dbeta(.5, .5)
#     }
#   }
# }', file=survsection_jags)
# 
# # bundle data to pass into JAGS
# survsection_data <- list(survtable=survtable, 
#                          sectiontable=sectiontable,
#                   firstdead=firstdead,
#                   firstpresent=firstpresent,
#                   n=nrow(survtable),
#                   np=ncol(survtable)-1)
# 
# # JAGS controls
# niter <- 100000
# ncores <- min(10, parallel::detectCores()-1)  # number of cores to use
# 
# {
#   tstart <- Sys.time()
#   print(tstart)
#   survsection_jags_out <- jagsUI::jags(model.file=survsection_jags, data=survsection_data,
#                                 parameters.to.save=c("p","survtable"),
#                                 n.chains=ncores, parallel=T, n.iter=niter,
#                                 n.burnin=niter/2, n.thin=niter/2000)
#   print(Sys.time() - tstart)
# }
# 
# # checking output & assessing convergence (all looks good)
# nbyname(survsection_jags_out)
# par(mfrow=c(2,2))
# plotRhats(survsection_jags_out)
# traceworstRhat(survsection_jags_out)
# 
# caterpillar(survsection_jags_out, p="p", row=1)
# caterpillar(survsection_jags_out, p="p", row=2)
# caterpillar(survsection_jags_out, p="p", row=3)
# 
# comparecat(list(as.data.frame(survsection_jags_out$sims.list$p[,1,]),
#                 as.data.frame(survsection_jags_out$sims.list$p[,2,]),
#                 as.data.frame(survsection_jags_out$sims.list$p[,3,])))
# 
# # plotting survival & survival probability
# par(mfrow=c(1,1))
# par(mar=c(6.1, 4.1, 4.1, 2.1))
# plot(NA, xlim=c(1,ncol(survtable)), ylim=0:1,
#      xlab="", ylab="Survival Probability", xaxt="n")
# for(i in 1:nrow(survtable)) {
#   jmeans <- survsection_jags_out$mean$survtable[i,] + runif(ncol(survsection_jags_out$mean$survtable), -0.01, 0.01)
#   jdates <- 1:ncol(survsection_jags_out$mean$survtable) + runif(ncol(survsection_jags_out$mean$survtable), -0.05, 0.05)
#   lines(jdates, jmeans, col=adjustcolor(1, alpha.f = .4))
#   points(jdates, jmeans, col=adjustcolor(1, alpha.f = .4))
# }
# axis(side=1, at=1:length(plotdates), plotdates, las=2)
# par(mar=parmar)
# 
# ## comparing DIC scores - it seems separating by section is better.  Weird!!
# survsection_jags_out$DIC  # 1263.714
# surv_jags_out$DIC         # 1623.812






## Bayesian survival model - next version
## Investigating whether individual-level variables (size, life history, section)
## affect survival probability, globally (not looking at interactions!!)

## this is just a list of variables
# tagging_loc <- telem3d$tagging_location[,1]
# # avgupriver
# # uprivercut
# # sectionmode
# 
# tagging_hab <- telem3d$tagging_habitat[,1]
# tagging_hab[tagging_hab=="lake"] <- "tributary"
# 
# life_hist <- telem3d$life_history[,1]
# life_hist[life_hist=="lake"] <- "tributary"
# 
# total_length <- telem3d$total_length[,1]
# # lengthcut

# skeleton("surv")
# specify model
surv_vbls_jags <- tempfile()
cat('model {
  for(i in 1:n) {  
    for(j in firstpresent[i]:firstdead[i]) {          # for each survey
      # survtable[i,j] ~ dbin(p[j-1], survtable[i,j-1])   # for each event present  
      survtable[i,j] ~ dbin(p[i,j], survtable[i,j-1])
      logit(p[i,j]) <- b0[j-1] 
      # + b_section[sectionmode[i]] 
      # + b_life[life_hist[i]] 
      # + b_length[lengthcut[i]] 
    }
  }
  
  for(j in 1:np) {
    b0[j] ~ dnorm(0, 0.1)
  }
  
  for(i_section in 1:(n_section-1)) {
    b_section[i_section] ~ dnorm(0, 0.1)
  }
  b_section[n_section] <- -sum(b_section[1:(n_section-1)])
  
  for(i_life in 1:(n_life-1)) {
    b_life[i_life] ~ dnorm(0, 0.1)
  }
  b_life[n_life] <- -sum(b_life[1:(n_life-1)])
  
  for(i_length in 1:(n_length-1)) {
    b_length[i_length] ~ dnorm(0, 0.1)
  }
  b_length[n_length] <- -sum(b_length[2:(n_length-1)])
}', file=surv_vbls_jags)

# bundle data to pass into JAGS
surv_vbls_data <- list(survtable=survtable,  
                  firstdead=firstdead,
                  firstpresent=firstpresent,
                  n=nrow(survtable),
                  np=ncol(survtable)-1,
                  sectionmode=as.numeric(as.factor(sectionmode)),
                  n_section=length(unique(sectionmode)),
                  # sectionmode=as.numeric(as.factor(paste(life_hist, sectionmode))),
                  # n_section=length(unique(paste(life_hist, sectionmode))),
                  life_hist=as.numeric(as.factor(life_hist)),
                  n_life=length(unique(life_hist)),
                  lengthcut=as.numeric(as.factor(lengthcut)),
                  n_length=length(levels(lengthcut)))
surv_vbls_data$survtable <- surv_vbls_data$survtable[!is.na(lengthcut),]
surv_vbls_data$firstdead <- surv_vbls_data$firstdead[!is.na(lengthcut)]
surv_vbls_data$firstpresent <- surv_vbls_data$firstpresent[!is.na(lengthcut)]
surv_vbls_data$sectionmode <- surv_vbls_data$sectionmode[!is.na(lengthcut)]
surv_vbls_data$life_hist <- surv_vbls_data$life_hist[!is.na(lengthcut)]
surv_vbls_data$lengthcut <- surv_vbls_data$lengthcut[!is.na(lengthcut)]
surv_vbls_data$n <- sum(!is.na(lengthcut))

# JAGS controls
niter <- 100*1000  # 100k takes 4.3 min
ncores <- min(10, parallel::detectCores()-1)  # number of cores to use

{
  tstart <- Sys.time()
  print(tstart)
  surv_vbls_jags_out <- jagsUI::jags(model.file=surv_vbls_jags, data=surv_vbls_data,
                                parameters.to.save=c("p","b0","b_section","b_life","b_length"), #"survtable",
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

# checking output & assessing convergence 
nbyname(surv_vbls_jags_out)     # how many parameters are returned
par(mfrow=c(2,2))
plotRhats(surv_vbls_jags_out)   # plotting Rhat for each parameter: near 1 is good
traceworstRhat(surv_vbls_jags_out)  # trace plots for worst parameters

par(mfrow=c(4,4))
tracedens_jags(surv_vbls_jags_out, p="b0")

par(mfrow=c(2,2))
caterpillar(surv_vbls_jags_out, p="b0")
abline(h=0, lty=3)
caterpillar(surv_vbls_jags_out, p="b_section")
abline(h=0, lty=3)
caterpillar(surv_vbls_jags_out, p="b_life")
abline(h=0, lty=3)
caterpillar(surv_vbls_jags_out, p="b_length")
abline(h=0, lty=3)

par(mfrow=c(2,2))
caterpillar(expit(surv_vbls_jags_out$sims.list$b0))
abline(h=0:1, lty=3)
caterpillar(exp(surv_vbls_jags_out$sims.list$b_section))
abline(h=0:1, lty=3)
caterpillar(exp(surv_vbls_jags_out$sims.list$b_life))
abline(h=0:1, lty=3)
caterpillar(exp(surv_vbls_jags_out$sims.list$b_length))
abline(h=0:1, lty=3)

## comparing DIC scores - it seems separating by section is better.  Weird!!
surv_jags_out$DIC         # 1623.812
survsection_jags_out$DIC  # 1263.714
surv_vbls_jags_out$DIC    # 1523.894, 1535.742 without length, 1563.367 with life x section interaction
# try running THIS model all the ways, REMOVING NA VALUES IN LENGTH:
# - taking out all explanatory variables              DIC = 1589.783 at 500k reps
# - life hist + section as main effects               DIC = 1524.214 at 500k reps
# - life hist + section as main effects, plus length  DIC = 1547.579 at 500k reps
# - life hist x section as interaction                DIC = 1539.344 at 500k reps
# - life hist x section as interaction, plus length   DIC = 1563.367 at 500k reps



# plotting survival & survival probability
par(mfrow=c(1,1))
par(mar=c(6.1, 4.1, 4.1, 2.1))
plot(NA, xlim=c(1,ncol(survtable)), ylim=0:1,
     xlab="", ylab="Survival Probability", xaxt="n")
for(i in 1:nrow(survtable)) {
  jmeans <- surv_vbls_jags_out$mean$survtable[i,] + runif(ncol(surv_vbls_jags_out$mean$survtable), -0.01, 0.01)
  jdates <- 1:ncol(surv_vbls_jags_out$mean$survtable) + runif(ncol(surv_vbls_jags_out$mean$survtable), -0.05, 0.05)
  lines(jdates, jmeans, col=adjustcolor(1, alpha.f = .4))
  points(jdates, jmeans, col=adjustcolor(1, alpha.f = .4))
}
axis(side=1, at=1:length(plotdates), plotdates, las=2)
par(mar=parmar)


### TRYING A MORE STRUCTURED APPROACH TO MODEL SELECTION
# bundle data to pass into JAGS
surv_vbls_data <- list(survtable=survtable,  
                       firstdead=firstdead,
                       firstpresent=firstpresent,
                       n=nrow(survtable),
                       np=ncol(survtable)-1,
                       sectionmode=as.numeric(as.factor(sectionmode)),
                       n_section=length(unique(sectionmode)),
                       lifesection=as.numeric(as.factor(paste(life_hist, sectionmode))),
                       n_lifesection=length(unique(paste(life_hist, sectionmode))),
                       life_hist=as.numeric(as.factor(life_hist)),
                       n_life=length(unique(life_hist)),
                       lengthcut=as.numeric(as.factor(lengthcut)),
                       n_length=length(levels(lengthcut)))
surv_vbls_data$survtable <- surv_vbls_data$survtable[!is.na(lengthcut),]
surv_vbls_data$firstdead <- surv_vbls_data$firstdead[!is.na(lengthcut)]
surv_vbls_data$firstpresent <- surv_vbls_data$firstpresent[!is.na(lengthcut)]
surv_vbls_data$sectionmode <- surv_vbls_data$sectionmode[!is.na(lengthcut)]
surv_vbls_data$lifesection <- surv_vbls_data$lifesection[!is.na(lengthcut)]
surv_vbls_data$life_hist <- surv_vbls_data$life_hist[!is.na(lengthcut)]
surv_vbls_data$lengthcut <- surv_vbls_data$lengthcut[!is.na(lengthcut)]
surv_vbls_data$n <- sum(!is.na(lengthcut))

# JAGS controls
niter <- 1000*1000  # 100k takes 4.3 min  - now takes 1 minute and change
ncores <- min(10, parallel::detectCores()-1)  # number of cores to use

jagsouts <- list()

## Model 1: baseline
surv_vbls_jags <- tempfile()
cat('model {
  for(i in 1:n) {  
    for(j in firstpresent[i]:firstdead[i]) {          # for each survey
      # survtable[i,j] ~ dbin(p[j-1], survtable[i,j-1])   # for each event present  
      survtable[i,j] ~ dbin(p[i,j], survtable[i,j-1])
      logit(p[i,j]) <- b0[j-1] 
      # + b_section[sectionmode[i]]
      # + b_life[life_hist[i]] 
      # + b_lifesection[lifesection[i]] 
      # + b_length[lengthcut[i]] 
    }
  }
  
  for(j in 1:np) {
    b0[j] ~ dnorm(0, 0.1)
  }
  
  # for(i_section in 1:(n_section-1)) {
  #   b_section[i_section] ~ dnorm(0, 0.1)
  # }
  # b_section[n_section] <- -sum(b_section[1:(n_section-1)])
  
  # for(i_life in 1:(n_life-1)) {
  #   b_life[i_life] ~ dnorm(0, 0.1)
  # }
  # b_life[n_life] <- -sum(b_life[1:(n_life-1)])
  
  # for(i_lifesection in 1:(n_lifesection-1)) {
  #   b_lifesection[i_lifesection] ~ dnorm(0, 0.1)
  # }
  # b_lifesection[n_lifesection] <- -sum(b_lifesection[1:(n_lifesection-1)])
  
  # for(i_length in 1:(n_length-1)) {
  #   b_length[i_length] ~ dnorm(0, 0.1)
  # }
  # b_length[n_length] <- -sum(b_length[2:(n_length-1)])
}', file=surv_vbls_jags)

{
  tstart <- Sys.time()
  print(tstart)
  surv_vbls_jags_out <- jagsUI::jags(model.file=surv_vbls_jags, data=surv_vbls_data,
                                     parameters.to.save=c("p","b0","b_section","b_life","b_lifesection","b_length"), #"survtable",
                                     n.chains=ncores, parallel=T, n.iter=niter,
                                     n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}
jagsouts[[1]] <- surv_vbls_jags_out



## Model 2: baseline + section + life
surv_vbls_jags <- tempfile()
cat('model {
  for(i in 1:n) {  
    for(j in firstpresent[i]:firstdead[i]) {          # for each survey
      # survtable[i,j] ~ dbin(p[j-1], survtable[i,j-1])   # for each event present  
      survtable[i,j] ~ dbin(p[i,j], survtable[i,j-1])
      logit(p[i,j]) <- b0[j-1] 
      + b_section[sectionmode[i]]
      + b_life[life_hist[i]]
      # + b_lifesection[lifesection[i]] 
      # + b_length[lengthcut[i]] 
    }
  }
  
  for(j in 1:np) {
    b0[j] ~ dnorm(0, 0.1)
  }
  
  for(i_section in 1:(n_section-1)) {
    b_section[i_section] ~ dnorm(0, 0.1)
  }
  b_section[n_section] <- -sum(b_section[1:(n_section-1)])
  
  for(i_life in 1:(n_life-1)) {
    b_life[i_life] ~ dnorm(0, 0.1)
  }
  b_life[n_life] <- -sum(b_life[1:(n_life-1)])
  
  # for(i_lifesection in 1:(n_lifesection-1)) {
  #   b_lifesection[i_lifesection] ~ dnorm(0, 0.1)
  # }
  # b_lifesection[n_lifesection] <- -sum(b_lifesection[1:(n_lifesection-1)])
  
  # for(i_length in 1:(n_length-1)) {
  #   b_length[i_length] ~ dnorm(0, 0.1)
  # }
  # b_length[n_length] <- -sum(b_length[2:(n_length-1)])
}', file=surv_vbls_jags)

{
  tstart <- Sys.time()
  print(tstart)
  surv_vbls_jags_out <- jagsUI::jags(model.file=surv_vbls_jags, data=surv_vbls_data,
                                     parameters.to.save=c("p","b0","b_section","b_life","b_lifesection","b_length"), #"survtable",
                                     n.chains=ncores, parallel=T, n.iter=niter,
                                     n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}
jagsouts[[2]] <- surv_vbls_jags_out





## Model 3: baseline + section + life + length
surv_vbls_jags <- tempfile()
cat('model {
  for(i in 1:n) {  
    for(j in firstpresent[i]:firstdead[i]) {          # for each survey
      # survtable[i,j] ~ dbin(p[j-1], survtable[i,j-1])   # for each event present  
      survtable[i,j] ~ dbin(p[i,j], survtable[i,j-1])
      logit(p[i,j]) <- b0[j-1] 
      + b_section[sectionmode[i]]
      + b_life[life_hist[i]]
      # + b_lifesection[lifesection[i]] 
      + b_length[lengthcut[i]]
    }
  }
  
  for(j in 1:np) {
    b0[j] ~ dnorm(0, 0.1)
  }
  
  for(i_section in 1:(n_section-1)) {
    b_section[i_section] ~ dnorm(0, 0.1)
  }
  b_section[n_section] <- -sum(b_section[1:(n_section-1)])
  
  for(i_life in 1:(n_life-1)) {
    b_life[i_life] ~ dnorm(0, 0.1)
  }
  b_life[n_life] <- -sum(b_life[1:(n_life-1)])
  
  # for(i_lifesection in 1:(n_lifesection-1)) {
  #   b_lifesection[i_lifesection] ~ dnorm(0, 0.1)
  # }
  # b_lifesection[n_lifesection] <- -sum(b_lifesection[1:(n_lifesection-1)])
  
  for(i_length in 1:(n_length-1)) {
    b_length[i_length] ~ dnorm(0, 0.1)
  }
  b_length[n_length] <- -sum(b_length[2:(n_length-1)])
}', file=surv_vbls_jags)

{
  tstart <- Sys.time()
  print(tstart)
  surv_vbls_jags_out <- jagsUI::jags(model.file=surv_vbls_jags, data=surv_vbls_data,
                                     parameters.to.save=c("p","b0","b_section","b_life","b_lifesection","b_length"), #"survtable",
                                     n.chains=ncores, parallel=T, n.iter=niter,
                                     n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}
jagsouts[[3]] <- surv_vbls_jags_out





## Model 4: baseline + section x life 
surv_vbls_jags <- tempfile()
cat('model {
  for(i in 1:n) {  
    for(j in firstpresent[i]:firstdead[i]) {          # for each survey
      # survtable[i,j] ~ dbin(p[j-1], survtable[i,j-1])   # for each event present  
      survtable[i,j] ~ dbin(p[i,j], survtable[i,j-1])
      logit(p[i,j]) <- b0[j-1] 
      # + b_section[sectionmode[i]]
      # + b_life[life_hist[i]]
      + b_lifesection[lifesection[i]]
      # + b_length[lengthcut[i]]
    }
  }
  
  for(j in 1:np) {
    b0[j] ~ dnorm(0, 0.1)
  }
  
  # for(i_section in 1:(n_section-1)) {
  #   b_section[i_section] ~ dnorm(0, 0.1)
  # }
  # b_section[n_section] <- -sum(b_section[1:(n_section-1)])
  
  # for(i_life in 1:(n_life-1)) {
  #   b_life[i_life] ~ dnorm(0, 0.1)
  # }
  # b_life[n_life] <- -sum(b_life[1:(n_life-1)])
  
  for(i_lifesection in 1:(n_lifesection-1)) {
    b_lifesection[i_lifesection] ~ dnorm(0, 0.1)
  }
  b_lifesection[n_lifesection] <- -sum(b_lifesection[1:(n_lifesection-1)])
  
  # for(i_length in 1:(n_length-1)) {
  #   b_length[i_length] ~ dnorm(0, 0.1)
  # }
  # b_length[n_length] <- -sum(b_length[2:(n_length-1)])
}', file=surv_vbls_jags)

{
  tstart <- Sys.time()
  print(tstart)
  surv_vbls_jags_out <- jagsUI::jags(model.file=surv_vbls_jags, data=surv_vbls_data,
                                     parameters.to.save=c("p","b0","b_section","b_life","b_lifesection","b_length"), #"survtable",
                                     n.chains=ncores, parallel=T, n.iter=niter,
                                     n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}
jagsouts[[4]] <- surv_vbls_jags_out





## Model 5: baseline + section x life + length
surv_vbls_jags <- tempfile()
cat('model {
  for(i in 1:n) {  
    for(j in firstpresent[i]:firstdead[i]) {          # for each survey
      # survtable[i,j] ~ dbin(p[j-1], survtable[i,j-1])   # for each event present  
      survtable[i,j] ~ dbin(p[i,j], survtable[i,j-1])
      logit(p[i,j]) <- b0[j-1] 
      # + b_section[sectionmode[i]]
      # + b_life[life_hist[i]]
      + b_lifesection[lifesection[i]]
      + b_length[lengthcut[i]]
    }
  }
  
  for(j in 1:np) {
    b0[j] ~ dnorm(0, 0.1)
  }
  
  # for(i_section in 1:(n_section-1)) {
  #   b_section[i_section] ~ dnorm(0, 0.1)
  # }
  # b_section[n_section] <- -sum(b_section[1:(n_section-1)])
  
  # for(i_life in 1:(n_life-1)) {
  #   b_life[i_life] ~ dnorm(0, 0.1)
  # }
  # b_life[n_life] <- -sum(b_life[1:(n_life-1)])
  
  for(i_lifesection in 1:(n_lifesection-1)) {
    b_lifesection[i_lifesection] ~ dnorm(0, 0.1)
  }
  b_lifesection[n_lifesection] <- -sum(b_lifesection[1:(n_lifesection-1)])
  
  for(i_length in 1:(n_length-1)) {
    b_length[i_length] ~ dnorm(0, 0.1)
  }
  b_length[n_length] <- -sum(b_length[2:(n_length-1)])
}', file=surv_vbls_jags)

{
  tstart <- Sys.time()
  print(tstart)
  surv_vbls_jags_out <- jagsUI::jags(model.file=surv_vbls_jags, data=surv_vbls_data,
                                     parameters.to.save=c("p","b0","b_section","b_life","b_lifesection","b_length"), #"survtable",
                                     n.chains=ncores, parallel=T, n.iter=niter,
                                     n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}
jagsouts[[5]] <- surv_vbls_jags_out


## simplifying model 2
simpsection <- ifelse(sectionmode=="Upper", 2, 1) # combining "Lower" and "Middle" 
simpsection <- simpsection[!is.na(lengthcut)]
surv_vbls_data$simpsection <- simpsection
surv_vbls_data$n_simpsection <- 2

## Model 6: baseline + section (2 levels) + life 
surv_vbls_jags <- tempfile()
cat('model {
  for(i in 1:n) {  
    for(j in firstpresent[i]:firstdead[i]) {          # for each survey
      # survtable[i,j] ~ dbin(p[j-1], survtable[i,j-1])   # for each event present  
      survtable[i,j] ~ dbin(p[i,j], survtable[i,j-1])
      logit(p[i,j]) <- b0[j-1] 
      # + b_section[sectionmode[i]]
      + b_life[life_hist[i]]
      # + b_lifesection[lifesection[i]]
      # + b_length[lengthcut[i]]
      + b_simpsection[simpsection[i]]
    }
  }
  
  for(j in 1:np) {
    b0[j] ~ dnorm(0, 0.1)
  }
  
  # for(i_section in 1:(n_section-1)) {
  #   b_section[i_section] ~ dnorm(0, 0.1)
  # }
  # b_section[n_section] <- -sum(b_section[1:(n_section-1)])
  
  for(i_life in 1:(n_life-1)) {
    b_life[i_life] ~ dnorm(0, 0.1)
  }
  b_life[n_life] <- -sum(b_life[1:(n_life-1)])
  
  # for(i_lifesection in 1:(n_lifesection-1)) {
  #   b_lifesection[i_lifesection] ~ dnorm(0, 0.1)
  # }
  # b_lifesection[n_lifesection] <- -sum(b_lifesection[1:(n_lifesection-1)])
  
  # for(i_length in 1:(n_length-1)) {
  #   b_length[i_length] ~ dnorm(0, 0.1)
  # }
  # b_length[n_length] <- -sum(b_length[2:(n_length-1)])
  
  for(i_simpsection in 1:(n_simpsection-1)) {
    b_simpsection[i_simpsection] ~ dnorm(0, 0.1)
  }
  b_simpsection[n_simpsection] <- -sum(b_simpsection[1:(n_simpsection-1)])
}', file=surv_vbls_jags)

{
  tstart <- Sys.time()
  print(tstart)
  surv_vbls_jags_out <- jagsUI::jags(model.file=surv_vbls_jags, data=surv_vbls_data,
                                     parameters.to.save=c("p","b0","b_section","b_life",
                                                          "b_lifesection","b_length","b_simpsection"), #"survtable",
                                     n.chains=ncores, parallel=T, n.iter=niter,
                                     n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}
jagsouts[[6]] <- surv_vbls_jags_out



## simplifying model 4
simplifesection <- ifelse(life_hist=="mainstem",
                          ifelse(sectionmode=="Upper", 2, 1), 3) 
table(simplifesection, sectionmode, life_hist)
simplifesection <- simplifesection[!is.na(lengthcut)]
surv_vbls_data$simplifesection <- simplifesection
surv_vbls_data$n_simplifesection <- 3

## Model 7: baseline + simplified section x life (3 levels)
surv_vbls_jags <- tempfile()
cat('model {
  for(i in 1:n) {  
    for(j in firstpresent[i]:firstdead[i]) {          # for each survey
      # survtable[i,j] ~ dbin(p[j-1], survtable[i,j-1])   # for each event present  
      survtable[i,j] ~ dbin(p[i,j], survtable[i,j-1])
      logit(p[i,j]) <- b0[j-1] 
      # + b_section[sectionmode[i]]
      # + b_life[life_hist[i]]
      # + b_lifesection[lifesection[i]]
      # + b_length[lengthcut[i]]
      # + b_simpsection[simpsection[i]]
      + b_simplifesection[simplifesection[i]]
    }
  }
  
  for(j in 1:np) {
    b0[j] ~ dnorm(0, 0.1)
  }
  
  # for(i_section in 1:(n_section-1)) {
  #   b_section[i_section] ~ dnorm(0, 0.1)
  # }
  # b_section[n_section] <- -sum(b_section[1:(n_section-1)])
  
  # for(i_life in 1:(n_life-1)) {
  #   b_life[i_life] ~ dnorm(0, 0.1)
  # }
  # b_life[n_life] <- -sum(b_life[1:(n_life-1)])
  
  # for(i_lifesection in 1:(n_lifesection-1)) {
  #   b_lifesection[i_lifesection] ~ dnorm(0, 0.1)
  # }
  # b_lifesection[n_lifesection] <- -sum(b_lifesection[1:(n_lifesection-1)])
  
  # for(i_length in 1:(n_length-1)) {
  #   b_length[i_length] ~ dnorm(0, 0.1)
  # }
  # b_length[n_length] <- -sum(b_length[2:(n_length-1)])
  
  # for(i_simpsection in 1:(n_simpsection-1)) {
  #   b_simpsection[i_simpsection] ~ dnorm(0, 0.1)
  # }
  # b_simpsection[n_simpsection] <- -sum(b_simpsection[1:(n_simpsection-1)])
  
  for(i_simplifesection in 1:(n_simplifesection-1)) {
    b_simplifesection[i_simplifesection] ~ dnorm(0, 0.1)
  }
  b_simplifesection[n_simplifesection] <- -sum(b_simplifesection[1:(n_simplifesection-1)])
}', file=surv_vbls_jags)

{
  tstart <- Sys.time()
  print(tstart)
  surv_vbls_jags_out <- jagsUI::jags(model.file=surv_vbls_jags, data=surv_vbls_data,
                                     parameters.to.save=c("p","b0","b_section","b_life",
                                                          "b_lifesection","b_length",
                                                          "b_simpsection","b_simplifesection"), #"survtable",
                                     n.chains=ncores, parallel=T, n.iter=niter,
                                     n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}
jagsouts[[7]] <- surv_vbls_jags_out



## comparing DIC scores across models!!
sapply(jagsouts, function(x) x$DIC)
# [1] 1592.138 1518.726 1521.173 1544.973 1566.863
# [1] 1592.145 1540.600 1540.535 1543.989 1556.535 1531.136 1536.323 ??????

sapply(jagsouts, function(x) max(unlist(x$Rhat), na.rm=T))
# [1] 1.001506 1.002097 1.002676 1.001322 1.001226 1.003476 1.002978

## comparing parameter inferences
comparecat(jagsouts, p="b0")

comparecat(jagsouts, p="b_section")
comparecat(jagsouts, p="b_simpsection")
comparecat(jagsouts, p="b_life")
comparecat(jagsouts, p="b_lifesection")
comparecat(jagsouts, p="b_simplifesection")


## choose a model
## re-run with full data, if we take length out
## maybe even re-run all models without length and choose again
## caterpillar plots of baseline probs and effects
## caterpillar plots of expit(b0) and exp(effects)
## compare expit(b0) to p from previous model
## put survtable back in and do time series plot

## think about post pred, loocv, etc
## plot model diagnostic stuff anyway

### actually - looks like lower/middle mainstem vs. upper mainstem vs. tribs (3 cats)


## ok, let's build some bones
surv_vbls_jags_out <- jagsouts[[2]]

caterpillar(surv_vbls_jags_out, p="b0", xax=rep("",length(plotdates)-1))
axis(side=1, at=0:surv_data$np, plotdates, las=2)
caterpillar(surv_vbls_jags_out, p="b_section", xax=levels(as.factor(sectionmode)))
abline(h=0:1, lty=3)
caterpillar(surv_vbls_jags_out, p="b_life", xax=levels(as.factor(life_hist)))
abline(h=0:1, lty=3)

caterpillar(surv_jags_out, p="p", xax=rep("", surv_data$np), col=2)  # 50 & 95% credible intervals for each p!
# mnDates1 <- as.character(mnDates)
# plotdates <- c("Capture 1", mnDates1[2:3],"Capture 2", mnDates1[4:5],"Capture 3",mnDates1[6:length(mnDates)])
axis(side=1, at=0:surv_data$np, plotdates, las=2)

caterpillar(exp(surv_vbls_jags_out$sims.list$b_section), xax=levels(as.factor(sectionmode)))
abline(h=0:1, lty=3)
caterpillar(exp(surv_vbls_jags_out$sims.list$b_life), xax=levels(as.factor(life_hist)))
abline(h=0:1, lty=3)

# caterpillar(expit(surv_vbls_jags_out$sims.list$b0), xax=rep("",length(plotdates)-1))
# axis(side=1, at=0:surv_data$np, plotdates, las=2)

xx1 <- as.data.frame(surv_jags_out$sims.list$p)
xx2 <- as.data.frame(expit(surv_vbls_jags_out$sims.list$b0))
colnames(xx1) <- plotdates[-1]
colnames(xx2) <- plotdates[-1]
comparecat(list(xx1, xx2), col=c(2,4))
