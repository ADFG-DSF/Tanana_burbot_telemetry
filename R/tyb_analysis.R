#### Analysis script for Tanana Burbot Telemetry
#### Matt Tyers, fall 2023


#### Data & setup


# load packages
library(riverdist)   # for spatial river network analysis (note: used v0.16.0)
library(tidyverse)   # for streamlined data manipulation
library(knitr)       # for table output to Markdown
library(jagsUI)      # for running Bayesian survival analysis 
library(jagshelper)  # for plotting & model diagnostics
# Note: JAGS can be downloaded here: https://sourceforge.net/projects/mcmc-jags/


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
# all.equal(telemdata$upriver, telemdata$upriver2)


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

## make a vector of mean flight dates associated with each number
## (easier for interpretation)
mnDates <- as.Date(round(tapply(telemdata$date, telemdata$flight_num, mean)), origin="1970-01-01")
mnDates_df <- data.frame(flight_num=as.numeric(names(mnDates)), mnDate=mnDates)
telemdata <- left_join(telemdata, mnDates_df)  # join these to original data

## make 3d list!!
telem3d <- df23d(df=telemdata, x1="unique_id_num", x2="flight_num", 
                 stopifmultiple = F, tiebreaker="power")  # wow this is slow






#### Analysis!

### length composition

# now we have a way of getting records for each individual fish
# check if we can just grab the first flight number (marking event)
colMeans(!is.na(telem3d$total_length))  # all good: all rows are not NA

# grab unique lengths
lengths <- telem3d$total_length[,1]
tagging_locations <- telem3d$tagging_location[,1]

# histogram form
hist(lengths, xlab="Total Length (mm)", main="")   # insert title in main=

# one-way ANOVA: any difference by tagging location?
boxplot(lengths~tagging_locations, ylab="Total Length (mm)", xlab="Tagging Location")
lm(lengths~tagging_locations) %>% anova
## yes!! evidence of difference in mean lengths by tagging area
## fish tagged lower down were bigger

# Multiple comparison (Tukey's HSD)
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
# network_dens <- with(telemdata1, makeriverdensity(seg=seg, vert=vert, survey=mnDate,
#                                        rivers=tyb_trim, bw=20*1000))
# # plotriverdensitypoints(network_dens)  # maybe more points than needed
# par(mfrow=c(4,4))
# plot(network_dens, scalebyN=F)
# par(mfrow=c(1,1))
# plot(network_dens, points=F, scalebyN = F)  # maybe I should try animating this

# kernel density (linear)
bandwidth <- 50  # kernel bandwidth, in km (trial and error)
lin_denses <- with(telemdata1, tapply(upriver, flight_num, density, bw=bandwidth))
maxy <- max(sapply(lin_denses, function(x) max(x$y)))
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


# all in one plot 
# upm <- telem3d$upriver[,-1]  # extracting just the matrix of upriver positions minus tagging
# plotx <- 1:ncol(upm)  # survey numbers
# plotx <- mnDates[-1]  # actual dates instead

# # v1: X is date, Y is position
# plot(NA, xlim=range(plotx), ylim=range(upm, na.rm=T),
#      ylab="Upriver position (km)", xlab="Survey Date")
# for(i in 1:nrow(upm)) {
#   points(plotx, upm[i,])
#   lines(plotx, upm[i,])
#   lines(plotx[!is.na(upm[i,])], upm[i,!is.na(upm[i,])], lty=3, col=adjustcolor(1,alpha.f=.3))
# }
# 
# # v2: X is position, Y is date
# plot(NA, ylim=c(max(plotx), min(plotx)), xlim=range(upm, na.rm=T),
#      xlab="Upriver position (km)", ylab="Survey Date")
# for(i in 1:nrow(upm)) {
#   points(upm[i,], plotx)
#   lines(upm[i,], plotx)
#   lines(upm[i,!is.na(upm[i,])], plotx[!is.na(upm[i,])], lty=3, col=adjustcolor(1,alpha.f=.3))
# }

# I wonder if it's more fair to give actual dates (including/highlighting additional fish tagged)
upm1 <- telem3d$upriver
upm1[telem3d$current_state=="dead"] <- NA   ### TAKING OUT MORTALITY
daym <- telem3d$date
parmar <- par("mar")  # storing original margin settings
par(mar=c(6.1, 4.1, 4.1, 2.1))  # tweaking margin for this plot
plot(NA, xlim=range(daym, na.rm=T), ylim=range(upm1, na.rm=T),
     ylab="Upriver position (km)", xaxt='n', xlab="")
for(i in 1:nrow(upm1)) {
  points(daym[i,], upm1[i,], col=adjustcolor(1, alpha.f=.5))
  lines(daym[i,], upm1[i,], col=adjustcolor(1, alpha.f=.5))
  lines(daym[i,][!is.na(upm1[i,])], upm1[i,!is.na(upm1[i,])], lty=3, col=adjustcolor(1,alpha.f=.3))
}
axis(side=1, at=as.numeric(mnDates[-1]), mnDates[-1], las=2)

par(mar=c(5.1, 6.1, 4.1, 2.1))  # tweaking margin for this plot
plot(NA, ylim=c(max(daym, na.rm=T), min(daym, na.rm=T)), xlim=range(upm1, na.rm=T),
     xlab="Upriver position (km)", yaxt='n', ylab="")
for(i in 1:nrow(upm1)) {
  points(upm1[i,], daym[i,], col=adjustcolor(1, alpha.f=.5))
  lines(upm1[i,], daym[i,], col=adjustcolor(1, alpha.f=.5))
  lines(upm1[i,!is.na(upm1[i,])], daym[i,][!is.na(upm1[i,])], lty=3, col=adjustcolor(1,alpha.f=.3))
}
axis(side=2, at=as.numeric(mnDates[-1]), mnDates[-1], las=2)
par(mar=parmar)  # resetting margins to default state




### estimate proportions within each river section, per survey
# note: combined Yukon with Lower!!

telemdata1$section[telemdata1$section=="Yukon"] <- "Lower"
sectiontable_alive <- with(telemdata1, table(flight_num, section))
mosaicplot(sectiontable_alive, col=grey.colors(3, rev=T), main="", xlab="Flight Number", ylab="")
chisq.test(sectiontable_alive)  # weak evidence of difference in proportions

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




## describe seasonal distributions and migrations (this will be the big one)
# - homerange by individual
hr <- homerange(unique=telem_alive$unique_id_num,
                survey=telem_alive$flight_num,
                seg=telem_alive$seg,
                vert=telem_alive$vert,
                rivers=tyb_trim)
hr_table <- hr$ranges
hr_table$range <- hr_table$range/1000
hist(hr_table$range, main="", xlab="Minimum homerange (km)")
hr_table[which.max(hr_table$range),]
par(mfrow=c(1,1))
plot(tyb_trim)
riverpoints(seg=telem_alive$seg[telem_alive$unique_id_num==115],
            vert=telem_alive$vert[telem_alive$unique_id_num==115],
            rivers=tyb_trim, pch=16)
#### seems to be a bug in homerange??  maybe it's fine
# - maybe it's just a bug in plot.homerange: tries to plot for indiv 28


# what if we did total distance INSTEAD of homerange?
# - for each individual: subset, sort by date, calculate distance sequentially
indiv <- sort(unique(telem_alive$unique_id_num))
nobs <- cumuldist <- ndays <- rep(NA, length(indiv))
for(i in 1:length(indiv)) {
  # print(i)
  d1 <- telem_alive[telem_alive$unique_id_num==indiv[i],]
  d2 <- d1[order(d1$date),]
  nobs[i] <- nrow(d2)
  ndays[i] <- d2$date %>% range %>% diff %>% as.numeric
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
dtab <- data.frame(nobs, ndays, cumuldist)
rownames(dtab) <- indiv
hist(dtab$cumuldist)
plot(nobs, cumuldist)
boxplot(cumuldist~nobs)
hist(cumuldist/(nobs-1))
plot(nobs, cumuldist/(nobs-1))
boxplot(cumuldist/(nobs-1) ~ nobs)  # seems consistent enough to use as metric
boxplot(cumuldist/ndays ~ nobs)  # less consistent, one big outlier
hist(cumuldist/ndays, breaks=10)

## looking at the top mover
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
dtab$homerange <- 0
for(i in 1:nrow(hr_table)) dtab$homerange[rownames(dtab)==hr_table$ID[i]] <- hr_table$range[i]


# - movement by individual
#    - ~size, ~stock, ~avg upriver position, ~life_history, ~tagging_habitat, ~section mode

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
par(mfrow=c(2,2))
boxplot(telem3d$total_length[,1] ~ sectionmode, xlab="River section", ylab="Total length (mm)")
boxplot(telem3d$total_length[,1] ~ telem3d$life_history[,1], xlab="Life history", ylab="Total length (mm)")
boxplot(telem3d$total_length[,1] ~ telem3d$tagging_habitat[,1], xlab="Tagging habitat", ylab="Total length (mm)")
mosaicplot(table(sectionmode, telem3d$life_history[,1]), xlab="River section", ylab="Life history", main="",col=T)

# - percent overlap between each pair of surveys (homerange) - might not do this





#### survival analysis thing? might be cool
# WRITE THIS UP

## NEW NEW STRATEGY:
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
# - 0=alive, 1=dead, NA=not found
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
nbyname(surv_jags_out)
par(mfrow=c(2,2))
plotRhats(surv_jags_out)
traceworstRhat(surv_jags_out)


# plot 50 & 95% credible intervals for each p!
par(mfrow=c(1,1))
parmar <- par("mar")  # storing original margin settings
par(mar=c(6.1, 4.1, 4.1, 2.1))
caterpillar(surv_jags_out, p="p", xax=rep("", surv_data$np))  # 50 & 95% credible intervals for each p!
mnDates1 <- as.character(mnDates)
plotdates <- c("Capture 1", mnDates1[2:3],"Capture 2", mnDates1[4:5],"Capture 3",mnDates1[6:length(mnDates)])
axis(side=1, at=0:surv_data$np, plotdates, las=2)


# plotting survival & survival probability
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




# for checking a sequence of logical vectors for each time step:
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
# plot 50 & 95% credible intervals for each p!
par(mfrow=c(1,1))
par(mar=c(6.1, 4.1, 4.1, 2.1))
caterpillar(surv_jags_out, p="p", xax=rep("", surv_data$np))  # 50 & 95% credible intervals for each p!
mnDates1 <- as.character(mnDates)
plotdates <- c("Capture 1", mnDates1[2:3],"Capture 2", mnDates1[4:5],"Capture 3",mnDates1[6:length(mnDates)])
axis(side=1, at=0:surv_data$np, plotdates, las=2)

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



### does survival/mortality differ by section?? harvest certainly does

# build section matrix from telem3d$section to correspond with survtable
# - insert columns for capture events (take entries from original column 1)
# - for NA, impute highest-frequency section for each row (individual)
# - combine Lower and Yukon

# recoding
sectiontable <- telem3d$section
sectiontable[sectiontable=="Lower"] <- 1
sectiontable[sectiontable=="Yukon"] <- 1
sectiontable[sectiontable=="Middle"] <- 2
sectiontable[sectiontable=="Upper"] <- 3
sectiontable <- matrix(as.numeric(sectiontable), nrow=nrow(sectiontable), ncol=ncol(sectiontable))

# adding columns for capture events
sectiontagged <- sectiontable[,1]
sectiontable[,1] <- NA
sectiontable <- cbind(sectiontable[,1:3], NA, sectiontable[,4:5], NA, sectiontable[,6:ncol(sectiontable)])
sectiontable[tagdatecut==1, 1] <- sectiontagged[tagdatecut==1]
sectiontable[tagdatecut==2, 4] <- sectiontagged[tagdatecut==2]
sectiontable[tagdatecut==3, 7] <- sectiontagged[tagdatecut==3]

# imputing missing values
themode <- function(x) as.numeric(names(sort(table(x), decreasing=T))[1])
for(i in 1:nrow(sectiontable)) {
  sectiontable[i, is.na(sectiontable[i,])] <- themode(sectiontable[i,])
}

## the next Bayesian model!!
survsection_jags <- tempfile()
cat('model {
  for(i in 1:n) {  
    for(j in firstpresent[i]:firstdead[i]) {          # for each survey
      survtable[i,j] ~ dbin(p[sectiontable[i,j], j-1], survtable[i,j-1])   # for each event present   
    }
  }
  for(i in 1:3) {
    for(j in 1:np) {
      p[i,j] ~ dbeta(.5, .5)
    }
  }
}', file=survsection_jags)

# bundle data to pass into JAGS
survsection_data <- list(survtable=survtable, 
                         sectiontable=sectiontable,
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
  survsection_jags_out <- jagsUI::jags(model.file=survsection_jags, data=survsection_data,
                                parameters.to.save=c("p","survtable"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

# checking output & assessing convergence (all looks good)
nbyname(survsection_jags_out)
par(mfrow=c(2,2))
plotRhats(survsection_jags_out)
traceworstRhat(survsection_jags_out)

caterpillar(survsection_jags_out, p="p", row=1)
caterpillar(survsection_jags_out, p="p", row=2)
caterpillar(survsection_jags_out, p="p", row=3)

comparecat(list(as.data.frame(survsection_jags_out$sims.list$p[,1,]),
                as.data.frame(survsection_jags_out$sims.list$p[,2,]),
                as.data.frame(survsection_jags_out$sims.list$p[,3,])))

# plotting survival & survival probability
par(mfrow=c(1,1))
par(mar=c(6.1, 4.1, 4.1, 2.1))
plot(NA, xlim=c(1,ncol(survtable)), ylim=0:1,
     xlab="", ylab="Survival Probability", xaxt="n")
for(i in 1:nrow(survtable)) {
  jmeans <- survsection_jags_out$mean$survtable[i,] + runif(ncol(survsection_jags_out$mean$survtable), -0.01, 0.01)
  jdates <- 1:ncol(survsection_jags_out$mean$survtable) + runif(ncol(survsection_jags_out$mean$survtable), -0.05, 0.05)
  lines(jdates, jmeans, col=adjustcolor(1, alpha.f = .4))
  points(jdates, jmeans, col=adjustcolor(1, alpha.f = .4))
}
axis(side=1, at=1:length(plotdates), plotdates, las=2)
par(mar=parmar)


survsection_jags_out$DIC  # 1263.714
surv_jags_out$DIC         # 1623.812
