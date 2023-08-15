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


# proportion



## run timing??



## estimate proportions within each river section, per survey



## identify probable spawning areas in the mainstem Tanana River during late January

# empirical k-functions 
# review this methodology!  it is very slow and doesn't seem to gain anything
telemdata1 <- telemdata %>% filter(flight_num>=1 & current_state!="dead")
par(mfrow=c(4,4))
with(telemdata1, kfunc(seg=seg, vert=vert, survey=mnDate,
      rivers=tyb_trim, envreps=100, maxdist=200000))

# kernel density (network)
# darn, this doesn't seem to gain anything either
network_dens <- with(telemdata1, makeriverdensity(seg=seg, vert=vert, survey=mnDate,
                                       rivers=tyb_trim, bw=20*1000))
# plotriverdensitypoints(network_dens)  # maybe more points than needed
par(mfrow=c(4,4))
plot(network_dens, scalebyN=F)
par(mfrow=c(1,1))
plot(network_dens, points=F, scalebyN = F)  # maybe I should try animating this

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

## I SHOULD BE CONSIDERING MORTALITY IN THIS
with(telemdata, table(current_state, flight_num))



## describe seasonal distributions and migrations (this will be the big one)



## survival analysis thing? might be cool

# define survival table
# - each row is an individual
# - each column is a flight event (not counting tagging)
# - 0=alive, 1=dead, NA=not found
survtable <- telem3d$current_state_num[,-1]            # take out tagging event
survtable <- survtable[rowMeans(is.na(survtable))<1,]  # take out indivs that were never found

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
      }
      if(survtable[i,j]==0) {
        if(is.na(firstalive[i])) firstalive[i] <- j
        survtable[i, firstalive[i]:j] <- 0 
      }
    }
    j <- j+1
    if(j>ncol(survtable)) founddead <- T
  }
}
firstdead[is.na(firstdead)] <- ncol(survtable)  # necessary for Bayesian model to function


# Bayesian model
# - Of interest is parameter vector p[], which represents the probability of
#   mortality for each sequential pair of surveys.
# - If an individual is not observed (has value NA) in a given survey, the unobserved
#   state will be treated as an unknown parameter

library(jagsUI)
library(jagshelper)
# skeleton("surv")
# specify model
surv_jags <- tempfile()
# cat('model {
#   for(i in 1:n) {                       # for each individual
#     for(j in 1:firstdead[i]) {          # for each survey
#       survtable[i,j] ~ dbern(p[j])   ### this does not account for prev time step!!
#     }
#   }
#   for(j in 1:np) {
#     p[j] ~ dbeta(.1,.1)
#   }
# }', file=surv_jags)
cat('model {
  for(i in 1:n) {                       # for each individual
    survtable[i,1] ~ dbern(p[1])  # first observation
  }
  for(i in firstmult:n) {  
    for(j in 2:firstdead[i]) {          # for each survey
      survtable[i,j] ~ dbin(p[j], survtable[i,j-1])  ## what if i make 0=dead, 1=alive?   
    }
  }
  for(j in 1:np) {
    p[j] ~ dbeta(5,5)
  }
}', file=surv_jags)

# bundle data to pass into JAGS
# surv_data <- list(survtable=survtable,
#                   firstdead=firstdead,
#                   n=nrow(survtable),
#                   np=ncol(survtable))
surv_data <- list(survtable=1-survtable[order(firstdead),],  ## what if i make 0=dead, 1=alive? 
                  firstdead=firstdead[order(firstdead)],
                  n=nrow(survtable),
                  np=ncol(survtable))
surv_data$firstmult <- which.max(surv_data$firstdead>1)

# JAGS controls
niter <- 10000
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
# nbyname(surv_jags_out)
# par(mfrow=c(2,2))
# plotRhats(surv_jags_out)
# traceworstRhat(surv_jags_out)


# plot 50 & 95% credible intervals for each p!
par(mfrow=c(1,1))
par(mar=c(6.1, 4.1, 4.1, 2.1))
caterpillar(surv_jags_out, p="p", xax=rep("", surv_data$np))  # 50 & 95% credible intervals for each p!
axis(side=1, at=1:surv_data$np, mnDates[-1], las=2)
par(mar=parmar)

# ppost <- jags_df(surv_jags_out, p="p")
pmed <- surv_jags_out$q50$p
pse <- surv_jags_out$sd$p
plo <- surv_jags_out$q2.5$p
phi <- surv_jags_out$q97.5$p
ndigits <- 3
ci95 <- paste0("(",round(plo, ndigits)," - ", round(phi, ndigits),")")
ptable <- data.frame(p_est=round(pmed, ndigits),
                     SE=round(pse, ndigits),
                     CI95=ci95)
rownames(ptable) <- mnDates[-1]
ptable   # print to console
knitr::kable(ptable)  # print to markdown


# I'm so curious how this compares to raw estimate (without modeling unknown states)
raw_n <- sum(!is.na(survtable[,1]))
raw_p <- mean(survtable[,1], na.rm=T)
for(j in 2:ncol(survtable)) {
  raw_n[j] <- sum(survtable[,j-1]==0 & !is.na(survtable[,j]), na.rm=T)
  raw_p[j] <- sum(survtable[,j-1]==0 & survtable[,j]==1, na.rm=T)/raw_n[j]
}
points(raw_p)  # always lower... did i do it right??
# table(survtable[,2], survtable[,3], useNA="always")   ### starting to check....
# 
# for(j in 2:ncol(survtable)) {
#   temptable <- table(survtable[,j-1], survtable[,j], useNA="always")
#   raw_n[j] <- sum(temptable[1,1:2])
#   raw_p[j] <- temptable[1,2]/raw_n[j]
# }

# plo <- raw_p - 2*sqrt(raw_p*(1-raw_p)/(raw_n-1))
# phi <- raw_p + 2*sqrt(raw_p*(1-raw_p)/(raw_n-1))
# points(plo, pch="-")
# points(phi, pch="-")

# table: 
# - confirmed alive at j-1
# - confirmed newly dead at j
# - unknown at j-1 with %% avg mort
# - unknown at j with %% avg mort
