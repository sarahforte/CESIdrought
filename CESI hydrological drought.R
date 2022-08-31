#######################################################################
#
#    -###- CANADIAN ENVIRONMENTAL SUSTAINABILITY INDICATORS -###-
#              -- Water Quantity Indicator Calculator --
#
#  This project contains scripts used to automate the calculation
#          of CESI drought/dry year indicator.
#
#  Environment and Climate Change Canada
#  Created: August 2022
#
#######################################################################

# These scripts calculate the presence/absence of hydrological drought
# during the summer at hydrometric station based on a 30 year reference
# period.
# Different methodologies are used for perennial and intermittent
# watercourses, but both involve defining a station-specific threshold
# and comparing the yearly summer data to the threshold.
# Kriged maps and trends over 50 years are calculated with the results.
# The specific steps are:
#    1) Define variables
#    2) Create tables to compile results
#    3) Retrieve hydrometric data from hydat database
#    4) Check if there is sufficient data for calculations
#    5) Identify summer period
#    6A) For perennial watercourses:
#       i. calculate rolling average
#      ii. define threshold by fitting curve to 30-year reference data
#     iii. evaluate each year against threshold
#    6B) For intermittent watercourses:
#       i. calculate dry spell durations
#      ii. define threshold from cumulative distribution of dry spell duration
#     iii. evaluate each year against threshold
#    7) Identify 50-year trends using logistic regression
#    8) Save output as csv files
#    9) Krige results into pan-Canadian maps
#   10) Determine 50-year trends with logistical regression

# Assuming there are three sub-folders in working directory: "Dependencies"
# (with files for station list [RHBN_U.csv] and country boundaries
# [CanadaBound.shp]), "Output" (in which to save the output), and "R" (in which
# 2 script files are located [CESI hydrological drought.R-THIS ONE] and
# [CESI hydrological drought-functions.R])
#######################################################################

# load necessary libraries
library(tidyhydat) # for interacting with hydat database
library(lfstat)    # for fitting curves to low flow statistics
library(sp) # to create spatial objects
library(rgdal) # read vector maps into spatial objects
library(automap) # to model variogram
library(gstat)  # for Kriging
library(raster) # processing spatial raster data. !!!overwrites dplyr::select!!!
library(dplyr)     # for data tidying and data tables

source('./R/CESI hydrological drought-functions.R')

######### STEP 1) Define variables #########
stations <- read.csv("./Dependencies/RHBN_U.csv", header = TRUE) #list of hydrometric stations for calculations
yrs.of.int <- c(1970:2019) #period of interest
yrs.of.ref <- c(1981:2010) #30-year reference period
n.day <- 7  #number of days for rolling average
return.freq <- 5 #return frequency in years for n.day average
perc <- 0.9 #percentile for threshold dry spell duration
summer.beg <- "05-01" #first day of summer period
summer.end <- "10-01" #last day of summer period
map.years <- c(2000:2019) #range of years for output maps


######### STEP 2) Create tables #########
threshold.per <- data.frame(STATION_NUMBER=NA, censored=NA, freq0=NA,
                    best.dist=NA, R2=NA,param1=NA, param2=NA, param3=NA,
                    ret.freq.flow=NA)
threshold.int <- data.frame(STATION_NUMBER=NA, duration=NA)

drought.per <- data.frame(STATION_NUMBER=NA, Year=NA, win_min=NA,
                  d.less.ret.flow=NA, flag=NA)
drought.int <- data.frame(STATION_NUMBER=NA, Year=NA, dr_events=NA, flag=NA)
drought <- data.frame(matrix(ncol=2+length(yrs.of.int), nrow=nrow(stations)))
colnames(drought) <- c("STATION_NUMBER", "type", yrs.of.int)
drought$STATION_NUMBER <- stations$STATION_NUMBER

log.reg <- data.frame(matrix(ncol=6, nrow=nrow(stations)))
colnames(log.reg) <- c("STATION_NUMBER", "slope", "intercept", "R2", "p", "significant")
log.reg$STATION_NUMBER <- stations$STATION_NUMBER


######### STEP 3) Retrieve hydrometric data #########
for (j in 1:length(stations$STATION_NUMBER)) { #start loop that cycles through each station
  print(paste0(j, " - ", stations$STATION_NUMBER[j]))
  flow.daily <- hy_daily_flows(stations[j,1])
  flow.daily$Year <- as.numeric(format(flow.daily$Date, "%Y")) #re-arrange year to faciliate sorting
  flow.daily <- filter(flow.daily,flow.daily$Year %in% yrs.of.int) #only keep data from years of interest
  valid_year <- unique(flow.daily$Year) #note which years have data
  len.v.year <- length(valid_year)
  flow.daily$M_D <- substr(flow.daily$Date, 6, 10) #also reformat months and days for use later
  flow.daily$Month <- as.numeric(format(flow.daily$Date, "%m"))


######### STEP 4) Check if there is sufficient data #########
  # requirements for calculating thresholds are at least 150 days of data in at
  # least 20 years during the reference period
  ref <- flow.daily %>% filter(Year %in% yrs.of.ref)
  summ <- ref %>% group_by(Year) %>% summarise(count=length(Value[!is.na(Value)]))
  if (sum(summ$count>=150)>=20) { #use if loop to stop calculations in case of insufficient data
    ref_years <-  summ %>% filter(summ$count>=150) %>% pull(Year) #note which of the reference years are available for that station

######### STEP 5) Identify summer period #########
    # 3 criteria are used to identify summer:
    #     1) from May 1 to October 1
    #     2) ice-free period during this window
    #     3) period after freshet during this window
    # Start by creating tables to compile summer flow and rolling averages
    days <- format(seq.Date(as.Date("05-01",format="%m-%d"), as.Date("10-01",
              format="%m-%d"), by = 1), format="%m-%d") #create vector with all possible days of summer (155)
    flow.summer <- data.frame(matrix(ncol=155,nrow=len.v.year))
    colnames(flow.summer) <- c("Year", days)
    flow.summer[,1] <- valid_year

    max_w_per_year <- as.numeric(as.Date(summer.end, format="%m-%d") -
                        as.Date(summer.beg, format="%m-%d")) - n.day + 2
    windows <- format(seq.Date(as.Date(summer.beg, format="%m-%d")+(n.day-1),
                 as.Date(summer.end, format="%m-%d"), by = 1), format = "%m-%d")
    averages <- data.frame(matrix(ncol=max_w_per_year+2,nrow=len.v.year))
    colnames(averages) <- c("Year", windows, "ann.min")
    averages[,1] <- valid_year

    # Calculate median freshet date from 30-year median hydrograph as the date
    # with max flow between 1 march and 31 july
    ## This step occurs before loop because data will be used in loop
    median.flow <- data.frame(Date=format(seq.Date(as.Date("03-01",format="%m-%d"),
                     as.Date("07-31",format="%m-%d"), by=1), format="%m-%d"),
                     flow=NA)
    for (b in 1:nrow(median.flow)) {
      median.flow$flow[b] <- median(flow.daily$Value[flow.daily$M_D==median.flow$Date[b]],
                                na.rm = TRUE)
    }
    med.freshet <- median.flow %>% filter(flow==max(median.flow$flow, na.rm=TRUE)) %>%
                     pull(Date)

    # Find summer period for all years
    for (k in 1:len.v.year) { #cycle through all the years for that station
      summer <- merge(days, flow.daily[(flow.daily$Year==valid_year[k] &
                  flow.daily$M_D>="05-01" & flow.daily$M_D<="10-01"),c(7,4,5,8)],
                  by.x="x", by.y="M_D", all.x=TRUE) #isolate data from all summer days (May 1-Oct 1)
      summer <- ice.free(summer) #only keep flow values for ice-free period
      summer <- after.freshet(summer, med.freshet) #only keep flow values after freshet
      flow.summer[k,2:155] <- summer$Value


######### STEP 6A) For perennial watercourses: #########
#----- i. calculate rolling average for summer windows  -----
      for (l in 1:max_w_per_year) {
        if (all(!is.na(flow.summer[k,l:(l+n.day-1)]))) { # only calculate averages when there are values for all the days in the window
          averages[k,l+1] <- mean(as.numeric(flow.summer[k,l:(l+n.day-1)]))
        }
      }
      #annual minimum window
      if (sum(!is.na(averages[k,-1]))>0) { # only calculate for those years with at least one average value during summer
        averages$ann.min[k] <- min(averages[k,-1],na.rm=TRUE)
      }
    } #this closes the loop cycling through all the years for that station

#----- ii. define threshold by fitting curve to 30-year reference data -----
    curve.data <- ret.freq.curve(return.freq, averages[,c(1,ncol(averages))], ref_years)
    threshold.per[nrow(threshold.per)+1,] <- c(stations$STATION_NUMBER[j], curve.data)

#----- iii. evaluate which years have average minimum flows below the threshold -----
    if (!is.na(curve.data[8])) { #only evaluate if there is a threshold
      # identify stations with thresholds as perennial
      drought$type[j] <- "perennial"
      # number of days below
      d.less.ret.flow <- rowSums(averages[,c(-1,-ncol(averages))]<as.numeric(curve.data[8]),
                           na.rm=TRUE)
      # flag for at least one period below return flow
      less.ret.flow <- ifelse(d.less.ret.flow>0, 1, 0)

      # Add station data to both detailed and general compilation tables
      drought.per[(nrow(drought.per)+1):(nrow(drought.per)+len.v.year),] <-
        data.frame(rep(stations$STATION_NUMBER[j],len.v.year), valid_year,
                   averages$ann.min, d.less.ret.flow, less.ret.flow)
      drought[j,which(names(drought) %in% valid_year)] <- less.ret.flow
    }


######### STEP 6B) For intermittent watercourses: #########
    # Identify hydrometric stations on intermittent watercourses as those that:
    #     1) have more than 25% annual n-day minima as 0 during 30 reference
    #          years (censored=TRUE)
    #     2) have a return flow frequency of less than 0.001 m3/s (detection limit)
    #     3) do not have a calculated return frequency flow.
    if (curve.data[1]==TRUE & (is.na(curve.data[8]) | curve.data[8]<0.001)) {
      # flag station as intermittent
      drought$type[j] <- "intermittent"

#----- i. calculate dry spell durations  -----
      no.flow.spells <- data.frame(matrix(ncol=79,nrow=len.v.year)) #number of columns is max possible number of dry spells (total summer days/2) + 1
      no.flow.spells[,1] <- valid_year

      for (k in 1:len.v.year) { #cycle through all the years for that station
        counter <- ifelse(!is.na(flow.summer$'05-01'[k]) &
                      flow.summer$'05-01'[k]==0, 1, 0)
        column <- 2
        if (sum(!is.na(flow.summer[k,]))>1) { #only calculate for years with data
          for (t in 2:ncol(flow.summer)) { #cycle through each of the columns with flow data
            if (!is.na(flow.summer[k,t])){ #for those columns with a value
              if (flow.summer[k,t]==0){ #if the flow value is zero, increase the counter by one
                counter <- counter+1
              } else { #still for cases when column has flow value
                if (flow.summer[k,t]>0 & counter!=0) { #if value is above zero and counter is not zero, that denotes the end of a dry spell
                  no.flow.spells[k,column] <- counter #note the length of the dry spell
                  column <- column+1  #increase the column number not to overwrite for the next noting
                  counter <- 0 #set the counter back to zero
                }
              }
              # case when summer ends in drought or NA
              if ((flow.summer[j,ncol(flow.summer)]==0 |
                   is.na(flow.summer[j,ncol(flow.summer)])) & counter>0) {
                no.flow.spells[k,column] <- counter
              }
            }
          }
        }
      }

#----- ii. define threshold from cumulative distribution of dry spell duration  -----
      dry.spells <- no.flow.spells %>% filter(no.flow.spells$X1 %in% ref_years &
                     !is.na(no.flow.spells$X2)) %>% dplyr::select(2:79) %>%
                     apply(1,function(x) na.omit(x)) %>%  unlist()
      # There needs to be a minimum of 11 dry spells during the 30 year reference
      # period to calculate a threshold
      if (length(dry.spells)>10) {
        Th <- quantile(dry.spells, probs=perc, na.rm=TRUE)
        threshold.int[nrow(threshold.int)+1,] <- c(stations$STATION_NUMBER[j], Th)

#----- iii. evaluate each year against threshold  -----
        events <- rowSums(no.flow.spells[,-1]>Th, na.rm=TRUE)
        # flag for at least one period below return flow
        dry <- ifelse(events>0, 1, 0)
        # Add station data to both detailed and general compilation tables
        drought.int[(nrow(drought.int)+1):(nrow(drought.int)+len.v.year),] <-
          data.frame(rep(stations$STATION_NUMBER[j],len.v.year), valid_year, events, dry)
        drought[j,which(names(drought) %in% valid_year)] <- dry
      } #close loop for intermittent stations with threshold
    } # close loop for all intermittent stations
    print(drought$type[j])

######### STEP 7) Identify 50-year trends #########
    # Since results are a binary presence/absence, use logistic regression to test
    # for trends - whether the presence of droughts is increasing or decreasing
    # check if there is sufficient data to calculate a trend
    # minimum requirements: some data 1970-1975, >=30 points & no gap over 10 years
    data <- drought[j,(which(!is.na(drought[j,3:52]))+2)] #use the presence/absence data for station
    goodyears <- as.numeric(names(data)) #find which years have data
    gap.check <- na.omit(goodyears - lag(goodyears)) #check for holes in data continuity
    if (all(any(goodyears %in% c(1970:1975)), (length(goodyears) >= 30),
            (max(gap.check) <= 11))){
      # calculate regression as a generalized linear model of binomial family
      logistic <- glm(as.numeric(data) ~ goodyears,data,family=binomial)
      # Extract model results
      coeff <- logistic$coefficients[c(2,1)]
      ll.null <- logistic$null.deviance/-2
      ll.proposed <- logistic$deviance/-2
      R2 <- (ll.null-ll.proposed)/ll.null
      p <- 1-pchisq(2*(ll.proposed-ll.null),df=(length(logistic$coefficients)-1))
      significant <- case_when(p<=0.05 ~ TRUE,
                               p>0.05 ~ FALSE,
                               is.na(p) ~ NA)
      # Add station data to both detailed and general compilation tables
      log.reg[j,2:6] <- c(coeff, R2, p, significant)
    }
  } #close loop for calculations on those stations which have sufficient data
} #close loop cycling through all the stations


######### STEP 8) Save output as csv files #########
# delete the first row of the threshold tables and two detailed drought tables
# as they will be NA
threshold.per <- threshold.per[-1,]
threshold.int <- threshold.int[-1,]
drought.per <- drought.per[-1,]
drought.int <- drought.int[-1,]
# harmonize TRUE/FALSE entries
threshold.per$censored[threshold.per$censored==0] <- FALSE
threshold.per$censored[threshold.per$censored==1] <- TRUE

# then write all the files, omitting the row numbers
write.csv(threshold.per, "./Output/CESIdrought-perennial thresholds.csv", row.names = FALSE)
write.csv(threshold.int, "./Output/CESIdrought-intermittent thresholds.csv", row.names = FALSE)
write.csv(drought.per, "./Output/CESIdrought-perennial drought.csv", row.names = FALSE)
write.csv(drought.int, "./Output/CESIdrought-intermittent drought.csv", row.names = FALSE)
write.csv(drought, "./Output/CESIdrought-drought presence-absence.csv", row.names = FALSE)
write.csv(log.reg, "./Output/CESIdrought-logistic.regression.csv", row.names = FALSE)


######### STEP 9) Krige results into pan-Canadian maps #########
pan.can.krige(drought, map.years, "./Dependencies", "./Output/CESI drought 2000-2019.pdf")


