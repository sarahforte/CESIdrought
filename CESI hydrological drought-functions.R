# functions to accompany CESI hydrological drought script

### Function to calculate ice-free period using "B" flag in hydat database which
### indicates ice influence
##### Input:  data = dataframe with at minimum 3 columns: x (MM-DD), Month 
#####                       Symbol (site condition code) 

ice.free <- function(data) {
  # First B day
  firstbday <- data %>% filter(Month > 8) %>%  filter(Symbol == "B")  %>% head(1)
  firstbday.date <- ifelse(nrow(firstbday)>0, firstbday$x, NA)
  if (!is.na(firstbday.date)) {
    FU <- which(as.Date(data$x,format="%m-%d")==as.Date(firstbday.date,format="%m-%d"))
    data$Value[FU:154] <- NA
  }

  # Last B day
  lastbday <- data %>% filter(Month < 6) %>% filter(Symbol == "B") %>% tail(1)
  lastbday.date <- ifelse(nrow(lastbday)>0, lastbday$x, NA)
  if (!is.na(lastbday.date)) {
    BU <- which(as.Date(data$x,format="%m-%d")==as.Date(lastbday.date,format="%m-%d"))
    data$Value[1:BU] <- NA
  }
  return(data)
}
#-------------------------------------------------------------------------------

### Function to calculate date of freshet as date of max flow in 6 week window
### centered on date of maximum flow from 30 year median hydrograph
##### Input:  data = dataframe with at minimum 2 columns: x (MM-DD), Value (flow), 
#####         date.med.fresh = numeric value of median freshet date

after.freshet <- function(data, date.med.fresh) {
  window.start <- which(data$x==format((as.Date(date.med.fresh, format="%m-%d")-21),
                    format="%m-%d"))
  window.end <- which(data$x==format((as.Date(date.med.fresh, format="%m-%d")+21),
                  format="%m-%d"))
  if (length(window.start)>0 & length(window.end)>0) {
    FM <- which(data$Value==max(data$Value[window.start:window.end], na.rm=TRUE))
    if (!is.na(FM[1])) {
      data$Value[1:FM[1]] <- NA
    }
  }
  return(data)
}
#-------------------------------------------------------------------------------

### Function to fit a curve to the lower end of a non-exceedence distribution
### to find a return frequency flow
##### Input:  ret.freq = numeric value of drought return frequency
#####         data = dataframe with at minimum 2 columns: Year, 
#####                         ann.min (minimum annual 7-day average flow)
#####         years = vector of numeric values of years from which to calculate return frequency flow

ret.freq.curve <- function(ret.freq, data, years) {
  # Curves for 7 different distribution types are tested to keep the best fit
  dis <- c("pe3", # Pearson type III
           "wei", # Weibull
           "ln3", # lognormal
           "nor", # normal
           "exp", # exponential
           "gevR", # generalized extreme value
           "gum") # Gumbell

  ann.min.ref <- data$ann.min[which(data$Year %in% years)] %>% na.omit()
  if (sum(ann.min.ref!=0)>(0.75*length(ann.min.ref))) { #only calculate curve if more than 75% of annual minima values are non-zero
    fit <- suppressWarnings(evquantile(evfit(ann.min.ref, distribution = dis, zeta=0,
             extreme = "minimum"), ret.freq))
    # note for output if the frequency of zero values
    output <- c(fit$is.censored, fit$freq.zeros, rep(NA, 6))

    #determine best fitting function based on largest R^2 value
    ## estimates are only made on censored values, so create dataframe with
    ## censored values and estimates
    cens <- data.frame(values=fit$values[fit$values!=0],fit$estimates)
    cens <- cens[order(cens$values),]
    ## calculate half of the closest even number of censored values since we are
    ## only interested in the fit of the smaller values (lowest half)
    half <- ceiling(nrow(cens)/2)
    ## calculate R2
    R2 <- cor(cens$values[1:half],cens[1:half,2:8])^2
    #### This will be empty when the lowest half values are identical, in which
    #### a curve should not be calculated
    if (any(!is.na(R2))) { #only for those data for which there is at least one curve
      best.dist <- which(R2==max(R2)) #find which distribution curve has the highest R2
      # if best distribution gives return frequency <0.001, when there are no 0
      # values in the minimum data, try next best option
      if (fit$is.censored==FALSE & fit$T_Years_Event[best.dist]<0.001) {
        best.dist2 <- which(R2==sort(R2,partial=6)[6])
        if (fit$T_Years_Event[best.dist2]>0.001) {
          best.dist <- best.dist2
        }
      }
      # set to 0 negative return frequency flows, an artifact of calculations and not possible
      ret.flow <- ifelse(fit$T_Years_Event[best.dist]<0,  0,
                    fit$T_Years_Event[best.dist])
      # note for output details of curve
      output[3:8] <- c(dis[best.dist], max(R2),
                unlist(fit$parameters[best.dist])[1],
                unlist(fit$parameters[best.dist])[2],
                unlist(fit$parameters[best.dist])[3], ret.flow)
    }

  } else { #for stations with too many zeros
    output <- c(TRUE, ">25%", rep(NA, 6))
  }
  return(output)
}
#-------------------------------------------------------------------------------

### Function to calculate drought events
##### Input:  data = dataframe with 4 columns: day, month, year and flow
#####         dr_thresh = numeric value of drought threshold
#####         stn = alphanumeric station name

drought_days <- function(data, dr_thresh, stn) {
  # transform the data into object to use lfstat package
  data.lf <- createlfobj(data, hyearstart=1)
  flowunit(data.lf) <- "m^3/s"
  # lfstat function
  dr <- find_droughts(data.lf, threshold=dr_thresh)
  if (max(dr$event.no)>0) { # if there are drought events
    ## pool adjacent events with less than 5 days between and 
    ## Vabove / min(Vbelowi, Vbelowj)>=0.1 where, Vabove is the volume
    ## exceeding the threshold, Vbelowi and Vbelowj are the volumes 
    ## below the threshold during events i and j.
    if (max(dr$event.no, na.rm=TRUE)>1) { dr <- pool_ic(dr, tmin=5, ratio=0.1) }
    
    summary2 <- possibly(base::summary, otherwise = data.frame(event.no=NULL))
    dur<-summary2(dr) #summary function that eliminates droughts less than 5 days
    
    
    ## write messages about pooling/elimination
    if (!is.null(attributes(dur)$deficit$n.pooled)) {
      if (attributes(dur)$deficit$n.pooled>0) {
        message(paste0(stn," has ", attributes(dur)$deficit$n.pooled,
                " pooled event(s) in ",valid_year[k]))
      }
    }
    if (attributes(dur)$deficit$omitted>0) {
      message(paste0(stn," has ", attributes(dur)$deficit$omitted,
              " eliminated event(s) in ",valid_year[k]))
    }
    
    ## if all droughts were eliminated, note that there were none
    if (nrow(dur)==0) { 
      output <- data.frame(dr_events=0, dr_max_dur=0, dr_days=0)
    } else { # note results
      output <- data.frame(dr_events=nrow(dur), dr_max_dur=max(dur$duration), 
                           dr_days=sum(dur$duration))
    }
  } else { ## for those years without drought
    output <- data.frame(dr_events=0, dr_max_dur=0, dr_days=0)
  }
  return(output)
}
#-------------------------------------------------------------------------------

### Function to trends using 3 tests, when possible a hurdle test for those series
### with 3 or more zeros, a negative binomial test otherwise, and finally a
### stationarity check using all zeros or Wald-Wolfowitz test
##### Input:  data = dataframe with at minimum 2 columns: year and dr_days (drought days)
#####         dr_thresh = numeric value of drought threshold
#####         stn = alphanumeric station name

################################################################################
################################################################################
#' hurdle_test
#'
#' @description
#' hurdle_test performs a hurdle test and returns multiple test statistics given
#' a variable y and a variable x. Hurdle test can be used to fit a linear model
#' when a data set contains an excess amount of zeros.
#'
#' @param y_var a vector of numbers representing the y variable
#' @param x_var a vector of numbers representing the x variable
#' @param zero_fam a character string representing the model family for the zero
#' portion of the data. Defaulted to "binomial" indicating binomial family
#' @param count_dist a character string representing the count distribution.
#' Defaulted to 'negbin', the negative binomial distribution
#' @param zero_dist a character string representing the zero distribution.
#' Defaulted to "binomial", the binomial distribution
#' @param link a character string representing the link for binomial zero hurdle
#' model. Defaulted to 'logit', the logarithm function
#'
#' @return a dataframe of the form || slope || intercept || ConfidenceLevel || Test type,
#' or a character string if the test is not applicable;
#'
#' @details
#' The hurdle test deals with zero and non-zero portions separatly using two
#' different models: the negative binomial distribution for the non-zero portion
#' and the binomial distribution for the zero portion. Confidence intervals are
#' created to report the level of confidence in the trend detected. A generalized
#' linear model (glm) will be then established between the two distributions with
#' a logarithm link function. The slope and intercept of the resulting model will
#' be recorded and returned.
#' Things to take notice of:
#' 1. y and x variables need to have the same length
#' 2. Although the detail section only describes the condition with a negative
#'    binomial count distribution and a binomial zero distribution, this function
#'    can be used for other distributions as well.
#' 3. Do note that this function does not check for the accuracy of the fitted
#'    model, neither does it provide any recommendations on what distribution
#'    fits the best. In order to use other families, it is highly suggested to
#'    look at the concept of general linear models beforehand.
#' 4. The program will return an error if the y variable does not contain any zeros.
#'    In order to calculate the trend using this function, make sure that there
#'    is at least one zero value in the y variable. If the error message "Data
#'    cannot be fitted using hurdle model" appears, user might want to consider
#'    using negative binomial model (provided as the negbin() function in the
#'    same package) for cases with excess amount of zeros, or Mann-Kendall Test
#'    for cases with limited ties and zeros.
#'
#' @export
hurdle_test <- function(y_var, x_var, zero_fam="binomial", count_dist="negbin",
                        zero_dist="binomial", link="logit"){
  # Are we confident there is a trend?
  # Count portion of hurdle
  y_count <- y_var[y_var>0]
  x_count <- x_var[y_var>0]
  model.count <- tryCatch(MASS::glm.nb(y_count~x_count),
                          error=function(e){return("A")})
  low.c<- ifelse(is.character(model.count), NA, exp(confint(profile(model.count), level=0.9))[2,1])
  hi.c <- ifelse(is.character(model.count), NA, exp(confint(profile(model.count), level=0.9))[2,2])
  # Zero portion of hurdle
  zero_indicator <- !I(y_var==0)
  model.zero <- tryCatch(glm(zero_indicator~x_var, family=zero_fam),
                         error=function(e){return("A")})
  low.z<-ifelse(is.character(model.zero), NA, exp(confint(profile(model.zero), level=0.9))[2,1])
  hi.z<-ifelse(is.character(model.zero), NA, exp(confint(profile(model.zero), level=0.9))[2,2])
  if(!is.na(low.c)&!is.na(hi.c)&!is.na(low.z)&!is.na(hi.z)){
    # combined confidence test @ 70% and 90% confidence
    pass <- ifelse(all(low.c*low.z<=1, hi.c*hi.z>=1), "Maybe?", "Confident")
    if (pass == "Maybe?"){
      low.c<- exp(confint(profile(model.count), level=0.7))[2,1]
      hi.c <- exp(confint(profile(model.count), level=0.7))[2,2]
      low.z<-exp(confint(profile(model.zero), level=0.7))[2,1]
      hi.z<-exp(confint(profile(model.zero), level=0.7))[2,2]
      pass <- ifelse(all(low.c*low.z<=1, hi.c*hi.z>=1), "Uncertain", "Likely")
    }
    # Get slope and intercept from hurdle
    model <- hurdle(y_var~x_var, dist=count_dist, zero.dist=zero_dist, link=link)
    fitted <- unname(model$fitted.values)
    slope <- round((fitted[length(fitted)] - fitted[1])/
                     (max(x_var) - min(x_var)),2)
    intercept <- round((fitted[1]-min(x_var)*((fitted[length(fitted)] - fitted[1])/
                                                (max(x_var) - min(x_var)))),2)
    years.for.trend <- sum(!is.na(y_var))
    output <- data.frame(slope=slope, intercept=intercept, CATTrend=pass,
                         years.for.trend=years.for.trend, test="Hurdle")
  }else{ # data can't be fitted
    output <- data.frame(slope=NA, intercept=NA, CATTrend=NA, years.for.trend=NA,
                         test="Hurdle")
  }
  return(output)
}
################################################################################
################################################################################
#' negbin_test
#'
#' @description
#' negbin fits a negative binomial model and returns multiple test statistics
#' given a variable y and a variable x.
#'
#' @param y_var a vector of numbers representing the y variable
#' @param x_var a vector of numbers representing the x variable
#'
#' @return a dataframe of the form || slope || intercept || ConfidenceLevel || Test type
#'
#' @details
#' The negbin function fits a negative binomial model to a sets of x and y variables
#' Confidence intervals are created to measure the level of confidence in the trend.
#' The slope and intercept of the resulting model are recorded and returned.
#' Things to take notice:
#' 1. y and x variables need to have the same length
#' 2. Theoretically saying, negbin is the most versatile among the three trend tests.
#'    The preferred condition to use negbin function is when there is neither an
#'    excess amount of zeros, nor can the trend be calculated by Mann-Kendall test.
#'    For those conditions, consider using the two other functions instead.
#'
#' @export
negbin_test <- function(y_var, x_var){
  model <- tryCatch(glm.nb(y_var~x_var), error=function(e){return("A")})
  low<-ifelse(is.character(model), NA, exp(confint(profile(model), level=0.9))[2,1])
  hi<-ifelse(is.character(model), NA, exp(confint(profile(model), level=0.9))[2,2])
  if (!is.na(low)&!is.na(hi)){
    # confidence test @ 70% and 90% confidence
    pass <- ifelse(all(low<=1, hi>=1), "Maybe?", "Confident")
    if (pass=="Maybe?"){
      low <- exp(confint(profile(model), level=0.7))[2,1]
      hi <-exp(confint(profile(model), level=0.7))[2,2]
      pass <- ifelse(all(low<=1, hi>=1),"Uncertain", "Likely")
    }
    # Get slope from negative binomial
    fitted <- unname(model$fitted.values)
    slope <- round((fitted[length(fitted)] - fitted[1])/
                     (max(x_var) - min(x_var)),2)
    intercept <- round((fitted[1]-min(x_var)*((fitted[length(fitted)] - fitted[1])/
                        (max(x_var) - min(x_var)))),2)
    years.for.trend <- sum(!is.na(y_var))
    output <- data.frame(slope=slope, intercept=intercept, CATTrend=pass,
                         years.for.trend=years.for.trend, test="Negative Binomial")
  }else{ # data can't be fitted
    output <- data.frame(slope=NA, intercept=NA, CATTrend=NA, years.for.trend=NA,
                         test="Negative Binomial")
  }
  return(output)
}
################################################################################
# CESI project-Stationarity Test
# This part of the code provides a function that checks if a set of data is
# stationary based on if it is all zeros or the Wald-Wolfowitz test.
################################################################################
stationarity_test <- function(var){
  if (sum(var!=0)==0) { # tests don't work for only zeros, so single those out
    output <- data.frame(slope=0, intercept=mean(var, na.rm=TRUE), CATTrend="Confident",
                         years.for.trend=length(!is.na(var)), test="all zeros")
  } else {
    w_w <- ww.test(var)
    if (w_w$p.value>0.05 & !is.na(w_w$p.value)) { # only continue if null hypothesis is not rejected
      output <- data.frame(slope=0, intercept=mean(var, na.rm=TRUE), CATTrend="Confident",
                           years.for.trend=length(!is.na(var)), test="W-W")
    } else { 
      output <- data.frame(slope=NA, intercept=NA, CATTrend="Uncertain", 
                           years.for.trend=NA, test="W-W")
    }
  }
  return(output)
}
################################################################################
################################################################################

identify_trends <- function(data) {
  # test if it is possible to use hurdle test
  hurdle <- FALSE #default to False
  if (sum(data$dr_days==0) >= 3) {
    model <- tryCatch(hurdle(data$dr_days~data$year, data, dist="negbin",
                      zero.dist = "negbin"), error=function(e){return("A")},
                      warning=function(w){return("B")})
    if(!is.character(model2)){
      hurdlechk <- data.frame(hurdlechk=hurdletest(model)[2,4])
      if(!is.nan(hurdlechk[1,1])){ hurdle <- ifelse(hurdlechk>0.1, TRUE, FALSE) }
    } else {
      hurdlechk <- data.frame(hurdlechk=NA)
    }
  }
  
  if(hurdle){
    #Apply the hurdle model
    message("Hurdle test")
    output <- hurdle_test(data$dr_days,data$year)
    output <- merge(hurdlechk,output)
    
  } else {
    #Apply the negative binomial model
    message("Negative Binomial test")
    output <-negbin_test(data$dr_days,data$year)
    output <- merge(hurdlechk,output)
  }
  
  # test for stationarity on subset of stations
  if (output$CATTrend=="Uncertain" | is.na(output$CATTrend)) {
    message("testing for stationarity")
    output <- stationarity_test(data$dr_days)
    output <- merge(hurdlechk,output)
  }
  output$mapslope <- case_when(grepl("Likely", output$CATTrend)    ~ output$slope,
                               grepl("Confident", output$CATTrend) ~ output$slope)
  return(output)
}
#-------------------------------------------------------------------------------

### Function to Krige data using 5 isotropic variograms for 5 combined ecozones.
### stationarity check using all zeros or Wald-Wolfowitz test
##### Input:  data = dataframe with at minimum 2 columns: STATION_NUMBER and mapslope
#####         out = name of output file

krige_data <- function(data.t, output.name) {
  # read in shapefile with canadian boundary to set canvas boundary
  Canada <- readOGR(dsn="./Dependencies", "CanadaBound")
  crs <- Canada@proj4string
  # set projection
  crs_wgs  <- CRS( "+init=epsg:4326")
  # read in shapefiles with kriging zones which are grouped ecozones
  Kzones <- readOGR(dsn="./Dependencies", "5 merged ecozones")
  Kzones <- spTransform(Kzones, CRSobj = crs)
                      
  # get station coordinates
  stn_data <- as.data.frame(hy_stations() %>% dplyr::select(STATION_NUMBER, LATITUDE, LONGITUDE))
  stn_data <- stn_data[stn_data$STATION_NUMBER %in% data.t$STATION_NUMBER,]
  dr_data <- left_join(stn_data, data.t, by="STATION_NUMBER") #combine data with station coordinates
  coordinates(dr_data) <- ~ LONGITUDE + LATITUDE #define which columns are coordinates
  proj4string(dr_data) <- crs_wgs #note projection to be used for maps
  dr_data <- spTransform(dr_data, CRSobj = crs)

  # create grid with 10km cells covering the whole country onto which data
  # should be interpolated
  Cangrid <- makegrid(Canada, cellsize=10000)
  Cangrid <- SpatialPoints(Cangrid, proj4string=crs)
  grids <- list()

  # variogram model curve types
  Mtype <- c("Sph", "Exp")
  Mcol <- c("red", "green")

  # create pdf file to populate with variograms
  pdf(file = paste0("./Output/", output.name, ".pdf"), width = 8.5, height = 11)
  par(mfrow = c(5,4))
  # add map of zones to show which are for which variogram
  # transform zone outline and create labels
  par(mar = c(0, 0, 0, 0))
  Zoutline <- spTransform(Kzones, CRSobj = crs)
  labels <- gCentroid(Zoutline, byid=TRUE)
  labels@coords[1,2] <- 100000
  ID <- 1:length(labels)
  plot(Zoutline, col=rainbow(5)[ID], border="black", lwd = 0.4, lty=2)
  plot(Canada, add=TRUE, col=NA, border = "grey45")
  text(labels, ID, cex=1)
  title(main="Grouped ecozones", line=-1)
  

  # work on each kriging zone separately
  for (z in 1:length(Kzones)) {
    message(paste0("kriging zone ", z))
    ## create buffer around zone
    Szone <- gSimplify(Kzones[z,1], tol=5000)
    Bzone <- gBuffer(Szone, width=300000) # buffer zone is 300 km
    Mzone <- gBuffer(Szone, width=25000) # overlap zone is 50 km, 25 on each side
    ## extract stations in buffered zone
    Zstns <- dr_data[Bzone,]
  
    ## find average distance between stations in zone
    dist <- as.data.frame(as.matrix(dist(Zstns@coords)))
    dist[upper.tri(dist, diag=TRUE)] <- NA 
    MDstns <- mean(apply(dist[-1,],1,FUN=min, na.rm=TRUE))
  
    ## create variogram
    var <- variogram(Zstns@data$mapslope~1, Zstns, width=(MDstns/3), cutoff=500000) #calculates sample variogram values with 500km max range and lag distance proportional to the mean min distance between points
  
    ## two options for fitting variogram model, try automap version first, and if
    ## it creates output with nugget of 0, try gstat version.
    ### automap fitting - entirely independent including creation of new variogram,
    ### including determining which model of exponential or spherical fits best
    out <- autofitVariogram(Zstns@data$mapslope~1, Zstns, model=c("Exp", "Sph"), width= MDstns, cutoff=500000)# the best fitting model from exponential or spherical is chosen
    model <- out$var_model
    if (model$psill[1]==0) {
      ### gstat fitting - works best with start values, so calculate these firt
      ### create new variogram
      v <- variogram(Zstns@data$mapslope~1, Zstns)
      ### note data from variogram as start values 
      nugget0 <- min(v$gamma)
      psill0 <- ((max(v$gamma)+median(v$gamma))/2)-nugget0
      range0 <- 0.1*sqrt((bbox(Zstns)[1,1]-bbox(Zstns)[1,2])^2+(bbox(Zstns)[2,1]-bbox(Zstns)[2,2])^2)
      ### find best model between exponential and spherical
      vmods <- list()
      SSErr <- rep(NA, length(Mtype))
      for (k in 1:length(Mtype)) {
        vmods[[k]] <- fit.variogram(v, vgm(psill0, model=Mtype[k], range0, nugget0))
        SSErr[k] <- attr(vmods[[k]], "SSErr")
      }
      k <- which(SSErr==min(SSErr))
      ### if alternate model is better, use it
      if (vmods[[k]]$psill[1] > model$psill[1]) { 
        message("'gstat' variogram fitting has better results that 'automap' fitting")
        model <- vmods[[k]]
      }
    }
  
    ## sort variogram data to plot, and change units from m to km
    dist <- var$dist/1000
    vario <- var$gamma
    ## plot variogram first
    par(mar = c(4, 4, 2, 1))
    plot(dist, vario, xlab="", ylab="semivariance", xlim=c(0,500), ylim=c(0,
      (ceiling(max(var$gamma, na.rm=TRUE)*100))/100), las=1, col="blue")
    title(main=paste0("Variogram zone ",z), line=0.5)
    title(xlab="distance [km]", line=2)
    ## sort model data to plot
    dist0 <- seq(0, max(dist), length.out=100)
    nugget <- model$psill[1]
    sill <- sum(model$psill[1:2])
    range <- model$range[2]/1000
    k <- case_when(model$model[2]=="Exp" ~ 2, #exponential
                   model$model[2]=="Sph" ~ 1) #spherical
    MV <- case_when(k==2 ~ (sill-nugget)*(1-exp(-dist0/range))+nugget, #exponential
                    k==1 ~ ifelse(dist0>range,sill,(sill-nugget)*(1.5*dist0/range-0.5*(dist0/range)^3)+nugget)) #spherical
    rangeC <- case_when(k==2 ~ 3*range, #exponential
                        k==1 ~ range) #spherical
    ## plot model and write its parameters
    lines(dist0, MV, col=Mcol[k])
    mtext(paste0("model: ", model$model[2], "   "), side=1, line=-4, adj=1, cex=0.4, col=Mcol[k])
    mtext(paste0("nugget: ", round(nugget, 3), "   "), side=1, line=-3, adj=1, cex=0.4, col=Mcol[k])
    mtext(paste0("sill: ", round(sill, 3), "   "), side=1, line=-2, adj=1, cex=0.4, col=Mcol[k])
    mtext(paste0("range: ", round(rangeC, 0), "   "), side=1, line=-1, adj=1, cex=0.4, col=Mcol[k])
  
  
    ## krige
    OK <- krige(formula=Zstns@data$mapslope~1, locations=Zstns, newdata=Cangrid,
              model=model, maxdist=500000, nmin=8)
    ## transform results into dataframe so they can be transformed into a raster
    OKfr <- OK %>% as.data.frame %>% dplyr::select(X = x1, Y = x2, Z = var1.pred)
    OKraster <- rasterFromXYZ(OKfr, res=c(10000, 10000), crs = Canada@proj4string)
    ## only keep the part of the raster than is in the zone
    grids[[z]] <- mask(OKraster, Mzone)
  }
  
  # combine grids from all the zones, using the mean for the 50 km (25 km on each side) overlap
  grids$fun <- mean
  OKmosaic <- do.call(mosaic, grids)
  output <- mask(OKmosaic, Canada)
  
  # plot the krigged results
  # find raster min max
  zmax <-setMinMax(output)@data@max
  zmin <-setMinMax(output)@data@min
  par(mar = c(0, 0, 0, 0))
  plot(Canada, col=NA, border = NA)
  mtext(paste0("data min: ", round(min(data.t$mapslope, na.rm=NA),2) , " max: ", round(max(data.t$mapslope, na.rm=NA),2)), 1, -4, adj=0, cex=0.4)
  mtext(paste0("grid min: ", round(zmin,2) , " max: ", round(zmax,2)), 1, -3, adj=0, cex=0.4)
  
  #close pdf file
  dev.off()
  
  # end
  return(output)
}


#--------------------------------------------------
### Function to create map using CESI colours of output
##### Input:  input = gridded data to map
#####         echelle = number, limit to use for colour scale
#####         divisions = number, of divisions for colorbar labels in legend
#####         leg.title = character string with legend title
#####         out = character string of name of output file
#####         type = character string for type of output file, either "PNG" or "PDF"

map_results <- function(input, echelle, divisions, leg.title, output.name, type) {
  
  # import objects to map
  Canada <- readOGR(dsn="./Dependencies", "CanadaBound")
  crs <- Canada@proj4string
  crs_wgs  <- CRS( "+init=epsg:4326")
  Lakes <- readOGR(dsn="./Dependencies", "CANwaterbodies")
  
  # north arrow and position
  n_arrow <- readPNG("./Dependencies/esri6north.png")
  centre_n_arrow <- c(2056967,-417531)
  
  # scale bar position
  loc_scale=c(2465000, -500000)
  
  # prepare colour pallet - CESI colours
  orange <- rgb(161,83,34, maxColorValue=255)
  white <- rgb(247,247,247, maxColorValue=255)
  blue <- rgb(15,76,106, maxColorValue=255)
  blue_lakes <- rgb(150,198,250, maxColorValue=255)
  grey_background <- rgb(164,163,173, maxColorValue=255)
  grey_borders <- rgb(130,130,130, maxColorValue=255)
  colours <- colorRampPalette(c(blue,white,orange))
  
    # find raster min max
  zmax <-setMinMax(input)@data@max
  zmin <-setMinMax(input)@data@min
  
  # open file depending on type
  if (type=="PNG") {
    png(file=paste0("./Output/", output.name, ".png"), width=5, height=5,
        units="in", res=600)
  } else {
    pdf(file=paste0("./Output/", output.name, ".pdf"), width=5, height=5)
  }
  par(mar = c(0, 0, 0, 0))
  
  #plot
  plot(Canada, col=grey_background, border = NA)
  plot(surface, add=TRUE, col=colours(101), zlim=c(-echelle,echelle), axes=FALSE,
        legend.only=TRUE, axis.args=list(at=seq(-echelle,echelle,(echelle*2/divisions)),
        labels=seq(-echelle,echelle, (echelle*2/divisions)), cex.axis=0.6),
        smallplot=c(0.73,0.75, 0.57,0.82), legend.args=list(text=leg.title,
        side=3, line=0.25, cex=0.5))
  plot(surface, add=TRUE, col=orange, zlim=c(echelle,zmax), legend=FALSE, axes=FALSE)
  plot(surface, add=TRUE, col=blue, zlim=c(zmin,-echelle), legend=FALSE, axes=FALSE)
  plot(surface, add=TRUE, col=colours(1000), zlim=c(-echelle,echelle), legend=FALSE)
  plot(Canada, add=TRUE, col=NA, border=grey_borders, lwd = 0.3)
  plot(Lakes, add=TRUE, col=blue_lakes, border=grey_borders, lwd=0.1)
  legend("topright", inset = c(0.16, 0.43), cex = 0.6, bty="n",
         border= grey_borders, fill = grey_background, legend="No data")
  
  # calculations to find angle of north arrow
  dim_can_lat <- c(41.91, 83.11)
  dim_can_fig <- extent(Canada)[3:4]
  loc_n_fig <- ((dim_can_fig[2]-dim_can_fig[1])/(dim_can_lat[2]-dim_can_lat[1])*(90-dim_can_lat[2]))+dim_can_fig[2]
  angle <- atan((centre_n_arrow[1]-mean(extent(Canada)[1:2]))/(loc_n_fig-centre_n_arrow[2]))*pi
  # calculations to find corners of north arrow
  size_n_arrow <- c(150000, dim(n_arrow)[1]*150000/dim(n_arrow)[2])
  x_left <- centre_n_arrow[1]-size_n_arrow[1]/2
  x_right <- centre_n_arrow[1]+size_n_arrow[1]/2
  y_bottom <- centre_n_arrow[2]-size_n_arrow[2]/2
  y_top <- centre_n_arrow[2]+size_n_arrow[2]/2
  # plot north arrow
  rasterImage(n_arrow, x_left, y_bottom, x_right, y_top, angle=angle)
  
  # plot scale
  raster::scalebar(d=500000, xy=loc_scale, type="bar", divs=2, lonlat=FALSE,
                   label=c(0,250,500), adj=c(0.5,-0.75), cex=0.5)
  mtext("km ", side=1, line=-4.65, adj=1, cex=0.5)
  
  #close file
  dev.off()
}