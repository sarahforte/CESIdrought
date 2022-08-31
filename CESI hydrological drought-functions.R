# functions to accompany CESI hydrological drought script

### Function to calculate ice-free period using "B" flag in hydat database which
### indicates ice influence
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
ret.freq.curve <- function(ret.freq, data, years) {

# library(lfstat)    # for fitting curves to low flow statistics

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

### Function to Krige pan-canadian data
### Current colour scale is for presence/absence data but that could easily be modified
### Output is pdf file with both semi-variogram and kriged map for each year
pan.can.krige <- function(data, years, can.bound.path, output.name) {

#  library(sp) # to create spatial objects
#  library(rgdal) # read vector maps into spatial objects
#  library(automap) # to model variogram
#  library(gstat) # for Kriging
#  library(raster) # processing spatial raster data. !!!overwrites dplyr::select!!!

  # read in shapefile with canadian boundary "CanadaBound"
  Canada <- readOGR(dsn = can.bound.path, "CanadaBound")
  crs <- Canada@proj4string
  crs_wgs  <- CRS( "+init=epsg:4326")
  # create grid with 10km cells covering the whole country onto which data
  # should be interpolated
  Cangrid <- makegrid(Canada, cellsize=10000)
  Cangrid <- SpatialPoints(Cangrid, proj4string=crs)

  # prepare colour pallet with two CESI colours
  orange <- rgb(161,83,34, maxColorValue = 255)
  blue <- rgb(15,76,106, maxColorValue = 255)
  colours <- colorRampPalette(c(blue,orange))

  # get station coordinates
  stn_data <- as.data.frame(hy_stations() %>% dplyr::select(STATION_NUMBER, LATITUDE, LONGITUDE))
  stn_data <- stn_data[stn_data$STATION_NUMBER %in% data$STATION_NUMBER,]

  # create pdf file to populate with variograms and maps
  pdf(file = output.name, width = 8.5, height = 11)
  par(mfrow = c(5,4))

  # work on data one year at a time
  for (a in 1:length(years)) { #cycle trough all the years to be mapped
    print(years[a])
    dr_data <- drought[,(colnames(drought)=="STATION_NUMBER" | colnames(drought)==map.years[a])] #extract data for 1 year
    dr_data <- left_join(stn_data, dr_data, by="STATION_NUMBER") #combine data with station coordinates
    dr_data <- na.omit(dr_data)
    coordinates(dr_data) <- ~ LONGITUDE + LATITUDE #define which columns are coordinates
    proj4string(dr_data) <- crs_wgs #note projection to be used for maps
    dr_data <- spTransform(dr_data, CRSobj = crs)

    # generate variogram, model variogram and plot both
    var <- variogram(dr_data@data[,2]~1, dr_data, cutoff=500000, width=25000) #calculates sample variogram values with 500km max range and 30km lag distance
    model <- autofitVariogram(dr_data@data[,2]~1, dr_data, model=c("Exp", "Sph",
               "Gau"), cutoff=500000) # the best fitting model from exponential, spherical or gaussian is chosen

    # sort data to plot, and change units from m to km
    dist <- var$dist/1000
    vario <- var$gamma
    nugget <- model$var_model$psill[1]
    sill <- sum(model$var_model$psill[1:2])
    range <- model$var_model$range[2]/1000
    dist0 <- seq(0, max(dist), length.out=100)
    mod_var <- case_when(model$var_model$model[2]=="Exp" ~ (sill-nugget)*(1-exp(-dist0/range))+nugget, #exponential
                         model$var_model$model[2]=="Sph" ~ ifelse(dist0>range,
                           sill,(sill-nugget)*(1.5*dist0/range-0.5*(dist0/range)^3)+nugget), #spherical
                         model$var_model$model[2]=="Gau" ~ (sill-nugget)*(1-exp(-(dist0/range)^2))+nugget) # gaussian
    rangeC <- case_when(model$var_model$model[2]=="Exp" ~ 3*range, #exponential
                        model$var_model$model[2]=="Sph" ~ range, #spherical
                        model$var_model$model[2]=="Gau" ~ 2*range) # gaussian

    # plot output
    par(mar = c(4, 4, 2, 1))
    plot(dist, vario, xlab="", ylab="semivariance", xlim=c(0,500), ylim=c(0,
          (ceiling(max(var$gamma, na.rm=TRUE)*100))/100), las=1, col="blue")
    lines(dist0, mod_var, col="blue")
    title(main=paste0(years[a], " semivariogram"), cex=0.7)
    title(xlab="distance [km]", line=2)
    mtext(paste0("model: ", model$var_model$model[2], "   "), side=1, line=-6, adj=1, cex=0.7)
    mtext(paste0("nugget: ", round(nugget, 2), "   "), side=1, line=-5, adj=1, cex=0.7)
    mtext(paste0("sill: ", round(sill, 2), "   "), side=1, line=-4, adj=1, cex=0.7)
    mtext(paste0("range: ", round(rangeC, 0), "   "), side=1, line=-3, adj=1, cex=0.7)

    max_int_d <- ifelse(rangeC>500, 500000, ceiling(rangeC/100)*100000) #maximum interpolation distance, 500km if range is over 500km, or next nearest 100 if range below 500km
    # Krige data using ordinary kriging and the variogram model created
    OK <- krige(formula=dr_data@data[,2]~1, locations=dr_data, newdata=Cangrid,
              model=model$var_model, maxdist=max_int_d, nmin=8, nmax=20)
    # transform results into dataframe so they can be transformed into a raster
    OKfr <- OK %>% as.data.frame %>% dplyr::select(X = x1, Y = x2, Z = var1.pred)
    OKraster <- rasterFromXYZ(OKfr, res=c(10000, 10000), crs = Canada@proj4string)
    # only keep the part of the raster than is in canadian boundaries
    OKraster <- mask(OKraster, Canada)

    # plot output
    #par(mar = c(0, 0, 0, 0), oma = c(0.5, 0.75, 0.5, 0.75))
    par(mar = c(0, 0, 0, 0))
    plot(Canada, col="grey85", border = NA)
    plot(OKraster, add=TRUE, col=colours(2), zlim=c(0,1), legend=FALSE)
    plot(Canada, add=TRUE, col=NA, border = "white", lwd = 0.3)
    legend("topright", inset=0.13, legend = c("abscence", "presence"),
         col = NA, pch = NA, cex=0.7, fill = c(blue, orange),border = NA,
         box.lty=0)
    title(main=paste0(years[a], " Drought"), line = -2, cex=0.7)
  } #end loop cycling through years to map
  dev.off() #close pdf opened to save plots
} #close function
