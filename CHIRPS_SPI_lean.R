
# Now, need to create RasterBrick and calculate SPI
## spi_pentad from spi_computation.R

require(raster)
library(raster)

setwd("/Users/nathansuberi/Desktop/WRI_Programming/chirps_hist")
files <- list.files(pattern="*.tif$")
fnames <- list.files(pattern="*.tif$")
rnames <- list.files()
files

files[1]
ref <- raster(files[1])
ref

two_years <- files[1:144]

s <- stack(two_years)
brick <- brick(s)
ncell(brick) # should be 14400000
length(brick[10000]) # 144 observations, 1 for each time slot
brick[1100000] # random, finding non-zero entries

## How to create a timeseries object that includes the pentads?
#time_series = ts(data = NA, start=1981, frequency=72)
#head(time_series)
#as.Date(time(time_series)[150], origin="1981-01-01")

# Example from SO: https://stackoverflow.com/questions/29202021/r-how-to-extract-dates-from-a-time-series
#x = seq(1,768)
#myts <- ts(x, start=1982, frequency=24)
#time(myts)
#library(zoo)
#as.yearmon(time(myts))

#### Handler

getSPI <- function(cell, brick, startYear, startMonth, scale) {
  #format to time series
  timeseries <- ts(t(as.matrix((brick[cell]))), start=c(startYear, startMonth), frequency=72)
  return(spi_pentad(timeseries, scale, na.rm=TRUE))
}

getSPI(1100000, brick, 1981, 1, 6) # scale = 6, means this is a monthly SPI

#### Main SPI-pentad method
install.packages("parallel")
install.packages("lmomco")
library(lmomco)
library(parallel)

spi_pentad <- function(data, scale, kernel=list(type='rectangular',shift=0),
                       distribution='Gamma', fit='ub-pwm', na.rm=FALSE, 
                       ref.start=NULL, ref.end=NULL, x=FALSE, params=NULL, ...) {
  
  # Frequency, number of periods before looping
  # For pentads, this is 72. 6 pentads per month * 12 months = 72.
  scale <- as.numeric(scale)
  # print(paste("scale is: ", scale))
  
  # Decide whether to allow times series with any missing data
  na.rm <- as.logical(na.rm)
  # print(paste("Allow missing data: ", na.rm))
  
  # NOT SURE WHAT THIS IS
  x <- as.logical(x)
  
  if (sum(is.na(data))>0 & na.rm==FALSE) {
    stop('Error: Data must not contain NAs')
  }
  if (distribution!='Gamma') {
    stop('Gamma distributions are used for Standard Precipitation Indices"')
  }
  
  ### Convert the data to a time series if it is not already
  # if (!is.ts(data)) {
  #  data <- ts(as.matrix(data), frequency = 72)
  #} else {
  #  data <- ts(as.matrix(data), frequency=frequency(data), start=start(data))
  #}
  
  # m = 1, only precipitation data
  m <- ncol(data)
  # 72 pentads per year
  fr <- frequency(data) 
  # print(paste("dims of data: ", m, ", num periods: ", fr))
  
  # create coefficient matrix for Gamma distribution coefficients
  coef <- array(NA,c(2,m,fr),list(par=c('alpha','beta'),"precipitation",NULL))
  
  std <- data*NA
  
  ## Apply kernel to aggregate data
  
  # Cumulative series (acu)
  acu <- data[,1]
  
  # Create a kernal @ scale of time series
  # returns a vector of same length as time series
  wgt <- kern(scale,"rectangular",0)
  
  # https://stat.ethz.ch/R-manual/R-devel/library/stats/html/embed.html
  # "embeds the time series in a low-dimensional Euclidean Space." Thanks.
  acu[scale:length(acu)] <- rowSums(embed(acu, scale)*wgt, na.rm=FALSE)
  acu[1:{scale-1}] <- NA
  
  # Loop through the pentads, 1-72
  for (pentad in (1:fr)) {
    
    # Filter pentad, excluding NAs
    ps <- which( cycle(acu) == pentad )
    
    # 1 is considered NA here because of the steps above...
    # STILL NOT SURE WHY THOSE VALS ARE ERASED
    ps <- ps[ !is.na(acu[ps]) ]
    
    # Pentadly series, sorted
    pentad_ser <- sort(acu[ps])
    
    # If no data for that pentad, no std dev
    if (length(pentad_ser)==0) {
      std[ps] <- NA
      next()
    }
    
    # if the data has NA sd, no std dev
    # looking at a list of precipitation observations,
    # in a given brick cell, 
    # in a given pentad
    if (is.na(sd(pentad_ser,na.rm=TRUE)) | (sd(pentad_ser, na.rm=TRUE) == 0)) {
      std[ps] <- NA
      next()
    }
    
    # Probability of pentadly precipitation = 0 (pze)
    zeros <- sum(pentad_ser==0)
    pze <- sum(pentad_ser==0)/length(pentad_ser)
    
    
    # Adding nmom=2 as an option made all the conditions below clear
    pwm <- pwm.ub( pentad_ser[pentad_ser>0] , nmom=2)
    lmom <- pwm2lmom( pwm )
    if (!are.lmom.valid(lmom) | is.na(sum(lmom[[1]])) | is.nan(sum(lmom[[1]]))) {
      next()
    }					
    
    gampar <- pargam(lmom)
    # Compute standardized values
    # What is diff between qnorm and pnorm?
    std[ps,1] <- qnorm( cdfgam( acu[ps], gampar ) )
    std[ps,1] <- qnorm( pze + (1-pze) * pnorm(std[ps,1] ) )
    # fill in coefficient matrix from above
    coef[,1,pentad] <- gampar$para
    
  } # next c (month)
  
  colnames(std) <- "SPI"
  
  # Create z, which returns the stds
  z <- list(call=match.call(expand.dots=FALSE),
            fitted=std,coefficients=coef,scale=scale,kernel=list(type="rectangular",
                                                              shift=0),values=kern(scale,"rectangular",0),
            distribution="Gamma",fit="ub-pwm",na.action=FALSE)
  
  
  # Here it is again... why return the original data?
  if (x) z$data <- data
  
  # Not considering reference windows
  # if (!is.null(ref.start)) z$ref.period <- rbind(ref.start,ref.end)
  
  class(z) <- 'spi_pentad'
  return(z)
}

## Helper function, smooth over data to 
# incorporate past events in the SPI measurement

kern <- function(scale, type='rectangular', shift=0) {
  if(type!='rectangular' & type!='triangular' & type!='circular' & type!='gaussian') {
    stop('type must be one of: rectangular, triangular, circular, gaussian')
  }
  #
  s <- scale
  
  # Don't respond to shift
  #h <- shift
  #if(h>=s) {
  #  stop('Parameter "shift" must be lower than "scale"')
  #}
  #if(h<0) {
  #  stop('Parameter "shift" cannot have a negative value')
  #}
  
  if(s<=0) {
    stop('Parameter "scale" must be higher than zero')
  }
  
  # Kernal functions for incorporating past data-points in 
  # time-series evaluation
  if(type=='rectangular' | s==1) k <- rep(1,s)
  if(type=='triangular') k <- s:1
  if(type=='circular') k <- (s^2+(1-(1:s)^2))
  if(type=='gaussian') k <- (1/0.4)*1/sqrt(2*pi*1^2)*
    exp(-(seq(0,-3,-3/(s-1))-0)^2/2*1^2)
  
  # Not clear how this works, and it locks my computer
  # if(h) k <- c(k[(h+1):2],k[1:(s-h)])
  
  #return(k/sum(k))
  return(k/mean(k))
}

###### Get to work!

# Get answers, but don't think they're right
one_call <- getSPI(1100000, brick, 1981, 1, 72)

# loop through all cells in the brick, calculate SPI values for that cell
# in all raster layers in the Brick. 
# The beauty of the Brick.
l <- lapply( seq(ncell(brick)), function(cell){
  
  if( !is.na(brick[cell][1]) ) {
    return(c(getSPI(cell, brick, 1981, 1, 72)$fitted))
  } else {
    #filler for NA cells in the first of the timeseries.
    #assume: if NA in the first cell then NA in the rest (need to check this)
    
    # Reminder that each value corresponds to a file
    return(rep(NA, length(files)))
  }
})

if(require(rgdal)){
  brick_file <- writeRaster(brick, filename="raster_brick.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
}

# This is no longer a RasterBrick.
# Next - see how fast RasterBrick computation can be done on AWS
brick <- raster("raster_brick.tif")


#### CONTINUE FROM HERE once l finishes

coordref <- CRS(brick@crs@projargs)
orig <- origin(brick)

setwd("/Users/nathansuberi/Desktop/WRI_Programming/chirps_hist/spi")

# write them spi rasters to 'spi' folder
if(dir=="current") {
  if("spi" %in% rnames) stop("spi folder already exists, either remove or move it")
  dir.create("spi")
} else {
  if("spi" %in% rnames) stop("spi folder already exists, either remove or move it")
  dir.create(paste0(dir, "/spi"))
}

# Remember this, from above???
# ref <- raster(files[1])

# why offset by timescale here? BECAUSE you need at least # of TIMESCALE observations 
# in order for the kernel to work
# check that coordref lays 2000 x 7200 cells2
lapply(seq(from=72+1, to=length(files)), function(layer) {
  data <- sapply(l, function(x) {
    x[layer]
  })
  r <- ref
  values(r) <- data
  crs(r) <- coordref
  origin(r) <- orig
  fileName <- paste0(dir,"/spi/", "spi_", fnames[layer])
  writeRaster(r, fileName, format="GTiff")
})

