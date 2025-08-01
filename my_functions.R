my_func_get_single_ppt_data <- function(file) {
  
  # get data
  raw <- read_csv(paste0("~/OneDrive - University of Toronto/Projects/Yas accent/data_analysis/gam/", file))
  # select channel
  df <- raw %>%
    # select channel
    select(c(time, all_of(chan), participant, item)) %>%
    rowwise() %>%
    mutate(uV = mean(c_across(all_of(chan)), na.rm = TRUE)) %>%
    ungroup() %>%
    select(c(time, participant, item, uV)) %>%
    # split into conditions
    separate(col = item, into = c("condition", "item"), sep="/") %>%
    mutate(participant = as.factor(participant),
           item = as.factor(item),
           condition = as.factor(condition),
           time = as.numeric(time)) %>%
    mutate(condition = as.factor(condition)) %>%
    droplevels()
  
  return(df)
  
}

my_func_extract_gam_measures <- function(dat, search_min, search_max) {
  
  # subset search data
  sdat <- dat[dat$time>=search_min & dat$time<=search_max, ]
  
  # get derivative and search for peak (derivative=0) of a negativity (previous derivative value < 0, which means the actual ERP waveform is decreasing before this point)
  drv <- data.frame(diff(dat$fit)/diff(dat$time))  # derivative
  colnames(drv) <- 'dYdX'
  drv$time <- rowMeans(embed(dat$time,2)) # center the X values for plotting
  drv$dYdX.next <- c(drv$dYdX[2:nrow(drv)],NA)
  drv$dYdX.prev <- c(NA,drv$dYdX[1:(nrow(drv)-1)])
  
  # MMN peak: going down (<0) then going up (>0)
  drv$local_peak <- ((drv$dYdX < 0 & drv$dYdX.next > 0) | (drv$dYdX.next > 0 & drv$dYdX == 0 & drv$dYdX.prev < 0))
  
  # if at least one local peak in the search time window
  if (sum(drv[drv$time>=search_min & drv$time<=search_max, ]$local_peak, na.rm=TRUE) >= 1) {
    hasPeak = TRUE
    # get all peak times
    all_peak_times <- drv[which(drv$local_peak & drv$time>=search_min & drv$time<=search_max), "time"]
    # initialize peak height with some larger value
    peak_height <- Inf
    # loop over local peak times
    for (peak_ind in 1:length(all_peak_times)) {
      # get the two fitted data points centering the local peak
      peakdat = dat[dat$time >= floor(all_peak_times[peak_ind]) & dat$time <= ceiling(all_peak_times[peak_ind]), ]
      # if the current height is smaller than the original peak height
      if ( min(peakdat$fit) < peak_height) {
        # update peak height
        peak_height <- min(peakdat$fit)
        # update peak time
        peak_time <- all_peak_times[peak_ind]
        # update se
        peak_se <- peakdat[which.min(peakdat$fit),]$se.fit
        # update NMP (original code doesn't have 1.96 factor)
        NMP <- peak_height / (1.96*peak_se) # relative peak measure (if < 1 then 95%CI overlaps with 0 at point of peak)
      }
    }
  } else { # if no local peak
    hasPeak <- FALSE
    # get general peak in search span
    subdat <- dat[dat$time>=search_min & dat$time<=search_max, ] # subset data
    # find peak
    peak_height <- min(subdat$fit)
    peak_index <- which.min(subdat$fit)
    # get time
    peak_time <- subdat[peak_index, "time"] # first time value with peak value
    peak_se <- subdat[peak_index, ]$se.fit
    NMP <- peak_height / (1.96*peak_se)
  }
  
  # if we are looking for a valley but the value is positive, then no correct positivity/negativity
  if (peak_height >= 0) {
    hasPeak <- FALSE
    peak_height <- NA
    peak_time <- NA
    peak_se <- NA
    NMP <- NA
  } else if ((peak_time == min(sdat$time))) { # if peak time is the first point, there is so no real minimum
    hasPeak <- FALSE
    peak_height <- NA
    peak_time <- NA
    peak_se <- NA
    NMP <- NA
  }
  
  # get area and fractional are latency
  if (is.na(peak_time)) {
    area <- NA
    half_area_latency <- NA
    firsttime <- NA
    lasttime <- NA
  } else {
    # initialize are
    area <- 0
    start = round(peak_time) # start time to take integral from
    # firsttime = search_min
    # lasttime = search_max
    # area to the right from the peak
    for (i in start:search_max) {
      val = sdat[sdat$time == i,]$fit
      
      # if derivative <=0
      if (val <= 0) {
        area = area + abs(val)
      } else { # end of peak, so stop going in this direction
        break
      }
      lasttime = i
    }
    # area to the left from the peak
    beforestart = start-1
    if (beforestart >= search_min) {
      for (j in beforestart:search_min) { # to the left from the peak
        val = sdat[sdat$time == j,]$fit
        # if derivative <=0
        if (val <= 0) {
          area = area + abs(val)
        } else { # end of peak, so stop going in this direction
          break
        }
        firsttime = j
      }
    }
    
    # get half area latency
    halfarea <- 0
    for (k in firsttime:lasttime) {
      val = sdat[sdat$time == k, ]$fit
      halfarea = halfarea + abs(val)
      if (halfarea >= 0.5 * area) {
        half_area_latency = k
        break
      }
    }
  } # get area and fractional are latency end
  
  return(list(
    sdat = sdat,
    drv = drv,
    hasPeak = hasPeak, 
    area = area, 
    peak_height = peak_height, 
    peak_se = peak_se, 
    NMP = NMP, 
    peak_time = peak_time,
    firsttime = firsttime,
    lasttime = lasttime,
    half_area_latency = half_area_latency))
  
}

# function compute p-value matrix
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}




# x is a matrix containing the data
# method : correlation method. "pearson"" or "spearman"" is supported
# removeTriangle : remove upper or lower triangle
# results :  if "html" or "latex"
# the results will be displayed in html or latex format
corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                    result=c("none", "html", "latex")){
  #Compute correlation matrix
  require(Hmisc)
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value 
  
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
  
  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  if (result[1]=="none") return(Rnew)
  else{
    if(result[1]=="html") print(xtable(Rnew), type="html")
    else print(xtable(Rnew), type="latex") 
  }
} 

