# Martin Holdrege

# Functions for use in final project (part 3) for Geospatial Analysis (Spring '19)

# stack_NA ---------------------------------------------------------------

stack_NA <- function(stack) {
  # args:
  #   stack --raster stack
  # returns:
  #   stack w/ negatives returned as NA
  stack2 <- calc(stack, fun = function(x){
    ifelse(x < 0, NA, x)
  })
  stack2
}

# mean daily for a month, to total monthly ----------------------------------

mean2total_prcp <- function(x, month) {
  # args:
  #   x--vector of monthly precip values, measured as units/day
  #   month (numeric vector, of length x)
  # returns:
  #   monthly precip (units/month) (does not account for leap years)
  stopifnot(length(x) == length(month))
  days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  daysx <- days[month] # ndays in each given month
  out <- x*daysx
  out
}

# month_check -------------------------------------------------------------

month_check <- function(x) {
  # args:
  #   x--vector (1 per month of year)
  # returns:
  #   Null if all true (check to be used in other functions)
  stopifnot(length(x) == 12,
            is.vector(x),
            is.numeric(x))
}


# precip concentration index ----------------------------------------------

PCI <- function(x, na.rm = FALSE) {
  # args:
  #   x, where x is length 12 vector of monthly sums of precip
  #   na.rm--remove NAs?
  # returns:
  #   Precip concentration index, see Fatichi et al 2012
  month_check(x)
  
  p_yr <- sum(x, na.rm = na.rm) # annual precip
  p_m2 <- x^2 # squared monthly precip
  pci <- (100*sum(p_m2, na.rm = na.rm))/p_yr^2
  pci
}


# annual CV -----------------------------------------------------------------

CV_annual <- function(x, na.rm = FALSE) {
  # args:
  #   x, where x is length 12 vector of monthly sums of precip
  #   na.rm--remove NAs?
  # returns:
  #   coefficient of variation of monthly precip
  month_check(x)
  
  xbar <- mean(x, na.rm = na.rm)
  sd <- sd(x, na.rm = na.rm)
  cv <- sd/xbar
  cv
}

# slope --------------------------------------------------------------------

slope <- function(y, x, na.rm) {
  # args:
  #   y--response vector
  #   x--predictor vector
  #   na.rm --not instituted, catches na.rm from stackApply
  # returns:
  #   slope, change in y with change in x 
  if(all(is.na(y))){
    return(NA)
  } # lm throws error if all NAs 
  
  stopifnot(is.numeric(y),
            is.numeric(x))
  
  beta <- lm(y ~ x)$coefficients[2] 
  beta
}


# slope from each pixel in stack ------------------------------------------

stack_slope <- function(x, time, indices = NULL) {
  # args:
  #  x raster stack 
  #  numerical vector same length as layers in x
  #  indices as passed to stackApply, by default lm applied to whole stack
  # returns:
  #   raster slope y ~ time, for each pixel
  
  if (is.null(indices)){
    indices = rep(1, x@file@nbands)
  }
  
  stackApply(x, indices = indices, fun = function(y, na.rm){
    slope(y = y, x = time, na.rm = na.rm)
  })
}

# sum or NA ---------------------------------------------------------------

sum_na <- function(x, na.rm) {
  # args:
  #  x--numeric vector
  # returns:
  #  sum of x (or NA, not 0 if all NAs)
  if(all(is.na(x))){
    return(NA)
  }
  sum(x, na.rm = na.rm)
}

# world map ---------------------------------------------------------------

plot_precip <- function(r, title = "", midpoint = NULL, pal = "RdBu", legend_title = "",
                        world = wrld2) {
  # args:
  #   r--raster
  #   title--plot title
  #   midpoint--midpoint of diverging color scale (mean used by default)
  #   pal--color pallette
  #   world -- world map shapefile
  # returns:
  #   global map
  
  if(is.null(midpoint)){
    midpoint = mean(r@data@values, na.rm = TRUE)
  }
  
  tm_style(style = "white", 
           main.title = title,
           main.title.size = 1)  +
    tm_shape(r) +
    tm_raster(title = legend_title, style = "cont", palette = pal, midpoint = midpoint, 
              colorNA = NULL) +
    tm_shape(world) +
    tm_polygons(alpha = 0, border.apha = 1) +
    tm_legend(text.size = 0.8,
              title.size = 1,
              position = c("left", "bottom"),
              bg.color = "white",
              bg.alpha = 0.7,
              frame = "gray50",
              height = 0.6)
}


# histogram function -----------------------------------------------------

g_hist1 <- function(df, title = NULL, xlab = NULL){
  # args:
  #   df--dataframe with a column labeled "x"
  # returns:
  #   histogram with mean and median marked
  
  mn <- mean(df$x, na.rm = TRUE)
  md <- median(df$x, na.rm = TRUE)
  
  g1 <- ggplot(df) +
    geom_histogram(aes(x = x), bins = 75) + 
    theme_classic() +
    labs(title = title, x = xlab) +
    geom_vline(aes(xintercept = mn, color = "Mean")) +
    geom_vline(aes(xintercept = md, color = "Median"), linetype = "dashed") +
    scale_color_manual(name = "statistics", 
                       values = c(Median = "blue", Mean = "red")) +
    theme(legend.position = c(0.8, 0.5))
  g1
}


# Moran's I at multiple lags ----------------------------------------------

moran_cgm <- function(r, k = 4, order = 8) {
  # args:
  #   r--raster
  #   k--number of nearest neighbors (for knn)
  #   order--number of lags calculated
  # returns:
  #   object of class "spcor"
  rp <- rasterToPoints(r, spatial = TRUE) # cell centroids
  
  knn <- knn2nb(knearneigh(rp, k = k)) # neighborhood
  
  cgm <- sp.correlogram(knn, rp$index_1, order = order, method = "I", 
                        zero.policy = TRUE)
  cgm
}

