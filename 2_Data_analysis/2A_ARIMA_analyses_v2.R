#2A_ARIMA_analyses_v2
#Author: Mary Lofton
#Date: 27NOV20

# load packages
pacman::p_load(tidyverse, lubridate, forecast, tseries, utils)

# get data
mydata <- read_csv("./2_Data_analysis/FP_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_Depth_m)
mydata <- mydata[,c(2,1,3:21)]

# convert to timeseries
max_biomass_2016 <- ts(data = mydata$Max_biomass_ugL[1:23], frequency = 52, start = c(2016, 17))
max_biomass_2017 <- ts(data = mydata$Max_biomass_ugL[25:46], frequency = 52, start = c(2017, 18))
max_biomass_2018 <- ts(data = mydata$Max_biomass_ugL[48:67], frequency = 52, start = c(2018, 19))
max_biomass_2019 <- ts(data = mydata$Max_biomass_ugL[75:104], frequency = 52, start = c(2019, 16))

peak_depth_2016 <- ts(data = mydata$Peak_depth_m[1:23], frequency = 52, start = c(2016, 17))
peak_depth_2017 <- ts(data = mydata$Peak_depth_m[25:46], frequency = 52, start = c(2017, 18))
peak_depth_2018 <- ts(data = mydata$Peak_depth_m[48:67], frequency = 52, start = c(2018, 19))
peak_depth_2019 <- ts(data = mydata$Peak_depth_m[75:104], frequency = 52, start = c(2019, 16))

peak_magnitude_2016 <- ts(data = mydata$Peak_magnitude_ugL[1:23], frequency = 52, start = c(2016, 17))
peak_magnitude_2017 <- ts(data = mydata$Peak_magnitude_ugL[25:46], frequency = 52, start = c(2017, 18))
peak_magnitude_2018 <- ts(data = mydata$Peak_magnitude_ugL[48:67], frequency = 52, start = c(2018, 19))
peak_magnitude_2019 <- ts(data = mydata$Peak_magnitude_ugL[75:104], frequency = 52, start = c(2019, 16))

peak_width_2016 <- ts(data = mydata$Peak_width_m[1:23], frequency = 52, start = c(2016, 17))
peak_width_2017 <- ts(data = mydata$Peak_width_m[25:46], frequency = 52, start = c(2017, 18))
peak_width_2018 <- ts(data = mydata$Peak_width_m[48:67], frequency = 52, start = c(2018, 19))
peak_width_2019 <- ts(data = mydata$Peak_width_m[75:104], frequency = 52, start = c(2019, 16))


# A. Model form selection
# 
# 1. Evaluate stationarity (of full timeseries and each year alone)

# The basic stationarity diagnostics are the following:
# 
#   1. Plot your data. 
ts.plot(max_biomass_2016, ylab = "Max. biomass (ug/L)")
ts.plot(max_biomass_2017, ylab = "Max. biomass (ug/L)")
ts.plot(max_biomass_2018, ylab = "Max. biomass (ug/L)")
ts.plot(max_biomass_2019, ylab = "Max. biomass (ug/L)")

ts.plot(peak_depth_2016, ylab = "Peak depth (m)")
ts.plot(peak_depth_2017, ylab = "Peak depth (m)")
ts.plot(peak_depth_2018, ylab = "Peak depth (m)")
ts.plot(peak_depth_2019, ylab = "Peak depth (m)")

ts.plot(peak_magnitude_2016, ylab = "Peak magnitude (ug/L)")
ts.plot(peak_magnitude_2017, ylab = "Peak magnitude (ug/L)")
ts.plot(peak_magnitude_2018, ylab = "Peak magnitude (ug/L)")
ts.plot(peak_magnitude_2019, ylab = "Peak magnitude (ug/L)")

ts.plot(peak_width_2016, ylab = "Peak width (m)")
ts.plot(peak_width_2017, ylab = "Peak width (m)")
ts.plot(peak_width_2018, ylab = "Peak width (m)")
ts.plot(peak_width_2019, ylab = "Peak width (m)")

# LOOK FOR: 
#     A. An increasing trend
# - max biomass and peak magnitude increase over the season
# - peak depth may decrease over the season?
# - peak width doesn't look like it increases or decreases
#     B. A non-zero level (if no trend)
# - peak width looks like it might have a level ~3
#     C. Strange shocks or steps in your data (indicating something dramatic changed like the data collection methodology)
# - something funky about peak depth at the end of 2017; it's very static

#   2. Apply stationarity tests
stationarity.tests <- matrix(NA, nrow = 16, ncol = 4)
all.ts <- list(max_biomass_2016,max_biomass_2017,max_biomass_2018,max_biomass_2019,
               peak_depth_2016,peak_depth_2017,peak_depth_2018,peak_depth_2019,
               peak_magnitude_2016,peak_magnitude_2017,peak_magnitude_2018,peak_magnitude_2019,
               peak_width_2016,peak_width_2017,peak_width_2018,peak_width_2019)

for (i in 1:16){
  adf <- adf.test(all.ts[[i]])
  kpss <- kpss.test(all.ts[[i]])
  adf_diffs <- ndiffs(all.ts[[i]], test = "adf")
  kpss_diffs <- ndiffs(all.ts[[i]], test = "kpss")
  stationarity.tests[i,1] <- adf$p.value
  stationarity.tests[i,2] <- adf_diffs
  stationarity.tests[i,3] <- kpss$p.value
  stationarity.tests[i,4] <- kpss_diffs
}

stationarity.table <- data.frame(adf = stationarity.tests[,1], adf_diffs = stationarity.tests[,2],
                                 kpss = stationarity.tests[,3],kpss_diffs = stationarity.tests[,4])
rownames(stationarity.table) <- c("Max_biomass_ugL_2016",
                                  "Max_biomass_ugL_2017",
                                  "Max_biomass_ugL_2018",
                                  "Max_biomass_ugL_2019",
                                  "Peak_depth_m_2016",
                                  "Peak_depth_m_2017",
                                  "Peak_depth_m_2018",
                                  "Peak_depth_m_2019",
                                  "Peak_magnitude_ugL_2016",
                                  "Peak_magnitude_ugL_2017",
                                  "Peak_magnitude_ugL_2018",
                                  "Peak_magnitude_ugL_2019",
                                  "Peak_width_m_2016",
                                  "Peak_width_m_2017",
                                  "Peak_width_m_2018",
                                  "Peak_width_m_2019")


#     A. adf.test() p-value should be less than 0.05 (reject null)
#     B. kpss.test() p-value should be greater than 0.05 (do not reject null)
#   3. If stationarity tests are failed, then try differencing to correct
#     Try ndiffs() in the forecast package or manually try different differences.

# 2. Selection of the differencing level (d) – to fix stationarity problems
#     See point 3 in stationarity diagnostics above.
# - overall, max_biomass and peak_magnitude should probably all be differenced by 1
# - peak_depth only needs a diff in 2016 and it could be either 1 or 2 depending on test
# - peak_width doesn't need to be differenced

# 3. Selection of the AR level (p)
for (i in 1:4){
  fit <- auto.arima(all.ts[[i]])
  print(fit)
  checkresiduals(fit)
  fc <- forecast(fit, h = 4)
  plot(fc)
}
# 4. Selection of the MA level (q)
# B. Parameter estimation
# 
# C. Model checking


##NOTES:

# use the following to get you started on an algorithm of all possible combinations
# https://stat.ethz.ch/R-manual/R-devel/library/utils/html/combn.html

# do this for auto.arima() to see which is best. save the aicc

# that's it.