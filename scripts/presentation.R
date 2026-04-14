
library(TSA)
library(fUnitRoots)
library(forecast)
library(CombMSC)
library(lmtest)
library(tseries)

# CLEARING WORKSPACE

# clearing all the previous stuff
rm(list = ls())
#--------------------------------------------------------------------------------

# LOADING LIBRARIES

# loading the TSA package (i need this for my time series functions)

library(TSA)
library(tseries)
library(fUnitRoots)
library(lmtest)
library(forecast)
library(fpp2)  # a10 dataset is from this package
#--------------------------------------------------------------------------------

# HELPER FUNCTIONS FOR MODEL COMPARISON AND RESIDUAL ANALYSIS
# These functions will help evaluate AIC/BIC and analyze residuals of ARIMA/GARCH models
#--------------------------------------------------------------------------------

# function to sort models based on aic or bic values
sort.score <- function(x, score = c("bic", "aic")) {
  if (score == "aic") {
    x[with(x, order(AIC)), ]
  } else if (score == "bic") {
    x[with(x, order(BIC)), ]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
}

# function to analyze residuals from ARIMA, GARCH, or hybrid models
residual.analysis <- function(model, std = TRUE, start = 2, class = c("ARIMA", "GARCH", "ARMA-GARCH", "fGARCH")[1]) {
  library(TSA)
  library(FitAR)
  
  if (class == "ARIMA") {
    if (std == TRUE) {
      res.model = rstandard(model)
    } else {
      res.model = residuals(model)
    }
  } else if (class == "GARCH") {
    res.model = model$residuals[start:model$n.used]
  } else if (class == "ARMA-GARCH") {
    res.model = model@fit$residuals
  } else if (class == "fGARCH") {
    res.model = model@residuals
  } else {
    stop("The argument 'class' must be either 'ARIMA' or 'GARCH'")
  }
  
  par(mfrow = c(3, 2))
  plot(res.model, type = 'o', ylab = 'Standardised residuals', main = "Time series plot of standardised residuals")
  abline(h = 0)
  hist(res.model, main = "Histogram of standardised residuals")
  qqnorm(res.model, main = "QQ plot of standardised residuals")
  qqline(res.model, col = 2)
  seasonal_acf(res.model, main = "ACF of standardised residuals")
  print(shapiro.test(res.model))
  LBQPlot(res.model, lag.max = 30, StartLag = 1, k = 0, SquaredQ = FALSE)
  par(mfrow = c(1, 1))
}
#--------------------------------------------------------------------------------

##################################################
# Following functions are developed by           #
# MATH1318 students                              #
# Le Van Tra Tran and Tin Trung Pham             #
# in 2024. We thank them for their contribution! #
##################################################

helper <- function(class = c("acf", "pacf"), ...) {
  params <- match.call(expand.dots = TRUE)
  params <- as.list(params)[-1]
  
  if (class == "acf") {
    acf_values <- do.call(acf, c(params, list(plot = FALSE)))
  } else if (class == "pacf") {
    acf_values <- do.call(pacf, c(params, list(plot = FALSE)))
  }
  
  acf_data <- data.frame(
    Lag = as.numeric(acf_values$lag),
    ACF = as.numeric(acf_values$acf)
  )
  
  seasonal_lags <- acf_data$Lag %% 1 == 0
  
  if (class == "acf") {
    do.call(acf, c(params, list(plot = TRUE)))
  } else if (class == "pacf") {
    do.call(pacf, c(params, list(plot = TRUE)))
  }
  
  for (i in which(seasonal_lags)) {
    segments(x0 = acf_data$Lag[i], y0 = 0, x1 = acf_data$Lag[i], y1 = acf_data$ACF[i], col = "red")
  }
}

# seasonal_acf ------------------------------------------------------------
seasonal_acf <- function(...) {
  helper(class = "acf", ...)
}

# seasonal_pacf -----------------------------------------------------------
seasonal_pacf <- function(...) {
  helper(class = "pacf", ...)
}
#--------------------------------------------------------------------------------

# SETTING WORKING DIRECTORY AND LOADING DATA

# setting the working directory to the folder where my data file is stored
setwd("/Users/kamleshhhkale/Documents/RMIT/SEM 3/TIME SERIES ANALYSIS/ASS-FINAL")

# checking if i have set the correct working directory
getwd()

# loading the a10 dataset (monthly anti-diabetic drug sales in AUD)
data(a10)

# checking the first few observations to see the format
print(head(a10))

# checking the summary statistics and structure of the raw data
print(summary(a10))
str(a10)

# checking for any missing values in the data
na_indices <- which(is.na(a10))
cat("missing values are found at indices:", na_indices, "\n")
#--------------------------------------------------------------------------------

# 'a10' is already a ts object with start = Jan 1992 and end = Jul 2008, frequency = 12
a10_ts <- a10
a10_ts
plot(a10_ts)

# Trend: ✅ seeing a clear upward trend in anti-diabetic drug sales over time
# Seasonality: ✅ strong seasonality, with visible peaks around mid-year (likely winter effect)
# Changing variance: ✅ slight increase in variance over time, especially after 2002
# Behaviour: ⚠️ since there is strong seasonality, it's hard to clearly tell AR or MA behaviour right now
# Change Point: ❌ no obvious or sudden change point is visible — looks fairly continuous throughout
#--------------------------------------------------------------------------------

# INITIAL PLOT AND SUMMARY STATS

# checking descriptive statistics of the time series data
print(summary(a10_ts))
str(a10_ts)

# plotting a time series plot to see the overall pattern of anti-diabetic drug sales
plot(a10_ts, type = "o", main = "Fig 1. Time Series Plot of Monthly Anti-Diabetic Drug Sales", 
     ylab = "Sales (in Million AUD)", cex.main = 1.8)


# adding seasonal markers to the time series plot
points(time(a10_ts), a10_ts, pch = as.vector(season(a10_ts)), col = "blue")

# plotting acf and pacf side by side to check autocorrelation patterns
par(mfrow = c(1, 2))
seasonal_acf(a10_ts, lag.max = 64, main = "Fig 2. ACF plot of a10 series")
seasonal_pacf(a10_ts, lag.max = 64, main = "Fig 3. PACF plot of a10 series")
par(mfrow = c(1, 1))

# checking stationarity using the augmented dickey-fuller (ADF) test
adf.test(a10_ts)

# plotting a qq plot to check normality of a10 series
qqnorm(y = a10_ts, main = "Fig 4. QQ plot of a10 series", cex.main = 1.8)
qqline(y = a10_ts, col = 2, lwd = 1, lty = 2)

# performing shapiro-wilk test to check normality
shapiro.test(a10_ts)
#--------------------------------------------------------------------------------

# LOG TRANSFORMATION TO STABILIZE VARIANCE

# applying log transformation to the a10_ts series before differencing
model_log <- log(a10_ts)

# plotting the log-transformed series to see stabilized variance
plot(model_log, type = "o", ylab = "Log(Sales in AUD)", 
     main = "Fig 5. Time Series Plot of Log-Transformed a10", cex.main = 1.8)
#--------------------------------------------------------------------------------

# MODEL SPECIFICATION USING RESIDUAL APPROACH

# we're following the residual approach (not classical) to specify the SARIMA model
# first, we fit a plain model with only seasonal differencing (D = 1)
# this helps us see if the seasonal trend can be removed

# fitting arima model with no AR or MA terms, but with seasonal differencing on log-transformed data
model_log_000_010 <- Arima(model_log, order = c(0, 0, 0), seasonal = list(order = c(0, 1, 0), period = 12))

# extracting residuals
res_log_000_010 <- residuals(model_log_000_010)

# plotting time series of residuals
par(mfrow = c(1, 1))
plot(res_log_000_010, xlab = "Time", ylab = "Residuals", 
     main = "Fig 6. Time Series Plot of Residuals (Log + D = 1)")

# plotting seasonal acf and pacf of the residuals
par(mfrow = c(1, 2))
seasonal_acf(res_log_000_010, lag.max = 48, main = "Fig 7. ACF of Residuals (Log + D = 1)")
seasonal_pacf(res_log_000_010, lag.max = 48, main = "Fig 8. PACF of Residuals (Log + D = 1)")
par(mfrow = c(1, 1))
#--------------------------------------------------------------------------------

# ADDING SEASONAL ARMA COMPONENT (P = 3, Q = 2)

# based on the ACF and PACF plots, we now add the seasonal component (P=3, D=1, Q=2)
# non-seasonal p, d, q are still 0 for now — we’ll tune them next

model_log_000_312 <- Arima(model_log, order = c(0, 0, 0), 
                           seasonal = list(order = c(3, 1, 2), period = 12))

# extracting residuals
res_log_000_312 <- residuals(model_log_000_312)

# plotting residuals and their acf/pacf
par(mfrow = c(1, 1))
plot(res_log_000_312, xlab = "Time", ylab = "Residuals", 
     main = "Fig 9. Residuals from SARIMA(0,0,0)(3,1,2)[12] on Log a10")

par(mfrow = c(1, 2))
seasonal_acf(res_log_000_312, lag.max = 48, main = "Fig 10. ACF of Residuals (SARIMA(0,0,0)(3,1,2))")
seasonal_pacf(res_log_000_312, lag.max = 48, main = "Fig 11. PACF of Residuals (SARIMA(0,0,0)(3,1,2))")
par(mfrow = c(1, 1)) 
#--------------------------------------------------------------------------------
 
# ADF TEST ON RESIDUALS FROM model_log_000_312

adf.test(res_log_000_312)
#--------------------------------------------------------------------------------

# APPLYING FIRST DIFFERENCING TO REMOVE REMAINING NON-STATIONARITY

# Based on ADF test results (p ≈ 0.066), we now apply first-order differencing (d = 1)
# while retaining the seasonal SARIMA component (P=3, D=1, Q=2)

model_log_010_312 <- Arima(model_log, 
                           order = c(0, 1, 0), 
                           seasonal = list(order = c(3, 1, 2), period = 12))

# extracting residuals
res_log_010_312 <- residuals(model_log_010_312)

# plotting residuals and their acf/pacf
par(mfrow = c(1, 1))
plot(res_log_010_312, xlab = "Time", ylab = "Residuals", 
     main = "Fig 12. Residuals from SARIMA(0,1,0)(3,1,2)[12] on Log a10")

par(mfrow = c(1, 2))
seasonal_acf(res_log_010_312, lag.max = 48, main = "Fig 13. ACF of Residuals (SARIMA(0,1,0)(3,1,2))")
seasonal_pacf(res_log_010_312, lag.max = 48, main = "Fig 14. PACF of Residuals (SARIMA(0,1,0)(3,1,2))")
par(mfrow = c(1, 1))

# ADF TEST ON FINAL RESIDUALS
adf.test(res_log_010_312)
#--------------------------------------------------------------------------------

# FITTING FINAL SARIMA MODEL (3,1,2)(3,1,2)[12] AFTER EXAMINING ACF/PACF

# Based on residual patterns and confirmed stationarity,
# we now fit the full SARIMA model including both seasonal and non-seasonal terms

model_log_312_312 <- Arima(model_log, 
                           order = c(3, 1, 2), 
                           seasonal = list(order = c(3, 1, 2), period = 12))

# extracting residuals
res_log_312_312 <- residuals(model_log_312_312)

# plotting time series of residuals
par(mfrow = c(1, 1))
plot(res_log_312_312, xlab = "Time", ylab = "Residuals", 
     main = "Fig 15. Residuals from SARIMA(3,1,2)(3,1,2)[12]")

# plotting acf and pacf to check if remaining autocorrelations are removed
par(mfrow = c(1, 2))
seasonal_acf(res_log_312_312, lag.max = 48, main = "Fig 16. ACF of Residuals (SARIMA(3,1,2)(3,1,2))")
seasonal_pacf(res_log_312_312, lag.max = 48, main = "Fig 17. PACF of Residuals (SARIMA(3,1,2)(3,1,2))")
par(mfrow = c(1, 1))

# running adf test to reconfirm residual stationarity
adf.test(res_log_312_312)
#--------------------------------------------------------------------------------

# EXTENDED AUTOCORRELATION FUNCTION (EACF)

# we're applying EACF on the residuals of model_log_010_312
# this model had d = 1 and seasonal (3,1,2), but p and q were still 0
# since there was leftover autocorrelation in the residuals (as seen in ACF/PACF),
# we now use EACF to explore good combinations of (p,q) for the non-seasonal part

eacf(res_log_010_312)

# The tentative models are specified as 
# SARIMA(3,1,4)x(3,1,2)_12
# SARIMA(3,1,5)x(3,1,2)_12
# SARIMA(4,1,4)x(3,1,2)_12  
# SARIMA(4,1,5)x(3,1,2)_12
# SARIMA(3,1,2)x(3,1,2)_12

#--------------------------------------------------------------------------------
# BIC TABLE FOR MODEL SELECTION (USING armasubsets)

# using armasubsets to visually explore BIC scores for different (p,q) combinations
# we apply this on the residuals of model_log_010_312, which had leftover signal
# this helps narrow down promising ARMA(p,q) structures for the non-seasonal part

par(mfrow = c(1, 1))
bic_table <- armasubsets(y = res_log_010_312, nar = 10, nma = 10, y.name = "p", ar.method = "ols")

# plotting the BIC table
plot(bic_table)

# The tentative models are specified as 
# SARIMA(4,1,6)x(3,1,2)_12
# SARIMA(4,1,7)x(3,1,2)_12
# SARIMA(1,1,8)x(3,1,2)_12
#--------------------------------------------------------------------------------

# Fitting the 8 SARIMA models using ML
m3_314 <- Arima(model_log, order=c(3,1,4), seasonal=list(order=c(3,1,2), period=12), method="ML")
m3_315 <- Arima(model_log, order=c(3,1,5), seasonal=list(order=c(3,1,2), period=12), method="ML")
m3_414 <- Arima(model_log, order=c(4,1,4), seasonal=list(order=c(3,1,2), period=12), method="ML")
m3_415 <- Arima(model_log, order=c(4,1,5), seasonal=list(order=c(3,1,2), period=12), method="ML")
m3_312 <- Arima(model_log, order=c(3,1,2), seasonal=list(order=c(3,1,2), period=12), method="ML")
m3_416 <- Arima(model_log, order=c(4,1,6), seasonal=list(order=c(3,1,2), period=12), method="ML")
m3_417 <- Arima(model_log, order=c(4,1,7), seasonal=list(order=c(3,1,2), period=12), method="ML")
m3_118 <- Arima(model_log, order=c(1,1,8), seasonal=list(order=c(3,1,2), period=12), method="ML")

# AIC and BIC comparison
sc.AIC <- AIC(m3_314, m3_315, m3_414, m3_415, m3_312, m3_416, m3_417, m3_118)
sc.BIC <- BIC(m3_314, m3_315, m3_414, m3_415, m3_312, m3_416, m3_417, m3_118)

sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "bic")

# Accuracy measures
Sm3_314 <- accuracy(m3_314)[1:7]
Sm3_315 <- accuracy(m3_315)[1:7]
Sm3_414 <- accuracy(m3_414)[1:7]
Sm3_415 <- accuracy(m3_415)[1:7]
Sm3_312 <- accuracy(m3_312)[1:7]
Sm3_416 <- accuracy(m3_416)[1:7]
Sm3_417 <- accuracy(m3_417)[1:7]
Sm3_118 <- accuracy(m3_118)[1:7]

df.Smodels <- data.frame(
  rbind(Sm3_314, Sm3_315, Sm3_414, Sm3_415, 
        Sm3_312, Sm3_416, Sm3_417, Sm3_118)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE", "ACF1")
rownames(df.Smodels) <- c("SARIMA(3,1,4)x(3,1,2)_12", "SARIMA(3,1,5)x(3,1,2)_12",
                          "SARIMA(4,1,4)x(3,1,2)_12", "SARIMA(4,1,5)x(3,1,2)_12",
                          "SARIMA(3,1,2)x(3,1,2)_12", "SARIMA(4,1,6)x(3,1,2)_12",
                          "SARIMA(4,1,7)x(3,1,2)_12", "SARIMA(1,1,8)x(3,1,2)_12")
round(df.Smodels, 3)
