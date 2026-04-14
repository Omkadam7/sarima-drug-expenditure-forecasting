# CLEARING WORKSPACE

# clearing all the previous stuff
rm(list = ls())
#--------------------------------------------------------------------------------

# LOADING LIBRARIES

library(fpp2)        # h02 dataset (PBS corticosteroids under ATC code H02)
library(tseries)     # adf.test(), shapiro.test()
library(TSA)         # LBQPlot() and time series diagnostics
library(lmtest)      # for coeftest()
#--------------------------------------------------------------------------------

# HELPER FUNCTIONS FOR MODEL COMPARISON AND RESIDUAL ANALYSIS

# These functions will help evaluate AIC/BIC and analyze residuals of ARIMA/GARCH models
 
# function to sort models based on aic or bic values
sort.score <- function(x, score = c("bic", "aic")){
  if (score == "aic"){
    x[with(x, order(AIC)),]
  } else if (score == "bic") {
    x[with(x, order(BIC)),]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
}

# function to analyze residuals from ARIMA, GARCH, or hybrid models
residual.analysis <- function(model, std = TRUE, class = "ARIMA") {
  library(TSA)
  
  if (class == "ARIMA") {
    res.model <- if (std) rstandard(model) else residuals(model)
  } else {
    stop("Only ARIMA models supported.")
  }
  
  par(mfrow = c(3, 2))  # 3 rows, 2 columns layout
  
  # Plot 1: Time series
  plot(res.model, type = 'o', ylab = 'Standardised residuals',
       main = "Time series plot of standardised residuals")
  abline(h = 0)
  
  # Plot 2: Histogram
  hist(res.model, main = "Histogram of standardised residuals")
  
  # Plot 3: QQ plot
  qqnorm(res.model, main = "QQ plot of standardised residuals")
  qqline(res.model, col = 2)
  
  # Plot 4: ACF
  seasonal_acf(res.model, main = "ACF of standardised residuals")
  
  # Plot 5: Shapiro-Wilk result (printed in console)
  print(shapiro.test(res.model))
  
  # Plot 6: Ljung-Box p-values as base R line plot
  lbq_pvals <- sapply(1:20, function(lag) {
    Box.test(res.model, lag = lag, type = "Ljung-Box")$p.value
  })
  plot(1:20, lbq_pvals, type = "b", pch = 19,
       main = "Ljung-Box P-values by Lag",
       xlab = "Lag", ylab = "P-value")
  abline(h = 0.05, col = "red", lty = 2)
  
  par(mfrow = c(1, 1))  # reset layout
}

##################################################
# Following functions are developed by           #  
# MATH1318 students                              #
# Le Van Tra Tran and Tin Trung Pham             #
# in 2024. WE thank them for their contribution! #
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

# seasonal_acf
seasonal_acf <- function(...) {helper(class = "acf", ...)}

# seasonal_pacf
seasonal_pacf <- function(...) {helper(class = "pacf", ...)}
#--------------------------------------------------------------------------------

# SETTING WORKING DIRECTORY AND LOADING DATA

# The dataset represents monthly government expenditure under Australia’s Pharmaceutical 
# Benefits Scheme (PBS) for drugs falling under ATC Code H02 (Corticosteroids). The data 
# spans from July 1991 to June 2008, measured in millions of Australian dollars (AUD).

# setting the working directory to the folder where my r script file is stored
setwd("/Users/kamleshhhkale/Documents/RMIT/SEM 3/TIME SERIES ANALYSIS/ASS-FINAL")

# checking if i have set the correct working directory
getwd()
  
# loading the h02 dataset
data(h02)

# 'h02' is already a ts object from July 1991 to June 2008, measured in millions of AUD
drug.ts <- h02 # main time series used in the analysis
#--------------------------------------------------------------------------------

# INITIAL PLOT AND SUMMARY STATS

# checking the data to see the format
drug.ts

# checking descriptive statistics of the time series data
print(summary(drug.ts))
str(drug.ts)

# checking for any missing values in the data
na_indices <- which(is.na(drug.ts))
cat("missing values are found at indices:", na_indices, "\n")

plot(drug.ts)

# Trend: ✅ seeing a clear upward trend in monthly PBS expenditure on corticosteroids (ATC H02)
# Seasonality: ✅ strong seasonal pattern, expenditure peaks mid-year (possibly winter-driven)
# Changing variance: ✅ slight increase in variance over time,
# Behaviour: ⚠️ since there is strong seasonality, it's hard to clearly tell AR or MA behaviour right now
# Change Point: ❌ no obvious or sudden change point is visible — looks fairly continuous throughout

# plotting a time series plot to see the overall pattern of monthly government drug expenditure
plot(drug.ts, type = "o", main = "Fig 1. Time Series Plot of Monthly Government Drug Expenditure", 
     ylab = "Sales (in Million AUD)", cex.main = 1.8)
# adding seasonal markers to the time series plot
points(time(drug.ts), drug.ts, pch = as.vector(season(drug.ts)), col = "blue")

y_t <- drug.ts[-1]           # current values
y_t_minus_1 <- drug.ts[-length(drug.ts)]  # lagged values

# plotting the scatterplot
plot(y_t_minus_1, y_t,
     main = "Fig 1.1. Scatter Plot of Consecutive Time Points",
     xlab = expression(Y[t-1]),
     ylab = expression(Y[t]),
     pch = 19, col = "steelblue")

# adding a trend line
abline(lm(y_t ~ y_t_minus_1), col = "red", lwd = 2)
#--------------------------------------------------------------------------------

# ACF & PACF TO CHECK FOR SEASONALITY
par(mfrow = c(1, 2))
seasonal_acf(drug.ts, lag.max = 64, main = "Fig 2. ACF plot of drug series")
seasonal_pacf(drug.ts, lag.max = 64, main = "Fig 3. PACF plot of drug series")
par(mfrow = c(1, 1))

# seasonality and existence of trend are obvious from the ACF and PACF plots
#--------------------------------------------------------------------------------

# CHECKING STATIONARITY AND NORMALITY
adf.test(drug.ts)
pp.test(drug.ts)

qqnorm(y = drug.ts, main = "Fig 4. QQ plot of drug series", cex.main = 1.8)
qqline(y = drug.ts, col = 2, lwd = 1, lty = 2)
shapiro.test(drug.ts)
#--------------------------------------------------------------------------------

# MODEL SPECIFICATION USING RESIDUAL APPROACH

# we're following the residual approach (not classical) to specify the SARIMA model

# first, we fit a plain model with only seasonal differencing (D = 1)
# this helps us see if the seasonal trend can be removed
# fitting arima model with no AR or MA terms, but with seasonal differencing raw data

# m1
m1.drug <- Arima(drug.ts, order = c(0, 0, 0), seasonal = list(order = c(0, 1, 0), period = 12))
res.m1 <- residuals(m1.drug)
par(mfrow = c(1, 1))
plot(res.m1, xlab = "Time", ylab = "Residuals", main = "Fig 5. Time Series Plot of Residuals (m1)")
par(mfrow = c(1, 2))
seasonal_acf(res.m1, lag.max = 48, main = "Fig 6. ACF of Residuals (m1)")
seasonal_pacf(res.m1, lag.max = 48, main = "Fig 7. PACF of Residuals (m1)")
par(mfrow = c(1, 1))

# based on residual ACF (Fig 6), we observe one significant red bar → Capital Q = 1
# based on residual PACF (Fig 7), we observe two significant red bars → Capital P = 2

# also, after applying seasonal differencing (D = 1), the significant spike at lag 1 in PACF disappeared
# this suggests that the seasonal differencing has already helped capture some of the seasonal structure

# so, now we will add the SARMA(2,1) component and see if we get rid of the seasonal component completely
#--------------------------------------------------------------------------------

# m2
m2.drug = Arima(drug.ts,order=c(0,0,0),seasonal=list(order=c(2,1,1), period=12))
res.m2 <- residuals(m2.drug)
par(mfrow = c(1, 1))
plot(res.m2, xlab = "Time", ylab = "Residuals", main = "Fig 8. Time Series Plot of Residuals (m2)")
par(mfrow = c(1, 2))
seasonal_acf(res.m2, lag.max = 48, main = "Fig 9. ACF of Residuals (m2)")
seasonal_pacf(res.m2, lag.max = 48, main = "Fig 10. PACF of Residuals (m2)")
par(mfrow = c(1, 1))

# from the ACF (Fig 9), there’s just one near miss — no clear seasonal pattern left overall
# but the PACF (Fig 10) still shows 2 strong seasonal red bars, meaning some seasonal effect remains

# instead of increasing P or Q again (which will just make the model more complex),
# we’ll now try adding an ordinary difference (d = 1) in m3
# this might help reduce the trend and remove the leftover seasonal spikes in PACF
#--------------------------------------------------------------------------------

# m3
m3.drug = Arima(drug.ts,order=c(0,1,0),seasonal=list(order=c(2,1,1), period=12))
res.m3 <- residuals(m3.drug)
par(mfrow = c(1, 1))
plot(res.m3, xlab = "Time", ylab = "Residuals", main = "Fig 11. Time Series Plot of Residuals (m3)")
par(mfrow = c(1, 2))
seasonal_acf(res.m3, lag.max = 48, main = "Fig 12. ACF of Residuals (m3)")
seasonal_pacf(res.m3, lag.max = 48, main = "Fig 13. PACF of Residuals (m3)")
par(mfrow = c(1, 1))

# after applying ordinary differencing (d = 1) in m3, we finally got rid of the seasonal lags
# both ACF (Fig 12) and PACF (Fig 13) no longer show any strong seasonal red bars — great!

# now, we’ll move on to specifying the ARIMA (non-seasonal) component
# in ACF (Fig 12), we see 3 significant bars before the first red seasonal lag — so we'll set q = 3
# in PACF (Fig 13), we see 4 significant bars before the first red seasonal lag — so we'll set p = 4
#--------------------------------------------------------------------------------

# m4
m4.drug = Arima(drug.ts,order=c(4,1,3),seasonal=list(order=c(2,1,1), period=12), method = "ML")
res.m4 <- residuals(m4.drug)
par(mfrow = c(1, 1))
plot(res.m4, xlab = "Time", ylab = "Residuals", main = "Fig 14. Time Series Plot of Residuals (m4)")
par(mfrow = c(1, 2))
seasonal_acf(res.m4, lag.max = 48, main = "Fig 15. ACF of Residuals (m4)")
seasonal_pacf(res.m4, lag.max = 48, main = "Fig 16. PACF of Residuals (m4)")
par(mfrow = c(1, 1))

# after fitting m4 with SARIMA(4,1,3)×(2,1,1)[12], the ACF (Fig 15) and PACF (Fig 16) look pretty good
# there are no strong seasonal spikes left — just a few small ones here and there, which we can treat as 
# random noise or mild outliers

# if we want to explore further, we can push the AR (p) value up to 6 since we saw about 6 smaller black 
# bars before the seasonal cut-off and ACF still shows up to 5 slight spikes — so MA (q) could go up to 5 too

# instead of guessing, we’ll now use the EACF and BIC plot to get a clearer idea of which AR and MA 
# combinations are worth trying next
#--------------------------------------------------------------------------------

# MODEL SPECIFICATION: EACF

# now we'll use the EACF (Extended Autocorrelation Function) to explore possible AR and MA orders
# we are using residuals of m3 because there is still leftover signal in them

eacf(res.m3)

# from the EACF table, we found the top-left 'o' at p = 0 and q = 5 — that’s our starting point
# now looking at its neighbors:
# now looking at its neighbors:
# → (p = 1, q = 5) — (one step below)
# → (p = 1, q = 4) — (diagonally left-down)
# → (p = 0, q = 6) — (just one step right)

# so the tentative models we’ll check are:
# SARIMA(0,1,5)x(2,1,1)[12]
# SARIMA(1,1,5)x(2,1,1)[12] 
# SARIMA(1,1,4)x(2,1,1)[12] 
# SARIMA(0,1,6)x(2,1,1)[12] 
#--------------------------------------------------------------------------------

# MODEL SPECIFICATION: BIC TABLE

# now let’s check the BIC table to see which AR and MA combinations are more promising

par(mfrow = c(1, 1))
bic_table = armasubsets(y = res.m3, nar = 8, nma = 8, y.name = "p", ar.method = "ols")
plot(bic_table)

# from the BIC plot, we can see two of the darkest cells (lowest BIC values) at:
# → p = 4, q = 1
# → p = 4, q = 2

# so the tentative models we’ll check are:
# SARIMA(4,1,1)x(2,1,1)[12]  
# SARIMA(4,1,2)x(2,1,1)[12]  
#--------------------------------------------------------------------------------

# FINAL SET OF POSSIBLE MODELS

# 1. SARIMA(4,1,3)×(2,1,1)[12]    (m4)
# 2. SARIMA(0,1,5)x(2,1,1)[12]    (EACF)
# 3. SARIMA(1,1,5)x(2,1,1)[12]    (EACF)
# 4. SARIMA(1,1,4)x(2,1,1)[12]    (EACF)
# 5. SARIMA(0,1,6)x(2,1,1)[12]    (EACF)
# 6. SARIMA(4,1,1)x(2,1,1)[12]    (BIC TABLE)
# 7. SARIMA(4,1,2)x(2,1,1)[12]    (BIC TABLE)
#--------------------------------------------------------------------------------

# FITTING SARIMA MODELS (ML AND CSS)

# 1. SARIMA(4,1,3)x(2,1,1)[12]
m_413.drugML <- Arima(drug.ts, order = c(4, 1, 3), seasonal = list(order = c(2, 1, 1), period = 12), method = "ML")
coeftest(m_413.drugML)
residual.analysis(model = m_413.drugML)

m_413.drugCSS <- Arima(drug.ts, order = c(4, 1, 3), seasonal = list(order = c(2, 1, 1), period = 12), method = "CSS")
coeftest(m_413.drugCSS)
residual.analysis(model = m_413.drugCSS)
#----------------------------------------


# 2. SARIMA(0,1,5)x(2,1,1)[12]
m_015.drugML <- Arima(drug.ts, order = c(0, 1, 5), seasonal = list(order = c(2, 1, 1), period = 12), method = "ML")
coeftest(m_015.drugML)
residual.analysis(model = m_015.drugML)

m_015.drugCSS <- Arima(drug.ts, order = c(0, 1, 5), seasonal = list(order = c(2, 1, 1), period = 12), method = "CSS")
coeftest(m_015.drugCSS)
residual.analysis(model = m_015.drugCSS)
#----------------------------------------


# 3. SARIMA(1,1,5)x(2,1,1)[12]
m_115.drugML <- Arima(drug.ts, order = c(1, 1, 5), seasonal = list(order = c(2, 1, 1), period = 12), method = "ML")
coeftest(m_115.drugML)
residual.analysis(model = m_115.drugML)

m_115.drugCSS <- Arima(drug.ts, order = c(1, 1, 5), seasonal = list(order = c(2, 1, 1), period = 12), method = "CSS")
coeftest(m_115.drugCSS)
residual.analysis(model = m_115.drugCSS)
#----------------------------------------


# 4. SARIMA(1,1,4)x(2,1,1)[12]
m_114.drugML <- Arima(drug.ts, order = c(1, 1, 4), seasonal = list(order = c(2, 1, 1), period = 12), method = "ML")
coeftest(m_114.drugML)
residual.analysis(model = m_114.drugML)

m_114.drugCSS <- Arima(drug.ts, order = c(1, 1, 4), seasonal = list(order = c(2, 1, 1), period = 12), method = "CSS")
coeftest(m_114.drugCSS)
residual.analysis(model = m_114.drugCSS)
#----------------------------------------


# 5. SARIMA(0,1,6)x(2,1,1)[12]
m_016.drugML <- Arima(drug.ts, order = c(0, 1, 6), seasonal = list(order = c(2, 1, 1), period = 12), method = "ML")
coeftest(m_016.drugML)
residual.analysis(model = m_016.drugML)

m_016.drugCSS <- Arima(drug.ts, order = c(0, 1, 6), seasonal = list(order = c(2, 1, 1), period = 12), method = "CSS")
coeftest(m_016.drugCSS)
residual.analysis(model = m_016.drugCSS)
#----------------------------------------


# 6. SARIMA(4,1,1)x(2,1,1)[12]
m_411.drugML <- Arima(drug.ts, order = c(4, 1, 1), seasonal = list(order = c(2, 1, 1), period = 12), method = "ML")
coeftest(m_411.drugML)
residual.analysis(model = m_411.drugML)

m_411.drugCSS <- Arima(drug.ts, order = c(4, 1, 1), seasonal = list(order = c(2, 1, 1), period = 12), method = "CSS")
coeftest(m_411.drugCSS)
residual.analysis(model = m_411.drugCSS)
#----------------------------------------


# 7. SARIMA(4,1,2)x(2,1,1)[12]
m_412.drugML <- Arima(drug.ts, order = c(4, 1, 2), seasonal = list(order = c(2, 1, 1), period = 12), method = "ML")
coeftest(m_412.drugML)
residual.analysis(model = m_412.drugML)

m_412.drugCSS <- Arima(drug.ts, order = c(4, 1, 2), seasonal = list(order = c(2, 1, 1), period = 12), method = "CSS")
coeftest(m_412.drugCSS)
residual.analysis(model = m_412.drugCSS)
#--------------------------------------------------------------------------------

# AIC AND BIC COMPARISON 

# storing AIC values for each SARIMA model
drug.aic <- AIC(m_413.drugML, m_015.drugML, m_115.drugML, m_114.drugML,
                m_016.drugML, m_411.drugML, m_412.drugML)

# storing BIC values for each SARIMA model
drug.bic <- BIC(m_413.drugML, m_015.drugML, m_115.drugML, m_114.drugML,
                m_016.drugML, m_411.drugML, m_412.drugML)

# sorting based on AIC and BIC
sort.score(drug.aic, score = "aic")
sort.score(drug.bic, score = "bic")
#--------------------------------------------------------------------------------

# ACCURACY METRICS FOR EACH MODEL

acc_413 <- accuracy(m_413.drugML)[1:7]
acc_015 <- accuracy(m_015.drugML)[1:7]
acc_115 <- accuracy(m_115.drugML)[1:7]
acc_114 <- accuracy(m_114.drugML)[1:7]
acc_016 <- accuracy(m_016.drugML)[1:7]
acc_411 <- accuracy(m_411.drugML)[1:7]
acc_412 <- accuracy(m_412.drugML)[1:7]

# combine all into a dataframe
drug_models_accuracy <- data.frame(
  rbind(acc_413, acc_015, acc_115, acc_114, acc_016, acc_411, acc_412)
)

# setting column and row names
colnames(drug_models_accuracy) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE", "ACF1")
rownames(drug_models_accuracy) <- c("SARIMA(4,1,3)x(2,1,1)[12]",
                                    "SARIMA(0,1,5)x(2,1,1)[12]",
                                    "SARIMA(1,1,5)x(2,1,1)[12]",
                                    "SARIMA(1,1,4)x(2,1,1)[12]",
                                    "SARIMA(0,1,6)x(2,1,1)[12]",
                                    "SARIMA(4,1,1)x(2,1,1)[12]",
                                    "SARIMA(4,1,2)x(2,1,1)[12]")

# view rounded accuracy metrics
round(drug_models_accuracy, 3)
#--------------------------------------------------------------------------------

# OVER PARAMETERISATION 

# Best model identified based on AIC, BIC, and forecast accuracy:
# --> SARIMA(4,1,1)×(2,1,1)[12]

# To test whether adding more parameters improves the model,
# we will check over-parameterised versions:
#   1. SARIMA(5,1,1)×(2,1,1)[12]
#   2. SARIMA(4,1,2)×(2,1,1)[12] --> already tested earlier

# 8. SARIMA(5,1,1)×(2,1,1)[12] 
m_511.drugML <- Arima(drug.ts, order = c(5, 1, 1), seasonal = list(order = c(2, 1, 1), period = 12), method = "ML")
coeftest(m_511.drugML)
residual.analysis(model = m_511.drugML)

m_511.drugCSS <- Arima(drug.ts, order = c(5, 1, 1), seasonal = list(order = c(2, 1, 1), period = 12), method = "CSS")
coeftest(m_511.drugCSS)
residual.analysis(model = m_511.drugCSS)
#--------------------------------------------------------------------------------

# AIC & BIC COMPARISON INCLUDING SARIMA(5,1,1)×(2,1,1)[12]

# storing AIC values for all 8 models
drug.aic <- AIC(m_413.drugML, m_015.drugML, m_115.drugML, m_114.drugML,
                m_016.drugML, m_411.drugML, m_412.drugML, m_511.drugML)

# storing BIC values for all 8 models
drug.bic <- BIC(m_413.drugML, m_015.drugML, m_115.drugML, m_114.drugML,
                m_016.drugML, m_411.drugML, m_412.drugML, m_511.drugML)

# sorting by AIC and BIC
sort.score(drug.aic, score = "aic")
sort.score(drug.bic, score = "bic")
#--------------------------------------------------------------------------------

# ACCURACY METRICS FOR EACH MODEL

acc_413 <- accuracy(m_413.drugML)[1:7]
acc_015 <- accuracy(m_015.drugML)[1:7]
acc_115 <- accuracy(m_115.drugML)[1:7]
acc_114 <- accuracy(m_114.drugML)[1:7]
acc_016 <- accuracy(m_016.drugML)[1:7]
acc_411 <- accuracy(m_411.drugML)[1:7]
acc_412 <- accuracy(m_412.drugML)[1:7]
acc_511 <- accuracy(m_511.drugML)[1:7]

# combining all accuracy metrics
drug_models_accuracy <- data.frame(
  rbind(acc_413, acc_015, acc_115, acc_114, acc_016, acc_411, acc_412, acc_511)
)

# setting column and row names
colnames(drug_models_accuracy) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE", "ACF1")
rownames(drug_models_accuracy) <- c("SARIMA(4,1,3)x(2,1,1)[12]",
                                    "SARIMA(0,1,5)x(2,1,1)[12]",
                                    "SARIMA(1,1,5)x(2,1,1)[12]",
                                    "SARIMA(1,1,4)x(2,1,1)[12]",
                                    "SARIMA(0,1,6)x(2,1,1)[12]",
                                    "SARIMA(4,1,1)x(2,1,1)[12]",
                                    "SARIMA(4,1,2)x(2,1,1)[12]",
                                    "SARIMA(5,1,1)x(2,1,1)[12]")

# display rounded accuracy table
round(drug_models_accuracy, 3)
#--------------------------------------------------------------------------------

# FINAL FORECAST USING BEST MODEL: SARIMA(4,1,1)×(2,1,1)[12]

# fitting the best model one last time to ensure it's fresh in memory
best.model <- Arima(drug.ts, order = c(4, 1, 1), seasonal = list(order = c(2, 1, 1), period = 12), method = "ML")

# generating forecast for the next 10 months
drug.forecast <- forecast(best.model, h = 10)

# plotting forecast results
autoplot(drug.forecast) +
  ggtitle("Fig 27. Forecast of Monthly PBS Corticosteroid Expenditure (Next 10 Months)") +
  xlab("Time") +
  ylab("Expenditure (in Million AUD)") +
  theme_minimal()

# printing forecast values
print(drug.forecast)                   
#--------------------------------------------------------------------------------
