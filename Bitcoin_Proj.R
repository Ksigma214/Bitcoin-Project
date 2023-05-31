#%============================================================================%=
####                                  Setup                                 ####
#%============================================================================%=

# install packages and load libraries
rm(list=ls())

pkgs <- c('TSA', 'tseries' ,'fUnitRoots', 'FSAdata','forecast', 'dynlm','xtable','stargazer','png','magick',
          'lmtest','dplyr','kableExtra','AID','nortest','xts','uroot','vars','ggplot2','moments','tsBSS',
          'urca','MASS','knitr','sandwich','grangers','lubridate','FitAR','gridExtra')

new.pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(new.pkgs)) install.packages(new.pkgs)
invisible(lapply(pkgs, require, character.only = T))

sort.score <- function(x, score = c("bic", "aic")){
  if (score == "aic"){
    x[with(x, order(AIC)),]
  } else if (score == "bic") {
    x[with(x, order(BIC)),]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
}
residual.analysis <- function(model, std = TRUE, start = 2) {
  library(TSA)
  library(FitAR)
  if (std == TRUE) {
    res.model <- rstandard(model)
  } else {
    res.model <- residuals(model)
  }
  par(mfrow = c(3, 2))
  # Time series plot of standardised residuals
  plot(res.model, type = 'o', ylab = 'Standardised residuals', main = "Time series plot of standardised residuals")
  abline(h = 0)
  # Histogram of standardised residuals with fitted normal curve
  hist(res.model, main = "Histogram of standardised residuals", prob = TRUE)
  curve(dnorm(x, mean = mean(res.model), sd = sd(res.model)), add = TRUE, col = 'red')
  # PACF of standardised residuals
  pacf(res.model, lag.max = 20, main = "PACF of standardised residuals")
  # ACF of standardised residuals
  acf(res.model, lag.max = 20, main = "ACF of standardised residuals")
  # QQ plot of standardised residuals
  qqnorm(res.model, main = "QQ plot of standardised residuals")
  qqline(res.model, col = 2)
  # Shapiro-Wilk normality test
  print(shapiro.test(res.model))
  # Komogrov test
  print(ks.test(res.model, rnorm))
  # jarque.bera Test
  print(jarque.bera.test(res.model))
  # Ljung-Box Q-Test
  k <- 0
  LBQPlot(res.model, lag.max = length(model$residuals) - 50, StartLag = k + 1, k = 0, SquaredQ = FALSE)
  par(mfrow = c(1, 1))
}

setwd("D:/桌面/M1 APE/2nd Semester/Metrics2/Bitcoin Project")
list.files() 

#%============================================================================%=
                      #### PART I: Univariate analysis ####   
#%============================================================================%=

#%============================================================================%=
                  ####    (a) Data presentation      ####
#%============================================================================%=

#%--------------- Data Cleaning and Defining ---------------#%
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")

# S&P
sp <- read.csv2("S&P 500 Historical Data_monthly.csv",header = TRUE, sep = ",",dec = ".",nrows = 152)
sp$Date <- as.Date(sp$Date, format = "%m/%d/%Y")
sp <- arrange(sp, Date)
sp_xts <- xts(sp[,"Price"], order.by = sp$Date)
sp_price <- ts(sp_xts, frequency = 12)
lsp500 <- log(sp_price)
dlsp500 <- diff(lsp500)
# Bitcoin
bc <- read.csv2("Bitcoin Historical Data_monthly.csv",header = TRUE, sep = ",",dec = ".",nrows = 152)
bc$Date <- as.Date(bc$Date, format = "%Y/%m/%d")
bc <- arrange(bc, Date)
bc_xts <- xts(bc[,2], order.by = bc$Date)
bc_price <- ts(bc_xts, frequency = 12)
lbitcoin <- log(bc_price)
dlbitcoin <- diff(lbitcoin)
# Gold
xu <- read.csv2("XAU_USD Historical Data_monthly.csv",header = TRUE, sep = ",",dec = ".",nrows = 152)
xu$Date <- as.Date(xu$Date, format = "%m/%d/%Y")
xu <- arrange(xu, Date)
xu_xts <- xts(xu[,2], order.by = xu$Date)
xu_price <- ts(xu_xts,frequency = 12)
lgold <- log(xu_price)
dlgold <- diff(lgold)

#% ----------------- Data Presentation --------------------%#

# The graph for the raw data
par(mfrow=c(1,3))
plot(sp_price,main = "Raw - S&P 500")
plot(bc_price,main = "Raw - Bitcoin")
plot(xu_price,main = "Raw - Gold")

# summary table data
## sp
summary(sp_price)
skewness(sp_price)
kurtosis(sp_price)
## bc
summary(bc_price)
skewness(bc_price)
kurtosis(bc_price)
## xu
summary(xu_price)
skewness(xu_price)
kurtosis(xu_price)

#####    We will do (b),(c) and (d) for 3 time series separately
#%============================================================================%=
             ####     (b) Unit Root Test -  S&P 500     #### 
#%============================================================================%=

# Decomposition of the raw data - S&P 500
sp_decomp_gr <- decompose(log(sp_price))
plot(sp_decomp_gr)

# ACF,PACF for log S&P 500 and diff log S&P 500
par(mfrow=c(1,2))
Pacf(lsp500,lag.max = 20) 
Acf(lsp500,lag.max = 20) 

#%-------------------For Log S&P 500----------------------%#
# To select the optimal number of lags :

# For log data, use the information criteria 
VARselect(lsp500, lag.max=24) # monthly data, we use lag.max=24 # 1 or 2 or 3

# ADF-test
summary(ur.df(lsp500, type="trend", lags=1)) 
summary(ur.df(lsp500, type="trend", lags=2))
summary(ur.df(lsp500, type="trend", lags=3))

# ERS:
# ERS (DF-GLS) Lag=1 : 
summary(ur.ers(lsp500, model="trend", type="DF-GLS", lag.max=1)) 
# ERS (P-test) Lag=1 :
summary(ur.ers(lsp500, model="trend", type="P-test", lag.max=1)) 
# ERS (DF-GLS) Lag=2 : 
summary(ur.ers(lsp500, model="trend", type="DF-GLS", lag.max=2)) 
# ERS (P-test) Lag=2 :
summary(ur.ers(lsp500, model="trend", type="P-test", lag.max=2)) 
# ERS (DF-GLS) Lag=3 : 
summary(ur.ers(lsp500, model="trend", type="DF-GLS", lag.max=3)) 
# ERS (P-test) Lag=3 :
summary(ur.ers(lsp500, model="trend", type="P-test", lag.max=3)) 
# PP:
summary(ur.pp(lsp500, model="trend", type="Z-tau", use.lag=1)) 
summary(ur.pp(lsp500, model="trend", type="Z-tau", use.lag=2)) 
summary(ur.pp(lsp500, model="trend", type="Z-tau", use.lag=3))

# KPSS:
summary(ur.kpss(lsp500, type="tau", lags="short")) #Lags = 4
summary(ur.kpss(lsp500, type="tau", lags="long")) #Lags = 13

#%-------------------For Diff Log S&P 500----------------------%#
par(mfrow=c(1,1))
plot(dlsp500)
# For diff log data, use the information criteria (works for all of these tests except KPSS)
VARselect(dlsp500, lag.max=24) # monthly data, we use lag.max=24 # 1 or 2
# ADF-test
summary(ur.df(dlsp500, type="drift", lags=1)) 
summary(ur.df(dlsp500, type="drift", lags=2))
# 2) PP test 
#- H0: the time series data contains a unit root, indicating that the series is non-stationary
summary(ur.pp(dlsp500, model="constant", type="Z-tau", use.lag=1)) 
summary(ur.pp(dlsp500, model="constant", type="Z-tau", use.lag=2)) 
# 3) ERS test
# ERS (DF-GLS) : 
summary(ur.ers(dlsp500, model="constant", type="DF-GLS", lag.max=1)) 
summary(ur.ers(dlsp500, model="constant", type="DF-GLS", lag.max=2)) 
# ERS (P-test) :
summary(ur.ers(dlsp500, model="constant", type="P-test", lag.max=1)) 
summary(ur.ers(dlsp500, model="constant", type="P-test", lag.max=2)) 
# 4) KPSS test
summary(ur.kpss(dlsp500, type="mu", lags="short")) #4 Lags
summary(ur.kpss(dlsp500, type="mu", lags="long")) #13 lags

#%============================================================================%=
 ####     (c) Identification of the ARMA or ARIMA process  -  S&P 500     #### 
#%============================================================================%=

par(mfrow=c(1,2))
Pacf(dlsp500,lag.max = 20) 
Acf(dlsp500,lag.max = 20)
# arima(1,0,0) for dlsp500
arima_dlsp500_10 = Arima(dlsp500, order=c(1,0,0))
coeftest(arima_dlsp500_10)
# arima(0,0,1) for dlsp500
arima_dlsp500_01 = Arima(dlsp500, order=c(0,0,1))
coeftest(arima_dlsp500_01)
# arima(1,0,1) for dlsp500
arima_dlsp500_11 = Arima(dlsp500, order=c(1,0,1))
coeftest(arima_dlsp500_11)
# Candidates are arima(1,1,0) and arima(0,1,1)
residual.analysis(arima_dlsp500_10,std=F)
residual.analysis(arima_dlsp500_01,std=F)
# Then we choose between AR(1) and MA(1)
sort.score(AIC(arima_dlsp500_01,arima_dlsp500_10),score = "aic")
sort.score(BIC(arima_dlsp500_01,arima_dlsp500_10),score = "bic")
auto.arima(lsp500,
           stepwise=FALSE, approximation=FALSE) # Arima(0,1,1)


#%============================================================================%=
  ####     (d) Forecasts: in-sample and out-of-sample -  S&P 500     #### 
#%============================================================================%=

# In-sample forecasting (= predictions = fitted values)
## for dlsp500
f_dlsp500_is = arima_dlsp500_01$fitted
# Out-of-sample forecasting (= real life)
## Short horizon
f_dlsp500_oos = forecast(object=dlsp500, model=arima_dlsp500_01, h=6) 
## check 6 months
f_dlsp500_oos = forecast(object=dlsp500, model=arima_dlsp500_01, h=6) 
dlsp500_s <- window(dlsp500,end=2022.9)
f_dlsp500_oos_s = forecast(object=dlsp500_s, model=arima_dlsp500_01, h=6)
## Long horizon
f_dlsp500_oos = forecast(object=dlsp500, model=arima_dlsp500_01, h=12) 
## check 12 months
f_dlsp500_oos = forecast(object=dlsp500, model=arima_dlsp500_01, h=12) 
dlsp500_l <- window(dlsp500, end=2022.4)
f_dlsp500_oos_l = forecast(object=dlsp500_l, model=arima_dlsp500_01, h=12)
## Graph for dlsp
par(mfrow=c(5,1))
plot(cbind(arima_dlsp500_01$x,f_dlsp500_is), plot.type="s", col=c("black","red"), ylab = "diff-log S&P500",main="In-sample forecast for log S&P500 ARIMA(0,1,1)")
legend(x=0, y=-0.2, legend=c("Data","In-sample forecast"),
       col=c("black","red"),lwd=2)
plot(f_dlsp500_oos,main="Out-of-sample forecast 6 months for log S&P500 ARIMA(0,1,1)", plot.type="s", ylab = "diff-log S&P500")
plot(f_dlsp500_oos_s,main="Check out-of-sample forecast 6 months for log S&P500 ARIMA(0,1,1)", plot.type="s", ylab = "diff-log S&P500")
lines(dlsp500, col="black")
plot(f_dlsp500_oos,main="Out-of-sample forecast 12 months for log S&P500 ARIMA(0,1,1)", plot.type="s", ylab = "diff-log S&P500")
plot(f_dlsp500_oos_l,main="Check Out-of-sample forecast 12 months for log S&P500 ARIMA(0,1,1)", plot.type="s", ylab = "diff-log S&P500")
lines(dlsp500, col="black")


#%============================================================================%=
             ####     (b) Unit Root Test -  Bitcoin     #### 
#%============================================================================%=

# Decomposition of the raw data - Bitcoin
bc_decomp_gr <- decompose(log(bc_price))
plot(bc_decomp_gr)
# ACF,PACF for log Bitcoin and diff log Bitcoin
par(mfrow=c(1,2))
Pacf(lbitcoin,lag.max = 20) 
Acf(lbitcoin,lag.max = 20)
#%-------------------For Log Bitcoin----------------------%#
VARselect(lbitcoin, lag.max=24) # monthly data, we use lag.max=24 # 1 or 2

# ADF-test
summary(ur.df(lbitcoin, type="trend", lags=1)) 
summary(ur.df(lbitcoin, type="trend", lags=2))

# ERS:
# ERS (DF-GLS) Lag=1 : 
summary(ur.ers(lbitcoin, model="trend", type="DF-GLS", lag.max=1)) 
# ERS (P-test) Lag=1 :
summary(ur.ers(lbitcoin, model="trend", type="P-test", lag.max=1)) 
# ERS (DF-GLS) Lag=2 : 
summary(ur.ers(lbitcoin, model="trend", type="DF-GLS", lag.max=2)) 
# ERS (P-test) Lag=2 :
summary(ur.ers(lbitcoin, model="trend", type="P-test", lag.max=2)) 
# PP:
summary(ur.pp(lbitcoin, model="trend", type="Z-tau", use.lag=1)) 
summary(ur.pp(lbitcoin, model="trend", type="Z-tau", use.lag=2)) 
# KPSS:
summary(ur.kpss(lbitcoin, type="tau", lags="short")) #Lags = 4
summary(ur.kpss(lbitcoin, type="tau", lags="long")) #Lags = 13

par(mfrow=c(1,1))
plot(dlbitcoin)
#%-------------------For Diff Log Bitcoin----------------------%#
VARselect(lgold, lag.max=24) # monthly data, we use lag.max=24 # 1 lag with all criteria

# ADF-test:
summary(ur.df(dlbitcoin, type="drift", lags=1)) 

# 2) PP test - H0: the time series data contains a unit root, indicating that the series is non-stationary
summary(ur.pp(dlbitcoin, model="constant", type="Z-tau", use.lag=1))
# 3) ERS test
# ERS (DF-GLS) : 
summary(ur.ers(dlbitcoin, model="constant", type="DF-GLS", lag.max=1)) 
# ERS (P-test) :
summary(ur.ers(dlbitcoin, model="constant", type="P-test", lag.max=1)) 
# 4) KPSS test
summary(ur.kpss(dlbitcoin, type="mu", lags="short")) #4 Lags
summary(ur.kpss(dlbitcoin, type="mu", lags="long")) #13 lags

#%============================================================================%=
####     (c) Identification of the ARMA or ARIMA process  -  Bitcoin     #### 
#%============================================================================%=
par(mfrow=c(1,2))
Pacf(dlbitcoin,lag.max = 20) 
Acf(dlbitcoin,lag.max = 20)
# arima(0,0,8) for dlbitcoin
arima_dlbitcoin_08 = Arima(dlbitcoin, order=c(0,0,8))
coeftest(arima_dlbitcoin_08)
# arima(1,0,0) for dlbitcoin
arima_dlbitcoin_10 = Arima(dlbitcoin, order=c(1,0,0))
coeftest(arima_dlbitcoin_10)
# arima(1,0,8) for dlbitcoin
arima_dlbitcoin_18 = Arima(dlbitcoin, order=c(1,0,8))
coeftest(arima_dlbitcoin_18)
# arima(0,0,1) for dlbitcoin
arima_dlbitcoin_01 = Arima(dlbitcoin, order=c(0,0,1))
coeftest(arima_dlbitcoin_01)
# arima(1,0,1) for dlbitcoin
arima_dlbitcoin_11 = Arima(dlbitcoin, order=c(1,0,1))
coeftest(arima_dlbitcoin_11)
# Candidates are ARIMA(1,1,0), ARIMA(0,1,1) and ARIMA(0,1,8)
residual.analysis(arima_dlbitcoin_10,std=F)
residual.analysis(arima_dlbitcoin_08,std=F)
residual.analysis(arima_dlbitcoin_01,std=F)
## Then we choose between AR(1), MA(1), MA(8)
sort.score(AIC(arima_dlbitcoin_08,arima_dlbitcoin_10,arima_dlbitcoin_01),score = "aic")
sort.score(BIC(arima_dlbitcoin_08,arima_dlbitcoin_10,arima_dlbitcoin_01),score = "bic")
auto.arima(dlbitcoin,stepwise=FALSE, approximation=FALSE) 

#%============================================================================%=
####     (d) Forecasts: in-sample and out-of-sample -  Bitcoin     #### 
#%============================================================================%=

# In-sample forecasting (= predictions = fitted values)
## for dlbitcoin
f_dlbitcoin_is = arima_dlbitcoin_10$fitted
# Out-of-sample forecasting (= real life)
## Short horizon
f_dlbitcoin_oos = forecast(object=dlbitcoin, model=arima_dlbitcoin_10, h=6) 
## check 6 months
f_dlbitcoin_oos = forecast(object=dlbitcoin, model=arima_dlbitcoin_10, h=6) 
dlbitcoin_s <- window(dlbitcoin,end=2022.9)
f_dlbitcoin_oos_s = forecast(object=dlbitcoin_s, model=arima_dlbitcoin_10, h=6)
## Long horizon
f_dlbitcoin_oos = forecast(object=dlbitcoin, model=arima_dlbitcoin_10, h=12) 
## check 12 months
f_dlbitcoin_oos = forecast(object=dlbitcoin, model=arima_dlbitcoin_10, h=12) 
dlbitcoin_l <- window(dlbitcoin, end=2022.4)
f_dlbitcoin_oos_l = forecast(object=dlbitcoin_l, model=arima_dlbitcoin_10, h=12)
## Graph for Bitcoin
par(mfrow=c(5,1))
plot(cbind(arima_dlbitcoin_10$x,f_dlbitcoin_is), plot.type="s", col=c("black","red"), ylab = "diff-log Bitcoin",main="In-sample forecast for log Bitcoin ARIMA(1,1,0)")
legend(x=0, y=-0.2, legend=c("Data","In-sample forecast"),
       col=c("black","red"),lwd=2)
plot(f_dlbitcoin_oos,main="Out-of-sample forecast 6 months for log Bitcoin ARIMA(1,1,0)", plot.type="s", ylab = "diff-log Bitcoin")
plot(f_dlbitcoin_oos_s,main="Check out-of-sample forecast 6 months for log Bitcoin ARIMA(1,1,0)", plot.type="s", ylab = "diff-log Bitcoin")
lines(dlbitcoin, col="black")
plot(f_dlbitcoin_oos,main="Out-of-sample forecast 12 months for log Bitcoin ARIMA(1,1,0)", plot.type="s", ylab = "diff-log Bitcoin")
plot(f_dlbitcoin_oos_l,main="Check Out-of-sample forecast 12 months for log Bitcoin ARIMA(1,1,0)", plot.type="s", ylab = "diff-log Bitcoin")
lines(dlbitcoin, col="black") 

#%============================================================================%= 
              ####     (b) Unit Root Test -  Gold    #### 
#%============================================================================%=

# Decomposition of raw gold data
xu_decomp_gr <- decompose(log(xu_price))
plot(xu_decomp_gr)
# ACF,PACF for log Bitcoin and diff log Bitcoin
plot(xu_decomp_gr)
par(mfrow=c(1,2))
Pacf(lgold,lag.max = 20)
Acf(lgold,lag.max = 20) 
#%-------------------For Log Gold----------------------%#
VARselect(dlbitcoin, lag.max=24) # monthly data, we use lag.max=24 # 1 lag with all criteria

# We justify the usage of trend as type in tests for Log gold
par(mfrow=c(1,1))
time <- 1:length(lgold)
# Fit a linear regression model with the time trend
modeltest <- lm(lgold ~ time)
summary(modeltest)
# Get the fitted values from the model
fitted_values <- fitted(modeltest)
# Plot the original time series and the fitted values
plot(time, lgold, type = "l", col = "blue", xlab = "Time", ylab = "Time Series Data")
lines(time, fitted_values, col = "red")
# Add a legend
legend("topleft", legend = c("Original", "Fitted"), col = c("blue", "red"),lty=1)

# ADF-test
summary(ur.df(lgold, type="trend", lags=1)) 

# ERS:
# ERS (DF-GLS) Lag=1 : 
summary(ur.ers(lgold, model="trend", type="DF-GLS", lag.max=1)) 
# ERS (P-test) Lag=1 :
summary(ur.ers(lgold, model="trend", type="P-test", lag.max=1)) 
# PP:
summary(ur.pp(lgold, model="trend", type="Z-tau", use.lag=1)) 
# KPSS:
summary(ur.kpss(lgold, type="tau", lags="short")) #Lags = 4
summary(ur.kpss(lgold, type="tau", lags="long")) #Lags = 13
# KPSS:
summary(ur.kpss(lgold, type="mu", lags="short")) #Lags = 4
summary(ur.kpss(lgold, type="mu", lags="long")) #Lags = 13

#%-------------------For Diff Log Gold----------------------%#
par(mfrow=c(1,1))
plot(dlgold)
VARselect(dlgold, lag.max=24) # monthly data, we use lag.max=24 # 1 lag with all criteria
# ADF-test:
summary(ur.df(dlgold, type="drift", lags=1)) # we reject H0 (UR) at 1% level, for whichever value of alpha.
# 2) PP test - H0: the time series data contains a unit root, indicating that the series is non-stationary
summary(ur.pp(dlgold, model="constant", type="Z-tau", use.lag=1)) 
# 3) ERS test
# ERS (DF-GLS) : 
summary(ur.ers(dlgold, model="constant", type="DF-GLS", lag.max=1)) 
# ERS (P-test) :
summary(ur.ers(dlgold, model="constant", type="P-test", lag.max=1))
# 4) KPSS test
summary(ur.kpss(dlgold, type="mu", lags="short")) #4 Lags
summary(ur.kpss(dlgold, type="mu", lags="long")) #13 lags

#%============================================================================%=
####     (c) Identification of the ARMA or ARIMA process  -  Gold     #### 
#%============================================================================%=
par(mfrow=c(1,2))
Pacf(dlgold,lag.max = 20)
Acf(dlgold,lag.max = 20)
# arima(0,0,0) for gold
arima_dlgold = Arima(dlgold, order=c(0,0,0))
coeftest(arima_dlgold)
auto.arima(lgold,stepwise=FALSE, approximation=FALSE) # Arima(0,1,0)
par(mar = c(2, 4, 3.5, 2))
residual.analysis(arima_dlgold,std=F)

auto.arima(dlgold,stepwise=FALSE, approximation=FALSE) # Arima(0,1,0)

#%============================================================================%=
  ####     (d) Forecasts: in-sample and out-of-sample -  Gold     #### 
#%============================================================================%=

# In-sample forecasting (= predictions = fitted values)
## for dlgold
f_dlgold_is = arima_dlgold$fitted
# Out-of-sample forecasting (= real life)
## Short horizon
f_dlgold_oos = forecast(object=dlgold, model=arima_dlgold, h=6) 
## check 6 months
f_dlgold_oos = forecast(object=dlgold, model=arima_dlgold, h=6) 
dlgold_s <- window(dlgold,end=2022.9)
f_dlgold_oos_s = forecast(object=dlgold_s, model=arima_dlgold, h=6)
## Long horizon
f_dlgold_oos = forecast(object=dlgold, model=arima_dlgold, h=12) 
## check 12 months
f_dlgold_oos = forecast(object=dlgold, model=arima_dlgold, h=12) 
dlgold_l <- window(dlgold, end=2022.4)
f_dlgold_oos_l = forecast(object=dlgold_l, model=arima_dlgold, h=12)
## Graph for Gold
par(mfrow=c(5,1))
plot(cbind(arima_dlgold$x,f_dlgold_is), plot.type="s", col=c("black","red"), ylab = "diff-log Gold",main="In-sample forecast for log Gold ARIMA(0,1,0)")
legend(x=0, y=-0.2, legend=c("Data","In-sample forecast"),
       col=c("black","red"),lwd=2)
plot(f_dlgold_oos,main="Out-of-sample forecast 6 months for log Gold ARIMA(0,1,0)", plot.type="s", ylab = "diff-log Gold")
plot(f_dlgold_oos_s,main="Check out-of-sample forecast 6 months for log Gold ARIMA(0,1,0)", plot.type="s", ylab = "diff-log Gold")
lines(dlgold, col="black")
plot(f_dlgold_oos,main="Out-of-sample forecast 12 months for log Gold ARIMA(0,1,0)", plot.type="s", ylab = "diff-log Gold")
plot(f_dlgold_oos_l,main="Check Out-of-sample forecast 12 months for log Gold ARIMA(0,1,0)", plot.type="s", ylab = "diff-log Gold")
lines(dlgold, col="black")


#%============================================================================%=
   #### PART II:  Multivariate analysis of the stationary variables ####
#%============================================================================%=

# Import data
vardata <-  read.csv2("VAR.csv",dec=".", sep=",")
vardata$Date <- as.Date(vardata$Date, format = "%m/%d/%Y")
vardata_xts <- xts(vardata[,-1], order.by = vardata$Date)
vardata = ts(vardata[,-1],start=c(2010,10), end=c(2023,5), frequency=12)

sp <- ts(vardata_xts[,"SP"], start=c(2010,10), frequency = 12)
dlsp500 <-  diff(log(sp))

btc <- ts(vardata[,"BTC"], start=c(2010,10), frequency = 12)
dlbitcoin = diff(log(btc))

xu <- ts(vardata[,"XU"], start=c(2010,10), frequency = 12)
dlgold = diff(log(xu))

par(mfrow=c(2,1))
par(mar = c(1, 4, 1, 2) + 0.1)
# The two surveys are advanced indicators of the industrial production.
plot(vardata[,c("SP","BTC","XU")], plot.type="s", col=c("black","red","green"),ylab="Price")
legend(x=2010, y=65000, legend=c("S&P 500","Bitcoin","Gold"),
       col=c("black","red","green"),lwd=2)
# To see it more clearly:
plot(scale(vardata[,c("SP","BTC","XU")]), plot.type="s", col=c("black","red","green"),ylab="Price_Scaled")
legend(x=2010, y=4, legend=c("S&P 500","Bitcoin","Gold"),
       col=c("black","red","green"),lwd=2)

#==============================================================================#
           # Selection of the lag order of the multivariate model #
#==============================================================================#

# According to PART I, we will include dlsp,dlbitcoin,dlgold in the VAR model.
dataVAR_balanced = cbind(dlgold,dlsp500,dlbitcoin)
VARselect(dataVAR_balanced, lag.max=24) # monthly data, we use lag.max=24
                                        # 1 lag with all criteria

#==============================================================================#
           # Example of estimation of the multivariate model #
#==============================================================================#

# Significance of the last lag and the constant
VAR1 = VAR(dataVAR_balanced, p=1)
summary(VAR1)
res_VAR1 <- residuals(VAR1)

# Plot inverse roots
par(mfrow = c(1, 1))
par(mar=c(2,4,3,2))
plot(Re(roots(VAR1)), Im(roots(VAR1)), xlim = c(-1, 1), ylim = c(-1, 1), xlab = "Real", ylab = "Imaginary", main = "Inverse Roots of AR Characteristic Polynomial", type = "p", col = "blue")
abline(h = 0, v = 0, lty = 2, col = "gray")
points(cos(seq(0, 2 * pi, length.out = 100)), sin(seq(0, 2 * pi, length.out = 100)), type = "l", col = "gray")

# Non-autocorrelated residuals
serial.test(VAR1, lags.bg = 4, type = "BG")
serial.test(VAR1, lags.bg = 8, type = "BG")
serial.test(VAR1, lags.bg = 12, type = "BG")
serial.test(VAR1, lags.bg = 24, type = "BG")
serial.test(VAR1, lags.bg = 30, type = "BG")

serial.test(VAR1, lags.pt=4, type="PT.asymptotic") 
serial.test(VAR1, lags.pt=8, type="PT.asymptotic")
serial.test(VAR1, lags.pt=12, type="PT.asymptotic") 
serial.test(VAR1, lags.pt=24, type="PT.asymptotic") 
serial.test(VAR1, lags.pt=30, type="PT.asymptotic") 

serial.test(VAR1, lags.pt=4, type="PT.adjusted") 
serial.test(VAR1, lags.pt=8, type="PT.adjusted") 
serial.test(VAR1, lags.pt=12, type="PT.adjusted") 
serial.test(VAR1, lags.pt=24, type="PT.adjusted") 
serial.test(VAR1, lags.pt=30, type="PT.adjusted") 

lbtest(res_VAR1,4, type = "linear")
lbtest(res_VAR1,8, type = "linear")
lbtest(res_VAR1,12, type = "linear")
lbtest(res_VAR1,24, type = "linear")
lbtest(res_VAR1,30, type = "linear")
# The residuals are not autocorrelated 

# Normally distributed residuals?
normality.test(VAR1) # the residuals are not normally distributed

#==============================================================================# 
                       #### Granger non-causality Test ####
#==============================================================================#

# The Granger causality tests presented in the course are all in the framework of ***Gaussian*** VAR models.
# We have seen that the residuals of our VAR(1) are not normally distributed.
# Thus, we should actually rely on ***bootstrap*** so as to work with the most appropriate critical values.

# Granger causality ("standard" or instantaneous) between two subvectors

# Between (dlbitcoin,dlgold) and dlsp. Here, the 2nd subvector includes only 1 variable:
# Normality 
causality(VAR1, cause="dlsp500") 
# Bootstrap
causality(VAR1, cause="dlsp500", boot=TRUE, boot.runs=1000) 
# Normality
causality(VAR1, cause=c("dlbitcoin","dlgold")) 
# Bootstrap
causality(VAR1, cause=c("dlbitcoin","dlgold"), boot=TRUE, boot.runs=1000) 

# Between (dlsp,dlgold) and dlbitcoin. Here, the 2nd subvector includes only 1 variable.
# Normality
causality(VAR1, cause="dlbitcoin") 
# Bootstrap
causality(VAR1, cause="dlbitcoin", boot=TRUE, boot.runs=1000) 
# Normality
causality(VAR1, cause=c("dlsp500","dlgold"))  
# Bootstrap
causality(VAR1, cause=c("dlsp500","dlgold"), boot=TRUE, boot.runs=1000) 

# Between (dlsp,dlbitcoin) and dlgold. Here, the 2nd subvector includes only 1 variable.
# Normality
causality(VAR1, cause="dlgold")
# Bootstap
causality(VAR1, cause="dlgold", boot=TRUE, boot.runs=1000)
# Normality
causality(VAR1, cause=c("dlsp500","dlbitcoin")) 
# Bootstrap
causality(VAR1, cause=c("dlsp500","dlbitcoin"), boot=TRUE, boot.runs=1000) 

# Partial Granger Test:
# "Standard" - Under the assumption of normality (likely irrelevant in our case), for a given equation of the VAR:
# If p=1, test that *the* coefficient of the variable x is null in the equation of the variable y, we can directly read it with the result of the individual t-test.
summary(VAR1)

# "Instantaneous"
gic_eq_dlgold = dynlm(dlgold ~ dlsp500 + dlbitcoin + L(dataVAR_balanced,1))
summary(gic_eq_dlgold) 

gic_eq_dlsp500 = dynlm(dlsp500 ~ dlbitcoin + dlgold + L(dataVAR_balanced,1))
summary(gic_eq_dlsp500)         

gic_eq_dlbitcoin = dynlm(dlbitcoin ~ dlsp500 + dlgold + L(dataVAR_balanced,1))
summary(gic_eq_dlbitcoin)  



# Inference on conditional Granger-causality
# x <- dlsp
Granger.inference.conditional(dlsp500, dlbitcoin, dlgold, ic.chosen = c("AIC","HQ","SC","FPE"),
                              max.lag = min(4, length(dlsp500) - 1), plot = T, type.chosen = "const",
                              p1 = 1, p2 = 1, nboots = 1000, conf = 0.95, bp = NULL,
                              ts_boot = 1)
Granger.inference.conditional(dlsp500, dlbitcoin, dlgold, ic.chosen = c("AIC","HQ","SC","FPE"),
                              max.lag = min(4, length(dlsp500) - 1), plot = T, type.chosen = "const",
                              p1 = 1, p2 = 1, nboots = 1000, conf = 0.90, bp = NULL,
                              ts_boot = 1)
Granger.inference.conditional(dlsp500, dlgold, dlbitcoin, ic.chosen = c("AIC","HQ","SC","FPE"),
                              max.lag = min(4, length(dlsp500) - 1), plot = T, type.chosen = "const",
                              p1 = 1, p2 = 1, nboots = 1000, conf = 0.90, bp = NULL,
                              ts_boot = 1)

# x <- dlbitcoin
Granger.inference.conditional(dlbitcoin, dlsp500, dlgold, ic.chosen = c("AIC","HQ","SC","FPE"),
                              max.lag = min(4, length(dlbitcoin) - 1), plot = T, type.chosen = "const",
                              p1 = 1, p2 = 1, nboots = 1000, conf = 0.90, bp = NULL,
                              ts_boot = 1)
Granger.inference.conditional(dlbitcoin, dlgold, dlsp500, ic.chosen = c("AIC","HQ","SC","FPE"),
                              max.lag = min(4, length(dlbitcoin) - 1), plot = T, type.chosen = "const",
                              p1 = 1, p2 = 1, nboots = 1000, conf = 0.90, bp = NULL,
                              ts_boot = 1)
# x <- dlgold
Granger.inference.conditional(dlgold, dlsp500, dlbitcoin, ic.chosen = c("AIC","HQ","SC","FPE"),
                              max.lag = min(4, length(dlgold) - 1), plot = T, type.chosen = "const",
                              p1 = 1, p2 = 1, nboots = 1000, conf = 0.90, bp = NULL,
                              ts_boot = 1)
Granger.inference.conditional(dlgold, dlbitcoin, dlsp500, ic.chosen = c("AIC","HQ","SC","FPE"),
                              max.lag = min(4, length(dlgold) - 1), plot = T, type.chosen = "const",
                              p1 = 1, p2 = 1, nboots = 1000, conf = 0.90, bp = NULL,
                              ts_boot = 1)


#==============================================================================#
                               #### Forecast ####
#==============================================================================#
# In-sample
par(mfrow=c(3,1))
par(mar = c(2, 3, 2, 2))
fv = ts(fitted(VAR1), end=c(2023,5), frequency=12)
plot(dlsp500,ylab="S&P 500",main="In-Sample Forecasting - S&P 500");lines(fv[,"dlsp500"], col="red")
plot(dlbitcoin,ylab="Bitcoin",main="In-Sample Forecasting - BTC"); lines(fv[,"dlbitcoin"], col="red")
plot(dlgold,ylab="Gold",main="In-Sample Forecasting - Gold"); lines(fv[,"dlgold"], col="red")

# Out-of-sample:

## Short horizon
forecast_VAR1 = forecast(VAR1, h=6)
plot(forecast_VAR1,main="Out-of-sample Forecasting - VAR(1)")
## Check 6 months
dataVAR_2022m9 = window(dataVAR_balanced,end=c(2022,9))
VAR1_2022m9 = VAR(dataVAR_2022m9, p=1,type = "none")
forecast_VAR1 = forecast(VAR1_2022m9, h=8)
par(mfrow=c(3,1))
par(mar = c(2, 3, 2, 2))
plot(forecast_VAR1$forecast$dlsp500,main="Check Forecasts from VAR(1) - S&P 500"); lines(dataVAR_balanced[,"dlsp500"])
plot(forecast_VAR1$forecast$dlbitcoin,main="Check Forecasts from VAR(1) -Bitcoin"); lines(dataVAR_balanced[,"dlbitcoin"])
plot(forecast_VAR1$forecast$dlgold,main="Check Forecasts from VAR(1) -Gold"); lines(dataVAR_balanced[,"dlgold"])

## long horizon
forecast_VAR1 = forecast(VAR1, h=12)
plot(forecast_VAR1,main="Out-of-sample Forecasting - VAR(1)")
## Check 6 months
dataVAR_2022m3 = window(dataVAR_balanced,end=c(2022,3))
VAR1_2022m3 = VAR(dataVAR_2022m3, p=1,type = "none")
forecast_VAR1 = forecast(VAR1_2022m3, h=14)
par(mfrow=c(3,1))
par(mar = c(2, 3, 2, 2))
plot(forecast_VAR1$forecast$dlsp500,main="Check Forecasts from VAR(1) - S&P 500"); lines(dataVAR_balanced[,"dlsp500"])
plot(forecast_VAR1$forecast$dlbitcoin,main="Check Forecasts from VAR(1) -Bitcoin"); lines(dataVAR_balanced[,"dlbitcoin"])
plot(forecast_VAR1$forecast$dlgold,main="Check Forecasts from VAR(1) -Gold"); lines(dataVAR_balanced[,"dlgold"])

#==============================================================================#
                      #### Impulse Response Function ####
#==============================================================================#

dataVAR_balanced = cbind(dlgold,dlsp500,dlbitcoin)
dataVAR_for_oirf = dataVAR_balanced[,c("dlgold","dlsp500","dlbitcoin")]
VAR1nc_for_oirf = VAR(dataVAR_for_oirf, p=1)

# Standard OIRFs
oirf = irf(VAR1nc_for_oirf, ortho="TRUE", n.ahead=5*12, boot=TRUE, ci=0.95, runs=1000)

# Cumulative OIRFs
oirf_cum = irf(VAR1nc_for_oirf, ortho="TRUE", cumulative=TRUE, n.ahead=5*12, boot=TRUE, ci=0.95, runs=1000)

# Forecast Error Variance Decomposition
fevd <- fevd(VAR1,n.ahead =20)
plot(fevd)

#==============================================================================# 
                             #### Appendix ####
#==============================================================================#
par(mfrow=c(3,3))
# Response of DLGOLD to DLGOLD
impulse = "dlgold"; response = "dlgold" 
plot.ts(cbind(oirf$irf[[impulse]][,response], 
              oirf$Lower[[impulse]][,response], 
              oirf$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"), 
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Response of DLGOLD to DLSP
impulse = "dlsp500"; response = "dlgold"
plot.ts(cbind(oirf$irf[[impulse]][,response],
              oirf$Lower[[impulse]][,response],
              oirf$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Response of DLGOLD to dlbitcoin 
impulse = "dlbitcoin"; response = "dlgold" 
plot.ts(cbind(oirf$irf[[impulse]][,response], 
              oirf$Lower[[impulse]][,response], 
              oirf$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"), 
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")


# Response of DLSP500 to DLGOLD
impulse = "dlgold"; response = "dlsp500"
plot.ts(cbind(oirf$irf[[impulse]][,response],
              oirf$Lower[[impulse]][,response],
              oirf$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")


# Response of DLSP500 to DLSP500
impulse = "dlsp500"; response = "dlsp500"
plot.ts(cbind(oirf$irf[[impulse]][,response],
              oirf$Lower[[impulse]][,response],
              oirf$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Response of DLSP500 to DLBITCOIN
impulse = "dlbitcoin"; response = "dlsp500"
plot.ts(cbind(oirf$irf[[impulse]][,response],
              oirf$Lower[[impulse]][,response],
              oirf$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Response of DLBITCOIN to DLGOLD
impulse = "dlsp500"; response = "dlgold"
plot.ts(cbind(oirf$irf[[impulse]][,response],
              oirf$Lower[[impulse]][,response],
              oirf$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Response of DLBITCOIN to DLSP500
impulse = "dlsp500"; response = "dlbitcoin"
plot.ts(cbind(oirf$irf[[impulse]][,response],
              oirf$Lower[[impulse]][,response],
              oirf$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Response of DLBITCOIN to DLBITCOIN
impulse = "dlbitcoin"; response = "dlbitcoin"
plot.ts(cbind(oirf$irf[[impulse]][,response],
              oirf$Lower[[impulse]][,response],
              oirf$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")



par(mfrow=c(3,3))
# Cumulated Response of DLGOLD to DLGOLD
impulse = "dlgold"; response = "dlgold" 
plot.ts(cbind(oirf_cum$irf[[impulse]][,response], 
              oirf_cum$Lower[[impulse]][,response], 
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"), 
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLGOLD to DLSP500
impulse = "dlgold"; response = "dlsp500"
plot.ts(cbind(oirf_cum$irf[[impulse]][,response],
              oirf_cum$Lower[[impulse]][,response],
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLGOLD to dlbitcoin 
impulse = "dlbitcoin"; response = "dlgold" 
plot.ts(cbind(oirf_cum$irf[[impulse]][,response], 
              oirf_cum$Lower[[impulse]][,response], 
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"), 
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLSP500 to DLGOLD
impulse = "dlsp500"; response = "dlgold"
plot.ts(cbind(oirf_cum$irf[[impulse]][,response],
              oirf_cum$Lower[[impulse]][,response],
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLSP500 to DLSP500
impulse = "dlsp500"; response = "dlsp500"
plot.ts(cbind(oirf_cum$irf[[impulse]][,response],
              oirf_cum$Lower[[impulse]][,response],
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLSP500 to DLBITCOIN
impulse = "dlsp500"; response = "dlbitcoin"
plot.ts(cbind(oirf_cum$irf[[impulse]][,response],
              oirf_cum$Lower[[impulse]][,response],
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLBITCOIN to DLGOLD
impulse = "dlbitcoin"; response = "dlgold"
plot.ts(cbind(oirf_cum$irf[[impulse]][,response],
              oirf_cum$Lower[[impulse]][,response],
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLBITCOIN to DLSP500
impulse = "dlbitcoin"; response = "dlsp500"
plot.ts(cbind(oirf_cum$irf[[impulse]][,response],
              oirf_cum$Lower[[impulse]][,response],
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLBITCOIN to DLBITCOIN
impulse = "dlbitcoin"; response = "dlbitcoin"
plot.ts(cbind(oirf_cum$irf[[impulse]][,response],
              oirf_cum$Lower[[impulse]][,response],
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

par(mfrow=c(3,3))

# Cumulated Response of DLGOLD to DLGOLD
impulse = "dlgold"; response = "dlgold" 
plot.ts(cbind(oirf_cum$irf[[impulse]][,response], 
              oirf_cum$Lower[[impulse]][,response], 
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"), 
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLGOLD to DLSP
impulse = "dlsp500"; response = "dlgold"
plot.ts(cbind(oirf_cum$irf[[impulse]][,response],
              oirf_cum$Lower[[impulse]][,response],
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLGOLD to dlbitcoin 
impulse = "dlbitcoin"; response = "dlgold" 
plot.ts(cbind(oirf_cum$irf[[impulse]][,response], 
              oirf_cum$Lower[[impulse]][,response], 
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"), 
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLSP500 to DLGOLD
impulse = "dlgold"; response = "dlsp500"
plot.ts(cbind(oirf_cum$irf[[impulse]][,response],
              oirf_cum$Lower[[impulse]][,response],
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLSP500 to DLSP500
impulse = "dlsp500"; response = "dlsp500"
plot.ts(cbind(oirf_cum$irf[[impulse]][,response],
              oirf_cum$Lower[[impulse]][,response],
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLSP500 to DLBITCOIN
impulse = "dlbitcoin"; response = "dlsp500"
plot.ts(cbind(oirf_cum$irf[[impulse]][,response],
              oirf_cum$Lower[[impulse]][,response],
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLBITCOIN to DLGOLD
impulse = "dlsp500"; response = "dlgold"
plot.ts(cbind(oirf_cum$irf[[impulse]][,response],
              oirf_cum$Lower[[impulse]][,response],
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLBITCOIN to DLSP500
impulse = "dlsp500"; response = "dlbitcoin"
plot.ts(cbind(oirf_cum$irf[[impulse]][,response],
              oirf_cum$Lower[[impulse]][,response],
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

# Cumulated Response of DLBITCOIN to DLBITCOIN
impulse = "dlbitcoin"; response = "dlbitcoin"
plot.ts(cbind(oirf_cum$irf[[impulse]][,response],
              oirf_cum$Lower[[impulse]][,response],
              oirf_cum$Upper[[impulse]][,response]),
        plot.type="s", col=c("black","red","red"),
        lty=c(1,2,2), xlab="", ylab="")
abline(h=0, col="red")

