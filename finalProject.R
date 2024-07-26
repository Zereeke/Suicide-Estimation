#5586782_IB9KW0

#Uncomment to install the packages (Ctrl+Shift+C for Windows, Cmd+Shift+C for Mac)

# install.packages("tidyr")
# install.packages("ggplot2")
# install.packages("lubridate")
# install.packages("dplyr")
# install.packages("tidyverse")
# install.packages("data.table")
# install.packages("stringr")
# install.packages("tseries")
# install.packages("forecast")
# install.packages("lmtest")
# install.packages("car")
# install.packages("nortsTest")

library(tseries)
library(tidyr)
library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyverse)
library(data.table)
library(stringr)
library(forecast)
library(lmtest)
library(car)
library(nortsTest)
library(grid)

##Data Loading and Preprocessing

##Loading and Preprocessing Google Trends Data
trendsData <- read.csv("trendsData.csv",
                       skip=2,
                       na.strings="<1", header=TRUE)
names(trendsData) <- c("Date", "Depression", "Anxiety", "Stress", "Disorder", "Suicide")
trendsData$Date <- paste(trendsData$Date, "01", sep="-")
trendsData$Date <- as.Date(trendsData$Date, format = "%Y-%m-%d") #changing to date type
trendsData$Depression <- as.numeric(trendsData$Depression)      #changing to numeric type
trendsData$Anxiety <- as.numeric(trendsData$Anxiety)
trendsData$Stress <- as.numeric(trendsData$Stress)
trendsData$Disorder <- as.numeric(trendsData$Disorder)
trendsData$Suicide <- as.numeric(trendsData$Suicide)
trendsData[is.na(trendsData)] <- 0       #if na is there, replace it with 0                

##Loading and Preprocessing CDC Data
cdcData <- read.table("CDC_Data.txt", header = TRUE, sep = "\t", quote = "\"", nrows = 204)
suicideData <- subset(cdcData, select= c("Month.Code", "Deaths")) 

names(suicideData) <- c("Date", "SuicideNumbers")
suicideData$Date <- paste(suicideData$Date, "01", sep="/")
suicideData$Date <- as.Date(suicideData$Date, format = "%Y/%m/%d")  #changing to date type

##Loading and adding the Population Data
populationData <- read.csv("Population.csv",
                           na.strings="<1", header = TRUE)
suicideData$Population <- c(populationData$POPTHM) #extracting population data

##Calculating the Suicide Rate
suicideData$SuicideRate <- (suicideData$SuicideNumbers/suicideData$Population) *100000

#Merging the two data sets by Date
Final_Data <- merge(trendsData, suicideData, by.x =c("Date"), by.y=c("Date"))

#Go to the end to see the analysis steps

##Testing correlations of two datasets
correlationSteps <- function(Series, SeriesExog) {
  pearson <- cor.test(Series, SeriesExog)
  correlationPearson <- pearson$estimate #getting pearson's r
  pValuePearson <- pearson$p.value
  pearson$parameter <- pearson$parameter
  
  shapiroSeries <- shapiro.test(Series)
  wValueSeries <- shapiroSeries$statistic #getting shapiro's w
  pValueSeries <- shapiroSeries$p.value
  
  shapiroSeriesExog <- shapiro.test(SeriesExog)
  wValueSeriesExog <- shapiroSeriesExog$statistic #getting shapiro's w
  pValueSeriesExog <- shapiroSeriesExog$p.value
  
  kendall <- cor.test(Series, SeriesExog, method = "kendall")
  correlationKendall <- kendall$estimate #getting kendall's tau
  pValueKendall <- kendall$p.value
  sampleSize <- length(Series)
  
  correlationData <- data.frame (
    testType = c("Pearson Correlation (r, p, df)", "Kendall Correlation(tau, p, sample size)"),
    testStatistic = c(correlationPearson, correlationKendall),
    pValue = c(pValuePearson, pValueKendall),
    df_SampleSize = c(pearson$parameter, sampleSize)
    )
  
  shapiroData <- data.frame(
    testType = c("Base Shapiro-Wilk Test(w, p)", "Exogenous Shapiro-Wilk Test(w, p) "),
    testStatistic = c(wValueSeries, wValueSeriesExog),
    pvalue = c(pValueSeries, pValueSeriesExog)
  )
  return(list(correlationData, shapiroData))
}

#Plotting the time series and correlation graphs 
plotData <- function(Data, xAxis, yAxis, title = "Your Visualisation", type="Time Series") {
  if (type == "Time Series"){
    p <- ggplot(Data, aes_string(x = xAxis, y = yAxis)) +
      geom_line(color = "darkblue", size = 1) +
      labs(title = title,
           y = yAxis,
           x = xAxis) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=75, hjust=1))+
      scale_x_date(date_breaks = "1 year", date_labels = "%b %Y")
  }
  else if (type == "Correlation"){
    Pearson <- cor.test(Data[[xAxis]], Data[[yAxis]])
    Kendall <- cor.test(Data[[xAxis]], Data[[yAxis]], method="kendall")
    p <- ggplot(Data, aes_string(x = xAxis, y = yAxis)) +
      geom_point(color = "seagreen", size = 1) +
      labs(title = title,
           y = yAxis,
           x = xAxis) +
      annotate("text", x=80, y =5, label = paste("Pearson's Correlation Coefficient:", 
                                                 Pearson$estimate,"\nKendall's Correlation Coefficient:", 
                                                 Kendall$estimate), size=3) +
      theme_minimal()
  }
  return(p)
}


#Getting the value of d: stationarity test using Augmented Dickey-Fuller (ADF) test
arimaStationary <- function(Series, SeriesExog = NULL){ #Getting the base and exogenous series
  if (length(SeriesExog) == 0 || is.null(SeriesExog)){
    #Conduct ADF test for the base series
    stationary_test <- adf.test(Series, alternative = "stationary")
    d <- 0
    #Increment differencing until the series becomes stationary with p > 0.05
    while(stationary_test$p.value > 0.05){ 
      Series <- diff(Series)  # Difference the series
      stationary_test <- adf.test(Series, alternative = "stationary")
      d <- d + 1 #Increment d
    }
  }
  else{
    #Conduct ADF tests for both base and exogenous variables
    stationary_test_1 <- adf.test(Series, alternative = "stationary")
    stationary_test_2 <- adf.test(SeriesExog, alternative = "stationary")
    d <- 0
    #Increment differencing until both series become stationary
    while(stationary_test_1$p.value > 0.05 || stationary_test_2$p.value > 0.05){
      if(stationary_test_1$p.value > 0.05 && stationary_test_2$p.value < 0.05) {
        Series <- diff(Series) # Difference the base variable
      }
      else if(stationary_test_1$p.value < 0.05 && stationary_test_2$p.value > 0.05){
        SeriesExog <- diff(SeriesExog) # Difference the exogenous variable
      }
      else {
        Series <- diff(Series) # Difference both variables
        SeriesExog <- diff(SeriesExog)
      }
      #Re-conduct ADF tests for both variables
      stationary_test_1 <- adf.test(Series, alternative = "stationary")
      stationary_test_2 <- adf.test(SeriesExog, alternative = "stationary")
      d <- d + 1 #Increment d
    }
  }
  result <- list(Series, SeriesExog, d) #Return the differenced series and d value
  return(result)
}

#Finding p,d, and q parameters for ARIMA modeling, these parameters were compared with auto.arima() and these showed smaller errors
orderExtraction <- function(Series, SeriesExog =NULL){
  sTest <- arimaStationary(Series, SeriesExog) # Stationarity test using arimaStationary function
  
  # Initializing variables
  best_aic <- Inf
  min_rmse <- Inf
  order <- c(0, 0, 0)
  p <- 0:8
  q <- 0:8
  d <- sTest[[3]] #Getting the value of d from the stationarity test
  
  # Generating all possible combinations of p, d, and q
  pdq_combination <- expand.grid(p, d, q)
  colnames(pdq_combination) <- c("p", "d", "q")
  
  # Main loop to iterate through each combination
  if (is.null(SeriesExog) || length(SeriesExog) == 0) { # If no exogenous series provided
    for (i in 1:nrow(pdq_combination)){
      pdq <- pdq_combination[i,]
        tryCatch(
          {
          model = arima(Series, order = c(pdq$p, pdq$d, pdq$q)) # Trying to fit ARIMA model for the current combination
          fitted_values <- fitted(model)
          mse <- sqrt(mean(( Series - fitted_values)^2)) # Calculating root mean square error
          aic <- AIC(model) # Calculating AIC (Akaike Information Criterion)
          residuals <- as.numeric(abs(model$residuals))
          
          if((Box.test(residuals, lag = 12, type = "Ljung-Box"))$p.value>0.01){ # Checking for Ljung-Box test p-value
            if ((mse < min_rmse) && (aic < best_aic)) {
              min_rmse <- mse
              best_aic <- aic
              order <- pdq # Updating order if criteria met
            }
          }
        },
        warning = function(w){
          # Ignoring warnings
        },
        error = function(e){
          # Ignoring errors
        }
        )
      }
    }
  else {
    # If exogenous series provided
    for (i in 1:nrow(pdq_combination)){
      pdq <- pdq_combination[i,]
        tryCatch({
          # Trying to fit ARIMA model with exogenous variables
          model = arima(Series, order = c(pdq$p,pdq$d, pdq$q), xreg = SeriesExog)
          fitted_values <- fitted(model)
          mse <- sqrt(mean(( Series - fitted_values)^2)) # Calculating root mean square error
          aic <- AIC(model) # Calculating AIC (Akaike Information Criterion)
          residuals <- as.numeric(abs(model$residuals))
          
          # Checking for Ljung-Box test p-value
          if((Box.test(residuals, lag = 12, type = "Ljung-Box"))$p.value>0.01){
            if ((mse < min_rmse) && (aic < best_aic)){
              min_rmse <- mse
              best_aic <- aic
              order <- pdq # Updating order if criteria met
            }
          }
        },
        warning = function(w) {
          # Ignoring warnings
        },
        error = function(e){
          # Ignoring errors
        }
        )
    }
  }
  #Returning the minimum RMSE and the optimal order
  return(list(min_rmse, order))
}

#Building Arima Model
arimaModel <- function(Series, SeriesExog = NULL, order = NULL) {
  # If order parameter is not provided, extract optimal orders using orderExtraction function
  if (is.null(order)) {
    getOrder <- orderExtraction(Series)[[2]] # Extract optimal order for Series
    getOrderExog <- orderExtraction(Series, SeriesExog)[[2]] # Extract optimal order for Series with exogenous variables
  }
  else{
    getOrderExog <- order # Use provided order for Series with exogenous variables
    getOrder <- order # Use provided order for Series
  }
  
  # Fit ARIMA model based on whether exogenous variables are provided or not
  if (length(SeriesExog) == 0 || is.null(SeriesExog)){
    Model <- arima(Series, order= c(getOrder$p, getOrder$d, getOrder$q)) # Fit ARIMA model without exogenous variables
  }
  else {
    Model <- arima(Series, order=c(getOrderExog$p,getOrderExog$d, getOrderExog$q),xreg=SeriesExog) # Fit ARIMA model with exogenous variables
  }
  
  # Return the fitted ARIMA model
  return(Model)
}

#Nowcasting using the entire data set and exogenous variables
nowcastData <- function(Series, SeriesExog =NULL, order = NULL) {
  # Define start and end dates for the data
  start_date <- as.Date("2004-01-01", format="%Y-%m-%d")
  end_date <- as.Date("2020-12-01", format="%Y-%m-%d")
  
  # Create a sequence of all dates in the range
  all_dates <- seq(start_date, end_date, by = "month")
  
  # Fit an ARIMA model
  Model <- arimaModel(Series, SeriesExog, order)
  
  # Compute residuals and mean
  residuals <- abs(Model$residuals)
  mean <- mean(residuals)
  
  # Compute fitted values with mean added to address non-zero mean issue
  fitted_values <- fitted(Model) + mean 
  
  # Compute upper and lower bounds for 95% and 80% prediction intervals
  upper_95 <- fitted_values + 1.96*sqrt(Model$sigma2)
  lower_95 <- fitted_values - 1.96*sqrt(Model$sigma2)
  upper_80 <- fitted_values + 1.28*sqrt(Model$sigma2)
  lower_80 <- fitted_values - 1.28*sqrt(Model$sigma2)
  
  # Create a data frame for comparison
  comparison <- data.frame(
    Date = all_dates,
    Actual = Series,
    Nowcast_Values = as.numeric(fitted_values),
    Nowcast_95_Upper = as.numeric(upper_95),
    Nowcast_95_Lower = as.numeric(lower_95),
    Nowcast_80_Upper = as.numeric(upper_80),
    Nowcast_80_Lower = as.numeric(lower_80)
  )
  return(comparison)
}

#Compare actual and nowcasted values
plotNowcast <- function(Series, SeriesExog=NULL, order=NULL) {
  # Compute nowcasted values
  nowcasted <- nowcastData(Series, SeriesExog = SeriesExog, order = order)
  improvementResiduals <- round(as.numeric(improvArima(Series, SeriesExog = SeriesExog, order = order)[[4]]),2)
  # Extract names of series
  Series_name <- sub(".+\\$", "", deparse(substitute(Series)))
  SeriesExog_name <- sub(".+\\$", "", deparse(substitute(SeriesExog)))
  
  # Create plot based on whether exogenous variables are provided or not because of the title difference
  if (length(SeriesExog) == 0 || is.null(SeriesExog)){
    p <- ggplot(nowcasted, aes(x = Date)) +
      geom_ribbon(aes(ymin = Nowcast_95_Lower, ymax = Nowcast_95_Upper, fill = "95% Prediction Interval"), alpha = 0.7) +
      geom_ribbon(aes(ymin = Nowcast_80_Lower, ymax = Nowcast_80_Upper, fill = "80% Prediction Interval"),  alpha = 0.5) +
      geom_line(aes(y = Nowcast_Values, color = "Estimated"), size = 1) +
      geom_point(aes(y = Actual, color = "Observed"), size = 1.2) +
      labs(title = (paste("Comparison of Observed and Estimated Suicide Rates with Historical Data of", Series_name)),
           y = "Suicide Rate(per 100,000 population)",
           x = "Years")
  }
  else {
    p <- ggplot(nowcasted, aes(x = Date)) +
      geom_ribbon(aes(ymin = Nowcast_95_Lower, ymax = Nowcast_95_Upper, fill = "95% Prediction Interval"), alpha = 0.7) +
      geom_ribbon(aes(ymin = Nowcast_80_Lower, ymax = Nowcast_80_Upper, fill = "80% Prediction Interval"),  alpha = 0.5) +
      geom_line(aes(y = Nowcast_Values, color = "Estimated"), size = 1) +
      geom_point(aes(y = Actual, color = "Observed"), size = 1.2) +
      labs(title = (paste("Comparison of Observed and Estimated Suicide Rates with Prediction Intervals Modelled with\n", 
                          Series_name, "and", SeriesExog_name, "(Mean Absolute Error Reduction:", improvementResiduals,"%)")),
           y = "Suicide Rate(per 100,000 population)",
           x = "Years") 
  }
  
  # Customize plot aesthetics
  p <- p + scale_color_manual(name ="Values", values = c("Observed"="firebrick", "Estimated"="skyblue")) +
    scale_fill_manual(name ="Intervals", values = c("95% Prediction Interval"="pink", "80% Prediction Interval"="yellow")) +
    theme_minimal()+ #editing the aesthetics of the plot
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x=element_text(angle=75, hjust=1),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.spacing.x = unit(0.5, 'cm'),
          legend.box = "horizontal", 
          panel.background = element_rect(fill = "beige"))+
    scale_x_date(date_breaks = "1 year", date_labels = "%Y")
      
  return(p)
}

#Calculating Residuals
residualCheck <- function(Series, SeriesExog =NULL, order=NULL) {
  model = arimaModel(Series, SeriesExog, order = order)  
  residuals = round(as.numeric(abs(model$residuals)), 4) #getting the residuals
  return(residuals)
}

#Checking for Improvement 
improvArima <- function(Series, SeriesExog, order=NULL){
  start_date <- as.Date("2004-01-01", format="%Y-%m-%d")
  end_date <- as.Date("2020-12-01", format="%Y-%m-%d")
  all_dates <- seq(start_date, end_date, by = "month")
  
  advancedFitMAE <- 0
  advancedFitAbsErrs <- NULL
  improvement <- 0
  baseFitAbsErrs <- residualCheck(Series)
  baseFitMAE <- mean(baseFitAbsErrs)
  
  advancedFitAbsErrs <- residualCheck(Series, SeriesExog, order = order)
  advancedFitMAE <- mean(advancedFitAbsErrs)
  
  improvement <- (1 - (advancedFitMAE / baseFitMAE)) * 100
  
  dataResiduals <- data.frame(
    Date = all_dates,
    baseFitResiduals = baseFitAbsErrs,
    advancedFitResiduals = advancedFitAbsErrs
  )
  return(list(dataResiduals, baseFitMAE, advancedFitMAE, improvement))
}

#Comparing the error of residuals
errorCheck <- function(Series, SeriesExog, order=NULL){
  Series_name <- sub(".+\\$", "", deparse(substitute(Series))) #extracting the name of the series
  SeriesExog_name <- sub(".+\\$", "", deparse(substitute(SeriesExog)))
  data <- improvArima(Series, SeriesExog, order = order)
  Residuals <- data[[1]] #Getting values for the plot
  baseFitMAE <- round(data[[2]],4)
  advancedFitMAE <- round(data[[3]],4)
  improvement <- round(data[[4]],2)
  longResiduals <- melt(setDT(Residuals), id.vars = c("Date"), variable.name = "Residuals") #converting the Residual's horizontal into vertical to plot the bar graph
  p <- ggplot(longResiduals, aes(x = Date, y = value, fill = Residuals)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = (paste("Error Comaprison between Base Model and Exogenous Model of\n", Series_name, "and", SeriesExog_name)),
         x = "Years",
         y = "Absolute Residual Errors")+
    annotate("text", x=as.Date("2005-01-01"), y =650, #positioning the text
             label = paste("Base Model MAE:", baseFitMAE
                           ,"\nExogenous Model MAE:", advancedFitMAE
                           ,"\nImprovement %:", improvement), hjust=0, size=2.8)+
    theme_minimal() + #editing the aesthetics of the plot
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x=element_text(angle=75, hjust=1),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.spacing.x = unit(0.5, 'cm'),
          legend.box = "horizontal", 
          panel.background = element_rect(fill = "beige"))+
    scale_x_date(date_breaks = "1 year", date_labels = "%Y")
  return(p)
}

#Diagnostics:

#Checking Heteroscedasticity and autocorrelation
diagnosticsCheck <- function(residuals, lag =12, alpha = 0.05) {
  arch_test <- arch.test(residuals, arch= "Lm", lag = lag) #ARCH LM Test
  arch_test_statistic <- arch_test$statistic
  arch_testp <- arch_test$p.value
  
  lb_test <- Box.test(residuals, lag = lag, type = "Ljung-Box") #Ljung Box Test
  lb_test_statistic <- lb_test$statistic
  lb_testp <- lb_test$p.value
  
  
  diagnosticData <- data.frame(
    testType = c("ARCH (Lm) Test", "Jlung Box Test"),
    testStatistic = c(arch_test_statistic, lb_test_statistic),
    pValue = c(arch_testp, lb_testp)
  )
  return(diagnosticData)
}

##Analysis Starts here : Copy the code step by step into the console, or run it at once to see the entire output

#Step 1: Basic evaluation
# Checking if the suicide rate has normal distribution or not

print("Step1: Plotting the histogram of base time series, summary of the dataset, and zscore")
print(ggplot(data=Final_Data, aes(x=Final_Data$SuicideRate)) + geom_histogram())
print(summary(Final_Data$SuicideRate))

#Checking outliers using zscore
print(head(scale(Final_Data$SuicideRate)))


# Step 2: Correlation Tests of all search terms
print("Step2: Pearson's correlation, shapiro tests, and Kendall's correlation of base time series and exogenous time series")
correlationData <- list(Depression = correlationSteps(Final_Data$SuicideRate, Final_Data$Depression),
                        Anxiety = correlationSteps(Final_Data$SuicideRate, Final_Data$Anxiety),
                        Stress = correlationSteps(Final_Data$SuicideRate, Final_Data$Stress),
                        Disorder = correlationSteps(Final_Data$SuicideRate, Final_Data$Disorder),
                        Suicide = correlationSteps(Final_Data$SuicideRate, Final_Data$Suicide))
print(correlationData)

#Step 3: Getting order of base model and exogenous models
print("Step3: Extracting order of base model and exogenous models")
orderData <- list(historicalSuicideRate = orderExtraction(Final_Data$SuicideRate)[[2]], #getting the 2nd item(order) of the returned list as 1st contains the RMSE
                  suicideRateXDepression = orderExtraction(Final_Data$SuicideRate, Final_Data$Depression)[[2]],
                  suicideRateXAnxiety = orderExtraction(Final_Data$SuicideRate, Final_Data$Anxiety)[[2]],
                  suicideRateXStress = orderExtraction(Final_Data$SuicideRate, Final_Data$Stress)[[2]],
                  suicideRateXDisorder = orderExtraction(Final_Data$SuicideRate, Final_Data$Disorder)[[2]],
                  suicideRateXSuicide = orderExtraction(Final_Data$SuicideRate, Final_Data$Suicide)[[2]])

model <- names(orderData)
p <- numeric(length(orderData))
d <- numeric(length(orderData))
q <- numeric(length(orderData))
i = 1
for(order in orderData) {
  p[i] <- order$p
  d[i] <- order$d
  q[i] <- order$q
  i <- i+1
}

arimaOrder <- data.frame(model = model, p = p, d = d, q = q)
print(arimaOrder)

print("Step4: Internally arima models have been created")
## Step4: Creating models, this step can be skipped as models are created internally in the functions used
arimaModels <- list(historicalSuicideRate = arimaModel(Final_Data$SuicideRate, order = orderData[[1]]),
                      suicideRateXDepression = arimaModel(Final_Data$SuicideRate, Final_Data$Depression, order = orderData[[2]]),
                      suicideRateXAnxiety = arimaModel(Final_Data$SuicideRate, Final_Data$Anxiety, order = orderData[[3]]),
                      suicideRateXStress = arimaModel(Final_Data$SuicideRate, Final_Data$Stress, order = orderData[[4]]),
                      suicideRateXDisorder = arimaModel(Final_Data$SuicideRate, Final_Data$Disorder, order = orderData[[5]]),
                      suicideRateXSuicide = arimaModel(Final_Data$SuicideRate, Final_Data$Suicide, order = orderData[[6]]))

print("Step5: Estimating suicide rates using all the models and the comparison with actual values")
##Step 5: Estimating Suicide Rates using all the models
nowcastedData <- list(historicalSuicideRate = nowcastData(Final_Data$SuicideRate, order = orderData[[1]]),
                        suicideRateXDepression = nowcastData(Final_Data$SuicideRate, Final_Data$Depression, order = orderData[[2]]),
                        suicideRateXAnxiety = nowcastData(Final_Data$SuicideRate, Final_Data$Anxiety, order = orderData[[3]]),
                        suicideRateXStress = nowcastData(Final_Data$SuicideRate, Final_Data$Stress, order = orderData[[4]]),
                        suicideRateXDisorder = nowcastData(Final_Data$SuicideRate, Final_Data$Disorder, order = orderData[[5]]),
                        suicideRateXSuicide = nowcastData(Final_Data$SuicideRate, Final_Data$Suicide, order = orderData[[6]]))


print(nowcastedData)



#While this is highly unlikely, but if the graphs are not plotted with this error: Error in UseMethod("depth") : no applicable method for 'depth' applied to an object of class "NULL"
##run this in the console: dev.off()

print("Step6: Plotting the comparison between observed and nowcasted values for each model: look at the plots")
##Step 6: Plotting the comparison between observed and nowcasted values for each model
plotNowcastedData <- list(historicalSuicideRate = plotNowcast(Final_Data$SuicideRate, order = orderData[[1]]),
                      suicideRateXDepression = plotNowcast(Final_Data$SuicideRate, Final_Data$Depression, order = orderData[[2]]),
                      suicideRateXAnxiety = plotNowcast(Final_Data$SuicideRate, Final_Data$Anxiety, order = orderData[[3]]),
                      suicideRateXStress = plotNowcast(Final_Data$SuicideRate, Final_Data$Stress, order = orderData[[4]]),
                      suicideRateXDisorder = plotNowcast(Final_Data$SuicideRate, Final_Data$Disorder, order = orderData[[5]]),
                      suicideRateXSuicide = plotNowcast(Final_Data$SuicideRate, Final_Data$Suicide, order = orderData[[6]]))

#Plotting combined graphs
edittedxAxis <- c(1,2,4)

for (i in edittedxAxis) { #editing x axis
  plotNowcastedData[[i]] <- plotNowcastedData[[i]] +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(size=10))
}

edittedyAxis <- c(1,3,5) #editing y axis
for (i in edittedyAxis) {
  plotNowcastedData[[i]] <- plotNowcastedData[[i]] +
    theme(axis.title.y = element_blank(),
          plot.title = element_text(size=10))
}

grid.newpage()
grid.draw(rbind(ggplotGrob(plotNowcastedData[[1]]),
                ggplotGrob(plotNowcastedData[[2]]),
                ggplotGrob(plotNowcastedData[[3]])))


grid.newpage()
grid.draw(rbind(ggplotGrob(plotNowcastedData[[1]]),
                ggplotGrob(plotNowcastedData[[4]]),
                ggplotGrob(plotNowcastedData[[5]]))) #ignoring suicide search term as it only showed 0.55% improvement in MAE


# Step 7: Comparing the MAE and residual errors
print("Step7a: Calculating residual errors for each model for the next step")
# Calculate residual errors for each model
residualErrors <- list(historicalSuicideRate = residualCheck(Final_Data$SuicideRate, order = orderData[[1]]),
  suicideRateXDepression = residualCheck(Final_Data$SuicideRate, Final_Data$Depression, order = orderData[[2]]),
  suicideRateXAnxiety = residualCheck(Final_Data$SuicideRate, Final_Data$Anxiety, order = orderData[[3]]),
  suicideRateXStress = residualCheck(Final_Data$SuicideRate, Final_Data$Stress, order = orderData[[4]]),
  suicideRateXDisorder = residualCheck(Final_Data$SuicideRate, Final_Data$Disorder, order = orderData[[5]]),
  suicideRateXSuicide = residualCheck(Final_Data$SuicideRate, Final_Data$Suicide, order = orderData[[6]]))


#Plotting the normal distribution and MAE improvement for each
print("Step7b: Plotting the normal distribution and MAE improvement for each model: look at the plots")
plotImprovementMAE <- list(suicideRateXDepression = errorCheck(Final_Data$SuicideRate, Final_Data$Depression, order = orderData[[2]]),
                      suicideRateXAnxiety = errorCheck(Final_Data$SuicideRate, Final_Data$Anxiety, order = orderData[[3]]),
                      suicideRateXStress = errorCheck(Final_Data$SuicideRate, Final_Data$Stress, order = orderData[[4]]),
                      suicideRateXDisorder = errorCheck(Final_Data$SuicideRate, Final_Data$Disorder, order = orderData[[5]]),
                      suicideRateXSuicide = errorCheck(Final_Data$SuicideRate, Final_Data$Suicide, order = orderData[[6]]))

#Plotting combined graphs
edittedxAxis <- c(1,3,4)

for (i in edittedxAxis) { #Editing the x axis
  plotImprovementMAE[[i]] <- plotImprovementMAE[[i]] +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(size=10))
}

edittedyAxis <- c(2,3,5) #Editing the y axis
for (i in edittedyAxis) {
  plotImprovementMAE[[i]] <- plotImprovementMAE[[i]] +
    theme(axis.title.y = element_blank(),
          plot.title = element_text(size=10))
}

grid.newpage()
grid.draw(rbind(ggplotGrob(plotImprovementMAE[[1]]),
                ggplotGrob(plotImprovementMAE[[2]])))


grid.newpage()
grid.draw(rbind(ggplotGrob(plotImprovementMAE[[3]]),
                ggplotGrob(plotImprovementMAE[[4]]),
                ggplotGrob(plotImprovementMAE[[5]])))


# Step 8: Testing Assumptions: heteroskedasticity and autocorrelation in the residual errors
print("Step8: Model Diagnostics: Checking heteroskedasticity and autocorrelation in the residual errors")
diagnosticsCheck <- list(historicalSuicideRate = diagnosticsCheck(residualErrors[[1]]),
  suicideRateXDepression = diagnosticsCheck(residualErrors[[2]]),
  suicideRateXAnxiety = diagnosticsCheck(residualErrors[[3]]),
  suicideRateXStress = diagnosticsCheck(residualErrors[[4]]),
  suicideRateXDisorder = diagnosticsCheck(residualErrors[[5]]),
  suicideRateXSuicide = diagnosticsCheck(residualErrors[[6]]))
print(diagnosticsCheck)

print("Ignore the warnings in this context as even if the pvalue is smaller than what's given in the output, the analysis procedure won't change for the arima model")


