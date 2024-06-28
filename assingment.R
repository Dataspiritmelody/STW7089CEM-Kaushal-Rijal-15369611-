
#Libraries used in the assignment

install.packages('matlib')#library for Linear Algebra and Multivariate Statistics
library(matlib)

install.packages("ggplot2") #Library for creating time series plot
library(ggplot2)

install.packages("rsample")
library(rsample)

install.packages("readxl")
library(readxl)


#importing data 

gene_data <- (read_excel("C:\\Users\\dipak\\OneDrive\\Desktop\\aasigmrnt\\gene_data_1718251756549.xlsx"))
#gene_data$x1 

#separating time data
time <- gene_data$`Time (min)`
Time <- as.matrix(time)
colnames(Time) <- c("Time (min)")

#separating X data for task one 
genes_1 <- as.matrix(gene_data[2:6])

# x data for task two 
genes_2 <- as.matrix(gene_data[,c(2,4,5,6)])

#given y = x2 (y_data)
Y <- as.matrix(gene_data[,"x2"])
colnames(Y) <- c("y1")

#Task 1: Preliminary data analysis    
#--------------------------------------------------------------------------------
# facet wrap 
#Time series plots (of all gene expression)
 
X_genes_1 <- ts(genes_1) #ts function create time series object
plot(X_genes_1, main = "Time Series Plot: All genes", xlab = "Time (min)", ylab ="x_genes")


#Task 1: Distribution of Each Genes
density(genes_1)
#All Input Signals
par(mfrow=c(1,2)) #plots 2x2=4 plots in same display
#Density curve of all genes 
plot(density(genes_1), type = "l", xlab = "X genes", main = "Density curve of all Genes")
#Creating histogram of X Signal
hist(genes_1, freq = FALSE, xlab = "X genes ", main = "Histogram and density curve of all Genes") #plot Histogram of X
#Plot density in histogram
lines(density(genes_1), col = "black", lwd = 1) #Plot density

par(mfrow=c(2,2)) #plots 2x2=4 plots in same display

par(mfrow=c(2,3))
#Creating histogram of X1 Gene
hist(genes_1[ ,"x1"], freq = FALSE, xlab = "x1", main = "Histogram & Density Curve of x1 Gene")
#Plot density in histogram
lines(density(genes_1[,"x1"]), lwd = 1, xlab = "x1")

#Creating histogram of X2 Gene
hist(genes_1[,"x2"], freq = FALSE, xlab = "x2", main = "Histogram & density curve of x2 Gene")
#Plot density in histogram
lines(density(genes_1[,"x2"]), lwd = 1, xlab = "x2")

#Creating histogram of X3 Gene
hist(genes_1[,"x3"], freq = FALSE, xlab = "x3", main = "Histogram & density curve of x3 Gene")
#Plot density in histogram
lines(density(genes_1[,"x3"]), lwd = 1, xlab = "x3")

#Creating histogram of X4 Gene
hist(genes_1[,"x4"], freq = FALSE, xlab = "x4", main = "Histogram & density curve of x4 Gene")
#Plot density in histogram
lines(density(genes_1[,"x4"]), lwd = 1, xlab = "x4")

#Creating histogram of X5 Gene
hist(genes_1[,"x5"], freq = FALSE, xlab = "x5", main = "Histogram & density curve of x5 Gene")
#Plot density in histogram
lines(density(genes_1[,"x5"]), lwd = 1, xlab = "x5")

#Task 1: Correlation and scatter plots to examine their dependencies 
#--------------------------------------------------------------------------------
par(mfrow=c(2,2)) #plots 2x2=4 plots in same display
#Creating Plot
plot(genes_1[,"x1"], genes_1[,"x2"], pch = 1, col="blue", main = "Correlation between x1 and x2", xlab = "x1" ,ylab="x2")
#regression line
abline(lm(genes_1[,"x2"] ~ genes_1[,"x1"]), col = "red", lwd = 1)

#Creating Plot
plot(genes_1[,"x2"], genes_1[,"x3"], pch = 1, col="blue", main = "Correlation between x2 and x3", xlab = "x2" , ylab= "x3")
#regression line
abline(lm(genes_1[,"x3"] ~ genes_1[,"x2"]), col = "red", lwd = 1)

#Creating Plot
plot(genes_1[,"x3"], genes_1[,"x4"], pch = 1, col="blue", main = "Correlation between x3 and x4", xlab = "x3" , ylab = "x4")
#regression line
abline(lm(genes_1[,"x4"] ~ genes_1[,"x3"]), col = "red", lwd = 1)

#Creating Plot
plot(genes_1[,"x4"], genes_1[,"x5"], pch = 1, col="blue", main = "Correlation between x4 and x5", xlab = "x4" , ylab ="x5")
#regression line
abline(lm(genes_1[,"x5"] ~ genes_1[,"x4"]), col = "red", lwd = 1)


#Task 2:  Regression – modeling the relationship between EEG signals
#--------------------------------------------------------------------------------
onesMatrix <- matrix(1 , length(genes_2)/4,1) # Creating matrix of ones
onesMatrix

# Creating Data for model 1 from given equation
XData_model1 <- cbind(onesMatrix, (genes_2[,"x4"]), (genes_2[,"x3"])^2)
XData_model1
# Creating Data for model 2 from given equation
XData_model2 <- cbind(onesMatrix, (genes_2[,"x4"]), (genes_2[,"x3"])^2 ,genes_2[,"x5"])
XData_model2
# Creating Data for model 3 from given equation
XData_model3 <- cbind(onesMatrix, (genes_2[,"x3"]), (genes_2[,"x4"]), (genes_2[,"x5"])^3)
XData_model3
# Creating Data for model 4 from given equation
XData_model4 <- cbind(onesMatrix, (genes_2[,"x4"]), (genes_2[,"x3"])^2, (genes_2[,"x5"])^3)
XData_model4
# Creating Data for model 5 from given equation
XData_model5 <- cbind(onesMatrix, (genes_2[,"x4"]), (genes_2[,"x1"])^2, (genes_2[,"x3"])^2)
XData_model5


#Task 2.1: 
#--------------------------------------------------------------------------------
#Calculation the least square (thetaHat)
#for model1
model1_Thetahat <- solve(t(XData_model1) %*% XData_model1) %*% t(XData_model1) %*% Y
model1_Thetahat

# Model 2
model2_Thetahat <- solve(t(XData_model2) %*% XData_model2) %*% t(XData_model2) %*% Y
model2_Thetahat

# Model 3
model3_Thetahat <- solve(t(XData_model3) %*% XData_model3) %*% t(XData_model3) %*% Y
model3_Thetahat

# Model 4
model4_Thetahat <- solve(t(XData_model4) %*% XData_model4) %*% t(XData_model4) %*% Y
model4_Thetahat

# Model 5
model5_Thetahat <- solve(t(XData_model5) %*% XData_model5) %*% t(XData_model5) %*% Y
model5_Thetahat


#Task 2.2: Model Residual Error
#--------------------------------------------------------------------------------
# Calculating y-hat  y^ = x*theta
# Model 1
model1_YHat <- XData_model1 %*% model1_Thetahat
model1_YHat

# Model 2
model2_YHat <- XData_model2 %*% model2_Thetahat
model2_YHat

# Model 3
model3_YHat <- XData_model3 %*% model3_Thetahat
model3_YHat

# Model 4
model4_YHat <- XData_model4 %*% model4_Thetahat
model4_YHat

# Model 5
model5_YHat <- XData_model5 %*% model5_Thetahat
model5_YHat


# Calculating RSS
# Model 1
RSS_model1 <- sum((Y-model1_YHat)^2)
RSS_model1

# Model 2
RSS_model2 <- sum((Y-model2_YHat)^2)
RSS_model2

# Model 3
RSS_model3 <- sum((Y-model3_YHat)^2)
RSS_model3

# Model 4
RSS_model4 <- sum((Y-model4_YHat)^2)
RSS_model4

# Model 5
RSS_model5 <- sum((Y-model5_YHat)^2)
RSS_model5


#-------------------------------------------------------------------------------- 
#Task 2.3: Calculating Likelihood and Variance for all models
#--------------------------------------------------------------------------------
n <- length(Y)
n 
#Calculating length of Y
# Calculating Variance for
# Model 1
VAR_model1 <- RSS_model1/(n-1)
VAR_model1

# Model 2
VAR_model2 <- RSS_model2/(n-1)
VAR_model2

# Model 3
VAR_model3 <- RSS_model3/(n-1)
VAR_model3

# Model 4
VAR_model4 <- RSS_model4/(n-1)
VAR_model4

# Model 5
VAR_model5 <- RSS_model5/(n-1)
VAR_model5


# Calculating likelihood for
# Model 1
Likelihood_model1 <- -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model1))-(1/(2*VAR_model1))*RSS_model1
Likelihood_model1

# Model 2
Likelihood_model2 <- -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model2))-(1/(2*VAR_model2))*RSS_model2
Likelihood_model2

# Model 3
Likelihood_model3 <- -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model3))-(1/(2*VAR_model3))*RSS_model3
Likelihood_model3

# Model 4
Likelihood_model4 <- -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model4))-(1/(2*VAR_model4))*RSS_model4
Likelihood_model4

# Model 5
Likelihood_model5 <- -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model5))-(1/(2*VAR_model5))*RSS_model5
Likelihood_model5


#-------------------------------------------------------------------------------- 
#Task 2.4: Calculating AIC and BIC for all models
#--------------------------------------------------------------------------------
# Calculating AIC
# Model 1
AIC_model1 <- 2*(length(model1_Thetahat))-2*Likelihood_model1
AIC_model1

# Model 2
AIC_model2 <- 2*(length(model2_Thetahat))-2*Likelihood_model2
AIC_model2

# Model 3
AIC_model3 <- 2*(length(model3_Thetahat))-2*Likelihood_model3
AIC_model3

# Model 4
AIC_model4 <- 2*(length(model4_Thetahat))-2*Likelihood_model4
AIC_model4

# Model 5
AIC_model5 <- 2*(length(model1_Thetahat))-2*Likelihood_model5
AIC_model5

# Calculating BIC
# Model 1
BIC_model1 <- length(model1_Thetahat)*log(n)-2*Likelihood_model1
BIC_model1

# Model 2
BIC_model2 <- length(model2_Thetahat)*log(n)-2*Likelihood_model2
BIC_model2

# Model 3
BIC_model3 <- length(model3_Thetahat)*log(n)-2*Likelihood_model3
BIC_model3

# Model 4
BIC_model4 <- length(model4_Thetahat)*log(n)-2*Likelihood_model4
BIC_model4

# Model 5
BIC_model5 <- length(model5_Thetahat)*log(n)-2*Likelihood_model5
BIC_model5


#-------------------------------------------------------------------------------- 
#Task 2.5: Calculating Error for all models and Plotting Q-Q plot with Q-Q line for them
#--------------------------------------------------------------------------------

par(mfrow = c(3, 2))
# Model 1
Error_model1 <- Y - model1_YHat #Error
qqnorm(Error_model1, col = "#336600", main = "Q-Q plot of Model 1") # Plots Graph
qqline(Error_model1, col = "#e60000", lwd = 1) # Adds Q-Q line on graph

# Model 2
Error_model2 <- Y - model2_YHat #Error
qqnorm(Error_model2, col = "#336600", main = "Q-Q plot of Model 2") # Plots Graph
qqline(Error_model2, col = "#e60000", lwd = 1) # Adds Q-Q line on graph

# Model 3
Error_model3 <- Y - model3_YHat #Error
qqnorm(Error_model3, col = "#336600", main = "Q-Q plot of Model 3") # Plots Graph
qqline(Error_model3, col = "#e60000", lwd = 1) # Adds Q-Q line on graph

# Model 4
Error_model4 <- Y - model4_YHat #Error
qqnorm(Error_model4, col = "#336600", main = "Q-Q plot of Model 4") # Plots Graph
qqline(Error_model4, col = "#e60000", lwd = 1) # Adds Q-Q line on graph

# Model 5
Error_model5 <- Y - model5_YHat #Error
qqnorm(Error_model5, col = "#336600", main = "Q-Q plot of Model 5") # Plots Graph
qqline(Error_model5, col = "#e60000", lwd = 1) # Adds Q-Q line on graph



#-------------------------------------------------------------------------------- 
#Task 2.7: Splitting input and output data sets into training and testing sets in ratio of 7:3
#--------------------------------------------------------------------------------

# Splitting the data (Training Data set)
XSplit <- initial_split(data = as.data.frame(genes_2), prop = .7)
YSplit <- initial_split(data = as.data.frame(Y), prop = .7)
# Training Data set
# Y Data
Y_Training_Set <- training(YSplit) # Y Training data set
Y_Training_Data <- as.matrix(Y_Training_Set) # Y Training data
# X Data
X_Training_Set <- training(XSplit) # X Training data set
X_Training_Data <- as.matrix(X_Training_Set) # X Training data

# Testing Data set
# Y Data
Y_Testing_Set <- testing(YSplit) # Y Testing data set
Y_Testing_Data <- as.matrix(Y_Testing_Set) # Y Testing data
# X Data
X_Testing_Set <- testing(XSplit) # X Testing data set
X_Testing_Data <- as.matrix(X_Testing_Set) # X Testing data

length(X_Training_Set$x1)
length(X_Training_Set[,"x1"])
str(X_Training_Data)
# Selecting Model 5, estimating model parameters using training data
TrainingOneXMatrix <- matrix(1, length(X_Training_Set$x1), 1) # ones matrix for training set
TrainingXModel <- cbind(TrainingOneXMatrix, X_Training_Set[, "x4"], (X_Training_Set[, "x1"])^2, (X_Training_Set[, "x3"])^2) #Training Model
TrainingThetaHat <- solve(t(TrainingXModel) %*% TrainingXModel) %*% t(TrainingXModel) %*% Y_Training_Data
TrainingThetaHat



# Computing output/prediction of Model2 using testing data set
TestingYHat <- X_Testing_Data %*% TrainingThetaHat
TestingYHat
RSStesting <- sum((Y_Testing_Set - TestingYHat)^2)
RSStesting
?t.test()
mean(Y_Training_Data)
var(Y_Training_Data)
#used t-test as sample size is less than 30 i.e 4 and varriance is know (calculated)
t.test(Y_Training_Data, mu = 
         , alternative = "two.sided", conf.level = 0.95)
C_I1 <- 1.064562
C_I2 <- 1.178227
#meu <- 0.1720075
?abline()
# With 95% of confidence interval, predicting the model and plotting them with testing data and error bars
par(mfrow = c(1, 1))
TrainingDensity <- density(Y_Training_Data) # Density of training data of output signal
TrainingDensity
plot(TrainingDensity, col="#336600", lwd = 2, main="Distribution of Output  Training Data")
abline(v = C_I1, col = "#e60000", lty=2)
abline(v = C_I2, col = "#e60000", lty=2)
#abline(v = meu, col = "#1a1a1a", lty=2)

residual <- ((Y_Testing_Set - TestingYHat)) # Calculating Error
residual

# plotting Error Bars
# Calculating Standard Deviation (Sigma)
Sigma <- sqrt(VAR_model5) #Variance of model2 from task 2.3
Sigma
XData_model5 #Data model 2 from task 2.1

dataFrame <- data.frame(xAxis = XData_model5,yAxis = Y)
dataFrame
ggplot(dataFrame) +
  geom_bar( aes(x=xAxis.1, y=y1), stat="identity", fill="#336600", alpha=0.7) +
  geom_errorbar( aes(x=xAxis.1, ymin=y1-Sigma, ymax=y1+Sigma), width=0.4, colour="#e60000", alpha=0.9, linewidth=1) +
  labs(title="Error Bar (Model 5 - x1)", x="Model 5 - x1", y = "Output Data")

ggplot(dataFrame) +
  geom_bar( aes(x=xAxis.2, y=y1), stat="identity", fill="#336600", alpha=0.7) +
  geom_errorbar( aes(x=xAxis.2, ymin=y1-Sigma, ymax=y1+Sigma), width=0.4, colour="#e60000", alpha=0.9, linewidth=1) +
  labs(title="Error Bar (Model 5 - x3)", x="Model 5 - x3", y = "Output  Data")

ggplot(dataFrame) +
  geom_bar( aes(x=xAxis.3, y=y1), stat="identity", fill="#336600", alpha=0.7) +
  geom_errorbar( aes(x=xAxis.3, ymin=y1-Sigma, ymax=y1+Sigma), width=0.4, colour="#e60000", alpha=0.9, linewidth=1) +
  labs(title="Error Bar (Model 5 - x4)", x="Model 5 - x4", y = "Output Data")

ggplot(dataFrame) +
  geom_bar( aes(x=xAxis.4, y=y1), stat="identity", fill="#336600", alpha=0.7) +
  geom_errorbar( aes(x=xAxis.4, ymin=y1-Sigma, ymax=y1+Sigma), width=0.4, colour="#e60000", alpha=0.9, linewidth=1) +
  labs(title="Error Bar (Model 5 - x5)", x="Model 5 - x5", y = "Output Data")



#-------------------------------------------------------------------------------- 
#Task 3: Approximate Bayesian Computation (ABC)
#--------------------------------------------------------------------------------
array1 <- 0
array2 <- 0
f_value <- 0
s_value <- 0

# Model 5 thetahat values from Task 2.1
ThetaBias <-  1.2951518 # chosen parameter
ThetaA <- 0.8312983 # chosen parameter
ThetaB <-  0.5385828 # set constant
ThetaC <- 0.1096679 # set constant
Epsilon <- RSS_model5 * 2 ## fixing value of epsilon, RSS_model2 from task 2.2
num <- 100 # number of iteration
##Calculating Y-hat for performing rejection ABC
counter <- 0
for (i in 1:num) {
  range1 <- runif(1, -0.483065688, 0.483065688) # calculating the range
  range1
  range2 <- runif(1, -0.143578928, 0.143578928)
  range2
  NewThetahat <- matrix(c(range1, range2, ThetaB, ThetaC))
  NewYHat <- XData_model2 %*% NewThetahat ## New Y hat and model2_Thetahat from task2.1
  NewRSS <- sum((Y - NewYHat)^2)
  NewRSS
  if (NewRSS > Epsilon){ #Performing rejection ABC
    array1[i] = range1
    array2[i] = range2
    counter = counter + 1
    Fvalue = matrix(array1)
    Svalue = matrix(array2)
  }
}

# Plotting the graph
plot(Fvalue, Svalue, col = c("#ff1a1a", "#3366ff"), main = "Joint and Marginal Posterior Distribution Model 5")













