#-------------------------------------------------------------------#
# This is a modification of the original script by Samuels-Crow.
# Three covariates below are removed from the original model:
# 1) VPD (min)
# 2) VPD (max)
# 3) deep soil moisture
#
# Thus, daily VPD range (VPDrng) is also removed.
# 
# In addition, shallow soil moisture is renamed as SWC
#-------------------------------------------------------------------#

### This script file will allow the user to run the SAM model to evaluate controls on NEE or ET for 
### either the US-Vcp or US-Mpj site. The user should update the file to reflect the working directory,
### input data files, and flux of interest. Items that the user should change are indicated by “XXXX”.

#Part 1: Load data, construct data list for JAGS model, generate initial values ----

setwd("C:/Users/Owner/Dropbox/Samuelscrow_2020_JGR") # set to your working directory
library(rjags) # install rjags if necessary

dataIN = read.csv("./DataFiles/Mpj_covariate.csv") #Choose the covariate file of interest. This should be a file that has all covariates for all time periods to ensure that 6-months of precipitation is supplied to the model – not just growing season precipitation
YIN    = read.csv("./DataFiles/Mpj_YVar.csv") #Choose the response variable file of interest. This file includes growing-season response (Y) variables of interest along with indices that link the Y variables back to covariates of interest.

# The covariate data file (for dataIN) should include (1) the serial date, and (2) the following covariates (as columns in the data file):
# (a) precipitation
# (b) air temperature
# (c) VPD (average)
# (d) VPD (min) <-- removed
# (e) VPD (max) <-- removed
# (f) PAR
# (g) shallow soil moisture <-- renamed as SWC
# (h) deep soil moisture <-- removed

# The response variable data file (for YIN) should include (as columns):
# (a) ET
# (b) NEE
# (c) index to link the growing season Y variables back to the appropriate row in the covariate data set

# Nparms is the number of driving variables included to calculate main effects. 

# Below, define the range of VPD in any given day (β6 in Tables S1a, S1b, S4) and assign column for Y values
#-------------------------------------------------------------------#
# do not compute the V1 (develop_mol)
# V1 = dataIN[, 9] - dataIN[, 8] ## Calculate the difference between minimum and maximum daily VPD (XX1 = column number for VPD_max, XX2 = column number for VPDmin)
Y  = YIN[, 1] # Change "XX3" to the column for the response variable of interest (ET or NEE)
#-------------------------------------------------------------------#

Nj <- 10  # what is this?

# write.csv(data.frame(j = c(1:10),
#                      ID1 = rep(1:4, times = 4:1),
#                      ID2 = c(2:5, 3:5, 4:5, 5)),
#           "./DataFiles/InteractionIND.csv",
#           row.names = FALSE)
jIND <- read.csv("./DataFiles/InteractionIND.csv") ## This file (see E.1, below) provides indices to calculate interactions between covariates
#-------------------------------------------------------------------#

#-------------------------------------------------------------------#
Nstart = 265  ## Choose the starting index. This is an index for a row in the Y data file. The value in the indexed row should be greater than 1 to accommodate calculation of antecedent values.
Nend   = nrow(YIN)  ## Choose the ending index. 
Yday = YIN[, 3] ## Choose column in YIN that provides indices linking response variables with covariates
#-------------------------------------------------------------------#

#-------------------------------------------------------------------#
# Data read into model -- covariates are centered and standardized in the data list (IMPORTANT: replace name of covariate in quotes with the appropriate column number in the dataIN file):

# Changes for the "data" variable
# Nparms 7 --> 5
# removed: Sdeep and V_rng
# renamed: Sshall --> SWC
# PAR: 7th column of dataIN --> 6th column of dataIN
data = list(Nstart = Nstart, Nend = Nend, Nlag = 7, NlagP = 9, Nparms = 5, Yday = Yday, ID1 = jIND[, 2], ID2 = jIND[, 3], 
            Y = Y,
            VPD_F = (dataIN[, 3]-mean(dataIN[, 3], na.rm = TRUE))/sd(dataIN[, 3], na.rm = TRUE),
            TA_F = (dataIN[, 2]-mean(dataIN[, 2], na.rm = TRUE))/sd(dataIN[, 2], na.rm = TRUE),
            P_F = (dataIN[, 4]-mean(dataIN[, 4], na.rm = TRUE))/sd(dataIN[, 4], na.rm = TRUE), 
            PAR = (dataIN[, 6]-mean(dataIN[, 6], na.rm = TRUE))/sd(dataIN[, 6], na.rm = TRUE),
            SWC = (dataIN[, 5]-mean(dataIN[, 5], na.rm = TRUE))/sd(dataIN[, 5], na.rm = TRUE), 
            P1 = c(6, 13, 20, 27, 55, 83, 111, 139, 167),  # indexes for NlagP
            P2 = c(0, 7, 14, 21, 28, 56, 84, 112, 140))
#-------------------------------------------------------------------#

# Initial values are estimated using a linear model.
# As in the data list (above), covariates are centered and standardized.
# Replace name of covariate in quotes with the appropriate column number in the dataIN file.

# changes for the "X1" variable
# removed: SWC-Deep and delta VPD
# renamed(only a comment in the tail): SWC-Shallow --> SWC
# PAR: 7th column of dataIN --> 6th column of dataIN
X1  = cbind((dataIN[Yday[Nstart:Nend], 3]-mean(dataIN[Yday[Nstart:Nend], 3], na.rm = TRUE))/sd(dataIN[Yday[Nstart:Nend], 3], na.rm = TRUE),  # VPDavg
            (dataIN[Yday[Nstart:Nend], 2]-mean(dataIN[Yday[Nstart:Nend], 2], na.rm = TRUE))/sd(dataIN[Yday[Nstart:Nend], 2], na.rm = TRUE),  # Tair
            (dataIN[Yday[Nstart:Nend], 4]-mean(dataIN[Yday[Nstart:Nend], 4], na.rm = TRUE))/sd(dataIN[Yday[Nstart:Nend], 4], na.rm = TRUE),  # precipitation
            (dataIN[Yday[Nstart:Nend], 6]-mean(dataIN[Yday[Nstart:Nend], 6], na.rm = TRUE))/sd(dataIN[Yday[Nstart:Nend], 6], na.rm = TRUE),  # PAR
            (dataIN[Yday[Nstart:Nend], 5]-mean(dataIN[Yday[Nstart:Nend], 5], na.rm = TRUE))/sd(dataIN[Yday[Nstart:Nend], 5], na.rm = TRUE))  # SWC
            
X1a = cbind(X1[,1]^2,  # VPDavg ^ 2
            X1[,2]^2)  # Tair ^ 2
X2  = cbind(X1[,1]*X1[,2],  # VPDavg x Tair
            X1[,1]*X1[,3],  # VPDavg x precipitation 
            X1[,1]*X1[,4],  # VPDavg x PAR 
            X1[,1]*X1[,5],  # VPDavg x SWC-shallow 
            X1[,2]*X1[,3],  # Tair x precipitation
            X1[,2]*X1[,4],  # Tair x PAR
            X1[,2]*X1[,5],  # Tair x SWC-shallow
            X1[,3]*X1[,4],  # precipitation x PAR
            X1[,3]*X1[,5],  # precipitation x SWC-shallow
            X1[,4]*X1[,5])  # PAR x SWC-shallow

# changes for the "fit" variable
# removed: X1[,6], X1[,7]
# renamed(only a comment in the tail): SWC-Shallow --> SWC
# PAR: 7th column of dataIN --> 6th column of dataIN
fit <- lm(Y[Nstart:Nend] ~ X1[,1] + X1[,2] + X1[,3] + X1[,4] + X1[,5] + X1a[,1] + X1a[,2] + X2[,1] + X2[,2]  + X2[,3]  + X2[,4]  + X2[,5]  + X2[,6]  + X2[,7]  + X2[,8]  + X2[,9]  + X2[,10])

beta0  = fit$coefficients[1]  # Intercept
# changes for the "beta1" variable
# fit$coefficients[2:8] --> fit$coefficients[2:6]
beta1  = fit$coefficients[2:6]  # linear terms

# changes for the "beta1a" variable
# fit$coefficients[9:10] --> fit$coefficients[7:8]
beta1a = fit$coefficients[7:8]  # quadratic terms


# interaction terms
# changes for the "beta2" variable
# [11:14] --> [9:12]
# [15:17] --> [13:15]
# [18:19] --> [16:17]
# [20] --> [18]
beta2 <- matrix(data = 0, nrow = 10, ncol = 1)

beta2[1:4,1]   = as.numeric(fit$coefficients[9:12])
beta2[5:7,1]   = as.numeric(fit$coefficients[13:15])
beta2[8:9,1]   = as.numeric(fit$coefficients[16:17])
beta2[10,1]    = as.numeric(fit$coefficients[18])

# changes for the "inits" variable
# (for each list) 
# remove sigSd, sigVrng, muSd, muVrng
# rename Ss as SWC
inits = list(list(beta0 = beta0, beta1 = beta1, beta1a = beta1a, beta2 = beta2, sig = 1, sigV = 1, sigT = 1, sigP = 1, sigPAR = 1, sigSWC = 1, muV=-1, muT=0.1, muP=0.3, muPAR = 0.5, muSWC = -0.1),
             list(beta0 = beta0/10, beta1 = beta1/10, beta1a = beta1a/10, beta2 = beta2/10, sig = 1/2, sigV = 1/2, sigT = 1/2, sigP = 1/2, sigPAR = 1/2, sigSWC = 1/2, muV=-1/2, muT=0.01, muP=0.03, muPAR = 0.05, muSWC = -0.01),
             list(beta0 = beta0*10, beta1 = beta1*10, beta1a = beta1a*10, beta2 = beta2*10, sig = 2, sigV = 2, sigT = 2, sigP = 2, sigPAR = 2, sigSWC = 2, muV=-2, muT=1, muP=3, muPAR = 5, muSWC = -1))

#-------------------------------------------------------------------#

#Part 2: Initialize JAGS Model ----

n.adapt = 100 # adjust this number (and n.iter) as appropriate 
n.iter = 1000
n.chains = 3

start.time <- Sys.time()
jm1.b=jags.model("./ModelCode/JAGSModel_SAM_SamuelsCrowetal_JGR_2020.R",
                 data=data,
                 n.chains = n.chains,
                 n.adapt = n.adapt,
                 inits = inits)
Sys.time() - start.time  # 7.526699 mins
#-------------------------------------------------------------------#

#Part 3: Run JAGS Model

load.module("dic")

### Choose the parameters to monitor. For this analysis, we allowed the model to converge while monitoring variables included in "zc1"
### below before monitoring dYdX (zc1dYdX), Xant (zc1X), or Y (zc1Y)

n.iter = 40000
start.time <- Sys.time()
# changes for the "zc1" variable
# remove wSd, wVrng
# rename wSs as wSWC
zc1 = coda.samples(jm1.b,
                   variable.names = c("deviance", "beta0", "beta1", "beta1a", "beta2", "wT", "wV", "wP", "wSWC", "wP.weekly", "wP.monthly", "sig"),
                   n.iter = n.iter,
                   thin = 40,
                   n.adapt=1000)
Sys.time() - start.time

### Evaluate convergence before monitoring the variables below

n.iter = 4000
start.time <- Sys.time()
# changes for the "zc1dYdX" variable
# remove dYdSd
# rename dYdSs as dYdSWC
zc1dYdX = coda.samples(jm1.b,variable.names = c("dYdVPD","dYdT","dYdSWC","dYdP"),n.iter)
Sys.time() - start.time

start.time <- Sys.time()
# changes for the "zc1X" variable
# remove Sdeep_ant and Vrngant
# rename Sshall_ant as SWC_ant
zc1X = coda.samples(jm1.b,variable.names = c("VPDant", "TAant", "PPTant", "SWC_ant", "PAR"),n.iter)
Sys.time() - start.time

n.iter = 5000

start.time <- Sys.time()
zc1Y = coda.samples(jm1.b,variable.names=c("Y.rep"),n.iter=n.iter)
Sys.time() - start.time