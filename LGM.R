getwd()
setwd("C:/Users/adetu/Desktop")

library(foreign) 
library(dplyr)
library(tidyr)
library(naniar)
library(ggplot2)
library(haven)
library(lavaan)
library(semPlot)
library(knitr)

ECLS_data <- read_sav("ECLS1.sav", encoding = NULL, user_na = FALSE, col_select = NULL, skip = 0, n_max = Inf, .name_repair = "universal")
ECLS_data %>% replace_with_na(replace = list(P1CURMAR = c(-1,-2,-5,-9,-8,-7), 
                                 P4CURMAR = c(-1,-2,-5,-9,-8,-7),
                                 P7CURMAR = c(-1,-2,-5,-9,-8,-7),
                                 P9CURMAR = c(-1,-2,-5,-9,-8,-7)))

ECLS_data <- na.omit(ECLS_data)

#Latent Factors

lf <- c("X1MTHETK5", "X4MTHETK5", "X7MTHETK5", "X9MTHETK5")

summary(ECLS_data[lf])
lf.cov <- cov(ECLS_data[lf])
lf.cov 
lf.mean <- colMeans(ECLS_data[lf])
lf.mean

#Latent Model Design

time <- c("X1MTHETK5", "X4MTHETK5", "X7MTHETK5", "X9MTHETK5", "P1EXPECT", "X12SESL", "X_CHSEX_R", "X2DISABL", "X1HPARNT", "P4HLPHWK", "P7HWKCMP", "T2PARIN", "X7PWKMEM", "X7TWKMEM")

ECLS.cov <- cov(ECLS_data[time])
ECLS.cov
ECLS.mean <- colMeans(ECLS_data[time])
ECLS.mean
names(ECLS.mean) <- colnames(ECLS.cov) <- rownames(ECLS.cov) <- c("X1MTHETK5", "X4MTHETK5", "X7MTHETK5", "X9MTHETK5", "P1EXPECT", "X12SESL", "X_CHSEX_R", "X2DISABL", "X1HPARNT", "P4HLPHWK", "P7HWKCMP", "T2PARIN", "X7PWKMEM", "X7TWKMEM")

# Basic latent curve model specification

library(lavaan)

# Hypothesis Question 1

ECLS.model1 <- '
                #regressions
                X1MTHETK5 + X4MTHETK5 + X7MTHETK5 + X9MTHETK5 ~  P1EXPECT + X12SESL + X_CHSEX_R + X2DISABL
                
                # intercept
                i =~ 1*X1MTHETK5 + 1*X4MTHETK5 + 1*X7MTHETK5 + 1*X9MTHETK5
                i~~0*i
                
                # slope
                s =~ 0*X1MTHETK5 + 1*X4MTHETK5 + 2*X7MTHETK5 + 3*X9MTHETK5
                s~0*1
                s~~0*i
                
                # residual variances
                X1MTHETK5~~r*X1MTHETK5
                X4MTHETK5~~r*X4MTHETK5
                X7MTHETK5~~r*X7MTHETK5
                X9MTHETK5~~r*X9MTHETK5
                
                '

ECLS.fit1 <- growth(ECLS.model1, sample.cov=ECLS.cov, sample.mean=ECLS.mean, sample.nobs = 2921)

semPaths(ECLS.fit1)

summary(ECLS.fit1, standardized = TRUE)

AIC(ECLS.fit1)

fitMeasures(ECLS.fit1, c("cfi","rmsea","srmr"))

#Hypothesis Question 2

ECLS.model2 <- '
                #regressions
                X1MTHETK5 + X4MTHETK5 + X7MTHETK5 + X9MTHETK5 ~  X_CHSEX_R 
                
                # intercept
                i =~ 1*X1MTHETK5 + 1*X4MTHETK5 + 1*X7MTHETK5 + 1*X9MTHETK5
                i~~0*i
                
                # slope
                s =~ 0*X1MTHETK5 + 1*X4MTHETK5 + 2*X7MTHETK5 + 3*X9MTHETK5
                s~0*1
                s~~0*i
                
                # residual variances
                X1MTHETK5~~r*X1MTHETK5
                X4MTHETK5~~r*X4MTHETK5
                X7MTHETK5~~r*X7MTHETK5
                X9MTHETK5~~r*X9MTHETK5
                '

ECLS.fit2 <- growth(ECLS.model2, sample.cov=ECLS.cov, sample.mean=ECLS.mean, sample.nobs=2921)

semPaths(ECLS.fit2)

summary(ECLS.fit2, standardized = TRUE)

AIC(ECLS.fit2)

fitMeasures(ECLS.fit2, c("cfi","rmsea","srmr"))

#Hypothesis Question 3

ECLS.model3 <- '
                #regressions
                X1MTHETK5 + X4MTHETK5 + X7MTHETK5 + X9MTHETK5 ~ X1HPARNT + P4HLPHWK + P7HWKCMP + T2PARIN + P1EXPECT
                
                # intercept
                i =~ 1*X1MTHETK5 + 1*X4MTHETK5 + 1*X7MTHETK5 + 1*X9MTHETK5
                i~~0*i
                
                # slope
                s =~ 0*X1MTHETK5 + 1*X4MTHETK5 + 2*X7MTHETK5 + 3*X9MTHETK5
                s~0*1
                s~~0*i
                
                # residual variances
                X1MTHETK5~~r*X1MTHETK5
                X4MTHETK5~~r*X4MTHETK5
                X7MTHETK5~~r*X7MTHETK5
                X9MTHETK5~~r*X9MTHETK5
                
                '
ECLS.fit3 <- growth(ECLS.model3, sample.cov=ECLS.cov, sample.mean=ECLS.mean, sample.nobs=2921)

semPaths(ECLS.fit3)

summary(ECLS.fit3, standardized = TRUE)

AIC(ECLS.fit3)

fitMeasures(ECLS.fit3, c("cfi","rmsea","srmr"))

# Hypothesis Question 4

ECLS.model4 <- '
                #regressions
                X1MTHETK5 + X4MTHETK5 + X7MTHETK5 + X9MTHETK5 ~  X7PWKMEM + X7TWKMEM
                
                # intercept
                i =~ 1*X1MTHETK5 + 1*X4MTHETK5 + 1*X7MTHETK5 + 1*X9MTHETK5
                i~~0*i
                
                # slope
                s =~ 0*X1MTHETK5 + 1*X4MTHETK5 + 2*X7MTHETK5 + 3*X9MTHETK5
                s~0*1
                s~~0*i
                
                # residual variances
                X1MTHETK5~~r*X1MTHETK5
                X4MTHETK5~~r*X4MTHETK5
                X7MTHETK5~~r*X7MTHETK5
                X9MTHETK5~~r*X9MTHETK5
                
                '
ECLS.fit4 <- growth(ECLS.model4, sample.cov=ECLS.cov, sample.mean=ECLS.mean, sample.nobs=2921)

semPaths(ECLS.fit4)

summary(ECLS.fit4, standardized = TRUE)

AIC(ECLS.fit4)

fitMeasures(ECLS.fit4, c("cfi","rmsea","srmr"))

#END