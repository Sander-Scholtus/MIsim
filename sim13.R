
rm(list = ls())
gc()

library(openxlsx)
library(sampling)
library(MASS)
library(foreach)
library(doParallel)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('sim_functions.R')

#####

mu.x <- (55 + 15)/2
sigma.x <- sqrt((55 - 15)^2 / 12)
mu.y <- 5 + 0.33 * mu.x
sigma.y <- sqrt((0.33^2) * sigma.x^2 + 200)
mu.yalt <- 0.33 * mu.x


resMIa <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 2000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '5 + 0.33 * pop.x + pop.eps'
                        , model.logit.pi = '1.53 + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x)'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim13a_HD_20230305'
                        , seed = 20230305
)

#####

resMIb <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 4000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '5 + 0.33 * pop.x + pop.eps'
                        , model.logit.pi = '-0.46 + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x)'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim13b_HD_20230305'
                        , seed = 20230305
)

#####

resMIc <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 8000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '5 + 0.33 * pop.x + pop.eps'
                        , model.logit.pi = '-1.53 + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x)'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim13c_HD_20230305'
                        , seed = 20230305
)

#####

resMId <- runSimulation(S = 100000L
                        , S_at_a_time = 5000L
                        , N = 16000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '5 + 0.33 * pop.x + pop.eps'
                        , model.logit.pi = '-2.38 + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x)'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim13d_HD_20230305'
                        , seed = 20230305
)

#####

resMIe <- runSimulation(S = 100000L
                        , S_at_a_time = 5000L
                        , N = 32000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '5 + 0.33 * pop.x + pop.eps'
                        , model.logit.pi = '-3.14 + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x)'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim13e_HD_20230305'
                        , seed = 20230305
)
