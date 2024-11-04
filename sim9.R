
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
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '5 + 0.33 * pop.x + pop.eps'
                        , model.logit.pi = '1.53 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim9a_HD_20230305'
                        , seed = 20230305
)

#####

resMIb <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 4000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '5 + 0.33 * pop.x + pop.eps'
                        , model.logit.pi = '-0.46 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim9b_HD_20230305'
                        , seed = 20230305
)

#####

resMIc <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 8000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '5 + 0.33 * pop.x + pop.eps'
                        , model.logit.pi = '-1.53 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim9c_HD_20230305'
                        , seed = 20230305
)

#####

resMId <- runSimulation(S = 100000L
                        , S_at_a_time = 5000L
                        , N = 16000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '5 + 0.33 * pop.x + pop.eps'
                        , model.logit.pi = '-2.38 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim9d_HD_20230305'
                        , seed = 20230305
)

#####

resMIe <- runSimulation(S = 100000L
                        , S_at_a_time = 5000L
                        , N = 32000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '5 + 0.33 * pop.x + pop.eps'
                        , model.logit.pi = '-3.14 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim9e_HD_20230305'
                        , seed = 20230305
)

#####

resMIf <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 2000L
                        , cfs = 1L
                        , use.intercept = FALSE
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '0.33 * pop.x + pop.eps'
                        , model.logit.pi = '1.53 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.yalt)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim9f_HD_no_icpt_20230305'
                        , seed = 20230305
)

#####

resMIg <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 4000L
                        , cfs = 1L
                        , use.intercept = FALSE
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '0.33 * pop.x + pop.eps'
                        , model.logit.pi = '-0.46 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.yalt)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim9g_HD_no_icpt_20230305'
                        , seed = 20230305
)

#####

resMIh <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 8000L
                        , cfs = 1L
                        , use.intercept = FALSE
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '0.33 * pop.x + pop.eps'
                        , model.logit.pi = '-1.53 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.yalt)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim9h_HD_no_icpt_20230305'
                        , seed = 20230305
)

#####

resMIi <- runSimulation(S = 100000L
                        , S_at_a_time = 5000L
                        , N = 16000L
                        , cfs = 1L
                        , use.intercept = FALSE
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '0.33 * pop.x + pop.eps'
                        , model.logit.pi = '-2.38 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.yalt)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim9i_HD_no_icpt_20230305'
                        , seed = 20230305
)

#####

resMIj <- runSimulation(S = 100000L
                        , S_at_a_time = 5000L
                        , N = 32000L
                        , cfs = 1L
                        , use.intercept = FALSE
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '0.33 * pop.x + pop.eps'
                        , model.logit.pi = '-3.14 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.yalt)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim9j_HD_no_icpt_20230305'
                        , seed = 20230305
)


#####

resMIk <- runSimulation(S = 100000L
                        , S_at_a_time = 5000L
                        , N = 80000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                        , model.y = '5 + 0.33 * pop.x + pop.eps'
                        , model.logit.pi = '-4.10 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                        , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim9k_HD_20230305'
                        , seed = 20230305
)


##### Method of Kim & Yang

resMIa.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 2000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                           , model.y = '5 + 0.33 * pop.x + pop.eps'
                           , model.logit.pi = '1.53 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                           , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim9a_KY_20230305'
                           , seed = 20230305
)

#####

resMIb.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 4000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                           , model.y = '5 + 0.33 * pop.x + pop.eps'
                           , model.logit.pi = '-0.46 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                           , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim9b_KY_20230305'
                           , seed = 20230305
)

#####

resMIc.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 8000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                           , model.y = '5 + 0.33 * pop.x + pop.eps'
                           , model.logit.pi = '-1.53 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                           , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim9c_KY_20230305'
                           , seed = 20230305
)

#####

resMId.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 5000L
                           , N = 16000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                           , model.y = '5 + 0.33 * pop.x + pop.eps'
                           , model.logit.pi = '-2.38 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                           , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim9d_KY_20230305'
                           , seed = 20230305
)

#####

resMIe.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 5000L
                           , N = 32000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                           , model.y = '5 + 0.33 * pop.x + pop.eps'
                           , model.logit.pi = '-3.14 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                           , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim9e_KY_20230305'
                           , seed = 20230305
)

#####

resMIf.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 2000L
                           , cfs = 1L
                           , use.intercept = FALSE
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                           , model.y = '0.33 * pop.x + pop.eps'
                           , model.logit.pi = '1.53 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.yalt)/sigma.y)'
                           , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , writeExcel = TRUE
                           , nameExcel = 'sim9f_KY_no_icpt_20230305'
                           , seed = 20230305
)

#####

resMIg.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 4000L
                           , cfs = 1L
                           , use.intercept = FALSE
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                           , model.y = '0.33 * pop.x + pop.eps'
                           , model.logit.pi = '-0.46 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.yalt)/sigma.y)'
                           , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , writeExcel = TRUE
                           , nameExcel = 'sim9g_KY_no_icpt_20230305'
                           , seed = 20230305
)

#####

resMIh.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 8000L
                           , cfs = 1L
                           , use.intercept = FALSE
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                           , model.y = '0.33 * pop.x + pop.eps'
                           , model.logit.pi = '-1.53 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.yalt)/sigma.y)'
                           , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , writeExcel = TRUE
                           , nameExcel = 'sim9h_KY_no_icpt_20230305'
                           , seed = 20230305
)

#####

resMIi.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 5000L
                           , N = 16000L
                           , cfs = 1L
                           , use.intercept = FALSE
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                           , model.y = '0.33 * pop.x + pop.eps'
                           , model.logit.pi = '-2.38 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.yalt)/sigma.y)'
                           , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , writeExcel = TRUE
                           , nameExcel = 'sim9i_KY_no_icpt_20230305'
                           , seed = 20230305
)

#####

resMIj.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 5000L
                           , N = 32000L
                           , cfs = 1L
                           , use.intercept = FALSE
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                           , model.y = '0.33 * pop.x + pop.eps'
                           , model.logit.pi = '-3.14 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.yalt)/sigma.y)'
                           , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , writeExcel = TRUE
                           , nameExcel = 'sim9j_KY_no_icpt_20230305'
                           , seed = 20230305
)

#####

resMIk.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 5000L
                           , N = 80000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(200))'
                           , model.y = '5 + 0.33 * pop.x + pop.eps'
                           , model.logit.pi = '-4.10 - 0.33 * pop.u + 0.6 * ((pop.y[ ,s] - mu.y)/sigma.y)'
                           , model.logit.phi = '0.5 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , writeExcel = TRUE
                           , nameExcel = 'sim9k_KY_20230305'
                           , seed = 20230305
)


