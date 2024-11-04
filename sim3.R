
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
mu.y <- 5 + mu.x
sigma.y <- sqrt(sigma.x^2 + 126)
mu.y.alt <- 5 + mu.x/2
sigma.y.alt <- sqrt((sigma.x / 2)^2 + 126)


resMIa <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 2000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                        , model.y = '5 + pop.x + pop.eps'
                        , model.logit.pi = '1.77 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                        , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim3a_HD_20230305'
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
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                        , model.y = '5 + pop.x + pop.eps'
                        , model.logit.pi = '-0.20 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                        , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim3b_HD_20230305'
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
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                        , model.y = '5 + pop.x + pop.eps'
                        , model.logit.pi = '-1.27 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                        , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim3c_HD_20230305'
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
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                        , model.y = '5 + pop.x + pop.eps'
                        , model.logit.pi = '-2.11 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                        , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim3d_HD_20230305'
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
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                        , model.y = '5 + pop.x + pop.eps'
                        , model.logit.pi = '-2.87 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                        , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim3e_HD_20230305'
                        , seed = 20230305
)

#####

resMIf <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 2000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                        , model.y = '5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '1.77 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y.alt)/sigma.y.alt) * sqrt(1.25) - 0.5)'
                        , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim3f_HD_20230305'
                        , seed = 20230305
)

#####

resMIg <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 2000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 15, max = 55)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                        , model.y = '5 + pop.x + pop.eps'
                        , model.logit.pi = '1.77 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                        , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) - 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim3g_HD_20230305'
                        , seed = 20230305
)



##### Using draws from a normal distribution instead of hot deck imputations

resMIa.ND <- runSimulation(S = 100000L
                           , S_at_a_time = 10000L
                           , N = 2000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                           , model.y = '5 + pop.x + pop.eps'
                           , model.logit.pi = '1.77 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                           , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'normal'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = TRUE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim3a_ND_20230305'
                           , seed = 20230305
)

#####

resMIb.ND <- runSimulation(S = 100000L
                           , S_at_a_time = 10000L
                           , N = 4000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                           , model.y = '5 + pop.x + pop.eps'
                           , model.logit.pi = '-0.20 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                           , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'normal'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = TRUE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim3b_ND_20230305'
                           , seed = 20230305
)

#####

resMIc.ND <- runSimulation(S = 100000L
                           , S_at_a_time = 10000L
                           , N = 8000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                           , model.y = '5 + pop.x + pop.eps'
                           , model.logit.pi = '-1.27 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                           , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'normal'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = TRUE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim3c_ND_20230305'
                           , seed = 20230305
)

#####

resMId.ND <- runSimulation(S = 100000L
                           , S_at_a_time = 5000L
                           , N = 16000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                           , model.y = '5 + pop.x + pop.eps'
                           , model.logit.pi = '-2.11 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                           , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'normal'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = TRUE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim3d_ND_20230305'
                           , seed = 20230305
)

#####

resMIe.ND <- runSimulation(S = 100000L
                           , S_at_a_time = 5000L
                           , N = 32000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                           , model.y = '5 + pop.x + pop.eps'
                           , model.logit.pi = '-2.87 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                           , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'normal'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = TRUE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim3e_ND_20230305'
                           , seed = 20230305
)

#####

resMIf.ND <- runSimulation(S = 100000L
                           , S_at_a_time = 10000L
                           , N = 2000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                           , model.y = '5 + 0.5 * pop.x + pop.eps'
                           , model.logit.pi = '1.77 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y.alt)/sigma.y.alt) * sqrt(1.25) - 0.5)'
                           , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'normal'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = TRUE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim3f_ND_20230305'
                           , seed = 20230305
)

#####

resMIg.ND <- runSimulation(S = 100000L
                           , S_at_a_time = 10000L
                           , N = 2000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                           , model.y = '5 + pop.x + pop.eps'
                           , model.logit.pi = '1.77 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                           , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) - 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'normal'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = TRUE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim3g_ND_20230305'
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
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                           , model.y = '5 + pop.x + pop.eps'
                           , model.logit.pi = '1.77 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                           , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim3a_KY_20230305'
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
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                           , model.y = '5 + pop.x + pop.eps'
                           , model.logit.pi = '-0.20 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                           , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim3b_KY_20230305'
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
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                           , model.y = '5 + pop.x + pop.eps'
                           , model.logit.pi = '-1.27 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                           , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim3c_KY_20230305'
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
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                           , model.y = '5 + pop.x + pop.eps'
                           , model.logit.pi = '-2.11 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                           , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim3d_KY_20230305'
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
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                           , model.y = '5 + pop.x + pop.eps'
                           , model.logit.pi = '-2.87 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                           , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim3e_KY_20230305'
                           , seed = 20230305
)

#####

resMIf.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 2000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                           , model.y = '5 + 0.5 * pop.x + pop.eps'
                           , model.logit.pi = '1.77 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y.alt)/sigma.y.alt) * sqrt(1.25) - 0.5)'
                           , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim3f_KY_20230305'
                           , seed = 20230305
)

#####

resMIg.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 2000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 15, max = 55)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = sqrt(126))'
                           , model.y = '5 + pop.x + pop.eps'
                           , model.logit.pi = '1.77 - 0.33 * pop.u + 0.5 * (((pop.y[ ,s] - mu.y)/sigma.y) * sqrt(1.25) - 0.5)'
                           , model.logit.phi = '1 + 0.5 * ((pop.x - mu.x)/sigma.x) - 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim3g_KY_20230305'
                           , seed = 20230305
)


