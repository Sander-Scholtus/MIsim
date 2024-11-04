
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

# original scenario from Kim & Yang's paper
resMIa <- runSimulation(S = 100000L
                        , S_at_a_time = 5000L
                        , N = 50000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '-3.66 - 0.33 * pop.u + 0.1 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim10a_HD_20230305'
                        , seed = 20230305
)


resMIb <- runSimulation(S = 100000L
                        , S_at_a_time = 5000L
                        , N = 50000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '-3.66 - 0.33 * pop.u + 0.5 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim10b_HD_20230305'
                        , seed = 20230305
)



##### Method of Kim & Yang

resMIa.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 5000L
                           , N = 50000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                           , model.logit.pi = '-3.66 - 0.33 * pop.u + 0.1 * pop.y[ ,s]'
                           , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim10a_KY_20230305'
                           , seed = 20230305
)

#####

resMIb.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 5000L
                           , N = 50000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                           , model.logit.pi = '-3.66 - 0.33 * pop.u + 0.5 * pop.y[ ,s]'
                           , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim10b_KY_20230305'
                           , seed = 20230305
)
