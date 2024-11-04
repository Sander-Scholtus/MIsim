
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
                        , model.logit.pi = '-4 - 0.5 * pop.x'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = FALSE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim11a_HD_20230305'
                        , seed = 20230305
)


resMIb <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 1800L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '0.25 - 0.5 * pop.x'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = FALSE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim11b_HD_20230305'
                        , seed = 20230305
)


resMIc <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 1600L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '0.5 - 0.5 * pop.x'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = FALSE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim11c_HD_20230305'
                        , seed = 20230305
)

resMId <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 1400L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '1 - 0.5 * pop.x'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = FALSE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim11d_HD_20230305'
                        , seed = 20230305
)
