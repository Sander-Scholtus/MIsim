
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

resMIa1 <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 10000L
                        , cfs = 2L
                        , use.intercept = TRUE
                        , m = 100L
                        , model.x = 'runif(n = N, min = 2, max = 5)'
                        , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                        , model.y = '(4/5) * pop.x + pop.eps'
                        , model.pi = '200 * pop.x / (3.5 * N)'
                        , model.logit.phi = '8 - 2 * pop.x'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = FALSE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim8a_HD_20230305'
                        , seed = 20230305
)

resMIa0 <- runSimulation(S = 100000L
                         , S_at_a_time = 10000L
                         , N = 10000L
                         , cfs = 1L
                         , use.intercept = FALSE
                         , m = 100L
                         , model.x = 'runif(n = N, min = 2, max = 5)'
                         , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                         , model.y = '(4/5) * pop.x + pop.eps'
                         , model.pi = '200 * pop.x / (3.5 * N)'
                         , model.logit.phi = '8 - 2 * pop.x'
                         , sigma.est = 'Hajekplus'
                         , draw.u = 'hotdeck'
                         , informative = FALSE
                         , runKY = TRUE
                         , fpcKY = TRUE
                         , approxKY = TRUE
                         , runTradMI = TRUE
                         , writeExcel = TRUE
                         , nameExcel = 'sim8a_HD_no_icpt_20230305'
                         , seed = 20230305
)


resMIb1 <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 10000L
                        , cfs = 2L
                        , use.intercept = TRUE
                        , m = 100L
                        , model.x = 'runif(n = N, min = 2, max = 5)'
                        , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                        , model.y = '(4/5) * pop.x + pop.eps'
                        , model.pi = '200 * pop.y[ ,s] / (2.8 * N)'
                        , model.logit.phi = '8 - 2 * pop.x'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim8b_HD_20230305'
                        , seed = 20230305
)

resMIb0 <- runSimulation(S = 100000L
                         , S_at_a_time = 10000L
                         , N = 10000L
                         , cfs = 1L
                         , use.intercept = FALSE
                         , m = 100L
                         , model.x = 'runif(n = N, min = 2, max = 5)'
                         , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                         , model.y = '(4/5) * pop.x + pop.eps'
                         , model.pi = '200 * pop.y[ ,s] / (2.8 * N)'
                         , model.logit.phi = '8 - 2 * pop.x'
                         , sigma.est = 'Hajekplus'
                         , draw.u = 'hotdeck'
                         , informative = TRUE
                         , runKY = TRUE
                         , fpcKY = TRUE
                         , approxKY = TRUE
                         , runTradMI = TRUE
                         , writeExcel = TRUE
                         , nameExcel = 'sim8b_HD_no_icpt_20230305'
                         , seed = 20230305
)


resMIc1 <- runSimulation(S = 100000L
                         , S_at_a_time = 5000L
                         , N = 20000L
                         , cfs = 2L
                         , use.intercept = TRUE
                         , m = 100L
                         , model.x = 'runif(n = N, min = 2, max = 5)'
                         , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                         , model.y = '(4/5) * pop.x + pop.eps'
                         , model.pi = '400 * pop.y[ ,s] / (2.8 * N)'
                         , model.logit.phi = '8 - 2 * pop.x'
                         , sigma.est = 'Hajekplus'
                         , draw.u = 'hotdeck'
                         , informative = TRUE
                         , runKY = TRUE
                         , fpcKY = TRUE
                         , approxKY = TRUE
                         , runTradMI = TRUE
                         , writeExcel = TRUE
                         , nameExcel = 'sim8c_HD_20230305'
                         , seed = 20230305
)

resMIc0 <- runSimulation(S = 100000L
                         , S_at_a_time = 5000L
                         , N = 20000L
                         , cfs = 1L
                         , use.intercept = FALSE
                         , m = 100L
                         , model.x = 'runif(n = N, min = 2, max = 5)'
                         , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                         , model.y = '(4/5) * pop.x + pop.eps'
                         , model.pi = '400 * pop.y[ ,s] / (2.8 * N)'
                         , model.logit.phi = '8 - 2 * pop.x'
                         , sigma.est = 'Hajekplus'
                         , draw.u = 'hotdeck'
                         , informative = TRUE
                         , runKY = TRUE
                         , fpcKY = TRUE
                         , approxKY = TRUE
                         , runTradMI = TRUE
                         , writeExcel = TRUE
                         , nameExcel = 'sim8c_HD_no_icpt_20230305'
                         , seed = 20230305
)


##### Methode Kim & Yang

resMIa1.KY <- runSimulation(S = 10000L
                         , S_at_a_time = 10000L
                         , N = 10000L
                         , cfs = 2L
                         , use.intercept = TRUE
                         , m = 100L
                         , model.x = 'runif(n = N, min = 2, max = 5)'
                         , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                         , model.y = '(4/5) * pop.x + pop.eps'
                         , model.pi = '200 * pop.x / (3.5 * N)'
                         , model.logit.phi = '8 - 2 * pop.x'
                         , sigma.est = 'Hajekplus'
                         , draw.u = 'hotdeck'
                         , informative = FALSE
                         , runKY = TRUE
                         , fpcKY = TRUE
                         , approxKY = FALSE
                         , writeExcel = TRUE
                         , nameExcel = 'sim8a_KY_20230305'
                         , seed = 20230305
)

resMIa0.KY <- runSimulation(S = 10000L
                         , S_at_a_time = 10000L
                         , N = 10000L
                         , cfs = 1L
                         , use.intercept = FALSE
                         , m = 100L
                         , model.x = 'runif(n = N, min = 2, max = 5)'
                         , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                         , model.y = '(4/5) * pop.x + pop.eps'
                         , model.pi = '200 * pop.x / (3.5 * N)'
                         , model.logit.phi = '8 - 2 * pop.x'
                         , sigma.est = 'Hajekplus'
                         , draw.u = 'hotdeck'
                         , informative = FALSE
                         , runKY = TRUE
                         , fpcKY = TRUE
                         , approxKY = FALSE
                         , writeExcel = TRUE
                         , nameExcel = 'sim8a_KY_no_icpt_20230305'
                         , seed = 20230305
)


resMIb1.KY <- runSimulation(S = 10000L
                         , S_at_a_time = 10000L
                         , N = 10000L
                         , cfs = 2L
                         , use.intercept = TRUE
                         , m = 100L
                         , model.x = 'runif(n = N, min = 2, max = 5)'
                         , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                         , model.y = '(4/5) * pop.x + pop.eps'
                         , model.pi = '200 * pop.y[ ,s] / (2.8 * N)'
                         , model.logit.phi = '8 - 2 * pop.x'
                         , sigma.est = 'Hajekplus'
                         , draw.u = 'hotdeck'
                         , informative = TRUE
                         , runKY = TRUE
                         , fpcKY = TRUE
                         , approxKY = FALSE
                         , writeExcel = TRUE
                         , nameExcel = 'sim8b_KY_20230305'
                         , seed = 20230305
)

resMIb0.KY <- runSimulation(S = 10000L
                         , S_at_a_time = 10000L
                         , N = 10000L
                         , cfs = 1L
                         , use.intercept = FALSE
                         , m = 100L
                         , model.x = 'runif(n = N, min = 2, max = 5)'
                         , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                         , model.y = '(4/5) * pop.x + pop.eps'
                         , model.pi = '200 * pop.y[ ,s] / (2.8 * N)'
                         , model.logit.phi = '8 - 2 * pop.x'
                         , sigma.est = 'Hajekplus'
                         , draw.u = 'hotdeck'
                         , informative = TRUE
                         , runKY = TRUE
                         , fpcKY = TRUE
                         , approxKY = FALSE
                         , writeExcel = TRUE
                         , nameExcel = 'sim8b_KY_no_icpt_20230305'
                         , seed = 20230305
)
