require(Rcpp)
require(magrittr)
require(coda)
require(rbenchmark)
library(here)

sourceCpp(here("extra", "ExactHMC.cpp"))

exactHMC <- function(n,
                     initial_position,
                     constraint_direc,
                     constraint_bound,
                     precisionMat,
                     mu,
                     precision = TRUE,
                     total_time = pi / 2,
                     seed = 666) {
  set.seed(seed)
  cholesky <- Cholesky(precisionMat)
  results <- matrix(nrow = n, ncol = ncol(constraint_direc))
  whitened_constraints <- WhitenConstraints(constraint_direc,
                                            constraint_bound,
                                            cholesky,
                                            mu,
                                            precision)
  sample <- initial_position
  for (i in 1:n) {
    initial_momentum <- rnorm(ncol(constraint_direc))
    sample <- GenerateSample(
      sample,
      initial_momentum,
      whitened_constraints$direc,
      whitened_constraints$direc_rownorm_sq,
      whitened_constraints$bound,
      cholesky,
      mu,
      total_time,
      precision
    )
    results[i, ] <- sample
  }
  return(results)
}


#' Get the minimal ESS
#'
#' @param res an array containing MCMC samples (each row is an iteration)
#'
#' @return the minimal ESS across all dimensions
#'
#' @examples
getMinESS <- function(res) {
  return(mcmc(res, thin = 1) %>% effectiveSize() %>% min())
}


benchMarkTMVN <-
  function(nEXACT = 100,
           nTN = 100,
           nHZZ = 100,
           nNUTS = 100,
           repeatTimes = 1,
           forOneSampleFlg = T,
           returnSamplesFlg = F,
           dimension,
           meanVec,
           precMat,
           covMat,
           pInitial,
           constraintDirec,
           constraintBdry,
           lbTN,
           ubTN,
           constraintHZZ) {
    timingCols = c("test",
                   "replications",
                   "elapsed",
                   "relative",
                   "user.self",
                   "sys.self")
    
    timeList <- list()
    samplesList <- list()
    numTest <- 1
    ## ExactHMC ####
    if (nEXACT > 0) {
      timeList[[numTest]] <- benchmark(
        "exactHMC" = {
          samplesEXACT <- exactHMC(
            n = nEXACT,
            initial_position = pInitial,
            constraint_direc = constraintDirec,
            #diag(dimension),
            constraint_bound = constraintBdry,
            #rep(0, dimension),
            precisionMat = precMat,
            mu = meanVec,
            precision = TRUE,
            total_time = pi / 2,
            seed = 666
          )
        },
        replications = repeatTimes,
        columns = timingCols,
        order = NULL
      )
      samplesList <-
        c(samplesList, samplesEXACT = list(samplesEXACT))
      numTest <- numTest + 1
    }
    ## TN ####
    if (nTN > 0) {
      set.seed(666)
      timeList[[numTest]] <- benchmark(
        "TN" = {
          samplesTN <- TruncatedNormal::rtmvnorm(
            n = nTN,
            mu = meanVec,
            sigma = covMat,
            lb = lbTN,
            ub = ubTN
          )
        },
        replications = repeatTimes,
        columns = timingCols,
        order = NULL
      )
      samplesList <- c(samplesList, samplesTN = list(samplesTN))
      numTest <- numTest + 1
    }
    ## HZZ and NUTS (one e-sample) ####
    if (forOneSampleFlg) {
      if (nHZZ > 0) {
        set.seed(666)
        engine <- createEngine(
          dimension = dimension,
          mask = rep(1, dimension),
          observed = rep(1, dimension),
          parameterSign = constraintHZZ,
          flags = 128,
          info = 1,
          seed = 666
        )
        
        samplesHZZ  <-  array(0, c(nHZZ, dimension))
        
        timeList[[numTest]] <- benchmark(
          "HZZ" = {
            setMean(sexp = engine$engine, mean = meanVec)
            setPrecision(sexp = engine$engine, precision = precMat)
            
            HZZtime <-
              sqrt(2) / sqrt(min(mgcv::slanczos(
                A = precMat, k = 1, kl = 1
              )[['values']]))
            
            for (i in 1:nHZZ) {
              momentum <- drawMomentum(dimension)
              samplesHZZ[i,] <- getSample(
                position = pInitial,
                momentum = momentum,
                t = HZZtime,
                nutsFlg = F,
                engine = engine
              )
            }
          },
          replications = repeatTimes,
          columns = timingCols,
          order = NULL
        )
        samplesList <- c(samplesList, samplesHZZ = list(samplesHZZ))
        numTest <- numTest + 1
      }
      
      if (nNUTS > 0) {
        set.seed(666)
        baseStep <-
          0.1 / sqrt(min(mgcv::slanczos(
            A = precMat, k = 1, kl = 1
          )[['values']]))
        engine <- createNutsEngine(
          dimension = dimension,
          mask = rep(1, dimension),
          observed = rep(1, dimension),
          parameterSign = constraintHZZ,
          flags = 128,
          info = 1,
          seed = 666,
          randomFlg = T,
          stepSize = baseStep,
          mean = meanVec,
          precision = precMat
        )
        
        samplesNUTS  <-  array(0, c(nNUTS, dimension))
        
        timeList[[numTest]] <- benchmark(
          "NUTS" = {
            baseStep <-
              0.1 / sqrt(min(mgcv::slanczos(
                A = precMat, k = 1, kl = 1
              )[['values']]))
            for (i in 1:nNUTS) {
              momentum <- drawMomentum(dimension)
              samplesNUTS[i,] <- getSample(
                position = pInitial,
                momentum = momentum,
                t = baseStep,
                nutsFlg = T,
                engine = engine
              )
            }
          },
          replications = repeatTimes,
          columns = timingCols,
          order = NULL
        )
        samplesList <-
          c(samplesList, samplesNUTS = list(samplesNUTS))
        numTest <- numTest + 1
      }
    } else {
      ## HZZ and NUTS (multiple e-sample) ####
      if (nHZZ > 0) {
        set.seed(666)
        timeList[[numTest]] <- benchmark(
          "HZZ" = {
            samplesHZZ <- rcmg(
              n = nHZZ,
              mean = meanVec,
              prec = precMat,
              constraits = constraintHZZ,
              t = 1,
              burnin = 0,
              cppFlg = T,
              nutsFlg = F
            )
          },
          replications = repeatTimes,
          columns = timingCols,
          order = NULL
        )
        samplesList <- c(samplesList, samplesHZZ = list(samplesHZZ))
        numTest <- numTest + 1
      }
      
      if (nNUTS > 0) {
        set.seed(666)
        timeList[[numTest]] <- benchmark(
          "NUTS" = {
            samplesNUTS <- rcmg(
              n = nNUTS,
              mean = meanVec,
              prec = precMat,
              constraits = constraintHZZ,
              t = 1,
              burnin = 0,
              cppFlg = T,
              nutsFlg = T
            )
          },
          replications = repeatTimes,
          columns = timingCols,
          order = NULL
        )
        samplesList <-
          c(samplesList, samplesNUTS = list(samplesNUTS))
        numTest <- numTest + 1
      }
    }
    ## return list ####
    if (returnSamplesFlg) {
      return(list(timeList = timeList,
                  samplesList = samplesList))
    } else{
      return(list(timeList = timeList,
                  samplesList = NA))
    }
  }
