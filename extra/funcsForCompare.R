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
                     jittering = T,
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
    if (jittering) {
      total_time <- runif(1, pi / 8, pi / 2)
    }
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
  function(nEHMC = 100,
           nZHMC = 100,
           nZNUTS = 100,
           nMT = 100,
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
    if (nEHMC > 0) {
      timeList[[numTest]] <- benchmark(
        "EHMC" = {
          samplesEHMC <- exactHMC(
            n = nEHMC,
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
        c(samplesList, samplesEHMC = list(samplesEHMC))
      numTest <- numTest + 1
    }

    ## HZZ and NUTS (one e-sample) ####
    if (forOneSampleFlg) {
      if (nZHMC > 0) {
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
        
        samplesZHMC  <-  array(0, c(nZHMC, dimension))
        
        timeList[[numTest]] <- benchmark(
          "ZHMC" = {
            setMean(sexp = engine$engine, mean = meanVec)
            setPrecision(sexp = engine$engine, precision = precMat)
            
            HZZtime <-
              sqrt(2) / sqrt(min(mgcv::slanczos(
                A = precMat, k = 1, kl = 1
              )[['values']]))
            
            for (i in 1:nZHMC) {
              momentum <- drawMomentum(dimension)
              samplesZHMC[i,] <- getSample(
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
        samplesList <- c(samplesList, samplesZHMC = list(samplesZHMC))
        numTest <- numTest + 1
      }
      
      if (nZNUTS > 0) {
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
        
        samplesZNUTS  <-  array(0, c(nZNUTS, dimension))
        
        timeList[[numTest]] <- benchmark(
          "ZNUTS" = {
            baseStep <-
              0.1 / sqrt(min(mgcv::slanczos(
                A = precMat, k = 1, kl = 1
              )[['values']]))
            for (i in 1:nZNUTS) {
              momentum <- drawMomentum(dimension)
              samplesZNUTS[i,] <- getSample(
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
          c(samplesList, samplesZNUTS = list(samplesZNUTS))
        numTest <- numTest + 1
      }
    } else {
      ## HZZ and NUTS (multiple e-sample) ####
      if (nZHMC > 0) {
        set.seed(666)
        timeList[[numTest]] <- benchmark(
          "ZHMC" = {
            samplesZHMC <- rcmg(
              n = nZHMC,
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
        samplesList <- c(samplesList, samplesZHMC = list(samplesZHMC))
        numTest <- numTest + 1
      }
      
      if (nZNUTS > 0) {
        set.seed(666)
        timeList[[numTest]] <- benchmark(
          "ZNUTS" = {
            samplesZNUTS <- rcmg(
              n = nZNUTS,
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
          c(samplesList, samplesZNUTS = list(samplesZNUTS))
        numTest <- numTest + 1
      }
    }
    
    ## TN ####
    if (nMT > 0) {
      set.seed(666)
      timeList[[numTest]] <- benchmark(
        "MT" = {
          samplesMT <- TruncatedNormal::rtmvnorm(
            n = nMT,
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
      samplesList <- c(samplesList, samplesMT = list(samplesMT))
      numTest <- numTest + 1
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
