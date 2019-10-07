library(rbenchmark)
library(fOptions)

source("src/fun-generateAntitheticGeometricBrownianMotions.R")
source("src/fun-getAsianOptionPriceMCadv.R")
source("src/fun-generateGeometricBrownianMotionsEnd.R")
source("src/fun-generateGeometricBrownianMotions.R")
source("src/fun-calculatePriceAsianGeometricDiscreteCall.R")
source("src/fun-verifyMCParameters.R")
source("src/fun-transformUniformToNormalByBoxMuller.R")

S0    <- 50
r     <- 0.02
sigma <- 0.3
TTM   <- 1.5
K     <- 50
n     <- 100
M     <- 10000
type  <- "call"

confLevel <- 0.95

simulateAll <- function(S0, r, sigma, TTM, K, n, M){
  brownianMotions           <- 
    generateAntitheticGeometricBrownianMotions(S0, r, sigma, TTM, K, n, M)
  
  return(0)
}

brownianMotions        <- 
  generateAntitheticGeometricBrownianMotions(S0, r, sigma, TTM, K, n, M)

payoffFunction         <- function(x) ifelse(mean(x) > K, mean(x) - K, 0)

auxilaryPayoffFunction <- function(x) ifelse(exp(mean(log(x))) > K,
                                             exp(mean(log(x))) - K, 0)

restAll <- function(S0, r, sigma, TTM, K, n, M, brownianMotions){
  discountedPayoffs         <- 
    apply(brownianMotions,
          1,
          FUN = function(x) payoffFunction(x)) * exp(-r * TTM)
  
  auxilaryDiscountedPayoffs <- 
    apply(brownianMotions,
          1,
          FUN = function(x) auxilaryPayoffFunction(x))  * exp(-r * TTM)
  
  auxilaryPrice             <- 
    calculatePriceAsianGeometricDiscreteCall(S0, r, sigma, TTM, K, n)
  
  theta           <- 
    -cov(discountedPayoffs, auxilaryDiscountedPayoffs) /
    var(auxilaryDiscountedPayoffs)
  
  modifiedPayoffs <- 
    discountedPayoffs + theta * (auxilaryDiscountedPayoffs - auxilaryPrice)
  price           <- mean(modifiedPayoffs)
  
  confInterval <- price + c(-1, 1) * qnorm(confLevel + (1-confLevel) / 2) *
    sd(modifiedPayoffs) / sqrt(M)
  
  return(0)
}

restWithout <- function(S0, r, sigma, TTM, K, n, M, brownianMotions){
  discountedPayoffs         <- 
    apply(brownianMotions,
          1,
          FUN = function(x) payoffFunction(x)) * exp(-r * TTM)
 
   price           <- mean(discountedPayoffs) 
  
  confInterval <- price + c(-1, 1) * qnorm(confLevel + (1 - confLevel) / 2) *
    sd(discountedPayoffs) / sqrt(M)
  
  return(0)
}

generateGeometricBrownianMotionsForLoop <- function(S0, r, sigma, TTM, K, n, M){
  deltaT            <- TTM / M
  discountedPayoffs <- numeric(M)
  for (i in 1:M)
    discountedPayoffs[i] <-
      payoffFunction(S0 * exp(cumsum((r - 0.5 * sigma ^ 2) * deltaT + sqrt(deltaT) *
                                       sigma * rnorm(n)))) * exp(-r * TTM)
  
  return(0)
}

restWithout(S0, r, sigma, TTM, K, n, M, brownianMotions)

benchmark1 <- 
  benchmark(t1 <- verifyMCParameters(S0, r, sigma, TTM, K, n, M, confLevel),
            t2 <- simulateAll(S0, r, sigma, TTM, K, n, M),
            t3 <- restAll(S0, r, sigma, TTM, K, n, M, brownianMotions))

benchmark_rest <- 
  benchmark(t1 <- restAll(S0, r, sigma, TTM, K, n, M, brownianMotions),
            t2 <- restWithout(S0, r, sigma, TTM, K, n, M, brownianMotions))

benchmark_for <- 
  benchmark(
    t1 <- 
      generateGeometricBrownianMotions(S0, r, sigma, TTM, K, n, M),
    t2 <-
      generateAntitheticGeometricBrownianMotions(S0, r, sigma, TTM, K, n, M),
    t3 <- 
      generateGeometricBrownianMotionsForLoop(S0, r, sigma, TTM, K, n, M))


benchmark_halton <- 
  benchmark(
    t1 <- generateGeometricBrownianMotionsEnd(S0, r, sigma, TTM, 1, M),
    t2 <- generateGeometricBrownianMotionsForLoop(S0, r, sigma, TTM, K, 1, M),
    t3 <- generateGeometricBrownianMotions(S0, r, sigma, TTM, K, 1, M))

saveRDS(benchmark1, file = "data/benchmark1.Rds")
saveRDS(benchmark_rest, file = "data/benchmark_rest.Rds")
saveRDS(benchmark_for, file = "data/benchmark_for.Rds")
saveRDS(benchmark_halton, file = "data/benchmark_halton.Rds")
