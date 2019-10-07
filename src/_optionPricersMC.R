getAsianOptionPriceMC <- function(S0, r, sigma, TTM, K, n, M = 10000, type){
  if(type == "call"){
    payoffFunction <- function(x) ifelse(mean(x) > K, mean(x) - K, 0)
  } else if( type == "put"){
    payoffFunction <- function(x) ifelse(K > mean(x), K - mean(x), 0)
  } else {
    stop("Wrong option type!")
  }
  
  deltaT            <- TTM / n
  discountedPayoffs <- numeric(M)
  
  for (i in 1:M)
    discountedPayoffs[i] <- payoffFunction(S0 * exp(cumsum((r - 0.5 * sigma ^ 2) * deltaT +
                                                             sqrt(deltaT) * sigma * rnorm(n)))) * 
    exp(-r * TTM)
  
  price             <- mean(discountedPayoffs) 
  confInterval      <- price + c(-1, 1) * 1.96 * sd(discountedPayoffs) / sqrt(M)
  
  result            <- list(price, confInterval)
  names(result)     <- c("price", "confInterval95")
  return(result)
}

getAsianOptionPriceMCAdv <- function(S0, r, sigma, TTM, K, n, M = NULL, type){
  # verify type
  if (type == "call") {
    payoffFunction         <- function(x) ifelse(mean(x) > K, mean(x) - K, 0)
    auxilaryPayoffFunction <- function(x) ifelse(exp(mean(log(x))) > K, exp(mean(log(x))) - K, 0)
  } else if (type == "put") {
    payoffFunction         <- function(x) ifelse(K > mean(x), K - mean(x), 0)
    auxilaryPayoffFunction <- function(x) ifelse(K > exp(mean(log(x))), K - exp(mean(log(x))), 0)
  } else {
    stop("Wrong option type!")
  }

  if(M = NULL){
    warning("The number of Monte Carlo trajectories simulation (M) has not been specified.\n
           The default value of 10000 is used.\n")
    M <- 10000
  }
  
  # vanilla case
  if (n == 1) {
    price             <- mean(discountedPayoffs) 
    confInterval      <- price + c(-1, 1) * 1.96 * sd(discountedPayoffs) / sqrt(M)
  } else{ # stricte asian option
    brownianMotions           <- generateAntitheticGeometricBrownianMotions(S0, r, sigma, TTM, K, n, M)
    discountedPayoffs         <- apply(brownianMotions,
                                       1,
                                       FUN = function(x) payoffFunction(x)) * exp(-r * TTM)
    auxilaryDiscountedPayoffs <- apply(brownianMotions,
                                       1,
                                       FUN = function(x) auxilaryPayoffFunction(x))  * exp(-r * TTM)
    auxilaryPrice             <- calculatePriceAsianGeometricDiscreteCall(S0, r, sigma, TTM, K, n)
    
    theta           <- -cov(discountedPayoffs, auxilaryDiscountedPayoffs) / var(auxilaryDiscountedPayoffs)
    modifiedPayoffs <- discountedPayoffs + theta * (auxilaryDiscountedPayoffs - auxilaryPrice)
    price           <- mean(modifiedPayoffs)
    
    confInterval <- price + c(-1, 1) * 1.96 * sd(modifiedPayoffs) / sqrt(M)
  }
 
  result            <- list(price, confInterval)
  names(result)     <- c("price", "confInterval95")
  return(result)
}
