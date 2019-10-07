getAsianOptionPriceMCAdv <- 
  function(S0, r, sigma, TTM, K, n, M = NULL, type, confLevel = 0.95){
    
  # Verifying type.
  if (type == "call") {
    payoffFunction         <- 
      function(x) ifelse(mean(x) > K, mean(x) - K, 0)
    
    auxilaryPayoffFunction <- 
      function(x) ifelse(exp(mean(log(x))) > K, exp(mean(log(x))) - K, 0)
    
  } else if (type == "put") {
    
    payoffFunction         <- 
      function(x) ifelse(K > mean(x), K - mean(x), 0)
    
    auxilaryPayoffFunction <- 
      function(x) ifelse(K > exp(mean(log(x))), K - exp(mean(log(x))), 0)
    
  } else {
    stop("Wrong option type. Please set type as a 'call' or 'put'.")
  }
  
  # Verifying M and use default if needed.
  if(is.null(M)){
    warning("The number of Monte Carlo trajectories simulation (M) has not been
specified.\n The default value of 10 000 is used.\n")
    M <- 10000
  }
  
  # Verifying parameters.
  verify <- verifyMCParameters(S0, r, sigma, TTM, K, n, M, confLevel)
  if (verify$assert == 1) {
    stop("Execution has been haulted. \n")
  } else{
    n         <- verify$n
    M         <- verify$M
    confLevel <- verify$confLevel
  }
    
  if (n == 1) { # Plain vanilla option.
    warning("The considered option is vanilla, so the Halton series has been
            used to improve precision.\n")
    
    # Redefining payoff functions.
    if (type == "call") {
      payoffFunction         <- function(x) ifelse(x > K, x - K, 0)
    } else if (type == "put") {
      payoffFunction         <- function(x) ifelse(K > x, K - x, 0)
    }
    
    brownianMotions   <- 
      generateGeometricBrownianMotionsEnd(S0, r, sigma, TTM, M)
    
    discountedPayoffs <- 
      payoffFunction(brownianMotions) * exp(-r * TTM)
    
    price             <- 
      mean(discountedPayoffs) 
    
    confInterval      <- 
      price + c(-1, 1) * qnorm(confLevel + (1-confLevel) / 2) *
      sd(discountedPayoffs) / sqrt(M)
    
  } else{ # Stricte Asian option.
    brownianMotions           <- 
      generateAntitheticGeometricBrownianMotions(S0, r, sigma, TTM, K, n, M)
    
    discountedPayoffs         <- 
      apply(brownianMotions,
            1,
            FUN = function(x) payoffFunction(x)) * exp(-r * TTM)
    
    auxilaryDiscountedPayoffs <- 
      apply(brownianMotions,
            1,
            FUN = function(x) auxilaryPayoffFunction(x))  * exp(-r * TTM)
    
    if(type == "call"){
      auxilaryPrice <- 
        calculatePriceAsianGeometricDiscreteCall(S0, r, sigma, TTM, K, n)
    } else {
      auxilaryPrice <- 
        calculatePriceAsianGeometricDiscretePut(S0, r, sigma, TTM, K, n)
    }
   
    theta           <- 
      -cov(discountedPayoffs, auxilaryDiscountedPayoffs) /
      var(auxilaryDiscountedPayoffs)
    
    modifiedPayoffs <- 
      discountedPayoffs + theta * (auxilaryDiscountedPayoffs - auxilaryPrice)
    
    price           <- 
      mean(modifiedPayoffs)
    
    confInterval <- 
      price + c(-1, 1) * qnorm(confLevel + (1-confLevel) / 2) *
      sd(modifiedPayoffs) / sqrt(M)
  }
  
  result            <- list(price, confInterval)
  names(result)     <- 
    c("price", paste0("confInterval", format(100*confLevel, digits = 4), "%"))
  result
}
