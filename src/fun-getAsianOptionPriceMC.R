getAsianOptionPriceMC <- function(S0, r, sigma, TTM, K, n, M = NULL, type, confLevel = 0.95){
  # Verifying type.
  if(type == "call"){
    payoffFunction <- function(x) ifelse(mean(x) > K, mean(x) - K, 0)
  } else if( type == "put"){
    payoffFunction <- function(x) ifelse(K > mean(x), K - mean(x), 0)
  } else {
    stop("Wrong option type. Please set type as a 'call' or 'put'.")
  }
  
  # Verifying M and use default if needed.
  if(is.null(M)){
    warning("The number of Monte Carlo trajectories simulation (M) has not
been specified.\n The default value of 10 000 is used.\n")
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
    
  # Predefined values.
  deltaT            <- TTM / n
  discountedPayoffs <- numeric(M)
  
  # Running the Monte-Carlo engine.
  for (i in 1:M)
    discountedPayoffs[i] <- 
    payoffFunction(S0 * exp(cumsum((r - 0.5 * sigma ^ 2) * deltaT +
                                     sqrt(deltaT) * sigma * rnorm(n)))) * 
    exp(-r * TTM)
  
  # Calculating price and confidence interval.
  price             <- mean(discountedPayoffs) 
  confInterval      <- price + c(-1, 1) * qnorm(confLevel + (1-confLevel) / 2) *
    sd(discountedPayoffs) / sqrt(M)
  
  result            <- list(price, confInterval)
  names(result)     <- 
    c("price", paste0("confInterval", format(100*confLevel, digits = 4), "%"))
  
  result
}
