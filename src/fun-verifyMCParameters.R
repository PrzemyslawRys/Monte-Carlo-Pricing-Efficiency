verifyMCParameters <- function(S0, r, sigma, TTM, K, n, M, confLevel){
  assertVar <- 0
  if (S0 <= 0) {
    message("Initial underlying price need to be positive, due to the
         Black-Scholes model assumptions.\n
         Please call function again, with the proper value.\n")
    assertVar <- 1
  }

  if (r <= 0)
    warning("The interest rate (r) is negative.\n.")
  
  if (sigma < 0) {
    message("The underlying price volatility (sigma) need to be positive. \n
         Please call function again, with the proper value.\n")
    assertVar <- 1
  }
  
  if (sigma == 0)
    warning("The underlying price is deterministic, due to the zero volatility
            (sigma).")
  
  if (TTM <= 0) {
    message("Time to maturity (TTM) need to be positive. \n
         Please call function again, with the proper value.\n")
    assertVar <- 1
  }
    
  if (n != round(n)) {
    warning(paste0("The number of periods for average purposes has been
                   rounded to", round(n), ".\n"))
    n <- round(n)
  }
  
  if (n < 1) {
    message("The number of periods for average purposes (n) need to higher
or equal to 1. \n
         Please call function again, with the proper value.\n")
    assertVar <- 1
  }
  
  
  if (M != round(M)) {
    warning(paste0("The number of Monte Carlo iterations (M) has been rounded to ", round(M), ".\n"))
    M <- round(M)
  }
  if(M <= 1){
    message("The number of Monte Carlo iterations (M) need to be higher or equal to 1. \n
         Please call function again, with the proper value.\n")
    assertVar <- 1
  }
    
  if (confLevel < 0.5 | confLevel > 1) {
    warning("The confidence level should be in (0.5, 1). The default value of 95% has been used. \n")
    confLevel <- 0.95
  }
  
  results        <- list(assertVar, n, M, confLevel)
  names(results) <- c("assert", "n", "M", "confLevel")
  return(results)
}