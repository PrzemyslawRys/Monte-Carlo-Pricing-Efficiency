calculatePriceAsianGeometricDiscretePut <- function(S0, r, sigma, TTM, K, n){
  b1 <- (log(S0 / K) + (r - 0.5 * sigma ^ 2) * TTM * (n + 1) / (2 * n) +
           sigma ^ 2 * TTM * (n + 1) * (2 * n + 1) / (6 * n ^ 2)) /
    ((sigma / n) * sqrt(TTM * (n + 1) * (2 * n + 1) / 6))
  
  b2 <- (log(S0 / K) + (r - 0.5 * sigma ^ 2) * TTM * (n + 1) / (2 * n)) /
    ((sigma / n) * sqrt(TTM * (n + 1) * (2 * n + 1) / 6))
  
  exp(-r * TTM) * K * pnorm(-b2) - 
           exp(-r * TTM + (r - 0.5 * sigma ^ 2) * TTM * (n + 1 ) / (2 * n) + 
      sigma ^ 2 * TTM * (n + 1) * (2 * n + 1)/ (12 * n  ^ 2)) * S0 * pnorm(-b1)
  
}