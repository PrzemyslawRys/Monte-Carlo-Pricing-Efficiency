generateAntitheticGeometricBrownianMotions <- function(S0, r, sigma, TTM, K, n, M){
  deltaT <- TTM / n
  normM  <- matrix(rnorm(M*n), M, n)
  
  t(apply(rbind(normM, -normM),
          1,
          function(x) S0 * exp(cumsum((r - 0.5 * sigma ^ 2) *
                                        deltaT + sqrt(deltaT) * sigma * x))))

}
