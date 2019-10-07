generateGeometricBrownianMotionsEnd <- function(S0, r, sigma, TTM, n, M){
  unifQMC   <- runif.halton(ceiling(M/2), 2, 1)
    as.numeric(
      apply(unifQMC,
        1,
        FUN = function(x) transformUniformToNormalByBoxMuller(x[1], x[2])))[1:M]
}
