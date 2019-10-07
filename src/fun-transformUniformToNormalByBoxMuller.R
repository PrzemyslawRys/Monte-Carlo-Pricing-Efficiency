transformUniformToNormalByBoxMuller <- function(u1, u2){
  r     <- sqrt(-2 * log(u1))
  theta <- 2 * pi * u2
  
  c(r * cos(theta), r * sin(theta))
}