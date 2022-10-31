library("Matrix")
library("huge")
library("testit")
library("gdata")
set.seed(111)

e_step <- function(omega, pi, v0, v1) {
  a <- dnorm(omega, 0, v1) * pi
  b <- dnorm(omega, 0, v0) * (1 - pi)
  p <- a / (a + b)
  d <- ((1 - p) / v0 ^ 2) + (p / v1 ^ 2)
  return(list("p" = p, "d" = d))
}

cm_step <-
  function(alpha,
           beta,
           lambda,
           v0,
           v1,
           s,
           pi,
           omega,
           n,
           estep) {
    p <- nrow(omega)
    pm1 <- p - 1
    sum_delta <- sum(upperTriangle(estep$p))
    pp2 <- p * (p - 1) / 2
    new_pi <- (alpha + sum_delta - 1) / (alpha + beta + pp2 - 2)
    
    for (t in 1:p) {
      print(omega)
      # Swap the columns
      omega[, c(t, p)] <- omega[, c(p, t)]
      estep$d[, c(t, p)] <- estep$d[, c(p, t)]
      s[, c(t, p)] <- s[, c(p, t)]
      
      # Swap the rows
      omega[c(t, p),] <- omega[c(p, t),]
      estep$d[c(t, p),] <- estep$d[c(p, t),]
      s[c(t, p),] <- s[c(p, t),]
      
      # Get relevant values
      s_22 <- s[p, p]
      s_12 <- s[1:pm1, p]
      d_12 <- diag(p - 1) * estep$d[1:pm1, p]
      omega_11_inv <- solve(omega[1:pm1, 1:pm1])
      c_inv <- solve((s_22 + lambda) * omega_11_inv + d_12)
      
      # Update SNS part
      omega_12_lp1 <- -c_inv %*% s_12
      omega[1:pm1, p] <- omega_12_lp1
      omega[p, 1:pm1] <- t(omega_12_lp1)
      
      # Update exp bit
      v <- t(omega_12_lp1) %*% omega_11_inv %*% omega_12_lp1
      omega[p, p] <- v + n/(lambda + s_22)
      
      # Swap back everything
      omega[, c(t, p)] <- omega[, c(p, t)]
      estep$d[, c(t, p)] <- estep$d[, c(p, t)]
      s[, c(t, p)] <- s[, c(p, t)]
      omega[c(t, p),] <- omega[c(p, t),]
      estep$d[c(t, p),] <- estep$d[c(p, t),]
      s[c(t, p),] <- s[c(p, t),]
    }
    return(list("omega" = omega, "pi" = new_pi))
  }

log_lik <- function(alpha, beta, lambda, pi, omega, n, estep) {
  sns_part <- -sum(upperTriangle(omega ^ 2 * estep$d)) / 2
  exp_part <- -sum(lambda * diag(omega)) / 2
  ber_part <-
    sum(upperTriangle(log(1 - pi) + (estep$p * log(pi / (
      1 - pi
    )))))
  misc_part <-
    (n / 2 * log(det(omega))) - (sum(diag(s %*% omega)) / 2)
  + ((alpha - 1) * log(pi)) + ((beta - 1) * log(1 - pi))
  
  return(sns_part + exp_part + ber_part + misc_part)
}

v0 <- 0.1
v1 <- 10
alpha <- 1
beta <- 1
lambda <- 1

d <- 6
n <- 200
steps <- 30

random <- huge.generator(n = n, d = d, prob = 1/d)
x <- random$data
write.csv(x, "x.csv", row.names = FALSE)
s <- t(x) %*% x

omega <- diag(d)
pi <- 0.5
log_liks <- c()

for (t in 1:steps) {
  estep <- e_step(omega, pi, v0, v1)
  cmstep <-
    cm_step(alpha, beta, lambda, v0, v1, s, pi, omega, n, estep)
  omega <- cmstep$omega
  pi <- cmstep$pi
  l <- log_lik(alpha, beta, lambda, pi, omega, n, estep)
  log_liks <- append(log_liks, l)
}
