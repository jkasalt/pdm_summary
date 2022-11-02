library("Matrix")
library("huge")
library("testit")
library("gdata")
library("matrixcalc")
set.seed(111)

e_step <- function(omega, pi, v0, v1) {
  a <- dnorm(omega, 0, v1) * pi
  b <- dnorm(omega, 0, v0) * (1 - pi)
  p <- a / (a + b)
  d <- (1 - p) / (v0 ^ 2) + p / (v1 ^ 2)
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
      omega[p, p] <- v + n / (lambda + s_22)
      
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

build_path <- function(p, len = 20) {
  tt <- seq(1, 0, length.out = len)
  output <- c()
  for (t in tt) {
    g <- Matrix(p >= t, sparse = TRUE)
    output <- append(output, g)
  }
  return(output)
}

bmg_roc <- function(p, theta, len = 200) {
  path <- build_path(p, len)
  huge.roc(path, theta)
}

prior_sns <- function(p, lambda, v0, v1, pi) {
  pp2 <- p * (p - 1) / 2
  repeat {
    # Build omega matrix using prior
    ## Build delta
    delta <- rbinom(pp2, 1, pi)
    
    ## Build omega
    omega <- matrix(0, nrow = p, ncol = p)
    ### Set upper triangle
    omega[upper.tri(omega)] <-
      delta * rnorm(pp2, 0, v1) ^ 2 + (1 - delta) * rnorm(pp2, 0, v0) ^ 2
    ### Set lower triangle
    omega <- omega + t(omega)
    ### Set diagonal
    diag(omega) <- rexp(p, lambda / 2)
    
    if (is.positive.definite(omega)) {
      break
    }
  }
  return(omega)
}

ecm <- function(x) {
  v0 <- 0.01
  v1 <- 10
  a <- 1
  b <- 1
  lambda <- 1
  p <- ncol(x)
  pp2 <- p * (p - 1) / 2
  steps <- 1000
  
  pi <- 0.5
  # omega <- prior_sns(p, lambda, v0, v1, pi)
  omega <- diag(p)
  
  # Do the actual algo
  log_liks <- c()
  for (t in 1:steps) {
    assert("omega must be positive definite",
           is.positive.definite(omega))
    estep <- e_step(omega, pi, v0, v1)
    cmstep <-
      cm_step(a, b, lambda, v0, v1, s, pi, omega, n, estep)
    omega <- cmstep$omega
    pi <- cmstep$pi
    l <- log_lik(a, b, lambda, pi, omega, n, estep)
    print(l)
    log_liks <- append(log_liks, l)
    
    # Return early if relative increase in likelihood is too low
    relative_increase <-
      log_liks[length(log_liks) - 1] / log_liks[length(log_liks)]
    if (relative_increase <= 1 + 1e-3) {
      break
    }
  }
  return(list(
    "omega" = omega,
    "prob" = estep$p,
    "log_liks" = log_liks
  ))
}

d <- 50
n <- 2000
random <- huge.generator(n = n, d = d)
x <- random$data
write.csv(x, "x.csv", row.names = FALSE)
s <- t(x) %*% x
random_ecm <- ecm(x)
#
# cluster <- huge.generator(n, d, "cluster")
# cluster_ecm <- ecm(cluster$data, omega)
