library("Matrix")
library("huge")
library("testit")
library("gdata")
library("matrixcalc")
set.seed(111)

e_step <- function(omega, pi, v0, v1) {
  a <- dnorm(omega, 0, v1, log = TRUE) + log(pi)
  b <- dnorm(omega, 0, v0, log = TRUE) + log(1 - pi)
  m <- pmax(a, b)
  p <- exp(a - m) / (exp(a - m) + exp(b - m))
  diag(p) <- 0
  d <- (1 - p) / (v0 ^ 2) + p / (v1 ^ 2)
  diag(d) <- 0
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

li_omega <- function(omega, lambda, S, Ed, N) {
  M <- nrow(omega)
  for (i in 1:M) {
    v <- N / (lambda + S[i, i])
    u <-
      -solve((lambda + S[i, i]) * solve(omega[-i, -i])
             + diag(Ed[i, -i])) %*% S[i, -i]
    omega[i, -i] <- u
    omega[-i, i] <- omega[i, -i]
    omega[i, i] <- v + t(u) %*% solve(omega[-i, -i]) %*% u
  }
  return(omega)
}

log_lik <- function(alpha, beta, lambda, pi, omega, n, estep, s) {
  eps <- .Machine$double.eps ^ 0.5
  p <- nrow(omega)
  sns_part <- -sum(upperTriangle(omega ^ 2 * estep$d)) / 2
  exp_part <- -lambda / 2 * sum(diag(omega))
  ber_part <-
    log(pi / (1 - pi)) * sum(upperTriangle(estep$p)) + p * (p - 1) / 2 * log(1 - pi)
  misc_part <-
    (n / 2 * log(det(omega))) - (sum(diag(s %*% omega)) / 2)
  + ((alpha - 1) * log(pi + eps)) + ((beta - 1) * log(1 - pi + eps))
  
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
  return(huge.roc(path, theta, verbose=FALSE))
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
      delta * rnorm(pp2, 0, v1) + (1 - delta) * rnorm(pp2, 0, v0)
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

select_v0 <- function(graph, n_folds=5) {
  v0s <- seq(0.005, 0.05, length.out = 40)
  n <- nrow(graph$data)
  folds = c(1, 1:n_folds * n / 5)
  vals_v0s <- c()
  for (v0 in v0s) {
    print(v0)
    vals <- c()
    for (i in 1:n_folds) {
      x <- graph$data[-folds[i]:-folds[i + 1], ]
      e <- ecm(x, cov(x), v0)
      val <- bmg_roc(e$prob, graph$theta)$F1
      vals <- append(vals, val)
    }
    vals_v0s <- append(vals_v0s, mean(vals))
  }
  return(v0s[which.max(vals_v0s)])
}

ecm <- function(x,
                sigmahat,
                v0,
                tol = 1e-3,
                maxiter = 1000) {
  v1 <- 100
  ridge <- 0.5
  a <- 1
  b <- 1
  lambda <- 1
  p <- ncol(x)
  pp2 <- p * (p - 1) / 2
  steps <- 60
  s <- t(x) %*% x
  
  pi <- 0.5
  # omega <- prior_sns(p, lambda, v0, v1, pi)
  # omega <- diag(p)
  omega <-
    as(nearPD(solve(sigmahat + ridge * diag(p)))$mat, "matrix")
  
  log_liks <- c()
  last_lik <- NULL
  for (t in 1:steps) {
    assert("omega must be positive definite",
           is.positive.definite(omega))
    estep <- e_step(omega, pi, v0, v1)
    cmstep <-
      cm_step(a, b, lambda, v0, v1, s, pi, omega, n, estep)
    #omega <- cmstep$omega
    omega <- li_omega(omega, lambda, s, estep$d, n)
    pi <- cmstep$pi
    l <- log_lik(a, b, lambda, pi, omega, n, estep, s)
    log_liks <- append(log_liks, l)
    
    # Return early if relative increase in likelihood is too low
    if (is.null(last_lik) || abs(last_lik - l) >= tol) {
      last_lik <- l
    } else {
      break
    }
  }
  return(list(
    "omega" = omega,
    "pi" = pi,
    "prob" = estep$p,
    "log_liks" = log_liks
  ))
}

f1s_emgs <- c()
f1s_huge <- c()
n <- 200
for(d in c(25, 50, 100, 200)) {
  graph <- huge.generator(n = n, d = d)
  write.csv(graph$data, "x.csv", row.names = FALSE)
  v0 <- select_v0(graph)
  cm <- ecm(graph$data, graph$sigmahat, v0)
  f1s_emgs <- append(f1s_emgs, bmg_roc(cm$prob, graph$theta)$F1)
  h <- huge(graph$data)
  f1s_huge <- append(f1s_huge, huge.roc(h$path, graph$theta)$F1)
}

compare <- function(graph, v0) {
  cm <- ecm(graph$data, graph$sigmahat, v0)
  bmg_roc(cm$prob, graph$theta)
  h <- huge(graph$data)
  huge.roc(h$path, graph$theta, verbose = FALSE)
}