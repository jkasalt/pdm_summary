library(huge)
library(Matrix)
library(MASS)
library(gdata)
library(testit)
library(matrixcalc)
library(LaplacesDemon)

K <- 4

hyperparams <- list(
  "v0" = 0.08,
  "v1" = 10,
  "lambda" = rep(2, K),
  "n" = c(50, 50, 50, 50),
  "p" = 20
)

## INPUTS:
## y is a 3D array that contains K  matrices of size n_k \times p
## Omega is a 3D array that contains K matrices of size p \times p
## RETURNS:
## todo!
bmg.ecm <- function(y, Omega, theta, Sigma, params, maxiter = 30) {
  S <- list()
  for (k in 1:K) {
    S[[k]] <- t(y[[k]]) %*% y[[k]]
  }

  log_liks <- c()
  # Run ECM
  for (t in 1:maxiter) {
    updated_vals <- bmg.mstep(Omega, theta, Sigma, S, params)
    Omega <- updated_vals$Omega
    theta <- updated_vals$theta
    
    # Get posterior likelihood
    log_lik <- bmg.posterior_likelihood(y, Omega, theta, Sigma, params)
    log_liks <- append(log_liks, log_lik)
  }

  # Get posterior edge inclusion probabilities
  estep <- bmg.estep(Omega, theta, params)

  # TODO: q is not posterior edge inclusion!
  return(list("theta" = theta, "Omega" = Omega, "prob" = estep$prob, "log_liks" = log_liks))
}

bmg.estep <- function(Omega, theta, params) {
  v0 <- params$v0
  v1 <- params$v1

  # TODO: maybe use log to make computation more stable
  a1 <- dnorm(Omega, mean = 0, sd = v1) * pnorm(theta)
  a0 <- dnorm(Omega, mean = 0, sd = v0) * (1 - pnorm(theta))
  pdeltais1 <- a1 / (a1 + a0)
  print(sum(pdeltais1 > 0.5))
  
  d <- (1 - pdeltais1) / v0^2 + pdeltais1 / v1^2

  mr1_log <- dnorm(theta, log = TRUE) - pnorm(theta, log.p = TRUE)
  mr0_log <- dnorm(theta, log = TRUE) - pnorm(theta, log.p = TRUE, lower.tail = F)

  q <- theta - exp(mr0_log) + pdeltais1 * (exp(mr1_log) + exp(mr0_log))
  return(list("d" = d, "q" = q, "prob" = pdeltais1))
}

bmg.posterior_likelihood <- function(y, Omega, theta, Sigma, params) {
  p <- params$p
  n <- params$n
  lambda <- params$lambda
  Sigma_inv <- solve(Sigma)
  estep <- bmg.estep(Omega, theta, params)
  const <- log(2 * pi) * K * p * (p-1) / 2 - log(det(Sigma)) * p * (p-1) / 4 
  
  # Sum the terms that depend only on k
  sum_k <- 0
  for (k in 1:K) {
    S_k <- t(y[[k]]) %*% y[[k]]
    Omega_k <- Omega[,,k]
    lambda_k <- lambda[k]
    this_k <- 0.5 * (n[k] * log(det(Omega_k)) - sum(diag(S_k %*% Omega_k)) - lambda_k * sum(diag(Omega_k)))
    sum_k <- sum_k + this_k
  }
  
  # Sum the terms that depend only on i,j
  sum_ij <- 0
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      this_ij <- - 0.5 * t(theta[i,j,]) %*% Sigma_inv %*% theta[i,j,]
      sum_ij <- sum_ij + this_ij
    }
  }
  
  # Sum the terms that depend on both i,j, and k
  omega_bit <- -0.5 * sum(Omega^2 * estep$d)
  theta_bit <- - 0.5 * sum(theta^2 - 2 * theta * estep$q)
  sum_ijk <-omega_bit + theta_bit
  
  return(sum_k + sum_ij + sum_ijk)
}

#' Calculates the next iterates for omega and theta
#'
#' @param Omega is a 3D array that contains K matrices of size p \times p.
#' @param theta is a 3D array that contains K matrices of size K \times K.
#' @param Sigma is a K \times K matrix.
#' @param S is a 3D array that contains K matrices of size p \times p that is obtained by t(Y) %*% Y
#' @param params are the parameters
#' 
#' @return A list containing the new values for theta and Omega, Omega is modified in-place
bmg.mstep <- function(Omega, theta, Sigma, S, params) {
  estep <- bmg.estep(Omega, theta, params)
  d <- estep$d
  q <- estep$q
  p <- params$p
  lambda <- params$lambda
  Sigma_inv <- solve(Sigma)
  new_theta <- array(dim = c(p, p, K))

  for (k in 1:K) {
    Omega_k <- Omega[, , k]
    # Compute new theta
    # TODO: check if there is a faster way to do this e.g. tensors?
    for (i in 1:p) {
      for (j in 1:p) {
        new_theta[i, j, k] <-
          (q[i, j, k] - Sigma_inv[k, -k] %*% theta[i, j, -k])  / (1 + Sigma_inv[k, k])
      }
    }

    # Compute new omega
    d_k <- d[, , k]
    S_k <- S[[k]]
    lambda_k <- lambda[k]
    pp2 <- p * (p - 1) / 2
    pm1 <- p - 1
    for (t in 1:p) {
      # Swap the columns
      Omega_k[, c(t, p)] <- Omega_k[, c(p, t)]
      d_k[, c(t, p)] <- d_k[, c(p, t)]
      S_k[, c(t, p)] <- S_k[, c(p, t)]

      # Swap the rows
      Omega_k[c(t, p), ] <- Omega_k[c(p, t), ]
      d_k[c(t, p), ] <- d_k[c(p, t), ]
      S_k[c(t, p), ] <- S_k[c(p, t), ]

      # Get relevant values
      S_k_22 <- S_k[p, p]
      S_k_12 <- S_k[1:pm1, p]
      d_k_12 <- diag(p - 1) * d_k[1:pm1, p]
      Omega_k_11_inv <- solve(Omega_k[1:pm1, 1:pm1])
      c_inv <- solve((S_k_22 + lambda_k) * Omega_k_11_inv + d_k_12)

      # Update SNS part
      Omega_k_12_lp1 <- -c_inv %*% S_k_12
      Omega_k[1:pm1, p] <- Omega_k_12_lp1
      Omega_k[p, 1:pm1] <- t(Omega_k_12_lp1)

      # Update exp part
      v <- t(Omega_k_12_lp1) %*% Omega_k_11_inv %*% Omega_k_12_lp1
      Omega_k[p, p] <- v + n[k] / (lambda_k + S_k_22)

      # Swap back everything
      Omega_k[, c(t, p)] <- Omega_k[, c(p, t)]
      d_k[, c(t, p)] <- d_k[, c(p, t)]
      S_k[, c(t, p)] <- S_k[, c(p, t)]
      Omega_k[c(t, p), ] <- Omega_k[c(p, t), ]
      d_k[c(t, p), ] <- d_k[c(p, t), ]
      S_k[c(t, p), ] <- S_k[c(p, t), ]
    }
    Omega[, , k] <- Omega_k
  }
  return(list("theta" = new_theta, "Omega" = Omega))
}

pos_def <- function(omega, eps=1e-4) {
  res <- TRUE
  eigenvalues <- eigen(omega)$values
  for (v in eigenvalues) {
    if (v < eps) {
      res <- FALSE
    }
  }
  return(res)
}

gen <- function(K=1, prob=0.1, which="random", params) {
  res <- list()
  p <- params$p
  n <- params$n
  for (k in 1:K) {
    graph <- huge.generator(n = n[k], d = p, graph = which)
    num_edges <- sum(graph$theta == TRUE) / 2
    while (TRUE) {
      g <- graph$theta
      g_omega <- graph$omega
      num_swap <- rbinom(1, size = num_edges, prob = prob)
      for (s in 1:num_swap) {
        assert(sum(g) == sum(graph$theta))
        edges_mask <- g & upper.tri(g)
        non_edges_mask <- !g & upper.tri(g)
        edges <- which(edges_mask, arr.ind = TRUE)
        non_edges <- which(non_edges_mask, arr.ind = TRUE)
        
        # Pick random edges
        swapping_off <- sample(1:nrow(edges), 1)
        
        # Pick random non-edges
        swapping_on <- sample(1:nrow(non_edges), 1)
        
        # Execute swap
        pos_off <- edges[swapping_off,]
        pos_on <- non_edges[swapping_on,]
        
        g[pos_off[1], pos_off[2]] <- 0
        g[pos_on[1], pos_on[2]] <- 1
        
        tmp <- g_omega[pos_off[1], pos_off[2]]
        g_omega[pos_off[1], pos_off[2]] <-
          g_omega[pos_on[1], pos_on[2]]
        g_omega[pos_on[1], pos_on[2]] <- tmp
      }
      g[lower.tri(g)] <- t(g)[lower.tri(g)]
      g_omega[lower.tri(g_omega)] <- t(g_omega)[lower.tri(g_omega)]
      g <- drop0(g)
      if (pos_def(g_omega)) {
        data <-
          mvrnorm(
            n = params$n[k],
            mu = rep(0, params$p),
            Sigma = solve(g_omega)
          )
        res[[k]] <- list("theta" = g, "omega" = g_omega, "data" = data)
        break
      }
    }
  }
  return(res)
}

n <- hyperparams$n
p <- hyperparams$p


g <- gen(K, params = hyperparams)
Sigma <- array(rinvwishart(K, diag(K)), dim = c(K, K))
theta_0 <- array(dim = c(p,p,K))
for (i in (1:(p-1))) {
  for (j in ((i+1):p)) {
    theta_0[i,j,] <- mvrnorm(n = 1, rep(-4, K), Sigma)
    theta_0[j,i,] <- theta_0[i,j,]
  }
}

for (i in 1:p) {
  theta_0[i,i,] <- 0
}

y <- lapply(g, function(gg) gg$data)
Omega_0 <- rWishart(K, p, diag(p))
Om0 <- Omega_0
cm <- bmg.ecm(y, Omega_0, theta_0, Sigma, hyperparams)
cm.estep <- bmg.estep(cm$Omega, cm$theta, hyperparams)

# Plotting
num <- 1
xlim <- c(min(cm$Omega[,,num]), max(cm$Omega[,,num]))
ylim <- c(min(g[[num]]$omega), max(g[[num]]$omega))
true_positive_mask <- cm.estep$prob[,,num] > 0.5 & g[[num]]$omega > 1e-4
false_positive_mask <- cm.estep$prob[,,num] > 0.5 & g[[num]]$omega < 1e-4
true_negative_mask <- cm.estep$prob[,,num] < 0.5 & g[[num]]$omega < 1e-4
false_negative_mask <- cm.estep$prob[,,num] < 0.5 & g[[num]]$omega > 1e-4

mask <- true_positive_mask
plot(cm$Omega[,,num][mask], g[[num]]$omega[mask], xlim=xlim, ylim=ylim, col = "red")
mask <- true_negative_mask
points(cm$Omega[,,num][mask], g[[num]]$omega[mask], col = "orange")
mask <- false_positive_mask
points(cm$Omega[,,num][mask], g[[num]]$omega[mask], col = "blue")
mask <- false_negative_mask
points(cm$Omega[,,num][mask], g[[num]]$omega[mask], col = "purple")
legend(x="topleft", legend=c("True +", "True -", "False +", "False -"), fill = c("red", "orange", "blue", "purple"))
