library(huge)
library(Matrix)
library(MASS)
library(gdata)
library(testit)
library(matrixcalc)
library(LaplacesDemon)
library(foreach)
library(parallel)
set.seed(111)
source("ecm.r")

## INPUTS: y is a 3D array that contains K matrices of size n_k \times p Omega
## is a 3D array that contains K matrices of size p \times p RETURNS: todo!
bmg.ecm <- function(y, Omega, theta, Sigma, params, maxiter = 1000) {
    S <- list()
    K <- dim(Omega)[3]
    if (is.na(K)) {
        K <- 1
        dim(Omega) <- c(dim(Omega)[1], dim(Omega)[2], 1)
        dim(theta) <- c(dim(theta)[1], dim(theta)[2], 1)
        y <- list(y)
    }
    for (k in 1:K) {
        S[[k]] <- t(y[[k]]) %*% y[[k]]
    }

    log_liks <- c()
    # Run ECM
    for (t in 1:maxiter) {
        estep <- bmg.estep(Omega, theta, params)
        updated_vals <- bmg.mstep(Omega, theta, Sigma, S, params)
        Omega <- updated_vals$Omega
        theta <- updated_vals$theta
        Sigma <- updated_vals$Sigma

        # Get posterior likelihood
        log_lik <- bmg.posterior_likelihood(y, Omega, theta, Sigma, estep, params)
        log_liks <- append(log_liks, log_lik)

        if (length(log_liks) >= 2 && abs(log_liks[t] - log_liks[t - 1]) < 0.01) {
            break
        }
    }

    # Get posterior edge inclusion probabilities
    estep <- bmg.estep(Omega, theta, params)

    # TODO: q is not posterior edge inclusion!
    return(list(theta = theta, Omega = Omega, Sigma = Sigma, prob = estep$prob, log_liks = log_liks))
}

bmg.estep <- function(Omega, theta, params) {
    v0 <- params$v0
    v1 <- params$v1

    a1 <- dnorm(Omega, mean = 0, sd = v1, log = T) + pnorm(theta, log.p = T)
    a0 <- dnorm(Omega, mean = 0, sd = v0, log = T) + pnorm(theta, log.p = T, lower.tail = F)
    m <- pmax(a1, a0)
    pdeltais1 <- exp(a1 - m)/(exp(a1 - m) + exp(a0 - m))
    # print(sum(pdeltais1 > 0.5))

    d <- (1 - pdeltais1)/v0^2 + pdeltais1/v1^2

    mr1_log <- dnorm(theta, log = TRUE) - pnorm(theta, log.p = TRUE)
    mr0_log <- dnorm(theta, log = TRUE) - pnorm(theta, log.p = TRUE, lower.tail = F)

    q <- theta - exp(mr0_log) + pdeltais1 * (exp(mr1_log) + exp(mr0_log))
    return(list(d = d, q = q, prob = pdeltais1))
}

bmg.posterior_likelihood <- function(y, Omega, theta, Sigma, estep, params) {
    p <- params$p
    n <- params$n
    K <- dim(Omega)[3]
    lambda <- params$lambda
    Sigma_inv <- solve(Sigma)

    # Sum the terms that depend only on k
    sum_k <- 0
    for (k in 1:K) {
        S_k <- t(y[[k]]) %*% y[[k]]
        Omega_k <- Omega[, , k]
        lambda_k <- lambda[k]
        this_k <- 0.5 * (n[k] * log(det(Omega_k)) - sum(diag(S_k %*% Omega_k)) - lambda_k *
            sum(diag(Omega_k)))
        sum_k <- sum_k + this_k
    }

    # Sum the terms that depend only on i,j
    theta_bar <- params$theta_bar
    sum_ij <- 0
    for (i in 1:(p - 1)) {
        for (j in (i + 1):p) {
            this_ij <- t(theta[i, j, ] - theta_bar) %*% Sigma_inv %*% (theta[i, j,
                ] - theta_bar)
            sum_ij <- sum_ij + this_ij
        }
    }

    # Sum the terms that depend on both i,j, and k
    omega_bit <- -0.5 * sum(Omega^2 * estep$d)
    theta_bit <- sum(theta * estep$q) - 0.5 * sum(theta^2)
    sum_ijk <- omega_bit + theta_bit

    # Sigma inverse Wishart bit
    sigma_bit <- -0.5 * (params$nu + K + 1) * log(det(Sigma)) - 0.5 * p * (p - 1)/4 *
        log(det(Sigma)) - 0.5 * sum(diag(params$Psi %*% Sigma_inv))

    sum_k + sum_ij + sum_ijk + sigma_bit
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
    n <- params$n
    K <- dim(Omega)[3]
    lambda <- params$lambda
    theta_bar <- params$theta_bar
    Sigma_inv <- solve(Sigma)
    new_theta <- array(dim = c(p, p, K))

    for (k in 1:K) {
        Omega_k <- Omega[, , k]
        for (i in 1:p) {
            for (j in 1:p) {
                new_theta[i, j, k] <- (q[i, j, k] + Sigma_inv[k, ] %*% rep(theta_bar,
                  K) - Sigma_inv[k, -k] %*% theta[i, j, -k])/(1 + Sigma_inv[k, k])
            }
        }
        diag(new_theta[, , k]) <- 0

        # Compute new omega
        d_k <- d[, , k]
        S_k <- S[[k]]
        lambda_k <- lambda[k]
        pp2 <- p * (p - 1)/2
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
            Omega_k[p, p] <- v + n[k]/(lambda_k + S_k_22)

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
    # Update sigma
    # Theta <- new_theta
    # mask <- apply(Theta, 3, upper.tri)
    # dim(Theta) <- c(p^2, K)
    # Theta <- matrix(Theta[mask], ncol = K)
    # Psi <- params$Psi
    # nu <- params$nu
    # TT <- t(Theta - theta_bar) %*% (Theta - theta_bar)
    # Sigma <- (TT + Psi)/(pp2 + nu + K + 1)

    return(list(theta = new_theta, Omega = Omega, Sigma = Sigma))
}

bmg.sel_v0 <- function(graph) {
    v0s <- seq(1e-3, 0.24, length.out=30)
    K <- length(graph)
    params <- params(1)
    init_val <- init_mg(params, 1)
    sapply(1:K, function(k) {
      g <- graph[[k]]
      vals <- sapply(v0s, function(v0) {
          params$v0 <- v0
          fit <- bmg.ecm(g$data, init_val$Omega_0[,,1], init_val$theta_0[,,1],
                         init_val$Sigma, params)
          # Using extended BIC
          om <- fit$Omega
          om[fit$prob < 0.5] <- 0.0
          diag(om[,,1]) <- diag(fit$Omega[,,1])
          estep <- bmg.estep(om, fit$theta, params)
          om <- om[,,1]
          s <- t(g$data) %*% g$data
          l <- params$n/2 * log(det(om)) - 1/2 * sum(diag(s %*% om))
          num_edges <- sum(fit$prob > 0.5)
          gamma <- 0.0
          
          # AIC
          val <- -2 * l + num_edges * 2 + 4 * num_edges * gamma * log(params$p)
          # print(c(v0, val))
          val
      })
      v0s[which.min(vals)]
    })
}

pos_def <- function(omega, eps = 0.0001) {
    eigenvalues <- eigen(omega)$values
    for (v in eigenvalues) {
        if (v < eps) {
            return(FALSE)
        }
    }
    TRUE
}

gen <- function(K = 1, prob = 0.1, graph_kind = "random", params) {
    res <- list()
    p <- params$p
    n <- params$n
    graph <- huge.generator(n = n[1], d = p, graph = graph_kind, verbose = F)
    num_edges <- sum(graph$theta == TRUE)/2
    for (k in 1:K) {
        while (TRUE) {
            g <- graph$theta
            g_omega <- graph$omega
            num_swap <- rbinom(1, size = num_edges, prob = prob)
            if (num_swap != 0) {
                for (s in 1:num_swap) {
                  # assert(sum(g) == sum(graph$theta))
                  edges_mask <- g & upper.tri(g)
                  non_edges_mask <- !g & upper.tri(g)
                  edges <- which(edges_mask, arr.ind = TRUE)
                  non_edges <- which(non_edges_mask, arr.ind = TRUE)

                  # Pick random edges
                  swapping_off <- sample(1:nrow(edges), 1)

                  # Pick random non-edges
                  swapping_on <- sample(1:nrow(non_edges), 1)

                  # Execute swap
                  pos_off <- edges[swapping_off, ]
                  pos_on <- non_edges[swapping_on, ]

                  g[pos_off[1], pos_off[2]] <- 0
                  g[pos_on[1], pos_on[2]] <- 1

                  g[pos_off[2], pos_off[1]] <- 0
                  g[pos_on[2], pos_on[1]] <- 1

                  tmp <- g_omega[pos_off[1], pos_off[2]]
                  g_omega[pos_off[1], pos_off[2]] <- g_omega[pos_on[1], pos_on[2]]
                  g_omega[pos_on[1], pos_on[2]] <- tmp
                }
            }
            g_omega[lower.tri(g_omega)] <- t(g_omega)[lower.tri(g_omega)]
            g <- drop0(g)
            if (pos_def(g_omega)) {
                data <- mvrnorm(n = params$n[k], mu = rep(0, params$p), Sigma = solve(g_omega))
                res[[k]] <- list(theta = g, omega = g_omega, data = data)
                break
            }
        }
    }
    return(res)
}

### MAIN ###
bmg.main <- function(prob=0.1) {
    K <- 5
    params <- params(K)
    n <- params$n
    p <- params$p
    g <- gen(K, params = params, prob = prob, graph_kind = "random")
  
    v0s <- bmg.sel_v0(g)
    params$v0 <- v0s
    print(c("v0:", params$v0))
      
    Sigma <- array(rep(0.9, K^2), dim = c(K, K)) + 0.1 * diag(K)
    theta_0 <- array(dim = c(p, p, K))
    for (i in (1:(p - 1))) {
        for (j in ((i + 1):p)) {
            theta_0[i, j, ] <- mvrnorm(n = 1, rep(params$theta_bar, K), Sigma)
            theta_0[j, i, ] <- theta_0[i, j, ]
        }
    }

    for (i in 1:p) {
        theta_0[i, i, ] <- 0
    }

    y <- lapply(g, function(gg) gg$data)
    Omega_0 <- rWishart(K, p, diag(p))
    cm <- bmg.ecm(y, Omega_0, theta_0, Sigma, params)
    
    sg <- lapply(1:K, function(k) {
        params <- params(1)
        params$v0 <- v0s[k]
        sgk <- bmg.ecm(y[[k]], Omega_0[, , k], theta_0[,,k], Sigma = array(0.0021, dim=c(1,1)), params)
        against.plot(cm, sgk, g, num = k)
    })

    bmg.plot(cm, g, params)

    print(sapply(1:K, function(i) f1(g[[i]], cm$prob[, , i])))

    return(cm)
}

### EXPERIMENTS ###

params <- function(K) {
    list(v0 = rep(0.05, K),
         v1 = rep(10, K),
         lambda = rep(2, K),
         n = rep(50, K),
         p = 20,
         Psi = diag(K) + matrix(1, ncol=K, nrow=K),
         nu = K,
         theta_bar = -1
    )
}

init_mg <- function(params, K) {
    p <- params$p
    Sigma <- array(0.01, dim=c(K,K)) + diag(K) * 0.001
    theta_0 <- array(dim = c(p, p, K))
    for (i in (1:(p - 1))) {
        for (j in ((i + 1):p)) {
            theta_0[i, j, ] <- mvrnorm(n = 1, rep(params$theta_bar, K), Sigma)
            theta_0[j, i, ] <- theta_0[i, j, ]
        }
    }

    for (i in 1:p) {
        theta_0[i, i, ] <- 0
    }

    Omega_0 <- rWishart(K, p, diag(p))

    return(list(Omega_0 = Omega_0, Sigma = Sigma, theta_0 = theta_0))
}

exp1 <- function(which = "scale-free") {
    K <- 5
    cl <- makeCluster(detectCores() - 1)
    clusterExport(cl, c("params", "huge", "huge.select", "f1", "bmg.ecm", "init_mg",
        "gen", "huge.generator", "which", "drop0", "pos_def", "mvrnorm", "bmg.estep",
        "bmg.mstep", "bmg.posterior_likelihood", "bmg.sel_v0"))
    res <- parLapply(cl, 1:10, function(i) {
        params <- params(K)
        g <- gen(K = K, prob = 0.1, graph_kind = which, params = params)
        # params$v0 <- bmg.sel_v0(g)
        mb <- mean(sapply(g, function(gg) {
            h <- huge.select(huge(gg$data, verbose = F), verbose = F)
            max(sapply(h$path, function(p) f1(gg, p)))
        }))

        mg5 <- {
            init_val <- init_mg(params, K)
            y <- lapply(g, function(gg) gg$data)
            fit <- bmg.ecm(y, init_val$Omega, init_val$theta, init_val$Sigma, params)
            mean(sapply(1:length(g), function(i) f1(g[[i]], fit$prob[, , i])))
        }

        mg1 <- mean(sapply(g, function(gg) {
            K <- 1
            params <- params(K)
            # params$v0 <- bmg.sel_v0(list(gg))
            init_val <- init_mg(params, K)
            fit <- bmg.ecm(gg$data, init_val$Omega[, , 1], init_val$theta[, , 1],
                init_val$Sigma, params)
            f1(gg, fit$prob[, , 1])
        }))
        return(list(mb = mb, mg1 = mg1, mg5 = mg5))
    })
    print(c(mean(sapply(res, function(e) e$mb)), mean(sapply(res, function(e) e$mg1)),
        mean(sapply(res, function(e) e$mg5))))
    stopCluster(cl)
    res
}

exp2 <- function(which) {
     # Create multi-threading cluster
     cl <- makeCluster(detectCores() - 1)
     clusterExport(cl, c("init_mg", "mvrnorm", "bmg.ecm", "bmg.estep", "bmg.mstep",
     "f1", "bmg.posterior_likelihood"))
    res <- lapply(c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), function(prob) {
        print(prob)
        # Generate K-graphs `num_samples` times with that prob setting
        K <- 5
        params <- params(K)
        num_sample <- 20
        ggg <- lapply(1:num_sample, function(i) {
            gen(K = K, graph_kind = which, prob = prob, params = params)
        })

        # Do multi-graph inference on each of those
        fit_info <- lapply(ggg, function(g) {
            K <- length(g)
            init_val <- init_mg(params = params, K = K)
            y <- lapply(g, function(gg) gg$data)
            fit <- bmg.ecm(y, init_val$Omega_0, init_val$theta_0, init_val$Sigma,
                params)
            f1s <- sapply(1:K, function(i) f1(g[[i]], fit$prob[, , i]))
            list(f1s = f1s, Sigma = fit$Sigma)
        })
        f1s <- sapply(fit_info, function(f) f$f1s)
        # Return average f1 score with standard error
        means <- rowMeans(f1s)
        sds <- sqrt(rowMeans((f1s - rowMeans(f1s))^2))
        ses <- qnorm(0.975, mean=0, sd=num_sample/(num_sample-1) * sds / sqrt(num_sample))
        list(means = means, ses = ses, Sigmas = lapply(fit_info, function(f) f$Sigma))
    })
    stopCluster(cl)
    res
}


exp3 <- function(graph_kind) {
    num_rep <- 20
    cl <- makeCluster(detectCores() - 1)
    clusterExport(cl, c("init_mg", "mvrnorm", "bmg.ecm", "bmg.estep", "bmg.mstep",
        "f1", "bmg.posterior_likelihood", "gen", "huge.generator", "params", "which",
        "drop0", "pos_def", "mvrnorm", "triu"))
    parLapply(cl, 1:num_rep, function(i) {
        print(i)
        # Generate the graph now
        g <- gen(K = 10, prob = 0.1, graph_kind = graph_kind, params = params(10))
        sapply(c(1, 2, 5, 10), function(K) {
            pars <- params(K)
            init_val <- init_mg(pars, K)
            y <- lapply(g[1:K], function(gg) gg$data)
            fit <- bmg.ecm(y, init_val$Omega_0, init_val$theta_0, init_val$Sigma,
                pars)
            mean(sapply(1:K, function(i) f1(g[[i]], fit$prob[, , i])))
        })
    })
}

ecm.avg_f1 <- function() {
    sapply(seq(0.001, 0.12, length.out = 50), function(v0) {
        res <- sapply(1:15, function(i) {
            g <- huge.generator(n = 50, d = 20, graph = "random", verbose = F)
            init_val <- init_mg(params(1), 1)
            fit <- ecm(g$data, init_val$Omega[, , 1], v0 = v0, v1 = 10)
            f1(g, fit$prob)
        })
        mean(res)
    })
}

### Plotting ###
bmg.plot <- function(cm, graph, params, num = 1, xlab = "Estimated Omega entries",
    ylab = "True Omega entries") {
    cm.estep <- bmg.estep(cm$Omega, cm$theta, params)
    Omega_hat <- cm$Omega[, , num]
    # Get masks for TP, FP etc
    true_positive_mask <- upper.tri(Omega_hat) & cm.estep$prob[, , num] > 0.5 & graph[[num]]$omega >
        0.05
    false_positive_mask <- upper.tri(Omega_hat) & cm.estep$prob[, , num] > 0.5 &
        graph[[num]]$omega < 0.05
    true_negative_mask <- upper.tri(Omega_hat) & cm.estep$prob[, , num] < 0.5 & graph[[num]]$omega <
        0.05
    false_negative_mask <- upper.tri(Omega_hat) & cm.estep$prob[, , num] < 0.5 &
        graph[[num]]$omega > 0.05

    # Jitter the values
    p <- params$p
    Omega_hat <- Omega_hat + matrix(rnorm(p^2, sd = 0.001), ncol = p)
    Omega_true <- graph[[num]]$omega + matrix(rnorm(p^2, sd = 0.001), ncol = p)

    # Get xlim and ylim
    entries_hat <- Omega_hat[upper.tri(Omega_hat)]
    entries_true <- Omega_true[upper.tri(Omega_hat)]
    xlim <- c(min(entries_hat), max(entries_hat))
    ylim <- c(min(entries_true), max(entries_true))

    # Do the plotting
    mask <- true_positive_mask
    plot(cm$Omega[, , num][mask], Omega_true[mask], xlim = xlim, ylim = ylim, col = "red",
        xlab = xlab, ylab = ylab)
    mask <- true_negative_mask
    points(cm$Omega[, , num][mask], Omega_true[mask], col = "orange")
    mask <- false_positive_mask
    points(cm$Omega[, , num][mask], Omega_true[mask], col = "blue")
    mask <- false_negative_mask
    points(cm$Omega[, , num][mask], Omega_true[mask], col = "purple")
    legend(x = "topleft", legend = c("True +", "True -", "False +", "False -"), fill = c("red",
        "orange", "blue", "purple"))
}

ecm.plot <- function(omega, ecm, xlab = "Estimated Omega entries", ylab = "True Omega entries") {
    omega_hat <- ecm$omega
    # Get masks
    true_positive_mask <- upper.tri(omega_hat) & ecm$prob > 0.5 & omega > 0.0001
    false_positive_mask <- upper.tri(omega_hat) & ecm$prob > 0.5 & omega < 0.0001
    true_negative_mask <- upper.tri(omega_hat) & ecm$prob < 0.5 & omega < 0.0001
    false_negative_mask <- upper.tri(omega_hat) & ecm$prob < 0.5 & omega > 0.0001

    # Jitter the values
    p <- ncol(omega)
    omega_hat <- omega_hat + matrix(rnorm(p^2, sd = 0.001), ncol = p)
    omega_true <- omega + matrix(rnorm(p^2, sd = 0.001), ncol = p)

    entries <- omega_hat[upper.tri(omega_hat)]
    true_entries <- omega_true[upper.tri(omega)]
    xlim <- c(min(entries), max(entries))
    ylim <- c(min(true_entries), max(true_entries))
    mask <- true_positive_mask
    plot(omega_hat[mask], omega_true[mask], xlim = xlim, ylim = ylim, col = "red",
        xlab = xlab, ylab = ylab)
    mask <- true_negative_mask
    points(omega_hat[mask], omega_true[mask], col = "orange")
    mask <- false_positive_mask
    points(omega_hat[mask], omega_true[mask], col = "blue")
    mask <- false_negative_mask
    points(omega_hat[mask], omega_true[mask], col = "purple")
    legend(x = "topleft", legend = c("True +", "True -", "False +", "False -"), fill = c("red",
        "orange", "blue", "purple"))
}

against.plot <- function(mg, sg, graph, num = 1) {
    omg <- mg$Omega[, , num]
    osg <- sg$Omega[, , 1]

    pos <- upper.tri(omg) & graph[[num]]$omega > 0.001
    neg <- upper.tri(omg) & !graph[[num]]$omega > 0.001

    xlim <- c(min(omg[upper.tri(omg)]), max(omg[upper.tri(omg)]))
    ylim <- c(min(osg[upper.tri(osg)]), max(osg[upper.tri(osg)]))

    plot(omg[neg], osg[neg], xlim = xlim, ylim = ylim, col = "blue", xlab = "Multi-graph Omega entries",
        ylab = "Single-graph Omega entries")
    points(omg[pos], osg[pos], xlim = xlim, ylim = ylim, col = "red")
    title(c("Single-graph vs Multi-graph;", num))
}


### Scoring ###
f1 <- function(graph, prob) {
    truth <- abs(graph$omega) > 0.1
    declared <- prob > 0.5
    TP <- sum(truth & declared)
    FP <- sum(!truth & declared)
    FN <- sum(truth & !declared)

    precision <- TP/(TP + FP)
    recall <- TP/(TP + FN)

    if (is.na(precision) || is.na(recall)) {
        return(0)
    }

    2 * precision * recall/(precision + recall)
}
