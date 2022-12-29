library("Matrix")
library("huge")
library("testit")
library("gdata")
library("matrixcalc")
# set.seed(111)

e_step <- function(omega, pi, v0, v1) {
  a <- dnorm(omega, 0, v1, log = TRUE) + log(pi)
  b <- dnorm(omega, 0, v0, log = TRUE) + log(1 - pi)
  m <- pmax(a, b)
  p <- exp(a - m) / (exp(a - m) + exp(b - m))
  diag(p) <- 0
  d <- (1 - p) / (v0^2) + p / (v1^2)
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
      omega[c(t, p), ] <- omega[c(p, t), ]
      estep$d[c(t, p), ] <- estep$d[c(p, t), ]
      s[c(t, p), ] <- s[c(p, t), ]

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

      # Update exp part
      v <- t(omega_12_lp1) %*% omega_11_inv %*% omega_12_lp1
      omega[p, p] <- v + n / (lambda + s_22)

      # Swap back everything
      omega[, c(t, p)] <- omega[, c(p, t)]
      estep$d[, c(t, p)] <- estep$d[, c(p, t)]
      s[, c(t, p)] <- s[, c(p, t)]
      omega[c(t, p), ] <- omega[c(p, t), ]
      estep$d[c(t, p), ] <- estep$d[c(p, t), ]
      s[c(t, p), ] <- s[c(p, t), ]
    }
    return(list("omega" = omega, "pi" = new_pi))
  }

log_lik <- function(alpha, beta, lambda, pi, omega, n, estep, s) {
  eps <- .Machine$double.eps^0.5
  p <- nrow(omega)
  sns_part <- -sum(upperTriangle(omega^2 * estep$d)) / 2
  exp_part <- -lambda / 2 * sum(diag(omega))
  ber_part <-
    log(pi / (1 - pi) + eps) * sum(upperTriangle(estep$p)) + p * (p - 1) / 2 * log(1 - pi + eps)
  misc_part <-
    (n / 2 * log(det(omega))) - (sum(diag(s %*% omega)) / 2)
  +((alpha - 1) * log(pi + eps)) + ((beta - 1) * log(1 - pi + eps))

  return(sns_part + exp_part + ber_part + misc_part)
}

build_path <- function(p, len = 20) {
  tt <- seq(0.8, 0, length.out = len)
  output <- c()
  for (t in tt) {
    g <- Matrix(p >= t, sparse = TRUE)
    output <- append(output, g)
  }
  return(output)
}

bmg_roc <- function(p, theta, len = 200) {
  path <- build_path(p, len)
  return(my_roc(path, theta, verbose = FALSE))
}

select_v0 <- function(graph, alpha = 1, beta = 1, lambda = 2, pi = 0.5, v1 = 100) {
  v0s <- seq(1e-3, 0.12, length.out = 50)
  # v0s <- exp(seq(log(1e-5), log(0.3), length.out = 150))
  n <- nrow(graph$data)
  p <- ncol(graph$data)
  pp2 <- p * (p - 1) / 2
  vals_v0s <- c()
  x <- graph$data
  s <- t(x) %*% x
  omega <- rinvwishart(p, diag(p))
  # omega <- as(solve(nearPD(s)$mat), "matrix")
  # as(nearPD(solve(cov(x) + ridge * diag(p)))$mat, "matrix")
  aucs <- c()
  f1s <- c()
  tr_aic <- c()
  logdet_aic <- c()
  k_aic <- c()
  fb_corrs <- c()
  omega_diffs <- c()
  last_omega <- NULL
  omegas <- list()
  for (v0 in v0s) {
    print(v0)
    e <- ecm(x, omega, v0, pi, alpha, beta, lambda, v1)
    k <- sum(e$prob > 0.5)
    print(k)
    om <- e$omega
    om[e$prob <= 0.5] <- 0.0
    diag(om) <- diag(e$omega)
    val <-
      sum(diag(s %*% om)) - n * determinant(om, logarithm = TRUE)$modulus[1] + 2 * k
    vals_v0s <- append(vals_v0s, val)
    i <- length(vals_v0s)

    omegas[[i]] <- e$omega

    # Will be plotted
    fb_corrs <-
      append(fb_corrs, log(val * 1 / ((pp2 + 1) * choose(pp2, k))))
    roc <- bmg_roc(e$prob, graph$theta)
    aucs <- append(aucs, roc$AUC)
    f1s <- append(f1s, max(roc$F1))
    tr_aic <- append(tr_aic, sum(diag(s %*% e$omega)))
    logdet_aic <- append(logdet_aic, -n * log(det(e$omega)))
    k_aic <- append(k_aic, k)

    # Build omega_diffs vector
    # if (!is.null(last_omega)) {
    #   l2_diff <- sqrt(sum((e$omega - last_omega)^2))
    #   omega_diffs <- append(omega_diffs, l2_diff)
    # }
    # last_omega <- e$omega


    # Warm-start next iteration
    # omega <- e$omega
  }
  v0_oracle <- v0s[which.max(aucs)]
  # plot(v0s, aucs, type = "l")
  # abline(v = v0_oracle, col = "red")
  #
  # plot(v0s, f1s, type = "l")
  # abline(v = v0_oracle, col = "red")
  #
  # plot(v0s, vals_v0s, type = "l")
  # abline(v = v0_oracle, col = "red")

  # plot(v0s, tr_aic, type = "l", col = 3)
  # abline(v = v0_oracle, col = "red")

  # plot(v0s, logdet_aic, type = "l", col = 4)
  # abline(v = v0_oracle, col = "red")

  plot(v0s, k_aic / p^2, type = "l", col = 5)
  abline(v = v0_oracle, col = "red")

  # plot(v0s, fb_corrs, type = "l", col = 6)
  # abline(v = v0_oracle, col = "red")

  # plot(v0s[1:length(omega_diffs)], omega_diffs, type = "l", col=7)
  # abline(v = v0_oracle, col = "red")

  # for (i in 1:p) {
  #   for (j in 1:p) {
  #     for (k in 1:length(v0s)) {
  #
  #     }
  #   }
  # }

  return(list(
    "v0s" = v0s,
    "oracle" = v0_oracle,
    "aucs" = aucs,
    "ks" = k_aic,
    "f1s" = f1s,
    # "omega_diffs" = omega_diffs,
    "vals" = vals_v0s
  ))
}

ecm <- function(x,
                omega,
                v0 = 0.05,
                pi = 0.5,
                a = 1,
                b = 1,
                lambda = 1,
                v1 = 100,
                tol = 1e-3,
                maxiter = 1000) {
  p <- ncol(x)
  s <- t(x) %*% x

  # omega <- prior_sns(p, lambda, v0, v1, pi)
  # omega <- diag(p)

  log_liks <- c()
  last_lik <- NULL
  for (t in 1:maxiter) {
    assert(
      "omega must be positive definite",
      is.positive.definite(omega)
    )
    estep <- e_step(omega, pi, v0, v1)
    cmstep <-
      cm_step(a, b, lambda, v0, v1, s, pi, omega, n, estep)
    # omega <- cmstep$omega
    omega <- cmstep$omega
    pi <- cmstep$pi
    l <- log_lik(a, b, lambda, pi, omega, n, estep, s)
    log_liks <- append(log_liks, l)

    # Return early if relative increase in likelihood is too low
    if (is.null(last_lik) || abs(last_lik - l) >= tol) {
      last_lik <- l
    } else {
      return(list(
        "omega" = omega,
        "pi" = pi,
        "prob" = estep$p,
        "log_liks" = log_liks
      ))
    }
  }
  return(list(
    "omega" = omega,
    "pi" = pi,
    "prob" = estep$p,
    "log_liks" = log_liks
  ))
}

# Copy-pasted from
# https://github.com/HMJiangGatech/huge/blob/master/R/huge.roc.R
# So that it doesn't do a plot each time
my_roc <- function(path, theta, verbose = TRUE) {
  gcinfo(verbose = FALSE)
  ROC <- list()

  theta <- as.matrix(theta)
  d <- ncol(theta)
  pos.total <- sum(theta != 0)
  neg.total <- d * (d - 1) - pos.total

  if (verbose) {
    cat("Computing F1 scores, false positive rates and true positive rates....")
  }
  ROC$tp <- rep(0, length(path))
  ROC$fp <- rep(0, length(path))
  ROC$F1 <- rep(0, length(path))
  for (r in 1:length(path)) {
    tmp <- as.matrix(path[[r]])
    tp.all <- (theta != 0) * (tmp != 0)
    diag(tp.all) <- 0
    ROC$tp[r] <- sum(tp.all != 0) / pos.total
    fp.all <- (theta == 0) * (tmp != 0)
    diag(fp.all) <- 0
    ROC$fp[r] <- sum(fp.all != 0) / neg.total

    fn <- 1 - ROC$tp[r]
    precision <- ROC$tp[r] / (ROC$tp[r] + ROC$fp[r])
    recall <- ROC$tp[r] / (ROC$tp[r] + fn)
    ROC$F1[r] <- 2 * precision * recall / (precision + recall)
    if (is.na(ROC$F1[r])) {
      ROC$F1[r] <- 0
    }
  }
  if (verbose) {
    cat("done.\n")
  }

  rm(precision, recall, tp.all, fp.all, path, theta, fn)
  gc()

  ord.fp <- order(ROC$fp)

  tmp1 <- ROC$fp[ord.fp]
  tmp2 <- ROC$tp[ord.fp]
  par(mfrow = c(1, 1))
  ROC$AUC <- sum(diff(tmp1) * (tmp2[-1] + tmp2[-length(tmp2)])) / 2

  rm(ord.fp, tmp1, tmp2)
  gc()
  class(ROC) <- "roc"
  return(ROC)
}

plot.roc <- function(roc) {
  ord.fp <- order(roc$fp)
  tmp1 <- roc$fp[ord.fp]
  tmp2 <- roc$tp[ord.fp]
  par(mfrow = c(1, 1))
  plot(
    tmp1,
    tmp2,
    type = "b",
    main = "ROC Curve",
    xlab = "False Postive Rate",
    ylab = "True Postive Rate",
    ylim = c(0, 1)
  )
}

ecm.main <- function() {
  a <- 1
  b <- 9
  
  
  f1s_emgs <- c()
  f1s_huge <- c()
  
  aucs <- list()
  v0_infos <- list()
  
  num_rep <- 1
  n <- 100
  
  for (rep in 1:num_rep) {
    for (d in c(30)) {
      graph <- huge.generator(n = n, d = d, graph = "hub")
      x <- graph$data
      omega <- precision_mat(graph)
      write.csv(graph$data, "x.csv", row.names = FALSE)
      v0_info <- select_v0(graph, b = b)
      cm <- ecm(x = x, omega = omega, v0 = v0_info$oracle, b = b)
      f1s_emgs <-
        append(f1s_emgs, max(bmg_roc(cm$prob, graph$theta)$F1))
      h <- huge(graph$data)
      h <- huge.select(h)
      f1s_huge <-
        append(f1s_huge, max(my_roc(h$path, graph$theta)$F1))
    }
    ecm.plot(graph$omega, cm$omega, v0 = v0_info$oracle, v1 = 100)
    v0_infos[[rep]] <- v0_info
    aucs[[rep]] <- v0_info$aucs
    r <- Reduce("+", aucs) / num_rep
  }
}