# the main part
# given delta_j and simulated draws, calculate the market share
share <- function(delta, mu) {
    numer <- exp(delta + mu)
    denom <- 1 + colSums(numer) # for each column sum up the rows
    sij <- numer / denom
    sj <- rowSums(sij) / 200 # for each row sum up the columns
    return(sj)
}

# given the market share sj, back out delta_j, and the number of iterations
blp_contraction <- function(sj, mu, delta_init) {
    J <- length(delta_init) # get the number of products
    delta_before <- rep(0, J)
    delta_next <- delta_init
    iter <- 0
    while (max(abs(delta_next - delta_before)) > 1e-7 && iter < 500) {
        delta_before <- delta_next
        iter <- iter + 1
        delta_next <- delta_before + log(sj) - log(share(delta_before, mu))
    }
    delta <- delta_next
    flag_cv <- (iter == 500 || max(is.na(delta)) == 1)
    return(list(delta, flag_cv, iter))
}

blp_moment_condition <- function(theta, data, var_iv_new) {
    delta_new <- c()
    sigma <- theta
    # flag_cv_new <- c()
    # iter_new <- c()
    for (i in unique(data$year)) {
        J <- data[year == i, .N]
        sj <- data[year == i, mktshr]
        mu <- sigma * outer(data[year == i, size], draw)
        delta_init <- data[year == i, log(mktshr) - log(shr_0)]
        contraction_result <- blp_contraction(sj, mu, delta_init)
        delta_market <- contraction_result[[1]]
        delta_new <- c(delta_new, delta_market)
        # flag_cv_new <- c(flag_cv_new, contraction_result[[2]])
        # iter_new <- c(iter_new, contraction_result[[3]])
    }
    data$y <- delta_new
    iv_formula <- as.formula(paste("y ~", paste(var_exo, collapse = " + "), "|", paste(var_end, collapse = " + "), "~", paste(var_iv_new, collapse = " + ")))
    residual <- feols(iv_formula, data = data)$residuals
    var_xz <- c(var_exo, var_iv_new)
    Z <- as.matrix(data[, ..var_xz])
    g <- residual * Z
    return(g)
}

# gmm estimation (concentrate out beta)
blp_moment_condition_2 <- function(theta, data, var_iv_new) {
    delta_new <- c()
    sigma <- theta
    # flag_cv_new <- c()
    # iter_new <- c()
    for (i in unique(data$year)) {
        J <- data[year == i, .N]
        sj <- data[year == i, mktshr]
        mu <- sigma * outer(data[year == i, size], draw)
        delta_init <- data[year == i, log(mktshr) - log(shr_0)]
        contraction_result <- blp_contraction(sj, mu, delta_init)
        delta_market <- contraction_result[[1]]
        delta_new <- c(delta_new, delta_market)
        # flag_cv_new <- c(flag_cv_new, contraction_result[[2]])
        # iter_new <- c(iter_new, contraction_result[[3]])
    }
    data$y <- delta_new
    iv_formula <- as.formula(paste("y ~", paste(var_exo, collapse = " + "), "|", paste(var_end, collapse = " + "), "~", paste(var_iv_new, collapse = " + ")))
    residual <- feols(iv_formula, data = data)$residuals
    var_xz <- c(var_exo, var_iv_new)
    Z <- as.matrix(data[, ..var_xz])
    W <- solve(t(Z) %*% Z)
    g <- t(Z) %*% residual
    gmm_obj <- t(g) %*% W %*% g
    return(gmm_obj)
}

blp_intermediate <- function(theta, data, var_iv_new) {
    delta_new <- c()
    sigma <- theta
    mu_new <- c()
    # flag_cv_new <- c()
    # iter_new <- c()
    for (i in unique(data$year)) {
        J <- data[year == i, .N]
        sj <- data[year == i, mktshr]
        mu <- sigma * outer(data[year == i, size], draw)
        delta_init <- data[year == i, log(mktshr) - log(shr_0)]
        contraction_result <- blp_contraction(sj, mu, delta_init)
        delta_market <- contraction_result[[1]]
        delta_new <- c(delta_new, delta_market)
        mu_new <- rbind(mu_new, mu)
        # flag_cv_new <- c(flag_cv_new, contraction_result[[2]])
        # iter_new <- c(iter_new, contraction_result[[3]])
    }
    data$y <- delta_new
    iv_formula <- as.formula(paste("y ~ ", paste(var_exo, collapse = " + "), "|", paste(var_end, collapse = " + "), "~", paste(var_iv_new, collapse = " + ")))
    reg_iv <- feols(iv_formula, data = data, cluster = "modelid")
    return(list(delta_new, reg_iv, mu_new))
}
