rm(list = ls())
pacman::p_load(data.table, AER, gmm)
# the random coefficient is on the size
dt <- readRDS("Data/carpanel_q3.rds")
# here we need to create dummies for the door
dr_dummy <- model.matrix(~ dt$dr - 1) # -1 removes the intercept, so you get dummies for all categories
colnames(dr_dummy) <- c("door2", "door3", "door4", "door5")
dt <- cbind(dt, dr_dummy[, -1]) # remove the first column, which is the base category

var_exo <- c("dpm", "door3", "door4", "door5", "at", "ps", "air", "drv", "wt", "hp2wt", "hp", "euro", "japan", "wb")
var_end <- c("p", "size")
var_iv <- readRDS("Data/IV_lst_nested.rds")
# recast the data type
var_x <- c(var_exo, var_end)
dt[, (var_x) := lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols = var_x]

# GH2020 (modified) instruments for variable with a random coefficient, here size
# for each market, calculate the difference matrix
dt_new <- data.table()
for (i in unique(dt$year)) {
    dt_year <- dt[year == i]
    size <- dt_year$size
    diff <- abs(outer(size, size, "-"))
    diff_range <- range(diff)
    # for each percentile
    for (j in 1:4) {
        dt_year[, paste0("iv", j) := (diff <= quantile(diff_range, j / 5)) %*% size]
    }
    dt_new <- rbind(dt_new, dt_year)
}
dt <- cbind(dt, dt_new[, .(iv1, iv2, iv3, iv4)])
saveRDS(dt, "Data/carpanel_q4.rds")
var_iv <- c(var_iv, paste0("iv", 1:4)) # update the iv list

# the main part
# given delta_j and simulated draws, calculate the market share
share <- function(delta, mu) {
    numer <- exp(delta + mu)
    denom <- 1 + colSums(numer)
    sij <- numer / denom
    sj <- rowSums(sij) / 200
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

# given the delta_j, calculate the residual from gmm estimation with IV
# moment_condition <- function(beta, data) {
#     var_x <- c(var_exo, var_end)
#     residual <- data$y - data[, ..var_x] %*% beta
#     g <- c()
#     for (i in 1:4) {
#         g <- cbind(g, residual * dt[, paste0("iv", i)])
#     }
#     return(g)
# }


blp_moment_condition <- function(theta, data) {
    delta_new <- c()
    sigma <- theta[1]
    beta <- theta[2:length(theta)]
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
    residual <- as.matrix(data$y) - as.matrix(data[, ..var_x]) %*% beta
    g <- c()
    for (i in var_iv) {
        g <- cbind(g, residual * data[, ..i])
    }
    return(g)
}

# draw 200 normal random draws to approximate the integral
n <- 200
# draw in the normal(0, 1)
draw <- rnorm(n)
# check the variables
var_x
var_iv

# gmm estimation (simultaneously)
reggmm <- gmm(blp_moment_condition, x = dt, t0 = c(1, rep(1, length(var_x))))
summary(reggmm)

blp_moment_condition_2 <- function(theta, data) {
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
    iv_formula <- as.formula(paste("y ~", paste(c(var_exo, var_end), collapse = " + "), "|", paste(c(var_exo, var_iv), collapse = " + ")))
    residual <- ivreg(iv_formula, data = data)$residuals # use iv regression or gmm regression whatever you like
    g <- c()
    for (i in var_iv) {
        g <- cbind(g, residual * data[, ..i])
    }
    return(g)
}
# gmm estimation (concentrate out beta)
reggmm2 <- gmm(blp_moment_condition_2, x = dt, t0 = 1, optfct = "optim", method = "BFGS")
summary(reggmm2)
sigma <- reggmm2$coefficients
blp_intermediate <- function(theta, data) {
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
    iv_formula <- as.formula(paste("y ~", paste(c(var_exo, var_end), collapse = " + "), "|", paste(c(var_exo, var_iv), collapse = " + ")))
    iv_reg <- ivreg(iv_formula, data = data) # use iv regression or gmm regression whatever you like
    return(list(delta_new, iv_reg))
}
blp_result <- blp_intermediate(sigma, dt)
ivreg<-blp_result[[2]]
beta <- ivreg$coefficients
delta <- blp_result[[1]]
save(dt,reggmm2, ivreg, sigma, beta, delta, file = "Data/blp_results.rda")
