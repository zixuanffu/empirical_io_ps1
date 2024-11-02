rm(list = ls())
pacman::p_load(data.table, fixest, gmm)
# the random coefficient is on the size
dt <- readRDS("Data/carpanel_q3.rds")

var_exo <- c("dpm", "door3", "door4", "door5", "at", "ps", "air", "drv", "wt", "hp2wt", "hp", "euro", "japan", "wb")
var_end <- c("p")


# instr class 1: local differentiation

char_lst <- c("const", "disp", "size", var_exo)
for (x in char_lst) {
    sd_x <- sd(dt[, get(x)])
    for (i in unique(dt$year)) {
        diff_x <- abs(outer(dt[year == i, get(x)], dt[year == i, get(x)], "-"))
        diff_x <- abs(diff_x) < 0.5 * sd_x
        iv_x <- colSums(diff_x)
        dt[year == i, paste0(x, "_local") := iv_x - 1] # exclude the own product
    }
}

# instr class 2: quadratic differentiation
for (x in char_lst) {
    for (i in unique(dt$year)) {
        diff_x <- outer(dt[year == i, get(x)], dt[year == i, get(x)], "-")
        diff_x <- diff_x^2
        iv_x <- colSums(diff_x)
        dt[year == i, paste0(x, "_quad") := iv_x]
    }
}

saveRDS(dt, "Data/carpanel_q4.rds")

IV_lst <- c(paste(char_lst, "_local", sep = ""), paste(char_lst, "_quad", sep = ""))
saveRDS(IV_lst, "Data/iv_q4.rds")


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


# load the data
dt <- readRDS("Data/carpanel_q4.rds")

# draw 500 normal random draws to approximate the integral
n <- 500
# draw in the normal(0, 1)
draw <- rnorm(n)
# check the variables
var_x <- c(var_exo, var_end)
var_iv_blp <- c("const_rival", "const_ownothers", "dpm_rival", "dpm_ownothers", "hp2wt_rival", "hp2wt_ownothers", "size_rival", "size_ownothers", "air_rival", "air_ownothers")
var_iv_nl <- c("const_rival", "const_ownothers", "dpm_rival", "dpm_ownothers", "hp2wt_rival", "hp2wt_ownothers", "air_rival", "air_ownothers", "const_rival_g", "hp2wt_rival_g", "dpm_rival_g", "air_rival_g")

# examine the variance of the characteristics over the year
market_var <- data.table("modelid" = unique(dt$modelid))
for (x in c("disp", "size", var_exo)) {
    for (y in unique(dt$modelid)) {
        market_var[modelid == y, paste0(x, "_sd") := sd(dt[modelid == y, get(x)])]
    }
}

var_iv_local <- c(var_iv_nl, "size_local", "disp_local", "wb_local", "hp_local", "hp2wt_local")

# gmm estimation (simultaneously)

blp_moment_condition <- function(theta, data, var_iv_new) {
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
    var_xz <- c(var_exo, var_iv_new)
    Z <- as.matrix(data[, ..var_xz])
    W <- solve(t(Z) %*% Z)
    g <- t(Z) %*% residual
    gmm_obj <- t(g) %*% W %*% g
    return(gmm_obj)
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


# plot
sigma_list <- seq(0, 4, 0.1)
# local diff iv
g_list <- c()
for (i in sigma_list) {
    g_list <- c(g_list, blp_moment_condition_2(i, dt, var_iv_local))
}
pdf("Results/Figures/gmm_obj_iv_local.pdf")
plot(sigma_list, g_list, type = "l", xlab = "sigma", ylab = "gmm objective function")
dev.off()
# nested logit iv
g_list <- c()
for (i in sigma_list) {
    g_list <- c(g_list, blp_moment_condition_2(i, dt, var_iv_nl))
}
pdf("Results/Figures/gmm_obj_iv_nl.pdf")
plot(sigma_list, g_list, type = "l", xlab = "sigma", ylab = "gmm objective function")
dev.off()
# blp iv
g_list <- c()
for (i in sigma_list) {
    g_list <- c(g_list, blp_moment_condition_2(i, dt, var_iv_blp))
}
pdf("Results/Figures/gmm_obj_iv_blp.pdf")
plot(sigma_list, g_list, type = "l", xlab = "sigma", ylab = "gmm objective function")
dev.off()



# to store intermediate results
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

# gmm estimation
sigma_opt <- optim(2, blp_moment_condition_2, data = dt, var_iv_new = var_iv_local, method = "BFGS")
sigma <- sigma_opt$par
blp_result <- blp_intermediate(sigma, dt, var_iv_nl)
ivreg <- blp_result[[2]]
etable(ivreg,
    cluster = "modelid", fitstat = ~ n + ar2 + ivf1 + ivf1.p + sargan + sargan.p,
    file = "Results/Tables/rand_coef.tex", tex = TRUE, replace = TRUE
)
summary(ivreg, cluster = "modelid")
fitstat(ivreg, ~ n + ar2 + ivf1 + ivf1.p + sargan + sargan.p)
beta <- ivreg$coefficients
delta <- blp_result[[1]]
mu <- blp_result[[3]]
save(dt, delta, sigma, ivreg, beta, mu, file = "Data/blp_results.rda")
