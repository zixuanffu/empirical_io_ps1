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

source("Code/Q4_function.R")
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

# # examine the variance of the characteristics over the year
# market_var <- data.table("modelid" = unique(dt$modelid))
# for (x in c("disp", "size", var_exo)) {
#     for (y in unique(dt$modelid)) {
#         market_var[modelid == y, paste0(x, "_sd") := sd(dt[modelid == y, get(x)])]
#     }
# }

var_iv_local <- c(var_iv_nl, "size_local", "wb_local", "hp_local", "hp2wt_local")
var_iv_quad <- c(var_iv_nl, "wb_quad", "hp_quad", "size_quad", "hp2wt_local")

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

# quadratic iv
g_list <- c()
for (i in sigma_list) {
    g_list <- c(g_list, blp_moment_condition_2(i, dt, var_iv_quad))
}
pdf("Results/Figures/gmm_obj_iv_quad.pdf")
plot(sigma_list, g_list, type = "l", xlab = "sigma", ylab = "gmm objective function")
dev.off()

# gmm estimation
library(gmm)
?gmm
gmm_moment <- function(theta, data, var_iv_new = var_iv_local) {
    return(blp_moment_condition(theta = theta, data = data, var_iv_new = var_iv_local))
}

sigma_opt <- optim(2, blp_moment_condition_2, data = dt, var_iv_new = var_iv_local, method = "BFGS")
sigma <- sigma_opt$par
blp_result <- blp_intermediate(sigma_gmm$coefficients, dt, var_iv_local)
ivreg <- blp_result[[2]]
summary(ivreg, cluster = "modelid")
etable(ivreg,
    cluster = "modelid", fitstat = ~ n + ar2 + ivf1 + ivf1.p + sargan + sargan.p,
    file = "Results/Tables/rand_coef.tex", tex = TRUE, replace = TRUE
)
fitstat(ivreg, ~ n + ar2 + ivf1 + ivf1.p + sargan + sargan.p)
beta <- ivreg$coefficients
delta <- blp_result[[1]]
mu <- blp_result[[3]]
save(dt, delta, sigma, ivreg, beta, mu, var_exo, var_end, file = "Data/blp_results.rda")

sigma_gmm <- gmm(gmm_moment, x = dt, t0 = 1)
sigma_gmm_v <- sigma_gmm$coefficients
sigma_gmm_vcov <- sigma_gmm$vcov
blp_result <- blp_intermediate(sigma_gmm_v, dt, var_iv_local)
ivreg <- blp_result[[2]]
summary(ivreg, cluster = "modelid")
etable(ivreg,
    cluster = "modelid", fitstat = ~ n + ar2 + ivf1 + ivf1.p + sargan + sargan.p,
    file = "Results/Tables/rand_coef_gmm.tex", tex = TRUE, replace = TRUE
)
beta <- ivreg$coefficients
delta <- blp_result[[1]]
mu <- blp_result[[3]]
save(dt, delta, sigma, ivreg, beta, mu, var_exo, var_end, file = "Data/blp_results_gmm.rda")
