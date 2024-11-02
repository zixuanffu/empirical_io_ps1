rm(list = ls())
pacman::p_load(data.table, fixest, gmm)
# the random coefficient is on the size
dt <- readRDS("Data/carpanel_q3.rds")

var_exo <- c("dpm", "door3", "door4", "door5", "at", "ps", "air", "drv", "wt", "hp2wt", "hp", "euro", "japan", "wb", "size")
var_end <- c("p")


# instr class 1: local differentiation

char_lst <- c("const", "disp", "size", var_exo)
for (x in char_lst) {
    sd_x <- sd(dt[, get(x)])
    for (i in unique(dt$year)) {
        diff_x <- abs(outer(dt[year == i, get(x)], dt[year == i, get(x)], "-"))
        diff_x <- abs(diff_x) < (0.5 * sd_x)
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

# instr class 3: lecture notes
for (x in char_lst) {
    for (i in unique(dt$year)) {
        diff_x <- outer(dt[year == i, get(x)], dt[year == i, get(x)], "-")
        cut1 <- max(abs(diff_x)) / 5
        diff_x <- (abs(diff_x) <= cut1)
        diff_x <- diff_x - diag(diag(diff_x))
        iv_x <- diff_x %*% dt[year == i, get(x)]
        dt[year == i, paste0(x, "_instr") := iv_x]
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
# var_iv_local <- c("size_local", "wb_local", "hp_local", "hp2wt_local", "disp_local")
# var_iv_quad <- c("wb_quad", "hp_quad", "size_quad", "hp2wt_local", "disp_quad")
var_iv_instr <- c("size_instr", "hp_instr", "hp2wt_instr", "dpm_instr", "air_instr", "wb_instr")

# plot
sigma_list <- seq(0, 4, 0.1)
# # local diff iv
# g_list <- c()
# for (i in sigma_list) {
#     g_list <- c(g_list, blp_moment_condition_2(i, dt, var_iv_local))
# }
# pdf("Results/Figures/gmm_obj_iv_local.pdf")
# plot(sigma_list, g_list, type = "l", xlab = "sigma", ylab = "gmm objective function")
# dev.off()
# # nested logit iv
# g_list <- c()
# for (i in sigma_list) {
#     g_list <- c(g_list, blp_moment_condition_2(i, dt, var_iv_nl))
# }
# pdf("Results/Figures/gmm_obj_iv_nl.pdf")
# plot(sigma_list, g_list, type = "l", xlab = "sigma", ylab = "gmm objective function")
# dev.off()
# # blp iv
# g_list <- c()
# for (i in sigma_list) {
#     g_list <- c(g_list, blp_moment_condition_2(i, dt, var_iv_blp))
# }
# pdf("Results/Figures/gmm_obj_iv_blp.pdf")
# plot(sigma_list, g_list, type = "l", xlab = "sigma", ylab = "gmm objective function")
# dev.off()

# # quadratic iv
# g_list <- c()
# for (i in sigma_list) {
#     g_list <- c(g_list, blp_moment_condition_2(i, dt, var_iv_quad))
# }
# pdf("Results/Figures/gmm_obj_iv_quad.pdf")
# plot(sigma_list, g_list, type = "l", xlab = "sigma", ylab = "gmm objective function")
# dev.off()

# instr iv
g_list <- c()
for (i in sigma_list) {
    g_list <- c(g_list, blp_moment_condition_2(i, dt, var_iv_instr))
}
pdf("Results/Figures/gmm_obj_iv_instr.pdf")
plot(sigma_list, g_list, type = "l", xlab = "sigma", ylab = "gmm objective function")
dev.off()


# optim
sigma_opt <- optim(2, blp_moment_condition_2, data = dt, var_iv_new = var_iv_instr)
sigma <- sigma_opt$par
blp_result <- blp_intermediate(sigma, dt, var_iv_instr)
ivreg <- blp_result[[2]]
summary(ivreg, cluster = "modelid")
etable(ivreg,
    cluster = "modelid", fitstat = ~ n + ar2 + ivf1 + ivf1.p + sargan + sargan.p,
    file = "Results/Tables/rand_coef.tex", tex = TRUE, replace = TRUE
)
beta <- ivreg$coefficients
delta <- blp_result[[1]]
mu <- blp_result[[3]]
save(dt, delta, sigma, ivreg, beta, mu, var_exo, var_end, file = "Data/blp_results.rda")


gmm_moment <- function(theta, data, var_iv_new = var_iv_instr) {
    return(blp_moment_condition(theta = theta, data = data, var_iv_new = var_iv_new))
}
sigma_gmm <- gmm(gmm_moment, x = dt, t0 = 1)
sigma <- sigma_gmm$coefficients
sigma_vcov <- sqrt(sigma_gmm$vcov)
blp_result <- blp_intermediate(sigma, dt, var_iv_instr)
ivreg <- blp_result[[2]]
summary(ivreg, cluster = "modelid")
etable(ivreg,
    cluster = "modelid", fitstat = ~ n + ar2 + ivf1 + ivf1.p + sargan + sargan.p,
    file = "Results/Tables/rand_coef_gmm.tex", tex = TRUE, replace = TRUE
)
beta <- ivreg$coefficients
delta <- blp_result[[1]]
mu <- blp_result[[3]]
save(dt, delta, sigma, sigma_vcov, ivreg, beta, mu, var_exo, var_end, file = "Data/blp_results_gmm.rda")
