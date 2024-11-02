rm(list = ls())
pacman::p_load(data.table, xtable)
source("Code/Q4_function.R")
load("Data/blp_results.rda")

summary(ivreg, cluster = "modelid")
alpha <- ivreg$coefficients["fit_p"]
# well, a bit too small
sj <- share(delta, mu)
sj_true <- dt$mktshr

# solve by optim
marginal_cost_1 <- function(c, data) {
    p <- data$p
    sj <- data$mktshr
    brand <- model.matrix(~ data$firmids - 1)
    ownership <- brand %*% t(brand)
    omega <- (alpha) * (diag(sj) - sj %*% t(sj)) * ownership
    eq <- omega %*% (p - c) + sj
    res <- sum(eq) * 10^6
    return(res)
}
# mc <- optim(rep(0, nrow(dt[year == 1977])), marginal_cost_1, data = dt[year == 1977], control = list(
#     reltol = 1e-8, # Set relative tolerance
#     abstol = 1e-10, # Set absolute tolerance
#     maxit = 10000 # Set maximum number of iterations
# ))

# marginal_cost_1(mc$par, data = dt[year == 1977])

# solve directly
marginal_cost_2 <- function(data) {
    year <- unique(data$year)
    for (i in year) {
        p <- data[year == i]$p
        sj <- data[year == i]$mktshr
        brand <- model.matrix(~ data[year == i]$firmids - 1)
        ownership <- brand %*% t(brand)
        omega <- (alpha) * (diag(sj) - sj %*% t(sj)) * ownership
        omega_inv <- solve(omega)
        data[year == i, `:=`(mc = p + omega_inv %*% sj, mk = -omega_inv %*% sj, elas = omega_inv %*% sj * p)]
    }
    return(data)
}
dt_postmerger <- marginal_cost_2(data = dt)
# merger simulation
dt_postmerger[firmids == 17, firmids := 7]
price <- function(data) {
    year <- unique(data$year)
    for (i in year) {
        c <- data[year == i]$mc
        sj <- data[year == i]$mktshr
        brand <- model.matrix(~ data[year == i]$firmids - 1)
        ownership <- brand %*% t(brand)
        omega <- (alpha) * (diag(sj) - sj %*% t(sj)) * ownership
        omega_inv <- solve(omega)
        data[year == i, `:=`(p_new = c - omega_inv %*% sj)]
    }
    return(data)
}
dt_postmerger <- price(data = dt_postmerger)

dt_post <- dt_postmerger[, lapply(.SD, mean), by = .(year), .SDcols = c("p", "mc", "mk", "elas", "p_new")]

print(xtable(dt_post), type = "latex", floating = FALSE, file = "Results/Tables/post_merger_mean.tex")

# consumer surplus in monetary terms
var_x <- c("const", "p_new", var_exo)
delta_new <- as.matrix(dt_postmerger[, ..var_x]) %*% ivreg$coefficients
dt_postmerger[, mktshr_new := share(as.vector(delta_new), as.matrix(mu))]
dt_postmerger[, shr_0_new := 1 - sum(mktshr_new), by = year]

# pre-merger
dt_cs <- unique(dt_postmerger[, .(year, nb_hh, shr_0, shr_0_new)])
inclu_util <- -log(dt_cs$shr_0)
nb_hh <- dt_cs$nb_hh
cs <- nb_hh * inclu_util / (-alpha)
# post-merger
inclu_util_new <- -log(dt_cs$shr_0_new)
nb_hh <- dt_cs$nb_hh
cs_new <- nb_hh * inclu_util_new / (-alpha)

dt_cs[, cs_pre := cs]
dt_cs[, cs_post := cs_new]
dt_cs[, diff := cs_new - cs]

print(xtable(dt_cs), floating = FALSE, type = "latex", file = "Results/Tables/post_merger_cs.tex")
