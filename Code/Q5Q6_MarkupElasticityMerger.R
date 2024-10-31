rm(list = ls())
pacman::p_load(data.table)
load("Data/blp_results.rda")

summary(ivreg, cluster = "modelid")
alpha <- ivreg$coefficients["fit_p"]
alpha <- -0.3 # this should be the right value but how???
mu <- sigma * outer(dt[, size], draw)
sj <- share(delta, mu)
sj_true <- dt$mktshr
# calculate the marginal cost
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
mc <- optim(rep(0, nrow(dt[year == 1977])), marginal_cost_1, data = dt[year == 1977], control = list(
    reltol = 1e-8, # Set relative tolerance
    abstol = 1e-10, # Set absolute tolerance
    maxit = 10000 # Set maximum number of iterations
))

marginal_cost_1(mc$par, data = dt[year == 1977])

# solve directly
marginal_cost_2 <- function(data) {
    p <- data$p
    sj <- data$mktshr
    brand <- model.matrix(~ data$firmids - 1)
    ownership <- brand %*% t(brand)
    omega <- (alpha) * (diag(sj) - sj %*% t(sj)) * ownership
    omega_inv <- solve(omega)
    mc <- p + omega_inv %*% sj
    return(mc)
}
mc <- marginal_cost_2(data = dt[year == 1977])

# markup
markup <- function(data) {
    p <- data$p
    sj <- data$mktshr
    brand <- model.matrix(~ data$firmids - 1)
    ownership <- brand %*% t(brand)
    omega <- (alpha) * (diag(sj) - sj %*% t(sj)) * ownership
    omega_inv <- solve(omega)
    mk <- -omega_inv %*% sj
    return(mk)
}
mk <- markup(data = dt[year == 1977])

elasticity <- function(data) {
    p <- data$p
    sj <- data$mktshr
    omega <- (alpha) * (diag(sj) - sj %*% t(sj))
    omega_inv <- solve(omega)
    elas <- omega_inv %*% sj * p
    return(elas)
}

elas <- elasticity(data = dt[year == 1977])

# merger simulation
dt_merger <- copy(dt)
dt_merger$c <- mc
dt_merger[firmids == 17, firmids := 7]
price <- function(data) {
    c <- data$c
    sj <- data$mktshr
    brand <- model.matrix(~ data$firmids - 1)
    ownership <- brand %*% t(brand)
    omega <- (alpha) * (diag(sj) - sj %*% t(sj)) * ownership
    p_new <- c - solve(omega) %*% sj
    return(p_new)
}

# consumer surplus in monetary terms

# # pre-merger
# inclu_util <- -log(s0)
# cs <- num_population * inclu_util / (-alpha)
# # post-merger
# inclu_util <- -log(1 - sum(sj))
# cs <- num_population * inclu_util / (-alpha)
