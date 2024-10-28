rm(list = ls())
pacman::p_load(data.table, AER, gmm)
load("Data/blp_results.rda")

summary(ivreg)
alpha <- ivreg$coefficients["p"]

# calculate the marginal cost
marginal_cost_1 <- function(c, data) {
    p <- data$p
    sj <- data$mktshr
    brand <- model.matrix(~ data$firmids - 1)
    ownership <- brand %*% t(brand)
    omega <- (-alpha) * (diag(sj) - sj %*% t(sj)) * ownership
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
    omega <- (-alpha) * (diag(sj) - sj %*% t(sj)) * ownership
    mc <- p + solve(omega) %*% sj
    return(mc)
}
mc <- marginal_cost_2(data = dt[year == 1977])

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
# pre-merger
inclu_util <- -log(s0)
cs <- num_population * inclu_util / (-alpha)
# post-merger
inclu_util <- -log(1 - sum(sj))
cs <- num_population * inclu_util / (-alpha)
