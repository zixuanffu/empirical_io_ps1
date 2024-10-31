rm(list = ls())
pacman::p_load(data.table, fixest)
dt <- readRDS("Data/carpanel_q2.rdS")

# nested based on size
dt[, cat := as.factor(cat)]
# in addition to the previous instrument for price, we will construct new instruments for market share in the nest
dt[, shr_g := sum(mktshr), by = .(year, cat)]
dt[, mktshr_g := mktshr / shr_g]


var_exo <- c("const", "dpm", "door3", "door4", "door5", "at", "ps", "air", "drv", "wt", "hp2wt", "hp", "euro", "japan", "size", "wb")
var_end <- c("p", "log(mktshr_g)")

# instr class 1: the number of products in the nest
dt[, num_products_own_g := .N, by = .(year, cat, firmids)]
dt[, num_products_rival_g := .N - num_products_own_g, by = .(year, cat)]


# instr class 2: the characteristic of other firms in the group
char_lst <- var_exo
for (iv in char_lst) {
    dt[, paste0(iv, "_own_g") := sum(get(iv)), by = .(year, firmids, cat)]
    dt[, paste0(iv, "_rival_g") := sum(get(iv)) - get(paste0(iv, "_own_g")), by = .(year, cat)]
}
saveRDS(dt, "Data/carpanel_q3.rds")
IV_lst <- c(paste(char_lst, "rival_g", sep = "_"), "num_products_rival_g")
saveRDS(IV_lst, "Data/iv_q3.rds")

# estimation: ols
ols_formula1 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end), collapse = " + ")))
reg_ols1 <- feols(ols_formula1, data = dt)
summary(reg_ols1, cluster = "modelid")
ols_formula2 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end), collapse = " + "), "|firmids"))
ols_formula2
reg_ols2 <- feols(ols_formula2, data = dt)
summary(reg_ols2, cluster = "modelid")

# estimation: iv
var_iv <- c("dpm_rival", "hp2wt_rival", "air_rival", "air_ownothers", "num_products_rival_g", "hp2wt_rival_g", "dpm_rival_g", "air_rival_g")
iv_formula1 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(var_exo, collapse = " + "), "|", paste(var_end, collapse = " + "), "~", paste(var_iv, collapse = " + ")))
reg_iv1 <- feols(iv_formula1, data = dt)
summary(reg_iv1, cluster = "modelid")
iv_formula2 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(var_exo, collapse = " + "), "|firmids|", paste(var_end, collapse = " + "), "~", paste(var_iv, collapse = " + ")))
reg_iv2 <- feols(iv_formula2, data = dt)
summary(reg_iv2, cluster = "modelid")

# export
etable(reg_ols1, reg_ols2, reg_iv1, reg_iv2,
    cluster = "modelid",
    headers = list(Estimator = list("OLS" = 2, "IV" = 2)),
    fitstat = ~ n + ar2 + ivwald1 + ivwald1.p + sargan + sargan.p,
    title = "Logit",
    file = "Results/Tables/nested_logit.tex", tex = TRUE, replace = TRUE
)
# appendix
# diagnostic test: https://cran.r-project.org/web/packages/ivreg/vignettes/Diagnostics-for-2SLS-Regression.html
