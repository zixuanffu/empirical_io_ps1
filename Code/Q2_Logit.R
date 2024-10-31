rm(list = ls())
pacman::p_load(data.table, fixest)
dt <- readRDS("Data/carpanel.rds")
# we have 5 years, which corresponds to 5 markets
unique(dt$year)
# calculate outside option share s0
dt[, shr_0 := 1 - sum(mktshr), by = year]
dt[, const := 1]
dt[, door2 := (dr == 2)]
dt[, door3 := (dr == 3)]
dt[, door4 := (dr == 4)]
dt[, door5 := (dr == 5)]
# set categorical variables
dt[, dr := as.factor(dr)]
dt[, firmids := as.factor(firmids)]
dt[, name := as.factor(name)]
dt[, modelid := as.numeric(name)]
# the variables to be included in the regression
var_exo <- c("const", "dpm", "door3", "door4", "door5", "at", "ps", "air", "drv", "wt", "hp2wt", "hp", "euro", "japan", "size", "wb")
var_end <- c("p")

# we construct a set of instruments following blp1995 (hp2wt, air, dpm, size)
char_lst <- var_exo
for (iv in char_lst) {
    dt[, paste0(iv, "_own") := sum(get(iv)), by = .(firmids, year)] # sum of x_jt,k for all products that belong to the same firm in the same year
    dt[, paste0(iv, "_ownothers") := get(paste0(iv, "_own")) - get(iv)] # sum of x_jt,k for other products that belong to the same firm in the same year
    dt[, paste0(iv, "_rival") := sum(get(iv)) - get(paste0(iv, "_own")), by = .(year)] # sum of x_jt,k for all products that belong to the rival firms
}
IV_lst <- c(paste(char_lst, "rival", sep = "_"), paste(char_lst, "ownothers", sep = "_"))
saveRDS(IV_lst, "Data/iv_q2.rds")
saveRDS(dt, "Data/carpanel_q2.rds")

# assume E(\xi_j|p_j,x_j)=0
# ols regression using feols
ols_formula1 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end), collapse = " + ")))
reg_ols1 <- feols(ols_formula1, data = dt)
summary(reg_ols1, cluster = "modelid")
ols_formula2 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end), collapse = " + "), "|firmids"))
reg_ols2 <- feols(ols_formula2, data = dt)
summary(reg_ols2, cluster = "modelid")

# assume p is correlated with unoberved product char, that is E(\xi_j|p_j,x_j) \neq 0
var_iv <- c("dpm_rival", "hp2wt_rival", "air_rival", "air_ownothers")
iv_formula1 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(var_exo, collapse = " + "), "|", paste(var_end, collapse = " + "), "~", paste(var_iv, collapse = " + ")))
reg_iv1 <- feols(iv_formula1, data = dt, cluster = "modelid")
summary(reg_iv1, cluster = "modelid")
iv_formula2 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(var_exo, collapse = " + "), "|firmids|", paste(var_end, collapse = " + "), "~", paste(var_iv, collapse = " + ")))
reg_iv2 <- feols(iv_formula2, data = dt, cluster = "modelid")
summary(reg_iv2, cluster = "modelid")

etable(reg_ols1, reg_ols2, reg_iv1, reg_iv2,
    cluster = "modelid",
    headers = list(Estimator = list("OLS" = 2, "IV" = 2)),
    fitstat = ~ n + ar2 + ivwald1 + ivwald1.p + sargan + sargan.p,
    title = "Logit",
    file = "Results/Tables/logit.tex", tex = TRUE, replace = TRUE
)
