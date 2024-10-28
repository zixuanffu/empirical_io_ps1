rm(list = ls())
pacman::p_load(data.table, plm, AER)
dt <- readRDS("Data/carpanel.rds")
# we have 5 years, which corresponds to 5 markets
unique(dt$year)
# calculate outside option share s0
dt[, shr_0 := 1 - sum(mktshr), by = year]
dt[, const := 1]
# set categorical variables
dt[, dr := as.factor(dr)]
dt[, firmids := as.factor(firmids)]
saveRDS(dt, "Data/carpanel_q2.rds")
# the variables to be included in the regression
var_exo <- c("dpm", "dr", "at", "ps", "air", "drv", "wt", "hp2wt", "hp", "euro", "japan", "size", "wb")
var_end <- c("p")
# assume E(\xi_j|p_j,x_j)=0
ols_formula1 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end), collapse = " + ")))
# ols regression using lm and plm (compare the resutls)
summary(lm(ols_formula1, data = dt))
summary(plm(ols_formula1, data = dt, model = "pooling"))
ols_formula2 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end, "firmids"), collapse = " + ")))
summary(lm(ols_formula2, data = dt))
summary(plm(ols_formula2, data = dt, model = "pooling"))

# assume p is correlated with unoberved product char, that is E(\xi_j|p_j,x_j) \neq 0

# we construct a set of instruments following blp1995 hp2wt, air, dpm, size
char_lst <- c("const", "hp2wt", "air", "dpm", "size")
for (iv in char_lst) {
    dt[, paste0(iv, "_own") := sum(get(iv)), by = .(firmids, year)] # sum of hp2wt for other products that belong to the same firm in the same year
    dt[, paste0(iv, "_ownothers") := get(iv) - get(paste0(iv, "_own"))] # products by the same firm
    dt[, paste0(iv, "_rival") := sum(get(iv)) - get(paste0(iv, "_own")), by = .(year)] # products by the rival firms
}

IV_lst <- c(paste(char_lst, "rival", sep = "_"), paste(char_lst, "ownothers", sep = "_"))
saveRDS(IV_lst, "Data/IV_lst.rds")
iv_formula1 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end), collapse = " + "), "|.-", var_end, "+", paste(IV_lst, collapse = " + ")))
iv_formula1_aer <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end), collapse = " + "), "|", paste(c(var_exo, IV_lst), collapse = " + ")))
# iv regression with aer and plm
summary(ivreg(iv_formula1_aer, data = dt))
summary(plm(iv_formula1, data = dt, model = "pooling"))
iv_formula2 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end, "firmids"), collapse = " + "), "|.-", var_end, "+", paste(IV_lst, collapse = " + ")))
iv_formula2_aer <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end, "firmids"), collapse = " + "), "|", paste(c(var_exo, IV_lst, "firmids"), collapse = " + ")))
summary(ivreg(iv_formula2_aer, data = dt))
summary(plm(iv_formula2, data = dt, model = "pooling"))
