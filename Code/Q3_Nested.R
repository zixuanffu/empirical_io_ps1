rm(list = ls())
pacman::p_load(data.table, plm, AER)
dt <- readRDS("Data/carpanel_q2.rdS")

# nested based on size
dt[, cat := as.factor(cat)]
# in addition to the previous instrument for price, we will construct new instruments for market share in the nest
dt[, shr_g := sum(mktshr), by = .(year, cat)]
dt[, mktshr_g := mktshr / shr_g]

# attempt 0: ols
var_exo <- c("dpm", "dr", "at", "ps", "air", "drv", "wt", "hp2wt", "hp", "euro", "japan", "size", "wb")
var_end <- c("p", "log(mktshr_g)")
ols_formula1 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end), collapse = " + ")))
summary(plm(ols_formula1, data = dt, model = "pooling"))
ols_formula2 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end, "firmids"), collapse = " + ")))
summary(plm(ols_formula2, data = dt, model = "pooling"))

# attempt 1: the number of products in the nest
dt[, num_products_g := .N, by = .(year, cat)] # not a lot of variation



# attempt 2: the characteristic of other firms in the group
char_lst <- c("const", "hp2wt", "air", "dpm", "size")
for (iv in char_lst) {
    dt[, paste0(iv, "_own") := sum(get(iv)), by = .(firmids, year)] # sum of hp2wt for other products that belong to the same firm in the same year
    dt[, paste0(iv, "_ownothers") := get(iv) - get(paste0(iv, "_own"))] # products by the same firm
    dt[, paste0(iv, "_rival") := sum(get(iv)) - get(paste0(iv, "_own")), by = .(year)] # products by the rival firms
}

for (iv in char_lst) {
    dt[, paste0(iv, "_own_g") := sum(get(iv)), by = .(year, firmids, cat)]
    dt[, paste0(iv, "_rival_g") := sum(get(iv)) - get(paste0(iv, "_own_g")), by = .(year, cat)]
}
saveRDS(dt, "Data/carpanel_q3.rds")
IV_lst <- c(paste(char_lst, "rival", sep = "_"), paste(char_lst, "ownothers", sep = "_"))
IV_lst <- c(IV_lst, paste(char_lst, "rival_g", sep = "_"))
saveRDS(IV_lst, "Data/IV_lst_nested.rds")
# iv regression with plm
iv_formula1 <- as.formula(
    paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end), collapse = " + "), "|.-", paste(c(var_end), collapse = " - "), "+", paste(IV_lst, collapse = " + "))
)
iv_reg1 <- plm(iv_formula1, data = dt, model = "pooling")
summary(iv_reg1)
# iv regression with aer
iv_formula1_aer <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end), collapse = " + "), "|", paste(c(var_exo, IV_lst), collapse = " + ")))
iv_reg1_aer <- ivreg(iv_formula1_aer, data = dt)
summary(iv_reg1_aer, diagnostics = TRUE)

# with brand fixed effect firmids
iv_formula2 <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end, "firmids"), collapse = " + "), "|.-", paste(c(var_end), collapse = " - "), "+", paste(IV_lst, collapse = " + ")))
iv_reg2 <- plm(iv_formula2, data = dt, model = "pooling")
summary(iv_reg2)

iv_formula2_aer <- as.formula(paste("log(mktshr)-log(shr_0) ~", paste(c(var_exo, var_end, "firmids"), collapse = " + "), "|", paste(c(var_exo, IV_lst, "firmids"), collapse = " + ")))
iv_reg2_aer <- ivreg(iv_formula2_aer, data = dt)
summary(iv_reg2_aer, diagnostics = TRUE)

# appendix
# diagnostic test: https://cran.r-project.org/web/packages/ivreg/vignettes/Diagnostics-for-2SLS-Regression.html
