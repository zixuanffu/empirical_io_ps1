rm(list = ls())
pacman::p_load(data.table)
# construct market share
car <- fread("Data/carpanel.csv") # this is aggregate data not individual data
hh <- fread("Data/UShouseholds.csv")
car <- merge(car, hh, by = "year")
car[, mktshr := q / (nb_hh * 1000)]

# construct dolar per miles
gas <- fread("Data/gasprice.csv")
car_gas <- merge(car, gas, by = "year") # mpg is miles per gallon
cpi <- fread("Data/USCPI.csv")
car_gas <- merge(car_gas, cpi, by = "year")
car_gas[, p := p * 100 / cpi]
car_gas[, gasprice := gasprice * 100 / cpi]
car_gas[, dpm := gasprice / mpg]
car_gas[, p := p / 1000]
# save
saveRDS(car_gas, "Data/carpanel.rds")
