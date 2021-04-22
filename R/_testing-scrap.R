rm(list = ls())

# require(data.table)
# require(fixest)
# require(tidyverse)

# formula = as.formula("incearn_ln ~ reform_math + x1 + x2")
# data = mutate(bacondecomp::math_reform,
#   x1 = rnorm(n()),
#   x2 = rnorm(n()))
# id_var = "state"
# time_var = "class"

# _______________________________________________________________
## ----
# _______________________________________________________________

source("R/bacon.R")
test_data <- data.table::data.table(bacondecomp::math_reform)

bacon(
  incearn_ln ~ reform_math,
  data = test_data,
  id_var = "state",
  time_var = "class"
)

# _______________________________________________________________
# installing package from gh:
# devtools::install_github("darinchristensen/bacondecomp", force = TRUE)
# 
# test_data <- data.table::data.table(bacondecomp::math_reform)
# 
# bacondecomp::bacon(
#   incearn_ln ~ reform_math,
#   data = test_data,
#   id_var = "state",
#   time_var = "class"
# )


# profvis::profvis(
#   bacondecomp::bacon(
#     incearn_ln ~ reform_math,
#     data = bacondecomp::math_reform,
#     id_var = "state",
#     time_var = "class"
#   )
# )

system.time(
  bacondecomp::bacon(
    incearn_ln ~ reform_math,
    data = bacondecomp::math_reform,
    id_var = "state",
    time_var = "class"
  )
)

# _______________________________________________________________

test_data <- bacondecomp::castle %>% data.table()
system.time(
  bacon(
    l_homicide ~ post + l_pop + l_income,
    dt = test_data,
    id_var = "state",
    time_var = "year"
  )
)

system.time(
  bacondecomp::bacon(
    l_homicide ~ post + l_pop + l_income,
    data = bacondecomp::castle,
    id_var = "state",
    time_var = "year"
  )
)


