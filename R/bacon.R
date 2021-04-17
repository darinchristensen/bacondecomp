#' Goodman-Bacon Decomposition
#'
#' bacon() is a function that performs the Goodman-Bacon decomposition for
#'  differences-in-differences with variation in treatment timing (with or
#'  without time-varying covariates).
#'
#' @param formula an object of class "formula": a symbolic
#'  representation of the model to be fitted. Must be  of the form y ~ D + controls,
#'  where y is the outcome variable,  D is the binary
#'  treatment indicator, and `controls` can be any additional control variables. Do not
#'  include the fixed effects in the formula. If using `.` notation must be of
#'  the form y ~ D + . - FE1 - FE2
#' @param dt a data.table containing the variables in the model.
#' @param id_var character, the name of id variable for units.
#' @param time_var character, the name of time variable.
#' @param quietly logical, default = FALSE, if set to TRUE then bacon() does not
#'  print the summary of estimates/weights by type (e.g. Treated vs Untreated)
#'
#' @return If control variables are included in the formula, then an object of
#'  class "list" with three elements:
#'  \item{Omega}{a number between 0 and 1, the weight of the within timing group
#'   coefficient}
#'  \item{beta_hat_w}{a number, the within timing group coefficient}
#'  \item{two_by_twos}{a data.frame with the covariate adjusted 2x2 estimates
#'   and weights}
#'
#' If not control variables are included then only the two_by_twos data.frame
#'  is returned.
#'
#' @examples
#' # Castle Doctrine (Uncontrolled)
#' df_bacon <- bacon(l_homicide ~ post,
#'                   data = bacondecomp::castle,
#'                   id_var = "state",
#'                   time_var = "year")
#'
#' # Castle Doctrine (Controlled)
#' ret_bacon <- bacon(l_homicide ~ post + l_pop + l_income,
#'                    data = bacondecomp::castle,
#'                    id_var = "state",
#'                    time_var = "year")
#'
#' @import stats
#' @import tidyverse
#' @import fixest
#' @import data.table
#'
#' @export
bacon <- function(formula,
  dt,
  id_var,
  time_var,
  quietly = FALSE) {
  
  # loading required packages:
  pkgs <- c("tidyverse", "fixest", "data.table")
  suppressPackageStartupMessages(sapply(pkgs, require, character.only = TRUE))
  
  # Evaluate formula in data environment
  formula <- formula(terms (formula, data = dt))
  dt <- copy(dt)
  
  # Unpack variable names and rename variables
  vars <- unpack_variable_names(formula)
  outcome_var <- vars$outcome_var
  treated_var <- vars$treated_var
  control_vars <- vars$control_vars
  dt <- rename_vars(dt, id_var, time_var, outcome_var, treated_var)
  
  # Check for NA observations
  nas <- sum(is.na(dt[, c("id", "time", "outcome", "treated"), with = FALSE]))
  
  if (length(control_vars > 0)) {
    nas_control <- sum(is.na(dt[, control_vars, with = FALSE]))
    nas <- nas + nas_control
    rm(nas_control)
  }
  
  if (nas > 0) {
    stop("NA observations")
  }
  rm(nas)
  
  # Create 2x2 grid of treatment groups
  treatment_group_calc <- create_treatment_groups(dt, control_vars, return_merged_df = TRUE)
  two_by_twos <- treatment_group_calc$two_by_twos
  dt <- treatment_group_calc$dt
  rm(treatment_group_calc)
  
  two_by_twos[, c("estimate", "weight") := NA_real_]
  
  # _______________________________________________________________
  # Uncontrolled ----
  if (length(control_vars) == 0) {
    # Iterate through treatment group dyads
    for (i in 1:nrow(two_by_twos)) {
      data1 <- subset_data(dt, two_by_twos$treated[i], two_by_twos$untreated[i])
      
      # Calculate estimate and weight
      weight <- calculate_weights(
        data = data1,
        treated_group = two_by_twos$treated[i],
        untreated_group = two_by_twos$untreated[i]
      )
      
      estimate <- feols(outcome ~ treated | time + id, data = data1)$coefficients
      
      two_by_twos$estimate[i] <- estimate
      two_by_twos$weight[i] <- weight
    }
    
    # Rescale weights to sum to 1
    two_by_twos <- scale_weights(two_by_twos)
    
    if (quietly == FALSE) {
      print_summary(two_by_twos)
    }
    
    return(two_by_twos)
    
  } else if (length(control_vars) > 0) {
    # _______________________________________________________________
    # Controled ----
    # Predict Treatment and calulate demeaned residuals
    control_formula <- str_c("treated ~ ", 
      str_flatten(control_vars, " + "), " | time + id") %>%
      as.formula()
    
    dt <- run_fwl(dt, control_formula)
    
    # Calculate within treatment group estimate and its weight
    Omega <- calculate_Omega(dt)
    beta_hat_w <- calculate_beta_hat_w(dt)
    
    # Collapse controls and predicted treatment to treatment group/year level
    r_collapse_x_p <- collapse_x_p(dt, control_vars)
    dt <- r_collapse_x_p$data
    g_control_formula <- r_collapse_x_p$g_control_formula
    rm(r_collapse_x_p)
    
    # Iterate through treatment group dyads
    for (i in 1:nrow(two_by_twos)) {
      data1 <- dt[treat_time == two_by_twos$treated[i] | 
          treat_time ==  two_by_twos$untreated[i]]
      
      # Calculate between group estimate and weight
      weight_est <- calc_controlled_beta_weights(data1, g_control_formula)
      s_kl <- weight_est$s_kl
      beta_hat_d_bkl <- weight_est$beta_hat_d_bkl
      rm(weight_est)
      
      two_by_twos$weight[i] <- s_kl
      two_by_twos$estimate[i] <- beta_hat_d_bkl
    }
    
    # Rescale weights to sum to 1
    two_by_twos <- scale_weights(two_by_twos)
    
    if (quietly == FALSE) {
      print_summary(two_by_twos)
    }
    
    r_list <- list(
      "beta_hat_w" = beta_hat_w,
      "Omega" = Omega,
      "two_by_twos" = two_by_twos
    )
    return(r_list)
  }
}


#' Unpack Variable Names from Formula
#'
#' @param formula formula call from bacon()
#'
#' @return A list with 3 elements:
#' \item{outcome_var}{name of the outcome variable from the LHS of formula}
#' \item{treated_var}{name of the treatment indicator variable, the first
#'  variable from the RHS of formula}
#' \item{control_vars}{All variables on RHS of formula except the treated_var.
#'  Only need to know if length(control_vars) > 0}
#'
#' @noRd
unpack_variable_names <- function(formula) {
  outcome_var <- as.character(formula)[2]
  right_side_vars <- as.character(formula)[3]
  right_side_vars <- strsplit(right_side_vars, " \\+ ")[[1]]
  treated_var <- right_side_vars[1]
  control_vars <- right_side_vars[-1]
  r_list <-
    list(
      outcome_var = outcome_var,
      treated_var = treated_var,
      control_vars = control_vars
    )
  return(r_list)
}


#' Rename Variables
#'
#' Rename id, time, outcome, and treatment indicator variables
#'
#' @param data a data.frame
#' @param id_var character, name of id variable
#' @param time_var character, name of time variable
#' @param outcome_var character, name of outcome variable
#' @param treated_var character, name of binary treatment variable
#'
#' @return A data.frame with renamed columns
#'
#' @noRd
rename_vars <-
  function(data,
    id_var,
    time_var,
    outcome_var,
    treated_var) {
    setnames(
      data,
      c(id_var, time_var, outcome_var, treated_var),
      c("id", "time", "outcome", "treated")
    )
    return(data)
  }


#' Create Grid of Treatment Groups
#'
#' @param dt a data.table used to create groups - MUST obey naming convention
#' used in `bacon()`. i.e. columns are ["id", "time", "outcome", "treated"]
#'
#' @param return_merged_df Defaults to `FALSE` whether to return merged data
#' as well as grid of treatment groups.
#' @param control_vars list of control variables
#'
#' @return If return_merge_df is `TRUE` then a list with two elements:
#' \item{two_by_twos}{data.frame with all treatment group dyads and NA
#'  weight/estimate columns}
#' \item{data}{updated data.frame with the new variable treat_time}
#'
#' If return_merge_df is `FALSE` then only the two_by_twos data.frame is
#'  returned
#'
#' @noRd
create_treatment_groups <-
  function(dt, control_vars, return_merged_df = FALSE) {
    df_treat <- dt[treated == 1, c("id", "time"), with = FALSE]
    df_treat <- df_treat[, list(time = min(time)), by = id]
    setnames(df_treat, "time", "treat_time")
    
    setkey(dt, id)
    setkey(df_treat, id)
    dt <- merge(dt, df_treat, all.x = TRUE)
    dt[is.na(treat_time), treat_time := 99999]
    
    # Check for weakly increasing treatment
    inc <-
      dt[(time >= treat_time &
          treated == 0) | (time < treat_time & treated == 1)] %>% nrow()
    if (inc > 0) {
      stop("Treatment not weakly increasing with time")
    }
    rm(inc)
    
    # inc <- ifelse(data$treat_time == 99999, 1,
    #   ifelse(data$time >= data$treat_time & data$treated == 1, 1,
    #     ifelse(data$time < data$treat_time & data$treated == 0,
    #       1, 0)))
    ttimes <- unique(dt$treat_time)
    min_time <- min(dt$time)
    two_by_twos <-
      expand_grid(treated = ttimes, untreated = ttimes) %>% data.table()
    two_by_twos <-
      two_by_twos[treated != untreated][treated != 99999][treated != min_time]
    
    if (length(control_vars) == 0) {
      # Create data.frame of all posible 2x2 estimates. Dyads may appear twice as
      # treatment groups can play the roll of both earlier and later treated
      two_by_twos[, type := case_when(
        untreated == 99999 ~ "Treated vs. Untreated",
        untreated == min_time ~ "Later vs Always Treated",
        treated > untreated ~ "Later vs. Earlier Treated",
        TRUE ~ "Earlier vs. Later Treated"
      )]
      
    } else if (length(control_vars) > 0) {
      # In the controlled decomposition, each dyad only appears once because we
      # do not make the distinction between earlier vs later treated
      two_by_twos <-
        two_by_twos[treated < untreated | untreated == 99999]
      
      two_by_twos[, type := case_when(
        untreated == 99999 ~ "Treated vs. Untreated",
        untreated == min_time |
          treated == min_time ~ "Later vs Always Treated",
        TRUE ~ "Both Treated"
      )]
    }
    
    # Whether or not to return the merged data too.
    if (return_merged_df == TRUE) {
      return_data <- list("two_by_twos" = two_by_twos, "dt" = dt)
    } else {
      return_data <- two_by_twos
    }
    return(return_data)
  }


#' Subset Data
#'
#' In the uncontrolled decomposition, subset_data() subsets the data to only
#'  two treatment groups, and subsets the time to either before the later group
#'  has been treated or after the earlier group has been treated.
#'
#' @param data a data.frame
#' @param treated_group integer, initial treatment time of group acting as
#'  treated
#' @param untreated_group interger, initial treatment time of group acting as
#'  untreated
#'
#' @return subsetted data.frame
#'
#' @noRd
subset_data <- function(data, treated_group, untreated_group) {
  data <- data[treat_time %in% c(treated_group, untreated_group)]
  if (treated_group < untreated_group) {
    data <- data[time < untreated_group]
  } else if (treated_group > untreated_group) {
    data <- data[time >= untreated_group]
  }
  return(data)
}


#' Calculate Weights for 2x2 Grid (uncontrolled)
#'
#'  Calculated weights for the uncontrolled decomposition using:
#'  n_u - observations in untreated group,
#'  n_k - observations in earlier treated group,
#'  n_l - observations in later treated group,
#'  D_k - proportion of time the earlier treated group was treated,
#'  D_l - proportion of time the later treated group was treated.
#'
#' @param data a data.frame
#' @param treated_group the identifier of the treated group
#' @param untreated_group the identifier of the untreated group
#'
#' @return Scalar weight for 2x2 grid
#'
#' @noRd
calculate_weights <-
  function(data, treated_group, untreated_group) {
    if (untreated_group == 99999) {
      # Treated vs untreated
      n_u <- sum(data$treat_time == untreated_group)
      n_k <- sum(data$treat_time == treated_group)
      n_ku <- n_k / (n_k + n_u)
      D_k <- mean(data[treat_time == treated_group]$treated)
      V_ku <- n_ku * (1 - n_ku) * D_k * (1 - D_k)
      weight1 <- (n_k + n_u) ^ 2 * V_ku
    } else if (treated_group < untreated_group) {
      # early vs late (before late is treated)
      n_k <- sum(data$treat_time == treated_group)
      n_l <- sum(data$treat_time == untreated_group)
      n_kl <- n_k / (n_k + n_l)
      D_k <- mean(data[treat_time == treated_group]$treated)
      D_l <- mean(data[treat_time == untreated_group]$treated)
      V_kl <-
        n_kl * (1 - n_kl) * (D_k - D_l) / (1 - D_l) * (1 - D_k) / (1 - D_l)
      weight1 <- ((n_k + n_l) * (1 - D_l)) ^ 2 * V_kl
    } else if (treated_group > untreated_group) {
      # late vs early (after early is treated)
      n_k <- sum(data$treat_time == untreated_group)
      n_l <- sum(data$treat_time == treated_group)
      n_kl <- n_k / (n_k + n_l)
      D_k <- mean(data[treat_time == untreated_group]$treated)
      D_l <- mean(data[treat_time == treated_group]$treated)
      V_kl <- n_kl * (1 - n_kl) * (D_l / D_k) * (D_k - D_l) / (D_k)
      weight1 <- ((n_k + n_l) * D_k) ^ 2 * V_kl
    }
    return(weight1)
  }


#' Scale Weights
#'
#' scale_weights scales the two_by_two weights so that they sum to 1
#'
#' @param two_by_twos data.frame containing 2x2 weights
#'
#' @return two_by_twos data.frame with rescaled weight column
#'
#' @noRd
scale_weights <- function(two_by_twos) {
  two_by_twos[, weight := weight / sum(weight)]
  return(two_by_twos)
}

#' Run Frisch-Waugh-Lowell Regression
#'
#' Predicts treatment using time-varing covariates and calculates predicted
#'  values/residuals
#'
#' @param data a data.frame
#' @param control_formula a fomula
#'
#' @return updated data.frame with predictions/residuals (and their
#'  transformations) from the FWL regession
#'
#' @noRd
run_fwl <- function(data, control_formula) {
  fit_fwl <- feols(control_formula, data = data)
  data$p <- predict(fit_fwl)
  data$d <- data$d_it <- fit_fwl$residuals
  
  # demean residulas
  data[, d_i_bar := mean(d), by = id]
  data[, d_t_bar := mean(d_it), by = id]
  data[, d_bar_bar := mean(d_it)]
  
  data[, d_it_til := d_it - d_i_bar - d_t_bar + d_bar_bar]
  
  data[, d_kt_bar := mean(d_it), by = list(treat_time, time)]
  data[, d_k_bar := mean(d_it), by = treat_time]
  
  data[, d_ikt_til := d_it - d_i_bar - (d_kt_bar - d_k_bar)]
  data[, d_kt_til := (d_kt_bar - d_k_bar) - (d_t_bar - d_bar_bar)]
  
  return(data)
}


#' Calculate Omega
#'
#' Calculate the weight for the within estimate in the controlled decompositon
#'
#' @param data a data frame with the columns d_ikt_til and d_it_til
#'
#' return Omega
#'
#' @noRd
calculate_Omega <- function(data) {
  # TODO Test that within + between = 1
  N <- nrow(data)
  Vd_w <- var(data$d_ikt_til) * (N - 1) / N
  V_d <- var(data$d_it_til) * (N - 1) / N
  Omega <- Vd_w / V_d
  return(Omega)
}

#' Calculate 1 - Omega
#'
#' Calculate the weight for the between estimates. This is used for testing
#'  purposes only.
#'
#' @param data a data frame with the columns d_kt_til and d_it_til
#'
#' @return 1 - Omega
#'
#' @noRd
calculate_one_minus_Omega <- function(data) {
  N <- nrow(data)
  Vd_b <- var(data$d_kt_til) * (N - 1) / N
  V_d <- var(data$d_it_til) * (N - 1) / N
  one_minus_Omega <- Vd_b / V_d
  return(one_minus_Omega)
}


#' Calculate Within Estimate
#'
#' @param data a data.frame with the columns outcome and d_ikt_til
#'
#' @return beta_hat_w, the within estimate
#'
#' @noRd
calculate_beta_hat_w <- function(data) {
  N <- nrow(data)
  C <- cov(data$outcome, data$d_ikt_til) * (N - 1) / N
  Vd_w <- var(data$d_ikt_til) * (N - 1) / N
  beta_hat_w <- C / Vd_w
  return(beta_hat_w)
}


#' Calculate Between Estimate
#'
#' Calculate the overall between estimate (not decomposed by treatement group
#'  dyads). This is used for testing purposes only.
#'
#' @param data a data.frame with the columns outcome and d_kt_til
#'
#' @return beta_hat_w, the between estimate
#'
#' @noRd
calculate_beta_hat_b <- function(data) {
  N <- nrow(data)
  C <- cov(data$outcome, data$d_kt_til) * (N - 1) / N
  Vd_b <- var(data$d_kt_til) * (N - 1) / N
  beta_hat_b <- C / Vd_b
  return(beta_hat_b)
}


#' Collapse Control Variable and Predicted Treatment
#'
#' Collapse controll variables and predicted treatment to treatment group/time
#'  level.
#'
#' @param data a data.frame
#' @param control_formula, a formula
#'
#' @return a list of two elements:
#' \item{data}{updated data frame with new treatment group/time level variables}
#' \item{g_control_formula}{the control formula but with group level variable
#'  names}
#'
#' @noRd
collapse_x_p <- function(data, control_vars) {
  # Group level Xs
  data[, str_c("g_", control_vars) := lapply(.SD, mean), 
    .SDcols = control_vars, 
    by = list(treat_time, time)]
  
  # Group level p
  data[, g_p := mean(p), by = list(treat_time, time)]
  
  g_control_formula <- str_c("g_", control_vars) %>%
    str_flatten(" + ") %>%
    str_c("~ ", .) %>%
    as.formula()
  
  r_list <- list(data = data, g_control_formula = g_control_formula)
  
  return(r_list)
}


#' Calculate Treatment Variance
#'
#' @param data a data.frame
#'
#' @return a list of two items:
#' \item{data}{a data.frame with the new demeaned treatment variable}
#' \item{VD}{the varince of the demeaned treatment variable}
#'
#' @noRd
calc_VD <- function(data) {
  fit_D <- feols(treated ~ 1 | id + time, data = data)
  data$Dtilde <- fit_D$residuals
  N <- nrow(data)
  VD <- var(data$Dtilde) * (N - 1) / N
  r_list <- list(data = data, VD = VD)
  return(r_list)
}


#' Partial Out FE
#'
#' Partial out time and unit fixed effects from group level control and
#'  and predicted treatment variables
#'
#' @param data, a data.frame with the variables in the g_control_formula and
#'  p
#' @param g_control_formula control formula with group level variable names
#'
#' @return updated data.frame with new partialled out group level variables
#'
#' @noRd
partial_group_x <- function(data, g_control_formula) {
  g_vars <- as.character(g_control_formula)[2] %>% str_split(" \\+ ") %>% unlist()
  
  for (v in g_vars) {
    g_formula <- str_c(v, "~ 1 | time + id") %>% as.formula()
    fit_g <- feols(g_formula, data = data)
    data[, str_c("p_", v) := fit_g$residuals]
  }
  return(data)
}


#' Calculate pgjtile and Rsq
#'
#' @param data, a data.frame with the variables in g_control_formula and
#'  Dtilde
#' @param g_control_formula control formula with group level variable names
#'
#' @return a list with two elements:
#' \item{data}{updated data.frame with new variable pgjtilde}
#' \item{Rsq}{Rsq from the regression to generate pgjtilde}
#'
#' @noRd
calc_pgjtile <- function(data, g_control_formula) {
  p_g_vars <- as.character(g_control_formula)[2] %>% 
    str_split(" \\+ ") %>% 
    unlist() %>%
    str_c("p_", .)
  
  p_g_control_formula <- str_c("Dtilde~", str_flatten(p_g_vars, "+")) %>%
    as.formula()
  
  fit_pgj <- lm(p_g_control_formula, data = data)
  data$pgjtilde <- predict(fit_pgj)
  Rsq <- summary(fit_pgj)$r.squared
  r_list <- list(data = data, Rsq = Rsq)
  return(r_list)
}


#' Calculate Vdp
#'
#' @param data a data.frame with the variables g_p, id, time, and ptilde
#'
#' @return list with two elements
#' \item{data}{data.frame with new variables ptilde and dp}
#' \item{Vdp}{Variance of dp}
#'
#' @noRd
calc_Vdp <- function(data) {
  N <- nrow(data)
  fit_p <- feols(g_p ~ 1 | id + time, data = data)
  data$ptilde <- fit_p$residuals
  data$dp <- data$pgjtilde - data$ptilde
  Vdp <- var(data$dp) * (N - 1) / N
  r_list <- list(data = data, Vdp = Vdp)
  return(r_list)
}


#' Calculate BD
#'
#' @param data a data.frame with all variables in the g_control_formula,
#'  outcome, treated, id, and time
#' @param g_control_formula control formula with group level variable names
#' @return BD
#'
#' @noRd
calc_BD <- function(data, g_control_formula) {
  BD_formula <- str_c(
    "outcome ~ treated + ", 
    as.character(g_control_formula)[2],
    " | time + id"
  ) %>%
    as.formula()
  
  fit_BD <- feols(BD_formula, data = data)
  BD <- fit_BD$coefficients["treated"]
  return(BD)
}


#' Calculate Bb
#'
#' @param data a data.frame with the variables outcome, dp, time, and id
#'
#' @return Bb
#'
#' @noRd
calc_Bb <- function(data) {
  fit_Bb <- feols(outcome ~ dp | time + id, data = data)
  Bb <- fit_Bb$coefficients["dp"]
  return(Bb)
}


#' Calculate Dyad Covariate Adjusted Between Beta Hat
#'
#' @param Rsq R squared from the regression in calc_pgjtilde
#' @param VD variance of the demeaned treatment indicator (in the dyad)
#' @param BD coeffieicient on treatment from a regression of the outcome on
#'  treatment and group level control variables with unite and time FE
#' @param Vdp variance of dp (pgjtilde - ptilde)
#' @param Bb coefficient on dp from a regression of the outcome on dp with unit
#'  and time fixed effects
#'
#' @return dyad's covariate adjusted between estimate
#'
#' @noRd
calculate_beta_hat_d_bkl <- function(Rsq, VD, BD, Vdp, Bb) {
  beta_hat_d_bkl <-
    ((1 - Rsq) * VD * BD + Vdp * Bb) / ((1 - Rsq) * VD + Vdp)
  return(beta_hat_d_bkl)
}

#' Calculate Dyad Weight (Covariate Adjusted)
#'
#' @param N observations in dyad
#' @param Rsq R squared from the regression in calc_pgjtilde
#' @param VD variance of the demeaned treatment indicator (in the dyad)
#' @param Vdp variance of dp (pgjtilde - ptilde)
#'
#' @return the weight given to the dyad's estimate
#'
#' @noRd
calculate_s_kl <- function(N, Rsq, VD, Vdp) {
  s_kl <- N ^ 2 * ((1 - Rsq) * VD + Vdp)
  return(s_kl)
}


#' Calculate Covariate Adjusted Betas and Weights
#'
#' Wrapper function for finding the dyad beta/weight in the controled
#'  decomposition
#'
#' @param data a data.frame
#' @param g_control_formula formula with group/timelevel control variables
#'
#' @return a list with two elements
#' \item{s_kl}{the weight given to the dyad kl's estimate}
#' \item{beta_hat_d_bkl}{the covariate adjusted 2x2 estimate for the dyad kl}
#'
#' @noRd
calc_controlled_beta_weights <- function(data, g_control_formula) {
  r_calc_VD <- calc_VD(data)
  VD <- r_calc_VD$VD
  data <- r_calc_VD$data
  
  data <- partial_group_x(data, g_control_formula)
  
  r_calc_pgjtilde <- calc_pgjtile(data, g_control_formula)
  Rsq <- r_calc_pgjtilde$Rsq
  data <- r_calc_pgjtilde$data
  
  r_calc_Vdp <- calc_Vdp(data)
  Vdp <- r_calc_Vdp$Vdp
  data <- r_calc_Vdp$data
  
  BD <- calc_BD(data, g_control_formula)
  Bb <- calc_Bb(data)
  
  N <- nrow(data)
  
  s_kl <- calculate_s_kl(N, Rsq, VD, Vdp)
  beta_hat_d_bkl <- calculate_beta_hat_d_bkl(Rsq, VD, BD, Vdp, Bb)
  
  r_list <- list(s_kl = s_kl, beta_hat_d_bkl = beta_hat_d_bkl)
  return(r_list)
}

print_summary <- function(two_by_twos, return_df = FALSE) {
  sum_df <- two_by_twos[, list(weight = sum(weight),
    avg_est = weighted.mean(estimate, weight)),
    by = type]
  
  sum_df <- mutate_if(sum_df, is.numeric, round, 5)
  print(sum_df)
  
  if (return_df == TRUE) {
    return(sum_df)
  }
}
