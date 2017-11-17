library(causaldrf)

hi_sample <- function(N){
  X1 <- rexp(N)
  X2 <- rexp(N)
  T <- rexp(N, X1 + X2)
  gps <- (X1 + X2) * exp(-(X1 + X2) * T)
  Y <- T + gps + rnorm(N)
  hi_data <- data.frame(cbind(X1, X2, T, gps, Y))
  return(hi_data)
}
hi_sim_data <- hi_sample(1000)

add_spl_estimate <- add_spl_est(Y = Y,
                                treat = T,
                                treat_formula = T ~ X1 + X2,
                                data = hi_sim_data,
                                grid_val = quantile(hi_sim_data$T,
                                                    probs = seq(0, .95, by = 0.01)),
                                knot_num = 3,
                                treat_mod = "Gamma",
                                link_function = "inverse")

gam_estimate <- gam_est(Y = Y,
                        treat = T,
                        treat_formula = T ~ X1 + X2,
                        data = hi_sim_data,
                        grid_val = quantile(hi_sim_data$T,
                                            probs = seq(0, .95, by = 0.01)),
                        treat_mod = "Gamma",
                        link_function = "inverse")


hi_estimate <- hi_est(Y = Y,
                      treat = T,
                      treat_formula = T ~ X1 + X2,
                      outcome_formula = Y ~ T + I(T^2) +
                        gps + I(gps^2) + T * gps,
                      data = hi_sim_data,
                      grid_val = quantile(hi_sim_data$T,
                                          probs = seq(0, .95, by = 0.01)),
                      treat_mod = "Gamma",
                      link_function = "inverse")


iptw_estimate <- iptw_est(Y = Y,
                          treat = T,
                          treat_formula = T ~ X1 + X2,
                          numerator_formula = T ~ 1,
                          data = hi_sim_data,
                          degree = 2,
                          treat_mod = "Gamma",
                          link_function = "inverse")


data("nmes_data")
nmes_data <- data.table(nmes_data)

pf_estimate <- reg_est(Y = TOTALEXP,
                       treat = packyears,
                       covar_formula = ~ 1,
                       data = full_data_orig,
                       degree = 2,
                       wt = full_data_orig$HSQACCWT,
                       method = "same")

reg_estimate <- reg_est(Y = TOTALEXP,
                        treat = packyears,
                        covar_formula = ~ LASTAGE + LASTAGE2 +
                          AGESMOKE + AGESMOKE2 + MALE + beltuse +
                          educate + marital + POVSTALB + RACE3,
                        covar_lin_formula = ~ 1,
                        covar_sq_formula = ~ 1,
                        data = full_data_orig,
                        degree = 2,
                        wt = full_data_orig$HSQACCWT,
                        method = "different")

spline_estimate <- prop_spline_est(Y = TOTALEXP,
                                   treat = packyears,
                                   covar_formula = ~ LASTAGE + LASTAGE2 +
                                     AGESMOKE + AGESMOKE2 + MALE + beltuse +
                                     educate + marital + POVSTALB + RACE3,
                                   covar_lin_formula = ~ 1,
                                   covar_sq_formula = ~ 1,
                                   data = full_data_orig,
                                   e_treat_1 = full_data_orig$est_treat,
                                   degree = 2,
                                   wt = full_data_orig$HSQACCWT,
                                   method = "different",
                                   spline_df = 5,
                                   spline_const = 4,
                                   spline_linear = 4,
                                   spline_quad = 4)

ivd_estimate <- prop_spline_est(Y = TOTALEXP,
                                treat = packyears,
                                covar_formula = ~ 1,
                                covar_lin_formula = ~ 1,
                                covar_sq_formula = ~ 1,
                                data = full_data_orig,
                                e_treat_1 = full_data_orig$est_treat,
                                degree = 2,
                                wt = full_data_orig$HSQACCWT,
                                method = "different",
                                spline_df = 5,
                                spline_const = 4,
                                spline_linear = 4,
                                spline_quad = 4)


bart_estimate <- bart_est(Y = iqsb.36,
                          treat = ncdctt,
                          outcome_formula = iqsb.36 ~ ncdctt + bw +
                            female + mom.lths +
                            site1 + site7 + momblack +
                            workdur.imp,
                          data = full_data_orig,
                          grid_val = grid_treat)

iw_estimate <- iw_est(Y = iqsb.36,
                      treat = ncdctt,
                      treat_formula = ncdctt ~ bw + female + mom.lths +
                        site1 + site7 + momblack +
                        workdur.imp,
                      data = full_data_orig,
                      grid_val = grid_treat,
                      bandw = 2 * bw.SJ(full_data_orig$ncdctt),
                      treat_mod = "Normal")

nw_estimate <- nw_est(Y = iqsb.36,
                      treat = ncdctt,
                      treat_formula = ncdctt ~ bw + female + mom.lths +
                        site1 + site7 + momblack +
                        workdur.imp,
                      data = full_data_orig,
                      grid_val = grid_treat,
                      bandw = 2 * bw.SJ(full_data_orig$ncdctt),
                      treat_mod = "Normal")

t_mod_list <- t_mod(treat = T,
                    treat_formula = T ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
                    data = example_data,
                    treat_mod = "Normal")
cond_exp_data <- t_mod_list$T_data
full_data <- cbind(example_data, cond_exp_data)
prop_spline_list <- prop_spline_est(Y = Y,
                                    treat = T,
                                    covar_formula = ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
                                    covar_lin_formula = ~ 1,
                                    covar_sq_formula = ~ 1,
                                    data = example_data,
                                    e_treat_1 = full_data$est_treat,
                                    degree = 1,
                                    wt = NULL,
                                    method = "different",
                                    spline_df = 5,
                                    spline_const = 4,
                                    spline_linear = 4,
                                    spline_quad = 4)
sample_index <- sample(1:1000, 100)
plot(example_data$T[sample_index],
     example_data$Y[sample_index],
     xlab = "T",
     ylab = "Y",
     main = "propensity spline estimate")
abline(prop_spline_list$param[1],
       prop_spline_list$param[2],
       lty = 2,
       col = "blue",
       lwd = 2)
legend('bottomright',
       "propensity spline estimate",
       lty = 2,
       bty = 'Y',
       cex = 1,
       col = "blue",
       lwd = 2)


bart_list <- bart_est(Y = Y,
        treat = A,
        outcome_formula = Y ~ A + W1 + W2 + W3 + W4 + W5 + W6 + W7 + 
                            W8 + W9 + W10 + W11 + W12 + W13 + 
                            W14 + W15 + W16 + W17 + W18 + W19 + 
                            W20 + W21 + W22 + W22 + W23 + W24 + 
                            W25,
        data = dt1,
        grid_val = c(0,1))
sample_index <- sample(1:400, 100)
plot(dt1[sample_index,A],
     dt1[sample_index,Y],
     xlab = "A",
     ylab = "Y",
     main = "bart estimate")
lines(c(1,2),bart_list$param, lty=2,lwd=2, col="blue")
