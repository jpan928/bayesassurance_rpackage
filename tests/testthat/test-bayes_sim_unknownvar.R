library(testthat)
library(ggplot2)
library(bayesassurance)


## bayes_sim_unknownvar() checks
test_that("Correct assurance produced", {
  n <- 100
  p <- 1
  sigsq <- 4
  epsilon <- 10e-7
  a_sig_d <- (sigsq / epsilon) + 2
  b_sig_d <- sigsq * (a_sig_d - 1)
  a_sig_a <- -p / 2
  b_sig_a <- 0
  
  set.seed(1234)
  out <- bayesassurance::bayes_sim_unknownvar(n = n, p = 1, 
         u = 1, C = 0, R = 70,
         Xn = NULL, Vn = NULL, Vbeta_d = 1e-8, 
         Vbeta_a_inv = 0, mu_beta_d = 0.25,
         mu_beta_a = 0, a_sig_a = a_sig_a, b_sig_a = b_sig_a, 
         a_sig_d = a_sig_d, b_sig_d = b_sig_d, 
         alt = "two.sided", alpha = 0.05, 
         mc_iter = 1000)

  expect_equal(out$assur_val, "Assurance: 0.214")
})
