####### bayesassurance Replication Script ########
## Examples shown in order of paper

## install and load the required packages
# install.packages("ggplot2")
# install.packages("MASS")
# install.packages("pbapply")
# install.packages("dplyr")
# install.packages("plot3D")
# install.packages("plotly")
# install.paackages("latex2exp")
# library(ggplot2)
# library(MASS)
# library(pbapply)
# library(dplyr)
# library(plot3D)
# library(plotly)
# library(latex2exp)


library(bayesassurance)
## Page 5
## Reproduces Figure 1
n <- seq(100, 250, 10)
n_a <- 10
n_d <- 10
theta_0 <- 0.15
theta_1 <- 0.25
sigsq <- 0.30

out <- assurance_nd_na(n = n, n_a = n_a,  n_d = n_d, theta_0 = theta_0, theta_1 = theta_1, sigsq = sigsq,
                          alt = "greater", alpha = 0.05)

head(out$assurance_table)
out$assurance_plot





## Page 6
## Reproduces Figure 2b
n <- seq(10, 250, 5)
n_a <- 1e-8 # large precision parameter in analysis stage signifies weak analysis prior
n_d <- 1e+8
theta_0 <- 0.15
theta_1 <- 0.25
sigsq <- 0.104

out <- assurance_nd_na(n = n, n_a = n_a,  n_d = n_d, theta_0 = theta_0, theta_1 = theta_1, sigsq = sigsq,
                       alt = "greater", alpha = 0.05)

head(out$assurance_table)
out$assurance_plot





## Page 8
pwr_freq(n = 20, theta_0 = 0.15, theta_1 = 0.35, sigsq = 0.30, alt = "greater", alpha = 0.05)





## Page 8
## Reproduces Figure 2a
n <- seq(10, 250, 5)
out <- pwr_freq(n = n, theta_0 = 0.15, theta_1 = 0.25, sigsq = 0.104, alt = "greater", alpha = 0.05)

out$pwr_table
out$pwr_plot





## Page 11
## Example 1 of Paper
n <- seq(100, 300, 10)

assur_vals <- bayesassurance::bayes_sim(n, p = 1, u = 1,
                                        C = 0.15, Xn = NULL, Vbeta_d = 0, Vbeta_a_inv = 0,
                                        Vn = NULL, sigsq = 0.265, mu_beta_d = 0.25, mu_beta_a = 0,
                                        alt = "greater", alpha = 0.05, mc_iter = 5000)
assur_vals$assurance_table
assur_vals$assurance_plot





## Page 13
## Example 2 of Paper (Entries taken from O'Hagan/Stevens 2001)
## Note the progress bar doesn't show up for scalar entries of n, plan to update
## this in future edits

n <- 285
p <- 4
K <- 20000 # threshold unit cost
C <- 0
u <- as.matrix(c(-K, 1, K, -1))
sigsq <- 4.04^2

## Assign mean parameters to analysis and design stage priors
mu_beta_d <- as.matrix(c(5, 6000, 6.5, 7200))
mu_beta_a <- as.matrix(rep(0, p))

## Assign correlation matrices (specified in paper)
## to analysis and design stage priors
Vbeta_a_inv <- matrix(rep(0, p^2), nrow = p, ncol = p)
Vbeta_d <- (1 / sigsq) * matrix(c(4, 0, 3, 0, 0, 10^7, 0,
                                     0, 3, 0, 4, 0, 0, 0, 0, 10^7), nrow = 4, ncol = 4)

tau1 <- tau2 <- 8700
sig <- sqrt(sigsq)
Vn <- matrix(0, nrow = n*p, ncol = n*p)
Vn[1:n, 1:n] <- diag(n)
Vn[(2*n - (n-1)):(2*n), (2*n - (n-1)):(2*n)] <- (tau1 / sig)^2 * diag(n)
Vn[(3*n - (n-1)):(3*n), (3*n - (n-1)):(3*n)] <- diag(n)
Vn[(4*n - (n-1)):(4*n), (4*n - (n-1)):(4*n)] <- (tau2 / sig)^2 * diag(n)

set.seed(10)
assur_val <- bayes_sim(n = 285, p = 4,  u = as.matrix(c(-K, 1, K, -1)),
                           C = 0, Xn = NULL, Vbeta_d = Vbeta_d, Vbeta_a_inv = Vbeta_a_inv,
                           Vn = Vn, sigsq = 4.04^2, mu_beta_d = as.matrix(c(5, 6000, 6.5, 7200)),
                           mu_beta_a = as.matrix(rep(0, p)), alt = "greater", alpha = 0.05, mc_iter = 500)

assur_val





## Page 16
## Example 3 of Paper
## Reproduces Figure 3
n <- seq(10, 100, 5)
ids <- c(1,2)
Vbeta_a_inv <- matrix(rep(0, 16), nrow = 4, ncol = 4)
sigsq = 100
Vbeta_d <- (1 / sigsq) * matrix(c(4, 0, 3, 0, 0, 6, 0, 0, 3, 0, 4, 0, 0, 0, 0, 6),
                                   nrow = 4, ncol = 4)

set.seed(12)
assur_out <- bayes_sim(n = n, p = NULL, u = c(1, -1, 1, -1), C = 0, Xn = NULL,
                          Vbeta_d = Vbeta_d, Vbeta_a_inv = Vbeta_a_inv,
                          Vn = NULL, sigsq = 100,
                          mu_beta_d = as.matrix(c(5, 6.5, 62, 84)),
                          mu_beta_a = as.matrix(rep(0, 4)), mc_iter = 5000,
                          alt = "two.sided", alpha = 0.05, longitudinal = TRUE, ids = ids,
                          from = 10, to = 120)

head(assur_out$assurance_table)
assur_out$assurance_plot





## Page 18
## Example 4 of Paper
## Note this code block takes some time to run as we are randomly generating values
## for two sets of unknown parameters.
n <- 285
p <- 4
K <- 20000
C <- 0
u <- as.matrix(c(-K, 1, K, -1))
Vbeta_a_inv <- matrix(rep(0, p^2), nrow = p, ncol = p)
sigsq <- 4.04^2
tau1 <- tau2 <- 8700

## manually-assigned correlation matrices
Vn <- matrix(0, nrow = n*p, ncol = n*p)
Vn[1:n, 1:n] <- diag(n)
Vn[(2*n - (n-1)):(2*n), (2*n - (n-1)):(2*n)] <- (tau1 / sig)^2 * diag(n)
Vn[(3*n - (n-1)):(3*n), (3*n - (n-1)):(3*n)] <- diag(n)
Vn[(4*n - (n-1)):(4*n), (4*n - (n-1)):(4*n)] <- (tau2 / sig)^2 * diag(n)

Vbeta_d <- (1 / sigsq) * matrix(c(4, 0, 3, 0, 0, 10^7, 0, 0, 3, 0, 4, 0,
                                     0, 0, 0, 10^7),
                                   nrow = 4, ncol = 4)

mu_beta_d <- as.matrix(c(5, 6000, 6.5, 7200))
mu_beta_a <- as.matrix(rep(0, p))
alpha <- 0.05
epsilon <- 10e-7
a_sig_d <- (sigsq / epsilon) + 2
b_sig_d <- sigsq * (a_sig_d - 1)
a_sig_a <- -p / 2
b_sig_a <- 0
R <- 150

set.seed(60)
assur_out <- bayes_sim_unknownvar(n = n, p = 4, u = as.matrix(c(-K, 1, K, -1)),
                                     C = 0, R = 150,Xn = NULL, Vn = Vn, Vbeta_d = Vbeta_d,
                                     Vbeta_a_inv = Vbeta_a_inv, mu_beta_d = mu_beta_d,
                                     mu_beta_a = mu_beta_a, a_sig_a = a_sig_a,
                                     b_sig_a = b_sig_a, a_sig_d = a_sig_d,b_sig_d = b_sig_d,
                                     alt = "two.sided", alpha = 0.05, mc_iter = 5000)
assur_out$assur_val





## Page 20
## Example 5 of Paper
## Reproduces Figure 4
## Note this code block takes some time to run as the code is checking across all
## combinations of n1 and n2 to produce contour plot.
n1 <- seq(20, 75, 5)
n2 <- seq(50, 160, 10)

set.seed(3)
assur_out <- bayes_sim_unbalanced(n1 = n1, n2 = n2, repeats = 1, u = c(1, -1),
                                     C = 0, Xn = NULL, Vbeta_d = matrix(c(50, 0, 0, 10),nrow = 2, ncol = 2),
                                     Vbeta_a_inv = matrix(rep(0, 4), nrow = 2, ncol = 2),
                                     Vn = NULL, sigsq = 100,  mu_beta_d = c(1.17, 1.25),
                                     mu_beta_a = c(0, 0), alt = "two.sided", alpha = 0.05, mc_iter = 5000,
                                     surface_plot = TRUE)

head(assur_out$assurance_table)
assur_out$contourplot





## Page 21
## Example 6 of Paper
## Reproduces Figure 5
n1 <- c(4, 5, 15, 25, 30, 100, 200)
n2 <- c(8, 10, 20, 40, 50, 200, 250)

mu_beta_d <- as.matrix(c(5, 6000, 6.5, 7200))
mu_beta_a <- as.matrix(rep(0, 4))
K = 20000 # threshold unit cost
C <- 0
u <- as.matrix(c(-K, 1, K, -1))
sigsq <- 4.04^2
Vbeta_a_inv <- matrix(rep(0, 16), nrow = 4, ncol = 4)
Vbeta_d <- (1 / sigsq) * matrix(c(4, 0, 3, 0, 0, 10^7, 0, 0,
                                     3, 0, 4, 0, 0, 0, 0, 10^7),nrow = 4, ncol = 4)

set.seed(3)
assur_out <- bayes_sim_unbalanced(n1 = n1, n2 = n2, repeats = 2,
                                     u = as.matrix(c(-K, 1, K, -1)), C = 0, Xn = NULL,
                                     Vbeta_d = Vbeta_d, Vbeta_a_inv = Vbeta_a_inv,
                                     Vn = NULL, sigsq = 4.04^2,
                                     mu_beta_d = as.matrix(c(5, 6000, 6.5, 7200)),
                                     mu_beta_a = as.matrix(rep(0, 4)),
                                     alt = "greater", alpha = 0.05, mc_iter = 5000,
                                     surface_plot = TRUE)

assur_out$assurance_table
assur_out$contourplot





## Page 24
## Example 7 of Paper
## Reproduces Figure 6

n <- seq(20, 145, 5)
set.seed(10)
out <- bayes_adcock(n = n, d = 0.20, mu_beta_a = 0.64, mu_beta_d = 0.9,
                       n_a = 20, n_d = 10, sig_sq = 0.265,
                       alpha = 0.05, mc_iter = 10000)
head(out$assurance_table)
out$assurance_plot





## Page 25
## Reproduces Figure 7
## Example 8 of Paper
## Please refer to fig7_replication.Rmd for step-by-step walk-through





## Page 28
## Example 9 of Paper
## Reproduces Figure 8
n <- seq(800, 2000, 50)

set.seed(30)
out <- bayes_sim_betabin(n1 = n, n2 = n, p1 = 0.25, p2 = 0.2,
                            alpha_1 = 0.5, beta_1 = 0.5, alpha_2 = 0.5,
                            beta_2 = 0.5, sig_level = 0.05, alt = "two.sided")

out$assurance_table
out$assurance_plot





## Page 29
## Example 10 of Paper
## Reproduces Figure 9
## Please refer to fig9_replication.Rmd for step-by-step walk-through





## Page 31
## Example 11 of Paper
## Reproduces Figure 10
## Please refer to fig10_replication.Rmd for step-by-step walk-through




## Page 35
## Example 12 of Paper
## Reproduces Figure 11
n <- seq(100, 1200, 10)

set.seed(3)
out <- bayes_goal_func(n, Xn = NULL, K = 1, pi = 0.5, u = 1, sigsq = 1,
                          beta_0 = 0.5, beta_1 = 0.6)

head(out$rc_table)
out$rc_plot

## adding labels to ggplot
out$rc_plot + ggplot2::geom_hline(yintercept = 0.9283, linetype = "dashed",
  color = "red", size = 0.45) + geom_vline(xintercept = 857,
  linetype = "dashed", color = "red", size = 0.45) +
  scale_y_continuous(breaks = sort(c(seq(0.7, 1, 0.1), 0.9283))) +
  scale_x_continuous(breaks = sort(c(seq(200, 1250, 250), 857)))





## Page 36
## Reproduces Figure 12
library(ggplot2)
library(latex2exp)

# Function that returns the minimum rate of correct classification r_star
rstar <- function(n_b, sig_sq, K, pi, delta){
  val1 <- (sqrt(sig_sq) * log((K*pi)/(1-pi)))/(sqrt(n_b)*delta) + (sqrt(n_b)*delta)/(2*sqrt(sig_sq))
  val2 <- (sqrt(sig_sq) * log((K*pi)/(1-pi)))/(sqrt(n_b)*delta) - (sqrt(n_b)*delta)/(2*sqrt(sig_sq))
  r_star <- K*pi*pnorm(val1) + (1-pi)*(1-pnorm(val2))
  return(r_star)
}

n <- seq(300, 10000, 100)
r_star <- rstar(n_b = n, sig_sq = 1, K = 1, pi = 0.5, delta = 0.05)
r_star1 <- rstar(n_b = n, sig_sq = 1, K = 1, pi = 0.5, delta = 0.10)
r_star2 <- rstar(n_b = n, sig_sq = 1, K = 1, pi = 0.5, delta = 0.01)
r_star3 <- rstar(n_b = n, sig_sq = 1, K = 1, pi = 0.5, delta = 0.03)
df_bayes <- as.data.frame(cbind(r_star, n))
df_bayes1 <- as.data.frame(cbind(r_star1, n))
df_bayes2 <- as.data.frame(cbind(r_star2, n))
df_bayes3 <- as.data.frame(cbind(r_star3, n))

## RUN THE FOLLOWING SEGMENTS ONE AT A TIME (separated by spaces)
p <- ggplot(df_bayes, aes(x = n, y = r_star)) + geom_line(lwd=1.1)

p1 <- p + xlab(TeX("Sample Size $n_b$")) + ylab("r*") +
  ggtitle(TeX("Bayesian Goal Function Curve for Varying $\\delta$"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(x=2000, y=0.95, label="delta = 0.10", size = 3.2)

p2 <- p1 + geom_line(data = df_bayes1, aes(x = n, y = r_star1), lwd = 1.1) +
  geom_text(x=2700, y=0.85, label="delta = 0.05", size = 3.2)

p_bayes <- p2 + geom_line(data = df_bayes3, aes(x = n, y = r_star3), lwd = 1.1) +
  geom_text(x=4600, y=0.80, label="delta = 0.03", size = 3.2)

p_bayes + geom_hline(yintercept = .9283, linetype = "dashed", color = "red", size = 0.45) +
  geom_vline(xintercept = 857, linetype = "dashed", color = "red", size = 0.45) +
  geom_vline(xintercept = 3426, linetype = "dashed", color = "red", size = 0.45) +
  geom_vline(xintercept = 9512, linetype = "dashed", color = "red", size = 0.45) +
  scale_y_continuous(breaks = sort(c(round(seq(0.6, 1.0, 0.1), 2), .9283)))




## Page 39
## Reproduces Figure 13
nf_formula <- function(alpha, beta, sig_sq, delta){
  n_f <- ((qnorm(alpha) + qnorm(beta))^2)*(sig_sq/delta)^2
  return(n_f)
}

# (power = 1 - beta_vals)
beta_vals <- seq(0.01, 0.93, 0.01)

# identifies specific sample sizes
# in frequentist setting (denoted as n_f) that meet
# the power criteria listed above
n_f <- sapply(beta_vals, nf_formula, alpha = 0.05, sig_sq = 1, delta = 0.1)

# computes corresponding rate of correct classifications (r*) for each
# of the sample sizes in n_f
r_star <- bayesassurance::bayes_goal_func(n = n_f, sigsq = 1, K = 1, pi = 0.5,
                                          u = 1, beta_0 = 0.5,beta_1 = 0.6)

# r* vs beta
library(ggplot2)
library(latex2exp)

betas <- seq(0.01, 0.93, 0.01)
n_f <- sapply(betas, nf_formula, alpha = 0.05, sig_sq = 1, delta = 0.1)
r_star <- bayesassurance::bayes_goal_func(n = n_f, sigsq = 1, u = 1, K = 1, pi = 0.5, beta_0 = 0.5,
                                          beta_1 = 0.6)
df <- as.data.frame(cbind(n_f, r_star$rc_table$`Rate of Correct Classification`, betas))
colnames(df) <- c("n_f", "r_star", "betas")
p <- ggplot(df, aes(x = 1-betas, y = r_star)) + geom_line(lwd=1.2)
r_star_beta_curve <- p + xlab(TeX("Power: $1-\\beta$")) + ylab("Rate of Correct Classification: r*") +
  labs(title = "r* and Power Values that Yield Equal Sample Sizes \n in Bayesian and Frequentist Setting ")+
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(x = 0.84, y = 0.94, label = "(0.90, 0.9283)", size = 3.2) +
  geom_text(x = 0.91, y = 0.97, label = "(0.95, 0.95)", size = 3.2) +
  geom_text(x = 0.68, y = 0.72, label = "y = x", size = 3.2) +
  geom_point(aes(x=0.90, y=0.9283), colour="blue") +
  geom_point(aes(x=0.95, y=0.95), colour="blue") +
  geom_abline(slope=1, intercept = 0, linetype = "dashed") + xlim(c(0.4, 1)) + ylim(c(0.4, 1))
r_star_beta_curve

df[c(5, 17, 35, 49, 72, 81, 87),]




## Page 40
## Reproduces Figure 14
library(bayesassurance)

out <- pwr_curve(n = seq(10, 200, 10), n_a = 1e-8, n_d = 1e+8,
                    sigsq = 0.104, theta_0 = 0.15,theta_1 = 0.25,
                    alt = "greater", alpha = 0.05,
                    bayes_sim = TRUE, mc_iter = 5000)

out$power_table
out$assurance_table
out$plot





## Page 42
## Reproduces Figure 15
library(bayesassurance)

out <- pwr_curve(n = seq(10, 200, 10), n_a = 1e-8, n_d = 1e-8,
             sigsq = 0.104, theta_0 = 0.15,theta_1 = 0.25,
             alt = "greater", alpha = 0.05,
             bayes_sim = TRUE, mc_iter = 5000)

out$plot




## Page 42
## Reproduces standard design matrix
n <- c(1,3,5,8)
gen_Xn(n = n)





## Page 44
## Reproduces longitudinal design matrix
ids <- c(1,2,3,4)
gen_Xn_longitudinal(ids, from = 1, to = 10, num_repeated_measures = 4)
