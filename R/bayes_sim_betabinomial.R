#' Bayesian Assurance Computation in the Beta-Binomial Setting
#'
#' Returns the Bayesian assurance corresponding to a hypothesis test for
#' difference in two independent proportions.
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline aes theme element_text xlab ylab ggtitle scale_x_continuous scale_y_continuous
#' @importFrom rlang .data
#' @importFrom pbapply pbsapply
#' @param n1 sample size of first group
#' @param n2 sample size of second group
#' @param p1 proportion of successes in first group. Takes on a NULL (default) assignment
#' if unknown.
#' @param p2 proportion of successes in second group. Takes on a NULL (default) assignment
#' if unknown.
#' @param alpha_1,beta_1 shape parameters for the distribution of p1 if p1 is not known:
#' \eqn{p1 ~ Beta(\alpha_1, \beta_1)}
#' @param alpha_2,beta_2 shape parameters for the distribution of p2 if p2 is not known:
#' \eqn{p2 ~ Beta(\alpha_2, \beta_2)}
#' @param sig_level significance level
#' @param alt a character string specifying the alternative hypothesis, must be one of
#' "two.sided" (default), "greater" or "less".
#' @param mc_iter number of MC samples evaluated under the analysis objective
#' @returns approximate Bayesian assurance of independent two-sample proportion test
#' @examples
#'
#' # Simulation that overlays frequentist power with Bayesian assurance
#' # for independent two-sample proportion test.
#'
#' library(ggplot2)
#'
#' #########################################################
#' # alpha1 = 0.5, beta1 = 0.5, alpha2 = 0.5, beta2 = 0.5 ##
#' #########################################################
#' n <- seq(40, 1000, 10)
#'
#' power_vals <- sapply(n, function(i) bayesassurance::propdiffCI_classic(i, p1 = 0.25, p2 = 0.2, alpha_1 = NULL,
#' beta_1 = NULL, alpha_2 = NULL, beta_2 = NULL, sig_level = 0.05))
#'
#' df1 <- as.data.frame(cbind(n, power_vals))
#'
#' assur_vals <- sapply(n, function(i) bayesassurance::bayes_sim_betabin(n1 = i, n2 = i, p1 = 0.25, p2 = 0.2, alpha_1 = 0.5, beta_1 = 0.5,
#'                                         alpha_2 = 0.5, beta_2 = 0.5, sig_level = 0.05, alt = "two.sided"))
#'
#' df2 <- as.data.frame(cbind(n, assur_vals))
#'
#' p1 <- ggplot(df1, alpha = 0.5, aes(x = n, y = power_vals, color="Frequentist Power"))
#' p1 <- p1 + geom_line(alpha = 0.5, aes(x = n, y = power_vals, color="Frequentist Power"), lwd = 1.2)
#' p2 <- p1 + geom_point(data = df2, alpha = 0.5, aes(x = n, y = assur_vals, color = "Bayesian Assurance"),
#' lwd = 1.2) + ylab("Power/Assurance") + xlab(~ paste("Sample Size n = ", "n"[1], " = ", "n"[2])) +
#' ggtitle("Power/Assurance Curves of Difference in Proportions")
#' p2 <- p2 + geom_label(aes(525, 0.6,
#'                          label = "~p[1]-p[2] == 0.05"), parse = TRUE,
#'                          color = "black", size = 3)
#'
#' @export
#'

bayes_sim_betabin <- function(n1, n2, p1, p2, alpha_1, alpha_2, beta_1, beta_2, sig_level, alt, mc_iter = 5000){
    count <- 0

    # set.seed(1)
    if(is.null(p1) == TRUE & is.null(p2) == TRUE){
      p1 <- rbeta(n=1, alpha_1, beta_1)
      p2 <- rbeta(n=1, alpha_2, beta_2)
    }else if(is.null(p1) == TRUE & is.null(p2) == FALSE){
      p1 <- rbeta(n=1, alpha_1, beta_1)
    }else if(is.null(p1) == FALSE & is.null(p2) == TRUE){
      p2 <- rbeta(n=1, alpha_2, beta_2)
    }

    MC_samp <- function(n1 = n1, n2 = n2){
      for(i in 1:mc_iter){

        # x1 and x2 are observed frequency values, simulate by drawing from the binomial distribution
        x1 <- rbinom(n=1, size=n1, prob=p1)
        x2 <- rbinom(n=1, size=n2, prob=p2)

        # estimating difference of proportions: p = p1 - p2
        p_post_mean <- ((alpha_1 + x1) / (alpha_1 + beta_1 + n1)) - ((alpha_2 + x2) / (alpha_2 + beta_2 + n2))
        p_post_var <- (((alpha_1 + x1) * (beta_1 + n1 - x1)) / ((alpha_1 + beta_1 + n1)^2 * (alpha_1 + beta_1 + n1 + 1))) +
          (((alpha_2 + x2) * (beta_2 + n2 - x2)) / ((alpha_2 + beta_2 + n2)^2 * (alpha_2 + beta_2 + n2 + 1)))

        # critical z-value of normal distribution
        if(alt == "two.sided"){
          z <- qnorm(1 - sig_level / 2)
        }else if(alt == "greater" | "less"){
          z <- qnorm(1 - sig_level)
        }

        # Condition to check if 0 falls in 100(1-sig_level)% Credible Interval
        if(alt == "two.sided"){
          # lower two-sided confidence bound
          lb <- p_post_mean - z * (sqrt(p_post_var))
          # upper two-sided confidence bound
          ub <- p_post_mean + z * (sqrt(p_post_var))
          if(0 < lb | 0 > ub){
            count <- count + 1
          }else{
            count <- count
          }
        }else if(alt == "less"){ # Condition to check if 0 falls above UB of 100(1-sig_level)% Credible Interval
          # lower two-sided confidence bound
          lb <- -Inf
          # upper two-sided confidence bound
          ub <- p_post_mean + z * (sqrt(p_post_var))
          if(0 > ub){
            count <- count + 1
          }else{
            count <- count
          }
        }else if(alt == "greater"){
          # lower two-sided confidence bound
          lb <- p_post_mean - z * (sqrt(p_post_var))
          # upper two-sided confidence bound
          ub <- Inf
          if(0 < lb){
            count <- count + 1
          }else{
            count <- count
          }
        }
      }
      assurance <- count / mc_iter
      return(assurance)
    }

    # objects returned if n1 and n2 are vectors
    if(length(n1) > 1 & length(n2) > 1){
      assurance_vals <- pbapply::pbmapply(MC_samp, n1 = n1, n2 = n2)
      assur_tab <- as.data.frame(cbind(n1, n2, assurance_vals))
      colnames(assur_tab) <- c("n1", "n2", "Assurance")

      # returns an assurance plot only if vectors n1 and n2 are identical, enabling
      # the plot to focus on a single n on the x-axis
      if(identical(n1, n2)){
        assur_plot <- ggplot2::ggplot(assur_tab, alpha = 0.5, aes(x = .data$`n1`,
                                                                  y = .data$Assurance)) +
          ggplot2::geom_point(aes(x = .data$`n1`, y = .data$Assurance), lwd = 1.2) +
          ggplot2::ggtitle("Assurance Values") +
          ggplot2::xlab("Sample Size n = n1 = n2") + ggplot2::ylab("Assurance")
        assur_plot <- structure(assur_plot, class = "ggplot")
        return(list(assurance_table = assur_tab, assurance_plot = assur_plot, mc_samples = mc_iter))
      }else{
        return(list(assurance_table = assur_tab, mc_samples = mc_iter))
      }
    }

    # objects returned when a scalar is passed in for n1 and n2
    return(list(assur_val = paste0("Assurance: ", round(assurance, 3)), mc_samples = mc_iter))

}





