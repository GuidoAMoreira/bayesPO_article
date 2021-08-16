library(bayesPO)
library(Rcpp)
gr <- (sqrt(5) + 1) / 2

#### Generating data function ####
cppFunction('NumericVector log1pexp_C(NumericVector x) {return log1p(exp(x));}')

genLogitLink <- function(b, d, lambdaStar){
  nb <- length(b)
  nd <- length(d)
  
  ## Applying the Lewis-Shedler '79 algorithm for occurrences
  total <- rpois(1, lambdaStar)
  position_total <- cbind(runif(total), runif(total))
  x_total <- matrix(c( # Covariates
    rep(1, total), rnorm(total * (nb - 1))
  ), ncol = nb)
  thinningProbabilities <- - log1pexp_C(-x_total %*% b)
  thinning <- log(runif(total)) <= thinningProbabilities
  n_occurrences <- sum(thinning)
  position_occurrences <- position_total[thinning, ]
  x_occurrences <- x_total[thinning, ]
  
  ## Thinning the occurrences for the observations
  w_occurrences <- matrix(c( # Covariates
    rep(1, n_occurrences), rnorm(n_occurrences * (nd - 1))
  ), ncol = nd)
  thinningProbabilities <- - log1pexp_C(-w_occurrences %*% d)
  thinning <- log(runif(n_occurrences)) <= thinningProbabilities
  # n_observed <- sum(thinning)
  position_observed <- position_occurrences[thinning, ]
  x_observed <- x_occurrences[thinning, ]
  w_observed <- w_occurrences[thinning, ]
  # n_unobserved <- sum(!thinning)
  position_unobserved <- position_occurrences[!thinning, ]
  x_unobserved <- x_occurrences[!thinning, ]
  w_unobserved <- w_occurrences[!thinning, ]
  
  list(position_observed = position_observed, x_observed = x_observed[,-1],
       w_observed = w_observed[,-1], position_unobserved = position_unobserved,
       x_unobserved = x_unobserved[,-1], w_unobserved = w_unobserved[,-1])
}

#### Simulate ####
set.seed(123)
b <- c(-1, 2); d <- c(3, 4); ls <- 1e4
po_data <- genLogitLink(b, d, ls)
print(paste0("po points: ", nrow(po_data$position_observed)))
print(paste0("unobserved points: ", nrow(po_data$position_unobserved)))

#### Estimate ####
# Need a matrix of background variables
nb <- length(b); nd <- length(d)
n_background <- 1e6
x_background <- matrix(rnorm(n_background * (nb - 1)), ncol = nb - 1)
w_background <- matrix(rnorm(n_background * (nd - 1)), ncol = nd - 1)

# Model building
po_model <- bayesPO_model(po = cbind(po_data$x_observed, po_data$w_observed),
                          intensitySelection = 1, observabilitySelection = 2,
                          intensityLink = "logit", observabilityLink = "logit",
                          initial_values = 4, joint_prior = prior(
                            NormalPrior(rep(0, 2), 10 * diag(2)),
                            NormalPrior(rep(0, 2), 10 * diag(2)),
                            GammaPrior(1e-4, 1e-4)
                          ))

estim = fit_bayesPO(object = po_model, background = cbind(x_background, w_background),
                    area = 1, mcmc_setup = list(
                      iter = 5e2, thin = 1, burnin = 0
                    ))

#### Convergence diagnostics ####
library(bayesplot)
color_scheme_set("green")
color_scheme_set("gray")
estimarray <- as.array(estim)
g <- mcmc_trace(estimarray)
g
ggplot2::ggsave(paste0("Aux/traceplots.pdf"), plot = g,
       device = "pdf", width = 8, height = 8 / gr)
# log Posterior seems to have converged

estim_plus = fit_bayesPO(estim, background = cbind(x_background, w_background),
                         mcmc_setup = list(
                           iter = 1e3, thin = 1, burnin = 0
                         ))
