library(bayesPO)
library(Rcpp)
gr <- (sqrt(5) + 1) / 2

#### Generating data function ####
cppFunction('NumericVector log1pexp_C(NumericVector x) {return log1p(exp(x));}')

genLogLink <- function(b, d){
  nb <- length(b)
  nd <- length(d)
  quitealot <- qnorm(1 - (1e-9) ^ (1/(nb - 1)))
  max_points <- sum(b[-1]) * quitealot + b[1] ## Should be an upper limit to the intensity
  
  ## Applying the Lewis-Shedler '79 algorithm for occurrences
  total <- rpois(1, exp(max_points))
  position_total <- cbind(runif(total), runif(total))
  x_total <- matrix(c( # Covariates
    rep(1, total), rnorm(total * (nb - 1))
  ), ncol = nb)
  thinningProbabilities <- x_total %*% b - max_points
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
b <- c(5, 2); d <- c(2, -2)
po_data <- genLogLink(b, d)
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
              initial_values = 1, joint_prior = prior(
                NormalPrior(rep(0, 2), 10 * diag(2)),
                NormalPrior(rep(0, 2), 10 * diag(2)),
                GammaPrior(1e-10, 1e-10)
              ))

estim = fit_bayesPO(object = po_model, background = cbind(x_background, w_background),
            area = 1, mcmc_setup = list(
              n_iter = 1e4, thin = 10, burnin = 2e5
            ))
# Saving results
saveRDS(estim, "Aux/loglin_single.rds")

#### Convergence diagnostics ####
library(bayesplot)
color_scheme_set("red")
estimarray <- as.array(estim)
mcmc_trace(estimarray)
# log Posterior seems to have converged

#### Plotting the intensity ####
library(ggplot2)
n_graph <- 500; xs <- seq(-1, 4, len = n_graph)
estim_intensities <- apply(as.matrix(estim), 1, \(params) exp(- log1pexp_C(-(
  params[1] + params[2] * xs
))) * params[5])
trueIntensity <- exp(cbind(rep(1, n_graph), xs) %*% b)
ggplot(
  data.frame(Covariate = xs,
             Intensity = trueIntensity),
  aes(Covariate, Intensity)
) + theme_bw() + geom_line() +
  geom_line(data = data.frame(Covariate = xs,
                              Intensity = apply(estim_intensities, 1, quantile, 0.025)),
            color = "gray") +
  geom_line(data = data.frame(Covariate = xs,
                              Intensity = rowMeans(estim_intensities)),
            color = "gray") +
  geom_line(data = data.frame(Covariate = xs,
                              Intensity = apply(estim_intensities, 1, quantile, 0.975)),
            color = "gray")

#### Do a lot of them ####
library(foreach)
library(doParallel)
library(ggplot2)
registerDoParallel(1)

seeds <- 1:20
b <- c(5, 2); d <- c(2, -2)
nb <- length(b); nd <- length(d)
n_graph <- 500; xs <- seq(-1, 4, len = n_graph)
ns <- rep(0, 20)

#manyIntensities <- foreach(i = 1:20) %dopar% {
manyIntensities <- foreach(i = c(17)) %dopar% {
  set.seed(seeds[i])
  po_data <- genLogLink(b, d)
  print(paste0("po points: ", nrow(po_data$position_observed)))
  print(paste0("unobserved points: ", nrow(po_data$position_unobserved)))
  
  # Need a matrix of background variables
  n_background <- 1e6
  x_background <- matrix(rnorm(n_background * (nb - 1)), ncol = nb - 1)
  w_background <- matrix(rnorm(n_background * (nd - 1)), ncol = nd - 1)
  
  # Model building
  po_model <- bayesPO_model(po = cbind(po_data$x_observed, po_data$w_observed),
                            intensitySelection = 1, observabilitySelection = 2,
                            intensityLink = "logit", observabilityLink = "logit",
                            initial_values = 1, joint_prior = prior(
                              NormalPrior(rep(0, 2), 10 * diag(2)),
                              NormalPrior(rep(0, 2), 10 * diag(2)),
                              GammaPrior(1e-20, 1e-20)
                            ))
  
  estim = fit_bayesPO(object = po_model, background = cbind(x_background, w_background),
                      area = 1, mcmc_setup = list(
                        iter = 1e4, thin = 10, burnin = 3e5
                      ))
  # Saving results
  nmb <- substr(format(i/100, nsmall = 2), 3,4)
  saveRDS(estim, paste0("Aux/loglin_multi",nmb,".rds"))
  
  estim_intensities <- apply(as.matrix(estim), 1, \(params) exp(- log1pexp_C(-(
    params[1] + params[2] * xs
  ))) * params[5])
  
  rowMeans(estim_intensities)
}

n_graph <- 500; xs <- seq(-1, 4, len = n_graph)
trueIntensity <- exp(cbind(rep(1, n_graph), xs) %*% b)
g <- ggplot(data.frame(Covariate = xs, Intensity = trueIntensity),
  aes(Covariate, Intensity)) + theme_bw() + geom_line()
for (i in 1:10)
  g <- g + geom_line(data = data.frame(Covariate = xs,
                                       Intensity = manyIntensities[[i]]),
                     color = "darkgrey")
g
ggsave("Aux/loglin_intensities.eps", plot = g,
       device = "eps", width = 8, height = 8 / gr)

for (i in 1:20) {
  set.seed(seeds[i])
  po_data <- genLogLink(b, d)
  ns[i] <- nrow(po_data$position_observed)
}

#### Recovering ####
library(ggplot2)

b <- c(5, 2); d <- c(2, -2)
n_graph <- 500; xs <- seq(-1, 3, len = n_graph)
trueIntensity <- exp(cbind(rep(1, n_graph), xs) %*% b)
g <- ggplot(data.frame(Covariate = xs, Intensity = trueIntensity),
            aes(Covariate, Intensity)) + theme_bw() + geom_line(size = 4) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))
grb <- ggplot(data.frame(Covariate = xs), aes(Covariate)) + theme_bw() +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))
rbMat <- matrix(0, 20, n_graph)
for (i in 1:20){
  nmb <- substr(format(i/100, nsmall = 2), 3, 4)
  estim <- readRDS(paste0("Aux/loglin_multi",nmb,".rds"))
  estim_intensities <- rowMeans(apply(as.matrix(estim), 1, \(params) exp(- log1pexp_C(-(
    params[1] + params[2] * xs
  ))) * params[5]))
  rb <- (estim_intensities - trueIntensity) / trueIntensity
  
  g <- g + geom_line(data = data.frame(Covariate = xs,
                                       Intensity = estim_intensities),
                     color = "darkgrey")
  grb <- grb + geom_line(aes(y = relBias), color = "darkgrey",
                         data = data.frame(Covariate = xs, relBias = rb))
  rbMat[i, ] <- rb
}
grb <- grb +
  geom_line(aes(y = rbInterval),
            data = data.frame(Covariate = xs, rbInterval = apply(rbMat, 2, quantile, 0.125))) +
  geom_line(aes(y = rbInterval),
            data = data.frame(Covariate = xs, rbInterval = apply(rbMat, 2, quantile, 0.875)))
g
grb
ggsave("Aux/loglin_intensities.eps", plot = g,
       device = "eps", width = 8, height = 8 / gr)
grb <- grb + labs(y = "Relative Bias")
ggsave(paste0("Aux/loglin_intensities_relBias.eps"), plot = grb,
       device = "eps", width = 8, height = 8 / gr)
