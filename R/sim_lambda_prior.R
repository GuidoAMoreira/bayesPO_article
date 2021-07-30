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

#### Do a lot of them ####
library(ggplot2)
theme_set(theme_bw())

set.seed(123)
b <- c(-1, 2); d <- c(3, 4); ls <- 1e4
nb <- length(b); nd <- length(d)
n_graph <- 500; xs <- seq(-3, 3, len = n_graph)
ns <- matrix(0, 4, 10)
results <- matrix(0, 1e4, 10)
lambdas <- matrix(0, 4, 10)

for(j in 1:4){
  for(i in 1:10){
    print(paste("Starting chain",i))
    po_data <- genLogitLink(b, d, ls)
    print(paste0("po points: ", nrow(po_data$position_observed)))
    print(paste0("unobserved points: ", nrow(po_data$position_unobserved)))
    ns[1, i] <- nrow(po_data$position_observed)
    
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
                                GammaPrior(ls * 10^(-j - 5), ls * 10^(-j - 5))
                              ))
    
    estim = fit_bayesPO(object = po_model, background = cbind(x_background, w_background),
                        area = 1, mcmc_setup = list(
                          n_iter = 2e4, thin = 2, burnin = 1e4
                        ))
    # Saving results
    nmb <- substr(format(i/100, nsmall = 2), 3, 4)
    saveRDS(estim, paste0("Aux/lambdaprior_multi",j,"_",nmb,".rds"))
    #estim <- readRDS(paste0("Aux/clogloglin_multi",nmb,".rds"))
    
    estim_intensities <- apply(as.matrix(estim), 1, \(params) exp(- log1pexp_C(-(
      params[1] + params[2] * xs
    ))) * params[5])
    
    results[, i] <- rowMeans(estim_intensities)
    lambdas[j, i] <- estim$parameters[5]
  }
  
  n_graph <- 500; xs <- seq(-3, 3, len = n_graph)
  trueIntensity <- exp(- log1pexp_C(-cbind(rep(1, n_graph), xs) %*% b)) * ls
  g <- ggplot(data.frame(Covariate = xs, Intensity = trueIntensity),
              aes(Covariate, Intensity)) + theme_bw() + geom_line()
  for (i in 1:10)
    g <- g + geom_line(data = data.frame(Covariate = xs,
                                         Intensity = results[, i]),
                       color = "darkgrey")
  g
  ggsave(paste0("Aux/mlambdaprior",j,"IntM.eps"), plot = g,
         device = "eps", width = 8, height = 8 / gr)
}

g <- ggplot(NULL) + geom_vline(xintercept = ls) + xlim(7500, 12500) + xlab(expression(lambda*"*")) +
  geom_point(aes(x = lambda, y = -1), data = data.frame(lambda = lambdas[1, ])) +
  geom_point(aes(x = lambda, y = -2), data = data.frame(lambda = lambdas[2, ])) +
  geom_point(aes(x = lambda, y = -3), data = data.frame(lambda = lambdas[3, ])) +
  geom_point(aes(x = lambda, y = -4), data = data.frame(lambda = lambdas[4, ])) +
  scale_y_continuous(name = "Prior variance",
                   breaks = -1:-4,
                   labels = c("100","1,000","10,000","100,000"))
g
ggsave(paste0("Aux/lambdapriorParameter.eps"),
       plot = g, device = "eps", width = 8, height = 8 / gr)