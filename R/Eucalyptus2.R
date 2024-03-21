library(dplyr)
library(bayesPO)
library(bayesplot)
color_scheme_set("green")


#### Data loading ####
load("Aux/Eucalyptus/Quad100m.RData", verbose = TRUE)
load("Aux/Eucalyptus/Eucalyptus sparsifolia Atlas 2012.RData", verbose = TRUE)
load("Aux/Eucalyptus/Euc Spars Soil.RData", verbose = TRUE)

#### Data preparation ####
quad <- quad %>% mutate(
  soil5 = (soil == 5) * 1, soil7 = (soil == 7) * 1, soil8 = (soil == 8) * 1,
  soil11 = (soil == 11) * 1, soil12 = (soil == 12) * 1, FC2 = FC^2, MNT2 = MNT^2,
  MXT2 = MXT^2, Rain2 = Rain^2, FC_MNT = FC * MNT, FC_MXT = FC * MXT,
  FC_Rain = FC * Rain, MNT_MXT = MNT * MXT, MNT_Rain = MNT * Rain,
  MXT_Rain = MXT * Rain, D.Main2 = D.Main^2, D.Urb2 = D.Urb^2,
  D.Main_D.Urb = D.Main * D.Urb, soil = as.numeric(soil))
means <- colMeans(quad[-grep("soil", names(quad))])
sds <- apply(quad[-grep("soil", names(quad))], 2, sd)
for (i in (1:ncol(quad))[-grep("soil", names(quad))]) {
  findStats <- which(names(means) == names(quad)[i])
  quad[, i] <- (quad[, i] - means[findStats]) / sds[findStats]
}
quad <- as.matrix(quad)

po_mat <- cbind(D.Main, D.Urb, FC, MNT, MXT, Rain,
                sp.soil == 5, sp.soil == 7, sp.soil == 8,
                sp.soil == 11, sp.soil == 12,
                D.Main^2, D.Urb^2, FC^2, MNT^2, MXT^2, Rain^2,
                D.Main * D.Urb, FC * MNT, FC * MXT, FC * Rain,
                MNT * MXT, MNT * Rain, MXT * Rain)
colnames(po_mat) <- c(
  "D.Main", "D.Urb", "FC", "MNT", "MXT", "Rain",
  "soil5", "soil7", "soil8", "soil11", "soil12",
  "D.Main2", "D.Urb2", "FC2", "MNT2", "MXT2", "Rain2",
  "D.Main_D.Urb", "FC_MNT", "FC_MXT", "FC_Rain",
  "MNT_MXT", "MNT_Rain", "MXT_Rain"
)
for (i in (1:ncol(po_mat))[-grep("soil", colnames(po_mat))]) {
  findStats <- which(names(means) == colnames(po_mat)[i])
  po_mat[, i] <- (po_mat[, i] - means[findStats]) / sds[findStats]
}

#### Fitting model ####

# Model definition
po <- bayesPO_model(
  po = po_mat,
  intensitySelection = colnames(po_mat)[-c(1, 2, 12, 13, 18)],
  observabilitySelection = colnames(po_mat)[c(1, 2, 12, 13, 18)],
  intensityLink = "logit", observabilityLink = "logit",
  initial_values = 3, joint_prior = prior(
    NormalPrior(rep(0, 20), 10 * diag(20)),
    NormalPrior(rep(0, 6), 10 * diag(6)),
    GammaPrior(1e-4, 1e-4)
  )
)

# Model fit
# High thin value chosen due tu high autocorrelation in some individual chains
model_fit <- fit_bayesPO(
  object = po, background = quad, area = nrow(quad),
  mcmc_setup = list(iter = 1e6, burnin = 1e5, thin = 100)
)

saveRDS(model_fit, "Aux/Eucalyptus/model_fit")
# model_fit <- readRDS("Aux/Eucalyptus/model_fit")

# Diagnosing
model_fit
mcmc_trace(model_fit)
mcmc_dens(model_fit)

model_fit$area

#### Lattice prediction ####
library(raster)
library(doParallel)
library(foreach)
registerDoParallel(10)
model_fit <- readRDS("Aux/Eucalyptus/model_fit")

sampleXp <- function(covs, params, sampleFrom) {
  intensityVars <- model_fit@original@iSelectedColumns
  observabilityVars <- model_fit@original@oSelectedColumns
  delta0_position <- which(names(params) == "Observability_Intercept")
  covs[is.na(covs)] <- 0
  calYsize <- params$n_Xprime
  covsInt <- do.call(c, lapply(intensityVars, function(n) which(n == names(covs))))
  paramsInt <- do.call(c, lapply(intensityVars, function(n) which(n == names(params))))
  covsObs <- do.call(c, lapply(observabilityVars, function(n) which(n == names(covs))))
  paramsObs <- do.call(c, lapply(observabilityVars, function(n) which(n == names(params))))
  
  accepted <- NULL
  while (length(accepted) < calYsize) {
    candidate <- sample(sampleFrom, 1)
    
    if (log(runif(1)) <=
        -log1p(exp(as.matrix(-cbind(1, covs[candidate, covsInt])) %*%
                   t(as.matrix(params[, c(1, paramsInt)])))) -
        log1p(exp(as.matrix(cbind(1, covs[candidate, covsObs])) %*%
                  t(as.matrix(params[, c(delta0_position, paramsObs)])))) )
      accepted <- c(accepted, candidate)
  }
  accepted
}

coords <- cbind(as.matrix(quad[, 1:2]), 1); colnames(coords) <- LETTERS[24:26]
df_fit <- as.data.frame(model_fit)
tempRaster <- rasterFromXYZ(coords)
sampleable <- which(!is.na(values(tempRaster)))

set.seed(123)
counts <- foreach(i = 1:nrow(df_fit)) %dopar% {
#for (i in 1:nrow(df_fit)) {
  sample <- 1:length(sampleable) %in% sampleXp(quad, df_fit[i, ], 1:length(sampleable))
  values(tempRaster)[sampleable] <- sample
  values(raster::aggregate(tempRaster, fact = 25, fun = max))
}
fullCounts25 <- rep(0, length(counts[[1]]))
for (i in 1:nrow(df_fit)) fullCounts25 <- fullCounts25 + counts[[i]]

saveRDS(fullCounts25, "Aux/Eucalyptus/counts25.rds")
# fullCounts25 <- readRDS("Aux/Eucalyptus/counts25.rds")

predicted <- raster::aggregate(tempRaster, fact = 25, fun = max)
values(predicted) <- fullCounts25 / nrow(df_fit)
plot(predicted)


mat <- cbind(coords, 1)
temp <- rowSums(as.data.frame(mini))
for (i in 1:nrow(mat))
  if (is.na(temp[i])) mat[i, 3] <- NA
colnames(mat) <- LETTERS[24:26]

mat10 <- raster::aggregate(rasterFromXYZ(mat), fact = 10, fun = max)
nas <- !is.na(mat10[1:length(mat10)])
values(mat10)[nas] <- fullCounts10[nas] / nrow(df_fit)
plot(mat10)

