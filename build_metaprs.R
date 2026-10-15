#!/usr/bin/env Rscript
# metaPRS construction from eligible, calculated PGS values.
# Input: IID, y (HF=1/control=0), AGE, SEX (0/1), PC1-PC10 and PGS columns.
# Run: Rscript build_metaprs.R data/analysis_data.csv results

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
})

safe_numeric <- function(x) {
    if (is.numeric(x) || is.integer(x)) 
        return(as.numeric(x))
    suppressWarnings(as.numeric(trimws(as.character(x))))
}

stratified_split <- function(y, fraction, seed) {
    set.seed(seed)
    out <- rep("test", length(y))
    for (cls in sort(unique(y))) {
        idx <- which(y == cls)
        n_train <- floor(length(idx) * fraction)
        n_train <- max(1L, min(length(idx) - 1L, n_train))
        out[sample(idx, n_train, replace = FALSE)] <- "training"
    }
    out
}

make_foldid <- function(y, seed) {
  set.seed(seed)
  k <- min(10L, min(table(y)))
  foldid <- integer(length(y))
  for (cls in sort(unique(y))) {
    idx <- which(y == cls)
    foldid[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  foldid
}

prepare_covariates <- function(data, training, variables) {
  X <- as.matrix(data[, ..variables])
  for (j in seq_len(ncol(X))) {
    x <- safe_numeric(X[, j])
    mu <- mean(x[training])
    sig <- sd(x[training])
    if (!is.finite(sig) || sig <= 0) { mu <- 0; sig <- 1 }
    X[, j] <- (x - mu) / sig
  }
  storage.mode(X) <- "double"
  X
}

screen_one_pgs <- function(y, pgs_z, cov_matrix) {
    ok <- is.finite(y) & is.finite(pgs_z) & apply(is.finite(cov_matrix), 1, all)
    empty <- c(beta = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_, n = sum(ok), converged = 0)
    if (sum(ok) < 20 || length(unique(y[ok])) != 2 || length(unique(pgs_z[ok])) < 2) 
        return(empty)
    X <- cbind(`(Intercept)` = 1, cov_matrix[ok, , drop = FALSE], PGS_z = pgs_z[ok])
    fit <- tryCatch(suppressWarnings(glm.fit(x = X, y = y[ok], family = binomial(), control = glm.control(maxit = 100, 
        epsilon = 1e-08))), error = function(e) NULL)
    if (is.null(fit)) 
        return(empty)
    pidx <- ncol(X)
    beta <- fit$coefficients[pidx]
    empty["beta"] <- beta
    empty["converged"] <- as.numeric(isTRUE(fit$converged))
    if (!is.finite(beta) || fit$rank == 0) 
        return(empty)
    pivot <- fit$qr$pivot[seq_len(fit$rank)]
    if (!(pidx %in% pivot)) 
        return(empty)
    R <- qr.R(fit$qr)[seq_len(fit$rank), seq_len(fit$rank), drop = FALSE]
    cov_rank <- tryCatch(chol2inv(R), error = function(e) NULL)
    if (is.null(cov_rank)) 
        return(empty)
    pos <- match(pidx, pivot)
    se <- sqrt(cov_rank[pos, pos])
    empty["se"] <- se
    if (!is.finite(se) || se <= 0) 
        return(empty)
    zv <- beta/se
    empty["z"] <- zv
    empty["p"] <- 2 * pnorm(-abs(zv))
    empty
}

greedy_prune <- function(ranked, correlations) {
  kept <- character()
  for (candidate in ranked) {
    r <- abs(correlations[candidate, kept])
    if (!any(is.finite(r) & r > 0.50)) kept <- c(kept, candidate)
  }
  kept
}

build_metaprs <- function(data,
                          pgs_columns = sort(grep("^PGS[0-9]+$", names(data), value = TRUE))) {
  split_seed <- fold_seed <- 20260624
  dt <- as.data.table(copy(data))
  dt[, y := safe_numeric(y)]
  covariates <- c("AGE", "SEX", paste0("PC", 1:10))
  dt[, dataset := stratified_split(y, 0.80, split_seed)]
  training <- which(dt$dataset == "training")
  controls <- which(dt$dataset == "training" & dt$y == 0)
  y_training <- dt$y[training]
  foldid <- make_foldid(y_training, fold_seed)
  C <- prepare_covariates(dt, training, covariates)
  normalization <- vector("list", length(pgs_columns))
  P <- matrix(NA_real_, nrow(dt), length(pgs_columns),
              dimnames = list(NULL, pgs_columns))
  for (j in seq_along(pgs_columns)) {
    variable <- pgs_columns[j]
    value <- safe_numeric(dt[[variable]])
    median_control <- median(value[controls], na.rm = TRUE)
    if (is.finite(median_control)) {
      value[!is.finite(value)] <- median_control
      mean_control <- mean(value[controls])
      sd_control <- sd(value[controls])
    } else {
      mean_control <- sd_control <- NA_real_
    }
    valid <- is.finite(mean_control) && is.finite(sd_control) && sd_control > 0
    normalization[[j]] <- data.table(
      variable = variable, training_control_median = median_control,
      training_control_mean = mean_control, training_control_sd = sd_control,
      valid = valid
    )
    if (valid) P[, j] <- (value - mean_control) / sd_control
  }
  normalization <- rbindlist(normalization)
  valid_pgs <- normalization[valid == TRUE, variable]
  P <- P[, valid_pgs, drop = FALSE]
  screening <- rbindlist(lapply(valid_pgs, function(variable) {
    fit <- screen_one_pgs(y_training, P[training, variable],
                          C[training, , drop = FALSE])
    data.table(variable = variable, beta = fit["beta"],
               standard_error = fit["se"], z_value = fit["z"],
               p_value = fit["p"], n_analysis = as.integer(fit["n"]),
               converged = as.logical(fit["converged"]))
  }))
  screening[, FDR_BH := p.adjust(p_value, method = "BH")]
  screening[, selected_FDR005 := is.finite(FDR_BH) & FDR_BH < 0.05]
  screening[, abs_z_internal := abs(z_value)]
  setorder(screening, FDR_BH, p_value, -abs_z_internal, variable)
  selected <- screening[selected_FDR005 == TRUE, variable]
  correlations <- cor(P[controls, selected, drop = FALSE],
                      use = "pairwise.complete.obs", method = "pearson")
  candidates <- greedy_prune(selected, correlations)
  X <- cbind(P[training, candidates, drop = FALSE], C[training, , drop = FALSE])
  penalty <- c(rep(1, length(candidates)), rep(0, ncol(C)))
  alphas <- seq(0, 1, by = 0.1)
  fits <- vector("list", length(alphas))
  cv_summary <- vector("list", length(alphas))
  for (i in seq_along(alphas)) {
    set.seed(fold_seed)
    fit <- cv.glmnet(
      x = X, y = y_training, family = "binomial", alpha = alphas[i],
      foldid = foldid, type.measure = "auc", standardize = FALSE,
      intercept = TRUE, nlambda = 100L, penalty.factor = penalty
    )
    index <- which.min(abs(log(fit$lambda) - log(fit$lambda.1se)))
    fits[[i]] <- fit
    cv_summary[[i]] <- data.table(
      alpha = alphas[i], lambda_1SE = fit$lambda.1se,
      mean_CV_AUC_at_lambda_1SE = fit$cvm[index]
    )
  }
  cv_summary <- rbindlist(cv_summary)
  setorder(cv_summary, -mean_CV_AUC_at_lambda_1SE, alpha)
  alpha <- cv_summary$alpha[1]
  fit <- fits[[match(alpha, alphas)]]
  lambda <- fit$lambda.1se
  coefficient_matrix <- as.matrix(coef(fit, s = lambda))
  coefficients <- data.table(
    variable = rownames(coefficient_matrix),
    coefficient = as.numeric(coefficient_matrix[, 1])
  )
  coefficients[, term_type := fifelse(variable %in% candidates, "PGS",
    fifelse(variable == "(Intercept)", "Intercept", "Adjustment covariate"))]
  selected_coefficients <- coefficients[term_type == "PGS" & coefficient != 0]
  raw_score <- if (nrow(selected_coefficients)) {
    as.vector(P[, selected_coefficients$variable, drop = FALSE] %*%
                selected_coefficients$coefficient)
  } else rep(0, nrow(dt))
  score_mean <- mean(raw_score[controls])
  score_sd <- sd(raw_score[controls])
  scores <- dt[, .(IID, y, dataset)]
  scores[, metaPRS_score := raw_score]
  scores[, metaPRS_score_z := (raw_score - score_mean) / score_sd]

  list(
    candidates = candidates, cv_summary = cv_summary,
    selected_alpha = alpha, selected_lambda = lambda,
    coefficients = coefficients, scores = scores,
    pgs_normalization = normalization,
    final_score_control_mean = score_mean, final_score_control_sd = score_sd,
    fitted_model = fit
  )
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  input_file <- if (length(args) >= 1L) args[1] else "data/analysis_data.csv"
  output_dir <- if (length(args) >= 2L) args[2] else "results"
  data <- fread(input_file, colClasses = c(IID = "character"))
  result <- build_metaprs(data)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, file.path(output_dir, "metaprs_model.rds"))
  fwrite(result$coefficients, file.path(output_dir, "metaprs_coefficients.csv"))
  fwrite(result$scores, file.path(output_dir, "metaprs_scores.csv"))
}
