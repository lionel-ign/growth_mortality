## model fit
fit_model_tree <- function(data, espar, indice) {
  
  # --- 1. Définir les covariables ---
  Xvar <- str_split_1(indice, pattern = "\\+|\\*")
  
  # --- 2. Filtrer espèce + NAs ---
  sub_dat <- data |>
    dplyr::filter(espar.x == espar) |>
    dplyr::select(all_of(Xvar), mort, ser_86, incref, npp) |>
    tidyr::drop_na(all_of(Xvar))
  
  # if no dead tree pass
  if (sum(sub_dat$mort, na.rm=TRUE) == 0) return(NULL)

  # --- 3. Préparer les variables ---
  sub_dat <- sub_dat |>
    dplyr::mutate(
      across(all_of(Xvar) & where(is.numeric), ~ as.numeric(scale(.x))),
      ser_id  = as.integer(as.factor(ser_86)),
      npp_id  = as.integer(as.factor(interaction(ser_86, npp))),
      increff = as.integer(as.factor(incref))
    )
  
  # --- 4. Formule ---
  ff <- as.formula(paste0(
    "mort ~ 1 + ", indice, " +
    f(ser_id,  model = 'iid',
      hyper = list(prec = list(prior = 'pc.prec', param = c(0.5, 0.01)))) +
    f(npp_id,  model = 'iid',
      hyper = list(prec = list(prior = 'pc.prec', param = c(0.5, 0.01)))) +
    f(increff, model = 'rw2',
      hyper = list(prec = list(prior = 'pc.prec', param = c(0.5, 0.01))))"
  ))
  
  # --- 5. Fit INLA ---
  ii <- INLA::inla(
    ff,
    data             = sub_dat,
    family           = "binomial",
    control.fixed    = list(mean = 0, prec = 4),
    control.family   = list(link = "logit"),
    control.predictor = list(compute = TRUE, link = 1),
    control.inla     = list(int.strategy = "eb"),
    control.compute  = list(config = TRUE)
  )
  
  # add espar and indice info in model object
  ii$.meta <- list(espar = espar, indice = indice)
  
  return(ii)
}

## extract model summary
extract_summary <- function(model, espar, indice) {
  if (is.null(model)) return(NULL)
  
  model$summary.fixed |>
    tibble::rownames_to_column("term") |>
    dplyr::mutate(espar = espar, indice = indice)
}

## extract model fitted values
extract_fit <- function(data, model, espar, indice) {
  if (is.null(model)) return(NULL)
  
  # --- 2. Variables du modèle ---
  Xvar <- str_split_1(indice, pattern = "\\+|\\*")
  
  # --- 3. Préparer les données ---
  sub_dat <- data |>
    dplyr::filter(espar.x == espar) |>
    dplyr::select(all_of(Xvar), mort, ser_86, incref, npp,
                  gtot, cdim) |>
    tidyr::drop_na(all_of(Xvar))
  
  # Séparer obs connues (mort non-NA) des prédictions futures
  is_obs   <- !is.na(sub_dat$mort)
  sub_obs  <- sub_dat[is_obs, ]
  proba    <- model$summary.fitted.values[is_obs, "mean"]
  
  # --- 4. Seuils par classe de taille ---
  cdiam <- tryCatch(
    cut(sub_obs$gtot, breaks = quantile(sub_obs$gtot, probs = seq(0, 1, by = 0.2))),
    error = function(e) factor(sub_obs$cdim)
  )
  
  t_l <- sapply(levels(cdiam), function(x) {
    tryCatch(
      pROC::coords(
        pROC::roc(sub_obs$mort[cdiam == x], proba[cdiam == x], quiet = TRUE),
        "best"
      )$threshold[1],
      error = function(e) 0.5
    )
  })
  
  # --- 5. Prédictions binaires ---
  sub_obs$pred <- as.integer(proba > unlist(t_l[as.numeric(cdiam)]))
  
  # --- 6. Métriques de classification ---
  cc_res <- sub_obs |>
    dplyr::mutate(class = dplyr::case_when(
      pred == 0 & mort == 0 ~ "TN",
      pred == 1 & mort == 0 ~ "FP",
      pred == 0 & mort == 1 ~ "FN",
      pred == 1 & mort == 1 ~ "TP"
    )) |>
    dplyr::count(class) |>
    dplyr::filter(!is.na(class)) |>
    tidyr::complete(class = c("TN", "TP", "FN", "FP"), fill = list(n = 0)) |>
    tidyr::pivot_wider(names_from = class, values_from = n, values_fill = 0) |>
    dplyr::mutate(
      tot         = FN + FP + TN + TP,
      accuracy    = (TN + TP) / tot,
      sensitivity = TP / (TP + FN),
      specificity = TN / (TN + FP),
      thresh      = mean(t_l)
    )
  
  # --- 7. Métriques supplémentaires ---
  cc_res$coef_gr   <- model$summary.fixed[2, "mean"]
  cc_res$coef_si   <- model$summary.fixed[3, "mean"]
  cc_res$coef_co   <- model$summary.fixed[4, "mean"]
  cc_res$mort_avg    <- arm::invlogit(model$summary.fixed["(Intercept)", "mean"])
  cc_res$prauc  <- MLmetrics::PRAUC(proba, sub_obs$mort)
  cc_res$p_dead <- mean(sub_obs$mort, na.rm = TRUE)
  cc_res$n_dead <- sum(sub_obs$mort, na.rm = TRUE)
  cc_res$n_tot  <- sum(!is.na(sub_obs$mort))
  cc_res$auc <- MLmetrics::AUC(proba, sub_obs$mort)
  cc_res$espar  <- espar
  cc_res$indice <- indice
  cc_res$mcc <- mltools::mcc(preds = sub_obs$pred,
                             actuals = sub_obs$mort)
  
  return(cc_res)
}

## extract preds
extract_pred <- function(data, model, espar, indice,
                         nSim = 200, chunk_size = 100){
  if (is.null(model)) return(NULL)
  inla.setOption(num.threads = "1:1")
  Xvar <- str_split_1(indice, pattern = "\\+|\\*")
  
  sub_dat <- data |>
    dplyr::filter(espar.x == espar) |>
    dplyr::select(all_of(Xvar), mort, ser_86, incref, npp, gtot, cdim, w, poids) |>
    tidyr::drop_na(all_of(Xvar))
  
  is_obs <- !is.na(sub_dat$mort)
  proba  <- model$summary.fitted.values[, "mean"]
  nObs   <- nrow(sub_dat)
  
  # --- Seuils par classe de taille ---
  cdiam <- tryCatch(
    base::cut(sub_dat$gtot,
              breaks = quantile(sub_dat$gtot, probs = seq(0, 1, by = 0.2)),
              include.lowest = TRUE),
    error = function(e) factor(sub_dat$cdim)
  )
  
  t_l <- sapply(levels(cdiam), function(x) {
    tryCatch(
      pROC::coords(
        pROC::roc(sub_dat$mort[cdiam == x & is_obs],
                  proba[cdiam == x & is_obs], quiet = TRUE),
        "best"
      )$threshold[1],
      error = function(e) 0.5
    )
  })
  
  thresh_vec <- t_l[as.numeric(cdiam)]  # vecteur de seuils aligné sur sub_dat
  
  # --- Traitement par chunks ---
  chunks    <- split(seq_len(nSim), ceiling(seq_len(nSim) / chunk_size))
  results   <- vector("list", length(chunks))
  
  for (i in seq_along(chunks)) {
    n_chunk <- length(chunks[[i]])
    
    samples <- inla.posterior.sample(
      n_chunk, model,
      selection         = list(Predictor = 1:nObs),
      use.improved.mean = FALSE,
      skew.corr         = FALSE
    )
    
    eta <- sapply(samples, function(s)
      s$latent[grep("^Predictor", rownames(s$latent)), ]
    )
    rm(samples); gc()
    
    p   <- arm::invlogit(eta)
    sim <- matrix(as.integer(p >= thresh_vec), nrow = nObs, ncol = n_chunk)
    rm(p, eta); gc()
    
    # Agréger directement — ne jamais stocker sim en entier longtemps
    ba_mat <- sim * sub_dat$gtot   # ba des arbres morts simulés
    
    results[[i]] <- as_tibble(ba_mat) |>
      mutate(npp  = sub_dat$npp,
             year = sub_dat$incref + 2005,
             w  = sub_dat$w,
             poid = sub_dat$poids) |>
      pivot_longer(starts_with("V")) |>
      group_by(npp, year, name, poid) |>
      summarise(ba_dead = sum(value * w) / sum(w), .groups = "drop") |>
      group_by(year, name) |>
      summarise(ba_dead = sum(ba_dead * poid) / sum(poid), .groups = "drop")
    
    rm(ba_mat, sim); gc()
    #message(glue::glue("  chunk {i}/{length(chunks)} OK"))
  }
  
  # --- Combiner les chunks et résumer ---
  bad <- dplyr::bind_rows(results) |>
    group_by(year) |>
    summarise(
      ba_mea = mean(ba_dead,                      na.rm = TRUE),
      ba_low = quantile(ba_dead, probs = 0.05,    na.rm = TRUE),
      ba_hig = quantile(ba_dead, probs = 0.95,    na.rm = TRUE)
    ) |>
    mutate(type = "pred", espar = espar, indice = indice)
  
  return(bad)
}

## compute bias between obs and pred
get_pred_bias <- function(all_pred, mort_trend){
  ## data wraggling for the ocre output
  mort_trend %>%
    filter(visite == "1") %>%
    select(espar.espar1, year, moy_ha) %>%
    rename(espar = espar.espar1,
           ba_pred = moy_ha) -> p_mod
  
  ## join with preds and compute mean bias per species
  inner_join(all_pred, p_mod, by = c("espar", "year")) %>%
    mutate(bias = ba_pred - ba_mea) %>%
    group_by(espar) %>%
    summarise(bias_med = median(bias)) -> espar_bias
  
  return(espar_bias)
}