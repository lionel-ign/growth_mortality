## boxplot dead/alive per growth metrics
plot_growth_status <- function(data, output_dir = "figures"){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  data %>%
    filter(type != "pred") %>%
    select(type, abg, abg2, abg2c, abg_diff, abg_pl, ir_trend,
           ir_sd, relbai, relbac) %>%
    pivot_longer(abg:relbac) %>%
    mutate(name = factor(name, labels = c("BAI1", "dBAI",
                                          "BAI1 plot", "BAI2",
                                          "BAI2 corr", "RI var",
                                          "RI trend",
                                          "relBAIc", "relBAIi"))) %>%
    ggplot(aes(x=type, y=value)) +
    geom_boxplot(outliers = FALSE) +
    facet_wrap(vars(name), scales = "free") +
    labs(x = "Tree status",
         y = "Growth metric") -> gg
  
  path <- file.path(output_dir, "boxplot_metrics.png")
  ggsave(path, gg, width = 8, height = 6, dpi = 150)
  return(path)  # tar_target(format = "file") attend un chemin
  
}

## plot growth metrics correlations
plot_growth_corr <- function(data, output_dir = "figures"){
  data %>%
    filter(type != "pred") %>%
    select(abg, abg2, abg2c, abg_diff, abg_pl, ir_trend,
           ir_sd, relbai, relbac) -> tt
  
  corr <- round(cor(tt, 
                    use = "pairwise.complete.obs"), 2)
  
  gg <- ggcorrplot(corr, legend.title = "Correlation\ncoefficient",
                   hc.order = TRUE, type = "lower", lab = TRUE,
                   insig = "blank",
                   outline.color = "white",
                   ggtheme = ggplot2::theme_gray,
                   colors = c("#6D9EC1", "white", "#E46726")
  )
  
  path <- file.path(output_dir, "growth_correlation.png")
  ggsave(path, gg, width = 8, height = 6, dpi = 150)
  return(path)  # tar_target(format = "file") attend un chemin
}

## plot trait correlation
plot_trait_corr <- function(traits, all_fit, output_dir = "figures"){
  
  all_fit %>%
    group_by(espar) %>%
    select(espar, n_dead, n_tot) %>%
    slice_head(n=1) %>%
    mutate(p_dead = n_dead / n_tot) %>%
    inner_join(traits) -> tt
  
  corr <- round(cor(tt[,c(3:4,7:16)], 
                    use = "pairwise.complete.obs"), 2)
  
  pval <- cor_pmat(tt[,c(3:4,7:16)])
  
  gg <- ggcorrplot(corr, legend.title = "Correlation\ncoefficient",
             hc.order = FALSE, type = "lower",
             p.mat = pval, lab = TRUE,
             insig = "blank",
             outline.color = "white",
             ggtheme = ggplot2::theme_gray,
             colors = c("#6D9EC1", "white", "#E46726")
  )
  
  path <- file.path(output_dir, "trait_corrmat.png")
  ggsave(path, gg, width = 8, height = 6, dpi = 150)
  return(path)  # tar_target(format = "file") attend un chemin
}


## Plot des métriques de fit (prauc, mort avg, coef, bias) par espèce 
plot_fitespar_metrics <- function(all_fit, all_summaries, sp_bias, lbls, output_dir = "figures") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # espar - species link
  tt <- read.csv("../IFN_stuff/data/pommier/espar_species.csv")
  
  # all_fit %>%
  #   inner_join(lbls, by = "indice") %>%
  #   filter(indiced %in% c("BAI1", "BAI2", "BAI2c", "dBAI", "relBAI",
  #                         "BAIc", "BAI1plot", "BAI1+\nDOM")) %>%
  #   select(espar, indice, prauc, mcc) %>%
  #   tidyr::pivot_longer(cols = c(mcc, prauc),
  #                       names_to = "metric", values_to = "value") %>%
  #   group_by(espar, metric) %>%
  #   summarise(mi = min(value),
  #             ma = max(value)) %>%
  #   inner_join(tt, by = "espar") -> dd
  
  all_summaries %>%
    filter(term == "abg_pl",
           indice == "abg_pl+gtot_0+ba_tot") %>%
    select(espar, mean, `0.025quant`, `0.975quant`) %>%
    rename(coef = mean) -> gcoef
  
  all_summaries %>%
    filter(term == "(Intercept)",
           indice == "abg_pl+gtot_0+ba_tot") %>%
    mutate(lci = arm::invlogit(mean - 1.96 * sd),
           uci = arm::invlogit(mean + 1.96 * sd),
           intercept = arm::invlogit(mean)) %>%
    select(espar, intercept, lci, uci) -> icoef
  
  all_fit %>%
    inner_join(lbls, by = "indice") %>%
    filter(indiced %in% c("BAI1plot")) %>%
    select(espar, indice, prauc, mcc) %>%
    inner_join(tt, by = "espar") %>%
    inner_join(sp_bias, by = "espar") %>%
    inner_join(gcoef, by = "espar") %>%
    inner_join(icoef, by = "espar") -> dd
    
  
  g1 <- ggplot(dd,
                aes(x = fct_reorder(species, prauc),
                    y = prauc)) +
    geom_hline(yintercept = mean(dd$prauc),
               linetype = "dashed") +
        geom_point() +
        coord_flip() +
        labs(x = "", y = "Performance score",
             subtitle = "prAUC") +
        theme_bw() +
        theme(text = element_text(size = 20))
  
  g2 <- ggplot(dd,
               aes(x = fct_reorder(species, mcc),
                   y = mcc)) +
    geom_hline(yintercept = mean(dd$mcc),
               linetype = "dashed") +
    geom_point() +
    coord_flip() +
    labs(x = "", y = "Performance score",
         subtitle = "MCC") +
    theme_bw() +
    theme(text = element_text(size = 20))
  
  g3 <- ggplot(dd,
               aes(x = fct_reorder(species, bias_med),
                   y = bias_med * 1e4)) +
    geom_hline(yintercept = mean(dd$bias_med * 1e4),
               linetype = "dashed") +
    geom_point() +
    coord_flip() +
    labs(x = "", y = "Median dead BA bias (cm²/ha)",
         subtitle = "Bias") +
    theme_bw() +
    theme(text = element_text(size = 20))
  
  g4 <- ggplot(dd,
               aes(x = fct_reorder(species, coef),
                   y = coef,
                   ymin = `0.025quant`,
                   ymax = `0.975quant`)) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = mean(dd$coef),
               linetype = "dashed") +
    geom_linerange() +
    geom_point() +
    coord_flip() +
    labs(x = "", y = "Growth regression\ncoefficient",
         subtitle = "Coefficient") +
    theme_bw() +
    theme(text = element_text(size = 20))
  
  # g4 <- ggplot(dd,
  #              aes(x = fct_reorder(species, intercept),
  #                  y = intercept,
  #                  ymin = lci,
  #                  ymax = uci)) +
  #   geom_hline(yintercept = mean(dd$intercept),
  #              linetype = "dashed") +
  #   geom_linerange() +
  #   geom_point() +
  #   coord_flip() +
  #   labs(x = "", y = "Predicted average\nmortality probability",
  #        subtitle = "Mortality probability") +
  #   theme_bw() +
  #   theme(text = element_text(size = 20))
  
  
  p <- patchwork::wrap_plots(g1, g2, axes = "collect")
  p2 <- patchwork::wrap_plots(g3, g4)
  pf <- patchwork::wrap_plots(p, p2, nrow = 2)

  path <- file.path(output_dir, "fit_espar_metrics.png")
  ggsave(path, pf, width = 12, height = 12, dpi = 150)
  return(path)  # tar_target(format = "file") attend un chemin
}

## plot des métriques de fit, comparison boxplot per metrics
plot_fit_metrics <- function(all_fit, lbls, output_dir = "figures"){
  
  all_fit %>%
    select(indice, prauc, mcc) %>%
    inner_join(lbls, by = "indice") -> dat
  
  ## Define your two groups of indices
  group_a <- c("BAI1", "BAI2", "BAI2c", "dBAI", "BAI1plot",
               "relBAI", "BAIc", "RItrend",
               "RIsd", "BAI1+\nDOM")
  group_b <- c("relBAI", "relBAC", "RItrend",
               "RIsd")
  
  make_plot <- function(data) {
    g1 <- ggplot(data, aes(x = fct_reorder(indiced, prauc, .fun = mean), y = prauc)) +
      stat_summary(fun.data = "mean_cl_boot") +
      theme_bw() +
      theme(text = element_text(size = 20)) +
      coord_flip() +
      labs(x = "Growth metric", y = "Performance", subtitle = "prAUC")
    
    g2 <- ggplot(data, aes(x = fct_reorder(indiced, mcc, .fun = mean), y = mcc)) +
      stat_summary(fun.data = "mean_cl_boot") +
      theme_bw() +
      theme(text = element_text(size = 20)) +
      coord_flip() +
      labs(x = "Growth metric", y = "Performance", subtitle = "MCC")
    
    
    
    patchwork::wrap_plots(g1, g2, ncol = 2, axes = "collect")
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  path_a <- file.path(output_dir, "fit_metrics_group_a.png")
  path_b <- file.path(output_dir, "fit_metrics_group_b.png")
  
  ggsave(path_a, make_plot(filter(dat, indiced %in% group_a)), width = 12, height = 8, dpi = 150)
  ggsave(path_b, make_plot(filter(dat, indiced %in% group_b)), width = 12, height = 8, dpi = 150)
  
  c(path_a, path_b)  # targets tracks both files
}

## plot des effets de la croissance
plot_abgint_effects <- function(all_summaries, output_dir = "figures") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --- 1. Extraire les indices concernés ---
  indices_abg  <- c("abg_pl+gtot_0+ba_tot")
  
  # --- 2. Effets simples (abg et abg2, sans interaction) ---
  simple_effects <- all_summaries |>
    dplyr::filter(
      indice %in% indices_abg,
      term %in% c("abg", "abg2", "abg_pl")
    ) |>
    dplyr::mutate(label = term)  # un seul point par espèce × indice
  
 
  
  # --- 4. Formater les effets simples pour bind ---
  simple_fmt <- simple_effects |>
    dplyr::mutate(
      q025  = `0.025quant`,
      q975  = `0.975quant`
    ) |>
    dplyr::select(espar, indice, mean, q025, q975, label)
  
  # --- 5. Combiner ---
  lbls <- read.csv("~/LIF/IFN_stuff/data/pommier/espar_species.csv")
  plot_dat <- dplyr::bind_rows(simple_fmt) |>
    inner_join(lbls, by = "espar")
  
  # --- 6. Figure ---
  p <- plot_dat |>
    ggplot(aes(x = mean, y = reorder(species, mean))) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(xmin = q025, xmax = q975),
                   width = 0) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    theme_bw() +
    labs(x = "Growth metric coefficient",
         y = "Species")
  
  path <- file.path(output_dir, "abg_effects.png")
  ggsave(path, p, width = 11, height = 9, dpi = 150)
  return(path)
}

## Plot des effets fixes issus de all_summary
plot_fixed_effects <- function(all_summaries, output_dir = "figures") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lbls <- read.csv("~/LIF/IFN_stuff/data/pommier/espar_species.csv")
  
  p <- all_summaries |>
    inner_join(lbls, by = "espar") |>
    dplyr::filter(!(term %in% c("(Intercept)", "ba_tot", "gtot", "gtot_0"))) |>
    ggplot(aes(x = mean, y = reorder(species, mean))) +
    geom_point(size = 2) +
    geom_errorbar(aes(xmin = `0.025quant`, xmax = `0.975quant`),
                   height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~term, scales = "free_x") +
    labs(x = "Effet (logit)", y = "",
         title = "") +
    theme_bw()
  
  path <- file.path(output_dir, "fixed_effects.png")
  ggsave(path, p, width = 14, height = 8, dpi = 150)
  return(path)
}

## plot des courbes mortalité/croissance prédite dans le modèle avec dom
plot_resp_dom <- function(data, all_summaries, traits, output_dir = "figures"){
  # ── 1. Extract coefficients, keyed by term name ───────────────────────────
  coefs <- all_summaries %>%
    filter(indice == "abg+gtot_0+ba_tot+lib") %>%
    select(term, mean, espar) %>%
    group_by(espar) %>%
    summarise(
      coef = list(setNames(mean, term)),  # named numeric vector → safe lookup
      .groups = "drop"
    )
  
  # ── 2. Per-species gradient ranges + scaling parameters ───────────────────
  gradients <- data %>%
    group_by(espar = espar.x) %>%
    mutate(  
      growth = as.numeric(scale(abg))
    ) %>%
    summarise(
      growth_min = quantile(growth, probs = 0.05),
      growth_max = quantile(growth, probs = 0.95),
      growth_mean = mean(abg),
      growth_sd   = sd(abg),
      .groups = "drop"
    )
  lbls <- read.csv("~/LIF/IFN_stuff/data/pommier/espar_species.csv")
  # ── 3. Build prediction grid, compute and back-transform ──────────────────
  pred_data <- gradients %>%
    inner_join(coefs, by = "espar") %>%
    rowwise() %>%
    mutate(
      grid = list({
        
        # Back-transform quantiles to original scale to define log-spaced seq
        growth_orig_min <- growth_min * growth_sd + growth_mean
        growth_orig_max <- growth_max * growth_sd + growth_mean
        
        # Generate grid in original scale, log-spaced for even visual spacing
        g <- expand.grid(
          growth_orig = seq(growth_orig_min, growth_orig_max, length.out = 50),
          lib = c("dom", "sup")
        )
        
        # Re-scale to model scale for prediction
        g$growth <- (g$growth_orig - growth_mean) / growth_sd
        
        b0     <- coef[["(Intercept)"]]
        b_abg  <- coef[["abg"]]
        b_size <- coef[["gtot_0"]]
        b_lib  <- coef[["libsup"]]
        
        g$pred <- arm::invlogit(b0 + b_abg * g$growth + ifelse(g$lib == "sup", b_lib, 0))
        
        g
      })
    ) %>%
    ungroup() %>%
    select(espar, grid) %>%
    tidyr::unnest(grid) %>%
    inner_join(lbls, by = "espar")
  
  # ── 4. Plot on original log scale ─────────────────────────────────────────
  p1 <- ggplot(pred_data, aes(x = growth_orig, y = pred, color = lib)) +
    geom_path() +
    facet_wrap(~species, scales = "free") + 
    labs(
      x     = "Growth",
      y     = "Probability of mortality"
    )
  
  
  # ── 5. Compute area between curves per species ────────────────────────────
  traits$Shade_tolerance[traits$espar == "51"] <- 2.1
  area_between <- pred_data %>%
    select(espar, species, growth_orig, lib, pred) %>%
    pivot_wider(names_from = lib, values_from = pred) %>%
    group_by(espar, species) %>%
    arrange(growth_orig, .by_group = TRUE) %>%
    summarise(
      area = sum(diff(growth_orig) * (abs(dom - sup)[-1] + abs(dom - sup)[-n()]) / 2),
      area_scaled = area / (max(growth_orig) - min(growth_orig)),
      .groups = "drop"
    ) %>%
    inner_join(select(traits, espar, Shade_tolerance), by = "espar")
  
  ## correlation coefficient
  rr <- cor.test(area_between$Shade_tolerance, area_between$area_scaled,
            method = "spearman")

  p2 <- ggplot(area_between, aes(x=Shade_tolerance, y=area_scaled, labels = species)) +
    ggrepel::geom_label_repel() +
    scale_y_log10() +
    annotate(geom = "text",
             x = min(area_between$Shade_tolerance),
             y = max(area_between$area_scaled),
             label = paste0("R = ", round(rr$estimate, 2)),
             hjust = 0,
             vjust = 1) +
    labs(x = "Shade tolerance",
         y = "Area between regression curves")
  
  # ── 5. Save ───────────────────────────────────────────────────────────────
  path_a <- file.path(output_dir, "response_curve_dom.png")
  path_b <- file.path(output_dir, "shade_area.png")
  
  ggsave(path_a, p1, width = 12, height = 8, dpi = 150)
  ggsave(path_b, p2, width = 12, height = 8, dpi = 150)
  
  c(path_a, path_b)  # targets tracks both files
  
  
}

## Plot comparatif prauc vs p_dead (performance relative à la prévalence)
plot_prauc_vs_prevalence <- function(all_fit, output_dir = "figures") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  p <- all_fit |>
    dplyr::filter(indice == "abg+gtot_0+ba_tot") |>
    ggplot(aes(x = p_dead, y = prauc, label = espar)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 15) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey50") +
    labs(x = "Mortality prevalence", y = "prAUC") +
    theme_bw()
  
  path <- file.path(output_dir, "prauc_vs_prevalence.png")
  ggsave(path, p, width = 10, height = 7, dpi = 150)
  return(path)
}

## plot des corrélations avec les traits
plot_corr_traits <- function(all_fit, sp_bias, traits, output_dir = "figures"){

  all_fit %>%
    filter(indice == "abg_pl+gtot_0+ba_tot") %>%
    inner_join(sp_bias, by = "espar") %>%
    select(espar, indice, mcc, mort_avg, contains("coef"), prauc,
           bias_med, n_dead, n_tot) %>%
    mutate(p_dead = n_dead / n_tot) %>%
    pivot_longer(mcc:bias_med) %>%
    inner_join(traits, by = c("espar" = "espar")) %>%
    select(espar, indice, name, value, contains("tolerance"),
           Wood_density, growth, Slenderness, height,
           age, p_dead, n_tot) %>%
    pivot_longer(5:16,
                 names_to = "trait", values_to = "tr_value") %>%
    nest(data = -c(name, trait)) %>%
    mutate(cc = map(data, ~cor.test(.x$value, .x$tr_value,
                                    use = "complete.obs",
                                    method = "spearman")),
           cct = map(cc, tidy)) %>%
    unnest(cct) %>%
    select(-data, -cc) %>%
    arrange(desc(estimate), p.value) %>%
    mutate(sig = ifelse(p.value < 0.1, "sig.", "not sig.")) %>%
    filter(name != "mcc") -> rr
  
  rr2 <- rr
  rr2$trait <- factor(rr2$trait, labels = c("Mean age",
                                          "Drought intolerance",
                                          "Drought tolerance",
                                          "Frost intolerance",
                                          "Mean growth",
                                          "Mean height",
                                          "Nb. tree",
                                          "Per. dead",
                                          "Shade tolerance",
                                          "Slenderness",
                                          "Waterlogging tolerance",
                                          "Mean wood density"))
  
  rr2$name <- factor(rr2$name,
                    levels = c("mort_avg", "coef_gr", "coef_si","coef_co", "prauc", "bias_med"),
                    labels = c( "Avg. mort.", "Coef. Gr.", "Coef. Size",
                                "Coef. Comp.","prAUC", "Bias"))
  
  grp_a <- c("Drought intolerance",
             "Drought tolerance",
             "Frost intolerance",
             "Shade tolerance",
             "Waterlogging tolerance",
             "Slenderness",
             "Mean age",
             "Mean wood density")
  grp_b <- c("Mean age",
             "Mean growth",
             "Mean height",
             "Mean wood density",
             "Slenderness")
  grp_c <- c("Per. dead", "Nb. tree")
  
  p1_data <- droplevels(subset(rr2, trait %in% grp_a))
  p1 <- ggplot(p1_data, aes(x = name, y = trait, fill = estimate)) +
    geom_raster() +
    geom_tile(
      data  = dplyr::filter(p1_data, sig == "sig."),  # uniquement les cases sig.
      color = "black",   # couleur du cadre
      fill  = NA,        # ne pas écraser le fill du raster
      linewidth = 1.2    # épaisseur du cadre
    ) +
    scale_fill_gradient2() +
    theme_bw() +
    theme(text = element_text(size = 16)) +
    labs(x = "Model metrics",
         y = "Species traits",
         fill = "Correlation\nestimates")
  
  p2_data <- droplevels(subset(rr2, trait %in% grp_b))
  p2 <- ggplot(p2_data, aes(x = name, y = trait, fill = estimate)) +
    geom_raster() +
    geom_tile(
      data  = dplyr::filter(p2_data, sig == "sig."),  # uniquement les cases sig.
      color = "black",   # couleur du cadre
      fill  = NA,        # ne pas écraser le fill du raster
      linewidth = 1.2    # épaisseur du cadre
    ) +
    scale_fill_gradient2() +
    theme_bw() +
    theme(text = element_text(size = 16)) +
    labs(x = "Model metrics",
         y = "Species characteristics",
         fill = "Correlation\nestimates")
  
  pp <- patchwork::wrap_plots(p1, p2, axes = "collect")
  
  p3_data <- droplevels(subset(rr2, trait %in% grp_c))
  p3 <- ggplot(p3_data, aes(x = name, y = trait, fill = estimate)) +
    geom_raster() +
    geom_tile(
      data  = dplyr::filter(p3_data, sig == "sig."),  # uniquement les cases sig.
      color = "black",   # couleur du cadre
      fill  = NA,        # ne pas écraser le fill du raster
      linewidth = 1.2    # épaisseur du cadre
    ) +
    scale_fill_gradient2() +
    theme_bw() +
    labs(x = "Model metrics",
         y = "Species characteristics",
         fill = "Correlation\nestimates")
  
  ## also plot the sig. correlations
  rr %>%
    filter(sig == "sig.") %>%
    mutate(pl = map2(name, trait, plot_corr_sig,
                     all_fit = all_fit,
                     traits = traits,
                     sp_bias = sp_bias)) -> tt
  
  
  
  path_a <- file.path(output_dir, "correlation_traits.png")
  path_b <- file.path(output_dir, "correlation_size.png")
  path_c <- file.path(output_dir, "correlation_biplot.png")
  
  ggsave(path_a, p1, width = 12, height = 11, dpi = 150)
  ggsave(path_b, p3, width = 12, height = 8, dpi = 150)
  ggsave(path_c, patchwork::wrap_plots(tt$pl,
                                       ncol = ceiling(nrow(tt) / 4)))
  
  c(path_a, path_b, path_c)  # targets tracks both files
}

## plot selected correlations
plot_corr_sig <- function(metric, trait,
                          all_fit, traits, sp_bias){
  
  all_fit %>%
    inner_join(sp_bias, by = "espar") %>%
    filter(indice == "abg_pl+gtot_0+ba_tot") %>%
    select(espar, all_of(metric), p_dead, n_tot) -> d1
  
  if (trait %in% c("n_tot", "p_dead")) {
    g <- ggplot(data = d1, aes(x = .data[[trait]], y = .data[[metric]])) +
      geom_point()
  } else {
    traits %>%
      select(espar, species, all_of(trait)) %>%
      inner_join(d1, by = "espar") -> d2
    
    g <- ggplot(data = d2, aes(x = .data[[trait]], y = .data[[metric]])) +
      geom_point()
  }
  
  return(g)
}


## compute anova 
plot_aov_traits <- function(all_fit, sp_bias, traits, output_dir = "figures"){
  
  ## run anovas for different model metrics
  all_fit %>%
    filter(indice == "abg+gtot_0+ba_tot") %>%
    inner_join(sp_bias, by = "espar") %>%
    select(espar, indice, accuracy, sensitivity, coef, prauc, bias_med) %>%
    pivot_longer(accuracy:bias_med) %>%
    inner_join(tr_tab, by = c("espar" = "espar")) %>%
    select(name, value, species, Wood.structure.types..WST.) %>%
    nest(data = -name) %>%
    mutate(ao = map(data, ~ aov(value ~ Wood.structure.types..WST., .x))) %>%
    mutate(aot = map(ao, tidy)) %>%
    unnest(aot) %>%
    filter(term != "Residuals") -> res_aov

  ## identify model metrics to check
  tgt <- res_aov$name[res_aov$p.value < 0.05]  
  
  ## plot for these tgt
  all_fit %>%
    filter(indice == "abg+gtot_0+ba_tot") %>%
    inner_join(sp_bias, by = "espar") %>%
    select(espar, !!(tgt)) %>%
    pivot_longer(-espar) %>%
    inner_join(tr_tab, by = c("espar" = "espar")) %>%
    mutate(namef = factor(name, labels = c("Bias", "Coefficient"))) %>%
    ggplot(aes(x=Wood.structure.types..WST.,
               y = value)) +
    stat_summary(fun.data = "mean_cl_boot", color = "red") +
    geom_jitter(width = 0.1) +
    facet_wrap(vars(namef), scales = "free") + 
    theme_classic() +
    labs(x = "Wood structure type",
         y = "Model metric") -> pp
  
  path <- file.path(output_dir, "traits_aov.png")
  ggsave(path, pp, width = 10, height = 7, dpi = 150)
  return(path)
}

## plot pred vs obs mort trends
plot_mort_trends <- function(mort_trend, all_pred, output_dir = "figures"){
  ## data wraggling for the ocre output
  mort_trend %>%
    filter(visite == "1") %>%
    select(year, espar.espar1, tot, moy_ha, ic_moy_ha) %>%
    mutate(ba_low = moy_ha-ic_moy_ha,
           ba_hig = moy_ha+ic_moy_ha,
           type = "past",
           indice = "obs") %>%
    select(espar.espar1, year, moy_ha,
           ba_low, ba_hig, type, indice) %>%
    rename(espar = espar.espar1,
           ba_mea = moy_ha) -> p_mod
  
  ## data wraggling for the model outputs
  lbls <- read.csv("~/LIF/IFN_stuff/data/pommier/espar_species.csv")
  all_pred %>%
    mutate(type = "pred") %>%
    bind_rows(p_mod) %>%
    inner_join(lbls,
               by = "espar") -> p
  
  ## dev version en base 100
  p_base100 <- p |>
    group_by(espar, indice, type) |>
    arrange(year) |>
    mutate(
      ref     = first(ba_mea),          # valeur de référence = première année
      ba_mea  = ba_mea  / ref,
      ba_low  = ba_low  / ref,    # IC bas divisé par la même référence
      ba_hig  = ba_hig  / ref
    ) |>
    ungroup() 
  
  # 1. Isolate observed (past) data — one row per year/espar
  past_dat <- p_base100 %>%
    filter(type == "past") %>%
    select(year, espar, type, past_ba_mea = ba_mea, past_ba_low = ba_low, past_ba_hig = ba_hig)
  
  # 2. Keep predicted data long, with indice preserved
  pred_dat <- p_base100 %>%
    filter(type == "pred") %>%
    select(year, espar, type, indice, pred_ba_mea = ba_mea, pred_ba_low = ba_low, pred_ba_hig = ba_hig)
  
  # 3. Join predictions to observations by year + espar
  overlap_summary <- pred_dat %>%
    left_join(past_dat, by = c("year", "espar")) %>%
    mutate(
      obs_se  = (past_ba_hig - past_ba_low) / (2 * 1.96),
      pred_se = (pred_ba_hig - pred_ba_low) / (2 * 1.96),
      se_diff = sqrt(obs_se^2 + pred_se^2),
      z_diff  = abs(past_ba_mea - pred_ba_mea) / se_diff,
      aligned = z_diff <= 1.96
    )
  
  # 4. Summarize per espar AND indice (i.e. per model, per species)
  overlap_counts <- overlap_summary %>%
    group_by(espar, indice) %>%
    summarise(
      n_years      = n(),
      n_aligned    = sum(aligned, na.rm = TRUE),
      prop_aligned = n_aligned / n_years,
      mean_z       = mean(z_diff, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(overlap = ifelse(prop_aligned >= 0.5, "yes", "no")) %>%
    group_by(espar) %>%
    slice_max(prop_aligned, n = 1)
  
  
  pred_dat %>%
    select(year, espar, type, indice, ba_mea = pred_ba_mea,
           ba_low = pred_ba_low, ba_hig = pred_ba_hig) %>%
    bind_rows(select(past_dat, year, espar, type, ba_mea = past_ba_mea,
                     ba_low = past_ba_low, ba_hig = past_ba_hig)) %>%
    inner_join(lbls, by = "espar") %>%
    mutate(type2 = case_when(type == "past" ~ "observed",
                             type == "pred" ~ "predicted"),
           species_short = gsub("\\b([A-Za-z]{3})[A-Za-z]{2,}\\b", "\\1.", species),
           ba_low = ifelse(ba_low < 0, 0, ba_low)) -> p_base100f
  
   # Build a lookup of species_short -> overlap (one value per species)
  overlap_counts %>%
    inner_join(lbls, by = "espar") %>%
    arrange(species) -> strip_lookup
  
  # Map overlap to colors
  strip_colors <- ifelse(strip_lookup$overlap == "yes", "blue", "orange")  
  
  gg <- ggplot(p_base100f, aes(x = year + 2.5, y = ba_mea,
                              color = type2, fill = type2)) +
    geom_ribbon(aes(ymin = ba_low, ymax = ba_hig), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    facet_wrap2(~species_short, scales = "free_y",
                strip = strip_themed(
                  background_x = elem_list_rect(fill = strip_colors)
                )) +
    labs(
      x     = "Year",
      y     = "Mortality (base 1 in first year)",
      color = NULL, fill = NULL
    ) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 12, color = "white"),
          axis.title = element_text(size = 20))
  
  path <- file.path(output_dir, "mort_trend.png")
  ggsave(path, gg, width = 10, height = 7, dpi = 150)
  return(path)
  
}

## plot spatial anomalies
plot_ser_anom <- function(data, model, espar, indice,
                          ser_path = "../IFN_stuff/data/geodata/ser_27572.gpkg") {
  if (is.null(model)) return(NULL)
  
  Xvar <- str_split_1(indice, pattern = "\\+|\\*")
  
  sub_dat <- data |>
    dplyr::filter(espar.x == espar) |>
    dplyr::select(all_of(Xvar), mort, ser_86, incref, npp, gtot, cdim, w) |>
    tidyr::drop_na(all_of(Xvar)) 
  sub_dat$ser_id <- as.integer(as.factor(sub_dat$ser_86))
  
  ser_id <- sub_dat |>
    dplyr::select(ser_86, ser_id) |>
    dplyr::group_by(ser_86) |>
    dplyr::slice_head(n = 1)
  
  # Random effects
  ser_val <- model$summary.random[["ser_id"]]
  ser_val$prob <- plogis(ser_val$mean + model$summary.fixed$mean[1])
  
  ss <- dplyr::inner_join(ser_id, ser_val[, c("ID", "prob")],
                          by = c("ser_id" = "ID"))
  
  ser_sf <- sf::st_read(ser_path, quiet = TRUE)
  
  ser_anom <- ss |>
    dplyr::right_join(ser_sf, by = c("ser_86" = "codeser")) |>
    sf::st_as_sf()
  
  # Species label for the panel title
  lbls <- read.csv("~/LIF/IFN_stuff/data/pommier/espar_species.csv")
  sp_label <- lbls$species[lbls$espar == model$.meta$espar]
  
  p <- ggplot(ser_anom, aes(fill = prob, geometry = geom)) +
    geom_sf(color = "grey80", linewidth = 0.1) +
    scale_fill_gradient2(
      midpoint = mean(ser_anom$prob, na.rm = TRUE),
      low      = scales::muted("blue"),
      high     = scales::muted("red"),
      name     = "Mortality\nprobability",
      na.value = "grey90"
    ) +
    labs(title = sp_label) +    # species name as panel title
    theme_void() +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.position = "bottom",
      legend.key.width = unit(1.2, "cm"),
      legend.title    = element_text(size = 7),
      legend.text     = element_text(size = 6)
    )
  
  attr(p, "espar") <- espar
  return(p)
}

## plot temporal anomalies
plot_incref_anom <- function(data, model, espar, indice){
  if (is.null(model)) return(NULL)
  
  Xvar <- str_split_1(indice, pattern = "\\+|\\*")
  
  sub_dat <- data |>
    dplyr::filter(espar.x == espar) |>
    dplyr::select(all_of(Xvar), mort, ser_86, incref, npp, gtot, cdim, w) |>
    tidyr::drop_na(all_of(Xvar)) 
  sub_dat$incref_id <- as.integer(as.factor(sub_dat$incref))
  
  incref_id <- sub_dat |>
    dplyr::select(incref, incref_id) |>
    dplyr::group_by(incref) |>
    dplyr::slice_head(n = 1)
  
  # Random effects
  sim <- inla.posterior.sample(n = 250, model)
  
  intercept_samp <- sapply(
    sim,
    function(x) x$latent[rownames(x$latent) == "(Intercept):1", 1]
  )
  
  # RW2 samples
  rw2_names <- paste0("increff:", rownames(model$summary.random$increff))
  rw2_samp <- sapply(
    sim,
    function(x) x$latent[rw2_names, 1]
  )
  
  # Linear predictor
  eta <- sweep(rw2_samp, 2, intercept_samp, "+")
  
  # Inverse-logit
  prob <- plogis(eta)
  
  # Summaries
  prob_mean <- apply(prob, 1, mean)
  prob_ci <- t(apply(prob, 1, quantile, probs = c(0.1, 0.9)))
  
  ## put together in a ofd
  incref_val <- data.frame(ID = as.integer(rownames(model$summary.random$increff)),
                           prob = prob_mean,
                           prob_low = prob_ci[,1],
                           prob_hig = prob_ci[,2])
  
  ## add the ser id
  ss <- inner_join(incref_id,
                   incref_val[,c("ID", "prob",
                                 "prob_low", "prob_hig")],
                   by = c("incref_id" = "ID")) %>%
    filter(incref < 12)
  
  # Species label for the panel title
  lbls <- read.csv("~/LIF/IFN_stuff/data/pommier/espar_species.csv")
  sp_label <- lbls$species[lbls$espar == model$.meta$espar]
  
  p <- ggplot(ss, aes(x=incref+2005,
                      y = prob,
                      ymin = prob_low,
                      ymax = prob_hig)) +
    geom_ribbon(alpha = 0.1) +
    geom_path() +
    labs(title = sp_label,
         x = "Year",
         y = "Mortality probability") +
    theme(text = element_text(size = 15))
  
  attr(p, "espar") <- espar
  return(p)
}

## plot size-growth mortality probability
plot_size_growth <- function(data, all_summaries, output_dir = "figures") {
  
  # ── 1. Extract coefficients, keyed by term name ───────────────────────────
  coefs <- all_summaries %>%
    filter(indice == "abg+gtot_0+ba_tot") %>%
    select(term, mean, espar) %>%
    group_by(espar) %>%
    summarise(
      coef = list(setNames(mean, term)),  # named numeric vector → safe lookup
      .groups = "drop"
    )
  
  # ── 2. Per-species gradient ranges + scaling parameters ───────────────────
  gradients <- data %>%
    group_by(espar = espar.x) %>%
    mutate(
      size   = as.numeric(scale(gtot_0)),  # same scaling as model fitting
      growth = as.numeric(scale(abg))
    ) %>%
    summarise(
      size_min   = quantile(size,   probs = 0.05),
      size_max   = quantile(size,   probs = 0.95),
      growth_min = quantile(growth, probs = 0.05),
      growth_max = quantile(growth, probs = 0.95),
      size_mean   = mean(gtot_0),  # scaling params on original scale
      size_sd     = sd(gtot_0),
      growth_mean = mean(abg),
      growth_sd   = sd(abg),
      .groups = "drop"
    )
  
  # ── 3. Build prediction grid, compute and back-transform ──────────────────
  pred_data <- gradients %>%
    inner_join(coefs, by = "espar") %>%
    rowwise() %>%
    mutate(
      grid = list({
        
        # Back-transform quantiles to original scale to define log-spaced seq
        size_orig_min   <- size_min   * size_sd   + size_mean
        size_orig_max   <- size_max   * size_sd   + size_mean
        growth_orig_min <- growth_min * growth_sd + growth_mean
        growth_orig_max <- growth_max * growth_sd + growth_mean
        
        # Generate grid in original scale, log-spaced for even visual spacing
        g <- expand.grid(
          size_orig   = exp(seq(log(size_orig_min),   log(size_orig_max),   length.out = 30)),
          growth_orig = exp(seq(log(growth_orig_min), log(growth_orig_max), length.out = 30))
        )
        
        # Re-scale to model scale for prediction
        g$size   <- (g$size_orig   - size_mean)   / size_sd
        g$growth <- (g$growth_orig - growth_mean) / growth_sd
        
        b0     <- coef[["(Intercept)"]]
        b_abg  <- coef[["abg"]]
        b_size <- coef[["gtot_0"]]
        
        g$pred <- arm::invlogit(b0 + b_abg * g$growth + b_size * g$size)
        
        g
      })
    ) %>%
    ungroup() %>%
    select(espar, grid) %>%
    tidyr::unnest(grid)
  
  # ── 4. Plot on original log scale ─────────────────────────────────────────
  p <- ggplot(pred_data, aes(x = size_orig, y = growth_orig, fill = pred)) +
    geom_tile() +
    facet_wrap(~espar, scales = "free") +
    scale_x_log10() +
    scale_y_log10() +
    scale_fill_distiller(
      palette = "RdBu",
      name   = "P(death)"
    ) +
    labs(
      x     = "Basal area (log scale)",
      y     = "Growth (log scale)"
    ) +
    theme_bw()
  
  # ── 5. Save ───────────────────────────────────────────────────────────────
  out_path <- file.path(output_dir, "size_growth_predictions.png")
  ggsave(out_path, p, width = 12, height = 8, dpi = 150)
  return(out_path)
}

## a function to combine the ser anomalies map
combine_ser_maps <- function(..., ncol = 4, filter_espar = NULL, output_dir = "figures") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  plots <- purrr::compact(list(...))  # drops NULLs (non-matching indices)
  # Filter by espar if requested
  if (!is.null(filter_espar)) {
    plots <- purrr::keep(plots, ~attr(.x, "espar") %in% filter_espar)
  }
  
  # Sort by espar attribute
  espar_order <- purrr::map_chr(plots, ~attr(.x, "espar"))
  plots <- plots[order(espar_order)]
  
  if (length(plots) == 0) {
    warning("No plots to combine.")
    return(NULL)
  }
  
  final <- patchwork::wrap_plots(plots, ncol = ncol)
  
  path <- file.path(output_dir, paste0("ser_map",
                                       paste(filter_espar, collapse = "-"),
                                       ".png"))
  ggsave(path, final,
         width  = 4 * ncol,
         height = ceiling(length(plots) / ncol) * 5,
         dpi    = 150)
  return(path)
}

## a function to synthetize anomalies ser maps
synth_ser_map <- function(..., filter_espar = NULL, output_dir = "figures"){
  plots <- purrr::compact(list(...))  # drops NULLs (non-matching indices)
  # Filter by espar if requested
  if (!is.null(filter_espar)) {
    plots <- purrr::keep(plots, ~attr(.x, "espar") %in% filter_espar)
  }
  
  # extract plot data
  plot_dat <- purrr::map(plots, ~ {
    dat <-.x$data
    dat$espar <- attr(.x, "espar")
    dat
  }) |>
    list_rbind()
  
  # count per ser the nb of species with mortality prob below the average
  plot_dat %>%
    group_by(espar) %>%
    mutate(avg_prob  = mean(prob, na.rm = TRUE),
           is_neg = prob > avg_prob) %>%
    ungroup() %>%
    group_by(ser_86, geom) %>%
    summarise(nb_neg = sum(is_neg, na.rm = TRUE),
              nb_sp = sum(!is.na(prob))) %>%
    mutate(p_sp = nb_neg / nb_sp) -> syn
  
  final <- ggplot() +
    geom_sf(data = syn, aes(fill = p_sp * 100, geometry = geom)) +
    scale_fill_gradient(low = "grey90",
                        high = scales::muted("red"),
                        name = "Percent of species") +
    theme_void()
  
  path <- file.path(output_dir, "ser_synth.png")
  ggsave(path, final, dpi = 150)
  return(path)
  
}

## a function to combine incref anomalies
combine_incref_anom <- function(..., ncol = 4, output_dir = "figures"){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  plots <- purrr::compact(list(...))  # drops NULLs (non-matching indices)

  # Sort by espar attribute
  espar_order <- purrr::map_chr(plots, ~attr(.x, "espar"))
  plots <- plots[order(espar_order)]
  
  if (length(plots) == 0) {
    warning("No plots to combine.")
    return(NULL)
  }
  
  final <- patchwork::wrap_plots(plots, ncol = ncol,
                                 axes = "collect",
                                 axes_title = "collect") 
  
  path <- file.path(output_dir, paste0("incref_anom.png"))
  ggsave(path, final,
         width  = 4 * ncol,
         height = ceiling(length(plots) / ncol) * 5,
         dpi    = 150)
  return(path)
}

## helper function to compute integrate
int_f <- function(x, mu1, mu2, sd1, sd2) {
  f1 <- dnorm(x, mean=mu1, sd=sd1)
  f2 <- dnorm(x, mean=mu2, sd=sd2)
  pmin(f1, f2)
}

## identify sp with deviating mortality preds
get_sp_deviation <- function(all_pred, mort_trend){
  ## data wraggling for the ocre output
  mort_trend %>%
    filter(visite == "1") %>%
    select(year, espar.espar1, tot, moy_ha, ic_moy_ha) %>%
    mutate(ba_low = moy_ha-ic_moy_ha,
           ba_hig = moy_ha+ic_moy_ha,
           type = "past") %>%
    select(espar.espar1, year, moy_ha,
           ba_low, ba_hig, type) %>%
    rename(espar = espar.espar1,
           ba_mea = moy_ha) -> p_mod
  
  ## data wraggling for the model outputs
  lbls <- read.csv("~/LIF/IFN_stuff/data/pommier/espar_species.csv")
  all_pred %>%
    select(-indice) %>%
    mutate(type = "pred") %>%
    bind_rows(p_mod) %>%
    inner_join(lbls,
               by = "espar") -> p
  
  ## version en base 100
  p_base100 <- p |>
    group_by(espar, type) |>
    arrange(year) |>
    mutate(
      ref     = first(ba_mea),          # valeur de référence = première année
      ba_mea  = ba_mea  / ref,
      ba_low  = ba_low  / ref,    # IC bas divisé par la même référence
      ba_hig  = ba_hig  / ref,
      sd = (ba_mea - ba_low) / 1.96
    ) |>
    ungroup() 
  
  ##check overlap between confidence intervals
  overlap_summary <- p_base100 %>%
    select(year, espar, type, ba_mea, ba_low, ba_hig) %>%
    pivot_wider(
      names_from = type,
      values_from = c(ba_mea, ba_low, ba_hig),
      names_glue = "{type}_{.value}"
    ) %>%
    mutate(
      obs_se  = (past_ba_hig  - past_ba_low)  / (2 * 1.96),
      pred_se = (pred_ba_hig - pred_ba_low) / (2 * 1.96),
      se_diff = sqrt(obs_se^2 + pred_se^2),
      z_diff  = abs(past_ba_mea - pred_ba_mea) / se_diff,
      aligned = z_diff <= 1.96   # means are within ~1 combined SE of 95% CI width
    )
  
  overlap_counts <- overlap_summary %>%
    group_by(espar) %>%
    summarise(
      n_years      = n(),
      n_aligned    = sum(aligned, na.rm = TRUE),
      prop_aligned = n_aligned / n_years,
      mean_z       = mean(z_diff, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(prop_aligned >= 0.5)
  
  return(overlap_counts$espar)
}
  
## classify nowcasted trends
get_class_trends <- function(all_pred, mort_trend, sp_dev, output_dir = "tables"){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  lbls <- read.csv("~/LIF/IFN_stuff/data/pommier/espar_species.csv") 

  ## linear trend in observed data
  mort_trend %>%
    filter(visite == "1") %>%
    select(espar.espar1, year, moy_ha) %>%
    nest(data = -espar.espar1) %>%
    mutate(model = map(data, ~ lm(moy_ha~year, .x))) %>%
    mutate(modelt = map(model, tidy)) %>%
    unnest(modelt) %>%
    filter(term == "year") %>%
    mutate(class = case_when(p.value > 0.05 ~ "stable",
                             p.value <= 0.05 & estimate > 0 ~ "increasing",
                             p.value <= 0.05 & estimate < 0 ~ "decreasing",
                             .default = NA)) %>%
    inner_join(lbls, by = c("espar.espar1" = "espar")) %>%
    filter(espar.espar1 %in% sp_dev) %>%
    select(species, class) -> tab1
  
  ## linear trends in predictions
  all_pred %>%
    nest(data = -espar) %>%
    mutate(model = map(data, ~ lm(ba_mea~year, .x))) %>%
    mutate(modelt = map(model, tidy)) %>%
    unnest(modelt) %>%
    filter(term == "year") %>%
    mutate(class = case_when(p.value > 0.05 ~ "stable",
                             p.value <= 0.05 & estimate > 0 ~ "increasing",
                             p.value <= 0.05 & estimate < 0 ~ "decreasing",
                             .default = NA)) %>%
    inner_join(lbls, by = "espar") %>%
    filter(espar %in% sp_dev) %>%
    select(species, class) -> tab2
  

  tabb <- inner_join(tab1, tab2, by = "species") %>%
    rename(`Trend obs. 2019-2010` = class.x,
           `Trend pred. 2024-2010` = class.y)
  
  path <- file.path(output_dir, "trend_class.csv")
  write.csv(tabb, path, row.names = FALSE)
  return(path)
  
}

## make summary stats table
get_summary_table <- function(data, output_dir = "tables"){
  
  lbls <- read.csv("~/LIF/IFN_stuff/data/pommier/espar_species.csv") 
  
  data %>%
    group_by(espar.x) %>%
    summarise(nb_plot = length(unique(npp)),
              nb_alive = sum(type == "alive"),
              nb_dead = sum(type == "dead")) %>%
    mutate(prop_dead = round((nb_dead / (nb_dead + nb_alive)) * 100, 2)) %>%
    inner_join(lbls, by = c("espar.x" = "espar")) %>%
    arrange(espar.x) %>%
    select(species, nb_plot, nb_alive, nb_dead, prop_dead) -> out
  
  
  path <- file.path(output_dir, "summary_tab.csv")
  write.csv2(out, path, row.names = FALSE)
  return(path)
  
}

## get bias between observed and predicted trends
get_bias_trends <- function(all_pred, mort_trend){
  ## data wraggling for the ocre output
  mort_trend %>%
    filter(visite == "1") %>%
    rename(espar = espar.espar1,
           ba_obs = moy_ha) %>%
    select(espar, year, ba_obs)-> p_mod
  
  ## data wraggling for the model outputs
  all_pred %>%
    select(espar, year, ba_mea) %>%
    inner_join(p_mod, by = c("espar", "year")) %>%
    mutate(bias = ba_mea - ba_obs) %>%
    group_by(espar) %>%
    summarise(bias_med = median(bias)) -> tt

  return(tt)
}

