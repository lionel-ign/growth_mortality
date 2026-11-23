####################################
# Script to reproduce the analysis #
#      Author: Lionel Hertzog      #
####################################

# This script uses targets to run the whole analysis reported
# in the manuscript

# load libraries
library(targets)
library(tarchetypes)
library(dplyr)
library(crew)

# define to crew, one fast and one slow for the heavy computation
controllers <- list(
  crew_controller_local(name = "fast", workers = 4, seconds_idle = 30),
  crew_controller_local(name = "slow", workers = 1, seconds_idle = 60)
)

# set options to run the pipeline
tar_option_set(
  packages = c("tibble", "INLA", "tidyverse", "sf", "broom",
               "data.table", "inventR", "pROC", "MLmetrics",
               "ggcorrplot", "ggh4x", "tidyr"),
  format = "qs",
  error = "null",
  controller = crew_controller_group(controllers)
)

# read the script files
tar_source()

# Définir les combinaisons espèce × indice
tgt_ess <-  c("03", "02", "09", "61", "62", "52", "51",
              "10", "64", "11", "05", "17C", "12V", 
              "53CO", "06", "54", "14", "13G", "24",
              "15S", "22M", "63", "21C", "57A", "73", "20P")

combinations <- tidyr::expand_grid(
  espar = tgt_ess,
  indice = c("abg+gtot_0+ba_tot",
             "abg+gtot_0+ba_tot+lib",
             "abg2+gtot+ba_tot",
             "abg2c+gtot+ba_tot",
             "abg_pl+gtot_0+ba_tot",
             "abg_diff+gtot+ba_tot",
             "ir_trend+gtot_0+ba_tot",
             "ir_sd+gtot_0+ba_tot",
             "relbai+gtot_0+ba_tot",
             "relbac+gtot_0+ba_tot")
)

## define nicer labelling for indices
lbls <- combinations |>
  group_by(indice) |>
  slice_head(n=1) |>
  mutate(indiced = case_when(
    indice == "abg+gtot_0+ba_tot+lib" ~ "BAI1+\nDOM",
    indice == "abg+gtot_0+ba_tot"     ~ "BAI1",
    indice == "abg2+gtot+ba_tot"      ~ "BAI2",
    indice == "abg2c+gtot+ba_tot"     ~ "BAI2c",
    indice == "abg_diff+gtot+ba_tot"  ~ "dBAI",
    indice == "abg_pl+gtot_0+ba_tot"  ~ "BAI1plot",
    indice == "ir_trend+gtot_0+ba_tot" ~ "RItrend",
    indice == "ir_sd+gtot_0+ba_tot" ~ "RIsd",
    indice == "relbai+gtot_0+ba_tot" ~ "relBAI",
    indice == "relbac+gtot_0+ba_tot" ~ "BAIc",
    .default = indice
  )) |>
  select(indice, indiced)

## read-in the data
here::i_am("_targets.R")
data <- arrow::read_feather("data/tree_df.feather")
traits <- read.csv("data/tree_traits.csv")
mort_trend <- read.csv("data/mort_trend.csv")



# dispatch the model and post-processing per species and indice
mapped <- tar_map(
  values = combinations,
  names  = c(espar, indice),
  
  tar_target(
    name    = model,
    command = fit_model_tree(data, espar, indice)
  ), # fit the models
  tar_target(
    name  = fit,
    command = extract_fit(data, model, espar, indice)
    ),  # extract fitness metrics
  tar_target(
    name    = summary,
    command = extract_summary(model, espar, indice)
  ),    # extract model summaries
  tar_target(
    name = pred,
  command = if (indice %in% c("abg_pl+gtot_0+ba_tot")) {
    extract_pred(data, model, espar, indice, nSim = 1000)
  } else {
    NULL
  },   #extract predictions for the selected indice
  resources = tar_resources(
    crew = tar_resources_crew(controller = "slow")
  )),
  tar_target(
    name = ser_map,
    command = if (indice == "abg_pl+gtot_0+ba_tot") {
      plot_ser_anom(data, model, espar, indice)
    } else {
      NULL
    }), # get map of mortality anomalies by ser
  tar_target(
    name = incref_anom,
    command = if (indice == "abg_pl+gtot_0+ba_tot") {
      plot_incref_anom(data, model, espar, indice)
    } else {
      NULL
    }) # get temporal anomalies
)


## get data, actual mapped run and concatenate outputs
list(
  # 1. create data, for information only, the files are saved in the github data/ subfolder
  # tar_target(name = data, command = get_data_tree(tgt_ess)),
  # tar_target(name = traits, command = get_data_traits(tgt_ess)),
  # tar_target(name = mort_trend, command = get_data_trend(tgt_ess)),
  
  # 2. run models
  mapped,  # le tar_map stocké dans une variable
  
  # 3. combine extracted infos from models
  tar_combine(
    name    = all_summaries,
    mapped[["summary"]],   
    command = dplyr::bind_rows(!!!.x)
  ),
  tar_combine(
    name    = all_fit,
    mapped[["fit"]],   
    command = dplyr::bind_rows(!!!.x)
  ),
  tar_combine(
    name    = all_pred,
    mapped[["pred"]],  
    command = dplyr::bind_rows(!!!.x)
  ),
  tar_target(
    name = espar_bias,
    command = get_pred_bias(all_pred, mort_trend)
  ), # compute bias between prediction and observation
  # 4. plot extracted infos
  tar_target(
    name = plot_growth,
    command = plot_growth_status(data),
    format = "file"
  ),
  tar_target(
    name = growth_corr,
    command = plot_growth_corr(data),
    format = "file"
  ), # correlation matrix between growth indices
  tar_target(
    name = plot_trait,
    command = plot_trait_corr(traits, all_fit),
    format = "file"
  ), # correlation matrix between traits
  tar_target(
    name = plot_espar_metrics,
    command = plot_fitespar_metrics(all_fit, all_summaries, sp_bias, lbls),
    format = "file"
  ), # dotplot of model fitness per species
  tar_target(
    name = plot_metrics,
    command = plot_fit_metrics(all_fit, lbls),
    format = "file"
  ), # plot fitness metrics aggregated
  tar_target(
    name = plot_dom,
    command = plot_resp_dom(data, all_summaries, traits),
    format = "file"
  ), # plot effect of tree dominance
  tar_target(
    name = plot_effects,
    command = plot_abgint_effects(all_summaries),
    format = "file"
  ), 
  tar_target(
    name = plot_corr_trait,
    command = plot_corr_traits(all_fit, sp_bias, traits),
    format = "file"
  ), # correlation matrix between model fitness and species traits
  tar_target(
    name = plot_trend,
    command = plot_mort_trends(mort_trend, all_pred),
    format = "file"
  ),
  # Combine all ser_map panels into one figure — format = "file"
  tar_combine(
    ser_maps_grp1,
    mapped[["ser_map"]],
    command = combine_ser_maps(!!!.x, ncol = 3,
                               filter_espar = c("02", "03", "05", "06",
                                                "09", "10","11", "12V", "17C")),
    format  = "file",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slow")
    )
  ),
  tar_combine(
    ser_maps_grp2,
    mapped[["ser_map"]],
    command = combine_ser_maps(!!!.x, ncol = 3,
                               filter_espar = c("51", "53CO", "52", "62", "57A",
                                                "64", "61","54", "73", "63")),
    format  = "file",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slow")
    )
  ),
  tar_combine(
    ser_maps_grp3,
    mapped[["ser_map"]],
    command = combine_ser_maps(!!!.x, ncol = 4,
                               filter_espar = c("13G", "14", "15S", "20P",
                                                "21C", "22M", "24")),
    format  = "file",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slow")
    )
  ),
  tar_combine(
    ser_synth,
    mapped[["ser_map"]],
    command = synth_ser_map(!!!.x),
    format  = "file",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slow")
    )
  ),
  # Combine all incref_anom panels into one figure — format = "file"
  tar_combine(
    incref_anom_grp,
    mapped[["incref_anom"]],
    command = combine_incref_anom(!!!.x, ncol = 5),
    format  = "file",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slow")
    )
  ),
  tar_target(
    name = class_trend,
    command = get_class_trends(all_pred, mort_trend, sp_dev),
    format = "file"
  ), # classify mortality trends
  tar_target(
    name = sp_dev,
    command = get_sp_deviation(all_pred, mort_trend)
  ), 
  tar_target(
    name = sp_bias,
    command = get_bias_trends(all_pred, mort_trend)
  ), # get bias between predicted and observed mortality
  tar_target(
    name = summary_tab,
    command = get_summary_table(data),
    format = "file"
  )
)