## Funciton to get the tree-level data

## species making up x% of the wood ressource
sel_sp <- function(thresh = 0.8){
  require(data.table)
  require(inventR)
  require(dplyr)
  
  dat <- ocre_nm(ans = 2010:2019,
                 us_nm = "1",
                 vars = "V",
                 vent_a = "espar")
  
  
  dat |>
    filter(visite == "1") |>
    select(espar.espar1, tot) |>
    arrange(desc(tot)) |>
    mutate(cumm = cumsum(tot),
           cump = cumm / last(cumm)) |>
    filter(cump <= thresh) -> dd
  
  return(as.character(dd$espar.espar1))
}

## gather tree-level data
get_data_tree <- function(sel){
  ## connect to db
  conn <- inventR::connect_db()
  
  ## grab growth and veget5 infos
  g3a <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "g3arbre"))
  g3a5 <- tbl(conn, dbplyr::in_schema("prod_exp", "g3arbre"))
  g3aa <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "g3agedom"))
  g3f <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "g3foret"))
  e2p <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "e2point"))
  e1p <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "e1point"))
  
  
  ## grab the tree-level data
  g3a %>%
    filter(espar %in% sel) %>%
    full_join(g3aa, by = c("npp", "a", "incref")) %>%
    group_by(npp) %>%
    mutate(ba_tot = sum(gtot * wac)) %>%
    ungroup() %>%
    filter(simplif == "0", ir5 > 0, incref > 4) %>%
    mutate(type = case_when(veget5 %in% c("0") ~ "alive",
                            veget5 %in% c("M") ~ "dead",
                            incref > 14 ~ "pred",
                            .default = NA)) %>%
    filter(!is.na(type)) %>%
    select(incref, npp, a, ba_tot, c13.x, espar.x, ori,
           htot, gtot, ir5, abg, type, age13, lib,
           w, wac, contains("ai")) %>%
    mutate(ai2 = ai2 - ai1,
           ai3 = ai3 - ai2,
           ai4 = ai4 - ai3,
           ai5 = ir5 - ai4) %>%
    left_join(g3a5, by = c("npp", "a")) %>%
    mutate(abg2 = ((c135 ** 2 - c13.x ** 2) / (4*pi)),
           gtot_0 = gtot - (abg * 5)) %>%
    mutate(abg2 = ifelse(abg2 <= 0, NA, abg2),
           gtot_0 = ifelse(gtot_0 <= 0, NA, gtot_0),
           abg2c = ifelse(type == "dead", (abg2) / 0.5, abg2),
           relbai = abg / gtot, # wunder 2008
           relbac = ((gtot / gtot_0) ** 0.2) - 1) %>% # hülseman 2018
    mutate(abg_diff = (abg2-abg) / abg,
           lib = case_when(lib == "0" ~ "sup",
                           lib %in% c("1", "2") ~ "dom",
                           .default = NA)) %>%
    group_by(npp, incref) %>%
    mutate(abg_pl = abg / mean(abg, na.rm=TRUE)) %>%
    select(npp, incref, a, espar.x, ori, lib, htot, age13, w, wac, type, ba_tot, gtot, gtot_0, abg, abg2,
           abg_diff, abg_pl, contains("ai"), ir5, abg2c, relbai, relbac) %>%
    collect() %>%
    ungroup() -> dd
  
  ## compute trend and var in ai
  # Assume your data frame is called df, with columns ai1..ai5
  mat <- as.matrix(dd[!is.na(dd$ai5), c("ai1","ai2","ai3","ai4","ai5")])
  
  ## --- 1. Linear trend (slope of growth over the 5 years) ---
  
  x <- 1:5
  xbar <- mean(x)
  Sxx <- sum((x - xbar)^2)
  
  # Slope for each row using the OLS formula:
  # slope = sum((x - xbar) * y) / Sxx   (the cross term with ybar vanishes)
  df <- dd[!is.na(dd$ai5),c("npp", "a")]
  df$ir_trend <- as.numeric(mat %*% (x - xbar)) / Sxx
  
  # Optional: intercept too, if you want predicted values / fitted line
  # ybar <- rowMeans(mat)
  # df$intercept <- ybar - df$trend * xbar
  
  ## --- 2. Variability (CV) of growth across the 5 years ---
  
  row_mean <- rowMeans(mat)
  df$ir_sd <- sqrt(rowSums((mat - row_mean)^2) / (ncol(mat) - 1))
    
  dd <- left_join(dd, df, by = c("npp", "a"))
  
  dd %>%
    group_by(npp, incref) %>%
    mutate(ba_larger = sapply(gtot, function(x) sum(gtot[gtot >= x]) - x)) -> dd
  
  
  ## add plot-level infos
  g3f %>%
    inner_join(e2p, by = "npp") %>%
    inner_join(select(e1p, npp, idp), by = "npp") %>%
    select(npp, idp, greco, ser_86, pro_psg, sfo_nm, sver, poids) %>%
    filter(npp %in% local(dd$npp)) %>%
    collect() %>%
    inner_join(dd, by = "npp") -> dat
  
  ## create response variables
  dat$mort <- case_when(dat$type == "alive" ~ 0,
                        dat$type == "dead" ~ 1,
                        dat$type == "pred" ~ NA)
  
  dat <- mutate(dat, diam = 200 * sqrt((gtot / pi)),
                cdim = case_when(diam < 22.5 ~ "pb",
                                 diam < 37.5 ~ "mb",
                                 .default = "gb"))
  
  return(dat)
}

## function to save tree-level data
save_data_tree <- function(data){
  
  ## replace w with wac and drop and give dummy plot weights
  data %>%
    mutate(w = wac,
           poids = 2) %>%
    select(-wac, -npp) -> save_dat
  
  path <- file.path("data", "tree_df.parquet")
  arrow::write_parquet(save_dat, path)
  return(path)
}

## generate mort trends from ocre
get_data_trend <- function(sel){
  
  
  yy <- data.table(year = 2010:2019,
                   veget5 = c("veget5:veget56:(M)","veget5:veget56:(M)",
                              "veget5:veget57:(M)","veget5:veget57:(M)",
                              rep("veget5:veget59:(M)", 6)))
  p_trend <- yy[,ocre_nm(ans = year, vars = "GTOT",
                         moy_ha = TRUE,
                         vent_a = c(veget5,
                                    paste0("espar:espar1:(",paste(sel, collapse = ","),")"))),
                by = c("year", "veget5")]
  
  return(p_trend)
}

## gather traits data
get_data_traits <- function(sel){
  ## litt traits
  tr_tab <- read.csv("~/LIF/Mortalite/data/preprocessed/tr_tab.csv")
  xdm_tab <- read.csv("~/LIF/Mortalite/data/xdm_species_density_data.csv")
  
  ## connect to db
  conn <- inventR::connect_db()
  
  ## grab growth and veget5 infos
  g3a <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "g3arbre"))
  g3aa <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "g3agedom"))
  
  ## compute db traits, average height, average age,
  ## average growth, average height/dbh
  g3a %>%
    left_join(g3aa, by = c("npp", "a", "espar")) %>%
    select(npp, a, espar, d13, ir5, htot, age13) %>%
    group_by(espar) %>%
    summarise(growth = mean(ir5 * 1000, na.rm = TRUE),
              height = mean(htot, na.rm = TRUE),
              age = mean(age13, na.rm = TRUE),
              slender = mean(htot / d13, na.rm = TRUE)) %>%
    collect() -> tr_db
  
  ## join the tables
  tr_tab %>%
    full_join(select(xdm_tab,Species, Mean), by = c("species" = "Species")) %>%
    full_join(tr_db, by = "espar") %>%
    filter(espar %in% sel) -> tr_tabf
  
  names(tr_tabf)[c(4, 5, 9, 13)] <- c("Frost_intolerance",
                                     "Drought_intolerance",
                                     "Wood_density",
                                     "Slenderness")

  return(tr_tabf)
}


