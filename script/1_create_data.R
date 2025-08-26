## Create data script

library(tidyverse)
library(inventR)
## connect to db
conn <- inventR::connect_db()

## grab growth and veget5 infos
g3a <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "g3arbre"))
g3a5 <- tbl(conn, dbplyr::in_schema("prod_exp", "g3arbre"))
g3aa <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "g3agedom"))
g3f <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "g3foret"))
e2p <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "e2point"))

## dev compute linear slope
tt <- as.data.frame(matrix(1:5 + rnorm(100), ncol=5))
tt %>%
  rowwise() %>%
  mutate(sd_y = sd(c_across(V1:V5)),
         sd_x = sd(1:5),
         r = cov(c_across(V1:V5), 1:5) / (sd_x*sd_y),
         b = (r*sd_y)/sd_x)

## grab the tree-level data
g3a %>%
  full_join(g3aa, by = c("npp", "a", "incref")) %>%
  filter(simplif == "0", ir5 > 0, incref > 4, incref < 15) %>%
  mutate(type = case_when(veget5 %in% c("0", "1") ~ "alive",
                          veget5 %in% c("5", "M", "2") ~ "dead",
                          .default = NA)) %>%
  filter(!is.na(type)) %>%
  select(incref, npp, a, c13.x, ess, lib, htot, gtot, ir5, abg, type, age13, lib, wac) %>%
  inner_join(g3a5, by = c("npp", "a")) %>%
  mutate(abg2 = 0.2 * ((c135 ** 2 - c13.x ** 2) / (4*pi)),
         gtot_0 = gtot - (abg * 5)) %>%
  mutate(abg2 = ifelse(abg2 <= 0, NA, abg2),
         gtot_0 = ifelse(gtot_0 <= 0, NA, gtot_0)) %>%
  mutate(abgg = abg / gtot_0,
         abgg2 = abg2 / gtot,
         abg_diff = (abg2-abg) / abg) %>%
  group_by(npp) %>%
  mutate(abg_pl = abg / mean(abg, na.rm=TRUE),
         abg2_pl = abg2 / mean(abg2, na.rm = TRUE)) %>%
  select(npp, a, ess, lib, htot, age13, wac, type, gtot, gtot_0, abg, abg2,
         abgg, abgg2, abg_diff, abg_pl, abg2_pl, starts_with("ir")) %>%
  collect() %>%
  ungroup() %>%
  rowwise() %>%
  mutate(sd_y = sd(c_across(starts_with("ir"))),
         sd_x = sd(1:5),
         r = cov(c_across(starts_with("ir")), 1:5) / (sd_x*sd_y),
         ir_trend = (r*sd_y)/sd_x,
         ir_cv = mean(c_across(starts_with("ir")))/sd_y) %>%
  ungroup() -> dd

## add plot-level infos
g3f %>%
  inner_join(e2p, by = "npp") %>%
  select(npp, greco, ser_86, pro_nm, pro_psg, sfo_nm, sver) %>%
  filter(npp %in% local(dd$npp)) %>%
  collect() %>%
  inner_join(dd, by = "npp") -> dat

## remove uneeded columns
dat <- select(dat, -ir1, -ir2, -ir3, -ir4, -ir5, -starts_with("sd"), -r)

## save this
write.csv(dat, "~/LIF/Mortalite/data/maindat.csv", row.names = FALSE)


## species selection
## select the 15 most common species
sel_ess <- ocre_nm(ans = 2005:2024, vars = "V", vent_a = "ess")[]
sel_ess <- arrange(sel_ess, desc(tot))  
target_ess <- data.frame(ess = sel_ess$ess[c(1:15)])
## now with name labs
target_ess <- inner_join(target_ess, ListDataMod("ess"), by = c("ess" = "mode"))
## add wd infos
ww <- readxl::read_excel("LIF/Mortalite/data/WoodDesigns_XDM_pourHenri.xlsx")
names(ww) <- c("ess", "type")
target_ess$definition[c(1, 2, 4, 5, 11, 12, 13, 14, 15)] <- c("Quercus robur",
                                                              "Quercus petraea",
                                                              "Abies alba", "Picea abies",
                                                              "Quercus pubescens",
                                                              "Fraxinus excelsior",
                                                              "Betula pubescens",
                                                              "Pinus nigra var. corsicana",
                                                              "Pinus nigra subsp. nigra")
target_ess <- inner_join(target_ess, ww, by = c("definition" = "ess"))
target_ess$type[target_ess$definition == "Fagus sylvatica"] <- "RP"

write.csv(target_ess, "~/LIF/Mortalite/data/sp_sel.csv", row.names = FALSE)
