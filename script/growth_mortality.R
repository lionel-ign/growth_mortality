# dead tree growing
library(tidyverse)
library(patchwork)
library(broom)
library(transport)
library(inventR)
## connect to db
conn <- inventR::connect_db()

## grab growth and veget5 infos
g3a <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "g3arbre"))
g3a5 <- tbl(conn, dbplyr::in_schema("prod_exp", "g3arbre"))
g3aa <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "g3agedom"))
g3f <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "g3foret"))
e2p <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "e2point"))

## grab the data
g3a %>%
  full_join(g3aa, by = c("npp", "a", "incref")) %>%
  filter(simplif == "0", ir5 > 0, incref > 4, incref < 15) %>%
  mutate(type = case_when(veget5 %in% c("0", "1") ~ "alive",
                          veget5 %in% c("5", "M", "2") ~ "dead",
                          .default = NA)) %>%
  filter(!is.na(type)) %>%
  select(incref, npp, a, c13.x, ess, htot, gtot, ir5, abg, type, age13, lib, wac) %>%
  inner_join(g3a5, by = c("npp", "a")) %>%
  mutate(abg2 = 0.2 * ((c135 ** 2 - c13.x ** 2) / (4*pi)),
         gtot_0 = gtot - (abg * 5)) %>%
  mutate(abg2 = ifelse(abg2 < 0, NA, abg2),
         gtot_0 = ifelse(gtot_0 <= 0, NA, gtot_0)) %>%
  mutate(abgg = abg / gtot_0,
         abgg2 = abg2 / gtot) %>%
  collect() -> dd

g3f %>%
  inner_join(e2p, by = "npp") %>%
  select(npp, greco, ser_86, pro_nm, pro_psg, sfo_nm, sver) %>%
  filter(npp %in% local(dd$npp)) %>%
  collect() %>%
  inner_join(dd, by = "npp") -> dd

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

## descriptive stats
### cumulative prop of species
dd %>%
  count(ess) %>%
  arrange(desc(n)) %>%
  mutate(cu = cumsum(n) / sum(n),
         rank = 1:n()) -> rr
ggplot(rr, aes(x=rank, y=cu)) +
  geom_line() +
  geom_hline(yintercept = 0.9) +
  geom_vline(xintercept = 23)
### only keep species making up 90% of the set
ddf <- filter(dd, ess %in% target_ess$ess)

### compare distrib
g1 <- ggplot(ddf, aes(x=gtot, fill=type)) +
  geom_density(alpha=0.5) +
  scale_x_log10(expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(x = "Basal area (m²)")

g2 <- ggplot(ddf, aes(x=htot, fill=type)) +
  geom_density(alpha=0.5) +
  scale_x_log10(expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(x = "Hauteur (m)")

g3 <- ggplot(ddf, aes(x=abg, fill=type)) +
  geom_density(alpha=0.5) +
  scale_x_log10(expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(x = "Basal area\nincrement (m²/year)")

g4 <- ggplot(ddf, aes(x=abgg * 100, fill=type)) +
  geom_density(alpha=0.5) +
  scale_x_log10(expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(x = "Rate of basal\narea increment (%/year)")

g5 <- ggplot(ddf, aes(x=abg_prop * 100, fill=type)) +
  geom_density(alpha=0.5) +
  scale_x_continuous(expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(x = "Rate of basal\narea increment (%/year)")

ga <- patchwork::wrap_plots(g1, g2, g3, g4, ncol=2, guides = "collect")

ggsave("LIF/Mortalite/figures/descri.png", ga)

## impact of access to light
ddf %>%
  count(type, lib) %>%
  group_by(lib) %>%
  mutate(prop = n / sum(n)) %>%
  filter(type == "dead") %>%
  ungroup() %>%
  mutate(ll = c("no access to light",
                "some access to light",
                "full access to light")) -> ll
gl <- ggplot(ll, aes(x=prop*100, y=fct_reorder(ll, prop))) +
  geom_bar(stat="identity", width = 0.5) +
  theme_classic() +
  scale_x_continuous(expand = expansion()) +
  labs(x = "Percent of dead trees per light access class",
       y = "")

## impact of property
ddf %>%
  mutate(pro = case_when(pro_nm %in% c("1O", "1S", "2O", "CO", "2S", "2M") ~ "public",
                          pro_psg == "1" ~ "private_psg",
                          pro_psg == "2" ~ "privat_no_psg")) %>%
  count(type, pro)   %>%
  group_by(pro) %>%
  mutate(prop = n / sum(n)) %>%
  filter(type == "dead") -> ddp

gp <- ggplot(ddp, aes(x=prop*100, y=fct_reorder(pro, prop))) +
  geom_bar(stat="identity", width = 0.5) +
  theme_classic() +
  scale_x_continuous(expand = expansion()) +
  labs(x = "Percent of dead trees per property class",
       y = "")

glp <- gl+gp+theme(element_text(size=20))
ggsave("LIF/Mortalite/figures/prop_prolib.png", glp, width = 8, height = 6)

## plot per species the dg ~ g curves
dd %>%
  select(ess, type, abg, abg2, gtot_0, gtot) %>%
  rename(abg_t0 = abg,
         gtot_t1 = gtot,
         abg_t1 = abg2,
         gtot_t0 = gtot_0) %>%
  pivot_longer(3:6, names_to = c(".value", "time"),
               names_sep = "_") -> dde

for(ess in target_ess$ess){
  ## remove outliers in g
  dde %>%
    filter(ess == ess) %>%
    filter(gtot < quantile(gtot, probs = 0.9, na.rm = TRUE),
           abg < quantile(abg, probs = 0.9, na.rm = TRUE)) -> ddf
  g <- ggplot(ddf[ddf$ess == ess,], aes(x = gtot * 1e4, y=abg * 1e4, color=type,
                                        linetype = time)) +
    stat_smooth(linewidth = 1.25) +
    labs(title = paste("Species: ", target_ess$libelle[target_ess$ess == ess]),
         x = "Basal area (cm²)",
         y = "Basal area increment (cm²/year)")
  
  ggsave(paste0("~/LIF/Mortalite/figures/ess/gtot_abg_abs_", ess, ".png"), g)
}

## explore imbalance between live and dead trees
ee <- ListDataMod("ess")
dd %>%
  count(ess, type) %>%
  pivot_wider(names_from = type,
              values_from = n) %>%
  filter(!is.na(dead),
         alive + dead >= 1000) %>%
  mutate(prop_dead = dead / (dead + alive)) %>%
  inner_join(ee, by = c("ess" = "mode")) %>%
  ggplot(aes(x=fct_reorder(libelle, prop_dead), y=prop_dead)) +
  geom_point() +
  coord_flip()

## now make a logistic model to predict death from basal area growth

## new dev with bootstrapp and reweighting
foo <- function(data, model){
  ## create a binary response
  data$mort <- ifelse(data$type == "alive", 0, 1)
  
  ## create size index
  data <- mutate(data, class = case_when(gtot < (pi * (0.225 / 2) ** 2) ~ "pb",
                                         gtot >= (pi * (0.225 / 2) ** 2) & gtot < (pi * (0.375 / 2) ** 2) ~ "mb",
                                         gtot >= (pi * (0.375 / 2) ** 2) ~ "gb"))
  
  ## scale the variables
  data <- mutate(data, across(gtot:abgg2, scale))
  


  set.seed(123)
  train_index <- caret::createDataPartition(data$mort, p = 0.8, list = FALSE)
  train <- data[train_index, ]
  test <- data[-train_index, ]
  
  n_boot <- 100
  metrics <- data.frame(accuracy = numeric(7*n_boot),
                        specificity = numeric(7*n_boot),
                        sensitivity = numeric(7*n_boot),
                        alpha = numeric(7*n_boot),
                        beta = numeric(7*n_boot),
                        gamma = numeric(7*n_boot),
                        class = rep(c("all", "gb", "mb", "pb",
                                      "lib_0", "lib_1", "lib_2"), n_boot),
                        id_boot = rep(1:n_boot, each = 7))
  
  for (i in 1:n_boot) {
    # Bootstrap sample
    boot_indices <- sample(nrow(train), replace = TRUE)
    boot_train <- train[boot_indices, ]
    
    ## get weights
    class_weights <- table(boot_train$mort)
    inv_weights <- ifelse(boot_train$mort == 1,
                          1 / class_weights[2],  # minority (dead)
                          1 / class_weights[1])
    
    # Fit weighted model
    if(model == "abg*gtot"){
      m <- glm(mort ~ gtot_0 * abg, train, family = "binomial")
    } else if(model == "abg"){
      m <- glm(mort ~ abg, train, family = "binomial")
    } else if(model == "abg+gtot"){
      m <- glm(mort ~ gtot_0 + abg, boot_train, family = binomial(link="probit"),
               weights = inv_weights)
    } else if(model == "abg_prop"){
      m <- glm(mort ~ abg_prop, train, family = "binomial")
    } else if(model == "abgg"){
      m <- glm(mort ~ abgg, train, family = "binomial")
    } else if(model == "abg2"){
      m <- glm(mort ~ abg2, train, family = "binomial")
    } else if(model == "abg2+gtot"){
      m <- glm(mort ~ gtot + abg2, boot_train, family = binomial(link = "probit"),
               weights = inv_weights)
    } else if(model == "abgg2"){
      m <- glm(mort ~ abgg2, train, family = "binomial")
    } else if(model == "abgt"){
      m <- glm(mort ~ abgt, train, family = "binomial")
    } else if(model == "abgt+gtot"){
      m <- glm(mort ~ abgt + gtot_0, train, family = "binomial")
    }

    
    # Predict on fixed test set
    p <- predict(m, newdata = test, type = "response")
    pred_class <- rbinom(nrow(test), size = 1, prob = p)
    
    # Confusion matrix and metrics for all and per size classes
    cml <- list("all" = table(Pred = pred_class, True = test$mort),
                "gb" = table(Pred = pred_class[test$class == "gb"],
                             True = test$mort[test$class == "gb"]),
                "mb" = table(Pred = pred_class[test$class == "mb"],
                             True = test$mort[test$class == "mb"]),
                "pb" = table(Pred = pred_class[test$class == "pb"],
                             True = test$mort[test$class == "pb"]),
                "lib_0" = table(Pred = pred_class[test$lib == "0"],
                                True = test$mort[test$lib == "0"]),
                "lib_1" = table(Pred = pred_class[test$lib == "1"],
                                True = test$mort[test$lib == "1"]),
                "lib_2" = table(Pred = pred_class[test$lib == "2"],
                                True = test$mort[test$lib == "2"]))
    
    # Defensive code in case a class is missing
    TP <- sapply(cml, function(x) tryCatch(x["1", "1"], error = function(e) 0))
    TN <- sapply(cml, function(x) tryCatch(x["0", "0"], error = function(e) 0))
    FP <- sapply(cml, function(x) tryCatch(x["1", "0"], error = function(e) 0))
    FN <- sapply(cml, function(x) tryCatch(x["0", "1"], error = function(e) 0))
    
    metrics$accuracy[metrics$id_boot == i] <- (TP + TN) / sapply(cml, sum)
    metrics$sensitivity[metrics$id_boot == i] <- TP / (TP + FN)
    metrics$specificity[metrics$id_boot == i] <- TN / (TN + FP)
    metrics$alpha[metrics$id_boot == i] <- coef(m)["(Intercept)"]
    metrics$beta[metrics$id_boot == i] <- coef(m)[grep("abg",attr(coef(m), "names"))]
    metrics$gamma[metrics$id_boot == i] <- coef(m)[grep("gtot",attr(coef(m), "names"))]
  }
  ## get summary stats
  metrics %>%
    pivot_longer(accuracy:gamma) %>%
    group_by(name, class) %>%
    summarise(me = mean(value),
              lo = quantile(value, probs = 0.05, na.rm = TRUE),
              hi = quantile(value, probs = 0.95, na.rm = TRUE)) -> out
  
  return(out)
}


## do this across species
dd %>%
  mutate(abgt = abg + abg2) %>%
  filter(ess %in% target_ess$ess) %>%
  select(ess, type,lib, gtot, gtot_0, abg, abg2, abgt, abgg, abgg2) -> df
  
ddd <- rbind(df, df)
ddd$model <- rep(c("abg+gtot", 
                   "abg2+gtot"), each = nrow(df))
ddd %>%
  nest(data = -c(ess, model)) %>%
  mutate(aa = map2(data, model, function(x, y) foo(x, y))) %>%
  unnest(aa) -> res

res <- inner_join(res, target_ess, by = c("ess"))
res$class <- factor(res$class, labels = c("all tree","large tree", "shaded",
                                          "codominant", "dominant", "medium tree",
                                          "small tree"))
## first fig tree size
fit <- filter(res, name %in% c("accuracy", "sensitivity", "specificity"),
              class %in% c("all tree", "large tree", "medium tree", "small tree"))
fit$model <- factor(fit$model, labels = c("n-10 <-> n-5",
                                          "n-5 <-> n"))
gb <- ggplot(fit, aes(x=class, y=me, color=model)) +
  stat_summary(fun.data = "mean_cl_boot") +
  geom_jitter(alpha = 0.25, width = 0.1) +
  facet_wrap(vars(name), nrow = 3) +
  geom_hline(yintercept = 0.5) +
  coord_flip() +
  theme(text = element_text(size = 20))  +
  labs(x = "Tree size",
       y = "Model fitness")

ggsave("~/LIF/Mortalite/figures/mode_size.png", gb, height = 8, width = 6)

# now with tree access to light
fit <- filter(res, name %in% c("accuracy", "sensitivity", "specificity"),
              class %in% c("all tree", "dominant", "codominant", "shaded"))
fit$model <- factor(fit$model, labels = c("n-10 <-> n-5",
                                          "n-5 <-> n"))
gb <- ggplot(fit, aes(x=class, y=me, color=model)) +
  stat_summary(fun.data = "mean_cl_boot") +
  geom_jitter(alpha = 0.25, width = 0.1) +
  facet_wrap(vars(name), nrow = 3) +
  geom_hline(yintercept = 0.5) +
  coord_flip() +
  theme(text = element_text(size = 20))  +
  labs(x = "Tree size",
       y = "Model fitness")

ggsave("~/LIF/Mortalite/figures/mode_light.png", gb, height = 8, width = 6)

# now with species ring types
fit <- filter(res, name %in% c("accuracy", "sensitivity", "specificity"),
              class == "all tree")
fit$model <- factor(fit$model, labels = c("n-10 <-> n-5",
                                          "n-5 <-> n"))
fit$type <- factor(fit$type, labels = c("coniferous", "diffuse\nporous", "ring\nporous"))
gb <- ggplot(fit, aes(x=type, y=me, color=model)) +
  stat_summary(fun.data = "mean_cl_boot") +
  geom_jitter(alpha = 0.25, width = 0.1) +
  facet_wrap(vars(name), nrow = 3) +
  geom_hline(yintercept = 0.5) +
  coord_flip() +
  theme(text = element_text(size = 20))  +
  labs(x = "Wood structure",
       y = "Model fitness")

ggsave("~/LIF/Mortalite/figures/mode_woodstr.png", gb, height = 8, width = 6)



ga <- ggplot(fit, aes(x=fct_reorder(libelle, me), y=me,
                ymin=lo, ymax=hi, color = class)) +
  geom_linerange(position = position_dodge2(width=0.5)) +
  geom_point(position = position_dodge2(width=0.5)) +
  facet_wrap(vars(name, model), scales="free_x", ncol = 2) +
  coord_flip() +
  labs(x = "",
       y = "Model fitness") +
  geom_hline(yintercept = 0.5) +
  scale_color_discrete(name = "Tree size")

ggsave("LIF/Mortalite/figures/acc_spec_sens_size.png", ga)

gb <- ggplot(fit, aes(x=class, y=me, color=model)) +
  stat_summary(fun.data = "mean_cl_boot") +
  geom_jitter(alpha = 0.1) +
  facet_wrap(vars(name))                                                                                           
                                                                                              
          
ggplot(fit, aes(x=type, y=me, color=model)) +
  facet_wrap(vars(name)) +
  geom_hline(yintercept = 0.5) +
  theme(axis.text.x = element_text(angle = 45)) +  geom_jitter(width = 0.1) +
  
  scale_color_discrete(name = "Growth period") +
  labs(x = "Tree size", y = "Model fitness")                                                                                     


ggsave("LIF/Mortalite/figures/acc_spec_sens_jitter.png", gb)

ggplot(subset(fit, class %in% c("shaded", "codominant", "dominant") & name == "accuracy"),
       aes(x=fct_reorder(libelle, me), y=me,
           ymin=lo, ymax=hi, color = model, shape = class)) +
  geom_linerange(position = position_dodge2(width=0.5)) +
  geom_point(position = position_dodge2(width=0.5)) +
  coord_flip() +
  labs(x = "",
       y = "Model fitness") +
  geom_hline(yintercept = 0.5) +
  scale_color_discrete(name = "Tree size")

ggplot(subset(fit, class %in% c("shaded", "codominant", "dominant") & name == "specificity"),
       aes(x=fct_reorder(libelle, me), y=me,
           ymin=lo, ymax=hi, color = model, shape = class)) +
  geom_linerange(position = position_dodge2(width=0.5)) +
  geom_point(position = position_dodge2(width=0.5)) +
  coord_flip() +
  labs(x = "",
       y = "Model fitness") +
  geom_hline(yintercept = 0.5) +
  scale_color_discrete(name = "Tree size")

## now plot the model coefs
res %>%
  mutate(model = ifelse(model == "abg+gtot", "m1", "m2")) %>%
  filter(name %in% c("alpha", "beta")) %>%
  select(ess, type, libelle, model, name, me) %>%
  pivot_wider(names_from = c(model, name),
              values_from = me) -> cc

ggplot(cc, aes(x=m1_alpha, y=m1_beta,
               xend = m2_alpha, yend = m2_beta,
               color = libelle)) +
  geom_segment()

## now plot the predicted curves
res %>%
  select(class, ess, libelle, model, name, me) %>%
  filter(name %in% c("alpha", "beta", "gamma"),
         class == "all tree") %>%
  pivot_wider(names_from = name,
              values_from = me) %>%
  rowwise() %>%
  mutate(xrange = list(seq(-2, 2, length.out = 50)),
         pred = list(arm::invlogit(alpha + beta * xrange))) %>%
  unnest(c(xrange, pred)) -> pp

pp$model <- factor(pp$model, labels = c("n-10 <-> n-5",
                                        "n-5 <-> n"))

gc <- ggplot(pp, aes(x=xrange, y=pred, color=model)) +
  geom_path() +
  facet_wrap(vars(libelle)) +
  labs(x = "Basal area growth (standardized)",
       y = "P(Death)") +
  scale_color_discrete(name = "Growth period")

ggsave("LIF/Mortalite/figures/pred_curves.png", gc)


df %>%
  group_by(ess) %>%
  summarise(gavg = mean(gtot),
            xmin = quantile(abg, probs = 0.1),
            xmax = quantile(abg, probs = 0.9)) %>%
  inner_join(pre, by = "ess") %>%
  rowwise() %>%
  mutate(xrange = list(seq(xmin, xmax, length.out = 50)),
         pred = list(arm::invlogit(alpha + gamma * gavg + beta * xrange))) -> tt


gse <- ggplot(subset(res, term == "sensitivity"), aes(x=fct_reorder(libelle, me), y=me,
                                                  ymin=lo, ymax=hi, color = model)) +
  geom_linerange(position = position_dodge2(width=0.5)) +
  geom_point(position = position_dodge2(width=0.5)) +
  facet_grid(vars(size), scales="free") +
  coord_flip() +
  labs(x = "",
       y = "Sensitivity")

gsp <- ggplot(subset(res, term == "specificity"), aes(x=fct_reorder(libelle, me), y=me,
                                                      ymin=lo, ymax=hi, color = model)) +
  geom_linerange(position = position_dodge2(width=0.5)) +
  geom_point(position = position_dodge2(width=0.5)) +
  facet_grid(vars(size), scales="free") +
  coord_flip() +
  labs(x = "",
       y = "Specificity")


ggsave("~/LIF/Mortalite/figures/model_res.png", gr, width=6, height=8)



ddf %>%
  filter(ess == "02") %>%
  select(gtot, abg, lib, type) %>%
  mutate(libf = factor(lib, labels = c("no", "little", "full"))) -> d02


# rebalance

# fit a model
m0 <- glm(mort ~ gtot_std + libf, undersampled, family = "binomial")
p_m0 <- rbinom(n = nrow(test), size = 1,
             p = predict(m0, newdata=test, type="response"))
caret::confusionMatrix(factor(p_m0),
                       factor(test$mort))

# then another one with growth
m1 <- glm(mort ~ gtot_std + libf + abg_std, undersampled, family = "binomial")
p_m1 <- rbinom(n = nrow(test), size = 1,
               p = predict(m1, newdata=test, type="response"))
caret::confusionMatrix(factor(p_m1),
                       factor(test$mort))

# a more complex model with interaction with growth
m2 <- glm(mort ~ gtot_std + libf * abg_std, undersampled, family = "binomial")
p_m2 <- rbinom(n = nrow(test), size = 1,
               p = predict(m2, newdata=test, type="response"))
caret::confusionMatrix(factor(p_m2),
                       factor(test$mort))

## now rf models
undersampled$mortf <- factor(undersampled$mort)
r0 <- randomForest(mortf ~ gtot_std + libf, data = undersampled,
                   importance = TRUE)
caret::confusionMatrix(predict(r0, newdata = test),
                       factor(test$mort))

## with growth
r1 <- randomForest(mortf ~ gtot_std + libf + abg_std, data = undersampled,
                   importance = TRUE)
caret::confusionMatrix(predict(r1, newdata = test),
                       factor(test$mort))


m_base <- glm(mort ~ abg_std, dat_train, family = "binomial")
pp <- rbinom(n = nrow(dat_test), size = 1,
             p = predict(m_base, newdata=dat_test, type="response"))
caret::confusionMatrix(factor(pp),
                       factor(dat_test$mort))
## 93% accuracy, but super low specificity

## rebalance dataset
ddf <- select(ddf, ess, gtot, abg, type)
### try for one species
d02 <- filter(ddf, ess == "02")
ii <- pdfos(as.data.frame(d02[,-1]), numInstances = 2000, classAttr = "type")
nd <- rbind(d02[,-1], ii)
plotComparison(as.data.frame(d02[,-1]), nd, attrs = c("gtot", "abg"),
               classAttr = "type")

## now fit a model with gtot as covariable
m_g <- glm(mort ~ abgg, dat_train, family = "binomial")
## predict with constant gtot
dat_g <- dat_test
dat_g$gtot_std <- 0
pp <- rbinom(n = nrow(dat_test), size = 1,
             p = predict(m_g, newdata=dat_test, type="response"))
caret::confusionMatrix(factor(pp),
                       factor(dat_test$mort))



g1 <- ggplot(dds, aes(x=abg, fill=type)) +
  geom_density(alpha=0.25) +
  facet_wrap(vars(libelle)) +
  scale_x_log10() +
  labs(x = "Basal area increment (m²/year)")

g2 <- ggplot(dds, aes(x=gtot, fill=type)) +
  geom_density(alpha=0.25) +
  facet_wrap(vars(libelle)) +
  scale_x_log10() +
  labs(x = "Basal area (m²)")

g3 <- ggplot(dds, aes(x=(abg/gtot)*100, fill=type)) +
  geom_density(alpha=0.25) +
  facet_wrap(vars(libelle)) +
  scale_x_log10() +
  labs(x = "Rate of basal area increment (%/year)")

ga <- wrap_plots(list(g2, g1, g3), nrow=3) + plot_layout(guides = "collect")
ggsave("LIF/Mortalite/figures/dens_ess.png", ga, width=8, height=16)


g1 <- ggplot(dd, aes(x=abg, fill=type)) +
  geom_density(alpha=0.25) +
  facet_wrap(vars(greco)) +
  scale_x_log10() +
  labs(x = "Basal area increment (m²/year)")

g2 <- ggplot(dd, aes(x=gtot, fill=type)) +
  geom_density(alpha=0.25) +
  facet_wrap(vars(greco)) +
  scale_x_log10() +
  labs(x = "Basal area (m²)")

g3 <- ggplot(dd, aes(x=(abg/gtot)*100, fill=type)) +
  geom_density(alpha=0.25) +
  facet_wrap(vars(greco)) +
  scale_x_log10() +
  labs(x = "Rate of basal area increment (%/year)")

ga <- wrap_plots(list(g2, g1, g3), nrow=3) + plot_layout(guides = "collect")
ggsave("LIF/Mortalite/figures/dens_greco.png", ga, width=8, height=16)




## plot gtot vs growth per type
for(ess in sel_ess$ess[1:23]){
  ## remove outliers in g
  dd %>%
    filter(ess == ess) %>%
    filter(gtot < quantile(gtot, probs = 0.9),
           abg < quantile(abg, probs = 0.9)) -> ddf
  g <- ggplot(ddf[ddf$ess == ess,], aes(x = gtot, y=abg, color=type,
                                         fill = type, alpha = type)) +
    geom_jitter() +
    stat_smooth(linewidth = 1.25) +
    scale_alpha_manual(values = c("dead" = 0.1, "alive" = 0.01)) +
    labs(title = paste("Species: ", ess))
  
  ggsave(paste0("~/LIF/Mortalite/figures/ess/gtot_abg_abs_", ess, ".png"), g)
}
ggplot(ddf[ddf$ess == "02",], aes(x = gtot, y=abg, color=type,
                                  fill = type, alpha = type)) +
  geom_jitter() +
  stat_smooth(linewidth = 1.25) +
  scale_x_log10() +
  scale_y_log10() +
  scale_alpha_manual(values = c("dead" = 0.1, "alive" = 0.01))


## analysis integral between dead and alive trees per size classes
### dev first for one species
d02 <- filter(dd, ess == "02")
d02 %>%
  select(gtot, abg, type) %>%
  nest(data = -type) %>%
  mutate(mo = map(data, ~ gam(abg ~ s(gtot, bs = "cs"), data = .x, method = "REML"))) %>%
  mutate(pred = map2(mo, data, ~ predict(.x,
                                  newdata = data.frame(gtot = seq(min(.y$gtot),
                                                                  max(.y$gtot),
                                                                  length.out = 100)))),
         gtot = map(data, ~ seq(min(.x$gtot),
                                max(.x$gtot),
                                length.out = 100))) %>%
  unnest(c(pred, gtot)) %>%
  select(type, pred, gtot) -> pp

ggplot(pp, aes(x = gtot, y = pred, color = type)) +
  geom_path() +
  scale_x_log10() +
  scale_y_log10()


get_area <- function(dat){
  x1 <- dat$gtot[dat$type =="alive"]
  x2 <- dat$gtot[dat$type == "dead"]
  y1 <- dat$pred[dat$type =="alive"]
  y2 <- dat$pred[dat$type == "dead"]
  
  # Create a common x grid (e.g., 100 points between min and max of overlap)
  x_common <- seq(max(min(x1), min(x2)), min(max(x1), max(x2)), length.out = 100)
  
  # Interpolate y1 and y2 on the common x
  y1_interp <- approx(x1, y1, xout = x_common)$y
  y2_interp <- approx(x2, y2, xout = x_common)$y
  
  # Compute the vertical difference
  diff_y <- y1_interp - y2_interp
  
  # Compute area between the curves
  area <- trapz(x_common, diff_y)
  
  return(area)
}

## do this per size classes
df <- data.frame(class = c("PB", "MB", "GB"),
                 limit_sup = c(pi * (0.225 / 2) ** 2,
                               pi * (0.375 / 2) ** 2,
                               max(d02$gtot)),
                 limit_inf = c(0,
                               pi * (0.225 / 2) ** 2,
                               pi * (0.375 / 2) ** 2),
                 area = NA)
for(i in 1:3){
  df$area[i] <- get_area(pp[pp$gtot < df$limit_sup[i] &
                            pp$gtot >= df$limit_inf[i],])
}


## estimate coeffs of growth ~ gtot model
foo_lm <- function(data){
  ## std the vars
  #data <- mutate(data, across(gtot:abg2, scale))
  m1 <- lm(log(abg+1e-6) ~ type * log(gtot_0), data)
  m2 <- lm(log(abg2+1e-6) ~ type * log(gtot), data)
  
  ## bring things together
  out <- data.frame(coef = c(coef(m1), coef(m2)),
             low = c(confint(m1)[,1], confint(m2)[,1]),
             hig = c(confint(m1)[,2], confint(m2)[,2]),
             name = rep(c("int", "int_mort", "slp", "slp_mort"), 2),
             model = rep(c("initial", "revisit"), each = 4))
  
  return(out)
}


df %>%
  select(ess, type, gtot:abg2) %>%
  filter(ess != "06") %>%
  nest(data = -ess) %>%
  mutate(mm = map(data, ~foo_lm(.x))) %>%
  unnest(mm) -> rr

ggplot(rr, aes(x=fct_reorder(ess, coef), y=coef, ymin=low, ymax=hig,
               color=model)) +
  geom_linerange() +
  geom_point() +
  facet_wrap(vars(name), scales="free") +
  coord_flip()

ww <- readxl::read_excel("LIF/Mortalite/data/WoodDesigns_XDM_pourHenri.xlsx")
names(ww) <- c("ess", "type")
target_ess$definition[c(2, 6, 7, 10, 11, 12, 13, 14)] <- c("Quercus pubescens",
                                                           "Quercus petraea",
                                                           "Quercus robur",
                                                           "Fraxinus excelsior",
                                                           "Picea abies",
                                                           "Abies alba",
                                                           "Betula pubescens",
                                                           "Acer campestre")
rr %>%
  inner_join(target_ess, by = "ess") %>%
  inner_join(ww, by = c("definition" = "ess")) -> rrw

ggplot(rrw, aes(x=type, y=coef, color=model)) +
  geom_jitter(width=0.1)+
  facet_wrap(vars(name), scales="free")

rrw %>%
  select(type, definition, coef, name, model) %>%
  pivot_wider(names_from = name, values_from = coef) %>%
  mutate(int_m = int+int_mort,
         slp_m = slp+slp_mort) -> rrww

ggplot(rrww, aes(x=type, y=int_m, color=model)) +
  geom_jitter(width=0.1)

ggplot(rrww, aes(x=int_m, y=slp_m, color=type)) +
  geom_jitter(width=0.1) +
  facet_wrap(vars(model))


## try out estimation of monod function on the relation dg~g
dd %>%
  filter(ess == "02",
         type == "alive",
         gtot_0 >= quantile(gtot_0, probs = 0.1, na.rm=TRUE),
         gtot_0 <= quantile(gtot_0, probs = 0.9, na.rm=TRUE)) %>%
  select(abg, gtot_0, ess, abg2, gtot) %>%
  mutate(abgs = scale(abg)[,1],
         gtots = scale(gtot_0)[,1]) -> tt

LLmonod <- function(x, y, a = 0.2, b = 1, sigma = 0.1){
  -sum(stats::dnorm(y, mean = a * (x / (x + b)), sd = sigma), log = TRUE)
}

mm <- mle2(LLmonod,
           start = list(a = 0.2, b = -0.5 , sigma = 1),
           method = "L-BFGS-B",
           optimizer = "optim",
           data = list(x = tt$gtot_0, y=tt$abg),
           lower = c(sigma = 1e-6))


prior1 <- prior(normal(0, 2), nlpar = "alpha") +
  prior(normal(0, 2), nlpar = "beta")

mn2 <- nls(abg ~ SSmicmen(gtot_0, 0.02, 0.1), data = tt)
mn2 <- nls(abg2 ~ a * (gtot /(gtot + b)), data = na.omit(tt), start = list( a = 0.2, b = 1))

ggplot(tt, aes(x=gtot_0, y=abg)) +
  geom_smooth(method="nls", 
              formula=y~ a * (x / (x + b)), # this is an nls argument
              method.args = list(start=c(a = 0.2, b = 1)), # this too
              se=FALSE) +
  stat_smooth(color="red")

## a helper function
foo_nls <- function(data){
  data <- na.omit(data)
  start_pars <- list(a = max(data$abg),
                     b = median(data$gtot_0))
  m <- tryCatch(nlsLM(abg ~ a * (gtot_0 / (gtot_0 + b)), data = data, start = start_pars),
                error = function(e) NA)
  if(class(m) == "nls"){
    cc_int <- tryCatch(confint(m),
                       error = function(e) as.data.frame(matrix(NA, ncol=2, nrow=2)))
    r2 <- cor(data$abg, predict(m)) ** 2
    tt <- as.data.frame(cbind(coef(m), cc_int))
    names(tt) <- c("est", "low", "hig")
    tt$coef <- c("saturation", "half-rate")
    tt$r2 <- r2
  }
  else NA 
  
  return(tt)
}

## go through the species
dd %>%
  filter(ess %in% target_ess$ess) %>%
  select(ess, type, gtot_0, abg) %>%
  nest(data = -c(ess, type)) %>%
  mutate(nlm = map(data, ~foo_nls(.x))) %>%
  unnest(nlm) -> mm

mm %>%
  filter(!is.na(est)) %>%
  select(ess, type, est:r2) %>%
  inner_join(ListDataMod("ess"), by = c("ess" = "mode")) -> md

gd <- ggplot(md, aes(x=fct_reorder(libelle, est), y=est, ymin=low, ymax=hig, color = type,
               group = type)) +
  geom_linerange(position = position_dodge2(width=0.25)) +
  geom_point(position = position_dodge2(width=0.25)) +
  facet_wrap(vars(coef), scales="free") +
  coord_flip() +
  scale_y_log10() +
  labs(x = "",
       y = "Estimate (95% CI)")

ggsave("~/LIF/Mortalite/figures/mm_fit.png", gd)

## now with wood str
ww <- readxl::read_xlsx("~/LIF/Mortalite/data/WoodDesigns_XDM_pourHenri.xlsx")

## compute deltaG along the gtot gradient



## try a linear log-log model
foo_ll <- function(data){
  data <- na.omit(data)
  m <- tryCatch(lm(log(abg) ~ log(gtot_0), data = data),
                error = function(e) NA)
  if(class(m) == "lm"){
    cc_int <- tryCatch(confint(m),
                       error = function(e) as.data.frame(matrix(NA, ncol=2, nrow=2)))
    r2 <- cor(data$abg, predict(m)) ** 2
    tt <- as.data.frame(cbind(coef(m), cc_int))
    names(tt) <- c("est", "low", "hig")
    tt$coef <- c("intercept", "slope")
    tt$r2 <- r2
  }
  else NA 
  
  return(tt)
}

## go through the species
dd %>%
  filter(ess %in% target_ess$ess) %>%
  select(ess, type, gtot_0, abg) %>%
  nest(data = -c(ess, type)) %>%
  mutate(nlm = map(data, ~foo_ll(.x))) %>%
  unnest(nlm) -> mm

mm %>%
  select(ess, type, est:r2) %>%
  inner_join(ListDataMod("ess"), by = c("ess" = "mode")) -> md

ggplot(md, aes(x=fct_reorder(libelle, est), y=est, ymin=low, ymax=hig, color = type)) +
  geom_linerange() +
  geom_point() +
  facet_wrap(vars(coef), scales="free") +
  coord_flip()


nls(abg ~ a * (gtot / (gtot + b)), data = tt, start = list(a = 0.01, b = 0.1))


newd <- data.frame(gtot_0 = seq(min(tt$gtot_0),
                                max(tt$gtot_0),
                                length.out = 50))
newd$abg <- predict(mn, newdata = newd)

ggplot(tt, aes(x=gtot_0, y=abg)) +
  geom_point(alpha=0.01) +
  geom_path(data = newd, color="red")

## get nrmse
nrmse <- sqrt(mean((tt$abg - predict(mn))**2)) / mean(tt$abg)

## short-term forecasting of tree survival / death
## for now growing stock of tree that survived given low sensitivity of models

## derive ocre estimation
## separate call per veget5 definition
yy <- data.table(year = 2010:2019)

p1 <- yy[1:2,ocre_nm(ans = year, var = "V",
        vent_a = c(paste0("ess:ess:(",paste(target_ess$ess, collapse = ","),"):TRUE"),
                   "veget5:veget56:0:FALSE"),
        moy_ha = TRUE), by = year]
p2 <- yy[3:4,ocre_nm(ans = year, var = "V",
                     vent_a = c(paste0("ess:ess:(",paste(target_ess$ess, collapse = ","),"):TRUE"),
                                "veget5:veget57:0:FALSE"),
                     moy_ha = TRUE), by = year]
p3 <- yy[5:10,ocre_nm(ans = year, var = "V",
                       vent_a = c(paste0("ess:ess:(",paste(target_ess$ess, collapse = ","),"):TRUE"),
                                  "veget5:veget59:(0, 1):FALSE"),
                       moy_ha = TRUE), by = year]
## combine
pa <- rbind(p1, p2, p3)
pa %>%
  filter(visite == "1") %>%
  select(year,ess.ess, moy_ha, ic_moy_ha) -> paf


## fit the logistic model to all species
foo_model <- function(data){
  ## create a binary response
  data$mort <- ifelse(data$type == "alive", 0, 1)
  
  ## scale the variables
  data$gtot_0 <- scale(data$gtot_0)[,1]
  data$abg <- scale(data$abg)[,1]
  ## get weights
  class_weights <- table(data$mort)
  inv_weights <- ifelse(data$mort == 1,
                          1 / class_weights[2],  # minority (dead)
                          1 / class_weights[1])
    
  # Fit weighted model
  m <- glm(mort ~ gtot_0 + abg, data, family = binomial(link="probit"),
               weights = inv_weights)
  
  return(m)
}

dd %>%
  select(ess, type, gtot_0, abg) %>%
  filter(ess %in% target_ess$ess) %>%
  nest(data = -ess) %>%
  mutate(model = map(data, ~foo_model(.x))) -> m_ess

## grab the tree-level data from after the 5-year period
g3a %>%
  filter(incref >= 15,
         ess %in% target_ess$ess) %>%
  select(npp, incref, w, wac, ess, abg, gtot, v) %>%
  mutate(gtot_0 = gtot - (abg * 5)) %>%
  collect() -> new_tree

## grab scale infos
dd %>%
  filter(ess %in% target_ess$ess) %>%
  group_by(ess) %>%
  summarise(abg_m = mean(abg, na.rm = TRUE),
            abg_sd = sd(abg, na.rm = TRUE),
            gtot_m = mean(gtot, na.rm = TRUE),
            gtot_sd = sd(gtot, na.rm = TRUE)) -> ss_ess

new_tree %>%
  inner_join(ss_ess, by = "ess") %>%
  mutate(abg = (abg-abg_m)/abg_sd,
         gtot_0 = (gtot_0 - gtot_m) / gtot_sd) %>%
  nest(new_tree = -ess) %>%
  inner_join(m_ess, by = "ess") %>%
  mutate(pred = map2(model, new_tree, ~ predict(.x, .y, type = "response"))) %>%
  unnest(c(new_tree, pred)) %>%
  mutate(surv = ifelse(pred < 0.5, 1, 0)) %>%
  group_by(ess, incref, npp) %>%
  summarise(v_pred = sum(wac * v * ifelse(pred < 0.7, 1, 0))) %>%
  inner_join(nppp, by = "npp") %>%
  group_by(ess, incref) %>%
  summarise(v = sum(v_pred * poids) / sum(poids)) %>%
  mutate(type = "pred") -> pp
  
g3a %>%
  filter(incref > 4, incref < 15, ess %in% target_ess$ess,
         veget5 %in% c("0", "1")) %>%
  group_by(ess, incref, npp) %>%
  summarise(mean_v = sum(w * v)) %>%
  collect() %>%
  inner_join(nppp, by = "npp") %>%
  group_by(ess, incref) %>%
  summarise(v = sum(mean_v * poids) / sum(poids)) %>%
  mutate(type = "obs") -> oo

aa <- rbind(oo, pp)
gp <- ggplot(aa, aes(x=incref+2005, y=v, color=type, group=type)) +
  geom_path() +
  facet_wrap(vars(ess), scales = "free") +
  labs(x = "Year",
       y = "Growing stock of surviving trees")
  
## plot of G, G-5

## dev one species
g3a %>%
  filter(ess == "02", incref > 4,
         incref < 15, simplif == "0") %>%
  select(ess, npp, a, gtot, abg, c13, veget5) %>%
  inner_join(select(g3a5, npp, a, c135), by = c("npp", "a")) %>%
  mutate(gtot_0 = gtot - (abg * 5),
         gtot_2 = pi * (c135 / (2*pi)) ** 2,
         type = case_when(veget5 %in% c("0", "1") ~ "alive",
                          veget5 %in% c("5", "M", "2") ~ "dead",
                          .default = NA)) %>%
  filter(!is.na(type)) %>%
  rename(gtot_1 = gtot) %>%
  collect() -> ee

ee <- filter(ee, gtot_0 > quantile(ee$gtot_0, probs = 0.05),
             gtot_0 < quantile(ee$gtot_0, probs = 0.95))

ggplot(ee, aes(x = gtot_0, xend = gtot_1, y = gtot_1, yend = gtot_2,
               color = type)) +
  geom_segment(arrow = arrow(length = unit(0.1, "inches"))) +
  scale_x_log10() +
  scale_y_log10() 

## tryout new index
g3a %>%
  filter(incref > 4,
         incref < 15, simplif == "0",
         gtot > 0) %>%
  select(ess, npp, a, gtot, abg, c13, veget5) %>%
  inner_join(select(g3a5, npp, a, c135), by = c("npp", "a")) %>%
  mutate(gtot_0 = gtot - (abg * 5),
         gtot_2 = pi * (c135 / (2*pi)) ** 2,
         abg2 = 0.2 * ((c135 ** 2 - c13 ** 2) / (4*pi)),
         type = case_when(veget5 %in% c("0", "1") ~ "alive",
                          veget5 %in% c("5", "M", "2") ~ "dead",
                          .default = NA)) %>%
  filter(gtot_0 > 0) %>%
  group_by(npp, type) %>%
  summarise(abg_1 = mean(abg, na.rm = TRUE),
            abg_2 = mean(abg2, na.rm = TRUE),
            abgg_1 = mean(abg / gtot_0, na.rm = TRUE),
            abgg_2 = mean(abg2 / gtot, na.rm = TRUE)) %>%
  filter(!is.na(type)) %>%
  pivot_wider(names_from = type, values_from = c(abg_1, abg_2, abgg_1, abgg_2), values_fill = 0) %>%
  collect() -> rr

rr %>%
  filter(abg_1_alive > 0, abg_2_alive > 0) %>%
  summarise(r1 = abg_1_dead / abg_1_alive,
            r2 = abg_2_dead / abg_2_alive,
            rr1 = abgg_1_dead / abgg_1_alive,
            rr2 = abgg_2_dead / abgg_2_alive) -> rrr

rrr %>%
  pivot_longer(r1:rr2) %>%
  group_by(name) %>%
  summarise(prop = sum(value > 0 & value < 1, na.rm = TRUE) / sum(value > 0, na.rm = TRUE)) %>%
  mutate(type = c("BAI", "BAI", "BAI/BA", "BAI/BA"),
         time = c("n-10", "n-5", "n-10", "n-5")) %>%
  ggplot(aes(x=time, y=prop, fill=type)) +
  geom_bar(stat = "identity", position = position_dodge(),
           width=0.25) +
  labs(x = "Time of measurement",
       y = "Proportion of plots where growth dead < growth alive") -> g1

ggsave("~/LIF/Mortalite/figures/ratio_growth.png", g1)

g3f %>%
  inner_join(e2p, by = "npp") %>%
  select(npp, pro_psg, esspre_fr, sver) %>%
  filter(npp %in% local(rrr$npp)) %>%
  collect() %>%
  inner_join(rrr, by = "npp") %>%
  pivot_longer(r1:rr2) -> rf

lab <- data.frame(name = c("r1", "r2", "rr1", "rr2"),
                  type = c("BAI", "BAI", "BAI/BA", "BAI/BA"),
                  time = c("n-10", "n-5", "n-10", "n-5"))

rf %>%
  group_by(name, pro_psg) %>%
  summarise(prop = sum(value > 0 & value < 1, na.rm = TRUE) / sum(value > 0, na.rm = TRUE)) %>%
  inner_join(lab, by = "name") %>%
  inner_join(ListDataMod("pro_psg"), by = c("pro_psg" = "mode")) %>%
  ggplot(aes(x=time, y=prop, fill=type)) +
  geom_bar(stat = "identity", position = position_dodge(),
           width=0.25) +
  coord_flip() +
  facet_grid(vars(libelle), labeller = labeller(libelle = label_wrap_gen(20))) +
  labs(x = "Time of measurement",
       y = "Proportion of plots where growth dead < growth alive") -> g1

rf %>%
  group_by(name, esspre_fr) %>%
  summarise(prop = sum(value > 0 & value < 1, na.rm = TRUE) / sum(value > 0, na.rm = TRUE)) %>%
  inner_join(lab, by = "name") %>%
  ggplot(aes(x=time, y=prop, fill=type)) +
  geom_bar(stat = "identity", position = position_dodge(),
           width=0.25) +
  facet_grid(vars(esspre_fr)) +
  labs(x = "Time of measurement",
       y = "Proportion of plots where growth dead < growth alive") -> g1

rf %>%
  filter(sver %in% c("3", "4", "5", "6")) %>%
  group_by(name, sver) %>%
  summarise(prop = sum(value > 0 & value < 1, na.rm = TRUE) / sum(value > 0, na.rm = TRUE)) %>%
  inner_join(lab, by = "name") %>%
  inner_join(ListDataMod("sver"), by = c("sver" = "mode")) %>%
  ggplot(aes(x=time, y=prop, fill=libelle)) +
  geom_bar(stat = "identity", position = position_dodge(),
           width=0.25) +
  facet_wrap(vars(type)) +
  labs(x = "Time of measurement",
       y = "Proportion of plots where growth dead < growth alive") -> g1


## dev look at model missclassification
dd %>%
  filter(!is.na(gtot_0)) %>%
  filter(ess %in% target_ess$ess) %>%
  group_by(ess) %>%
  mutate(mort = ifelse(type == "dead", 1, 0),
         w = case_when(mort == 1 ~ 1 / sum(mort == 1),
                       mort == 0 ~ 1 / sum(mort == 0))) %>%
  select(npp, ess, greco, incref, a, gtot, gtot_0, abg, abg2, lib, mort, w) %>%
  nest(data = -ess) %>%
  mutate(m = map(data, ~glm(mort ~ gtot_0 + abg, data = .x,
                            family = "binomial", weights = w)),
         p = map(m, ~predict(.x, type = "response")),
         pred_class = map(p, ~rbinom(length(.x), size = 1, prob = .x))) %>%
  unnest(c(data,pred_class)) %>%
  mutate(type = case_when(mort == 0 & pred_class == 0 ~ "TN",
                          mort == 0 & pred_class == 1 ~ "FP",
                          mort == 1 & pred_class == 0 ~ "FN",
                          mort == 1 & pred_class == 1 ~ "TP")) -> pp


count(pp, type, ess, greco)
ggplot(pp, aes(x=greco, fill=type)) +
  geom_bar(position = "dodge") +
  facet_wrap(vars(type), scales="free")

pp %>%
  count(greco, type) %>%
  pivot_wider(names_from = type, 
              values_from = n, 
              values_fill = 0) %>%
  mutate(acc = (TN+TP)/(TN+TP+FP+FN),
            sen = TN / (TN + FP),
            spe = TP / (TP + FN))

count(pp, incref, type) %>%
  pivot_wider(names_from = type, 
              values_from = n, 
              values_fill = 0) %>%
  mutate(acc = (TN+TP)/(TN+TP+FP+FN),
         sen = TN / (TN + FP),
         spe = TP / (TP + FN)) %>%
  pivot_longer(acc:spe) %>%
  ggplot(aes(x=incref, y=value, color=name)) +
  geom_path()


pp %>% ungroup() %>% count(greco, type) %>%
  pivot_wider(names_from = type, 
              values_from = n,
              values_fill = 0) %>%
  mutate(accuracy = (TN+TP) / (TN+TP+FP+FN),
         sensitivity = TP / (TP + FN),
         specificity = TN / (TN + FP)) %>%
  pivot_longer(accuracy:specificity) %>%
  inner_join(ListDataMod("greco"), by = c("greco"="mode")) %>%
  ggplot(aes(x=fct_reorder(libelle, value), y=value, color = name)) +
  geom_point(position = position_dodge2(width = 0.25), size = 1.5) +
  coord_flip() +
  theme(text = element_text(size = 20)) +
  labs(x = "",
       y = "Model fitness") -> gg

ggsave("~/LIF/Mortalite/figures/fit_greco.png", gg)

count(pp, type)%>% 
  pivot_wider(names_from = type, 
                           values_from = n,
                           values_fill = 0) %>%
  mutate(accuracy = (TN+TP) / (TN+TP+FP+FN),
         sensitivity = TP / (TP + FN),
         specificity = TN / (TN + FP)) %>%
  pivot_longer(accuracy:specificity) %>%
  inner_join(ListDataMod("ess"), by = c("ess"="mode")) %>%
  ggplot(aes(x=fct_reorder(libelle, value), y=value, color = name)) +
  geom_point(position = position_dodge2(width = 0.25), size = 1.5) +
  theme(text = element_text(size = 20)) +
  coord_flip() +
  labs(x = "",
       y = "Model fitness") -> ge

ggsave("~/LIF/Mortalite/figures/fit_ess.png", ge)

## look at temporal changes
pp %>% ungroup() %>% count(incref, type) %>%
  pivot_wider(names_from = type, 
              values_from = n,
              values_fill = 0) %>%
  mutate(accuracy = (TN+TP) / (TN+TP+FP+FN),
         sensitivity = TP / (TP + FN),
         specificity = TN / (TN + FP)) %>%
  pivot_longer(accuracy:specificity) %>%
  ggplot(aes(x=incref+2005, y=value, color=name)) +
  geom_path() +
  scale_x_continuous(breaks = seq(2010, 2020, by = 2))

## now for some other plot characteristics
dd %>%
  select(npp, a, pro_psg, sver) %>%
  inner_join(pp, by = c("npp", "a")) %>%
  filter(!is.na(sver), !(sver %in% c("0", "X"))) %>%
  count(sver, type) %>%
  pivot_wider(names_from = type, 
                                values_from = n,
                                values_fill = 0) %>%
  mutate(accuracy = (TN+TP) / (TN+TP+FP+FN),
         sensitivity = TP / (TP + FN),
         specificity = TN / (TN + FP)) %>%
  pivot_longer(accuracy:specificity) %>%
  inner_join(ListDataMod("sver"), by = c("sver" = "mode")) %>%
  ggplot(aes(x=fct_reorder(libelle, value), y=value, color = name)) +
  geom_point(position = position_dodge2(width = 0.25), size = 1.5) +
  coord_flip() +
  theme(text = element_text(size = 20)) +
  labs(x = "",
       y = "Model fitness")

dd %>%
  select(npp, a, pro_psg, sver) %>%
  inner_join(pp, by = c("npp", "a")) %>%
  filter(!is.na(pro_psg)) %>%
  count(pro_psg, type) %>%
  pivot_wider(names_from = type, 
              values_from = n,
              values_fill = 0) %>%
  mutate(accuracy = (TN+TP) / (TN+TP+FP+FN),
         sensitivity = TP / (TP + FN),
         specificity = TN / (TN + FP)) %>%
  pivot_longer(accuracy:specificity) %>%
  inner_join(ListDataMod("pro_psg"), by = c("pro_psg" = "mode")) %>%
  ggplot(aes(x=fct_reorder(libelle, value), y=value, color = name)) +
  geom_point(position = position_dodge2(width = 0.25), size = 1.5) +
  coord_flip() +
  theme(text = element_text(size = 20)) +
  labs(x = "",
       y = "Model fitness")

pp %>%
  inner_join(subset(mm, lag == 1), by = "npp") %>%
  ggplot(aes(x=tmean_winter, y=psum_winter, color=type)) +
  geom_point()


## some graphs of mortality dynamics
yy <- data.table(year = 2008:2024)

## first mortality of first visit
m1 <- yy[,ocre_nm(ans=year, vent_a = c("RAKVEGET:RAKVEGET:1", "datemort", "veget"),
                  vars = "V", vifs = FALSE, morts = TRUE, moy_ha = TRUE),
         by = year]

## then mortality of second visit
m2 <- yy[3:13,ocre_nm(ans=year, vars = "vmort_an", 
                      moy_ha = TRUE), by = year]

## some cleaning
m1 %>%
  filter(datemort.datemort %in% c("1", "2"),
         veget.dom8m1 %in% c("5")) %>%
  select(year, datemort.datemort, moy_ha, ic_moy_ha) %>%
  mutate(datemort.datemort = ifelse(datemort.datemort == "1",
                                    "1st visit -5year", "1st visit +5year")) -> m1d

m2 %>%
  filter(veget5.veget56 %in% c("M")) %>%
  group_by(year) %>%
  summarise(moy_ha = sum(moy_ha),
            ic_moy_ha = sqrt(sum(ic_moy_ha ** 2))) %>%
  mutate(datemort.datemort = "2nd visit -5year",
         year = year+5) -> m2d

ma <- rbind(m1d, m2d[,names(m1d)])

gt <- ggplot(ma, aes(x=year, y=moy_ha, ymin=moy_ha-ic_moy_ha, ymax = moy_ha+ic_moy_ha)) +
  geom_ribbon(aes(fill = datemort.datemort), alpha = 0.1) +
  geom_path(aes(color= datemort.datemort)) +
  scale_color_discrete(name = "type") +
  scale_fill_discrete(name = "type") +
  labs(y = "Volume of dead trees (m^3/ha)")

ggsave("LIF/Mortalite/figures/mort_trend.png", gt)
