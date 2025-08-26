## data analysis

library(tidyverse)
library(broom)

## load data
dat <- read.csv("~/LIF/Mortalite/data/maindat.csv")
sp_sel <- read.csv("~/LIF/Mortalite/data/sp_sel.csv")

## only keep selected species
dat <- filter(dat, ess %in% sp_sel$ess)

## helper function to fit a logistic model per species and indice
## and get accuracy metrics
foo_log <- function(data, indice){
  ## create the response variable
  data$mort <- ifelse(data$type == "dead", 1, 0)
  ## subset to the variable of interest
  Xvar <- strsplit(indice, "\\+")[[1]]
  sub_dat <- select(data, all_of(Xvar), mort)
  ## remove the NAs
  sub_dat <- na.omit(sub_dat)
  ## std the predictors
  sub_dat <- mutate(sub_dat, across(Xvar, scale))
  ## train and test partition
  set.seed(123)
  train_index <- tryCatch(caret::createDataPartition(sub_dat$mort, p = 0.8, list = FALSE),
                          error = function(e) sample(1:nrow(sub_dat),
                                                     ceiling(0.8*nrow(sub_dat))))
  train <- sub_dat[train_index, ]
  test <- sub_dat[-train_index, ]
  ## make weights
  class_weights <- table(train$mort)
  inv_weights <- ifelse(train$mort == 1,
                        1 / class_weights[2],  # minority (dead)
                        1 / class_weights[1])
  ## create formula
  ff <- as.formula(paste0("mort~", indice))
  ## fit model
  m <- glm(ff, data = train, weights = inv_weights,
           family = "binomial")
  ## make prediction
  pred <- rbinom(nrow(test), 1, predict(m, newdata = test, type = "response"))
  
  ## derive accuracy metrics
  cml <- table(Pred=pred, True = test$mort)
  # Defensive code in case a class is missing
  TP <- tryCatch(cml["1", "1"], error = function(e) 0)
  TN <- tryCatch(cml["0", "0"], error = function(e) 0)
  FP <- tryCatch(cml["1", "0"], error = function(e) 0)
  FN <- tryCatch(cml["0", "1"], error = function(e) 0)
  
  out <- data.frame(accuracy = (TP + TN) / sum(cml),
                    sensitivity = TP / (TP + FN),
                    specificity = TN / (TN + FP))
  
  return(out)
}

## make model across indicators
idl <- c("abg+gtot_0", "abg2+gtot", "abgg", "abgg2", "abg_diff+gtot",
         "abg_pl+gtot_0", "abg2_pl+gtot", "ir_trend+gtot", "ir_cv+gtot")

ex_dat <- dat[rep(1:nrow(dat), each = 9),]
ex_dat$indice <- rep(idl, times = nrow(dat))

ex_dat %>%
  nest(data=-c(ess, indice)) %>%
  mutate(res = map2(data, indice, ~foo_log(.x, .y))) %>%
  unnest(res) -> res_log

res_log %>%
  select(-data) %>%
  pivot_longer(accuracy:specificity) -> dd

g1 <- ggplot(dd, aes(x=fct_reorder(indice, value), y=value)) +
  stat_summary(fun.data = "mean_cl_boot", color="red") +
  geom_jitter(width = 0.25, alpha=0.25) +
  facet_wrap(vars(name), scales = "free_x") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  coord_flip() +
  labs(x = "Predictors",
       y = "Model precision")

ggsave("~/LIF/Mortalite/figures/model_res_9ind.png", g1)

## get top ind per species and the associated values
res_log %>%
  group_by(ess) %>%
  slice_max(accuracy, n = 1) %>%
  arrange(desc(accuracy)) %>%
  select(-data) -> tt
  
write.csv(tt, "~/LIF/Mortalite/table/top_ind.csv", row.names=FALSE)


## helper functions to get regression curves
foo_reg <- function(data, indice){
  ## create the response variable
  data$mort <- ifelse(data$type == "dead", 1, 0)
  ## subset to the variable of interest
  Xvar <- strsplit(indice, "\\+")[[1]]
  sub_dat <- select(data, Xvar, mort)
  ## remove the NAs
  sub_dat <- na.omit(sub_dat)
  ## std the predictors
  sub_dat <- mutate(sub_dat, across(Xvar, ~scale(.x)[,1]))
  ## make weights
  class_weights <- table(sub_dat$mort)
  inv_weights <- ifelse(sub_dat$mort == 1,
                        1 / class_weights[2],  # minority (dead)
                        1 / class_weights[1])
  ## create formula
  ff <- as.formula(paste0("mort~", indice))
  ## fit model
  m <- glm(ff, data = sub_dat, weights = inv_weights,
           family = "binomial")
  #print(coef(m))
  
  ## make predictions
  newdat <- data.frame(gtot = 0,
                       abg = seq(quantile(sub_dat[,1, drop=TRUE], probs = 0.05),
                                 quantile(sub_dat[,1, drop=TRUE], probs = 0.95),
                                 length.out = 50))
  names(newdat) <- Xvar
  preds <- predict(m, newdata = newdat, type="link",se.fit=TRUE)
  out <- data.frame(lwr = m$family$linkinv(preds$fit - 1.64 * preds$se.fit),
                    upr = m$family$linkinv(preds$fit + 1.64 * preds$se.fit),
                    fit = m$family$linkinv(preds$fit),
                    abg = newdat[,2])
  
  return(out)
}

## make regression curves from the abg*+gtot model
ex_dat %>%
  filter(indice %in% c("abg+gtot_0",
                       "abg2+gtot")) %>%
  nest(data = -c(ess, indice)) %>%
  mutate(rr = map2(data, indice, ~foo_reg(.x, .y))) %>%
  unnest(rr) -> res_reg

ggplot(res_reg, aes(x=abg, y=fit, color=indice)) +
  geom_path() +
  facet_wrap(vars(ess))
