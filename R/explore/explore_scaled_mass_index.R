# see https://doi.org/10.1016/j.anbehav.2012.12.002 for a similar analysis that
# tried to use the SMI but ended up simply using the raw body mass data instead

# take the first observation from each bird in the dataset 
# (i.e., so that each data point is independent for the SMA regression)
cap_05_09_std_pca_1 <- 
  cap_05_09_std_pca %>% 
  group_by(capture_id) %>% 
  arrange(capture_id, gpsdt_std) %>% 
  slice(1) %>% 
  ungroup() %>% 
  na.omit() %>% 
  dplyr::select(capture_id, gpsdt_std, weight, culmen, totalHead, tarsus, wing, structure_pc1)

cor.test(cap_05_09_std_pca_1$weight, cap_05_09_std_pca_1$wing)
cor.test(cap_05_09_std_pca_1$weight, cap_05_09_std_pca_1$tarsus)

# estimate the slope of the log-weight and log-wing SMA regression
weight_wing_sma_slope <- coef(sma(log(cap_05_09_std_pca_1$weight) ~ log(cap_05_09_std_pca_1$wing), method = "MA"))[2]
weight_tarsus_sma_slope <- coef(sma(log(cap_05_09_std_pca_1$weight) ~ log(cap_05_09_std_pca_1$tarsus)))[2]

# calculate the average wing size
avg_wing <- mean(cap_05_09_std_pca_1$wing)
avg_tarsus <- mean(cap_05_09_std_pca_1$tarsus)

# calculate scaled mass index for full dataset
cap_05_09_std_pca_smi <-
  cap_05_09_std_pca %>% 
  mutate(wing_10 = wing/10) %>%
  ungroup() %>% 
  mutate(smi_wing = weight * (avg_wing/wing)^weight_wing_sma_slope,
         smi_tarsus = weight * (avg_tarsus/tarsus)^weight_tarsus_sma_slope)

# check if smi and wing length are not correlated
test <- 
  cap_05_09_std_pca_smi %>% 
  group_by(capture_id) %>% 
  arrange(capture_id, gpsdt_std) %>% 
  slice(1) %>% 
  ungroup() %>% 
  na.omit()

cor.test(test$weight, test$smi_wing)
cor.test(test$tarsus, test$smi_tarsus)

mod_smi <-
  lmer(smi_wing ~ date_deviance + first_date + last_date +
         (1 | capture_id),
       data = cap_05_09_std_pca_smi)

# Derive confidence intervals of effect sizes from parametric bootstrapping
tidy_mod_smi <-
  tidy(mod_smi, conf.int = TRUE, conf.method = "boot", nsim = 1000)

# run rptR to obtain repeatabilities of random effects
rpt_mod_smi <-
  rpt(smi_wing ~ date_deviance + first_date + last_date +
        (1 | capture_id),
      grname = c("capture_id", "Fixed"),
      data = cap_05_09_std_pca_smi,
      datatype = "Gaussian",
      nboot = 1000, npermut = 1000, ratio = TRUE,
      adjusted = TRUE, ncores = 4, parallel = TRUE)

# run partR2 on each model to obtain marginal R2, parameter estimates, and beta
# weights
R2m_mod_smi <-
  partR2(mod_smi,
         partvars = c("date_deviance",
                      "first_date",
                      "last_date"),
         R2_type = "marginal",
         nboot = 1000,
         CI = 0.95,
         max_level = 1)

R2c_mod_smi <-
  partR2(mod_smi,
         partvars = c("date_deviance",
                      "first_date",
                      "last_date"),
         R2_type = "conditional",
         nboot = 1000,
         CI = 0.95,
         max_level = 1)

# males typically decrease their weight over time
summary(glht(mod_smi))
plot(allEffects(mod_smi))

#### Table of effect sizes (van de Pol method) ----
# Retrieve sample sizes
sample_sizes <-
  cap_05_09_std_pca %>% 
  ungroup() %>% 
  summarise(Year = n_distinct(year_),
            Individual = n_distinct(capture_id),
            Observations = nrow(.))

sample_sizes <- 
  as.data.frame(t(as.data.frame(sample_sizes))) %>%
  rownames_to_column("term") %>% 
  rename(estimate = V1) %>% 
  mutate(stat = "n")

# # dataset summary
# cap_05_09_std_pca %>% 
#   ungroup() %>% 
#   summarise(max_weight = max(weight, na.rm = TRUE),
#             min_weight = min(weight, na.rm = TRUE),
#             mean_weight = mean(weight, na.rm = TRUE),
#             sd_weight = sd(weight, na.rm = TRUE)) %>% 
#   t()


# clean model component names
mod_comp_names <- 
  data.frame(comp_name = c("Within ind. temporal change",
                           "Between ind. effect of season (first measure)",
                           "Between ind. effect of season (last measure)",
                           "Year 2006",
                           "Year 2007",
                           "Year 2008",
                           "Year 2009",
                           "Total Marginal \U1D479\U00B2",
                           "Within ind. temporal change",
                           "Between ind. effect of season (first measure)",
                           "Between ind. effect of season (last measure)",
                           "Year",
                           "Total Conditional \U1D479\U00B2",
                           "Individual",
                           "Residual",
                           "Individual",
                           "Residual",
                           "Years",
                           "Individuals",
                           "Observations"))

# Fixed effect sizes (non-standardized)
fixefTable <- 
  stats_mod_weight$tidy %>% 
  dplyr::filter(effect == "fixed") %>% 
  dplyr::select(term, estimate, conf.low, conf.high) %>% 
  as.data.frame() %>% 
  mutate(stat = "fixed")

# Fixed effect sizes (standardized)
fixef_bw_Table <- 
  stats_mod_weight$partR2m$BW %>% 
  as.data.frame() %>% 
  mutate(stat = "fixed_bw") %>% 
  rename(conf.low = CI_lower,
         conf.high = CI_upper)

# Semi-partial R2 estimates
R2Table <- 
  bind_rows(stats_mod_weight$partR2m$R2,
            stats_mod_weight$partR2c$R2[1,]) %>% 
  dplyr::select(term, estimate, CI_lower, CI_upper) %>% 
  as.data.frame() %>% 
  mutate(stat = "partR2") %>% 
  rename(conf.low = CI_lower,
         conf.high = CI_upper)

# Random effects variances
ranefTable <- 
  stats_mod_weight$tidy %>% 
  dplyr::filter(effect == "ran_pars") %>% 
  dplyr::select(group, estimate, conf.low, conf.high) %>% 
  as.data.frame() %>% 
  mutate(stat = "rand") %>% 
  rename(term = group) %>% 
  mutate(estimate = estimate^2,
         conf.high = conf.high^2,
         conf.low = conf.low^2)

# Adjusted repeatabilities
coefRptTable <- 
  stats_mod_weight$rptR$R_boot %>% 
  dplyr::select(-Fixed) %>% 
  mutate(residual = 1 - rowSums(.)) %>% 
  apply(., 2, 
        function(x) c(mean (x), quantile (x, prob = c(0.025, 0.975)))) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("term") %>% 
  rename(estimate = V1,
         conf.low = `2.5%`,
         conf.high = `97.5%`) %>% 
  mutate(stat = "RptR")

# Store all parameters into a single table and clean it up
allCoefs_mod <- 
  bind_rows(fixef_bw_Table,
            R2Table,
            ranefTable, 
            coefRptTable, 
            sample_sizes) %>% 
  bind_cols(.,
            mod_comp_names) %>%
  mutate(coefString = ifelse(!is.na(conf.low),
                             paste0("[", 
                                    round(conf.low, 2), ", ", 
                                    round(conf.high, 2), "]"),
                             NA),
         effect = c(rep("Fixed effects \U1D6FD (standardized)", nrow(fixef_bw_Table)),
                    rep("Partitioned \U1D479\U00B2", nrow(R2Table)),
                    rep("Random effects \U1D70E\U00B2", nrow(ranefTable)),
                    rep("Adjusted repeatability \U1D45F", nrow(coefRptTable)),
                    rep("Sample sizes \U1D45B", nrow(sample_sizes)))) %>%
  dplyr::select(effect, everything())

# # re-organize model components for table
# allCoefs_mod <-
#   allCoefs_mod[c(1:4, 8, 6:9, 16, 15, 10:14, 17:28), ]
# 
# allCoefs_mod %>% 
#   round(estimate)

# draw gt table
mod_weight_wing_table <- 
  allCoefs_mod %>% 
  dplyr::select(effect, comp_name, estimate, coefString) %>% 
  gt(rowname_col = "row",
     groupname_col = "effect") %>% 
  cols_label(comp_name = html("<i>male Pectoral Sandpiper body mass dynamics</i>"),
             estimate = "Mean estimate",
             coefString = "95% confidence interval") %>% 
  fmt_number(columns = c(estimate),
             rows = 1:19,
             decimals = 2,
             use_seps = FALSE) %>% 
  fmt_number(columns = c(estimate),
             rows = 20:22,
             decimals = 0,
             use_seps = FALSE) %>% 
  sub_missing(columns = 1:4,
              missing_text = "") %>% 
  cols_align(align = "left",
             columns = c(comp_name)) %>% 
  tab_options(row_group.font.weight = "bold",
              row_group.background.color = brewer.pal(9,"Greys")[3],
              table.font.size = 12,
              data_row.padding = 3,
              row_group.padding = 4,
              summary_row.padding = 2,
              column_labels.font.size = 14,
              row_group.font.size = 12,
              table.width = pct(60))

mod_weight_wing_table

#     Peig & Green
#    "New perspectives for estimating body condition from mass/length data:
#     the scaled mass index as an alternative method"
#     Oikos 118: 1883-1891, 2009

cap_05_09_body_1 <- 
  cap_05_09_std_pca %>% 
  group_by(capture_id) %>% 
  arrange(capture_id, gpsdt_std) %>% 
  slice(1) %>% 
  ungroup() %>% 
  na.omit() %>% 
  dplyr::select(capture_id, gpsdt_std, weight, culmen, totalHead, tarsus,  wing, structure_pc1)

b.msa.ols_1 <- coef(sma(log(cap_05_09_body_1$weight) ~ log(cap_05_09_body_1$wing) ))[2]
SMI.ols_1 <- cap_05_09_body_1$weight * (mean(cap_05_09_body_1$wing) / cap_05_09_body_1$wing) ^ b.msa.ols_1
cap_05_09_body_1 <- 
  bind_cols(cap_05_09_body_1, SMI.ols_1) %>% 
  rename(SMI_ols_1 = `...9`)

cap_05_09_body_all <- 
  cap_05_09_std_pca %>% 
  group_by(capture_id) %>% 
  arrange(capture_id, gpsdt_std) %>% 
  # slice(1) %>% 
  ungroup() %>% 
  na.omit() %>% 
  dplyr::select(capture_id, gpsdt_std, weight, culmen, totalHead, tarsus,  wing, structure_pc1)

b.msa.ols_all <- coef(sma(log(cap_05_09_body_all$weight) ~ log(cap_05_09_body_all$wing)))[2]
SMI.ols_all <- cap_05_09_body_all$weight * (mean(cap_05_09_body_all$wing) / cap_05_09_body_all$wing) ^ b.msa.ols_all
cap_05_09_body_all <- 
  bind_cols(cap_05_09_body_all, SMI.ols_all) %>% 
  rename(SMI_ols_all = `...9`) %>% 
  mutate(SMI_ols_1 = weight * (143.7568/wing)^b.msa.ols_1)

cap_05_09_body_all %>% 
  group_by(capture_id) %>% 
  arrange(capture_id, gpsdt_std) %>% 
  slice(c(1, 2)) %>% 
  group_by(capture_id) %>% 
  mutate(rank = rank(gpsdt_std)) %>% 
  pivot_wider(id_cols = capture_id, names_from = rank, values_from = c(SMI_ols_1, SMI_ols_all)) %>% 
  mutate(SMI_ols_1_diff = SMI_ols_1_1 - SMI_ols_1_2,
         SMI_ols_all_diff = SMI_ols_all_1 - SMI_ols_all_2) %>% 
  ggplot() +
  geom_point(aes(x = SMI_ols_1_diff, y = SMI_ols_all_diff))
  
logM.ols <- lm(log(cap_05_09_body$weight) ~ log(cap_05_09_body$wing))
logM.rob <- rlm(log(cap_05_09_body$weight) ~ log(cap_05_09_body$wing), method = "M")
b.msa.ols <- coef(sma(log(cap_05_09_body$weight) ~ log(cap_05_09_body$wing)))[2]
b.msa.rob <- coef(sma(log(cap_05_09_body$weight) ~ log(cap_05_09_body$wing), robust = T))[2]
SMI.ols <- cap_05_09_body$weight * (mean(cap_05_09_body$wing) / cap_05_09_body$wing) ^ b.msa.ols
SMI.rob <- cap_05_09_body$weight * (mean(cap_05_09_body$wing) / cap_05_09_body$wing) ^ b.msa.rob
res <- data.frame(SMI.ols, SMI.rob, cap_05_09_body$wing, cap_05_09_body$weight)

library(dplyr)

pred.DT <- 
  data.frame(x = seq(min(cap_05_09_body$wing), max(cap_05_09_body$wing), length = 100)) %>%
  mutate(y.ols = predict(logM.ols, newdata = .) %>% exp(),
         y.rob = predict(logM.rob, newdata = .) %>% exp())


pred.DT <-
  data.table(x = seq(min(x), max(x), length = 100)) %>%
  .[, y.ols := predict(logM.ols, newdata = .) %>% exp] %>%
  .[, y.rob := predict(logM.rob, newdata = .) %>% exp]

scaledMassIndex <-
  function(x, y, x.0 = mean(x)) {
    require(smatr)
    require(magrittr)
    require(MASS)
    require(data.table)
    logM.ols <- lm(log(y) ~ log(x))
    logM.rob <- rlm(log(y) ~ log(x), method = "M")
    b.msa.ols <- coef(sma(log(y) ~ log(x)))[2]
    b.msa.rob <- coef(sma(log(y) ~ log(x), robust = T))[2]
    SMI.ols <- y * (x.0 / x) ^ b.msa.ols
    SMI.rob <- y * (x.0 / x) ^ b.msa.rob
    res <- data.frame(SMI.ols, SMI.rob, x, y)
    pred.DT <-
      data.table(x = seq(min(x), max(x), length = 100)) %>%
      .[, y.ols := predict(logM.ols, newdata = .) %>% exp] %>%
      .[, y.rob := predict(logM.rob, newdata = .) %>% exp]
    attr(res, "b.msa") <- c(ols = b.msa.ols, rob = b.msa.rob)
    return(res)
  }

library(magrittr)
library(data.table)
library(ggplot2)
library(ggpubr)
set.seed(123)
dt <-
  data.table(`Body length` = seq(1, 5, 0.5)) %>%
  .[, `Body weight` := 2 * `Body length` ^ (3 + rnorm(length(`Body length`), sd = 0.001))] %>%
  rbind(., data.table(`Body length` = 6, `Body weight` = 200)) %>%
  .[order(`Body length`), ]

x = dt$`Body length`
y = dt$`Body weight`
x.0 = mean(dt$`Body length`)

dt.SMI <-
  scaledMassIndex(dt$`Body length`, dt$`Body weight`, x.0 = 2)
summary(dt.SMI)

cor.test(dt.SMI$SMI.ols, dt.SMI$x)

g1 <-
  ggplot(dt, aes(`Body length`, `Body weight`)) +
  geom_point() +
  geom_line(aes(x, y, color = Method),
            data = attr(dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")) +
  scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
g2 <-
  ggplot(dt, aes(log10(`Body length`), log10(`Body weight`))) +
  geom_point() +
  geom_line(
    aes(log10(x), log10(y), color = Method),
    data = attr(dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")
  ) +
  scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
g3 <-
  ggplot(melt(dt.SMI, "x", c("SMI.ols", "SMI.rob"), "Method"),
         aes(x, value)) +
  geom_point(aes(color = Method)) +
  xlab("Body length") +
  ylab("SMI (at body length = 2)") +
  scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
ggarrange(g1,
          g2,
          g3,
          ncol = 1,
          nrow = 3,
          align = "hv")