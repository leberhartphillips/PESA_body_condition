source("R/DAT/wrangle_datasets.R")
source("R/project/project_plotting.R")
cap_all_phenology <-
  cap_all %>% 
  mutate(capture_id = paste(ID, year_, sep = "_")) %>% 
  filter(!is.na(gpsdt)) %>% 
  mutate(gpsdt = as.POSIXct(gpsdt, origin = "1970-01-01")) %>% 
  mutate(gpsdt_year_start = ymd_hms(paste0(year(gpsdt), "-01-01 00:00:00"))) %>% 
  mutate(gpsdt_since_year_start = as.numeric(gpsdt - gpsdt_year_start)) %>% 
  arrange(ID, gpsdt) %>% 
  group_by(ID) %>% 
  slice(1) %>% 
  arrange(year_) %>% 
  group_by(year_) %>% 
  summarize(mean_date = mean(gpsdt_since_year_start),
            sd_date = sd(gpsdt_since_year_start))

# assess the distribution of repeated measures
cap_05_09 %>% 
  dplyr::select(year_, ID, n_by_id) %>% 
  distinct() %>% 
  count(n_by_id)

# assess variation in rate of body mass change between first and second captures
cap_05_09_delta %>% 
  dplyr::select(year_, capture_id, gpsdt, delta_weight, delta_time) %>% 
  na.omit() %>% 
  mutate(mass_change = delta_weight/as.numeric(delta_time)) %>% 
  arrange(desc(mass_change)) %>% 
  group_by(capture_id) %>% 
  slice(1) %>% 
  ggplot() +
  geom_histogram(aes(mass_change))

cap_05_09_delta %>% 
  dplyr::select(year_, capture_id, gpsdt, delta_weight, delta_time) %>% 
  na.omit() %>% 
  mutate(mass_change = delta_weight/as.numeric(delta_time)) %>% 
  mutate(delta_time = as.numeric(delta_time)/60/60/24) %>% 
  arrange(desc(mass_change)) %>% 
  group_by(capture_id) %>% 
  slice(1) %>% 
  sumtable(vars = c("mass_change", "delta_time", "delta_weight"))

# weird outlier with a change of 7.5g in 2.5 hours!! (maybe incorrect?)
cap_all %>% 
  filter(ID == "162116388") %>% 
  mutate(gpsdt = as.POSIXct(gpsdt, origin = "1970-01-01"))

# standardize capture date according to the year
cap_05_09_std <- 
  cap_05_09 %>% 
  filter(n_by_id > 1, !is.na(weight)) %>% 
  mutate(capture_id = paste(ID, year_, sep = "_")) %>% 
  dplyr::select(-capdt, -ID) %>% 
  mutate(gpsdt_year_start = ymd_hms(paste0(year(gpsdt), "-01-01 00:00:00"))) %>% 
  mutate(gpsdt_since_year_start = as.numeric(gpsdt - gpsdt_year_start)) %>% 
  left_join(cap_all_phenology, by = "year_") %>%
  mutate(gpsdt_std = gpsdt_since_year_start - mean_date) %>% 
  group_by(capture_id) %>% 
  arrange(capture_id, gpsdt_since_year_start) %>% 
  mutate(measure_deviation = gpsdt_std - gpsdt_std[which.min(gpsdt_std)],
         first_measure = gpsdt_std[which.min(gpsdt_std)],
         last_measure = gpsdt_std[which.max(gpsdt_std)]) %>% 
  dplyr::select(year_, capture_id, gpsdt, measure_deviation, first_measure, 
                last_measure, weight, culmen, totalHead, tarsus,  wing, n_by_id)

# PCA of static body structural measurements (culmen, totalHead, tarsus, wing)
static_measures_pca <-
  cap_05_09_std %>% 
  group_by(capture_id) %>% 
  arrange(gpsdt) %>% 
  slice(1) %>% 
  ungroup() %>% 
  na.omit() %>%
  dplyr::select(culmen, totalHead, tarsus,  wing) %>%
  princomp()

# check the PCA results
summary(static_measures_pca)
biplot(static_measures_pca, cex = 0.7)

# bind PC1 to the original dataframe
cap_05_09_std_pca <- 
  cap_05_09_std %>% 
  group_by(capture_id) %>% 
  arrange(capture_id, gpsdt) %>% 
  slice(1) %>% 
  ungroup() %>% 
  na.omit() %>% 
  dplyr::select(capture_id, culmen, totalHead, tarsus,  wing) %>%
  bind_cols(., static_measures_pca$scores[, 1]) %>% 
  rename(structure_pc1 = `...6`) %>% 
  left_join(cap_05_09_std %>% 
              dplyr::select(-c(culmen, totalHead, tarsus,  wing)), ., 
            by = "capture_id", multiple = "all") %>% 
  na.omit()

ggplot(cap_05_09_std_pca, 
       aes(x = structure_pc1, y = weight)) + 
  geom_point() + 
  labs(x = "PC1", y = "Body Mass") + 
  theme_minimal()

ggplot(cap_05_09_std_pca, 
       aes(x = wing, y = weight)) + 
  geom_point() + 
  labs(x = "wing length", y = "Body Mass") + 
  theme_minimal()

# check 
cap_05_09_std_pca %>% 
  arrange(capture_id, gpsdt) %>% 
  slice(1) %>% 
  ungroup() %>% 
  na.omit() %>% 
  dplyr::select(weight, culmen, totalHead, tarsus,  wing, structure_pc1) %>% 
  cor() %>% 
  corrplot(type = "upper", method = "number", tl.srt = 45)

# model
# mod_weight_PC1 <-
#   lmer(weight ~ structure_pc1 + measure_deviation + first_measure +
#          (1 | capture_id),
#        data = cap_05_09_std_pca)

mod_weight_wing <-
  lmer(weight ~ wing + measure_deviation + first_measure + last_measure +
         (1 | capture_id),
       data = cap_05_09_std_pca)

# tidy_mod_weight_PC1 <-
#   tidy(mod_weight, conf.int = TRUE, conf.method = "boot", nsim = 1000)

tidy_mod_weight_wing <-
  tidy(mod_weight_wing, conf.int = TRUE, conf.method = "boot", nsim = 1000)

# run rptR to obtain repeatabilities of random effects
rpt_mod_weight_wing <-
  rpt(weight ~ wing + measure_deviation + first_measure + last_measure +
        (1 | capture_id),
      grname = c("capture_id", "Fixed"),
      data = cap_05_09_std_pca,
      datatype = "Gaussian",
      nboot = 1000, npermut = 1000, ratio = TRUE,
      adjusted = TRUE, ncores = 4, parallel = TRUE)

# run partR2 on each model to obtain marginal R2, parameter estimates, and beta
# weights
R2m_mod_weight_wing <-
  partR2(mod_weight_wing,
         partvars = c("wing",
                      "measure_deviation",
                      "first_measure",
                      "last_measure"),
         R2_type = "marginal",
         nboot = 1000,
         CI = 0.95,
         max_level = 1)

R2c_mod_weight_wing <-
  partR2(mod_weight_wing,
         partvars = c("wing",
                      "measure_deviation",
                      "first_measure",
                      "last_measure"),
         R2_type = "conditional",
         nboot = 1000,
         CI = 0.95,
         max_level = 1)

stats_mod_weight_wing <- 
  list(mod = mod_weight_wing,
       tidy = tidy_mod_weight_wing,
       rptR = rpt_mod_weight_wing,
       partR2m = R2m_mod_weight_wing,
       partR2c = R2c_mod_weight_wing)

save(stats_mod_weight_wing,
     file = "R/output/stats_mod_weight_wing.rds")

# quick visualizations of model

plot(allEffects(mod_weight_wing))

# Get the residuals
resid <- residuals(mod_weight)

# Plot the residuals vs fitted values
ggplot(data.frame(resid = resid, fitted = fitted(mod_weight)), aes(x = fitted, y = resid)) +
  geom_point() +
  labs(x = "Fitted Values", y = "Residuals") +
  ggtitle("Residuals vs Fitted Values")

# Check for normality
qqPlot(resid)

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

# dataset summary
cap_05_09_std_pca %>% 
  ungroup() %>% 
  summarise(max_weight = max(weight, na.rm = TRUE),
            min_weight = min(weight, na.rm = TRUE),
            mean_weight = mean(weight, na.rm = TRUE),
            sd_weight = sd(weight, na.rm = TRUE)) %>% 
  t()


# clean model component names
mod_comp_names <- 
  data.frame(comp_name = c("Wing length",
                           "Within ind. temporal change",
                           "Between ind. effect of season (first measure)",
                           "Between ind. effect of season (last measure)",
                           "Total Marginal \U1D479\U00B2",
                           "Wing length",
                           "Within ind. temporal change",
                           "Between ind. effect of season (first measure)",
                           "Between ind. effect of season (last measure)",
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
  stats_mod_weight_wing$tidy %>% 
  dplyr::filter(effect == "fixed") %>% 
  dplyr::select(term, estimate, conf.low, conf.high) %>% 
  as.data.frame() %>% 
  mutate(stat = "fixed")

# Fixed effect sizes (standardized)
fixef_bw_Table <- 
  stats_mod_weight_wing$partR2m$BW %>% 
  as.data.frame() %>% 
  mutate(stat = "fixed_bw") %>% 
  rename(conf.low = CI_lower,
         conf.high = CI_upper)

# Semi-partial R2 estimates
R2Table <- 
  bind_rows(stats_mod_weight_wing$partR2m$R2,
            stats_mod_weight_wing$partR2c$R2[1,]) %>% 
  dplyr::select(term, estimate, CI_lower, CI_upper) %>% 
  as.data.frame() %>% 
  mutate(stat = "partR2") %>% 
  rename(conf.low = CI_lower,
         conf.high = CI_upper)

# Random effects variances
ranefTable <- 
  stats_mod_weight_wing$tidy %>% 
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
  stats_mod_weight_wing$rptR$R_boot %>% 
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
#   allCoefs_mod[c(5, 1:4, 6:9, 16, 15, 10:14, 17:28), ]
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
             rows = 1:14,
             decimals = 2,
             use_seps = FALSE) %>% 
  fmt_number(columns = c(estimate),
             rows = 15:17,
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

# export table to disk
mod_weight_wing_table %>%
  gtsave("mod_weight_wing_table.rtf", path = "R/products/tables/rtf/")

mod_weight_wing_table %>%
  gtsave("mod_weight_wing_table.png", path = "R/products/tables/png/")

#### Forest plot of results ----
# Standardized fixed effects
mod_weight_wing_forest_plot_fixef <-
  allCoefs_mod %>%
  filter(str_detect(effect, "Fixed") & 
           term != "(Intercept)") %>%
  mutate(comp_name = fct_relevel(comp_name,
                                 "Between ind. effect of season (last measure)",
                                 "Between ind. effect of season (first measure)", 
                                 "Within ind. temporal change",
                                 "Wing length")) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbarh(aes(xmin = conf.low,
                     xmax = conf.high,
                     y = comp_name),
                 alpha = 1, color = col_all, 
                 size = 0.5,
                 height = 0) +
  geom_point(aes(y = comp_name, x = estimate),
             size = 3, shape = 21, 
             fill = "#ECEFF4", col = col_all, 
             alpha = 1, stroke = 0.5) +
  luke_theme +
  theme(axis.title.x = element_text(size = 10, hjust = 0.5),
        plot.title = element_text(face = 'italic', hjust = 0.5)) +
  ylab("Fixed effects") +
  xlab(expression(italic(paste("              Standardized effect size (", beta,")" %+-% "95% CI", sep = "")))) #+
  # ggtitle('male Pectoral Sandpiper body mass dynamics')

# Semi-partial R2 estimates
mod_weight_wing_forest_plot_partR2 <-
  allCoefs_mod %>%
  filter(str_detect(effect, "Partitioned") & str_detect(comp_name, "Conditional", negate = TRUE)) %>%
  mutate(comp_name = fct_relevel(comp_name,
                                 "Between ind. effect of season (last measure)",
                                 "Between ind. effect of season (first measure)", 
                                 "Within ind. temporal change",
                                 "Wing length",
                                 "Total Marginal \U1D479\U00B2")) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbarh(aes(xmin = conf.low,
                     xmax = conf.high,
                     y = comp_name),
                 alpha = 1, color = col_all, 
                 size = 0.5,
                 height = 0) +
  geom_point(aes(y = comp_name, x = estimate),
             size = 3, shape = 21, 
             fill = "#ECEFF4", col = col_all, 
             alpha = 1, stroke = 0.5) +
  luke_theme +
  theme(axis.title.x = element_text(size = 10, hjust = 0.5)) +
  scale_y_discrete(labels = c("Between ind. effect of season (last measure)" = expression("Between ind. effect of season (last measure)"),
                              "Between ind. effect of season (first measure)" = expression("Between ind. effect of season (first measure)"),
                              "Within ind. temporal change" = expression("Within ind. temporal change"),
                              "Wing length" = expression("Wing length"),
                              "Total Marginal \U1D479\U00B2" = expression(paste("Total marginal ", italic("R"), ''^{2}, sep = "")))) +
  ylab(expression(paste("Semi-partial ", italic("R"),''^{2}, sep = ""))) +
  xlab(expression(italic(paste("               Variance explained (R", ''^{2}, ")" %+-% "95% CI", sep = ""))))

# Random effect variances
mod_weight_wing_forest_plot_randef <-
  allCoefs_mod %>%
  filter(str_detect(effect, "Random")) %>%
  mutate(comp_name = fct_relevel(comp_name,
                                 "Residual",
                                 "Individual")) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbarh(aes(xmin = conf.low,
                     xmax = conf.high,
                     y = comp_name),
                 alpha = 1, color = col_all, 
                 size = 0.5,
                 height = 0) +
  geom_point(aes(y = comp_name, x = estimate),
             size = 3, shape = 21, 
             fill = "#ECEFF4", col = col_all, 
             alpha = 1, stroke = 0.5) +
  luke_theme +
  theme(axis.title.x = element_text(size = 10, hjust = 0.5)) +
  ylab("Random\neffects") +
  xlab(expression(italic(paste("Variance (", sigma, ''^{2}, ")" %+-% "95% CI", sep = ""))))

# Adjusted repeatabilities
mod_weight_wing_forest_plot_rptR <-
  allCoefs_mod %>%
  filter(str_detect(effect, "repeat")) %>%
  mutate(comp_name = fct_relevel(comp_name,
                                 "Residual",
                                 "Individual")) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbarh(aes(xmin = conf.low,
                     xmax = conf.high,
                     y = comp_name),
                 alpha = 1, color = col_all, 
                 size = 0.5,
                 height = 0) +
  geom_point(aes(y = comp_name, x = estimate),
             size = 3, shape = 21, 
             fill = "#ECEFF4", col = col_all, 
             alpha = 1, stroke = 0.5) +
  luke_theme +
  theme(axis.title.x = element_text(size = 10, hjust = 0.5)) +
  ylab("Intra-class\ncorrelation") +
  xlab(expression(italic(paste("              Adjusted repeatability (r)" %+-% "95% CI", sep = ""))))

# Patchwork plot
mod_weight_wing_forest_plot_combo <-
  (mod_weight_wing_forest_plot_fixef / mod_weight_wing_forest_plot_partR2 / 
     # mod_weight_wing_forest_plot_randef /
     mod_weight_wing_forest_plot_rptR) + 
  plot_annotation(tag_levels = 'A', title = 'male pectoral sandpiper body mass dynamics', 
                  theme = theme(plot.title = element_text(face = 'italic', hjust = 0.85))) +
  plot_layout(heights = unit(c(4.5, 4, 
                               # 2.5,
                               2.5), c('cm', 'cm', 
                                       # 'cm',
                                       'cm')))

mod_weight_wing_forest_plot_combo

# export plot to disk
ggsave(plot = mod_weight_wing_forest_plot_combo,
       filename = "R/products/figures/jpg/mod_weight_wing_forest.jpg",
       width = 6.5,
       height = 7, units = "in")

ggsave(plot = mod_weight_wing_forest_plot_combo,
       filename = "R/products/figures/svg/mod_weight_wing_forest.svg",
       width = 6.5,
       height = 7, units = "in")

load("R/output/stats_mod_weight_wing.rds")
# extract fitted values of chick weight v egg volume model
mod_weight_wing_fits <- 
  as.data.frame(effect(term = "wing", mod = stats_mod_weight_wing$mod, 
                       xlevels = list(wing = seq(min(cap_05_09_std_pca[, "wing"], na.rm = TRUE),
                                                   max(cap_05_09_std_pca[, "wing"], na.rm = TRUE), 0.01))))
# model summary a diagnostics
summary(stats_mod_weight_wing$mod)
plot(allEffects(stats_mod_weight_wing$mod))
coefplot2(stats_mod_weight_wing$mod)
summary(glht(stats_mod_weight_wing$mod))

#### Plot: weight v wing length ----
wing_weight_plot <-
  ggplot() +
  # geom_errorbarh(data = cap_05_09_std_pca,
  #                aes(y = weight, x = wing, 
  #                    xmin = weight - sd_egg_volume, 
  #                    xmax = weight + sd_egg_volume), 
  #                alpha = 0.3, size = 0.5, linetype = "solid",
  #                color = brewer.pal(8, "Set1")[c(2)]) +
  geom_errorbar(data = cap_05_09_std_pca %>% 
                  group_by(capture_id, wing) %>% 
                  summarise(mean_weight = mean(weight),
                            max_weight = max(weight),
                            min_weight = min(weight)),
                aes(y = mean_weight, x = wing,
                    ymin = min_weight,
                    ymax = max_weight),
                alpha = 0.3, size = 0.5, width = 0, linetype = "solid",
                color = brewer.pal(8, "Set1")[c(2)]) +
  geom_point(data = 
               cap_05_09_std_pca %>% 
               group_by(capture_id, wing) %>% 
               summarise(mean_weight = mean(weight),
                         sd_weight = sd(weight)),
             aes(x = wing, y = mean_weight),
             alpha = 0.4,
             shape = 19, #21, 
             color = brewer.pal(8, "Set1")[c(2)]) +
  geom_line(data = mod_weight_wing_fits, aes(x = wing, y = fit),
            lwd = 0.5) +
  geom_ribbon(data = mod_weight_wing_fits, aes(x = wing, 
                                               ymax = upper, ymin = lower),
              lwd = 1, alpha = 0.25) +
  luke_theme +
  theme(panel.border = element_blank(),
        plot.margin = margin(0, 0, 0, 0.5, "cm"),
        axis.title.y = element_text(vjust = 5)) +
  scale_y_continuous(limits = c(min(cap_05_09_std_pca$weight, na.rm = TRUE), 
                                max(cap_05_09_std_pca$weight, na.rm = TRUE) * 1.05)) +
  ylab("body mass (g; mean and range)") +
  xlab("wing length (mm)")

ggsave(plot = wing_weight_plot,
       filename = "R/products/figures/svg/wing_weight_plot.svg",
       width = 5.29*2,
       height = 5.29*2, units = "cm")

wing_dist_plot <-
  cap_05_09_std_pca %>% 
  ggplot() + 
  geom_boxplot(aes(x = wing, y = 31), 
               fill = brewer.pal(8, "Set1")[c(2)], 
               color = brewer.pal(8, "Set1")[c(2)],
               width = 2, alpha = 0.5) +
  geom_jitter(aes(x = wing, y = 28), 
              fill = brewer.pal(8, "Set1")[c(2)], 
              color = brewer.pal(8, "Set1")[c(2)],
              height = 1, alpha = 0.5) +
  geom_histogram(alpha = 0.5, aes(wing), 
                 fill = brewer.pal(8, "Set1")[c(2)], 
                 binwidth = 1) +
  luke_theme +
  theme(axis.title.x = element_text(hjust = 0.4),
        panel.border = element_blank(),
        plot.margin = margin(0, 0, 0, 0.5, "cm")) +
  ylab("Number of males                        ") +
  xlab("wing length (mm)") +
  # scale_x_reverse(limits = c(25/10, 20/10),
  #                 breaks = c(seq(from = 2, to = 2.4, by = 0.1))) +
  scale_y_continuous(limits = c(0, 32),
                     breaks = c(0, 5, 10, 15, 20, 25), position = "left")

weight_dist_plot <-
  cap_05_09_std_pca %>% 
  ggplot() + 
  geom_boxplot(aes(x = wing, y = 31), 
               fill = brewer.pal(8, "Set1")[c(2)], 
               color = brewer.pal(8, "Set1")[c(2)],
               width = 2, alpha = 0.5) +
  geom_jitter(aes(x = wing, y = 28), 
              fill = brewer.pal(8, "Set1")[c(2)], 
              color = brewer.pal(8, "Set1")[c(2)],
              height = 1, alpha = 0.5) +
  geom_histogram(alpha = 0.5, aes(wing), 
                 fill = brewer.pal(8, "Set1")[c(2)], 
                 binwidth = 1) +
  luke_theme +
  theme(axis.title.x = element_text(hjust = 0.4),
        panel.border = element_blank(),
        plot.margin = margin(0, 0, 0, 0.5, "cm")) +
  ylab("Number of males                        ") +
  xlab("wing length (mm)") +
  # scale_x_reverse(limits = c(25/10, 20/10),
  #                 breaks = c(seq(from = 2, to = 2.4, by = 0.1))) +
  scale_y_continuous(limits = c(0, 32),
                     breaks = c(0, 5, 10, 15, 20, 25), position = "left")
