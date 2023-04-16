# 
# cap_all <- 
#   suppressWarnings({
#     dbq(q = "SELECT * FROM PESAatBARROW.CAPTURES") %>% 
#       dplyr::select(year_, ID, start_capture_date_time, gps_date_time, weight, 
#                     culmen, totalHead, tarsus, wing, thr_width, thr_height) %>% 
#       rename(capdt = start_capture_date_time,
#              gpsdt = gps_date_time) %>% 
#       left_join(., dbq(q = "SELECT * FROM PESAatBARROW.SEX") %>% 
#                   dplyr::select(ID, sex), 
#                 by = "ID", multiple = "all") %>%
#       distinct() %>% 
#       filter(sex == 1) %>% 
#       group_by(ID, year_) %>% 
#       mutate(n_by_id = n()) %>% 
#       # use cap dt for NA obs in the gps dt
#       mutate(gpsdt = ifelse(is.na(gpsdt), capdt, gpsdt)) %>% 
#       mutate(gpsdt = as.POSIXct(gpsdt, origin = "1970-01-01")) %>% 
#       arrange(year_, ID, capdt) %>% 
#       ungroup()
#   })
# 
# cap_05_09 <-
#   cap_all %>% 
#   filter(year_ <= 2009 & year_ >= 2005)

load("/Users/luketheduke2/Documents/Academic/Postdoc_Bart/projects/PESA_body_condition/R/DAT/cap_05_09.rds")
load("/Users/luketheduke2/Documents/Academic/Postdoc_Bart/projects/PESA_body_condition/R/DAT/res.rds")
load("/Users/luketheduke2/Documents/Academic/Postdoc_Bart/projects/PESA_body_condition/R/DAT/cap.rds")

# extract dates of first resighting
first_ress_05_09 <- 
  PESAatBARROW.RESIGHTINGS %>% 
  dplyr::select(year_, UL, LL, UR, LR, gps_date_time) %>% 
  group_by(year_, UL, LL, UR, LR) %>% 
  mutate(n_obs = n()) %>% 
  group_by(year_, UL, LL, UR, LR) %>% 
  arrange(gps_date_time) %>% 
  slice(1) %>% 
  rename(first_gpsdt_res = gps_date_time) %>% 
  filter(year_ <= 2009 & year_ >= 2005) %>% 
  dplyr::select(-n_obs) %>% 
  ungroup()

# extract dates of last resighting
last_ress_05_09 <- 
  PESAatBARROW.RESIGHTINGS %>% 
  dplyr::select(year_, UL, LL, UR, LR, gps_date_time) %>% 
  group_by(year_, UL, LL, UR, LR) %>% 
  mutate(n_obs = n()) %>% 
  group_by(year_, UL, LL, UR, LR) %>% 
  arrange(desc(gps_date_time)) %>% 
  slice(1) %>% 
  rename(last_gpsdt_res = gps_date_time) %>% 
  filter(year_ <= 2009 & year_ >= 2005) %>% 
  dplyr::select(-n_obs) %>% 
  ungroup()

# get the combo-ring key
ring_combo_key <- 
  PESAatBARROW.CAPTURES %>% 
  dplyr::select(year_, ID, ul, ll, ur, lr) %>% 
  mutate(capture_id = paste(ID, year_, sep = "_")) %>% 
  dplyr::select(-ID) %>% 
  rename(UL = ul,
         LL = ll,
         UR = ur,
         LR = lr)

# extract dates of first capture
first_caps_05_09 <- 
  PESAatBARROW.CAPTURES %>% 
  dplyr::select(year_, ID, ul, ll, ur, lr, gps_date_time, start_capture_date_time) %>% 
  rename(capdt = start_capture_date_time,
         gpsdt = gps_date_time) %>% 
  mutate(gpsdt = ifelse(is.na(gpsdt), capdt, gpsdt)) %>% 
  mutate(gpsdt = as.POSIXct(gpsdt, origin = "1970-01-01")) %>% 
  rename(UL = ul,
         LL = ll,
         UR = ur,
         LR = lr) %>% 
  dplyr::select(-c(capdt)) %>% 
  group_by(year_, ID, UL, LL, UR, LR) %>% 
  mutate(n_obs = n()) %>% 
  arrange(gpsdt) %>% 
  slice(1) %>% 
  rename(first_gpsdt_cap = gpsdt) %>% 
  filter(year_ <= 2009 & year_ >= 2005) %>% 
  dplyr::select(-n_obs)

# extract dates of last capture
last_caps_05_09 <- 
  PESAatBARROW.CAPTURES %>% 
  dplyr::select(year_, ID, ul, ll, ur, lr, gps_date_time, start_capture_date_time) %>% 
  rename(capdt = start_capture_date_time,
         gpsdt = gps_date_time) %>% 
  mutate(gpsdt = ifelse(is.na(gpsdt), capdt, gpsdt)) %>% 
  mutate(gpsdt = as.POSIXct(gpsdt, origin = "1970-01-01")) %>% 
  rename(UL = ul,
         LL = ll,
         UR = ur,
         LR = lr) %>% 
  dplyr::select(-c(capdt)) %>% 
  group_by(year_, ID, UL, LL, UR, LR) %>% 
  mutate(n_obs = n()) %>% 
  arrange(desc(gpsdt)) %>% 
  slice(1) %>% 
  rename(last_gpsdt_cap = gpsdt) %>% 
  filter(year_ <= 2009 & year_ >= 2005) %>% 
  dplyr::select(-n_obs)

# join dataframes and calculate tenure and clean up
tenure_05_09 <- 
  left_join(first_ress_05_09, last_ress_05_09, multiple = "all") %>%
  left_join(., ring_combo_key, multiple = "all") %>% 
  left_join(., first_caps_05_09, multiple = "all") %>% 
  left_join(., last_caps_05_09, multiple = "all") %>% 
  group_by(capture_id) %>% 
  mutate(min_date = min(first_gpsdt_res, last_gpsdt_res, first_gpsdt_cap, last_gpsdt_cap, na.rm = TRUE),
         max_date = max(first_gpsdt_res, last_gpsdt_res, first_gpsdt_cap, last_gpsdt_cap, na.rm = TRUE)) %>% 
  mutate(tenure = as.numeric(max_date - min_date)) %>% 
  filter(tenure < 500) %>% 
  filter(!is.na(capture_id)) %>% 
  distinct() %>% 
  group_by(capture_id) %>% 
  mutate(n_obs = n()) %>% 
  arrange(desc(n_obs)) %>% 
  ungroup()

# plot of tenure distributions by year
tenure_05_09 %>% 
  ggplot() +
  geom_histogram(aes(tenure)) +
  facet_grid(year_ ~ .)

# join the tenure data with the data of repeated measures
cap_05_09_std_pca_pat_delta_res_ten <-
  tenure_05_09 %>% 
  select(-year_) %>% 
  left_join(cap_05_09_std_pca_pat_delta_res, ., by = "capture_id", multiple = "all") %>% 
  group_by(capture_id) %>% 
  mutate(n_obs = n()) %>% 
  ungroup() %>% 
  arrange(desc(n_obs))

# tenure distributions of the repeated measures dataset
cap_05_09_std_pca_pat_delta_res_ten %>%
  ggplot() +
  geom_histogram(aes(tenure)) +
  facet_grid(year_ ~ .)

# prepare the data that includes all individuals regardless of if they were sampled twice or not
# first extract the annual phenologies
cap_05_09_phenology <-
  cap_05_09 %>% 
  mutate(capture_id = paste(ID, year_, sep = "_")) %>% 
  filter(!is.na(gpsdt)) %>% 
  mutate(gpsdt = as.POSIXct(gpsdt, origin = "1970-01-01")) %>% 
  mutate(gpsdt_year_start = ymd_hms(paste0(year(gpsdt), "-01-01 00:00:00"))) %>% 
  mutate(gpsdt_since_year_start = as.numeric(gpsdt - gpsdt_year_start)) %>% 
  group_by(ID) %>% 
  arrange(ID, gpsdt) %>% 
  slice(1) %>% 
  arrange(year_) %>% 
  group_by(year_) %>% 
  summarize(mean_date = mean(gpsdt_since_year_start),
            sd_date = sd(gpsdt_since_year_start))

# standardize dates by year
cap_05_09_std_all <- 
  cap_05_09 %>% 
  filter(!is.na(weight)) %>%
  mutate(capture_id = paste(ID, year_, sep = "_")) %>% 
  dplyr::select(-capdt, -ID) %>% 
  mutate(gpsdt_year_start = ymd_hms(paste0(year(gpsdt), "-01-01 00:00:00"))) %>% 
  mutate(gpsdt_since_year_start = as.numeric(gpsdt - gpsdt_year_start)) %>% 
  left_join(cap_05_09_phenology, by = "year_") %>%
  mutate(gpsdt_std = gpsdt_since_year_start - mean_date) %>% 
  group_by(capture_id) %>% 
  arrange(capture_id, gpsdt_since_year_start) %>% 
  mutate(date_deviance = gpsdt_std - gpsdt_std[which.min(gpsdt_std)],
         first_date = gpsdt_std[which.min(gpsdt_std)],
         last_date = gpsdt_std[which.max(gpsdt_std)]) %>% 
  dplyr::select(year_, capture_id, gpsdt, gpsdt_std, date_deviance, first_date, 
                last_date, weight, culmen, totalHead, tarsus,  wing, n_by_id)

# run the PCA
static_measures_pca_all <-
  cap_05_09_std_all %>% 
  group_by(capture_id) %>% 
  arrange(gpsdt) %>% 
  slice(1) %>% 
  ungroup() %>% 
  na.omit() %>%
  dplyr::select(culmen, totalHead, tarsus, wing) %>%
  princomp()

# check the PCA results
summary(static_measures_pca_all)
biplot(static_measures_pca_all, cex = 0.7)

# bind PC1 to the original dataframe
cap_05_09_std_pca_all <- 
  cap_05_09_std_all %>% 
  group_by(capture_id) %>% 
  arrange(capture_id, gpsdt) %>% 
  slice(1) %>% 
  ungroup() %>% 
  na.omit() %>% 
  dplyr::select(capture_id, culmen, totalHead, tarsus,  wing) %>%
  bind_cols(., static_measures_pca_all$scores[, 1]) %>% 
  rename(structure_pc1 = `...6`) %>% 
  left_join(cap_05_09_std_all %>% 
              dplyr::select(-c(culmen, totalHead, tarsus,  wing)), ., 
            by = "capture_id", multiple = "all") %>% 
  na.omit() %>% 
  mutate(year_ = factor(year_, levels = c(2005, 2006, 2007, 2008, 2009))) %>% 
  mutate(log_weight = log(weight),
         log_wing = log(wing))

# join the paternity data 
cap_05_09_std_pca_pat_delta_all <-
  pat %>%
  mutate(capture_id = paste(IDfather, year_, sep = "_")) %>% 
  dplyr::select(-c(IDfather, year_)) %>% 
  left_join(cap_05_09_std_pca_all, ., by = "capture_id", multiple = "all") %>% 
  # assume that males with NA for paternity have 0
  mutate(no_pat = ifelse(is.na(N_young), 1, 0),
         pat = ifelse(!is.na(N_young), 1, 0),
         no_fem = ifelse(is.na(N_females), 1, 0),
         fem = ifelse(!is.na(N_females), 1, 0)) %>%
  group_by(capture_id) %>% 
  arrange(gpsdt) %>% 
  slice(1) %>%
  ungroup()

# residual index preparation
mod_all <- lm(weight ~ wing, data = cap_05_09_std_pca_pat_delta_all)
summary(glht(mod_all))
plot(allEffects(mod_all))

# check the residuals
ggplot(data.frame(resid = residuals(mod_all), fitted = fitted(mod_all)), aes(x = fitted, y = resid)) +
  geom_point() +
  labs(x = "Fitted Values", y = "Residuals") +
  ggtitle("Residuals vs Fitted Values")

# check that the row number is the same (for sanity!)
nrow(cap_05_09_std_pca_pat_delta_all)
length(residuals(mod_all))

# make a dataframe of the residuals and the unique capture events
cap_05_09_std_pca_pat_delta_res_all_ten <-
  data.frame(capture_id = cap_05_09_std_pca_pat_delta_all$capture_id, 
             gpsdt = cap_05_09_std_pca_pat_delta_all$gpsdt,
             mod_res = residuals(mod_all)) %>% 
  left_join(cap_05_09_std_pca_pat_delta_all, ., by = c("capture_id", "gpsdt"), multiple = "all") %>%
  distinct() %>% 
  left_join(tenure_05_09 %>% dplyr::select(-year_), by = "capture_id", multiple = "all") %>% 
  group_by(capture_id) %>% 
  mutate(n_obs = n()) %>% 
  filter(n_obs == 1)

# mixed model of residuals and tenure
mod_res_pat_bin_all_ten <-
  glmer(cbind(pat, no_pat) ~ mod_res + tenure + (1 | year_),
        data = cap_05_09_std_pca_pat_delta_res_all_ten, family = "binomial")
tbl_regression(mod_res_pat_bin_all_ten, intercept = TRUE)
ggeffect(mod_res_pat_bin_all_ten) %>% plot()

mod_res_ten_bin_all_ten <-
  lmer(tenure ~ weight + (1 | year_),
        data = cap_05_09_std_pca_pat_delta_res_all_ten)
tbl_regression(mod_res_ten_bin_all_ten, intercept = TRUE)
ggeffect(mod_res_ten_bin_all_ten) %>% plot()

ggplot(cap_05_09_std_pca_pat_delta_res_all_ten, 
       aes(y = pat, x = mod_res)) +
  geom_point() +
  geom_smooth(
    method = glm,
    method.args = list(family = "binomial")) +
  # facet_grid(year_ ~ .) +
  luke_theme +
  ylab("probability of siring at least one offspring") +
  xlab("relative body mass (residuals, g)")

ggplot(cap_05_09_std_pca_pat_delta_res_all_ten, 
       aes(y = pat, x = mod_res)) +
  geom_point() +
  geom_smooth(
    method = glm,
    method.args = list(family = "binomial")) +
  # facet_grid(year_ ~ .) +
  luke_theme +
  ylab("probability of siring at least one offspring") +
  xlab("relative body mass (residuals, g)")

ggplot(cap_05_09_std_pca_pat_delta_res_all_ten, 
       aes(y = tenure, x = weight)) +
  geom_point() +
  geom_smooth(
    method = lm) +
  facet_grid(year_ ~ .) +
  luke_theme +
  ylab("tenure") +
  xlab("body mass (g)")

ggplot(cap_05_09_std_pca_pat_delta_res_all_ten, 
       aes(y = pat, x = tenure)) +
  geom_point() +
  geom_smooth(
    method = glm,
    method.args = list(family = "binomial")) +
  # facet_grid(year_ ~ .) +
  luke_theme +
  ylab("probability of siring at least one offspring") +
  xlab("tenure (days)")

ggplot(cap_05_09_std_pca_pat_delta_res_all_ten, 
       aes(y = pat, x = mod_res)) +
  geom_point() +
  geom_smooth(
    method = glm,
    method.args = list(family = "binomial")) +
  facet_grid(year_ ~ .) +
  luke_theme +
  ylab("probability of siring at least one offspring") +
  xlab("relative body mass (residuals, g)")


