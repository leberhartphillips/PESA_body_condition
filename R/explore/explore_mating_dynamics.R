nests <- 
  dbq(q = "SELECT * FROM PESAatBARROW.NESTS") %>% 
  mutate(year_start = ymd_hms(paste0(year(laying_date), "-01-01 00:00:00"))) %>% 
  mutate(julian_laydate = as.numeric(laying_date - year_start))

nests %>% 
  ggplot() +
  geom_histogram(aes(x = julian_laydate), binwidth = 1) +
  facet_grid(year_ ~ .)

str(nests$laying_date)

cap_05_09_std_pca_delta <- 
  cap_05_09_std_pca %>% 
  group_by(capture_id) %>% 
  mutate(delta_weight = weight - lag(weight),
         delta_time = (gpsdt - lag(gpsdt))) %>% 
  mutate(delta_time = as.numeric(delta_time)/60/60/24) %>% 
  na.omit() %>% 
  mutate(mass_change = delta_weight/as.numeric(delta_time)) %>% 
  # mutate(delta_time_trans = log(mass_change + 1)) %>% 
  arrange(desc(mass_change)) %>% 
  group_by(capture_id) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(-c(gpsdt, n_by_id, culmen, totalHead, tarsus)) %>% 
  filter(mass_change<30)

cap_05_09_std_pca_smi_delta <- 
  cap_05_09_std_pca_smi %>% 
  group_by(capture_id) %>% 
  mutate(delta_weight_smi = smi_wing - lag(smi_wing),
         delta_time = (gpsdt - lag(gpsdt))) %>% 
  mutate(delta_time = as.numeric(delta_time)/60/60/24) %>% 
  na.omit() %>% 
  mutate(mass_change = delta_weight_smi/as.numeric(delta_time)) %>% 
  # mutate(delta_time_trans = log(mass_change + 1)) %>% 
  arrange(desc(mass_change)) %>% 
  group_by(capture_id) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(year_, capture_id, delta_time, delta_weight_smi, mass_change, first_date, last_date, smi_wing, wing) %>% 
  filter(mass_change < 30)

hist(cap_05_09_std_pca_smi_delta$mass_change)

cap_05_09_std_pca_smi_delta_pat <- 
  pat %>%
  mutate(capture_id = paste(IDfather, year_, sep = "_")) %>% 
  select(-c(IDfather, year_)) %>% 
  left_join(cap_05_09_std_pca_smi_delta, ., by = "capture_id")

#### do heavier males arrive earlier? ----
mod_arrival_wing_weight <- 
  lmer(first_date ~ smi_wing + (1 | year_),
        data = cap_05_09_std_pca_smi_delta_pat)

# yes, but not significant (p = 0.08)
summary(glht(mod_arrival_wing_weight))
plot(allEffects(mod_arrival_wing_weight))

#### do bigger males arrive earlier? ----
mod_arrival_wing_wing <- 
  lmer(first_date ~ wing + (1 | year_),
       data = cap_05_09_std_pca_smi_delta_pat)

# no, no relationship
summary(glht(mod_arrival_wing_wing))
plot(allEffects(mod_arrival_wing_wing))

#### is the change in body mass over season associated with arrival date? ----
mod_delta_weight <-
  lmer(mass_change ~ first_date + (1 | year_),
       data = cap_05_09_std_pca_smi_delta_pat)

# yes, but not significant (p = 0.08): earlier males increase weight, late males
# loose weight
summary(glht(mod_delta_weight))
plot(allEffects(mod_delta_weight))

#### Do earlier males get more ladies? ----
mod_n_fem_arrival <-
  glmer(N_females ~ first_date + (1 | year_),
       data = cap_05_09_std_pca_smi_delta_pat, family = poisson)

# no relationship
summary(glht(mod_n_fem_arrival))
plot(allEffects(mod_n_fem_arrival))

#### Do earlier males get more offspring? ----
mod_n_young_arrival <-
  glmer(N_young ~ first_date + (1 | year_),
       data = cap_05_09_std_pca_smi_delta_pat, family = poisson)

# no relationship
summary(glht(mod_n_young_arrival))
plot(allEffects(mod_n_young_arrival))

#### is the change in body mass associated with polygyny? ----
mod_n_female_arrival_mass_change <-
  glmer(N_females ~ mass_change + first_date + (1 | year_),
       data = cap_05_09_std_pca_smi_delta_pat, family = poisson)

# no relationship
summary(glht(mod_n_female_arrival_mass_change))
plot(allEffects(mod_n_female_arrival_mass_change))

#### is the change in body mass associated with paternity? ----
mod_n_young_arrival_mass_change <-
  glmer(N_young ~ mass_change + first_date + (1 | year_),
        data = cap_05_09_std_pca_smi_delta_pat, family = poisson)

# no relationship
summary(glht(mod_n_female_arrival_mass_change))
plot(allEffects(mod_n_female_arrival_mass_change))

#### is the change in body mass associated with polygyny? ----
mod_n_females_mass_change <-
  glmer(N_females ~ mass_change + (1 | year_),
        data = cap_05_09_std_pca_smi_delta_pat, family = poisson)

# no relationship
summary(glht(mod_n_females_mass_change))
plot(allEffects(mod_n_females_mass_change))

#### is the change in body mass associated with paternity? ----
mod_n_young_mass_change <-
  glmer(N_young ~ mass_change + (1 | year_),
        data = cap_05_09_std_pca_smi_delta_pat, family = poisson)

# no relationship
summary(glht(mod_n_young_mass_change))
plot(allEffects(mod_n_young_mass_change))

#### is the excess weight associated with polygyny? ----
mod_n_females_first_weight_smi <-
  glmer(N_females ~ smi_wing + (1 | year_),
        data = cap_05_09_std_pca_smi_delta_pat, family = poisson)

# no relationship
summary(glht(mod_n_females_first_weight_smi))
plot(allEffects(mod_n_females_first_weight_smi))

#### is the excess weight associated with paternity? ----
mod_n_young_first_weight_smi <-
  glmer(N_young ~ smi_wing + (1 | year_),
        data = cap_05_09_std_pca_smi_delta_pat, family = poisson)

# no relationship
summary(glht(mod_n_young_first_weight_smi))
plot(allEffects(mod_n_young_first_weight_smi))

#### is body size associated with polygyny? ----
mod_n_females_wing <-
  glmer(N_females ~ wing + (1 | year_),
        data = cap_05_09_std_pca_smi_delta_pat, family = poisson)

# no relationship
summary(glht(mod_n_females_wing))
plot(allEffects(mod_n_females_wing))

#### is body size associated with paternity? ----
mod_n_young_wing <-
  glmer(N_young ~ wing + (1 | year_),
        data = cap_05_09_std_pca_smi_delta_pat, family = poisson)

# yes, but not significant (p = 0.06): larger bodied males sire more offspring
summary(glht(mod_n_young_wing))
plot(allEffects(mod_n_young_wing))

#### 
# mod_n_young_mass_change <-
#   glmer(N_young ~ mass_change + (1 | year_),
#         data = cap_05_09_std_pca_delta_pat, family = poisson)
# 
# summary(glht(mod_n_young_mass_change))
# plot(allEffects(mod_n_young_mass_change))
# 
# mod_n_young_first_weight_mass_change <-
#   glmer(N_young ~ weight + mass_change + (1 | year_),
#       data = cap_05_09_std_pca_delta_pat, family = poisson)
# 
# summary(glht(mod_n_young_first_weight_mass_change))
# plot(allEffects(mod_n_young_first_weight_mass_change))
# summary(mod_n_young_first_weight_mass_change) # singular because year_explains zero variance, keep it in for completeness
# 
# ggplot(data=  cap_05_09_std_pca_delta_pat) +
#   geom_histogram(aes(weight))+
#   facet_grid(year_ ~ .)
# 
# ggplot(data= cap_05_09_std_pca_delta_pat) +
#   geom_point(aes(x = weight, y = mass_change))
# 
# cor.test(cap_05_09_std_pca_delta_pat$weight, cap_05_09_std_pca_delta_pat$mass_change)
