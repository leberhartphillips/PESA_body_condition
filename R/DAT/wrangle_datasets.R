
#! NOTE
  # Extract data from scidb; prepare data.

# load libraries
source(here::here("R/project/project_libraries.R"))

# query scidb for all male PESA morphometric data (i.e., all years)
cap_all <- 
  suppressWarnings({
    dbq(q = "SELECT * FROM PESAatBARROW.CAPTURES") %>% 
      dplyr::select(year_, ID, start_capture_date_time, gps_date_time, weight, 
                    culmen, totalHead, tarsus, wing, thr_width, thr_height) %>% 
      rename(capdt = start_capture_date_time,
             gpsdt = gps_date_time) %>% 
      left_join(., dbq(q = "SELECT * FROM PESAatBARROW.SEX") %>% 
                  dplyr::select(ID, sex), 
                by = "ID", multiple = "all") %>%
      distinct() %>% 
      filter(sex == 1) %>% 
      group_by(ID, year_) %>% 
      mutate(n_by_id = n()) %>% 
      # use cap dt for NA obs in the gps dt
      mutate(gpsdt = ifelse(is.na(gpsdt), capdt, gpsdt)) %>% 
      mutate(gpsdt = as.POSIXct(gpsdt, origin = "1970-01-01")) %>% 
      arrange(year_, ID, capdt) %>% 
      ungroup()
    })

# write.csv(cap_all, file = "R/data/cap_all.csv")

# query scidb for PESA morphometric data collected between 2005 and 2009 (i.e.,
# when repeated measures were taken)
cap_05_09 <-
  cap_all %>% 
  filter(year_ <= 2009 & year_ >= 2005)

# save(cap_05_09,
#      file = "R/DAT/cap_05_09.rds")

cap_05_09_delta <- 
  cap_05_09 %>% 
  filter(n_by_id > 1, !is.na(weight)) %>% 
  mutate(capture_id = paste(ID, year_, sep = "_")) %>% 
  group_by(capture_id) %>% 
  mutate(delta_weight = weight - lag(weight),
         delta_time = gpsdt - lag(gpsdt))

# write.csv(cap_05_09_delta, file = "R/data/cap_05_09_delta.csv")

# query scidb for paternity data
pat <- 
  suppressWarnings({
    dbq(q = "SELECT * FROM PESAatBARROW.PATERNITY") %>% 
      group_by(IDfather, year_) %>% 
      summarise(N_females = n_distinct(IDmother),
                N_young = n_distinct(IDchick)) %>% 
      filter(!is.na(IDfather)) %>% 
      mutate(N_females = as.integer(N_females),
             N_young = as.integer(N_young)) %>% 
      ungroup()
  })

# query scidb for tenure data
# ten <- 
#   dbq(q = "SELECT * FROM PESAatBARROW.CAPTURES") %>% 
#   dplyr::select(year_, ID, ul, ll, ur, lr, start_capture_date_time) %>% 
#   left_join(., dbq(q = "SELECT * FROM PESAatBARROW.SEX") %>% 
#               select(ID, sex), 
#             by = "ID", multiple = "all") %>%
#   distinct() %>% 
#   filter(sex == 1)
# 
# dbq(q = "SELECT * FROM PESAatBARROW.RESIGHTINGS") %>% 
#   dplyr::select(year_, UL, LL, UR, LR, gps_date_time) %>% 
#   
# 
# 
# ten <- dbq(q = "SELECT A.year_, A.ID, A.combo, B.lastCapt, A.firstCapt
# FROM
#   (
#     SELECT year_, ID, FUNCTIONS.combo(ul, ll, ur, lr) combo, min(start_capture_date_time) firstCapt FROM PESAatBARROW.CAPTURES
#     WHERE year_ between 2005 and 2009 AND
#     ID in (SELECT ID from PESAatBARROW.SEX where sex = 1)
#         GROUP BY year_, FUNCTIONS.combo(ul, ll, ur, lr)
# 
#         ) A
# 
#   LEFT JOIN
#   (
# 
#     SELECT year_, FUNCTIONS.combo(ul, ll, ur, lr) combo, max(gps_date_time) lastCapt FROM PESAatBARROW.RESIGHTINGS
#     WHERE year_ between 2005 and 2009
#       GROUP BY year_, FUNCTIONS.combo(ul, ll, ur, lr)
# 
#       ) B
#         ON A.combo = B.combo AND A.year_ = B.year_")
# 
#  ten[, tenureDays := difftime(lastCapt, firstCapt, units = "days")|>as.numeric()]

  # transform cap to wide format (for males captured at least twice)
  # wcap <- dcast(cap_05_09[cap_05_09$n_by_id > 1 & !is.na(cap_05_09$weight)], ID + year_ ~ capture_id, value.var = c("weight", "gpsdt"))
  # 
  # wcap[, delta_weight := weight_2 - weight_1]
  # wcap[, delta_time := difftime(gpsdt_2, gpsdt_1, units = "days")]

# #DATA: prepare tenures
#   # there are 3 males with negative tenures (resighted right after capture)
#   ten[tenureDays <= 0, tenureDays := NA]
#   
#   # we'll assume that birds never seen again had minimum tenure
#   ten[is.na(tenureDays), tenureDays := min(ten$tenureDays, na.rm = TRUE)]