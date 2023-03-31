
#! NOTE
  # Extract data from scidb; prepare data.

#! TOOLS, SETTINGS
  sapply(
  c('data.table',
  'dbo',
  'stringr'
  ,'here',
  'glue'),
  require, character.only = TRUE, quietly = TRUE)

#DATA: fetch captures, paternity, tenure from database
  cap <- dbq(q = "SELECT DISTINCT c.year_, c.ID, start_capture_date_time capdt, gps_date_time gpsdt, weight,
                  totalHead, tarsus, wing
                FROM PESAatBARROW.CAPTURES c,
                PESAatBARROW.SEX s
                  WHERE
                    c.year_ between 2005 and 2009 AND -- subset  for breeding study seasons
                    s.ID = c.ID AND
                    s.sex = 1  -- males
                    ")

  # Paternity 
  pat <- dbq(q = "SELECT year_,  COUNT(DISTINCT IDmother) AS N_females , COUNT(IDmother) AS N_young, IDfather ID 
                    FROM PESAatBARROW.PATERNITY
						          WHERE IDfather IS NOT NULL AND IDfather IN (SELECT ID FROM PESAatBARROW.SEX WHERE sex = 1)
						              GROUP BY IDfather, year_")
  pat[, N_females := as.integer(N_females)]
  pat[, N_young := as.integer(N_young)]

  # Tenure
  ten <- dbq(q = "SELECT A.year_, A.ID, A.combo, B.lastCapt, A.firstCapt
  FROM
    (
      SELECT year_, ID, FUNCTIONS.combo(ul, ll, ur, lr) combo, min(start_capture_date_time) firstCapt FROM PESAatBARROW.CAPTURES
      WHERE year_ between 2005 and 2009 AND
      ID in (SELECT ID from PESAatBARROW.SEX where sex = 1)
          GROUP BY year_, FUNCTIONS.combo(ul, ll, ur, lr)

          ) A

    LEFT JOIN
    (

      SELECT year_, FUNCTIONS.combo(ul, ll, ur, lr) combo, max(gps_date_time) lastCapt FROM PESAatBARROW.RESIGHTINGS
      WHERE year_ between 2005 and 2009
        GROUP BY year_, FUNCTIONS.combo(ul, ll, ur, lr)

        ) B
          ON A.combo = B.combo AND A.year_ = B.year_")

   ten[, tenureDays := difftime(lastCapt, firstCapt, units = "days")|>as.numeric()]


#DATA: prepare captures
  # n by ID & year
  cap[, n_by_id := .N, .(year_, ID)]
  # gps dt is missing for 3 individuals (we'll use capture dt)
  cap[is.na(gpsdt), gpsdt := capdt]

  # rank id by time of capture
  setorder(cap, year_, ID, capdt)
  cap[, capture_id := 1:.N, .(ID, year_)]


  # transform cap to wide format (for males captured at least twice)
  wcap <- dcast(cap[n_by_id > 1 & !is.na(weight)], ID + year_ ~ capture_id, value.var = c("weight", "gpsdt"))

  wcap[, delta_weight := weight_2 - weight_1]
  wcap[, delta_time := difftime(gpsdt_2, gpsdt_1, units = "days")]

#DATA: prepare tenures
  # there are 3 males with negative tenures (resighted right after capture)
  ten[tenureDays <= 0, tenureDays := NA]
  
  # we'll assume that birds never seen again had minimum tenure
  ten[is.na(tenureDays), tenureDays := min(ten$tenureDays, na.rm = TRUE)]