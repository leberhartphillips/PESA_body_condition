## install dependent packages

# a vector of all the packages needed in the project's scripts
packages_required_in_project <- 
  c("data.table",
    "dbo",
    "glue",
    "here",
    "lubridate",
    "stringr",
    "tidyverse",
    "vtable",
    "corrplot",
    "lme4",
    "car",
    "broom.mixed",
    "rptR",
    "partR2",
    "effects",
    "gt",
    "RColorBrewer",
    "patchwork"
  )

# of the required packages, check if some need to be installed
new.packages <- 
  packages_required_in_project[!(packages_required_in_project %in% 
                                   installed.packages()[,"Package"])]

# install all packages that are not locally available
if(length(new.packages)) install.packages(new.packages)

# load all the packages into the current R session
lapply(packages_required_in_project, require, character.only = TRUE)