# library imports ---------------------------------------------------------

# Set working directory to script location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(alfred)
library(tidyverse)
library(lubridate)


# import data -------------------------------------------------------------

gdp <- get_fred_series("GDPC1") #  Billions of Chained 2012 Dollars, Seasonally Adjusted
core_cpi <- get_fred_series("CPILFESL") # Consumer Price Index for All Urban Consumers: All Items Less Food & Energy
                                        # Index 1982-1984=100, Seasonally Adjusted 

save(gdp,core_cpi, file = "FREDdata.RData")


# load data ---------------------------------------------------------------

load("FREDdata.RData")

# transform data ----------------------------------------------------------

lgdp <- gdp %>% 
  mutate(qtr = quarter(date,with_year = TRUE),
         lgdp = log(GDPC1)) %>% 
  select(qtr, lgdp)

core_infl <- core_cpi %>% 
  mutate(qtr = quarter(date,with_year = TRUE)) %>% 
  group_by(qtr) %>% 
  summarise(avg_core_cpi = mean(CPILFESL)) %>% 
  mutate(core_infl = avg_core_cpi/lag(avg_core_cpi)-1) %>% 
  select(qtr, core_infl)

dta <- lgdp %>% 
  left_join(core_infl, by = "qtr") %>% 
  na.omit()


# export to m file --------------------------------------------------------


fname <- "estdata.m"

write("dta = [", fname)

write.table(dta, file = fname, append = TRUE, quote = FALSE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE)

write("];

y = dta(:,2);
Pi = dta(:,3);

clear dta;", fname, append = TRUE)