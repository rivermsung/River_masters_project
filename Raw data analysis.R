########################################
### Raw Data Input & Manipulation    ###
### Combining Env. variabls w/ RLC   ###
###       Output: Masterlist         ###
########################################

library(readxl)
library(dplyr)
library(writexl)

library(lubridate)
library(ggplot2)
library(hms)

### ----------------------------------------------------------
### 1. Calculate avg PAR for each RLC
### ----------------------------------------------------------

raw_PAM <- read_excel("2025.04_TI_PAM_R.xlsx", sheet = 1)

#Compute mean Ext. PAR per Unique ID
raw_PAM_updated <- raw_PAM %>%
  group_by(`Unique ID`) %>%
  mutate(`avg PAR` = mean(`Ext. PAR`, na.rm = TRUE)) %>%
  ungroup()

str(raw_PAM_updated)

### ----------------------------------------------------------
### 2. Visualizing PAR wrt Time of Day
### ----------------------------------------------------------
raw_PAR = raw_PAR[complete.cases(raw_PAR), ]
raw_PAR <- raw_PAM_updated %>%
  mutate(
    Time_only = hms::as_hms(format(lubridate::ymd_hms(Time), "%H:%M:%S"))
  )

ggplot(raw_PAR, aes(x = Time_only, y = `Ext. PAR`)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_line(alpha = 0.6) +
  scale_x_time(
    limits = c(as_hms("11:00:00"), as_hms("14:00:00")),
    labels = scales::time_format("%H:%M")
  ) +
  labs(
    x = "Time of Day",
    y = "External PAR (µmol m⁻² s⁻¹)",
    title = "External PAR Over the Course of the Day"
  ) +
  theme_minimal()






