library(dplyr)
library(lubridate)

# read in data
data <- read.csv('game_gameinfo_gamesummary_linescore.csv')
# initial checks
head(data)
names(data)

# initial cleanup; filter for correct years
data_filtered <- data %>%
  mutate(date = as.Date(game_date)) %>%
  filter(
    year(date) >= 2014,
    month(date) >= 7,
    season_type == "Regular Season"
  )

View(data_filtered)
# 4299 games
nrow(data_filtered)
