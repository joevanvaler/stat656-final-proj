library(dplyr)
library(lubridate)

# read in data
data <- read.csv('game_gameinfo_gamesummary_linescore.csv')


# initial cleanup; filter for correct years
data_filtered <- data %>%
  mutate(date = as.Date(game_date)) %>%
  filter(
    date <= as.Date("2022-07-01"),
    date >= as.Date("2014-07-01"),
    season_type == "Regular Season"
  )

View(head(data_filtered))
View(tail(data_filtered))
# 4299 games
nrow(data_filtered)
