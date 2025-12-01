# install.packages("sqldf")
rm(list=ls())
setwd("C:/Users/LENOVO/Desktop/STAT 656/NBA/csv")
library(sqldf)
library(dplyr)

game <- read.csv("game.csv", header = TRUE)
game_info <- read.csv("game_info.csv", header = TRUE)
game_summary <- read.csv("game_summary.csv", header = TRUE)
line_score <- read.csv("line_score.csv", header = TRUE)

# Remove duplicate rows
game_info_1row <- game_info %>% distinct(game_id, .keep_all = TRUE)
game_summary_1row <- game_summary %>% distinct(game_id, .keep_all = TRUE)
line_score_1row <- line_score %>% distinct(game_id, .keep_all = TRUE)

# Remove overlapping columns from line_score (keep only unique quarter data)
line_score_clean <- line_score_1row %>%
  select(-team_id_home, -team_id_away, 
         -team_abbreviation_home, -team_abbreviation_away,
         -pts_home, -pts_away)

# Remove overlapping columns from game_info
game_info_clean <- game_info_1row %>%
  select(-game_date)

# Remove overlapping columns from game_summary 
game_summary_clean <- game_summary_1row %>%
  select(-game_date_est, -game_sequence)

# Now merge
join1 <- game %>% left_join(game_info_clean, by = "game_id")
join2 <- join1 %>% left_join(game_summary_clean, by = "game_id")
join3 <- join2 %>% left_join(line_score_clean, by = "game_id")

readr::write_excel_csv(join3, "game_gameinfo_gamesummary_linescore.csv")

