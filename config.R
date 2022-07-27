## Common shared constants

# available egg grades and their ranges in gram
ordered.grades <- c("Jumbo", "XLarge", "Large")
grades <- c("Oversize","SJumbo", "Jumbo", "XLarge", "55", "Large", "Medum", "Undersize")
ranges.1 <- list(Medum=c(41.7, 50), Large=c(50, 54), '55'=c(54, 58.2), XLarge=c(58.2, 66.7), Jumbo=c(66.7, 68.3), SJumbo=c(68.3, 79))
ranges.2 <- list(Medum=c(41.7, 50), Large=c(50, 54), '55'=c(54, 58.2), XLarge=c(58.2, 66.7), Jumbo=c(66.7, 74), SJumbo=c(74, 79))

# egg supplier id for each trial
suppliers <- list(trial1=2287, trial2=2287, trial3=c(2288,2289), trial4=2290, trial5=2290, trial6=2290)

# supply dist by gram and grade
gram_path <- "./data/supply_egg_dist_per_gram_10Jun2022.csv"
grade_path <- "./data/supply_egg_dist_per_grade_10Jun2022.csv"

# distribution intermediate screenshots
snapshot_path <- "./data/distribution_snapshots.csv"

# trial-order conversion table path
trial2order_path <- "./data/trial_to_order_id_2.csv"

# track-trace path
track_trace_path <- "./data/KFP-MOBA-Track-Trace-10Jun2022.csv"

# lane-detail path
lane_detail_path <- "./data/Lane-Running-Details-10Jun2022-Ben.csv"