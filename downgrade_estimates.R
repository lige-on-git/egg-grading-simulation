## Step 2 of the trial data analysis for the conference paper
## Develop the mixture distribution of truncated normal distributions.

setwd("/media/lige/Samsung2TB/Projects of Lige/Egg Package Project/Egg_packaging_Lige/June_Factory_Trial")
source("./formal_analysis/utils.R")

track_trace <- fread(track_trace_path)
lane_detail <- fread(lane_detail_path)

## __________________________________________________
## Test with trial 1 (only accept a single downgrade)
trial.1 <- trial_stats_summary(1, track_trace)
get_simple_ratios(trial.1)
mu.supply <- trial.1[['supply_dist_est']][['mu']]
var.supply <- trial.1[['supply_dist_est']][['var']]

# By comparing A and B, give the following summary:
# Virtually no SJ -> J & no J -> XL downgrades. In these 2 grades most giveaway are contributed by in-grade variation
# Extremely large downgrade from XL -> 55/L

# C. calculate or simulate error of giveaway? Use combined error.

# initial analysis
ordered_num.total <- trial.1[['large_order']][['eggs']] + trial.1[['xlarge_order']][['eggs']] + trial.1[['jumbo_order']][['eggs']]
trial.1_results <- trial.run(ordered_num.total, mu.supply, var.supply, trial.1, list(Jumbo=c("SJumbo"), XLarge=c("Jumbo"), Large=c("XLarge")))

# find that the current initial egg number estimation will result in insufficient XL egg number.
# So, use the XL egg deficit number to point-estimate how many total egg number is required.

deficit_num <- -1*trial.1_results[['egg_number_left']][['XLarge']][['n_hat']]
grade_prop <- egg_grade_initial_proportion("XLarge", mu.supply, var.supply)
extra_total_num <- deficit_num/grade_prop

# second analysis after update the total egg number
ordered_num.total <- round(ordered_num.total + extra_total_num)  # use integer for binomial sampling
trial.1_results <- trial.run(ordered_num.total, mu.supply, var.supply, trial.1, list(Jumbo=c("SJumbo"), XLarge=c("Jumbo"), Large=c("XLarge")))

# manually calculate the final giveaway ratio error (follow formula to do error propagation)
weight_lowerbound <- 50*16105+58.2*37733+66.7*9524
e1 <- err_product(9445,40,69.08,0.01)[['dz']]
e2 <- err_product(79,15,75.44,0.02)[['dz']]
e3 <- err_product(37200,100,62.47,0.01)[['dz']]
e4 <- err_product(535,15,69.08,0.01)[['dz']]
e5 <- err_product(7222,28,55.74,0.01)[['dz']]
e6 <- err_product(8883,30,62.47,0.01)[['dz']]

sqrt(e1^2+e2^2+e3^2+e4^2+e5^2+e6^2)/weight_lowerbound  # error


## _____________________________________________
## Test with trial 2 (accept 2 downgrade grades)
trial.2 <- trial_stats_summary(2, track_trace)
get_simple_ratios(trial.2)
mu.supply2 <- trial.2[['supply_dist_est']][['mu']]
var.supply2 <- trial.2[['supply_dist_est']][['var']]  # same supply distribution as trial 1

# initial analysis
ordered_num.total2 <- trial.2[['large_order']][['eggs']] + trial.2[['xlarge_order']][['eggs']] + trial.2[['jumbo_order']][['eggs']]
trial.2_results <- trial.run(ordered_num.total2, mu.supply2, var.supply2, trial.2, list(Jumbo=c("SJumbo"), XLarge=c(), Large=c("XLarge","Jumbo")))

# insufficient L egg number
deficit_num <- -1*trial.2_results[['egg_number_left']][['Large']][['n_hat']]
grade_prop <- egg_grade_initial_proportion("Large", mu.supply2, var.supply2)
extra_total_num <- deficit_num/grade_prop

# second analysis after update the total egg number
ordered_num.total2 <- round(ordered_num.total2 + extra_total_num)  # use integer for binomial sampling
trial.2_results <- trial.run(ordered_num.total2, mu.supply2, var.supply2, trial.2, list(Jumbo=c("SJumbo"), XLarge=c(), Large=c("XLarge","Jumbo")))

# more L egg number
surplus_num <- +1*trial.2_results[['egg_number_left']][['Large']][['n_hat']]
grade_prop <- egg_grade_initial_proportion("Large", mu.supply2, var.supply2)
less_total_num <- surplus_num/grade_prop

# third analysis after update the total egg number
ordered_num.total2 <- round(ordered_num.total2 - less_total_num)  # use integer for binomial sampling
trial.2_results <- trial.run(ordered_num.total2, mu.supply2, var.supply2, trial.2, list(Jumbo=c("SJumbo"), XLarge=c(), Large=c("XLarge","Jumbo")))


## ________________________________________________
## Test with trial 4 (no downgrades accepted to XL)
trial.4 <- trial_stats_summary(4, track_trace, FALSE)  # display snapshots rather than final accumulated dist; but still use final accumulated dist in calculation
get_simple_ratios(trial.4)
mu.supply4 <- trial.4[['supply_dist_est']][['mu']]
var.supply4 <- trial.4[['supply_dist_est']][['var']]   # same supply distribution as trial 5 & 6

# initial analysis
ordered_num.total4 <- trial.4[['large_order']][['eggs']] + trial.4[['xlarge_order']][['eggs']] + trial.4[['jumbo_order']][['eggs']]
trial.4_results <- trial.run(ordered_num.total4, mu.supply4, var.supply4, trial.4, list(Jumbo=c("SJumbo"), XLarge=c(), Large=c()))

# still insufficient L egg number
deficit_num <- -1*trial.4_results[['egg_number_left']][['Large']][['n_hat']]
grade_prop <- egg_grade_initial_proportion("Large", mu.supply4, var.supply4)
extra_total_num <- deficit_num/grade_prop

# second analysis after update the total egg number
ordered_num.total4 <- round(ordered_num.total4 + extra_total_num)  # use integer for binomial sampling
trial.4_results <- trial.run(ordered_num.total4, mu.supply4, var.supply4, trial.4, list(Jumbo=c("SJumbo"), XLarge=c(), Large=c()))


## ________________________________________________
## Test with trial 5 (accept 2 downgrade grades)
trial.5 <- trial_stats_summary(5, track_trace)
get_simple_ratios(trial.5)
mu.supply5 <- trial.5[['supply_dist_est']][['mu']]
var.supply5 <- trial.5[['supply_dist_est']][['var']]  # same supply distribution as trial 4 & 6

# initial analysis
ordered_num.total5 <- trial.5[['large_order']][['eggs']] + trial.5[['xlarge_order']][['eggs']] + trial.5[['jumbo_order']][['eggs']]
trial.5_results <- trial.run(ordered_num.total5, mu.supply5, var.supply5, trial.5, list(Jumbo=c("SJumbo"), XLarge=c("Jumbo", "SJumbo"), Large=c("XLarge","Jumbo")))

# insufficient L egg number
deficit_num <- -1*trial.5_results[['egg_number_left']][['Large']][['n_hat']]
grade_prop <- egg_grade_initial_proportion("Large", mu.supply5, var.supply5)
extra_total_num <- deficit_num/grade_prop

# second analysis after update the total egg number
ordered_num.total5 <- round(ordered_num.total5 + extra_total_num)  # use integer for binomial sampling
trial.5_results <- trial.run(ordered_num.total5, mu.supply5, var.supply5, trial.5, list(Jumbo=c("SJumbo"), XLarge=c("Jumbo", "SJumbo"), Large=c("XLarge","Jumbo")))


## ________________________________________________
## Test with trial 6 (accept 2 downgrade grades)
trial.6 <- trial_stats_summary(6, track_trace)
get_simple_ratios(trial.6)
mu.supply6 <- trial.6[['supply_dist_est']][['mu']]
var.supply6 <- trial.6[['supply_dist_est']][['var']]  # same supply distribution as trial 4 & 5

# initial analysis
ordered_num.total6 <- trial.6[['large_order']][['eggs']] + trial.6[['xlarge_order']][['eggs']] + trial.6[['jumbo_order']][['eggs']]
trial.6_results <- trial.run(ordered_num.total6, mu.supply6, var.supply6, trial.6, list(Jumbo=c("SJumbo"), XLarge=c("Jumbo", "SJumbo"), Large=c("XLarge","Jumbo")))

# insufficient L egg number
deficit_num <- -1*trial.6_results[['egg_number_left']][['Large']][['n_hat']]
grade_prop <- egg_grade_initial_proportion("Large", mu.supply6, var.supply6)
extra_total_num <- deficit_num/grade_prop

# second analysis after update the total egg number
ordered_num.total6 <- round(ordered_num.total6 + extra_total_num)  # use integer for binomial sampling
trial.6_results <- trial.run(ordered_num.total6, mu.supply6, var.supply6, trial.6, list(Jumbo=c("SJumbo"), XLarge=c("Jumbo", "SJumbo"), Large=c("XLarge","Jumbo")))













# 
