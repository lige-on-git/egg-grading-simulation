## Step 2 of the trial data analysis for the conference paper
## Develop the mixture distribution of truncated normal distributions.

setwd("/media/lige/Samsung2TB/Projects of Lige/Egg Package Project/Egg_packaging_Lige/June_Factory_Trial")
source("./formal_analysis/utils.R")
library(plotly)  # very nice 3D plot - need to install packages openssl, httr, and plotly in order (need apt install some packages too)

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
rules.2 <- list(Jumbo=c("SJumbo"), XLarge=c(), Large=c("XLarge","Jumbo"))

# store means of simulated eggs packed in the Large order and the N2/N1 ratios
means.L.2 <- c()
ratios.21.2 <- c()

for(ratio in 0.01*10^(seq(0,2,0.1))){
  # ratio: N2/N1 - treat as a constant (between 0 and 1? since n(SJ)<n(J)<n(XL), I anticipate N1<N2)
  # initial analysis
  ordered_num.total2 <- trial.2[['large_order']][['eggs']] + trial.2[['xlarge_order']][['eggs']] + trial.2[['jumbo_order']][['eggs']]
  trial.2_results <- trial.run(ordered_num.total2, mu.supply2, var.supply2, trial.2, rules.2,
                               mode='exploration', ratios=c(0,ratio))  # only Large order accepts two downgraded grades
  
  # auto check which grade has insufficient supply
  grade.insufficient <- check_insufficient_eggs(trial.2_results)
  
  if(!is.null(grade.insufficient)){
    # offset insufficient egg number
    deficit_num <- -1*trial.2_results[['egg_number_left']][[grade.insufficient]][['n_hat']]
    grade_prop <- egg_grade_initial_proportion(grade.insufficient, mu.supply2, var.supply2)
    extra_total_num <- deficit_num/grade_prop
    
    # second analysis after update the total egg number
    ordered_num.total2 <- round(ordered_num.total2 + extra_total_num)  # use integer for binomial sampling
    trial.2_results <- trial.run(ordered_num.total2, mu.supply2, var.supply2, trial.2, rules.2,
                                 mode='exploration', ratios=c(0,ratio))
  }
  # run simulation to check if the current ratio is correct
  packed_eggs <- trial.reproduce.wrapper(trial.2_results, trial.2, downgrade_rules=rules.2)

  means.L.2 <- append(means.L.2, mean(packed_eggs[["Large"]]))
  ratios.21.2 <- append(ratios.21.2, ratio)
}

# check which ratio can result in the most accurate Large order mean
large.mean <- trial.2$large_order$mean
plot(ratios.21.2, means.L.2)
abline(h=large.mean)

means.L.2.copy <- means.L.2
means.L.2.copy <- abs(means.L.2.copy - large.mean)
best.ratio21.2 <- ratios.21.2[which.min(means.L.2.copy)]

# using the "best" N2/N1 ratio to estimate the "best" downgrade probability
ordered_num.total2 <- trial.2[['large_order']][['eggs']] + trial.2[['xlarge_order']][['eggs']] + trial.2[['jumbo_order']][['eggs']]
trial.2_results <- trial.run(ordered_num.total2, mu.supply2, var.supply2, trial.2, list(Jumbo=c("SJumbo"), XLarge=c(), Large=c("XLarge","Jumbo")),
                             mode='exploration', ratios=c(0,best.ratio21))  # only Large order accepts two downgraded grades;

grade.insufficient <- check_insufficient_eggs(trial.2_results)

deficit_num <- -1*trial.2_results[['egg_number_left']][[grade.insufficient]][['n_hat']]
grade_prop <- egg_grade_initial_proportion(grade.insufficient, mu.supply2, var.supply2)
extra_total_num <- deficit_num/grade_prop

ordered_num.total2 <- round(ordered_num.total2 + extra_total_num)  # use integer for binomial sampling
trial.2_results <- trial.run(ordered_num.total2, mu.supply2, var.supply2, trial.2, list(Jumbo=c("SJumbo"), XLarge=c(), Large=c("XLarge","Jumbo")),
                             mode='exploration', ratios=c(0,best.ratio21))


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
rules.5 <- list(Jumbo=c("SJumbo"), XLarge=c("Jumbo", "SJumbo"), Large=c("XLarge","Jumbo"))

# store means of simulated eggs packed in the Large order and the N2/N1 ratios
means.L.5 <- c()
ratios.21.5.XL <- c()  # ratio to downgrade to XLarge order
ratios.21.5.L <- c()   # ratio to downgrade to Large order

for(ratio1 in 0.01*10^(seq(0,2,0.1))){
  for(ratio2 in 0.01*10^(seq(0,2,0.1))){
    # ratio: N2/N1 - treat as a constant (between 0 and 1? since n(SJ)<n(J)<n(XL), I anticipate N1<N2)
    # initial analysis
    
    ordered_num.total5 <- trial.5[['large_order']][['eggs']] + trial.5[['xlarge_order']][['eggs']] + trial.5[['jumbo_order']][['eggs']]
    trial.5_results <- trial.run(ordered_num.total5, mu.supply5, var.supply5, trial.5, rules.5,
                                 mode='exploration', ratios=c(ratio1,ratio2))  # both XLarge (ratio1) and Large (ratio2) orders accept two downgraded grades)
                                
    # auto check which grade has insufficient supply
    grade.insufficient <- check_insufficient_eggs(trial.5_results)
    
    if(!is.null(grade.insufficient)){
      # offset insufficient egg number
      deficit_num <- -1*trial.5_results[['egg_number_left']][[grade.insufficient]][['n_hat']]
      grade_prop <- egg_grade_initial_proportion(grade.insufficient, mu.supply5, var.supply5)
      extra_total_num <- deficit_num/grade_prop
      
      # second analysis after update the total egg number
      
      ordered_num.total5 <- round(ordered_num.total5 + extra_total_num)  # use integer for binomial sampling
      trial.5_results <- trial.run(ordered_num.total5, mu.supply5, var.supply5, trial.5, rules.5,
                                   mode='exploration', ratios=c(ratio1,ratio2))
    }
    # run simulation to check if the current ratio is correct
    packed_eggs <- trial.reproduce.wrapper(trial.5_results, trial.5, downgrade_rules=rules.5)
    
    means.L.5 <- append(means.L.5, mean(packed_eggs[["Large"]]))
    ratios.21.5.XL <- append(ratios.21.5.XL, ratio1)
    ratios.21.5.L <- append(ratios.21.5.L, ratio2)
  }
}

# check which ratio can result in the most accurate Large order mean
plot_ly(x=ratios.21.5.L, y=ratios.21.5.XL, z=means.L.5, type="scatter3d", mode="markers", color=means.L.5)

large.mean <- trial.5$large_order$mean
plot(ratios.21.5.L, means.L.5, xlim=c(0, 0.05), ylim=c(57.3,57.5))
abline(h=large.mean) 
plot(ratios.21.5.XL, means.L.5, ylim=c(57.3,57.5))
abline(h=large.mean)

means.L.5.copy <- means.L.5
means.L.5.copy <- abs(means.L.5.copy - large.mean)
best.ratio21.5.XL <- ratios.21.5.XL[which.min(means.L.5.copy)]
best.ratio21.5.L <- ratios.21.5.L[which.min(means.L.5.copy)]

# using the "best" N2/N1 ratio to estimate the "best" downgrade probability
ordered_num.total5 <- trial.5[['large_order']][['eggs']] + trial.5[['xlarge_order']][['eggs']] + trial.5[['jumbo_order']][['eggs']]
trial.5_results <- trial.run(ordered_num.total5, mu.supply5, var.supply5, trial.5, rules.5,
                             mode='exploration', ratios=c(best.ratio21.5.XL,best.ratio21.5.L))  # both XLarge (ratio1) and Large (ratio2) orders accept two downgraded grades)

grade.insufficient <- check_insufficient_eggs(trial.5_results)

deficit_num <- -1*trial.5_results[['egg_number_left']][[grade.insufficient]][['n_hat']]
grade_prop <- egg_grade_initial_proportion(grade.insufficient, mu.supply5, var.supply5)
extra_total_num <- deficit_num/grade_prop

ordered_num.total5 <- round(ordered_num.total5 + extra_total_num)  # use integer for binomial sampling
trial.5_results <- trial.run(ordered_num.total5, mu.supply5, var.supply5, trial.5, rules.5,
                             mode='exploration', ratios=c(best.ratio21.5.XL,best.ratio21.5.L))



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


## ________________________________________________
## Significance test of two means with error bars
test_diff_using_err <- function(ratio1, ratio1.err, ratio2, ratio2.err){
  # use two-sample t test to test if two variables (with error bars) are significantly equal or not
  t <- abs(ratio1-ratio2)/sqrt((ratio1.err^2+ratio2.err^2)/2)
  df <- 1
  pt(t, df, lower.tail=FALSE)
}
test_diff_using_err(0.0844, 0.0001, 0.0842, 0.0001)
plot(seq(-0.0005,0.0005,0.00001), dnorm(seq(-0.0005,0.0005,0.00001), 0, 0.0001), type='l')
lines(seq(-0.0005,0.0005,0.00001), dnorm(seq(-0.0005,0.0005,0.00001), 0.0003, 0.0001))
