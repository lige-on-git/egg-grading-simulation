## Step 1 of the trial data analysis for the conference paper - Study supply distributions

## 1. Compare if grade count histogram can deduce a good estimation of mu and sigma square. We will later use grade counts to check
## whether egg distribution is stable in different time periods

setwd("/media/lige/Samsung2TB/Projects of Lige/Egg Package Project/Egg_packaging_Lige/June_Factory_Trial")
getwd()
source("./formal_analysis/config.R")
source("./formal_analysis/utils.R")

hist_gram <- fread(gram_path)
hist_grade <- fread(grade_path)

windarra_grade <- hist_grade[suppliername=="Windarra Farm" & (Description %in% grades) & (Eggs != 0), .(Description, Eggs, Weight)]
windarra_gram <- hist_gram[suppliername=="Windarra Farm" & (EggWeight>41.7) & (EggWeight<79), .(EggWeight, Eggs)]

windarra_grade[, sum(Eggs)]  # once include Blood, Crack and Dirt, the total egg numbers will mostly match (can't include Liquid)
windarra_gram[, sum(Eggs)]

find_mean(windarra_gram)     # overall sample mean
find_mean(windarra_grade)
windarra_grade[, .(Description, grade_avg = Weight/10/Eggs)]  # sample mean of each bin

estimate_var(windarra_gram)  # overall sample variance
estimate_var(windarra_grade, mode="mean")
estimate_var(windarra_grade, mode="uniform")

var_dist(windarra_grade, 200)  # sampling distribution of variance

# null hypothesis: variances are equal
# p-value ~0.09. not significant enough to reject. i.e variances are reasonably equal
# HOWEVER, CAN ONLY SHOW VARIANCE EQUALITY OF THE WINDARRA FARM. ENOUGH FOR OUR ANALYSIS. BUT NOT EQUAL IN GENERAL.
# (use gram distribution in general comparisons)
F_stat <- estimate_var(windarra_gram) / var_dist(windarra_grade, 200)[["var_hat"]]; F_stat
pf(F_stat, windarra_gram[, sum(Eggs)-nrow(windarra_gram)], windarra_grade[, sum(Eggs)-nrow(windarra_grade)], lower.tail = FALSE)



## 2. Chi-square test of the overall weight distribution of Windarra Farm: confirm normality
## Since the available data are already in the discrete form (histograms in gram and egg grade),
## the non-parametric Chi-square goodness of fit test is adequate.

# observed probability in bins
bins <- windarra_gram[, EggWeight]; bins
nums <- windarra_gram[, Eggs]; nums
total_num <- windarra_gram[, sum(Eggs)]; total_num
prob_actual <- windarra_gram[, Eggs / sum(Eggs)]; prob_actual

# expected probability in bins based on a normal distribution
mu_est <- find_mean(windarra_gram)
sigma_est <- sqrt(var_dist(windarra_grade, 30)[["var_hat"]])

cdf_lo <- pnorm(bins, mean = mu_est, sd = sigma_est, lower.tail = TRUE)
cdf_hi <- pnorm(bins+1, mean = mu_est, sd = sigma_est, lower.tail = TRUE)
prob_est <- cdf_hi - cdf_lo

# null hypothesis: data is sampled from a normal distribution.
# p-value ~0.24. not significant to reject. i.e it's reasonably to assume normality
chisq.test(prob_actual, prob_est)



## 3. Using three snapshots to find grade count distribution in the middle two periods. Check whether distribution is consistent.
snapshots_dist_summary()



## Trial 2 (Sky Farm Euroa) and 6 (Windarra Farm) distributions

windarra_dist <- hist_gram[suppliername=="Windarra Farm" & (EggWeight>41.7) & (EggWeight<79), .(EggWeight, Eggs)]; windarra_dist
skyfarm_dist <- hist_gram[suppliername=="Sky Farm Euroa" & (EggWeight>41.7) & (EggWeight<79), .(EggWeight, Eggs)]; skyfarm_dist

find_mean(windarra_dist)
find_mean(skyfarm_dist)
estimate_var(windarra_dist)
estimate_var(skyfarm_dist)

# when compare - do error propagation

# trial 6 and 2
62.0; 62.5
23.3; 21.9

# var by gram and by grade
23.3; 23.2+0.05

# var
22.2-0.1; 20.5+0.1; 21.7+/-0.1 (can use t test for these two mean values and there s.e.)

# mean
61.9; 62.6; 62.3 (need a method to estimate uncertainty of the means - use measurement error - for simplicity, we only compute uncertainty for grade distribution)


## Now, the job is to justify that giveaway of trial 2 and 6 have no significant difference 
## get the model build first. if can't justify using distribution, then maybe use error propagation
