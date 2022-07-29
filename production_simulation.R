## Step 3 of the trial data analysis for the conference paper
## Use estimated downgrade results to reproduce the trials

setwd("/media/lige/Samsung2TB/Projects of Lige/Egg Package Project/Egg_packaging_Lige/June_Factory_Trial")
source("./formal_analysis/utils.R")

# trial 1
packed_eggs <- trial.reproduce.wrapper(trial.1_results, trial.1, downgrade_rules=list(Jumbo=c("SJumbo"), XLarge=c("Jumbo"), Large=c("XLarge")))
ratios.1 <- ratio_samplings(times=30, trial.1_results, trial.1, downgrade_rules=list(Jumbo=c("SJumbo"), XLarge=c("Jumbo"), Large=c("XLarge")))
giveaway.ratios.1 <- ratios.1[['giveaway.ratios']]
downgrade.ratios.1 <- ratios.1[['downgrade.ratios']]

mean(giveaway.ratios.1)
sd(giveaway.ratios.1)
mean(downgrade.ratios.1)
sd(downgrade.ratios.1)


# trial 2
packed_eggs <- trial.reproduce.wrapper(trial.2_results, trial.2, downgrade_rules=list(Jumbo=c("SJumbo"), XLarge=c(), Large=c("XLarge","Jumbo")))
ratios.2 <- ratio_samplings(times=30, trial.2_results, trial.2, downgrade_rules=list(Jumbo=c("SJumbo"), XLarge=c(), Large=c("XLarge","Jumbo")))
giveaway.ratios.2 <- ratios.2[['giveaway.ratios']]
downgrade.ratios.2 <- ratios.2[['downgrade.ratios']]

mean(giveaway.ratios.2)
sd(giveaway.ratios.2)
mean(downgrade.ratios.2)
sd(downgrade.ratios.2)


# trial 4
packed_eggs <- trial.reproduce.wrapper(trial.4_results, trial.4, downgrade_rules=list(Jumbo=c("SJumbo"), XLarge=c(), Large=c()))
ratios.4 <- ratio_samplings(times=30, trial.4_results, trial.4, downgrade_rules=list(Jumbo=c("SJumbo"), XLarge=c(), Large=c()))
giveaway.ratios.4 <- ratios.4[['giveaway.ratios']]
downgrade.ratios.4 <- ratios.4[['downgrade.ratios']]

mean(giveaway.ratios.4)
sd(giveaway.ratios.4)
mean(downgrade.ratios.4)
sd(downgrade.ratios.4)


# trial 5
packed_eggs <- trial.reproduce.wrapper(trial.5_results, trial.5, downgrade_rules=list(Jumbo=c("SJumbo"), XLarge=c("Jumbo", "SJumbo"), Large=c("XLarge","Jumbo")))
ratios.5 <- ratio_samplings(times=30, trial.5_results, trial.5, downgrade_rules=list(Jumbo=c("SJumbo"), XLarge=c("Jumbo", "SJumbo"), Large=c("XLarge","Jumbo")))
giveaway.ratios.5 <- ratios.5[['giveaway.ratios']]
downgrade.ratios.5 <- ratios.5[['downgrade.ratios']]

mean(giveaway.ratios.5)
sd(giveaway.ratios.5)
mean(downgrade.ratios.5)
sd(downgrade.ratios.5)


# trial 6
packed_eggs <- trial.reproduce.wrapper(trial.6_results, trial.6, downgrade_rules=list(Jumbo=c("SJumbo"), XLarge=c("Jumbo", "SJumbo"), Large=c("XLarge","Jumbo")))
ratios.6 <- ratio_samplings(times=30, trial.6_results, trial.6, downgrade_rules=list(Jumbo=c("SJumbo"), XLarge=c("Jumbo", "SJumbo"), Large=c("XLarge","Jumbo")))
giveaway.ratios.6 <- ratios.6[['giveaway.ratios']]
downgrade.ratios.6 <- ratios.6[['downgrade.ratios']]

mean(giveaway.ratios.6)
sd(giveaway.ratios.6)
mean(downgrade.ratios.6)
sd(downgrade.ratios.6)


# to compare egg avg by ordered grades between simulations and actual records
mean(packed_eggs[["Jumbo"]])
mean(packed_eggs[["XLarge"]])
mean(packed_eggs[["Large"]])
mean(packed_eggs[["EOL"]])

trial.2$jumbo_order$mean
trial.2$xlarge_order$mean
trial.2$large_order$mean

# to compare estimated and simulated total consumed egg numbers (include EOL)
ordered_num.total6
length(packed_eggs[["Jumbo"]])+length(packed_eggs[["XLarge"]])+
length(packed_eggs[["Large"]])+length(packed_eggs[["EOL"]])

# check egg distributions in each bucket
hist(packed_eggs[["Jumbo"]])
abline(v=get_grade_bounds("Jumbo")[1], col='red')
abline(v=get_grade_bounds("Jumbo")[2], col='red')

hist(packed_eggs[["XLarge"]])
abline(v=get_grade_bounds("XLarge")[1], col='red')
abline(v=get_grade_bounds("XLarge")[2], col='red')

hist(packed_eggs[["Large"]])
abline(v=get_grade_bounds("Large")[1], col='red')
abline(v=get_grade_bounds("Large")[2], col='red')

hist(packed_eggs[["EOL"]])
