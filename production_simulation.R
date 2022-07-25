## Step 3 of the trial data analysis for the conference paper
## Use estimated downgrade results to reproduce the trials

setwd("/media/lige/Samsung2TB/Projects of Lige/Egg Package Project/Egg_packaging_Lige/June_Factory_Trial")
source("./formal_analysis/config.R")
source("./formal_analysis/utils.R")



times <- 100
g.ratios <- d.ratios <- c()
for(i in (1:times)){
  packed_eggs <- trial.reproduce.wrapper(trial_results=trial.1_results, trial_stats=trial.1, 
                                         downgrade_rules=list(Jumbo=c("SJumbo"), XLarge=c("Jumbo"), Large=c("XLarge")))
  ratio_stats <- get_giveaway_ratio(packed_eggs, trial.1)
  g.ratios <- append(g.ratios, ratio_stats[['overall_giveaway.ratio']])
  d.ratios <- append(d.ratios, ratio_stats[['overall_downgrade.ratio']])
}

sd(g.ratios)
mean(g.ratios)

sd(d.ratios)
mean(d.ratios)



mean(packed_eggs[["Jumbo"]])
mean(packed_eggs[["XLarge"]])
mean(packed_eggs[["Large"]])
mean(packed_eggs[["EOL"]])

length(packed_eggs[["Jumbo"]])+
length(packed_eggs[["XLarge"]])+
length(packed_eggs[["Large"]])+
length(packed_eggs[["EOL"]])

hist(packed_eggs[["EOL"]])





# 