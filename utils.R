## Note: when estimating sampling distribution of a parameter, the order of magnitude of s.e. / uncertainty
## doesn't depends on number of loops, but the sample size inside each loop. So, use the actual size.

## Shared functions
if(!require("data.table")){
  install.packages('data.table')
}
library('data.table')
library('truncnorm')  # https://cran.r-project.org/web/packages/truncnorm/truncnorm.pdf
source("./formal_analysis/config.R")

## _____________________________________________________________________________
## Data Retrieval and Exploration
trial2order <- function(trial_num){
  # return a vector of orders given a trial number
  fread(trial2order_path)[TrialNumber==trial_num, OrderId]
}

id2grade <- function(order_id){
  # return a grade string of "L", "XL", or "J" given an order id
  fread(trial2order_path)[OrderId==order_id, OrderGrade]
}

order_egg_num <- function(order_id, track_trace_table){
  # return how many eggs are packed into an order
  track_trace_table[OrderId==order_id, Eggs][1]
}

order_egg_weight <- function(order_id, track_trace_table){
  # return the total egg weight (gram) packed into an order
  track_trace_table[OrderId==order_id, OrderWeightGram][1]
}

find_mean <- function(histogram_table){
  # histogram_table: either a data.table histogram with two columns - x-axis (gram) and frequency
  #                  or a data.table bar chart with three columns - descriptions, and the corresponding total number of weight
  # return the sample mean as a point estimate
  if (length(histogram_table)==3){
    histogram_table[, sum(Weight/10) / sum(Eggs)]
  }
  else if (length(histogram_table)==2){
    histogram_table[, sum(EggWeight*Eggs) / sum(Eggs)]
  }
  else {
    print("the input data has incorrect structure")
  }
}

estimate_var <- function(histogram_table, mode="mean"){
  # histogram_table: either a data.table histogram with two columns - x-axis (gram) and frequency
  #                  or a data.table bar chart with three columns - descriptions, and the corresponding total number of weight
  # mode: "mean" - assume samples are concentrated at the mean value of bins
  #       "uniform" - assume samples are uniformly distributed within bins (can only choose this if input a bar chart)
  # return an estimated sample variance of normal distribution as a point estimate
  if (length(histogram_table)==2){
    mean <- find_mean(histogram_table)
    return(histogram_table[, sum((EggWeight - mean)^2*Eggs) / (sum(Eggs)-1)])
  }
  else if (length(histogram_table)==3 & mode=="mean"){
    mean <- find_mean(histogram_table)
    return(histogram_table[, sum((Weight/10/Eggs - mean)^2*Eggs) / (sum(Eggs)-1)])
  }
  else if (length(histogram_table)==3 & mode=="uniform"){
    # uniformly sample for egg number times for each bin/grade 
    pool <- c()
    for(grade in names(ranges.1)){
      lo <- ranges.1[[grade]][1]
      hi <- ranges.1[[grade]][2]
      # multinomial sampling of egg number in each bin/grade (instead of using the original egg number in each bin)
      histogram_table[, Eggs_sampled := rmultinom(1, size = histogram_table[, sum(Eggs)], prob = histogram_table[, Eggs / sum(Eggs)])]
      num <- histogram_table[Description==grade, Eggs_sampled]
      samples <- runif(num, lo, hi)
      pool <- append(pool, samples)
    }
    # remove the sampled egg number column
    histogram_table[, Eggs_sampled := NULL]
    return(var(pool))
  }
  else if (mode!="mean" & mode!="uniform"){
    print("the choice of estimation mode is invalid")
  }
  else {
    print("the input data has incorrect structure")
  }
}

var_dist <- function(histogram_table, times=100){
  # histogram_table: a data.table bar chart with three columns - descriptions, and the corresponding total number of weight
  # times: sampling times to estimate the sampling distribution
  # return the mean and s.e. of the sampling distribution of the sample variance of the total supply
  sample_var <- c()
  for(i in 1:times) {
    variance <- estimate_var(histogram_table, mode="uniform")
    sample_var <- append(sample_var, variance)
  }
  return(list(var_hat=mean(sample_var), var_se=sqrt(var(sample_var))))
}

fit_normal_dist <- function(supply_id){
  # fit the overall normal dist parameters based on the gram dist
  # return estimated mu and var
  hist_gram <- fread(gram_path)
  gram_dist <- hist_gram[SupplyId==supply_id & (EggWeight>41.7) & (EggWeight<79), .(EggWeight, Eggs)]
  return(list(mu.overall=find_mean(gram_dist), var.overall=estimate_var(gram_dist)))
}

get_grade_bounds <- function(grade_name){
  # grade_name: select one from c('Large', 'XLarge', 'Jumbo', 'SJumbo)
  # return the lower and higher bounds of a grade
  if (grade_name=="Large"){
    lower_bound <- ranges.2[["Large"]][1]
    higher_bound <- ranges.2[["55"]][2]
  }
  else{
    lower_bound <- ranges.2[[grade_name]][1]
    higher_bound <- ranges.2[[grade_name]][2]
  }
  return(c(lower_bound, higher_bound))
}

egg_grade_initial_proportion <- function(grade_name, mu, var){
  # grade_name: select one from c('Large', 'XLarge', 'Jumbo', 'SJumbo)
  # given a fitted normal distribution (mu & var)
  # return the proportion of a egg grade (theoretical proportions based on an estimated normal dist)
  lower_bound <- get_grade_bounds(grade_name)[1]
  higher_bound <- get_grade_bounds(grade_name)[2]
  return(pnorm(higher_bound,mu,sqrt(var)) - pnorm(lower_bound,mu,sqrt(var)))
}

egg_grade_initial_number <- function(num, grade_name, mu, var, times=100){
  # num: total egg number
  # grade_name: select one from c('Large', 'XLarge', 'Jumbo', 'SJumbo)
  # times: sampling times to estimate the sampling distribution
  # given a fitted normal distribution (mu & var), return an sampling dist of egg numbers of a grade
  grade_prop <- egg_grade_initial_proportion(grade_name, mu, var)
  binom_samples <- rbinom(n=times, size=num, prob=grade_prop)
  return(list(n_hat=mean(binom_samples), n_se=sqrt(var(binom_samples))))
}

egg_grade_sampler <- function(n, mu, var, grade_name){
  # grade_name: select one from c('Large', 'XLarge', 'Jumbo', 'SJumbo)
  # n: simulate how many samples
  # mu: mean of the overall normal distribution
  # var: variance of the overall normal distribution
  # return a vector of values sampled from truncated normal distribution
  lower_bound <- get_grade_bounds(grade_name)[1]
  higher_bound <- get_grade_bounds(grade_name)[2]
  return(rtruncnorm(n, lower_bound, higher_bound, mu, sqrt(var)))
}

grade_mean_dist <- function(grade_name, supply_id, times=100){
  # grade_name: select one from c('Large', 'XLarge', 'Jumbo')
  # supply_id: the supplier id whose egg distribution is used to fit normal distribution
  # times: sampling times to estimate the sampling distribution
  # return the mean and s.e. of the sampling distribution of the sample mean of a grade of supply
  
  # fit normal parameters using gram distribution (here we don't consider measurement errors of mu and var)
  normal_fit <- fit_normal_dist(supply_id)
  mu.overall <- normal_fit[['mu.overall']]
  var.overall <- normal_fit[['var.overall']]
  
  # simulate egg numbers of the input grade using grade distribution (multinomial random sampling)
  hist_grade <- fread(grade_path)
  grade_dist <- hist_grade[SupplyId==supply_id & (Description %in% grades) & (Eggs != 0), .(Description, Eggs, Weight)]
  sample_means <- c()
  
  for(i in 1:times){
    grade_dist[, Eggs_sampled := rmultinom(1, size = grade_dist[, sum(Eggs)], prob = grade_dist[, Eggs / sum(Eggs)])]
    if (grade_name=="Large"){
      # number to randomly sample (aka sample size): s.e. of estimated mean depended on this number
      num2sample <- sum(grade_dist[Description=='55' | Description=="Large", Eggs_sampled])
    }
    else{
      num2sample <- grade_dist[Description==grade_name, Eggs_sampled]
    }
    
    # sample eggs in that grade, and estimate avg weight (truncated normal sampling)
    truncate_samples <-  egg_grade_sampler(num2sample, mu.overall, var.overall, grade_name)
    sample_means <- append(sample_means, mean(truncate_samples))
    
    grade_dist[, Eggs_sampled := NULL]  # remove the temporary multinomial samples
  }
  return(list(mu_hat=mean(sample_means), mu_se=sqrt(var(sample_means))))
}

## _____________________________________________________________________________
## Error Propagation
sum_weight_err <- function(n){
  # n is number of eggs
  # suppose dx=0.1g is the error of an egg weight
  # return the measurement error of the total weight
  dx <- 0.2
  error <- sqrt(dx^2 * n)
}

err_sum <- function(a, x, dx, b, y, dy, add=TRUE){
  # z=ax+by (add=TRUE) or z=ax-by (add=FALSE) 
  # dx is the error of x, and dy is the error of y
  # return z +/- dz
  error <- sqrt((a*dx)^2 + (b*dy)^2)
  if(add){
    list(z=a*x+b*y, dz=error)
  }
  else{
    list(z=a*x-b*y, dz=error)
  }
}

err_product <- function(x, dx, y, dy, product=TRUE){
  # z=xy (product=TRUE) or x/y (product=FALSE)
  # dx is the error of x, and dy is the error of y
  # return z +/- dz
  error <- sqrt((dx/x)^2 + (dy/y)^2)
  if(product){
    list(z=(x*y), dz=error*(x*y))
  }
  else{
    list(z=(x/y), dz=error*(x/y))
  }
}

## _____________________________________________________________________________
## Trial Summary
snapshots_dist_summary <- function(){
  # return distribution stats of three snapshot period (for trial 4, 5, and 5)
  
  snapshots <- fread(snapshot_path)  # Eggs (doze); Weight (kg)
  # accumulated stats
  period.hist0 <- snapshots[Accumulated_egg_number=="6.7k", .(Description, Eggs, Weight)]
  period.hist1 <- snapshots[Accumulated_egg_number=="13.3k", .(Description, Eggs, Weight)]
  period.hist2 <- snapshots[Accumulated_egg_number=="19.6k", .(Description, Eggs, Weight)]
  period.hist3 <- snapshots[Accumulated_egg_number=="25.4k", .(Description, Eggs, Weight)]
  
  # individual stats (subtraction)
  period.hist3[, Eggs := Eggs - period.hist2[, Eggs]][, Weight := Weight - period.hist2[, Weight]]
  period.hist2[, Eggs := Eggs - period.hist1[, Eggs]][, Weight := Weight - period.hist1[, Weight]]
  period.hist1[, Eggs := Eggs - period.hist0[, Eggs]][, Weight := Weight - period.hist0[, Weight]]
  
  # convert to default format (Weight in 0.1g)
  period.hist1[, c("Eggs", "Weight") := .(Eggs*12, Weight*1000*10)]
  period.hist2[, c("Eggs", "Weight") := .(Eggs*12, Weight*1000*10)]
  period.hist3[, c("Eggs", "Weight") := .(Eggs*12, Weight*1000*10)]
  
  # remove invalid rows
  period.hist1 <- period.hist1[Eggs != 0]
  period.hist2 <- period.hist2[Eggs != 0]
  period.hist3 <- period.hist3[Eggs != 0]
  
  dist.mean.1 <- find_mean(period.hist1)  # snapshot distribution means 
  dist.mean.2 <- find_mean(period.hist2)
  dist.mean.3 <- find_mean(period.hist3)
  
  dist.var.1 <- var_dist(period.hist1)  # snapshot distribution vars
  dist.var.2 <- var_dist(period.hist2)
  dist.var.3 <- var_dist(period.hist3) 
  
  # convert to natural format and add grade avg (Weight in 1g)
  period.hist1[, c("Weight", "grade_avg") := .(Weight/10, Weight/10/Eggs)]
  period.hist2[, c("Weight", "grade_avg") := .(Weight/10, Weight/10/Eggs)]
  period.hist3[, c("Weight", "grade_avg") := .(Weight/10, Weight/10/Eggs)]
  
  return(list(snapshot1=list(total_mean=dist.mean.1, total_var=dist.var.1, dist=period.hist1),
              snapshot2=list(total_mean=dist.mean.2, total_var=dist.var.2, dist=period.hist2),
              snapshot3=list(total_mean=dist.mean.3, total_var=dist.var.3, dist=period.hist3)))
}

trial_stats_summary <- function(trial_num, track_trace_table, accum_dist=TRUE, times=100){
  # trial_num: an int in c(1:6) which represents one of the six factory trials.
  # track_trace_table: the recorded track trace data for all the trials.
  # times: sampling times to estimate the sampling distribution
  # accum_dist: if TRUE present the overall supply distribution of a supplier;
  # if FALSE present an intermediate supply distribution of a supplier (only available for trial 4, 5, and 6).
  # return the summary data of a trial
  if(!(trial_num %in% c(1:6))){
    print("Trial number not supported. Must be in c(1:6).")
    return(NULL)
  }
  trial.ids <- trial2order(trial_num)  # a numeric vector of 3 ids in the order of L, XL, and J
  supplyId <- suppliers[[paste0("trial",trial_num)]]
  
  # A. actual packed eggs
  trial.L.number <- order_egg_num(order_id=trial.ids[1], track_trace_table)
  trial.L.weight <- order_egg_weight(order_id=trial.ids[1], track_trace_table)
  trial.XL.number <- order_egg_num(order_id=trial.ids[2], track_trace_table)
  trial.XL.weight <- order_egg_weight(order_id=trial.ids[2], track_trace_table)
  trial.J.number <- order_egg_num(order_id=trial.ids[3], track_trace_table)
  trial.J.weight <- order_egg_weight(order_id=trial.ids[3], track_trace_table)
  
  # A. actual packed egg avg weight by order (each order has a grade)
  L.avg <- trial.L.weight / trial.L.number
  XL.avg <- trial.XL.weight / trial.XL.number
  J.avg <- trial.J.weight / trial.J.number
  
  # B. recorded avg weight of a supplied grade (values reported here only subject to measurement errors)
  if(accum_dist){
    hist_grade <- fread(grade_path)
    hist_grade <- hist_grade[SupplyId==supplyId & (Description %in% grades) & (Eggs != 0), .(Description, Eggs, Weight)]
    hist_grade <- hist_grade[, .(Description, Eggs, Weight = Weight/10, grade_avg = Weight/10/Eggs)]  # sample mean of each bin
  }
  else if(trial_num %in% c(4,5,6)){
    hist_grade <- snapshots_dist_summary()[[paste0("snapshot",trial_num-3)]][["dist"]]
  }
  else{
    print("accum_dist=FALSE doesn't support the input trial number")
    return(NULL)
  }
  
  large_num <- hist_grade[Description=="55",Eggs] + hist_grade[Description=="Large",Eggs]
  large_weight <- hist_grade[Description=="55",Weight] + hist_grade[Description=="Large",Weight]
  large_avg <- large_weight / large_num  # sample mean of L and 55 bins combined
  
  large_grade <- data.table(Description = "L.combined", Eggs = large_num, Weight = large_weight, grade_avg=large_avg)
  hist_grade <- rbindlist(list(hist_grade, large_grade))  # integrate into hist_grade table
  
  # B. modeled avg weight of a supplied grade - use repeated sampling to estimate s.e. of sample mean
  # (sample means reported here only subject to the randomness of a normal distribution)
  J.mean.dist <- grade_mean_dist('Jumbo', supplyId, times)
  XL.mean.dist <- grade_mean_dist('XLarge', supplyId, times)
  L.mean.dist <- grade_mean_dist('Large', supplyId, times)  # 55 & L combined
  SJ.mean.dist <- grade_mean_dist('SJumbo', supplyId, times)
  
  # fit normal parameters using gram distribution (here we don't consider measurement errors of mu and var)
  normal_fit <- fit_normal_dist(supplyId)
  mu.overall <- normal_fit[['mu.overall']]
  var.overall <- normal_fit[['var.overall']]
  
  return(list(large_order=list(eggs=trial.L.number, weight=trial.L.weight, mean=L.avg),
         xlarge_order=list(eggs=trial.XL.number, weight=trial.XL.weight, mean=XL.avg),
         jumbo_order=list(eggs=trial.J.number, weight=trial.J.weight, mean=J.avg),
         supply_record=hist_grade,
         supply_simu_L=L.mean.dist, supply_simu_XL=XL.mean.dist,
         supply_simu_J=J.mean.dist, supply_simu_SJ=SJ.mean.dist,
         supply_dist_est=list(mu=mu.overall, var=var.overall)))
}

## _____________________________________________________________________________
## Downgrade Analysis and Control
grade2entry <- function(grade_name){
  # grade_name: select one from c('Large', 'XLarge', 'Jumbo', 'SJumbo')
  # change the grade name to entry name to access simulated average weight of that grade in <trial_stats>
  if(grade_name=="SJumbo"){
    return("supply_simu_SJ")
  }
  else if(grade_name=="Jumbo"){
    return("supply_simu_J")
  }
  else if(grade_name=="XLarge"){
    return("supply_simu_XL")
  }
  else if(grade_name=="Large"){
    return("supply_simu_L")
  }
  else{
    print("input grade invalid")
    return(NULL)
  }
}

grade2entry.2 <- function(grade_name){
  # grade_name: select one from c('Large', 'XLarge', 'Jumbo')
  # change the grade name to entry name to access the real order info in <trial_stats>
  if(grade_name=="Jumbo"){
    return("jumbo_order")
  }
  else if(grade_name=="XLarge"){
    return("xlarge_order")
  }
  else if(grade_name=="Large"){
    return("large_order")
  }
  else{
    print("input grade invalid")
    return(NULL)
  }
}

update_egg_numbers <- function(grade_name, removed_num, removed_num_err, egg_number_list){
  # grade_name: select one from c('Large', 'XLarge', 'Jumbo', 'SJumbo')
  # removed_num: egg number to be removed from the original list of numbers
  # removed_num_err: error to propagate
  # egg_number_list: the original list of egg numbers
  egg_number <- egg_number_list[[grade_name]]
  new_number <- err_sum(1,egg_number[['n_hat']],egg_number[['n_se']], 
                        1,removed_num,removed_num_err,add=FALSE)
  egg_number_list[[grade_name]][['n_hat']] <- new_number[['z']]
  egg_number_list[[grade_name]][['n_se']] <- new_number[['dz']]
  
  return(egg_number_list)
}

grade_avg_stats <- function(grade_name, trial_stats){
  # given a grade name and a trial summary stats (returned by trial_stats_summary() function)
  # return the average egg weight and error in that grade
  if(grade_name=='effective_grade'){
    grade_avg <- trial_stats[['effective_grade']][['mu_hat']]
    grade_err <- trial_stats[['effective_grade']][['mu_se']]
  }
  else if(grade_name %in% c('Large', 'XLarge', 'Jumbo', 'SJumbo')){
    grade_avg <- trial_stats[[grade2entry(grade_name)]][['mu_hat']]
    grade_err <- trial_stats[[grade2entry(grade_name)]][['mu_se']]
  }
  else{
    print("input grade invalid")
    return(NULL)
  }
  return(list(egg_avg=grade_avg, avg_err=grade_err))
}

packed_egg_nums <- function(ordered_grade, downgraded_grade, packed_num, packed_weight, packed_weight.se, trial_stats){
  # ordered_grade: egg grade to be packed in the order; select one from c('Large', 'XLarge', 'Jumbo')
  # downgraded_grade: egg grade accepted to downgrade to packed_grade; select one from c('XLarge', 'Jumbo', 'SJumbo')
  # packed_num: total number of eggs packed in this order
  # packed_weight: total weight of eggs packed in this order (in gram)
  # packed_weight.se: measurement error of packed_weight (in gram)
  # trial_stats: the summary stats of the trial
  # return packed egg numbers of the ordered grade and the downgraded grade and their errors
  
  ordered_avg <- grade_avg_stats(ordered_grade, trial_stats)[['egg_avg']]
  ordered_avg_err <- grade_avg_stats(ordered_grade, trial_stats)[['avg_err']]
  
  downgraded_avg <- grade_avg_stats(downgraded_grade, trial_stats)[['egg_avg']]
  downgraded_avg_err <- grade_avg_stats(downgraded_grade, trial_stats)[['avg_err']]
  
  if(downgraded_avg <= ordered_avg){
    print("downgrade grade cannot be equal or lower to the ordered grade")
    return(NULL)
  }
  else{
    weight_diff <- err_sum(1,downgraded_avg,downgraded_avg_err,
                           1,ordered_avg,ordered_avg_err, add=FALSE)
    avg_weight_diff <- weight_diff[['z']]
    avg_weight_diff.se <- weight_diff[['dz']]
  }
  
  # refer to the formulas to calculate the packed egg numbers (solving two liner equations)
  downgraded_numerator <- err_sum(1,packed_weight,packed_weight.se,
                                  packed_num,ordered_avg,ordered_avg_err, add=FALSE)
  downgraded_num <- err_product(downgraded_numerator[['z']],downgraded_numerator[['dz']],
                                avg_weight_diff,avg_weight_diff.se, product=FALSE)    # e.g. n_sj- in my notes
  
  order_grade_numerator <- err_sum(packed_num,downgraded_avg,downgraded_avg_err,
                                   1,packed_weight,packed_weight.se, add=FALSE)
  order_grade_num <- err_product(order_grade_numerator[['z']],order_grade_numerator[['dz']],
                                 avg_weight_diff,avg_weight_diff.se, product=FALSE)   # e.g. n_j in my notes
  
  return(list(downgraded_num=downgraded_num, order_grade_num=order_grade_num))
}

accept_no_downgrade <- function(ordered_grade, packed_weight, packed_weight.se, trial_stats, egg_number_list){
  # The function to control the analysis if no grade is accepted to downgrade to the ordered_grade
  # simply calculate how many eggs of the current ordered_grade need to be packed
  packed_curr_grade_num <- err_product(packed_weight, packed_weight.se, trial_stats[[grade2entry(ordered_grade)]][['mu_hat']],
                                       trial_stats[[grade2entry(ordered_grade)]][['mu_se']], product=FALSE)
  
  # still record same-grade packing ratio. can't assume 100% packing rate
  same_grade_ratio <- err_product(packed_curr_grade_num[['z']], packed_curr_grade_num[['dz']],
                                 egg_number_list[[ordered_grade]][['n_hat']], egg_number_list[[ordered_grade]][['n_se']], product=FALSE)
  # calculate same-grade ratio before updated egg numbers
  egg_number_list <- update_egg_numbers(ordered_grade, packed_curr_grade_num[['z']], packed_curr_grade_num[['dz']], egg_number_list)
  
  return(list(same_grade_ratio=same_grade_ratio, egg_number_left=egg_number_list))
}

accept_one_downgrade <- function(ordered_grade, downgraded_grade, packed_num, packed_weight, packed_weight.se, trial_stats, egg_number_list){
  # The function to control the analysis if a grade only accept downgrade from one other single grade
  # parameters mostly same as function <packed_egg_nums>
  # egg_number_list: an list that includes egg numbers left in the pool and their errors
  check_has_downgrade(ordered_grade, trial_stats)
  
  # (1) - packed egg number (including ordered_grade and downgraded downgraded_grade egg number) and error
  total_num_to_pack <- packed_egg_nums(ordered_grade, downgraded_grade, packed_num, packed_weight, packed_weight.se, trial_stats)
  downgraded_num <- total_num_to_pack[['downgraded_num']]  # n_sj-: number of downgraded SJ eggs
  packed_curr_grade_num <- total_num_to_pack[['order_grade_num']]  # n_j: number of packed J eggs
  
  # (2) - based on ratio of J and SJ eggs (mixture ratio), estimate # of generated SJ eggs and its downgrade ratio
  # calculate downgrade ratio and same-grade ratio before updated egg numbers
  downgrade_ratio <- err_product(downgraded_num[['z']], downgraded_num[['dz']], 
                                 egg_number_list[[downgraded_grade]][['n_hat']], egg_number_list[[downgraded_grade]][['n_se']], product=FALSE)
  
  same_grade_ratio <- err_product(packed_curr_grade_num[['z']], packed_curr_grade_num[['dz']],
                                 egg_number_list[[ordered_grade]][['n_hat']], egg_number_list[[ordered_grade]][['n_se']], product=FALSE)
  
  egg_number_list <- update_egg_numbers(downgraded_grade, downgraded_num[['z']], downgraded_num[['dz']], egg_number_list)
  egg_number_list <- update_egg_numbers(ordered_grade, packed_curr_grade_num[['z']], packed_curr_grade_num[['dz']], egg_number_list)
  
  return(list(downgrade_ratio=downgrade_ratio, same_grade_ratio=same_grade_ratio, egg_number_left=egg_number_list))
}

accept_two_downgrades <- function(ordered_grade, downgraded_grades, packed_num, packed_weight, packed_weight.se, trial_stats, egg_number_list, 
                                  mode='probability-ratio', ratio=NULL){
  # The function to control the analysis if a grade accept downgrade from two other grades (see notes for algorithm details)
  # parameters exactly same as function <accept_one_downgrade>,
  # except that the downgraded grade is now an effective grade
  # downgraded_grades: a vector of two grades selected from c('XLarge', 'Jumbo', 'SJumbo'). Any combination accepted
  # mode: one of c('probability-ratio','exploration') - see function <add_effective_avg> for details
  # ratio: N2/N1 - treat as a constant
  if(length(downgraded_grades) != 2){
    print("downgrade grade number must be two")
    return(NULL)
  }
  check_has_downgrade(ordered_grade, trial_stats)
  
  # add effective grade average
  trial_stats <- add_effective_avg(downgraded_grades, egg_number_list, trial_stats, mode, ratio)
  
  # packed egg numbers
  total_num_to_pack <- packed_egg_nums(ordered_grade, 'effective_grade', packed_num, packed_weight, packed_weight.se, trial_stats)
  effective_num <- total_num_to_pack[['downgraded_num']]           # number of the combined downgraded eggs
  packed_curr_grade_num <- total_num_to_pack[['order_grade_num']]  # number of the eggs that are packed in their own order
  
  # split the combined downgraded egg numbers
  if(mode=='probability-ratio'){
    split_results <- split_comb_downgrade(effective_num, egg_number_list, 
                                          ratio=two_grades_ratio(downgraded_grades[2], downgraded_grades[1], egg_number_list)[['z']])
  }
  else if(mode=='exploration' & is.numeric(ratio)){
    split_results <- split_comb_downgrade(effective_num, egg_number_list, ratio)
  }
  else{
    print("Invalid N2/N1 ratio mode or ratio data type")
    return(NULL)
  }
  n1 <- split_results[['n1_result']][['z']]
  n1.err <- split_results[['n1_result']][['dz']]
  n2 <- split_results[['n2_result']][['z']]
  n2.err <- split_results[['n2_result']][['dz']]
  
  # N1 + N2
  left_num.1 <- egg_number_list[[downgraded_grades[1]]][['n_hat']]
  left_num.1.err <- egg_number_list[[downgraded_grades[1]]][['n_se']]
  left_num.2 <- egg_number_list[[downgraded_grades[2]]][['n_hat']]
  left_num.2.err <- egg_number_list[[downgraded_grades[2]]][['n_se']]
  sum12 <- err_sum(1,left_num.1,left_num.1.err, 1,left_num.2,left_num.2.err, add=TRUE)
  
  # calculate downgrade ratio and same-grade ratio
  downgrade_ratio <- err_product(effective_num[['z']], effective_num[['dz']], 
                                 sum12[['z']], sum12[['dz']], product=FALSE)
  
  same_grade_ratio <- err_product(packed_curr_grade_num[['z']], packed_curr_grade_num[['dz']],
                                  egg_number_list[[ordered_grade]][['n_hat']], egg_number_list[[ordered_grade]][['n_se']], product=FALSE)
  
  # updated egg numbers
  egg_number_list <- update_egg_numbers(downgraded_grades[1], n1, n1.err, egg_number_list)
  egg_number_list <- update_egg_numbers(downgraded_grades[2], n2, n2.err, egg_number_list)
  egg_number_list <- update_egg_numbers(ordered_grade, packed_curr_grade_num[['z']], packed_curr_grade_num[['dz']], egg_number_list)
  
  # remove effective grade average
  trial_stats <- remove_effective_avg(trial_stats)
  
  return(list(downgrade_ratio=downgrade_ratio, same_grade_ratio=same_grade_ratio, egg_number_left=egg_number_list))
}

no_downgrade_wrapper <- function(ordered_grade, trial_stats, egg_number_left){
  # if one order doesn't accept any downgrade grades
  order_num <- trial_stats[[grade2entry.2(ordered_grade)]][['eggs']]
  order_weight <- trial_stats[[grade2entry.2(ordered_grade)]][['weight']]
  order_weight.err <- sum_weight_err(order_num)
  
  order_update <- accept_no_downgrade(ordered_grade, order_weight, order_weight.err, trial_stats, egg_number_left)
  same_grade_ratio <- order_update[['same_grade_ratio']]
  egg_number_left <- order_update[['egg_number_left']]
  
  return(list(downgrade_ratio=NULL, same_grade_ratio=same_grade_ratio, egg_number_left=egg_number_left))
}

one_downgrade_wrapper <- function(ordered_grade, downgraded_grade, trial_stats, egg_number_left){
  # if one order only accept one downgrade grades
  order_num <- trial_stats[[grade2entry.2(ordered_grade)]][['eggs']]
  order_weight <- trial_stats[[grade2entry.2(ordered_grade)]][['weight']]
  order_weight.err <- sum_weight_err(order_num)
  
  order_update <- accept_one_downgrade(ordered_grade, downgraded_grade, order_num, order_weight, order_weight.err, trial_stats, egg_number_left)
  downgrade_ratio <- order_update[['downgrade_ratio']]
  same_grade_ratio <- order_update[['same_grade_ratio']]
  egg_number_left <- order_update[['egg_number_left']]
  
  return(list(downgrade_ratio=downgrade_ratio, same_grade_ratio=same_grade_ratio, egg_number_left=egg_number_left))
}

two_downgrade_wrapper <- function(ordered_grade, downgraded_grades, trial_stats, egg_number_left, mode='probability-ratio', ratio=NULL){
  # if one order accept two downgrade grades
  # mode: one of c('probability-ratio','exploration') - see function <add_effective_avg> for details
  # ratio: N2/N1 - treat as a constant
  order_num <- trial_stats[[grade2entry.2(ordered_grade)]][['eggs']]
  order_weight <- trial_stats[[grade2entry.2(ordered_grade)]][['weight']]
  order_weight.err <- sum_weight_err(order_num)
  
  order_update <- accept_two_downgrades(ordered_grade, downgraded_grades, order_num, order_weight, order_weight.err, trial_stats, egg_number_left, mode, ratio)
  downgrade_ratio <- order_update[['downgrade_ratio']]
  same_grade_ratio <- order_update[['same_grade_ratio']]
  egg_number_left <- order_update[['egg_number_left']]
  
  return(list(downgrade_ratio=downgrade_ratio, same_grade_ratio=same_grade_ratio, egg_number_left=egg_number_left))
}

trial.run <- function(total_ordered_num, mu.overall, var.overall, trial_stats, downgrade_rules, mode='probability-ratio', ratios=c(0,0)){
  ## a simple wrapper function to run factory trial analysis for one time
  ## there are 3 orders in each trial, follow the J-XL-L order
  ## downgrade_rules: the accepted downgrade rules for 3 orders. Three rules available: none, one, or two allowed downgrades.
  ## e.g. trial 3, downgrade_rules <- list(Jumbo=c("SJumbo"), XLarge=c(), Large=c("XLarge", "Jumbo")), i.e the Jumbo order accept one,
  ## the XLarge order accept none, and the Large order accept two downgraded grades.
  ## mode: one of c('probability-ratio','exploration') - see function <add_effective_avg> for details - available when there are two downgraded grades
  ## ratios: if c(ratio1, ratio2), both XLarge and Large orders accept two downgraded grades; - use if mode=="exploration"
  ## if c(0, ratio), only Large order accepts two downgraded grades;
  ## if c(ratio, 0), only XLarge order accepts two downgraded grades.
  
  sjumbo_num_dist <- egg_grade_initial_number(total_ordered_num, "SJumbo", mu.overall, var.overall)
  jumbo_num_dist <- egg_grade_initial_number(total_ordered_num, "Jumbo", mu.overall, var.overall)
  xlarge_num_dist <- egg_grade_initial_number(total_ordered_num, "XLarge", mu.overall, var.overall)
  large_num_dist <- egg_grade_initial_number(total_ordered_num, "Large", mu.overall, var.overall)
  
  egg_number_left <- list(SJumbo=sjumbo_num_dist, Jumbo=jumbo_num_dist, XLarge=xlarge_num_dist, Large=large_num_dist)
  
  # Jumbo order
  if(length(downgrade_rules[['Jumbo']])==1){
    jumbo_order_update <- one_downgrade_wrapper('Jumbo', downgrade_rules[['Jumbo']], trial_stats, egg_number_left)
  }
  else if(length(downgrade_rules[['Jumbo']])==0){
    jumbo_order_update <- one_downgrade_wrapper('Jumbo', downgrade_rules[['Jumbo']], trial_stats, egg_number_left)
  }
  egg_number_left <- jumbo_order_update[['egg_number_left']]
  
  # XLarge order
  if(length(downgrade_rules[['XLarge']])==2 & mode=='probability-ratio'){
    xlarge_order_update <- two_downgrade_wrapper('XLarge', downgrade_rules[['XLarge']], trial_stats, egg_number_left, mode, ratio=NULL)
  }
  else if(length(downgrade_rules[['XLarge']])==2 & mode=='exploration'){
    if(ratios[1] == 0){
      xlarge_order_update <- two_downgrade_wrapper('XLarge', downgrade_rules[['XLarge']], trial_stats, egg_number_left, mode, ratio=NULL)  # convert 0 to NULL
    }
    else{
      xlarge_order_update <- two_downgrade_wrapper('XLarge', downgrade_rules[['XLarge']], trial_stats, egg_number_left, mode, ratio=ratios[1])
    }
  }
  else if(length(downgrade_rules[['XLarge']])==1){
    xlarge_order_update <- one_downgrade_wrapper('XLarge', downgrade_rules[['XLarge']], trial_stats, egg_number_left)
  }
  else if(length(downgrade_rules[['XLarge']])==0){
    xlarge_order_update <- no_downgrade_wrapper('XLarge', trial_stats, egg_number_left)
  }
  egg_number_left <- xlarge_order_update[['egg_number_left']]
  
  # Large order
  if(length(downgrade_rules[['Large']])==2 & mode=='probability-ratio'){
    large_order_update <- two_downgrade_wrapper('Large', downgrade_rules[['Large']], trial_stats, egg_number_left, mode, ratio=NULL)
  }
  else if(length(downgrade_rules[['Large']])==2 & mode=='exploration'){
    if(ratios[2] == 0){
      large_order_update <- two_downgrade_wrapper('Large', downgrade_rules[['Large']], trial_stats, egg_number_left, mode, ratio=NULL)  # convert 0 to NULL
    }
    else{
      large_order_update <- two_downgrade_wrapper('Large', downgrade_rules[['Large']], trial_stats, egg_number_left, mode, ratio=ratios[2])
    }
  }
  else if(length(downgrade_rules[['Large']])==1){
    large_order_update <- one_downgrade_wrapper('Large', downgrade_rules[['Large']], trial_stats, egg_number_left)
  }
  else if(length(downgrade_rules[['Large']])==0){
    large_order_update <- no_downgrade_wrapper('Large', trial_stats, egg_number_left)
  }
  egg_number_left <- large_order_update[['egg_number_left']]
  
  return(list(egg_number_left=egg_number_left, 
              sj_downgrade_ratio=jumbo_order_update[['downgrade_ratio']], j_same_grade_ratio=jumbo_order_update[['same_grade_ratio']],
              j_downgrade_ratio=xlarge_order_update[['downgrade_ratio']], xl_same_grade_ratio=xlarge_order_update[['same_grade_ratio']],
              xl_downgrade_ratio=large_order_update[['downgrade_ratio']], l_same_grade_ratio=large_order_update[['same_grade_ratio']]))
}

check_has_downgrade <- function(ordered_grade, trial_stats){
  # check if a grade has accepted downgrades in reality
  grade_avg <- trial_stats[[grade2entry(ordered_grade)]][['mu_hat']]
  packed_avg <- trial_stats[[paste0(tolower(ordered_grade),'_order')]][['mean']]
  if(grade_avg>=packed_avg){
    print("*********")
    print(paste0("There shouldn't have any downgrades to ",ordered_grade," grade."))
    print("Please carefully check and adjust downgrade rules and make change accordingly.")
    print("*********")
  }
}

split_comb_downgrade <- function(effective_num, egg_number_list, ratio){
  # split the effective downgraded number to numbers in two grades (n1- and n2-)
  # ratio: N2/N1 - treat as a constant
  
  # n1-
  n1_result <- list(z=effective_num[['z']]/(1+ratio), dz=effective_num[['dz']]/(1+ratio))
  
  # n2-
  n2_result <- err_sum(1,effective_num[['z']],effective_num[['dz']], 1,n1_result[['z']],n1_result[['dz']], add=FALSE)
  
  return(list(n1_result=n1_result, n2_result=n2_result))
}

two_grades_ratio <- function(grade.1, grade.2, egg_number_list){
  # calculate N1/N2 and its error, where N1 and N2 are the remaining egg numbers of grade 1 and 2
  N1 <- egg_number_list[[grade.1]][['n_hat']]
  N1.err <- egg_number_list[[grade.1]][['n_se']]
  N2 <- egg_number_list[[grade.2]][['n_hat']]
  N2.err <- egg_number_list[[grade.2]][['n_se']]
  return(err_product(N1, N1.err, N2, N2.err, product=FALSE))
}

add_effective_avg <- function(downgraded_grades, egg_number_list, trial_stats, mode='probability-ratio', ratio=NULL){
  # downgraded_grades: a vector of two grades selected from c('XLarge', 'Jumbo', 'SJumbo'). Any combination accepted
  # Add the grade avg and error of the effective grade combined by the two downgraded grades to <trial_stats>.
  # mode: one of c('probability-ratio','exploration')
  # <probability-ratio>: n2-/n2-=P2/P1=N2/N1 - no need to consider error of the ratio
  # <exploration>: any constant ratio of N2/N1 accepted - no need to consider error of the ratio
  # ratio: N2/N1 - treat as a constant
  
  if(length(downgraded_grades)==2){
    grade.1 <- downgraded_grades[1]
    grade.2 <- downgraded_grades[2]
  }
  else{
    print("To add effective average must only have two downgraded grades")
    return(NULL)
  }
  
  # get w1 and w2
  grade.1.avg <- trial_stats[[grade2entry(grade.1)]][['mu_hat']]
  grade.1.avg.err <- trial_stats[[grade2entry(grade.1)]][['mu_se']]
  grade.2.avg <- trial_stats[[grade2entry(grade.2)]][['mu_hat']]
  grade.2.avg.err <- trial_stats[[grade2entry(grade.2)]][['mu_se']]
  
  # N2/N1
  if(mode=='probability-ratio'){
    ratio <- two_grades_ratio(grade.2, grade.1, egg_number_list)[['z']]
  }
  else if(mode=='exploration' & !is.numeric(ratio)){
    print("In the ratio exploration mode the input <ratio> must be assigned a valid number")
    return(NULL)
  }
  else if(mode!='exploration'){
    print("Invalid mode when assigning N2/N1 ratio")
    return(NULL)
  }
  
  # (N2/N1)*w2
  right <- ratio*grade.2.avg
  right.err <- ratio*grade.2.avg.err
  
  # w1 + (N2/N1)*w2
  eff.avg <- err_sum(1,grade.1.avg,grade.1.avg.err, 1,right,right.err, add=TRUE)[['z']]
  eff.avg.err <- err_sum(1,grade.1.avg,grade.1.avg.err, 1,right,right.err, add=TRUE)[['dz']]
  
  # update
  trial_stats[['effective_grade']] <- list(mu_hat=eff.avg, mu_se=eff.avg.err)
  return(trial_stats)
}

remove_effective_avg <- function(trial_stats){
  # remove the temporary effective grade from <trial_stats>
  trial_stats[['effective_grade']] <- NULL
  return(trial_stats)
}

check_insufficient_eggs <- function(egg_number_list){
  # check which grade has the most insufficient supplied egg number
  left_num_by_grade <- c()
  for(grade in ordered.grades){
    egg_left <- egg_number_list[['egg_number_left']][[grade]][['n_hat']]
    egg_left.error <- egg_number_list[['egg_number_left']][[grade]][['n_se']]
    left_num_by_grade <- append(left_num_by_grade, egg_left+egg_left.error)
  }
  if(min(left_num_by_grade)<0){
    return(ordered.grades[which.min(left_num_by_grade)])
  }
  # return NULL if no insufficient grade
}

## _____________________________________________________________________________
## Trial Reproduction and Subsequent Analysis
trial.reproduce <- function(total_ordered_nums, mu.overall, var.overall, downgrade_rules, 
                            downgrade_ratios, downgrade_ratios.err, same_grade_ratio, same_grade_ratio.err){
  # total_ordered_nums: total ordered egg numbers of three orders - e.g. list(Jumbo=9524, XLarge=37733, Large=16105)
  # randomly sample eggs from a normal distribution with parameters mu.overall and var.overall
  # downgrade_rules: e.g. list(Jumbo=c("SJumbo"), XLarge=c(), Large=c("XLarge", "Jumbo"))
  # downgrade_ratios & same_grade_ratio: must match downgrade_rules, e.g. list(Jumbo=0.42, XLarge=0.0, Large=0.18)
  # downgrade_ratios.err & same_grade_ratio.err: must match downgrade_ratios, e.g. list(Jumbo=0.04, XLarge=0.0, Large=0.02)
  # i.e. 0.42 probability for a SJumbo egg to downgrade when passing the Jumbo production lanes
  
  jumbo_num <- total_ordered_nums[['Jumbo']]
  xlarge_num <- total_ordered_nums[['XLarge']]
  large_num <- total_ordered_nums[['Large']]
  EOL_num <- 0  # end-of-lane
  
  # pack actual eggs inside orders
  jumbo_eggs <- xlarge_eggs <- large_eggs <- EOL_eggs <- c()
  
  while(TRUE){
    if(jumbo_num==0 & xlarge_num==0 & large_num==0){
      break
    }
    egg <- rnorm(1, mu.overall, sqrt(var.overall))
    # an egg passing through Jumbo, XLarge, and Large lanes in order
    if(jumbo_num>0){
      results <- passing_production_lane("Jumbo", jumbo_num, egg, downgrade_rules, 
                                         downgrade_ratios, downgrade_ratios.err, same_grade_ratio, same_grade_ratio.err)
      jumbo_num <- results[1]
      packed <- results[2]
      if(packed){
        jumbo_eggs <- append(jumbo_eggs, egg)
        next
      }
    }
    if(xlarge_num>0){
      results <- passing_production_lane("XLarge", xlarge_num, egg, downgrade_rules, 
                                         downgrade_ratios, downgrade_ratios.err, same_grade_ratio, same_grade_ratio.err)
      xlarge_num <- results[1]
      packed <- results[2]
      if(packed){
        xlarge_eggs <- append(xlarge_eggs, egg)
        next
      }
    }
    if(large_num>0){
      results <- passing_production_lane("Large", large_num, egg, downgrade_rules,
                                         downgrade_ratios, downgrade_ratios.err, same_grade_ratio, same_grade_ratio.err)
      large_num <- results[1]
      packed <- results[2]
      if(packed){
        large_eggs <- append(large_eggs, egg)
        next
      }
    }
    EOL_num <- EOL_num + 1
    EOL_eggs <- append(EOL_eggs, egg)
  }
  return(list(Jumbo=jumbo_eggs, XLarge=xlarge_eggs, Large=large_eggs, EOL=EOL_eggs))
}

passing_production_lane <- function(order_grade, left_num, egg.weight, downgrade_rules, 
                                    downgrade_ratios, downgrade_ratios.err, same_grade_ratio, same_grade_ratio.err){
  # represents the production lane of any of c("Jumbo", "XLarge", "Large") ordered grade
  # left_num: how many ordered number of this grade still required to complete
  accepted.grades <- downgrade_rules[[order_grade]]
  downgrade.ratio <- downgrade_ratios[[order_grade]]   # probability of a higher grade egg to downgrade
  downgrade.ratio.err <- downgrade_ratios.err[[order_grade]]
  same_grade.ratio <- same_grade_ratio[[order_grade]]  # probability of a current grade egg to get packed
  same_grade.ratio.err <- same_grade_ratio.err[[order_grade]]
  egg.grade <- weight2grade(egg.weight)
  
  packed <- FALSE
  if(egg.grade==order_grade){
    # sample same grade ratio by considering its error - assume normal distribution
    same_grade.ratio <- rnorm(1, mean=same_grade.ratio, sd=same_grade.ratio.err)
    if(runif(1)<=same_grade.ratio){
      left_num <- left_num - 1
      packed <- TRUE
    }
  }
  else if(egg.grade %in% accepted.grades){
    # sample downgrade ratio only if downgrade accepted
    downgrade.ratio <- rnorm(1, mean=downgrade.ratio, sd=downgrade.ratio.err)
    if(runif(1)<=downgrade.ratio){
      left_num <- left_num - 1
      packed <- TRUE
    }
  }
  return(c(left_num, packed))
}

weight2grade <- function(weight){
  # convert egg weight to grade
  if(weight>get_grade_bounds("SJumbo")[1] & weight<get_grade_bounds("SJumbo")[2]){
    return("SJumbo")
  }
  else if(weight>get_grade_bounds("Jumbo")[1] & weight<get_grade_bounds("Jumbo")[2]){
    return("Jumbo")
  }
  else if(weight>get_grade_bounds("XLarge")[1] & weight<get_grade_bounds("XLarge")[2]){
    return("XLarge")
  }
  else if(weight>get_grade_bounds("Large")[1] & weight<get_grade_bounds("Large")[2]){
    return("Large")
  }
  else{
    return("Others")
  }
}

trial.reproduce.wrapper <- function(trial_results, trial_stats, downgrade_rules){
  # a simple wrapper to run trial reproduction with minimum input parameters
  
  downgrade_ratios <- list(Jumbo=trial_results[['sj_downgrade_ratio']][['z']], 
                           XLarge=trial_results[['j_downgrade_ratio']][['z']], Large=trial_results[['xl_downgrade_ratio']][['z']])
  downgrade_ratios.err <- list(Jumbo=trial_results[['sj_downgrade_ratio']][['dz']], 
                               XLarge=trial_results[['j_downgrade_ratio']][['dz']], Large=trial_results[['xl_downgrade_ratio']][['dz']])
  
  same_grade_ratios <- list(Jumbo=trial_results[['j_same_grade_ratio']][['z']], 
                            XLarge=trial_results[['xl_same_grade_ratio']][['z']], Large=trial_results[['l_same_grade_ratio']][['z']])
  same_grade_ratios.err <- list(Jumbo=trial_results[['j_same_grade_ratio']][['dz']], 
                                XLarge=trial_results[['xl_same_grade_ratio']][['dz']], Large=trial_results[['l_same_grade_ratio']][['dz']])
  
  mu <- trial_stats[['supply_dist_est']][['mu']]
  var <- trial_stats[['supply_dist_est']][['var']]
  total_ordered_nums <- list(Jumbo=trial_stats[['jumbo_order']][['eggs']],
                             XLarge=trial_stats[['xlarge_order']][['eggs']], Large=trial_stats[['large_order']][['eggs']])
  
  packed_eggs <- trial.reproduce(total_ordered_nums, mu, var, downgrade_rules, 
                                 downgrade_ratios, downgrade_ratios.err, same_grade_ratios, same_grade_ratios.err)
}

get_simple_ratios <- function(trial_stats){
  # use trial stats to compute the simple giveaway and downgrade ratios
  # of the overall process and of each ordered grade (J, XL, L)
  giveaway.numerator <- downgrade.numerator <- denominator <- 0     # to compute the overall downgrade and giveaway ratio
  giveaway.numerator.err <- downgrade.numerator.err <- 0
  stats <- list()
  
  for(grade in ordered.grades){
    stats[[grade]][['lower_bound']] <- get_grade_bounds(grade)[1]
    stats[[grade]][['packed_avg']] <- trial_stats[[grade2entry.2(grade)]][['mean']]
    stats[[grade]][['packed_num']] <- trial_stats[[grade2entry.2(grade)]][['eggs']]
    stats[[grade]][['supply_grade_avg']] <- trial_stats[[grade2entry(grade)]][['mu_hat']]
    
    packed_avg <- stats[[grade]][['packed_avg']]  # avg weight packed into an ordered grade
    supply_grade_avg <- stats[[grade]][['supply_grade_avg']]  # avg weight of a supplied grade
    packed_avg.err <- sum_weight_err(stats[[grade]][['packed_num']]) / stats[[grade]][['packed_num']]  # measurement error
    supply_grade_avg.err <- trial_stats[[grade2entry(grade)]][['mu_se']]  # modeled error
    
    if(packed_avg>=supply_grade_avg){
      giveaway_weight <- packed_avg - stats[[grade]][['lower_bound']]
      giveaway_weight.err <- err_sum(1,packed_avg,packed_avg.err, 1,stats[[grade]][['lower_bound']],0.01, add=FALSE)[['dz']]
      downgrade_weight <- packed_avg - supply_grade_avg
      downgrade_weight.err <- err_sum(1,packed_avg,packed_avg.err, 1,supply_grade_avg,supply_grade_avg.err, add=FALSE)[['dz']]
      stats[[grade]][['downgrade_ratio']] <- downgrade_weight / stats[[grade]][['lower_bound']]
      stats[[grade]][['downgrade_ratio_err']] <- downgrade_weight.err / stats[[grade]][['lower_bound']]
      stats[[grade]][['giveaway_ratio']] <- giveaway_weight / stats[[grade]][['lower_bound']]
      stats[[grade]][['giveaway_ratio_err']] <- giveaway_weight.err / stats[[grade]][['lower_bound']]
    }
    else{  # don't allow "upgrade"; set downgrade to 0 if negative
      downgrade_weight <- 0
      giveaway_weight <- supply_grade_avg - stats[[grade]][['lower_bound']]
      giveaway_weight.err <- supply_grade_avg.err
      stats[[grade]][['downgrade_ratio']] <- 0
      stats[[grade]][['downgrade_ratio_err']] <- 0
      stats[[grade]][['giveaway_ratio']] <- giveaway_weight / stats[[grade]][['lower_bound']]
      stats[[grade]][['giveaway_ratio_err']] <- giveaway_weight.err / stats[[grade]][['lower_bound']]
    }
    giveaway.numerator <- err_sum(1,giveaway.numerator,giveaway.numerator.err, 
                                  stats[[grade]][['packed_num']],giveaway_weight,giveaway_weight.err, add=TRUE)[['z']]
    giveaway.numerator.err <- err_sum(1,giveaway.numerator,giveaway.numerator.err, 
                                      stats[[grade]][['packed_num']],giveaway_weight,giveaway_weight.err, add=TRUE)[['dz']]
    downgrade.numerator <- err_sum(1,downgrade.numerator,downgrade.numerator.err, 
                                   stats[[grade]][['packed_num']],downgrade_weight,downgrade_weight.err, add=TRUE)[['z']]
    downgrade.numerator.err <- err_sum(1,downgrade.numerator,downgrade.numerator.err, 
                                       stats[[grade]][['packed_num']],downgrade_weight,downgrade_weight.err, add=TRUE)[['dz']]
    denominator <- denominator + stats[[grade]][['packed_num']] * stats[[grade]][['lower_bound']]
  }
  return(list(ratios_by_grade=stats,
              overall_giveaway.ratio=giveaway.numerator/denominator,
              overall_giveaway.ratio.err=giveaway.numerator.err/denominator,
              overall_downgrade.ratio=downgrade.numerator/denominator,
              overall_downgrade.ratio.err=downgrade.numerator.err/denominator))
}

get_simu_ratios <- function(packed_eggs, trial_stats){
  # for a simulation row data collection stored in <packed_eggs>,
  # we return giveaway and downgrade ratios of the overall process and of each ordered grade (J, XL, L)
  grades <- names(packed_eggs)
  giveaway.numerator <- downgrade.numerator <- denominator <- 0     # to compute the overall downgrade and giveaway ratio
  stats <- list()
  
  for(grade in grades[1:(length(grades)-1)]){
    stats[[grade]][['lower_bound']] <- get_grade_bounds(grade)[1]
    stats[[grade]][['packed_avg']] <- mean(packed_eggs[[grade]])
    stats[[grade]][['packed_num']] <- length(packed_eggs[[grade]])  # involve randomness
    supply_grade_avg <- trial_stats[[grade2entry(grade)]]
    stats[[grade]][['supply_grade_avg']] <- rnorm(1, mean=supply_grade_avg[['mu_hat']], sd=supply_grade_avg[['mu_se']])  # add randomness
    
    stats[[grade]][['downgrade_ratio']] <- (stats[[grade]][['packed_avg']] - stats[[grade]][['supply_grade_avg']]) / stats[[grade]][['lower_bound']]
    stats[[grade]][['giveaway_ratio']] <- (stats[[grade]][['packed_avg']] - stats[[grade]][['lower_bound']]) / stats[[grade]][['lower_bound']]
    
    giveaway.numerator <- giveaway.numerator + stats[[grade]][['packed_num']] * (stats[[grade]][['packed_avg']] - stats[[grade]][['lower_bound']])
    downgrade.numerator <- downgrade.numerator + stats[[grade]][['packed_num']] * (stats[[grade]][['packed_avg']] - stats[[grade]][['supply_grade_avg']])
    denominator <- denominator + stats[[grade]][['packed_num']] * stats[[grade]][['lower_bound']]
  }
  return(list(ratios_by_grade=stats,
              overall_giveaway.ratio=giveaway.numerator/denominator,
              overall_downgrade.ratio=downgrade.numerator/denominator))
}

ratio_samplings <- function(times, trial_results, trial_stats, downgrade_rules){
  # to run <get_simu_ratios> multiple times
  # return giveaway.ratio and downgrade.ratio in two vectors
  g.ratios <- d.ratios <- c()
  for(i in (1:times)){
    packed_eggs <- trial.reproduce.wrapper(trial_results, trial_stats, downgrade_rules)
    ratio_stats <- get_simu_ratios(packed_eggs, trial_stats)
    g.ratios <- append(g.ratios, ratio_stats[['overall_giveaway.ratio']])
    d.ratios <- append(d.ratios, ratio_stats[['overall_downgrade.ratio']])
  }
  return(list(giveaway.ratios=g.ratios, downgrade.ratios=d.ratios))
}
