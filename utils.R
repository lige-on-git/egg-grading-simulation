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
## Data Retrival and Exploration
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
  dx <- 0.1
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
  egg_number_list <- update_egg_numbers(ordered_grade, packed_curr_grade_num[['z']], packed_curr_grade_num[['dz']], egg_number_list)
  
  # still record same-grade packing ratio. can't assume 100% packing rate
  same_grade_ratio <- err_product(packed_curr_grade_num[['z']], packed_curr_grade_num[['dz']],
                                 egg_number_list[[ordered_grade]][['n_hat']], egg_number_list[[ordered_grade]][['n_se']], product=FALSE)
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

accept_two_downgrades <- function(ordered_grade, downgraded_grades, packed_num, packed_weight, packed_weight.se, trial_stats, egg_number_list){
  # The function to control the analysis if a grade accept downgrade from two other grades (see notes for algorithm details)
  # parameters exactly same as function <accept_one_downgrade>,
  # except that the downgraded grade is now an effective grade
  # downgraded_grades: a vector of two grades selected from c('XLarge', 'Jumbo', 'SJumbo'). Any combination accepted.
  if(length(downgraded_grades) != 2){
    print("downgrade grade number must be two")
    return(NULL)
  }
  check_has_downgrade(ordered_grade, trial_stats)
  
  # add effective grade average
  trial_stats <- add_effective_avg(downgraded_grades, egg_number_list, trial_stats)
  
  # packed egg numbers
  total_num_to_pack <- packed_egg_nums(ordered_grade, 'effective_grade', packed_num, packed_weight, packed_weight.se, trial_stats)
  effective_num <- total_num_to_pack[['downgraded_num']]  # number of combined downgraded eggs
  packed_curr_grade_num <- total_num_to_pack[['order_grade_num']]
  
  # split the combined downgraded egg numbers
  split_results <- split_comb_downgrade(downgraded_grades[1], downgraded_grades[2], effective_num, egg_number_list)
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

two_downgrade_wrapper <- function(ordered_grade, downgraded_grades, trial_stats, egg_number_left){
  # if one order accept two downgrade grades
  order_num <- trial_stats[[grade2entry.2(ordered_grade)]][['eggs']]
  order_weight <- trial_stats[[grade2entry.2(ordered_grade)]][['weight']]
  order_weight.err <- sum_weight_err(order_num)
  
  order_update <- accept_two_downgrades(ordered_grade, downgraded_grades, order_num, order_weight, order_weight.err, trial_stats, egg_number_left)
  downgrade_ratio <- order_update[['downgrade_ratio']]
  same_grade_ratio <- order_update[['same_grade_ratio']]
  egg_number_left <- order_update[['egg_number_left']]
  
  return(list(downgrade_ratio=downgrade_ratio, same_grade_ratio=same_grade_ratio, egg_number_left=egg_number_left))
}

trial.run <- function(total_ordered_num, mu.overall, var.overall, trial_stats, downgrade_rules){
  ## a simple wrapper function to run factory trial analysis for one time
  ## there are 3 orders in each trial, follow the J-XL-L order
  ## downgrade_rules: the accepted downgrade rules for 3 orders. Three rules available: none, one, or two allowed downgrades.
  ## e.g. trial 3, downgrade_rules <- list(Jumbo=c("SJumbo"), XLarge=c(), Large=c("XLarge", "Jumbo")), i.e the Jumbo order accept one,
  ## the XLarge order accept none, and the Large order accept two downgraded grades.
  
  sjumbo_num_dist <- egg_grade_initial_number(total_ordered_num, "SJumbo", mu.overall, var.overall)
  jumbo_num_dist <- egg_grade_initial_number(total_ordered_num, "Jumbo", mu.overall, var.overall)
  xlarge_num_dist <- egg_grade_initial_number(total_ordered_num, "XLarge", mu.overall, var.overall)
  large_num_dist <- egg_grade_initial_number(total_ordered_num, "Large", mu.overall, var.overall)
  
  egg_number_left <- list(SJumbo=sjumbo_num_dist, Jumbo=jumbo_num_dist, XLarge=xlarge_num_dist, Large=large_num_dist)
  
  # Jumbo order
  if(length(downgrade_rules[['Jumbo']])==2){
    jumbo_order_update <- two_downgrade_wrapper('Jumbo', downgrade_rules[['Jumbo']], trial_stats, egg_number_left)
  }
  else if(length(downgrade_rules[['Jumbo']])==1){
    jumbo_order_update <- one_downgrade_wrapper('Jumbo', downgrade_rules[['Jumbo']], trial_stats, egg_number_left)
  }
  else if(length(downgrade_rules[['Jumbo']])==0){
    jumbo_order_update <- one_downgrade_wrapper('Jumbo', downgrade_rules[['Jumbo']], trial_stats, egg_number_left)
  }
  egg_number_left <- jumbo_order_update[['egg_number_left']]
  
  # XLarge order
  if(length(downgrade_rules[['XLarge']])==2){
    xlarge_order_update <- two_downgrade_wrapper('XLarge', downgrade_rules[['XLarge']], trial_stats, egg_number_left)
  }
  else if(length(downgrade_rules[['XLarge']])==1){
    xlarge_order_update <- one_downgrade_wrapper('XLarge', downgrade_rules[['XLarge']], trial_stats, egg_number_left)
  }
  else if(length(downgrade_rules[['XLarge']])==0){
    xlarge_order_update <- no_downgrade_wrapper('XLarge', trial_stats, egg_number_left)
  }
  egg_number_left <- xlarge_order_update[['egg_number_left']]
  
  # Large order
  if(length(downgrade_rules[['Large']])==2){
    large_order_update <- two_downgrade_wrapper('Large', downgrade_rules[['Large']], trial_stats, egg_number_left)
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

split_comb_downgrade <- function(grade.1, grade.2, effective_num, egg_number_list){
  # split the effective downgraded number to numbers in two grades (n1- and n2-)
  # N2/N1
  ratio_results <- two_grades_ratio(grade.2, grade.1, egg_number_list)
  ratio <- ratio_results[['z']]
  ratio.err <- ratio_results[['dz']]
  
  # n1-
  n1_result <- err_product(effective_num[['z']], effective_num[['dz']], 1+ratio, ratio.err, product=FALSE)
  
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

add_effective_avg <- function(downgraded_grades, egg_number_list, trial_stats){
  # downgraded_grades: a vector of two grades selected from c('XLarge', 'Jumbo', 'SJumbo'). Any combination accepted
  # Add the grade avg and error of the effective grade combined by the two downgraded grades to <trial_stats>.
  # Algorithm details see my notebook.
  grade.1 <- downgraded_grades[1]
  grade.2 <- downgraded_grades[2]
  
  # get w1 and w2
  grade.1.avg <- trial_stats[[grade2entry(grade.1)]][['mu_hat']]
  grade.1.avg.err <- trial_stats[[grade2entry(grade.1)]][['mu_se']]
  grade.2.avg <- trial_stats[[grade2entry(grade.2)]][['mu_hat']]
  grade.2.avg.err <- trial_stats[[grade2entry(grade.2)]][['mu_se']]
  
  # N2/N1
  ratio_results <- two_grades_ratio(grade.2, grade.1, egg_number_list)
  ratio <- ratio_results[['z']]
  ratio.err <- ratio_results[['dz']]
  
  # (N2/N1)*w2
  right_result <- err_product(ratio, ratio.err, grade.2.avg, grade.2.avg.err, product=TRUE)
  right <- right_result[['z']]
  right.err <- right_result[['dz']]
  
  # w1 + (N2/N1)*w2
  eff.avg_result <- err_sum(1,grade.1.avg,grade.1.avg.err, 1,right,right.err, add=TRUE)
  eff.avg <- eff.avg_result[['z']]
  eff.avg.err <- eff.avg_result[['dz']]
  
  # update
  trial_stats[['effective_grade']] <- list(mu_hat=eff.avg, mu_se=eff.avg.err)
  return(trial_stats)
}

remove_effective_avg <- function(trial_stats){
  # remove the temporary effective grade from <trial_stats>
  trial_stats[['effective_grade']] <- NULL
  return(trial_stats)
}

## _____________________________________________________________________________
## Trial Reproduction and Subsequent Analysis
