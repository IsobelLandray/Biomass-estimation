
#THE DATASET AND PLOTS OUTPUTTED ALL SAY _PLUSLIMIT AS I HAVE NOT ONLY PUT A PRIOR ON THE PIVOT BUT ALSO 
#RESTRICTED IT TO BE AT LEAST 18.5 AS THIS HELPED WITH SOME OF THE CONVERGENCE - WENT FROM 44 GROUPS
#NOT CONVERGING TO 16 GROUPS NOT CONVERGING
#CAN SEE IN WORD DOC: SPLITNORMMODELS_CONVERGENCE.docx that 

library(rstan)
library(tidyverse)
library(splines)
library(cowplot)
library(sn)
library(e1071)

###################################
#DO IN LOOP TO RECORD FOR ALL GROUPS:
DF<-read.csv("~/OneDrive/Documents/UAH/Weightofnations/Code/Draft code Sep2025/Output/df_bmipopheight_2022_18plus.csv")
DF_ADULT_noNA<-DF %>% filter(!is.na(`prev_less18.5`))
DF_ADULT_noNA$MEAN_midpoints<-((18.5+10)/2*DF_ADULT_noNA$prev_less18.5)+((20+18.5)/2*DF_ADULT_noNA$prev_18.5to20)+((25+20)/2*DF_ADULT_noNA$prev_20to25)+((30+25)/2*DF_ADULT_noNA$prev_25to30)+((35+30)/2*DF_ADULT_noNA$prev_30to35)+((40+35)/2*DF_ADULT_noNA$prev_35to40)+((50+40)/2*DF_ADULT_noNA$prev_over40)
DF_ADULT_noNA$SD_midpoints <- sqrt(
  (( (18.5+10)/2 - DF_ADULT_noNA$MEAN_midpoints )^2 * DF_ADULT_noNA$prev_less18.5) +
    (( (20+18.5)/2 - DF_ADULT_noNA$MEAN_midpoints )^2 * DF_ADULT_noNA$prev_18.5to20) +
    (( (25+20)/2 - DF_ADULT_noNA$MEAN_midpoints )^2 * DF_ADULT_noNA$prev_20to25) +
    (( (30+25)/2 - DF_ADULT_noNA$MEAN_midpoints )^2 * DF_ADULT_noNA$prev_25to30) +
    (( (35+30)/2 - DF_ADULT_noNA$MEAN_midpoints )^2 * DF_ADULT_noNA$prev_30to35) +
    (( (40+35)/2 - DF_ADULT_noNA$MEAN_midpoints )^2 * DF_ADULT_noNA$prev_35to40) +
    (( (50+40)/2 - DF_ADULT_noNA$MEAN_midpoints )^2 * DF_ADULT_noNA$prev_over40))

DF_ADULT_noNA$skewness<-NA
DF_ADULT_noNA$omega<-NA
DF_ADULT_noNA$alpha<-NA
DF_ADULT_noNA$xi<-NA
for(i in 1:nrow(DF_ADULT_noNA)){
  
  row<-as.list(DF_ADULT_noNA[i,])
  
  # Create synthetic data from bin midpoints and prevalence
  bmi_values <- c((18.5+10)/2, (20+18.5)/2, (25+20)/2, (30+25)/2, (35+30)/2, (40+35)/2, (50+40)/2)
  weights <- c(row$prev_less18.5, row$prev_18.5to20, row$prev_20to25, row$prev_25to30, row$prev_30to35, row$prev_35to40, row$prev_over40)
  
  # Repeat midpoints according to prevalence (scaled to sample size)
  synthetic_data <- rep(bmi_values, times = round(weights * row$total_population *100))
  gamma1 <- skewness(synthetic_data, type=1)
  delta <- sign(gamma1) * sqrt( (pi/2) * (abs(gamma1))^(2/3) / ( (abs(gamma1))^(2/3) + ((4 - pi)/2)^(2/3) ) )
  if(abs(delta)>=1){
    delta<-sign(delta)*0.99999 #restrict delta to be between -1 and 1
  }
  alpha <- delta / sqrt(1 - delta^2)
  omega <- row$SD_midpoints / sqrt(1 - (2 * delta^2 / pi))
  xi <- row$MEAN_midpoints - omega * delta * sqrt(2 / pi)
  
  
  DF_ADULT_noNA[i,"skewness"]<-gamma1
  DF_ADULT_noNA[i,"delta"]<-delta
  DF_ADULT_noNA[i,"alpha"]<-alpha
  DF_ADULT_noNA[i,"omega"]<-omega
  DF_ADULT_noNA[i,"xi"]<-xi
  
}

bin_edges <- c(10, 18.5, 20, 25, 30, 35, 40, 50)


results <- data.frame(
  index = integer(),
  name = character(),
  sex = character(),
  agegroup = character(),
  sigma1 = numeric(),
  sigma2 = numeric(),
  pivot = numeric(),
  total_abs_diff_bayes = numeric(),
  total_abs_diff_skew = numeric(),
  stringsAsFactors = FALSE
)

for (i in nrow(DF_ADULT_noNA)) {
  
  # Observed proportions
  proportions <- c(
    DF_ADULT_noNA[i, "prev_less18.5"][[1]],
    DF_ADULT_noNA[i, "prev_18.5to20"][[1]],
    DF_ADULT_noNA[i, "prev_20to25"][[1]],
    DF_ADULT_noNA[i, "prev_25to30"][[1]],
    DF_ADULT_noNA[i, "prev_30to35"][[1]],
    DF_ADULT_noNA[i, "prev_35to40"][[1]],
    DF_ADULT_noNA[i, "prev_over40"][[1]]
  )
  
  # Skip if missing data
  if (any(is.na(proportions))) next
  
  # Bin structure
  bmi_bins <- data.frame(
    mid = c((18.5 + 10) / 2, 19.25, 22.5, 27.5, 32.5, 37.5, 45),
    width = c(8.5, 1.5, 5, 5, 5, 5, 10),
    prevalence = proportions
  )
  
  bmi_bins$density <- bmi_bins$prevalence / bmi_bins$width
  sorted_indices <- order(bmi_bins$density, decreasing = TRUE)
  max_index <- sorted_indices[1]
  second_max_index <- sorted_indices[2]
  
  # Find index of highest-density bin
  max_index <- which.max(bmi_bins$density)
  
  highest_density_bin_lower <- bmi_bins[max_index,]$mid-bmi_bins[max_index,]$width/2
  highest_density_bin_upper <- bmi_bins[max_index,]$mid+bmi_bins[max_index,]$width/2
  midpoint_density<-bmi_bins[max_index,]$mid
  
  # Convert to counts
  N <- 100000
  y <- round(proportions * N)
  
  stan_data <- list(
    K = length(y),
    y = y,
    N = N,
    bin_edges = bin_edges,
    pivot_prior1=highest_density_bin_lower,
    pivot_prior2=highest_density_bin_upper
  )
  
  # Fit Bayesian two-half-normal model
  fit_n <- tryCatch({
    stan(
      file = "halfnormals_scaled_pivotprior.stan",
      data = stan_data,
      iter = 4000,
      chains = 4,
      seed = 123,
      refresh = 0,
      control = list(adapt_delta = 0.99, max_treedepth = 15)
    )
  }, error = function(e) {
    cat("Skipping index", i, "- half-normal model failed:", conditionMessage(e), "\n")
    return(structure(NULL, class = "try-error"))
  })
  
  if (inherits(fit_n, "try-error")) next
  
  
  posterior <- rstan::extract(fit_n)
  
  # Posterior means of parameters
  sigma1_fit <- mean(posterior$sigma1)
  sigma2_fit <- mean(posterior$sigma2)
  pivot_fit <- mean(posterior$pivot)
  p_fit <- colMeans(posterior$p)
  summary_fit <- summary(fit_n)$summary
  rhat_sigma1 <- summary_fit["sigma1", "Rhat"]
  rhat_sigma2 <- summary_fit["sigma2", "Rhat"]
  rhat_pivot  <- summary_fit["pivot", "Rhat"]
  
  
  #NORMAL (for comparison - of two half normal vs normal)
  stan_data_normal <- list(
    K = length(y),
    y = y,
    N = N,
    bin_edges = bin_edges,
    mean_prior=midpoint_density
  )
  
  fit_n_normal <- tryCatch({
    stan(
      file = "normal_forcompar.stan",
      data = stan_data_normal,
      iter = 4000,
      chains = 4,
      seed = 123,
      refresh = 0,
      control = list(adapt_delta = 0.99, max_treedepth = 15)
    )
  }, error = function(e) {
    cat("Skipping index", i, "- normal model failed:", conditionMessage(e), "\n")
    return(structure(NULL, class = "try-error"))
  })
  
  if (inherits(fit_n_normal, "try-error")) next
  
  
  posterior_normal <- rstan::extract(fit_n_normal)
  sigma_fit_normal <- mean(posterior_normal$sigma)
  mu_fit_normal <- mean(posterior_normal$mu)
  p_fit_normal <- colMeans(posterior_normal$p)
  summary_fit_normal <- summary(fit_n_normal)$summary
  rhat_sigma_normal <- summary_fit_normal["sigma", "Rhat"]
  rhat_mu_normal <- summary_fit_normal["mu", "Rhat"]

  # Extract skew-normal parameters
  xi <- DF_ADULT_noNA[i, "xi"][[1]]
  omega <- DF_ADULT_noNA[i, "omega"][[1]]
  alpha <- DF_ADULT_noNA[i, "alpha"][[1]]
  
  # Compute predicted proportions from skew-normal
  predicted_skew <- numeric(length(bin_edges) - 1)
  for (k in seq_along(predicted_skew)) {
    predicted_skew[k] <- psn(bin_edges[k + 1], xi, omega, alpha) -
      psn(bin_edges[k], xi, omega, alpha)
  }
  
  # Compute total absolute differences
  total_abs_diff_bayes <- sum(abs(p_fit - proportions))
  total_abs_diff_bayes_normal <- sum(abs(p_fit_normal - proportions))
  total_abs_diff_skew <- sum(abs(predicted_skew - proportions))
  
  # Store results
  results <- rbind(
    results,
    data.frame(
      index = i,
      name = DF_ADULT_noNA[i, "name"][[1]],
      sex = DF_ADULT_noNA[i, "SEX"][[1]],
      agegroup = DF_ADULT_noNA[i, "heightagegroups"][[1]],
      sigma1 = sigma1_fit,
      sigma2 = sigma2_fit,
      pivot = pivot_fit,
      sigma_normal = sigma_fit_normal,
      mu_normal = mu_fit_normal,
      total_abs_diff_bayes = total_abs_diff_bayes,
      total_abs_diff_bayes_normal = total_abs_diff_bayes_normal,
      total_abs_diff_skew = total_abs_diff_skew,
      rhat_pivot = rhat_pivot,
      rhat_sigma1 = rhat_sigma1,
      rhat_sigma2 = rhat_sigma2
    )
  )
  
  cat("Finished:", i, "of", nrow(DF_ADULT_noNA), "\n")
}

boxplot(results$rhat_pivot)
boxplot(results$rhat_sigma1)
boxplot(results$rhat_sigma2)

results %>% filter(rhat_pivot>1.05)

results %>% filter(rhat_pivot>1.05) %>% summarise(max=max(pivot),min=min(pivot), sum=sum(sex))

plot(results$rhat_pivot, results$pivot)

mean(results$total_abs_diff_bayes)
mean(results$total_abs_diff_bayes_normal)
mean(results$total_abs_diff_skew)
ggplot(results,aes(x=total_abs_diff_bayes,y=total_abs_diff_skew))+geom_point()+
  geom_abline(intercept=0,slope=1, color="red") #Bayes two half-normal model does much better than skew-normal from synthetic data

ggplot(results,aes(x=total_abs_diff_bayes,y=total_abs_diff_bayes_normal))+geom_point()+
  geom_abline(intercept=0,slope=1, color="red") #Bayes two half-normal model does much better than Bayes normal model

ggplot(results,aes(x=total_abs_diff_skew,y=total_abs_diff_bayes_normal))+geom_point()+
  geom_abline(intercept=0,slope=1, color="red") #skew-normal from synthetic data and Bayesian normal model do similarly well

mean(results$total_abs_diff_bayes<=results$total_abs_diff_skew)
max(-results$total_abs_diff_bayes+results$total_abs_diff_skew)

boxplot(results$total_abs_diff_bayes)

results$bayesvsskew<-results$total_abs_diff_bayes-results$total_abs_diff_skew

results %>% arrange(desc(bayesvsskew)) #Iceland, sex=0, age=60-69, index=1160 Bayes: 0.24, skew: 0.047

results %>% arrange(desc(total_abs_diff_bayes)) 

results_ranked <- results %>%
  arrange(desc(total_abs_diff_bayes)) %>%
  mutate(rank_diff = row_number())

plot(results_ranked$rank_diff,results_ranked$rhat_pivot)
plot(results_ranked$total_abs_diff_bayes,results_ranked$rhat_pivot)


saveRDS(results,"twohalfnorms_results_newmid_rules_modeshift_2022_31Oct_pivotprior_pluslimit_plusnormal_updatedpriors.RDS")

#MAKE PLOTS BY COUNTRY:
DF_ADULT_noNA$name <- if_else(DF_ADULT_noNA$name=="C\x99te d'Ivoire", "Côte d'Ivoire", DF_ADULT_noNA$name)

plot_list_country<-list()
for(j in 1:length(unique(as.factor(DF_ADULT_noNA$name)))){
  DF_ADULT_noNA_country<-DF_ADULT_noNA %>% filter(as.numeric(as.factor(name))==j)
  results_country <- results %>% filter(as.numeric(as.factor(name))==j)
  plot_list<-list()
  for(i in results_country$index){
    
    proportions <- c(
      DF_ADULT_noNA[i, "prev_less18.5"][[1]],
      DF_ADULT_noNA[i, "prev_18.5to20"][[1]],
      DF_ADULT_noNA[i, "prev_20to25"][[1]],
      DF_ADULT_noNA[i, "prev_25to30"][[1]],
      DF_ADULT_noNA[i, "prev_30to35"][[1]],
      DF_ADULT_noNA[i, "prev_35to40"][[1]],
      DF_ADULT_noNA[i, "prev_over40"][[1]]
    )
    
    country <- DF_ADULT_noNA[i, "name"][[1]]
    sex <- DF_ADULT_noNA[i, "SEX"][[1]]
    age <- DF_ADULT_noNA[i, "heightagegroups"][[1]]
    
    title_text <- paste(country, sex, age, sep = " - ")
    
    if (any(is.na(proportions))) next
    bmi_bins <- data.frame(
      mid = c((18.5 + 10) / 2, 19.25, 22.5, 27.5, 32.5, 37.5, 45),
      width = c(8.5, 1.5, 5, 5, 5, 5, 10),
      prevalence = proportions
    )
    
    pivot_t <- results[i,"pivot"]
    sigma1_t <- results[i,"sigma1"]
    sigma2_t <- results[i,"sigma2"]
    two_half_normal <- function(x, pivot, sigma1, sigma2) {
      ifelse(
        x <= pivot,
        2/(sigma1+sigma2) * dnorm((x-pivot)/sigma1),
        2/(sigma1+sigma2) * dnorm((x-pivot)/sigma2)
      )
    }
    
    mu_norm <- results[i,"mu_normal"]
    sigma_norm <- results[i,"sigma_normal"]
    
    p <- ggplot(bmi_bins, aes(x = mid, y = prevalence / width)) +
      geom_col(aes(width = width), fill = "steelblue", alpha = 0.7) +
      stat_function(
        fun = local({
          pivot <- pivot_t; sigma1 <- sigma1_t; sigma2 <- sigma2_t
          function(x) two_half_normal(x, pivot, sigma1, sigma2)
        }),
        color = "purple", size = 1.2, n = 500
      ) +
      stat_function(
        fun = local({
          mu <- mu_norm; sigma <- sigma_norm
          function(x) dnorm(x, mu, sigma)
        }),
        color = "red", size = 1.2, n = 500
      ) +
      stat_function(
        fun = local({
          xi <- DF_ADULT_noNA[i,"xi"][[1]]
          omega <- DF_ADULT_noNA[i,"omega"][[1]]
          alpha <- DF_ADULT_noNA[i,"alpha"][[1]]
          function(x) dsn(x, xi = xi, omega = omega, alpha = alpha)
        }),
        color = "darkgreen", size = 1.2, n = 500
      ) +
      labs(
        title= title_text,
        x = "BMI",
        y = "Density"
      ) +
      theme_minimal()
    plot_list[[length(plot_list)+1]]<-p
  }
  plot_list_country[[length(plot_list_country)+1]]<-plot_grid(plotlist=plot_list,ncol=2)
  print(country)
}


# Open a PDF device
pdf("country_plots_2022_31Oct_pivotprior_pluslimit_plusnormalposterior_updatedpriors.pdf", width = 8, height = 10) 
plot.new()
text(0.5, 0.6, "Green = Skewed normal distribution from synthetic data", cex = 1.2)
text(0.5, 0.4, "Purple = Bayesian multinomial model with two half-normal distributions", cex = 1.2)
text(0.5, 0.2, "Red = Bayesian multinomial model with one normal distribution", cex = 1.2)
for (pp in 1:length(plot_list_country)) {
  print(plot_list_country[[pp]])
}
dev.off()
