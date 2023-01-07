# Mixed effects binomial in R for lesion paper

library(lme4)
library(stringr)

# Read the data frame.

fileIn <- '/Users/frederictheunissen/Code/zebra-finch-memory-lesions/data/behavior/fig4cCount.csv'
post_counts <- read.csv(fileIn)

# make new colums for the count of interruption
ints <- str_split(post_counts$Counts, ',')
n <- length(ints)

for (i in 1:n) {
  post_counts$ints[i] <- strtoi(substr(ints[[i]][1], start=2, stop=nchar(ints[[i]][1])))
  post_counts$tots[i] <- strtoi(substr(ints[[i]][2], start=2, stop=nchar(ints[[i]][2])-1))
  if (post_counts$Treatment[i] == 'NCM') {
    post_counts$Treats[i] = 'NCM'
  }
  else {
    post_counts$Treats[i] = 'HVC+Ctrl'
    
  }
}

for (treat in c('CTRL', 'HVC', 'NCM')) {
  treat_model <- glmer('cbind(ints, tots-ints) ~ Rewarded +((1 + Rewarded)|Subject)', data = post_counts, subset = (post_counts$Treatment == treat), family = binomial)
  treat_modelS <- summary(treat_model)
  print(sprintf('Treatment: %s', treat))
  print(sprintf('\t Log2(OR)= %.2f +- %.3f Z=%.3f p = %.4f', 
                (1/log(2))*treat_modelS$coefficients[2,1],
                (1/log(2))*treat_modelS$coefficients[2,2],
                treat_modelS$coefficients[2,3],
                treat_modelS$coefficients[2,4]))
}


# Now for HVC and Control
treat_model <- glmer('cbind(ints, tots-ints) ~ Rewarded +((1 + Rewarded)|Subject)', data = post_counts, subset = (post_counts$Treats == "HVC+Ctrl"), family = binomial)
treat_modelS <- summary(treat_model)
print(sprintf('Treatment: %s', "HVC+Ctrl"))
print(sprintf('\t Log2(OR)= %.2f +- %.3f Z=%.3f p = %.4f', 
              (1/log(2))*treat_modelS$coefficients[2,1],
              (1/log(2))*treat_modelS$coefficients[2,2],
              treat_modelS$coefficients[2,3],
              treat_modelS$coefficients[2,4]))

null_model <- glmer('cbind(ints, tots-ints) ~ Rewarded + ((1 + Rewarded)|Subject)', data = post_counts, family = binomial)
full_model <- glmer('cbind(ints, tots-ints) ~ Treatment*Rewarded + ((1 + Rewarded)|Subject)', data = post_counts, family = binomial)

# Comparing the full model to the null model
print('null vs full')
print("---------")
print(anova(null_model, full_model))

full_model <- glmer('cbind(ints, tots-ints) ~ Treats*Rewarded + ((1 + Rewarded)|Subject)', data = post_counts, family = binomial)
full_modelS <- summary(full_model)
print(sprintf('\t Log2(OR)= %.2f +- %.3f Z=%.3f p = %.4f', 
              (1/log(2))*full_modelS$coefficients[4,1],
              (1/log(2))*full_modelS$coefficients[4,2],
              full_modelS$coefficients[4,3],
              full_modelS$coefficients[4,4]))

# Coefficients for the slopes and intercepts for Table 4.
print(sprintf('Intercept of Control: %.2f +- %.2f', 
              -(1/log(2))*sum_full_model$coefficients['RewardedTrue', 'Estimate'],
              (1/log(2))*sum_full_model$coefficients['RewardedTrue', 'Std. Error']))

print(sprintf('Slope of Control: %.2f +- %.2f', 
              -(1/log(2))*sum_full_model$coefficients['k:RewardedTrue', 'Estimate'],
              (1/log(2))*sum_full_model$coefficients['k:RewardedTrue', 'Std. Error']))

k_model <- glmer('cbind(ints, tots-ints) ~ k*Rewarded + ((1 + Rewarded*k)|Subject)', data = learning_curve_counts, family = binomial, subset = (learning_curve_counts$Treatment == 'NCM'))
sum_k_model <- summary(k_model)
print(sprintf('Intercept of NCM: %.2f +- %.2f', 
              -(1/log(2))*sum_k_model$coefficients['RewardedTrue', 'Estimate'],
              (1/log(2))*sum_k_model$coefficients['RewardedTrue', 'Std. Error']))

print(sprintf('Slope of NCM: %.2f +- %.2f', 
              -(1/log(2))*sum_k_model$coefficients['k:RewardedTrue', 'Estimate'],
              (1/log(2))*sum_k_model$coefficients['k:RewardedTrue', 'Std. Error']))

# Does one get the same result from the full model?
print(sprintf('Intercept of NCM (from Full Model): %.2f +- %.2f', 
              -(1/log(2))*(sum_full_model$coefficients['RewardedTrue', 'Estimate'] + sum_full_model$coefficients['TreatsNCM:RewardedTrue', 'Estimate']) ,
              (1/log(2))*sum_full_model$coefficients['TreatsNCM:RewardedTrue', 'Std. Error']))

print(sprintf('Slope of NCM: %.2f +- %.2f', 
              -(1/log(2))*(sum_full_model$coefficients['k:RewardedTrue', 'Estimate'] + sum_full_model$coefficients['TreatsNCM:k:RewardedTrue', 'Estimate']),
              (1/log(2))*sum_full_model$coefficients['TreatsNCM:k:RewardedTrue', 'Std. Error']))



