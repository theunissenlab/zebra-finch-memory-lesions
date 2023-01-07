# Mixed effects binomial in R for lesion paper

library(lme4)
library(stringr)

# Read the data frame.

fileIn <- '/Users/frederictheunissen/Code/zebra-finch-memory-lesions/data/behavior/before_d1Count.csv'
learning_curve_counts <- read.csv(fileIn)

# make new colums for the count of interruption
ints <- str_split(learning_curve_counts$Counts, ',')
n <- length(ints)

for (i in 1:n) {
  learning_curve_counts$ints[i] <- strtoi(substr(ints[[i]][1], start=2, stop=nchar(ints[[i]][1])))
  learning_curve_counts$tots[i] <- strtoi(substr(ints[[i]][2], start=2, stop=nchar(ints[[i]][2])-1))
  if (learning_curve_counts$Treatment[i] == 'NCM') {
    learning_curve_counts$Treats[i] = 'NCM'
  }
  else {
    learning_curve_counts$Treats[i] = 'HVC+Ctrl'
    
  }
}



k_model <- glmer('cbind(ints, tots-ints) ~ k*Rewarded + ((1 + Rewarded*k)|Subject)', data = learning_curve_counts, family = binomial)
full_model <- glmer('cbind(ints, tots-ints) ~ Treats*k*Rewarded + ((1 + Rewarded*k)|Subject)', data = learning_curve_counts, family = binomial)
anova(k_model, full_model)

#
coefFull <- full_model@beta
kval <- 0:max(learning_curve_counts$k)
logOddnoreCtrl <- coefFull[1] + kval*coefFull[3] 
logOddreCtrl <- coefFull[1] + kval*coefFull[3] + coefFull[4] + kval*coefFull[7]
log2OddsRatCtrl <- (1/log(2))*(logOddnoreCtrl-logOddreCtrl)
log2OddsRatCtrl <- (1/log(2))*(-coefFull[4] - kval*coefFull[7])

logOddnoreNCM <- coefFull[1] + coefFull[2]+ kval*(coefFull[3] + coefFull[5])
logOddreNCM<- coefFull[1] + coefFull[2] + kval*(coefFull[3] + coefFull[5]) +
  (coefFull[4] + coefFull[6]) + kval*(coefFull[7] + coefFull[8])
log2OddsRatNCM <- (1/log(2))*(logOddnoreNCM-logOddreNCM)
log2OddsRatNCM <-  (1/log(2))*(-coefFull[4] - coefFull[6]) + kval*(-coefFull[7] - coefFull[8])


plot(kval, log2OddsRatCtrl, col = 'black', type = 'l' , ylim = c(-1,5))
lines(kval, log2OddsRatNCM, col = 'pink')


sum_full_model = summary(full_model)

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



