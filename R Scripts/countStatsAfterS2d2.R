# Mixed effects binomial in R for lesion paper

library(lme4)
library(stringr)

# Read the data frame.

fileIn <- '/Users/frederictheunissen/Code/zebra-finch-memory-lesions/data/behavior/afterS2d2Count.csv'
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



null_model <- glmer('cbind(ints, tots-ints) ~ Rewarded + ((1 + Rewarded)|Subject)', data = learning_curve_counts, family = binomial)
full_model <- glmer('cbind(ints, tots-ints) ~ Treats*Rewarded + ((1 + Rewarded)|Subject)', data = learning_curve_counts, family = binomial)
anova(null_model, full_model)

#
coefFull <- full_model@beta

logOddnoreCtrl <- coefFull[1] 
logOddreCtrl <- coefFull[1] + coefFull[3] 
log2OddsRatCtrl <- (1/log(2))*(logOddnoreCtrl-logOddreCtrl)
log2OddsRatCtrl <- (1/log(2))*(-coefFull[3])

logOddnoreNCM <- coefFull[1] + coefFull[2] 
logOddreNCM<- coefFull[1] + coefFull[2] + (coefFull[3] + coefFull[4]) 
log2OddsRatNCM <- (1/log(2))*(logOddnoreNCM-logOddreNCM)
log2OddsRatNCM <-  (1/log(2))*(-coefFull[3] - coefFull[4]) 


# Odds and difference in Odds ratios
sum_full_model <- summary(full_model)
print(sprintf('Control+HVC: Log2(OR) =  %.2f', 
              -(1/log(2))*sum_full_model$coefficients['RewardedTrue', 'Estimate']))

print(sprintf('NCM: log2(OR) = %.2f', 
              -(1/log(2))*(sum_full_model$coefficients['RewardedTrue', 'Estimate']+ sum_full_model$coefficients['TreatsNCM:RewardedTrue', 'Estimate'])))
  
print(sprintf('Difference in Log2(OR) = %.2f +- %.2f', 
             -(1/log(2))*sum_full_model$coefficients['TreatsNCM:RewardedTrue', 'Estimate'],
             (1/log(2))*sum_full_model$coefficients['TreatsNCM:RewardedTrue', 'Std. Error']))

print(sprintf('Wald test for Difference Z = %.2f p = %.3f', 
             -(1/log(2))*sum_full_model$coefficients['TreatsNCM:RewardedTrue', 'z value'],
             (1/log(2))*sum_full_model$coefficients['TreatsNCM:RewardedTrue', 'Pr(>|z|)']))
