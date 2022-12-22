# Mixed effects binomial in R for lesion paper

library(lme4)
library(stringr)

# Read the data frame.

fileIn <- '/Users/frederictheunissen/Code/zebra-finch-memory-lesions/data/behavior/fig5bCount.csv'
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
k_model <- glmer('cbind(ints, tots-ints) ~ k*Rewarded + ((1 + Rewarded*k)|Subject)', data = learning_curve_counts, family = binomial)
T_model <- glmer('cbind(ints, tots-ints) ~ Treatment*Rewarded + ((1 + Rewarded)|Subject)', data = learning_curve_counts, family = binomial)
Tandk_model <- glmer('cbind(ints, tots-ints) ~ (Treatment + k)*Rewarded + ((1 + Rewarded*k)|Subject)', data = learning_curve_counts, family = binomial)
full_model <- glmer('cbind(ints, tots-ints) ~ Treatment*k*Rewarded + ((1 + Rewarded*k)|Subject)', data = learning_curve_counts, family = binomial)

# Comparing the k model to the full model
print('k vs full')
print("---------")
print("We are testing if a model of the learning curves that includes a lesion treatment (3 lines) is a better fit than one without (1 line)")
print(anova(k_model, full_model))

print('k vs T+k')
print("--------")
print("We are testing if a model of the learning curves that includes a lesion treatment (3 parallel lines) is a better fit than one without (1 line)")
print(anova(k_model, Tandk_model))

print('T vs full')
print("--------")
print("We are testing if a model of the learning curves that includes informative trials (3 lines with slopes) is a better fit than one without slope (3 horizontal lines)")
print(anova(T_model, full_model))

print('Null vs T')
print('---------')
print("We are testing if a model of the learning curves that includes Treatment (3 horizontal lines) is a better fit than one without treatment (1 horizontal line)")
print(anova(null_model, T_model))

print('T+k vs full')
print('---------')
print("We are testing if a model of the learning curves that includes Treatment and k (3 arbitrary lines) is a better fit than one without 3 parallel lines")
anova(Tandk_model, full_model)

# Make log(OR) plots from the coefficients of the Full model
coefFull <- full_model@beta
kval <- 0:max(learning_curve_counts$k)

# Control
logOddnoreCtrl <- coefFull[1] + kval*coefFull[4] 
logOddreCtrl <- coefFull[1] + kval*coefFull[4] + coefFull[5] + kval*coefFull[10]
log2OddsRatCtrl <- (1/log(2))*(logOddnoreCtrl-logOddreCtrl)

# HVC 
logOddnoreHVC <- coefFull[1] + coefFull[2]+ kval*(coefFull[4] + coefFull[6])
logOddreHVC <- coefFull[1] + coefFull[2] + kval*(coefFull[4] + coefFull[6]) +
  (coefFull[5] + coefFull[8]) + kval*(coefFull[10] + coefFull[11])
log2OddsRatHVC <- (1/log(2))*(logOddnoreHVC-logOddreHVC)

# NCM
logOddnoreNCM <- coefFull[1] + coefFull[3]+ kval*(coefFull[4] + coefFull[7])
logOddreNCM <- coefFull[1] + coefFull[3] + kval*(coefFull[4] + coefFull[7]) +
  (coefFull[5] + coefFull[9]) + kval*(coefFull[10] + coefFull[12])
log2OddsRatNCM <- (1/log(2))*(logOddnoreNCM-logOddreNCM)


plot(kval, log2OddsRatCtrl, col = 'black', type = 'l' , ylim = c(-1,5))
lines(kval, log2OddsRatHVC, col = 'pink')
lines(kval, log2OddsRatNCM, col = 'green')

# Repeat for pair-wise comparisons

# Control vs NCM
full_modelCtrlNCM <- glmer('cbind(ints, tots-ints) ~ Treatment*k*Rewarded + ((1 + k*Rewarded)|Subject)', 
                           data = learning_curve_counts, family = binomial,
                           subset = (learning_curve_counts$Treatment != 'HVC'))

coefFullCtrlNCM <- full_modelCtrlNCM@beta

logOddnoreCtrl <- coefFullCtrlNCM[1] + kval*coefFullCtrlNCM[3] 
logOddreCtrl <- coefFullCtrlNCM[1] + kval*coefFullCtrlNCM[3] + coefFullCtrlNCM[4] + kval*coefFullCtrlNCM[7]
log2OddsRatCtrl <- (1/log(2))*(logOddnoreCtrl-logOddreCtrl)
log2OddsRatCtrl <- (1/log(2))*(-coefFullCtrlNCM[4] - kval*coefFullCtrlNCM[7])

logOddnoreNCM <- coefFullCtrlNCM[1] + coefFullCtrlNCM[2]+ kval*(coefFullCtrlNCM[3] + coefFullCtrlNCM[5])
logOddreNCM <- coefFullCtrlNCM[1] + coefFullCtrlNCM[2] + kval*(coefFullCtrlNCM[3] + coefFullCtrlNCM[5]) +
  (coefFullCtrlNCM[4] + coefFullCtrlNCM[6]) + kval*(coefFullCtrlNCM[7] + coefFullCtrlNCM[8])
log2OddsRatNCM <- (1/log(2))*(logOddnoreNCM-logOddreNCM)
log2OddsRatNCM <-  (1/log(2))*(-coefFullCtrlNCM[4] - coefFullCtrlNCM[6]) + kval*(-coefFullCtrlNCM[7] - coefFullCtrlNCM[8])

plot(kval, log2OddsRatCtrl, col = 'black', type = 'l' , ylim = c(-1,5))
lines(kval, log2OddsRatNCM, col = 'green')

sum_full_modelCtrlNCM <- summary(full_modelCtrlNCM)
print(sprintf('Intercept difference (Ctrl - NCM) = %.2f +- %.3f', (1/log(2))*coefFullCtrlNCM[6], (1/log(2))*sum_full_modelCtrlNCM$coefficients[6,'Std. Error'] ) )
print(sprintf('                      Z = %.3f, p=%.4f', sum_full_modelCtrlNCM$coefficients[6,'z value'], sum_full_modelCtrlNCM$coefficients[6,'Pr(>|z|)']))

print(sprintf('Slope difference (Ctrl - NCM) = %.2f +- %.3f', (1/log(2))*coefFullCtrlNCM[8], (1/log(2))*sum_full_modelCtrlNCM$coefficients[8,'Std. Error'] ) )
print(sprintf('                      Z = %.3f, p=%.4f', sum_full_modelCtrlNCM$coefficients[8,'z value'], sum_full_modelCtrlNCM$coefficients[8,'Pr(>|z|)']))

# HVC vs NCM
full_modelHVCNCM <- glmer('cbind(ints, tots-ints) ~ Treatment*k*Rewarded + ((1 + Rewarded*k)|Subject)', 
                           data = learning_curve_counts, family = binomial,
                           subset = (learning_curve_counts$Treatment != 'CTRL'))

coefFullHVCNCM <- full_modelHVCNCM@beta
# HVC
logOddnoreHVC <- coefFullHVCNCM[1] + kval*coefFullHVCNCM[3] 
logOddreHVC <- coefFullHVCNCM[1] + kval*coefFullHVCNCM[3] + coefFullHVCNCM[4] + kval*coefFullHVCNCM[7]
log2OddsRatHVC <- (1/log(2))*(logOddnoreHVC-logOddreHVC)
log2OddsRatHVC <- (1/log(2))*(-coefFullHVCNCM[4] - kval*coefFullHVCNCM[7])

#NCM
logOddnoreNCM <- coefFullCtrlNCM[1] + coefFullCtrlNCM[2]+ kval*(coefFullCtrlNCM[3] + coefFullCtrlNCM[5])
logOddreNCM <- coefFullCtrlNCM[1] + coefFullCtrlNCM[2] + kval*(coefFullCtrlNCM[3] + coefFullCtrlNCM[5]) +
  (coefFullCtrlNCM[4] + coefFullCtrlNCM[6]) + kval*(coefFullCtrlNCM[7] + coefFullCtrlNCM[8])
log2OddsRatNCM <- (1/log(2))*(logOddnoreNCM-logOddreNCM)
log2OddsRatNCM <- (1/log(2))*(-coefFullCtrlNCM[4] - coefFullCtrlNCM[6]) + kval*(-coefFullCtrlNCM[7] - coefFullCtrlNCM[8])

plot(kval, log2OddsRatHVC, col = 'pink', type = 'l' , ylim = c(-1,5))
lines(kval, log2OddsRatNCM, col = 'green')

sum_full_modelHVCNCM <- summary(full_modelHVCNCM)
print(sprintf('Intercept difference (HVC - NCM) = %.2f +- %.3f', (1/log(2))*coefFullHVCNCM[6], (1/log(2))*sum_full_modelHVCNCM$coefficients[6,'Std. Error'] ) )
print(sprintf('                      Z = %.3f, p=%.4f', sum_full_modelHVCNCM$coefficients[6,'z value'], sum_full_modelHVCNCM$coefficients[6,'Pr(>|z|)']))

print(sprintf('Slope difference (HVC - NCM) = %.2f +- %.3f', (1/log(2))*coefFullHVCNCM[8], (1/log(2))*sum_full_modelHVCNCM$coefficients[8,'Std. Error'] ) )
print(sprintf('                      Z = %.3f, p=%.4f', sum_full_modelHVCNCM$coefficients[8,'z value'], sum_full_modelHVCNCM$coefficients[8,'Pr(>|z|)']))

# Ctrl vs HVC
full_modelCtrlHVC <- glmer('cbind(ints, tots-ints) ~ Treatment*k*Rewarded + ((1 + k*Rewarded)|Subject)', 
                          data = learning_curve_counts, family = binomial,
                          subset = (learning_curve_counts$Treatment != 'NCM'))

coefFullCtrlHVC <- full_modelCtrlHVC@beta

logOddnoreCtrl <- coefFullCtrlHVC[1] + kval*coefFullCtrlHVC[3] 
logOddreCtrl <- coefFullCtrlHVC[1] + kval*coefFullCtrlHVC[3] + coefFullCtrlHVC[4] + kval*coefFullCtrlHVC[7]
log2OddsRatCtrl <- (1/log(2))*(logOddnoreCtrl-logOddreCtrl)
log2OddsRatCtrl <- (1/log(2))*(-coefFullCtrlHVC[4] - kval*coefFullCtrlHVC[7])

logOddnoreHVC <- coefFullCtrlHVC[1] + coefFullCtrlHVC[2]+ kval*(coefFullCtrlHVC[3] + coefFullCtrlHVC[5])
logOddreHVC <- coefFullCtrlHVC[1] + coefFullCtrlHVC[2] + kval*(coefFullCtrlHVC[3] + coefFullCtrlHVC[5]) +
  (coefFullCtrlHVC[4] + coefFullCtrlHVC[6]) + kval*(coefFullCtrlHVC[7] + coefFullCtrlHVC[8])
log2OddsRatHVC <- (1/log(2))*(logOddnoreHVC-logOddreHVC)
log2OddsRatHVC <-  (1/log(2))*(-coefFullCtrlHVC[4] - coefFullCtrlHVC[6]) + kval*(-coefFullCtrlHVC[7] - coefFullCtrlHVC[8])

plot(kval, log2OddsRatCtrl, col = 'black', type = 'l' , ylim = c(-1,5))
lines(kval, log2OddsRatHVC, col = 'pink')

sum_full_modelCtrlHVC <- summary(full_modelCtrlHVC)
print(sprintf('Intercept difference (Ctrl - HVC) = %.2f +- %.3f', (1/log(2))*coefFullCtrlHVC[6], (1/log(2))*sum_full_modelCtrlHVC$coefficients[6,'Std. Error'] ) )
print(sprintf('                      Z = %.3f, p=%.4f', sum_full_modelCtrlHVC$coefficients[6,'z value'], sum_full_modelCtrlHVC$coefficients[6,'Pr(>|z|)']))

print(sprintf('Slope difference (Ctrl - HVC) = %.2f +- %.3f', (1/log(2))*coefFullCtrlHVC[8], (1/log(2))*sum_full_modelCtrlHVC$coefficients[8,'Std. Error'] ) )
print(sprintf('                      Z = %.3f, p=%.4f', sum_full_modelCtrlHVC$coefficients[8,'z value'], sum_full_modelCtrlHVC$coefficients[8,'Pr(>|z|)']))


# Repeat the model comparison after combining HVC and control
k_model <- glmer('cbind(ints, tots-ints) ~ k*Rewarded + ((1 + Rewarded*k)|Subject)', data = learning_curve_counts, family = binomial)
full_model <- glmer('cbind(ints, tots-ints) ~ Treats*k*Rewarded + ((1 + Rewarded*k)|Subject)', data = learning_curve_counts, family = binomial)
anova(k_model, full_model)

#
coefFull <- full_model@beta
logOddnoreCtrl <- coefFull[1] + kval*coefFull[3] 
logOddreCtrl <- coefFull[1] + kval*coefFull[3] + coefFull[4] + kval*coefFull[7]
log2OddsRatCtrl <- (1/log(2))*(logOddnoreCtrl-logOddreCtrl)
log2OddsRatCtrl <- (1/log(2))*(-coefFull[4] - kval*coefFull[7])

logOddnoreNCM <- coefFull[1] + coefFull[2]+ kval*(coefFull[3] + coefFull[5])
logOddreNCM<- coefFull[1] + coefFull[2] + kval*(coefFull[3] + coefFull[5]) +
  (coefFull[4] + coefFull[6]) + kval*(coefFull[7] + coefFull[8])
log2OddsRatNCM <- (1/log(2))*(logOddnoreHVC-logOddreHVC)
log2OddsRatNCM <-  (1/log(2))*(-coefFull[4] - coefFull[6]) + kval*(-coefFull[7] - coefFull[8])


sum_full_model <- summary(full_model)

print(sprintf('Intercept difference (NCM - HVC+Ctrl) = %.2f +- %.3f', (1/log(2))*coefFull[6], (1/log(2))*sum_full_model$coefficients[6,'Std. Error'] ) )
print(sprintf('                      Z = %.3f, p=%.4f', sum_full_model$coefficients[6,'z value'], sum_full_model$coefficients[6,'Pr(>|z|)']))

