# Mixed effects binomial in R for lesion paper

library(lme4)
library(stringr)

# Read the data frame.

fileIn <- '/Users/frederictheunissen/Code/zebra-finch-memory-lesions/data/behavior/callCountS2.csv'
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


base_model <- glmer('cbind(ints, tots-ints) ~ Treatment*Rewarded + ((1 + CallType*Rewarded)|Subject)', data = post_counts, family = binomial)
callType_model <- glmer('cbind(ints, tots-ints) ~ Treatment*Rewarded + CallType*Rewarded + ((1 + CallType*Rewarded)|Subject)', data = post_counts, family = binomial)
full_model <- glmer('cbind(ints, tots-ints) ~ Treatment*Rewarded*CallType + ((1 + CallType*Rewarded)|Subject)', data = post_counts, family = binomial)

# This is significant because call type changes interruption rates
anova(base_model, full_model)

# Note that none of the CallType fixed effect parameters are significant in the full model.
summary(full_model)


for (ctype in c('SO', 'DC')) {
  
  print(sprintf("===================== %s ===============================", ctype))
  
  for (treat in c('CTRL', 'HVC', 'NCM')) {
    treat_model <- glmer('cbind(ints, tots-ints) ~ Rewarded +((1 + Rewarded)|Subject)',
                         data = post_counts,
                         subset = ((post_counts$Treatment == treat)&(post_counts$CallType == ctype)),
                         family = binomial)
    
    treat_modelS <- summary(treat_model)
    print(sprintf('Treatment: %s', treat))
    print(sprintf('\t Log2(OR)= %.2f +- %.3f Z=%.3f p = %.4f', 
                  (1/log(2))*treat_modelS$coefficients[2,1],
                  (1/log(2))*treat_modelS$coefficients[2,2],
                  treat_modelS$coefficients[2,3],
                  treat_modelS$coefficients[2,4]))
  }
  
  full_model <- glmer('cbind(ints, tots-ints) ~ Treats*Rewarded + ((1 + Rewarded)|Subject)', 
                      data = post_counts, 
                      subset = (post_counts$CallType == ctype),
                      family = binomial)
  full_modelS <- summary(full_model)
  print(sprintf('Difference NCM vc HVC+Control: Log2(OR)= %.2f +- %.3f Z=%.3f p = %.4f', 
                (1/log(2))*full_modelS$coefficients[4,1],
                (1/log(2))*full_modelS$coefficients[4,2],
                full_modelS$coefficients[4,3],
                full_modelS$coefficients[4,4]))
}
