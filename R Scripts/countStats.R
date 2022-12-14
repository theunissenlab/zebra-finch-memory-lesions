# Mixed effects binomial in R for lesion paper

library(lme4)
library(stringr)

# Read the data frame.

fileIn <- '/Users/frederictheunissen/Code/zebra-finch-memory-lesions/data/behavior/fig5Count.csv'
learning_curve_counts <- read.csv(fileIn)

# make new colums for the count of interruption
ints <- str_split(learning_curve_counts$Counts, ',')
n <- length(ints)

for (i in 1:n) {
  learning_curve_counts$ints[i] <- strtoi(substr(ints[[i]][1], start=2, stop=nchar(ints[[i]][1])))
  learning_curve_counts$tots[i] <- strtoi(substr(ints[[i]][2], start=2, stop=nchar(ints[[i]][2])-1))
}

null_model <- glmer('cbind(ints, tots-ints) ~ Rewarded + (1 + Rewarded|Subject)', data = learning_curve_counts, family = binomial)
k_model <- glmer('cbind(ints, tots-ints) ~ k*Rewarded + (1 + Rewarded|Subject)', data = learning_curve_counts, family = binomial)
T_model <- glmer('cbind(ints, tots-ints) ~ Treatment*Rewarded + (1 + Rewarded|Subject)', data = learning_curve_counts, family = binomial)
Tandk_model <- glmer('cbind(ints, tots-ints) ~ (Treatment + k)*Rewarded + (1 + Rewarded|Subject)', data = learning_curve_counts, family = binomial)
full_model <- glmer('cbind(ints, tots-ints) ~ Treatment*k*Rewarded + (1 + Rewarded|Subject)', data = learning_curve_counts, family = binomial)

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


