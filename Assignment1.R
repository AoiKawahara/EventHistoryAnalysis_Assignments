#### libraries ####
library(haven)
library(survival)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gmodels)
library(reshape2)

#### setup ####
setwd('/Users/aoikawahara/Library/CloudStorage/OneDrive-KULeuven/Leuven/04_Event History Analysis/Assignment/Workplace')
getwd()

data <- read_sav("FFS_ThirdBirthPPF_Assignment.sav")
CrossTable(data$EXP_LINE,data$BIRTH3)

#### 1a ####
# Excel file

#### 1b ####
data <- data %>%
  group_by(EXP_LINE) %>%
  mutate(ht_obs = mean(BIRTH3, na.rm = TRUE)) %>%
  ungroup()
ggplot(data, aes(EXP_LINE,ht_obs)) + geom_line() +  theme_bw()

data$oddsht_obs <- data$ht_obs/(1-data$ht_obs)
data$logitht_obs <- log(data$ht_obs/(1-data$ht_obs))


#### 1c #####
# Excel file

#### 2a ####
# constant specification
glm1 = glm(BIRTH3 ~ 1, data=data, family=binomial(link = "logit"))
summary(glm1)

data$predl.glm1 = predict.glm(glm1)
data$predo.glm1 = exp(data$predl.glm1)
data$predp.glm1 = data$predo.glm1/(1+data$predo.glm1)

ggplot(data, aes(as.numeric(EXP_LINE), predp.glm1)) + geom_line() + theme_bw()


# linear specification
glm2 = glm(BIRTH3 ~ 1 + EXP_LINE, data=data, family=binomial(link = "logit"))
summary(glm2)

data$predl.glm2 = predict.glm(glm2)
data$predo.glm2 = exp(data$predl.glm2)
data$predp.glm2 = data$predo.glm2/(1+data$predo.glm2)

ggplot(data, aes(x = EXP_LINE)) +
  geom_line(aes(y = ht_obs, group = 1), color = "blue") +
  geom_line(aes(y = predp.glm2, group = 1), color = "darkolivegreen") +
  theme_minimal()


# quadratic specification
data$t2 <- as.numeric(data$EXP_LINE)^2

glm3 = glm(BIRTH3 ~ 1 + EXP_LINE + t2, data=data, family=binomial(link = "logit"))
summary(glm3)

data$predl.glm3 = predict.glm(glm3)
data$predo.glm3 = exp(data$predl.glm3)
data$predp.glm3 = data$predo.glm3/(1+data$predo.glm3)

ggplot(data, aes(x = EXP_LINE)) +
  geom_line(aes(y = ht_obs, group = 1), color = "blue") +
  geom_line(aes(y = predp.glm3, group = 1), color = "darkolivegreen") +
  theme_minimal()


# cubic specification
data$t3 <- as.numeric(data$EXP_LINE)^3

glm4 = glm(BIRTH3 ~ 1 + EXP_LINE + t2 + t3, data=data, family=binomial(link = "logit"))
summary(glm4)

data$predl.glm4 = predict.glm(glm4)
data$predo.glm4 = exp(data$predl.glm4)
data$predp.glm4 = data$predo.glm4/(1+data$predo.glm4)

ggplot(data, aes(x = EXP_LINE)) +
  geom_line(aes(y = ht_obs, group = 1), color = "blue") +
  geom_line(aes(y = predp.glm4, group = 1), color = "darkolivegreen") +
  theme_minimal()


# 4th-order polynomial
data$t4 <- as.numeric(data$EXP_LINE)^4

glm5 = glm(BIRTH3 ~ 1 + EXP_LINE + t2 + t3 + t4, data=data, family=binomial(link = "logit"))
summary(glm5)

data$predl.glm5 = predict.glm(glm5)
data$predo.glm5 = exp(data$predl.glm5)
data$predp.glm5 = data$predo.glm5/(1+data$predo.glm5)

ggplot(data, aes(x = EXP_LINE)) +
  geom_line(aes(y = ht_obs, group = 1), color = "blue") +
  geom_line(aes(y = predp.glm5, group = 1), color = "darkolivegreen") +
  theme_minimal()


# general specification
glm6 = glm(BIRTH3 ~ as.factor(EXP_LINE), data=data, family=binomial(link = "logit"))
summary(glm6)

data$predl.glm6 = predict.glm(glm6)
data$predo.glm6 = exp(data$predl.glm6)
data$predp.glm6 = data$predo.glm6/(1+data$predo.glm6)

ggplot(data, aes(x = EXP_LINE)) +
  geom_line(aes(y = ht_obs, group = 1), color = "blue") +
  geom_line(aes(y = predp.glm6, group = 1), color = "darkolivegreen") +
  theme_minimal()

# step function - 分岐点調整中
crosstab = table(data$EXP_LINE,data$EXP_CAT4)
print(crosstab)

glm7 = glm(BIRTH3 ~ 1 + t4 + as.factor(EXP_CAT4) , data=data, family=binomial(link = "logit"))
summary(glm7)

data$predl.glm7 = predict.glm(glm7)
data$predo.glm7 = exp(data$predl.glm7) 
data$predp.glm7 = data$predo.glm7/(1+data$predo.glm7)

ggplot(data, aes(x = EXP_LINE)) +
  geom_line(aes(y = ht_obs, group = 1), color = "blue") +  # Connected line for ht_obs
  geom_line(aes(y = predp.glm7, group = 1), color = "red") +   # Connected line for ht_cte
  theme_minimal()

# plot all the estimated models
ggplot(data, aes(as.numeric(EXP_LINE), predp.glm1)) + geom_line() + theme_bw() + 
  geom_line(aes(y = predp.glm2), color = "blue") + 
  geom_line(aes(y = predp.glm3), color = "red") + 
  geom_line(aes(y = predp.glm4), color = "darkolivegreen") +
  geom_line(aes(y = predp.glm5), color = "orange") + 
  geom_line(aes(y = predp.glm6), color = "black")


#### 2b ####

# constant
test_glm1glm6 <- anova(glm1,glm6, test = "LRT")
print(test_glm1glm6)

# linear
test_glm1glm2 <- anova(glm1,glm2, test = "LRT")
print(test_glm1glm2)

test_glm2glm6 <- anova(glm2,glm6, test = "LRT")
print(test_glm2glm6)

# quadratic
test_glm2glm3 <- anova(glm2,glm3, test = "LRT")
print(test_glm2glm3)

test_glm3glm6 <- anova(glm3,glm6, test = "LRT")
print(test_glm3glm6)

# cubic
test_glm3glm4 <- anova(glm3,glm4, test = "LRT")
print(test_glm3glm4)

test_glm4glm6 <- anova(glm4,glm6, test = "LRT")
print(test_glm4glm6)

# 4th order
test_glm4glm5 <- anova(glm4,glm5, test = "LRT")
print(test_glm4glm5)

test_glm5glm6 <- anova(glm5,glm6, test = "LRT")
print(test_glm5glm6)

#### 2c ####
# Calculate mean of BIRTH3 by exposure
data <- data %>%
  group_by(AGEKID2) %>%
  mutate(ht_agekid2 = mean(BIRTH3, na.rm = TRUE)) %>%
  ungroup()

ggplot(data, aes(AGEKID2,ht_agekid2)) + geom_line() +  theme_bw()

data$a2 = data$AGEKID2^2
data$a3 = data$AGEKID2^3
data$a4 = data$AGEKID2^4

# 4th-order polynomial model
glm5age = glm(BIRTH3 ~ 1 + EXP_LINE + t2 + t3 + t4 + AGEKID2 + a2 + a3 + a4, data=data, family=binomial(link = "logit"))
summary(glm5age)

data$predl.glm5age = predict.glm(glm5age)
data$predo.glm5age = exp(data$predl.glm5age)
data$predp.glm5age = data$predo.glm5age/(1+data$predo.glm5age)

test_glm5glm5age <- anova(glm5,glm5age, test = "LRT")
print(test_glm5glm5age)

ggplot(data, aes(as.numeric(EXP_LINE), predp.glm5age, linetype=as.factor(AGEKID2))) + geom_line() + theme_bw()

#### 2d ####
data$ASTATUS4 <- factor(data$ASTATUS4,levels=c(1,2,3,4,5), 
                            labels=c("Education","Full-time employed","Part-time employed","Unemployed","Other"))

data.substatus <- subset(data, as.numeric(ASTATUS4) != "NA")
table(data.substatus$ASTATUS4, useNA="always")

glm5sta = glm(BIRTH3 ~ 1 + EXP_LINE + t2 + t3 + t4 + ASTATUS4 , data=data.substatus, family=binomial(link = "logit"))
summary(glm5sta)

# likelihood ratio test


#### 2e ####

