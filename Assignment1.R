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

data.glm <- data %>%
  group_by(EXP_LINE) %>%
  summarize(predp.glm1 = first(predp.glm1))

ggplot(data, aes(as.numeric(EXP_LINE), predp.glm1)) +
  geom_line(aes(y = ht_obs, group = 1), color = "blue") +
  geom_line(aes(y = predp.glm1, group = 1), color = "darkolivegreen") +
  theme_bw()


# linear specification
glm2 = glm(BIRTH3 ~ 1 + EXP_LINE, data=data, family=binomial(link = "logit"))
summary(glm2)

data$predl.glm2 = predict.glm(glm2)
data$predo.glm2 = exp(data$predl.glm2)
data$predp.glm2 = data$predo.glm2/(1+data$predo.glm2)

data.glm <- data %>%
  group_by(EXP_LINE) %>%
  summarize(predp.glm2 = first(predp.glm2))


ggplot(data, aes(x = EXP_LINE)) +
  geom_line(aes(y = ht_obs, group = 1), color = "blue") +
  geom_line(aes(y = predp.glm2, group = 1), color = "darkolivegreen") +
  theme_minimal()


# quadratic specification
glm3 = glm(BIRTH3 ~ 1 + EXP_LINE + EXP_QUAD, data=data, family=binomial(link = "logit"))
summary(glm3)

data$predl.glm3 = predict.glm(glm3)
data$predo.glm3 = exp(data$predl.glm3)
data$predp.glm3 = data$predo.glm3/(1+data$predo.glm3)

data.glm <- data %>%
  group_by(EXP_LINE) %>%
  summarize(predp.glm3 = first(predp.glm3))

ggplot(data, aes(x = EXP_LINE)) +
  geom_line(aes(y = ht_obs, group = 1), color = "blue") +
  geom_line(aes(y = predp.glm3, group = 1), color = "darkolivegreen") +
  theme_minimal()


# cubic specification
data$EXP_CUBE <- as.numeric(data$EXP_LINE)^3

glm4 = glm(BIRTH3 ~ 1 + EXP_LINE + EXP_QUAD + EXP_CUBE, data=data, family=binomial(link = "logit"))
summary(glm4)

data$predl.glm4 = predict.glm(glm4)
data$predo.glm4 = exp(data$predl.glm4)
data$predp.glm4 = data$predo.glm4/(1+data$predo.glm4)

data.glm <- data %>%
  group_by(EXP_LINE) %>%
  summarize(predp.glm4 = first(predp.glm4))

ggplot(data, aes(x = EXP_LINE)) +
  geom_line(aes(y = ht_obs, group = 1), color = "blue") +
  geom_line(aes(y = predp.glm4, group = 1), color = "darkolivegreen") +
  theme_minimal()


# 4th-order polynomial
data$EXP_4TH <- as.numeric(data$EXP_LINE)^4

glm5 = glm(BIRTH3 ~ 1 + EXP_LINE + EXP_QUAD + EXP_CUBE + EXP_4TH, data=data, family=binomial(link = "logit"))
summary(glm5)

data$predl.glm5 = predict.glm(glm5)
data$predo.glm5 = exp(data$predl.glm5)
data$predp.glm5 = data$predo.glm5/(1+data$predo.glm5)

data.glm <- data %>%
  group_by(EXP_LINE) %>%
  summarize(predp.glm5 = first(predp.glm5))

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

data.glm <- data %>%
  group_by(EXP_LINE) %>%
  summarize(predp.glm6 = first(predp.glm6))

ggplot(data, aes(x = EXP_LINE)) +
  geom_line(aes(y = predp.glm6, group = 1), color = "darkolivegreen") +
  theme_minimal()


# step function
crosstab = table(data$EXP_LINE,data$EXP_CAT4)
print(crosstab)

glm7 = glm(BIRTH3 ~ 1 + as.factor(EXP_CAT4) + EXP_LINE + EXP_QUAD + EXP_CUBE + EXP_4TH, data=data, family=binomial(link = "logit"))
summary(glm7)

data$predl.glm7 = predict.glm(glm7)
data$predo.glm7 = exp(data$predl.glm7) 
data$predp.glm7 = data$predo.glm7/(1+data$predo.glm7)

data.glm <- data %>%
  group_by(EXP_LINE) %>%
  summarize(predp.glm7 = first(predp.glm7))

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

# step function

test_glm6glm7 <- anova(glm6, glm7, test = "LRT")
print(test_glm6glm7)

test_glm5glm7 <- anova(glm5, glm7, test = "LRT")
print(test_glm5glm7)

#### 2c ####
table(data$AGEKID2, useNA="always")

# Calculate mean of BIRTH3 by exposure
data <- data %>%
  group_by(AGEKID2) %>%
  mutate(ht_agekid2 = mean(BIRTH3, na.rm = TRUE)) %>%
  ungroup()

ggplot(data, aes(AGEKID2,ht_agekid2)) +
  geom_line() +
  theme_bw()

data$a2 = data$AGEKID2^2
data$a3 = data$AGEKID2^3
data$a4 = data$AGEKID2^4


# 4th-order polynomial model
# age as linear
glm5age1 = glm(BIRTH3 ~ 1 + EXP_LINE + EXP_QUAD + EXP_CUBE + EXP_4TH + AGEKID2, data=data, family=binomial(link = "logit"))
summary(glm5age1)

data$predl.glm5age1 = predict.glm(glm5age1)
data$predo.glm5age1 = exp(data$predl.glm5age1)
data$predp.glm5age1 = data$predo.glm5age1/(1+data$predo.glm5age1)

test_glm5glm5age1 <- anova(glm5,glm5age1, test = "LRT")
print(test_glm5glm5age1)

ggplot(data, aes(as.numeric(EXP_LINE), predp.glm5age1, color=as.factor(AGEKID2))) +
  geom_line() +
  theme_bw()

# age as quadratic
glm5age2 = glm(BIRTH3 ~ 1 + EXP_LINE + EXP_QUAD + EXP_CUBE + EXP_4TH + AGEKID2 + a2, data=data, family=binomial(link = "logit"))
summary(glm5age2)

data$predl.glm5age2 = predict.glm(glm5age2)
data$predo.glm5age2 = exp(data$predl.glm5age2)
data$predp.glm5age2 = data$predo.glm5age2/(1+data$predo.glm5age2)

test_glm5age1glm5age2 <- anova(glm5age1,glm5age2, test = "LRT")
print(test_glm5age1glm5age2)

ggplot(data, aes(as.numeric(EXP_LINE), predp.glm5age2, color=as.factor(AGEKID2))) +
  geom_line() +
  theme_bw()


# age as cubic
glm5age3 = glm(BIRTH3 ~ 1 + EXP_LINE + EXP_QUAD + EXP_CUBE + EXP_4TH + AGEKID2 + a2 + a3, data=data, family=binomial(link = "logit"))
summary(glm5age3)

data$predl.glm5age3 = predict.glm(glm5age3)
data$predo.glm5age3 = exp(data$predl.glm5age3)
data$predp.glm5age3 = data$predo.glm5age3/(1+data$predo.glm5age3)

test_glm5age1glm5age3 <- anova(glm5age1,glm5age3, test = "LRT")
print(test_glm5age1glm5age3)

# age as 4th-order
glm5age4 = glm(BIRTH3 ~ 1 + EXP_LINE + EXP_QUAD + EXP_CUBE + EXP_4TH + AGEKID2 + a2 + a3 + a4, data=data, family=binomial(link = "logit"))
summary(glm5age4)

data$predl.glm5age4 = predict.glm(glm5age4)
data$predo.glm5age4 = exp(data$predl.glm5age4)
data$predp.glm5age4 = data$predo.glm5age4/(1+data$predo.glm5age4)

test_glm5age1glm5age4 <- anova(glm5age1,glm5age4, test = "LRT")
print(test_glm5age1glm5age4)

# age as categorical (= general specification)
glm5age5 = glm(BIRTH3 ~ 1 + EXP_LINE + EXP_QUAD + EXP_CUBE + EXP_4TH + as.factor(AGEKID2), data=data, family=binomial(link = "logit"))
summary(glm5age5)

data$predl.glm5age5 = predict.glm(glm5age5)
data$predo.glm5age5 = exp(data$predl.glm5age5)
data$predp.glm5age5 = data$predo.glm5age5/(1+data$predo.glm5age5)

test_glm5age1glm5age5 <- anova(glm5age1,glm5age5, test = "LRT")
print(test_glm5age1glm5age5)

ggplot(data, aes(as.numeric(EXP_LINE), predp.glm5age5, color=as.factor(AGEKID2))) +
  geom_line() +
  ylim(0,0.40) +
  theme_bw()


test_glm5age1glm5age5 <- anova(glm5age1, glm5age5, test = "LRT")
print(test_glm5age1glm5age5)

## Suggest the linear model for the effect of age

#### 2d ####
data$ASTATUS4 <- factor(data$ASTATUS4,levels=c(1,2,3,4,5), 
                            labels=c("Education","Full-time employed","Part-time employed","Unemployed","Other"))

data.substatus <- subset(data, as.numeric(ASTATUS4) != "NA")
table(data.substatus$ASTATUS4, useNA="always")

glm5sta = glm(BIRTH3 ~ 1 + EXP_LINE + EXP_QUAD + EXP_CUBE + EXP_4TH + ASTATUS4 , data=data.substatus, family=binomial(link = "logit"))
summary(glm5sta)

data.substatus$predl.glm5sta = predict.glm(glm5sta)
data.substatus$predo.glm5sta = exp(data.substatus$predl.glm5sta)
data.substatus$predp.glm5sta = data.substatus$predo.glm5sta/(1+data.substatus$predo.glm5sta)

ggplot(data.substatus, aes(as.numeric(EXP_LINE), predp.glm5sta, linetype=as.factor(ASTATUS4))) +
  geom_line() +
  theme_bw()


# Since LRT assumes the people in the models are the same,
# the number of cases must be the same across the models.


# 4th-order polynomial specification with subset data
glm5sub = glm(BIRTH3 ~ 1 + EXP_LINE + EXP_QUAD + EXP_CUBE + EXP_4TH, data=data.substatus, family=binomial(link = "logit"))
summary(glm5sub)

test_glm5subglm5sta <- anova(glm5sub,glm5sta, test = "LRT")
print(test_glm5subglm5sta)


#### 2e ####
# model with both age and status
glm5netfull = glm(BIRTH3 ~ 1 + EXP_LINE + EXP_QUAD + EXP_CUBE + EXP_4TH +
                 AGEKID2 + ASTATUS4,
                 data=data.substatus, family=binomial(link = "logit"))
summary(glm5netfull)

# 4th-order polynomial specification with subset data only with age
glm5netage = glm(BIRTH3 ~ 1 + EXP_LINE + EXP_QUAD + EXP_CUBE + EXP_4TH + AGEKID2,
                 data=data.substatus, family=binomial(link = "logit"))
summary(glm5netage)


# likelihood ratio test
test_glm5netfullglm5netage <- anova(glm5netage,glm5netfull, test = "LRT")
print(test_glm5netfullglm5netage)

test_glm5netfullglm5sta <- anova(glm5sta,glm5netfull, test = "LRT")
print(test_glm5netfullglm5sta)


#### 2f ####
glm5netfull$devres <- residuals(glm5netfull, "deviance")

plot(glm5netfull$devres,
     xlab = "Person-Period",
     ylab = "Deviance Residuals")

glm5netfull$ssdevres <- tapply(glm5netfull$devres^2, data.substatus$NR, sum)

plot(glm5netfull$ssdevres,
     xlab = "Identification Number",
     ylab = "Sum of Squared Residuals")
text(x = which(glm5netfull$ssdevres >= 8),
     y = glm5netfull$ssdevres[glm5netfull$ssdevres >= 8],
     labels = names(glm5netfull$ssdevres)[glm5netfull$ssdevres >= 8],
     pos = 3,
     cex = 0.8,
     col = "red")
glm5netfull$ssdevres

#### 2g ####
glm5netfull_cloglog = glm(BIRTH3 ~ 1 + EXP_LINE + EXP_QUAD + EXP_CUBE + EXP_4TH +
                            AGEKID2 + ASTATUS4,
                          data=data.substatus, family=binomial(link = "cloglog"))
summary(glm5netfull_cloglog)
