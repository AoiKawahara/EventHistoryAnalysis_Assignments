#### libraries ####
library(haven)
library(survival)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gmodels)
library(reshape2)
library(nnet)

#### setup ####
setwd('/Users/aoikawahara/Library/CloudStorage/OneDrive-KULeuven/Leuven/04_Event History Analysis/Assignment/Workplace')
getwd()

data <- read_sav("QASS General & Cause-specific Mortality.sav")
options(max.print=20000)

#### 1a ####
table(data$MONTH_CONT,data$DEATH)

#### 1b ####
data$ACTIVITY_fac <- factor(data$ACTIVITY,
                        levels=c(1,2,3,4,5),
                        labels=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled", "5.Unknown"))
table(data$ACTIVITY_fac, useNA="always")

KMact <- survfit(Surv(data$MONTH_CONT, data$DEATH) ~ 1 + ACTIVITY, data = data)
summary(KMact)

# survivor function
plot(KMact, xlab="Time", ylab="S(t)", conf.int=F, mark.time=F, col=c(2,3,4,5), lty=1, lwd=1, ylim = c(0.85,1))
legend(x = "bottomleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

# cumulative hazard function
plot(KMact, fun="cumhaz", xlab="Time", ylab="H(t)=-ln[S(t)]", conf.int=F, col=c(2,3,4,5), ylim=c(0,0.15))
legend(x = "topleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

# log cumulative hazard function
logcumhaz <- function(y) {
  to.return <- log(-log(y))
  return(to.return)
}

plot(KMact, fun=logcumhaz, , xlab="Time", ylab="Log[H(t)]", ylim=c(-8,-2), conf.int=F, mark.time=F, col=c(2,3,4,5))
legend(x = "bottomright",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

#### 1c ####
data.subacts <- subset(data, as.numeric(ACTIVITY) != "NA")
data.subacts$ACTIVITY

coxAct <- coxph(Surv(data$MONTH_CONT, data$DEATH) ~ 1 + as.factor(ACTIVITY), data = data)
summary(coxAct)

# recover cumulative hazard function
plot(survfit(coxAct, newdata = data.frame(ACTIVITY = 1)),
     conf.int = FALSE,
     col="red", fun="cumhaz", ylim=c(0,0.15))
lines(survfit(coxAct, newdata = data.frame(ACTIVITY = 2)),
      conf.int = FALSE,
      col="yellow", fun="cumhaz")
lines(survfit(coxAct, newdata = data.frame(ACTIVITY = 3)),
      conf.int = FALSE,
      col="lightgreen", fun="cumhaz")
lines(survfit(coxAct, newdata = data.frame(ACTIVITY = 4)),
      conf.int = FALSE,
      col="blue", fun="cumhaz")
legend(x = "topleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))

# log cumulative hazard function
plot(survfit(coxAct, newdata = data.frame(ACTIVITY = 1)),
     conf.int = FALSE,
     col="red", fun=logcumhaz, ylim=c(-10,-2))
lines(survfit(coxAct, newdata = data.frame(ACTIVITY = 2)),
      conf.int = FALSE,
      col="yellow", fun=logcumhaz)
lines(survfit(coxAct, newdata = data.frame(ACTIVITY = 3)),
      conf.int = FALSE,
      col="lightgreen", fun=logcumhaz)
lines(survfit(coxAct, newdata = data.frame(ACTIVITY = 4)),
      conf.int = FALSE,
      col="blue", fun=logcumhaz)
legend(x = "bottomright",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))

cox <- coxph(Surv(data.subacts$MONTH_CONT, data.subacts$DEATH) ~ 1, data = data.subacts)
summary(cox)

test_coxcoxAct <- anova(coxAct, cox, test = "LRT")
print(test_coxcoxAct)

#### 1d ####
data$CAUSE_fac <- factor(data$CAUSE,
                         levels=c(0,1,2,3,4),
                         labels=c("0.Censored", "1.Cancer", "2.Circular", "3.Respiratory", "4.Other"))
table(data$MONTH_CONT, data$CAUSE_fac)

#### 1e ####

#### CANCER ####
KMactCan <- survfit(Surv(data$MONTH_CONT, data$CANCER) ~ 1 + ACTIVITY, data = data)

plot(KMactCan, xlab="Time", ylab="S(t)", conf.int=F, mark.time=F, col=c(2,3,4,5), lty=1, lwd=1, ylim = c(0.96,1))
legend(x = "bottomleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactCan, fun="cumhaz", xlab="Time", ylab="H(t)=-ln[S(t)]", conf.int=F, col=c(2,3,4,5), ylim=c(0,0.04))
legend(x = "topleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactCan, fun=logcumhaz, , xlab="Time", ylab="Log[H(t)]", ylim=c(-10,-3), conf.int=F, mark.time=F, col=c(2,3,4,5))
legend(x = "bottomright",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

coxActCan<- coxph(Surv(MONTH_CONT, CANCER) ~ 1 + as.factor(ACTIVITY), data = data.subacts)
summary(coxActCan)

plot(survfit(coxActCan, newdata = data.frame(ACTIVITY = 1)),
     conf.int = FALSE,
     col="red", fun="cumhaz", ylim=c(0,0.05))
lines(survfit(coxActCan, newdata = data.frame(ACTIVITY = 2)),
      conf.int = FALSE,
      col="yellow", fun="cumhaz")
lines(survfit(coxActCan, newdata = data.frame(ACTIVITY = 3)),
      conf.int = FALSE,
      col="lightgreen", fun="cumhaz")
lines(survfit(coxActCan, newdata = data.frame(ACTIVITY = 4)),
      conf.int = FALSE,
      col="blue", fun="cumhaz")
legend(x = "topleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))

plot(survfit(coxActCan, newdata = data.frame(ACTIVITY = 1)),
     conf.int = FALSE,
     col="red", fun=logcumhaz, ylim=c(-10,-2))
lines(survfit(coxActCan, newdata = data.frame(ACTIVITY = 2)),
      conf.int = FALSE,
      col="yellow", fun=logcumhaz)
lines(survfit(coxActCan, newdata = data.frame(ACTIVITY = 3)),
      conf.int = FALSE,
      col="lightgreen", fun=logcumhaz)
lines(survfit(coxActCan, newdata = data.frame(ACTIVITY = 4)),
      conf.int = FALSE,
      col="blue", fun=logcumhaz)
legend(x = "bottomright",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))

coxCan <- coxph(Surv(MONTH_CONT, CANCER) ~ 1, data = data.subacts)
summary(coxCan)

test_coxCancoxActCan <- anova(coxActCan, coxCan, test = "LRT")
print(test_coxCancoxActCan)

#### CIRCULAR ####
KMactCir <- survfit(Surv(data$MONTH_CONT, data$CIRCULAR) ~ 1 + ACTIVITY, data = data)

plot(KMactCir, xlab="Time", ylab="S(t)", conf.int=F, mark.time=F, col=c(2,3,4,5), lty=1, lwd=1, ylim = c(0.95,1))
legend(x = "bottomleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactCir, fun="cumhaz", xlab="Time", ylab="H(t)=-ln[S(t)]", conf.int=F, col=c(2,3,4,5), ylim=c(0,0.055))
legend(x = "topleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactCir, fun=logcumhaz, , xlab="Time", ylab="Log[H(t)]", ylim=c(-10,-3), conf.int=F, mark.time=F, col=c(2,3,4,5))
legend(x = "bottomright",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

coxActCir <- coxph(Surv(MONTH_CONT, CIRCULAR) ~ 1 + as.factor(ACTIVITY), data = data.subacts)
summary(coxActCir)

plot(survfit(coxActCir, newdata = data.frame(ACTIVITY = 1)),
     conf.int = FALSE,
     col="red", fun="cumhaz", ylim=c(0,0.06))
lines(survfit(coxActCir, newdata = data.frame(ACTIVITY = 2)),
      conf.int = FALSE,
      col="yellow", fun="cumhaz")
lines(survfit(coxActCir, newdata = data.frame(ACTIVITY = 3)),
      conf.int = FALSE,
      col="lightgreen", fun="cumhaz")
lines(survfit(coxActCir, newdata = data.frame(ACTIVITY = 4)),
      conf.int = FALSE,
      col="blue", fun="cumhaz")
legend(x = "topleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))

plot(survfit(coxActCir, newdata = data.frame(ACTIVITY = 1)),
     conf.int = FALSE,
     col="red", fun=logcumhaz, ylim=c(-10,-2))
lines(survfit(coxActCir, newdata = data.frame(ACTIVITY = 2)),
      conf.int = FALSE,
      col="yellow", fun=logcumhaz)
lines(survfit(coxActCir, newdata = data.frame(ACTIVITY = 3)),
      conf.int = FALSE,
      col="lightgreen", fun=logcumhaz)
lines(survfit(coxActCir, newdata = data.frame(ACTIVITY = 4)),
      conf.int = FALSE,
      col="blue", fun=logcumhaz)
legend(x = "bottomright",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))

coxCir <- coxph(Surv(MONTH_CONT, CIRCULAR) ~ 1, data = data.subacts)
summary(coxCir)

test_coxCircoxActCir <- anova(coxCir, coxActCir, test = "LRT")
print(test_coxCircoxActCir)

#### RESPIRATORY ####
KMactRes <- survfit(Surv(data$MONTH_CONT, data$RESPIRATORY) ~ 1 + ACTIVITY, data = data)

plot(KMactRes, xlab="Time", ylab="S(t)", conf.int=F, mark.time=F, col=c(2,3,4,5), lty=1, lwd=1, ylim = c(0.98,1))
legend(x = "bottomleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactRes, fun="cumhaz", xlab="Time", ylab="H(t)=-ln[S(t)]", conf.int=F, col=c(2,3,4,5), ylim=c(0,0.025))
legend(x = "topleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactRes, fun=logcumhaz, , xlab="Time", ylab="Log[H(t)]", ylim=c(-10,-3), conf.int=F, mark.time=F, col=c(2,3,4,5))
legend(x = "bottomright",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

coxActRes <- coxph(Surv(MONTH_CONT, RESPIRATORY) ~ 1 + as.factor(ACTIVITY), data = data.subacts)
summary(coxActRes)

plot(survfit(coxActRes, newdata = data.frame(ACTIVITY = 1)),
     conf.int = FALSE,
     col="red", fun="cumhaz", ylim=c(0,0.025))
lines(survfit(coxActRes, newdata = data.frame(ACTIVITY = 2)),
      conf.int = FALSE,
      col="yellow", fun="cumhaz")
lines(survfit(coxActRes, newdata = data.frame(ACTIVITY = 3)),
      conf.int = FALSE,
      col="lightgreen", fun="cumhaz")
lines(survfit(coxActRes, newdata = data.frame(ACTIVITY = 4)),
      conf.int = FALSE,
      col="blue", fun="cumhaz")
legend(x = "topleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))

plot(survfit(coxActRes, newdata = data.frame(ACTIVITY = 1)),
     conf.int = FALSE,
     col="red", fun=logcumhaz, ylim=c(-11,-4))
lines(survfit(coxActRes, newdata = data.frame(ACTIVITY = 2)),
      conf.int = FALSE,
      col="yellow", fun=logcumhaz)
lines(survfit(coxActRes, newdata = data.frame(ACTIVITY = 3)),
      conf.int = FALSE,
      col="lightgreen", fun=logcumhaz)
lines(survfit(coxActRes, newdata = data.frame(ACTIVITY = 4)),
      conf.int = FALSE,
      col="blue", fun=logcumhaz)
legend(x = "bottomright",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))

coxRes <- coxph(Surv(MONTH_CONT, RESPIRATORY) ~ 1, data = data.subacts)
summary(coxRes)

test_coxRescoxActRes <- anova(coxRes, coxActRes, test = "LRT")
print(test_coxRescoxActRes)


#### OTHER ####
KMactOth <- survfit(Surv(data$MONTH_CONT, data$OTHER) ~ 1 + ACTIVITY, data = data)

plot(KMactOth, xlab="Time", ylab="S(t)", conf.int=F, mark.time=F, col=c(2,3,4,5), lty=1, lwd=1, ylim = c(0.965,1))
legend(x = "bottomleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactOth, fun="cumhaz", xlab="Time", ylab="H(t)=-ln[S(t)]", conf.int=F, col=c(2,3,4,5), ylim=c(0,0.035))
legend(x = "topleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactOth, fun=logcumhaz, , xlab="Time", ylab="Log[H(t)]", ylim=c(-10,-3), conf.int=F, mark.time=F, col=c(2,3,4,5))
legend(x = "bottomright",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c(2,3,4,5),
       lty=c(1,1))

coxActOth <- coxph(Surv(MONTH_CONT, OTHER) ~ 1 + as.factor(ACTIVITY), data = data.subacts)
summary(coxActOth)

plot(survfit(coxActOth, newdata = data.frame(ACTIVITY = 1)),
     conf.int = FALSE,
     col="red", fun="cumhaz", ylim=c(0,0.025))
lines(survfit(coxActOth, newdata = data.frame(ACTIVITY = 2)),
      conf.int = FALSE,
      col="yellow", fun="cumhaz")
lines(survfit(coxActOth, newdata = data.frame(ACTIVITY = 3)),
      conf.int = FALSE,
      col="lightgreen", fun="cumhaz")
lines(survfit(coxActOth, newdata = data.frame(ACTIVITY = 4)),
      conf.int = FALSE,
      col="blue", fun="cumhaz")
legend(x = "topleft",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))

plot(survfit(coxActOth, newdata = data.frame(ACTIVITY = 1)),
     conf.int = FALSE,
     col="red", fun=logcumhaz, ylim=c(-10,-3))
lines(survfit(coxActOth, newdata = data.frame(ACTIVITY = 2)),
      conf.int = FALSE,
      col="yellow", fun=logcumhaz)
lines(survfit(coxActOth, newdata = data.frame(ACTIVITY = 3)),
      conf.int = FALSE,
      col="lightgreen", fun=logcumhaz)
lines(survfit(coxActOth, newdata = data.frame(ACTIVITY = 4)),
      conf.int = FALSE,
      col="blue", fun=logcumhaz)
legend(x = "bottomright",
       legend=c("1.Employed", "2.Unemployed", "3.Retired", "4.Disabled"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))

coxOth <- coxph(Surv(MONTH_CONT, OTHER) ~ 1, data = data.subacts)
summary(coxOth)

test_coxOthcoxActOth <- anova(coxOth, coxActOth, test = "LRT")
print(test_coxOthcoxActOth)


#### 2a ####
data.ppf <- read_sav("QASS General & Cause-specific Mortality PPF.sav")

data.ppf$CAUSE_fac <- factor(data.ppf$CAUSE,
                         levels=c(0,1,2,3,4),
                         labels=c("0.Survival or emigration", "1.Cancer", "2.Circular", "3.Respiratory", "4.Other"))
table(data.ppf$EXPOSURE, data.ppf$CAUSE_fac)

#### 2b 2c ####
glmcan = glm(CANCER ~ 1 + EXPOSURE, data=data.ppf, family=binomial(link = "logit"))
summary(glmcan)
glmActcan = glm(CANCER ~ 1 + EXPOSURE + as.factor(ACTIVITY), data=data.ppf, family=binomial(link = "logit"))
summary(glmActcan)

test_glmcanglmActcan <- anova(glmcan,glmActcan, test = "LRT")
print(test_glmcanglmActcan)

coxCanPpf <- coxph(Surv(EXPOSURE, CANCER) ~ 1 + as.factor(ACTIVITY), data = data.ppf)
summary(coxCanPpf)

glmcir = glm(CIRCULAR ~ 1 + EXPOSURE, data=data.ppf, family=binomial(link = "logit"))
summary(glmcir)
glmActcir = glm(CIRCULAR ~ 1 + EXPOSURE + as.factor(ACTIVITY), data=data.ppf, family=binomial(link = "logit"))
summary(glmActcir)

test_glmcirglmActcir <- anova(glmcir,glmActcir, test = "LRT")
print(test_glmcirglmActcir)

coxCirPpf <- coxph(Surv(EXPOSURE, CIRCULAR) ~ 1 + as.factor(ACTIVITY), data = data.ppf)
summary(coxCirPpf)

glmres = glm(RESPIRATORY ~ 1 + EXPOSURE, data=data.ppf, family=binomial(link = "logit"))
summary(glmres)
glmActres = glm(RESPIRATORY ~ 1 + EXPOSURE + as.factor(ACTIVITY), data=data.ppf, family=binomial(link = "logit"))
summary(glmActres)

test_glmresglmActres <- anova(glmres,glmActres, test = "LRT")
print(test_glmresglmActres)

coxResPpf <- coxph(Surv(EXPOSURE, RESPIRATORY) ~ 1 + as.factor(ACTIVITY), data = data.ppf)
summary(coxResPpf)

glmoth = glm(OTHER ~ 1 + EXPOSURE, data=data.ppf, family=binomial(link = "logit"))
summary(glmoth)
glmActoth = glm(OTHER ~ 1 + EXPOSURE + as.factor(ACTIVITY), data=data.ppf, family=binomial(link = "logit"))
summary(glmActoth)

test_glmothglmActoth <- anova(glmoth,glmActoth, test = "LRT")
print(test_glmothglmActoth)

coxOthPpf <- coxph(Surv(EXPOSURE, OTHER) ~ 1 + as.factor(ACTIVITY), data = data.ppf)
summary(coxOthPpf)
