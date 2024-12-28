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

#### 1a ####
options(max.print=20000)
table <- table(data$MONTH_CONT,data$DEATH)

#### 1b ####
table(data$ACTIVITY, useNA="always")

KMact <- survfit(Surv(data$MONTH_CONT, data$DEATH) ~ 1 + ACTIVITY, data = data)
options(max.print=2000)
summary(KMact)

# survivor function
plot(KMact, xlab="Time", ylab="S(t)", conf.int=F, mark.time=F, col=c(2,3,4,5), lty=1, lwd=1, ylim = c(0.85,1))
legend(x = "bottomleft",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

# cumulative hazard function
plot(KMact, fun="cumhaz", xlab="Time", ylab="H(t)=-ln[S(t)]", conf.int=F, col=c(2,3,4,5), ylim=c(0,0.15))
legend(x = "topleft",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

# log cumulative hazard function
logcumhaz <- function(y) {
  to.return <- log(-log(y))
  return(to.return)
}

plot(KMact, fun=logcumhaz, , xlab="Time", ylab="Log[H(t)]", ylim=c(-8,-2), conf.int=F, mark.time=F, col=c(2,3,4,5))
legend(x = "bottomright",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

#### 1c ####
coxAct <- coxph(Surv(data$MONTH_CONT, data$DEATH) ~ ACTIVITY, data = data)
summary(coxAct)

# cumulative hazard function
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
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
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
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))

# likelihood ratio test
cox <- coxph(Surv(data$MONTH_CONT, data$DEATH) ~ 1, data = data)
test_coxcoxAct <- anova(cox,coxAct, test = "LRT")
print(test_coxcoxAct)

#### 1d ####
table(data$MONTH_CONT, data$CAUSE)

#### CANCER ####
KMactCan <- survfit(Surv(data$MONTH_CONT, data$CANCER) ~ 1 + ACTIVITY, data = data)

plot(KMactCan, xlab="Time", ylab="S(t)", conf.int=F, mark.time=F, col=c(2,3,4,5), lty=1, lwd=1, ylim = c(0.96,1))
legend(x = "bottomleft",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactCan, fun="cumhaz", xlab="Time", ylab="H(t)=-ln[S(t)]", conf.int=F, col=c(2,3,4,5), ylim=c(0,0.04))
legend(x = "topleft",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactCan, fun=logcumhaz, , xlab="Time", ylab="Log[H(t)]", ylim=c(-10,-3), conf.int=F, mark.time=F, col=c(2,3,4,5))
legend(x = "bottomright",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

coxActCan<- coxph(Surv(data$MONTH_CONT, data$CANCER) ~ ACTIVITY, data = data)
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
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
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
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))

#### CIRCULAR ####
KMactCir <- survfit(Surv(data$MONTH_CONT, data$CIRCULAR) ~ 1 + ACTIVITY, data = data)

plot(KMactCir, xlab="Time", ylab="S(t)", conf.int=F, mark.time=F, col=c(2,3,4,5), lty=1, lwd=1, ylim = c(0.95,1))
legend(x = "bottomleft",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactCir, fun="cumhaz", xlab="Time", ylab="H(t)=-ln[S(t)]", conf.int=F, col=c(2,3,4,5), ylim=c(0,0.055))
legend(x = "topleft",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactCir, fun=logcumhaz, , xlab="Time", ylab="Log[H(t)]", ylim=c(-10,-3), conf.int=F, mark.time=F, col=c(2,3,4,5))
legend(x = "bottomright",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

coxActCir <- coxph(Surv(data$MONTH_CONT, data$CIRCULAR) ~ ACTIVITY, data = data)
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
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
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
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))

#### RESPIRATORY ####
KMactRes <- survfit(Surv(data$MONTH_CONT, data$RESPIRATORY) ~ 1 + ACTIVITY, data = data)

plot(KMactRes, xlab="Time", ylab="S(t)", conf.int=F, mark.time=F, col=c(2,3,4,5), lty=1, lwd=1, ylim = c(0.98,1))
legend(x = "bottomleft",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactRes, fun="cumhaz", xlab="Time", ylab="H(t)=-ln[S(t)]", conf.int=F, col=c(2,3,4,5), ylim=c(0,0.025))
legend(x = "topleft",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactRes, fun=logcumhaz, , xlab="Time", ylab="Log[H(t)]", ylim=c(-10,-3), conf.int=F, mark.time=F, col=c(2,3,4,5))
legend(x = "bottomright",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

coxActRes <- coxph(Surv(data$MONTH_CONT, data$RESPIRATORY) ~ ACTIVITY, data = data)
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
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
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
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))

#### OTHER ####
KMactOth <- survfit(Surv(data$MONTH_CONT, data$OTHER) ~ 1 + ACTIVITY, data = data)

plot(KMactOth, xlab="Time", ylab="S(t)", conf.int=F, mark.time=F, col=c(2,3,4,5), lty=1, lwd=1, ylim = c(0.965,1))
legend(x = "bottomleft",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactOth, fun="cumhaz", xlab="Time", ylab="H(t)=-ln[S(t)]", conf.int=F, col=c(2,3,4,5), ylim=c(0,0.035))
legend(x = "topleft",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

plot(KMactOth, fun=logcumhaz, , xlab="Time", ylab="Log[H(t)]", ylim=c(-10,-3), conf.int=F, mark.time=F, col=c(2,3,4,5))
legend(x = "bottomright",
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c(2,3,4,5),
       lty=c(1,1))

coxActOth <- coxph(Surv(data$MONTH_CONT, data$OTHER) ~ ACTIVITY, data = data)
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
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
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
       legend=c("ACTIVITY = 1", "ACTIVITY = 2", "ACTIVITY = 3", "ACTIVITY = 4"),
       col=c("red", "yellow", "lightgreen", "blue"), lty=c(1,1))


#### 1f ####
coxCan <- coxph(Surv(data$MONTH_CONT, data$CANCER) ~ ACTIVITY, data = data)
test_coxCancoxActCan <- anova(coxCan,coxActCan, test = "LRT")
print(test_coxCancoxActCan)

coxCir <- coxph(Surv(data$MONTH_CONT, data$CIRCULAR) ~ ACTIVITY, data = data)
test_coxCircoxActCir <- anova(coxCir,coxActCir, test = "LRT")
print(test_coxCircoxActCir)

coxRes <- coxph(Surv(data$MONTH_CONT, data$RESPIRATORY) ~ ACTIVITY, data = data)
test_coxRescoxActRes <- anova(coxRes,coxActRes, test = "LRT")
print(test_coxRescoxActRes)

coxOth <- coxph(Surv(data$MONTH_CONT, data$OTHER) ~ ACTIVITY, data = data)
test_coxOthcoxActOth <- anova(coxOth,coxActOth, test = "LRT")
print(test_coxOthcoxActOth)





