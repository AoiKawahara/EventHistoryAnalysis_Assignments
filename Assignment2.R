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

data <- read_sav("QASS General & Cause-specific Mortality.sav")

#### 1a ####
table <- table(data$MONTH_CONT,data$DEATH)
table

#### 1b ####
# create actuarial life-table
data$CENSOR <- 1
data$CENSOR[data$DEATH == 1] <- 0
data2 <- select(data, MONTH_CONT, CENSOR)

interval<- rep(dim(data2)[1], 0)             # intervalというdata2の行数と同じ長さのベクトルを生成
for (i in 1:dim(data2)[1]){
  interval[i] = floor(data2$MONTH_CONT[i])}  # MONTH_COUNTの整数部分をintervalに代入

int.event <- rep(0, 35)
int.censor <- rep(0, 35)
for (i in 1:dim(data2)[1]){
  int <- interval[i]
  if (data2$CENSOR[i] == 0) int.event[int] <- int.event[int] + 1
  if (data2$CENSOR[i] == 1) int.censor[int] <- int.censor[int] + 1}  # int.eventとint.censorにintervalごとのeventとcensorの数を代入

interval <- c(1:35) # calculate widths of intervals (= 1 for all intervals)
interval.end <- c(2:36)
interval.w <- interval.end - interval

int.risk <- c(20000, rep(0, 34)) # calculate risk set
for (i in 2:35){
  int.risk[i] <- int.risk[i-1] - (int.event[i-1] + int.censor[i-1])}

phat = (int.event / int.risk)
discrete.h = (phat / interval.w)  # discrete-time hazard (= phat)
actuarial.h = (int.event / (int.risk - (int.censor/2) - (int.event/2))) / interval.w  # actuarial hazard

discrete.s <- c(1 - discrete.h[1], rep(0, 34))  # conditional probability of survival
actuarial.s <- c(1, rep(0,34))                  # conditional probability of survival based on actuarial hazard
for (i in 2:35){
  discrete.s[i] <- (1 - phat[i]) * discrete.s[i-1]
  actuarial.s[i] <- actuarial.s[i-1] * (1 - int.event[i] / (int.risk[i] - int.censor[i]/2))}	

int.table <- cbind(interval, interval.end, int.risk, int.event, int.censor, phat, discrete.s, discrete.h, actuarial.s, actuarial.h)
int.table

attach(data.frame(int.table))

# plot survivor function
s.hat.steps <- stepfun(interval, c(1, actuarial.s))
plot(s.hat.steps, do.points = FALSE, xlim = c(0, 40), 
     ylab = "Actuarial Survival", xlab = "Time", main = "")

#### 1c ####
KM1 <- survfit(Surv(data$MONTH_CONT, data$DEATH) ~ 1, data = data) 
options(max.print=20000)
summary(KM1)

mean(KM1$time)
quantile(KM1$time)

# survivor function S(t)
plot(KM1, conf.int=F, mark.time=T, col=1, lty=1, lwd=1, ylim = c(0.9,1),
     ylab = "S(t)", xlab = "Time", main = "")

# cumulative hazard function -ln(S(t))
plot(KM1, fun="cumhaz", ylim=c(0,0.1), conf.int=F)

#### 1d ####
logcumhaz <- function(y) {
  to.return <- log(-log(y))
  return(to.return)
}
plot(KM1, fun=logcumhaz, conf.int=F)  # log cumulative hazard ln(-ln(S(t)))

#### 1e ####
cox1 <- coxph(Surv(data$MONTH_CONT, data$DEATH) ~ 1, data = data)
summary(cox1)

cox1res <- resid(cox1, type="martingale")
summary(cox1res)

#### 2a ####
plotmartingale <- qplot(AGE, cox1res, data=data, color=factor(DEATH))
plotmartingale +
  scale_color_discrete(name="DEATH", label=c("No", "Yes")) +
  ylab("Martingale Residuals") +
  xlab("AGE") +
  theme_bw()

ggplot(data, aes(AGE, cox1res)) +
  geom_line() +
  theme_bw()

# simple regression of martingale residuals
resmodel = lm(cox1res ~ AGE, data = data)
summary(resmodel)

reslinear = predict(resmodel)

ggplot(data, aes(AGE, reslinear)) +
  geom_line() +
  theme_bw()

# loess smooth of martingale residuals
loessmodel = loess(cox1res ~ data$AGE)
summary(loessmodel)

loessvalues = predict(loessmodel)

ggplot(data, aes(AGE, loessvalues)) +
  geom_line() +
  theme_bw()

#### 2b ####
coxAge <- coxph(Surv(data$MONTH_CONT, data$DEATH) ~ 1 + AGE, data = data)
summary(coxAge)

test_cox1coxAge <- anova(cox1,coxAge, test = "LRT")
print(test_cox1coxAge)

#### 2c ####
coxAge.partial <- resid(coxAge, "schoenfeld")
summary(coxAge.partial)

data.partial <- subset(data, CENSOR!= 1)
data.partial <- data.frame(data.partial, coxAge.partial)
data.partial$RANK <- rank(data.partial$MONTH_CONT, na.last = TRUE, ties.method = "first")

loessmodel_coxAge = loess(coxAge.partial ~ RANK, data = data.partial)
summary(loessmodel_coxAge)
data.partial$loessvalues_coxAge = predict(loessmodel_coxAge)

ggplot(data.partial, aes(x = RANK)) +
  geom_point(aes(y = coxAge.partial), size = 0.5) +
  geom_line(aes(y = loessvalues_coxAge)) +
  theme_bw()

#### 2d ####
KMedu <- survfit(Surv(MONTH_CONT, DEATH) ~ 1 + EDUCATION, data = data)
summary(KMedu)

# survivor function
plot(KMedu, xlab="Time", ylab="S(t)", conf.int=F, mark.time=F, col=c(2,3,4,5,6), lty=1, lwd=1, ylim = c(0.85,1))
legend(x = "bottomleft",
       legend=c("EDUCATION = 1", "EDUCATION = 2", "EDUCATION = 3", "EDUCATION = 4", "EDUCATION = 5"),
       col=c(2,3,4,5,6),
       lty=c(1,1))

# cumulative hazard function
plot(KMedu, fun="cumhaz", xlab="Time", ylab="H(t)=-ln[S(t)]", ylim=c(0,0.12), col=c(2,3,4,5,6), conf.int=F)
legend(x = "topleft",
       legend=c("EDUCATION = 1", "EDUCATION = 2", "EDUCATION = 3", "EDUCATION = 4", "EDUCATION = 5"),
       col=c(2,3,4,5,6),
       lty=c(1,1))

# log cumulative hazard function
plot(KMedu, fun=logcumhaz, , xlab="Time", ylab="Log[H(t)]", ylim=c(-8,-2), conf.int=F, mark.time=F, col=c(2,3,4,5,6))
legend(x = "bottomright",
       legend=c("EDUCATION = 1", "EDUCATION = 2", "EDUCATION = 3", "EDUCATION = 4", "EDUCATION = 5"),
       col=c(2,3,4,5,6),
       lty=c(1,1))

#### 2e ####
coxEdu <- coxph(Surv(data$MONTH_CONT, data$DEATH) ~ 1 + EDUCATION, data = data)
summary(coxEdu)

test_cox1coxEdu <- anova(cox1,coxEdu, test = "LRT")
print(test_cox1coxEdu)

coxEdu.partial <- resid(coxEdu, "schoenfeld")
summary(coxEdu.partial)

data.partial.Edu <- subset(data, CENSOR!= 1)
data.partial.Edu <- data.frame(data.partial.Edu, coxEdu.partial)
data.partial.Edu$RANK <- rank(data.partial.Edu$MONTH_CONT, na.last = TRUE, ties.method = "first")

loessmodel_coxEdu = loess(coxEdu.partial ~ RANK, data = data.partial.Edu)
summary(loessmodel_coxEdu)
data.partial.Edu$loessvalues_coxEdu = predict(loessmodel_coxEdu)

ggplot(data.partial.Edu, aes(x = RANK)) +
  geom_point(aes(y = coxEdu.partial), size = 0.5) +
  geom_line(aes(y = loessvalues_coxEdu)) +
  theme_bw()


# cumulative hazard function
plot(survfit(coxEdu, newdata = data.frame(EDUCATION = 1)),
     conf.int = FALSE,
     col="red", fun="cumhaz", ylim=c(0,0.12))
lines(survfit(coxEdu, newdata = data.frame(EDUCATION = 2)),
      conf.int = FALSE,
      col="orange", fun="cumhaz")
lines(survfit(coxEdu, newdata = data.frame(EDUCATION = 3)),
      conf.int = FALSE,
      col="yellow", fun="cumhaz")
lines(survfit(coxEdu, newdata = data.frame(EDUCATION = 4)),
      conf.int = FALSE,
      col="lightgreen", fun="cumhaz")
lines(survfit(coxEdu, newdata = data.frame(EDUCATION = 5)),
      conf.int = FALSE,
      col="blue", fun="cumhaz")
legend(x = "topleft",
       legend=c("EDUCATION = 1", "EDUCATION = 2", "EDUCATION = 3", "EDUCATION = 4", "EDUCATION = 5"),
       col=c("red", "orange", "yellow", "lightgreen", "blue"), lty=c(1,1))

# log cumulative hazard function
plot(survfit(coxEdu, newdata = data.frame(EDUCATION = 1)),
     conf.int = FALSE,
     col="red", fun=logcumhaz, ylim=c(-10,-2))
lines(survfit(coxEdu, newdata = data.frame(EDUCATION = 2)),
      conf.int = FALSE,
      col="orange", fun=logcumhaz)
lines(survfit(coxEdu, newdata = data.frame(EDUCATION = 3)),
      conf.int = FALSE,
      col="yellow", fun=logcumhaz)
lines(survfit(coxEdu, newdata = data.frame(EDUCATION = 4)),
      conf.int = FALSE,
      col="lightgreen", fun=logcumhaz)
lines(survfit(coxEdu, newdata = data.frame(EDUCATION = 5)),
      conf.int = FALSE,
      col="blue", fun=logcumhaz)
legend(x = "bottomright",
       legend=c("EDUCATION = 1", "EDUCATION = 2", "EDUCATION = 3", "EDUCATION = 4", "EDUCATION = 5"),
       col=c("red", "orange", "yellow", "lightgreen", "blue"), lty=c(1,1))

#### 2f ####
coxAgeEdu <- coxph(Surv(data$MONTH_CONT, data$DEATH) ~ 1 + AGE + EDUCATION, data = data)
summary(coxAgeEdu)

# effect of AGE
test_coxEducoxAgeEdu <- anova(coxEdu, coxAgeEdu, test = "LRT")
print(test_coxEducoxAgeEdu)

# effect of EDUCATION
test_coxAgecoxAgeEdu <- anova(coxAge, coxAgeEdu, test = "LRT")
print(test_coxAgecoxAgeEdu)
