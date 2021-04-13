
# Necromass decomposition - Jana Englmeier
#
# This document reproduces tables and figures.
#

wd <- "J:/PhD_Jana/Statistik/decomposition_Jana/final_data_Jana/final_data"
setwd(wd)

# Packages
library("mgcv")
library("survival")
library("tramME")
library("coxme")
library("parallel")
library("sjPlot")
library("stargazer")

# Load data

## Carrion data - Overview
x2 <- read.csv("carrion.csv", sep=";", dec=",", fileEncoding="latin1")

# PlotID:                         unique identifier for sample site
# Landscape_type_2:               landscape category (a: near-natural, b: agricultural, c: urban)
# Plot_type_3:                    habitat category (a: forest, b: grassland, c: arable field, d: settlement)
# subject:                        rat: carrion with insect access; ratAA: carrion without insect access
# Insects:                        Insect access to carrion yes/ no
# start_day1:                     earliest day of full decomposition (when rat was last seen not fully decomposed)
# end_day1:                       latest day of full decomposition (when rat was decomposed)
# comments:                       comments on missing visits/ photos
# mummified_last_visit:           1: yes; 0: no
# Übergang.6b....6:               transition from mummified stage (6b) to full decomposition (6) happened; 1: yes; 0: no
# end.stage.6.oder.6b.erreicht:   carrion reached one of the final decomposition stages 1: yes; 0: no
# average_Temp_DL_exp:            average temperature measured by data logger on study site during experiment (May until August)
# average_Humidity_DL_exp:        average humidity measured by data logger on study sites during experiment (May until August)
# elevation:                      elevation of study site
# notes:                          comments on climate data
# average_Temp_DL:                average temperature measured by data logger on study sites during entire field season (March until September)
# av_prec_exp:                    modelled average precipitation on study sites during experiment (data from Deutscher Wetterdienst)
# av_temp_exp:                    modelled average temperature on study sites during experiment (data from Deutscher Wetterdienst)  
# Annual_P_LT:                    modelled mean-annual precipitation over past 30 years on study sites (long-term data) (data from Deuscher Wetterdienst)
# Annual_T2AVG_LT:                modelled mean-annual temperature over past 30 years on study sites (long-term data) (data from Deuscher Wetterdienst)





## Dung data - Overview
wisent <- read.csv("dung.csv", sep=";", dec=",")

# PlotID:                       unique identifier for sample site
# Landscape_type_2:             landscape category (a: near-natural, b: agricultural, c: urban)
# Plot_type_3:                  habitat category (a: forest, b: grassland, c: arable field, d: settlement)
# Date.start:                   start day of experiment
# Date.end:                     end day of experiment
# Duration:                     duration of experiment (dung exposure)
# Insects:                      insect access yes: 1; 0: no
# Ausgangsgewicht:              initital wet weight
# Ausgangsgewicht_trocken:      initial dry weight
# dung_wet:                     final wet weight
# dung_dry:                     final dry weight
# weight_loss:                  dry weight loss from start until end of experiment 
# weight_loss_percent:          dry weight loss in percent from start until end of experiment
# elevation:                    elevation of study site
# average_Temp_DL_exp:          average temperature measured by data logger on study site during experiment (May until August)
# average_Humidity_DL_exp:      average humidity measured by data logger on study sites during experiment (May until August)
# av_prec_exp:                  modelled average precipitation on study sites during experiment (data from Deutscher Wetterdienst)
# av_temp_exp:                  modelled average temperature on study sites during experiment (data from Deutscher Wetterdienst)  
# Annual_P_LT:                  modelled mean-annual precipitation over past 30 years on study sites (long-term data) (data from Deuscher Wetterdienst)
# Annual_T2AVG_LT:              modelled mean-annual temperature over past 30 years on study sites (long-term data) (data from Deuscher Wetterdienst)


## rescale predictor variables: elevation per 100 m
x2$elevation100 <- x2$elevation / 100
wisent$elevation100 <- wisent$elevation / 100



## Check for correlation between environmental variables using wisent-data set including all study sites

# elevation and local temperature (data logger) during field season
shapiro.test(wisent$elevation) # p < 0.05
cor.test(wisent$elevation, wisent$average_Temp_DL, method="spearman") # rho = -0.43, p < 0.05
  # --> elevation and average temperature measured by data logger correlate only moderately, hence both are included in models

# elevation and long-term precipitation
cor.test(wisent$elevation, wisent$Annual_P_LT, method="spearman") # rho = 0.74, p < 0.05

# elevation and long-term temperature
cor.test(wisent$elevation, wisent$Annual_T2AVG_LT, method="spearman") # rho = -0.84, p < 0.05

  # --> elevation highly correlates with long-term temperature and precipitation and is hence chosen as a surrogate for macroclimate



########################################
### carrion data preparation for modelling

x2$end_day1 <- as.numeric(x2$end_day1)


x2[["Insects"]] <- factor(x2[["Insects"]]) # Insects 1: insect access; Insects 0: no insect access
x2[["Landscape"]] <- factor(x2[["Landscape_type_2"]], 
                            levels = letters[1:3], 
                            labels = c("seminatural", "agriculture", "urban"))

x2[["Habitat"]] <- factor(x2[["Plot_type_3"]], 
                          levels = letters[1:4], 
                          labels = c("forest", "meadow", "arable field", "settlement"))

x2[["PlotID"]] <- factor(x2[["PlotID"]])

x2[["Status"]] <- x2[["end.stage.6.oder.6b.erreicht"]] == 1 # 1: reached final decomposition stage 6 or 6b; 0: did not reach final decomp. stage (Rechtszensierung)

# replace end_day1 = 1000 by start_day1 for rats that did not fully decompose; Status = 0
x2$end_day_gam = x2$end_day1
x2$end_day_gam[x2$Status == FALSE] = x2$start_day1[x2$Status==FALSE]



## generalized additive models for carrion decomposition (family = cox.ph)

# gam 1) including elevation as a surrogate for macroclimate 
# gam 2) including long-term temperature and precipitation as a surrogate for macroclimate


# 1) gam including elevation and local temperature measured by data logger ## 
gam_rat <- gam(end_day_gam~ Insects + Habitat + Landscape  + average_Temp_DL_exp +s(elevation100) + s(PlotID, bs="re"),
               family=cox.ph(), data=x2, weights=Status)
summary(gam_rat)

### just a check
tmp_gam <- gam(end_day_gam~ Insects + Habitat + Landscape  + average_Temp_DL_exp + elevation100 + s(PlotID, bs="re"),
               family=cox.ph(), data=x2, weights=Status)
tmp_ME <- CoxphME(Surv(end_day_gam, Status) ~ Insects + Habitat + Landscape  + average_Temp_DL_exp + elevation100
                  + (1 | PlotID), data = x2)
tmp_me <- coxme(Surv(end_day_gam, Status) ~ Insects + Habitat + Landscape  + average_Temp_DL_exp + elevation100
                + (1 | PlotID), data = x2)

### s(elevation) not needed -> linear Cox model
logLik(gam_rat)
logLik(tmp_gam)

logLik(tmp_ME)
logLik(tmp_me)

### <TH> "almost" equivalent. Note that fitting these models is really hard
summary(gam_rat)
summary(tmp_gam)
summary(tmp_ME)
summary(tmp_me)


# partial plot gam_rat
plot(gam_rat, shade=TRUE, cex.lab=1.2, ylab="Partial effect (elevation)", xlab="Elevation [100m]", main="Carrion (decay rate)", cex.main=1.3, select=1, cex.lab=1.2)
text(2, 3, "B)", cex=1.5)


# 2) gam including long-term temperature and precipitation as surrogates for macroclimate ##
gam_rat2 <- gam(end_day_gam~ Insects + Habitat + Landscape  + Annual_T2AVG_LT + Annual_P_LT  + s(PlotID, bs="re"),
                family=cox.ph(), data=x2, weights=Status)
summary(gam_rat2)


# compare AIC of gam with elevation and gam with long-term weather data to decide for surrogate for macroclimate
AIC(gam_rat) #2208.31 mit elevation
AIC(gam_rat2) # 2279.51 mit long-term temperature + precipitation




## mixed-effect parametric cox regression

start <- x2$start_day1 # start_day1: earliest potential date of full carrion decomposition
stop <- x2$end_day1 # end_day1: latest day of full carrion decomposition
stop[stop == 1000] <- Inf
x2$time <- with(x2, Surv(time = start, time2 = stop, type = "interval2"))

rats2 <- x2[, c("PlotID", "time", "Insects", 
                "Landscape", "Habitat", "elevation100", "average_Temp_DL_exp", "Annual_T2AVG_LT", "Annual_P_LT")]

#dump("rats2")

# parametric cox regression with elevation as surrogate for macroclimate
cm_elevation <- CoxphME(time ~ Insects + Landscape + Habitat  +  average_Temp_DL_exp + elevation100 +
                 (1 | PlotID), data = rats2)
summary(cm_elevation)
cf_elevation <- exp(coef(cm_elevation))
cbind(cf_elevation, exp(confint(cm_elevation))[names(cf_elevation),])

# parametric cox regression with long-term temperature and precipitation as surrogate for macroclimate
cm_LT_data <- CoxphME(time ~ Insects + Habitat + Landscape  + Annual_T2AVG_LT + Annual_P_LT +  
                 (1 | PlotID), data = rats2)
summary(cm_LT_data)
cf_LT <- exp(coef(cm_LT_data))
cbind(cf_LT, exp(confint(cm_LT_data))[names(cf_LT),])


## plot survival curve
### <TH> mark right-censoring times
### this is the Turnbull estimator 
fit <- survfit(time ~ Insects, data = rats2)
plot(fit, xlab = "Time after carrion exposure [days]", ylab = "Probability of full decomposition", 
     fun="event", lty=1:2, main="", cex.main=1.3, cex.lab=1.2, mark = "/")
legend("bottomright", lty=2:1, 
       legend=c("allowed", "excluded"), bty="n", title="   Insects", title.adj=0)

### marginal distribution functions
## A function to evaluate the joint cdf of the response and the random effects:
## Takes a vector of random effect and covariates values, evaluates the conditional
## distribution at these values and multiplies it with the pdf of the random effects
jointCDF <- function(re, nd, mod) {
  nd <- nd[rep(1, length(re)), ]
  nd$PlotID <- seq(nrow(nd)) ## to take vector-valued REs
  pr <- predict(mod, newdata = nd, ranef = re, type = "distribution") *
    dnorm(re, 0, sd = sqrt(varcov(mod)[[1]][1, 1]))
  c(pr)
}
## Marginalize the joint cdf by integrating out the random effects
## using adaptive quadratures
marginalCDF <- function(nd, mod) {
  nd$cdf <- integrate(jointCDF, lower = -Inf, upper = Inf, nd = nd, mod = mod)$value
  nd
}
## Set up the grid on which we evaluate the marginal distribution
nd <- x3[idx <- rep(1:2, each = 16),]
nd$ID <- factor(idx) 
nd$time <- 1:16 * 5 #### <TH> ~ between 5 and 80
## Calls marginalCDF on each row of nd
## (done in parallel to speed up computations)
mp <- parallel::mclapply(split(nd, seq(nrow(nd))),
                         marginalCDF, mod = cm4, mc.cores = 4)
mp <- do.call("rbind", mp)

xyplot(cdf ~ time | ID, data = mp, type = "l")






###################################################

#### Dung decomposition #####

# data preparation

### <TH> das sieht ziemlich konstant aus, dann macht es nicht viel
### Sinn, das als log-offset zu nehmen
wisent$Ausgangsgewicht_trocken <- as.numeric(wisent$Ausgangsgewicht_trocken)

wisent[["Insects"]] <- factor(wisent[["Insects"]]) # 1: insect access; 0: no insects
wisent[["Landscape"]] <- factor(wisent[["Landscape_type_2"]], 
                                levels = letters[1:3], 
                                labels = c("seminatural", "agriculture", "urban"))

wisent[["Habitat"]] <- factor(wisent[["Plot_type_3"]], 
                              levels = letters[1:4], 
                              labels = c("forest", "meadow", "arable field", "settlement"))

wisent[["PlotID"]] <- factor(wisent[["PlotID"]])


# exclude dung pats with increasing weight due to invertebrates inside the dung pat
wisent<-subset(wisent,weight_loss<100) 


## gams family = gaussian 
# 1) including elevation as a surrogate for macroclimate 
# 2) including long-term temp. and precipitation as a surrogate for macroclimate

# 1) gam including elevation and local temperature (data logger)
gam_dung <- gam(dung_dry ~ offset(log(Ausgangsgewicht_trocken)+log(I(Duration/29))) + Insects + 
                  Habitat +  Landscape  + average_Temp_DL_exp +
                  s(elevation100)  + s(PlotID,bs="re")
                ,family = gaussian(link="log"), method = "REML", data=wisent)
summary(gam_dung)


# gam without offset(Ausgangsgewicht)
gam_dung1 <- gam(dung_dry ~ offset(log(I(Duration/29))) + Insects + 
                   Habitat +  Landscape  + average_Temp_DL_exp +
                   s(elevation100)  + s(PlotID,bs="re")
                 ,family = gaussian(link="log"), method = "REML", data=wisent)
summary(gam_dung1)

# partial plot 
plot(gam_dung, shade=TRUE, cex.lab=1.2, ylab="Partial effect (elevation)", xlab="Elevation [100m]", main="Dung (final dry weight)", cex.main=1.3, select=1, cex.lab=1.2) 
text(2, 0.51, "A)", cex=1.5)


# 2) gam including long-term temperature and precipitation
# merge data frames wisent and climate to include long-term precipitation and temperature
gam_dung2 <- gam(dung_dry ~ offset(log(Ausgangsgewicht_trocken)+log(I(Duration/29))) + Insects + 
                   Habitat +  Landscape  + Annual_T2AVG_LT + Annual_P_LT +
                   s(PlotID,bs="re")
                 ,family = gaussian(link="log"), method = "REML", data=wisent)
summary(gam_dung2)



AIC(gam_dung) #  2714.377 mit elevation
AIC(gam_dung2) # 2775.368 mit long-term temp + precip


# produce output tables for manuscript
## gams including elevation as macroclimate surrogate
tab_model(gam_dung, gam_rat, transform = NULL, show.ci = FALSE, show.se=TRUE)

## gams including long-term temperature and precipitation data as macroclimate surrogate (for Appendix)
tab_model(gam_dung2, gam_rat2, transform = NULL, show.ci = FALSE, show.se=TRUE)
