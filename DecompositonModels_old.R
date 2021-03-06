
#wd <- "J:/PhD_Jana/Statistik/decomposition_Jana/final_data_Jana/final_data"
#setwd(wd)


###### Decomposition Models #######

library("mgcv")
library("survival")
library("tramME")
library("coxme")
library("parallel")


#### data sets ####

# Rats (carrion measurements, environmental parameters incl. elevation and data logger temperature on plot)
x2 <- read.csv("Torsten_modifiedJE.csv", sep=";", dec=",", fileEncoding="latin1") # Carrion dataset

# Dung (including measurements, environmental parameters incl. elevation and data logger temperature on plot)
wisent<-read.csv2("wisent2_JE_new.csv", header=TRUE, dec=".")


# klima- dataset:
# PlotID: Study site ID
# average_Temp_DL: local temperature measured by datalogger on study sites during experiment (also included in data sets x2 and wisent)
# elevation: elevation of study sites (also included in data sets x2 and wisent)
# av_prec_exp: modelled average precipitation during experiment (data from Deutscher Wetterdienst)
# av_temp_exp: modelled average temperature during experiment (data from Deutscher Wetterdienst)
# Annual_P_LT: modelled mean-annual precipitation over past 30 years on study sites (long-term data) (data from Deuscher Wetterdienst)
# Annual_T2AVG_LT: modelled mean-annual temperature over past 30 years on study sites (long-term data) (data from Deuscher Wetterdienst))
climate <- read.csv("klima.csv", row.names = 1L)

#### Correlation of environmental variables ####

# elevation and local temperature (data logger) during experiment
hist(climate$elevation)
hist(climate$average_Temp_DL)
shapiro.test(climate$elevation) # p < 0.05
cor.test(climate$elevation, climate$average_Temp_DL, method="spearman") # rho = -0.43, p < 0.05

# elevation and long-term precipitation
cor.test(climate$elevation, climate$Annual_P_LT, method="spearman") # rho = 0.74, p < 0.05

# elevation and long-term temperautre
cor.test(climate$elevation, climate$Annual_T2AVG_LT, method="spearman") # rho = -0.84, p < 0.05

# --> elevation highly correlates with long-term temperature and precipitation and is hence chosen as a surrogate for macroclimate




###############################
#### Carrion decomposition ####

### data preparation

x2$end_day1 <- as.numeric(x2$end_day1)
### <TH> das sieht als Zensierungszeiten komisch aus. Ist 1000 IMMER der 
### letzte Tag, an dem die Ratte beobachtet wurde OHNE Zersetzung?
#<JE> Wenn eines der finalen Zersetzungsstadien 6 oder 6b nicht erreicht wurde, sondern nur Stadium 4 oder 5, 
#haben wir der Ratte einen extrem hohen Wert gegeben, daf??r dass sie irgendwann in der Zukunft noch zersetzt wird
x2$end_day1[x2[["end.stage.6.oder.6b.erreicht"]] == 0] # 1: end stage (finales Zersetzungsstadium) erreicht; 0: nicht erreicht

### <TH> was bedeutet 0 und was 1? labels richtig setzen
x2[["Insects"]] <- factor(x2[["Insects"]]) # <JE> Insects 0: kein Insektenzugang; Insects 1: Insekten hatten Zugang zur Ratte
x2[["Landscape"]] <- factor(x2[["Landscape_type_2"]], 
                            levels = letters[1:3], 
                            labels = c("seminatural", "agriculture", "urban"))

x2[["Habitat"]] <- factor(x2[["Plot_type_3"]], 
                          levels = letters[1:4], 
                          labels = c("forest", "meadow", "arable field", "settlement"))

x2[["PlotID"]] <- factor(x2[["PlotID"]])

### <TH> das sollte ein logical sein, mit TRUE (= 1?) fuer events und FALSE
### fuer Rechtszensierung (= 0?) <JE> hier bin ich etwas verwirrt
x2[["Status"]] <- x2[["end.stage.6.oder.6b.erreicht"]] == 1

# downscaling of elevation to adjust temperature and elevation scales
x2$elevation100 <- x2$elevation / 100


#### generalized additive models for carrion decomposition 

## gams family = cox.ph() 
# 1) including elevation as a surrogate for macroclimate 
# 2) including long-term temperature and precipitation as a surrogate for macroclimate

# 1) gam including elevation and local temperature measured by data logger ## 
x2$end_day_gam = x2$end_day1
x2$end_day_gam[x2$Status == FALSE] = x2$start_day1[x2$Status==FALSE]

gam_rat <- gam(end_day_gam~ Insects + Habitat + Landscape  + average_Temp_DL +s(elevation100) + s(PlotID, bs="re"),
               family=cox.ph(), data=x2, weights=Status)
summary(gam_rat)

### just a check
tmp_gam <- gam(end_day_gam~ Insects + Habitat + Landscape  + average_Temp_DL + elevation100 + s(PlotID, bs="re"),
               family=cox.ph(), data=x2, weights=Status)
tmp_ME <- CoxphME(Surv(end_day_gam, Status) ~ Insects + Habitat + Landscape  + average_Temp_DL + elevation100
                  + (1 | PlotID), data = x2)
tmp_me <- coxme(Surv(end_day_gam, Status) ~ Insects + Habitat + Landscape  + average_Temp_DL + elevation100
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



# partial plot
plot(gam_rat, shade=TRUE, cex.lab=1.2, ylab="Partial effect (elevation)", xlab="Elevation [100m]", main="Carrion (decay rate)", cex.main=1.3, select=1, cex.lab=1.2)
text(2, 3, "B)", cex=1.5)


# 2) gam including long-term temperature and precipitation ##
# merge data frames x2 and climate to include long-term temperature and precipitation in model
x3 = merge(x2, climate, by="PlotID")

gam_rat2 <- gam(end_day_gam~ Insects + Habitat + Landscape  + Annual_T2AVG_LT + Annual_P_LT  + s(PlotID, bs="re"),
                family=cox.ph(), data=x3, weights=Status)
summary(gam_rat2)


AIC(gam_rat) #2225.31 mit elevation
AIC(gam_rat2) # 2296.973 mit long-term temperature + precipitation


# produce output tables

#library("sjPlot")
#library("stargazer")
#tab_model(gam_rat, gam_rat2)



## mixed-effect parametric cox regression ###
### <TH> was ist start_day1?

#<JE> start_day1 ist der fr??heste Tag an dem die Ratte zersetzt gewesen sein k??nnte; 
# Bsp: die Ratte war an Tag 27 noch nicht zersetzt, Tag 38 kommen wir wieder an den Plot und die Ratte war zersetzt,
# d.h. die Ratte muss irgendwann zwischen Tag 28 (start_day1) und 38 (end_day1) komplett zersetzt worden sein
# war die Ratte an unserem letzten Besuch immer noch nicht zersetzt, haben wir ihr den Wert 1000 als end_day1 zugeordnet
start <- x2$start_day1 
stop <- x2$end_day1
stop[stop == 1000] <- Inf
x3$time <- with(x3, Surv(time = start, time2 = stop, type = "interval2"))

rats2 <- x3[, c("PlotID", "time", "Insects", 
                "Landscape", "Habitat", "elevation100", "average_Temp_DL.x", "Annual_T2AVG_LT", "Annual_P_LT")]

#dump("rats2")

cm3 <- CoxphME(time ~ Insects + Landscape + Habitat  +  average_Temp_DL.x + elevation100 +
                 (1 | PlotID), data = rats2)
summary(cm3)
cf3 <- exp(coef(cm3))
cbind(cf3, exp(confint(cm3))[names(cf3),])

cm4 <- CoxphME(time ~ Insects + Habitat + Landscape  + Annual_T2AVG_LT + Annual_P_LT +  
                 (1 | PlotID), data = rats2)
summary(cm4)
cf4 <- exp(coef(cm4))
cbind(cf4, exp(confint(cm4))[names(cf4),])




# plot survival curve
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
sd(wisent$Ausgangsgewicht_trocken)


### <TH> labels?
wisent[["Insects"]] <- factor(wisent[["Insects"]]) # 1: insect access; 0: no insects
wisent[["Landscape"]] <- factor(wisent[["Landscape_type_2"]], 
                                levels = letters[1:3], 
                                labels = c("seminatural", "agriculture", "urban"))

wisent[["Habitat"]] <- factor(wisent[["Plot_type_3"]], 
                              levels = letters[1:4], 
                              labels = c("forest", "meadow", "arable field", "settlement"))

wisent[["PlotID"]] <- factor(wisent[["PLOTID"]])

### <TH> was macht das subset?
wisent<-subset(wisent,weight_loss<100) # nur Dung mit Gewichtsabnahme wird in Datensatz aufgenommen, eine Dungprobe hatte an Gewicht zugenommen (bspw. K??fer und W??rmer, die im Dung waren und mit gewogen wurden)
wisent$elevation100<-wisent$elevation/100



## gams family = gaussian 
# 1) including elevation as a surrogate for macroclimate 
# 2) including long-term temp. and precipitation as a surrogate for macroclimate

# 1) gam including elevation and local temperature (data logger)
### <TH> Duration is Beobachtungsdauer? <JE> Genau, die Dauer des Experiments. Beim Dung gibt es nur zwei Zeitpunkte, Tag an dem wir ihn ausgebracht haben und Tag an dem wir wieder eingesammelt haben
gam_dung <- gam(dung_dry ~ offset(log(Ausgangsgewicht_trocken)+log(I(Duration/29))) + Insects + 
                  Habitat +  Landscape  + average_Temp_DL +
                  s(elevation100)  + s(PlotID,bs="re")
                ,family = gaussian(link="log"), method = "REML", data=wisent)
summary(gam_dung)

# partial plot
plot(gam_dung, shade=TRUE, cex.lab=1.2, ylab="Partial effect (elevation)", xlab="Elevation [100m]", main="Dung (final dry weight)", cex.main=1.3, select=1, cex.lab=1.2) 
text(2, 0.51, "A)", cex=1.5)

# 2) gam including long-term temperature and precipitation
# merge data frames wisent and climate to include long-term precipitation and temperature
wisent2 = merge(wisent, climate, by="PlotID")
gam_dung2 <- gam(dung_dry ~ offset(log(Ausgangsgewicht_trocken)+log(I(Duration/29))) + Insects + 
                   Habitat +  Landscape  + Annual_T2AVG_LT + Annual_P_LT +
                   s(PlotID,bs="re")
                 ,family = gaussian(link="log"), method = "REML", data=wisent2)
summary(gam_dung2)



AIC(gam_dung) #  2714.377 mit elevation
AIC(gam_dung2) # 2775.368 mit long-term temp + precip

tab_model(gam_dung, gam_dung2)

