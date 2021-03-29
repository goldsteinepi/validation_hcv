#################
# Validation of Hep C testing and diagnosis
# Citation: Goldstein ND, Kahal D, Testa K, Burstyn I. Data quality in electronic health record research: an approach for validation and quantitative bias analysis for the data scientist. Manuscript in preparation.
# 4/3/19 -- Neal Goldstein
#################


### JOIN VALIDATION DATA ###

#read data, risk cohort, universal cohort, validation of risk cohort
load(file="Retrospective.Rdata")
load(file="Prospective.Rdata")
validation = read.csv("validation.csv", as.is=T, stringsAsFactors=F, na.strings="")

#merge
retro = merge(retro, validation, by.x="mrn", by.y="MRN", all.x=T, all.y=F)
retro$Determination[is.na(retro$Determination)] = "Negative"

#check whether anyone in the universal cohort had a risk factor for HCV screen; we infer this from a hep c test being ordered previously
#also check whether anyone had previously diagnosed hcv, as again would place individuals in the risk factor group
pro$riskfactor = NA
pro$hepc_old = NA
for (i in 1:nrow(pro))
{
  pro$riskfactor[i] = ifelse(length(retro$hepc_test_ordered[retro$mrn==pro$mrn[i]])>0, retro$hepc_test_ordered[retro$mrn==pro$mrn[i]], 0)
  pro$hepc_old[i] = ifelse(length(retro$hepc[retro$mrn==pro$mrn[i]])>0, retro$hepc[retro$mrn==pro$mrn[i]], 0)
}
rm(i)

#save analytic data
rm(demo,hepc,hepc_orders,hepc_results,hiv,hiv_orders,hiv_results,validation)
retro$mrn = NULL
retro$name = NULL
retro$dob = NULL
pro$mrn = NULL
pro$name = NULL
pro$dob = NULL
save.image(file="validation.RData")


### FUNCTIONS ###

library(psych) #describe, describeBy
library(gmodels) #CrossTable
library(epiR) #sensitivity, specificity, predictive values
library(blm) #expit
library(sandwich) #robust standard errors
library(reldist) #weighted quantiles


### LOAD DATA ###

load("validation.RData")


### DETERMINE INFECTION STATUS ###

#ignore HCV+ lab results for 2 individuals without corresponding HCV diagnostic code: we do not know the reason for the test
retro$Determination[which(retro$hepc==0 & (retro$Determination=="PCR Positive" | retro$Determination=="Cured"))] = "Negative"

#remove W=1 individuals from universal screening era (as they would have been screened under risk factor paradigm); that is, these individuals are W=0
pro = pro[pro$riskfactor==0 & pro$hepc_old==0, ]

#create a diagnostic variable for the universal screening group
pro$Determination = ifelse((!is.na(pro$hepc_test_pcr) & (pro$hepc_test_pcr==1)), 1, 0)

#recode race: 0=white, 1=non-white
retro$race = ifelse(retro$race==0, 0, 1)
pro$race = ifelse(pro$race==0, 0, 1)


### DESCRIPTIVES OF COHORT ###

describe(retro$last_visit_age); quantile(retro$last_visit_age, c(0.25, 0.75))
CrossTable(retro$sex)
CrossTable(retro$race)
CrossTable(retro$ethnicity)
CrossTable(retro$insurance==1)
CrossTable(retro$hepc_test_ordered)
CrossTable(retro$hepc_test_resulted)
CrossTable(retro$hepc)
CrossTable(retro$Determination=="PCR Positive" | retro$Determination=="Cured")

describe(pro$last_visit_age); quantile(pro$last_visit_age, c(0.25, 0.75))
CrossTable(pro$sex)
CrossTable(pro$race)
CrossTable(pro$ethnicity)
CrossTable(pro$insurance==1)
CrossTable(pro$hepc_test_ordered)
CrossTable(pro$hepc_test_resulted)
CrossTable(pro$Determination)

describe(pro$last_visit_age[pro$hepc_test_resulted==1]); quantile(pro$last_visit_age[pro$hepc_test_resulted==1], c(0.25, 0.75))
CrossTable(pro$sex[pro$hepc_test_resulted==1])
CrossTable(pro$race[pro$hepc_test_resulted==1])
CrossTable(pro$ethnicity[pro$hepc_test_resulted==1])
CrossTable(pro$insurance[pro$hepc_test_resulted==1]==1)
#CrossTable(pro$hepc_test_ordered[pro$hepc_test_resulted==1])
CrossTable(pro$Determination[pro$hepc_test_resulted==1])

describe(pro$last_visit_age[pro$hepc_test_resulted==0]); quantile(pro$last_visit_age[pro$hepc_test_resulted==0], c(0.25, 0.75))
CrossTable(pro$sex[pro$hepc_test_resulted==0])
CrossTable(pro$race[pro$hepc_test_resulted==0])
CrossTable(pro$ethnicity[pro$hepc_test_resulted==0])
CrossTable(pro$insurance[pro$hepc_test_resulted==0]==1)
CrossTable(pro$hepc_test_ordered[pro$hepc_test_resulted==0])
#CrossTable(pro$Determination[pro$hepc_test_resulted==0])


### VALIDATION SUB-COHORT ###

#start with risk factor: W=1
validation = retro[retro$hepc==1, c("redcap_id", "race", "Determination")]
validation$cohort = "risk factor"
validation$W = 1
validation$X = ifelse((validation$Determination=="PCR Positive" | validation$Determination=="Cured"), 1, 0)
validation$Determination = NULL

#combine with universal: W=0
validation = rbind(validation, data.frame("redcap_id"=pro$redcap_id[pro$hepc_test_resulted==1], "race"=pro$race[pro$hepc_test_resulted==1], "cohort"="universal", "W"=0, "X"=ifelse(is.na(pro$hepc_test_pcr[pro$hepc_test_resulted==1]), 0, 1), stringsAsFactors=F))

#estimates of accuracy from model
model1 = glm(X ~ W, family=binomial(link="logit"), data=validation)
model2 = glm(W ~ X, family=binomial(link="logit"), data=validation)
summary(model1)
summary(model2)

#accuracy parameters
PPV = expit(model1$coefficients[1] + model1$coefficients[2])
FOR = (expit(model1$coefficients[1]))
NPV = 1 - FOR
FDR = 1 - PPV

SN = expit(model2$coefficients[1] + model2$coefficients[2])
FPR = expit(model2$coefficients[1])
SP = 1 - FPR
FNR = 1 - SN

#select precision estimates obtained through bootstrapping
set.seed(777)
PPV_CI = NA
for (i in 1:1000){
  validation_sample = validation[sample(1:nrow(validation), nrow(validation), replace=T), ]
  boot_model = glm(X ~ W, family=binomial(link="logit"), data=validation_sample)
  PPV_CI = c(PPV_CI, expit(boot_model$coefficients[1] + boot_model$coefficients[2]))
}
rm(i, validation_sample, boot_model)
quantile(PPV_CI, probs=c(0.025,0.5,0.975), na.rm=T)

NPV_CI = NA
for (i in 1:1000){
  validation_sample = validation[sample(1:nrow(validation), nrow(validation), replace=T), ]
  boot_model = glm(X ~ W, family=binomial(link="logit"), data=validation_sample)
  NPV_CI = c(NPV_CI, 1-(expit(boot_model$coefficients[1])))
}
rm(i, validation_sample, boot_model)
quantile(NPV_CI, probs=c(0.025,0.5,0.975), na.rm=T)

SN_CI = NA
for (i in 1:1000){
  validation_sample = validation[sample(1:nrow(validation), nrow(validation), replace=T), ]
  boot_model = glm(W ~ X, family=binomial(link="logit"), data=validation_sample)
  SN_CI = c(SN_CI, expit(boot_model$coefficients[1] + boot_model$coefficients[2]))
}
rm(i, validation_sample, boot_model)
quantile(SN_CI, probs=c(0.025,0.5,0.975), na.rm=T)

SP_CI = NA
for (i in 1:1000){
  validation_sample = validation[sample(1:nrow(validation), nrow(validation), replace=T), ]
  boot_model = glm(W ~ X, family=binomial(link="logit"), data=validation_sample)
  SP_CI = c(SP_CI, 1-expit(boot_model$coefficients[1]))
}
rm(i, validation_sample, boot_model)
quantile(SP_CI, probs=c(0.025,0.5,0.975), na.rm=T)

#corroborate agains a 2x2 table
contigency = CrossTable(validation$W, validation$X, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)$t
contigency

#re-organize table and calculate accuracy
contigency = matrix(c(contigency[4],contigency[3],contigency[2],contigency[1]),nrow=2,ncol=2)
epi.tests(contigency)

# #exact logistic regression for PPV: see https://stats.idre.ucla.edu/r/dae/exact-logistic-regression/
# library(elrm) #exact logistic regression (if difficulty installing, see https://stackoverflow.com/questions/65307023/is-there-an-alternative-to-elrm-exact-logistic-regression-in-r-the-elrm-pack)
# x = xtabs(~X + W, data=validation)
# cdat = data.frame(W = c(0,1), X = x[2, ], ntrials = colSums(x))
# model1.elrm = elrm(formula=X/ntrials ~ W, interest=~W, r=2, dataset=cdat, iter=100000)
# summary(model1.elrm)
# model1.elrm$coeffs
# model1.elrm$coeffs.ci

# #Firth's Bias-Reduced Logistic Regression
# library(logistf)
# model1.firth = logistf(X ~ W, data=validation)
# PPV = expit(model1.firth$coefficients[1] + model1.firth$coefficients[2])
# PPV_lo = expit(confint(model1.firth)[1,1] + confint(model1.firth)[2,1])
# PPV_hi = expit(confint(model1.firth)[1,2] + confint(model1.firth)[2,2])

#amplified data for covariate stratification on ethnicity
set.seed(777)
validation_100x = validation[sample(1:nrow(validation), nrow(validation)*100, replace=T), ]

#2x2 table
CrossTable(validation_100x$W[ validation_100x$race==0], validation_100x$X[ validation_100x$race==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(validation_100x$W[ validation_100x$race==1], validation_100x$X[ validation_100x$race==1], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)

#estimates of accuracy from model
model3 = glm(X ~ W*race, family=binomial(link="logit"), data=validation_100x)
model4 = glm(W ~ X*race, family=binomial(link="logit"), data=validation_100x)
summary(model3)
summary(model4)

#accuracy parameters
PPV0 = expit(model3$coefficients[1] + model3$coefficients[2])
PPV1 = expit(model3$coefficients[1] + model3$coefficients[2] + model3$coefficients[3] + model3$coefficients[4])
NPV0 = 1 - (expit(model3$coefficients[1]))
NPV1 = 1 - (expit(model3$coefficients[1] + model3$coefficients[3]))

SN0 = expit(model4$coefficients[1] + model4$coefficients[2])
SN1 = expit(model4$coefficients[1] + model4$coefficients[2] + model4$coefficients[3] + model4$coefficients[4])
SP0 = 1 - (expit(model4$coefficients[1]))
SP1 = 1 - (expit(model4$coefficients[1] + model4$coefficients[3]))

#select precision estimates obtained through bootstrapping
PPV0_CI = NA
for (i in 1:1000){
  validation_100x_sample = validation_100x[sample(1:nrow(validation_100x), nrow(validation_100x), replace=T), ]
  boot_model = glm(X ~ W*race, family=binomial(link="logit"), data=validation_100x_sample)
  PPV0_CI = c(PPV0_CI, expit(boot_model$coefficients[1] + boot_model$coefficients[2]))
}
rm(i, validation_100x_sample, boot_model)
quantile(PPV0_CI, probs=c(0.025,0.5,0.975), na.rm=T)

PPV1_CI = NA
for (i in 1:1000){
  validation_100x_sample = validation_100x[sample(1:nrow(validation_100x), nrow(validation_100x), replace=T), ]
  boot_model = glm(X ~ W*race, family=binomial(link="logit"), data=validation_100x_sample)
  PPV1_CI = c(PPV1_CI, expit(boot_model$coefficients[1] + boot_model$coefficients[2] + boot_model$coefficients[3] + boot_model$coefficients[4]))
}
rm(i, validation_100x_sample, boot_model)
quantile(PPV1_CI, probs=c(0.025,0.5,0.975), na.rm=T)

NPV0_CI = NA
for (i in 1:1000){
  validation_100x_sample = validation_100x[sample(1:nrow(validation_100x), nrow(validation_100x), replace=T), ]
  boot_model = glm(X ~ W*race, family=binomial(link="logit"), data=validation_100x_sample)
  NPV0_CI = c(NPV0_CI, 1-(expit(boot_model$coefficients[1])))
}
rm(i, validation_100x_sample, boot_model)
quantile(NPV0_CI, probs=c(0.025,0.5,0.975), na.rm=T)

NPV1_CI = NA
for (i in 1:1000){
  validation_100x_sample = validation_100x[sample(1:nrow(validation_100x), nrow(validation_100x), replace=T), ]
  boot_model = glm(X ~ W*race, family=binomial(link="logit"), data=validation_100x_sample)
  NPV1_CI = c(NPV1_CI, 1-(expit(boot_model$coefficients[1] + boot_model$coefficients[3])))
}
rm(i, validation_100x_sample, boot_model)
quantile(NPV1_CI, probs=c(0.025,0.5,0.975), na.rm=T)

SN0_CI = NA
for (i in 1:1000){
  validation_100x_sample = validation_100x[sample(1:nrow(validation_100x), nrow(validation_100x), replace=T), ]
  boot_model = glm(W ~ X*race, family=binomial(link="logit"), data=validation_100x_sample)
  SN0_CI = c(SN0_CI, expit(boot_model$coefficients[1] + boot_model$coefficients[2]))
}
rm(i, validation_100x_sample, boot_model)
quantile(SN0_CI, probs=c(0.025,0.5,0.975), na.rm=T)

SN1_CI = NA
for (i in 1:1000){
  validation_100x_sample = validation_100x[sample(1:nrow(validation_100x), nrow(validation_100x), replace=T), ]
  boot_model = glm(W ~ X*race, family=binomial(link="logit"), data=validation_100x_sample)
  SN1_CI = c(SN1_CI, expit(boot_model$coefficients[1] + boot_model$coefficients[2] + boot_model$coefficients[3] + boot_model$coefficients[4]))
}
rm(i, validation_100x_sample, boot_model)
quantile(SN1_CI, probs=c(0.025,0.5,0.975), na.rm=T)

SP0_CI = NA
for (i in 1:1000){
  validation_100x_sample = validation_100x[sample(1:nrow(validation_100x), nrow(validation_100x), replace=T), ]
  boot_model = glm(W ~ X*race, family=binomial(link="logit"), data=validation_100x_sample)
  SP0_CI = c(SP0_CI, 1-expit(boot_model$coefficients[1]))
}
rm(i, validation_100x_sample, boot_model)
quantile(SP0_CI, probs=c(0.025,0.5,0.975), na.rm=T)

SP1_CI = NA
for (i in 1:1000){
  validation_100x_sample = validation_100x[sample(1:nrow(validation_100x), nrow(validation_100x), replace=T), ]
  boot_model = glm(W ~ X*race, family=binomial(link="logit"), data=validation_100x_sample)
  SP1_CI = c(SP1_CI, 1-(expit(boot_model$coefficients[1] + boot_model$coefficients[3])))
}
rm(i, validation_100x_sample, boot_model)
quantile(SP1_CI, probs=c(0.025,0.5,0.975), na.rm=T)


### QBA ###

## simple case of prevalance

#regression parameters
beta0 = summary(model1)$coefficients[1, 1]
beta0_se = summary(model1)$coefficients[1, 2]
beta1 = summary(model1)$coefficients[2, 1]
beta1_se = summary(model1)$coefficients[2, 2]

#generate 1,000 realizations of x
set.seed(777)
nsamples = 1000
x_matrix = matrix(NA, nrow=nrow(retro), ncol=nsamples)
sim_PPVs=NA; sim_NPVs=NA; sim_prev=NA
for (i in 1:nsamples) {
  #sample parameters
  beta_dot0 = rnorm(1, beta0, beta0_se)
  beta_dot1 = rnorm(1, beta1, beta1_se)
  
  #model of true prevalence
  logit.px = beta_dot0 + beta_dot1*retro$hepc
  px = expit(logit.px)
  
  #generate x
  x_matrix[, i] = rbinom(nrow(retro), 1, px)
  
  #save key estimates
  sim_PPVs = c(sim_PPVs, expit(beta_dot0 + beta_dot1))
  sim_NPVs = c(sim_NPVs, (1 - (expit(beta_dot0))))
  sim_prev = c(sim_prev, mean(px))
}
rm(i, logit.px, px, nsamples, beta_dot0, beta_dot1)

quantile(na.omit(sim_PPVs), probs=c(0.025,0.5,0.975))
quantile(na.omit(sim_NPVs), probs=c(0.025,0.5,0.975))
quantile(colMeans(x_matrix), probs=c(0.025,0.5,0.975))

#diagnostic plots
par(mfrow=c(2, 2))
plot(density(na.omit(sim_PPVs)), xlab="Positive Predictive Value", main="A)", col="#DDDDDD")
polygon(density(na.omit(sim_PPVs)), col="#DDDDDD", border=0)

plot(density(na.omit(sim_NPVs)), xlab="Negative Predictive Value", main="B)")
polygon(density(na.omit(sim_NPVs)), col="#DDDDDD", border=0)

plot(density(na.omit(sim_prev)), xlab="Prevalence", main="C)")
polygon(density(na.omit(sim_prev)), col="#DDDDDD", border=0)

## extension to risk factor of race

#regression parameters
beta00 = summary(model3)$coefficients[1, 1]
beta00_se = summary(model3)$coefficients[1, 2]
beta10 = summary(model3)$coefficients[2, 1]
beta10_se = summary(model3)$coefficients[2, 2]
beta01 = summary(model3)$coefficients[3, 1]
beta01_se = summary(model3)$coefficients[3, 2]
beta11 = summary(model3)$coefficients[4, 1]
beta11_se = summary(model3)$coefficients[4, 2]

#generate 1,000 realizations of x
set.seed(777)
nsamples = 1000
x_matrix = matrix(NA, nrow=nrow(retro), ncol=nsamples)
sim_PPV0s=NA; sim_PPV1s=NA; sim_NPV0s=NA; sim_NPV1s=NA; sim_prev=NA; sim_RRs=NA; sim_ll=NA;
for (i in 1:nsamples) {
  #sample parameters
  beta_dot0 = rnorm(1, beta00, beta00_se)
  beta_dot1 = rnorm(1, beta10, beta10_se)
  beta_dot2 = rnorm(1, beta01, beta01_se)
  beta_dot3 = rnorm(1, beta11, beta11_se)
  
  #model of true prevalence conditional on race
  logit.px = beta_dot0 + beta_dot1*retro$hepc + beta_dot2*retro$race + beta_dot3*retro$hepc*retro$race
  px = expit(logit.px)
  
  #generate x (NAs are expected due to missing race; excluded later in complete case analysis)
  x_matrix[, i] = rbinom(nrow(retro), 1, px)
  
  #estimate RR
  model5 = glm(x_matrix[, i] ~ retro$race, family=poisson())
  
  #sample RR from distribution
  sim_RRs = c(sim_RRs, exp(rnorm(1, summary(model5)$coefficients[2,1], summary(model5)$coefficients[2,2])))
  
  #save key estimates
  sim_PPV0s = c(sim_PPV0s, expit(beta_dot0 + beta_dot1))
  sim_PPV1s = c(sim_PPV1s, expit(beta_dot0 + beta_dot1 + beta_dot2 + beta_dot3))
  sim_NPV0s = c(sim_NPV0s, (1 - (expit(beta_dot0))))
  sim_NPV1s = c(sim_NPV1s, (1 - (expit(beta_dot0 + beta_dot2))))
  sim_prev = c(sim_prev, mean(px, na.rm=T))
  sim_ll = c(sim_ll, logLik(model5))
  
}
rm(i, logit.px, px, model5, nsamples, beta_dot0, beta_dot1, beta_dot2, beta_dot3)

quantile(sim_PPV0s, probs=c(0.025,0.5,0.975), na.rm=T)
quantile(sim_PPV1s, probs=c(0.025,0.5,0.975), na.rm=T)
quantile(sim_NPV0s, probs=c(0.025,0.5,0.975), na.rm=T)
quantile(sim_NPV1s, probs=c(0.025,0.5,0.975), na.rm=T)
quantile(sim_RRs, probs=c(0.025,0.5,0.975), na.rm=T)
quantile(sim_prev, probs=c(0.025,0.5,0.975), na.rm=T)

#compare weighted average using likelihood (LL is scaled first)
hist(sim_ll)
mean(sim_RRs,na.rm=T)
weighted.mean(sim_RRs, exp(scale(sim_ll)),na.rm=T)
wtd.quantile(sim_RRs, q=c(0.025,0.5,0.975), weight=exp(scale(sim_ll)))

# #plot log(rr) by log-likelihood
# plot(log(sim_RRs_sampled), sim_ll, xlab="log(RR)", ylab="LL")
# text(1,-380, paste("r =", round(cor(log(sim_RRs), sim_ll, use="complete.obs"),2)))

#diagnostic plots
par(mfrow=c(2, 2))
plot(density(na.omit(sim_PPV0s)), xlab="Positive Predictive Value", main="A)", col="#DDDDDD", xlim=c(0.3,0.8))
polygon(density(na.omit(sim_PPV0s)), col="#DDDDDD", border=0)
lines(density(na.omit(sim_PPV1s)), col="#888888")
polygon(density(na.omit(sim_PPV1s)), col="#888888", border=0)
legend(0.65,18, legend="White", fill="#DDDDDD", cex=1, bty="n", x.intersp=0.2)
legend(0.65,15.5, legend="Non-White", fill="#888888", cex=1, bty="n", x.intersp=0.2)

plot(density(na.omit(sim_NPV0s)), xlab="Negative Predictive Value", main="B)", col="#DDDDDD", xlim=c(0.95,1))
polygon(density(na.omit(sim_NPV0s)), col="#DDDDDD", border=0)
lines(density(na.omit(sim_NPV1s)), col="#888888")
polygon(density(na.omit(sim_NPV1s)), col="#888888", border=0)
legend(0.95,300, legend="White", fill="#DDDDDD", cex=1, bty="n", x.intersp=0.2)
legend(0.95,260, legend="Non-White", fill="#888888", cex=1, bty="n", x.intersp=0.2)

plot(density(na.omit(sim_prev)), xlab="Prevalence", main="C)", col="#DDDDDD")
polygon(density(na.omit(sim_prev)), col="#DDDDDD", border=0)

plot(density(na.omit(sim_RRs)), xlab="Relative Risk", main="D)", col="#DDDDDD")
polygon(density(na.omit(sim_RRs)), col="#DDDDDD", border=0)

#comparison
model6 = glm(hepc ~ race, data=retro, family=poisson())

#for RR must estimate robust standard errors as detailed in Zou, AJE 2004; implementation: https://stats.idre.ucla.edu/r/dae/poisson-regression/
model6.cov.matrix = vcovHC(model6, type="HC0")
model6.std.err = sqrt(diag(model6.cov.matrix))
model6.r.est = cbind(Estimate= coef(model6), "Robust SE" = model6.std.err, "Pr(>|z|)" = 2 * pnorm(abs(coef(model6)/model6.std.err), lower.tail=FALSE), LL = coef(model6) - 1.96 * model6.std.err, UL = coef(model6) + 1.96 * model6.std.err)
exp(model6.r.est[, 1])
exp(model6.r.est[, 4:5])


