#################
# Validation of Hep C testing and diagnosis
# Citation: Goldstein ND, Kahal D, Testa K. Accuracy of chronic Hepatitis C diagnosis in the electronic medical record: implications for bias correction approaches. Manuscript in preparation.
# 4/3/19 -- Neal Goldstein
#################

### FUNCTIONS ###

library(psych) #describe, describeBy
library(gmodels) #CrossTable
library(epiR) #sensitivity, specificity, predictive values


### JOIN VALIDATION DATA ###

#read data
load(file="Retrospective.Rdata")
validation = read.csv("validation.csv", as.is=T, stringsAsFactors=F, na.strings="")

#merge
retro = merge(retro, validation, by.x="mrn", by.y="MRN", all.x=T, all.y=F)
retro$Determination[is.na(retro$Determination)] = "Negative"

#save analytic data
retro$mrn = NULL
retro$name = NULL
retro$dob = NULL
save(retro, file="validation.RData")


### LOAD DATA ###

load("validation.RData")


### ANALYSIS ###

describe(retro$last_visit_age); IQR(retro$last_visit_age)
CrossTable(retro$sex)
CrossTable(retro$race==0)
CrossTable(retro$ethnicity)
CrossTable(retro$insurance==1)
describe(retro$n_visits); IQR(retro$n_visits)
CrossTable(retro$hepc_test_ordered)
CrossTable(retro$hepc_test_ab)
CrossTable(!is.na(retro$hepc_test_ab) | !is.na(retro$hepc_test_pcr))
CrossTable(retro$hepc_test_ab)
CrossTable(retro$hepc_test_pcr)
CrossTable(retro$hepc)
CrossTable(retro$Determination=="PCR Positive" | retro$Determination=="Cured")
CrossTable(retro$Determination=="PCR Positive" | retro$Determination=="Cured" | retro$Determination=="Probable" | retro$Determination=="Suspect")

#ICD accuracy           
#strictest definition based on documented labs or treatment
table(retro$Determination)
def1 = CrossTable(retro$hepc, (retro$Determination=="PCR Positive" | retro$Determination=="Cured"), prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
epi.tests(matrix(c(def1$t[4],def1$t[3],def1$t[2],def1$t[1]),nrow=2,ncol=2))

#relaxed definition incporating provider note (but w/out lab)
table(retro$Determination)
def2 = CrossTable(retro$hepc, (retro$Determination=="PCR Positive" | retro$Determination=="Cured" | retro$Determination=="Probable" | retro$Determination=="Suspect"), prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
epi.tests(matrix(c(def2$t[4],def2$t[3],def2$t[2],def2$t[1]),nrow=2,ncol=2))

#antibody accuracy
#strictest definition based on documented labs or treatment
table(retro$Determination)
def1 = CrossTable(retro$hepc_test_ab, (retro$Determination=="PCR Positive" | retro$Determination=="Cured"), prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
epi.tests(matrix(c(def1$t[4],def1$t[3],def1$t[2],def1$t[1]),nrow=2,ncol=2))

#relaxed definition incporating provider note (but w/out lab)
table(retro$Determination)
def2 = CrossTable(retro$hepc_test_ab, (retro$Determination=="PCR Positive" | retro$Determination=="Cured" | retro$Determination=="Probable" | retro$Determination=="Suspect"), prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
epi.tests(matrix(c(def2$t[4],def2$t[3],def2$t[2],def2$t[1]),nrow=2,ncol=2))

#explore predictors of misclassification using confirmed or suspected definition
emrdx = retro
emrdx$def1 = ifelse(emrdx$Determination=="PCR Positive" | emrdx$Determination=="Cured", 1, 0)
emrdx$def2 = ifelse(emrdx$Determination=="PCR Positive" | emrdx$Determination=="Cured" | emrdx$Determination=="Probable" | emrdx$Determination=="Suspect", 1, 0)
emrdx$misclass = ifelse(emrdx$hepc != emrdx$def2, 1, 0)

summary(glm(misclass ~ last_visit_age, data=emrdx, family=binomial()))
summary(glm(misclass ~ sex, data=emrdx, family=binomial()))
summary(glm(misclass ~ (race==0), data=emrdx, family=binomial()))
summary(glm(misclass ~ ethnicity, data=emrdx, family=binomial()))
summary(glm(misclass ~ (insurance==1), data=emrdx, family=binomial()))
summary(glm(misclass ~ n_visits, data=emrdx, family=binomial()))
summary(glm(misclass ~ hepc_test_ordered, data=emrdx, family=binomial()))

age_model = glm(misclass ~ last_visit_age, data=emrdx, family=binomial())
exp(coef(age_model))
exp(confint(age_model))

ethnicity_model = glm(misclass ~ ethnicity, data=emrdx, family=binomial())
1/exp(coef(ethnicity_model))
1/exp(confint(ethnicity_model))

#differential misclassification checks
0.96 - 0.96*.05
0.96 - 0.96*.1
0.96 - 0.96*.25
