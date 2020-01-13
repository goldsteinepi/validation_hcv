#################
# Validation of Hep C testing and diagnosis
# Citation: Goldstein ND, Kahal D, Testa K, Burstyn I. Accuracy of chronic Hepatitis C diagnosis in the electronic medical record: implications for recalling patients for treatment. Manuscript in preparation.
# 4/3/19 -- Neal Goldstein
#################

### FUNCTIONS ###

library(psych) #describe, describeBy
library(gmodels) #CrossTable
library(epiR) #sensitivity, specificity, predictive values


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


### LOAD DATA ###

load("validation.RData")


### DETERMINE INFECTION STATUS ###

#ignore HCV+ lab results for 2 individuals without corresponding HCV diagnostic code: we do not know the reason for the test
retro$Determination[which(retro$hepc==0 & (retro$Determination=="PCR Positive" | retro$Determination=="Cured"))] = "Negative"

#remove W=1 individuals from universal screening era (as they would have been screened under risk factor paradigm); that is, these individuals are W=0
pro = pro[pro$riskfactor==0 & pro$hepc_old==0, ]

#create a diagnostic variable for the universal screening group
pro$Determination = ifelse((!is.na(pro$hepc_test_pcr) & (pro$hepc_test_pcr==1)), 1, 0)


### ANALYSIS ###

describe(retro$last_visit_age); IQR(retro$last_visit_age)
CrossTable(retro$sex)
CrossTable(retro$race==0)
CrossTable(retro$ethnicity)
CrossTable(retro$insurance==1)
CrossTable(retro$hepc_test_ordered)
CrossTable(retro$hepc_test_resulted)
CrossTable(retro$hepc)
CrossTable(retro$Determination=="PCR Positive" | retro$Determination=="Cured")

describe(pro$last_visit_age); IQR(pro$last_visit_age)
CrossTable(pro$sex)
CrossTable(pro$race==0)
CrossTable(pro$ethnicity)
CrossTable(pro$insurance==1)
CrossTable(pro$hepc_test_ordered)
CrossTable(pro$hepc_test_resulted)
CrossTable(pro$Determination)

describe(pro$last_visit_age[pro$hepc_test_resulted==1]); IQR(pro$last_visit_age[pro$hepc_test_resulted==1])
CrossTable(pro$sex[pro$hepc_test_resulted==1])
CrossTable(pro$race[pro$hepc_test_resulted==1]==0)
CrossTable(pro$ethnicity[pro$hepc_test_resulted==1])
CrossTable(pro$insurance[pro$hepc_test_resulted==1]==1)
#CrossTable(pro$hepc_test_ordered[pro$hepc_test_resulted==1])
CrossTable(pro$Determination[pro$hepc_test_resulted==1])

describe(pro$last_visit_age[pro$hepc_test_resulted==0]); IQR(pro$last_visit_age[pro$hepc_test_resulted==0])
CrossTable(pro$sex[pro$hepc_test_resulted==0])
CrossTable(pro$race[pro$hepc_test_resulted==0]==0)
CrossTable(pro$ethnicity[pro$hepc_test_resulted==0])
CrossTable(pro$insurance[pro$hepc_test_resulted==0]==1)
CrossTable(pro$hepc_test_ordered[pro$hepc_test_resulted==0])
#CrossTable(pro$Determination[pro$hepc_test_resulted==0])

#compute validation metrics for risk-based cohort
riskcohort = CrossTable(retro$hepc, (retro$Determination=="PCR Positive" | retro$Determination=="Cured"), prop.r=F, prop.t=F, prop.chisq=F, chisq=T)$t
riskcohort

#re-organize table and calculate PPV
riskcohort = matrix(c(riskcohort[4],riskcohort[3],riskcohort[2],riskcohort[1]),nrow=2,ncol=2)
epi.tests(riskcohort)

#data for universal screening cohort
universalcohort = matrix(c(0,5,0,(341-5)),nrow=2,ncol=2)
epi.tests(universalcohort)
NPV = c(0.97,0.99,1.00)

#re-classify individuals in the risk-based cohort using the mean and bounds
riskcohort_mean = riskcohort
riskcohort_mean[2,1] = round(sum(riskcohort_mean[2,]) * (1-NPV[2])) # Pr(X=0|W=1)
riskcohort_mean[2,2] = round(sum(riskcohort_mean[2,]) * (NPV[2])) # Pr(X=0|W=1)
riskcohort_mean
epi.tests(riskcohort_mean)

riskcohort_lo = riskcohort
riskcohort_lo[2,1] = round(sum(riskcohort_lo[2,]) * (1-NPV[1])) # Pr(X=0|W=1)
riskcohort_lo[2,2] = round(sum(riskcohort_lo[2,]) * (NPV[1])) # Pr(X=0|W=1)
riskcohort_lo
epi.tests(riskcohort_lo)

riskcohort_hi = riskcohort
riskcohort_hi[2,1] = round(sum(riskcohort_hi[2,]) * (1-NPV[3])) # Pr(X=0|W=1)
riskcohort_hi[2,2] = round(sum(riskcohort_hi[2,]) * (NPV[3])) # Pr(X=0|W=1)
riskcohort_hi
epi.tests(riskcohort_hi)


