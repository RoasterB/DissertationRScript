library(tidyverse)
library(lubridate)
library(epitools)
library(finalfit)
library(broom)
library(yardstick)
library(naniar)
library(survival)
library(survminer)
library(survMisc)
library(forcats)
library(forestplot)
library(emmeans)
library(paletteer)

data <-read.csv('***.csv')
#Create factors and levels, approx age, age groups, deprivation levels,
#time to discharge(can't be 0 so add 1), urgency and discharged under 18yrs.
#If no discharge, add an end of study date

datafactors<-data%>%
  mutate(Gender = factor(Gender),
         Ethnicity = factor(Ethnicity, levels = c('White - British', 'White - Irish', 'White - Any other White background',
                                                  'Asian or Asian British - Bangladeshi','Asian or Asian British - Indian',
                                                  'Asian or Asian British - Pakistani', 'Asian or Asian British - Any other Asian background',
                                                  'Other Ethnic Groups - Chinese','Black or Black British - African',
                                                  'Black or Black British - Caribbean', 'Black or Black British - Any other Black background', 
                                                  'Mixed - White and Black Caribbean', 'Mixed - White and Asian',
                                                  'Mixed - White and Black African', 'Mixed - Any other mixed background',
                                                  'Other Ethnic Groups - Any other ethnic group' )),
         ApproxAge = ifelse(MonthRef > 06, YearRef-YearOfBirth, YearRef-YearOfBirth-1),
         AgeGroup= case_when(ApproxAge <13~ 'Below 13yrs',
                             ApproxAge <=15 & ApproxAge >=13 ~ '13-15yrs',
                             ApproxAge > 15 ~ '16-18yrs'),
         AgeGroup= factor(AgeGroup, levels = c('Below 13yrs','13-15yrs', '16-18yrs')),
         DepGrp = case_when(IDACIDecile %in% c( '10','9')~ 'Most Affluent',
                            IDACIDecile %in% c( '8','7','6')~ 'Fairly Affluent',
                            IDACIDecile %in% c('5','4','3')~ 'Fairly Deprived',
                            IDACIDecile %in% c('2','1')~'Most Deprived'),
         DepGrp = factor(DepGrp, levels = c('Most Affluent','Fairly Affluent','Fairly Deprived','Most Deprived' )),
         DateReferral=as.Date(DateReferral),
         X1stApptDate=as.Date(X1stApptDate),
         DateDischarge= as.Date(DateDischarge),
         TimeToDis = DateDischarge - X1stApptDate+1,
         TimeToDis = as.numeric(TimeToDis)
         ReferralSource= replace(ReferralSource, ReferralSource %in% c('Deliberate Self Harm/Out of Hours','A&E',
                                                                              'STAR Assessment Team/Out of hours',
                                                                       'Adult Mental Health'), 'Emergency'),
         Urgency = ifelse((LengthWait < 21|ReferralSource == 'Emergency', "1", "0"),
         YearRef = format(as.Date(DateReferral, format="%d/%m/%Y"),"%Y"),
         YearRef = as.numeric(YearRef),
         MonthRef = format(as.Date(DateReferral, format="%d/%m/%Y"),"%m"),
         YearDis = format(as.Date(DateDischarge, format="%d/%m/%Y"),"%Y"),
         YearDis = as.numeric(YearDis),
         DisAge = YearDis -YearOfBirth),
         Dischargedunder18 = ifelse(DisAge>= 18|Discharged==0, '0','1'),
         Dischargedunder18=as.numeric(Dischargedunder18),
         DateDischarge=if_else(is.na(DateDischarge),
                               ymd("2026-01-28"),DateDischarge))
           )
#Summary to find missing values and check variables
summary(datafactors)

#Create a missing values dataset
missingvalues<-datafactors%>% 
  filter(is.na(Ethnicity)|is.na(Gender)|is.na(IDACIDecile))
#Remove missing values
nomissing<-datafactors%>%
  filter(!is.na(Ethnicity),!is.na(Gender),!is.na(IDACIDecile))

#Ethnicity counts with urgent referral
table(df$Ethnicity, df$Urgency)

#Keep only those with >5 urgents
df<-nomissing%>%
  filter(Ethnicity %in%c('White - British',
                         'Asian or Asian British - Pakistani',
                         'Black or Black British - African'))%>%
  mutate(Ethnicity = factor(Ethnicity))

#Overview
summary(df)
#Proportions 
df%>%
  group_by(Gender) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))#repeated for each confounder

#Relationship between potential confounders and ethnicity

#Chisq if counts in cells of table >5 otherwise use Fishers

table(df$Ethnicity,df$Gender)
fisher.test(table(df$Ethnicity,df$Gender))
chisq.test(table(df$Ethnicity,df$Gender))#repeated for each confounder

#Visualisation
propethID<-df%>%
  group_by(Ethnicity, DepGrp) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))
ggplot(propethID,aes(x=Ethnicity,y=prop,fill = DepGrp))+
  geom_col(position = 'dodge')+
  scale_fill_paletteer_d("MexBrewer::Frida")+
  theme_bw()+
  labs(y='Proportion', 
       title = 'Proportion in Deprivation Groups by Ethnicity', fill= 'Deprivation Level')

#Relationship between potential confounders and urgency
#Chisq or Fishers as appropriate
table(df$Gender, df$Urgency)
chisq.test(table(df$Ethnicity, df$Urgency))
fisher.test(table(df$Ethnicity, df$Urgency))#repeat for other confounders

#Visualisation
propgenurg<-df%>%
  group_by(Gender,Urgency) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))
ggplot(propgenurg,aes(x=Gender,y=prop,fill = Urgency))+
  geom_col(position = 'dodge')+
  scale_fill_paletteer_d("Manu::Kereru", labels = c('Routine', 'Urgent'))+
  theme_bw()+
  labs(y='Proportion', 
       title = 'Proportion with Urgent Referral by Gender')+
  scale_x_discrete(label = c('Female', 'Male'))


# Ethnicity and urgency
table(df$Ethnicity, df$Urgency)
chisq.test(table(df$Ethnicity, df$Urgency))
fisher.test(table(df$Ethnicity, df$Urgency))

#Crude ORs
genderOR <- glm(Urgency ~ Gender, 
                data = df, 
                family = binomial)
tidy(genderOR,exponentiate = TRUE, conf.int=TRUE)#repeated for each confounder

ethOR <- glm(Urgency ~ Ethnicity, 
             data = ltdeth, 
             family = binomial)
tidy(ethOR,exponentiate = TRUE, conf.int=TRUE)

#Full logistic regression model
urgencyall_relationship <- glm(Urgency ~ Ethnicity+AgeGroup+Gender+DepGrp, 
                               data = df, 
                               family = binomial)

tidy(urgencyall_relationship,exponentiate = TRUE, conf.int=TRUE)

#Create tibble of data
ORdata <- tibble(mean  = c(0.85, 1.04, 2.35, 2.30, 3.62, 3.53, 
                           1.21, 1.18, 0.82, 0.83, 0.91, 0.92, 
                           1.61, 1.57, 4.03, 3.88),
                 lower = c(0.65, 0.79, 1.68, 1.63, 2.52, 2.43, 0.82, 0.80, 0.55, 0.55, 0.65, 0.65, 0.85, 0.82, 1.49, 1.41 ),
                 upper = c(1.11, 1.36, 3.33, 3.27, 5.26, 5.16, 1.79, 1.76, 1.23, 1.26, 1.30, 1.32, 2.84, 2.81, 10.10, 9.96),
                 labeltext = c("Gender- Male", "Gender- Male", "Age - 13-15yrs", "Age - 13-15yrs","        - 16-18yrs","        - 16-18yrs",
                               "Deprivation - Fairly affluent", "Deprivation - Fairly affluent","                - Fairly deprived","                - Fairly deprived","                - Most deprived", "                - Most deprived",
                               "Ethnicity- Asian Pakistani", "Ethnicity- Asian Pakistani", "Ethnicity- Black African", "Ethnicity- Black African"),
                 model = c("Crude", "Adjusted","Crude", "Adjusted","Crude", "Adjusted",
                           "Crude", "Adjusted","Crude", "Adjusted","Crude", "Adjusted",
                           "Crude", "Adjusted","Crude", "Adjusted"),
                 urgent = c("98", " ", "128", " ", "89", " ",
                            "60", " ", "50", " ","100", " ",
                            "14", " ", "7", " " ),
                 routine = c("718", " ", "744", " ", "336", " ",
                             "318", " ", "388", " ", "700", " ",
                             "60", " ", "12", " "),
                 OR = c("0.85", "1.04", "2.35", "2.30", "3.62", "3.53",
                        "1.21", "1.18", "0.82", "0.35", "0.91", "0.92",
                        "1.61", "1.57", "4.03", "3.88"))
#Create forestplot of both crude and adjusted

base_data %>%
  group_by(model) %>%
  forestplot(clip = c(0, 8),
             ci.vertices = TRUE,
             boxsize = .15,
             lineheight = "lines",
             legend_args = fpLegend(pos = list(x = .85, y = 0.4)),
             xlab="Odds Ratio")%>%
  fp_add_lines(h_2 = gpar(lty = 1, col ="steelblue"),
               h_3 = gpar(lty = 1,col ="steelblue"), 
               h_5= gpar(lty = 1,col ="steelblue"), 
               h_8 = gpar(lty = 1,col ="steelblue"))%>%
  fp_add_header("Variable") %>%
  fp_set_style(box = c("blue", "darkred") %>%
                 lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))%>%
  fp_decorate_graph(grid = structure(1))
# Predicted probs
predp<-emmeans(urgencyall_relationship,~ Ethnicity, type = 'response')
predp
#Survival analysis
table(nomissing$Ethnicity,nomissing$Dischargedunder18)

#Keep ethnicities with >10 discharges
survdata<-nomissing%>%
  filter(Ethnicity %in% c('White - British',
                          'Asian or Asian British - Pakistani',
                          'Black or Black British - African',
                          'Black or Black British - Caribbean'))%>%
  mutate(Ethnicity = factor(Ethnicity))

#Check variables
summary(survdata)

#remove anomalies
survdata<-survdata%>%
  filter(TimetoDis>0)

#Confounding- relationship with ethnicity repeat steps in log regression but including new ethnicity
#Create survival object
survival<-Surv(survdata$LengthActive, survdata$Dischargedunder18)
#KM Survival in confounders
#repeat for each confounder
kmfitg <- survfit(survival ~ Gender, data = limitsurv)
kmfitg
#KM curve
ggsurvplot(
  kmfitg,
  data = survdata,
  pval = TRUE,
  xlab='Time since first appointment(days)',
  ylab= ' Proportion not discharged',
  legend='right',
  legend.title = "Gender",
  legend.labs = c('Female',
                  'Male'),
  title = 'Kaplan-Meier curve for Time to Discharge by Gender')
#table of medians
summary(kmfitg)$table[ , "median"]  
#Survival for ethnicity
#KM as done for confounders

#Cox proportional hazards
#For each confounder

coxgen<-coxph(
  survival ~ Gender,
  data=survdata)
summary(coxgen)
#Test proprotional hazards
phtestgen<-cox.zph(coxgen)
phtestgen
#view reisdual plots
plot(phtestgen)
#Main explanatory variable fo Ethnicity
coxeth<-coxph(
  survival ~ Ethnicity,
  data=survdata)
summary(coxeth)
phtest<-cox.zph(coxeth)
phtest
plot(phtest)
#Model including all variables
coxmodel <- coxph(
  survival ~ Ethnicity + Gender + AgeGroup + DepGrp,
  data = survdata
)
summary(coxmodel)
finalphtest<-cox.zph(coxmodel)
plot(finalphtest)

#Visualisation of full cox for ethnicity

final_plot <- ggadjustedcurves(
  coxmodel,
  data = survdata,
  variable = "Ethnicity",        
  method = "average",            
  legend.title = "Ethnicity",
  xlab = "Time since first appointment (days)",
  ylab = "Proportion not discharged",
  title = 'Cox Proportional Hazards Model: Adjusted Survival Estimates by Ethnicity'
) + theme_minimal()
final_plot

