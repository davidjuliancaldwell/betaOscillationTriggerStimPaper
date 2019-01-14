# ------------------------------------------------------------------------
setwd('C:/Users/david/SharedCode/betaStimPaper')

library('Hmisc')
library('nlme')
library('ggplot2')
library('drc')
library('minpack.lm')
library('lmtest')
library('glmm')
library("lme4")
library('multcomp')
library('plyr')
library('here')
library('lmerTest')
library('sjPlot')
library('emmeans')

rootDir = here()

savePlot = 0
figWidth = 8 
figHeight = 6 

chanInt = 14
chanInt1 = paste0(7,chanInt)
chanInt2 = paste0(8,chanInt)

# ------------------------------------------------------------------------
data <- read.table(here("Experiment","BetaTriggeredStim","betaStim_outputTable_50.csv"),header=TRUE,sep = ",",stringsAsFactors=F,
                   colClasses=c("magnitude"="numeric","betaLabels"="factor","sid"="factor","numStims"="factor","stimLevel"="numeric","channel"="factor","subjectNum"="factor","phaseClass"="factor","setToDeliverPhase"="factor"))
data <- subset(data, magnitude<1500)
data <- subset(data, magnitude>25)

data <- subset(data,!is.nan(data$magnitude))
#data <- subset(data,data$numStims!='Null')
# rename for ease
data$numStims <- revalue(data$numStims, c("Test 1"="[1,2]","Test 2"="[3,4]","Test 3"="[5,inf)"))
#data$phaseClass <- revalue(data$phaseClass, c("90"=0,"270"=1))

data$percentDiff = 0
for (name in unique(data$sid)){
  for (chan in unique(data[data$sid == name,]$channel)){
    for (numStimTrial in unique(data$numStims)){
      numBase = nrow(data[data$sid == name & data$channel == chan & data$numStims == 'Base',])
      base = data[data$sid == name & data$channel == chan & data$numStims == 'Base',]$magnitude
      baseMean = mean(base)
      data[data$sid == name & data$channel == chan & data$numStims == 'Base',]$percentDiff = 100*(base - baseMean)/baseMean
      for (typePhase in unique(data$phaseClass)){
        percentDiff = 100*((data[data$sid == name & data$channel == chan & data$numStims == numStimTrial & data$phaseClass == typePhase,]$magnitude)-baseMean)/baseMean
        data[data$sid == name & data$channel == chan & data$numStims == numStimTrial & data$phaseClass == typePhase,]$percentDiff = percentDiff
      }
    }
  }
}

sapply(data,class)
#summaryData = ddply(data[data$numStims != "Base",] , .(sid,phaseClass,numStims,channel), function(x) mean(x[,"percentDiff"]))

dataNoBaseline = data[data$numStims != "Base",]
dataSubjOnly <- subset(data,data$sid=='0b5a2e' | data$sid=='0b5a2ePlayBack')
dataSubjChanOnly <- subset(dataSubjOnly,dataSubjOnly$channel == chanInt1 | dataSubjOnly$channel == chanInt2)
#summaryData = ddply(dataSubjOnly[dataSubjOnly$numStims != "Base",] , .(sid,phaseClass,numStims,channel,betaLabels), summarize, percentDiff = mean(percentDiff))
#summaryData = ddply(dataSubjOnly, .(sid,phaseClass,numStims,channel,betaLabels), summarize, percentDiff = mean(percentDiff))
summaryData = ddply(dataSubjOnly, .(sid,phaseClass,numStims,channel,betaLabels), summarize, meanMag = mean(magnitude), sdMag = sd(magnitude))

summaryDataChan = subset(summaryData,summaryData$chan == chanInt1 | summaryData$chan == chanInt2)
# ------------------------------------------------------------------------


# Change box plot colors by groups
# ggplot(summaryData, aes(x=numStims, y=percentDiff,fill=phaseClass)) +
#   geom_boxplot(notch=TRUE)
# Change the position
p<-ggplot(dataSubjChanOnly, aes(x=numStims, y=magnitude,fill=sid)) + theme_light(base_size = 18) +
  geom_boxplot(notch=TRUE,position=position_dodge(1)) +
  labs(x = 'Number of conditioning stimuli',colour = 'closed loop vs. control',title = 'Closed loop vs. control cortical evoked potentials', y = expression(paste("Voltage (",mu,"V)"))) +
  scale_fill_hue(name="Experimental\nCondition",
                      breaks=c("0b5a2e", "0b5a2ePlayBack"),
                      labels=c("closed-loop", "control")) + 
  ylim(0,max(dataSubjChanOnly$magnitude+20))
  p

  
  p2<-ggplot(dataSubjChanOnly, aes(x=numStims, y=magnitude,fill=sid)) + theme_light(base_size = 18) +
    geom_violin(position=position_dodge(1)) +
    labs(x = 'Number of conditioning stimuli',colour = 'closed loop vs. control',title = 'Closed loop vs. control cortical evoked potentials', y = expression(paste("Voltage (",mu,"V)"))) +
    scale_fill_hue(name="Experimental\nCondition",
                   breaks=c("0b5a2e", "0b5a2ePlayBack"),
                   labels=c("closed-loop", "control")) + 
    ylim(0,max(dataSubjChanOnly$magnitude+20))
  p2
  
  

# ------------------------------------------------------------------------


fit.lm    = lm(magnitude ~ numStims+sid + numStims*sid,data=dataSubjChanOnly)

summary(fit.lm)
plot(fit.lm)
summary(glht(fit.lm,linfct=mcp(sid="Tukey")))
emmeans(fit.lm, list(pairwise ~ numStims), adjust = "tukey")
emmeans(fit.lm, list(pairwise ~ sid), adjust = "tukey")

emm_s.t <- emmeans(fit.lm, pairwise ~ sid | numStims)
emm_s.t <- emmeans(fit.lm, pairwise ~ numStims | sid)


tab_model(fit.lm)
tab_model(
  m1, m2, 
  pred.labels = c("Intercept", "Age (Carer)", "Hours per Week", "Gender (Carer)",
                  "Education: middle (Carer)", "Education: high (Carer)", 
                  "Age (Older Person)"),
  dv.labels = c("First Model", "M2"),
  string.pred = "Coeffcient",
  string.ci = "Conf. Int (95%)",
  string.p = "P-Value"
)

summary(glht(fit.lm,linfct=mcp(sid="Tukey")))
summary(glht(fit.lm,linfct=mcp(numStims="Tukey")))

p <- ggplot(dataNoBaseline, aes(x=numStims, y=percentDiff, colour=phaseClass)) +
  geom_point(size=3) +
  geom_line(aes(y=predict(fm2), group=Subject, size="Subjects")) +
  geom_line(data=newdat, aes(y=predict(fm2, level=0, newdata=newdat), size="Population")) +
  scale_size_manual(name="Predictions", values=c("Subjects"=0.5, "Population"=3)) +
  theme_bw(base_size=22) 
print(p)

p2 <- ggplot(dataNoBaseline, aes(x=numStims, y=percentDiff,fill=phaseClass)) + theme_classic(base_size = 18) +
  geom_boxplot(binaxis='y',binwidth=2,stackdir='center', 
               position=position_dodge(0.8)) +
  geom_dotplot(data = dataNoBaseline, aes(y=predict(fit.lmm2,level=0,newdata = dataNoBaseline))) +
  labs(x = 'Number of conditioning stimuli',colour = 'delivered phase',title = 'Dose dependence as a function of phase of stimulation',y = 'Percent difference from baseline')+ 
  geom_hline(yintercept=0) 
p2

#fit.lm    = glm(magnitude ~ stimLevel + numStims + subjectNum + channel + phaseClass,data=data)
#fit.lm    = glm(magnitude ~ numStims + channel + phaseClass,data=data)
fit.lm2    = glm(magnitude ~ numStims+phaseClass+betaLabels,data=data)

summary(fit.lm2)
plot(fit.lm2)
summary(glht(fit.lm2,linfct=mcp(phaseClass="Tukey")))
summary(glht(fit.lm2,linfct=mcp(numStims="Tukey")))


#fit.lmm = lmer(magnitude ~ stimLevel + numStims + subjectNum + channel + phaseClass + (1|subjectNum) + (1|numStims) + (1|channel)+(1|stimLevel),data=data)
fit.lmm = lmer(percentDiff~numStims+phaseClass+betaLabels+channel+ (1| sid) ,data=dataNoBaseline)

summary(fit.lmm)
plot(fit.lmm)
#confint(fit.lmm,method="boot")
summary(glht(fit.lmm,linfct=mcp(numStims="Tukey")))
summary(glht(fit.lmm,linfct=mcp(betaLabels="Tukey")))
summary(glht(fit.lmm,linfct=mcp(phaseClass="Tukey")))

fit.lmm2 = lmer(percentDiff~numStims+phaseClass+betaLabels + numStims * phaseClass + (1 | sid/channel) ,data=dataNoBaseline)
fit.lmm2 = lmer(percentDiff~numStims+phaseClass+betaLabels  + (1 | sid/channel)  ,data=dataNoBaseline)

RIaS = unlist(ranef(fit.lmm2))
FixedEff = fixef(fit.lmm2)
summary(fit.lmm2)
#fit.lmm2 = lmer(percentDiff~numStims+phaseClass+betaLabels + (1 + numStims|sid/channel) + (phaseClass|sid/channel) ,data=dataNoBaseline)
plot(fit.lmm2)
qqnorm(resid(fit.lmm2))
qqline(resid(fit.lmm2))  #summary(fit.lmm2)
#confint(fit.lmm2,method="boot")
summary(glht(fit.lmm2,linfct=mcp(numStims="Tukey")))
summary(glht(fit.lmm2,linfct=mcp(betaLabels="Tukey")))
summary(glht(fit.lmm2,linfct=mcp(phaseClass="Tukey")))

#

#fit.lmm3 = lmer(percentDiff~numStims+phaseClass + betaLabels +  (1 | sid/channel) ,data=summaryData)
fit.lmm3 = lmer(percentDiff~numStims+phaseClass + betaLabels  + numStims*betaLabels + numStims*phaseClass + (1 | sid/channel) ,data=dataNoBaseline)

qqnorm(resid(fit.lmm3))
qqline(resid(fit.lmm3))  #summary(fit.lmm2)
summary(glht(fit.lmm3,linfct=mcp(numStims="Tukey")))
summary(glht(fit.lmm3,linfct=mcp(betaLabels="Tukey")))
summary(glht(fit.lmm3,linfct=mcp(phaseClass="Tukey")))

fit.lmm4 = lmer(percentDiff~numStims+phaseClass + betaLabels +  (1 | channel) ,data=dataSubjOnly)
summary(fit.lmm4)
qqnorm(resid(fit.lmm4))
qqline(resid(fit.lmm4))  #summary(fit.lmm2)
summary(glht(fit.lmm4,linfct=mcp(numStims="Tukey")))
summary(glht(fit.lmm4,linfct=mcp(betaLabels="Tukey")))
summary(glht(fit.lmm4,linfct=mcp(phaseClass="Tukey")))

# xtra ss test
anova(fit.lm, fit.lmm,fit.lmm2)
# Likelihood ratio test

lrtest(fit.lm,fit.lmm,fit.lmm2)

# aic
AIC(fit.lm,fit.lmm,fit.lmm2)

BIC(fit.lm,fit.lmm,fit.lmm2)
