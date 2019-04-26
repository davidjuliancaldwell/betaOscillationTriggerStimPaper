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
library('wesanderson')

rootDir = here()

savePlot = 0
figWidth = 8 
figHeight = 6 

chanInt = 14
chanInt1 = paste0(7,chanInt)
chanInt2 = paste0(8,chanInt)

# ------------------------------------------------------------------------
data <- read.table(here("data","output_table","betaStim_outputTable_50.csv"),header=TRUE,sep = ",",stringsAsFactors=F,
                   colClasses=c("magnitude"="numeric","betaLabels"="factor","sid"="factor","numStims"="factor","stimLevel"="numeric","channel"="factor","subjectNum"="factor","phaseClass"="factor","setToDeliverPhase"="factor"))
data <- subset(data, magnitude<1500)
data <- subset(data, magnitude>25)

data <- subset(data,!is.nan(data$magnitude))
#data <- subset(data,data$numStims!='Null')
# rename for ease
data$numStims <- revalue(data$numStims, c("Test 1"="[1,2]","Test 2"="[3,4]","Test 3"="[5,inf)"))
#data$phaseClass <- revalue(data$phaseClass, c("90"=0,"270"=1))

data$percentDiff = 0
data$absDiff = 0
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
        absDiff = data[data$sid == name & data$channel == chan & data$numStims == numStimTrial & data$phaseClass == typePhase,]$magnitude-baseMean
        data[data$sid == name & data$channel == chan & data$numStims == numStimTrial & data$phaseClass == typePhase,]$absDiff = absDiff
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
  scale_fill_manual(name="Experimental\nCondition",
                    breaks=c("0b5a2e", "0b5a2ePlayBack"),
                    labels=c("Closed-loop", "Control"),
                    values=wes_palette(n=2, name="GrandBudapest1")) +
  ylim(0,max(dataSubjChanOnly$magnitude+20)) 
  p

  figHeight = 6
  figWidth = 8
  ggsave(paste0("betaStim_control_subj_7.png"), units="in", width=figWidth, height=figHeight,dpi=600)
  ggsave(paste0("betaStim_control_subj_7.eps"), units="in", width=figWidth, height=figHeight, dpi=600, device=cairo_ps)
  
  
  p2<-ggplot(dataSubjChanOnly, aes(x=numStims, y=magnitude,fill=sid)) + theme_light(base_size = 18) +
    geom_violin(position=position_dodge(1)) +
    labs(x = 'Number of conditioning stimuli',colour = 'closed loop vs. control',title = 'Closed loop vs. control cortical evoked potentials', y = expression(paste("Voltage (",mu,"V)"))) +
    scale_fill_hue(name="Experimental\nCondition",
                   breaks=c("0b5a2e", "0b5a2ePlayBack"),
                   labels=c("Closed-loop", "Control")) + 
    ylim(0,max(dataSubjChanOnly$magnitude+20)) 
  p2
  
# ------------------------------------------------------------------------


fit.lm    = lm(magnitude ~ numStims+sid + numStims:sid,data=dataSubjChanOnly)

summary(fit.lm)
plot(fit.lm)
summary(glht(fit.lm,linfct=mcp(sid="Tukey")))
emmeans(fit.lm, list(pairwise ~ numStims), adjust = "tukey")
emmeans(fit.lm, list(pairwise ~ sid), adjust = "tukey")

emm_s.t <- emmeans(fit.lm, pairwise ~ sid | numStims)
emm_s.t <- emmeans(fit.lm, pairwise ~ numStims | sid)

anova(fit.lm)
tab_model(fit.lm)

summary(glht(fit.lm,linfct=mcp(sid="Tukey")))
summary(glht(fit.lm,linfct=mcp(numStims="Tukey")))
