# ------------------------------------------------------------------------

setwd('C:/Users/david/SharedCode/betaStimPaper')

library('Hmisc')
library('ggplot2')
library('glmm')
library("lme4")
library('multcomp')
library('plyr')
library('here')
library('lmerTest')
library('sjPlot')
library('emmeans')

savePlot = 0
figWidth = 8 
figHeight = 6 

# ------------------------------------------------------------------------

data <- read.table(here("data","output_table","betaStim_outputTable_50.csv"),header=TRUE,sep = ",",stringsAsFactors=F,
                   colClasses=c("magnitude"="numeric","betaLabels"="factor","sid"="factor","numStims"="factor","stimLevel"="numeric","channel"="factor","subjectNum"="factor","phaseClass"="factor","setToDeliverPhase"="factor"))
data <- subset(data, magnitude<1500)
data <- subset(data, magnitude>25)

data <- subset(data,!is.nan(data$magnitude))
data <- subset(data,data$sid!='702d24')
data <- subset(data,data$sid!='0b5a2ePlayBack')
data <- subset(data,data$numStims!='Null')
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

# confirm nothing below 150 uV 
summaryData = ddply(data, .(sid,phaseClass,numStims,channel,betaLabels), summarize, magnitude = mean(magnitude))

summaryData = ddply(data[data$numStims != "Base",] , .(sid,phaseClass,numStims,channel,betaLabels), summarize, percentDiff = mean(percentDiff))

dataNoBaseline = data[data$numStims != "Base",]
dataSubjOnly <- subset(data,data$sid=='0b5a2e')

# ------------------------------------------------------------------------

#data <- read.table(here("Experiment","BetaTriggeredStim","betaStim_outputTable.csv"),header=TRUE,sep = ",",stringsAsFactors=F)
ggplot(data, aes(x=magnitude)) + 
  geom_histogram(binwidth=100)

# # Change box plot colors by groups
# ggplot(data, aes(x=numStims, y=magnitude, fill=phaseClass)) +
#   geom_boxplot()
# Change the position
p<-ggplot(data, aes(x=numStims, y=magnitude, fill=phaseClass)) +
  geom_boxplot(position=position_dodge(1))
p

# Change box plot colors by groups
# ggplot(summaryData, aes(x=numStims, y=percentDiff,fill=phaseClass)) +
#   geom_boxplot(notch=TRUE)
# Change the position
p<-ggplot(summaryData, aes(x=numStims, y=percentDiff,fill=phaseClass)) +
  geom_boxplot(notch=TRUE,position=position_dodge(1)) +
  geom_hline(yintercept=0)
p

p2 <- ggplot(summaryData, aes(x=numStims, y=percentDiff,fill=phaseClass)) + theme_classic(base_size = 18) +
  geom_dotplot(binaxis='y',binwidth=2,stackdir='center', 
               position=position_dodge(0.8)) +
  geom_pointrange(mapping = aes(x = numStims, y = percentDiff,color=phaseClass),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median,
                  position=position_dodge(0.8),size=1.2,color="black",show.legend = FALSE) +  
  labs(x = 'Number of conditioning stimuli',colour = 'delivered phase',title = 'Dose dependence as a function of phase of stimulation',y = 'Percent difference from baseline')+ 
  geom_hline(yintercept=0) 
p2

pd1 = position_dodge(0.2)
pd2 = position_dodge(0.65)

p2 <- ggplot(summaryData, aes(x=numStims, y=percentDiff,color=phaseClass)) + theme_light(base_size = 18) +
  geom_point(position=position_jitterdodge(dodge.width=0.65, jitter.height=0, jitter.width=0.25),
             alpha=0.7) +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=0.05, position=pd1) +
  stat_summary(fun.y=mean, geom="point", size=2, position=pd1) +  
  labs(x = 'Number of conditioning stimuli',colour = 'delivered phase',title = 'Dose dependence as a function of phase of stimulation',y = 'Percent difference from baseline')+ 
  geom_hline(yintercept=0) 
p2

pd1 = position_dodge(0.2)
pd2 = position_dodge(0.65)

p2 <- ggplot(summaryData, aes(x=numStims, y=percentDiff,color=phaseClass)) + theme_light(base_size = 14) +
  geom_point(position=position_jitterdodge(dodge.width=0.65, jitter.height=0, jitter.width=0.25),
             alpha=0.7) +
  stat_summary(fun.data=median_hilow,fun.args=(conf.int =0.5), geom="errorbar", width=0.05, position=pd1) +
  stat_summary(fun.y=median, geom="point", size=2, position=pd1) +  
  labs(x = 'Number of conditioning stimuli',colour = 'delivered phase',title = 'Dose dependence as a function of phase of stimulation',y = 'Percent difference from baseline')+ 
  geom_hline(yintercept=0) +
  scale_color_hue(labels=c("depolarizing", "hyperpolarizing")) 
p2
figHeight = 4
figWidth = 8
ggsave(paste0("betaStim_dose_phase.png"), units="in", width=figWidth, height=figHeight,dpi=600)
ggsave(paste0("betaStim_dose_phase.eps"), units="in", width=figWidth, height=figHeight, dpi=600, device=cairo_ps)


p2 <- ggplot(summaryData, aes(x=numStims, y=percentDiff,fill=phaseClass)) + 
  geom_boxplot(mapping = aes(x = numStims, y = percentDiff,fill=phaseClass),
               position=position_dodge(0.8),notch=TRUE)  + 
  geom_dotplot(binaxis='y',binwidth=2,stackdir='center', 
               position=position_dodge(0.8))+
  labs(x = 'Number of conditioning stimuli',colour = 'delivered phase',title = 'Dose dependence as a function of phase of stimulation',y = 'Percent difference from baseline')

p2 
p2 + geom_hline(yintercept=0) + theme_classic()

# ------------------------------------------------------------------------

#
############ BEST ONE RIGHT NOW
fit.lmm3 = lme4::lmer(percentDiff~numStims+phaseClass + betaLabels + + numStims:betaLabels + numStims:phaseClass + (1 | sid/channel) ,data=summaryData)
#fit.lmm3 = lmerTest::lmer(percentDiff~numStims+phaseClass + betaLabels + (1 | sid/channel) ,data=dataNoBaseline)

#fit.lmm3 = lmerTest::lmer(percentDiff~numStims+phaseClass + betaLabels  + numStims:betaLabels + numStims:phaseClass + (1 | sid/channel) ,data=dataNoBaseline)
fit.lmm3 = lmerTest::lmer(absDiff~numStims+phaseClass + betaLabels  + numStims:betaLabels + numStims:phaseClass + (1 | sid/channel) ,data=dataNoBaseline)

RIaS = unlist(ranef(fit.lmm3))
FixedEff = fixef(fit.lmm3)
emm_s.t <- emmeans(fit.lmm3, pairwise ~ numStims | phaseClass)
emm_s.t <- emmeans(fit.lmm3, pairwise ~ numStims | betaLabels)
anova(fit.lmm3)


tab_model(fit.lmm3)

summary(fit.lmm3)

figHeight = 4
figWidth = 8
png("betaStim_residuals_allSubjs.png",width=figWidth,height=figHeight,units="in",res=600)
plot(fit.lmm3)
dev.off()

setEPS()
postscript("betaStim_residuals_allSubjs.eps",width=figWidth,height=figHeight)
 plot(fit.lmm3)
dev.off()


figHeight = 4
figWidth = 8
png("betaStim_qq_allSubjs.png",width=figWidth,height=figHeight,units="in",res=600)
qqPlot <- qqnorm(resid(fit.lmm3)) 
qqline(resid(fit.lmm3))  #summary(fit.lmm2)dev.off()
dev.off()

setEPS()
postscript("betaStim_qq_allSubjs.eps",width=figWidth,height=figHeight)
qqPlot <- qqnorm(resid(fit.lmm3)) 
qqline(resid(fit.lmm3)) 
dev.off()


summary(glht(fit.lmm3,linfct=mcp(numStims="Tukey")))
summary(glht(fit.lmm3,linfct=mcp(betaLabels="Tukey")))
summary(glht(fit.lmm3,linfct=mcp(phaseClass="Tukey")))

