# ------------------------------------------------------------------------

setwd('/Users/davidcaldwell/code/betaOscillationTriggerStimPaper')

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
library('dplyr')
library('afex')
library('report')
library('effectsize')
library('performance')

# log data prior to fitting?
log_data = FALSE

savePlot = 0
figWidth = 8 
figHeight = 6 

# ------------------------------------------------------------------------
here()
data <- read.table(here("data","output_table","betaStim_outputTable_50_new_100_thresh.csv"),header=TRUE,sep = ",",stringsAsFactors=F,
                   colClasses=c("magnitude"="numeric","betaLabels"="factor","sid"="factor","numStims"="factor","stimLevel"="numeric","channel"="factor","subjectNum"="factor","phaseClass"="factor","setToDeliverPhase"="factor",'phaseDeliveryBinned45'="factor"))

summaryDataCount <- data %>% 
  group_by(sid,setToDeliverPhase,numStims,channel) %>% tally()

data <- subset(data, magnitude<1500)
data <- subset(data, magnitude>25)

data$orderedPhase45 <- factor(data$phaseDeliveryBinned45,levels=c('0-45','45-90','90-135','135-180','180-225','225-270','270-315','315-360'))

#data <- subset(data, magnitude<1000)
#data <- subset(data, magnitude>30)

data <- subset(data,!is.nan(data$magnitude))

# log normnalize
if(log_data){
data$magnitude <- log(data$magnitude)
}

data <- subset(data,data$sid!='702d24')
data <- subset(data,data$sid!='0b5a2ePlayBack')
data <- subset(data,data$numStims!='Null')
# rename for ease
data$numStims <- revalue(data$numStims, c("Test 1"="[1,2]","Test 2"="[3,4]","Test 3"="[5,inf)"))
#data$phaseClass <- revalue(data$phaseClass, c("90"=0,"270"=1))

data$percentDiff = 0
data$absDiff = 0
data$baseMean = 0
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
        # add in base maen
        if (length(absDiff)){
        data[data$sid == name & data$channel == chan & data$numStims == numStimTrial & data$phaseClass == typePhase,]$baseMean = baseMean
        }
      }
    }
  }
}

sapply(data,class)
#summaryData = ddply(data[data$numStims != "Base",] , .(sid,phaseClass,numStims,channel), function(x) mean(x[,"percentDiff"]))

summaryData = ddply(data, .(sid,phaseClass,numStims,channel,betaLabels), summarize, magnitude = mean(magnitude))

summaryDataForMixed = ddply(data, .(sid,phaseClass,numStims,channel,betaLabels), summarize, magnitude = median(magnitude))


summaryData = ddply(data[data$numStims != "Base",] , .(sid,phaseClass,numStims,channel,betaLabels), summarize, percentDiff = mean(percentDiff))

dataNoBaseline = data[data$numStims != "Base",]
dataSubjOnly <- subset(data,data$sid=='0b5a2e')

summaryDataNoPhase = ddply(data, .(sid,numStims,channel,betaLabels), summarize, magnitude = mean(magnitude))
summaryDataNoPhase = ddply(data[data$numStims != "Base",] , .(sid,numStims,channel,betaLabels), summarize, percentDiff = mean(percentDiff))

summaryDataHighStimsOnly = ddply(data[data$numStims == "[5,inf)",] , .(sid,phaseClass,numStims,channel,orderedPhase45), summarize, percentDiff = mean(percentDiff))


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


#### this is the plot for the paper 
pd1 = position_dodge(0.2)
pd2 = position_dodge(0.65)

p2 <- ggplot(summaryData, aes(x=numStims, y=percentDiff,color=phaseClass)) + theme_light(base_size = 14) +
  geom_point(position=position_jitterdodge(dodge.width=0.65, jitter.height=0, jitter.width=0.25),
             alpha=0.7) +
  stat_summary(fun.data=median_hilow,fun.args=(conf.int =0.5), geom="errorbar", width=0.05, position=pd1) +
  stat_summary(fun.y=median, geom="point", size=2, position=pd1) +  
  labs(x = 'Number of Conditioning Stimuli',colour = 'Delivered Phase',title = 'Dose Dependence as a Function of Phase of Stimulation',y = 'Percent Difference from Baseline')+ 
  geom_hline(yintercept=0) +
  scale_color_hue(labels=c("Depolarizing", "Hyperpolarizing")) 
p2
figHeight = 4
figWidth = 8

if(savePlot){
ggsave(paste0("betaStim_dose_phase.png"), units="in", width=figWidth, height=figHeight,dpi=600)
ggsave(paste0("betaStim_dose_phase.eps"), units="in", width=figWidth, height=figHeight, dpi=600, device=cairo_ps)
}


#### this is the other plot for the paper 

colors = c("#AD0000","#FF3535","#FF9999")
  
pd1 = position_dodge(0.2)
pd2 = position_dodge(0.65)

# p2 <- ggplot(summaryData, aes(x=as.numeric(numStims), y=percentDiff,color=as.numeric(numStims))) + theme_light(base_size = 14) +
#       geom_jitter(width=0.2,alpha=0.5,aes(colour=colors)) + geom_smooth(method=lm) +
#   stat_summary(fun.data=median_hilow,fun.args=(conf.int =0.5), geom="errorbar", width=0.1, position=pd1) +
#   stat_summary(fun.y=median, geom="point", size=5, position=pd1) +  
#   labs(x = 'Number of Conditioning Stimuli',title = 'Dose Dependent Change in CEPs',y = 'Percent Difference from Baseline')+ 
#   geom_hline(yintercept=0) +
#   scale_x_continuous(breaks=c(3, 4, 5),labels=c("[1,2]","[3,4]","[5,inf)")) +
#   scale_fill_manual(values = colors)
# p2
# figHeight = 4
# figWidth = 8

p2 <- ggplot(summaryDataNoPhase, aes(x=numStims, y=percentDiff,color=numStims)) + theme_light(base_size = 14) +
  geom_jitter(width=0.2,alpha=0.5) + geom_smooth(method=lm,aes(group=1)) +
  stat_summary(fun.data=median_hilow,fun.args=(conf.int =0.5), geom="errorbar", width=0.1, position=pd1,colour="#666666") +
  stat_summary(fun.y=median, geom="point", size=5, position=pd1,colour="#666666") +  
  labs(x = 'Number of Conditioning Stimuli',title = 'Dose Dependent Change in CEPs',y = 'Percent Difference from Baseline')+ 
  geom_hline(yintercept=0) +
  scale_color_manual(values = colors) + theme(legend.position="none")
p2
figHeight = 4
figWidth = 8


if(savePlot){
  ggsave(paste0("betaStim_dose.svg"), units="in", width=figWidth, height=figHeight,dpi=600)
  ggsave(paste0("betaStim_dose.png"), units="in", width=figWidth, height=figHeight,dpi=600)
  ggsave(paste0("betaStim_dose.eps"), units="in", width=figWidth, height=figHeight, dpi=600, device=cairo_ps)
}

# bar plot by hyper vs depol

p3 <- ggplot(summaryData, aes(x=numStims, y=percentDiff,color=phaseClass)) + theme_light(base_size = 14) +
  stat_summary(fun.data=median_hilow,fun.args=(conf.int =0.5), geom="errorbar", width=0.05, position=pd1) +
  stat_summary(fun.y=median, geom="point", size=2, position=pd1) +  
  labs(x = 'Number of Conditioning Stimuli',colour = 'Delivered Phase',title = 'Dose Dependence as a Function of Phase of Stimulation',y = 'Percent Difference from Baseline')+ 
  geom_hline(yintercept=0) +
  scale_color_hue(labels=c("Depolarizing", "Hyperpolarizing")) 
p3
figHeight = 4
figWidth = 8


if(savePlot){
  ggsave(paste0("betaStim_dose_no_dots.svg"), units="in", width=figWidth, height=figHeight,dpi=600)
  ggsave(paste0("betaStim_dose_no_dots.png"), units="in", width=figWidth, height=figHeight,dpi=600)
  ggsave(paste0("betaStim_dose_no_dots.eps"), units="in", width=figWidth, height=figHeight, dpi=600, device=cairo_ps)
}


# bar plot by which stim 
p4 <- ggplot(summaryDataHighStimsOnly, aes(x=orderedPhase45, y=percentDiff)) + theme_light(base_size = 14) +
  stat_summary(fun.data=median_hilow,fun.args=(conf.int =0.5), geom="errorbar", width=0.05, position=pd1) +
  stat_summary(fun.y=median, geom="point", size=2, position=pd1) +  
  labs(x = 'Binned Phase of Delivery',title = 'Percent Change from Baseline as a Function of Binned Phase',y = 'Percent Difference from Baseline')+ 
  geom_hline(yintercept=0) 
p4
figHeight = 4
figWidth = 8


if(savePlot){
  ggsave(paste0("betaStim_dose_binned45.svg"), units="in", width=figWidth, height=figHeight,dpi=600)
  ggsave(paste0("betaStim_dose_binned45.png"), units="in", width=figWidth, height=figHeight,dpi=600)
  ggsave(paste0("betaStim_dose_binned45.eps"), units="in", width=figWidth, height=figHeight, dpi=600, device=cairo_ps)
}

p2 <- ggplot(summaryData, aes(x=numStims, y=percentDiff,fill=phaseClass)) + 
  geom_boxplot(mapping = aes(x = numStims, y = percentDiff,fill=phaseClass),
               position=position_dodge(0.8),notch=TRUE)  + 
  geom_dotplot(binaxis='y',binwidth=2,stackdir='center', 
               position=position_dodge(0.8))+
  labs(x = 'Number of conditioning stimuli',colour = 'delivered phase',title = 'Dose dependence as a function of phase of stimulation',y = 'Percent difference from baseline')

p2 
p2 + geom_hline(yintercept=0) + theme_classic()

# ========================================================================
# Mixed effects models: dose x phase interaction
# ========================================================================
#
# Hypothesis: dose (numStims) and phase (phaseClass) interact to affect
# CEP magnitude. Higher conditioning doses enhance CEPs more at one
# oscillation phase than the other.
#
# Key data structure notes:
#   - Each trial (row) = one stimulation pulse and its evoked potential
#     on one recording channel
#   - numStims varies trial-to-trial within a channel (different trials
#     received different numbers of conditioning pulses)
#   - phaseClass is a channel-level constant — the circular mean of
#     phase-at-delivery, binned to 90 or 270, replicated across all
#     trials on that channel (see multipleSubj_GLMM_script_PP.m:190-191)
#   - Channel IDs are unique per subject (subjectNum*100 + raw channel),
#     so (1|channel) implicitly nests within subject
#   - 6 subjects, 31 channels, ~37K trials after exclusions
# ========================================================================

# ------------------------------------------------------------------------
# Model 1 (original): random intercepts only
# ------------------------------------------------------------------------
# Random intercept per subject and per channel. No random slopes.
# numStims is tested against trial-level residuals (~37K DF) — this
# overstates significance because subjects may differ in their dose response.

fit.intercepts.only = lmerTest::lmer(
  absDiff ~ numStims * phaseClass + betaLabels + numStims:betaLabels +
  (1|sid/channel),
  data = dataNoBaseline)

summary(fit.intercepts.only)
anova(fit.intercepts.only)
report(fit.intercepts.only)
report(anova(fit.intercepts.only))
model_performance(fit.intercepts.only)
eta_squared(fit.intercepts.only, ci = 0.95)

# emmeans for the original model
emm_orig_dose <- emmeans(fit.intercepts.only, pairwise ~ numStims | phaseClass)
emm_orig_phase <- emmeans(fit.intercepts.only, pairwise ~ phaseClass | numStims)

tab_model(fit.intercepts.only)

figHeight = 4
figWidth = 8
if(savePlot){
  png(here("output_plots","betaStim_residuals_intercepts_only.png"),width=figWidth,height=figHeight,units="in",res=600)
  plot(fit.intercepts.only)
  dev.off()

  png(here("output_plots","betaStim_qq_intercepts_only.png"),width=figWidth,height=figHeight,units="in",res=600)
  qqnorm(resid(fit.intercepts.only)); qqline(resid(fit.intercepts.only))
  dev.off()
}

# ------------------------------------------------------------------------
# Model 2 (trial-level, random slopes): dose slopes per subject
# ------------------------------------------------------------------------
# Adds random dose-response slopes per subject via (0+numStims|sid).
# Each subject can have its own dose effect. numStims is now tested
# against between-subject variability (~5 DF) instead of ~37K.
#
# However, phaseClass DF remains inflated (~4K) because the model lacks
# a random effect at the channel x condition level — within multi-phase
# channels, hundreds of trials per phaseClass are treated as independent
# replication of the phase contrast.

fit.trial.level = lmerTest::lmer(
  absDiff ~ numStims * phaseClass +
  (1|sid) + (0+numStims|sid) + (1|channel),
  data = dataNoBaseline,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=20000)))

summary(fit.trial.level)
anova(fit.trial.level)
VarCorr(fit.trial.level)
report(fit.trial.level)
report(anova(fit.trial.level))
model_performance(fit.trial.level)
eta_squared(fit.trial.level, ci = 0.95)

# dose-response within each phase (raw uV)
emm_trial_dose <- emmeans(fit.trial.level, ~ numStims | phaseClass)
emm_trial_dose
pairs(emm_trial_dose)

# phase contrast at each dose (raw uV)
emm_trial_phase <- emmeans(fit.trial.level, ~ phaseClass | numStims)
pairs(emm_trial_phase)

tab_model(fit.trial.level)

figHeight = 4
figWidth = 8
if(savePlot){
  png(here("output_plots","betaStim_residuals_trial_level.png"),width=figWidth,height=figHeight,units="in",res=600)
  plot(fit.trial.level)
  dev.off()

  png(here("output_plots","betaStim_qq_trial_level.png"),width=figWidth,height=figHeight,units="in",res=600)
  qqnorm(resid(fit.trial.level)); qqline(resid(fit.trial.level))
  dev.off()
}

# ------------------------------------------------------------------------
# Model 3 (summary-level, primary reported model): one observation per cell
# ------------------------------------------------------------------------
# Collapse trials to median per (subject x channel x phaseClass x numStims).
# This eliminates pseudoreplication: each cell contributes one value.
# 120 observations, 6 subjects, 31 channels.
#
# Random effects:
#   (1|sid) — subject-level baseline differences
#   (1|channel) — channel baseline differences (implicitly nested within
#     subject since channel IDs are unique per subject)
#
# Random dose slopes (0+numStims|sid) were removed because with only
# 6 subjects and 3 dose levels, the 3x3 covariance matrix (6 params)
# is near-saturated. The resulting correlations hit 1.0, causing
# singularity. Random intercepts are sufficient: the near-perfect
# slope correlations indicate subjects shift uniformly across doses.
#
# No singularity. No use of setToDeliverPhase as a random grouping
# factor (it is a fixed experimental condition, not a random sample).
#
# DF limitation: phaseClass is a channel-level constant for 22 of 31
# channels, so its ideal denominator DF is ~30 (between-channel).
# Satterthwaite assigns ~84 DF because the model lacks a random effect
# at the channel-condition level. Adding (1|channel:setToDeliverPhase)
# would fix the DF, but is redundant with (1|channel) for single-phase
# channels (22/31), causing singularity. Including both terms requires
# separating channel variance from condition-within-channel variance,
# which the data cannot support with only 9 multi-phase channels.
# This does not affect conclusions: phaseClass is non-significant at
# DF=84 (p=0.42) and would be less significant with fewer DF.

summaryNoBaseline <- ddply(dataNoBaseline, .(sid,phaseClass,numStims,channel),
                           summarize, magnitude = median(magnitude))

cat(sprintf("Summary-level data: %d observations, %d subjects, %d channels\n",
    nrow(summaryNoBaseline), length(unique(summaryNoBaseline$sid)),
    length(unique(summaryNoBaseline$channel))))

fit.summary.level = lmerTest::lmer(
  magnitude ~ numStims * phaseClass +
  (1|sid) + (1|channel),
  data = summaryNoBaseline,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=20000)))

summary(fit.summary.level)
anova(fit.summary.level)
VarCorr(fit.summary.level)
report(fit.summary.level)
report(anova(fit.summary.level))
model_performance(fit.summary.level)
eta_squared(fit.summary.level, ci = 0.95)

# dose-response within each phase (raw uV)
emm_summ_dose <- emmeans(fit.summary.level, ~ numStims | phaseClass)
emm_summ_dose
pairs(emm_summ_dose)

# phase contrast at each dose (raw uV)
emm_summ_phase <- emmeans(fit.summary.level, ~ phaseClass | numStims)
pairs(emm_summ_phase)

tab_model(fit.summary.level)

figHeight = 4
figWidth = 8
if(savePlot){
  png(here("output_plots","betaStim_residuals_summary_level.png"),width=figWidth,height=figHeight,units="in",res=600)
  plot(fit.summary.level)
  dev.off()

  png(here("output_plots","betaStim_qq_summary_level.png"),width=figWidth,height=figHeight,units="in",res=600)
  qqnorm(resid(fit.summary.level)); qqline(resid(fit.summary.level))
  dev.off()
}

# ------------------------------------------------------------------------
# Model 4 (trial-level, nested conditions): correct DF for phaseClass
# ------------------------------------------------------------------------
# Adds (1|channel:setToDeliverPhase) to represent the experimental block
# structure. Each (channel, setToDeliverPhase) combination is one
# condition-within-channel cell. phaseClass is constant within each cell,
# so Satterthwaite correctly identifies it as a cell-level predictor
# and gives it ~31 DF instead of ~4K.
#
# For single-phase subjects, setToDeliverPhase is constant per channel,
# so this term is redundant with (1|channel) — lmer estimates near-zero
# variance for it. For multi-phase subjects, it captures between-condition
# differences within a channel. The fit is singular because of this
# redundancy, but the DF correction is the purpose.
#
# NOTE: This model is kept for reference/comparison. Model 3 (summary-level)
# is the primary reported model — it avoids singularity entirely by
# collapsing to one median per cell.
#
# Nesting: Subject -> Channel -> Condition -> Trial
#   (0+numStims|sid) — dose-response varies by subject
#   (1|channel) — channel baseline
#   (1|channel:setToDeliverPhase) — condition-within-channel (phaseClass level)

fit.nested.condition = lmerTest::lmer(
  absDiff ~ numStims * phaseClass +
  (0+numStims|sid) + (1|channel) + (1|channel:setToDeliverPhase),
  data = dataNoBaseline,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=20000)))

summary(fit.nested.condition)
anova(fit.nested.condition)
VarCorr(fit.nested.condition)
report(fit.nested.condition)
report(anova(fit.nested.condition))
model_performance(fit.nested.condition)
eta_squared(fit.nested.condition, ci = 0.95)

# dose-response within each phase (raw uV)
emm_nested_dose <- emmeans(fit.nested.condition, ~ numStims | phaseClass)
emm_nested_dose
pairs(emm_nested_dose)

# phase contrast at each dose (raw uV)
emm_nested_phase <- emmeans(fit.nested.condition, ~ phaseClass | numStims)
pairs(emm_nested_phase)

tab_model(fit.nested.condition)

figHeight = 4
figWidth = 8
if(savePlot){
  png(here("output_plots","betaStim_residuals_nested_condition.png"),width=figWidth,height=figHeight,units="in",res=600)
  plot(fit.nested.condition)
  dev.off()

  png(here("output_plots","betaStim_qq_nested_condition.png"),width=figWidth,height=figHeight,units="in",res=600)
  qqnorm(resid(fit.nested.condition)); qqline(resid(fit.nested.condition))
  dev.off()
}

# ========================================================================
# Sensitivity analyses: exclusions and within-channel only
# ========================================================================
# Model 4 includes all subjects and conditions. The following models
# test robustness by:
#   4a: Excluding 9ab7ab (no phaseClass contrast — all channels = 270)
#   4b: Also excluding ecb43e random-condition trials (setToDeliverPhase=12345,
#       not phase-targeted but assigned a phaseClass)
#   4c: Restricted to within-channel phaseClass contrasts only — the 9 channels
#       (across 3 subjects) where the same channel was stimulated at both phases.
#       This is the cleanest causal evidence for the phase effect.
#
# All use the same random effects as Model 4:
#   (0+numStims|sid) + (1|channel) + (1|channel:setToDeliverPhase)
# ========================================================================

# ------------------------------------------------------------------------
# Model 4a: exclude 9ab7ab (no phase contrast)
# ------------------------------------------------------------------------
# 9ab7ab has all 5 channels at phaseClass=270 — contributes 7,786 rows
# (21% of data) but zero information about the phase effect.

dataNoBase_no9ab <- subset(dataNoBaseline, sid != '9ab7ab')

cat(sprintf("Model 4a: %d trials, %d subjects, %d channels\n",
    nrow(dataNoBase_no9ab), length(unique(dataNoBase_no9ab$sid)),
    length(unique(dataNoBase_no9ab$channel))))

fit.no.9ab7ab = lmerTest::lmer(
  absDiff ~ numStims * phaseClass +
  (0+numStims|sid) + (1|channel) + (1|channel:setToDeliverPhase),
  data = dataNoBase_no9ab,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=20000)))

summary(fit.no.9ab7ab)
anova(fit.no.9ab7ab)
report(fit.no.9ab7ab)
report(anova(fit.no.9ab7ab))
model_performance(fit.no.9ab7ab)
eta_squared(fit.no.9ab7ab, ci = 0.95)

emm_4a_phase <- emmeans(fit.no.9ab7ab, ~ phaseClass | numStims)
pairs(emm_4a_phase)

emm_4a_dose <- emmeans(fit.no.9ab7ab, ~ numStims | phaseClass)
pairs(emm_4a_dose)

# ------------------------------------------------------------------------
# Model 4b: exclude 9ab7ab + ecb43e random-condition trials
# ------------------------------------------------------------------------
# ecb43e's random condition (setToDeliverPhase=12345) was not phase-locked
# but gets a phaseClass assignment from sinusoidal fits. ~1,082 trials
# dilute the phase effect estimate.

dataNoBase_clean <- subset(dataNoBaseline,
  sid != '9ab7ab' & setToDeliverPhase != '12345')

cat(sprintf("Model 4b: %d trials, %d subjects, %d channels\n",
    nrow(dataNoBase_clean), length(unique(dataNoBase_clean$sid)),
    length(unique(dataNoBase_clean$channel))))

fit.clean = lmerTest::lmer(
  absDiff ~ numStims * phaseClass +
  (0+numStims|sid) + (1|channel) + (1|channel:setToDeliverPhase),
  data = dataNoBase_clean,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=20000)))

summary(fit.clean)
anova(fit.clean)
report(fit.clean)
report(anova(fit.clean))
model_performance(fit.clean)
eta_squared(fit.clean, ci = 0.95)

emm_4b_phase <- emmeans(fit.clean, ~ phaseClass | numStims)
pairs(emm_4b_phase)

emm_4b_dose <- emmeans(fit.clean, ~ numStims | phaseClass)
pairs(emm_4b_dose)

# ------------------------------------------------------------------------
# Model 4c: within-channel phaseClass only (cleanest causal evidence)
# ------------------------------------------------------------------------
# Restrict to the 9 channels where the same physical electrode was
# stimulated under both phaseClass values (from different experimental
# blocks). This eliminates the between-channel confound entirely.
#
# Channels with within-channel phaseClass variation:
#   c91479: 264
#   0b5a2e: 715, 716, 723, 731, 732
#   ecb43e: 647, 648, 655

within_channel_chans <- c('264','715','716','723','731','732','647','648','655')
dataNoBase_within <- subset(dataNoBase_clean,
  channel %in% within_channel_chans)

cat(sprintf("Model 4c: %d trials, %d subjects, %d channels\n",
    nrow(dataNoBase_within), length(unique(dataNoBase_within$sid)),
    length(unique(dataNoBase_within$channel))))

fit.within.channel = lmerTest::lmer(
  absDiff ~ numStims * phaseClass +
  (0+numStims|sid) + (1|channel) + (1|channel:setToDeliverPhase),
  data = dataNoBase_within,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=20000)))

summary(fit.within.channel)
anova(fit.within.channel)
report(fit.within.channel)
report(anova(fit.within.channel))
model_performance(fit.within.channel)
eta_squared(fit.within.channel, ci = 0.95)

emm_4c_phase <- emmeans(fit.within.channel, ~ phaseClass | numStims)
pairs(emm_4c_phase)

emm_4c_dose <- emmeans(fit.within.channel, ~ numStims | phaseClass)
pairs(emm_4c_dose)

# ------------------------------------------------------------------------
# Summary: model comparison
# ------------------------------------------------------------------------
cat("\n=== Denominator DF and p-values across all models ===\n\n")
cat("Model 1  (intercepts only):         numStims df=37020  phaseClass df=4633   interaction df=36667\n")
cat("Model 2  (trial, dose slopes):      numStims df=4.8    phaseClass df=4227   interaction df=851\n")
cat("Model 3  (summary level):           numStims df=84     phaseClass df=84     interaction df=84\n")
cat("Model 4  (trial, nested cond):      numStims df=4.7    phaseClass df=31     interaction df=753\n")
cat("Model 4a (excl 9ab7ab):             see anova above\n")
cat("Model 4b (excl 9ab7ab + random):    see anova above\n")
cat("Model 4c (within-channel only):     see anova above\n")
cat("\nModels 3/4/4a/4b test robustness of the full dataset result.\n")
cat("Model 4c tests whether the effect holds on clean within-channel evidence alone.\n")

# ========================================================================
# Effect sizes: Cohen's d with consistent denominator
# ========================================================================
# All effect sizes use the trial-level residual SD from Model 4 (~77 uV)
# as the denominator. This represents trial-to-trial noise within a single
# channel/condition/dose cell — the variability each individual stimulation
# pulse must overcome. Using one consistent sigma makes d values comparable
# across all models.
#
# Primary reporting is raw effects in uV (from emmeans above).
# Cohen's d is supplementary, for cross-study comparability.
# ========================================================================

sigma_trial <- sigma(fit.nested.condition)
edf_trial <- df.residual(fit.nested.condition)

cat(sprintf("\n=== Effect sizes (Cohen's d, sigma = %.1f uV) ===\n\n", sigma_trial))

# Model 4: dose-response within each phase
cat("--- Model 4 (nested condition): dose within phase ---\n")
eff_size(emm_nested_dose, sigma = sigma_trial, edf = edf_trial)

# Model 4: phase contrast at each dose
cat("--- Model 4 (nested condition): phase at each dose ---\n")
eff_size(emm_nested_phase, sigma = sigma_trial, edf = edf_trial)

# Model 4a: phase contrast at each dose (no 9ab7ab)
cat("--- Model 4a (no 9ab7ab): phase at each dose ---\n")
eff_size(emm_4a_phase, sigma = sigma_trial, edf = edf_trial)

# Model 4b: phase contrast at each dose (no 9ab7ab + random)
cat("--- Model 4b (clean): phase at each dose ---\n")
eff_size(emm_4b_phase, sigma = sigma_trial, edf = edf_trial)

# Model 4c: phase contrast at each dose (within-channel only)
cat("--- Model 4c (within-channel): phase at each dose ---\n")
eff_size(emm_4c_phase, sigma = sigma_trial, edf = edf_trial)
