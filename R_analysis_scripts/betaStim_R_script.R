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
library('officer')
library('flextable')

# log data prior to fitting?
log_data = FALSE

savePlot = 1
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

# 702d24 previously excluded (1 channel only) — now included
data <- subset(data,data$sid!='0b5a2ePlayBack')
data <- subset(data,data$numStims!='Null')
# rename for ease
data$numStims <- revalue(data$numStims, c("Test 1"="[1,2]","Test 2"="[3,4]","Test 3"="[5,inf)"))
#data$phaseClass <- revalue(data$phaseClass, c("90"=0,"270"=1))

data$percentDiff = 0
data$absDiff = 0
data$baseMedian = 0
for (name in unique(data$sid)){
  for (chan in unique(data[data$sid == name,]$channel)){
    for (numStimTrial in unique(data$numStims)){
      numBase = nrow(data[data$sid == name & data$channel == chan & data$numStims == 'Base',])
      base = data[data$sid == name & data$channel == chan & data$numStims == 'Base',]$magnitude
      baseMedian = median(base)
      data[data$sid == name & data$channel == chan & data$numStims == 'Base',]$percentDiff = 100*(base - baseMedian)/baseMedian
      for (typePhase in unique(data$phaseClass)){
        percentDiff = 100*((data[data$sid == name & data$channel == chan & data$numStims == numStimTrial & data$phaseClass == typePhase,]$magnitude)-baseMedian)/baseMedian
        data[data$sid == name & data$channel == chan & data$numStims == numStimTrial & data$phaseClass == typePhase,]$percentDiff = percentDiff
        absDiff = data[data$sid == name & data$channel == chan & data$numStims == numStimTrial & data$phaseClass == typePhase,]$magnitude-baseMedian
        data[data$sid == name & data$channel == chan & data$numStims == numStimTrial & data$phaseClass == typePhase,]$absDiff = absDiff
        if (length(absDiff)){
        data[data$sid == name & data$channel == chan & data$numStims == numStimTrial & data$phaseClass == typePhase,]$baseMedian = baseMedian
        }
      }
    }
  }
}

sapply(data,class)
#summaryData = ddply(data[data$numStims != "Base",] , .(sid,phaseClass,numStims,channel), function(x) mean(x[,"percentDiff"]))

summaryData = ddply(data, .(sid,phaseClass,numStims,channel,betaLabels), summarize, magnitude = median(magnitude))

summaryDataForMixed = ddply(data, .(sid,phaseClass,numStims,channel,betaLabels), summarize, magnitude = median(magnitude))


summaryData = ddply(data[data$numStims != "Base",] , .(sid,phaseClass,numStims,channel,betaLabels), summarize, percentDiff = median(percentDiff))

dataNoBaseline = data[data$numStims != "Base",]
dataSubjOnly <- subset(data,data$sid=='0b5a2e')

summaryDataNoPhase = ddply(data, .(sid,numStims,channel,betaLabels), summarize, magnitude = median(magnitude))
summaryDataNoPhase = ddply(data[data$numStims != "Base",] , .(sid,numStims,channel,betaLabels), summarize, percentDiff = median(percentDiff))

summaryDataHighStimsOnly = ddply(data[data$numStims == "[5,inf)",] , .(sid,phaseClass,numStims,channel,orderedPhase45), summarize, percentDiff = median(percentDiff))


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
  #ggsave(paste0("betaStim_dose.svg"), units="in", width=figWidth, height=figHeight,dpi=600)
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
  #ggsave(paste0("betaStim_dose_no_dots.svg"), units="in", width=figWidth, height=figHeight,dpi=600)
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
  #ggsave(paste0("betaStim_dose_binned45.svg"), units="in", width=figWidth, height=figHeight,dpi=600)
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

# ========================================================================
# Summary-level models: one median per (subject x channel x phaseClass x
# numStims) cell. Eliminates pseudoreplication. Three complementary
# specifications that converge on the same result: dose-dependent
# enhancement selective to hyperpolarizing (270) channels.
#
# Data:
#   - summaryNB: 120 obs (no baseline), outcome = absDiff (median)
#   - summaryAll: 151 obs (with baseline), outcome = magnitude (median)
#   - summaryNB_ancova: 120 obs (no baseline), outcome = magnitude (median)
#     with baseline magnitude as covariate
#   - 6 subjects, 31 channels
#   - Channel IDs unique per subject, so (1|channel) nests within subject
# ========================================================================

# --- Prepare summary datasets ---
summaryNB <- ddply(dataNoBaseline, .(sid,phaseClass,numStims,channel),
                   summarize, absDiff = median(absDiff), magnitude = median(magnitude))

summaryAll <- ddply(data, .(sid,phaseClass,numStims,channel),
                    summarize, magnitude = median(magnitude))
summaryAll$numStims <- relevel(summaryAll$numStims, ref = "Base")
summaryAll$doseNum <- as.numeric(factor(summaryAll$numStims,
                      levels = c("Base","[1,2]","[3,4]","[5,inf)"))) - 1

# baseline covariate: one median per channel, pooled across phaseClass
basePerChan <- ddply(data[data$numStims == "Base",], .(sid, channel),
                     summarize, baselineMag = median(magnitude))
summaryNB_ancova <- merge(summaryNB, basePerChan, by = c("sid", "channel"))
summaryNB_ancova$baselineMag_c <- summaryNB_ancova$baselineMag - mean(summaryNB_ancova$baselineMag)
summaryNB_ancova$doseNum <- as.numeric(factor(summaryNB_ancova$numStims,
                            levels = c("[1,2]","[3,4]","[5,inf)"))) - 1

cat(sprintf("Summary data: %d obs (no baseline), %d obs (with baseline), %d subjects, %d channels\n",
    nrow(summaryNB), nrow(summaryAll), length(unique(summaryAll$sid)),
    length(unique(summaryAll$channel))))

# ------------------------------------------------------------------------
# Model 3a: absDiff, random intercepts only (original approach)
# ------------------------------------------------------------------------
# Outcome: median absolute difference from baseline (pre-subtracted per
# channel using pooled baseMedian). No baseline category in model.
# Random intercepts only — categorical dose slopes per subject cause
# singularity with 6 subjects (correlations hit 1.0).

fit.absDiff = lmerTest::lmer(
  absDiff ~ numStims * phaseClass +
  (1|sid) + (1|channel),
  data = summaryNB,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=20000)))

cat("\n=== Model 3a: absDiff, random intercepts only ===\n")
summary(fit.absDiff)
anova(fit.absDiff)
VarCorr(fit.absDiff)
report(fit.absDiff)
report(anova(fit.absDiff))
model_performance(fit.absDiff)
eta_squared(fit.absDiff, ci = 0.95)

emm_3a_dose <- emmeans(fit.absDiff, ~ numStims | phaseClass)
emm_3a_dose
pairs(emm_3a_dose)
confint(pairs(emm_3a_dose))

emm_3a_phase <- emmeans(fit.absDiff, ~ phaseClass | numStims)
pairs(emm_3a_phase)
confint(pairs(emm_3a_phase))

tab_model(fit.absDiff)

# emmip plot: absDiff emmeans (already baseline-subtracted)
emm_3a_df <- as.data.frame(emm_3a_dose)
emm_3a_df$phaseClass <- factor(emm_3a_df$phaseClass, levels = c("90", "270"))

pd <- position_dodge(0.15)
p_3a <- ggplot(emm_3a_df, aes(x = numStims, y = emmean, color = phaseClass, group = phaseClass)) +
  theme_light(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(position = pd, linewidth = 0.8) + geom_point(position = pd, size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = pd, linewidth = 0.8) +
  labs(x = "Number of Conditioning Stimuli",
       y = expression(paste("Absolute Difference from Baseline (", mu, "V)")),
       color = "Phase Class",
       title = "Model 3a: absDiff Emmeans (baseline pre-subtracted)") +
  scale_color_hue(labels = c("Depolarizing (90)", "Hyperpolarizing (270)"))
p_3a

figHeight = 4.5
figWidth = 6.5
if(savePlot){
  ggsave(here("output_plots","betaStim_model3a_emmip.png"), plot = p_3a,
         units = "in", width = figWidth, height = figHeight, dpi = 600)
  ggsave(here("output_plots","betaStim_model3a_emmip.eps"), plot = p_3a,
         units = "in", width = figWidth, height = figHeight, dpi = 600, device = cairo_ps)
}

# forest plot: Model 3a pairwise dose contrasts within each phase (Tukey CIs)
dose_3a_ci <- as.data.frame(confint(pairs(emm_3a_dose)))
dose_3a_ci$phaseClass <- factor(dose_3a_ci$phaseClass, levels = c("90", "270"))
dose_3a_ci$sig <- ifelse(dose_3a_ci$lower.CL > 0 | dose_3a_ci$upper.CL < 0, "sig", "ns")

pd_forest <- position_dodge(0.4)
p_3a_forest <- ggplot(dose_3a_ci, aes(x = contrast, y = estimate, color = phaseClass, shape = sig)) +
  theme_light(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(position = pd_forest, size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15, position = pd_forest, linewidth = 0.8) +
  scale_color_hue(labels = c("Depolarizing (90)", "Hyperpolarizing (270)")) +
  scale_shape_manual(values = c("ns" = 1, "sig" = 16), guide = "none") +
  labs(x = "Dose Contrast",
       y = expression(paste("Difference (", mu, "V) with 95% Tukey CI")),
       color = "Phase Class",
       title = "Model 3a: Pairwise Dose Contrasts Within Each Phase") +
  coord_flip()
p_3a_forest

figHeight = 4
figWidth = 7
if(savePlot){
  ggsave(here("output_plots","betaStim_model3a_forest.png"), plot = p_3a_forest,
         units = "in", width = figWidth, height = figHeight, dpi = 600)
  ggsave(here("output_plots","betaStim_model3a_forest.eps"), plot = p_3a_forest,
         units = "in", width = figWidth, height = figHeight, dpi = 600, device = cairo_ps)
}

# ------------------------------------------------------------------------
# Model 3b: magnitude with baseline category + random linear dose slope
# ------------------------------------------------------------------------
# Outcome: raw magnitude (median). Baseline included as a level of
# numStims (4 levels: Base, [1,2], [3,4], [5,inf)).
# Random linear dose slope (doseNum = 0,1,2,3) per subject — allows
# subjects to differ in overall dose sensitivity without the singularity
# caused by categorical dose slopes (4x4 covariance from 6 subjects).
# Fixed effects keep numStims categorical for non-linear dose patterns.

fit.modelD = lmerTest::lmer(
  magnitude ~ numStims * phaseClass +
  (1 + doseNum|sid) + (1|channel),
  data = summaryAll,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=50000)))

cat("\n=== Model 3b: magnitude + baseline category + random dose slope ===\n")
summary(fit.modelD)
anova(fit.modelD)
VarCorr(fit.modelD)
report(fit.modelD)
report(anova(fit.modelD))
model_performance(fit.modelD)
eta_squared(fit.modelD, ci = 0.95)

emm_3b_dose <- emmeans(fit.modelD, ~ numStims | phaseClass)
emm_3b_dose
contr_3b_base <- contrast(emm_3b_dose, method = "trt.vs.ctrl", ref = "Base")
contr_3b_base
confint(contr_3b_base)

emm_3b_phase <- emmeans(fit.modelD, ~ phaseClass | numStims)
pairs(emm_3b_phase)
confint(pairs(emm_3b_phase))

tab_model(fit.modelD)

# emmip plot: Model 3b contrasts vs baseline
contr_3b_ci <- as.data.frame(confint(contr_3b_base))
contr_3b_ci$dose <- sub(" - .*", "", contr_3b_ci$contrast)
contr_3b_ci$phaseClass <- as.character(contr_3b_ci$phaseClass)
ref_3b <- data.frame(phaseClass = as.character(unique(contr_3b_ci$phaseClass)),
  dose = "Base", estimate = 0, SE = 0, lower.CL = 0, upper.CL = 0)
plot_3b <- bind_rows(
  ref_3b %>% select(phaseClass, dose, estimate, SE, lower.CL, upper.CL),
  contr_3b_ci %>% select(phaseClass, dose, estimate, SE, lower.CL, upper.CL))
plot_3b$dose <- factor(plot_3b$dose, levels = c("Base", "[1,2]", "[3,4]", "[5,inf)"))
plot_3b$phaseClass <- factor(plot_3b$phaseClass, levels = c("90", "270"))

p_3b <- ggplot(plot_3b, aes(x = dose, y = estimate, color = phaseClass, group = phaseClass)) +
  theme_light(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(position = pd, linewidth = 0.8) + geom_point(position = pd, size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = pd, linewidth = 0.8) +
  labs(x = "Number of Conditioning Stimuli",
       y = expression(paste(Delta, " Magnitude from Baseline (", mu, "V)")),
       color = "Phase Class",
       title = "Model 3b: Contrasts vs Baseline (random dose slope)") +
  scale_color_hue(labels = c("Depolarizing (90)", "Hyperpolarizing (270)"))
p_3b

figHeight = 4.5
figWidth = 6.5
if(savePlot){
  ggsave(here("output_plots","betaStim_model3b_emmip.png"), plot = p_3b,
         units = "in", width = figWidth, height = figHeight, dpi = 600)
  ggsave(here("output_plots","betaStim_model3b_emmip.eps"), plot = p_3b,
         units = "in", width = figWidth, height = figHeight, dpi = 600, device = cairo_ps)
}

# forest plot: Model 3b dose vs baseline within each phase (Dunnett CIs)
contr_3b_forest <- as.data.frame(confint(contr_3b_base))
contr_3b_forest$phaseClass <- factor(contr_3b_forest$phaseClass, levels = c("90", "270"))
contr_3b_forest$sig <- ifelse(contr_3b_forest$lower.CL > 0 | contr_3b_forest$upper.CL < 0, "sig", "ns")

p_3b_forest <- ggplot(contr_3b_forest, aes(x = contrast, y = estimate, color = phaseClass, shape = sig)) +
  theme_light(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(position = pd_forest, size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15, position = pd_forest, linewidth = 0.8) +
  scale_color_hue(labels = c("Depolarizing (90)", "Hyperpolarizing (270)")) +
  scale_shape_manual(values = c("ns" = 1, "sig" = 16), guide = "none") +
  labs(x = "Dose vs Baseline",
       y = expression(paste("Difference (", mu, "V) with 95% Dunnett CI")),
       color = "Phase Class",
       title = "Model 3b: Dose vs Baseline Contrasts Within Each Phase") +
  coord_flip()
p_3b_forest

figHeight = 4
figWidth = 7
if(savePlot){
  ggsave(here("output_plots","betaStim_model3b_forest.png"), plot = p_3b_forest,
         units = "in", width = figWidth, height = figHeight, dpi = 600)
  ggsave(here("output_plots","betaStim_model3b_forest.eps"), plot = p_3b_forest,
         units = "in", width = figWidth, height = figHeight, dpi = 600, device = cairo_ps)
}

# ------------------------------------------------------------------------
# Model 3c: ANCOVA — baseline magnitude as covariate + random dose slope
# ------------------------------------------------------------------------
# Outcome: raw magnitude (median) on non-baseline trials only.
# Baseline magnitude (median per channel, pooled across phaseClass) is a
# fixed-effect covariate, centered at the grand mean (~310 uV).
# This absorbs between-channel magnitude differences more effectively
# than random intercepts alone (baseline coeff ~ 1.03, near 1:1).
# Channel random intercept captures residual channel variance after
# baseline adjustment. Random linear dose slope per subject as in 3b.

cat(sprintf("\n=== Model 3c: ANCOVA (baseline covariate, centered at %.1f uV) ===\n",
    mean(summaryNB_ancova$baselineMag)))

fit.ancova = lmerTest::lmer(
  magnitude ~ numStims * phaseClass + baselineMag_c +
  (1 + doseNum|sid) + (1|channel),
  data = summaryNB_ancova,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=50000)))

summary(fit.ancova)
anova(fit.ancova)
VarCorr(fit.ancova)
report(fit.ancova)
report(anova(fit.ancova))
model_performance(fit.ancova)
eta_squared(fit.ancova, ci = 0.95)

emm_3c_dose <- emmeans(fit.ancova, ~ numStims | phaseClass)
emm_3c_dose
pairs(emm_3c_dose)
confint(pairs(emm_3c_dose))

emm_3c_phase <- emmeans(fit.ancova, ~ phaseClass | numStims)
pairs(emm_3c_phase)
confint(pairs(emm_3c_phase))

tab_model(fit.ancova)

# emmip plot: ANCOVA contrasts vs [1,2]
contr_3c_ref <- contrast(emm_3c_dose, method = "trt.vs.ctrl", ref = 1)
contr_3c_ci <- as.data.frame(confint(contr_3c_ref))
contr_3c_ci$dose <- sub(" - .*", "", contr_3c_ci$contrast)
contr_3c_ci$phaseClass <- as.character(contr_3c_ci$phaseClass)
ref_3c <- data.frame(phaseClass = as.character(unique(contr_3c_ci$phaseClass)),
  dose = "[1,2]", estimate = 0, SE = 0, lower.CL = 0, upper.CL = 0)
plot_3c <- bind_rows(
  ref_3c %>% select(phaseClass, dose, estimate, SE, lower.CL, upper.CL),
  contr_3c_ci %>% select(phaseClass, dose, estimate, SE, lower.CL, upper.CL))
plot_3c$dose <- factor(plot_3c$dose, levels = c("[1,2]", "[3,4]", "[5,inf)"))
plot_3c$phaseClass <- factor(plot_3c$phaseClass, levels = c("90", "270"))

p_3c <- ggplot(plot_3c, aes(x = dose, y = estimate, color = phaseClass, group = phaseClass)) +
  theme_light(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(position = pd, linewidth = 0.8) + geom_point(position = pd, size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = pd, linewidth = 0.8) +
  labs(x = "Number of Conditioning Stimuli",
       y = expression(paste(Delta, " Magnitude from [1,2] (", mu, "V)")),
       color = "Phase Class",
       title = "Model 3c: ANCOVA Dose-Response (baseline covariate)") +
  scale_color_hue(labels = c("Depolarizing (90)", "Hyperpolarizing (270)"))
p_3c

figHeight = 4.5
figWidth = 6.5
if(savePlot){
  ggsave(here("output_plots","betaStim_model3c_emmip.png"), plot = p_3c,
         units = "in", width = figWidth, height = figHeight, dpi = 600)
  ggsave(here("output_plots","betaStim_model3c_emmip.eps"), plot = p_3c,
         units = "in", width = figWidth, height = figHeight, dpi = 600, device = cairo_ps)
}

# forest plot: Model 3c pairwise dose contrasts within each phase (Tukey CIs)
dose_3c_ci <- as.data.frame(confint(pairs(emm_3c_dose)))
dose_3c_ci$phaseClass <- factor(dose_3c_ci$phaseClass, levels = c("90", "270"))
dose_3c_ci$sig <- ifelse(dose_3c_ci$lower.CL > 0 | dose_3c_ci$upper.CL < 0, "sig", "ns")

p_3c_forest <- ggplot(dose_3c_ci, aes(x = contrast, y = estimate, color = phaseClass, shape = sig)) +
  theme_light(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(position = pd_forest, size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15, position = pd_forest, linewidth = 0.8) +
  scale_color_hue(labels = c("Depolarizing (90)", "Hyperpolarizing (270)")) +
  scale_shape_manual(values = c("ns" = 1, "sig" = 16), guide = "none") +
  labs(x = "Dose Contrast",
       y = expression(paste("Difference (", mu, "V) with 95% Tukey CI")),
       color = "Phase Class",
       title = "Model 3c: ANCOVA Pairwise Dose Contrasts Within Each Phase") +
  coord_flip()
p_3c_forest

figHeight = 4
figWidth = 7
if(savePlot){
  ggsave(here("output_plots","betaStim_model3c_forest.png"), plot = p_3c_forest,
         units = "in", width = figWidth, height = figHeight, dpi = 600)
  ggsave(here("output_plots","betaStim_model3c_forest.eps"), plot = p_3c_forest,
         units = "in", width = figWidth, height = figHeight, dpi = 600, device = cairo_ps)
}

# ========================================================================
# Model 3d: Ordinal dose (polynomial contrasts) + linear random slope
# ========================================================================
# Convert numStims to ordered factor. R auto-generates polynomial contrasts:
#   .L (linear): tests whether magnitude increases linearly with dose
#   .Q (quadratic): tests whether the dose-response curve accelerates/decelerates
# Random slope uses only the linear polynomial component extracted from
# contr.poly(3)[,1] = c(-0.707, 0, 0.707). This is a linear rescaling of
# doseNum = (0,1,2), specifically dose_linpoly = (doseNum - 1) / sqrt(2).
# Same random effect structure as Model 3c — same convergence properties.
# Fixed effects are reparameterized (polynomial vs treatment contrasts)
# but the model has identical df and likelihood as 3c.

summaryNB_ancova$numStims_ord <- ordered(
  summaryNB_ancova$numStims,
  levels = c("[1,2]", "[3,4]", "[5,inf)"))

# extract linear polynomial contrast as numeric column for random slope
poly_lin <- contr.poly(3)[, 1]
summaryNB_ancova$dose_linpoly <- poly_lin[as.numeric(summaryNB_ancova$numStims_ord)]

cat("\n=== Model 3d: Ordinal dose (polynomial contrasts) ===\n")
cat("Linear polynomial values:", poly_lin, "\n")

# --- 3d-i: afex::mixed with expand_re + per_parameter ---
# || gives uncorrelated .L and .Q random slopes. per_parameter tests
# .L and .Q fixed effects separately. If .Q random variance ~ 0,
# confirms linear-only slope is sufficient.
fit.ordinal_afex = afex::mixed(
  magnitude ~ numStims_ord * phaseClass + baselineMag_c +
  (numStims_ord || sid) + (1 | channel),
  data = summaryNB_ancova,
  expand_re = TRUE,
  per_parameter = "numStims_ord",
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

cat("\n--- 3d-i: afex (expand_re, .L+.Q uncorrelated random slopes) ---\n")
cat("Singular:", isSingular(fit.ordinal_afex$full_model), "\n")
summary(fit.ordinal_afex$full_model)
cat("\nafex ANOVA (per-parameter for ordered factor):\n")
print(anova(fit.ordinal_afex))
VarCorr(fit.ordinal_afex$full_model)

# --- 3d-ii: lmerTest with manual dose_linpoly (.L only random slope) ---
# Restricts random slope to linear component only. Compare to 3d-i
# to verify .Q random slope is negligible.
fit.ordinal_lmer = lmerTest::lmer(
  magnitude ~ numStims_ord * phaseClass + baselineMag_c +
  (1 + dose_linpoly | sid) + (1 | channel),
  data = summaryNB_ancova,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

cat("\n--- 3d-ii: lmerTest (.L only random slope) ---\n")
cat("Singular:", isSingular(fit.ordinal_lmer), "\n")
summary(fit.ordinal_lmer)
anova(fit.ordinal_lmer)

cat(sprintf("\n3d-i (afex, .L+.Q):    AIC=%.1f  BIC=%.1f  singular=%s\n",
    AIC(fit.ordinal_afex$full_model), BIC(fit.ordinal_afex$full_model),
    isSingular(fit.ordinal_afex$full_model)))
cat(sprintf("3d-ii (lmer, .L only): AIC=%.1f  BIC=%.1f  singular=%s\n",
    AIC(fit.ordinal_lmer), BIC(fit.ordinal_lmer), isSingular(fit.ordinal_lmer)))

model_performance(fit.ordinal_lmer)
eta_squared(fit.ordinal_lmer, ci = 0.95)

# emmeans: dose within phase
emm_3d_dose <- emmeans(fit.ordinal_lmer, ~ numStims_ord | phaseClass)
emm_3d_dose
pairs(emm_3d_dose)
confint(pairs(emm_3d_dose))

# phase contrast at each dose
emm_3d_phase <- emmeans(fit.ordinal_lmer, ~ phaseClass | numStims_ord)
pairs(emm_3d_phase)
confint(pairs(emm_3d_phase))

# emmip plot: ordinal model
emm_3d_df <- as.data.frame(emm_3d_dose)
emm_3d_df$phaseClass <- factor(emm_3d_df$phaseClass, levels = c("90", "270"))

p_3d <- ggplot(emm_3d_df, aes(x = numStims_ord, y = emmean, color = phaseClass, group = phaseClass)) +
  theme_light(base_size = 14) +
  geom_line(position = pd, linewidth = 0.8) + geom_point(position = pd, size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = pd, linewidth = 0.8) +
  labs(x = "Number of Conditioning Stimuli",
       y = expression(paste("Predicted Magnitude (", mu, "V)")),
       color = "Phase Class",
       title = "Model 3d: Ordinal Dose (polynomial contrasts)") +
  scale_color_hue(labels = c("Depolarizing (90)", "Hyperpolarizing (270)"))
p_3d

figHeight = 4.5
figWidth = 6.5
if(savePlot){
  ggsave(here("output_plots","betaStim_model3d_ordinal_emmip.png"), plot = p_3d,
         units = "in", width = figWidth, height = figHeight, dpi = 600)
  ggsave(here("output_plots","betaStim_model3d_ordinal_emmip.eps"), plot = p_3d,
         units = "in", width = figWidth, height = figHeight, dpi = 600, device = cairo_ps)
}

# ========================================================================
# Model 3e: Fully numeric dose (fixed + random)
# ========================================================================
# Both fixed and random effects use numeric doseNum (0, 1, 2).
# Most parsimonious — forces linear dose-response. Cannot detect
# non-linear (quadratic) patterns. Nested within 3c/3d.

cat("\n=== Model 3e: Fully numeric dose ===\n")

fit.numeric = lmerTest::lmer(
  magnitude ~ doseNum * phaseClass + baselineMag_c +
  (1 + doseNum | sid) + (1 | channel),
  data = summaryNB_ancova,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

cat("Singular:", isSingular(fit.numeric), "\n")
summary(fit.numeric)
anova(fit.numeric)
VarCorr(fit.numeric)
model_performance(fit.numeric)
eta_squared(fit.numeric, ci = 0.95)

# emmeans: evaluate at doseNum = 0, 1, 2
emm_3e_dose <- emmeans(fit.numeric, ~ doseNum | phaseClass,
                        at = list(doseNum = c(0, 1, 2)))
emm_3e_dose
pairs(emm_3e_dose)
confint(pairs(emm_3e_dose))

# phase contrast at each dose
emm_3e_phase <- emmeans(fit.numeric, ~ phaseClass | doseNum,
                         at = list(doseNum = c(0, 1, 2)))
pairs(emm_3e_phase)
confint(pairs(emm_3e_phase))

# emmip plot: numeric model
emm_3e_df <- as.data.frame(emm_3e_dose)
emm_3e_df$phaseClass <- factor(emm_3e_df$phaseClass, levels = c("90", "270"))
emm_3e_df$dose_label <- factor(
  c("[1,2]", "[3,4]", "[5,inf)")[emm_3e_df$doseNum + 1],
  levels = c("[1,2]", "[3,4]", "[5,inf)"))

p_3e <- ggplot(emm_3e_df, aes(x = dose_label, y = emmean, color = phaseClass, group = phaseClass)) +
  theme_light(base_size = 14) +
  geom_line(position = pd, linewidth = 0.8) + geom_point(position = pd, size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = pd, linewidth = 0.8) +
  labs(x = "Number of Conditioning Stimuli",
       y = expression(paste("Predicted Magnitude (", mu, "V)")),
       color = "Phase Class",
       title = "Model 3e: Numeric Dose (linear only)") +
  scale_color_hue(labels = c("Depolarizing (90)", "Hyperpolarizing (270)"))
p_3e

figHeight = 4.5
figWidth = 6.5
if(savePlot){
  ggsave(here("output_plots","betaStim_model3e_numeric_emmip.png"), plot = p_3e,
         units = "in", width = figWidth, height = figHeight, dpi = 600)
  ggsave(here("output_plots","betaStim_model3e_numeric_emmip.eps"), plot = p_3e,
         units = "in", width = figWidth, height = figHeight, dpi = 600, device = cairo_ps)
}

# ========================================================================
# Likelihood ratio test: numeric (3e) vs ordinal (3d)
# ========================================================================
# Model 3e is nested within 3d (drops the quadratic component).
# Significant LRT means the quadratic term improves fit.
cat("\n=== LRT: Model 3e (numeric) vs Model 3d (ordinal) ===\n")
print(anova(fit.numeric, fit.ordinal_lmer))

# ------------------------------------------------------------------------
# Model comparison summary
# ------------------------------------------------------------------------
cat("\n=== Summary-level model comparison ===\n")
cat("Model 3a (absDiff, intercepts):  AIC=", round(AIC(fit.absDiff),1), " BIC=", round(BIC(fit.absDiff),1), " singular=", isSingular(fit.absDiff), "\n")
cat("Model 3b (mag+baseline, slope):  AIC=", round(AIC(fit.modelD),1), " BIC=", round(BIC(fit.modelD),1), " singular=", isSingular(fit.modelD), "\n")
cat("Model 3c (ANCOVA, slope):        AIC=", round(AIC(fit.ancova),1), " BIC=", round(BIC(fit.ancova),1), " singular=", isSingular(fit.ancova), "\n")
cat("Model 3d (ordinal, slope):       AIC=", round(AIC(fit.ordinal_lmer),1), " BIC=", round(BIC(fit.ordinal_lmer),1), " singular=", isSingular(fit.ordinal_lmer), "\n")
cat("Model 3e (numeric, slope):       AIC=", round(AIC(fit.numeric),1), " BIC=", round(BIC(fit.numeric),1), " singular=", isSingular(fit.numeric), "\n")

# residual diagnostics for primary model (3c)
figHeight = 4
figWidth = 8
if(savePlot){
  png(here("output_plots","betaStim_residuals_summary_level.png"),width=figWidth,height=figHeight,units="in",res=600)
  plot(fit.ancova)
  dev.off()

  png(here("output_plots","betaStim_qq_summary_level.png"),width=figWidth,height=figHeight,units="in",res=600)
  qqnorm(resid(fit.ancova)); qqline(resid(fit.ancova))
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

# Model 3a (absDiff): dose-response within each phase
cat("--- Model 3a (absDiff): dose within phase ---\n")
eff_size(emm_3a_dose, sigma = sigma_trial, edf = edf_trial)

# Model 3a (absDiff): phase contrast at each dose
cat("--- Model 3a (absDiff): phase at each dose ---\n")
eff_size(emm_3a_phase, sigma = sigma_trial, edf = edf_trial)

# Model 3c (ANCOVA): dose-response within each phase
cat("--- Model 3c (ANCOVA): dose within phase ---\n")
eff_size(emm_3c_dose, sigma = sigma_trial, edf = edf_trial)

# Model 3c (ANCOVA): phase contrast at each dose
cat("--- Model 3c (ANCOVA): phase at each dose ---\n")
eff_size(emm_3c_phase, sigma = sigma_trial, edf = edf_trial)

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

# ========================================================================
# Residual diagnostics: skewness, kurtosis, normality tests
# ========================================================================
# Evaluate residual distributions for all summary-level models (3a-3e).
# Guidelines: |skewness| < 1 and |excess kurtosis| < 2 are acceptable
# for LME with these sample sizes. Shapiro-Wilk is conservative at n=120+.

residual_diagnostics <- function(fit, model_name) {
  r <- resid(fit)
  n <- length(r)
  sw <- shapiro.test(r)
  skew <- (sum((r - mean(r))^3) / n) / (sum((r - mean(r))^2) / n)^1.5
  kurt <- (sum((r - mean(r))^4) / n) / (sum((r - mean(r))^2) / n)^2 - 3

  cat(sprintf("\n--- %s ---\n", model_name))
  cat("N residuals:", n, "\n")
  cat("Shapiro-Wilk W =", round(sw$statistic, 4),
      ", p =", format(sw$p.value, digits = 3), "\n")
  cat("Skewness:", round(skew, 3), "\n")
  cat("Excess kurtosis:", round(kurt, 3), "\n")

  data.frame(
    Model = model_name, N = n,
    Shapiro_W = round(sw$statistic, 4),
    Shapiro_p = sw$p.value,
    Skewness = round(skew, 3),
    Kurtosis = round(kurt, 3),
    stringsAsFactors = FALSE
  )
}

cat("\n=== RESIDUAL DIAGNOSTICS (summary-level models) ===\n")
diag_list <- list(
  residual_diagnostics(fit.absDiff, "3a: absDiff"),
  residual_diagnostics(fit.modelD, "3b: Baseline category"),
  residual_diagnostics(fit.ancova, "3c: ANCOVA"),
  residual_diagnostics(fit.ordinal_lmer, "3d: Ordinal"),
  residual_diagnostics(fit.numeric, "3e: Numeric")
)
diag_df <- do.call(rbind, diag_list)
print(diag_df)

# --- ggplot QQ plots and residuals-vs-fitted for all summary models ---
if (savePlot) {
  summary_models <- list(
    list(fit = fit.absDiff, name = "3a_absDiff"),
    list(fit = fit.modelD, name = "3b_baseline"),
    list(fit = fit.ancova, name = "3c_ANCOVA"),
    list(fit = fit.ordinal_lmer, name = "3d_ordinal"),
    list(fit = fit.numeric, name = "3e_numeric")
  )

  figHeight <- 5
  figWidth <- 6

  for (m in summary_models) {
    r <- resid(m$fit)
    n <- length(r)
    skew_val <- (sum((r - mean(r))^3) / n) / (sum((r - mean(r))^2) / n)^1.5
    kurt_val <- (sum((r - mean(r))^4) / n) / (sum((r - mean(r))^2) / n)^2 - 3

    # QQ plot
    qq_df <- data.frame(residual = r)
    p_qq <- ggplot(qq_df, aes(sample = residual)) +
      stat_qq(alpha = 0.6) + stat_qq_line(color = "red", linewidth = 0.8) +
      theme_light(base_size = 14) +
      labs(title = paste("QQ Plot: Model", m$name),
           subtitle = sprintf("Skew = %.2f, Kurtosis = %.2f", skew_val, kurt_val),
           x = "Theoretical Quantiles", y = "Sample Quantiles")

    ggsave(here("output_plots", paste0("betaStim_qq_", m$name, ".png")),
           plot = p_qq, units = "in", width = figWidth, height = figHeight, dpi = 600)

    # Residuals vs fitted
    rf_df <- data.frame(fitted = fitted(m$fit), residual = r)
    p_rf <- ggplot(rf_df, aes(x = fitted, y = residual)) +
      geom_point(alpha = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_smooth(method = "loess", se = TRUE, color = "blue", linewidth = 0.8) +
      theme_light(base_size = 14) +
      labs(title = paste("Residuals vs Fitted: Model", m$name),
           x = "Fitted Values", y = "Residuals")

    ggsave(here("output_plots", paste0("betaStim_resid_vs_fitted_", m$name, ".png")),
           plot = p_rf, units = "in", width = figWidth, height = figHeight, dpi = 600)
  }
}

# ========================================================================
# Export manuscript-ready .docx tables (officer + flextable)
# ========================================================================
if (requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("flextable", quietly = TRUE)) {
  library(officer)
  library(flextable)

  fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))

  outputDir <- here("output_plots")
  doc <- read_docx()

  # ------------------------------------------------------------------
  # Table 1: Residual Diagnostics
  # ------------------------------------------------------------------
  diag_out <- diag_df
  diag_out$Shapiro_p <- sapply(diag_out$Shapiro_p, fmt_p)
  doc <- body_add_par(doc, "Table: Residual Diagnostics", style = "heading 2")
  ft <- flextable(diag_out) |> autofit() |>
    set_caption("Residual normality diagnostics for summary-level models (3a-3e). Skewness and excess kurtosis computed from standardized residuals; Shapiro-Wilk tests normality.")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # ------------------------------------------------------------------
  # Table 2: Model Comparison Summary
  # ------------------------------------------------------------------
  comparison_df <- data.frame(
    Model = c("3a: absDiff (intercepts)",
              "3b: Magnitude + baseline",
              "3c: ANCOVA (primary)",
              "3d: Ordinal dose",
              "3e: Numeric dose"),
    N = c(nrow(summaryNB), nrow(summaryAll), nrow(summaryNB_ancova),
          nrow(summaryNB_ancova), nrow(summaryNB_ancova)),
    AIC = round(c(AIC(fit.absDiff), AIC(fit.modelD), AIC(fit.ancova),
                   AIC(fit.ordinal_lmer), AIC(fit.numeric)), 1),
    BIC = round(c(BIC(fit.absDiff), BIC(fit.modelD), BIC(fit.ancova),
                   BIC(fit.ordinal_lmer), BIC(fit.numeric)), 1),
    Singular = c(isSingular(fit.absDiff), isSingular(fit.modelD),
                 isSingular(fit.ancova), isSingular(fit.ordinal_lmer),
                 isSingular(fit.numeric)),
    stringsAsFactors = FALSE
  )
  doc <- body_add_par(doc, "Table: Model Comparison", style = "heading 2")
  ft <- flextable(comparison_df) |> autofit() |>
    set_caption("Summary-level model comparison. All models use median per (subject x channel x phaseClass x dose) cell.")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # ------------------------------------------------------------------
  # Helper: add ANOVA + Fixed Effects + Random Effects for a model
  # ------------------------------------------------------------------
  add_model_tables <- function(doc, fit, model_label) {
    # Random Effects
    vc <- VarCorr(fit)
    ngrps <- summary(fit)$ngrps
    re_rows <- list()
    for (nm in names(vc)) {
      v <- vc[[nm]]
      ng <- as.character(ngrps[nm])
      if (ncol(v) == 1) {
        re_rows[[length(re_rows) + 1]] <- data.frame(
          Component = nm, Term = "(Intercept)",
          Variance = round(v[1, 1], 4), SD = round(sqrt(v[1, 1]), 4),
          Corr = "", Groups = ng, stringsAsFactors = FALSE)
      } else {
        corr_val <- sprintf("%.2f", attr(v, "correlation")[2, 1])
        re_rows[[length(re_rows) + 1]] <- data.frame(
          Component = c(nm, ""), Term = rownames(v),
          Variance = round(diag(v), 4), SD = round(sqrt(diag(v)), 4),
          Corr = c("", corr_val), Groups = c(ng, ""),
          stringsAsFactors = FALSE)
      }
    }
    re_df <- do.call(rbind, re_rows)
    re_df <- rbind(re_df, data.frame(
      Component = "Residual", Term = "",
      Variance = round(sigma(fit)^2, 4), SD = round(sigma(fit), 4),
      Corr = "", Groups = ""))

    doc <- body_add_par(doc, paste("Random Effects:", model_label), style = "heading 2")
    ft <- flextable(re_df) |> autofit() |>
      set_caption(paste("Random effects for", model_label))
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, sprintf("Singular: %s. AIC = %.1f", isSingular(fit), AIC(fit)))
    doc <- body_add_par(doc, "")

    # Type III ANOVA
    aov_tbl <- as.data.frame(anova(fit, type = 3))
    aov_tbl$Effect <- rownames(aov_tbl)
    aov_tbl <- aov_tbl[, c("Effect", "Sum Sq", "Mean Sq", "NumDF", "DenDF", "F value", "Pr(>F)")]
    aov_tbl$`Sum Sq` <- round(aov_tbl$`Sum Sq`, 3)
    aov_tbl$`Mean Sq` <- round(aov_tbl$`Mean Sq`, 3)
    aov_tbl$DenDF <- round(aov_tbl$DenDF, 1)
    aov_tbl$`F value` <- round(aov_tbl$`F value`, 2)
    aov_tbl$p <- sapply(aov_tbl$`Pr(>F)`, fmt_p)
    aov_tbl$`Pr(>F)` <- NULL

    doc <- body_add_par(doc, paste("Type III ANOVA:", model_label), style = "heading 2")
    ft <- flextable(aov_tbl) |> autofit() |>
      set_caption(paste("Type III ANOVA (Satterthwaite df) for", model_label))
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")

    # Fixed Effects
    fe_tbl <- as.data.frame(summary(fit)$coefficients)
    fe_tbl$Predictor <- rownames(fe_tbl)
    fe_tbl <- fe_tbl[, c("Predictor", "Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    fe_tbl$Estimate <- round(fe_tbl$Estimate, 3)
    fe_tbl$`Std. Error` <- round(fe_tbl$`Std. Error`, 3)
    fe_tbl$df <- round(fe_tbl$df, 1)
    fe_tbl$`t value` <- round(fe_tbl$`t value`, 2)
    fe_tbl$p <- sapply(fe_tbl$`Pr(>|t|)`, fmt_p)
    fe_tbl$`Pr(>|t|)` <- NULL

    doc <- body_add_par(doc, paste("Fixed Effects:", model_label), style = "heading 2")
    ft <- flextable(fe_tbl) |> autofit() |>
      set_caption(paste("Fixed effect coefficients for", model_label))
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")

    doc
  }

  # ------------------------------------------------------------------
  # Tables for primary models: 3a, 3c, 3e
  # ------------------------------------------------------------------
  doc <- add_model_tables(doc, fit.absDiff, "Model 3a (absDiff, intercepts only)")
  doc <- add_model_tables(doc, fit.ancova, "Model 3c (ANCOVA, primary)")
  doc <- add_model_tables(doc, fit.numeric, "Model 3e (numeric dose)")

  # ------------------------------------------------------------------
  # EMM Dose Contrasts: Model 3a
  # ------------------------------------------------------------------
  dose_3a_tbl <- as.data.frame(confint(pairs(emm_3a_dose)))
  dose_3a_tbl$estimate <- round(dose_3a_tbl$estimate, 3)
  dose_3a_tbl$SE <- round(dose_3a_tbl$SE, 3)
  dose_3a_tbl$df <- round(dose_3a_tbl$df, 1)
  dose_3a_tbl$lower.CL <- round(dose_3a_tbl$lower.CL, 3)
  dose_3a_tbl$upper.CL <- round(dose_3a_tbl$upper.CL, 3)
  names(dose_3a_tbl)[names(dose_3a_tbl) == "lower.CL"] <- "CI lower"
  names(dose_3a_tbl)[names(dose_3a_tbl) == "upper.CL"] <- "CI upper"
  doc <- body_add_par(doc, "EMM Dose Contrasts: Model 3a", style = "heading 2")
  ft <- flextable(dose_3a_tbl) |> autofit() |>
    set_caption("Pairwise dose contrasts within each phase class (Model 3a, Tukey-adjusted)")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # EMM Phase Contrasts: Model 3a
  phase_3a_tbl <- as.data.frame(confint(pairs(emm_3a_phase)))
  phase_3a_tbl$estimate <- round(phase_3a_tbl$estimate, 3)
  phase_3a_tbl$SE <- round(phase_3a_tbl$SE, 3)
  phase_3a_tbl$df <- round(phase_3a_tbl$df, 1)
  phase_3a_tbl$lower.CL <- round(phase_3a_tbl$lower.CL, 3)
  phase_3a_tbl$upper.CL <- round(phase_3a_tbl$upper.CL, 3)
  names(phase_3a_tbl)[names(phase_3a_tbl) == "lower.CL"] <- "CI lower"
  names(phase_3a_tbl)[names(phase_3a_tbl) == "upper.CL"] <- "CI upper"
  doc <- body_add_par(doc, "EMM Phase Contrasts: Model 3a", style = "heading 2")
  ft <- flextable(phase_3a_tbl) |> autofit() |>
    set_caption("Phase contrasts at each dose level (Model 3a)")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # ------------------------------------------------------------------
  # EMM Dose Contrasts: Model 3c (ANCOVA, primary)
  # ------------------------------------------------------------------
  dose_3c_tbl <- as.data.frame(confint(pairs(emm_3c_dose)))
  dose_3c_tbl$estimate <- round(dose_3c_tbl$estimate, 3)
  dose_3c_tbl$SE <- round(dose_3c_tbl$SE, 3)
  dose_3c_tbl$df <- round(dose_3c_tbl$df, 1)
  dose_3c_tbl$lower.CL <- round(dose_3c_tbl$lower.CL, 3)
  dose_3c_tbl$upper.CL <- round(dose_3c_tbl$upper.CL, 3)
  names(dose_3c_tbl)[names(dose_3c_tbl) == "lower.CL"] <- "CI lower"
  names(dose_3c_tbl)[names(dose_3c_tbl) == "upper.CL"] <- "CI upper"
  doc <- body_add_par(doc, "EMM Dose Contrasts: Model 3c (ANCOVA)", style = "heading 2")
  ft <- flextable(dose_3c_tbl) |> autofit() |>
    set_caption("Pairwise dose contrasts within each phase class (Model 3c ANCOVA, Tukey-adjusted)")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # EMM Phase Contrasts: Model 3c
  phase_3c_tbl <- as.data.frame(confint(pairs(emm_3c_phase)))
  phase_3c_tbl$estimate <- round(phase_3c_tbl$estimate, 3)
  phase_3c_tbl$SE <- round(phase_3c_tbl$SE, 3)
  phase_3c_tbl$df <- round(phase_3c_tbl$df, 1)
  phase_3c_tbl$lower.CL <- round(phase_3c_tbl$lower.CL, 3)
  phase_3c_tbl$upper.CL <- round(phase_3c_tbl$upper.CL, 3)
  names(phase_3c_tbl)[names(phase_3c_tbl) == "lower.CL"] <- "CI lower"
  names(phase_3c_tbl)[names(phase_3c_tbl) == "upper.CL"] <- "CI upper"
  doc <- body_add_par(doc, "EMM Phase Contrasts: Model 3c (ANCOVA)", style = "heading 2")
  ft <- flextable(phase_3c_tbl) |> autofit() |>
    set_caption("Phase contrasts at each dose level (Model 3c ANCOVA)")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # ------------------------------------------------------------------
  # Effect Sizes (Cohen's d): Models 3a and 3c
  # ------------------------------------------------------------------
  # Dose within phase
  es_3a_dose <- as.data.frame(confint(eff_size(emm_3a_dose, sigma = sigma_trial, edf = edf_trial)))
  es_3a_dose$effect.size <- round(es_3a_dose$effect.size, 3)
  es_3a_dose$SE <- round(es_3a_dose$SE, 3)
  es_3a_dose$lower.CL <- round(es_3a_dose$lower.CL, 3)
  es_3a_dose$upper.CL <- round(es_3a_dose$upper.CL, 3)
  es_3a_dose$df <- round(es_3a_dose$df, 1)
  names(es_3a_dose)[names(es_3a_dose) == "effect.size"] <- "Cohen's d"
  names(es_3a_dose)[names(es_3a_dose) == "lower.CL"] <- "CI lower"
  names(es_3a_dose)[names(es_3a_dose) == "upper.CL"] <- "CI upper"
  doc <- body_add_par(doc, "Effect Sizes: Model 3a Dose Contrasts", style = "heading 2")
  ft <- flextable(es_3a_dose) |> autofit() |>
    set_caption(sprintf("Cohen's d for dose contrasts (Model 3a). Denominator sigma = %.1f uV (trial-level residual SD from Model 4).", sigma_trial))
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  es_3c_dose <- as.data.frame(confint(eff_size(emm_3c_dose, sigma = sigma_trial, edf = edf_trial)))
  es_3c_dose$effect.size <- round(es_3c_dose$effect.size, 3)
  es_3c_dose$SE <- round(es_3c_dose$SE, 3)
  es_3c_dose$lower.CL <- round(es_3c_dose$lower.CL, 3)
  es_3c_dose$upper.CL <- round(es_3c_dose$upper.CL, 3)
  es_3c_dose$df <- round(es_3c_dose$df, 1)
  names(es_3c_dose)[names(es_3c_dose) == "effect.size"] <- "Cohen's d"
  names(es_3c_dose)[names(es_3c_dose) == "lower.CL"] <- "CI lower"
  names(es_3c_dose)[names(es_3c_dose) == "upper.CL"] <- "CI upper"
  doc <- body_add_par(doc, "Effect Sizes: Model 3c Dose Contrasts", style = "heading 2")
  ft <- flextable(es_3c_dose) |> autofit() |>
    set_caption(sprintf("Cohen's d for dose contrasts (Model 3c ANCOVA). Denominator sigma = %.1f uV.", sigma_trial))
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # ------------------------------------------------------------------
  # Save .docx
  # ------------------------------------------------------------------
  docx_path <- paste0(outputDir, "/betaStim_statistical_tables.docx")
  print(doc, target = docx_path)
  cat("Saved manuscript tables to:", docx_path, "\n")

} else {
  cat("Install officer and flextable packages for .docx export:\n")
  cat("  install.packages(c('officer', 'flextable'))\n")
}
