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
library('marginaleffects')
library('dplyr')
library('afex')
library('report')
library('effectsize')
library('performance')
library('car')
library('officer')
library('flextable')

# Global contrast coding. contr.sum gives orthogonal (sum-to-zero) contrasts
# for unordered factors so Type III ANOVAs on factors involved in interactions
# yield marginal main effects. contr.poly (R default) retained for ordered
# factors. Must be set BEFORE factors are created/used in model fits.
options(contrasts = c("contr.sum", "contr.poly"))

# log data prior to fitting?
log_data = FALSE

# Print HTML tab_model() summaries? Off by default (opens RStudio Viewer).
showTabModel = FALSE

# Run the phase-quality threshold sensitivity analysis?
# Fits 10 lmer models (5a and 5a-gf2 at 5 thresholds each) — adds ~30-60s.
# Opt-in for manuscript sensitivity tables; off for routine runs.
runPhaseSensitivity = FALSE

# ecb43e setToDeliverPhase value marking random-phase (non-targeted) trials.
# These are kept for pure-dose plots but excluded from phase-dependent models.
RANDOM_PHASE_MARKER <- "12345"
excludeRandomPhase = TRUE

# Phase-quality filters. Applied per-model inside each fit block so each
# model can use the threshold appropriate for its sample size and random
# effects structure. 5a-gf2 uses 0.2 (not 0.3) because 0.3 is singular in
# the good-fit subset. 5a-gf uses per-burst scope to match its predictor.
minPhaseVecLength_5a  = 0.3   # Model 5a channel r threshold
minPhaseVecLength_gf2 = 0.2   # Model 5a-gf2 channel r threshold
minBurstVecLength_gf  = 0     # Model 5a-gf per-burst r threshold
minGoodBetaPerBurst   = 1     # minimum R²>0.7 stims per burst

# Drop rows where `col` is NA or below `min_val`. If `min_val <= 0`,
# returns df unchanged. Logs before/after counts when `label` is given.
apply_min_filter <- function(df, col, min_val, label = NULL) {
  if (is.null(min_val) || min_val <= 0) return(df)
  keep <- !is.na(df[[col]]) & df[[col]] >= min_val
  if (!is.null(label)) {
    cat(sprintf("%s filter (%s >= %g): %d -> %d trials\n",
                label, col, min_val, nrow(df), sum(keep)))
  }
  df[keep, ]
}

# Build a phase-curve plot title and filename stem from the filter
# threshold. Kept in a single helper so title and filename never drift.
# model_name like "Model 5a" / "Model 5a-gf2" becomes stem "model5a" /
# "model5a_gf2" (spaces removed, hyphens → underscores).
phase_curve_label <- function(model_name, r_thresh, n_obs) {
  stem <- gsub("[^A-Za-z0-9_]", "",
               gsub("-", "_",
                    gsub(" ", "", tolower(model_name))))
  if (r_thresh > 0) {
    list(
      title = sprintf("%s: Phase-Response Curve by Dose (phaseVecLength >= %.2f, N = %d)",
                      model_name, r_thresh, n_obs),
      fname = sprintf("betaStim_%s_phase_curve_r%02d", stem, round(100 * r_thresh)))
  } else {
    list(
      title = sprintf("%s: Phase-Response Curve by Dose (no quality filter, N = %d)",
                      model_name, n_obs),
      fname = sprintf("betaStim_%s_phase_curve", stem))
  }
}

savePlot = 1
figWidth = 8
figHeight = 6

# ------------------------------------------------------------------------
here()
data <- read.table(here("data","output_table","betaStim_outputTable_50_new_100_thresh.csv"),header=TRUE,sep = ",",stringsAsFactors=F,
                   colClasses=c("magnitude"="numeric","betaLabels"="factor","sid"="factor","numStims"="factor","stimLevel"="numeric","channel"="factor","subjectNum"="factor","phaseClass"="factor","phaseDeg"="numeric","setToDeliverPhase"="factor",'phaseDeliveryBinned45'="factor","phaseVecLength"="numeric","phaseCircStd"="numeric","phaseOmnibusP"="numeric"))

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

# --- Circular phase decomposition (Fisher 1993) ---
# Convert channel-level circular mean phase (degrees) to sin/cos predictors.
# Standard approach for including angular variables in linear models:
#   sin(phase) captures the 90-270 axis (positive = toward 90 deg)
#   cos(phase) captures the 0-180 axis (positive = toward 0 deg)
data$phase_rad <- data$phaseDeg * pi / 180
data$sin_phase <- sin(data$phase_rad)
data$cos_phase <- cos(data$phase_rad)
# Rounded phase as condition identifier for Model 5 grouping.
# phaseDeg is constant within each (channel x condition) pair (MATLAB repmat),
# so rounding is just float-safety. Keeps conditions with different measured
# phases separate, unlike phaseClass which can merge two conditions into one bin.
data$phaseDeg_round <- round(data$phaseDeg, 1)

# --- Per-subject baseline CEP magnitude summary (descriptive table) ---
# One row per subject summarizing the baseline evoked potential magnitude
# distribution across all baseline probe trials pooled over channels.
# Sorted by subjectNum (manuscript numeric ID); subject column formatted as
# "Subject N (sid)" for cross-reference. Written to CSV + docx as
# subject-characteristics Table 1. Uses the post-filter `data` (25-1500 uV
# trial filter applied, 0b5a2ePlayBack and Null bursts already removed) so
# the numbers reflect trials that entered the analysis.
subject_cep_summary <- data |>
  dplyr::filter(numStims == "Base") |>
  dplyr::group_by(subjectNum, sid) |>
  dplyr::summarise(
    n_trials   = dplyr::n(),
    n_channels = dplyr::n_distinct(channel),
    median_mag = round(median(magnitude), 1),
    mad_mag    = round(mad(magnitude), 1),
    q25_mag    = round(quantile(magnitude, 0.25), 1),
    q75_mag    = round(quantile(magnitude, 0.75), 1),
    .groups    = "drop") |>
  dplyr::arrange(as.integer(as.character(subjectNum))) |>
  dplyr::mutate(subject_label = sprintf("Subject %s (%s)", subjectNum, sid)) |>
  dplyr::select(subject_label, n_trials, n_channels,
                median_mag, mad_mag, q25_mag, q75_mag)
cat("\n--- Per-subject baseline CEP magnitude (Base trials only) ---\n")
print(subject_cep_summary)
write.csv(subject_cep_summary,
          here("output_plots", "betaStim_subject_cep_summary.csv"),
          row.names = FALSE)

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

summaryData = ddply(data, .(sid,phaseClass,numStims,channel,betaLabels), summarize,
                    magnitude = median(magnitude), sin_phase = first(sin_phase), cos_phase = first(cos_phase))

summaryData = ddply(data[data$numStims != "Base",] , .(sid,phaseClass,numStims,channel,betaLabels), summarize,
                    percentDiff = median(percentDiff), sin_phase = first(sin_phase), cos_phase = first(cos_phase))

dataNoBaseline = data[data$numStims != "Base",]

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
  ggsave(here("output_plots","betaStim_dose_phase.png"), plot=p2,
         units="in", width=figWidth, height=figHeight, dpi=600)
  ggsave(here("output_plots","betaStim_dose_phase.eps"), plot=p2,
         units="in", width=figWidth, height=figHeight, dpi=600, device=cairo_ps)
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
  ggsave(here("output_plots","betaStim_dose.png"), plot=p2,
         units="in", width=figWidth, height=figHeight, dpi=600)
  ggsave(here("output_plots","betaStim_dose.eps"), plot=p2,
         units="in", width=figWidth, height=figHeight, dpi=600, device=cairo_ps)
}

# ------------------------------------------------------------------------
# Exclude ecb43e random-phase trials from phase-related analyses below.
# The pure-dose plot (betaStim_dose.png) above includes them because it
# doesn't use phase. Everything downstream that depends on phase (models
# 3a-3e, 5a, 5a-gf, 5a-gf2, 6, residual diagnostics, ecb43e random has
# phaseDeg = noisy circular mean → pollutes phase predictors).
if (excludeRandomPhase) {
  n_before <- nrow(data)
  data <- data[!(data$sid == "ecb43e" &
                 as.character(data$setToDeliverPhase) == RANDOM_PHASE_MARKER), ]
  n_after <- nrow(data)
  cat(sprintf("\nExcluded %d ecb43e random-phase trials from phase analyses (%d → %d rows)\n",
      n_before - n_after, n_before, n_after))

  # rebuild summary tables that include phaseClass / phaseDeg predictors so
  # downstream models and plots reflect the filter (match original defs at
  # lines 109-120)
  summaryData <- ddply(data[data$numStims != "Base",],
    .(sid, phaseClass, numStims, channel, betaLabels),
    summarize, percentDiff = median(percentDiff),
    sin_phase = first(sin_phase), cos_phase = first(cos_phase))
  summaryDataHighStimsOnly <- ddply(data[data$numStims == "[5,inf)",],
    .(sid, phaseClass, numStims, channel, orderedPhase45),
    summarize, percentDiff = median(percentDiff))
  dataNoBaseline <- data[data$numStims != "Base",]
}

# Phase-quality filtering is applied per-model below (see 5a and 5a-gf2
# fit blocks), not globally, so each model can use the threshold that
# fits its sample size and random-effects structure.

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
  ggsave(here("output_plots","betaStim_dose_no_dots.png"), plot=p3,
         units="in", width=figWidth, height=figHeight, dpi=600)
  ggsave(here("output_plots","betaStim_dose_no_dots.eps"), plot=p3,
         units="in", width=figWidth, height=figHeight, dpi=600, device=cairo_ps)
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
  ggsave(here("output_plots","betaStim_dose_binned45.png"), plot=p4,
         units="in", width=figWidth, height=figHeight, dpi=600)
  ggsave(here("output_plots","betaStim_dose_binned45.eps"), plot=p4,
         units="in", width=figWidth, height=figHeight, dpi=600, device=cairo_ps)
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

if (showTabModel) tab_model(fit.intercepts.only)

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

if (showTabModel) tab_model(fit.trial.level)

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

if (showTabModel) tab_model(fit.absDiff)

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

if (showTabModel) tab_model(fit.modelD)

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

if (showTabModel) tab_model(fit.ancova)

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

# residual diagnostics for Model 3c (sensitivity model; 5a/5a-gf2 diagnostics below)
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

# ========================================================================
# Model 5: Prepare summary datasets with phaseDeg_round grouping
# ========================================================================
# Group by (sid, phaseDeg_round, numStims, channel) instead of phaseClass.
# phaseDeg_round keeps conditions with different measured phases separate,
# even when both would be binned to the same phaseClass (e.g., both < 180).
# sin_phase/cos_phase are constant within each cell (same measured phase),
# so first() is safe here.

dataNB_m5 <- apply_min_filter(dataNoBaseline, "phaseVecLength",
  minPhaseVecLength_5a, label = "Model 5a")

summaryNB_m5 <- ddply(dataNB_m5, .(sid, phaseDeg_round, numStims, channel),
  summarize, magnitude = median(magnitude),
  sin_phase = first(sin_phase), cos_phase = first(cos_phase),
  betaLabels = first(betaLabels))

# baseline covariate: one median per channel, pooled across conditions (unchanged)
summaryNB_m5 <- merge(summaryNB_m5, basePerChan, by = c("sid", "channel"))
summaryNB_m5$baselineMag_c <- summaryNB_m5$baselineMag - mean(summaryNB_m5$baselineMag)
summaryNB_m5$doseNum <- as.numeric(factor(summaryNB_m5$numStims,
  levels = c("[1,2]","[3,4]","[5,inf)"))) - 1

# ordinal dose coding
summaryNB_m5$numStims_ord <- ordered(summaryNB_m5$numStims,
  levels = c("[1,2]", "[3,4]", "[5,inf)"))
summaryNB_m5$dose_linpoly <- poly_lin[as.numeric(summaryNB_m5$numStims_ord)]

cat(sprintf("\nModel 5 summary data: %d obs, %d subjects, %d channels\n",
    nrow(summaryNB_m5), length(unique(summaryNB_m5$sid)),
    length(unique(summaryNB_m5$channel))))

# ========================================================================
# Model 5: Continuous circular phase (sin/cos decomposition) — ANCOVA
# ========================================================================
# Replaces binary phaseClass (90/270) with sin(phase) and cos(phase), the
# standard approach for including angular predictors in linear models
# (Fisher 1993, "Statistical Analysis of Circular Data").
#
# Why sin/cos decomposition:
#   phaseClass bins all phases 0-180 as "90" and 181-360 as "270", discarding
#   the continuous phase information. Channels with delivered phases of 85 and
#   175 degrees are both labelled "90" despite being ~90 degrees apart.
#   The sin/cos decomposition preserves the full 360 degrees of phase:
#     sin(phase): captures the 90-270 deg axis (positive = toward 90)
#     cos(phase): captures the 0-180 deg axis (positive = toward 0)
#   A joint test of sin + cos = omnibus test for any phase effect direction.
#
# Why ANCOVA (baselineMag_c as covariate):
#   Baseline CEP magnitude varies ~10-fold across channels (~50-600+ uV) due
#   to electrode proximity to cortical generators. Including baselineMag_c
#   (grand-mean-centered baseline median per channel) absorbs this between-
#   channel variance (channel random SD drops from ~153 to ~10 uV in existing
#   models). Baseline trials are excluded from the modeled data — they serve
#   only as the covariate source — avoiding the confound that arises when
#   baseline is a dose level (baseline trials inherit phaseClass/sin_phase
#   labels even though no phase-targeted stimulation occurred).
#   Standard for pre/post designs (Senn 2006; Van Breukelen 2006).
#
# Hypothesis: interaction between delivered phase and conditioning dose —
# the phase effect on CEP magnitude depends on number of conditioning stimuli.
#
# betaLabels (1 = beta reference channel, 0 = EP channel) included as a
# fixed effect to control for systematic magnitude differences at the beta
# recording electrode.
# ========================================================================

# --- Model 5a: Ordinal dose x sin/cos + betaLabels (ANCOVA) ---
# Primary model. Uses polynomial contrasts (.L linear, .Q quadratic) for dose.
# REML for parameter estimation; ML comparisons via explicit REML=FALSE fits.
cat("\n=== Model 5a: Continuous phase (sin/cos) + ordinal dose + ANCOVA ===\n")

fit.sincos.ordinal = lmerTest::lmer(
  magnitude ~ numStims_ord * (sin_phase + cos_phase) + betaLabels + baselineMag_c +
  (1 + dose_linpoly | sid) + (1 | channel),
  data = summaryNB_m5,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

cat("Singular:", isSingular(fit.sincos.ordinal), "\n")
summary(fit.sincos.ordinal)
cat("\nType III ANOVA (Satterthwaite):\n")
print(anova(fit.sincos.ordinal))
VarCorr(fit.sincos.ordinal)
report(fit.sincos.ordinal)
model_performance(fit.sincos.ordinal)

# --- Joint hypothesis tests ---
# (1) Omnibus phase main effect: are sin_phase and cos_phase jointly zero?
#     With polynomial contrasts, main effects represent the average phase
#     effect across all dose levels.
#
# (2) Phase x dose interaction: do the dose:sin and dose:cos interactions
#     jointly contribute?
#
# Two approaches:
#   (a) LRT (likelihood ratio test): compare nested ML-fitted models.
#       Preferred for small samples — chi-squared approximation is more
#       reliable than the Wald test from linearHypothesis.
#   (b) Wald test (linearHypothesis): complementary check using the
#       coefficient covariance matrix. Can be liberal with few clusters
#       (6 subjects), so we report alongside LRT rather than in isolation.

# -- (a) LRT comparisons (REML=FALSE required for differing fixed effects) --
fit.sincos.ordinal.ml = update(fit.sincos.ordinal, REML = FALSE)

# Reduced: no phase terms at all
fit.sincos.nophase.ml = lmerTest::lmer(
  magnitude ~ numStims_ord + betaLabels + baselineMag_c +
  (1 + dose_linpoly | sid) + (1 | channel),
  data = summaryNB_m5, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

cat("\nLRT: full model vs no-phase model (omnibus phase test, 6 df):\n")
print(anova(fit.sincos.nophase.ml, fit.sincos.ordinal.ml))

# Reduced: phase main effects only, no interaction
fit.sincos.noint.ml = lmerTest::lmer(
  magnitude ~ numStims_ord + sin_phase + cos_phase + betaLabels + baselineMag_c +
  (1 + dose_linpoly | sid) + (1 | channel),
  data = summaryNB_m5, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

cat("\nLRT: full model vs main-effects-only (dose x phase interaction, 4 df):\n")
print(anova(fit.sincos.noint.ml, fit.sincos.ordinal.ml))

# Reduced: dose + interaction, no phase main effects
# (tests whether phase main effects add beyond the interaction)
cat("\nLRT: no-phase vs main-effects-only (phase main effect, 2 df):\n")
print(anova(fit.sincos.nophase.ml, fit.sincos.noint.ml))

# -- (b) Wald tests (linearHypothesis) — complementary to LRT --
# Average phase modulation across doses (2-df joint test of main effects).
# Note: Wald chi-squared can be liberal with few clusters; interpret alongside LRT.
cat("\nWald test: average phase modulation (sin_phase = cos_phase = 0):\n")
print(car::linearHypothesis(fit.sincos.ordinal,
  c("sin_phase = 0", "cos_phase = 0")))

# Phase x dose interaction (4-df joint test of all interaction terms)
cat("\nWald test: phase x dose interaction:\n")
print(car::linearHypothesis(fit.sincos.ordinal,
  c("numStims_ord.L:sin_phase = 0", "numStims_ord.Q:sin_phase = 0",
    "numStims_ord.L:cos_phase = 0", "numStims_ord.Q:cos_phase = 0")))

# --- emmeans at key phase angles ---
# For a single phase angle, at= gives one sin/cos value each — no factorial issue.

# --- emmeans at key phase angles ---
# All pairwise dose contrasts use Tukey adjustment (default for pairs()).
# Phase contrasts (90 vs 270) are single comparisons — no adjustment needed.
# Dose contrasts at phase = 90 deg (depolarizing: sin=1, cos=0)
cat("\n--- emmeans: dose contrasts at phase = 90 deg (Tukey-adjusted) ---\n")
# betaLabels omitted so emmeans marginalizes equally across non-beta / beta
# channel types. The fixed effect of betaLabels is non-significant (p=0.52
# in 5a, 0.70 in 5a-gf2), so marginalization shifts EMMs by <1.3 uV but
# gives a cleaner population-level interpretation for the docx tables.
emm_5a_90 <- emmeans(fit.sincos.ordinal, ~ numStims_ord,
  at = list(sin_phase = 1, cos_phase = 0, baselineMag_c = 0))
print(emm_5a_90)
cat("Pairwise dose contrasts at 90 deg:\n")
print(confint(pairs(emm_5a_90)))

# Dose contrasts at phase = 270 deg (hyperpolarizing: sin=-1, cos=0)
cat("\n--- emmeans: dose contrasts at phase = 270 deg (Tukey-adjusted) ---\n")
emm_5a_270 <- emmeans(fit.sincos.ordinal, ~ numStims_ord,
  at = list(sin_phase = -1, cos_phase = 0, baselineMag_c = 0))
print(emm_5a_270)
cat("Pairwise dose contrasts at 270 deg:\n")
print(confint(pairs(emm_5a_270)))

# Phase contrast (90 vs 270) at each dose level — cos=0 for both, no factorial issue
cat("\n--- emmeans: phase 90 vs 270 at each dose ---\n")
emm_5a_phase <- emmeans(fit.sincos.ordinal, ~ numStims_ord * sin_phase,
  at = list(sin_phase = c(1, -1), cos_phase = 0, betaLabels = "0", baselineMag_c = 0))
print(contrast(emm_5a_phase, method = "pairwise", by = "numStims_ord"))

# --- Phase-response curve via emmeans (every 45 deg) ---
# betaLabels omitted from at=... and weights = "proportional" used so the
# curve marginalizes over betaLabels using the OBSERVED proportions in
# summaryNB_m5 (rather than pinning at "0" or using emmeans' default
# equal-weight averaging). This matches the G-computation / AME semantics
# used by the marginaleffects dose-effect-vs-phase curve below and by the
# AME numeric docx table — giving the plot a population-level rather than
# non-beta-specific interpretation. baselineMag_c is grand-mean centered
# so pinning at 0 is equivalent to averaging (E[baselineMag_c] = 0).
# Faceted beta-vs-non-beta version preserved separately as p_5a_beta.
phase_vals <- seq(0, 315, by = 45)
emm_curve <- lapply(phase_vals, function(ph) {
  em <- emmeans(fit.sincos.ordinal, ~ numStims_ord,
    at = list(sin_phase = sin(ph * pi / 180), cos_phase = cos(ph * pi / 180),
              baselineMag_c = 0),
    weights = "proportional")
  df <- as.data.frame(em)
  df$phase_deg <- ph
  df
})
emm_curve_df <- do.call(rbind, emm_curve)

lbl_5a <- phase_curve_label("Model 5a", minPhaseVecLength_5a, nrow(summaryNB_m5))
p_5a <- ggplot(emm_curve_df, aes(x = phase_deg, y = emmean, color = numStims_ord)) +
  theme_light(base_size = 14) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = numStims_ord), alpha = 0.15, color = NA) +
  labs(x = "Delivered Phase (degrees)",
       y = expression(paste("Predicted Magnitude (", mu, "V)")),
       color = "Dose", fill = "Dose",
       title = lbl_5a$title,
       subtitle = "emmeans +/- 95% CI; baseline at grand mean; betaLabels marginalized over observed proportions") +
  scale_x_continuous(breaks = seq(0, 315, by = 45))
p_5a

if(savePlot){
  ggsave(here("output_plots", paste0(lbl_5a$fname, ".png")), plot = p_5a,
         units = "in", width = 7, height = 4.5, dpi = 600)
  ggsave(here("output_plots", paste0(lbl_5a$fname, ".eps")), plot = p_5a,
         units = "in", width = 7, height = 4.5, dpi = 600, device = cairo_ps)
}

# --- Dose-effect-vs-phase curve (marginaleffects) ---
# Complements p_5a. p_5a shows predicted MAGNITUDE at each phase × dose;
# this shows the pairwise dose CONTRAST (effect in uV) as a function of
# phase. Answers: "where around the circle is the dose effect largest?"
# — a phase-stratified decomposition of the AME numeric table.
#
# At each hypothetical phase angle, override sin_phase/cos_phase on the
# observed data and call avg_comparisons() — which then averages the
# pairwise dose contrast across the OBSERVED distribution of betaLabels,
# baselineMag_c, sid, and channel (rather than pinning them). Exactly
# parallel to the AME numeric block, just stratified by phase. Pinning
# would give a conditional-on-covariates curve; this gives a population-
# level curve, consistent with the AME table caption.
# re.form = NA ignores random effects for population-level prediction.
# CIs are pointwise (no adjustment across the 12 phase x 3 contrast = 36
# tests); the AME table gives the FWER-controlled single number, this
# curve decomposes it descriptively.
phase_vals_curve <- seq(0, 330, by = 30)
cmp_5a_curve_df <- do.call(rbind, lapply(phase_vals_curve, function(ph) {
  nd <- summaryNB_m5
  nd$sin_phase <- sin(ph * pi / 180)
  nd$cos_phase <- cos(ph * pi / 180)
  df <- as.data.frame(avg_comparisons(fit.sincos.ordinal,
    variables = list(numStims_ord = "pairwise"),
    newdata = nd, re.form = NA))
  df$phase_deg <- ph
  df
}))

p_5a_effect <- ggplot(cmp_5a_curve_df,
    aes(x = phase_deg, y = estimate, color = contrast, fill = contrast)) +
  theme_light(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  labs(x = "Delivered Phase (degrees)",
       y = expression(paste("Dose Effect (", mu, "V)")),
       color = "Pairwise contrast", fill = "Pairwise contrast",
       title = paste0(lbl_5a$title, " - Dose Effect vs Phase"),
       subtitle = "avg_comparisons() at each phase, averaged over observed betaLabels and baseline; pointwise 95% CI") +
  scale_x_continuous(breaks = seq(0, 360, by = 45))
p_5a_effect

if(savePlot){
  ggsave(here("output_plots", paste0(lbl_5a$fname, "_dose_effect_vs_phase.png")),
         plot = p_5a_effect, units = "in", width = 7, height = 4.5, dpi = 600)
  ggsave(here("output_plots", paste0(lbl_5a$fname, "_dose_effect_vs_phase.eps")),
         plot = p_5a_effect, units = "in", width = 7, height = 4.5, dpi = 600,
         device = cairo_ps)
}

# --- Supplementary: phase-response curves faceted by betaLabels ---
# Same reference grid as above but with betaLabels as a predictor dimension
# rather than pinned at "0". Generates one sub-panel per channel type (non-
# beta vs beta-trigger). Useful sanity check since the betaLabels fixed
# effect was non-significant — the two panels should look nearly identical.
emm_curve_beta_5a <- do.call(rbind, lapply(c("0", "1"), function(bl) {
  do.call(rbind, lapply(phase_vals, function(ph) {
    em <- emmeans(fit.sincos.ordinal, ~ numStims_ord,
      at = list(sin_phase = sin(ph * pi / 180), cos_phase = cos(ph * pi / 180),
                betaLabels = bl, baselineMag_c = 0))
    df <- as.data.frame(em)
    df$phase_deg <- ph
    df$betaLabels <- bl
    df
  }))
}))
emm_curve_beta_5a$channel_type <- factor(
  ifelse(emm_curve_beta_5a$betaLabels == "0", "Non-beta (EP) channels",
         "Beta trigger channels"),
  levels = c("Non-beta (EP) channels", "Beta trigger channels"))

p_5a_beta <- ggplot(emm_curve_beta_5a, aes(x = phase_deg, y = emmean, color = numStims_ord)) +
  theme_light(base_size = 14) +
  facet_wrap(~ channel_type, nrow = 1) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = numStims_ord),
              alpha = 0.15, color = NA) +
  labs(x = "Delivered Phase (degrees)",
       y = expression(paste("Predicted Magnitude (", mu, "V)")),
       color = "Dose", fill = "Dose",
       title = paste(lbl_5a$title, "- Beta vs Non-beta Channels"),
       subtitle = "emmeans +/- 95% CI; baseline at grand mean. betaLabels fixed effect p = 0.52 (ns)") +
  scale_x_continuous(breaks = seq(0, 315, by = 90)) +
  theme(strip.text = element_text(face = "bold"))

if (savePlot) {
  ggsave(here("output_plots", paste0(lbl_5a$fname, "_beta_vs_nonbeta.png")),
         plot = p_5a_beta, units = "in", width = 11, height = 4.5, dpi = 600)
  ggsave(here("output_plots", paste0(lbl_5a$fname, "_beta_vs_nonbeta.eps")),
         plot = p_5a_beta, units = "in", width = 11, height = 4.5, dpi = 600,
         device = cairo_ps)
}

# --- Total-variance effect sizes ---
# Standard (conditional) Cohen's d uses residual SD only. Marginal (unconditional)
# d uses total variance = sigma^2_subject + sigma^2_channel + sigma^2_residual.
# This answers: "how large is the effect relative to ALL variability in CEP
# magnitude across subjects, channels, and residual?"
# More conservative, more generalizable to new subjects/channels.
# Reference: Westfall, Kenny & Judd (2014), JEPG.

# For models with random slopes, the marginal variance depends on the
# predictor value: Var(Y|x) = Var_int + 2*x*Cov + x^2*Var_slope + ...
# Average over the observed dose distribution: E[x]=0 (polynomial centered),
# E[x^2] = mean(contr.poly(3)[,1]^2) = 1/3. So covariance term drops out
# and slope variance is weighted by 1/3, not 1.
vc_5a <- VarCorr(fit.sincos.ordinal)
var_sid_int <- attr(vc_5a$sid, "stddev")["(Intercept)"]^2
var_sid_slope <- attr(vc_5a$sid, "stddev")["dose_linpoly"]^2
var_channel <- attr(vc_5a$channel, "stddev")["(Intercept)"]^2
var_resid <- sigma(fit.sincos.ordinal)^2
ex2 <- mean(contr.poly(3)[,1]^2)  # = 1/3 for 3 dose levels
total_sd_5a <- sqrt(var_sid_int + ex2 * var_sid_slope + var_channel + var_resid)
resid_sd_5a <- sigma(fit.sincos.ordinal)

cat(sprintf("\nVariance components: total SD = %.1f uV, residual SD = %.1f uV\n",
    total_sd_5a, resid_sd_5a))

# Effect sizes for dose contrasts at phase=90 and phase=270
compute_effect_sizes <- function(emm_obj, total_sd, resid_sd, label) {
  contr <- pairs(emm_obj)
  contr_df <- as.data.frame(confint(contr))
  contr_df$d_total <- contr_df$estimate / total_sd
  contr_df$d_total_lower <- contr_df$lower.CL / total_sd
  contr_df$d_total_upper <- contr_df$upper.CL / total_sd
  contr_df$d_conditional <- contr_df$estimate / resid_sd
  cat(sprintf("\n--- Effect sizes: %s ---\n", label))
  cat(sprintf("d_total denominator (total SD): %.1f uV\n", total_sd))
  cat(sprintf("d_conditional denominator (residual SD): %.1f uV\n", resid_sd))
  print(contr_df[, c("contrast", "estimate", "d_total", "d_total_lower", "d_total_upper", "d_conditional")])
  return(contr_df)
}

es_90 <- compute_effect_sizes(emm_5a_90, total_sd_5a, resid_sd_5a, "dose at phase=90")
es_270 <- compute_effect_sizes(emm_5a_270, total_sd_5a, resid_sd_5a, "dose at phase=270")

# --- Average Marginal Effect of dose (marginaleffects) ---
# emm_5a_90 / emm_5a_270 evaluate dose contrasts at two specific phases
# (sin=+/-1, cos=0) with baselineMag_c and betaLabels pinned at 0. The AME
# below instead marginalizes across the OBSERVED joint distribution of
# sin_phase, cos_phase, baselineMag_c, and betaLabels in summaryNB_m5 —
# the phase-weighted population-level dose effect. These are complementary:
# the emmeans output answers "what is the dose effect at this phase?", the
# AME answers "what is the dose effect averaged across the phases this
# study actually sampled?".
#
# Multiplicity adjustment: single-step max-t via multcomp::glht — the closest
# available analog to the Tukey-Kramer adjustment that emmeans::pairs() uses
# by default for the neighboring EMM contrasts. marginaleffects 0.29 does
# not expose "tukey" directly; single-step is the same family of FWER
# control via the joint multivariate-t distribution of the Wald statistics.
# hypotheses() drops the contrast label column, so we merge it back from
# the unadjusted avg_comparisons object for a readable table.
#
# Inference uses fixed-effect covariance only (standard marginaleffects
# behavior for lmerMod — population-averaged fixed-effect uncertainty).
cat("\n--- AME: Model 5a dose pairwise contrasts, marginalized over observed covariates ---\n")
ame_5a_raw <- avg_comparisons(fit.sincos.ordinal,
  variables = list(numStims_ord = "pairwise"))
ame_5a_adj <- as.data.frame(hypotheses(ame_5a_raw, multcomp = "single-step"))
ame_5a_adj$contrast <- as.data.frame(ame_5a_raw)$contrast
ame_5a_dose <- ame_5a_adj[, c("term", "contrast",
  setdiff(names(ame_5a_adj), c("term", "contrast")))]
print(ame_5a_dose)

# --- Model 5b: Numeric dose x sin/cos (sensitivity) ---
cat("\n=== Model 5b: Continuous phase (sin/cos) + numeric dose + ANCOVA ===\n")

fit.sincos.numeric = lmerTest::lmer(
  magnitude ~ doseNum * (sin_phase + cos_phase) + betaLabels + baselineMag_c +
  (1 + doseNum | sid) + (1 | channel),
  data = summaryNB_m5,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

cat("Singular:", isSingular(fit.sincos.numeric), "\n")
summary(fit.sincos.numeric)
print(anova(fit.sincos.numeric))

# --- Model 5c: Categorical dose x sin/cos (sensitivity) ---
cat("\n=== Model 5c: Continuous phase (sin/cos) + categorical dose + ANCOVA ===\n")

fit.sincos.categ = lmerTest::lmer(
  magnitude ~ numStims * (sin_phase + cos_phase) + betaLabels + baselineMag_c +
  (1 + doseNum | sid) + (1 | channel),
  data = summaryNB_m5,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

cat("Singular:", isSingular(fit.sincos.categ), "\n")
summary(fit.sincos.categ)
print(anova(fit.sincos.categ))

# --- Model 5 comparison ---
cat("\n=== Model 5 variants comparison ===\n")
cat("Model 5a (ordinal):     AIC=", round(AIC(fit.sincos.ordinal),1), " BIC=", round(BIC(fit.sincos.ordinal),1), " singular=", isSingular(fit.sincos.ordinal), "\n")
cat("Model 5b (numeric):     AIC=", round(AIC(fit.sincos.numeric),1), " BIC=", round(BIC(fit.sincos.numeric),1), " singular=", isSingular(fit.sincos.numeric), "\n")
cat("Model 5c (categorical): AIC=", round(AIC(fit.sincos.categ),1), " BIC=", round(BIC(fit.sincos.categ),1), " singular=", isSingular(fit.sincos.categ), "\n")

# ========================================================================
# Model 6: Continuous dose (natural spline) x sin/cos phase — EXPLORATORY
# ========================================================================
# CAUTION: Models 6 and 6-gf are EXPLORATORY only. Each unique stim count
# becomes its own summary cell, producing many more observations (682-1280)
# than the binned models (96-153). With only 7 subjects, a random dose slope
# per subject is singular — the model cannot account for between-subject
# heterogeneity in dose-response. Without random slopes, within-subject
# dose observations are treated as if independent, making p-values and CIs
# anti-conservative (too liberal). Model 5a (binned dose, with proper random
# slopes) is the primary inferential model. Model 6 is retained for its
# descriptive value in visualizing the continuous dose-response shape.
#
# Replaces binned dose ([1,2], [3,4], [5,inf)) with the actual number of
# conditioning stims delivered per burst (nCondStims = 1, 2, 3, ...).
# Uses natural splines (ns) to capture potential nonlinear dose-response
# without imposing a parametric form.
#
# nCondStims is the TOTAL number of conditioning stims in the burst, not
# the number with good beta fits. It reflects the experimental manipulation
# (how much conditioning was delivered), regardless of measurement quality.
#
# nCondStims is sourced from *_burst_phase_precision.csv files, which have
# the raw burst stim count from bursts(4, burstId) in the stim table.
#
# The merge uses (probeSample, channel) as the key: probeSample is the TDT
# sample number uniquely identifying each probe event, and channel is the
# encoded channel ID (subjectNum*100 + chan). Together they give a unique
# (probe, channel) pair matching between the main data and precision CSVs.
cat("\n=== Model 6 [EXPLORATORY]: Continuous dose (spline) x sin/cos phase ===\n")
cat("NOTE: p-values are anti-conservative — no random dose slopes (singular at n=7 subjects).\n")

precision_files_m6 <- Sys.glob(here("data", "output_table", "*_burst_phase_precision.csv"))
if (length(precision_files_m6) > 0) {
  precision_combined <- do.call(rbind, lapply(precision_files_m6, read.csv))
}

if (length(precision_files_m6) > 0) {
  library(splines)

  precision_dose <- precision_combined[, c("probeSample", "channelEncoded", "nCondStims")]
  precision_dose$channelEncoded <- as.factor(precision_dose$channelEncoded)

  data_m6 <- merge(data, precision_dose,
    by.x = c("probeSample", "channel"),
    by.y = c("probeSample", "channelEncoded"),
    all.x = TRUE)

  # exclude baselines (nCondStims=0) and trials with no precision match
  data_m6_NB <- data_m6[!is.na(data_m6$nCondStims) & data_m6$nCondStims > 0, ]

  cat(sprintf("Merged nCondStims: %d non-baseline trials, range [%d, %d]\n",
      nrow(data_m6_NB), min(data_m6_NB$nCondStims), max(data_m6_NB$nCondStims)))
  cat("Distribution of nCondStims:\n")
  print(table(data_m6_NB$nCondStims))

  # summary: one median per (sid, channel, phaseDeg_round, nCondStims)
  # Using phaseDeg_round to keep conditions separate (same as Model 5).
  # Each unique stim count gets its own cell — more granular than binned dose.
  data_m6_NB$phaseDeg_round <- round(data_m6_NB$phaseDeg, 1)
  summaryNB_m6 <- ddply(data_m6_NB, .(sid, channel, phaseDeg_round, nCondStims),
    summarize, magnitude = median(magnitude),
    sin_phase = first(sin_phase), cos_phase = first(cos_phase),
    betaLabels = first(betaLabels))

  # baseline covariate (from unfiltered data, same as Model 5)
  summaryNB_m6 <- merge(summaryNB_m6, basePerChan, by = c("sid", "channel"))
  summaryNB_m6$baselineMag_c <- summaryNB_m6$baselineMag - mean(summaryNB_m6$baselineMag)

  cat(sprintf("Model 6 summary: %d obs, %d subjects, %d channels, nCondStims range [%d, %d]\n",
      nrow(summaryNB_m6), length(unique(summaryNB_m6$sid)),
      length(unique(summaryNB_m6$channel)),
      min(summaryNB_m6$nCondStims), max(summaryNB_m6$nCondStims)))

  # Natural spline with 3 df for dose — captures nonlinearity without overfitting.
  # Boundary knots at min/max of observed data, internal knots at quantiles.
  # Interaction with sin/cos tests whether the dose-response shape differs by phase.
  fit.spline = lmerTest::lmer(
    magnitude ~ ns(nCondStims, df = 3) * (sin_phase + cos_phase) + betaLabels + baselineMag_c +
    (1 + nCondStims | sid) + (1 | channel),
    data = summaryNB_m6,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

  cat("Singular:", isSingular(fit.spline), "\n")
  print(summary(fit.spline))
  cat("\nType III ANOVA (Satterthwaite):\n")
  print(anova(fit.spline))

  # --- LRT: phase effect ---
  fit.spline.ml <- update(fit.spline, REML = FALSE)
  fit.spline.nophase.ml <- lmerTest::lmer(
    magnitude ~ ns(nCondStims, df = 3) + betaLabels + baselineMag_c +
    (1 + nCondStims | sid) + (1 | channel),
    data = summaryNB_m6, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

  cat("\nLRT: full model vs no-phase (omnibus phase test):\n")
  print(anova(fit.spline.nophase.ml, fit.spline.ml))

  # --- Wald test: average phase modulation ---
  cat("\nWald test: average phase modulation:\n")
  print(car::linearHypothesis(fit.spline, c("sin_phase = 0", "cos_phase = 0")))

  # --- Dose-response curve via emmeans ---
  # emmeans handles ns() correctly — it reconstructs the spline basis at new
  # values using the same knots from the original model fit.
  # Evaluate at 90 deg (sin=1, cos=0) and 270 deg (sin=-1, cos=0).
  # Limit dose range to where data is dense (most data at nCondStims <= 20).
  max_dose_plot <- min(20, max(summaryNB_m6$nCondStims))
  spline_at <- list(
    nCondStims = seq(1, max_dose_plot, by = 0.5),
    sin_phase = c(1, -1),   # 90 deg and 270 deg
    cos_phase = 0,           # cos=0 for both
    betaLabels = "0",
    baselineMag_c = 0)

  emm_spline <- emmeans(fit.spline, ~ nCondStims * sin_phase, at = spline_at)
  emm_spline_df <- as.data.frame(emm_spline)
  emm_spline_df$phase_label <- ifelse(emm_spline_df$sin_phase > 0, "90 deg", "270 deg")

  p_m6 <- ggplot(emm_spline_df, aes(x = nCondStims, y = emmean, color = phase_label)) +
    theme_light(base_size = 14) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = phase_label), alpha = 0.15, color = NA) +
    geom_rug(data = summaryNB_m6[summaryNB_m6$nCondStims <= max_dose_plot, ],
             aes(x = nCondStims, y = NULL), inherit.aes = FALSE, sides = "b", alpha = 0.3) +
    labs(x = "Number of Conditioning Stimuli",
         y = expression(paste("Predicted Magnitude (", mu, "V)")),
         color = "Phase", fill = "Phase",
         title = "Model 6 [EXPLORATORY]: Spline Dose-Response at 90 vs 270 deg") +
    scale_x_continuous(breaks = seq(1, max_dose_plot, by = 2))
  p_m6

  if(savePlot){
    ggsave(here("output_plots","betaStim_model6_spline_dose_response.png"), plot = p_m6,
           units = "in", width = 7, height = 4.5, dpi = 600)
    ggsave(here("output_plots","betaStim_model6_spline_dose_response.eps"), plot = p_m6,
           units = "in", width = 7, height = 4.5, dpi = 600, device = cairo_ps)
  }

  cat(sprintf("\nModel 6 AIC=%.1f BIC=%.1f\n", AIC(fit.spline), BIC(fit.spline)))

  # ========================================================================
  # Model 6-gf: Good-fit version — continuous spline dose x per-burst phase
  # ========================================================================
  # Restricts to trials where nGoodBeta > 0 on that channel (same filter as
  # Model 5a-gf). Uses per-burst burstCircMean for phase (not channel-level
  # average). Aggregates by (sid, channel, nCondStims) with circular mean
  # of per-burst sin/cos. This avoids the inflated df of the global Model 6
  # (1280 obs with many 1-2 trial cells).
  cat("\n=== Model 6-gf [EXPLORATORY]: Good-fit spline dose x per-burst phase ===\n")
  cat("NOTE: p-values are anti-conservative — no random dose slopes (singular at n=7 subjects).\n")

  # merge nCondStims + nGoodBeta + burstCircMean together
  precision_gf_m6 <- precision_combined[, c("probeSample", "channelEncoded",
    "nCondStims", "nGoodBeta", "burstCircMean")]
  precision_gf_m6$channelEncoded <- as.factor(precision_gf_m6$channelEncoded)

  data_m6_gf <- merge(data, precision_gf_m6,
    by.x = c("probeSample", "channel"),
    by.y = c("probeSample", "channelEncoded"),
    all.x = TRUE)

  # filter: non-baseline, good beta fit, valid nCondStims
  data_m6_gf <- data_m6_gf[!is.na(data_m6_gf$nGoodBeta) & data_m6_gf$nGoodBeta > 0 &
    !is.na(data_m6_gf$nCondStims) & data_m6_gf$nCondStims > 0 &
    data_m6_gf$numStims != "Base", ]

  # per-burst sin/cos from burstCircMean (actual delivered phase per burst)
  data_m6_gf$sin_phase_burst <- sin(data_m6_gf$burstCircMean * pi / 180)
  data_m6_gf$cos_phase_burst <- cos(data_m6_gf$burstCircMean * pi / 180)

  cat(sprintf("Good-fit trials: %d, nCondStims range [%d, %d]\n",
      nrow(data_m6_gf), min(data_m6_gf$nCondStims), max(data_m6_gf$nCondStims)))

  # summary: one median per (sid, channel, nCondStims)
  # Phase via circular mean of per-burst values: mean(sin), mean(cos)
  summaryNB_m6_gf <- ddply(data_m6_gf, .(sid, channel, nCondStims),
    summarize, magnitude = median(magnitude),
    sin_phase = mean(sin_phase_burst), cos_phase = mean(cos_phase_burst),
    betaLabels = first(betaLabels))

  # baseline covariate (unfiltered, same as other models)
  summaryNB_m6_gf <- merge(summaryNB_m6_gf, basePerChan, by = c("sid", "channel"))
  summaryNB_m6_gf$baselineMag_c <- summaryNB_m6_gf$baselineMag - mean(summaryNB_m6_gf$baselineMag)

  cat(sprintf("Model 6-gf summary: %d obs, %d subjects, %d channels, nCondStims range [%d, %d]\n",
      nrow(summaryNB_m6_gf), length(unique(summaryNB_m6_gf$sid)),
      length(unique(summaryNB_m6_gf$channel)),
      min(summaryNB_m6_gf$nCondStims), max(summaryNB_m6_gf$nCondStims)))

  # fit: intercepts only (same rationale as Model 5a-gf — random slope
  # variance estimated at zero with filtered data)
  fit.spline.gf = lmerTest::lmer(
    magnitude ~ ns(nCondStims, df = 3) * (sin_phase + cos_phase) + betaLabels + baselineMag_c +
    (1 | sid) + (1 | channel),
    data = summaryNB_m6_gf,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

  cat("Singular:", isSingular(fit.spline.gf), "\n")
  print(summary(fit.spline.gf))
  cat("\nType III ANOVA (Satterthwaite):\n")
  print(anova(fit.spline.gf))

  # LRT: phase effect
  fit.spline.gf.ml <- update(fit.spline.gf, REML = FALSE)
  fit.spline.gf.nophase.ml <- lmerTest::lmer(
    magnitude ~ ns(nCondStims, df = 3) + betaLabels + baselineMag_c +
    (1 | sid) + (1 | channel),
    data = summaryNB_m6_gf, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

  cat("\nLRT (6-gf): full vs no-phase:\n")
  print(anova(fit.spline.gf.nophase.ml, fit.spline.gf.ml))

  # Wald: average phase modulation
  cat("\nWald test (6-gf): average phase modulation:\n")
  print(car::linearHypothesis(fit.spline.gf, c("sin_phase = 0", "cos_phase = 0")))

  # emmeans dose-response curve at 90 vs 270
  max_dose_gf <- min(20, max(summaryNB_m6_gf$nCondStims))
  spline_at_gf <- list(
    nCondStims = seq(1, max_dose_gf, by = 0.5),
    sin_phase = c(1, -1),
    cos_phase = 0,
    betaLabels = "0",
    baselineMag_c = 0)

  emm_spline_gf <- emmeans(fit.spline.gf, ~ nCondStims * sin_phase, at = spline_at_gf)
  emm_spline_gf_df <- as.data.frame(emm_spline_gf)
  emm_spline_gf_df$phase_label <- ifelse(emm_spline_gf_df$sin_phase > 0, "90 deg", "270 deg")

  p_m6_gf <- ggplot(emm_spline_gf_df, aes(x = nCondStims, y = emmean, color = phase_label)) +
    theme_light(base_size = 14) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = phase_label), alpha = 0.15, color = NA) +
    geom_rug(data = summaryNB_m6_gf[summaryNB_m6_gf$nCondStims <= max_dose_gf, ],
             aes(x = nCondStims, y = NULL), inherit.aes = FALSE, sides = "b", alpha = 0.3) +
    labs(x = "Number of Conditioning Stimuli",
         y = expression(paste("Predicted Magnitude (", mu, "V)")),
         color = "Phase", fill = "Phase",
         title = "Model 6-gf [EXPLORATORY]: Spline Dose-Response (good fits, per-burst phase)") +
    scale_x_continuous(breaks = seq(1, max_dose_gf, by = 2))
  p_m6_gf

  if(savePlot){
    ggsave(here("output_plots","betaStim_model6_gf_spline_dose_response.png"), plot = p_m6_gf,
           units = "in", width = 7, height = 4.5, dpi = 600)
    ggsave(here("output_plots","betaStim_model6_gf_spline_dose_response.eps"), plot = p_m6_gf,
           units = "in", width = 7, height = 4.5, dpi = 600, device = cairo_ps)
  }

  cat(sprintf("\nModel 6-gf AIC=%.1f BIC=%.1f\n", AIC(fit.spline.gf), BIC(fit.spline.gf)))
  cat(sprintf("Comparison: Model 6 AIC=%.1f (1280 obs) vs 6-gf AIC=%.1f (%d obs)\n",
      AIC(fit.spline), AIC(fit.spline.gf), nrow(summaryNB_m6_gf)))

} else {
  cat("No burst_phase_precision CSVs found — skipping Model 6.\n")
}

# ========================================================================
# Model 5a-gf: Secondary analysis — restrict to good beta fits
# ========================================================================
# Filter trials to those where at least one conditioning stim in the burst
# had a good beta-band sinusoidal fit (R^2 > 0.7, frequency 12-20 Hz) on
# ANY channel. If that yields too few subjects, fall back to beta reference
# channel data only (dropping betaLabels since all data is from one channel
# type per subject).
cat("\n=== Secondary analysis: good beta fit restriction ===\n")

if (exists("precision_combined")) {
  # precision_combined already loaded before Model 6

  # --- Merge per-channel nGoodBeta AND burstCircMean ---
  # burstCircMean: circular mean phase (degrees) from the specific conditioning
  # burst preceding each probe, computed by compute_burst_phase_precision.m.
  # Includes ONLY conditioning stims that passed R^2 > 0.7 AND frequency 12-20 Hz
  # (beta band). If a burst of 5 stims had 1 good fit, burstCircMean = that
  # single stim's phase. nGoodBeta counts how many stims passed the threshold.
  #
  # Dose labels ([1,2], [3,4], [5,inf)) reflect TOTAL conditioning stims
  # delivered in the burst, not the number with good fits. Dose = experimental
  # manipulation; good-fit filter = measurement quality. These are independent.
  precision_sub <- precision_combined[, c("probeSample", "channelEncoded",
    "nGoodBeta", "burstCircMean", "burstVecLength")]
  precision_sub$channelEncoded <- as.factor(precision_sub$channelEncoded)

  data_merged <- merge(data, precision_sub,
    by.x = c("probeSample", "channel"),
    by.y = c("probeSample", "channelEncoded"),
    all.x = TRUE)

  # Trial quality filters: minimum good-fit stims per burst and optional
  # per-burst vector length concentration
  dataGoodFit <- apply_min_filter(data_merged, "nGoodBeta", minGoodBetaPerBurst)
  dataGoodFit <- apply_min_filter(dataGoodFit, "burstVecLength",
    minBurstVecLength_gf, label = "Per-burst")
  dataGoodFit_NB <- dataGoodFit[dataGoodFit$numStims != "Base", ]
  nSubj_gf <- length(unique(dataGoodFit_NB$sid))

  cat(sprintf("Per-channel good-fit filter (nGoodBeta >= %d): %d non-baseline trials from %d subjects\n",
      minGoodBetaPerBurst, nrow(dataGoodFit_NB), nSubj_gf))

  # --- Fallback: if too few subjects, use beta reference channel only ---
  use_beta_only <- nSubj_gf < 2
  if (use_beta_only) {
    cat("Too few subjects with per-channel good fits; restricting to beta reference channel.\n")
    # Beta reference channels are marked by betaLabels == "1" in the CSV
    beta_chans <- unique(as.character(data$channel[data$betaLabels == "1"]))
    dataGoodFit <- data_merged[data_merged$channel %in% beta_chans, ]
    dataGoodFit$nGoodBeta <- NULL
    dataGoodFit$burstCircMean <- NULL

    precision_beta <- precision_combined[precision_combined$channelEncoded %in% as.numeric(beta_chans), ]
    precision_beta_probe <- precision_beta[, c("probeSample", "nGoodBeta",
      "burstCircMean", "burstVecLength")]
    names(precision_beta_probe)[2:4] <- c("nGoodBeta_beta", "burstCircMean_beta",
      "burstVecLength_beta")
    dataGoodFit <- merge(dataGoodFit, precision_beta_probe, by = "probeSample", all.x = TRUE)
    dataGoodFit <- dataGoodFit[!is.na(dataGoodFit$nGoodBeta_beta) &
      dataGoodFit$nGoodBeta_beta >= minGoodBetaPerBurst, ]
    dataGoodFit$burstCircMean <- dataGoodFit$burstCircMean_beta
    dataGoodFit$burstVecLength <- dataGoodFit$burstVecLength_beta
    if (minBurstVecLength_gf > 0) {
      dataGoodFit <- dataGoodFit[!is.na(dataGoodFit$burstVecLength) &
        dataGoodFit$burstVecLength >= minBurstVecLength_gf, ]
    }
    dataGoodFit_NB <- dataGoodFit[dataGoodFit$numStims != "Base", ]
    nSubj_gf <- length(unique(dataGoodFit_NB$sid))
    cat(sprintf("Beta-channel-only good-fit filter: %d non-baseline trials from %d subjects\n",
        nrow(dataGoodFit_NB), nSubj_gf))
  }

  cat(sprintf("Good-fit trials: %d of %d non-baseline trials (%.0f%%)\n",
      nrow(dataGoodFit_NB), nrow(dataNoBaseline),
      100 * nrow(dataGoodFit_NB) / nrow(dataNoBaseline)))

  # --- Compute per-trial sin/cos from burst-specific measured phase ---
  # burstCircMean is the actual delivered phase for each specific burst,
  # computed from only the good-fit conditioning stims (R^2 > 0.7, 12-20 Hz).
  # This replaces the channel-level average phase (phaseDeg) for the good-fit
  # analysis, giving a more precise per-trial phase estimate.
  dataGoodFit_NB$sin_phase_burst <- sin(dataGoodFit_NB$burstCircMean * pi / 180)
  dataGoodFit_NB$cos_phase_burst <- cos(dataGoodFit_NB$burstCircMean * pi / 180)

  # --- Aggregate to summary level ---
  # Group by (sid, numStims, channel). No phaseClass or phaseDeg_round needed
  # because the per-burst phase varies trial-to-trial within a channel.
  # Phase aggregation uses mean(sin) and mean(cos) — the standard circular
  # mean in Cartesian form (Fisher 1993). This gives the average delivered
  # phase direction across good-fit bursts within each cell.
  # EP magnitude uses median, consistent with all other summary models.
  summaryNB_gf <- ddply(dataGoodFit_NB, .(sid, numStims, channel),
    summarize, magnitude = median(magnitude),
    sin_phase = mean(sin_phase_burst), cos_phase = mean(cos_phase_burst),
    betaLabels = first(betaLabels))

  # baselines from UNFILTERED data — baseline probes have no conditioning stims
  # so nGoodBeta is always 0; filtering would remove all baselines
  # basePerChan already computed (~line 408): one median per (sid, channel)
  summaryNB_gf <- merge(summaryNB_gf, basePerChan, by = c("sid", "channel"))
  summaryNB_gf$baselineMag_c <- summaryNB_gf$baselineMag - mean(summaryNB_gf$baselineMag)
  summaryNB_gf$doseNum <- as.numeric(factor(summaryNB_gf$numStims,
    levels = c("[1,2]","[3,4]","[5,inf)"))) - 1

  # ordinal dose coding
  summaryNB_gf$numStims_ord <- ordered(summaryNB_gf$numStims,
    levels = c("[1,2]", "[3,4]", "[5,inf)"))
  summaryNB_gf$dose_linpoly <- poly_lin[as.numeric(summaryNB_gf$numStims_ord)]

  nSubj_gf <- length(unique(summaryNB_gf$sid))
  nChan_gf <- length(unique(summaryNB_gf$channel))
  cat(sprintf("Good-fit summary: %d obs, %d subjects, %d channels\n",
      nrow(summaryNB_gf), nSubj_gf, nChan_gf))

  # refit Model 5a — drop betaLabels if beta-only fallback (all channels are beta)
  # 5a-gf uses per-burst phase which collapses within each (sid, numStims,
  # channel) cell — ~96 obs total. Random dose slope is singular at this
  # sample size (verified 2026-04-10 after the random-slope investigation).
  # Intercepts-only. Note: the p-value on dose without a random slope is
  # anti-conservative because within-subject dose observations are treated as
  # more independent than they are. For the stable, fully-grouped version
  # (5a-gf2, 153 obs) we use the full random-slope structure.
  gf_formula <- if (use_beta_only) {
    magnitude ~ numStims_ord * (sin_phase + cos_phase) + baselineMag_c +
      (1 | sid) + (1 | channel)
  } else {
    magnitude ~ numStims_ord * (sin_phase + cos_phase) + betaLabels + baselineMag_c +
      (1 | sid) + (1 | channel)
  }

  fit.sincos.ordinal.gf = lmerTest::lmer(gf_formula,
    data = summaryNB_gf,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

  cat("\nModel 5a-gf (good beta fits only, per-burst phase):\n")
  cat("Singular:", isSingular(fit.sincos.ordinal.gf), "\n")
  print(summary(fit.sincos.ordinal.gf))
  cat("\nType III ANOVA (Satterthwaite):\n")
  print(anova(fit.sincos.ordinal.gf))

  # joint tests (ML)
  fit.sincos.ordinal.gf.ml = update(fit.sincos.ordinal.gf, REML = FALSE)
  nophase_gf_formula <- if (use_beta_only) {
    magnitude ~ numStims_ord + baselineMag_c +
      (1 | sid) + (1 | channel)
  } else {
    magnitude ~ numStims_ord + betaLabels + baselineMag_c +
      (1 | sid) + (1 | channel)
  }
  fit.sincos.nophase.gf.ml = lmerTest::lmer(nophase_gf_formula,
    data = summaryNB_gf, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

  cat("\nLRT (good-fit): full vs no-phase (omnibus phase test):\n")
  print(anova(fit.sincos.nophase.gf.ml, fit.sincos.ordinal.gf.ml))

  cat("\nWald test (good-fit): average phase modulation:\n")
  print(car::linearHypothesis(fit.sincos.ordinal.gf,
    c("sin_phase = 0", "cos_phase = 0")))

  cat(sprintf("\nComparison: Model 5a AIC=%.1f vs 5a-gf AIC=%.1f\n",
      AIC(fit.sincos.ordinal), AIC(fit.sincos.ordinal.gf)))

  # --- Phase-response curve for good-fit model ---
  at_base_gf <- list(betaLabels = "0", baselineMag_c = 0)
  if (use_beta_only) at_base_gf <- list(baselineMag_c = 0)

  emm_curve_gf <- lapply(phase_vals, function(ph) {
    at_args <- c(list(sin_phase = sin(ph * pi / 180),
                      cos_phase = cos(ph * pi / 180)), at_base_gf)
    em <- emmeans(fit.sincos.ordinal.gf, ~ numStims_ord, at = at_args)
    df <- as.data.frame(em)
    df$phase_deg <- ph
    df
  })
  emm_curve_gf_df <- do.call(rbind, emm_curve_gf)

  p_5a_gf <- ggplot(emm_curve_gf_df, aes(x = phase_deg, y = emmean, color = numStims_ord)) +
    theme_light(base_size = 14) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = numStims_ord), alpha = 0.15, color = NA) +
    labs(x = "Delivered Phase (degrees)",
         y = expression(paste("Predicted Magnitude (", mu, "V)")),
         color = "Dose", fill = "Dose",
         title = "Model 5a-gf: Phase-Response Curve (good beta fits only)") +
    scale_x_continuous(breaks = seq(0, 315, by = 45))
  p_5a_gf

  if(savePlot){
    ggsave(here("output_plots","betaStim_model5a_gf_phase_curve.png"), plot = p_5a_gf,
           units = "in", width = 7, height = 4.5, dpi = 600)
    ggsave(here("output_plots","betaStim_model5a_gf_phase_curve.eps"), plot = p_5a_gf,
           units = "in", width = 7, height = 4.5, dpi = 600, device = cairo_ps)
  }

  # ========================================================================
  # Model 5a-gf2: Good-fit trials, CHANNEL-LEVEL phase predictor
  # ========================================================================
  cat("\n=== Model 5a-gf2: Good-fit trials, channel-level phase ===\n")

  dataGoodFit_NB$phaseDeg_round <- round(dataGoodFit_NB$phaseDeg, 1)

  # Channel-level phase-quality filter for 5a-gf2 (5a-gf uses per-burst
  # phase and already gets a per-burst filter upstream via minBurstVecLength_gf)
  data_gf2 <- apply_min_filter(dataGoodFit_NB, "phaseVecLength",
    minPhaseVecLength_gf2, label = "Model 5a-gf2")

  summaryNB_gf2 <- ddply(data_gf2, .(sid, phaseDeg_round, numStims, channel),
    summarize, magnitude = median(magnitude),
    sin_phase = first(sin_phase), cos_phase = first(cos_phase),
    betaLabels = first(betaLabels))

  summaryNB_gf2 <- merge(summaryNB_gf2, basePerChan, by = c("sid", "channel"))
  summaryNB_gf2$baselineMag_c <- summaryNB_gf2$baselineMag - mean(summaryNB_gf2$baselineMag)
  summaryNB_gf2$doseNum <- as.numeric(factor(summaryNB_gf2$numStims,
    levels = c("[1,2]","[3,4]","[5,inf)"))) - 1
  summaryNB_gf2$numStims_ord <- ordered(summaryNB_gf2$numStims,
    levels = c("[1,2]", "[3,4]", "[5,inf)"))
  summaryNB_gf2$dose_linpoly <- poly_lin[as.numeric(summaryNB_gf2$numStims_ord)]

  nSubj_gf2 <- length(unique(summaryNB_gf2$sid))
  nChan_gf2 <- length(unique(summaryNB_gf2$channel))
  cat(sprintf("Good-fit summary (channel-level phase): %d obs, %d subjects, %d channels\n",
      nrow(summaryNB_gf2), nSubj_gf2, nChan_gf2))

  # Match Model 5a's random structure: random dose slope per subject.
  # Post 702d24 fix (2026-04-10) this is non-singular for 5a-gf2 as well.
  fit.sincos.ordinal.gf2 = lmerTest::lmer(
    magnitude ~ numStims_ord * (sin_phase + cos_phase) + betaLabels + baselineMag_c +
    (1 + dose_linpoly | sid) + (1 | channel),
    data = summaryNB_gf2,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

  cat("Singular:", isSingular(fit.sincos.ordinal.gf2), "\n")
  print(summary(fit.sincos.ordinal.gf2))
  cat("\nType III ANOVA (Satterthwaite):\n")
  print(anova(fit.sincos.ordinal.gf2))

  fit.sincos.ordinal.gf2.ml = update(fit.sincos.ordinal.gf2, REML = FALSE)
  fit.sincos.nophase.gf2.ml = lmerTest::lmer(
    magnitude ~ numStims_ord + betaLabels + baselineMag_c +
    (1 + dose_linpoly | sid) + (1 | channel),
    data = summaryNB_gf2, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

  cat("\nLRT (gf2): full vs no-phase (omnibus phase test):\n")
  print(anova(fit.sincos.nophase.gf2.ml, fit.sincos.ordinal.gf2.ml))

  cat("\nWald test (gf2): average phase modulation (sin=cos=0 at mean dose):\n")
  print(car::linearHypothesis(fit.sincos.ordinal.gf2,
    c("sin_phase = 0", "cos_phase = 0")))

  # All pairwise dose contrasts use Tukey adjustment (default for pairs()).
  cat("\n--- emmeans (gf2): dose contrasts at phase 90 and 270 (Tukey-adjusted) ---\n")
  # betaLabels omitted so emmeans marginalizes equally across non-beta / beta
  # channel types. See comment at emm_5a_90 above.
  emm_gf2_90 <- emmeans(fit.sincos.ordinal.gf2, ~ numStims_ord,
    at = list(sin_phase = 1, cos_phase = 0, baselineMag_c = 0))
  emm_gf2_270 <- emmeans(fit.sincos.ordinal.gf2, ~ numStims_ord,
    at = list(sin_phase = -1, cos_phase = 0, baselineMag_c = 0))
  cat("At phase=90:\n"); print(as.data.frame(emm_gf2_90))
  cat("At phase=270:\n"); print(as.data.frame(emm_gf2_270))

  # Phase contrast (90 vs 270) at each dose level — cos=0 for both, no factorial issue.
  # Parallels emm_5a_phase at line 1160 so the docx and console include both primary
  # models' 90-vs-270 phase contrasts for direct side-by-side interpretation.
  cat("\n--- emmeans (gf2): phase 90 vs 270 at each dose ---\n")
  emm_gf2_phase <- emmeans(fit.sincos.ordinal.gf2, ~ numStims_ord * sin_phase,
    at = list(sin_phase = c(1, -1), cos_phase = 0, baselineMag_c = 0))
  print(contrast(emm_gf2_phase, method = "pairwise", by = "numStims_ord"))

  # --- Variance components + total-SD effect sizes for Model 5a-gf2 ---
  # Mirrors the 5a computation at lines 1199-1231. Same dose-polynomial
  # convention (contr.poly(3)[,1]) so E[x^2] = 1/3.
  vc_gf2 <- VarCorr(fit.sincos.ordinal.gf2)
  var_sid_int_gf2   <- attr(vc_gf2$sid,     "stddev")["(Intercept)"]^2
  var_sid_slope_gf2 <- attr(vc_gf2$sid,     "stddev")["dose_linpoly"]^2
  var_channel_gf2   <- attr(vc_gf2$channel, "stddev")["(Intercept)"]^2
  var_resid_gf2     <- sigma(fit.sincos.ordinal.gf2)^2
  total_sd_gf2 <- sqrt(var_sid_int_gf2 + ex2 * var_sid_slope_gf2 +
                       var_channel_gf2 + var_resid_gf2)
  resid_sd_gf2 <- sigma(fit.sincos.ordinal.gf2)

  cat(sprintf("\n5a-gf2 variance components: total SD = %.1f uV, residual SD = %.1f uV\n",
      total_sd_gf2, resid_sd_gf2))

  es_gf2_90  <- compute_effect_sizes(emm_gf2_90,  total_sd_gf2, resid_sd_gf2,
                                      "5a-gf2 dose at phase=90")
  es_gf2_270 <- compute_effect_sizes(emm_gf2_270, total_sd_gf2, resid_sd_gf2,
                                      "5a-gf2 dose at phase=270")

  # --- Average Marginal Effect of dose (marginaleffects), Model 5a-gf2 ---
  # Parallels the Model 5a AME block — marginalizes dose contrasts over the
  # observed (sin_phase, cos_phase, baselineMag_c, betaLabels) distribution
  # in summaryNB_gf2. Single-step max-t multiplicity correction (see rationale
  # at ame_5a_dose).
  cat("\n--- AME: Model 5a-gf2 dose pairwise contrasts, marginalized over observed covariates ---\n")
  ame_gf2_raw <- avg_comparisons(fit.sincos.ordinal.gf2,
    variables = list(numStims_ord = "pairwise"))
  ame_gf2_adj <- as.data.frame(hypotheses(ame_gf2_raw, multcomp = "single-step"))
  ame_gf2_adj$contrast <- as.data.frame(ame_gf2_raw)$contrast
  ame_gf2_dose <- ame_gf2_adj[, c("term", "contrast",
    setdiff(names(ame_gf2_adj), c("term", "contrast")))]
  print(ame_gf2_dose)

  # betaLabels marginalized using observed proportions in summaryNB_gf2 —
  # matches 5a (see comment at emm_curve above) and the AME table semantics.
  emm_curve_gf2 <- lapply(phase_vals, function(ph) {
    em <- emmeans(fit.sincos.ordinal.gf2, ~ numStims_ord,
      at = list(sin_phase = sin(ph*pi/180), cos_phase = cos(ph*pi/180),
                baselineMag_c = 0),
      weights = "proportional")
    df <- as.data.frame(em)
    df$phase_deg <- ph
    df
  })
  emm_curve_gf2_df <- do.call(rbind, emm_curve_gf2)

  lbl_gf2 <- phase_curve_label("Model 5a-gf2", minPhaseVecLength_gf2, nrow(summaryNB_gf2))
  p_5a_gf2 <- ggplot(emm_curve_gf2_df, aes(x = phase_deg, y = emmean, color = numStims_ord)) +
    theme_light(base_size = 14) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = numStims_ord), alpha = 0.15, color = NA) +
    labs(x = "Delivered Phase (degrees)",
         y = expression(paste("Predicted Magnitude (", mu, "V)")),
         color = "Dose", fill = "Dose",
         title = lbl_gf2$title,
         subtitle = "emmeans +/- 95% CI; baseline at grand mean; betaLabels marginalized over observed proportions") +
    scale_x_continuous(breaks = seq(0, 315, by = 45))
  p_5a_gf2

  if(savePlot){
    ggsave(here("output_plots", paste0(lbl_gf2$fname, ".png")), plot = p_5a_gf2,
           units = "in", width = 7, height = 4.5, dpi = 600)
    ggsave(here("output_plots", paste0(lbl_gf2$fname, ".eps")), plot = p_5a_gf2,
           units = "in", width = 7, height = 4.5, dpi = 600, device = cairo_ps)
  }

  # --- Dose-effect-vs-phase curve (marginaleffects), Model 5a-gf2 ---
  # Parallel to p_5a_effect. At each phase angle, override sin_phase /
  # cos_phase on summaryNB_gf2 and average the pairwise dose contrast
  # across the observed (betaLabels, baselineMag_c, sid, channel)
  # distribution — same marginalization as the AME numeric block.
  cmp_gf2_curve_df <- do.call(rbind, lapply(phase_vals_curve, function(ph) {
    nd <- summaryNB_gf2
    nd$sin_phase <- sin(ph * pi / 180)
    nd$cos_phase <- cos(ph * pi / 180)
    df <- as.data.frame(avg_comparisons(fit.sincos.ordinal.gf2,
      variables = list(numStims_ord = "pairwise"),
      newdata = nd, re.form = NA))
    df$phase_deg <- ph
    df
  }))

  p_5a_gf2_effect <- ggplot(cmp_gf2_curve_df,
      aes(x = phase_deg, y = estimate, color = contrast, fill = contrast)) +
    theme_light(base_size = 14) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    labs(x = "Delivered Phase (degrees)",
         y = expression(paste("Dose Effect (", mu, "V)")),
         color = "Pairwise contrast", fill = "Pairwise contrast",
         title = paste0(lbl_gf2$title, " - Dose Effect vs Phase"),
         subtitle = "avg_comparisons() at each phase, averaged over observed betaLabels and baseline; pointwise 95% CI") +
    scale_x_continuous(breaks = seq(0, 360, by = 45))
  p_5a_gf2_effect

  if(savePlot){
    ggsave(here("output_plots", paste0(lbl_gf2$fname, "_dose_effect_vs_phase.png")),
           plot = p_5a_gf2_effect, units = "in", width = 7, height = 4.5, dpi = 600)
    ggsave(here("output_plots", paste0(lbl_gf2$fname, "_dose_effect_vs_phase.eps")),
           plot = p_5a_gf2_effect, units = "in", width = 7, height = 4.5, dpi = 600,
           device = cairo_ps)
  }

  # --- Supplementary: phase-response curves faceted by betaLabels (gf2) ---
  emm_curve_beta_gf2 <- do.call(rbind, lapply(c("0", "1"), function(bl) {
    do.call(rbind, lapply(phase_vals, function(ph) {
      em <- emmeans(fit.sincos.ordinal.gf2, ~ numStims_ord,
        at = list(sin_phase = sin(ph * pi / 180), cos_phase = cos(ph * pi / 180),
                  betaLabels = bl, baselineMag_c = 0))
      df <- as.data.frame(em)
      df$phase_deg <- ph
      df$betaLabels <- bl
      df
    }))
  }))
  emm_curve_beta_gf2$channel_type <- factor(
    ifelse(emm_curve_beta_gf2$betaLabels == "0", "Non-beta (EP) channels",
           "Beta trigger channels"),
    levels = c("Non-beta (EP) channels", "Beta trigger channels"))

  p_5a_gf2_beta <- ggplot(emm_curve_beta_gf2,
                          aes(x = phase_deg, y = emmean, color = numStims_ord)) +
    theme_light(base_size = 14) +
    facet_wrap(~ channel_type, nrow = 1) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = numStims_ord),
                alpha = 0.15, color = NA) +
    labs(x = "Delivered Phase (degrees)",
         y = expression(paste("Predicted Magnitude (", mu, "V)")),
         color = "Dose", fill = "Dose",
         title = paste(lbl_gf2$title, "- Beta vs Non-beta Channels"),
         subtitle = "emmeans +/- 95% CI; baseline at grand mean. betaLabels fixed effect p = 0.70 (ns)") +
    scale_x_continuous(breaks = seq(0, 315, by = 90)) +
    theme(strip.text = element_text(face = "bold"))

  if (savePlot) {
    ggsave(here("output_plots", paste0(lbl_gf2$fname, "_beta_vs_nonbeta.png")),
           plot = p_5a_gf2_beta, units = "in", width = 11, height = 4.5, dpi = 600)
    ggsave(here("output_plots", paste0(lbl_gf2$fname, "_beta_vs_nonbeta.eps")),
           plot = p_5a_gf2_beta, units = "in", width = 11, height = 4.5, dpi = 600,
           device = cairo_ps)
  }

  cat(sprintf("\nComparison: 5a AIC=%.1f, 5a-gf (per-burst) AIC=%.1f, 5a-gf2 (channel) AIC=%.1f\n",
      AIC(fit.sincos.ordinal), AIC(fit.sincos.ordinal.gf), AIC(fit.sincos.ordinal.gf2)))

  # ========================================================================
  # Raw summary data: percent diff from baseline per cell, colored by
  # 8 phase bins (45° each, wrapped so 0/90/180/270 fall at bin CENTERS).
  # Two variants built from the same template:
  #   - summaryNB_m5  → Model 5a  (all trials, phaseVecLength ≥ minPhaseVecLength_5a)
  #   - summaryNB_gf2 → Model 5a-gf2 (good-fit trials, phaseVecLength ≥ minPhaseVecLength_gf2)
  # x-axis: same 3-level dose binning used by the LME models.
  # ========================================================================

  # Helper: 8 wedges of 45° each. Breaks at [−22.5, 22.5, 67.5, ..., 337.5]
  # wrap back to the "0" label via cut()'s label trick (first and last
  # intervals both get "0" so phases near 0° and near 360° merge).
  bin8phase_fn <- function(phase_deg) {
    b <- cut(phase_deg %% 360,
             breaks = c(-1, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 361),
             labels = c("0", "45", "90", "135", "180", "225", "270", "315", "0"))
    factor(b, levels = c("0", "45", "90", "135", "180", "225", "270", "315"))
  }

  # ------------------------------------------------------------------------
  # Compute percentDiffBase + phaseBin on both summary tables
  # ------------------------------------------------------------------------
  summaryNB_gf2$percentDiffBase <- 100 *
    (summaryNB_gf2$magnitude - summaryNB_gf2$baselineMag) / summaryNB_gf2$baselineMag
  summaryNB_gf2$phaseBin <- bin8phase_fn(summaryNB_gf2$phaseDeg_round)

  summaryNB_m5$percentDiffBase <- 100 *
    (summaryNB_m5$magnitude - summaryNB_m5$baselineMag) / summaryNB_m5$baselineMag
  summaryNB_m5$phaseBin <- bin8phase_fn(summaryNB_m5$phaseDeg_round)

  # Dodge dots by phase bin so each phase has its own sub-column within each
  # dose bin. To give the dose groups visible breathing room (and make the
  # 8-bin phase sub-columns readable), we manually place the three dose
  # groups at widely-spaced continuous x positions (1, 3, 5) rather than at
  # the default discrete positions (1, 2, 3). Each group spans dodge_width,
  # leaving a clear gap between groups with a dashed vertical separator.
  make_raw_plot <- function(df, model_label, r_thresh, n_subj, n_chan) {
    df$numStims_x <- c(1, 3, 5)[as.numeric(df$numStims_ord)]
    dodge_width <- 1.4
    ggplot(df,
      aes(x = numStims_x, y = percentDiffBase, color = phaseBin,
          group = phaseBin)) +
      theme_light(base_size = 14) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = c(2, 4), linetype = "dashed",
                 color = "grey70", linewidth = 0.5) +
      geom_point(alpha = 0.5, size = 2,
                 position = position_jitterdodge(jitter.width = 0.18,
                                                 dodge.width = dodge_width)) +
      stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
                   geom = "errorbar", width = 0.35, linewidth = 0.8,
                   position = position_dodge(width = dodge_width)) +
      stat_summary(fun = median, geom = "point", size = 3, shape = 23,
                   fill = "white", stroke = 1,
                   position = position_dodge(width = dodge_width)) +
      labs(x = "Number of Conditioning Stimuli",
           y = "Percent Difference from Baseline",
           color = "Delivered phase\n(deg, 8 bins × 45°)",
           title = sprintf("%s: Dose × Phase (8 bins)", model_label),
           subtitle = sprintf("%d cells from %d subjects × %d channels × phase conditions (phaseVecLength ≥ %.2f)",
                              nrow(df), n_subj, n_chan, r_thresh)) +
      scale_x_continuous(breaks = c(1, 3, 5),
                         labels = c("[1,2]", "[3,4]", "[5,inf)"),
                         limits = c(0.1, 5.9)) +
      scale_color_viridis_d(option = "turbo", end = 0.95)
  }

  nSubj_m5  <- length(unique(summaryNB_m5$sid))
  nChan_m5  <- length(unique(summaryNB_m5$channel))

  p_5a_raw      <- make_raw_plot(summaryNB_m5,  "Model 5a",     minPhaseVecLength_5a,
                                 nSubj_m5, nChan_m5)
  p_5a_gf2_raw  <- make_raw_plot(summaryNB_gf2, "Model 5a-gf2", minPhaseVecLength_gf2,
                                 nSubj_gf2, nChan_gf2)

  if (savePlot) {
    ggsave(here("output_plots","betaStim_model5a_raw_8phase.png"), plot = p_5a_raw,
           units = "in", width = 9, height = 5, dpi = 600)
    ggsave(here("output_plots","betaStim_model5a_raw_8phase.eps"), plot = p_5a_raw,
           units = "in", width = 9, height = 5, dpi = 600, device = cairo_ps)
    ggsave(here("output_plots","betaStim_model5a_gf2_raw_8phase.png"), plot = p_5a_gf2_raw,
           units = "in", width = 9, height = 5, dpi = 600)
    ggsave(here("output_plots","betaStim_model5a_gf2_raw_8phase.eps"), plot = p_5a_gf2_raw,
           units = "in", width = 9, height = 5, dpi = 600, device = cairo_ps)
  }

} else {
  cat("No burst_phase_precision CSVs found — skipping secondary analysis.\n")
  cat("Run compute_burst_phase_precision.m for each subject first.\n")
}

# ========================================================================
# Phase quality sensitivity analysis: fit Model 5a and 5a-gf2 at a grid of
# minPhaseVecLength thresholds and report dose.L / sin / cos p-values.
# Useful for manuscript sensitivity tables.
# ========================================================================
phase_quality_sensitivity <- function() {
  thresholds <- c(0, 0.1, 0.2, 0.3, 0.4)
  results <- list()

  fit_one <- function(data_in, model_label, r_thresh) {
    d <- data_in
    if (r_thresh > 0) {
      d <- d[!is.na(d$phaseVecLength) & d$phaseVecLength >= r_thresh, ]
    }
    if (nrow(d) == 0) return(NULL)
    summ <- ddply(d, .(sid, phaseDeg_round, numStims, channel),
      summarize, magnitude = median(magnitude),
      sin_phase = first(sin_phase), cos_phase = first(cos_phase),
      betaLabels = first(betaLabels))
    summ <- merge(summ, basePerChan, by = c("sid", "channel"))
    summ$baselineMag_c <- summ$baselineMag - mean(summ$baselineMag)
    summ$numStims_ord <- ordered(summ$numStims,
      levels = c("[1,2]", "[3,4]", "[5,inf)"))
    summ <- summ[!is.na(summ$numStims_ord), ]
    summ$dose_linpoly <- poly_lin[as.numeric(summ$numStims_ord)]
    summ$betaLabels <- as.factor(summ$betaLabels)

    fit <- tryCatch(
      lmerTest::lmer(
        magnitude ~ numStims_ord * (sin_phase + cos_phase) +
          betaLabels + baselineMag_c +
          (1 + dose_linpoly | sid) + (1 | channel),
        data = summ,
        control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000))),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)

    coefs <- summary(fit)$coefficients
    list(
      model = model_label,
      r_thresh = r_thresh,
      n_obs = nrow(summ),
      n_ch = length(unique(summ$channel)),
      n_subj = length(unique(summ$sid)),
      singular = isSingular(fit),
      dose_L_est = coefs["numStims_ord.L", "Estimate"],
      dose_L_p   = coefs["numStims_ord.L", "Pr(>|t|)"],
      sin_est    = coefs["sin_phase", "Estimate"],
      sin_p      = coefs["sin_phase", "Pr(>|t|)"],
      cos_est    = coefs["cos_phase", "Estimate"],
      cos_p      = coefs["cos_phase", "Pr(>|t|)"])
  }

  cat("\n\n=== PHASE QUALITY SENSITIVITY ANALYSIS ===\n")

  # Model 5a (all trials)
  for (th in thresholds) {
    r <- fit_one(dataNoBaseline, "5a", th)
    if (!is.null(r)) results[[length(results) + 1]] <- r
  }

  # Model 5a-gf2 (good-fit trials, channel-level phase)
  if (exists("dataGoodFit_NB")) {
    for (th in thresholds) {
      r <- fit_one(dataGoodFit_NB, "5a-gf2", th)
      if (!is.null(r)) results[[length(results) + 1]] <- r
    }
  }

  df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      Model      = x$model,
      `r_thresh` = x$r_thresh,
      N          = x$n_obs,
      `N_chan`   = x$n_ch,
      `N_subj`   = x$n_subj,
      Singular   = x$singular,
      dose_L_est = round(x$dose_L_est, 2),
      dose_L_p   = round(x$dose_L_p, 4),
      sin_est    = round(x$sin_est, 2),
      sin_p      = round(x$sin_p, 4),
      cos_est    = round(x$cos_est, 2),
      cos_p      = round(x$cos_p, 4),
      check.names = FALSE, stringsAsFactors = FALSE)
  }))
  cat("\n")
  print(df, row.names = FALSE)

  # Save to CSV for manuscript
  out_csv <- here("output_plots", "betaStim_phase_quality_sensitivity.csv")
  write.csv(df, out_csv, row.names = FALSE)
  cat(sprintf("\nSaved sensitivity table: %s\n", out_csv))

  invisible(df)
}
if (runPhaseSensitivity) {
  phase_quality_sensitivity_df <- phase_quality_sensitivity()
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
# NOTE: This model is kept for reference/comparison. Models 5a / 5a-gf2
# (continuous circular phase, summary-level) are the primary reported
# models — they use continuous phase (sin/cos), proper random dose slopes,
# and avoid the singularity issues of the trial-level nested-condition fit.
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

if (showTabModel) tab_model(fit.nested.condition)

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
  residual_diagnostics(fit.numeric, "3e: Numeric"),
  residual_diagnostics(fit.sincos.ordinal, "5a: Sincos ordinal"),
  residual_diagnostics(fit.sincos.numeric, "5b: Sincos numeric"),
  residual_diagnostics(fit.sincos.categ, "5c: Sincos categorical")
)
if (exists("fit.sincos.ordinal.gf")) {
  diag_list <- c(diag_list,
    list(residual_diagnostics(fit.sincos.ordinal.gf, "5a-gf: good-fit per-burst")))
}
if (exists("fit.sincos.ordinal.gf2")) {
  diag_list <- c(diag_list,
    list(residual_diagnostics(fit.sincos.ordinal.gf2, "5a-gf2: good-fit channel")))
}
diag_df <- do.call(rbind, diag_list)
print(diag_df)

# --- ggplot QQ plots and residuals-vs-fitted for all summary models ---
if (savePlot) {
  summary_models <- list(
    list(fit = fit.absDiff, name = "3a_absDiff"),
    list(fit = fit.modelD, name = "3b_baseline"),
    list(fit = fit.ancova, name = "3c_ANCOVA"),
    list(fit = fit.ordinal_lmer, name = "3d_ordinal"),
    list(fit = fit.numeric, name = "3e_numeric"),
    list(fit = fit.sincos.ordinal, name = "5a_sincos_ordinal"),
    list(fit = fit.sincos.numeric, name = "5b_sincos_numeric"),
    list(fit = fit.sincos.categ, name = "5c_sincos_categorical")
  )
  if (exists("fit.sincos.ordinal.gf")) {
    summary_models <- c(summary_models,
      list(list(fit = fit.sincos.ordinal.gf, name = "5a_gf_sincos")))
  }
  if (exists("fit.sincos.ordinal.gf2")) {
    summary_models <- c(summary_models,
      list(list(fit = fit.sincos.ordinal.gf2, name = "5a_gf2_sincos")))
  }

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
  # Subject characteristics: baseline CEP magnitude per subject
  # ------------------------------------------------------------------
  if (exists("subject_cep_summary")) {
    scs_out <- subject_cep_summary
    names(scs_out) <- c("Subject", "N trials", "N channels",
                        "Median (uV)", "MAD (uV)",
                        "Q25 (uV)", "Q75 (uV)")
    doc <- body_add_par(doc, "Subject Characteristics: Baseline CEP Magnitude",
                        style = "heading 2")
    ft <- flextable(scs_out) |> autofit() |>
      set_caption("Per-subject baseline evoked-potential magnitude (Base probe trials only, pooled across channels). Subjects labeled as 'Subject N (sid)' where N is the manuscript numeric ID. Median and MAD are robust statistics; Q25/Q75 give the interquartile range. Reflects trials post 25-1500 uV filter used in all downstream models.")
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")
  }

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
  # Base comparison rows (always present)
  models_base <- c("3a: absDiff (intercepts, sensitivity)",
                   "3b: Magnitude + baseline (sensitivity)",
                   "3c: ANCOVA (binary phaseClass, sensitivity)",
                   "3d: Ordinal dose (sensitivity)",
                   "3e: Numeric dose (sensitivity)",
                   "5a: Sin/cos ordinal (ANCOVA) -- PRIMARY",
                   "5b: Sin/cos numeric",
                   "5c: Sin/cos categorical")
  n_base   <- c(nrow(summaryNB), nrow(summaryAll), nrow(summaryNB_ancova),
                nrow(summaryNB_ancova), nrow(summaryNB_ancova),
                nrow(summaryNB_m5), nrow(summaryNB_m5), nrow(summaryNB_m5))
  aic_base <- c(AIC(fit.absDiff), AIC(fit.modelD), AIC(fit.ancova),
                AIC(fit.ordinal_lmer), AIC(fit.numeric),
                AIC(fit.sincos.ordinal), AIC(fit.sincos.numeric), AIC(fit.sincos.categ))
  bic_base <- c(BIC(fit.absDiff), BIC(fit.modelD), BIC(fit.ancova),
                BIC(fit.ordinal_lmer), BIC(fit.numeric),
                BIC(fit.sincos.ordinal), BIC(fit.sincos.numeric), BIC(fit.sincos.categ))
  sing_base <- c(isSingular(fit.absDiff), isSingular(fit.modelD),
                 isSingular(fit.ancova), isSingular(fit.ordinal_lmer),
                 isSingular(fit.numeric),
                 isSingular(fit.sincos.ordinal), isSingular(fit.sincos.numeric),
                 isSingular(fit.sincos.categ))

  # 5a-gf2 added when precision CSVs were present (gated upstream)
  if (exists("fit.sincos.ordinal.gf2")) {
    models_base <- c(models_base,
                     "5a-gf2: Sin/cos ordinal, good-fit + channel phase -- PRIMARY companion")
    n_base   <- c(n_base, nrow(summaryNB_gf2))
    aic_base <- c(aic_base, AIC(fit.sincos.ordinal.gf2))
    bic_base <- c(bic_base, BIC(fit.sincos.ordinal.gf2))
    sing_base <- c(sing_base, isSingular(fit.sincos.ordinal.gf2))
  }

  comparison_df <- data.frame(
    Model    = models_base,
    N        = n_base,
    AIC      = round(aic_base, 1),
    BIC      = round(bic_base, 1),
    Singular = sing_base,
    stringsAsFactors = FALSE)
  doc <- body_add_par(doc, "Table: Model Comparison", style = "heading 2")
  ft <- flextable(comparison_df) |> autofit() |>
    set_caption("Summary-level model comparison. Models 5a and 5a-gf2 are PRIMARY (continuous circular phase); 3-series are sensitivity. All models use median per cell.")
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

  # ==================================================================
  # PRIMARY MODELS — 5a and 5a-gf2 (continuous circular phase)
  # ==================================================================
  # These are the primary inferential models for the manuscript. Model 5a
  # is the conservative primary (channel-level phaseVecLength >= 0.3, all
  # conditioned trials); Model 5a-gf2 adds the good-fit burst filter
  # (nGoodBeta >= 1) and relaxes phaseVecLength to 0.2 to stay non-singular.
  # See CLAUDE.md and statistical_audit.md Finding 23 for rationale.
  doc <- body_add_par(doc, "PRIMARY MODELS (continuous circular phase)",
                      style = "heading 1")
  doc <- body_add_par(doc,
    paste("Models 5a and 5a-gf2 are the primary inferential models for the",
          "manuscript. Both use continuous circular phase via",
          "sin(phaseDeg) + cos(phaseDeg), ordinal dose with a random linear",
          "slope per subject, baseline ANCOVA, and betaLabels as an additive",
          "covariate. 5a applies the channel-level phase-quality filter",
          "(phaseVecLength >= 0.3); 5a-gf2 additionally restricts to",
          "conditioned trials from bursts with confirmed beta (nGoodBeta >= 1)",
          "and relaxes the phase-quality threshold to 0.2. See CLAUDE.md and",
          "statistical_audit.md Finding 23 for rationale."),
    style = "Normal")

  doc <- add_model_tables(doc, fit.sincos.ordinal,
                          "Model 5a (PRIMARY — sin/cos ordinal, ANCOVA; phaseVecLength >= 0.3)")

  # Format an already-computed contrast / summary object (from contrast() or
  # pairs(..., by = ...)) into a docx table with CIs, t-ratios, and p-values.
  # Used for phase-at-each-dose contrasts where the contrast call already has
  # `by = "numStims_ord"` baked in and we don't want render_emm_contrasts' pairs()
  # step to re-pair the result.
  render_contrast_direct <- function(doc, contr_obj, heading_str, caption_str) {
    tbl <- as.data.frame(summary(contr_obj, infer = c(TRUE, TRUE)))
    tbl$estimate <- round(tbl$estimate, 3)
    tbl$SE       <- round(tbl$SE, 3)
    tbl$df       <- round(tbl$df, 1)
    tbl$lower.CL <- round(tbl$lower.CL, 3)
    tbl$upper.CL <- round(tbl$upper.CL, 3)
    if ("t.ratio" %in% names(tbl)) tbl$t.ratio <- round(tbl$t.ratio, 2)
    if ("z.ratio" %in% names(tbl)) tbl$z.ratio <- round(tbl$z.ratio, 2)
    tbl$p <- sapply(tbl$p.value, fmt_p)
    tbl$p.value <- NULL
    names(tbl)[names(tbl) == "lower.CL"] <- "CI lower"
    names(tbl)[names(tbl) == "upper.CL"] <- "CI upper"
    doc <- body_add_par(doc, heading_str, style = "heading 2")
    ft <- flextable(tbl) |> autofit() |> set_caption(caption_str)
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")
    doc
  }

  # AME (marginaleffects::avg_comparisons) output: different column names from
  # emmeans so it gets its own render helper. Marginalizes over observed
  # covariates rather than pinning them at a reference value.
  render_ame_contrasts <- function(doc, ame_obj, heading_str, caption_str) {
    tbl <- as.data.frame(ame_obj)
    keep <- intersect(c("contrast", "estimate", "std.error", "statistic",
                        "p.value", "conf.low", "conf.high"), names(tbl))
    tbl <- tbl[, keep, drop = FALSE]
    tbl$estimate   <- round(tbl$estimate,   3)
    tbl$std.error  <- round(tbl$std.error,  3)
    if ("statistic" %in% names(tbl))  tbl$statistic  <- round(tbl$statistic,  2)
    if ("conf.low"  %in% names(tbl))  tbl$conf.low   <- round(tbl$conf.low,   3)
    if ("conf.high" %in% names(tbl))  tbl$conf.high  <- round(tbl$conf.high,  3)
    if ("p.value"   %in% names(tbl))  tbl$p          <- sapply(tbl$p.value, fmt_p)
    tbl$p.value <- NULL
    names(tbl)[names(tbl) == "std.error"] <- "SE"
    names(tbl)[names(tbl) == "conf.low"]  <- "CI lower"
    names(tbl)[names(tbl) == "conf.high"] <- "CI upper"
    doc <- body_add_par(doc, heading_str, style = "heading 2")
    ft <- flextable(tbl) |> autofit() |> set_caption(caption_str)
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")
    doc
  }

  # EMM Dose Contrasts: pairs() with estimates, CIs, t-ratio, AND p-values.
  # infer = c(TRUE, TRUE) requests both confidence interval and hypothesis test.
  render_emm_contrasts <- function(doc, emm_obj, heading_str, caption_str) {
    tbl <- as.data.frame(summary(pairs(emm_obj), infer = c(TRUE, TRUE)))
    tbl$estimate <- round(tbl$estimate, 3)
    tbl$SE       <- round(tbl$SE, 3)
    tbl$df       <- round(tbl$df, 1)
    tbl$lower.CL <- round(tbl$lower.CL, 3)
    tbl$upper.CL <- round(tbl$upper.CL, 3)
    if ("t.ratio" %in% names(tbl)) tbl$t.ratio <- round(tbl$t.ratio, 2)
    if ("z.ratio" %in% names(tbl)) tbl$z.ratio <- round(tbl$z.ratio, 2)
    tbl$p <- sapply(tbl$p.value, fmt_p)
    tbl$p.value <- NULL
    names(tbl)[names(tbl) == "lower.CL"] <- "CI lower"
    names(tbl)[names(tbl) == "upper.CL"] <- "CI upper"
    doc <- body_add_par(doc, heading_str, style = "heading 2")
    ft <- flextable(tbl) |> autofit() |> set_caption(caption_str)
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")
    doc
  }

  doc <- render_emm_contrasts(doc, emm_5a_90,
    "EMM Dose Contrasts: Model 5a @ phase=90 deg",
    "Pairwise dose contrasts at phase=90 deg (Model 5a, Tukey-adjusted).")
  doc <- render_emm_contrasts(doc, emm_5a_270,
    "EMM Dose Contrasts: Model 5a @ phase=270 deg",
    "Pairwise dose contrasts at phase=270 deg (Model 5a, Tukey-adjusted).")

  # EMM Phase Contrast: 90 vs 270 at each dose (Model 5a). Single comparison
  # per dose, no family-wise adjustment applied.
  doc <- render_contrast_direct(doc,
    contrast(emm_5a_phase, method = "pairwise", by = "numStims_ord"),
    "EMM Phase Contrasts (90 vs 270): Model 5a",
    "Phase 90 vs 270 contrast (estimate = EMM at 90 minus EMM at 270, in uV) at each dose level, Model 5a. Single comparison per dose, unadjusted.")

  # AME Dose Contrasts: Model 5a, phase-marginalized.
  if (exists("ame_5a_dose")) {
    doc <- render_ame_contrasts(doc, ame_5a_dose,
      "AME Dose Contrasts (marginalized over observed phase): Model 5a",
      "Pairwise dose contrasts (Model 5a, single-step max-t FWER correction), marginalized over the observed joint distribution of sin_phase, cos_phase, baselineMag_c, and betaLabels. Complements the at-90-deg / at-270-deg EMM contrasts above, which pin phase at +/-1 and covariates at 0.")
  }

  # Effect sizes for Model 5a (already in the file further down, but repeat
  # here so primary-model effect sizes live together with primary tables)
  es_90_primary <- es_90[, c("contrast", "estimate", "d_total",
                              "d_total_lower", "d_total_upper", "d_conditional")]
  es_90_primary[, -1] <- round(es_90_primary[, -1], 3)
  names(es_90_primary) <- c("Contrast", "Estimate (uV)", "d_total",
                             "d_total lower", "d_total upper", "d_conditional")
  doc <- body_add_par(doc, "Effect Sizes: Model 5a Dose @ phase=90",
                      style = "heading 2")
  ft <- flextable(es_90_primary) |> autofit() |>
    set_caption(sprintf("Total-variance Cohen's d for dose contrasts at phase=90 (Model 5a). d_total SD = %.1f uV, d_conditional SD = %.1f uV.",
                        total_sd_5a, resid_sd_5a))
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  es_270_primary <- es_270[, c("contrast", "estimate", "d_total",
                                "d_total_lower", "d_total_upper", "d_conditional")]
  es_270_primary[, -1] <- round(es_270_primary[, -1], 3)
  names(es_270_primary) <- c("Contrast", "Estimate (uV)", "d_total",
                              "d_total lower", "d_total upper", "d_conditional")
  doc <- body_add_par(doc, "Effect Sizes: Model 5a Dose @ phase=270",
                      style = "heading 2")
  ft <- flextable(es_270_primary) |> autofit() |>
    set_caption(sprintf("Total-variance Cohen's d for dose contrasts at phase=270 (Model 5a). d_total SD = %.1f uV.",
                        total_sd_5a))
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # --- Model 5a-gf2 (good-fit companion primary) ---
  if (exists("fit.sincos.ordinal.gf2")) {
    doc <- add_model_tables(doc, fit.sincos.ordinal.gf2,
      "Model 5a-gf2 (PRIMARY companion — good-fit + channel-level phase, phaseVecLength >= 0.2)")

    if (exists("emm_gf2_90")) {
      doc <- render_emm_contrasts(doc, emm_gf2_90,
        "EMM Dose Contrasts: Model 5a-gf2 @ phase=90 deg",
        "Pairwise dose contrasts at phase=90 deg (Model 5a-gf2, Tukey-adjusted).")
      doc <- render_emm_contrasts(doc, emm_gf2_270,
        "EMM Dose Contrasts: Model 5a-gf2 @ phase=270 deg",
        "Pairwise dose contrasts at phase=270 deg (Model 5a-gf2, Tukey-adjusted).")
    }

    # EMM Phase Contrast: 90 vs 270 at each dose (Model 5a-gf2). Parallels 5a.
    if (exists("emm_gf2_phase")) {
      doc <- render_contrast_direct(doc,
        contrast(emm_gf2_phase, method = "pairwise", by = "numStims_ord"),
        "EMM Phase Contrasts (90 vs 270): Model 5a-gf2",
        "Phase 90 vs 270 contrast (estimate = EMM at 90 minus EMM at 270, in uV) at each dose level, Model 5a-gf2. Single comparison per dose, unadjusted.")
    }

    # AME Dose Contrasts: Model 5a-gf2, phase-marginalized.
    if (exists("ame_gf2_dose")) {
      doc <- render_ame_contrasts(doc, ame_gf2_dose,
        "AME Dose Contrasts (marginalized over observed phase): Model 5a-gf2",
        "Pairwise dose contrasts (Model 5a-gf2, single-step max-t FWER correction), marginalized over the observed joint distribution of sin_phase, cos_phase, baselineMag_c, and betaLabels. Complements the at-90-deg / at-270-deg EMM contrasts above, which pin phase at +/-1 and covariates at 0.")
    }

    if (exists("es_gf2_90")) {
      es_gf2_90_export <- es_gf2_90[, c("contrast", "estimate", "d_total",
                                          "d_total_lower", "d_total_upper", "d_conditional")]
      es_gf2_90_export[, -1] <- round(es_gf2_90_export[, -1], 3)
      names(es_gf2_90_export) <- c("Contrast", "Estimate (uV)", "d_total",
                                    "d_total lower", "d_total upper", "d_conditional")
      doc <- body_add_par(doc, "Effect Sizes: Model 5a-gf2 Dose @ phase=90",
                          style = "heading 2")
      ft <- flextable(es_gf2_90_export) |> autofit() |>
        set_caption(sprintf("Total-variance Cohen's d for dose contrasts at phase=90 (Model 5a-gf2). d_total SD = %.1f uV, d_conditional SD = %.1f uV.",
                            total_sd_gf2, resid_sd_gf2))
      doc <- body_add_flextable(doc, ft)
      doc <- body_add_par(doc, "")

      es_gf2_270_export <- es_gf2_270[, c("contrast", "estimate", "d_total",
                                            "d_total_lower", "d_total_upper", "d_conditional")]
      es_gf2_270_export[, -1] <- round(es_gf2_270_export[, -1], 3)
      names(es_gf2_270_export) <- c("Contrast", "Estimate (uV)", "d_total",
                                      "d_total lower", "d_total upper", "d_conditional")
      doc <- body_add_par(doc, "Effect Sizes: Model 5a-gf2 Dose @ phase=270",
                          style = "heading 2")
      ft <- flextable(es_gf2_270_export) |> autofit() |>
        set_caption(sprintf("Total-variance Cohen's d for dose contrasts at phase=270 (Model 5a-gf2). d_total SD = %.1f uV.",
                            total_sd_gf2))
      doc <- body_add_flextable(doc, ft)
      doc <- body_add_par(doc, "")
    }
  }

  # ==================================================================
  # SUPPORTING / SENSITIVITY MODELS — 3a, 3c, 3e (binary phaseClass)
  # ==================================================================
  # Retained as robustness checks. Binary 90/270 binning discards circular
  # information; primary inference is in Models 5a / 5a-gf2 above.
  doc <- body_add_par(doc, "SUPPORTING / SENSITIVITY MODELS (binary phaseClass)",
                      style = "heading 1")
  doc <- body_add_par(doc,
    paste("Models 3a, 3c, and 3e use binary phaseClass (90 / 270 bins).",
          "Retained as robustness checks — they converge on the dose",
          "effect seen in 5a/5a-gf2 but discard circular-phase information",
          "and conflate distinct measured phases at multi-phase channels.",
          "These are not the primary inferential models for the manuscript."),
    style = "Normal")

  doc <- add_model_tables(doc, fit.absDiff, "Model 3a (absDiff, intercepts only — sensitivity)")
  doc <- add_model_tables(doc, fit.ancova, "Model 3c (ANCOVA, binary phaseClass — sensitivity)")
  doc <- add_model_tables(doc, fit.numeric, "Model 3e (numeric dose, binary phaseClass — sensitivity)")

  # EMM contrasts: Models 3a and 3c. Uses the render_emm_contrasts helper
  # (defined in the PRIMARY MODELS block above) so 3-series tables include
  # p-values + CIs + t-ratios identically to the 5a/5a-gf2 tables.
  doc <- render_emm_contrasts(doc, emm_3a_dose,
    "EMM Dose Contrasts: Model 3a",
    "Pairwise dose contrasts within each phase class (Model 3a, Tukey-adjusted).")
  doc <- render_emm_contrasts(doc, emm_3a_phase,
    "EMM Phase Contrasts: Model 3a",
    "Phase contrasts at each dose level (Model 3a).")
  doc <- render_emm_contrasts(doc, emm_3c_dose,
    "EMM Dose Contrasts: Model 3c (ANCOVA)",
    "Pairwise dose contrasts within each phase class (Model 3c ANCOVA, Tukey-adjusted).")
  doc <- render_emm_contrasts(doc, emm_3c_phase,
    "EMM Phase Contrasts: Model 3c (ANCOVA)",
    "Phase contrasts at each dose level (Model 3c ANCOVA).")

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

  # Note: Model 5a and 5a-gf2 effect-size tables are placed in the
  # "PRIMARY MODELS" section at the top of this document, not here.

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

# ========================================================================
# Per-subject per-channel phase scatter plots
# ========================================================================
# For each subject: one image with channels as rows, doses as columns.
# Each dot = one probe trial. x = per-burst delivered phase (burstCircMean),
# y = baseline-normalized EP magnitude (magnitude - channel baseline).
# Organized into multi-phase subjects (c91479, 702d24, 0b5a2e) and
# single-phase subjects (d5cd55, 7dbdec, 9ab7ab). ecb43e (triple) separate.
# Loess smoother shows phase-response trend per channel per dose.

if (savePlot) {
  cat("\n========== Per-subject per-channel phase scatter plots ==========\n")

  multi_sids <- c("c91479", "702d24", "0b5a2e")
  single_sids <- c("d5cd55", "7dbdec", "9ab7ab")
  triple_sids <- c("ecb43e")

  if (exists("precision_combined")) {
    # precision_combined already loaded before Model 6
    precision_combined$channelEncoded <- as.factor(precision_combined$channelEncoded)

    # merge burstCircMean into trial-level data
    prec_cols <- precision_combined[, c("probeSample", "channelEncoded", "burstCircMean")]
    data_phase <- merge(data, prec_cols,
      by.x = c("probeSample", "channel"), by.y = c("probeSample", "channelEncoded"), all.x = TRUE)

    # use basePerChan already computed (line ~407): one median per (sid, channel)
    # rename to avoid collision with existing baseMedian column from earlier loop
    base_for_scatter <- basePerChan
    names(base_for_scatter)[names(base_for_scatter) == "baselineMag"] <- "base_scatter"
    data_phase <- merge(data_phase, base_for_scatter, by = c("sid", "channel"), all.x = TRUE)
    data_phase$diff_from_base <- data_phase$magnitude - data_phase$base_scatter

    # exclude baselines for plotting
    data_phase_NB <- data_phase[data_phase$numStims != "Base", ]
    data_phase_NB$numStims <- factor(data_phase_NB$numStims, levels = c("[1,2]", "[3,4]", "[5,inf)"))

    # function to make per-subject plot
    make_subj_plot <- function(sid_val, subj_data, title_suffix = "") {
      subj_data$chan_label <- paste0("Ch ", subj_data$channel)
      chan_order <- sort(unique(subj_data$chan_label))
      subj_data$chan_label <- factor(subj_data$chan_label, levels = chan_order)
      nChan <- length(chan_order)

      p <- ggplot(subj_data, aes(x = burstCircMean, y = diff_from_base)) +
        theme_light(base_size = 10) +
        facet_grid(chan_label ~ numStims, scales = "free_y") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_point(alpha = 0.3, size = 0.8) +
        geom_smooth(method = "loess", se = TRUE, color = "red", linewidth = 0.6, span = 0.75) +
        labs(x = "Per-Burst Delivered Phase (degrees)",
             y = expression(paste(Delta, " from Baseline (", mu, "V)")),
             title = paste0(sid_val, ": Baseline-Normalized EP vs Burst Phase", title_suffix)) +
        scale_x_continuous(breaks = seq(0, 315, by = 90))

      fig_height <- max(3, nChan * 2)
      ggsave(here("output_plots", sprintf("betaStim_%s_phase_scatter.png", sid_val)),
             plot = p, units = "in", width = 10, height = fig_height, dpi = 600)
      cat(sprintf("  %s: %d channels, %d probes with burst phase\n",
          sid_val, nChan, sum(!is.na(subj_data$burstCircMean))))
    }

    # --- Multi-phase subjects ---
    cat("\nMulti-phase subjects (2 target phases per channel):\n")
    for (sid_val in multi_sids) {
      subj <- data_phase_NB[data_phase_NB$sid == sid_val & !is.na(data_phase_NB$burstCircMean), ]
      if (nrow(subj) > 0) make_subj_plot(sid_val, subj, " [multi-phase]")
    }

    # --- Single-phase subjects ---
    cat("\nSingle-phase subjects (1 target phase per channel):\n")
    for (sid_val in single_sids) {
      subj <- data_phase_NB[data_phase_NB$sid == sid_val & !is.na(data_phase_NB$burstCircMean), ]
      if (nrow(subj) > 0) make_subj_plot(sid_val, subj, " [single-phase]")
    }

    # --- Triple-condition subject (ecb43e) ---
    cat("\nTriple-condition subject:\n")
    for (sid_val in triple_sids) {
      subj <- data_phase_NB[data_phase_NB$sid == sid_val & !is.na(data_phase_NB$burstCircMean), ]
      if (nrow(subj) > 0) make_subj_plot(sid_val, subj, " [targeted + random]")
    }

  } else {
    cat("No burst_phase_precision CSVs found — skipping per-subject scatter.\n")
  }
}
