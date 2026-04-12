# ------------------------------------------------------------------------
# Compare phase-targeted vs random stimulation for subject 6 (ecb43e)
#
# ecb43e had three experimental conditions:
#   setToDeliverPhase = 270 (targeted depolarizing)
#   setToDeliverPhase = 90  (targeted hyperpolarizing)
#   setToDeliverPhase = 12345 (random phase — not phase-locked)
#
# This script tests whether phase-targeted stimulation produces different
# evoked potentials than random-phase stimulation on the same channel.
# ------------------------------------------------------------------------

setwd('/Users/davidcaldwell/code/betaOscillationTriggerStimPaper')

library('Hmisc')
library('ggplot2')
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

# Print HTML tab_model() summaries? Off by default (opens RStudio Viewer).
showTabModel = FALSE

# channel 55 on subject 6 (beta channel, subjectNum=6, encoded as 6*100+55=655)
chanInt = 55
chanInt1 = paste0(6, chanInt)

# ------------------------------------------------------------------------
data <- read.table(here("data","output_table","betaStim_outputTable_50_new_100_thresh.csv"),header=TRUE,sep = ",",stringsAsFactors=F,
                   colClasses=c("magnitude"="numeric","betaLabels"="factor","sid"="factor","numStims"="factor","stimLevel"="numeric","channel"="factor","subjectNum"="factor","phaseClass"="factor","setToDeliverPhase"="factor"))
data <- subset(data, magnitude<1500)
data <- subset(data, magnitude>25)

data <- subset(data,!is.nan(data$magnitude))
data <- subset(data,data$numStims!='Null')
# rename for ease
data$numStims <- revalue(data$numStims, c("Test 1"="[1,2]","Test 2"="[3,4]","Test 3"="[5,inf)"))

data$percentDiff = 0
data$absDiff = 0
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
         }
    }
  }
}

sapply(data,class)

dataNoBaseline = data[data$numStims != "Base",]
dataSubjOnly <- subset(data, data$sid=='ecb43e')
dataSubjChanOnly <- subset(dataSubjOnly, dataSubjOnly$channel == chanInt1 & dataSubjOnly$numStims != 'Base')

summaryData = ddply(dataSubjOnly, .(sid,setToDeliverPhase,numStims,channel,betaLabels), summarize, percentDiff = median(percentDiff))
summaryDataChan = subset(summaryData, summaryData$channel == chanInt1)

# ------------------------------------------------------------------------
# Box plot: phase-targeted vs random by dose level
# ------------------------------------------------------------------------

p <- ggplot(dataSubjChanOnly, aes(x=numStims, y=absDiff, fill=setToDeliverPhase)) +
  theme_light(base_size = 18) +
  geom_boxplot(notch=TRUE, position=position_dodge(1)) +
  labs(x = 'Number of conditioning stimuli',
       title = 'Changes in EPs for hyperpolarizing,\ndepolarizing, and random phase stimulation',
       y = expression(paste("Absolute difference from baseline (", mu, "V)"))) +
  scale_fill_hue(name="Experimental\nCondition",
                 breaks=c("12345","270","90"),
                 labels=c("Random","Depolarizing","Hyperpolarizing"))
p

figHeight = 6
figWidth = 8
if(savePlot){
ggsave(here("output_plots","betaStim_control_subj_6.png"), units="in", width=figWidth, height=figHeight, dpi=600)
ggsave(here("output_plots","betaStim_control_subj_6.eps"), units="in", width=figWidth, height=figHeight, dpi=600, device=cairo_ps)
}

# ------------------------------------------------------------------------
# Linear model: condition x dose interaction
# ------------------------------------------------------------------------

fit.lm = lm(magnitude ~ numStims + setToDeliverPhase + numStims:setToDeliverPhase,
            data = dataSubjChanOnly)

summary(fit.lm)
anova(fit.lm)

emm_cond <- emmeans(fit.lm, pairwise ~ setToDeliverPhase | numStims)
emm_cond

emm_dose <- emmeans(fit.lm, pairwise ~ numStims | setToDeliverPhase)
emm_dose

if (showTabModel) tab_model(fit.lm)

# ------------------------------------------------------------------------
# Effect sizes (Cohen's d) via emmeans::eff_size
# ------------------------------------------------------------------------

residual_sd <- sigma(fit.lm)
residual_edf <- df.residual(fit.lm)

# condition contrasts at each dose level
emm_by_dose <- emmeans(fit.lm, ~ setToDeliverPhase | numStims)
es_cond <- eff_size(emm_by_dose, sigma = residual_sd, edf = residual_edf)
es_cond

# dose-response within each condition
emm_by_cond <- emmeans(fit.lm, ~ numStims | setToDeliverPhase)
es_dose <- eff_size(emm_by_cond, sigma = residual_sd, edf = residual_edf)
es_dose

# ------------------------------------------------------------------------
# Permutation test: phase-targeted vs random
# ------------------------------------------------------------------------

set.seed(42)
nPerm <- 10000

# compare each targeted condition (270, 90) against random (12345)
# at each dose level

dose_levels <- c("[1,2]", "[3,4]", "[5,inf)")
target_phases <- c("270", "90")

cat("\n=== Permutation test: targeted vs random (channel", chanInt1, ") ===\n\n")

for (target in target_phases) {
  target_label <- ifelse(target == "270", "Depolarizing", "Hyperpolarizing")
  cat(sprintf("--- %s (%s) vs Random ---\n", target_label, target))

  for (dose in dose_levels) {
    dSub <- dataSubjChanOnly[dataSubjChanOnly$numStims == dose &
              (dataSubjChanOnly$setToDeliverPhase == target |
               dataSubjChanOnly$setToDeliverPhase == "12345"),]

    dSub$cond <- ifelse(dSub$setToDeliverPhase == target, "targeted", "random")
    n_targ <- sum(dSub$cond == "targeted")
    n_rand <- sum(dSub$cond == "random")

    obs_stat <- median(dSub$magnitude[dSub$cond == "targeted"]) -
                median(dSub$magnitude[dSub$cond == "random"])

    perm_stats <- numeric(nPerm)
    for (p in 1:nPerm) {
      shuf <- sample(dSub$cond)
      perm_stats[p] <- median(dSub$magnitude[shuf == "targeted"]) -
                       median(dSub$magnitude[shuf == "random"])
    }
    p_val <- mean(abs(perm_stats) >= abs(obs_stat))

    cat(sprintf("  %-8s  diff = %6.1f uV   perm p = %.4f   (n_targ=%d, n_rand=%d)\n",
      dose, obs_stat, p_val, n_targ, n_rand))
  }
  cat("\n")
}
