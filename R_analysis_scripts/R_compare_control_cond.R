# ------------------------------------------------------------------------
setwd('/Users/davidcaldwell/code/betaOscillationTriggerStimPaper')

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
library('dplyr')
library('here')
library('lmerTest')
library('sjPlot')
library('emmeans')
library('wesanderson')
library(car)

# Global contrast coding. contr.sum gives orthogonal (sum-to-zero) contrasts
# for unordered factors so Type III ANOVAs on factors involved in interactions
# yield marginal main effects. contr.poly (R default) retained for ordered
# factors. Must be set BEFORE factors are created/used in model fits.
options(contrasts = c("contr.sum", "contr.poly"))

rootDir = here()

savePlot = 1
figWidth = 8
figHeight = 6

# Print HTML tab_model() summaries? Off by default (opens RStudio Viewer).
showTabModel = FALSE

# ------------------------------------------------------------------------
# Good-fit / burst-quality filters for the 0b5a2e CL-vs-PB paired analysis
# (dataCL_gf / dataPB_gf below and all downstream scatter / permutation /
# per-channel plots). All filters are applied conjunctively (AND) on each
# probe's preceding CL burst. PB probes inherit their matched CL burst's
# metrics (PB was asynchronous, so its own phase fits are uninterpretable).
#
# Defaults match Model 5a-gf in betaStim_R_script.R:
#   nGoodBeta >= 1, no other burst-quality constraints.
# Set any filter to 0 / Inf to disable.
# ------------------------------------------------------------------------
minGoodBetaPerBurst_clpb   <- 1     # min R^2>0.7 beta (12-20 Hz) stims per burst
minBurstVecLength_clpb     <- 0.5   # min circular vector length of good fits in burst (0-1)
maxBurstCircStd_clpb       <- Inf   # max circular std of good fits in burst (degrees)
minGoodFitFrac_clpb        <- 0     # min nGoodBeta / nCondStims ratio (0-1)
minNCondStimsPerBurst_clpb <- 0     # min total conditioning stims in burst

# Channel-level phase-quality filter (analog of minPhaseVecLength_gf2 in
# betaStim_R_script.R). Drops entire (channel × setToDeliverPhase) cells
# from the clpb analysis if the CL session's pooled phase circular vector
# length was below this threshold — i.e., the channel's local beta was
# not tightly phase-locked during CL targeting, so the "phase condition"
# label is meaningless. Applied BEFORE the paired merge, based on CL's
# phaseVecLength (PB's own r is uninterpretable — asynchronous delivery).
# Set 0 to disable. Default 0.2 matches Model 5a-gf2 in the main script.
minPhaseVecLength_clpb <- 0.2

chanInt = 14
chanInt1 = paste0(7,chanInt)
chanInt2 = paste0(8,chanInt)

# ------------------------------------------------------------------------
data <- read.table(here("data","output_table","betaStim_outputTable_50_new_100_thresh.csv"),header=TRUE,sep = ",",stringsAsFactors=F,
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
#summaryData = ddply(data[data$numStims != "Base",] , .(sid,phaseClass,numStims,channel), function(x) mean(x[,"percentDiff"]))

dataNoBaseline = data[data$numStims != "Base",]
dataSubjOnly <- subset(data,data$sid=='0b5a2e' | data$sid=='0b5a2ePlayBack')
dataSubjChanOnly <- subset(dataSubjOnly,dataSubjOnly$channel == chanInt1 | dataSubjOnly$channel == chanInt2)
#summaryData = ddply(dataSubjOnly[dataSubjOnly$numStims != "Base",] , .(sid,phaseClass,numStims,channel,betaLabels), summarize, percentDiff = mean(percentDiff))
#summaryData = ddply(dataSubjOnly, .(sid,phaseClass,numStims,channel,betaLabels), summarize, percentDiff = mean(percentDiff))
summaryData = ddply(dataSubjOnly, .(sid,phaseClass,numStims,channel,betaLabels), summarize, medianMag = median(magnitude), sdMag = sd(magnitude))

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
  ggsave(here("output_plots","betaStim_control_subj_7.png"), units="in", width=figWidth, height=figHeight,dpi=600)
  ggsave(here("output_plots","betaStim_control_subj_7.eps"), units="in", width=figWidth, height=figHeight, dpi=600, device=cairo_ps)
  
  
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
if (showTabModel) tab_model(fit.lm)

summary(glht(fit.lm,linfct=mcp(sid="Tukey")))
summary(glht(fit.lm,linfct=mcp(numStims="Tukey")))

# ------------------------------------------------------------------------
# Log-transformed linear model (trial-level)
# ------------------------------------------------------------------------
# Log reduces skew (2.5 → 0.3) and downweights extreme CL trials.
# CL condition has heavy right tail that grows with dose (sd 138→243 uV);
# PB is symmetric (sd ~65 uV constant). Log makes residuals more normal.

dataSubjChanOnly$logMag <- log(dataSubjChanOnly$magnitude)

fit.lm.log = lm(logMag ~ numStims + sid + numStims:sid, data = dataSubjChanOnly)

cat("\n=== Log-transformed linear model (channel 14) ===\n")
summary(fit.lm.log)
cat("\nANOVA (Type III):\n")
print(car::Anova(fit.lm.log, type = 3))

emm_log_cond <- emmeans(fit.lm.log, ~ sid | numStims)
cat("\nCL vs PB at each dose (log scale, as % change):\n")
contr_log <- as.data.frame(confint(pairs(emm_log_cond)))
contr_log$pct_change <- round((exp(contr_log$estimate) - 1) * 100, 1)
contr_log$pct_lo <- round((exp(contr_log$lower.CL) - 1) * 100, 1)
contr_log$pct_hi <- round((exp(contr_log$upper.CL) - 1) * 100, 1)
print(contr_log[, c("numStims", "estimate", "pct_change", "pct_lo", "pct_hi")])

emm_log_dose <- emmeans(fit.lm.log, ~ numStims | sid)
cat("\nDose contrasts within CL (log, % change):\n")
dose_log <- as.data.frame(confint(pairs(emm_log_dose)))
dose_log$pct <- round((exp(dose_log$estimate) - 1) * 100, 1)
print(dose_log[dose_log$sid == "0b5a2e", c("contrast", "pct")])
cat("\nDose contrasts within PB (log, % change):\n")
print(dose_log[dose_log$sid == "0b5a2ePlayBack", c("contrast", "pct")])

cat(sprintf("\nResidual diagnostics: raw skew=%.2f, log skew=%.2f\n",
  mean((resid(fit.lm)/sd(resid(fit.lm)))^3),
  mean((resid(fit.lm.log)/sd(resid(fit.lm.log)))^3)))

# ------------------------------------------------------------------------
# Effect sizes (Cohen's d) via emmeans::eff_size
# ------------------------------------------------------------------------

residual_sd <- sigma(fit.lm)
residual_edf <- df.residual(fit.lm)

# Cohen's d: closed-loop vs playback at each dose level
emm_by_dose <- emmeans(fit.lm, ~ sid | numStims)
es_cl_vs_pb <- eff_size(emm_by_dose, sigma = residual_sd, edf = residual_edf)
es_cl_vs_pb

# Cohen's d: dose-response within each condition
emm_by_sid <- emmeans(fit.lm, ~ numStims | sid)
es_dose_response <- eff_size(emm_by_sid, sigma = residual_sd, edf = residual_edf)
es_dose_response

# Cohen's d bar plot: closed-loop vs playback by dose level
es_df <- as.data.frame(es_cl_vs_pb)
es_df <- es_df[es_df$numStims != "Base",]
es_df$numStims <- factor(es_df$numStims, levels = c("Null","[1,2]","[3,4]","[5,inf)"))

p_cohend <- ggplot(es_df, aes(x = numStims, y = effect.size)) +
  theme_light(base_size = 14) +
  geom_col(fill = wes_palette(n=1, name="GrandBudapest1"), width = 0.6) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_hline(yintercept = c(0.2, 0.5, 0.8), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  annotate("text", x = 4.4, y = 0.2, label = "small", size = 3, color = "grey40") +
  annotate("text", x = 4.4, y = 0.5, label = "medium", size = 3, color = "grey40") +
  annotate("text", x = 4.4, y = 0.8, label = "large", size = 3, color = "grey40") +
  labs(x = "Number of Conditioning Stimuli",
       y = "Cohen's d (95% CI)",
       title = "Effect Size: Closed-loop vs Playback Control") +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 30, 5, 5))
p_cohend

figHeight = 4
figWidth = 6
if(savePlot){
ggsave(here("output_plots","betaStim_control_cohens_d.png"), plot = p_cohend,
       units = "in", width = figWidth, height = figHeight, dpi = 600)
ggsave(here("output_plots","betaStim_control_cohens_d.eps"), plot = p_cohend,
       units = "in", width = figWidth, height = figHeight, dpi = 600, device = cairo_ps)
}

# ------------------------------------------------------------------------
# Permutation test: closed-loop vs playback on channel 14
# ------------------------------------------------------------------------

set.seed(42)
nPerm <- 10000

dataSubjChanOnly$condition <- ifelse(dataSubjChanOnly$sid == "0b5a2e", "CL", "PB")
dataPerm <- subset(dataSubjChanOnly, numStims != "Null")

dose_levels <- c("Base", "[1,2]", "[3,4]", "[5,inf)")

# per-dose permutation test
perm_results <- data.frame(numStims = character(), obs_diff = numeric(),
                           perm_p = numeric(), stringsAsFactors = FALSE)

for (dose in dose_levels) {
  dSub <- dataPerm[dataPerm$numStims == dose,]
  obs_stat <- median(dSub$magnitude[dSub$condition == "CL"]) -
              median(dSub$magnitude[dSub$condition == "PB"])

  perm_stats <- numeric(nPerm)
  for (p in 1:nPerm) {
    shuf <- sample(dSub$condition)
    perm_stats[p] <- median(dSub$magnitude[shuf == "CL"]) -
                     median(dSub$magnitude[shuf == "PB"])
  }
  p_val <- mean(abs(perm_stats) >= abs(obs_stat))
  perm_results <- rbind(perm_results,
    data.frame(numStims = dose, obs_diff = obs_stat, perm_p = p_val))
}

perm_results

# dose-response interaction: does the dose effect differ between CL and PB?
# Permute dose labels (Base vs [5,inf)) WITHIN each condition separately.
# This respects the session structure — CL and PB trials are from different
# sessions and should not be shuffled across conditions. Dose labels ARE
# exchangeable within a session.
dCL <- dataPerm[dataPerm$condition == "CL" &
                (dataPerm$numStims == "Base" | dataPerm$numStims == "[5,inf)"),]
dPB <- dataPerm[dataPerm$condition == "PB" &
                (dataPerm$numStims == "Base" | dataPerm$numStims == "[5,inf)"),]

# observed dose effect within each condition (median scale)
obs_cl_dose <- median(dCL$magnitude[dCL$numStims == "[5,inf)"]) -
               median(dCL$magnitude[dCL$numStims == "Base"])
obs_pb_dose <- median(dPB$magnitude[dPB$numStims == "[5,inf)"]) -
               median(dPB$magnitude[dPB$numStims == "Base"])
obs_interaction <- obs_cl_dose - obs_pb_dose

perm_interactions <- numeric(nPerm)
for (p in 1:nPerm) {
  # shuffle dose labels within CL
  shuf_cl <- sample(dCL$numStims)
  perm_cl_dose <- median(dCL$magnitude[shuf_cl == "[5,inf)"]) -
                  median(dCL$magnitude[shuf_cl == "Base"])
  # shuffle dose labels within PB
  shuf_pb <- sample(dPB$numStims)
  perm_pb_dose <- median(dPB$magnitude[shuf_pb == "[5,inf)"]) -
                  median(dPB$magnitude[shuf_pb == "Base"])
  perm_interactions[p] <- perm_cl_dose - perm_pb_dose
}

interaction_p <- mean(abs(perm_interactions) >= abs(obs_interaction))

cat(sprintf("\nDose-response interaction — median, dose-permuted (channel %d):\n", chanInt))
cat(sprintf("  CL dose effect ([5,inf) - Base):  %+6.1f uV\n", obs_cl_dose))
cat(sprintf("  PB dose effect ([5,inf) - Base):  %+6.1f uV\n", obs_pb_dose))
cat(sprintf("  Interaction (CL - PB):            %+6.1f uV   perm p = %.4f (two-sided)\n", obs_interaction, interaction_p))

# same on log scale (robust to heavy tails)
obs_cl_dose_log <- median(log(dCL$magnitude[dCL$numStims == "[5,inf)"])) -
                   median(log(dCL$magnitude[dCL$numStims == "Base"]))
obs_pb_dose_log <- median(log(dPB$magnitude[dPB$numStims == "[5,inf)"])) -
                   median(log(dPB$magnitude[dPB$numStims == "Base"]))
obs_interaction_log <- obs_cl_dose_log - obs_pb_dose_log

perm_interactions_log <- numeric(nPerm)
for (p in 1:nPerm) {
  shuf_cl <- sample(dCL$numStims)
  perm_cl_log <- median(log(dCL$magnitude[shuf_cl == "[5,inf)"])) -
                 median(log(dCL$magnitude[shuf_cl == "Base"]))
  shuf_pb <- sample(dPB$numStims)
  perm_pb_log <- median(log(dPB$magnitude[shuf_pb == "[5,inf)"])) -
                 median(log(dPB$magnitude[shuf_pb == "Base"]))
  perm_interactions_log[p] <- perm_cl_log - perm_pb_log
}

interaction_p_log <- mean(abs(perm_interactions_log) >= abs(obs_interaction_log))

cat(sprintf("\nDose-response interaction — log, dose-permuted (channel %d):\n", chanInt))
cat(sprintf("  CL dose effect (log):  %+.4f (%+.1f%%)\n", obs_cl_dose_log, (exp(obs_cl_dose_log)-1)*100))
cat(sprintf("  PB dose effect (log):  %+.4f (%+.1f%%)\n", obs_pb_dose_log, (exp(obs_pb_dose_log)-1)*100))
cat(sprintf("  Interaction (log):     %+.4f   perm p = %.4f (two-sided)\n", obs_interaction_log, interaction_p_log))

# store per-dose null distributions for plotting
perm_null_dists <- list()
obs_diffs_vec <- c()
perm_ps_vec <- c()

for (dose in dose_levels) {
  dSub <- dataPerm[dataPerm$numStims == dose,]
  obs_stat <- median(dSub$magnitude[dSub$condition == "CL"]) -
              median(dSub$magnitude[dSub$condition == "PB"])
  obs_diffs_vec <- c(obs_diffs_vec, obs_stat)

  perm_stats <- numeric(nPerm)
  for (p in 1:nPerm) {
    shuf <- sample(dSub$condition)
    perm_stats[p] <- median(dSub$magnitude[shuf == "CL"]) -
                     median(dSub$magnitude[shuf == "PB"])
  }
  perm_null_dists[[dose]] <- perm_stats
  perm_ps_vec <- c(perm_ps_vec, mean(abs(perm_stats) >= abs(obs_stat)))
}

# ------------------------------------------------------------------------
# Plot: null distribution for the interaction test (two-sided)
# ------------------------------------------------------------------------

null_df <- data.frame(value = perm_interactions)

p_perm_interaction <- ggplot(null_df, aes(x = value)) +
  theme_light(base_size = 14) +
  geom_histogram(aes(fill = abs(value) >= abs(obs_interaction)),
                 bins = 60, color = "grey50", show.legend = FALSE) +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#C93312")) +
  geom_vline(xintercept = obs_interaction, color = "#C93312", linewidth = 1.2) +
  geom_vline(xintercept = -obs_interaction, color = "#C93312", linewidth = 1.2, linetype = "dashed") +
  annotate("text", x = obs_interaction + 3, y = Inf, vjust = 2, hjust = 0,
           label = sprintf("observed = %.1f uV\np = %.3f (two-sided)", obs_interaction, interaction_p),
           color = "#C93312", size = 4) +
  labs(x = expression(paste(Delta, " (CL-PB gap at [5,inf)) - (CL-PB gap at Base)  [",mu,"V]")),
       y = "Count",
       title = sprintf("Permutation Test: Dose-Response Interaction (Channel %d)", chanInt),
       subtitle = sprintf("%s permutations", formatC(nPerm, format="d", big.mark=",")))
p_perm_interaction

figHeight = 4
figWidth = 7
if(savePlot){
ggsave(here("output_plots","betaStim_perm_interaction_ch14.png"), plot = p_perm_interaction,
       units = "in", width = figWidth, height = figHeight, dpi = 600)
ggsave(here("output_plots","betaStim_perm_interaction_ch14.eps"), plot = p_perm_interaction,
       units = "in", width = figWidth, height = figHeight, dpi = 600, device = cairo_ps)
}

# ------------------------------------------------------------------------
# Plot: per-dose CL-PB differences against null distribution
# ------------------------------------------------------------------------

dose_df <- data.frame(
  numStims = factor(dose_levels, levels = dose_levels),
  obs_diff = obs_diffs_vec,
  null_lo = sapply(perm_null_dists, quantile, 0.025),
  null_hi = sapply(perm_null_dists, quantile, 0.975),
  perm_p = perm_ps_vec
)
dose_df$p_label <- ifelse(dose_df$perm_p < 0.0001, "p < 0.0001",
                   ifelse(dose_df$perm_p < 0.001, sprintf("p = %.4f", dose_df$perm_p),
                          sprintf("p = %.2f", dose_df$perm_p)))

p_perm_dose <- ggplot(dose_df, aes(x = numStims)) +
  theme_light(base_size = 14) +
  geom_crossbar(aes(y = 0, ymin = null_lo, ymax = null_hi),
                fill = "grey85", color = "grey60", width = 0.5, middle.linewidth = 0) +
  geom_point(aes(y = obs_diff), color = wes_palette(n=1, name="GrandBudapest1"), size = 4) +
  geom_segment(aes(xend = numStims, y = 0, yend = obs_diff),
               color = wes_palette(n=1, name="GrandBudapest1"), linewidth = 1) +
  geom_text(aes(y = obs_diff, label = p_label), vjust = -1, size = 3.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Number of Conditioning Stimuli",
       y = expression(paste("CL - Playback  [",mu,"V]")),
       title = sprintf("Closed-loop vs Playback by Dose (Channel %d)", chanInt),
       subtitle = "Grey bars = 95% of null distribution") +
  coord_cartesian(ylim = c(-60, 140))
p_perm_dose

figHeight = 4.5
figWidth = 6
if(savePlot){
ggsave(here("output_plots","betaStim_perm_dose_ch14.png"), plot = p_perm_dose,
       units = "in", width = figWidth, height = figHeight, dpi = 600)
ggsave(here("output_plots","betaStim_perm_dose_ch14.eps"), plot = p_perm_dose,
       units = "in", width = figWidth, height = figHeight, dpi = 600, device = cairo_ps)
}

data$phaseDeg <- as.numeric(as.character(data$phaseDeg))
data$phase_rad <- data$phaseDeg * pi / 180
data$sin_phase <- sin(data$phase_rad)
data$cos_phase <- cos(data$phase_rad)
data$phaseDeg_round <- round(data$phaseDeg, 1)

# ========================================================================
# 0b5a2e (CL) vs 0b5a2ePlayBack (PB): MATCHED probe comparison
# ========================================================================
cat("\n========== 0b5a2e: Matched CL vs PB ==========\n")

dataCL_raw <- data[data$sid == "0b5a2e", ]
dataPB_raw <- data[data$sid == "0b5a2ePlayBack", ]
dataCL_raw$channel_raw <- as.factor(as.numeric(as.character(dataCL_raw$channel)) %% 100)
dataPB_raw$channel_raw <- as.factor(as.numeric(as.character(dataPB_raw$channel)) %% 100)

dataCL_raw <- dataCL_raw[order(dataCL_raw$channel_raw, dataCL_raw$probeSample), ]
dataCL_raw$probeIdx <- unlist(tapply(dataCL_raw$probeSample, dataCL_raw$channel_raw,
  function(x) seq_along(x)))
dataPB_raw <- dataPB_raw[order(dataPB_raw$channel_raw, dataPB_raw$probeSample), ]
dataPB_raw$probeIdx <- unlist(tapply(dataPB_raw$probeSample, dataPB_raw$channel_raw,
  function(x) seq_along(x)))

cl_precision <- read.csv(here("data", "output_table", "0b5a2e_burst_phase_precision.csv"))
cl_precision$channelEncoded <- as.factor(cl_precision$channelEncoded)
# Pull in all burst-quality columns we may filter on. goodFitFrac is derived
# (nGoodBeta / nCondStims) rather than stored.
dataCL_raw <- merge(dataCL_raw,
  cl_precision[, c("probeSample", "channelEncoded", "nGoodBeta", "nCondStims",
                   "burstVecLength", "burstCircStd")],
  by.x = c("probeSample", "channel"), by.y = c("probeSample", "channelEncoded"),
  all.x = TRUE)
dataCL_raw$goodFitFrac <- ifelse(is.na(dataCL_raw$nCondStims) | dataCL_raw$nCondStims == 0,
                                 NA, dataCL_raw$nGoodBeta / dataCL_raw$nCondStims)

# Propagate the CL burst metrics to matched PB probes (PB's own fits are
# uninterpretable because delivery was asynchronous).
cl_filter <- dataCL_raw[, c("channel_raw", "probeIdx", "nGoodBeta", "nCondStims",
                            "burstVecLength", "burstCircStd", "goodFitFrac")]
names(cl_filter)[3:7] <- c("nGoodBeta_CL", "nCondStims_CL",
                           "burstVecLength_CL", "burstCircStd_CL", "goodFitFrac_CL")
dataPB_raw <- merge(dataPB_raw, cl_filter, by = c("channel_raw", "probeIdx"), all.x = TRUE)

# Exact sign-flip permutation: enumerate all 2^n ± sign vectors and return
# the observed statistic, the complete null distribution, and the exact
# two-sided p-value. Use when n <= 20 (2^20 ≈ 10^6 configurations is the
# practical limit for an in-memory matrix). For larger n, use Monte Carlo.
#
# Rationale: for small-n paired/aggregated tests (e.g., 8 channels, 13-16
# channel × condition cells), Monte Carlo with 10,000 draws just resamples
# the same 2^n configurations many times over. Enumerating gives an exact,
# reproducible p-value with no Monte Carlo error.
#
# Fast path: when stat_fn is base `mean`, use a single matrix multiply
# (sign_grid %*% x) / n to compute all 2^n means at once. For other
# statistics (e.g., median) fall back to row-wise apply().
exact_signflip <- function(x, stat_fn = mean) {
  n <- length(x)
  if (n < 2) stop("exact_signflip: need at least 2 units")
  if (n > 20) {
    stop(sprintf("exact_signflip: n=%d too large for enumeration (2^n = %d); use Monte Carlo",
                 n, 2^n))
  }
  obs <- stat_fn(x)
  sign_grid <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  if (identical(stat_fn, mean)) {
    perm_stats <- as.numeric(sign_grid %*% x) / n
  } else {
    perm_stats <- apply(sign_grid, 1, function(s) stat_fn(s * x))
  }
  list(
    obs = obs,
    perm_stats = perm_stats,
    p_two_sided = mean(abs(perm_stats) >= abs(obs)),
    n_perms = length(perm_stats)
  )
}

# Conjunctive filter using the configurable thresholds at the top of the
# script. apply_burst_filters keeps each row only if every ACTIVE filter
# passes (NAs fail, since a missing precision row means the probe wasn't
# scored and we can't verify its burst quality).
apply_burst_filters <- function(df, ng_col, bvl_col, bcs_col, gff_col, ncs_col) {
  keep <- !is.na(df[[ng_col]]) & df[[ng_col]] >= minGoodBetaPerBurst_clpb
  if (minBurstVecLength_clpb > 0) {
    keep <- keep & !is.na(df[[bvl_col]]) & df[[bvl_col]] >= minBurstVecLength_clpb
  }
  if (is.finite(maxBurstCircStd_clpb)) {
    keep <- keep & !is.na(df[[bcs_col]]) & df[[bcs_col]] <= maxBurstCircStd_clpb
  }
  if (minGoodFitFrac_clpb > 0) {
    keep <- keep & !is.na(df[[gff_col]]) & df[[gff_col]] >= minGoodFitFrac_clpb
  }
  if (minNCondStimsPerBurst_clpb > 0) {
    keep <- keep & !is.na(df[[ncs_col]]) & df[[ncs_col]] >= minNCondStimsPerBurst_clpb
  }
  df[keep, ]
}

dataCL_gf <- apply_burst_filters(dataCL_raw, "nGoodBeta", "burstVecLength",
                                 "burstCircStd", "goodFitFrac", "nCondStims")
dataPB_gf <- apply_burst_filters(dataPB_raw, "nGoodBeta_CL", "burstVecLength_CL",
                                 "burstCircStd_CL", "goodFitFrac_CL", "nCondStims_CL")

cat(sprintf("clpb burst filters: nGoodBeta>=%d, burstVecLength>=%.2f, burstCircStd<=%s, goodFitFrac>=%.2f, nCondStims>=%d\n",
    minGoodBetaPerBurst_clpb, minBurstVecLength_clpb,
    ifelse(is.finite(maxBurstCircStd_clpb), sprintf("%.1f", maxBurstCircStd_clpb), "Inf"),
    minGoodFitFrac_clpb, minNCondStimsPerBurst_clpb))
cat(sprintf("CL burst-filtered: %d of %d (%.0f%%)\n",
    nrow(dataCL_gf), nrow(dataCL_raw), 100*nrow(dataCL_gf)/nrow(dataCL_raw)))
cat(sprintf("PB burst-filtered: %d of %d (%.0f%%)\n",
    nrow(dataPB_gf), nrow(dataPB_raw), 100*nrow(dataPB_gf)/nrow(dataPB_raw)))

# Channel-level phase-quality filter: drop whole (channel × setToDeliverPhase)
# cells whose CL-session pooled phaseVecLength is below threshold. Applied
# symmetrically to both CL and PB via cell-key matching so pair alignment
# is preserved. PB's own phaseVecLength is NOT used (asynchronous delivery
# makes it uninterpretable).
if (minPhaseVecLength_clpb > 0) {
  cl_keep_cells <- unique(dataCL_raw[!is.na(dataCL_raw$phaseVecLength) &
    dataCL_raw$phaseVecLength >= minPhaseVecLength_clpb,
    c("channel_raw", "setToDeliverPhase")])
  keep_key <- paste(cl_keep_cells$channel_raw, cl_keep_cells$setToDeliverPhase, sep = "__")
  cl_key_vec <- paste(dataCL_gf$channel_raw, dataCL_gf$setToDeliverPhase, sep = "__")
  pb_key_vec <- paste(dataPB_gf$channel_raw, dataPB_gf$setToDeliverPhase, sep = "__")
  n_before_cl <- nrow(dataCL_gf)
  n_before_pb <- nrow(dataPB_gf)
  dataCL_gf <- dataCL_gf[cl_key_vec %in% keep_key, ]
  dataPB_gf <- dataPB_gf[pb_key_vec %in% keep_key, ]
  cat(sprintf("Channel phaseVecLength >= %.2f (CL-based): CL %d->%d, PB %d->%d (kept %d cells)\n",
      minPhaseVecLength_clpb, n_before_cl, nrow(dataCL_gf),
      n_before_pb, nrow(dataPB_gf), length(keep_key)))
}

# Each condition uses its OWN measured phase (not CL phase for both).
# CL phase = where stims were intentionally targeted.
# PB phase = where stims accidentally landed (asynchronous delivery).
# If phase effect is real, CL should show modulation; PB should not,
# because PB phases are random relative to the oscillation.

dataCL_gf <- dataCL_gf[dataCL_gf$magnitude > 25 & dataCL_gf$magnitude < 1500, ]
dataPB_gf <- dataPB_gf[dataPB_gf$magnitude > 25 & dataPB_gf$magnitude < 1500, ]

dataCL_gf$condition <- "CL"
dataPB_gf$condition <- "PB"
common_cols <- intersect(names(dataCL_gf), names(dataPB_gf))
dataBoth <- rbind(dataCL_gf[, common_cols], dataPB_gf[, common_cols])
dataBoth$condition <- factor(dataBoth$condition, levels = c("PB", "CL"))
dataBoth_NB <- dataBoth[dataBoth$numStims != "Null" & dataBoth$numStims != "Base", ]

summaryMatched <- ddply(dataBoth_NB, .(condition, channel_raw, phaseDeg_round, numStims),
  summarize, magnitude = median(magnitude),
  sin_phase = first(sin_phase), cos_phase = first(cos_phase))

# baselines from UNFILTERED data — baseline probes have nGoodBeta=0
# (no conditioning stims), so good-fit filter would remove them all
dataCL_base <- dataCL_raw[dataCL_raw$numStims == "Base" & dataCL_raw$magnitude > 25 & dataCL_raw$magnitude < 1500, ]
dataPB_base <- dataPB_raw[dataPB_raw$numStims == "Base" & dataPB_raw$magnitude > 25 & dataPB_raw$magnitude < 1500, ]
dataCL_base$condition <- "CL"; dataCL_base$channel_raw <- as.factor(as.numeric(as.character(dataCL_base$channel)) %% 100)
dataPB_base$condition <- "PB"; dataPB_base$channel_raw <- as.factor(as.numeric(as.character(dataPB_base$channel)) %% 100)
baseMatched <- rbind(
  ddply(dataCL_base, .(condition, channel_raw), summarize, baselineMag = median(magnitude)),
  ddply(dataPB_base, .(condition, channel_raw), summarize, baselineMag = median(magnitude)))
summaryMatched <- merge(summaryMatched, baseMatched, by = c("condition", "channel_raw"))
summaryMatched$baselineMag_c <- summaryMatched$baselineMag - mean(summaryMatched$baselineMag)
summaryMatched$numStims_ord <- ordered(summaryMatched$numStims,
  levels = c("[1,2]", "[3,4]", "[5,inf)"))

cat(sprintf("\nMatched summary: %d obs, %d channels, CL=%d, PB=%d\n",
    nrow(summaryMatched), length(unique(summaryMatched$channel_raw)),
    sum(summaryMatched$condition == "CL"), sum(summaryMatched$condition == "PB")))

# Each condition uses its own measured phase. condition ref = PB.
# sin_phase/cos_phase = phase effect at PB (accidental phases — should be null).
# conditionCL:sin_phase = CL-specific phase effect BEYOND PB.
# No dose:condition (assumes same average dose-response in both sessions).
cat("\n--- Combined CL+PB model (each condition's own phase) ---\n")

fit.clpb = lmerTest::lmer(
  magnitude ~ numStims_ord * (sin_phase + cos_phase) +
              condition * (sin_phase + cos_phase) +
              baselineMag_c +
  (1 | channel_raw),
  data = summaryMatched,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

cat("Singular:", isSingular(fit.clpb), "\n")
print(summary(fit.clpb))
cat("\nType III ANOVA:\n")
print(anova(fit.clpb))

# coefficient names: check actual names from the model
cat("\nCoefficient names:", paste(names(fixef(fit.clpb)), collapse=", "), "\n")
cnames <- names(fixef(fit.clpb))
# find the condition:phase interaction terms
cond_sin <- cnames[grepl("condition.*sin_phase|sin_phase.*condition", cnames)]
cond_cos <- cnames[grepl("condition.*cos_phase|cos_phase.*condition", cnames)]
cat("\nWald: CL-specific phase effect beyond PB:\n")
print(car::linearHypothesis(fit.clpb,
  c(paste0(cond_sin, " = 0"), paste0(cond_cos, " = 0"))))

cat("\nWald: phase at PB reference (should be ns if asynchronous):\n")
print(car::linearHypothesis(fit.clpb, c("sin_phase = 0", "cos_phase = 0")))

phase_vals_ctrl <- seq(0, 315, by = 45)
emm_curves <- lapply(c("CL", "PB"), function(cond) {
  do.call(rbind, lapply(phase_vals_ctrl, function(ph) {
    em <- emmeans(fit.clpb, ~ numStims_ord,
      at = list(sin_phase = sin(ph*pi/180), cos_phase = cos(ph*pi/180),
                condition = cond, baselineMag_c = 0))
    df <- as.data.frame(em)
    df$phase_deg <- ph; df$condition <- cond; df
  }))
})
emm_ctrl_df <- do.call(rbind, emm_curves)

p_clpb <- ggplot(emm_ctrl_df, aes(x = phase_deg, y = emmean, color = numStims_ord)) +
  theme_light(base_size = 14) + facet_wrap(~ condition) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = numStims_ord), alpha = 0.15, color = NA) +
  labs(x = "Delivered Phase (degrees)",
       y = expression(paste("Predicted Magnitude (", mu, "V)")),
       color = "Dose", fill = "Dose",
       title = "CL vs PB: Phase-Response (matched, ANCOVA, CL phase + good-fit)") +
  scale_x_continuous(breaks = seq(0, 315, by = 90))
p_clpb
if(savePlot){
  ggsave(here("output_plots","betaStim_clpb_matched_phase_curve.png"), plot = p_clpb,
         units = "in", width = 10, height = 4.5, dpi = 600)
}

# ========================================================================
# 0b5a2e CL vs PB: SENSITIVITY — per-burst phase, each condition's own
# ========================================================================
# Phase = circular mean of burstCircMean from good-fit bursts within each
# summary cell. CL uses CL's burstCircMean, PB uses PB's burstCircMean.
# PB's per-burst phase = where stims actually landed (accidental phases).
cat("\n--- CL vs PB sensitivity: per-burst phase (each condition's own) ---\n")

# CL per-burst phase (already have precision data)
cl_burst_cols <- cl_precision[, c("probeSample", "channelEncoded", "burstCircMean")]
cl_burst_cols$channelEncoded <- as.factor(cl_burst_cols$channelEncoded)
dataCL_gf_burst <- merge(dataCL_gf, cl_burst_cols,
  by.x = c("probeSample", "channel"), by.y = c("probeSample", "channelEncoded"), all.x = TRUE)
dataCL_gf_burst$sin_phase_burst <- sin(dataCL_gf_burst$burstCircMean * pi / 180)
dataCL_gf_burst$cos_phase_burst <- cos(dataCL_gf_burst$burstCircMean * pi / 180)

# PB per-burst phase from PB's own precision data
pb_precision <- read.csv(here("data", "output_table", "0b5a2ePlayback_burst_phase_precision.csv"))
pb_precision$channelEncoded <- as.factor(pb_precision$channelEncoded)
# PB probes need to match by probeIdx (not probeSample — different sessions)
# assign probeIdx to PB precision data
pb_precision <- pb_precision[order(pb_precision$channelEncoded, pb_precision$probeSample), ]
pb_precision$probeIdx <- unlist(tapply(pb_precision$probeSample, pb_precision$channelEncoded,
  function(x) seq_along(x)))
pb_precision$channel_raw <- as.factor(as.numeric(as.character(pb_precision$channelEncoded)) %% 100)

dataPB_gf_burst <- merge(dataPB_gf,
  pb_precision[, c("channel_raw", "probeIdx", "burstCircMean")],
  by = c("channel_raw", "probeIdx"), all.x = TRUE)
dataPB_gf_burst$sin_phase_burst <- sin(dataPB_gf_burst$burstCircMean * pi / 180)
dataPB_gf_burst$cos_phase_burst <- cos(dataPB_gf_burst$burstCircMean * pi / 180)

# summary: group by (channel_raw, phaseDeg_round, numStims) per condition
# CL uses CL's phaseDeg_round; PB uses PB's own phaseDeg_round
summCL_burst <- plyr::ddply(
  dataCL_gf_burst[dataCL_gf_burst$numStims != "Null" & dataCL_gf_burst$numStims != "Base", ],
  .(channel_raw, phaseDeg_round, numStims), summarize,
  magnitude = median(magnitude),
  sin_phase = mean(sin_phase_burst, na.rm = TRUE),
  cos_phase = mean(cos_phase_burst, na.rm = TRUE))
summCL_burst$condition <- "CL"

summPB_burst <- plyr::ddply(
  dataPB_gf_burst[dataPB_gf_burst$numStims != "Null" & dataPB_gf_burst$numStims != "Base", ],
  .(channel_raw, phaseDeg_round, numStims), summarize,
  magnitude = median(magnitude),
  sin_phase = mean(sin_phase_burst, na.rm = TRUE),
  cos_phase = mean(cos_phase_burst, na.rm = TRUE))
summPB_burst$condition <- "PB"

summBurst <- rbind(summCL_burst, summPB_burst)
summBurst$condition <- factor(summBurst$condition, levels = c("PB", "CL"))
summBurst <- merge(summBurst, baseMatched, by = c("condition", "channel_raw"))
summBurst$baselineMag_c <- summBurst$baselineMag - mean(summBurst$baselineMag)
summBurst$numStims_ord <- ordered(summBurst$numStims, levels = c("[1,2]", "[3,4]", "[5,inf)"))

cat(sprintf("Per-burst phase summary: %d obs, CL=%d, PB=%d\n",
    nrow(summBurst), sum(summBurst$condition == "CL"), sum(summBurst$condition == "PB")))

fit.clpb.burst = lmerTest::lmer(
  magnitude ~ numStims_ord * (sin_phase + cos_phase) +
              condition * (sin_phase + cos_phase) +
              baselineMag_c +
  (1 | channel_raw),
  data = summBurst,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

cat("Singular:", isSingular(fit.clpb.burst), "\n")
print(summary(fit.clpb.burst))
cat("\nType III ANOVA:\n")
print(anova(fit.clpb.burst))

cnames_burst <- names(fixef(fit.clpb.burst))
cb_sin <- cnames_burst[grepl("condition.*sin_phase|sin_phase.*condition", cnames_burst)]
cb_cos <- cnames_burst[grepl("condition.*cos_phase|cos_phase.*condition", cnames_burst)]
cat("\nWald: CL-specific phase (per-burst) beyond PB:\n")
print(car::linearHypothesis(fit.clpb.burst,
  c(paste0(cb_sin, " = 0"), paste0(cb_cos, " = 0"))))

emm_burst_curves <- lapply(c("CL", "PB"), function(cond) {
  do.call(rbind, lapply(phase_vals_ctrl, function(ph) {
    em <- emmeans(fit.clpb.burst, ~ numStims_ord,
      at = list(sin_phase = sin(ph*pi/180), cos_phase = cos(ph*pi/180),
                condition = cond, baselineMag_c = 0))
    df <- as.data.frame(em)
    df$phase_deg <- ph; df$condition <- cond; df
  }))
})
emm_burst_df <- do.call(rbind, emm_burst_curves)

p_clpb_burst <- ggplot(emm_burst_df, aes(x = phase_deg, y = emmean, color = numStims_ord)) +
  theme_light(base_size = 14) + facet_wrap(~ condition) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = numStims_ord), alpha = 0.15, color = NA) +
  labs(x = "Delivered Phase (degrees)",
       y = expression(paste("Predicted Magnitude (", mu, "V)")),
       color = "Dose", fill = "Dose",
       title = "CL vs PB: Per-burst Phase (each condition's own)") +
  scale_x_continuous(breaks = seq(0, 315, by = 90))
p_clpb_burst
if(savePlot){
  ggsave(here("output_plots","betaStim_clpb_burst_phase_curve.png"), plot = p_clpb_burst,
         units = "in", width = 10, height = 4.5, dpi = 600)
}

# ========================================================================
# Paired probe-level: CL - PB differences by dose and phase
# ========================================================================
# For each matched probe pair (same channel_raw, same probeIdx):
#   1. Normalize each by its own session's channel baseline
#   2. diff = normalized_CL - normalized_PB
#   3. Sign-flip permutation test per dose: is mean diff != 0?
#   4. Plot diff vs CL per-burst phase (burstCircMean) at each dose
cat("\n========== Paired probe-level: CL - PB differences ==========\n")

baseCL_chan <- plyr::ddply(dataCL_raw[dataCL_raw$numStims == "Base" &
  dataCL_raw$magnitude > 25 & dataCL_raw$magnitude < 1500, ],
  .(channel_raw), summarize, baseCL = median(magnitude))
basePB_chan <- plyr::ddply(dataPB_raw[dataPB_raw$numStims == "Base" &
  dataPB_raw$magnitude > 25 & dataPB_raw$magnitude < 1500, ],
  .(channel_raw), summarize, basePB = median(magnitude))

# CL trials with per-burst phase from precision data
dataCL_paired <- merge(dataCL_gf,
  cl_precision[, c("probeSample", "channelEncoded", "burstCircMean")],
  by.x = c("probeSample", "channel"), by.y = c("probeSample", "channelEncoded"), all.x = TRUE)
# Carry phaseVecLength (channel-level r) and setToDeliverPhase (intended
# target) forward so the forest-plot labels can show both alongside the
# measured phase.
dataCL_paired <- dataCL_paired[, c("channel_raw", "probeIdx", "magnitude", "numStims",
  "burstCircMean", "phaseDeg_round", "phaseVecLength", "setToDeliverPhase")]
names(dataCL_paired)[3] <- "mag_CL"
names(dataCL_paired)[5] <- "cl_burst_phase"

# PB trials (magnitude only — phase comes from CL for x-axis)
dataPB_paired <- dataPB_gf[dataPB_gf$magnitude > 25 & dataPB_gf$magnitude < 1500,
  c("channel_raw", "probeIdx", "magnitude")]
names(dataPB_paired)[3] <- "mag_PB"

# merge matched pairs
paired <- merge(dataCL_paired, dataPB_paired, by = c("channel_raw", "probeIdx"))
paired <- merge(paired, baseCL_chan, by = "channel_raw")
paired <- merge(paired, basePB_chan, by = "channel_raw")

# baseline-normalize each session, then difference
paired$norm_CL <- paired$mag_CL - paired$baseCL
paired$norm_PB <- paired$mag_PB - paired$basePB
paired$diff <- paired$norm_CL - paired$norm_PB

paired <- paired[paired$numStims != "Base" & paired$numStims != "Null", ]
paired <- paired[!is.na(paired$cl_burst_phase), ]  # need per-burst phase for x-axis

cat(sprintf("Matched pairs with per-burst phase: %d\n", nrow(paired)))

# --- Channel x condition aggregated sign-flip permutation ---
# Aggregate to (channel x phaseDeg_round) means to avoid pseudoreplication
# AND keep the two target conditions (90 vs 270) separate. 0b5a2e targets
# both phases, so each channel has two conditions with different delivered
# phases. Collapsing across conditions averages ~90 and ~270 to ~180.
# With 8 channels x 2 conditions = up to 16 units per dose.
# Holm correction across 3 dose-level tests.
set.seed(42)
# nPerm already defined above (line 213)
dose_levels_paired <- c("[1,2]", "[3,4]", "[5,inf)")

chan_cond_means <- plyr::ddply(paired, .(channel_raw, phaseDeg_round, numStims), summarize,
  mean_diff = mean(diff), median_diff = median(diff),
  cl_cond_phase = mean(cl_burst_phase, na.rm = TRUE),
  n_probes = length(diff))

# Exact sign-flip enumeration. n = 13-16 channel-condition cells per dose
# after filters (2^n in {8192, 65536}), so a full enumeration is fast and
# gives a reproducible p-value with no Monte Carlo error.
perm_chan <- data.frame(numStims = character(), nUnits = integer(),
  obs_mean = numeric(), perm_p = numeric(), stringsAsFactors = FALSE)

for (dose in dose_levels_paired) {
  dSub <- chan_cond_means[chan_cond_means$numStims == dose, ]
  ef <- exact_signflip(dSub$mean_diff, stat_fn = mean)
  perm_chan <- rbind(perm_chan,
    data.frame(numStims = dose, nUnits = length(dSub$mean_diff),
               obs_mean = round(ef$obs, 2),
               perm_p = round(ef$p_two_sided, 4)))
}
perm_chan$perm_p_holm <- p.adjust(perm_chan$perm_p, method = "holm")

cat("\nChannel x condition sign-flip permutation (Holm-corrected):\n")
print(perm_chan)

# --- Per channel x condition sign-flip tests ---
# Each (channel, phase condition) tested separately at each dose.
# Uses median as test statistic. Reported descriptively.
cat("\n--- Per channel x condition sign-flip tests ---\n")
perm_list <- list(); idx <- 0L

for (ch in sort(unique(paired$channel_raw))) {
  for (phrd in sort(unique(paired$phaseDeg_round[paired$channel_raw == ch]))) {
    for (dose in dose_levels_paired) {
      dSub <- paired[paired$channel_raw == ch & paired$phaseDeg_round == phrd &
                     paired$numStims == dose, ]
      if (nrow(dSub) < 5) next
      obs_med <- median(dSub$diff)
      n <- nrow(dSub)
      perm_meds <- replicate(nPerm, {
        signs <- sample(c(-1, 1), n, replace = TRUE)
        median(dSub$diff * signs)
      })
      idx <- idx + 1L
      perm_list[[idx]] <- data.frame(channel_raw = as.character(ch), phaseDeg_round = phrd,
                   numStims = dose, n_probes = n,
                   obs_median_diff = round(obs_med, 1),
                   perm_p = mean(abs(perm_meds) >= abs(obs_med)),
                   cl_phase = round(mean(dSub$cl_burst_phase, na.rm = TRUE), 1),
                   phaseVecLength = round(dSub$phaseVecLength[1], 2),
                   setToDeliverPhase = as.character(dSub$setToDeliverPhase[1]))
    }
  }
}
perm_perchan <- do.call(rbind, perm_list)

# --- BH FDR within dose ---
# Each dose is a distinct scientific question (does CL exceed PB at this dose?).
# ~13-16 channel x phase cells per dose form a natural family. Mirrors the
# conditioned-vs-baseline per-cell correction at line 1511. Only one subject
# (0b5a2e) contributes, so no subject-level nesting is needed.
perm_perchan$perm_q_bh <- NA_real_
for (dose in dose_levels_paired) {
  rows <- perm_perchan$numStims == dose
  if (sum(rows) > 0) {
    perm_perchan$perm_q_bh[rows] <- p.adjust(perm_perchan$perm_p[rows], method = "BH")
  }
}
perm_perchan$perm_q_bh <- round(perm_perchan$perm_q_bh, 4)

cat(sprintf("Per channel x condition tests: %d cells with >= 5 probes\n", nrow(perm_perchan)))
print(perm_perchan[order(perm_perchan$numStims, perm_perchan$cl_phase), ])

cat("\nUnits with CL > PB (median diff > 0) per dose:\n")
for (dose in dose_levels_paired) {
  sub <- perm_perchan[perm_perchan$numStims == dose, ]
  cat(sprintf("  %s: %d/%d positive (%.0f%%), sig at p<0.05 (uncorr): %d, BH q<0.05: %d\n",
      dose, sum(sub$obs_median_diff > 0), nrow(sub),
      100 * mean(sub$obs_median_diff > 0),
      sum(sub$perm_p < 0.05), sum(sub$perm_q_bh < 0.05, na.rm = TRUE)))
}

# --- Direction by measured-phase region (depol 45-135 vs hyperpol 225-315) ---
# Bin cells by where their channel-level measured phase landed: the 90°
# (depolarizing, 45-135) region vs the 270° (hyperpolarizing, 225-315) region.
# Cells outside both bands are tracked as "other" (near 0° or 180°).
# For each dose, report count positive / count negative in each region, median
# effect, and FDR-significant-cell counts. This complements the per-cell forest
# plot by summarizing the direction of modulation across cells within each
# physiologically-meaningful phase neighborhood.
cat("\nDirection by measured-phase half (hyperpol half 0-180 deg vs depol half 180-360 deg):\n")
# phaseDeg_round is the channel-level circular mean (already in [0, 360) from
# circ_mean + round), so a direct band check suffices.
# Convention (per Zanos et al., Curr Biol 2018, PIIS0960982218309084):
# 90 deg = hyperpolarizing phase of the beta oscillation (center of 0-180 half)
# 270 deg = depolarizing phase of the beta oscillation (center of 180-360 half).
perm_perchan$phase_region <- factor(
  ifelse(perm_perchan$phaseDeg_round >= 0 & perm_perchan$phaseDeg_round < 180,
         "hyperpol half (0-180)", "depol half (180-360)"),
  levels = c("hyperpol half (0-180)", "depol half (180-360)"))

clpb_region_direction <- plyr::ddply(
  perm_perchan, .(numStims, phase_region), summarize,
  n_cells        = length(obs_median_diff),
  n_positive     = sum(obs_median_diff > 0),
  n_negative     = sum(obs_median_diff < 0),
  pct_positive   = round(100 * mean(obs_median_diff > 0), 1),
  median_effect  = round(median(obs_median_diff), 1),
  mean_effect    = round(mean(obs_median_diff), 1),
  n_sig_uncorr   = sum(!is.na(perm_p) & perm_p < 0.05),
  n_sig_fdr_bh   = sum(!is.na(perm_q_bh) & perm_q_bh < 0.05))
clpb_region_direction$numStims <- factor(
  clpb_region_direction$numStims, levels = dose_levels_paired)
clpb_region_direction <- clpb_region_direction[
  order(clpb_region_direction$numStims, clpb_region_direction$phase_region), ]
print(clpb_region_direction, row.names = FALSE)

# --- Dot plot + loess: diff vs CL per-burst phase, faceted by dose ---
paired$numStims <- factor(paired$numStims, levels = dose_levels_paired)

p_paired_dots <- ggplot(paired, aes(x = cl_burst_phase, y = diff)) +
  theme_light(base_size = 14) +
  facet_wrap(~ numStims, ncol = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = TRUE, color = "red", linewidth = 0.8) +
  labs(x = "CL Per-Burst Delivered Phase (degrees)",
       y = expression(paste(Delta, " Baseline-Normalized: CL - PB (", mu, "V)")),
       title = "Paired Probe Differences vs CL Burst Phase by Dose") +
  scale_x_continuous(breaks = seq(0, 315, by = 90))
p_paired_dots

if(savePlot){
  ggsave(here("output_plots","betaStim_clpb_paired_diff_vs_phase.png"), plot = p_paired_dots,
         units = "in", width = 10, height = 4, dpi = 600)
}

# --- Density of differences per dose ---
p_paired_density <- ggplot(paired, aes(x = diff)) +
  theme_light(base_size = 14) +
  facet_wrap(~ numStims, ncol = 3) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(data = data.frame(numStims = factor(perm_chan$numStims, levels = dose_levels_paired),
                                obs_mean = perm_chan$obs_mean),
             aes(xintercept = obs_mean), color = "red", linewidth = 0.8) +
  labs(x = expression(paste(Delta, " Baseline-Normalized: CL - PB (", mu, "V)")),
       y = "Density",
       title = "Distribution of Paired Differences (red = observed mean)")
p_paired_density

if(savePlot){
  ggsave(here("output_plots","betaStim_clpb_paired_diff_density.png"), plot = p_paired_density,
         units = "in", width = 10, height = 4, dpi = 600)
}

# --- Per channel x condition: histogram + density, faceted by dose ---
paired$chan_cond_label <- paste0("Ch ", paired$channel_raw, " (", paired$phaseDeg_round, " deg)")
chan_cond_order <- unique(paired[order(paired$phaseDeg_round), "chan_cond_label"])
paired$chan_cond_label <- factor(paired$chan_cond_label, levels = chan_cond_order)

perm_perchan$chan_cond_label <- paste0("Ch ", perm_perchan$channel_raw,
  " (", perm_perchan$phaseDeg_round, " deg)")
perm_perchan$chan_cond_label <- factor(perm_perchan$chan_cond_label, levels = chan_cond_order)
perm_perchan$numStims_f <- factor(perm_perchan$numStims, levels = dose_levels_paired)
perm_perchan$p_label <- ifelse(perm_perchan$perm_p < 0.001, "p<0.001",
  sprintf("p=%.3f", perm_perchan$perm_p))

p_perchan <- ggplot(paired, aes(x = diff)) +
  theme_light(base_size = 9) +
  facet_grid(chan_cond_label ~ numStims, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_histogram(aes(y = after_stat(density)), bins = 25, fill = "steelblue", alpha = 0.4) +
  geom_density(color = "darkblue", linewidth = 0.4) +
  geom_vline(data = perm_perchan, aes(xintercept = obs_median_diff),
             color = "red", linewidth = 0.5) +
  geom_text(data = perm_perchan, aes(label = p_label, x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.5, size = 2, color = "red") +
  labs(x = expression(paste(Delta, " CL - PB (", mu, "V)")),
       y = "Density",
       title = "Per Channel x Condition: CL-PB Differences by Dose") +
  coord_cartesian(xlim = c(-500, 500))
p_perchan

if(savePlot){
  ggsave(here("output_plots","betaStim_clpb_perchan_diff.png"), plot = p_perchan,
         units = "in", width = 10, height = 16, dpi = 600)
}

# --- Scatter: median diff vs CL phase at [5,inf) ---
# One dot per (channel x condition). Shows phase-response of CL advantage.
# Use channel-level condition phase (phaseDeg_round) for x-axis, not per-burst
# cl_phase. phaseDeg_round is the stable circular mean; per-burst phases are
# noisy and can compress two well-separated conditions toward the middle.
perm_5inf <- perm_perchan[perm_perchan$numStims == "[5,inf)", ]
perm_5inf$sig <- perm_5inf$perm_p < 0.05

p_phase_diff <- ggplot(perm_5inf, aes(x = phaseDeg_round, y = obs_median_diff)) +
  theme_light(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "loess", se = TRUE, color = "grey70", linewidth = 0.5) +
  geom_point(aes(color = sig), size = 4) +
  geom_text(aes(label = channel_raw), vjust = -1, size = 3) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red"),
                     labels = c("ns", "p<0.05"), name = "") +
  labs(x = "CL Delivered Phase (degrees)",
       y = expression(paste("Median CL - PB Difference (", mu, "V)")),
       title = "Phase-Response of CL Advantage at [5,inf) Dose") +
  scale_x_continuous(breaks = seq(0, 315, by = 45))
p_phase_diff

if(savePlot){
  ggsave(here("output_plots","betaStim_clpb_phase_response_5inf.png"), plot = p_phase_diff,
         units = "in", width = 7, height = 4.5, dpi = 600)
}

# ========================================================================
# ALL-BURSTS version: same analysis without good-fit filter
# ========================================================================
# Includes ALL matched probes (not just nGoodBeta > 0). Uses channel-level
# phase (phaseDeg_round) for grouping and x-axis. Tests whether the phase-
# response pattern is visible even without restricting to good fits.
cat("\n========== All-bursts paired analysis (no good-fit filter) ==========\n")

# build paired data from ALL CL and PB trials (no nGoodBeta filter)
dataCL_all <- dataCL_raw[dataCL_raw$magnitude > 25 & dataCL_raw$magnitude < 1500, ]
dataPB_all <- dataPB_raw[dataPB_raw$magnitude > 25 & dataPB_raw$magnitude < 1500, ]

dataCL_all_p <- dataCL_all[, c("channel_raw", "probeIdx", "magnitude", "numStims", "phaseDeg_round")]
names(dataCL_all_p)[3] <- "mag_CL"
dataPB_all_p <- dataPB_all[, c("channel_raw", "probeIdx", "magnitude")]
names(dataPB_all_p)[3] <- "mag_PB"

paired_all <- merge(dataCL_all_p, dataPB_all_p, by = c("channel_raw", "probeIdx"))
paired_all <- merge(paired_all, baseCL_chan, by = "channel_raw")
paired_all <- merge(paired_all, basePB_chan, by = "channel_raw")
paired_all$norm_CL <- paired_all$mag_CL - paired_all$baseCL
paired_all$norm_PB <- paired_all$mag_PB - paired_all$basePB
paired_all$diff <- paired_all$norm_CL - paired_all$norm_PB
paired_all <- paired_all[paired_all$numStims != "Base" & paired_all$numStims != "Null", ]

cat(sprintf("All-bursts matched pairs: %d\n", nrow(paired_all)))

# per channel x condition sign-flip tests
perm_list_all <- list(); idx_all <- 0L
for (ch in sort(unique(paired_all$channel_raw))) {
  for (phrd in sort(unique(paired_all$phaseDeg_round[paired_all$channel_raw == ch]))) {
    for (dose in dose_levels_paired) {
      dSub <- paired_all[paired_all$channel_raw == ch & paired_all$phaseDeg_round == phrd &
                         paired_all$numStims == dose, ]
      if (nrow(dSub) < 5) next
      obs_med <- median(dSub$diff)
      n <- nrow(dSub)
      perm_meds <- replicate(nPerm, {
        signs <- sample(c(-1, 1), n, replace = TRUE)
        median(dSub$diff * signs)
      })
      idx_all <- idx_all + 1L
      perm_list_all[[idx_all]] <- data.frame(channel_raw = as.character(ch), phaseDeg_round = phrd,
                   numStims = dose, n_probes = n,
                   obs_median_diff = round(obs_med, 1),
                   perm_p = mean(abs(perm_meds) >= abs(obs_med)))
    }
  }
}
perm_perchan_all <- do.call(rbind, perm_list_all)

cat(sprintf("All-bursts per channel x condition: %d cells\n", nrow(perm_perchan_all)))

# scatter: all-bursts vs good-fit at [5,inf)
perm_5inf_all <- perm_perchan_all[perm_perchan_all$numStims == "[5,inf)", ]
perm_5inf_all$sig <- perm_5inf_all$perm_p < 0.05
perm_5inf_all$filter <- "All bursts"
perm_5inf_gf <- perm_5inf
perm_5inf_gf$filter <- "Good-fit only"
perm_5inf_both <- rbind(
  perm_5inf_all[, c("channel_raw", "phaseDeg_round", "obs_median_diff", "perm_p", "sig", "filter")],
  perm_5inf_gf[, c("channel_raw", "phaseDeg_round", "obs_median_diff", "perm_p", "sig", "filter")])

p_phase_both <- ggplot(perm_5inf_both, aes(x = phaseDeg_round, y = obs_median_diff)) +
  theme_light(base_size = 14) +
  facet_wrap(~ filter) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "loess", se = TRUE, color = "grey70", linewidth = 0.5) +
  geom_point(aes(color = sig), size = 4) +
  geom_text(aes(label = channel_raw), vjust = -1, size = 3) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red"),
                     labels = c("ns", "p<0.05"), name = "") +
  labs(x = "CL Delivered Phase (degrees)",
       y = expression(paste("Median CL - PB Difference (", mu, "V)")),
       title = "Phase-Response at [5,inf): Good-fit vs All Bursts") +
  scale_x_continuous(breaks = seq(0, 315, by = 45))
p_phase_both

if(savePlot){
  ggsave(here("output_plots","betaStim_clpb_phase_response_5inf_comparison.png"), plot = p_phase_both,
         units = "in", width = 12, height = 4.5, dpi = 600)
}

paired$numStims <- factor(paired$numStims, levels = dose_levels_paired)

# ========================================================================
# PRIMARY channel-level figure (Option A): 8 channels × 3 doses grid.
# All channels on one image; facet_grid(channel_raw ~ numStims). Each dot
# is a matched probe pair; blue horizontal line = per-cell median; dashed
# line at 0. Replaces the eight separate betaStim_clpb_ch*_diff_vs_phase
# images as the primary manuscript figure (those are retained as
# supplementary below).
# ========================================================================
cell_medians <- plyr::ddply(paired, .(channel_raw, numStims), summarize,
  median_diff = median(diff))

p_clpb_grid <- ggplot(paired, aes(x = cl_burst_phase, y = diff)) +
  theme_light(base_size = 11) +
  facet_grid(channel_raw ~ numStims, scales = "fixed", switch = "y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(alpha = 0.35, size = 0.8, color = "grey30") +
  geom_hline(data = cell_medians, aes(yintercept = median_diff),
             color = "#d62728", linewidth = 0.6) +
  geom_smooth(method = "loess", se = FALSE, span = 0.9,
              color = "#1f77b4", linewidth = 0.5) +
  labs(x = "CL Per-Burst Delivered Phase (degrees)",
       y = expression(paste("Channel    ",
                             Delta, " CL - PB (", mu, "V)")),
       title = "0b5a2e: Matched CL - PB Differences by Channel × Dose",
       subtitle = sprintf("%d probe pairs across 8 channels × 3 doses (red = per-cell median, blue = loess)",
                          nrow(paired))) +
  scale_x_continuous(breaks = c(0, 90, 180, 270), limits = c(0, 360)) +
  theme(strip.text.y.left = element_text(angle = 0),
        panel.spacing.x = unit(0.4, "lines"),
        panel.spacing.y = unit(0.2, "lines"))

if (savePlot) {
  ggsave(here("output_plots", "betaStim_clpb_grid_ch_x_dose.png"),
         plot = p_clpb_grid,
         units = "in", width = 10, height = 12, dpi = 600)
  ggsave(here("output_plots", "betaStim_clpb_grid_ch_x_dose.eps"),
         plot = p_clpb_grid,
         units = "in", width = 10, height = 12, dpi = 600, device = cairo_ps)
}

# ========================================================================
# PRIMARY inferential figure (Option B): forest plot of per-cell effect
# sizes with bootstrap 95% CIs, significance from the sign-flip perms.
# One row per (channel × phase condition), three columns (doses). Shows
# direction, magnitude, and significance in a single compact frame. Pairs
# directly with the matched-pair permutation tests (perm_perchan) that
# are the primary inferential claim for 0b5a2e.
# ========================================================================
set.seed(42)
nBoot <- 2000
perm_perchan$lo_ci <- NA_real_
perm_perchan$hi_ci <- NA_real_
for (i in seq_len(nrow(perm_perchan))) {
  row <- perm_perchan[i, ]
  d <- paired$diff[paired$channel_raw == row$channel_raw &
                   paired$phaseDeg_round == row$phaseDeg_round &
                   paired$numStims == row$numStims]
  if (length(d) < 3) next
  boot_meds <- replicate(nBoot, median(sample(d, length(d), replace = TRUE)))
  perm_perchan$lo_ci[i] <- quantile(boot_meds, 0.025, names = FALSE)
  perm_perchan$hi_ci[i] <- quantile(boot_meds, 0.975, names = FALSE)
}

perm_perchan$chan_cond_lbl <- sprintf("Ch %s @ %.0f° (r=%.2f) [tgt %s°]",
                                       as.character(perm_perchan$channel_raw),
                                       perm_perchan$phaseDeg_round,
                                       perm_perchan$phaseVecLength,
                                       perm_perchan$setToDeliverPhase)
# Sort rows by measured delivered phase (0° at top, 360° at bottom),
# channel ascending as secondary tie-breaker. rev() is needed because
# ggplot's discrete y-axis draws the FIRST level at the bottom and the
# LAST level at the top.
row_order <- unique(perm_perchan[order(perm_perchan$phaseDeg_round,
                                        perm_perchan$channel_raw),
                                  "chan_cond_lbl"])
perm_perchan$chan_cond_lbl <- factor(perm_perchan$chan_cond_lbl,
                                      levels = rev(row_order))
perm_perchan$sig_uncorr <- !is.na(perm_perchan$perm_p) &
                            perm_perchan$perm_p < 0.05
perm_perchan$sig_fdr    <- !is.na(perm_perchan$perm_q_bh) &
                            perm_perchan$perm_q_bh < 0.05
perm_perchan$numStims_f <- factor(perm_perchan$numStims, levels = dose_levels_paired)

# Three-tier significance based on BH FDR within dose. Dot color reflects
# significance only (no beta-channel override) so significance is readable
# on every row. Beta trigger channel (Ch 31) is flagged separately by (a)
# pink y-axis tick labels and (b) a magenta ring behind the dot, matching
# the conditioned-vs-baseline forest plot convention.
perm_perchan$color_cat <- factor(
  ifelse(perm_perchan$sig_fdr, "FDR q<0.05 (within dose)",
    ifelse(perm_perchan$sig_uncorr, "p<0.05 uncorr", "ns")),
  levels = c("ns", "p<0.05 uncorr", "FDR q<0.05 (within dose)"))

# axis.text.y color vector: element_text accepts a vector of per-tick
# colors, applied in factor-level order (which for our factor is
# bottom-to-top = descending measured phase).
y_label_colors <- ifelse(
  grepl("^Ch 31 ", levels(perm_perchan$chan_cond_lbl)),
  "deeppink3", "black")

# Subset used to draw a magenta ring behind beta-trigger (Ch 31) rows so the
# flag is visible independent of the significance-driven fill color. The
# ring is mapped to a named shape aesthetic so it appears in the legend.
beta_ring <- perm_perchan[as.character(perm_perchan$channel_raw) == "31", ]
beta_ring$ring_label <- "Beta trigger channel"

p_clpb_forest <- ggplot(perm_perchan,
  aes(y = chan_cond_lbl, x = obs_median_diff)) +
  theme_light(base_size = 12) +
  facet_wrap(~ numStims_f, nrow = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = lo_ci, xmax = hi_ci, color = color_cat),
                 height = 0.25, linewidth = 0.6, na.rm = TRUE) +
  # Magenta ring (open circle) behind Ch 31 dots. Draw first so the filled
  # significance dot sits on top. Mapped to a named shape scale so it
  # generates its own legend entry.
  geom_point(data = beta_ring, aes(shape = ring_label),
             size = 5, stroke = 1.2, color = "deeppink3", fill = NA) +
  geom_point(aes(color = color_cat), size = 2.8) +
  scale_color_manual(values = c("ns" = "grey45",
                                 "p<0.05 uncorr" = "#ff8c00",
                                 "FDR q<0.05 (within dose)" = "#d62728"),
                     drop = FALSE, name = "") +
  scale_shape_manual(values = c("Beta trigger channel" = 1), name = "") +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2,
                         override.aes = list(size = 5, stroke = 1.2,
                                             color = "deeppink3"))) +
  labs(x = expression(paste("Median Baseline-Normalized: CL - PB (", mu, "V)")),
       y = NULL,
       title = "0b5a2e: Paired CL - PB Effects by Channel × Phase × Dose",
       subtitle = "Each probe baseline-normalized to its own session. Bootstrap 95% CIs (2000 resamples); matched-pair sign-flip perms, BH FDR within dose.") +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(color = y_label_colors),
        legend.position = "bottom")

if (savePlot) {
  ggsave(here("output_plots", "betaStim_clpb_forest_paired_effects.png"),
         plot = p_clpb_forest,
         units = "in", width = 11, height = 6, dpi = 600)
  ggsave(here("output_plots", "betaStim_clpb_forest_paired_effects.eps"),
         plot = p_clpb_forest,
         units = "in", width = 11, height = 6, dpi = 600, device = cairo_ps)
}

# --- SUPPLEMENTARY: Per-channel individual probe scatter (eight images) ---
# One image per channel, each with 3 dose subpanels. Retained for
# supplementary detail; the primary channel-level figure is the grid
# (betaStim_clpb_grid_ch_x_dose.png) immediately above.
if(savePlot){
  for (ch in sort(unique(paired$channel_raw))) {
    ch_data <- paired[paired$channel_raw == ch, ]
    if (nrow(ch_data) < 10) next

    p_ch <- ggplot(ch_data, aes(x = cl_burst_phase, y = diff)) +
      theme_light(base_size = 14) +
      facet_wrap(~ numStims, ncol = 3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(alpha = 0.4, size = 1.5) +
      geom_smooth(method = "loess", se = TRUE, color = "red", linewidth = 0.8, span = 0.75) +
      labs(x = "CL Per-Burst Delivered Phase (degrees)",
           y = expression(paste(Delta, " CL - PB (", mu, "V)")),
           title = sprintf("Channel %s: Paired CL-PB Differences by Dose", as.character(ch))) +
      scale_x_continuous(breaks = seq(0, 315, by = 90))

    ggsave(here("output_plots", sprintf("betaStim_clpb_ch%s_diff_vs_phase.png", as.character(ch))),
           plot = p_ch, units = "in", width = 10, height = 4, dpi = 600)
  }
  cat(sprintf("Saved per-channel scatter plots for %d channels\n",
      length(unique(paired$channel_raw))))
}

# ========================================================================
# ecb43e: 3-condition (270, 90, random)
# ========================================================================
# 4 channels. condType ref = random. sin/cos at random = accidental phase
# clustering (should be null). condTypetargeted:sin/cos = targeted-specific
# phase effect beyond random.
# Only 270 has baselines — used for all conditions.
# Random: 4 channels x 1 phase x 3 dose = 12 cells
# Targeted: 4 channels x 2 conditions (270+90) x 3 dose = 24 cells
cat("\n========== ecb43e: targeted vs random ==========\n")

dataEC <- data[data$sid == "ecb43e" & data$numStims != "Null", ]
dataEC <- dataEC[dataEC$magnitude > 25 & dataEC$magnitude < 1500, ]
dataEC$condType <- factor(ifelse(dataEC$setToDeliverPhase == "12345", "random", "targeted"),
                          levels = c("random", "targeted"))

dataEC_NB <- dataEC[dataEC$numStims != "Base", ]
summaryEC <- ddply(dataEC_NB, .(condType, channel, phaseDeg_round, numStims),
  summarize, magnitude = median(magnitude),
  sin_phase = first(sin_phase), cos_phase = first(cos_phase))

baseEC <- ddply(dataEC[dataEC$numStims == "Base", ], .(channel),
  summarize, baselineMag = median(magnitude))
summaryEC <- merge(summaryEC, baseEC, by = "channel")
summaryEC$baselineMag_c <- summaryEC$baselineMag - mean(summaryEC$baselineMag)
summaryEC$numStims_ord <- ordered(summaryEC$numStims,
  levels = c("[1,2]", "[3,4]", "[5,inf)"))

cat(sprintf("ecb43e summary: %d obs, %d channels, targeted=%d, random=%d\n",
    nrow(summaryEC), length(unique(summaryEC$channel)),
    sum(summaryEC$condType == "targeted"), sum(summaryEC$condType == "random")))

cat("\n--- ecb43e combined model ---\n")

fit.ec = lmerTest::lmer(
  magnitude ~ numStims_ord * (sin_phase + cos_phase) +
              condType * (sin_phase + cos_phase) +
              baselineMag_c +
  (1 | channel),
  data = summaryEC,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

cat("Singular:", isSingular(fit.ec), "\n")
print(summary(fit.ec))
cat("\nType III ANOVA:\n")
print(anova(fit.ec))

cnames_ec <- names(fixef(fit.ec))
ct_sin <- cnames_ec[grepl("condType.*sin_phase|sin_phase.*condType", cnames_ec)]
ct_cos <- cnames_ec[grepl("condType.*cos_phase|cos_phase.*condType", cnames_ec)]
cat("\nWald: targeted-specific phase effect beyond random:\n")
print(car::linearHypothesis(fit.ec,
  c(paste0(ct_sin, " = 0"), paste0(ct_cos, " = 0"))))

cat("\nWald: phase at random reference (should be ns):\n")
print(car::linearHypothesis(fit.ec, c("sin_phase = 0", "cos_phase = 0")))

emm_ec_curves <- lapply(c("targeted", "random"), function(ct) {
  do.call(rbind, lapply(phase_vals_ctrl, function(ph) {
    em <- emmeans(fit.ec, ~ numStims_ord,
      at = list(sin_phase = sin(ph*pi/180), cos_phase = cos(ph*pi/180),
                condType = ct, baselineMag_c = 0))
    df <- as.data.frame(em)
    df$phase_deg <- ph; df$condType <- ct; df
  }))
})
emm_ec_df <- do.call(rbind, emm_ec_curves)

p_ec <- ggplot(emm_ec_df, aes(x = phase_deg, y = emmean, color = numStims_ord)) +
  theme_light(base_size = 14) + facet_wrap(~ condType) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = numStims_ord), alpha = 0.15, color = NA) +
  labs(x = "Delivered Phase (degrees)",
       y = expression(paste("Predicted Magnitude (", mu, "V)")),
       color = "Dose", fill = "Dose",
       title = "ecb43e: Targeted vs Random Phase-Response by Dose") +
  scale_x_continuous(breaks = seq(0, 315, by = 90))
p_ec
if(savePlot){
  ggsave(here("output_plots","betaStim_ecb43e_phase_curve.png"), plot = p_ec,
         units = "in", width = 10, height = 4.5, dpi = 600)
}

# ========================================================================
# Percent modulation from baseline — all subjects, per good channel
# ========================================================================
cat("\n========== Percent modulation from baseline (all subjects) ==========\n")

# use the full data (already filtered to magnitude 25-1500, no NaN)
dataAll <- data[data$numStims != "Null", ]
dataAll$numStims <- plyr::revalue(dataAll$numStims,
  c("Test 1"="[1,2]","Test 2"="[3,4]","Test 3"="[5,inf)"))

# baseline per (sid, channel) — pooled across conditions
baseAll <- plyr::ddply(dataAll[dataAll$numStims == "Base", ], .(sid, channel),
  summarize, baseMedian = median(magnitude))

# median magnitude per (sid, channel, numStims) — pooled across phase conditions
dataAll_NB <- dataAll[dataAll$numStims != "Base", ]
doseAll <- plyr::ddply(dataAll_NB, .(sid, channel, numStims),
  summarize, mag = median(magnitude))
doseAll <- merge(doseAll, baseAll, by = c("sid", "channel"))
doseAll$pctChange <- 100 * (doseAll$mag - doseAll$baseMedian) / doseAll$baseMedian

# per-subject summary: mean across channels for each dose
pct_subj <- plyr::ddply(doseAll, .(sid, numStims), summarize,
  nChannels = length(channel),
  meanBaseline = round(mean(baseMedian), 1),
  meanMag = round(mean(mag), 1),
  meanPctChange = round(mean(pctChange), 1),
  medianPctChange = round(median(pctChange), 1))

cat("\nPercent modulation per subject x dose:\n")
print(pct_subj)

# grand summary across subjects
pct_grand <- plyr::ddply(doseAll, .(numStims), summarize,
  nSubjects = length(unique(sid)),
  nChannels = length(channel),
  meanPctChange = round(mean(pctChange), 1),
  sdPctChange = round(sd(pctChange), 1),
  medianPctChange = round(median(pctChange), 1))

cat("\nGrand mean percent modulation by dose:\n")
print(pct_grand)

# per-channel detail table
pct_chan <- doseAll[, c("sid", "channel", "numStims", "baseMedian", "mag", "pctChange")]
pct_chan$baseMedian <- round(pct_chan$baseMedian, 1)
pct_chan$mag <- round(pct_chan$mag, 1)
pct_chan$pctChange <- round(pct_chan$pctChange, 1)
pct_chan <- pct_chan[order(pct_chan$sid, pct_chan$channel, pct_chan$numStims), ]

# ========================================================================
# Conditioned vs Baseline: per-cell permutation + within-subject FDR
# ========================================================================
# Companion to Model 5a-gf2. The LMM estimates the AVERAGE effect across
# channels. This block estimates how that average is DISTRIBUTED across
# channels — asking in how many cells the effect is individually detectable
# within each subject (channels are nested within subjects).
#
# Scope: all 7 main subjects (0b5a2ePlayBack and ecb43e random excluded).
# Filters: channel-level phaseVecLength >= 0.2 (matches 5a-gf2),
#          >= 10 baseline probes, >= 5 conditioned probes per cell.
# Test: two-sample label-shuffle permutation (median diff, 10k MC).
# Correction: BH FDR within (subject x dose). Channels are nested within
# subjects, so the FDR family is one subject's cells at one dose. A pooled
# FDR column (perm_q_pooled) is also reported as a sensitivity check.
# ========================================================================
cat("\n========== Conditioned vs Baseline per-cell permutation ==========\n")

minPhaseVecLength_cb    <- 0.2   # match 5a-gf2 channel r threshold
minGoodBetaPerBurst_cb  <- 1     # match 5a-gf2 good-fit filter: each conditioned
                                 # trial's preceding burst must have >= 1 stim
                                 # with R^2 > 0.7 AND frequency in 12-20 Hz.
                                 # Set 0 to disable for sensitivity checks.
min_n_cond_cb <- 5               # min conditioned probes per cell (post-filter)
min_n_base_cb <- 10              # min baseline probes per channel
nPerm_cb <- 10000
nBoot_cb <- 2000
set.seed(42)

# `data` already has magnitude bounds, ecb43e random excluded,
# 0b5a2ePlayBack excluded. Re-apply safety filters in case this block is
# run standalone.
dataCB <- data[data$magnitude >= 25 & data$magnitude <= 1500 &
               !is.na(data$magnitude) &
               data$sid != "0b5a2ePlayBack" &
               data$numStims != "Null" &
               !(data$sid == "ecb43e" &
                 as.character(data$setToDeliverPhase) == "12345"), ]

# Merge nGoodBeta from per-subject burst phase precision CSVs. The merge
# key is (probeSample, channel) -> (probeSample, channelEncoded). Matches
# how Model 5a-gf2 applies its good-fit filter. Baseline trials don't have
# preceding bursts, so they keep nGoodBeta = NA and are exempted from the
# filter below.
precision_files_cb <- Sys.glob(here("data", "output_table",
                                     "*_burst_phase_precision.csv"))
if (length(precision_files_cb) > 0) {
  precision_cb <- do.call(rbind, lapply(precision_files_cb, read.csv))
  precision_cb_sub <- precision_cb[, c("probeSample", "channelEncoded",
                                        "nGoodBeta")]
  # Drop 0b5a2ePlayBack rows — only the 7 main subjects enter this analysis.
  dataCB <- merge(dataCB, precision_cb_sub,
                  by.x = c("probeSample", "channel"),
                  by.y = c("probeSample", "channelEncoded"),
                  all.x = TRUE)
  n_before_gf <- sum(dataCB$numStims != "Base")
  # Apply good-fit filter to conditioned trials only. Trials with no
  # precision match (NA nGoodBeta) on a conditioned row are dropped —
  # matches 5a-gf2 behavior of requiring a confirmed burst fit.
  keep_trial <- dataCB$numStims == "Base" |
    (!is.na(dataCB$nGoodBeta) & dataCB$nGoodBeta >= minGoodBetaPerBurst_cb)
  dataCB <- dataCB[keep_trial, ]
  n_after_gf <- sum(dataCB$numStims != "Base")
  cat(sprintf("Good-fit filter (nGoodBeta >= %d): kept %d of %d conditioned trials (%.1f%%)\n",
      minGoodBetaPerBurst_cb, n_after_gf, n_before_gf,
      100 * n_after_gf / n_before_gf))
} else {
  warning("No burst_phase_precision CSVs found — good-fit filter not applied. Results will be non-comparable to 5a-gf2.")
  dataCB$nGoodBeta <- NA
}

# Channel × condition combos passing the phaseVecLength filter.
keep_chan_cond <- unique(dataCB[!is.na(dataCB$phaseVecLength) &
  dataCB$phaseVecLength >= minPhaseVecLength_cb,
  c("sid", "channel", "setToDeliverPhase")])
cat(sprintf("Kept %d channel x condition combos at phaseVecLength >= %.2f\n",
    nrow(keep_chan_cond), minPhaseVecLength_cb))

doses_cb <- c("[1,2]", "[3,4]", "[5,inf)")
cell_list <- list(); idx_cb <- 0L

for (sid_val in sort(unique(as.character(dataCB$sid)))) {
  dS <- dataCB[as.character(dataCB$sid) == sid_val, ]
  for (ch in sort(unique(as.character(dS$channel)))) {
    dCh <- dS[as.character(dS$channel) == ch, ]
    # Baselines are channel-level (not tied to phase condition).
    mBase <- dCh$magnitude[dCh$numStims == "Base"]
    if (length(mBase) < min_n_base_cb) next

    for (phrd in sort(unique(dCh$phaseDeg_round))) {
      if (is.na(phrd)) next
      setPh <- unique(as.character(
        dCh$setToDeliverPhase[dCh$phaseDeg_round == phrd]))
      passes <- any(keep_chan_cond$sid == sid_val &
                    as.character(keep_chan_cond$channel) == ch &
                    as.character(keep_chan_cond$setToDeliverPhase) %in% setPh)
      if (!passes) next

      for (dose in doses_cb) {
        mCond <- dCh$magnitude[dCh$phaseDeg_round == phrd &
                                as.character(dCh$numStims) == dose]
        if (length(mCond) < min_n_cond_cb) next

        obs_diff <- median(mCond) - median(mBase)
        all_mags <- c(mCond, mBase)
        is_cond  <- c(rep(TRUE, length(mCond)), rep(FALSE, length(mBase)))
        perm_diffs <- replicate(nPerm_cb, {
          shuf <- sample(is_cond)
          median(all_mags[shuf]) - median(all_mags[!shuf])
        })
        perm_p <- mean(abs(perm_diffs) >= abs(obs_diff))

        boot_diffs <- replicate(nBoot_cb, {
          bC <- sample(mCond, replace = TRUE)
          bB <- sample(mBase, replace = TRUE)
          median(bC) - median(bB)
        })

        idx_cb <- idx_cb + 1L
        cell_list[[idx_cb]] <- data.frame(
          sid = sid_val, channel_raw = ch,
          phaseDeg_round = phrd,
          setToDeliverPhase = paste(setPh, collapse = "/"),
          numStims = dose,
          n_cond = length(mCond), n_base = length(mBase),
          median_cond = round(median(mCond), 1),
          median_base = round(median(mBase), 1),
          obs_diff = round(obs_diff, 1),
          lo_ci = round(quantile(boot_diffs, 0.025, names = FALSE), 1),
          hi_ci = round(quantile(boot_diffs, 0.975, names = FALSE), 1),
          perm_p = round(perm_p, 4),
          stringsAsFactors = FALSE)
      }
    }
  }
}
cb_perchan <- do.call(rbind, cell_list)
cat(sprintf("Built %d cells across %d subjects\n",
    nrow(cb_perchan), length(unique(cb_perchan$sid))))

# --- FDR within (subject x dose): primary correction ---
# Family = one subject's cells at one dose. Matches the nested design.
cb_perchan$perm_q <- NA_real_
cb_perchan$n_tests_family <- NA_integer_
for (sid_val in unique(cb_perchan$sid)) {
  for (dose in doses_cb) {
    rows <- cb_perchan$sid == sid_val & cb_perchan$numStims == dose
    if (sum(rows) > 0) {
      cb_perchan$perm_q[rows] <- p.adjust(cb_perchan$perm_p[rows], method = "BH")
      cb_perchan$n_tests_family[rows] <- sum(rows)
    }
  }
}
cb_perchan$perm_q <- round(cb_perchan$perm_q, 4)

# --- Pooled FDR within dose: sensitivity check ---
cb_perchan$perm_q_pooled <- NA_real_
for (dose in doses_cb) {
  rows <- cb_perchan$numStims == dose
  if (sum(rows) > 0) {
    cb_perchan$perm_q_pooled[rows] <-
      p.adjust(cb_perchan$perm_p[rows], method = "BH")
  }
}
cb_perchan$perm_q_pooled <- round(cb_perchan$perm_q_pooled, 4)

cb_perchan$sig_uncorr    <- cb_perchan$perm_p < 0.05
cb_perchan$sig_fdr       <- !is.na(cb_perchan$perm_q) &
                            cb_perchan$perm_q < 0.05
cb_perchan$sig_fdr_pooled <- !is.na(cb_perchan$perm_q_pooled) &
                             cb_perchan$perm_q_pooled < 0.05

# --- Per-subject x dose: primary reporting unit ---
cb_subj_fdr <- plyr::ddply(cb_perchan, .(sid, numStims), summarize,
  n_cells          = length(perm_p),
  n_sig_uncorr     = sum(sig_uncorr),
  n_sig_fdr        = sum(sig_fdr),
  pct_sig_fdr      = round(100 * mean(sig_fdr), 1),
  any_sig_uncorr   = any(sig_uncorr),
  any_sig_fdr      = any(sig_fdr),
  median_effect    = round(median(obs_diff), 1))
cb_subj_fdr$numStims <- factor(cb_subj_fdr$numStims, levels = doses_cb)
cb_subj_fdr <- cb_subj_fdr[order(cb_subj_fdr$sid, cb_subj_fdr$numStims), ]

cat("\nPer-subject x dose (FDR within subject x dose):\n")
print(cb_subj_fdr)

# --- Across-subject summary per dose ---
# Headline statistic: median of per-subject fractions modulated, and count
# of subjects with >=1 FDR-sig cell. Treats each subject as one unit.
cb_subj_presence <- plyr::ddply(cb_subj_fdr, .(numStims), summarize,
  n_subjects_total      = length(any_sig_fdr),
  n_subj_any_fdr        = sum(any_sig_fdr),
  pct_subj_any_fdr      = round(100 * mean(any_sig_fdr), 1),
  n_subj_any_uncorr     = sum(any_sig_uncorr),
  pct_subj_any_uncorr   = round(100 * mean(any_sig_uncorr), 1),
  median_pct_sig_fdr    = round(median(pct_sig_fdr), 1),
  mean_pct_sig_fdr      = round(mean(pct_sig_fdr), 1))
cb_subj_presence$numStims <- factor(cb_subj_presence$numStims, levels = doses_cb)
cb_subj_presence <- cb_subj_presence[order(cb_subj_presence$numStims), ]

cat("\nAcross-subject summary per dose:\n")
print(cb_subj_presence)

# --- Per-dose summary: collapses across subjects, reports within-subject FDR ---
# One row per dose bin. Counts cells (not subjects) passing uncorrected and
# within-subject FDR thresholds. The within-subject FDR is the primary
# correction — families built per (subject x dose) — and this summary
# aggregates the per-family outcomes into one global count per dose.
cb_summary <- plyr::ddply(cb_perchan, .(numStims), summarize,
  n_cells         = length(perm_p),
  n_subjects      = length(unique(sid)),
  n_channels      = length(unique(paste(sid, channel_raw))),
  n_sig_uncorr    = sum(sig_uncorr),
  pct_sig_uncorr  = round(100 * mean(sig_uncorr), 1),
  n_sig_fdr       = sum(sig_fdr),
  pct_sig_fdr     = round(100 * mean(sig_fdr), 1),
  median_effect   = round(median(obs_diff), 1),
  mean_effect     = round(mean(obs_diff), 1))
cb_summary$numStims <- factor(cb_summary$numStims, levels = doses_cb)
cb_summary <- cb_summary[order(cb_summary$numStims), ]

cat("\nPer-dose summary (within-subject FDR, collapsed across subjects):\n")
print(cb_summary)

# --- Pooled (secondary) summary for sensitivity check ---
cb_summary_pooled <- plyr::ddply(cb_perchan, .(numStims), summarize,
  n_cells         = length(perm_p),
  n_sig_uncorr    = sum(sig_uncorr),
  pct_sig_uncorr  = round(100 * mean(sig_uncorr), 1),
  n_sig_fdr_pooled  = sum(sig_fdr_pooled),
  pct_sig_fdr_pooled = round(100 * mean(sig_fdr_pooled), 1))
cb_summary_pooled$numStims <- factor(cb_summary_pooled$numStims, levels = doses_cb)
cb_summary_pooled <- cb_summary_pooled[order(cb_summary_pooled$numStims), ]

cat("\nPooled-FDR sensitivity summary per dose:\n")
print(cb_summary_pooled)

# --- Direction by measured-phase half across all subjects ---
# Analogous to the clpb_region_direction block for CL vs PB, but across all
# 7 subjects (101 cells total). Splits measured delivered phase into two
# halves per Zanos et al. (Curr Biol 2018, PIIS0960982218309084) convention:
#   hyperpol (0-180 deg), depol (180-360 deg). 90 deg is center of hyperpol
#   half; 270 deg is center of depol half.
# Reports per (dose x half): cell counts, pos/neg breakdown, median/mean
# effect, and within-subject FDR-significant counts.
cat("\nCond vs Base direction by measured-phase half (all subjects):\n")
cb_perchan$phase_half <- factor(
  ifelse(cb_perchan$phaseDeg_round >= 0 & cb_perchan$phaseDeg_round < 180,
         "hyperpol half (0-180)", "depol half (180-360)"),
  levels = c("hyperpol half (0-180)", "depol half (180-360)"))

cb_halves <- plyr::ddply(cb_perchan, .(numStims, phase_half), summarize,
  n_cells       = length(obs_diff),
  n_subjects    = length(unique(sid)),
  n_channels    = length(unique(paste(sid, channel_raw))),
  n_positive    = sum(obs_diff > 0),
  n_negative    = sum(obs_diff < 0),
  pct_positive  = round(100 * mean(obs_diff > 0), 1),
  median_effect = round(median(obs_diff), 1),
  mean_effect   = round(mean(obs_diff), 1),
  n_sig_uncorr  = sum(sig_uncorr),
  n_sig_fdr     = sum(sig_fdr))
cb_halves$numStims <- factor(cb_halves$numStims, levels = doses_cb)
cb_halves <- cb_halves[order(cb_halves$numStims, cb_halves$phase_half), ]
print(cb_halves, row.names = FALSE)

# --- Overlap with 5a-gf2 channels (internal consistency check) ---
if (exists("summaryNB_gf2")) {
  gf2_keys <- unique(paste(summaryNB_gf2$sid, summaryNB_gf2$channel,
                           summaryNB_gf2$phaseDeg_round, sep = "|"))
  cb_perchan$in_5agf2 <- paste(cb_perchan$sid, cb_perchan$channel_raw,
                                cb_perchan$phaseDeg_round, sep = "|") %in% gf2_keys
  sig_hi <- cb_perchan[cb_perchan$sig_fdr & cb_perchan$numStims == "[5,inf)", ]
  cat(sprintf("\n[5,inf) FDR-sig cells also in 5a-gf2 good-fit subset: %d of %d\n",
      sum(sig_hi$in_5agf2), nrow(sig_hi)))
}

# --- Forest plot: all subjects, faceted by dose ---
# Colored by within-subject FDR (primary; channels nested within subjects).
# Pooled FDR is retained as a CSV column (perm_q_pooled) but not plotted —
# on this dataset it was near-identical to within-subject and added no info.
# Rows sorted by measured phase (0 at top), with subject as tie-breaker.
# Channel labels strip the subject-number prefix (e.g., 714 -> 14).
# Beta trigger channels highlighted via pink y-axis tick labels.

cb_perchan$numStims_f <- factor(cb_perchan$numStims, levels = doses_cb)
# Raw channel number (strip subjectNum*100 prefix)
cb_perchan$channel_disp <- as.integer(as.character(cb_perchan$channel_raw)) %% 100
# Coded sid -> "Subject N" label using the subjectNum column carried in `data`
sid_num_map <- setNames(
  as.integer(as.character(unique(data[, c("sid", "subjectNum")])$subjectNum)),
  as.character(unique(data[, c("sid", "subjectNum")])$sid))
cb_perchan$subj_label <- sprintf("Subject %d",
  sid_num_map[as.character(cb_perchan$sid)])
cb_perchan$cell_label <- sprintf("%s Ch%d @ %.0f\u00b0",
                                  cb_perchan$subj_label,
                                  cb_perchan$channel_disp,
                                  cb_perchan$phaseDeg_round)

# Sort by phase ascending (0 at top); ggplot draws first factor level at the
# bottom, so rev() puts the smallest phase at the top of the y-axis.
row_order <- unique(cb_perchan[
  order(cb_perchan$phaseDeg_round, as.character(cb_perchan$sid)),
  "cell_label"])
cb_perchan$cell_label <- factor(cb_perchan$cell_label, levels = rev(row_order))

# --- Flag beta-trigger channels for highlighting (pink override) ---
# Beta reference channels per subject (from CLAUDE.md table).
beta_ref_map <- data.frame(
  sid      = c("d5cd55","c91479","7dbdec","9ab7ab","702d24","ecb43e","0b5a2e"),
  beta_raw = c(53,       64,      4,       51,      5,       55,      31),
  stringsAsFactors = FALSE)
cb_perchan$is_beta_ref <- mapply(function(s, c_raw) {
  row <- beta_ref_map[beta_ref_map$sid == s, ]
  if (nrow(row) == 0) return(FALSE)
  c_raw == row$beta_raw
}, as.character(cb_perchan$sid), cb_perchan$channel_disp)

# --- Color category: within-subject FDR (primary) ---
# Significance drives the dot/error-bar color; beta-trigger rows are flagged
# only via the pink y-axis tick label (set below).
cb_perchan$color_cat <- factor(
  ifelse(cb_perchan$sig_fdr, "FDR q<0.05 (within subj)",
    ifelse(cb_perchan$sig_uncorr, "p<0.05 uncorr", "ns")),
  levels = c("ns", "p<0.05 uncorr", "FDR q<0.05 (within subj)"))

# --- Y-axis label colors: pink for beta-trigger rows, black otherwise ---
# element_text() accepts a vector of colors applied in factor-level order
# (for our factor that's bottom-to-top). Build from the sorted levels so the
# mapping matches the plotted y-axis ticks exactly.
is_beta_by_level <- sapply(levels(cb_perchan$cell_label), function(lbl) {
  any(cb_perchan$is_beta_ref & cb_perchan$cell_label == lbl)
})
y_label_colors_cb <- ifelse(is_beta_by_level, "deeppink3", "black")

p_cb_forest <- ggplot(cb_perchan,
    aes(y = cell_label, x = obs_diff, color = color_cat)) +
  theme_light(base_size = 17) +
  facet_wrap(~ numStims_f, nrow = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = lo_ci, xmax = hi_ci),
                 height = 0.28, linewidth = 0.7) +
  geom_point(size = 2.8) +
  scale_color_manual(
    values = c("ns" = "grey60",
               "p<0.05 uncorr" = "orange",
               "FDR q<0.05 (within subj)" = "red"),
    drop = FALSE, name = "") +
  labs(x = expression(paste("Median Conditioned - Baseline (", mu,
                            "V) with 95% bootstrap CI")),
       y = "",
       title = "Per-cell EP modulation from baseline across subjects",
       subtitle = sprintf(
         "%d cells, %d subjects. BH FDR within (subject x dose). Pink y-axis labels = beta trigger channel.",
         nrow(cb_perchan), length(unique(cb_perchan$sid)))) +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 15),
        axis.text.y  = element_text(size = 13, color = y_label_colors_cb),
        axis.text.x  = element_text(size = 15),
        axis.title.x = element_text(size = 16),
        plot.title   = element_text(size = 19, face = "bold"),
        plot.subtitle = element_text(size = 14),
        strip.text   = element_text(size = 16, face = "bold"))

if (savePlot) {
  ggsave(here("output_plots", "betaStim_cond_vs_base_forest.png"),
         plot = p_cb_forest,
         units = "in", width = 16, height = 16, dpi = 600)
  ggsave(here("output_plots", "betaStim_cond_vs_base_forest.eps"),
         plot = p_cb_forest,
         units = "in", width = 16, height = 16, dpi = 600, device = cairo_ps)
}

# ========================================================================
# Null-burst probes vs Baseline permutation test (within-subject control)
# ------------------------------------------------------------------------
# 0b5a2e has a "null burst" condition (sham bursts; nullType=3 in the
# MATLAB extraction). Probes following null bursts should produce EPs
# equivalent to true baselines (>2s-after-burst probes) — this is the
# within-subject validity check for the null-burst control.
#
# The main script (betaStim_R_script.R) and this script both filter out
# numStims == "Null" from their primary analyses; here we run the test
# on the Null rows before that filter would apply. The script's top-of-
# file loading deliberately keeps Null trials in `data` (the filter line
# is commented out).
#
# Permutation: for each channel with enough trials (n >= 10 per group),
# pool Null + Base probes, shuffle the labels, recompute the difference
# of medians under permutation. Aggregated test uses sign-flip on
# per-channel median differences (one unit per channel), matching the
# clpb permutation approach above.
#
# Scope: 0b5a2e only for now. 0b5a2ePlayBack and ecb43e can be added by
# extending null_burst_subjects below.
# ========================================================================
cat("\n========== Null vs Baseline permutation (0b5a2e) ==========\n")

null_burst_subjects <- c("0b5a2e")
set.seed(42)
nPerm_nb <- 10000
min_n_per_group <- 10

nb_perchan_list <- list()
nb_per_subj_list <- list()
idx_nb <- 0L

for (sid_val in null_burst_subjects) {
  dN_all <- data[data$sid == sid_val &
                 (data$numStims == "Null" | data$numStims == "Base") &
                 !is.na(data$magnitude), ]
  if (nrow(dN_all) == 0) {
    cat(sprintf("%s: no Null/Base rows in CSV — skipping\n", sid_val))
    next
  }

  chan_diffs <- list()
  for (ch in sort(as.character(unique(dN_all$channel)))) {
    dC <- dN_all[as.character(dN_all$channel) == ch, ]
    mNull <- dC$magnitude[dC$numStims == "Null"]
    mBase <- dC$magnitude[dC$numStims == "Base"]
    if (length(mNull) < min_n_per_group || length(mBase) < min_n_per_group) next

    obs_diff <- median(mNull) - median(mBase)

    # Pooled two-sample permutation: shuffle the Null/Base labels among
    # all trials in this channel, recompute the difference of medians.
    all_mags <- c(mNull, mBase)
    is_null <- c(rep(TRUE, length(mNull)), rep(FALSE, length(mBase)))
    perm_diffs <- replicate(nPerm_nb, {
      shuf <- sample(is_null)
      median(all_mags[shuf]) - median(all_mags[!shuf])
    })
    perm_p <- mean(abs(perm_diffs) >= abs(obs_diff))

    idx_nb <- idx_nb + 1L
    nb_perchan_list[[idx_nb]] <- data.frame(
      sid = sid_val,
      channel = ch,
      n_null = length(mNull),
      n_base = length(mBase),
      median_null = round(median(mNull), 1),
      median_base = round(median(mBase), 1),
      obs_diff = round(obs_diff, 1),
      perm_p = round(perm_p, 4),
      stringsAsFactors = FALSE)

    chan_diffs[[ch]] <- obs_diff
  }

  # --- Aggregated test: exact sign-flip on per-channel median diffs ---
  # Each channel contributes one observation (median Null - median Base).
  # Under H0 (null bursts have no residual effect), the sign of each
  # channel's diff is arbitrary, so sign-flip permutation is the correct
  # null distribution. Test statistic = mean of channel diffs.
  # With n = 8 channels, 2^8 = 256 sign configurations — enumerated exactly.
  if (length(chan_diffs) >= 2) {
    ch_vec <- unlist(chan_diffs)
    ef_nb <- exact_signflip(ch_vec, stat_fn = mean)
    nb_per_subj_list[[sid_val]] <- data.frame(
      sid = sid_val,
      nChannels = length(ch_vec),
      mean_chan_diff = round(ef_nb$obs, 1),
      perm_p_aggregate = round(ef_nb$p_two_sided, 4),
      n_perms_exact = ef_nb$n_perms,
      stringsAsFactors = FALSE)
  }
}

nb_perchan <- if (length(nb_perchan_list) > 0)
  do.call(rbind, nb_perchan_list) else data.frame()
nb_per_subj <- if (length(nb_per_subj_list) > 0)
  do.call(rbind, nb_per_subj_list) else data.frame()

cat(sprintf("\nPer-channel Null vs Base (two-sided, n>=%d per group, %d perms):\n",
    min_n_per_group, nPerm_nb))
print(nb_perchan, row.names = FALSE)

cat("\nAggregated Null vs Base per subject (EXACT sign-flip enumeration on per-channel median diffs):\n")
print(nb_per_subj, row.names = FALSE)

if (nrow(nb_perchan) > 0) {
  cat(sprintf("\nCells with |Null - Base| perm p < 0.05: %d of %d\n",
      sum(nb_perchan$perm_p < 0.05), nrow(nb_perchan)))
}

# ========================================================================
# Prepare summary tables for export and write all as CSVs
# ------------------------------------------------------------------------
# CSVs are the machine-readable source of truth for the manuscript tables;
# the .docx report below reuses the same frames built here. Naming
# convention: betaStim_<scope>_<table>.csv in output_plots/.
# ========================================================================

# --- ANOVA tables (moved out of the docx block so CSVs can share them) ---
anova_clpb_tbl <- as.data.frame(anova(fit.clpb))
anova_clpb_tbl$Term <- rownames(anova_clpb_tbl)
anova_clpb_tbl <- anova_clpb_tbl[, c("Term", "Sum Sq", "Mean Sq", "NumDF",
                                      "DenDF", "F value", "Pr(>F)")]
anova_clpb_tbl[, 2:7] <- round(anova_clpb_tbl[, 2:7], 4)

anova_ec_tbl <- as.data.frame(anova(fit.ec))
anova_ec_tbl$Term <- rownames(anova_ec_tbl)
anova_ec_tbl <- anova_ec_tbl[, c("Term", "Sum Sq", "Mean Sq", "NumDF",
                                  "DenDF", "F value", "Pr(>F)")]
anova_ec_tbl[, 2:7] <- round(anova_ec_tbl[, 2:7], 4)

# --- perm_perchan: drop plot-formatting columns for the export ---
perm_perchan_export <- perm_perchan[, c(
  "channel_raw", "setToDeliverPhase", "phaseDeg_round", "phaseVecLength",
  "numStims", "n_probes", "obs_median_diff", "lo_ci", "hi_ci",
  "perm_p", "perm_q_bh"
)]
perm_perchan_export <- perm_perchan_export[order(
  perm_perchan_export$numStims,
  as.numeric(as.character(perm_perchan_export$channel_raw)),
  perm_perchan_export$phaseDeg_round), ]
perm_perchan_export$lo_ci <- round(perm_perchan_export$lo_ci, 1)
perm_perchan_export$hi_ci <- round(perm_perchan_export$hi_ci, 1)

# --- CSV writer helper ---
write_summary_csv <- function(df, name) {
  path <- here("output_plots", paste0(name, ".csv"))
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("Wrote %s (%d rows)\n", basename(path), nrow(df)))
}

cat("\n========== Writing summary / permutation CSVs ==========\n")

# Percent modulation tables (were previously docx-only)
write_summary_csv(pct_subj,  "betaStim_within_subject_pct_per_subject")
write_summary_csv(pct_grand, "betaStim_within_subject_pct_grand")
write_summary_csv(pct_chan,  "betaStim_within_subject_pct_per_channel")

# ANOVA tables (were previously docx-only)
write_summary_csv(anova_clpb_tbl, "betaStim_clpb_anova")
write_summary_csv(anova_ec_tbl,   "betaStim_ecb43e_anova")

# clpb sign-flip permutation tables (were previously console-only)
write_summary_csv(perm_chan,           "betaStim_clpb_perm_chan_aggregate")
write_summary_csv(perm_perchan_export, "betaStim_clpb_perm_perchan_bycell")
if (exists("clpb_region_direction") && nrow(clpb_region_direction) > 0) {
  write_summary_csv(clpb_region_direction, "betaStim_clpb_region_direction")
}
if (exists("cb_halves") && nrow(cb_halves) > 0) {
  write_summary_csv(cb_halves, "betaStim_cond_vs_base_phase_halves")
}

# Null vs baseline permutation tables (were previously console-only)
if (nrow(nb_perchan) > 0)  write_summary_csv(nb_perchan,  "betaStim_null_vs_base_perchan")
if (nrow(nb_per_subj) > 0) write_summary_csv(nb_per_subj, "betaStim_null_vs_base_aggregate")

# Conditioned vs Baseline per-cell permutation + FDR tables
if (exists("cb_perchan") && nrow(cb_perchan) > 0) {
  write_summary_csv(cb_perchan,        "betaStim_cond_vs_base_perchan")
  write_summary_csv(cb_subj_fdr,       "betaStim_cond_vs_base_per_subject")
  write_summary_csv(cb_subj_presence,  "betaStim_cond_vs_base_subject_presence")
  write_summary_csv(cb_summary,        "betaStim_cond_vs_base_summary")
  write_summary_csv(cb_summary_pooled, "betaStim_cond_vs_base_pooled_summary")

  # ----------------------------------------------------------------------
  # Trigger-channel modulation ranking per subject
  # ----------------------------------------------------------------------
  # Quantifies how often the beta-trigger channel (betaLabels == 1, the
  # channel the closed-loop phase-triggering was locked to) shows the
  # largest dose-dependent CEP modulation within its subject. Aggregates
  # cb_perchan's obs_diff (median conditioned - median baseline, uV) across
  # phase cells to one value per (subject x channel x dose) via max, then
  # ranks channels within (subject x dose). Rank 1 = biggest modulation.
  # Subjects with n_channels == 1 (only 702d24 after filtering) are trivially
  # rank 1 and flagged for separate interpretation.
  cb_channel_max <- plyr::ddply(cb_perchan,
    .(sid, channel_raw, numStims, is_beta_ref), summarize,
    obs_diff_max = max(obs_diff))
  cb_channel_max$rank_within <- ave(cb_channel_max$obs_diff_max,
    cb_channel_max$sid, cb_channel_max$numStims,
    FUN = function(x) rank(-x, ties.method = "min"))
  cb_channel_max$n_channels_subj <- ave(cb_channel_max$channel_raw,
    cb_channel_max$sid, cb_channel_max$numStims, FUN = length)

  trig_rank <- cb_channel_max[cb_channel_max$is_beta_ref == TRUE,
    c("sid", "numStims", "rank_within", "n_channels_subj", "obs_diff_max")]
  trig_rank <- trig_rank[order(trig_rank$sid, trig_rank$numStims), ]
  names(trig_rank) <- c("sid", "numStims", "trigger_rank",
                        "n_channels", "trigger_obs_diff_uV")
  trig_rank$trigger_is_top <- trig_rank$trigger_rank == 1
  trig_rank$trivial_single_channel <- trig_rank$n_channels == 1
  cat("\n--- Trigger-channel modulation rank within each subject x dose ---\n")
  print(trig_rank)

  trig_summary <- plyr::ddply(trig_rank, .(numStims), summarize,
    n_subjects_total             = length(sid),
    n_subjects_multichannel      = sum(!trivial_single_channel),
    n_trigger_top_all            = sum(trigger_is_top),
    n_trigger_top_multichannel   = sum(trigger_is_top & !trivial_single_channel),
    median_trigger_rank          = median(trigger_rank),
    median_trigger_rank_mc       = median(trigger_rank[!trivial_single_channel]))
  cat("\n--- Trigger-channel top-rank summary per dose ---\n")
  print(trig_summary)

  write_summary_csv(trig_rank,    "betaStim_trigger_channel_rank")
  write_summary_csv(trig_summary, "betaStim_trigger_channel_rank_summary")
}

# ========================================================================
# Export all tables to .docx
# ========================================================================
if (require(officer) && require(flextable)) {
  doc <- read_docx()
  doc <- body_add_par(doc, "Within-Subject Phase Analysis", style = "heading 1")

  # --- CL vs PB ANOVA ---
  doc <- body_add_par(doc, "0b5a2e: CL vs PB Combined Model (ANOVA)", style = "heading 2")
  ft <- flextable(anova_clpb_tbl) |> autofit() |>
    set_caption("Type III ANOVA: CL+PB combined model. 96 obs, 8 channels. condition ref = PB.")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # --- ecb43e ANOVA ---
  doc <- body_add_par(doc, "ecb43e: Targeted vs Random Combined Model (ANOVA)", style = "heading 2")
  ft <- flextable(anova_ec_tbl) |> autofit() |>
    set_caption("Type III ANOVA: ecb43e combined model. 36 obs, 4 channels. condType ref = random.")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # --- clpb channel-aggregated sign-flip (primary inferential claim) ---
  doc <- body_add_par(doc,
    "0b5a2e: CL vs PB Channel-Aggregated Sign-Flip (Exact, Holm-corrected)",
    style = "heading 2")
  ft <- flextable(perm_chan) |> autofit() |>
    set_caption("Exact sign-flip permutation on per-channel-condition mean differences (2^n sign configurations enumerated); test statistic = mean of 13-16 channel x phase-condition mean diffs per dose. Holm step-down across 3 doses. Primary inferential test for the CL vs PB comparison at the aggregated level.")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # --- clpb per-cell sign-flip ---
  doc <- body_add_par(doc,
    "0b5a2e: CL vs PB Per Channel x Phase x Dose (Sign-Flip + Bootstrap CIs)",
    style = "heading 2")
  ft <- flextable(perm_perchan_export) |> autofit() |>
    set_caption("Per-cell sign-flip permutation on trial-level paired differences (10k Monte Carlo draws). Bootstrap 95% CIs from 2000 resamples on paired$diff. BH FDR within dose (perm_q_bh) — each dose level forms its own correction family of 13-16 cells. Supplementary detail to the channel-aggregated test above; rows sorted by dose then channel then delivered phase.")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # --- Null vs Baseline per channel ---
  if (nrow(nb_perchan) > 0) {
    doc <- body_add_par(doc,
      "0b5a2e: Null-Burst vs Baseline Per Channel (Two-Sample Permutation)",
      style = "heading 2")
    ft <- flextable(nb_perchan) |> autofit() |>
      set_caption("Per-channel pooled two-sample permutation: shuffle Null/Base labels within channel without replacement (10k MC iterations), recompute median(Null) - median(Base). Within-subject validity check for the null-burst control condition.")
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")
  }

  # --- Null vs Baseline aggregated ---
  if (nrow(nb_per_subj) > 0) {
    doc <- body_add_par(doc,
      "0b5a2e: Null-Burst vs Baseline Aggregated (Exact Sign-Flip)",
      style = "heading 2")
    ft <- flextable(nb_per_subj) |> autofit() |>
      set_caption("Aggregated exact sign-flip on per-channel median differences (n=8 channels, 2^8 = 256 sign configurations enumerated). Non-significant mean chan diff indicates null-burst probes are statistically indistinguishable from baseline probes, validating the null-burst control.")
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")
  }

  # --- Percent modulation: per subject ---
  doc <- body_add_par(doc, "Percent Modulation from Baseline (per subject)", style = "heading 2")
  ft <- flextable(pct_subj) |> autofit() |>
    set_caption("Mean percent change from baseline per subject. Baseline = median of Base trials per channel, pooled across conditions.")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # --- Percent modulation: grand summary ---
  doc <- body_add_par(doc, "Percent Modulation from Baseline (grand summary)", style = "heading 2")
  ft <- flextable(pct_grand) |> autofit() |>
    set_caption("Grand mean percent change across all subjects and channels.")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # --- Percent modulation: per channel detail ---
  doc <- body_add_par(doc, "Percent Modulation from Baseline (per channel)", style = "heading 2")
  ft <- flextable(pct_chan) |> autofit() |>
    set_caption("Percent change per channel per dose. baseMedian = channel baseline (uV), mag = conditioned median (uV).")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # --- Conditioned vs Baseline: per-cell permutation + FDR ---
  if (exists("cb_perchan") && nrow(cb_perchan) > 0) {
    doc <- body_add_par(doc,
      "Conditioned vs Baseline: Per-Cell Permutation + Within-Subject FDR",
      style = "heading 1")
    doc <- body_add_par(doc,
      paste("Per-cell two-sample label-shuffle permutation (10,000 MC",
            "iterations, median difference statistic) of conditioned probe",
            "magnitudes vs channel baseline magnitudes, at each dose separately.",
            "Cell = (subject x channel x phaseDeg_round).",
            "Inclusion: channel-level phaseVecLength >= 0.2 (matches 5a-gf2),",
            ">=10 baseline probes, >=5 conditioned probes.",
            "BH FDR correction applied within (subject x dose) family because",
            "channels are nested within subjects."),
      style = "Normal")

    # Per-dose cell-count summary (collapses subjects; uses within-subject FDR)
    doc <- body_add_par(doc, "Per-dose cell-count summary",
                        style = "heading 2")
    ft <- flextable(cb_summary) |> autofit() |>
      set_caption(paste("Cells across subjects by dose. n_sig_uncorr = cells with",
                        "raw permutation p < 0.05. n_sig_fdr = cells surviving",
                        "Benjamini-Hochberg FDR (q < 0.05) applied within each",
                        "(subject x dose) family. Percentages are over n_cells."))
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")

    # Across-subject summary (subject-level presence)
    doc <- body_add_par(doc, "Subject-level presence summary per dose",
                        style = "heading 2")
    ft <- flextable(cb_subj_presence) |> autofit() |>
      set_caption(paste("Subject-level summary. n_subj_any_fdr = count of",
                        "subjects with >=1 FDR-significant cell at each dose.",
                        "median_pct_sig_fdr = median across subjects of the",
                        "per-subject fraction modulated."))
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")

    # Per-subject breakdown
    doc <- body_add_par(doc, "Per-subject cell counts", style = "heading 2")
    ft <- flextable(cb_subj_fdr) |> autofit() |>
      set_caption(paste("Per-subject x dose: cells tested, cells significant",
                        "at uncorrected p<0.05 and at FDR q<0.05",
                        "(correction family = subject's cells at that dose)."))
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")

    # Pooled FDR sensitivity
    doc <- body_add_par(doc,
      "Pooled-FDR sensitivity (secondary)", style = "heading 2")
    ft <- flextable(cb_summary_pooled) |> autofit() |>
      set_caption(paste("Pooled BH FDR across all cells within each dose",
                        "(ignores subject clustering). Reported as a",
                        "sensitivity check alongside the primary",
                        "within-subject FDR results above."))
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")

    # Trigger-channel modulation rank
    if (exists("trig_rank") && nrow(trig_rank) > 0) {
      doc <- body_add_par(doc, "Trigger-channel modulation rank per subject",
                          style = "heading 2")
      doc <- body_add_par(doc,
        paste("For each subject x dose, channels were ranked by their maximum",
              "conditioned-minus-baseline median difference (uV) across phase",
              "cells. Rank 1 indicates the channel with the largest dose-dependent",
              "CEP modulation. The beta-trigger channel is the channel used as",
              "the phase reference for closed-loop stimulation (betaLabels == 1).",
              "Subjects with only one included channel are trivially rank 1 and",
              "flagged."),
        style = "Normal")

      tr_out <- trig_rank
      names(tr_out) <- c("Subject", "Dose", "Trigger rank",
                         "N channels", "Trigger obs_diff (uV)",
                         "Trigger at top?", "Trivial (1 channel)")
      ft <- flextable(tr_out) |> autofit() |>
        set_caption("Per (subject x dose) rank of the beta-trigger channel within the subject's included channels, based on the maximum observed conditioned-minus-baseline median (uV) across phase cells.")
      doc <- body_add_flextable(doc, ft)
      doc <- body_add_par(doc, "")

      ts_out <- trig_summary
      names(ts_out) <- c("Dose", "N subjects (total)",
                         "N subjects (multi-channel)",
                         "Trigger top-ranked (all)",
                         "Trigger top-ranked (multi-channel)",
                         "Median trigger rank (all)",
                         "Median trigger rank (multi-channel)")
      doc <- body_add_par(doc, "Trigger-channel top-rank summary per dose",
                          style = "heading 2")
      ft <- flextable(ts_out) |> autofit() |>
        set_caption("Count of subjects where the beta-trigger channel was the top-modulated channel (rank 1) at each dose level. 'Multi-channel' excludes subjects with only 1 included channel (trivially rank 1). Median rank summarizes where the trigger channel falls when not at the top.")
      doc <- body_add_flextable(doc, ft)
      doc <- body_add_par(doc, "")
    }
  }

  docx_path <- here("output_plots", "betaStim_within_subject_tables.docx")
  print(doc, target = docx_path)
  cat(sprintf("\nSaved within-subject tables to: %s\n", docx_path))
}
