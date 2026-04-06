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
tab_model(fit.lm)

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
