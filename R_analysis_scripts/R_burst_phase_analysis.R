# R_burst_phase_analysis.R
# Test whether delivered beta phase predicts CEP magnitude.
# Merges per-burst phase precision (from compute_burst_phase_precision.m)
# with EP magnitudes using probeSample as the join key.
#
# Analyses:
# 1. 270° bin: bursts delivered within ±45° of 270° produce larger CEPs?
# 2. Four bins centered on 0, 90, 180, 270 (±45° each)
# 3. Phase error: angular distance from target → CEP magnitude

setwd('/Users/davidcaldwell/code/betaOscillationTriggerStimPaper')
library('plyr')
library('ggplot2')
library('lme4')
library('lmerTest')
library('emmeans')
library('here')

sids_to_analyze <- c('d5cd55', 'c91479', '7dbdec', '9ab7ab', '702d24', 'ecb43e', '0b5a2e')

# load EP data
ep <- read.table("data/output_table/betaStim_outputTable_50_new_100_thresh.csv",
  header=TRUE, sep=",", stringsAsFactors=F,
  colClasses=c("magnitude"="numeric","betaLabels"="factor","sid"="factor",
    "numStims"="factor","stimLevel"="numeric","channel"="factor",
    "subjectNum"="factor","phaseClass"="factor","setToDeliverPhase"="factor",
    "probeSample"="numeric"))
ep <- subset(ep, magnitude < 1500 & magnitude > 25 & !is.nan(magnitude))
ep$numStims <- revalue(ep$numStims, c("Test 1"="[1,2]","Test 2"="[3,4]","Test 3"="[5,inf)"))

# merge across all subjects
all_merged <- data.frame()

for (sid in sids_to_analyze) {
  bp_file <- paste0("data/output_table/", sid, "_burst_phase_precision.csv")
  if (!file.exists(bp_file)) {
    cat(sprintf("Skipping %s: no burst precision file\n", sid))
    next
  }
  bp <- read.csv(bp_file)
  bp_kept <- bp[bp$probeClass == 0 | bp$probeClass == 1, ]
  ep_subj <- ep[ep$sid == sid, ]
  merged <- merge(ep_subj, bp_kept,
    by.x = c("channel", "probeSample"),
    by.y = c("channelEncoded", "probeSample"),
    suffixes = c(".ep", ".bp"))
  cat(sprintf("%s: EP=%d, BP=%d, merged=%d\n", sid, nrow(ep_subj), nrow(bp_kept), nrow(merged)))
  all_merged <- rbind(all_merged, merged)
}
cat(sprintf("\nTotal merged: %d rows, %d subjects\n", nrow(all_merged), length(unique(all_merged$sid))))

# ============================================================
# Filter: conditioned probes with >= 1 good beta fit
# ============================================================
cond_beta <- all_merged[all_merged$probeClass == 1 & all_merged$nGoodBeta > 0, ]

# phase bins centered on 0, 90, 180, 270 (±45° each)
# near0: 315-360 or 0-45
# near90: 45-135
# near180: 135-225
# near270: 225-315
cond_beta$phaseBin <- ifelse(
  cond_beta$burstCircMean >= 315 | cond_beta$burstCircMean < 45, "0",
  ifelse(cond_beta$burstCircMean >= 45 & cond_beta$burstCircMean < 135, "90",
  ifelse(cond_beta$burstCircMean >= 135 & cond_beta$burstCircMean < 225, "180", "270")))
cond_beta$phaseBin <- factor(cond_beta$phaseBin, levels = c("0", "90", "180", "270"))

# binary: 270 bin vs others
cond_beta$is270 <- cond_beta$phaseBin == "270"

# phase error: angular distance from target
cond_beta$phaseError <- pmin(
  abs(cond_beta$burstCircMean - cond_beta$targetPhase),
  360 - abs(cond_beta$burstCircMean - cond_beta$targetPhase))

cat(sprintf("Conditioned with beta fits: %d rows\n", nrow(cond_beta)))
cat("Phase bins:\n")
print(table(cond_beta$phaseBin))
cat(sprintf("\n270 bin: %d, others: %d\n\n", sum(cond_beta$is270), sum(!cond_beta$is270)))

set.seed(42)
nPerm <- 10000

# ============================================================
# ANALYSIS 1: 270° bin vs others
# ============================================================
cat("============================================================\n")
cat("ANALYSIS 1: 270 deg bin (225-315) vs others\n")
cat("============================================================\n\n")

cat("--- Per subject at [5,inf) ---\n\n")
for (sid_val in unique(cond_beta$sid)) {
  high <- cond_beta[cond_beta$sid == sid_val & cond_beta$numStims == "[5,inf)", ]
  n270 <- sum(high$is270)
  nOther <- nrow(high) - n270
  if (n270 >= 3 & nOther >= 3) {
    med270 <- median(high$magnitude[high$is270])
    medOther <- median(high$magnitude[!high$is270])
    obs <- med270 - medOther
    perms <- numeric(nPerm)
    for (p in 1:nPerm) {
      shuf <- sample(high$is270)
      perms[p] <- median(high$magnitude[shuf]) - median(high$magnitude[!shuf])
    }
    perm_p <- mean(abs(perms) >= abs(obs))
    cat(sprintf("  %s: 270 n=%d med=%.1f, other n=%d med=%.1f, diff=%+.1f, p=%.4f\n",
      sid_val, n270, med270, nOther, medOther, obs, perm_p))
  } else {
    cat(sprintf("  %s: insufficient (270=%d, other=%d)\n", sid_val, n270, nOther))
  }
}

cat("\n--- All subjects pooled by dose ---\n\n")
for (dose in c("[1,2]", "[3,4]", "[5,inf)")) {
  rows <- cond_beta[cond_beta$numStims == dose, ]
  n270 <- sum(rows$is270)
  nOther <- nrow(rows) - n270
  med270 <- median(rows$magnitude[rows$is270])
  medOther <- median(rows$magnitude[!rows$is270])
  obs <- med270 - medOther
  perms <- numeric(nPerm)
  for (p in 1:nPerm) {
    shuf <- sample(rows$is270)
    perms[p] <- median(rows$magnitude[shuf]) - median(rows$magnitude[!shuf])
  }
  perm_p <- mean(abs(perms) >= abs(obs))
  cat(sprintf("  %s: 270 n=%d med=%.1f, other n=%d med=%.1f, diff=%+.1f, p=%.4f\n",
    dose, n270, med270, nOther, medOther, obs, perm_p))
}

# ============================================================
# ANALYSIS 2: Four phase bins at [5,inf)
# ============================================================
cat("\n============================================================\n")
cat("ANALYSIS 2: Four phase bins (0, 90, 180, 270) at [5,inf)\n")
cat("============================================================\n\n")

high_all <- cond_beta[cond_beta$numStims == "[5,inf)", ]

cat("--- All subjects pooled ---\n")
for (q in levels(high_all$phaseBin)) {
  rows <- high_all[high_all$phaseBin == q, ]
  if (nrow(rows) > 0) {
    cat(sprintf("  %s deg: n=%d, median mag=%.1f\n", q, nrow(rows), median(rows$magnitude)))
  }
}

cat("\n--- Per subject ---\n\n")
for (sid_val in unique(high_all$sid)) {
  subj <- high_all[high_all$sid == sid_val, ]
  cat(sprintf("  %s:\n", sid_val))
  for (q in levels(subj$phaseBin)) {
    rows <- subj[subj$phaseBin == q, ]
    if (nrow(rows) > 0) {
      cat(sprintf("    %s deg: n=%d, median=%.1f\n", q, nrow(rows), median(rows$magnitude)))
    }
  }
}

# ============================================================
# ANALYSIS 3: Phase error
# ============================================================
cat("\n============================================================\n")
cat("ANALYSIS 3: Phase error (angular distance from target)\n")
cat("    negative rho = tighter phase -> larger CEP\n")
cat("============================================================\n\n")

cat("--- Per subject at [5,inf) ---\n\n")
for (sid_val in unique(cond_beta$sid)) {
  high <- cond_beta[cond_beta$sid == sid_val & cond_beta$numStims == "[5,inf)", ]
  if (nrow(high) >= 10 & sd(high$phaseError) > 0) {
    ct <- suppressWarnings(cor.test(high$phaseError, high$magnitude, method="spearman"))
    cat(sprintf("  %s: n=%d, rho=%+.3f, p=%.4f\n", sid_val, nrow(high), ct$estimate, ct$p.value))
  } else {
    cat(sprintf("  %s: n=%d, insufficient\n", sid_val, nrow(high)))
  }
}

cat("\n--- All subjects pooled by dose ---\n\n")
for (dose in c("[1,2]", "[3,4]", "[5,inf)")) {
  rows <- cond_beta[cond_beta$numStims == dose, ]
  ct <- suppressWarnings(cor.test(rows$phaseError, rows$magnitude, method="spearman"))
  cat(sprintf("  %s: n=%d, rho=%+.3f, p=%.4f\n", dose, nrow(rows), ct$estimate, ct$p.value))
}

# ============================================================
# Summary-level mixed model: 270 bin effect
# ============================================================
cat("\n============================================================\n")
cat("Summary-level mixed model: 270 bin effect\n")
cat("============================================================\n\n")

cond_beta_nb <- cond_beta[cond_beta$numStims != "Base", ]
cond_beta_nb$phaseDir <- ifelse(cond_beta_nb$is270, "270", "other")

summary_phase <- ddply(cond_beta_nb, .(sid, channel, numStims, phaseDir),
  summarize, magnitude = median(magnitude))
basePerChan <- ddply(all_merged[all_merged$numStims == "Base", ],
  .(sid, channel), summarize, baselineMag = median(magnitude))
summary_phase <- merge(summary_phase, basePerChan, by = c("sid", "channel"), all.x = TRUE)
summary_phase <- summary_phase[!is.na(summary_phase$baselineMag), ]
summary_phase$baselineMag_c <- summary_phase$baselineMag - mean(summary_phase$baselineMag)
summary_phase$channel <- factor(summary_phase$channel)
summary_phase$sid <- factor(summary_phase$sid)

cat(sprintf("Summary data: %d cells, %d subjects, %d channels\n",
  nrow(summary_phase), length(unique(summary_phase$sid)), length(unique(summary_phase$channel))))
print(table(summary_phase$phaseDir, summary_phase$numStims))

fit <- lmerTest::lmer(magnitude ~ numStims * phaseDir + baselineMag_c +
  (1|sid) + (1|channel),
  data = summary_phase,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=20000)))

cat(sprintf("\nSingular: %s\n", isSingular(fit)))
summary(fit)
cat("\nType III ANOVA:\n")
print(anova(fit))
cat("\n270 bin vs other at each dose:\n")
emm <- emmeans(fit, ~ phaseDir | numStims)
print(confint(pairs(emm)))

# ============================================================
# Plots
# ============================================================
emm_df <- as.data.frame(emm)
emm_df$phaseDir <- factor(emm_df$phaseDir, levels = c("other", "270"))

pd <- position_dodge(0.15)
p1 <- ggplot(emm_df, aes(x = numStims, y = emmean, color = phaseDir, group = phaseDir)) +
  theme_light(base_size = 14) +
  geom_line(position = pd, linewidth = 0.8) +
  geom_point(position = pd, size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = pd, linewidth = 0.8) +
  labs(x = "Number of Conditioning Stimuli",
       y = expression(paste("Predicted Magnitude (", mu, "V)")),
       color = "Delivered Phase",
       title = "CEP Magnitude: 270 deg bin vs Others (all subjects)") +
  scale_color_hue(labels = c("Other phases", "270 deg (225-315)"))
ggsave(here("output_plots","betaStim_burst_270bin_vs_other.png"), plot = p1,
       units = "in", width = 6.5, height = 4.5, dpi = 600)

high_all <- cond_beta[cond_beta$numStims == "[5,inf)", ]
p2 <- ggplot(high_all, aes(x = phaseError, y = magnitude)) +
  theme_light(base_size = 14) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ sid, scales = "free_y") +
  labs(x = "Phase Error (degrees from target)",
       y = expression(paste("Magnitude (", mu, "V)")),
       title = "Phase Error vs CEP Magnitude at [5,inf) by Subject")
ggsave(here("output_plots","betaStim_phase_error_scatter.png"), plot = p2,
       units = "in", width = 10, height = 8, dpi = 600)

cat("\nPlots saved.\n")
