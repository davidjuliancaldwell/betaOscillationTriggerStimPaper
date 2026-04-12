library("Hmisc"); library("ggplot2"); library("lme4"); library("plyr"); library("here"); library("lmerTest"); library("emmeans"); library("dplyr")

data <- read.table(here("data","output_table","betaStim_outputTable_50_new_100_thresh.csv"),header=TRUE,sep=",",stringsAsFactors=F,
  colClasses=c("magnitude"="numeric","betaLabels"="factor","sid"="factor","numStims"="factor","stimLevel"="numeric","channel"="factor","subjectNum"="factor","phaseClass"="factor","phaseDeg"="numeric","setToDeliverPhase"="factor"))
data <- subset(data, magnitude<1500 & magnitude>25)
data <- subset(data,!is.nan(data$magnitude))
data <- subset(data, sid!="702d24" & sid!="0b5a2ePlayBack" & numStims!="Null")
data$numStims <- revalue(data$numStims, c("Test 1"="[1,2]","Test 2"="[3,4]","Test 3"="[5,inf)"))

# circular phase decomposition
data$phase_rad <- data$phaseDeg * pi / 180
data$sin_phase <- sin(data$phase_rad)
data$cos_phase <- cos(data$phase_rad)

# compute absDiff per trial (one baseMedian per channel, pooled across phaseClass)
data$absDiff <- 0
for (name in unique(data$sid)){
  for (chan in unique(data[data$sid == name,]$channel)){
    base <- data[data$sid == name & data$channel == chan & data$numStims == "Base",]$magnitude
    baseMedian <- median(base)
    data[data$sid == name & data$channel == chan,]$absDiff <- data[data$sid == name & data$channel == chan,]$magnitude - baseMedian
  }
}
dataNoBaseline <- data[data$numStims != "Base",]

# =============================================
# 1. absDiff model (no baseline category)
# =============================================
summaryNB_abs <- ddply(dataNoBaseline, .(sid,phaseClass,numStims,channel), summarize, absDiff = median(absDiff))

fit.abs <- lmerTest::lmer(absDiff ~ numStims * phaseClass + (1|sid) + (1|channel),
  data = summaryNB_abs, control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=20000)))

cat("=== 1. absDiff model (no baseline category) ===\n")
cat("Singular:", isSingular(fit.abs), "\n")
print(anova(fit.abs))

emm_abs <- emmeans(fit.abs, ~ numStims | phaseClass)
cat("\nEmmeans (already change from baseline):\n")
print(emm_abs)
cat("\nPairwise:\n")
print(confint(pairs(emm_abs)))

emm_abs_phase <- emmeans(fit.abs, ~ phaseClass | numStims)
cat("\nPhase contrasts:\n")
print(confint(pairs(emm_abs_phase)))

# =============================================
# 2. Model D: magnitude with baseline + random dose slope
# =============================================
summaryAll <- ddply(data, .(sid,phaseClass,numStims,channel), summarize, magnitude = median(magnitude))
summaryAll$numStims <- relevel(summaryAll$numStims, ref = "Base")
summaryAll$doseNum <- as.numeric(factor(summaryAll$numStims, levels=c("Base","[1,2]","[3,4]","[5,inf)"))) - 1

fit.d <- lmerTest::lmer(magnitude ~ numStims * phaseClass + (1+doseNum|sid) + (1|channel),
  data = summaryAll, control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=50000)))

cat("\n\n=== 2. Model D (magnitude + baseline category + random slope) ===\n")
cat("Singular:", isSingular(fit.d), "\n")
print(anova(fit.d))

emm_d <- emmeans(fit.d, ~ numStims | phaseClass)
contr_d <- contrast(emm_d, method = "trt.vs.ctrl", ref = "Base")
cat("\nContrasts vs baseline:\n")
print(confint(contr_d))

emm_d_phase <- emmeans(fit.d, ~ phaseClass | numStims)
cat("\nPhase contrasts:\n")
print(confint(pairs(emm_d_phase)))

# =============================================
# 3. ANCOVA: baseline magnitude as covariate
# =============================================
summaryNB_mag <- ddply(dataNoBaseline, .(sid,phaseClass,numStims,channel), summarize, magnitude = median(magnitude))

# one baseline median per channel (pooled across phaseClass)
basePerChan <- ddply(data[data$numStims == "Base",], .(sid, channel), summarize, baselineMag = median(magnitude))
summaryNB_mag <- merge(summaryNB_mag, basePerChan, by = c("sid", "channel"))

# center baseline
summaryNB_mag$baselineMag_c <- summaryNB_mag$baselineMag - mean(summaryNB_mag$baselineMag)
summaryNB_mag$doseNum <- as.numeric(factor(summaryNB_mag$numStims, levels=c("[1,2]","[3,4]","[5,inf)"))) - 1

cat(sprintf("\n\n=== 3. ANCOVA (baseline as covariate, centered at %.1f uV) ===\n", mean(summaryNB_mag$baselineMag)))

fit.ancova <- lmerTest::lmer(magnitude ~ numStims * phaseClass + baselineMag_c +
  (1+doseNum|sid) + (1|channel),
  data = summaryNB_mag, control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=50000)))

cat("Singular:", isSingular(fit.ancova), "\n")
summary(fit.ancova)
print(anova(fit.ancova))

emm_anc <- emmeans(fit.ancova, ~ numStims | phaseClass)
cat("\nEmmeans (at mean baseline):\n")
print(emm_anc)
cat("\nPairwise dose contrasts:\n")
print(confint(pairs(emm_anc)))

emm_anc_phase <- emmeans(fit.ancova, ~ phaseClass | numStims)
cat("\nPhase contrasts:\n")
print(confint(pairs(emm_anc_phase)))

# =============================================
# PLOTS
# =============================================
pd <- position_dodge(0.15)

# --- Plot 1: absDiff emmeans ---
emm_abs_df <- as.data.frame(emm_abs)
emm_abs_df$phaseClass <- factor(emm_abs_df$phaseClass, levels = c("90","270"))

p1 <- ggplot(emm_abs_df, aes(x = numStims, y = emmean, color = phaseClass, group = phaseClass)) +
  theme_light(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(position = pd, linewidth = 0.8) + geom_point(position = pd, size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = pd, linewidth = 0.8) +
  labs(x = "Number of Conditioning Stimuli",
       y = expression(paste("Absolute Difference from Baseline (", mu, "V)")),
       color = "Phase Class",
       title = "absDiff Model: Emmeans (baseline pre-subtracted)") +
  scale_color_hue(labels = c("Depolarizing (90)", "Hyperpolarizing (270)"))
ggsave(here("output_plots","betaStim_emmip_absDiff.png"), plot = p1, units = "in", width = 6.5, height = 4.5, dpi = 600)

# --- Plot 2: Model D contrasts vs baseline ---
contr_d_ci <- as.data.frame(confint(contr_d))
contr_d_ci$dose <- sub(" - .*", "", contr_d_ci$contrast)
contr_d_ci$phaseClass <- as.character(contr_d_ci$phaseClass)
ref_d <- data.frame(phaseClass = as.character(unique(contr_d_ci$phaseClass)),
  dose = "Base", estimate = 0, SE = 0, lower.CL = 0, upper.CL = 0)
plot_d <- bind_rows(
  ref_d %>% select(phaseClass, dose, estimate, SE, lower.CL, upper.CL),
  contr_d_ci %>% select(phaseClass, dose, estimate, SE, lower.CL, upper.CL))
plot_d$dose <- factor(plot_d$dose, levels = c("Base","[1,2]","[3,4]","[5,inf)"))
plot_d$phaseClass <- factor(plot_d$phaseClass, levels = c("90","270"))

p2 <- ggplot(plot_d, aes(x = dose, y = estimate, color = phaseClass, group = phaseClass)) +
  theme_light(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(position = pd, linewidth = 0.8) + geom_point(position = pd, size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = pd, linewidth = 0.8) +
  labs(x = "Number of Conditioning Stimuli",
       y = expression(paste(Delta, " Magnitude from Baseline (", mu, "V)")),
       color = "Phase Class",
       title = "Model D: Contrasts vs Baseline (random dose slope)") +
  scale_color_hue(labels = c("Depolarizing (90)", "Hyperpolarizing (270)"))
ggsave(here("output_plots","betaStim_emmip_modelD.png"), plot = p2, units = "in", width = 6.5, height = 4.5, dpi = 600)

# --- Plot 3: ANCOVA contrasts vs [1,2] ---
contr_anc_ref <- contrast(emm_anc, method = "trt.vs.ctrl", ref = 1)
contr_anc_ci <- as.data.frame(confint(contr_anc_ref))
contr_anc_ci$dose <- sub(" - .*", "", contr_anc_ci$contrast)
contr_anc_ci$phaseClass <- as.character(contr_anc_ci$phaseClass)
ref_anc <- data.frame(phaseClass = as.character(unique(contr_anc_ci$phaseClass)),
  dose = "[1,2]", estimate = 0, SE = 0, lower.CL = 0, upper.CL = 0)
plot_anc <- bind_rows(
  ref_anc %>% select(phaseClass, dose, estimate, SE, lower.CL, upper.CL),
  contr_anc_ci %>% select(phaseClass, dose, estimate, SE, lower.CL, upper.CL))
plot_anc$dose <- factor(plot_anc$dose, levels = c("[1,2]","[3,4]","[5,inf)"))
plot_anc$phaseClass <- factor(plot_anc$phaseClass, levels = c("90","270"))

p3 <- ggplot(plot_anc, aes(x = dose, y = estimate, color = phaseClass, group = phaseClass)) +
  theme_light(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(position = pd, linewidth = 0.8) + geom_point(position = pd, size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = pd, linewidth = 0.8) +
  labs(x = "Number of Conditioning Stimuli",
       y = expression(paste(Delta, " Magnitude from [1,2] (", mu, "V)")),
       color = "Phase Class",
       title = "ANCOVA: Dose-Response (baseline covariate, random dose slope)") +
  scale_color_hue(labels = c("Depolarizing (90)", "Hyperpolarizing (270)"))
ggsave(here("output_plots","betaStim_emmip_ancova.png"), plot = p3, units = "in", width = 6.5, height = 4.5, dpi = 600)

cat("\nAll 3 plots saved.\n")

# =============================================
# 4. Sin/cos ANCOVA: continuous phase, ordinal dose
# =============================================
summaryNB_sincos <- ddply(dataNoBaseline, .(sid,phaseClass,numStims,channel), summarize,
  magnitude = median(magnitude), sin_phase = first(sin_phase), cos_phase = first(cos_phase),
  betaLabels = first(betaLabels))
basePerChan_sc <- ddply(data[data$numStims == "Base",], .(sid, channel), summarize, baselineMag = median(magnitude))
summaryNB_sincos <- merge(summaryNB_sincos, basePerChan_sc, by = c("sid", "channel"))
summaryNB_sincos$baselineMag_c <- summaryNB_sincos$baselineMag - mean(summaryNB_sincos$baselineMag)
summaryNB_sincos$doseNum <- as.numeric(factor(summaryNB_sincos$numStims,
  levels = c("[1,2]","[3,4]","[5,inf)"))) - 1
summaryNB_sincos$numStims_ord <- ordered(summaryNB_sincos$numStims,
  levels = c("[1,2]", "[3,4]", "[5,inf)"))
poly_lin <- contr.poly(3)[, 1]
summaryNB_sincos$dose_linpoly <- poly_lin[as.numeric(summaryNB_sincos$numStims_ord)]

fit.sincos <- lmerTest::lmer(
  magnitude ~ numStims_ord * (sin_phase + cos_phase) + betaLabels + baselineMag_c +
  (1 + dose_linpoly | sid) + (1 | channel),
  data = summaryNB_sincos,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

cat("\n\n=== 4. Sin/cos ANCOVA (ordinal dose, continuous phase) ===\n")
cat("Singular:", isSingular(fit.sincos), "\n")
print(anova(fit.sincos))

# LRT: phase effect (ML)
fit.sincos.ml <- update(fit.sincos, REML = FALSE)
fit.nophase.ml <- lmerTest::lmer(
  magnitude ~ numStims_ord + betaLabels + baselineMag_c +
  (1 + dose_linpoly | sid) + (1 | channel),
  data = summaryNB_sincos, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))
cat("\nLRT: omnibus phase test:\n")
print(anova(fit.nophase.ml, fit.sincos.ml))

# emmeans at phase=90 and 270
emm_90 <- emmeans(fit.sincos, ~ numStims_ord,
  at = list(sin_phase = 1, cos_phase = 0, betaLabels = "0", baselineMag_c = 0))
emm_270 <- emmeans(fit.sincos, ~ numStims_ord,
  at = list(sin_phase = -1, cos_phase = 0, betaLabels = "0", baselineMag_c = 0))
cat("\nEmmeans at phase=90:\n"); print(emm_90)
cat("\nEmmeans at phase=270:\n"); print(emm_270)
cat("\nDose contrasts at 90:\n"); print(confint(pairs(emm_90)))
cat("\nDose contrasts at 270:\n"); print(confint(pairs(emm_270)))

# phase-response curve
phase_vals <- seq(0, 315, by = 45)
emm_curve <- lapply(phase_vals, function(ph) {
  em <- emmeans(fit.sincos, ~ numStims_ord,
    at = list(sin_phase = sin(ph*pi/180), cos_phase = cos(ph*pi/180),
              betaLabels = "0", baselineMag_c = 0))
  df <- as.data.frame(em)
  df$phase_deg <- ph
  df
})
emm_curve_df <- do.call(rbind, emm_curve)

p4 <- ggplot(emm_curve_df, aes(x = phase_deg, y = emmean, color = numStims_ord)) +
  theme_light(base_size = 14) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = numStims_ord), alpha = 0.15, color = NA) +
  labs(x = "Delivered Phase (degrees)",
       y = expression(paste("Predicted Magnitude (", mu, "V)")),
       color = "Dose", fill = "Dose",
       title = "Sin/Cos ANCOVA: Phase-Response Curve (emmeans +/- 95% CI)") +
  scale_x_continuous(breaks = seq(0, 315, by = 45))
ggsave(here("output_plots","betaStim_sincos_phase_curve.png"), plot = p4,
       units = "in", width = 7, height = 4.5, dpi = 600)

cat("\nAll 4 plots saved.\n")
