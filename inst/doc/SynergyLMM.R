## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", 
                      fig.align = "center", 
                      out.width = "100%",
                      prompt = TRUE)


## ----setup--------------------------------------------------------------------
library(SynergyLMM)

## -----------------------------------------------------------------------------
data("grwth_data")

## -----------------------------------------------------------------------------
head(grwth_data)

## -----------------------------------------------------------------------------
unique(grwth_data$Treatment)

## ----fig.width=12, fig.height=8-----------------------------------------------
lmm_ex <- lmmModel(data = grwth_data, sample_id = "subject", time = "Time", 
                   treatment = "Treatment", tumor_vol = "TumorVolume",
                   trt_control = "Control", drug_a = "DrugA", 
                   drug_b = "DrugB", combination = "Combination")

## -----------------------------------------------------------------------------
lmmModel_estimates(lmm_ex)

## ----fig.width=10, fig.height=10----------------------------------------------
ranefDiagnostics(lmm_ex)

## ----fig.width=10, fig.height=14----------------------------------------------
residDiagnostics(lmm_ex)

## ----fig.width=12, fig.height=8-----------------------------------------------
lmm_ex_var <- lmmModel(data = grwth_data, sample_id = "subject", time = "Time", 
                   treatment = "Treatment", tumor_vol = "TumorVolume",
                   trt_control = "Control", drug_a = "DrugA", 
                   drug_b = "DrugB", combination = "Combination",
                   weights = nlme::varIdent(form = ~1|SampleID))

## -----------------------------------------------------------------------------
lmmModel_estimates(lmm_ex_var)

## ----fig.width=10, fig.height=10----------------------------------------------
ranefD <- ranefDiagnostics(lmm_ex_var, verbose = FALSE)

# We can access to individual results of the diagnostics:
ranefD$Normality

## ----fig.width=10, fig.height=14----------------------------------------------
residD <- residDiagnostics(lmm_ex_var, verbose = FALSE)
residD$Normality

## ----fig.width=10, fig.height=10----------------------------------------------
ObsvsPred(lmm_ex_var, nrow = 8, ncol = 4)

## ----fig.width=10, fig.height=8-----------------------------------------------
CookDistance(lmm_ex_var)

## ----fig.width=10, fig.height=8-----------------------------------------------
logLikSubjectDisplacements(lmm_ex_var, var_name = "SampleID")

## ----fig.width=12, fig.height=10----------------------------------------------
bliss <- lmmSynergy(lmm_ex_var, method = "Bliss")

## ----error=TRUE---------------------------------------------------------------
try({
bliss <- lmmSynergy(lmm_ex_var, method = "Bliss", robust = TRUE)
})

## ----fig.width=12, fig.height=10----------------------------------------------
bliss <- lmmSynergy(lmm_ex_var, method = "Bliss", robust = TRUE, min_time = 6)

## -----------------------------------------------------------------------------
bliss$Synergy

## ----fig.width=12, fig.height=10----------------------------------------------
hsa <- lmmSynergy(lmm_ex_var, method = "HSA", robust = TRUE, min_time = 6)

## ----fig.width=12, fig.height=10----------------------------------------------
ra <- lmmSynergy(lmm_ex_var, method = "RA", robust = TRUE, min_time = 6, ra_nsim = 1000)

## -----------------------------------------------------------------------------
PostHocPwr(lmm_ex_var, nsim = 100, method = "Bliss")

## -----------------------------------------------------------------------------
# Vector with the time points
days <- unique(grwth_data$Time)

# Model estimates
estimates <- lmmModel_estimates(lmm_ex_var)

## ----fig.width=10, fig.height=8-----------------------------------------------
PwrSampleSize(npg = 1:10,
              time = days,
              grwrControl = round(estimates$Control,3),
              grwrA = round(estimates$DrugA,3),
              grwrB = round(estimates$DrugB, 3),
              grwrComb = round(estimates$Combination, 3),
              sd_ranef = round(estimates$sd_ranef, 3),
              sgma = round(estimates$sd_resid, 3),
              method = "Bliss")

## -----------------------------------------------------------------------------
max_time <- list(seq(0,9,3), seq(0,12,3), seq(0,15,3), 
                 seq(0,18,3), seq(0,21,3), seq(0,24,3), 
                 seq(0,27,3), seq(0,30,3))

## -----------------------------------------------------------------------------
# We can calculate the average sample size dividing the number of subjects
# by the number of groups, in this case, 4 groups
(npg <- round(length(unique(grwth_data$subject))/4,0))

## ----fig.width=10, fig.height=8-----------------------------------------------
PwrTime(npg = npg,
        time = max_time,
        type = "max",
        grwrControl = round(estimates$Control,3),
              grwrA = round(estimates$DrugA,3),
              grwrB = round(estimates$DrugB, 3),
              grwrComb = round(estimates$Combination, 3),
              sd_ranef = round(estimates$sd_ranef, 3),
              sgma = round(estimates$sd_resid, 3),
              method = "Bliss")

## -----------------------------------------------------------------------------
freq_time <- list(seq(0,18,1), seq(0,18,3), seq(0,18,6), seq(0,18,9),seq(0,18,18))

## ----fig.width=10, fig.height=8-----------------------------------------------
PwrTime(npg = npg,
        time = freq_time,
        type = "freq",
        grwrControl = round(estimates$Control,3),
        grwrA = round(estimates$DrugA,3),
        grwrB = round(estimates$DrugB, 3),
        grwrComb = round(estimates$Combination, 3),
        sd_ranef = round(estimates$sd_ranef, 3),
        sgma = round(estimates$sd_resid, 3),
        method = "Bliss")

## -----------------------------------------------------------------------------
estimates

## ----fig.width=14, fig.height=8-----------------------------------------------
APrioriPwr(npg = npg, # Sample size per group, calculated above
           time = days, # Time points of measurements, calculated above
           # Model estimates:
           grwrControl = round(estimates$Control,3),
           grwrA = round(estimates$DrugA,3),
           grwrB = round(estimates$DrugB, 3),
           grwrComb = round(estimates$Combination, 3),
           sd_ranef = round(estimates$sd_ranef, 3),
           sgma = round(estimates$sd_resid, 3),
           sd_eval = seq(0.01, 0.1, 0.01),
           sgma_eval = seq(0.01, 1, 0.01),
           grwrComb_eval = seq(-0.05, 0.1, 0.001)
           )
           

