## PLAUSL SAP v0.1 - Aim 1 analyses (sections 7.1-7.6) and Aim 3 data check (8.0)
## Usage (study data): from the data folder (final_plausl_data/), or with PLAUSL_DATA set to it:
##   Rscript <path to>/sap_aim1.R > sap_analysis/sap_aim1_output.txt   (writes sap_analysis/sap_analysis_dataset.csv)
## Usage (shared de-identified data, e.g. the Zenodo files): set PLAUSL_SHARED to the folder holding
##   plausl_participants.csv and plausl_sl_trials.csv:  PLAUSL_SHARED=<folder> Rscript <path to>/sap_aim1.R
##   (writes sap_analysis_dataset.csv to the current folder; the Aim 3 REDCap check is skipped)
## Needs R >= 4.4 with readxl, readr, jsonlite, lme4, lmerTest, sandwich, dplyr, tidyr. Firth regression, AUC and LOOCV are coded here.
## SL composite (from 2026-10-16): per task mean of z(2AFC accuracy) and z(SICR bigram recall on target trials), z within
## session; every model with SL as predictor adjusts for SICR bigram recall on foil trials (foil1/foil2, see build_sl).
suppressMessages({library(readxl); library(lme4); library(lmerTest); library(sandwich); library(dplyr)})
set.seed(2026)
options(width = 200, warn = -1)
SHARED <- nzchar(Sys.getenv("PLAUSL_SHARED"))      # TRUE: run on the shared de-identified files
D <- normalizePath(Sys.getenv(if (SHARED) "PLAUSL_SHARED" else "PLAUSL_DATA", "."))   # data folder
# private participant-specific corrections (study data only; the shared files already have invalid scores removed)
CFG <- if (SHARED) NULL else jsonlite::fromJSON(file.path(D, "study_config.json"), simplifyVector = FALSE)
hr <- function(s) cat("\n", strrep("=", 100), "\n", s, "\n", strrep("=", 100), "\n", sep = "")
num <- function(v) suppressWarnings(as.numeric(v))
z <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
rmean <- function(...) { x <- cbind(...); r <- rowMeans(x, na.rm = TRUE); r[is.nan(r)] <- NA; r }

## ---------------------------------------------------------------- data
if (SHARED) {
  m <- read.csv(file.path(D, "plausl_participants.csv"))
  m <- m %>% mutate(record_id = participant, pid = participant, age = age_years,
                    gender = ifelse(sex == "male", "M", ifelse(sex == "female", "F", NA)))
} else {
  m <- read_excel(file.path(D, "PLAUSL Data - Clinical Trials with SL scores.xlsx"), sheet = "Master Sheet")
  dem <- read.csv(file.path(D, "clinical_trials_with_demographics.csv"))
  m <- m %>% left_join(dem %>% transmute(record_id, pid = demog_record_id, age = num(age_at_consent)), by = "record_id")
}
a <- m %>% filter(completion_status == "completed", !is.na(diagnosis_group)) %>% as.data.frame()
for (v in setdiff(names(a), c("record_id", "diagnosis_group", "completion_status", "gender", "pid", "hispanic/latinx",
                              "participant", "sex", "race", "ethnicity")))
  if (!is.numeric(a[[v]])) a[[v]] <- num(a[[v]])
a$asd  <- as.integer(a$diagnosis_group == "ASD")
a$male <- ifelse(a$gender == "M", 1, ifelse(a$gender == "F", 0, NA))
a$nviq <- a$kbit_ss
a$lang1 <- a$language_composite
a$read1 <- a$reading_composite
a$lang2 <- a$language_composite_2                   # T2 language: Picture Vocabulary + TROG-2 (TILLS not given at T2)
comp2 <- function(df, cols) { x <- sapply(df[cols], num); n <- rowSums(!is.na(x)); r <- rowMeans(x, na.rm = TRUE); r[n < 2] <- NA; r }
a$read2 <- comp2(a, c("wrmt_letterid_ss_2", "wrmt_id_ss_2", "wrmt_wa_ss_2", "orr_ss_2"))
inv <- a   # reading T1 without letter ID: drop scores the config lists as invalid
for (x in CFG$merge$invalid_scores) if (x[[2]] %in% names(inv)) inv[[x[[2]]]][inv$record_id == x[[1]]] <- NA  # no-op when SHARED
a$read1_noLID <- comp2(inv, c("wrmt_id_ss", "wrmt_wa_ss", "orr_ss"))
a$read2_noLID <- comp2(a, c("wrmt_id_ss_2", "wrmt_wa_ss_2", "orr_ss_2"))
a$disorder <- a$tills_languagedisorder
a$disorder[is.na(a$tills_coresum)] <- NA          # a flag without a core sum is not used

# SL composites, each measure z-scored within session across the analysis set.
# Primary ("bigram"): task composite = mean of z(2AFC accuracy) and z(SICR bigram recall, target trials).
# SICR bigram recall on foil trials is returned as a covariate (foil1/foil2: mean of the two tasks' z-scores,
# built like the composite so children missing one task's SICR are kept; ssl_bfoil/lsl_bfoil per task).
# Earlier versions, kept for sensitivity analyses: "tmf" = serial target - foil (SAP draft), "target" = serial target.
build_sl <- function(df, recall = c("bigram", "tmf", "target"), drop_versionE = FALSE) {
  recall <- match.arg(recall); out <- list()
  # children whose letter SICR run came from task version E (JGA shown in recall), by session
  vE_ids <- function(s) if (SHARED) df$pid[which(df[[sprintf("lsl_sicr_versionE_s%d", s)]] == 1)] else {
    vE <- read.csv(file.path(D, "plausl_sl_sources.csv")) %>% filter(used == "yes", grepl("version E", note), measure == "lsl_sicr")
    vE$participant[vE$session == s] }
  for (s in 1:2) {
    tc <- list(); fc <- list()
    for (t in c("ssl", "lsl")) {
      afc <- df[[sprintf("%s_2afc_acc_s%d", t, s)]]
      tg <- df[[sprintf("%s_sicr_serial_target_s%d", t, s)]]; fl <- df[[sprintf("%s_sicr_serial_foil_s%d", t, s)]]
      bt <- df[[sprintf("%s_sicr_bigram_target_s%d", t, s)]]; bf <- df[[sprintf("%s_sicr_bigram_foil_s%d", t, s)]]
      idx <- switch(recall, bigram = bt, tmf = tg - fl, target = tg)
      if (drop_versionE && t == "lsl") { drop <- df$pid %in% vE_ids(s); idx[drop] <- NA; bf[drop] <- NA }
      out[[sprintf("%s_afc%d", t, s)]] <- afc; out[[sprintf("%s_idx%d", t, s)]] <- idx; out[[sprintf("%s_bfoil%d", t, s)]] <- bf
      tc[[t]] <- rmean(z(afc), z(idx)); out[[sprintf("%s_comp%d", t, s)]] <- tc[[t]]
      fc[[t]] <- z(bf)
    }
    out[[sprintf("sl%d", s)]] <- rmean(tc$ssl, tc$lsl)
    out[[sprintf("foil%d", s)]] <- rmean(fc$ssl, fc$lsl)
  }
  as.data.frame(out)
}
a <- cbind(a, build_sl(a))
COV <- "+ foil1 + asd + male + nviq"          # covariates for SL_T1 models (bigram foil recall, group, sex, NVIQ)
COVg <- "+ foil1 + male + nviq"               # within-group models
fm <- function(...) as.formula(paste0(...))
cat("Analysis set:", nrow(a), "children (", sum(a$asd), "ASD,", sum(1 - a$asd), "non-ASD )\n")

## ---------------------------------------------------------------- helpers
# standardized effect of `focal` (focal and outcome z-scored within the model's complete cases)
fitlm <- function(fml, data, focal, label = "") {
  vars <- all.vars(fml); d <- data[complete.cases(data[vars]), vars, drop = FALSE]
  y <- vars[1]; d[[y]] <- z(d[[y]]); d[[focal]] <- z(d[[focal]])
  f <- lm(fml, d); n <- nrow(d); co <- summary(f)$coefficients
  if (!focal %in% rownames(co)) return(NULL)
  ci <- confint(f)[focal, ]
  red <- update(f, as.formula(paste(". ~ . -", focal)))
  f2 <- (summary(f)$r.squared - summary(red)$r.squared) / (1 - summary(f)$r.squared)
  hc <- tryCatch({ se <- sqrt(diag(vcovHC(f, type = "HC3")))[focal]; 2 * pt(-abs(coef(f)[focal] / se), f$df.residual) }, error = function(e) NA)
  data.frame(model = label, n = n, df = f$df.residual, beta = co[focal, 1], lo = ci[1], hi = ci[2], p = co[focal, 4],
             p_HC3 = hc, f2 = f2, row.names = NULL)
}
pr <- function(tab) { num_cols <- sapply(tab, is.numeric); tab[num_cols] <- lapply(tab[num_cols], round, 3); print(tab, row.names = FALSE) }

## ---------------------------------------------------------------- 7.1 learning + reliability
hr("7.1  Evidence of learning in each task (supporting analyses)")
tr <- read.csv(file.path(D, if (SHARED) "plausl_sl_trials.csv" else "plausl_sl_scores_by_trial.csv"))
tr <- tr %>% inner_join(a %>% select(pid, asd), by = c("participant" = "pid")) %>% mutate(asd_c = asd - 0.5)
cat("\n2AFC: logistic GLMM  score ~ group(effect-coded) + (1|child) + (1|trial); intercept = mean of groups vs chance\n")
res <- list()
for (ms in c("ssl_afc", "lsl_afc")) for (s in 1:2) {
  d <- tr %>% filter(measure == ms, session == s, !is.na(score))
  g <- glmer(score ~ asd_c + (1 | participant) + (1 | trial), d, family = binomial)
  co <- summary(g)$coefficients
  pm <- d %>% group_by(participant) %>% summarise(acc = mean(score))
  tt <- t.test(pm$acc, mu = .5)
  res[[length(res) + 1]] <- data.frame(task = ms, session = s, children = nrow(pm), trials = nrow(d),
    p_correct = plogis(co[1, 1]), p_vs_chance = co[1, 4], group_OR = exp(co[2, 1]), p_group = co[2, 4],
    ttest_mean = mean(pm$acc), ttest_p = tt$p.value, singular = isSingular(g))
}
pr(do.call(rbind, res))

cat("\nSICR: linear mixed model  serial/6 ~ trial type * group + (1 + trial type|child) + (1|item)\n")
res <- list()
for (ms in c("ssl_sicr", "lsl_sicr")) for (s in 1:2) for (outc in c("serial", "bigrams")) {
  d <- tr %>% filter(measure == ms, session == s, !is.na(tag)) %>%
    mutate(y = if (outc == "serial") serial / 6 else bigrams / 4, tag_c = tag - 0.5)
  g <- tryCatch(lmer(y ~ tag_c * asd_c + (1 + tag_c | participant) + (1 | key), d),
                error = function(e) lmer(y ~ tag_c * asd_c + (1 | participant) + (1 | key), d))
  if (isSingular(g)) g <- lmer(y ~ tag_c * asd_c + (1 | participant) + (1 | key), d)
  co <- summary(g)$coefficients
  res[[length(res) + 1]] <- data.frame(task = ms, session = s, score = outc, children = n_distinct(d$participant),
    target_minus_foil = co["tag_c", 1], p_type = co["tag_c", 5], typeXgroup = co["tag_c:asd_c", 1],
    p_typeXgroup = co["tag_c:asd_c", 5], singular = isSingular(g))
}
pr(do.call(rbind, res))

cat("\nReliability: split-half (mean over 1,000 random splits of trials, Spearman-Brown corrected) and test-retest (session 1 vs 2)\n")
sb <- function(r) 2 * r / (1 + r)
res <- list()
for (ms in c("ssl_afc", "lsl_afc", "ssl_sicr", "lsl_sicr")) for (s in 1:2) {
  d <- tr %>% filter(measure == ms, session == s)
  # SICR: bigram recall on target trials (the composite component)
  d <- if (grepl("afc", ms)) d %>% filter(!is.na(score)) %>% mutate(v = score, tag = 1) else d %>% filter(tag == 1) %>% mutate(v = bigrams / 4)
  items <- sort(unique(d$trial))
  rs <- replicate(1000, {
    h <- d %>% mutate(half = trial %in% sample(items, length(items) %/% 2))
    w <- h %>% group_by(participant, half) %>% summarise(v = mean(v), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = half, values_from = v, names_prefix = "h") %>% filter(!is.na(hTRUE), !is.na(hFALSE))
    if (nrow(w) > 3) cor(w$hTRUE, w$hFALSE) else NA })
  r <- mean(rs, na.rm = TRUE)
  res[[length(res) + 1]] <- data.frame(measure = ms, session = s, n = n_distinct(d$participant), mean_split_half_r = r, spearman_brown = sb(r))
}
pr(do.call(rbind, res))
rt <- do.call(rbind, lapply(c("ssl_afc", "lsl_afc", "ssl_idx", "lsl_idx", "ssl_comp", "lsl_comp", "sl"), function(v) {
  x <- a[[paste0(v, 1)]]; y <- a[[paste0(v, 2)]]; ok <- complete.cases(x, y)
  data.frame(measure = v, n = sum(ok), test_retest_r = if (sum(ok) > 3) cor(x[ok], y[ok]) else NA,
             p = if (sum(ok) > 3) cor.test(x[ok], y[ok])$p.value else NA)
}))
pr(rt)
cat(sprintf("\nComponent intercorrelations (T1): syllable 2AFC vs SICR r = %.2f; letter 2AFC vs SICR r = %.2f; syllable vs letter composite r = %.2f\n",
    cor(a$ssl_afc1, a$ssl_idx1, use = "pair"), cor(a$lsl_afc1, a$lsl_idx1, use = "pair"), cor(a$ssl_comp1, a$lsl_comp1, use = "pair")))

## ---------------------------------------------------------------- 7.2 H1A
hr("7.2  H1A: concurrent associations  Outcome_T1 ~ SL_T1 + bigram foil recall + group + sex + NVIQ  (beta = SD units)")
cat(sprintf("Children with an SL_T1 composite: %d; with the bigram foil covariate as well: %d\n", sum(!is.na(a$sl1)), sum(!is.na(a$sl1) & !is.na(a$foil1))))
conf <- list()
h1a <- rbind(
  fitlm(fm("lang1 ~ sl1", COV), a, "sl1", "Language T1, whole sample [confirmatory]"),
  fitlm(fm("read1 ~ sl1", COV), a, "sl1", "Reading T1, whole sample [confirmatory]"),
  fitlm(fm("lang1 ~ sl1", COVg), subset(a, asd == 1), "sl1", "Language T1, ASD"),
  fitlm(fm("lang1 ~ sl1", COVg), subset(a, asd == 0), "sl1", "Language T1, non-ASD"),
  fitlm(fm("read1 ~ sl1", COVg), subset(a, asd == 1), "sl1", "Reading T1, ASD"),
  fitlm(fm("read1 ~ sl1", COVg), subset(a, asd == 0), "sl1", "Reading T1, non-ASD"),
  fitlm(lang2 ~ sl2 + foil2 + asd + male + nviq, a, "sl2", "Language T2 composite ~ SL T2 [replication]"),
  fitlm(read2 ~ sl2 + foil2 + asd + male + nviq, a, "sl2", "Reading T2 ~ SL T2 [replication]"))
pr(h1a)
conf$lang_concurrent <- h1a$p[1]; conf$read_concurrent <- h1a$p[2]
cat("\nExploratory: single tasks and components (whole sample; SICR models adjust for that task's bigram foil recall)\n")
ex <- list()
for (o in c("lang1", "read1")) for (p in c("ssl_comp1", "lsl_comp1", "ssl_afc1", "ssl_idx1", "lsl_afc1", "lsl_idx1")) {
  foil <- if (grepl("afc", p)) "" else paste0("+ ", substr(p, 1, 3), "_bfoil1")
  ex[[length(ex) + 1]] <- fitlm(fm(o, " ~ ", p, foil, " + asd + male + nviq"), a, p, paste(o, "~", p, foil))
}
pr(do.call(rbind, ex))
cat("(ssl_idx1 / lsl_idx1 = SICR bigram recall on target trials)\n")

cat("\nLeave-one-out cross-validated prediction (1,000 outcome permutations)\n")
loocv <- function(fml, d) { yh <- sapply(seq_len(nrow(d)), function(i) predict(lm(fml, d[-i, ]), d[i, ])); y <- d[[all.vars(fml)[1]]]
  c(q2 = 1 - sum((y - yh)^2) / sum((y - mean(y))^2), rmse = sqrt(mean((y - yh)^2))) }
for (o in c("lang1", "read1")) {
  full <- fm(o, " ~ sl1 ", COV); base <- fm(o, " ~ foil1 + asd + male + nviq")
  d <- a[complete.cases(a[all.vars(full)]), all.vars(full)]
  obs <- loocv(full, d); b0 <- loocv(base, d)
  perm <- replicate(1000, { dd <- d; dd[[o]] <- sample(dd[[o]]); loocv(full, dd)["q2"] })
  slperm <- replicate(1000, { dd <- d; dd$sl1 <- sample(dd$sl1); loocv(full, dd)["q2"] })
  cat(sprintf("%s: n=%d  q2(full)=%.3f RMSE=%.2f | q2(covariates only)=%.3f RMSE=%.2f | p(outcome permutation)=%.3f | p(SL permutation)=%.3f\n",
      o, nrow(d), obs["q2"], obs["rmse"], b0["q2"], b0["rmse"], mean(perm >= obs["q2"]), mean(slperm >= obs["q2"])))
}

## ---------------------------------------------------------------- 7.3 H1B(a) Firth logistic
hr("7.3  H1B(a): disorder status (TILLS core < 24) ~ SL_T1 + bigram foil recall + sex + NVIQ, Firth penalized logistic regression")
firth <- function(fml, d, focal) {
  X <- model.matrix(fml, d); y <- model.response(model.frame(fml, d)); j <- which(colnames(X) == focal)
  pll <- function(b) { eta <- X %*% b; p <- plogis(eta); w <- as.vector(p * (1 - p))
    sum(y * eta - log1p(exp(eta))) + 0.5 * as.numeric(determinant(crossprod(X * sqrt(w)), logarithm = TRUE)$modulus) }
  fitb <- function(fix = NULL) { k <- ncol(X); st <- rep(0, k)
    if (is.null(fix)) { o <- optim(st, function(b) -pll(b), method = "BFGS", control = list(maxit = 2000)); list(b = o$par, ll = -o$value) }
    else { o <- optim(st[-j], function(bb) { b <- numeric(k); b[-j] <- bb; b[j] <- fix; -pll(b) }, method = "BFGS", control = list(maxit = 2000)); -o$value } }
  full <- fitb(); prof <- function(v) 2 * (full$ll - fitb(v)) - qchisq(.95, 1)
  bj <- full$b[j]; lo <- tryCatch(uniroot(prof, c(bj - 15, bj))$root, error = function(e) NA); hi <- tryCatch(uniroot(prof, c(bj, bj + 15))$root, error = function(e) NA)
  p <- pchisq(2 * (full$ll - fitb(0)), 1, lower.tail = FALSE)
  pr_hat <- plogis(X %*% full$b); auc <- mean(outer(pr_hat[y == 1], pr_hat[y == 0], ">") + 0.5 * outer(pr_hat[y == 1], pr_hat[y == 0], "=="))
  data.frame(n = nrow(X), events = sum(y), OR_perSD = exp(bj), lo = exp(lo), hi = exp(hi), p_PLR = p, AUC = auc)
}
dd <- a[complete.cases(a[c("disorder", "sl1", "foil1", "male", "nviq")]), ]; dd$slz <- z(dd$sl1); dd$foilz <- z(dd$foil1)
cat("Disorder by group in this model sample:\n"); print(table(group = ifelse(dd$asd == 1, "ASD", "non-ASD"), disorder = dd$disorder))
pr(rbind(cbind(model = "disorder ~ SL + foil + sex + NVIQ", firth(disorder ~ slz + foilz + male + nviq, dd, "slz")),
         cbind(model = "disorder ~ SL + foil", firth(disorder ~ slz + foilz, dd, "slz"))))
cat("(Within-group models not fitted: too few ASD children with a core score, per SAP 7.3)\n")

## ---------------------------------------------------------------- 7.4 H1B(b) growth
hr("7.4  H1B(b): growth  Outcome_T2 ~ Outcome_T1 + SL_T1 + bigram foil recall + group + sex + NVIQ, and cross-lagged comparison")
gr <- rbind(
  fitlm(fm("read2 ~ read1 + sl1", COV), a, "sl1", "Reading T2 ~ SL T1 [confirmatory]"),
  fitlm(fm("lang2 ~ lang1 + sl1", COV), a, "sl1", "Language T2 composite ~ SL T1 [confirmatory]"),
  fitlm(read2 ~ read1 + sl1 + foil1, a, "sl1", "Reading growth, minimal model (foil only)"),
  fitlm(lang2 ~ lang1 + sl1 + foil1, a, "sl1", "Language growth, minimal model (foil only)"),
  fitlm(lang2 ~ picvocab_ss + sl1 + foil1, a, "sl1", "PicVocab T2 ~ PicVocab T1 + SL T1 (same-test baseline)"))
pr(gr)
conf$read_growth <- gr$p[1]; conf$lang_growth <- gr$p[2]
cat("\nCross-lagged: standardized SL_T1 -> Outcome_T2 vs Outcome_T1 -> SL_T2 (same children; bootstrap 2,000)\n")
xl <- function(o1, o2, d) {
  vs <- c(o1, o2, "sl1", "sl2", "foil1", "foil2", "asd", "male", "nviq"); d <- d[complete.cases(d[vs]), vs]
  for (v in c(o1, o2, "sl1", "sl2")) d[[v]] <- z(d[[v]])
  est <- function(dd) { f1 <- lm(fm(o2, " ~ ", o1, " + sl1 + foil1 + asd + male + nviq"), dd)
    f2 <- lm(fm("sl2 ~ sl1 + ", o1, " + foil1 + foil2 + asd + male + nviq"), dd)
    c(coef(f1)["sl1"], coef(f2)[o1]) }
  e <- est(d); bs <- t(replicate(2000, tryCatch(est(d[sample(nrow(d), replace = TRUE), ]), error = function(err) c(NA, NA))))
  dif <- bs[, 1] - bs[, 2]
  cat(sprintf("%s: n=%d  SL->outcome beta=%.2f | outcome->SL beta=%.2f | difference=%.2f, 95%% bootstrap CI [%.2f, %.2f]\n",
      o2, nrow(d), e[1], e[2], e[1] - e[2], quantile(dif, .025, na.rm = TRUE), quantile(dif, .975, na.rm = TRUE)))
}
xl("read1", "read2", a); xl("lang1", "lang2", a)

hr("6.3  Holm correction across the four confirmatory tests (SL_T1 composite)")
pv <- unlist(conf); print(data.frame(test = names(pv), p = round(pv, 4), p_holm = round(p.adjust(pv, "holm"), 4)), row.names = FALSE)

## ---------------------------------------------------------------- 7.5 H1C
hr("7.5  H1C (exploratory): group x SL interactions")
ix <- function(fml, d, term, label) { vars <- all.vars(fml); d <- d[complete.cases(d[vars]), vars]
  f <- lm(fml, d); co <- summary(f)$coefficients
  if (!term %in% rownames(co)) return(data.frame(model = label, n = nrow(d), df = f$df.residual, b = NA, lo = NA, hi = NA, p = NA))
  ci <- confint(f)[term, ]; data.frame(model = label, n = nrow(d), df = f$df.residual, b = co[term, 1], lo = ci[1], hi = ci[2], p = co[term, 4]) }
pr(rbind(
  ix(read2 ~ asd * (sl1 + read1) + foil1 + male + nviq, a, "asd:sl1", "Reading T2 ~ group*(SL T1 + Reading T1) [SAP model]"),
  ix(lang2 ~ asd * (sl1 + lang1) + foil1 + male + nviq, a, "asd:sl1", "Language T2 ~ group*(SL T1 + Language T1) [SAP model]"),
  ix(sl2 ~ asd * (sl1 + read1) + foil1 + foil2 + male + nviq, a, "asd:read1", "SL T2 ~ group*(SL T1 + Reading T1): group x Reading"),
  ix(lang1 ~ asd * sl1 + foil1 + male + nviq, a, "asd:sl1", "Language T1 ~ group*SL T1 (concurrent, added)"),
  ix(read1 ~ asd * sl1 + foil1 + male + nviq, a, "asd:sl1", "Reading T1 ~ group*SL T1 (concurrent, added)")))
cat("b = difference in slope (outcome points per composite unit), ASD minus non-ASD\n")

## ---------------------------------------------------------------- 7.6 sensitivity
hr("7.6  Sensitivity analyses for the confirmatory models (SL_T1 composite beta)")
sens <- function(d, label, lang = "lang1", read = "read1", cov = COV) rbind(
  fitlm(fm(lang, " ~ sl1 ", cov), d, "sl1", paste(label, "| language T1")),
  fitlm(fm(read, " ~ sl1 ", cov), d, "sl1", paste(label, "| reading T1")),
  fitlm(fm("read2 ~ ", read, " + sl1 ", cov), d, "sl1", paste(label, "| reading growth")))
swap <- function(newsl) { d <- a; d[names(newsl)] <- newsl; d }
b_tar <- build_sl(a, "target"); b_tmf <- build_sl(a, "tmf"); b_vE <- build_sl(a, "bigram", drop_versionE = TRUE)
NOFOIL <- "+ asd + male + nviq"
out3 <- function(v, d) { s <- d[[v]]; any_out <- rep(FALSE, nrow(d))
  for (g in 0:1) { i <- d$asd == g; zz <- (s[i] - mean(s[i], na.rm = TRUE)) / sd(s[i], na.rm = TRUE); any_out[i] <- !is.na(zz) & abs(zz) > 3 }; any_out }
noout <- a[!(out3("sl1", a) | out3("lang1", a) | out3("read1", a)), ]
cat("Children >3 SD from their group mean on SL T1, language T1 or reading T1:", a$record_id[out3("sl1", a) | out3("lang1", a) | out3("read1", a)], "\n")
pr(rbind(
  sens(a, "Primary"),
  sens(a, "Separate syllable and letter bigram-foil covariates", cov = "+ ssl_bfoil1 + lsl_bfoil1 + asd + male + nviq"),
  sens(a, "No foil covariate", cov = NOFOIL),
  sens(swap(b_tar[c("sl1", "sl2")]), "Earlier composite: SICR serial target recall, no foil (PRS v1)", cov = NOFOIL),
  sens(swap(b_tmf[c("sl1", "sl2")]), "Earlier composite: SICR serial target - foil (SAP draft)", cov = NOFOIL),
  sens(swap(b_vE[c("sl1", "sl2", "foil1", "foil2")]), "Letter SICR version E runs excluded"),
  sens(subset(a, is.na(nviq) | nviq >= 70), "KBIT-2 < 70 excluded"),
  sens(a, "Reading without WRMT letter ID", read = "read1_noLID"),
  fitlm(fm("lang1 ~ sl1 ", COV, " + age"), a, "sl1", "Age added | language T1"),
  fitlm(fm("read1 ~ sl1 ", COV, " + age"), a, "sl1", "Age added | reading T1"),
  sens(noout, "Outliers (>3 SD within group) excluded"),
  sens(swap(data.frame(sl1 = a$ssl_comp1, foil1 = z(a$ssl_bfoil1))), "Syllable composite only (syllable foil)"),
  sens(swap(data.frame(sl1 = a$lsl_comp1, foil1 = z(a$lsl_bfoil1))), "Letter composite only (letter foil)")))
cat("\nSpearman correlations (unadjusted), SL T1 composite:\n")
for (o in c("lang1", "read1")) { ok <- complete.cases(a$sl1, a[[o]]); ct <- cor.test(a$sl1[ok], a[[o]][ok], method = "spearman", exact = FALSE)
  cat(sprintf("  %s: rho=%.2f p=%.3f n=%d\n", o, ct$estimate, ct$p.value, sum(ok))) }
cat("\nMultiple imputation (SAP 6.4) not run: the 'mice' package is not installed.\n")

## ---------------------------------------------------------------- Aim 3
hr("8.0  Aim 3 data availability")
if (SHARED) {
  cat("Skipped: the serial recall (PSRT/OSRT) data are not part of the shared dataset.\n")
} else {
  r1 <- readr::read_csv(file.path(dirname(D), "beh", "plausl_time1_assessment_scores.csv"), col_types = readr::cols(.default = "c"), progress = FALSE)
  cat("Children with a PSRT/OSRT administration date in REDCap T1:", r1$record_id[!is.na(r1$psrt_dateadmin) | !is.na(r1$osrt_dateadmin)], "\n")
  cat(CFG$analysis$aim3_note, "\n", sep = "")
}
cat("SAP rule (>= 10 per group with usable T1 data) not met -> no inferential Aim 3 analyses.\n")
write.csv(a, if (SHARED) "sap_analysis_dataset.csv" else file.path(D, "sap_analysis", "sap_analysis_dataset.csv"), row.names = FALSE)
