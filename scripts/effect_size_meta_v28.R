suppressMessages({ library(metafor); library(jsonlite) })

## ---------------------------------------------------------------------------
## v2.8 module-level meta-analysis (full transparency version)
##  - GSE155776, GSE2513 (RMA), GSE183153 (voom logCPM): independent Hedges' g
##  - GSE51995: paired SMCRP (raw-score SD denominator, r-based variance)
##  - Primary: REML + Knapp-Hartung; Sensitivity: REML+adhoc(mKH), DL, exclude-GSE51995
##  - BH across 19 modules; per-cohort g/se; pair_r; module-level LODO
## ---------------------------------------------------------------------------

T <- "/mnt/c/Users/tidus/PycharmProjects/PterygiumAnalysis/results/paper_a/tables"
OUT <- "/mnt/c/Users/tidus/PycharmProjects/PterygiumAnalysis/results/paper_a/tables"
COHORTS <- c("GSE155776", "GSE51995", "GSE2513", "GSE183153")

## ---- load scores ----
## NOTE on the UV_response columns (fixed 2026-10-18):
##   module_scores_GSE2513.csv / GSE183153.csv were rebuilt on 2026-08-12 by
##   rebuild_gsva_reprocessed.R and ALREADY contain UV_response_UP / UV_response_DN.
##   The older module_scores_UVsplit_*.csv files (2026-08-11) carry the SAME column
##   names but different values (up to 0.28 apart on individual samples). The previous
##   code did a bare cbind(), so R silently kept the first (BASE) copy for GSE2513 and
##   GSE183153 while GSE155776 / GSE51995 (which only have the combined UV_response
##   column) fell back to the UVsplit copy -- and the Python audit script
##   export_audit_cohort.py made the opposite choice. That is what made the sensitivity
##   tables disagree with the primary analysis.
##   Fix: take non-UV columns from BASE, then append only the UV columns that BASE does
##   NOT already provide, and assert there are no duplicate names.
load_scores <- function(c) {
  base <- read.csv(file.path(T, paste0("module_scores_", c, ".csv")), row.names = 1, check.names = FALSE)
  uv <- read.csv(file.path(T, paste0("module_scores_UVsplit_", c, ".csv")), row.names = 1, check.names = FALSE)
  base <- base[, setdiff(colnames(base), "UV_response"), drop = FALSE]
  uv <- uv[, setdiff(colnames(uv), colnames(base)), drop = FALSE]
  common_rn <- intersect(rownames(base), rownames(uv))
  m <- cbind(base[common_rn, , drop = FALSE], uv[common_rn, , drop = FALSE])
  if (any(duplicated(colnames(m)))) {
    stop("duplicate module columns after cbind for ", c, ": ",
         paste(colnames(m)[duplicated(colnames(m))], collapse = ", "))
  }
  m
}
scores <- lapply(setNames(COHORTS, COHORTS), load_scores)
modules <- setdiff(colnames(scores$GSE155776), "UV_response")
cat("modules:", length(modules), "\n")

## ---- group maps ----
gm <- function(c) {
  rn <- rownames(scores[[c]])
  if (c == "GSE155776") return(ifelse(grepl("^Pterygium", rn), "pterygium", "control"))
  if (c == "GSE51995") {
    sh <- read.csv(file.path(T, "sample_sheets/GSE51995_sample_sheet.csv"))
    sid <- as.character(unlist(sh[["sample_id"]]))
    grp <- as.character(unlist(sh[["group"]]))
    m <- setNames(grp, sid)
    g <- unname(m[rn]); g[is.na(g)] <- "exclude"; g[g == "normal_conjunctiva"] <- "control"; return(g)
  }
  if (c == "GSE2513") {
    # correct mapping from OSNS_FES: pterygium = GSM48029-48036; control = GSM48026-28, GSM48037
    num <- as.numeric(sub("GSM", "", rn))
    return(ifelse(num >= 48029 & num <= 48036, "pterygium",
           ifelse(num %in% c(48026, 48027, 48028, 48037), "control", "exclude")))
  }
  if (c == "GSE183153") {
    num <- as.numeric(sub("GSM", "", rn))
    g <- ifelse(num >= 5552155, "pterygium", ifelse(num <= 5552152, "control", "exclude"))
    return(g)
  }
  stop("unknown cohort")
}

## ---- per-cohort Hedges' g (independent) via metafor escalc ----
per_effect <- function(mod, c) {
  sc <- scores[[c]]; g <- gm(c)
  pt <- sc[g == "pterygium", mod]; ct <- sc[g == "control", mod]
  if (length(pt) < 2 || length(ct) < 2) return(c(g = NA, se = NA))
  if (sd(pt) == 0 || sd(ct) == 0 || !is.finite(sd(pt)) || !is.finite(sd(ct))) return(c(g = NA, se = NA))
  e <- tryCatch(escalc(measure = "SMD", m1i = mean(pt), sd1i = sd(pt), n1i = length(pt),
                       m2i = mean(ct), sd2i = sd(ct), n2i = length(ct)), error = function(e) NULL)
  if (is.null(e)) return(c(g = NA, se = NA))
  c(g = unlist(e$yi)[1], se = sqrt(unlist(e$vi)[1]))
}

## ---- GSE51995 paired SMCRP ----
sh <- read.csv(file.path(T, "sample_sheets/GSE51995_sample_sheet.csv"))
pts <- setNames(as.character(unlist(sh[["sample_id"]]))[unlist(sh[["group"]]) == "pterygium"],
                as.character(unlist(sh[["patient_id"]]))[unlist(sh[["group"]]) == "pterygium"])
cts <- setNames(as.character(unlist(sh[["sample_id"]]))[unlist(sh[["group"]]) == "normal_conjunctiva"],
                as.character(unlist(sh[["patient_id"]]))[unlist(sh[["group"]]) == "normal_conjunctiva"])
common <- intersect(names(pts), names(cts))

paired_effect <- function(mod, r_override = NULL) {
  sc <- scores$GSE51995
  x <- sapply(common, function(p) sc[pts[p], mod])
  y <- sapply(common, function(p) sc[cts[p], mod])
  if (any(!is.finite(x)) || any(!is.finite(y)) || sd(x) == 0 || sd(y) == 0)
    return(c(g = NA, se = NA, r = NA))
  r <- if (is.null(r_override)) cor(x, y) else r_override
  e <- tryCatch(escalc(measure = "SMCRP", m1i = mean(x), m2i = mean(y),
                       sd1i = sd(x), sd2i = sd(y), ni = length(x), ri = r),
                error = function(e) NULL)
  if (is.null(e)) return(c(g = NA, se = NA, r = r))
  c(g = unlist(e$yi)[1], se = sqrt(unlist(e$vi)[1]), r = r)
}

## ---- run all ----
results <- list()
for (mod in modules) {
  per <- lapply(COHORTS, function(c) per_effect(mod, c))
  names(per) <- COHORTS
  pp <- paired_effect(mod)
  per$GSE51995 <- c(pp[1], pp[2])

  gs <- sapply(COHORTS, function(c) unname(per[[c]]["g"]))
  vs <- sapply(COHORTS, function(c) unname(per[[c]]["se"])^2)
  if (any(!is.finite(gs)) || any(!is.finite(vs))) {
    cat("skip", mod, "due to NA:", gs, "\n")
    results[[mod]] <- NULL
    next
  }

  fit_kh  <- rma(yi = gs, vi = vs, method = "REML", test = "knha")
  fit_mkh <- rma(yi = gs, vi = vs, method = "REML", test = "adhoc")
  fit_dl  <- rma(yi = gs, vi = vs, method = "DL", test = "z")

  rem <- COHORTS[COHORTS != "GSE51995"]
  fit_x <- rma(yi = gs[rem], vi = vs[rem], method = "REML", test = "knha")

  results[[mod]] <- list(
    module = mod,
    per_cohort = per,
    pair_r = pp["r"],
    kh = fit_kh, mkh = fit_mkh, dl = fit_dl, x = fit_x
  )
}

## ---- assemble table ----
rows <- lapply(results, function(R) {
  f <- R$kh; m <- R$mkh; d <- R$dl; x <- R$x
  ci_kh <- predict(f); ci_x <- predict(x)
  data.frame(
    module = R$module,
    pair_r_GSE51995 = round(as.numeric(R$pair_r), 3),
    kh_g = f$b[1,1], kh_se = f$se, kh_ci_low = ci_kh$ci.lb, kh_ci_high = ci_kh$ci.ub,
    kh_p = f$pval, kh_tau2 = f$tau2, kh_I2 = f$I2/100, kh_Q = f$QE, kh_Qp = f$QEp, kh_k = f$k,
    mkh_g = m$b[1,1], mkh_p = m$pval, mkh_tau2 = m$tau2, mkh_I2 = m$I2/100,
    dl_g = d$b[1,1], dl_se = d$se, dl_ci_low = d$ci.lb, dl_ci_high = d$ci.ub, dl_p = d$pval,
    dl_tau2 = d$tau2, dl_I2 = d$I2/100, dl_Q = d$QE, dl_Qp = d$QEp,
    x_g = x$b[1,1], x_ci_low = ci_x$ci.lb, x_ci_high = ci_x$ci.ub, x_p = x$pval, x_I2 = x$I2/100,
    GSE155776_g = R$per_cohort$GSE155776["g"], GSE155776_se = R$per_cohort$GSE155776["se"],
    GSE51995_g = R$per_cohort$GSE51995["g"], GSE51995_se = R$per_cohort$GSE51995["se"],
    GSE2513_g = R$per_cohort$GSE2513["g"], GSE2513_se = R$per_cohort$GSE2513["se"],
    GSE183153_g = R$per_cohort$GSE183153["g"], GSE183153_se = R$per_cohort$GSE183153["se"]
  )
})
tab <- do.call(rbind, rows)
tab$q_BH19 <- p.adjust(tab$kh_p, method = "BH")
tab <- tab[order(tab$kh_g, decreasing = TRUE), ]
write.csv(tab, file.path(OUT, "module_effect_meta_v28.csv"), row.names = FALSE)

## ---- module-level LODO (drop each cohort, primary REML+KH) ----
lodo <- list()
for (mod in modules) {
  R <- results[[mod]]
  for (drop in COHORTS) {
    keep <- COHORTS[COHORTS != drop]
    gs <- sapply(keep, function(c) unname(R$per_cohort[[c]]["g"]))
    vs <- sapply(keep, function(c) unname(R$per_cohort[[c]]["se"])^2)
    fit <- tryCatch(rma(yi = gs, vi = vs, method = "REML", test = "knha"), error = function(e) NULL)
    if (!is.null(fit)) {
      ci <- predict(fit)
      lodo[[paste(mod, drop)]] <- data.frame(module = mod, dropped = drop,
        pooled_g = fit$b[1,1], ci_low = ci$ci.lb, ci_high = ci$ci.ub,
        tau2 = fit$tau2, I2 = fit$I2/100, p = fit$pval)
    }
  }
}
lodo_df <- do.call(rbind, lodo)
write.csv(lodo_df, file.path(OUT, "module_lodo_v28.csv"), row.names = FALSE)

## ---- sensitivity: r sweep for GSE51995 (SMCRP with r=0/0.25/0.5/0.75) ----
r_sweep <- list()
for (mod in modules) {
  for (rv in c(0, 0.25, 0.5, 0.75)) {
    pp <- paired_effect(mod, rv)
    gs <- c(sapply(c("GSE155776", "GSE2513", "GSE183153"), function(c) unname(per_effect(mod, c)["g"])), pp["g"])
    vs <- c(sapply(c("GSE155776", "GSE2513", "GSE183153"), function(c) unname(per_effect(mod, c)["se"])^2), pp["se"]^2)
    fit <- rma(yi = gs, vi = vs, method = "REML", test = "knha")
    r_sweep[[paste(mod, rv)]] <- data.frame(module = mod, r = rv,
      pooled_g = fit$b[1,1], p = fit$pval, I2 = fit$I2/100)
  }
}
r_sweep_df <- do.call(rbind, r_sweep)
write.csv(r_sweep_df, file.path(OUT, "module_sensitivity_rsweep_v28.csv"), row.names = FALSE)

cat("\n=== FINAL (REML+KH primary, sorted) ===\n")
for (i in seq_len(nrow(tab))) {
  flag <- ifelse(tab$kh_p[i] < 0.05, " *", "")
  cat(sprintf("%-22s g=%+.2f [%+.2f,%+.2f] I2=%3.0f%% p=%.2e q=%.3f%s\n",
      tab$module[i], tab$kh_g[i], tab$kh_ci_low[i], tab$kh_ci_high[i],
      tab$kh_I2[i]*100, tab$kh_p[i], tab$q_BH19[i], flag))
}
cat("\nWROTE module_effect_meta_v28.csv, module_lodo_v28.csv, module_sensitivity_rsweep_v28.csv\n")
