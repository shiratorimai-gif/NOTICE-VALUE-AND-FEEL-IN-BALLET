# ============================================================
# 02. Domain-level CFA model comparisons
#
# Prerequisite:
#   Run 01_efa_item_selection_and_diagnostics.R first so that
#   the cleaned analytic data frame `df_raw` is available.
# ============================================================

pkgs <- c("dplyr", "tibble", "lavaan", "writexl")
invisible(lapply(pkgs, function(p) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)))
invisible(lapply(pkgs, library, character.only = TRUE))
options(contrasts = c("contr.sum", "contr.poly"))

# ============================================================
# 1. 項目セット
# ============================================================

att_salience <- c("q6s2", "q6s3", "q6s4", "q6s6")
att_body <- c("q6s7", "q6s8", "q6s9", "q6s10", "q6s11", "q6s12")
att_comp <- c("q6s13", "q6s14", "q6s15", "q6s17")
attention_items_48 <- c(att_salience, att_body, att_comp)

eval_tech <- c("q7s1", "q7s2", "q7s3", "q7s5", "q7s6", "q7s7", "q7s8")
eval_meaning <- c("q7s9", "q7s10", "q7s11", "q7s12", "q7s13", "q7s14", "q7s15", "q7s16", "q7s17")
eval_sim <- c("q7s19", "q7s20", "q7s21")
eval_comp <- c("q7s22", "q7s23", "q7s24")
evaluation_core_items_22 <- c(eval_tech, eval_meaning, eval_sim, eval_comp)

all_cfa_items <- c(attention_items_48, evaluation_core_items_22)

setdiff(all_cfa_items, names(df_raw))
intersect(c("q6s1", "q6s5_reverse", "q6s16_reverse", "q7s4", "q7s18"), all_cfa_items)

# ============================================================
# 2. lavaan用データ
# ============================================================

df_cfa <- df_raw %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(all_cfa_items), ~ ordered(as.numeric(as.character(.)), levels = 1:6)))

# ============================================================
# 3. CFAモデル
# ============================================================

model_att_1 <- '
Attention =~ q6s2 + q6s3 + q6s4 + q6s6 + q6s7 + q6s8 + q6s9 + q6s10 + q6s11 + q6s12 + q6s13 + q6s14 + q6s15 + q6s17
'

model_att_3 <- '
Att_salience =~ q6s2 + q6s3 + q6s4 + q6s6
Att_body =~ q6s7 + q6s8 + q6s9 + q6s10 + q6s11 + q6s12
Att_comp =~ q6s13 + q6s14 + q6s15 + q6s17
'

model_eval_core_22_1 <- '
Evaluation_core =~ q7s1 + q7s2 + q7s3 + q7s5 + q7s6 + q7s7 + q7s8 + q7s9 + q7s10 + q7s11 + q7s12 + q7s13 + q7s14 + q7s15 + q7s16 + q7s17 + q7s19 + q7s20 + q7s21 + q7s22 + q7s23 + q7s24
'

model_eval_core_22_4 <- '
Eval_tech =~ q7s1 + q7s2 + q7s3 + q7s5 + q7s6 + q7s7 + q7s8
Eval_meaning =~ q7s9 + q7s10 + q7s11 + q7s12 + q7s13 + q7s14 + q7s15 + q7s16 + q7s17
Eval_simulation =~ q7s19 + q7s20 + q7s21
Eval_composition =~ q7s22 + q7s23 + q7s24
'

# ============================================================
# 4. 推定
# ============================================================

fit_cfa <- function(model, items) {
  lavaan::cfa(
    model = model,
    data = df_cfa,
    ordered = items,
    estimator = "WLSMV",
    parameterization = "theta",
    std.lv = TRUE,
    missing = "listwise"
  )
}

fit_att_1 <- fit_cfa(model_att_1, attention_items_48)
fit_att_3 <- fit_cfa(model_att_3, attention_items_48)

fit_eval_core_22_1 <- fit_cfa(model_eval_core_22_1, evaluation_core_items_22)
fit_eval_core_22_4 <- fit_cfa(model_eval_core_22_4, evaluation_core_items_22)

# ============================================================
# 5. fit指標・負荷量・因子間相関
# ============================================================

get_fit <- function(fit, model_name, domain, n_factor, n_item, comparison_set) {
  fm <- lavaan::fitMeasures(fit)
  pick <- function(x, fallback = NULL) {
    if (x %in% names(fm)) return(unname(fm[x]))
    if (!is.null(fallback) && fallback %in% names(fm)) return(unname(fm[fallback]))
    NA_real_
  }
  tibble::tibble(
    domain = domain,
    comparison_set = comparison_set,
    model = model_name,
    n_factor = n_factor,
    n_item = n_item,
    converged = lavaan::lavInspect(fit, "converged"),
    chisq = pick("chisq.scaled", "chisq"),
    df = pick("df.scaled", "df"),
    pvalue = pick("pvalue.scaled", "pvalue"),
    cfi = pick("cfi.scaled", "cfi"),
    tli = pick("tli.scaled", "tli"),
    rmsea = pick("rmsea.scaled", "rmsea"),
    srmr = pick("srmr"),
    aic = pick("aic")
  )
}

fit_summary <- dplyr::bind_rows(
  get_fit(fit_att_1, "Attention_1factor", "Attention", 1, length(attention_items_48), "Attention_14items"),
  get_fit(fit_att_3, "Attention_3factor", "Attention", 3, length(attention_items_48), "Attention_14items"),
  get_fit(fit_eval_core_22_1, "Evaluation_core_22_1factor", "Evaluation", 1, length(evaluation_core_items_22), "Evaluation_core_22items"),
  get_fit(fit_eval_core_22_4, "Evaluation_core_22_4factor", "Evaluation", 4, length(evaluation_core_items_22), "Evaluation_core_22items")
)

get_loadings <- function(fit, model_name, domain, comparison_set) {
  lavaan::parameterEstimates(fit, standardized = TRUE) %>%
    dplyr::filter(op == "=~") %>%
    dplyr::transmute(domain = domain, comparison_set = comparison_set, model = model_name, factor = lhs, item = rhs, loading = est, se = se, z = z, p = pvalue, std_loading = std.all)
}

loadings <- dplyr::bind_rows(
  get_loadings(fit_att_1, "Attention_1factor", "Attention", "Attention_14items"),
  get_loadings(fit_att_3, "Attention_3factor", "Attention", "Attention_14items"),
  get_loadings(fit_eval_core_22_1, "Evaluation_core_22_1factor", "Evaluation", "Evaluation_core_22items"),
  get_loadings(fit_eval_core_22_4, "Evaluation_core_22_4factor", "Evaluation", "Evaluation_core_22items")
)

get_factor_corr <- function(fit, model_name, domain, comparison_set) {
  lv <- lavaan::lavNames(fit, "lv")
  lavaan::parameterEstimates(fit, standardized = TRUE) %>%
    dplyr::filter(op == "~~", lhs %in% lv, rhs %in% lv, lhs != rhs) %>%
    dplyr::transmute(domain = domain, comparison_set = comparison_set, model = model_name, factor1 = lhs, factor2 = rhs, r = std.all, p = pvalue)
}

factor_correlations <- dplyr::bind_rows(
  get_factor_corr(fit_att_3, "Attention_3factor", "Attention", "Attention_14items"),
  get_factor_corr(fit_eval_core_22_4, "Evaluation_core_22_4factor", "Evaluation", "Evaluation_core_22items")
)

# ============================================================
# 6. モデル比較
# 同じ項目セットのモデルだけ比較する
# ============================================================

safe_lrt <- function(fit1, fit2, domain, comparison_set, comparison) {
  out <- tryCatch(
    as.data.frame(lavaan::lavTestLRT(fit1, fit2)),
    error = function(e) data.frame(error = e$message)
  )
  out$domain <- domain
  out$comparison_set <- comparison_set
  out$comparison <- comparison
  out
}

model_comparisons <- dplyr::bind_rows(
  safe_lrt(fit_att_1, fit_att_3, "Attention", "Attention_14items", "1factor_vs_3factor"),
  safe_lrt(fit_eval_core_22_1, fit_eval_core_22_4, "Evaluation", "Evaluation_core_22items", "1factor_vs_4factor")
)

low_loading_items <- loadings %>%
  dplyr::filter(std_loading < .40) %>%
  dplyr::arrange(domain, model, std_loading)

high_factor_correlations <- factor_correlations %>%
  dplyr::filter(abs(r) >= .85) %>%
  dplyr::arrange(domain, dplyr::desc(abs(r)))

# ============================================================
# 7. Excel出力
# ============================================================

model_syntax <- tibble::tibble(
  model = c(
    "Attention_1factor",
    "Attention_3factor",
    "Evaluation_core_22_1factor",
    "Evaluation_core_22_4factor"
  ),
  syntax = c(
    model_att_1,
    model_att_3,
    model_eval_core_22_1,
    model_eval_core_22_4
  )
)

writexl::write_xlsx(
  list(
    fit_summary = fit_summary,
    model_comparisons = model_comparisons,
    loadings = loadings,
    low_loading_items = low_loading_items,
    factor_correlations = factor_correlations,
    high_factor_correlations = high_factor_correlations,
    model_syntax = model_syntax
  ),
  path = "Study2_domain_CFA_48items_revised.xlsx"
)
