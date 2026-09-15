# ============================================================
# 01. EFA item selection and diagnostics
#
# Expected input:
#   data/analytic_data.xlsx
#
# 0. パッケージ読み込み
# ============================================================

packages <- c( "tidyverse", "readxl", "janitor", "psych", "GPArotation", "lavaan", "semTools", "effectsize", "emmeans", "writexl" ) 
installed <- rownames(installed.packages()) 
for (p in packages) { if (!p %in% installed) { install.packages(p) } } 
library(tidyverse)
library(readxl)
library(janitor)
library(psych)
library(GPArotation)
library(lavaan)
library(semTools)
library(effectsize)
library(emmeans)
library(writexl)


# ============================================================
# 1. データ読み込み
# ============================================================

data_path <- file.path("data", "analytic_data.xlsx")

df_raw <- read_excel(data_path) %>%
  clean_names()


# ============================================================
# 2. データ確認
# ============================================================

dim(df_raw)
names(df_raw)

# ============================================================
# 3. 列名一覧
# ============================================================

colnames_table <- data.frame(no = seq_along(names(df_raw)), name = names(df_raw))
View(colnames_table)


# ============================================================
# 4. 因子分析用53項目の抽出
# ============================================================

items_53_names <- c(paste0("q6s", 1:4), "q6s5_reverse", paste0("q6s", 6:15), "q6s16_reverse", "q6s17", paste0("q7s", 1:28), paste0("q8s", 1:8))
setdiff(items_53_names, names(df_raw))

items_53_907 <- df_raw %>% dplyr::select(dplyr::all_of(items_53_names)) %>% dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(as.character(.))))

dim(items_53_907)
colnames(items_53_907)

# ============================================================
# 5. 記述統計：53項目907名
# ============================================================

desc_53_907 <- psych::describe(items_53_907)

desc_table_53_907 <- desc_53_907 %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "item") %>%
  dplyr::select(item, n, mean, sd, median, min, max, skew, kurtosis)

View(desc_table_53_907)

# ============================================================
# 6. 床効果・天井効果：53項目907名
# ============================================================

floor_ceiling_53_907 <- data.frame(
  item = names(items_53_907),
  prop_1 = sapply(items_53_907, function(x) mean(x == 1, na.rm = TRUE)),
  prop_6 = sapply(items_53_907, function(x) mean(x == 6, na.rm = TRUE))) %>%
  dplyr::mutate(
    prop_1_percent = round(prop_1 * 100, 1),
    prop_6_percent = round(prop_6 * 100, 1))

View(floor_ceiling_53_907)

floor_ceiling_flag_53_907 <- floor_ceiling_53_907 %>%
  dplyr::filter(prop_1 >= .50 | prop_6 >= .50)

View(floor_ceiling_flag_53_907)

# ============================================================
# 7. 記述統計表の統合
# ============================================================

desc_final_53_907 <- desc_table_53_907 %>%
  dplyr::left_join(floor_ceiling_53_907, by = "item") %>%
  dplyr::select(item, n, mean, sd, median, min, max, skew, kurtosis, prop_1, prop_1_percent, prop_6, prop_6_percent)

View(desc_final_53_907)

# ============================================================
# 8. 記述統計表のExcel出力
# ============================================================

writexl::write_xlsx(
  list(
    descriptive_statistics_53_907 = desc_final_53_907,
    floor_ceiling_flags_53_907 = floor_ceiling_flag_53_907),
  path = "Study2_descriptive_statistics_53items_N907.xlsx")


# ============================================================
# 9. ポリコリック相関行列・KMO・Bartlett検定：53項目907件
# ============================================================

poly_53_907 <- psych::polychoric(items_53_907)
cor_53_907 <- poly_53_907$rho

dim(cor_53_907)
round(cor_53_907[1:5, 1:5], 2)

kmo_53_907 <- psych::KMO(cor_53_907)
kmo_53_907

bartlett_53_907 <- psych::cortest.bartlett(R = cor_53_907, n = nrow(items_53_907))
bartlett_53_907

# ============================================================
# 10. 平行分析：53項目907名
# ============================================================

set.seed(1234)
parallel_53_907 <- psych::fa.parallel(cor_53_907, n.obs = nrow(items_53_907), fa = "fa", fm = "minres", n.iter = 100, main = "Parallel Analysis: 53 items, N = 907")
parallel_53_907$nfact

# ============================================================
# 11. EFA診断関数
# ============================================================

make_loading_diagnostics <- function(efa_obj, cutoff = .32, gap_cutoff = .20) {
  loading_matrix <- as.data.frame(unclass(efa_obj$loadings))
  abs_matrix <- abs(as.matrix(loading_matrix))
  factor_names <- colnames(loading_matrix)
  primary_factor <- apply(abs_matrix, 1, function(x) factor_names[which.max(x)])
  primary_loading <- apply(as.matrix(loading_matrix), 1, function(x) x[which.max(abs(x))])
  abs_primary_loading <- apply(abs_matrix, 1, max)
  second_abs_loading <- apply(abs_matrix, 1, function(x) sort(x, decreasing = TRUE)[2])
  loading_gap <- abs_primary_loading - second_abs_loading
  n_loadings_over_32 <- apply(abs_matrix >= cutoff, 1, sum)
  summary_table <- loading_matrix %>% tibble::rownames_to_column(var = "item") %>% dplyr::mutate(primary_factor = primary_factor, primary_loading = primary_loading, abs_primary_loading = abs_primary_loading, second_abs_loading = second_abs_loading, loading_gap = loading_gap, n_loadings_over_32 = n_loadings_over_32, low_loading = abs_primary_loading < cutoff, cross_loading_32 = n_loadings_over_32 >= 2, cross_loading_problem = cross_loading_32 & loading_gap < gap_cutoff)
  factor_counts <- data.frame(primary_factor = factor_names) %>% dplyr::left_join(summary_table %>% dplyr::count(primary_factor, name = "n_primary_items"), by = "primary_factor") %>% dplyr::mutate(n_primary_items = ifelse(is.na(n_primary_items), 0, n_primary_items))
  list(loadings = loading_matrix %>% tibble::rownames_to_column(var = "item"), summary = summary_table, factor_counts = factor_counts, low_loading = summary_table %>% dplyr::filter(low_loading), cross_loading_32 = summary_table %>% dplyr::filter(cross_loading_32), cross_loading_problem = summary_table %>% dplyr::filter(cross_loading_problem))
}

# ============================================================
# 12. EFA：3〜7因子解
# ============================================================

efa_3_53_907 <- psych::fa(r = cor_53_907, nfactors = 3, n.obs = nrow(items_53_907), fm = "minres", rotate = "oblimin")
efa_4_53_907 <- psych::fa(r = cor_53_907, nfactors = 4, n.obs = nrow(items_53_907), fm = "minres", rotate = "oblimin")
efa_5_53_907 <- psych::fa(r = cor_53_907, nfactors = 5, n.obs = nrow(items_53_907), fm = "minres", rotate = "oblimin")
efa_6_53_907 <- psych::fa(r = cor_53_907, nfactors = 6, n.obs = nrow(items_53_907), fm = "minres", rotate = "oblimin")
efa_7_53_907 <- psych::fa(r = cor_53_907, nfactors = 7, n.obs = nrow(items_53_907), fm = "minres", rotate = "oblimin")

# ============================================================
# 13. 診断表作成：3〜7因子解
# ============================================================

diag_3_53_907 <- make_loading_diagnostics(efa_3_53_907)
diag_4_53_907 <- make_loading_diagnostics(efa_4_53_907)
diag_5_53_907 <- make_loading_diagnostics(efa_5_53_907)
diag_6_53_907 <- make_loading_diagnostics(efa_6_53_907)
diag_7_53_907 <- make_loading_diagnostics(efa_7_53_907)

# ============================================================
# 14. 因子数ごとの適合度・診断比較
# ============================================================

efa_compare_53_907 <- data.frame(
  n_factors = 3:7,
  TLI = c(efa_3_53_907$TLI, efa_4_53_907$TLI, efa_5_53_907$TLI, efa_6_53_907$TLI, efa_7_53_907$TLI),
  RMSEA = c(efa_3_53_907$RMSEA[1], efa_4_53_907$RMSEA[1], efa_5_53_907$RMSEA[1], efa_6_53_907$RMSEA[1], efa_7_53_907$RMSEA[1]),
  RMSEA_lower = c(efa_3_53_907$RMSEA[2], efa_4_53_907$RMSEA[2], efa_5_53_907$RMSEA[2], efa_6_53_907$RMSEA[2], efa_7_53_907$RMSEA[2]),
  RMSEA_upper = c(efa_3_53_907$RMSEA[3], efa_4_53_907$RMSEA[3], efa_5_53_907$RMSEA[3], efa_6_53_907$RMSEA[3], efa_7_53_907$RMSEA[3]),
  BIC = c(efa_3_53_907$BIC, efa_4_53_907$BIC, efa_5_53_907$BIC, efa_6_53_907$BIC, efa_7_53_907$BIC),
  RMSR = c(efa_3_53_907$rms, efa_4_53_907$rms, efa_5_53_907$rms, efa_6_53_907$rms, efa_7_53_907$rms),
  n_low_loading = c(nrow(diag_3_53_907$low_loading), nrow(diag_4_53_907$low_loading), nrow(diag_5_53_907$low_loading), nrow(diag_6_53_907$low_loading), nrow(diag_7_53_907$low_loading)),
  n_cross_loading_32 = c(nrow(diag_3_53_907$cross_loading_32), nrow(diag_4_53_907$cross_loading_32), nrow(diag_5_53_907$cross_loading_32), nrow(diag_6_53_907$cross_loading_32), nrow(diag_7_53_907$cross_loading_32)),
  n_cross_loading_problem = c(nrow(diag_3_53_907$cross_loading_problem), nrow(diag_4_53_907$cross_loading_problem), nrow(diag_5_53_907$cross_loading_problem), nrow(diag_6_53_907$cross_loading_problem), nrow(diag_7_53_907$cross_loading_problem)),
  min_primary_items_per_factor = c(min(diag_3_53_907$factor_counts$n_primary_items), min(diag_4_53_907$factor_counts$n_primary_items), min(diag_5_53_907$factor_counts$n_primary_items), min(diag_6_53_907$factor_counts$n_primary_items), min(diag_7_53_907$factor_counts$n_primary_items)),
  n_empty_factors = c(sum(diag_3_53_907$factor_counts$n_primary_items == 0), sum(diag_4_53_907$factor_counts$n_primary_items == 0), sum(diag_5_53_907$factor_counts$n_primary_items == 0), sum(diag_6_53_907$factor_counts$n_primary_items == 0), sum(diag_7_53_907$factor_counts$n_primary_items == 0)))

View(efa_compare_53_907)

# ============================================================
# 15. 因子ごとの主所属項目数
# ============================================================

View(diag_3_53_907$factor_counts)
View(diag_4_53_907$factor_counts)
View(diag_5_53_907$factor_counts)
View(diag_6_53_907$factor_counts)
View(diag_7_53_907$factor_counts)

# ============================================================
# 16. EFA結果のExcel出力：53項目907名
# ============================================================

factor_cor_3_53_907 <- as.data.frame(round(efa_3_53_907$Phi, 3)) %>% tibble::rownames_to_column(var = "factor")
factor_cor_4_53_907 <- as.data.frame(round(efa_4_53_907$Phi, 3)) %>% tibble::rownames_to_column(var = "factor")
factor_cor_5_53_907 <- as.data.frame(round(efa_5_53_907$Phi, 3)) %>% tibble::rownames_to_column(var = "factor")
factor_cor_6_53_907 <- as.data.frame(round(efa_6_53_907$Phi, 3)) %>% tibble::rownames_to_column(var = "factor")
factor_cor_7_53_907 <- as.data.frame(round(efa_7_53_907$Phi, 3)) %>% tibble::rownames_to_column(var = "factor")

writexl::write_xlsx(
  list(
    desc_53_907 = desc_final_53_907,
    efa_compare_53_907 = efa_compare_53_907,
    counts_3 = diag_3_53_907$factor_counts, counts_4 = diag_4_53_907$factor_counts, counts_5 = diag_5_53_907$factor_counts, counts_6 = diag_6_53_907$factor_counts, counts_7 = diag_7_53_907$factor_counts,
    loadings_3 = diag_3_53_907$loadings, loadings_4 = diag_4_53_907$loadings, loadings_5 = diag_5_53_907$loadings, loadings_6 = diag_6_53_907$loadings, loadings_7 = diag_7_53_907$loadings,
    summary_3 = diag_3_53_907$summary, summary_4 = diag_4_53_907$summary, summary_5 = diag_5_53_907$summary, summary_6 = diag_6_53_907$summary, summary_7 = diag_7_53_907$summary,
    cross_problem_3 = diag_3_53_907$cross_loading_problem, cross_problem_4 = diag_4_53_907$cross_loading_problem, cross_problem_5 = diag_5_53_907$cross_loading_problem, cross_problem_6 = diag_6_53_907$cross_loading_problem, cross_problem_7 = diag_7_53_907$cross_loading_problem,
    factor_cor_3 = factor_cor_3_53_907, factor_cor_4 = factor_cor_4_53_907, factor_cor_5 = factor_cor_5_53_907, factor_cor_6 = factor_cor_6_53_907, factor_cor_7 = factor_cor_7_53_907
  ),
  path = "Study2_EFA_53items_N907_first_round.xlsx")

# ============================================================
# 17. 51項目版：q6s16_reverse と q7s18 を削除
# ============================================================

items_51_907 <- items_53_907 %>% dplyr::select(-q6s16_reverse, -q7s18)
dim(items_51_907)
colnames(items_51_907)

# ============================================================
# 18. 51項目版：相関行列・KMO・Bartlett
# ============================================================

poly_51_907 <- psych::polychoric(items_51_907)
cor_51_907 <- poly_51_907$rho

dim(cor_51_907)

kmo_51_907 <- psych::KMO(cor_51_907)
kmo_51_907

bartlett_51_907 <- psych::cortest.bartlett(R = cor_51_907, n = nrow(items_51_907))
bartlett_51_907

# ============================================================
# 19. 51項目版：平行分析
# ============================================================

set.seed(1234)
parallel_51_907 <- psych::fa.parallel(cor_51_907, n.obs = nrow(items_51_907), fa = "fa", fm = "minres", n.iter = 100, main = "Parallel Analysis: 51 items, N = 907")
parallel_51_907$nfact

# ============================================================
# 20. 51項目版：EFA 3〜7因子解
# ============================================================

efa_3_51_907 <- psych::fa(r = cor_51_907, nfactors = 3, n.obs = nrow(items_51_907), fm = "minres", rotate = "oblimin")
efa_4_51_907 <- psych::fa(r = cor_51_907, nfactors = 4, n.obs = nrow(items_51_907), fm = "minres", rotate = "oblimin")
efa_5_51_907 <- psych::fa(r = cor_51_907, nfactors = 5, n.obs = nrow(items_51_907), fm = "minres", rotate = "oblimin")
efa_6_51_907 <- psych::fa(r = cor_51_907, nfactors = 6, n.obs = nrow(items_51_907), fm = "minres", rotate = "oblimin")
efa_7_51_907 <- psych::fa(r = cor_51_907, nfactors = 7, n.obs = nrow(items_51_907), fm = "minres", rotate = "oblimin")

# ============================================================
# 21. 51項目版：診断表作成
# ============================================================

diag_3_51_907 <- make_loading_diagnostics(efa_3_51_907)
diag_4_51_907 <- make_loading_diagnostics(efa_4_51_907)
diag_5_51_907 <- make_loading_diagnostics(efa_5_51_907)
diag_6_51_907 <- make_loading_diagnostics(efa_6_51_907)
diag_7_51_907 <- make_loading_diagnostics(efa_7_51_907)

# ============================================================
# 22. 51項目版：因子数比較表
# ============================================================

efa_compare_51_907 <- data.frame(
  n_factors = 3:7,
  TLI = c(efa_3_51_907$TLI, efa_4_51_907$TLI, efa_5_51_907$TLI, efa_6_51_907$TLI, efa_7_51_907$TLI),
  RMSEA = c(efa_3_51_907$RMSEA[1], efa_4_51_907$RMSEA[1], efa_5_51_907$RMSEA[1], efa_6_51_907$RMSEA[1], efa_7_51_907$RMSEA[1]),
  RMSEA_lower = c(efa_3_51_907$RMSEA[2], efa_4_51_907$RMSEA[2], efa_5_51_907$RMSEA[2], efa_6_51_907$RMSEA[2], efa_7_51_907$RMSEA[2]),
  RMSEA_upper = c(efa_3_51_907$RMSEA[3], efa_4_51_907$RMSEA[3], efa_5_51_907$RMSEA[3], efa_6_51_907$RMSEA[3], efa_7_51_907$RMSEA[3]),
  BIC = c(efa_3_51_907$BIC, efa_4_51_907$BIC, efa_5_51_907$BIC, efa_6_51_907$BIC, efa_7_51_907$BIC),
  RMSR = c(efa_3_51_907$rms, efa_4_51_907$rms, efa_5_51_907$rms, efa_6_51_907$rms, efa_7_51_907$rms),
  n_low_loading = c(nrow(diag_3_51_907$low_loading), nrow(diag_4_51_907$low_loading), nrow(diag_5_51_907$low_loading), nrow(diag_6_51_907$low_loading), nrow(diag_7_51_907$low_loading)),
  n_cross_loading_32 = c(nrow(diag_3_51_907$cross_loading_32), nrow(diag_4_51_907$cross_loading_32), nrow(diag_5_51_907$cross_loading_32), nrow(diag_6_51_907$cross_loading_32), nrow(diag_7_51_907$cross_loading_32)),
  n_cross_loading_problem = c(nrow(diag_3_51_907$cross_loading_problem), nrow(diag_4_51_907$cross_loading_problem), nrow(diag_5_51_907$cross_loading_problem), nrow(diag_6_51_907$cross_loading_problem), nrow(diag_7_51_907$cross_loading_problem)),
  min_primary_items_per_factor = c(min(diag_3_51_907$factor_counts$n_primary_items), min(diag_4_51_907$factor_counts$n_primary_items), min(diag_5_51_907$factor_counts$n_primary_items), min(diag_6_51_907$factor_counts$n_primary_items), min(diag_7_51_907$factor_counts$n_primary_items)),
  n_empty_factors = c(sum(diag_3_51_907$factor_counts$n_primary_items == 0), sum(diag_4_51_907$factor_counts$n_primary_items == 0), sum(diag_5_51_907$factor_counts$n_primary_items == 0), sum(diag_6_51_907$factor_counts$n_primary_items == 0), sum(diag_7_51_907$factor_counts$n_primary_items == 0)))

View(efa_compare_51_907)

# ============================================================
# 23. 51項目版：Excel出力
# ============================================================

writexl::write_xlsx(
  list(
    efa_compare_51_907 = efa_compare_51_907,
    counts_3 = diag_3_51_907$factor_counts, counts_4 = diag_4_51_907$factor_counts, counts_5 = diag_5_51_907$factor_counts, counts_6 = diag_6_51_907$factor_counts, counts_7 = diag_7_51_907$factor_counts,
    loadings_3 = diag_3_51_907$loadings, loadings_4 = diag_4_51_907$loadings, loadings_5 = diag_5_51_907$loadings, loadings_6 = diag_6_51_907$loadings, loadings_7 = diag_7_51_907$loadings,
    summary_3 = diag_3_51_907$summary, summary_4 = diag_4_51_907$summary, summary_5 = diag_5_51_907$summary, summary_6 = diag_6_51_907$summary, summary_7 = diag_7_51_907$summary,
    cross_problem_3 = diag_3_51_907$cross_loading_problem, cross_problem_4 = diag_4_51_907$cross_loading_problem, cross_problem_5 = diag_5_51_907$cross_loading_problem, cross_problem_6 = diag_6_51_907$cross_loading_problem, cross_problem_7 = diag_7_51_907$cross_loading_problem),
  path = "Study2_EFA_51items_N907_after_removing_q6s16reverse_q7s18.xlsx")

# ============================================================
# 24. 49項目版：q6s1 と q7s4 を削除
# ============================================================

items_49_907 <- items_51_907 %>% dplyr::select(-q6s1, -q7s4)
dim(items_49_907)

# ============================================================
# 25. 49項目版：相関行列・KMO・Bartlett
# ============================================================

poly_49_907 <- psych::polychoric(items_49_907)
cor_49_907 <- poly_49_907$rho
kmo_49_907 <- psych::KMO(cor_49_907)
bartlett_49_907 <- psych::cortest.bartlett(R = cor_49_907, n = nrow(items_49_907))

dim(cor_49_907)
kmo_49_907
bartlett_49_907

# ============================================================
# 26. 49項目版：平行分析
# ============================================================

set.seed(1234)
parallel_49_907 <- psych::fa.parallel(cor_49_907, n.obs = nrow(items_49_907), fa = "fa", fm = "minres", n.iter = 100, main = "Parallel Analysis: 49 items, N = 907")
parallel_49_907$nfact

# ============================================================
# 27. 49項目版：EFA 3〜7因子解
# ============================================================

efa_3_49_907 <- psych::fa(r = cor_49_907, nfactors = 3, n.obs = nrow(items_49_907), fm = "minres", rotate = "oblimin")
efa_4_49_907 <- psych::fa(r = cor_49_907, nfactors = 4, n.obs = nrow(items_49_907), fm = "minres", rotate = "oblimin")
efa_5_49_907 <- psych::fa(r = cor_49_907, nfactors = 5, n.obs = nrow(items_49_907), fm = "minres", rotate = "oblimin")
efa_6_49_907 <- psych::fa(r = cor_49_907, nfactors = 6, n.obs = nrow(items_49_907), fm = "minres", rotate = "oblimin")
efa_7_49_907 <- psych::fa(r = cor_49_907, nfactors = 7, n.obs = nrow(items_49_907), fm = "minres", rotate = "oblimin")

# ============================================================
# 28. 49項目版：診断表
# ============================================================

diag_3_49_907 <- make_loading_diagnostics(efa_3_49_907)
diag_4_49_907 <- make_loading_diagnostics(efa_4_49_907)
diag_5_49_907 <- make_loading_diagnostics(efa_5_49_907)
diag_6_49_907 <- make_loading_diagnostics(efa_6_49_907)
diag_7_49_907 <- make_loading_diagnostics(efa_7_49_907)

# ============================================================
# 29. 49項目版：因子数比較表
# ============================================================

efa_compare_49_907 <- data.frame(
  n_factors = 3:7,
  TLI = c(efa_3_49_907$TLI, efa_4_49_907$TLI, efa_5_49_907$TLI, efa_6_49_907$TLI, efa_7_49_907$TLI),
  RMSEA = c(efa_3_49_907$RMSEA[1], efa_4_49_907$RMSEA[1], efa_5_49_907$RMSEA[1], efa_6_49_907$RMSEA[1], efa_7_49_907$RMSEA[1]),
  RMSEA_lower = c(efa_3_49_907$RMSEA[2], efa_4_49_907$RMSEA[2], efa_5_49_907$RMSEA[2], efa_6_49_907$RMSEA[2], efa_7_49_907$RMSEA[2]),
  RMSEA_upper = c(efa_3_49_907$RMSEA[3], efa_4_49_907$RMSEA[3], efa_5_49_907$RMSEA[3], efa_6_49_907$RMSEA[3], efa_7_49_907$RMSEA[3]),
  BIC = c(efa_3_49_907$BIC, efa_4_49_907$BIC, efa_5_49_907$BIC, efa_6_49_907$BIC, efa_7_49_907$BIC),
  RMSR = c(efa_3_49_907$rms, efa_4_49_907$rms, efa_5_49_907$rms, efa_6_49_907$rms, efa_7_49_907$rms),
  n_low_loading = c(nrow(diag_3_49_907$low_loading), nrow(diag_4_49_907$low_loading), nrow(diag_5_49_907$low_loading), nrow(diag_6_49_907$low_loading), nrow(diag_7_49_907$low_loading)),
  n_cross_loading_32 = c(nrow(diag_3_49_907$cross_loading_32), nrow(diag_4_49_907$cross_loading_32), nrow(diag_5_49_907$cross_loading_32), nrow(diag_6_49_907$cross_loading_32), nrow(diag_7_49_907$cross_loading_32)),
  n_cross_loading_problem = c(nrow(diag_3_49_907$cross_loading_problem), nrow(diag_4_49_907$cross_loading_problem), nrow(diag_5_49_907$cross_loading_problem), nrow(diag_6_49_907$cross_loading_problem), nrow(diag_7_49_907$cross_loading_problem)),
  min_primary_items_per_factor = c(min(diag_3_49_907$factor_counts$n_primary_items), min(diag_4_49_907$factor_counts$n_primary_items), min(diag_5_49_907$factor_counts$n_primary_items), min(diag_6_49_907$factor_counts$n_primary_items), min(diag_7_49_907$factor_counts$n_primary_items)),
  n_empty_factors = c(sum(diag_3_49_907$factor_counts$n_primary_items == 0), sum(diag_4_49_907$factor_counts$n_primary_items == 0), sum(diag_5_49_907$factor_counts$n_primary_items == 0), sum(diag_6_49_907$factor_counts$n_primary_items == 0), sum(diag_7_49_907$factor_counts$n_primary_items == 0)))

View(efa_compare_49_907)

# ============================================================
# 30. 49項目版：Excel出力
# ============================================================

writexl::write_xlsx(
  list(
    efa_compare_49_907 = efa_compare_49_907,
    counts_3 = diag_3_49_907$factor_counts, counts_4 = diag_4_49_907$factor_counts, counts_5 = diag_5_49_907$factor_counts, counts_6 = diag_6_49_907$factor_counts, counts_7 = diag_7_49_907$factor_counts,
    loadings_3 = diag_3_49_907$loadings, loadings_4 = diag_4_49_907$loadings, loadings_5 = diag_5_49_907$loadings, loadings_6 = diag_6_49_907$loadings, loadings_7 = diag_7_49_907$loadings,
    summary_3 = diag_3_49_907$summary, summary_4 = diag_4_49_907$summary, summary_5 = diag_5_49_907$summary, summary_6 = diag_6_49_907$summary, summary_7 = diag_7_49_907$summary,
    cross_problem_3 = diag_3_49_907$cross_loading_problem, cross_problem_4 = diag_4_49_907$cross_loading_problem, cross_problem_5 = diag_5_49_907$cross_loading_problem, cross_problem_6 = diag_6_49_907$cross_loading_problem, cross_problem_7 = diag_7_49_907$cross_loading_problem
  ),
  path = "Study2_EFA_49items_N907_after_removing_q6s1_q7s4.xlsx")

# ============================================================
# 31. 49項目3因子解：因子間相関・説明率
# ============================================================

factor_cor_3_49_907 <- as.data.frame(round(efa_3_49_907$Phi, 3)) %>% tibble::rownames_to_column(var = "factor")
variance_3_49_907 <- as.data.frame(efa_3_49_907$Vaccounted) %>% tibble::rownames_to_column(var = "index")

View(factor_cor_3_49_907)
View(variance_3_49_907)

round(efa_3_49_907$Phi, 2)
efa_3_49_907$Vaccounted

# ============================================================
# 32. 49項目3因子解：因子ごとの項目リスト
# ============================================================

items_by_factor_3_49 <- split(diag_3_49_907$summary$item, diag_3_49_907$summary$primary_factor)
items_by_factor_3_49

factor_item_table_3_49 <- diag_3_49_907$summary %>% dplyr::select(item, primary_factor, primary_loading, abs_primary_loading, second_abs_loading, loading_gap)
View(factor_item_table_3_49)

# ============================================================
# 33. 49項目3因子解：信頼性分析
# ============================================================

alpha_factor_1_49 <- psych::alpha(items_49_907[, items_by_factor_3_49[[1]]], check.keys = FALSE)
alpha_factor_2_49 <- psych::alpha(items_49_907[, items_by_factor_3_49[[2]]], check.keys = FALSE)
alpha_factor_3_49 <- psych::alpha(items_49_907[, items_by_factor_3_49[[3]]], check.keys = FALSE)

alpha_factor_1_49
alpha_factor_2_49
alpha_factor_3_49

# q6s5_reverse が入っている因子を確認
names(items_by_factor_3_49)[sapply(items_by_factor_3_49, function(x) "q6s5_reverse" %in% x)]

# q6s5_reverse を含む因子で check.keys = TRUE
q6s5_factor_name <- names(items_by_factor_3_49)[sapply(items_by_factor_3_49, function(x) "q6s5_reverse" %in% x)]
alpha_q6s5_check <- psych::alpha(items_49_907[, items_by_factor_3_49[[q6s5_factor_name]]], check.keys = TRUE)

alpha_q6s5_check
alpha_q6s5_check$keys

# ============================================================
# 34. 信頼性結果の整理・出力
# ============================================================

make_alpha_output <- function(alpha_obj, scale_name) {
  list(
    total = alpha_obj$total %>% as.data.frame() %>% tibble::rownames_to_column("index") %>% dplyr::mutate(scale = scale_name),
    item_stats = alpha_obj$item.stats %>% as.data.frame() %>% tibble::rownames_to_column("item") %>% dplyr::mutate(scale = scale_name),
    alpha_drop = alpha_obj$alpha.drop %>% as.data.frame() %>% tibble::rownames_to_column("item") %>% dplyr::mutate(scale = scale_name))}

alpha_out_1_49 <- make_alpha_output(alpha_factor_1_49, names(items_by_factor_3_49)[1])
alpha_out_2_49 <- make_alpha_output(alpha_factor_2_49, names(items_by_factor_3_49)[2])
alpha_out_3_49 <- make_alpha_output(alpha_factor_3_49, names(items_by_factor_3_49)[3])

writexl::write_xlsx(
  list(
    factor_cor_3_49 = factor_cor_3_49_907,
    variance_3_49 = variance_3_49_907,
    factor_item_table_3_49 = factor_item_table_3_49,
    alpha_total_1 = alpha_out_1_49$total,
    alpha_item_1 = alpha_out_1_49$item_stats,
    alpha_drop_1 = alpha_out_1_49$alpha_drop,
    alpha_total_2 = alpha_out_2_49$total,
    alpha_item_2 = alpha_out_2_49$item_stats,
    alpha_drop_2 = alpha_out_2_49$alpha_drop,
    alpha_total_3 = alpha_out_3_49$total,
    alpha_item_3 = alpha_out_3_49$item_stats,
    alpha_drop_3 = alpha_out_3_49$alpha_drop),
  path = "Study2_EFA_49items_3factor_reliability_N907.xlsx")

# ============================================================
# 35. 48項目版：q6s5_reverse を削除
# ============================================================

items_48_907 <- items_49_907 %>% dplyr::select(-q6s5_reverse)
dim(items_48_907)

# ============================================================
# 36. 48項目版：相関行列・KMO・Bartlett
# ============================================================

poly_48_907 <- psych::polychoric(items_48_907)
cor_48_907 <- poly_48_907$rho
kmo_48_907 <- psych::KMO(cor_48_907)
bartlett_48_907 <- psych::cortest.bartlett(R = cor_48_907, n = nrow(items_48_907))

dim(cor_48_907)
kmo_48_907
bartlett_48_907

# ============================================================
# 37. 48項目版：平行分析
# ============================================================

set.seed(1234)
parallel_48_907 <- psych::fa.parallel(cor_48_907, n.obs = nrow(items_48_907), fa = "fa", fm = "minres", n.iter = 100, main = "Parallel Analysis: 48 items, N = 907")
parallel_48_907$nfact

# ============================================================
# 38. 48項目版：EFA 3〜7因子解
# ============================================================

efa_3_48_907 <- psych::fa(r = cor_48_907, nfactors = 3, n.obs = nrow(items_48_907), fm = "minres", rotate = "oblimin")
efa_4_48_907 <- psych::fa(r = cor_48_907, nfactors = 4, n.obs = nrow(items_48_907), fm = "minres", rotate = "oblimin")
efa_5_48_907 <- psych::fa(r = cor_48_907, nfactors = 5, n.obs = nrow(items_48_907), fm = "minres", rotate = "oblimin")
efa_6_48_907 <- psych::fa(r = cor_48_907, nfactors = 6, n.obs = nrow(items_48_907), fm = "minres", rotate = "oblimin")
efa_7_48_907 <- psych::fa(r = cor_48_907, nfactors = 7, n.obs = nrow(items_48_907), fm = "minres", rotate = "oblimin")

# ============================================================
# 39. 48項目版：診断表
# ============================================================

diag_3_48_907 <- make_loading_diagnostics(efa_3_48_907)
diag_4_48_907 <- make_loading_diagnostics(efa_4_48_907)
diag_5_48_907 <- make_loading_diagnostics(efa_5_48_907)
diag_6_48_907 <- make_loading_diagnostics(efa_6_48_907)
diag_7_48_907 <- make_loading_diagnostics(efa_7_48_907)

# ============================================================
# 40. 48項目版：因子数比較表
# ============================================================

efa_compare_48_907 <- data.frame(
  n_factors = 3:7,
  TLI = c(efa_3_48_907$TLI, efa_4_48_907$TLI, efa_5_48_907$TLI, efa_6_48_907$TLI, efa_7_48_907$TLI),
  RMSEA = c(efa_3_48_907$RMSEA[1], efa_4_48_907$RMSEA[1], efa_5_48_907$RMSEA[1], efa_6_48_907$RMSEA[1], efa_7_48_907$RMSEA[1]),
  RMSEA_lower = c(efa_3_48_907$RMSEA[2], efa_4_48_907$RMSEA[2], efa_5_48_907$RMSEA[2], efa_6_48_907$RMSEA[2], efa_7_48_907$RMSEA[2]),
  RMSEA_upper = c(efa_3_48_907$RMSEA[3], efa_4_48_907$RMSEA[3], efa_5_48_907$RMSEA[3], efa_6_48_907$RMSEA[3], efa_7_48_907$RMSEA[3]),
  BIC = c(efa_3_48_907$BIC, efa_4_48_907$BIC, efa_5_48_907$BIC, efa_6_48_907$BIC, efa_7_48_907$BIC),
  RMSR = c(efa_3_48_907$rms, efa_4_48_907$rms, efa_5_48_907$rms, efa_6_48_907$rms, efa_7_48_907$rms),
  n_low_loading = c(nrow(diag_3_48_907$low_loading), nrow(diag_4_48_907$low_loading), nrow(diag_5_48_907$low_loading), nrow(diag_6_48_907$low_loading), nrow(diag_7_48_907$low_loading)),
  n_cross_loading_32 = c(nrow(diag_3_48_907$cross_loading_32), nrow(diag_4_48_907$cross_loading_32), nrow(diag_5_48_907$cross_loading_32), nrow(diag_6_48_907$cross_loading_32), nrow(diag_7_48_907$cross_loading_32)),
  n_cross_loading_problem = c(nrow(diag_3_48_907$cross_loading_problem), nrow(diag_4_48_907$cross_loading_problem), nrow(diag_5_48_907$cross_loading_problem), nrow(diag_6_48_907$cross_loading_problem), nrow(diag_7_48_907$cross_loading_problem)),
  min_primary_items_per_factor = c(min(diag_3_48_907$factor_counts$n_primary_items), min(diag_4_48_907$factor_counts$n_primary_items), min(diag_5_48_907$factor_counts$n_primary_items), min(diag_6_48_907$factor_counts$n_primary_items), min(diag_7_48_907$factor_counts$n_primary_items)),
  n_empty_factors = c(sum(diag_3_48_907$factor_counts$n_primary_items == 0), sum(diag_4_48_907$factor_counts$n_primary_items == 0), sum(diag_5_48_907$factor_counts$n_primary_items == 0), sum(diag_6_48_907$factor_counts$n_primary_items == 0), sum(diag_7_48_907$factor_counts$n_primary_items == 0))
)

View(efa_compare_48_907)

# ============================================================
# 41. 48項目3因子解：因子間相関・説明率・項目表
# ============================================================

factor_cor_3_48_907 <- as.data.frame(round(efa_3_48_907$Phi, 3)) %>% tibble::rownames_to_column(var = "factor")
variance_3_48_907 <- as.data.frame(efa_3_48_907$Vaccounted) %>% tibble::rownames_to_column(var = "index")
factor_item_table_3_48 <- diag_3_48_907$summary %>% dplyr::select(item, primary_factor, primary_loading, abs_primary_loading, second_abs_loading, loading_gap)

View(factor_cor_3_48_907)
View(variance_3_48_907)
View(factor_item_table_3_48)

# ============================================================
# 42. 48項目3因子解：信頼性分析
# ============================================================

items_by_factor_3_48 <- split(diag_3_48_907$summary$item, diag_3_48_907$summary$primary_factor)

alpha_factor_1_48 <- psych::alpha(items_48_907[, items_by_factor_3_48[[1]]], check.keys = FALSE)
alpha_factor_2_48 <- psych::alpha(items_48_907[, items_by_factor_3_48[[2]]], check.keys = FALSE)
alpha_factor_3_48 <- psych::alpha(items_48_907[, items_by_factor_3_48[[3]]], check.keys = FALSE)

alpha_factor_1_48
alpha_factor_2_48
alpha_factor_3_48

# ============================================================
# 43. 48項目版：Excel出力
# ============================================================

alpha_out_1_48 <- make_alpha_output(alpha_factor_1_48, names(items_by_factor_3_48)[1])
alpha_out_2_48 <- make_alpha_output(alpha_factor_2_48, names(items_by_factor_3_48)[2])
alpha_out_3_48 <- make_alpha_output(alpha_factor_3_48, names(items_by_factor_3_48)[3])

writexl::write_xlsx(
  list(
    efa_compare_48_907 = efa_compare_48_907,
    factor_cor_3_48 = factor_cor_3_48_907,
    variance_3_48 = variance_3_48_907,
    factor_item_table_3_48 = factor_item_table_3_48,
    counts_3 = diag_3_48_907$factor_counts,
    counts_4 = diag_4_48_907$factor_counts,
    counts_5 = diag_5_48_907$factor_counts,
    counts_6 = diag_6_48_907$factor_counts,
    counts_7 = diag_7_48_907$factor_counts,
    loadings_3 = diag_3_48_907$loadings,
    summary_3 = diag_3_48_907$summary,
    summary_4 = diag_4_48_907$summary,
    summary_5 = diag_5_48_907$summary,
    summary_6 = diag_6_48_907$summary,
    summary_7 = diag_7_48_907$summary,
    cross_problem_3 = diag_3_48_907$cross_loading_problem,
    cross_problem_4 = diag_4_48_907$cross_loading_problem,
    cross_problem_5 = diag_5_48_907$cross_loading_problem,
    cross_problem_6 = diag_6_48_907$cross_loading_problem,
    cross_problem_7 = diag_7_48_907$cross_loading_problem,
    alpha_total_1 = alpha_out_1_48$total,
    alpha_item_1 = alpha_out_1_48$item_stats,
    alpha_drop_1 = alpha_out_1_48$alpha_drop,
    alpha_total_2 = alpha_out_2_48$total,
    alpha_item_2 = alpha_out_2_48$item_stats,
    alpha_drop_2 = alpha_out_2_48$alpha_drop,
    alpha_total_3 = alpha_out_3_48$total,
    alpha_item_3 = alpha_out_3_48$item_stats,
    alpha_drop_3 = alpha_out_3_48$alpha_drop
  ),
  path = "Study2_EFA_48items_N907_after_removing_q6s5reverse.xlsx")

# ============================================================
# 44. McDonald's omega：48項目3因子解
# ============================================================

items_by_factor_3_48 <- split(diag_3_48_907$summary$item, diag_3_48_907$summary$primary_factor)
items_by_factor_3_48

omega_factor_1_48 <- psych::omega(items_48_907[, items_by_factor_3_48[[1]]], nfactors = 1, fm = "minres", poly = TRUE, plot = FALSE)
omega_factor_2_48 <- psych::omega(items_48_907[, items_by_factor_3_48[[2]]], nfactors = 1, fm = "minres", poly = TRUE, plot = FALSE)
omega_factor_3_48 <- psych::omega(items_48_907[, items_by_factor_3_48[[3]]], nfactors = 1, fm = "minres", poly = TRUE, plot = FALSE)

omega_factor_1_48$omega.tot
omega_factor_2_48$omega.tot
omega_factor_3_48$omega.tot

omega_summary_48 <- data.frame(
  factor = names(items_by_factor_3_48),
  n_items = sapply(items_by_factor_3_48, length),
  omega_total = c(omega_factor_1_48$omega.tot, omega_factor_2_48$omega.tot, omega_factor_3_48$omega.tot))

View(omega_summary_48)

writexl::write_xlsx(
  list(
    omega_summary_48 = omega_summary_48
  ),
  path = "Study2_omega_48items_3factor_N907.xlsx")

# ============================================================
# 44-2. ω算出用データの確認
# ============================================================

omega_check <- function(dat) {
  r <- psych::polychoric(dat)$rho
  data.frame(n_items = ncol(dat), min_eigen = min(eigen(r)$values), max_eigen = max(eigen(r)$values))
}

omega_check_1 <- omega_check(items_48_907[, items_by_factor_3_48[[1]]])
omega_check_2 <- omega_check(items_48_907[, items_by_factor_3_48[[2]]])
omega_check_3 <- omega_check(items_48_907[, items_by_factor_3_48[[3]]])

omega_check_1
omega_check_2
omega_check_3

# ============================================================
# 追加：KMO・Bartlett・平行分析結果をExcelへ保存
# ============================================================

make_diag_row <- function(stage, n_items, dat, kmo, bartlett, parallel) {
  tibble::tibble(
    stage = stage,
    n = nrow(dat),
    n_items = n_items,
    KMO_overall = unname(kmo$MSA),
    Bartlett_chisq = unname(bartlett$chisq),
    Bartlett_df = unname(bartlett$df),
    Bartlett_p = unname(bartlett$p.value),
    parallel_nfact = unname(parallel$nfact)
  )
}

make_item_msa <- function(stage, kmo) {
  tibble::tibble(
    stage = stage,
    item = names(kmo$MSAi),
    MSA = unname(kmo$MSAi)
  )
}

efa_diagnostics_summary <- dplyr::bind_rows(
  make_diag_row("53_items", 53, items_53_907, kmo_53_907, bartlett_53_907, parallel_53_907),
  make_diag_row("51_items", 51, items_51_907, kmo_51_907, bartlett_51_907, parallel_51_907),
  make_diag_row("49_items", 49, items_49_907, kmo_49_907, bartlett_49_907, parallel_49_907),
  make_diag_row("48_items", 48, items_48_907, kmo_48_907, bartlett_48_907, parallel_48_907)
)

efa_item_MSA <- dplyr::bind_rows(
  make_item_msa("53_items", kmo_53_907),
  make_item_msa("51_items", kmo_51_907),
  make_item_msa("49_items", kmo_49_907),
  make_item_msa("48_items", kmo_48_907)
)

writexl::write_xlsx(
  list(
    diagnostics_summary = efa_diagnostics_summary,
    item_MSA = efa_item_MSA,
    KMO_53_items = make_item_msa("53_items", kmo_53_907),
    polychoric_53 = as.data.frame(cor_53_907),
    polychoric_48 = as.data.frame(cor_48_907)
  ),
  path = "Study2_EFA_KMO_Bartlett_parallel_53to48.xlsx"
)
