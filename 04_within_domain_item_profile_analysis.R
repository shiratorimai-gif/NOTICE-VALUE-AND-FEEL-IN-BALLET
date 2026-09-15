# ============================================================
# 04. Within-domain item-profile analysis
#
# Prerequisite:
#   Run 01_efa_item_selection_and_diagnostics.R first so that
#   the cleaned analytic data frame `df_raw` is available.
#
# 0. Load packages
# ============================================================

pkgs <- c("dplyr", "tidyr", "tibble", "afex", "car", "writexl")
invisible(lapply(pkgs, function(p) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)))
invisible(lapply(pkgs, library, character.only = TRUE))
afex::afex_options(type = 3, correction_aov = "GG", es_aov = "pes")
options(contrasts = c("contr.sum", "contr.poly"))

# ============================================================
# 1. Create recent viewing-experience and lesson-experience groups
# ============================================================

group_col <- if ("x6group" %in% names(df_raw)) "x6group" else if ("6group" %in% names(df_raw)) "6group" else names(df_raw)[7]

df_prof <- df_raw %>%
  dplyr::mutate(
    row_id = dplyr::row_number(),
    group6_num = as.numeric(.data[[group_col]]),
    viewing_group = dplyr::case_when(group6_num %in% c(1, 2) ~ "high", group6_num %in% c(3, 4) ~ "moderate", group6_num %in% c(5, 6) ~ "low", TRUE ~ NA_character_),
    lesson_group = dplyr::case_when(group6_num %in% c(2, 4, 6) ~ "lesson_yes", group6_num %in% c(1, 3, 5) ~ "lesson_no", TRUE ~ NA_character_)
  )

df_prof$viewing_group <- factor(df_prof$viewing_group, levels = c("low", "moderate", "high"))
df_prof$lesson_group <- factor(df_prof$lesson_group, levels = c("lesson_no", "lesson_yes"))

# ============================================================
# 2. Define the final 48-item set
# ============================================================

attention_items_48 <- c(paste0("q6s", 2:4), paste0("q6s", 6:15), "q6s17")
evaluation_items_48 <- setdiff(paste0("q7s", 1:28), c("q7s4", "q7s18"))
emotion_items_48 <- paste0("q8s", 1:8)

domain_map_48 <- dplyr::bind_rows(
  data.frame(item = attention_items_48, domain = "Attention"),
  data.frame(item = evaluation_items_48, domain = "Evaluation"),
  data.frame(item = emotion_items_48, domain = "Emotion")
)

all_items_48 <- domain_map_48$item
setdiff(all_items_48, names(df_prof))

# ============================================================
# 3. Create long-format data with within-domain centered scores
# ============================================================

item_long_48 <- df_prof %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(all_items_48), ~ as.numeric(as.character(.)))) %>%
  dplyr::select(row_id, viewing_group, lesson_group, dplyr::all_of(all_items_48)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(all_items_48), names_to = "item", values_to = "score") %>%
  dplyr::left_join(domain_map_48, by = "item") %>%
  dplyr::group_by(row_id, domain) %>%
  dplyr::mutate(domain_mean = mean(score, na.rm = TRUE), centered_score = score - domain_mean) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(row_id = factor(row_id), item = factor(item))

# ============================================================
# 4. Run the item-profile analysis separately for each domain
# Item is the within-participant factor
# ============================================================

run_profile_aov <- function(domain_name) {
  dat <- item_long_48 %>%
    dplyr::filter(domain == domain_name) %>%
    droplevels()
  
  model <- afex::aov_ez(
    id = "row_id",
    dv = "centered_score",
    data = dat,
    within = "item",
    between = c("viewing_group", "lesson_group"),
    type = 3
  )
  
  out <- as.data.frame(model$anova_table) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::mutate(domain = domain_name)
  
  list(model = model, anova = out)
}

prof_att <- run_profile_aov("Attention")
prof_eval <- run_profile_aov("Evaluation")
prof_emo <- run_profile_aov("Emotion")

profile_anova_48 <- dplyr::bind_rows(prof_att$anova, prof_eval$anova, prof_emo$anova)

View(profile_anova_48)
profile_anova_48

# ============================================================
# 5. Extract the focal interaction terms
# ============================================================

profile_key_terms_48 <- profile_anova_48 %>%
  dplyr::filter(term %in% c(
    "viewing_group:item",
    "lesson_group:item",
    "viewing_group:lesson_group:item"
  ))

View(profile_key_terms_48)
profile_key_terms_48

# ============================================================
# 6. Item-wise 2 x 3 ANOVAs of centered scores: main effect of lesson experience
# p_Holm_domain adjusts p values separately within each domain
# ============================================================

run_item_anova_centered <- function(dat) {
  dat %>%
    dplyr::group_by(domain, item) %>%
    dplyr::group_modify(~{
      m <- lm(centered_score ~ viewing_group * lesson_group, data = .x)
      a <- car::Anova(m, type = 3) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("term")
      pcol <- grep("^Pr", names(a), value = TRUE)[1]
      fcol <- grep("^F", names(a), value = TRUE)[1]
      err_ss <- sum(resid(m)^2, na.rm = TRUE)
      a %>%
        dplyr::filter(term == "lesson_group") %>%
        dplyr::transmute(
          term = term,
          df = Df,
          F = .data[[fcol]],
          p = .data[[pcol]],
          partial_eta2 = `Sum Sq` / (`Sum Sq` + err_ss)
        )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(domain) %>%
    dplyr::mutate(p_Holm_domain = p.adjust(p, method = "holm")) %>%
    dplyr::ungroup()
}

lesson_item_tests_48 <- run_item_anova_centered(item_long_48)

sig_counts_centered_48 <- lesson_item_tests_48 %>%
  dplyr::group_by(domain) %>%
  dplyr::summarise(
    n_sig_domain_Holm = sum(p_Holm_domain < .05),
    .groups = "drop"
  )

sig_lesson_items_domain_Holm <- lesson_item_tests_48 %>%
  dplyr::filter(p_Holm_domain < .05) %>%
  dplyr::arrange(dplyr::desc(partial_eta2))

View(sig_counts_centered_48)
View(sig_lesson_items_domain_Holm)

# ============================================================
# 7. Export the results to Excel
# ============================================================

writexl::write_xlsx(
  list(
    profile_anova_48 = profile_anova_48,
    profile_key_terms_48 = profile_key_terms_48,
    lesson_item_tests_48 = lesson_item_tests_48,
    sig_counts_centered_48 = sig_counts_centered_48,
    sig_lesson_items_domain_Holm = sig_lesson_items_domain_Holm
  ),
  path = "Study2_centered_profile_rmANOVA_48items.xlsx"
)
