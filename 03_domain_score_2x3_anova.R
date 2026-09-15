# ============================================================
# 04. Domain score comparisons: 2 x 3 ANOVA
#
# Prerequisite:
#   Run 01_efa_item_selection_and_diagnostics.R first so that
#   the cleaned analytic data frame `df_raw` is available.
# ============================================================


# ============================================================
# 1. Packages and analysis options
# ============================================================

pkgs <- c("dplyr", "tibble", "car", "emmeans", "effectsize", "writexl")
invisible(lapply(pkgs, function(p) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)))
invisible(lapply(pkgs, library, character.only = TRUE))
options(contrasts = c("contr.sum", "contr.poly"))


# ============================================================
# 2. Create recent viewing-experience and lesson-experience groups
# ============================================================

group_col <- if ("x6group" %in% names(df_raw)) "x6group" else if ("6group" %in% names(df_raw)) "6group" else names(df_raw)[7]

df_score_base <- df_raw %>%
  dplyr::mutate(
    group6_num = as.numeric(.data[[group_col]]),
    lesson_group = dplyr::case_when(
      group6_num %in% c(2, 4, 6) ~ "lesson_yes",
      group6_num %in% c(1, 3, 5) ~ "lesson_no",
      TRUE ~ NA_character_
    ),
    viewing_group = dplyr::case_when(
      group6_num %in% c(1, 2) ~ "high",
      group6_num %in% c(3, 4) ~ "moderate",
      group6_num %in% c(5, 6) ~ "low",
      TRUE ~ NA_character_
    )
  )

df_score_base$lesson_group <- factor(
  df_score_base$lesson_group,
  levels = c("lesson_no", "lesson_yes")
)
df_score_base$viewing_group <- factor(
  df_score_base$viewing_group,
  levels = c("low", "moderate", "high")
)

table(
  df_score_base$viewing_group,
  df_score_base$lesson_group,
  useNA = "ifany"
)


# ============================================================
# 3. Calculate scores for the three domains using the final 48 items
# ============================================================

attention_items_48 <- c(
  paste0("q6s", 2:4),
  paste0("q6s", 6:15),
  "q6s17"
)
evaluation_items_48 <- setdiff(
  paste0("q7s", 1:28),
  c("q7s4", "q7s18")
)
emotion_items_48 <- paste0("q8s", 1:8)

score_items_48 <- c(
  attention_items_48,
  evaluation_items_48,
  emotion_items_48
)
setdiff(score_items_48, names(df_score_base))

score_data_48 <- df_score_base %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(score_items_48),
      ~ as.numeric(as.character(.))
    )
  ) %>%
  dplyr::mutate(
    evaluation_score = rowMeans(
      dplyr::select(., dplyr::all_of(evaluation_items_48)),
      na.rm = TRUE
    ),
    attention_score = rowMeans(
      dplyr::select(., dplyr::all_of(attention_items_48)),
      na.rm = TRUE
    ),
    emotion_score = rowMeans(
      dplyr::select(., dplyr::all_of(emotion_items_48)),
      na.rm = TRUE
    )
  )

score_desc_48 <- score_data_48 %>%
  dplyr::group_by(viewing_group, lesson_group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    evaluation_m = mean(evaluation_score, na.rm = TRUE),
    evaluation_sd = sd(evaluation_score, na.rm = TRUE),
    attention_m = mean(attention_score, na.rm = TRUE),
    attention_sd = sd(attention_score, na.rm = TRUE),
    emotion_m = mean(emotion_score, na.rm = TRUE),
    emotion_sd = sd(emotion_score, na.rm = TRUE),
    .groups = "drop"
  )

View(score_desc_48)
score_desc_48


# ============================================================
# 4. Define the 2 x 3 ANOVA function
# ============================================================

run_score_anova <- function(dat, outcome) {
  f <- as.formula(
    paste0(outcome, " ~ viewing_group * lesson_group")
  )
  m <- lm(f, data = dat)

  anova_tbl <- car::Anova(m, type = 3) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::mutate(outcome = outcome)

  p_col <- grep("^Pr", names(anova_tbl), value = TRUE)[1]
  anova_tbl <- anova_tbl %>%
    dplyr::mutate(p_value = .data[[p_col]])

  eta_tbl <- effectsize::eta_squared(m, partial = TRUE) %>%
    as.data.frame() %>%
    dplyr::mutate(outcome = outcome)

  emm_viewing <- emmeans::emmeans(
    m,
    pairwise ~ viewing_group,
    adjust = "holm"
  )$contrasts %>%
    as.data.frame() %>%
    dplyr::mutate(outcome = outcome)

  emm_lesson <- emmeans::emmeans(
    m,
    pairwise ~ lesson_group,
    adjust = "holm"
  )$contrasts %>%
    as.data.frame() %>%
    dplyr::mutate(outcome = outcome)

  simple_viewing_by_lesson <- emmeans::emmeans(
    m,
    pairwise ~ viewing_group | lesson_group,
    adjust = "holm"
  )$contrasts %>%
    as.data.frame() %>%
    dplyr::mutate(outcome = outcome)

  simple_lesson_by_viewing <- emmeans::emmeans(
    m,
    pairwise ~ lesson_group | viewing_group,
    adjust = "holm"
  )$contrasts %>%
    as.data.frame() %>%
    dplyr::mutate(outcome = outcome)

  list(
    model = m,
    anova = anova_tbl,
    eta = eta_tbl,
    emm_viewing = emm_viewing,
    emm_lesson = emm_lesson,
    simple_viewing_by_lesson = simple_viewing_by_lesson,
    simple_lesson_by_viewing = simple_lesson_by_viewing
  )
}


# ============================================================
# 5. Run the analyses for Evaluation, Attention, and Emotion
# ============================================================

res_eval_48 <- run_score_anova(score_data_48, "evaluation_score")
res_att_48 <- run_score_anova(score_data_48, "attention_score")
res_emo_48 <- run_score_anova(score_data_48, "emotion_score")

effects_for_holm_adjustment <- c(
  "viewing_group",
  "lesson_group",
  "viewing_group:lesson_group"
)

anova_scores_48 <- dplyr::bind_rows(
  res_eval_48$anova,
  res_att_48$anova,
  res_emo_48$anova
) %>%
  dplyr::group_by(term) %>%
  dplyr::mutate(
    p_holm_across_domains = if (
      dplyr::first(term) %in% effects_for_holm_adjustment
    ) {
      p.adjust(p_value, method = "holm")
    } else {
      NA_real_
    }
  ) %>%
  dplyr::ungroup()
eta_scores_48 <- dplyr::bind_rows(
  res_eval_48$eta,
  res_att_48$eta,
  res_emo_48$eta
)
emm_viewing_48 <- dplyr::bind_rows(
  res_eval_48$emm_viewing,
  res_att_48$emm_viewing,
  res_emo_48$emm_viewing
)
emm_lesson_48 <- dplyr::bind_rows(
  res_eval_48$emm_lesson,
  res_att_48$emm_lesson,
  res_emo_48$emm_lesson
)
simple_viewing_by_lesson_48 <- dplyr::bind_rows(
  res_eval_48$simple_viewing_by_lesson,
  res_att_48$simple_viewing_by_lesson,
  res_emo_48$simple_viewing_by_lesson
)
simple_lesson_by_viewing_48 <- dplyr::bind_rows(
  res_eval_48$simple_lesson_by_viewing,
  res_att_48$simple_lesson_by_viewing,
  res_emo_48$simple_lesson_by_viewing
)

View(anova_scores_48)
View(eta_scores_48)
View(emm_viewing_48)
View(emm_lesson_48)
View(simple_viewing_by_lesson_48)
View(simple_lesson_by_viewing_48)


# ============================================================
# 6. Export the results to Excel
# ============================================================

writexl::write_xlsx(
  list(
    score_desc_48 = score_desc_48,
    anova_scores_48 = anova_scores_48,
    eta_scores_48 = eta_scores_48,
    emm_viewing_48 = emm_viewing_48,
    emm_lesson_48 = emm_lesson_48,
    simple_viewing_by_lesson_48 = simple_viewing_by_lesson_48,
    simple_lesson_by_viewing_48 = simple_lesson_by_viewing_48
  ),
  path = "Study2_scale_score_comparison_48items_3factor.xlsx"
)
