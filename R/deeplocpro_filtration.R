library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

deeplocpro_results <- read.table("data/raw/detected_proteins_deeploc_prob.txt", h=T, sep ="\t")

deeplocpro <- deeplocpro_results |>
  mutate(
    across(starts_with("prob_"), as.numeric)
  ) |> 
  tidyr::separate(gene_callers_id, into = c("c_prefix", "contig_num", "gene_callers_id"), sep = "_", remove = FALSE) %>%
  select(-c(c_prefix, contig_num)) %>%
  mutate(gene_callers_id = as.character(gene_callers_id)) 


# ── 2. compute confidence metrics ────────────────────────────────────────────
prob_cols <- c(names(deeplocpro)[3:8])

deeplocpro <- deeplocpro |>
  rowwise() |>
  mutate(
    # max probability = confidence in top prediction
    prob_max    = max(c_across(all_of(prob_cols)), na.rm = TRUE),

    # second-highest probability
    prob_2nd    = sort(c_across(all_of(prob_cols)), decreasing = TRUE)[2],

    # margin = difference between top and second prediction
    # large margin = unambiguous; small margin = ambiguous between two classes
    margin      = prob_max - prob_2nd,

    # entropy = uncertainty across all 6 classes
    # low entropy = confident; high entropy = uncertain
    # H = -sum(p * log(p)), max H for 6 classes = log(6) = 1.79
    entropy     = -sum(
      ifelse(c_across(all_of(prob_cols)) > 0,
             c_across(all_of(prob_cols)) * log(c_across(all_of(prob_cols))),
             0),
      na.rm = TRUE
    ),
    entropy_norm = entropy / log(6),   # 0 = perfectly confident, 1 = uniform

    # confidence tier
    deeplocpro_confidence = case_when(
      prob_max >= 0.9  & margin >= 0.5  ~ "High",       # clear winner
      prob_max >= 0.7  & margin >= 0.3  ~ "Medium",     # probable but some uncertainty
      prob_max >= 0.5  & margin >= 0.1  ~ "Low",        # weakly predicted
      TRUE                              ~ "Ambiguous"   # two or more classes equally likely
    )
  ) |>
  ungroup() |> 
  dplyr::rename(deeplocpro_localization = Localization)

deeplocpro_filt<- deeplocpro |> filter(deeplocpro_confidence %in% c("High", "Medium"))

