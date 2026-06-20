library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

singalp_results <- read.table("data/raw/detected_proteins_signalp_prob.txt", h=T, sep ="\t")

signalp <- singalp_results |>
  mutate(
    gene_callers_id = gene_callers_id,
    # rename columns to clean names
    prob_OTHER      = OTHER,
    prob_SP         = SP.Sec.SPI.,
    prob_LIPO       = LIPO.Sec.SPII.,
    prob_TAT        = TAT.Tat.SPI.,
    prob_TATLIPO    = TATLIPO.Tat.SPII.,
    prob_PILIN      = PILIN.Sec.SPIII.
  ) |>
  select(gene_callers_id, Prediction, CS_Position,
         prob_OTHER, prob_SP, prob_LIPO, prob_TAT, prob_TATLIPO, prob_PILIN) |> 
  tidyr::separate(gene_callers_id, into = c("c_prefix", "contig_num", "gene_callers_id"), sep = "_", remove = FALSE) %>%
  select(-c(c_prefix, contig_num)) %>%
  mutate(gene_callers_id = as.character(gene_callers_id)) 


signalp_sig_cols <- c("prob_OTHER","prob_SP","prob_LIPO",
                      "prob_TAT","prob_TATLIPO","prob_PILIN")

signalp <- signalp |>
  rowwise() |>
  mutate(
    # probability of the winning prediction
    prob_max     = max(c_across(all_of(signalp_sig_cols)), na.rm = TRUE),

    # second highest probability
    prob_2nd     = sort(c_across(all_of(signalp_sig_cols)),
                        decreasing = TRUE)[2],

    # margin: how decisively the top class beats the runner-up
    # large margin = unambiguous; small margin = borderline call
    margin       = prob_max - prob_2nd,

    # normalised entropy: 0 = perfectly confident, 1 = uniform across all classes
    entropy_norm = {
      p <- c_across(all_of(signalp_sig_cols))
      p <- p[p > 0]
      -sum(p * log(p)) / log(length(signalp_sig_cols))
    },

    # confidence tier — informed by SignalP-6 paper thresholds
    # (Teufel et al. 2022 Nature Biotechnology)
    SignalP_confidence = case_when(
      Prediction == "OTHER"                        ~ NA_character_, # no signal predicted
      prob_max >= 0.9  & margin >= 0.5             ~ "High",
      prob_max >= 0.7  & margin >= 0.3             ~ "Medium",
      prob_max >= 0.5  & margin >= 0.1             ~ "Low",
      TRUE                                         ~ "Ambiguous"
    )
  ) |>
  ungroup() |> 
  rename(SignalP_Prediction = Prediction)

signalp_filt <- signalp |>
  filter(
    SignalP_Prediction != "OTHER",                           # remove non-signal proteins
    SignalP_confidence %in% c("High", "Medium")      # keep confident predictions
  )