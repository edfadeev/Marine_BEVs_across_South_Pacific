# ── OM beta-barrel family signatures ─────────────────────────────────────────
# Used in TIER 4 to resolve SignalP SP ambiguity (beta-barrel vs periplasmic)
# and in TIER 10 when SignalP SP is absent
#
# TIER 4  SignalP SP + TCDB/IPR/Pfam/NCBIfam OM family → Outer Membrane
#         SignalP SP alone is AMBIGUOUS in Gram-negatives:
#           - Classical secreted proteins → Periplasmic
#           - Beta-barrel OM proteins     → Outer Membrane (BAM complex insertion)
#         TCDB 1.B.* / IPR / Pfam resolve this ambiguity.
#         Phobius cannot detect beta-barrels so this tier is necessary.
#
# TIER 5  SignalP SP + Phobius SP (agree, no TM, no OM annotation) → Periplasmic
#         Two tools agree on Sec export; no OM family match → stays periplasmic
#
# TIER 6  SignalP SP alone (no TM, no OM annotation) → Periplasmic
#         Classical signal peptide; Phobius silent; no OM family → periplasmic
#
# TIER 7  Phobius SP alone (no SignalP hit, no TM) → Periplasmic
#
# TIER 8  Phobius NON_CYTOPLASMIC only (no SP, no TM) → Periplasmic
#
# TIER 9  SignalP PILIN → Extracellular
#         Pilin-type Sec signal → pilus surface assembly
#
# TIER 10 TCDB / IPR / Pfam / NCBIfam / blastp OM annotation alone
#         (no SignalP SP detected — e.g. signal already cleaved in DB entry,
#          or organism uses non-canonical signal)
#         → Outer Membrane
#
# TIER 11 TCDB IM family (2.A / 3.A / 1.A etc.) without Phobius TM
#         → Cytoplasmic Membrane
#         (Phobius may miss short/atypical TM helices)
#
# TIER 12 TCDB 1.C.* pore-forming toxins → Extracellular
#
# TIER 13 GO term evidence (no sequence/structure signal from above)
#         Cellular component GO terms — used as fallback only
#         ~10–15% FP rate from inferred GO annotations
#
# TIER 14 DeepLocPro High confidence (prob_max ≥ 0.9, margin ≥ 0.5)
#         ML model; ~15–20% error rate; cannot detect beta-barrels
#
# TIER 15 DeepLocPro Medium confidence (prob_max ≥ 0.7, margin ≥ 0.3)
#         Lower reliability — flag downstream
#
# TIER 16 Unassigned — no reliable evidence
#
# Input:  protein_annotations (17,835 × 20) — already in session
#         deeplocpro_filt     — gene_callers_id in full contig_gene format
#         signalp_filt        — gene_callers_id in full contig_gene format
#         phobius_prediction_df
#         GO_prediction_df
#
# Output: localization_df saved to data/protein_localization.tsv
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)

# ── helper: check pipe-delimited field for membership in a set ──────────────
any_match_fixed <- function(x, set) {
  # x   : character vector (pipe-delimited accessions)
  # set : character vector of exact tokens to match
  sapply(str_split(replace_na(x, ""), "\\|"),
         function(tokens) any(trimws(tokens) %in% set))
}

any_match_regex <- function(x, patterns) {
  # x        : character vector
  # patterns : regex patterns (any match → TRUE)
  sapply(replace_na(x, ""),
         function(s) any(sapply(patterns,
                                function(p) grepl(p, s, ignore.case = TRUE,
                                                  perl = TRUE))))
}

# =============================================================================
# REFERENCE SETS
# =============================================================================

# ── Outer Membrane beta-barrel signatures ────────────────────────────────────
OM_InterPro <- c(
  "IPR000531",  # TonB-dependent receptor plug domain
  "IPR010916",  # TonB-dependent receptor beta-barrel
  "IPR012910",  # TonB-dependent receptor (all)
  "IPR039426",  # B12
  "IPR036942",  # copper receptor
  "IPR002368",  # OmpA family
  "IPR023614",  # Porin superfamily
  "IPR001702",  # Porin
  "IPR013793",  # Porin
  "IPR021570"   # LamB-type porin
)

OM_NCBIfam <- c(
  "TIGR01352",  # TonB-dependent receptor family
  "TIGR01776",  # TonB-dependent siderophore receptor
  "TIGR01782",  # TonB-dependent receptor
  "TIGR01783",  # TonB-dependent siderophore receptor
  "TIGR04057",  # SusC/RagA
  "TIGR04056",  # SusC
  "TIGR01779",  # BtuB — cobalamin TonB receptor
  "TIGR01785",  # Haem TonB receptor
  "TIGR01786"  # Haem TonB receptor 2
)

OM_Pfam <- c(
  "PF00593",  # TonB-dependent receptor
  "PF07715",  # TonB-dependent receptor plug domain
  "PF13715",  # TonB-dependent receptor 3
  "PF01103",  # Omp85 superfamily
  "PF02321",  # Outer membrane efflux protein
  "PF00267",  # Gram-negative porin
  "PF07980"   # SusD family
)

# blastp annotation keywords that unambiguously identify OM proteins
OM_blastp_patterns <- c(
  r"(tonb.dependent receptor)",
  r"(outer membrane receptor)",
  r"(\bporin\b)",
  r"(ompa.like)",
  r"(ompa family)",
  r"(outer membrane beta.barrel)",
  r"(\bsusc\b)",
  r"(\braga\b)"
)

# ── GO cellular component terms per compartment ───────────────────────────────
GO_OM    <- c("GO:0009279",  # cell outer membrane
              "GO:0019867",  # outer membrane
              "GO:0015288",  # porin activity
              "GO:0046930")  # pore complex

GO_IM    <- c("GO:0005886",  # plasma membrane
              "GO:0005887")  # integral component of plasma membrane

GO_Peri  <- c("GO:0042597",  # periplasmic space
              "GO:0030288")  # outer membrane-bounded periplasmic space

GO_Cyto  <- c("GO:0005737",  # cytoplasm
              "GO:0005829",  # cytosol
              "GO:0009295")  # nucleoid

GO_Extra <- c("GO:0005576"  # extracellular region
              )  


# =============================================================================
# BUILD WORKING TABLE
# =============================================================================

localization_df <- protein_annotations_df |>
  mutate(

    # =========================================================================
    # FEATURE FLAGS
    # One logical flag per evidence type — used in the cascade below
    # =========================================================================

    # ── Phobius ───────────────────────────────────────────────────────────────
    # Detects alpha-helical TM topology only.
    # Cannot detect beta-barrel OM topology (model trained on alpha-helices).

    # Alpha-helical TM helix — strongest physical evidence for Inner Membrane
    Phobius_TM      = !is.na(Phobius_acc) &
                      str_detect(Phobius_acc, "\\bTRANSMEMBRANE\\b"),

    # Full signal peptide — exclude sub-region tokens (N/H/C region only)
    Phobius_SP      = !is.na(Phobius_acc) &
                      str_detect(Phobius_acc, "\\bSIGNAL_PEPTIDE\\b") &
                     !str_detect(Phobius_acc, "SIGNAL_PEPTIDE_[NHC]_REGION"),

    # Non-cytoplasmic domain without TM helix — surface/periplasm exposed
    Phobius_nonCyto = !is.na(Phobius_acc) &
                      str_detect(Phobius_acc, "NON_CYTOPLASMIC_DOMAIN") &
                     !str_detect(Phobius_acc, "\\bTRANSMEMBRANE\\b"),

    # ── SignalP-6 — High or Medium confidence only ────────────────────────────
    # LIPO (Sec/SPII): lipoprotein signal peptide → OM lipoprotein via Cys+1 rule
    # Most reliable SignalP class; near-zero FP rate
    SP_LIPO  = !is.na(SignalP_Prediction) &
               SignalP_Prediction == "LIPO" &
               SignalP_confidence %in% c("High", "Medium"),

    # TAT / TATLIPO: twin-arginine translocase
    # Substrate folds in cytoplasm THEN exported → ends up in periplasm
    # (or OM if beta-barrel — resolved by TIER 4 if needed)
    SP_TAT   = !is.na(SignalP_Prediction) &
               SignalP_Prediction %in% c("TAT", "TATLIPO") &
               SignalP_confidence %in% c("High", "Medium"),

    # SP (Sec/SPI): classical signal peptide → AMBIGUOUS destination
    # Periplasmic secreted proteins AND OM beta-barrel proteins both use Sec
    # → destination resolved by TCDB/IPR/Pfam in TIER 4
    SP_Sec   = !is.na(SignalP_Prediction) &
               SignalP_Prediction == "SP" &
               SignalP_confidence %in% c("High", "Medium"),

    # PILIN: type IV pilin Sec signal → surface pilus
    SP_PILIN = !is.na(SignalP_Prediction) &
               SignalP_Prediction == "PILIN" &
               SignalP_confidence %in% c("High", "Medium"),

    # ── DeepLocPro ────────────────────────────────────────────────────────────
    # Lowest priority: ML model; ~15–20% error; blind to beta-barrels
    DL_high   = !is.na(deeplocpro_confidence) &
                deeplocpro_confidence == "High",
    DL_medium = !is.na(deeplocpro_confidence) &
                deeplocpro_confidence == "Medium",

    # ── TCDB family class ─────────────────────────────────────────────────────
    # 1.B.* = outer membrane beta-barrel channels (all OM)
    TCDB_is_OM    = !is.na(TCDB_acc) &
                    str_detect(TCDB_acc, r"(^1\.B\.)"),

    # 2.A / 3.A / 1.A / 4.A / 5.A / 9.A = transporters (IM)
    TCDB_is_IM    = !is.na(TCDB_acc) &
                    str_detect(TCDB_acc, r"(^(2\.|3\.|1\.A\.|4\.|5\.|9\.))"),

    # 1.C.* = pore-forming toxins (extracellular/secreted)
    TCDB_is_toxin = !is.na(TCDB_acc) &
                    str_detect(TCDB_acc, r"(^1\.C\.)"),

    # ── Functional annotation → OM beta-barrel families ───────────────────────
    # Used to RESOLVE the SignalP SP ambiguity (TIER 4)
    # and as independent OM evidence when SP absent (TIER 10)
    Annot_OM_IPR    = any_match_fixed(InterPro_acc, OM_InterPro),
    Annot_OM_NCB    = any_match_fixed(NCBIfam_acc,  OM_NCBIfam),
    Annot_OM_Pfam   = any_match_fixed(Pfam_acc,     OM_Pfam),
    Annot_OM_blastp = any_match_regex(blastp_ann,   OM_blastp_patterns),
    # combined OM annotation flag (any source)
    Annot_OM        = Annot_OM_IPR | Annot_OM_NCB | Annot_OM_Pfam |
                      Annot_OM_blastp | TCDB_is_OM,

    # ── GO cellular component flags ───────────────────────────────────────────
    GO_flag_OM    = any_match_fixed(GO_terms, GO_OM),
    GO_flag_IM    = any_match_fixed(GO_terms, GO_IM),
    GO_flag_Peri  = any_match_fixed(GO_terms, GO_Peri),
    GO_flag_Cyto  = any_match_fixed(GO_terms, GO_Cyto),
    GO_flag_Extra = any_match_fixed(GO_terms, GO_Extra),

    # =========================================================================
    # CASCADE — final localization assignment
    # =========================================================================

    Localization_final = case_when(

      # ── TIER 1: SignalP LIPO → Outer Membrane ────────────────────────────────
      # Cys+1 biochemical rule — Sec/SPII lipoprotein signal
      # Highest confidence of all SignalP classes (~1% FP rate)
      # Overrides all functional annotation and predictive tools
      SP_LIPO                                        ~ "Outer Membrane",

      # ── TIER 2: SignalP TAT → Periplasmic ────────────────────────────────────
      # Twin-arginine signal — mechanistically distinct from Sec (~2% FP rate)
      # Substrate is fully folded before export → ends up in periplasm
      # Exception: if protein also has OM annotation → still periplasmic
      # (TAT-exported beta-barrels are extremely rare)
      SP_TAT                                         ~ "Periplasmic",

      # ── TIER 3: Phobius TM → Cytoplasmic Membrane ────────────────────────────
      # Alpha-helical TM helix — strongest physical evidence for IM (~3% FP rate)
      # Applied AFTER LIPO/TAT to avoid overriding lipoproteins with TM-like regions
      # NOT applied to OM beta-barrels (Phobius cannot detect them anyway)
      Phobius_TM                                     ~ "Cytoplasmic Membrane",

      # ── TIER 4: SignalP SP + OM annotation → Outer Membrane ──────────────────
      # SignalP SP alone is AMBIGUOUS:
      #   - Periplasmic proteins use Sec (SP present, stay in periplasm)
      #   - OM beta-barrels use Sec for IM crossing THEN BAM complex for OM insertion
      # TCDB 1.B.* / IPR / Pfam / NCBIfam resolve this ambiguity
      # This tier places functional annotation AFTER Phobius TM (tier 3) because:
      #   - If Phobius detects TM → definitely IM, regardless of OM family hit
      #   - OM family annotation WITHOUT SP is lower confidence (tier 10)
      SP_Sec & Annot_OM                              ~ "Outer Membrane",

      # ── TIER 5: SignalP SP + Phobius SP agree, no OM annotation → Periplasmic ─
      # Two tools agree on Sec export; no OM family match; no TM detected
      SP_Sec & Phobius_SP & !Annot_OM                ~ "Periplasmic",

      # ── TIER 6: SignalP SP alone, no TM, no OM annotation → Periplasmic ──────
      # Classical Sec signal; no TM; not an OM beta-barrel family
      SP_Sec & !Phobius_TM & !Annot_OM               ~ "Periplasmic",

      # ── TIER 7: Phobius SP alone → Periplasmic ───────────────────────────────
      # Phobius detects signal peptide; SignalP did not call SP
      # (may be below SignalP confidence threshold)
      Phobius_SP & !Phobius_TM & !SP_Sec             ~ "Periplasmic",

      # ── TIER 8: Phobius NON_CYTOPLASMIC only → Periplasmic ───────────────────
      # Surface-exposed domain without TM helix
      Phobius_nonCyto & !SP_Sec & !SP_LIPO & !SP_TAT ~ "Periplasmic",

      # ── TIER 9: SignalP PILIN → Extracellular ────────────────────────────────
      SP_PILIN                                        ~ "Extracellular",

      # ── TIER 10: OM annotation without SignalP SP ─────────────────────────────
      # Protein identified as OM beta-barrel family (IPR/Pfam/NCBIfam/TCDB/blastp)
      # but SignalP did not detect SP — may be:
      #   - Signal already cleaved in the DB reference sequence
      #   - Non-canonical signal below SignalP threshold
      #   - Protein genuinely lacks cleavable signal (rare in Gram-negative OM)
      # Lower confidence than TIER 4 because SP corroboration is absent
      Annot_OM                                        ~ "Outer Membrane",

      # ── TIER 11: TCDB IM family without Phobius TM ────────────────────────────
      # Transporter family assignment (2.A / 3.A / 1.A) without Phobius TM
      # Phobius may miss short or atypical TM helices
      TCDB_is_IM & !Phobius_TM & !SP_Sec             ~ "Cytoplasmic Membrane",

      # ── TIER 12: TCDB 1.C.* pore-forming toxins → Extracellular ──────────────
      TCDB_is_toxin                                   ~ "Extracellular",

      # ── TIER 13: GO term evidence ─────────────────────────────────────────────
      # Cellular component GO annotations used as fallback only
      # ~10–15% error rate from computationally inferred GO annotations
      # Priority within tier: OM > Peri > Extra > IM > Cyto
      # Cytoplasm excluded if conflicting with membrane/secretion GO
      GO_flag_OM   & !GO_flag_Cyto                    ~ "Outer Membrane",
      GO_flag_Peri & !GO_flag_Cyto                    ~ "Periplasmic",
      GO_flag_Extra                                   ~ "Extracellular",
      GO_flag_IM   & !GO_flag_Cyto                    ~ "Cytoplasmic Membrane",
      GO_flag_Cyto                                    ~ "Cytoplasmic",

      # ── TIER 14: DeepLocPro High confidence ───────────────────────────────────
      # All sequence-based and annotation evidence silent above
      # prob_max ≥ 0.9 & margin ≥ 0.5
      DL_high                                         ~ deeplocpro_localization,

      # ── TIER 15: DeepLocPro Medium confidence ─────────────────────────────────
      # prob_max ≥ 0.7 & margin ≥ 0.3 — flag as lower reliability
      DL_medium                                       ~ deeplocpro_localization,

      # ── TIER 16: fallback ─────────────────────────────────────────────────────
      TRUE                                            ~ "Unassigned"
    ),

    # ── source tracking ───────────────────────────────────────────────────────
    Localization_source = case_when(
      SP_LIPO                                         ~ "SignalP_LIPO",
      SP_TAT                                          ~ "SignalP_TAT",
      Phobius_TM                                      ~ "Phobius_TM",
      SP_Sec & Annot_OM                               ~ "SignalP_SP+Annot_OM",
      SP_Sec & Phobius_SP & !Annot_OM                 ~ "SignalP_SP+Phobius_SP",
      SP_Sec & !Phobius_TM & !Annot_OM                ~ "SignalP_SP",
      Phobius_SP & !Phobius_TM & !SP_Sec              ~ "Phobius_SP",
      Phobius_nonCyto & !SP_Sec & !SP_LIPO & !SP_TAT  ~ "Phobius_nonCyto",
      SP_PILIN                                        ~ "SignalP_PILIN",
      Annot_OM                                        ~ "Annot_OM_only",
      TCDB_is_IM & !Phobius_TM & !SP_Sec              ~ "TCDB_IM",
      TCDB_is_toxin                                   ~ "TCDB_toxin",
      GO_flag_OM   & !GO_flag_Cyto                    ~ "GO_OM",
      GO_flag_Peri & !GO_flag_Cyto                    ~ "GO_Peri",
      GO_flag_Extra                                   ~ "GO_Extra",
      GO_flag_IM   & !GO_flag_Cyto                    ~ "GO_IM",
      GO_flag_Cyto                                    ~ "GO_Cyto",
      DL_high                                         ~ "DeepLocPro_High",
      DL_medium                                       ~ "DeepLocPro_Medium",
      TRUE                                            ~ "Unassigned"
    ),

    # ── reliability score 0–4 ─────────────────────────────────────────────────
    # 4 = gold: two independent biochemical evidence types agree
    # 3 = strong: single unambiguous biochemical/physical predictor
    # 2 = moderate: one reliable topology or signal predictor
    # 1 = weak: GO or DeepLocPro High only
    # 0 = unreliable: DeepLocPro Medium or Unassigned
    reliability_score = case_when(
      # gold: SP corroborates OM annotation, or two signaling methods agree
      SP_LIPO & Annot_OM                              ~ 4L,
      SP_LIPO & !is.na(deeplocpro_localization) &
        deeplocpro_localization == "Outer Membrane"   ~ 4L,
      SP_Sec  & Annot_OM & Phobius_SP                 ~ 4L,
      SP_Sec  & Annot_OM & DL_high &
        deeplocpro_localization == "Outer Membrane"   ~ 4L,
      SP_TAT  & DL_high &
        deeplocpro_localization == "Periplasmic"      ~ 4L,
      Phobius_TM & TCDB_is_IM                         ~ 4L,
      # strong: single unambiguous predictor
      SP_LIPO                                         ~ 3L,
      SP_TAT                                          ~ 3L,
      SP_Sec & Annot_OM                               ~ 3L,
      SP_Sec & Phobius_SP                             ~ 3L,
      Phobius_TM & DL_high &
        deeplocpro_localization == "Cytoplasmic Membrane" ~ 3L,
      # moderate: one reliable predictor
      Phobius_TM                                      ~ 2L,
      TCDB_is_IM                                      ~ 2L,
      SP_Sec                                          ~ 2L,
      Phobius_SP                                      ~ 2L,
      Annot_OM                                        ~ 2L,
      Phobius_nonCyto                                 ~ 2L,
      # weak: GO or DeepLocPro High only
      GO_flag_OM | GO_flag_Peri | GO_flag_IM |
        GO_flag_Extra | GO_flag_Cyto                  ~ 1L,
      DL_high                                         ~ 1L,
      TRUE                                            ~ 0L
    ),

    # ── conflict flag ─────────────────────────────────────────────────────────
    # TRUE = DeepLocPro prediction disagrees with sequence/annotation-based call
    # useful for QC — investigate high-conflict proteins
    tools_conflict = case_when(
      Localization_source %in% c(
        "SignalP_LIPO", "SignalP_TAT", "Phobius_TM",
        "SignalP_SP+Annot_OM", "SignalP_SP+Phobius_SP", "SignalP_SP",
        "Annot_OM_only", "TCDB_IM"
      ) &
        !is.na(deeplocpro_localization) &
        Localization_final != deeplocpro_localization   ~ TRUE,
      !is.na(deeplocpro_localization)                   ~ FALSE,
      TRUE                                              ~ NA
    )
  )

# =============================================================================
# SUMMARY
# =============================================================================

cat("── Proteins per localization class ─────────────────────────────────\n")
localization_df |>
  dplyr::count(Localization_final, sort = TRUE) |>
  mutate(pct = round(100 * n / sum(n), 1)) 

cat("\n── Final localization × source ──────────────────────────────────────\n")
localization_df |>
  dplyr::count(Localization_source, Localization_final) |>
  arrange(Localization_source, Localization_final)

cat("\n── Reliability score distribution ───────────────────────────────────\n")
localization_df |>
  dplyr::count(reliability_score, Localization_final) |>
  tidyr::pivot_wider(
    names_from   = reliability_score,
    values_from  = n,
    values_fill  = 0,
    names_prefix = "score_"
  ) |>
  arrange(Localization_final) 

cat("\n── TonB receptor check (SP + OM annotation) ─────────────────────────\n")
localization_df |>
  filter(Annot_OM) |>
  dplyr::count(SP_Sec, Localization_final, Localization_source) 

cat("\n── DeepLocPro conflicts with sequence/annotation evidence ───────────\n")
localization_df |>
  filter(tools_conflict == TRUE) |>
  dplyr::count(Localization_source, deeplocpro_localization,
               Localization_final, sort = TRUE) 

# =============================================================================
# SAVE
# =============================================================================

localization_df |> 
  mutate (Localization_final=case_when(reliability_score==0|Localization_final=="" ~ "Unassigned", TRUE ~ Localization_final)) |> 
  select(gene_callers_id,
         Localization_final,
         Localization_source,
         reliability_score,
         tools_conflict,
         # evidence flags (for audit trail)
         SP_LIPO, SP_TAT, SP_Sec, SP_PILIN,
         Phobius_TM, Phobius_SP, Phobius_nonCyto,
         Annot_OM, Annot_OM_IPR, Annot_OM_NCB, Annot_OM_Pfam,
         Annot_OM_blastp, TCDB_is_OM, TCDB_is_IM, TCDB_is_toxin,
         GO_flag_OM, GO_flag_Peri, GO_flag_IM, GO_flag_Cyto, GO_flag_Extra,
         deeplocpro_localization, deeplocpro_confidence,
         SignalP_Prediction, SignalP_confidence) |>
  readr::write_tsv("data/derived/protein_localization.tsv")

cat("\nSaved: data/protein_localization.tsv\n")
cat(sprintf("Total proteins : %d\n", nrow(localization_df)))
cat(sprintf("Unassigned     : %d (%.1f%%)\n",
            sum(localization_df$Localization_final == "Unassigned"),
            100 * mean(localization_df$Localization_final == "Unassigned")))
