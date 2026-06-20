# =============================================================================
# process_NTA_results.R
# =============================================================================
# Sourced by Fig_3_S3-FCM_NTA.R. Produces:
#   SummaryData_df — particles concentration per sample (cells L-1)
#   SummaryData_df_no_dens - particles concentration per sample (cells L-1) w/o density gradient
# =============================================================================
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)   
library(ggplot2)   

# load custom ggplot theme and colour palettes
source("R/utils.R")

# -----------------------------------------------------------------------------
# 0. Parameters
# -----------------------------------------------------------------------------
PATH_DENS <- "data/raw/NTA/DensityGradient/"
PATH_NO_DENS <- "data/raw/NTA/No_densityGradient/"

# concentration factor calculation — documented step by step
CF_vivaflow  <- 60/0.5    # 60 L -> 0.5 L   = 120
CF_vivaspin  <- 0.1/0.003             # 100 mL -> 3 mL  = 33.3
CF_dens_grad    <- 0.003 / 0.001              # 3 mL -> 1 mL    = 3
CF_NTA_dilution <- 0.0002/0.0005               # 0.2 mL diluted into 0.5 mL before NTA measurement

CONC_FACTOR_DENS    <- CF_vivaflow * CF_vivaspin * CF_dens_grad* CF_NTA_dilution   # 4,800
CONC_FACTOR_NO_DENS <- CF_vivaflow                                    # 200

cat(sprintf("CONC_FACTOR_DENS:    %g\n", CONC_FACTOR_DENS))
cat(sprintf("CONC_FACTOR_NO_DENS: %g\n", CONC_FACTOR_NO_DENS))

KS_P_THRESH    <- 0.05   # KS p-value threshold for replicate similarity
MAX_REPLICATES <- 3      # maximum replicates to keep per station

# size cutoffs:
size_min_use <- 20   # lower size cutoff
size_max_use <- 250   # upper size cutoff

# -----------------------------------------------------------------------------
# 1. Station metadata
# -----------------------------------------------------------------------------
station_IDs <- read.table("data/samples_metadata.txt", header = TRUE) |>
  mutate(
    Region = factor(Region, levels = c("WEST", "GYRE", "TRAN", "UP")),
    Station_ID = factor(Station_ID, levels = c(
      "SO289_44", "SO289_43", "SO289_41", "SO289_39", "SO289_37", "SO289_34",
      "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20",
      "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9",  "SO289_6",
      "SO289_3",  "SO289_1"
    ))
  ) |>
  select(Station_name, Station_ID, Region) |>
  distinct()

# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------

# parse filename label into metadata list
# filenames: StationID_Media_Dye_Dilution_Method_FlowRate Date Stamp
parse_label <- function(basename_str, has_fraction = TRUE) {
  label <- str_split_1(basename_str, "_")
  stamp <- if (!is.na(label[6])) str_split_1(label[6], " ") else rep(NA_character_, 3)
  list(
    StationID    = label[1],
    Station_name = if (has_fraction) str_split_1(label[1], "-")[1] else label[1],
    Fraction     = if (has_fraction) str_split_1(label[1], "-")[2] else NA_character_,
    Media        = label[2],
    Dye          = label[3],
    Dilution     = label[4],
    Method       = label[5],
    Flow_rate    = stamp[1],
    Date         = stamp[2],
    Stamp        = stamp[3]
  )
}

# assign replicate labels A/B/C/D/E within each group
assign_replicates <- function(df, group_vars) {
  df |>
    group_by(across(all_of(group_vars))) |>
    mutate(Replicate = factor(
      Tech_run,
      levels  = unique(Tech_run),
      labels  = c("A", "B", "C", "D", "E", "F")[seq_along(unique(Tech_run))]
    )) |>
    ungroup()
}

# load all ParticleData CSV files from a path into one tidy dataframe
load_particle_data <- function(path, has_fraction = TRUE) {
  files <- list.files(path, pattern = "_ParticleData.csv",
                      full.names = TRUE, recursive = TRUE)
  map_dfr(files, function(f) {
    meta <- parse_label(basename(f), has_fraction)
    read.csv(f, header = TRUE, fileEncoding = "utf8") |>
      filter(Included.in.distribution. == "True") |>
      select(Particle.ID, Size.nm) |>
      mutate(across(everything(), as.numeric)) |>
      mutate(!!!meta, Tech_run = paste0(meta$Date, meta$Stamp))
  })
}

# load all Summary CSV files from a path into one tidy dataframe
load_summary_data <- function(path, has_fraction = TRUE,
                              conc_factor = 1, type_label = NULL) {
  files <- list.files(path, pattern = "_Summary.csv",
                      full.names = TRUE, recursive = FALSE)

  result <- purrr::map_dfr(files, function(f) {
    meta <- parse_label(basename(f), has_fraction)
    read.csv(f, header = TRUE, fileEncoding = "utf8", skip = 93, nrows = 1000) |>
      dplyr::select(Bin.centre..nm., Concentration..particles...ml.) |>
      dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) |>
      dplyr::mutate(
        !!!meta,
        Tech_run = paste0(meta$Date, meta$Stamp),
        Conc     = Concentration..particles...ml. / conc_factor
      )
  }) |>
    dplyr::left_join(station_IDs, by = "Station_name")

  # Assign replicates 
  result |>
    assign_replicates(c("Region", "Station_ID", "Fraction"))
}

# KS test: select most similar replicates per station (up to MAX_REPLICATES)
# handles stations with < 2 replicates gracefully
run_ks_tests <- function(df,
                         size_col = "Size.nm",
                         group_vars = "Station_ID",
                         rep_col   = "Replicate",
                         max_rep   = MAX_REPLICATES,
                         p_thresh  = KS_P_THRESH) {
  rep_levels <- c("A", "B", "C", "D", "E")
  rep_pairs  <- combn(rep_levels, 2, FUN = paste, collapse = "_", simplify = TRUE)

  groups <- df |>
    group_by(across(all_of(group_vars))) |>
    group_split(.keep = TRUE)

  map_dfr(groups, function(sub) {
    key <- sub |>
      distinct(across(all_of(group_vars)))

    avail <- intersect(rep_levels, as.character(unique(sub[[rep_col]])))

    # No or single replicate available — return first available (or NA if none)
    if (length(avail) < 2) {
      return(bind_cols(
        key,
        tibble(replicate = factor(ifelse(length(avail) == 0, NA, avail[1]),
                                  levels = rep_levels))
      ))
    }

    valid_pairs <- rep_pairs[vapply(strsplit(rep_pairs, "_"),
                                    function(rs) all(rs %in% avail),
                                    logical(1))]

    if (length(valid_pairs) == 0) {
      return(bind_cols(
        key,
        tibble(replicate = factor(avail[1], levels = rep_levels))
      ))
    }

    ks_res <- map_dfr(valid_pairs, function(p) {
      reps <- strsplit(p, "_")[[1]]
      x <- sub |> filter(.data[[rep_col]] == reps[1]) |> pull(all_of(size_col))
      y <- sub |> filter(.data[[rep_col]] == reps[2]) |> pull(all_of(size_col))
      if (length(x) < 2 || length(y) < 2) return(tibble())
      k <- ks.test(x, y)
      tibble(pair = p, D = as.numeric(k$statistic), p = as.numeric(k$p.value))
    })

    replicates <- ks_res |>
      filter(p > p_thresh) |>
      separate(pair, into = c("Rep_1", "Rep_2"), sep = "_", remove = TRUE) |>
      select(Rep_1, Rep_2) |>
      pivot_longer(everything(), values_to = "replicate") |>
      select(replicate) |>
      distinct() |>
      mutate(replicate = factor(replicate, levels = rep_levels))

    # If no pair passes threshold, fall back to first available replicate
    if (nrow(replicates) == 0) {
      return(bind_cols(
        key,
        tibble(replicate = factor(avail[1], levels = rep_levels))
      ))
    }

    bind_cols(key, replicates |> slice_head(n = max_rep))
  })
}

# assign Fraction -> Type label
label_type <- function(df) {
  df |> mutate(Type = factor(case_when(
    Fraction %in% c("DEFG","HI") ~ "BEVs",
    Fraction == "JP"             ~ "High-density",
    Fraction == "ABC"            ~ "Soluble",
    TRUE                         ~ "Total"
  ), levels = c("Soluble","BEVs","High-density","Total")))
}

# data-driven size cutoffs:
#   MIN = 2.5th percentile of particle size distribution
#   MAX = 97.5th percentile of particle size distribution
# Robust to instrument noise and sparse tails; no parametric assumptions.
compute_size_cutoffs <- function(df, size_col = "Size.nm",
                                 lower_q = 0.025,
                                 upper_q = 0.975) {
  sizes <- df[[size_col]]
  sizes <- sizes[!is.na(sizes)]
  tibble(
    size_min = quantile(sizes, lower_q, na.rm = TRUE),
    size_max = quantile(sizes, upper_q, na.rm = TRUE),
    lower_q  = lower_q,
    upper_q  = upper_q
  )
}


# Grouped version: computes cutoffs per group (e.g., per Type, Region, etc.)
compute_size_cutoffs_by <- function(df,
                                    group_vars = character(),
                                    size_col   = "Size.nm",
                                    lower_q    = 0.025,
                                    upper_q    = 0.975) {
  if (length(group_vars) == 0) {
    return(compute_size_cutoffs(df, size_col, lower_q, upper_q))
  }

  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars)), .drop = FALSE) |>
    dplyr::group_modify(~ compute_size_cutoffs(.x, size_col, lower_q, upper_q)) |>
    dplyr::ungroup()
}

# -----------------------------------------------------------------------------
# 3. Load raw particle data
# -----------------------------------------------------------------------------
raw_dens <- load_particle_data(PATH_DENS, has_fraction = TRUE) |>
  left_join(station_IDs, by = "Station_name") |>
  label_type() |>
  assign_replicates(c("Region", "Station_ID", "Fraction"))  # replicate within each Fraction

# no-density: Fraction = "Total", then label_type() assigns Type = "Total"
raw_no_dens <- load_particle_data(PATH_NO_DENS, has_fraction = FALSE) |>
  left_join(station_IDs, by = "Station_name") |>
  mutate(Fraction = "Total") |>
  label_type() |>                                            # Type = "Total"
  filter(Station_name %in% raw_dens$Station_name) |>
  droplevels() |>
  assign_replicates(c("Region", "Station_ID", "Fraction"))  # replicate within Fraction

raw_NTA <- bind_rows(raw_dens, raw_no_dens)

# -----------------------------------------------------------------------------
# 4. KS test replicate selection — per Station_ID × Fraction
# -----------------------------------------------------------------------------
# Run at Fraction level (DEFG / HI / ABC / JP / Total) — finest grouping.
# Each fraction gets its own independent replicate selection.
# Consolidation of DEFG+HI → Type="BEVs" happens AFTER replicate selection (step 6).

ks_tests <- run_ks_tests(
  raw_NTA,
  group_vars = c("Station_ID", "Fraction"),   # one KS test per fraction
  rep_col    = "Replicate"
)

# filter raw particle data to best replicates — join on Station_ID + Fraction
raw_NTA_best <- raw_NTA |>
  left_join(ks_tests,
            by           = c("Station_ID", "Fraction"),
            relationship = "many-to-many") |>
  filter(Replicate == replicate)

# -----------------------------------------------------------------------------
# 5. Size cutoffs
# -----------------------------------------------------------------------------
tidy_data_filt_type <- raw_NTA_best |>
  dplyr::filter(Size.nm >= size_min_use, Size.nm <= size_max_use)

# diagnostic plot: size distribution per Region × Fraction (before consolidation)
p_cutoff_diag_type <- raw_NTA_best |>
  ggplot(aes(x = Size.nm)) +
  geom_histogram(binwidth = 2, fill = "grey70") +
  geom_vline(xintercept = size_min_use, colour = "red",  linetype = "dashed") +
  geom_vline(xintercept = size_max_use, colour = "blue", linetype = "dashed") +
  facet_grid(Region ~ Fraction, scales = "free_y") +
  labs(x = "Size (nm)", y = "Particle count",
       title = "NTA size distribution — per Region × Fraction (best replicates)") +
  theme_EF

# -----------------------------------------------------------------------------
# 6. Summary concentration data (binned, dilution-corrected)
# -----------------------------------------------------------------------------

# ── 6a. density gradient fractions ───────────────────────────────────────────
# Load → size filter → replicate selection (on Fraction) →
# sum bins within Fraction per replicate → THEN consolidate DEFG+HI → Type

SummaryData_df_frac <- load_summary_data(
  PATH_DENS,
  has_fraction = TRUE,
  conc_factor  = CONC_FACTOR_DENS
) |>
  label_type() |>
  dplyr::filter(Bin.centre..nm. >= size_min_use,
                Bin.centre..nm. <= size_max_use) |>
  # replicate selection at Fraction level
  left_join(ks_tests,
            by           = c("Station_ID", "Fraction"),
            relationship = "many-to-many") |>
  filter(Replicate == replicate) |>
  select(-replicate) |>
  # sum bins within each Fraction × Replicate
  dplyr::group_by(Region, Station_ID, Type, Fraction, Replicate, Bin.centre..nm.) |>
  dplyr::summarise(Part.conc = sum(Conc), .groups = "drop")

# ── 6b. consolidate DEFG + HI → Type = "BEVs" ────────────────────────────────
# Each of DEFG and HI has its own best replicate (selected independently above).
# To get a combined BEVs concentration: sum Part.conc across DEFG and HI
# within the same Station × Replicate × Bin.
# Replicate label after consolidation = replicate pair used (e.g. "A+A", "A+B")

SummaryData_df_BEVs <- SummaryData_df_frac |>
  filter(Type == "BEVs") |>                                 # DEFG + HI only
  dplyr::group_by(Region, Station_ID, Type, Replicate, Bin.centre..nm.) |>
  dplyr::summarise(Part.conc = sum(Part.conc), .groups = "drop") |>
  mutate(Fraction = "DEFG+HI")                             # consolidated label

# non-BEVs fractions kept at Fraction level (ABC, JP)
SummaryData_df_other <- SummaryData_df_frac |>
  filter(Type != "BEVs")                                   # ABC, JP

SummaryData_df_dens <- bind_rows(SummaryData_df_BEVs, SummaryData_df_other) |>
  mutate(
    Fraction = factor(Fraction,
                      levels = c("DEFG+HI", "DEFG", "HI", "ABC", "JP")),
    Type     = factor(Type,
                      levels = c("Soluble","BEVs","High-density","Total"))
  )

# ── 6c. no density gradient (Total) ──────────────────────────────────────────
SummaryData_df_no_dens <- load_summary_data(
  PATH_NO_DENS,
  has_fraction = FALSE,
  conc_factor  = CONC_FACTOR_NO_DENS
) |>
  mutate(Fraction = "Total") |>
  label_type() |>
  filter(Station_name %in% unique(raw_dens$Station_name)) |>
  dplyr::filter(Bin.centre..nm. >= size_min_use,
                Bin.centre..nm. <= size_max_use) |>
  left_join(ks_tests,
            by           = c("Station_ID", "Fraction"),
            relationship = "many-to-many") |>
  droplevels() |> 
  filter(Replicate == replicate) |>
  select(-replicate) |>
  dplyr::group_by(Region, Station_ID, Type, Fraction, Replicate, Bin.centre..nm.) |>
  dplyr::summarise(Part.conc = sum(Conc), .groups = "drop")
  

# ── 6d. combine all ───────────────────────────────────────────────────────────
SummaryData_df <- bind_rows(SummaryData_df_dens, SummaryData_df_no_dens) |>
  mutate(
    Fraction   = factor(Fraction,
                        levels = c("DEFG+HI","DEFG","HI","ABC","JP","Total")),
    Type       = factor(Type,
                        levels = c("Soluble","BEVs","High-density","Total")),
    Station_ID = factor(Station_ID, levels = levels(station_IDs$Station_ID)),
    Region     = factor(Region,     levels = c("WEST","GYRE","TRAN","UP"))
  )

cat(sprintf(
  "SummaryData_df built: %d rows\n  Fractions: %s\n  Types    : %s\n",
  nrow(SummaryData_df),
  paste(levels(droplevels(SummaryData_df$Fraction)), collapse = " | "),
  paste(levels(droplevels(SummaryData_df$Type)),     collapse = " | ")
))

# -----------------------------------------------------------------------------
# 7. Save
# -----------------------------------------------------------------------------
NTA_dfs <- list(
  SummaryData_df     = SummaryData_df,        # binned conc — Fraction + Type
  tidy_data_filt_type = tidy_data_filt_type,  # particle-level, size-filtered
  ks_tests           = ks_tests,              # replicate lookup (Station×Fraction)
  raw_NTA_best       = raw_NTA_best           # best-replicate particle data
)
saveRDS(NTA_dfs, file = "data/derived/NTA_summary_data.rds")





# ── selected replicates only — use tidy_data_filt_type directly ──────────────
# tidy_data_filt_type already contains only KS-selected replicates + size filter
# avoids many-to-many join inflation from the previous approach

plot_df <- tidy_data_filt_type |>
  mutate(
    Fraction   = factor(Fraction,
                        levels = c("DEFG","HI","ABC","JP","Total")),
    Region     = factor(Region,
                        levels = c("WEST","GYRE","TRAN","UP")),
    Station_ID = factor(Station_ID,
                        levels = levels(station_IDs$Station_ID))
  )

frac_colours <- c(
  "DEFG"  = "#008585",
  "HI"    = "#74a892",
  "ABC"   = "#e5c185",
  "JP"    = "#c7522a",
  "Total" = "#6a0572"
)

frac_labels <- c(
  "DEFG"  = "DEFG (BEVs)",
  "HI"    = "HI (BEVs)",
  "ABC"   = "ABC (Soluble)",
  "JP"    = "JP (High-density)",
  "Total" = "Total"
)

stations_ordered <- levels(plot_df$Station_ID) |>
  intersect(unique(as.character(plot_df$Station_ID)))

pdf("./Figures/NTA_size_distributions_per_station.pdf",
    width = 10, height = 7)

purrr::walk(stations_ordered, function(stn) {

  df_stn <- plot_df |> filter(Station_ID == stn)
  if (nrow(df_stn) == 0) return(invisible(NULL))

  reg           <- unique(as.character(df_stn$Region))
  fracs_present <- levels(droplevels(df_stn$Fraction))

  p <- df_stn |>
    ggplot(aes(x        = Size.nm,
               colour   = Fraction,
               fill     = Fraction,
               group    = interaction(Fraction, Replicate),
               linetype = Replicate)) +
    geom_density(alpha = 0.15, linewidth = 0.7, adjust = 1.5) +
    geom_vline(xintercept = c(size_min_use, size_max_use),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    scale_colour_manual(values = frac_colours, breaks = fracs_present) +
    scale_fill_manual(  values = frac_colours, breaks = fracs_present) +
    scale_x_continuous(limits = c(size_min_use, size_max_use),
                       breaks = c(50, 100, 150, 200, 250)) +
    scale_linetype_discrete(name = "Replicate") +
    facet_wrap(~ Fraction, nrow = 1,
               labeller = as_labeller(frac_labels)) +
    labs(
      x        = "Particle size (nm)",
      y        = "Density",
      colour   = "Fraction",
      fill     = "Fraction",
      title    = sprintf("%s  |  Region: %s", stn, reg),
      subtitle = "KS-selected replicates only | line type = replicate"
    ) +
    theme_EF +
    theme(
      strip.text      = element_text(face = "bold", size = 9),
      legend.position = "bottom",
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
      plot.title      = element_text(face = "bold", size = 11),
      plot.subtitle   = element_text(size = 8, colour = "grey40")
    )

  print(p)
})

dev.off()

cat(sprintf("Saved: NTA_size_distributions_per_station.pdf — %d pages\n",
            length(stations_ordered)))

cat("process_NTA_results.R done\n")


#print session info and clean the workspace
sessionInfo()
rm(list = ls())

#unload all libraries loaded by this script
lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
