# =============================================================================
# process_FCM_results.R
# =============================================================================
# Sourced by Fig_3_S3-FCM_NTA.R. Produces:
#   counts_all — cell concentration per sample (cells L-1)
# =============================================================================

library(flowWorkspace)
library(ncdfFlow)
library(flowAI)
library(ggcyto)
library(tidyr)
library(dplyr)

# load custom ggplot theme and colour palettes
source("R/utils.R")

# -----------------------------------------------------------------------------
# 0. Parameters and metadata
# -----------------------------------------------------------------------------
PATH_FCM <- "data/raw/FCM"
set.seed(42)   # reproducible gating

meta_data <- read.table(file.path(PATH_FCM, "FCM_sample_metadata.txt"), header = TRUE) |> 
  mutate(
    Region     = factor(Region,     levels = c("WEST","GYRE","TRAN","UP")),
    Station_ID = factor(Station_ID, levels = c("SO289_44", "SO289_43", "SO289_41",  "SO289_39", "SO289_37", "SO289_34",
                                                    "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20", 
                                                    "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9", "SO289_6", 
                                                    "SO289_3", "SO289_1"))) 

SYBR_CHANNEL <- "Sybr Green 1-A"
FSC_A        <- "FSC-A"
FSC_H        <- "FSC-H"
SSC_A        <- "SSC-A"

# fixed lower bound for SYBR signal — below this = instrument noise
# set based on blank (no-stain) sample or negative control
# adjust after inspecting blank FCS file
SYBR_NOISE_THRESH <- 0.5   # on logicle-transformed scale (0–5)
SSC_NOISE_THRESH  <- 0.5   # minimum SSC for bacteria (removes instrument noise)

# flowClust settings for SYBR gate
K_CLUSTERS   <- 2      # HNA (high nucleic acid) + LNA (low nucleic acid)
Z_CUTOFF     <- 0.5    # more stringent than default 0.1 — reduces borderline events
NU_EST       <- FALSE  # use Gaussian mixture (not t-distribution)

# ── 2. load and QC (same as before) ──────────────────────────────────────────
fcsFiles <- list.files(file.path(PATH_FCM, "FCS_files"),
                       full.names = TRUE, recursive = TRUE, pattern = "\\.fcs$")
fs       <- load_cytoset_from_fcs(fcsFiles)
fs_qc    <- flow_auto_qc(fs)

# ── 3. logicle transform ──────────────────────────────────────────────────────
chnls <- c(FSC_A, FSC_H, SSC_A, SYBR_CHANNEL)
trans <- estimateLogicle(fs_qc[[1]], channels = chnls)

fs_qc_tran <- flowWorkspace::transform(fs_qc, trans)
gs          <- GatingSet(fs_qc_tran)

# ── 4. DEBRIS GATE — fixed rectangle based on SYBR noise threshold ────────────
# Key change: instead of per-sample quantile-based gate (which adapts to debris)
# use a FIXED lower bound on SYBR and SSC that excludes instrument noise.
# Bacteria: SSC > SYBR_NOISE_THRESH AND SYBR > SYBR_NOISE_THRESH
# This is more reproducible across samples.

debris_gate <- flowCore::rectangleGate(
  filterId = "nonDebris",
  "SSC-A"           = c(SSC_NOISE_THRESH, Inf),
  "Sybr Green 1-A"  = c(SYBR_NOISE_THRESH, Inf)
)

gs_pop_add(gs, debris_gate, parent = "root", name = "nonDebris")
recompute(gs)

# ── 5. SYBR GATE — per-sample flowClust K=2 (HNA + LNA) ─────────────────────
# Key change: K=2 (not K=3) — HNA and LNA are the two real bacterial clusters
# z.cutoff=0.5 (not 0.1) — excludes borderline events more aggressively
# Gate computed separately per sample — avoids pooled gate mis-assignment

sybr_gates_per_sample <- lapply(sampleNames(gs), function(sn) {
  ff <- gh_pop_get_data(gs[[sn]], "nonDebris")

  tryCatch({
    g <- openCyto::gate_flowclust_2d(
      ff,
      xChannel = SYBR_CHANNEL,
      yChannel = SSC_A,
      K        = K_CLUSTERS,
      trans    = 0,
      z.cutoff = Z_CUTOFF
    )
    g@filterId <- "Population"
    g
  }, error = function(e) {
    # fallback: simple rectangle gate if flowClust fails
    message(sprintf("flowClust failed for %s — using rectangle fallback: %s", sn, e$message))
    flowCore::rectangleGate(
      filterId = "Population",
      "Sybr Green 1-A" = c(SYBR_NOISE_THRESH + 0.3, Inf),
      "SSC-A"          = c(SSC_NOISE_THRESH,          Inf)
    )
  })
})
names(sybr_gates_per_sample) <- sampleNames(gs)

gs_pop_add(gs, sybr_gates_per_sample, parent = "nonDebris", name = "Population")
recompute(gs)

# ── 7. QC plots — one panel per sample ───────────────────────────────────────
# plot all samples to visually inspect gates — key for identifying outliers

p_debris <- ggcyto(gs,
                   aes(x = `SSC-A`, y = `Sybr Green 1-A`),
                   subset = "root") +
  geom_hex(bins = 128) +
  geom_gate("nonDebris") +
  geom_stats(gate = "nonDebris", type = c("percent","count")) +
  facet_wrap(~ name, ncol = 5) +
  labs(title = "Debris gate — fixed SYBR + SSC threshold") +
  theme_EF +
  theme(strip.text = element_text(size = 6))

p_sybr <- ggcyto(gs,
                 aes(x = `Sybr Green 1-A`, y = `SSC-A`),
                 subset = "nonDebris") +
  geom_hex(bins = 128) +
  geom_gate("Population") +
  geom_stats(gate = "Population", type = c("percent","count")) +
  facet_wrap(~ name, ncol = 5) +
  labs(title = "SYBR gate — per-sample flowClust K=2") +
  theme_EF +
  theme(strip.text = element_text(size = 6))

ggsave("./Figures/FCM_gating_debris_all.pdf",   p_debris,   width = 14, height = 10, units = "in")
ggsave("./Figures/FCM_gating_sybr_all.pdf",      p_sybr,     width = 14, height = 10, units = "in")

# ── 8. extract counts and compute concentrations ──────────────────────────────
flow_rate    <- read.csv(file.path(PATH_FCM, "flow_rates.csv"), header = TRUE)
flow_rate.lm <- lm(Volume ~ Flowrate, data = flow_rate)

meta_data <- meta_data |> 
  mutate(
    Volume     = predict(flow_rate.lm, newdata = pick(everything()))
  )

counts_all <- gs_pop_get_stats(gs, c("root","nonDebris","Population"),
                               "count") |>
  as.data.frame() |>
  pivot_wider(names_from = pop, values_from = count) |>
  mutate(
    Prop_nonDebris  = nonDebris / root,
    Prop_SYBR       = Population / nonDebris,
    Prop_SYBR_root  = Population / root
  ) |>
  separate(sample, into = c("Cruise","x","Sample_ID"), sep = "_") |>
  mutate(Sample_ID = gsub("\\.fcs$", "", Sample_ID)) |>
  select(-x) |>
  left_join(meta_data, by = "Sample_ID") |>
  mutate(
    Cell_conc = 1e6 * Dilution_factor * Population / Volume
  ) |>
  dplyr::filter(!is.na(Station_ID))

# ── 9. QC: flag samples with unusual cell counts ──────────────────────────────
# flag samples where:
#   - Prop_SYBR > 0.95 (gate too inclusive)
#   - Prop_SYBR < 0.05 (gate too exclusive or sample issue)
#   - Cell_conc > 3 SD from regional mean (outlier)

counts_all <- counts_all |>
  group_by(Region) |>
  mutate(
    log_conc     = log10(Cell_conc),
    conc_mean    = mean(log_conc, na.rm = TRUE),
    conc_sd      = sd(log_conc,   na.rm = TRUE),
    conc_z       = (log_conc - conc_mean) / conc_sd,
    QC_flag      = case_when(
      Prop_SYBR > 0.95    ~ "Gate too inclusive (Prop_SYBR > 95%)",
      Prop_SYBR < 0.05    ~ "Gate too exclusive (Prop_SYBR < 5%)",
      abs(conc_z) > 3     ~ "Outlier (>3 SD from regional mean)",
      TRUE                ~ "OK"
    )
  ) |>
  ungroup() |>
  select(-log_conc, -conc_mean, -conc_sd, -conc_z)

# report QC flags
counts_all |>
  dplyr::count(QC_flag, sort = TRUE) |>
  print()

write.table(counts_all, "data/derived/FCM_cell_counts.txt")
cat(sprintf("FCM processing done: %d samples | %d stations\n",
            nrow(counts_all), length(unique(counts_all$Station_ID))))


#print session info and clean the workspace
sessionInfo()
rm(list = ls())

#unload all libraries loaded by this script
lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
