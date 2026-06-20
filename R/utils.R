
require(ggplot2)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9",
                "#009E73", "#F0E442", "#0072B2", 
                "#D55E00", "#CC79A7")

# 21 colours
tol21rainbow <- c("#771155", "#AA4488","#CC99BB","#114477",
                  "#4477AA","#117744","#117777","#88CCAA",
                  "#77CCCC","#00ffff","#44AA77","#44AAAA",
                  "#777711","#AAAA44","#DDDD77","#774411",
                  "#AA7744","#DDAA77","#771122","#AA4455","#DD7788")

# theme
theme_EF <- theme_bw() +
  theme(panel.grid.major    = element_blank(),
        panel.grid.minor    = element_blank(),
        strip.background    = element_blank(),
        title               = element_text(size = 25, face = "bold"),
        legend.title        = element_text(size = 20),
        legend.text         = element_text(size = 18),
        strip.text          = element_text(size = 22, face = "bold"),
        axis.title          = element_text(size = 20, face = "bold"),
        axis.text           = element_text(size = 18, color = "black"),
        axis.line           = element_line(colour = "black", linewidth = 0.8),
        strip.placement     = "outside",
        strip.background.x  = element_blank(),
        panel.spacing.x     = unit(5, "pt"))

# ── helpers ───────────────────────────────────────────────────────────────────
# calculate NSAF (normalized spectral abundance factor) from raw spectral counts and protein lengths
add_nsaf_protein_abund <- function(protein_abund,
                                   protein_metadata,
                                   id_col = "gene_callers_id",
                                   length_col = "Length") {
  rn <- rownames(protein_abund)
  cn <- colnames(protein_abund)
  if (is.null(rn)) stop("protein_abund must have rownames (gene_callers_id)")
  mat <- as.matrix(protein_abund)
  meta_ids <- protein_metadata[[id_col]]
  prot_len <- as.numeric(protein_metadata[[length_col]][match(rn, meta_ids)])
  missing_len <- is.na(prot_len) | prot_len == 0
  if (any(missing_len)) {
    warning(sprintf("%d proteins have missing or zero length; their NSAF will be NA", sum(missing_len)))
    prot_len[missing_len] <- NA_real_
  }
  mat_per_len <- sweep(mat, 1, prot_len, "/")
  col_sums <- colSums(mat_per_len, na.rm = TRUE)
  col_sums[col_sums == 0] <- NA_real_
  mat_nsaf <- sweep(mat_per_len, 2, col_sums, "/")
  if (any(missing_len)) mat_nsaf[missing_len, ] <- NA_real_
  dimnames(mat_nsaf) <- list(rn, cn)
  if (is.data.frame(protein_abund)) {
    out <- as.data.frame(mat_nsaf)
    rownames(out) <- rn
    names(out) <- cn
    out
  } else {
    mat_nsaf
  }
}


# rename "Abundance F# Sample" -> Sample_ID (from your snippet)
rename_abundance_cols <- function(df, samples_df,
                                  prefix = "^Abundance\\s+",
                                  suffix = "\\s+Sample$") {
  ab_raw  <- grep(paste0(prefix, "F\\d+", suffix), names(df), value = TRUE)
  if (!length(ab_raw)) return(list(df = df, sample_cols = character()))
  file_id <- sub(paste0(prefix, "(F\\d+)", suffix), "\\1", ab_raw)
  missing <- setdiff(file_id, samples_df$File.ID)
  if (length(missing)) warning("File.IDs not in samples_df: ", paste(missing, collapse = ", "))
  map   <- setNames(samples_df$Sample_ID, samples_df$File.ID)
  new   <- unname(map[file_id])
  ok    <- !is.na(new)
  df2   <- df |> dplyr::rename(!!!setNames(ab_raw[ok], new[ok]))
  list(df = df2, sample_cols = new[ok])
}

run_source_if_missing <- function(target, script,
                                   local    = TRUE,
                                   chdir    = FALSE,
                                   encoding = getOption("encoding")) {
  if (!file.exists(target)) {
    message("Results not found: ", target, " — sourcing ", script)
    source(script, local = local, chdir = chdir, encoding = encoding)
  } else {
    message("Results exists: ", target, " — skipping ", script)
  }
  invisible(file.exists(target))
}


# standard error
se <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x) / length(x))
}

# =============================================================================
# resolve_function()
# Assigns a human-readable Function string to each protein row by selecting
# the first available annotation source in priority order:
#   1. NCBIfam_ann    — curated HMM family; most specific
#   2. TC_description — TCDB family description; best for transporters
#   3. COG20_ann      — COG curated function; underscores + tags cleaned
#   4. InterPro_ann   — domain-level annotation
#   5. Pfam_ann       — domain family fallback
#   6. blastp_ann     — sequence similarity (cleaned of organism tags)
#   7. "hypothetical protein"
#
# Requires columns: NCBIfam_ann, TC_description, COG20_ann (optional — if
# absent the step is silently skipped), InterPro_ann, Pfam_ann, blastp_ann
# =============================================================================

# ── internal helper: clean COG20_ann strings ──────────────────────────────────
# COG20_ann format: "Function_name_with_underscores_(GeneAbbr)_(PDB:XXXX)_(PUBMED:NNNNN)"
# → replace underscores with spaces
# → strip (PDB:...) and (PUBMED:...) tags
# → strip trailing single gene-name token in parentheses e.g. "(TufA)"
.clean_cog20 <- function(x) {
  x |>
    gsub("_", " ", x = _) |>
    gsub("\\s*\\(PDB:[^)]+\\)",    "", x = _) |>
    gsub("\\s*\\(PUBMED:[^)]+\\)", "", x = _) |>
    # trailing short gene abbreviation e.g. "(TufA)", "(HisJ)", "(RpsJ)"
    gsub("\\s*\\([A-Za-z][A-Za-z0-9]{1,8}\\)\\s*$", "", x = _) |>
    trimws()
}

resolve_function <- function(df) {

  has_value <- function(x) !is.na(x) & trimws(x) != "" & trimws(x) != "-"

  # COG20_ann is optional — present only after merge with cog20_df
  has_cog20 <- "COG20_ann" %in% names(df)

  df |>
    mutate(
      # clean COG20_ann if column exists
      COG20_ann_clean = if (has_cog20)
        dplyr::if_else(has_value(COG20_ann), .clean_cog20(COG20_ann), NA_character_)
      else
        NA_character_,
      
        Function_raw = case_when(
        # 1. NCBIfam — most specific curated HMM
        has_value(NCBIfam_ann)    ~ NCBIfam_ann,

        # 2. TCDB — best for transporters/channels
        has_value(TC_description) ~ TC_description,

        # 3. COG20 — curated function, cleaned
        !is.na(COG20_ann_clean)   ~ COG20_ann_clean,

        # 4. InterPro — domain-level, reliable
        has_value(InterPro_ann)   ~ InterPro_ann,

        # 5. Pfam — domain family fallback
        has_value(Pfam_ann)       ~ Pfam_ann,

        # 6. blastp — least specific; clean organism tags
        has_value(blastp_ann)     ~
          gsub("MULTISPECIES:\\s*", "",
               gsub("\\s*\\[.*?\\]", "", blastp_ann)),

        # 7. nothing available
        TRUE                      ~ "hypothetical protein"
      ),

      # keep only first pipe-delimited domain; guard against leading "|"
      Function = {
        f <- trimws(gsub("\\|.*", "", Function_raw))
        dplyr::if_else(f == "" | is.na(f), Function_raw, f)
      },
      Function = trimws(Function),
      Function = dplyr::if_else(
        is.na(Function) | Function == "" | Function == "-",
        "hypothetical protein",
        Function
      )
    ) |>
    select(-Function_raw, -COG20_ann_clean)
}

# =============================================================================
# categorize_single_loc()
# Like categorize_single() but uses BOTH Function and Localization_final.
#
# x        : character vector of Function strings (lowercased internally)
# loc      : character vector of Localization_final values (same length as x)
#            e.g. "Outer Membrane" | "Periplasmic" | "Cytoplasmic Membrane" |
#                 "Cytoplasmic" | "Extracellular" | "Unassigned"
# patterns : named list — each element is either:
#            (a) a plain regex string — matched against Function only
#            (b) a list(func = "regex", loc = "regex") — matched against
#                both Function AND Localization_final (both must match)
#
# Returns  : character vector of category labels; "Other" if no match
#
# Usage:
#   mutate(Category = categorize_single_loc(Function, Localization_final,
#                                           patterns_cyano_loc))
# =============================================================================
categorize_single_loc <- function(x, loc, patterns) {
  x_clean   <- tolower(trimws(
    ifelse(is.na(x) | trimws(x) == "" | trimws(x) == "-",
           "hypothetical protein", x)
  ))
  loc_clean <- tolower(trimws(replace_na(loc, "unassigned")))

  purrr::map2_chr(x_clean, loc_clean, function(s, l) {
    for (nm in names(patterns)) {
      pat <- patterns[[nm]]
      if (is.list(pat)) {
        # both Function AND Localization must match
        if (str_detect(s, pat$func) &&
            str_detect(l, pat$loc))  return(nm)
      } else {
        # Function only (backward compatible with patterns_cyano)
        if (str_detect(s, pat))      return(nm)
      }
    }
    "Other"
  })
}

patterns_bacteria_loc <- list(

  # ── Outer membrane porins / receptors ────────────────────────────────────
  "Iron uptake porins" = list(
    func = "(?i)(iron uptake porin|\\bsoma\\b|soma porin)",
    loc  = "outer membrane|unassigned"
  ),

  "Carbohydrate-selective porins" = list(
    func = "(?i)(carbohydrate[[:punct:][:space:]_-]*selective porin|carbohydrate porin|\\boprb\\b)",
    loc  = "outer membrane|unassigned"
  ),

  # Generic porins and porin-domain families (TonB-dependent receptors excluded)
  "Outer membrane porins (other)" = list(
    func = paste0(
      "(?i)(",
      "outer membrane (porin|protein insertion porin)|",
      "porin domain, gram[- ][[:alpha:]]+ type|",
      "the general bacterial porin \\(gbp\\) family|",
      "the cyanobacterial porin \\(cbp\\) family|",
      "the putative bacterial porin \\(pbp\\) family|",
      "the rhodobacter porca porin \\(rpp\\) family|",
      "the treponema porin major surface protein \\(tp[- ]msp\\) family|",
      "the omp[a]-ompf porin \\(oop\\) family|",
      "\\bomp[a-z]\\b|",
      "putative beta[- ]barrel porin[- ]2|\\bbbp2\\b|ompl[- ]like|",
      "phosphate[- ]selective porin|\\bphoe\\b|porin, alpha proteobacteria type|",
      "\\bporin\\b",
      ")"
    ),
    loc  = "outer membrane|periplasmic|unassigned"
  ),

  # ── TonB-dependent receptors and related systems ─────────────────────────
  "TonB-dependent receptor (TBDR)" = list(
    func = paste0(
      "(?i)(",
      "tonb[[:punct:][:space:]_-]*dependent[[:punct:][:space:]_-]*(receptor|outer membrane receptor|siderophore receptor|heme|hemoglobin|transferrin|lactoferrin)|",
      "outer membrane receptor (protein, )?fe|",
      "outer membrane receptor for (fe3\\+[^,]*|fe3\\+-dicitrate|monomeric catechols|ferrienterochelin|colicins|ferric coprogen|ferric[- ]rhodotorulic acid)|",
      "outer membrane cobalamin receptor protein btu[b]|\\bbtu[b]\\b|",
      "\\bcir[a]\\b|\\bfep[a]\\b|\\bfec[a]\\b|",
      "\\blbtu\\b|lbt[u]\\s+family|",
      "the outer membrane receptor \\(omr\\) family",
      ")"
    ),
    loc  = "outer membrane|periplasmic|unassigned"
  ),

  "SusC (polysaccharide utilization loci)" = list(
    func = "(?i)(\\bsusc\\b|\\braga\\b|susc[[:punct:][:space:]_-]*raga|raga[[:punct:][:space:]_-]*susc|susc[[:punct:][:space:]_-]*(family|subfamily)|raga[[:punct:][:space:]_-]*(family|subfamily)|tonb[[:punct:][:space:]_-]*linked outer membrane protein|tonb[[:punct:][:space:]_-]*dependent outer membrane receptor,? signature region)",
    loc  = "outer membrane|periplasmic|unassigned"
  ),

  "SusD (polysaccharide utilization loci)" = list(
    func = "(?i)(\\bsusd\\b|\\bragb\\b|ragb[[:punct:][:space:]_-]*susd|susd[[:punct:][:space:]_-]*like|ragb[[:punct:][:space:]_-]*like|ragb/susd( domain)?|susd[- ]like 2|ragb/susd domain)",
    loc  = "outer membrane|periplasmic|unassigned"
  ),

  "TonB system (TonB/ExbB/ExbD/Tol)" = list(
    func = paste0(
      "(?i)(",
      "(the )?tonb[- ]exbb[- ]exbd/tola[- ]tolq[- ]tolr.*(family|energizer|stabilizer)|",
      "the h\\+[- ]*or na\\+[- ]*translocating bacterial flagellar motor/exbbd outer membrane transport energizer \\(mot/exb\\) superfamily|",
      "biopolymer transport protein exbb/tolq|",
      "\\bexbb\\b|\\bexbd\\b|\\btol[aqr]\\b|",
      "\\btonb\\b(?![[:space:]-]*dependent[[:space:]-]*receptor)",
      ")"
    ),
    loc  = "outer membrane|periplasmic|cytoplasmic membrane|unassigned"
  ),

  # ── Outer membrane assembly ───────────────────────────────────────────────
  "BAM complex / OM assembly" = list(
    func = "(?i)(\\bbam[a-e]\\b|bam complex|outer membrane protein assembly factor|outer membrane protein assembly factor bama)",
    loc  = "outer membrane|periplasmic"
  ),

  # ── Photosystems and light harvesting ────────────────────────────────────
  "Photosystem I (PSI)" = list(
    func = "(?i)(photosystem i( core)? protein|\\bpsa[ablcdefgiklm]\\b|ycf4|psa[- ]linked|psa[a-z])",
    loc  = "thylakoid membrane|cytoplasmic membrane|unassigned"
  ),

  "Photosystem II (PSII)" = list(
    func = "(?i)(photosystem ii|psii|\\bpsb[a-z]\\b|d1 protein|d2 protein|cp43|cp47|psb[oqu]27?)",
    loc  = "thylakoid membrane|cytoplasmic membrane|unassigned"
  ),

  "Phycobilisome and phycobiliproteins" = list(
    func = "(?i)(phycobilisome|phycocyanin|phycoerythrin|allophycocyanin|\\bcpc[a-z]\\b|\\bcpe[a-z]\\b|\\bapc[a-z]\\b|apce[- ]cpcd fusion|orange carotenoid[- ]binding protein)",
    loc  = "thylakoid membrane|periplasmic|unassigned"
  ),

  "Anoxygenic reaction centers (PRC/BRC)" = list(
    func = "(?i)(photosynthetic reaction center|prc family|bacterial reaction center|puf[lm]|puh[a-z])",
    loc  = "cytoplasmic membrane|unassigned"
  ),

  "Rhodopsin" =
    "(?i)(rhodopsin|ion[[:punct:][:space:]_-]*translocating microbial rhodopsin|bacteriorhodopsin)",

  # ── Solute binding proteins — BROAD, order‑agnostic ──────────────────────
  "Solute binding proteins (SBP)" = list(
    func = paste0(
      "(?i)(",
      "(solute|substrate)[[:punct:][:space:]_-]*binding protein|",
      "periplasmic.*binding protein|",
      "\\bznua\\b|znua[[:punct:][:space:]_-]*like|\\bpsts\\b|",
      "spermidine|putrescine|",
      "phosphate[[:punct:][:space:]_-]*binding|iron[[:punct:][:space:]_-]*binding|molybdate|",
      "amino[[:punct:][:space:]_-]*acid[[:punct:][:space:]_-]*binding|",
      "oligopeptide[[:punct:][:space:]_-]*binding|",
      "bacterial extracellular solute[[:punct:][:space:]_-]*binding|",
      "ssua/thi5[[:punct:][:space:]_-]*like|",
      "(?=.*\\b(periplasmic|extracytoplasmic)\\b)(?=.*\\b(component|receptor)\\b)(?=.*\\b(abc|abc[[:punct:][:space:]_-]*type|trap|ttt|tctc|tricarboxylate|transport)\\b).+|",
      "taxi family|",
      "lipid[[:punct:][:space:]_-]*binding transport protein|\\btim44\\b|",
      "polyisoprenoid[- ]binding periplasmic protein|\\bycei\\b",
      ")"
    ),
    loc  = ".*"
  ),

  # ── Membrane transporters (split by energy requirement) ───────────────────
  "Membrane transporter: ATP-dependent (ABC)" = list(
    func = "(?i)(atp[[:punct:][:space:]_-]*binding[[:punct:][:space:]_-]*cassette|\\babc\\b|\\bmalk\\b).*(transporter|permease|protein|atpase([[:space:]_-]*component)?|atp[[:punct:][:space:]_-]*binding|superfamily|family)",
    loc  = "cytoplasmic membrane|outer membrane|periplasmic|unassigned"
  ),

  "Membrane transporter: secondary (MFS/SSS/RND/Amt/TRAP)" = list(
    func = paste0(
      "(?i)(",
      "\\bmfs\\b|mfs transporter|major facilitator|",
      "cation.*exchanger|outer membrane channel|monovalent cation:proton antiporter|",
      "membrane protein insertase yidc|cereblon.*\\bcrbn\\b.*family|",
      "\\bsss\\b|solute[[:punct:][:space:]_:-]*sodium[[:punct:][:space:]_:-]*symporter|sodium[[:punct:][:space:]_:-]*solute[[:punct:][:space:]_:-]*symporter|vc_2705|",
      "\\brnd (superfamily|transporter)\\b|\\bamt\\b|ammonium transporter|",
      "trap[[:punct:][:space:]_-]*type(?!.*(receptor|component))|",
      "tripartite.*periplasmic transporter(?!.*(receptor|component))|",
      "the[[:space:]]+tripartite[[:space:]]+atp[- ]independent periplasmic transporter \\(trap[- ]t\\) family|",
      "the tricarboxylate transporter \\(ttt\\) family|",
      "tripartite[- ]type tricarboxylate transporter|",
      "\\btct[ab]\\b|tricarboxylate transporter|",
      "large permease component|",
      "sodium/solute symporter|",
      "the putative peptide transporter carbon starvation csta \\(csta\\) family|\\bcsta\\b|",
      "na\\+(or h\\+)/acetate symporter actp|\\bactp\\b|",
      "phosphate/sulfate permease",
      ")"
    ),
    loc  = "cytoplasmic membrane|outer membrane|periplasmic|unassigned"
  ),

  # ── Secretion systems (beyond Sec/Tat already sprinkled elsewhere) ───────
  "Secretion systems (T2SS/T3SS/T4SS/T6SS)" = list(
    func = paste0(
      "(?i)(",
      "type ii secretion system|\\bgsp[cd]\\b|secretin|",
      "type iii secretion system|injectisome|\\bsct[a-z]\\b|",
      "type iv secretion system|\\bvirb\\d+\\b|trb|",
      "type vi secretion system|\\btss[a-z]\\b|\\bhcp\\b|\\bvgrg\\b",
      ")"
    ),
    loc  = "outer membrane|periplasmic|cytoplasmic membrane|unassigned"
  ),

  # ── Oxidative stress / redox ─────────────────────────────────────────────
  "Superoxide dismutase (NiSOD)" =
    "(?i)(\\bnisod\\b|\\bsodn\\b|\\bni[[:punct:][:space:]_-]*sod\\b|superoxide dismutase[[:punct:][:space:]_-]*ni\\b)",

  "Superoxide dismutases" =
    "(?i)(\\bsod[amfg]?\\b|superoxide dismutase)",

  "Ferritin / DPS" =
    "(?i)(ferritin|\\bdps\\b|bacterioferritin|ferritin/dps protein domain|dna[- ]binding ferritin[- ]like)",

  "Thioredoxin / glutaredoxin" =
    "(?i)(thioredoxin|glutaredoxin|glutathione peroxidase|peptide[[:punct:][:space:]_-]*methionine.*oxide reductase|\\bmsra\\b|thioredoxin reductase)",

  "Alkyl hydroperoxide / peroxiredoxin" =
    "(?i)(alkyl hydroperoxide|peroxiredoxin|redoxin|thiol[[:punct:][:space:]_-]*specific antioxidant|ahp[c])",

  "Flavodoxin / ferredoxin" =
    "(?i)(flavodoxin|ferredoxin|ferredoxin[[:punct:][:space:]_-]*nadp[+\\-]?[[:punct:][:space:]_-]*reductase|\\bpeth\\b|pet[h])",

  "Oxidoreductases (other)" =
    paste0("(?i)((\\boxidoreductase\\b|\\bshort[[:punct:][:space:]_-]*chain[[:punct:][:space:]_-]*dehydrogenase\\b|",
           "\\bdihydrolipoyl\\b|\\bdehydrogenase\\b|\\breductase\\b)",
           "(?!.*\\b(abc|ferredoxin[[:punct:][:space:]_-]*nadp|fnr|nitrit|sulfi[tte]|nitra[tte]|ribonucleotid|ribonucleosid|ketol[[:punct:][:space:]_-]*acid|geranylgeranyl|glucose[[:punct:][:space:]_-]*6[[:punct:][:space:]_-]*phosphate)\\b))"),

  # ── Nitrogen & sulfur metabolism ─────────────────────────────────────────
  "Nitrogen regulation (PII)" =
    "(?i)(nitrogen regulatory protein pii|\\bpii\\b|\\bglnb\\b|\\bntrc\\b|\\bnifh\\b|ntc[a])",

  "Urease / urea transport" =
    "(?i)(\\burease\\b|\\burt[abcdf]\\b|urea[[:punct:][:space:]_-]*abc transporter|urea abc transporter (substrate|atp)[- ]binding subunit|\\burt[e]\\b|urt[d])",

  "Nitrate / nitrite transport" =
    "(?i)(\\bnrt[abcd]\\b|\\bnark\\b|nitrate[[:punct:][:space:]_-]*nitrite mfs transporter|nitrate[[:punct:][:space:]_-]*abc)",

  "Sulfur metabolism" =
    "(?i)(sulfate adenylyltransferase|adenylyl[[:punct:][:space:]_-]*sulfate kinase|sulfate transporter|\\bcys[ijmk]\\b|\\bsulp\\b|sulfite reductase|o[[:punct:][:space:]_-]*acetylserine[[:punct:][:space:]_-]*sulfhydrylase)",

  # ── Proteostasis ─────────────────────────────────────────────────────────
  "Chaperones (GroEL/DnaK)" =
    "(?i)(chaperonin|\\bgroel\\b|\\bgroes\\b|\\bdnaj\\b|\\bdnak\\b|\\bcpn60\\b|hsp60|hsp70|\\bhsp70\\b|molecular chaperone|small heat shock|hsp20|ibp[a-z]|hsp90|molecular chaperone, hsp90 family|trigger factor)",

  "Protein folding (PPIase / cyclophilin)" =
    "(?i)(peptidyl[[:punct:][:space:]_-]*prolyl|cyclophilin|\\bppiase\\b|fkbp)",

  "Protease / AAA-ATPase" =
    "(?i)(protease|peptidase|\\bfts[h]\\b|\\blon\\b|\\bm24\\b|\\bclpb\\b|\\bclpa\\b|\\bclpp\\b|\\bclpc\\b|aaa[[:punct:][:space:]_-]*type.*core|aatpase.*aaa|endopeptidase la|atp[- ]dependent zinc metalloprotease fts[h]|atp[- ]dependent zn protease|do family serine endopeptidase|signal peptide peptidase)",

  # ── Gene expression ──────────────────────────────────────────────────────
  "Ribosomal proteins" =
    "(?i)(([35]0s ribosomal|\\bribosomal protein\\b|large ribosomal subunit|small ribosomal subunit|\\bpsrp\\b|\\bribosomal protein psrp\\b|\\brps\\d+\\b|\\brpl\\d+\\b|l\\d+p?/l\\d+e?))",

  "Translation factors" =
    "(?i)(elongation factor|initiation factor|translation elongation factor|\\bef[[:punct:][:space:]_-]*tu\\b|\\bif[[:punct:][:space:]_-]*[123]\\b|\\blepg?\\b|\\blepa\\b|translation termination factor rho|rho factor|translation elongation factor ts)",

  "Aminoacyl-tRNA synthetase" =
    "(?i)(trna ligase|trna synthetase|aminoacyl[[:punct:][:space:]_-]*trna|\\bgatb\\b|aspartate[[:punct:][:space:]_-]*trna ligase|leucine[[:punct:][:space:]_-]*trna ligase|arginine[[:punct:][:space:]_-]*trna ligase|cysteine[[:punct:][:space:]_-]*trna ligase|proline[[:punct:][:space:]_-]*trna ligase|serine[[:punct:][:space:]_-]*trna ligase|alanine[[:punct:][:space:]_-]*trna ligase)",

  "Two-component signalling" =
    "(?i)(histidine kinase|response regulator|\\bhamp\\b|omp[r]|nar[l]|fixj|signal transduction histidine kinase|signal transduction response regulator|sensor histidine kinase)",

  "RNA polymerase / processing" =
    "(?i)(dna[[:punct:][:space:]_-]*directed rna polymerase|\\brpo[abcd]\\b|polyribonucleotide|\\bpnpase\\b|\\brnase\\b|ribonuclease|\\brnase j\\b|beta[[:punct:][:space:]_-]*casp ribonuclease|rna polymerase rpb1|rna polymerase rpb2|sigma[[:punct:][:space:]_-]*70.*rna polymerase sigma|sigma factor|dna[- ]directed rna polymerase subunit (alpha|beta|beta''|gamma))",

  "RNA-binding / processing" =
    "(?i)(rna recognition motif|s1 domain|\\brrm\\b|rna[[:punct:][:space:]_-]*binding protein|au[[:punct:][:space:]_-]*1/ribonuclease|nuclear mrna exporter|nuclear pore|single[- ]stranded dna[- ]binding protein|ybab/ebfc family|ylqd)",

  "Transcription regulators" =
    "(?i)(transcription (regulator|factor|repressor|activator|termination)|\\babrb[[:punct:][:space:]_-]*like\\b|\\bntca\\b|\\bluxr\\b|\\bnrdr\\b|\\bnusa\\b|\\bnusg\\b|lysr family|omp[r] family|nar[l]/fixj|winged[- ]helix|wht[h]|anti[[:punct:][:space:]_-]*sigma factor antagonist|abrb[- ]like)",

  # ── Central metabolism ───────────────────────────────────────────────────
  "Carbohydrate metabolism" =
    "(?i)(enolase|pyruvate kinase|pfkb|aldolase|class [i|ii][[:punct:][:space:]_-]*fructose[[:punct:][:space:]_-]*bisphosphate aldolase|carbohydrate kinase|glycoside hydrolase|glycogen.*phosphorylase|alpha[[:punct:][:space:]_-]*glucan.*phosphorylase|rhamnonate aldolase|phosphoglucose isomerase|phosphoglycerate kinase|fructose[- ]1,6[- ]bisphosphatase|sedoheptulose[- ]1,7[- ]bisphosphatase|transketolase|glucose[- ]6[- ]phosphate dehydrogenase)",

  "Amino acid metabolism" =
    "(?i)(glutamate[[:punct:][:space:]_-]*5[[:punct:][:space:]_-]*semialdehyde|aspartate.*dehydrogenase|argininosuccinate|branched[[:punct:][:space:]_-]*chain[[:punct:][:space:]_-]*amino[[:punct:][:space:]_-]*acid|glutamine synthetase|glutamine[- ]hydrolyzing carbamoyl[- ]phosphate synthase|carbamoyl phosphate|acetylglutamate|agmatinase|4[[:punct:][:space:]_-]*hydroxy[[:punct:][:space:]_-]*tetrahydrodipicolinate|delta[[:punct:][:space:]_-]*aminolevulinic acid|anthranilate|acetolactate synthase|ketol[[:punct:][:space:]_-]*acid reductoisomerase|glutamyl[[:punct:][:space:]_-]*trna reductase|asparagine synthase|aminotransferase\\b|branched[- ]chain[- ]amino[- ]acid transaminase|cysteine synthase|o[- ]acetylhomoserine aminocarboxylpropyltransferase)",

  "Nucleotide metabolism" =
    "(?i)(dtmp kinase|ump kinase|ribonucleoside[[:punct:][:space:]_-]*triphosphate reductase|adenine phosphoribosyl|thymidylate synthase|adenosylhomocysteinase|adenylosuccinate lyase|nucleoside[- ]diphosphate kinase|phosphoribosyl[[:punct:][:space:]_-]*atp|ribose[- ]phosphate diphosphokinase|nrd[abj]|purm[[:punct:][:space:]_-]*like|puros|pseudouridine synthase|cyclic nucleotide[[:punct:][:space:]_-]*binding|methionine adenosyltransferase)",

  "Lipid / fatty acid metabolism" =
    "(?i)(3[[:punct:][:space:]_-]*oxoacyl|enoyl[[:punct:][:space:]_-]*acyl|acetyl[[:punct:][:space:]_-]*coa carboxylase|acyl[[:punct:][:space:]_-]*?acp|acyltransferase|fatty acid|beta[[:punct:][:space:]_-]*ketoacyl[[:punct:][:space:]_-]*acp synthase|phosphate acyltransferase plsx|acyl[[:punct:][:space:]_-]*acp[[:punct:][:space:]_-]*.*udp[[:punct:][:space:]_-]*n[[:punct:][:space:]_-]*acetylglucosamine|sterol carrier protein|\\bscp2\\b|enoyl[- ]coa hydratase|acyl[- ]coa reductase|acryloyl[- ]coa reductase|acetoacetate decarboxylase|2[- ]methylcitrate synthase|citrate synthase ii|enoyl[- ]coa hydratase/carnithine racemase)",

  "Cofactor / vitamin biosynthesis" =
    "(?i)(thiazole synthase|lumazine synthase|chorismate|1[[:punct:][:space:]_-]*deoxy[[:punct:][:space:]_-]*d[[:punct:][:space:]_-]*xylulose|cobalamin|thiamine|folate|biotin|riboflavin|6,7[[:punct:][:space:]_-]*dimethyl[[:punct:][:space:]_-]*8[[:punct:][:space:]_-]*ribityllumazine|magnesium[[:punct:][:space:]_-]*protoporphyrin|delta[[:punct:][:space:]_-]*aminolevulinic acid|precorrin|bifunctional hydroxymethylpyrimidine|\\bnmt1[[:punct:][:space:]_-]*like\\b|pyridoxal phosphate|\\bpdxj\\b|2[[:punct:][:space:]_-]*c[[:punct:][:space:]_-]*methyl[[:punct:][:space:]_-]*d[[:punct:][:space:]_-]*erythritol|cyanase|amp[[:punct:][:space:]_-]*dependent synthetase|coenzyme q[[:punct:][:space:]_-]*binding|inositol monophosphatase|inorganic pyrophosphatase|\\bthij\\b|\\bpfpi\\b|heme oxygenase|\\bmapeg\\b|alkaline phosphatase|histidine phosphatase|broad specificity phosphatase|pho[e])",

  "Carotenoid / pigment biosynthesis" =
    "(?i)(\\bcrt[a-z]\\b|phytoene synthase|lycopene cyclase|zeaxanthin|carotenoid|phycobilin lyase|heme oxygenase ho1|orange carotenoid[- ]binding protein|phytoene/squalene synthetase|magnesium chelatase)",

  # ── Cell structure / surface / misc ──────────────────────────────────────
  "Cell division / cytoskeleton" =
    "(?i)(\\bftsz\\b|\\bmine\\b|\\bmreb\\b|mincd|cell division|septal|topological specificity factor mine|\\bmrec\\b|bactofilin|tubulin/ftsz|tubulin family|actin[- ]related protein|mreb/mrl|mreb)",

  "DNA topology / replication" =
    "(?i)(dna gyrase|topoisomerase|dna replication|dna repair|recombinase|\\bdnab\\b|dna polymerase|chromosomal replication initiator|type i dna topoisomerase|uvr[ab]|excinuclease|replicative dna helicase|dna polymerase (i|iii subunit beta))",

  "DNA-binding (HU / histone-like)" =
    "(?i)(histone[[:punct:][:space:]_-]*like|\\bhu family\\b|dna[[:punct:][:space:]_-]*binding protein|\\bihf\\b|\\bhup\\b|histone[[:punct:][:space:]_-]*like dna[[:punct:][:space:]_-]*binding|h4[- ]like protein|h2a)",

  "Stomatin / Band 7 / SPFH" =
    "(?i)(band 7|stomatin|\\bspfh\\b|stomatin/podocin|prohibitin|\\bhfl[ck]\\b|flotillin|reggie|regulator of protease activity hfl[ck])",

  "Pentapeptide repeat" =
    "(?i)(pentapeptide repeat)",

  "Tetratricopeptide repeat (TPR)" =
    "(?i)(tetratricopeptide repeat|\\btpr repeat\\b|sh3 domain)",

  "Methyltransferase" =
    "(?i)(methyltransferase)",

  "RTX toxin" =
    "(?i)(pore[[:punct:][:space:]_-]*forming rtx toxin|rtx[[:punct:][:space:]_-]*toxin|rtx[- ]related|ca2\\+[[:punct:][:space:]_-]*binding protein)",

  "Glycan / starch enzymes" =
    "(?i)(glycosyl hydrolase|glucan branching enzyme|glycogen.*starch synthase|fructose[[:punct:][:space:]_-]*1,6[[:punct:][:space:]_-]*bisphosphatase|sedoheptulose[[:punct:][:space:]_-]*1,7[[:punct:][:space:]_-]*bisphosphatase|age family epimerase|nad[[:punct:][:space:]_-]*dependent epimerase/dehydratase|dihydroxy[[:punct:][:space:]_-]*acid.*dehydratase)",

  # ── Protein secretion / translocation ────────────────────────────────────
  "Protein secretion / Sec-Tat" =
    "(?i)(preprotein translocase|\\bseca\\b|\\bsecy\\b|yajc|tat pathway|twin arginine translocase|tat[a-c])",

  "Ion channels / mechanosensitive" =
    "(?i)(mechanosensitive ion channel|\\bmscs\\b|chloride channel|regulator of k\\+[[:punct:][:space:]_-]*conductance|small gtp[[:punct:][:space:]_-]*binding protein|cobw/hypb/ureg|small gtp[- ]binding protein domain)",

  # ── Additional systems ───────────────────────────────────────────────────
  "Cyanophycin metabolism" =
    "(?i)(cyanophycin|\\bcph[ab]\\b|cyanophycin synthetase|cyanophycinase)",

  "Gas vesicles" =
    "(?i)(\\bgvp[a-z]\\b|gas[[:punct:][:space:]_-]*vesicle)",

  "Motility / pili / competence" =
    "(?i)(\\bpil[a-z]\\b|type iv pil(us)?|\\bpilq\\b|\\bcomea\\b|\\bcomec\\b|competence|tad[abcde]?|flagellin|flgl|fimv|type ii secretion system secretin|gspd)",

  "CRISPR-Cas" =
    "(?i)(\\bcas\\d*\\b|crispr|\\bcsy\\b|\\bcsm\\b|\\bcmr\\b)",

  "LPS / lipid A biosynthesis" =
    "(?i)(\\blpx[abdkhl]\\b|\\bwaa\\w\\b|\\bkdt?a\\b|kdo transferase|lipid a|lps biosynthesis|kds[a-z]?|3[[:punct:][:space:]_-]*deoxy[[:punct:][:space:]_-]*8[[:punct:][:space:]_-]*phosphooctulonate synthase|2[[:punct:][:space:]_-]*keto[[:punct:][:space:]_-]*3[[:punct:][:space:]_-]*deoxy[[:punct:][:space:]_-]*octonate|kdo synthase)",

  "Peptidoglycan / cell wall" =
    "(?i)(\\bmur[acdefgl]\\b|\\bddl\\b|\\bpbp\\b|penicillin[[:punct:][:space:]_-]*binding protein|murein|carboxypeptidase|omp[a]\\b|peptidoglycan[- ]associated lipoprotein|\\bpal\\b|outer membrane protein ompa)",

  "Surface / repeat proteins" =
    "(?i)(bacterial surface antigen|\\bd15\\b|beta strand repeat|vhs domain|eve domain|extracellular substrate binding[[:punct:][:space:]_-]*like|\\bgrrp\\b|circularly permuted atp[[:punct:][:space:]_-]*grasp|carbon[[:punct:][:space:]_-]*nitrogen hydrolase|nucleotidyl transferase domain|2fe[[:punct:][:space:]_-]*2s iron[[:punct:][:space:]_-]*sulfur|beta[[:punct:][:space:]_-]*lactamase superfamily|hsp90/cdc37|apoptosis cell death regulator|protein continuous vascular ring|predicted lipid[- ]binding transport protein, tim44 family|ragb/susd domain|t9ss type a sorting domain)"
)

bacteria_categories_list <- list(
  "Genome organization and gene expression" = c(
    "DNA topology / replication",
    "DNA-binding (HU / histone-like)",
    "Transcription regulators",
    "Two-component signalling",
    "RNA polymerase / processing",
    "RNA-binding / processing",
    "Translation factors",
    "Aminoacyl-tRNA synthetase",
    "Ribosomal proteins",
    "CRISPR-Cas"
  ),

  "Proteostasis (folding, chaperones, proteolysis)" = c(
    "Protein folding (PPIase / cyclophilin)",
    "Chaperones (GroEL/DnaK)",
    "Protease / AAA-ATPase"
  ),

  "Transport and uptake" = c(
    "Outer membrane porins (other)",
    "Carbohydrate-selective porins",
    "Iron uptake porins",
    "Solute binding proteins (SBP)",
    "TonB-dependent receptor (TBDR)",
    "SusC (polysaccharide utilization loci)",
    "SusD (polysaccharide utilization loci)",
    "TonB system (TonB/ExbB/ExbD/Tol)",
    "Membrane transporter (ABC)",
    "Membrane transporter (MFS/SSS/RND/Amt/TRAP)",
    "Nitrate / nitrite transport",
    "Urease / urea transport",
    "Nitrogen regulation (PII)"
  ),

  "Phototrophy and energy transduction" = c(
    "Photosystem I (PSI)",
    "Photosystem II (PSII)",
    "Phycobilisome and phycobiliproteins",
    "Anoxygenic reaction centers (PRC/BRC)",
    "Rhodopsin",
    "Respiratory chain / NDH",
    "V-type H(+)-translocating pyrophosphatase (V-PPase)",
    "ATP synthase (membrane)",
    "ATP synthase (cytoplasmic)",
    "Flavodoxin / ferredoxin"
  ),

  "Carbon fixation and core metabolism" = c(
    "Carbon fixation (RuBisCO / Calvin cycle)",
    "Carbohydrate metabolism",
    "Glycan / starch enzymes",
    "Lipid / fatty acid metabolism",
    "Amino acid metabolism",
    "Nucleotide metabolism",
    "Cofactor / vitamin biosynthesis",
    "Carotenoid / pigment biosynthesis",
    "Bacterial microcompartments",
    "Cyanophycin metabolism"
  ),

  "Catch-all / unassigned" = c(
    "Other",
    "Alkyl hydroperoxide / peroxiredoxin",
    "Thioredoxin / glutaredoxin",
    "Superoxide dismutase (NiSOD)",
    "Superoxide dismutases",
    "Ferritin / DPS",
    "Oxidoreductases (other)",
    "Hypothetical / unknown",
    "Methyltransferase",
    "Cell division / cytoskeleton",
    "Stomatin / Band 7 / SPFH",
    "Surface / repeat proteins",
    "Pentapeptide repeat",
    "Tetratricopeptide repeat (TPR)",
    "Motility / pili / competence",
    "Protein secretion / Sec-Tat",
    "Secretion systems (T2SS/T3SS/T4SS/T6SS)",
    "Ion channels / mechanosensitive",
    "BAM complex / OM assembly",
    "LPS / lipid A biosynthesis",
    "Peptidoglycan / cell wall",
    "Gas vesicles"
  )
)



