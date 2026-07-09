library(dplyr)
library(purrr)  
library(readr)
library(tidyr)
library(phyloseq)
library(readxl)
library(ggplot2)
library(openxlsx)
library(tidyverse)

#set working directory
setwd("/Users/kristinasorokolet/Library/CloudStorage/Box-Box/Wright Lab Operations/Student Folders/Kris/KS_2/Proteomics/GLOBAL/RANDOM_BINS/males_CTABS/MS9_1P_75_unclus/KEGG_new/")

metadata_file <- "metadata.csv"

# -------------------------------
# Read metadata
# -------------------------------
meta <- read_csv(metadata_file, show_col_types = FALSE)

meta <- meta %>%
  mutate(group = recode(group,
                        "group1" = "dss_hfid",
                        "group2" = "dss",
                        "group3" = "control",
                        "group4" = "hfid"))


library(tidyverse)
library(KEGGREST)
library(pheatmap)
library(readr)

###################################################################
#CONVERT EGGNOG ANNOTATION FILE TO EXCEL FOR EASY VIEW
####################################################################
eggnog_files <- list.files(
  pattern = "emapper\\.(annotations|hits|seed_orthologs)$",
  full.names = TRUE
)

for (f in eggnog_files) {
  message("Processing: ", basename(f))
  
  # Read all lines
  lines <- readLines(f)
  
  # Remove meta lines starting with double ##
  lines <- lines[!grepl("^##", lines)]
  
  # Find header line (starts with #query), fallback to first line
  header_line <- lines[grepl("^#query", lines)][1]
  if (is.na(header_line)) {
    header_line <- lines[1]  # fallback
    message("⚠️ No #query header found, using first line as header")
  }
  
  # Remove leading # if present and split by tab
  header <- strsplit(sub("^#", "", header_line), "\t")[[1]]
  
  # Identify data lines
  data_start <- which(lines == header_line)[1] + 1
  if (data_start > length(lines)) {
    message("⚠️ No data found for this file, skipping")
    next
  }
  
  data_lines <- lines[data_start:length(lines)]
  
  # Read data
  df <- read.table(
    text = paste(data_lines, collapse = "\n"),
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = ""
  )
  
  # Assign column names
  colnames(df) <- header
  
  # Write Excel
  out_xlsx <- paste0(basename(f), ".xlsx")
  write.xlsx(df, file = out_xlsx, overwrite = TRUE)
}

message("✅ All eggNOG files converted to Excel (headers preserved)")

###################################################################################
#MERGE EGGNOG ANNOTATED KEGG PATHWAYS INTO THE PROCESSED PROTEIN CROSSTAB
###################################################################################

library(dplyr)
library(readr)
library(stringr)

crosstab <- read_tsv("imputed_data.tsv")

library(readr)

# Read all lines, skipping lines that start with ##
all_lines <- readLines("combined_fasta_males.emapper.annotations")
data_lines <- all_lines[!grepl("^##", all_lines)]

# Find the header (first line starting with #query or first line)
header_line <- data_lines[grepl("^#query", data_lines)][1]
if (is.na(header_line)) header_line <- data_lines[1]

# Remove leading # if present and split by tab
header <- strsplit(sub("^#", "", header_line), "\t")[[1]]

# Read the data, skipping everything above header_line
start_row <- which(data_lines == header_line)[1] + 1
eggnog <- read_tsv(
  file = I(paste(data_lines[start_row:length(data_lines)], collapse = "\n")),
  col_names = header,
  col_types = cols(.default = "c")
)

# Rename first column to Protein
colnames(eggnog)[1] <- "Protein"

head(eggnog)

eggnog_kegg <- eggnog %>%
  select(
    Protein,
    KEGG_Pathway,
    EC
  )

colnames(eggnog) <- c(
  "Protein",
  "Seed_ortholog",
  "Evalue",
  "Score",
  "EggNOG_OGs",
  "Max_annot_lvl",
  "COG_category",
  "Description",
  "Preferred_name",
  "GOs",
  "EC",
  "KEGG_ko",
  "KEGG_Pathway",
  "KEGG_Module",
  "KEGG_Reaction",
  "KEGG_rclass",
  "BRITE",
  "KEGG_TC",
  "CAZy",
  "BiGG_Reaction",
  "PFAMs"
)[seq_len(ncol(eggnog))]

eggnog_kegg <- eggnog %>%
  mutate(
    KEGG_Pathway = str_remove_all(KEGG_Pathway, "map[0-9]{5},?"),
    KEGG_Pathway = str_remove_all(KEGG_Pathway, ",$"),
    KEGG_Pathway = na_if(KEGG_Pathway, "")
  ) %>%
  filter(!is.na(KEGG_Pathway)) %>%
  select(Protein, KEGG_Pathway)

crosstab_kegg <- crosstab %>%
  left_join(eggnog_kegg, by = "Protein")

crosstab_kegg %>%
  summarise(
    total_proteins = n(),
    annotated = sum(
      !is.na(KEGG_Pathway) & KEGG_Pathway != "-"
    ),
    unannotated = sum(
      is.na(KEGG_Pathway) | KEGG_Pathway == "-"
    )
  )

# Save as TSV
write_tsv(crosstab_kegg, "imputed_w_taxa_kegg.tsv")

# Save as CSV
write_csv(crosstab_kegg, "imputed_w_taxa_kegg.csv")

####################################################################################################################
#COLLAPSE PROTEINS INTO KEGG PATHWAYS 
####################################################################################################################

library(dplyr)
library(tidyr)
library(KEGGREST)
library(stringr)

# -----------------------------
# Read data
# -----------------------------
crosstab <- read_csv("imputed_w_taxa_kegg.csv")

# If your data is log2-transformed, convert back to linear
df_linear <- crosstab %>%
  mutate(across(starts_with("LI"), ~ 2^. - 1))

sample_cols <- grep("^LI", colnames(df_linear), value = TRUE)

# -----------------------------
# Compute per-protein relative abundance
# -----------------------------
df_long <- df_linear %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  group_by(SampleID) %>%
  mutate(RelAbundance = Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

# -----------------------------
# Split multi-KO proteins and divide abundance
# -----------------------------
df_long <- df_long %>%
  separate_rows(KEGG_Pathway, sep = ",\\s*") %>%
  filter(!is.na(KEGG_Pathway) & KEGG_Pathway != "") %>%
  group_by(Protein, SampleID) %>%
  mutate(RelAbundance = RelAbundance / n()) %>%  # divide protein among KOs
  ungroup()

# -----------------------------
# Remove unwanted KOs (ambiguous or host associated, misannotated)
# -----------------------------
kos_to_remove <- c(
  "ko04212", "ko04626", "ko05165", "ko04930", "ko05203", 
  "ko00195", "ko05010", "ko05215", "ko05200", "ko04114", 
  "ko05152", "ko04940", "ko04113", "ko04213", "ko04217", 
  "ko0100", "ko04214", "ko05418", "ko05134", "ko04016",
  "ko05206", "ko05016", "ko04013", "ko05230", "ko04112", 
  "ko0464", "ko4141", "ko01100", "ko04146", "ko03013",
  "ko04964"
)

df_long <- df_long %>%
  filter(!KEGG_Pathway %in% kos_to_remove)

# -----------------------------
# Aggregate per KEGG pathway
# -----------------------------
df_kegg_rel <- df_long %>%
  group_by(KEGG_Pathway, SampleID) %>%
  summarise(RelAbundance = sum(RelAbundance, na.rm = TRUE), .groups = "drop")

# -----------------------------
# Pivot back to wide format for downstream use
# -----------------------------
df_kegg_rel_clean_filtered <- df_kegg_rel %>%
  pivot_wider(names_from = SampleID, values_from = RelAbundance, values_fill = 0)

# -----------------------------
# Add KEGG pathway names
# -----------------------------
all_kegg <- keggList("pathway", "ko")  # fetch all KO → pathway names

core_kos <- df_kegg_rel_clean_filtered$KEGG_Pathway
core_kos <- core_kos[core_kos %in% names(all_kegg)]

kegg_info_df <- data.frame(
  KEGG_Pathway = core_kos,
  pathway_name = all_kegg[core_kos],
  stringsAsFactors = FALSE
)

df_kegg_rel_clean_filtered <- df_kegg_rel_clean_filtered %>%
  left_join(kegg_info_df, by = "KEGG_Pathway")

# -----------------------------
# Fill missing pathways manually if needed (search on kegg online database)
# -----------------------------
manual_paths <- data.frame(
  KEGG_Pathway = c("ko00072", "ko00471", "ko00473", "ko01130"),
  pathway_name = c("Synthesis and degradation of ketone bodies",
                   "D-Glutamine and D-Glutamate metabolism",
                   "D-Alanine metabolism",
                   "Ubiquinone and other terpenoid quinone biosynthesis"),
  stringsAsFactors = FALSE
)

df_kegg_rel_clean_filtered <- df_kegg_rel_clean_filtered %>%
  rows_update(manual_paths, by = "KEGG_Pathway")

write_xlsx(
  df_kegg_rel_clean_filtered,
  "df_kegg_rel_collapsed.xlsx"
)

########################################################################################################################
#HEATMAP OF ALL KEGG PATHWAYS IDENTIFIED
########################################################################################################################
library(dplyr)
library(pheatmap)
library(viridis)

# -------------------------
# with clustered columns
# -------------------------

# Select sample columns (starting with LI)
sample_cols <- colnames(df_kegg_rel_clean_filtered)[grepl("^LI", colnames(df_kegg_rel_clean_filtered))]

# Combine pathway name + KO ID for row labels
df_labels <- df_kegg_rel_clean_filtered %>%
  mutate(row_label = paste0(pathway_name, " (", KEGG_Pathway, ")"))

# Create numeric matrix (only keeps row_label and sample abundance, makes row label rowname)
mat <- df_labels %>%
  select(row_label, all_of(sample_cols)) %>%
  column_to_rownames(var = "row_label") %>%
  as.data.frame()

# Ensure numeric matrix
mat[] <- lapply(mat, as.numeric)
mat <- as.matrix(mat)

# Row-wise Z-score (canonical)
mat_z <- t(scale(t(mat)))


# Match metadata
meta_li <- meta[meta$sample_id %in% colnames(mat), , drop = FALSE]
meta_li <- meta_li[match(colnames(mat), meta_li$sample_id), , drop = FALSE]

# Annotation dataframe for pheatmap
anno <- data.frame(Group = meta_li$group)
rownames(anno) <- meta_li$sample_id

group_colors <- list(
  Group = c(
    control  = "#0072B2",  # strong blue
    hfid     = "#E69F00",  # vivid orange
    dss      = "#CC79A7",  # magenta/purple
    dss_hfid = "#009E73"   # green
  )
)

png(
  filename = "kegg_pathway_heatmap_z-score_imputed_col.png",
  width = 3000,
  height = 7000,
  res = 300
)

pheatmap(
  mat_z,
  scale = "none",          
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 0,
  treeheight_col = 10,
  
  annotation_col = anno,
  annotation_colors = group_colors,
  
  border_color = NA,
  na_col = "black",
  color = viridis::plasma (100),
  
  main = "KEGG Pathways - males"
)

dev.off()

# -------------------------
#with unclustered columns 
# -------------------------
library(dplyr)
library(pheatmap)
library(tibble)
library(viridis)

# -----------------------------
# Select sample columns
# -----------------------------
sample_cols <- colnames(df_kegg_rel_clean_filtered)[grepl("^LI", colnames(df_kegg_rel_clean_filtered))]

# -----------------------------
# Create row labels
# -----------------------------
df_labels <- df_kegg_rel_clean_filtered %>%
  mutate(row_label = paste0(pathway_name, " (", KEGG_Pathway, ")"))

# -----------------------------
# Create numeric matrix
# -----------------------------
mat <- df_labels %>%
  select(row_label, all_of(sample_cols)) %>%
  column_to_rownames(var = "row_label") %>%
  as.data.frame()

mat[] <- lapply(mat, as.numeric)
mat <- as.matrix(mat)

# -----------------------------
# Row-wise Z-score
# -----------------------------
mat_z <- t(scale(t(mat)))

# -----------------------------
# Match metadata to matrix
# -----------------------------
meta_li <- meta[meta$sample_id %in% colnames(mat_z), , drop = FALSE]
meta_li <- meta_li[match(colnames(mat_z), meta_li$sample_id), , drop = FALSE]

anno <- data.frame(Group = meta_li$group)
rownames(anno) <- meta_li$sample_id

# -----------------------------
# Enforce group order
# -----------------------------
group_levels <- c("control", "hfid", "dss", "dss_hfid")
anno$Group <- factor(anno$Group, levels = group_levels)

# Order samples
ord <- order(anno$Group)

mat_z_ord <- mat_z[, ord]
anno_ord  <- anno[ord, , drop = FALSE]

# -----------------------------
# check this
# -----------------------------
stopifnot(all(colnames(mat_z_ord) == rownames(anno_ord)))

# -----------------------------
# group colors
# -----------------------------
group_colors <- list(
  Group = c(
    control  = "#0072B2",
    hfid     = "#E69F00",
    dss      = "#CC79A7",
    dss_hfid = "#009E73"
  )
)

# -----------------------------
# Plot
# -----------------------------
png(
  filename = "kegg_pathway_heatmap_z-score_imputed.png",
  width = 3000,
  height = 7000,
  res = 300
)

pheatmap(
  mat_z_ord,
  scale = "none",          
  cluster_rows = TRUE,
  cluster_cols = FALSE,  
  treeheight_row = 0,
  
  annotation_col = anno_ord,
  annotation_colors = group_colors,
  
  border_color = NA,
  na_col = "black",
  color = viridis::plasma(100),
  
  main = "KEGG Pathways - males"
)

dev.off()


############################################################################################################
# BROAD CATEGORIES KEGG HEATMAP
############################################################################################################

library(dplyr)
library(stringr)
library(KEGGREST)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(grid)

# --- Extract KEGG IDs ---
df_labels <- data.frame(
  rowname = rownames(mat_z),
  stringsAsFactors = FALSE
) %>%
  mutate(KEGG_ID = str_extract(rowname, "ko\\d{5}"))

# --- Fetch KEGG categories ---
get_category <- function(pid) {
  info <- tryCatch(keggGet(pid)[[1]], error = function(e) NULL)
  if (!is.null(info$CLASS)) info$CLASS else NA
}

valid_ids <- unique(df_labels$KEGG_ID)
valid_ids <- valid_ids[!is.na(valid_ids) & valid_ids != "-"]

pathway_categories <- data.frame(
  KEGG_ID = valid_ids,
  Category = sapply(valid_ids, get_category),
  stringsAsFactors = FALSE
)

# ---  Merge categories back to the original tables ---
df_labels <- df_labels %>%
  left_join(pathway_categories, by = "KEGG_ID")

# --- Parse KEGG hierarchy ---
df_labels <- df_labels %>%
  mutate(
    Level3 = trimws(sapply(strsplit(Category, ";"), `[`, 3))
  )

library(openxlsx)

# Save current df_labels
write.xlsx(df_labels, file = "df_labels.xlsx", rowNames = FALSE)
#now manually edit missing assignments and ensure broad categories are assigned properly

# -------------------------
#now use edited label file, upload it to break kegg pathways into categories
# -------------------------

library(ComplexHeatmap)
library(circlize)
library(viridis)
library(grid)
library(openxlsx)
library(dplyr)
library(stringr)

df_labels <- read.xlsx("df_labels_for_manual_edit.xlsx") %>%
  as.data.frame(stringsAsFactors = FALSE)

colnames(df_labels)[1] <- "rowname"

# If any rows have zero variance, scale() produces NA
# You said NA is fine, so we leave them as-is

# --- Keep only matching rownames ---
common_rows <- intersect(df_labels$rowname, rownames(mat_z))

# --- Subset BOTH objects ---
mat_z<- mat_z[common_rows, , drop = FALSE]
df_labels  <- df_labels[df_labels$rowname %in% common_rows, ]

# --- enforce identical order ---
df_labels <- df_labels[match(rownames(mat_z), df_labels$rowname), ]

# ---Assign row_split by KEGG Level2 ---
row_splits <- df_labels$Level2
row_splits[is.na(row_splits)] <- "Other"

# --- Keep row order as in mat_scaled ---
# This preserves original Z-scores

# --- Clean rownames for display ---
rownames(mat_z) <- gsub("^.*\\|\\s*", "", rownames(mat_z))

# --- Column annotation (groups) ---
meta_li <- meta[match(colnames(mat_z), meta$sample_id), , drop = FALSE]

# Define colors 
group_colors <- list(
  Group = c(
    control  = "#0072B2",  # strong blue
    hfid     = "#E69F00",  # vivid orange
    dss      = "#CC79A7",  # magenta/purple
    dss_hfid = "#009E73"   # green
  )
)

meta_li$group <- factor(meta_li$group, levels = names(group_colors$Group))

top_anno <- HeatmapAnnotation(
  Group = meta_li$group,
  col = list(Group = group_colors$Group)
)

# -------------------------
# now plot
# -------------------------

png(
  filename = "kegg_pathway_heatmap_z-score_imputed_categ_4-20_test.png",
  width = 3600,
  height = 5000,
  res = 300
)

ht <- Heatmap(
  mat_z,
  name = "Z-score",
  col = viridis::plasma (100),
  
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  column_title = "KEGG pathways - males",
  column_title_gp = gpar(
    fontface = "bold",
    fontsize = 16
  ),
  
  row_split = row_splits,
  
  top_annotation = top_anno,
  
  row_title_rot = 0,  
  
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 8)
)

draw(
  ht,
  heatmap_legend_side = "left",
  annotation_legend_side = "left",
  padding = unit(c(5, 5, 5, 40), "mm"),
  merge_legends = TRUE
)

dev.off()

############################################################################################################
# based on the same categories, make separate heatmaps for each broad category of KEGG pathways
############################################################################################################

#-----------------
# column clustered version first
#-----------------

library(openxlsx)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(grid)

# -----------------------------
# Load annotation file
# -----------------------------
df_labels <- read.xlsx("df_labels_for_manual_edit.xlsx") %>%
  as.data.frame(stringsAsFactors = FALSE)

colnames(df_labels)[1] <- "rowname"

# -----------------------------
# Match annotation to matrix
# -----------------------------
common_rows <- intersect(df_labels$rowname, rownames(mat_z))

mat_z <- mat_z[common_rows, , drop = FALSE]
df_labels <- df_labels[df_labels$rowname %in% common_rows, ]
df_labels <- df_labels[match(rownames(mat_z), df_labels$rowname), ]

# -----------------------------
# Define KEGG categories
# -----------------------------
row_splits <- df_labels$Level2
row_splits[is.na(row_splits)] <- "Other"

split_list <- split(seq_len(nrow(mat_z)), row_splits)

# -----------------------------
# Metadata + annotation
# -----------------------------
meta_li <- meta[match(colnames(mat_z), meta$sample_id), , drop = FALSE]

group_colors <- list(
  Group = c(
    control  = "#0072B2",
    hfid     = "#E69F00",
    dss      = "#CC79A7",
    dss_hfid = "#009E73"
  )
)

meta_li$group <- factor(meta_li$group, levels = names(group_colors$Group))

top_anno <- HeatmapAnnotation(
  Group = meta_li$group,
  col = list(Group = group_colors$Group)
)

# -----------------------------
# Output directory
# -----------------------------
dir.create("KEGG_category_heatmaps", showWarnings = FALSE)

# -----------------------------
# Loop over categories
# -----------------------------
for (cat in names(split_list)) {
  
  idx <- split_list[[cat]]
  mat_sub <- mat_z[idx, , drop = FALSE]
  
  if (nrow(mat_sub) < 3) next
  
  nrow_sub <- nrow(mat_sub)
  ncol_sub <- ncol(mat_sub)
  
  # heatmap geometry
  cell_size <- 4
  
  ht_width  <- unit(ncol_sub * cell_size, "mm")
  ht_height <- unit(nrow_sub * cell_size, "mm")
  
  # safe device sizing (prevents cropping)
  png(
    filename = paste0(
      "KEGG_category_heatmaps/heatmap_",
      gsub("[^A-Za-z0-9]", "_", cat),
      ".png"
    ),
    width = (ncol_sub * 70) + 1500,
    height = (nrow_sub * 70) + 1500,
    res = 300
  )
  
  ht <- Heatmap(
    mat_sub,
    name = "Z-score",
    col = viridis::plasma (100),
    
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = FALSE,
    
    top_annotation = top_anno,
    
    column_title = cat,
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    
    width = ht_width,
    height = ht_height,
    
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8)
  )
  
  draw(
    ht,
    heatmap_legend_side = "left",
    annotation_legend_side = "left"
  )
  
  dev.off()
}

############################################################################################################
# broad categories separate heatmaps - not column clustered
############################################################################################################

# -----------------------------
# Define group order + colors
# -----------------------------
group_levels <- c("control", "hfid", "dss", "dss_hfid")

group_colors <- list(
  Group = c(
    control  = "#0072B2",
    hfid     = "#E69F00",
    dss      = "#CC79A7",
    dss_hfid = "#009E73"
  )
)

# -----------------------------
# Align anno to mat_z FIRST
# -----------------------------
common_samples <- intersect(colnames(mat_z), rownames(anno))

mat_z <- mat_z[, common_samples]
anno  <- anno[common_samples, , drop = FALSE]

# strict check
stopifnot(all(colnames(mat_z) == rownames(anno)))

# -----------------------------
# Enforce group order
# -----------------------------
anno$Group <- factor(anno$Group, levels = group_levels)

# -----------------------------
# Order samples
# -----------------------------
ord <- order(anno$Group)

mat_z_ord <- mat_z[, ord]
anno_ord  <- anno[ord, , drop = FALSE]

# -----------------------------
# Create annotation
# -----------------------------
top_anno <- HeatmapAnnotation(
  Group = anno_ord$Group,
  col = group_colors
)

# -----------------------------
# Output directory
# -----------------------------
dir.create("KEGG_category_heatmaps_unclus", showWarnings = FALSE)

# -----------------------------
# Loop over KEGG categories
# -----------------------------
for (cat in names(split_list)) {
  
  idx <- split_list[[cat]]
  mat_sub <- mat_z_ord[idx, , drop = FALSE]
  
  if (nrow(mat_sub) < 3) next
  
  nrow_sub <- nrow(mat_sub)
  ncol_sub <- ncol(mat_sub)
  
  # geometry
  cell_size <- 4
  
  ht_width  <- unit(ncol_sub * cell_size, "mm")
  ht_height <- unit(nrow_sub * cell_size, "mm")
  
  png(
    filename = paste0(
      "KEGG_category_heatmaps_unclus/heatmap_",
      gsub("[^A-Za-z0-9]", "_", cat),
      ".png"
    ),
    width  = (ncol_sub * 80) + 1200,
    height = (nrow_sub * 80) + 1200,
    res = 300
  )
  
  ht <- Heatmap(
    mat_sub,
    name = "Z-score",
    col = viridis::plasma(100),
    
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    
    top_annotation = top_anno,
    
    column_title = cat,
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    
    width = ht_width,
    height = ht_height,
    
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8)
  )
  
  draw(
    ht,
    heatmap_legend_side = "left",
    annotation_legend_side = "left"
  )
  
  dev.off()
}

#############################################################################
# make KEGG pathway heatmap for only select pathways of interest, type kos
#############################################################################

library(dplyr)
library(pheatmap)
library(viridis)
library(tibble)

# Define pathways you WANT to show
kos_to_keep <- c(
  "ko01200",
  "ko00220",
  "ko00910",
  "ko00250",
  "ko00471"
  
)

df_kegg_rel_clean_filtered_select <- df_kegg_rel_clean_filtered %>%
  filter(KEGG_Pathway %in% kos_to_keep)

sample_cols <- colnames(df_kegg_rel_clean_filtered_select)[grepl("^LI", colnames(df_kegg_rel_clean_filtered_select))]

df_labels <- df_kegg_rel_clean_filtered_select %>%
  mutate(row_label = paste0(pathway_name, " (", KEGG_Pathway, ")"))

# Combine pathway name + KO ID for row labels
df_labels <- df_kegg_rel_clean_filtered_select %>%
  mutate(row_label = paste0(pathway_name, " (", KEGG_Pathway, ")"))

# Create numeric matrix
mat <- df_labels %>%
  select(row_label, all_of(sample_cols)) %>%
  column_to_rownames(var = "row_label") %>%
  as.data.frame()

# Force numeric conversion and replace Inf/NaN with NA
mat[] <- lapply(mat, as.numeric)
mat <- as.matrix(mat)
mat[!is.finite(mat)] <- NA

# Match metadata
meta_li <- meta[meta$sample_id %in% colnames(mat), , drop = FALSE]
meta_li <- meta_li[match(colnames(mat), meta_li$sample_id), , drop = FALSE]

# Annotation dataframe for pheatmap
anno <- data.frame(Group = meta_li$group)
rownames(anno) <- meta_li$sample_id

# Optional: assign custom colors for each group
group_colors <- list(
  Group = c(
    control  = "#0072B2",
    hfid     = "#E69F00",
    dss      = "#CC79A7",
    dss_hfid = "#009E73"
  )
)

#not enforce group order 
anno$Group <- as.character(anno$Group)

# Enforce desired group order
anno$Group <- factor(
  anno$Group,
  levels = c("control", "hfid", "dss", "dss_hfid")
)

# Reorder samples by group
ord <- order(anno$Group)
mat_ord <- mat[, ord]
anno_ord <- anno[ord, , drop = FALSE]


png(
  filename = "males_GDH_pathways_imp_clus.png",
  width = 3500,
  height = 3000,
  res = 300
)

pheatmap(
  mat_ord,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 0,
  treeheight_col = 8,
  annotation_col = anno_ord,
  annotation_colors = group_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  na_col = "black",
  color = viridis(100, option = "plasma"),
  main = "NADP-GDH pathways (males)",
  name = "Z-score",
  fontsize_row = 7.5,          # shrink row names
  fontsize_col = 7,          # shrink column names
  legend = TRUE,
  legend_labels = NULL,      
  cellwidth = 10,            # shrink cells horizontally
  cellheight = 10            # shrink cells vertically
)

dev.off()



####################################################################################################################################
# VOLCANO PLOTS OF PROTEINS IN SELECT KOs COLORED BY TAXA
####################################################################################################################################

####################################################################################################################################
# DSS-HFID VS HFID COLORED VOLCANO PLOT COMPARISON (ko00471|ko00460|ko00400|ko00410|ko00473|ko00450|ko00290|ko00340|ko00280|ko01230|ko00260|ko00330|ko00380|ko00300|ko00310|ko00250|ko00220|ko00350|ko00430|ko00270|ko00480|ko00360)
####################################################################################################################################

library(dplyr)
library(tibble)

crosstab <- read_csv("imputed_w_taxa_kegg.csv")

# Identify sample columns (all LI samples)
sample_cols <- colnames(crosstab)[grepl("^LI", colnames(crosstab))]

expr_matrix <- crosstab %>%
  select(Protein, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), ~ as.numeric(.))) %>%
  column_to_rownames("Protein") %>%
  as.matrix()

meta <- meta %>%
  filter(sample_id %in% colnames(expr_matrix)) %>%
  arrange(match(sample_id, colnames(expr_matrix)))

library(limma)

meta$group <- factor(meta$group, levels = c("hfid", "dss_hfid", "control", "dss"))

design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)

fit <- lmFit(expr_matrix, design)

contrast_matrix <- makeContrasts(
  dss_hfid_vs_hfid = dss_hfid - hfid,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

library(dplyr)
library(tibble)
library(limma)
library(ggplot2)

# -----------------------------
# limma result
# -----------------------------
res <- topTable(fit2, coef = "dss_hfid_vs_hfid", number = Inf)
res$Protein <- rownames(res)

# -----------------------------
# annotation join
# -----------------------------
annot <- crosstab %>%
  select(Protein, KEGG_Pathway, Taxa) %>%
  distinct() %>%
  mutate(
    Taxa = gsub("\\[|\\]", "", Taxa),
    Taxa = trimws(Taxa)
  )

# -----------------------------
# merge + filter ko of interest + clean stats
# -----------------------------
res_kegg <- res %>%
  left_join(annot, by = "Protein") %>%
  filter(grepl(
    "ko00471|ko00460|ko00400|ko00410|ko00473|ko00450|ko00290|ko00340|ko00280|ko01230|ko00260|ko00330|ko00380|ko00300|ko00310|ko00250|ko00220|ko00350|ko00430|ko00270|ko00480|ko00360",
    KEGG_Pathway
  )) %>%
  mutate(
    logFC = as.numeric(logFC),
    adj.P.Val = as.numeric(adj.P.Val),
    negLog10P = -log10(adj.P.Val),
    
    sig = ifelse(adj.P.Val < 0.05 & abs(logFC) >= 2, "sig", "ns"),
    
    # IMPORTANT: correct biological interpretation
    direction = case_when(
      logFC >=  2 ~ "DSS_HFID",  # right side
      logFC <= -2 ~ "HFID",      # left side
      TRUE ~ "NS"
    )
  ) %>%
  mutate(
    Taxa = gsub("\\[|\\]", "", Taxa),
    Taxa = trimws(Taxa),
    Taxa = ifelse(is.na(Taxa) | Taxa == "" | Taxa == "NA", NA, Taxa)
  ) %>%
  filter(!is.na(Taxa))

# -----------------------------
# select representative taxa
# -----------------------------
dss_taxa <- res_kegg %>%
  filter(sig == "sig", direction == "DSS_HFID") %>%
  arrange(desc(logFC)) %>%
  distinct(Taxa, .keep_all = TRUE) %>%
  slice_head(n = 10) %>%
  pull(Taxa)

hfid_taxa <- res_kegg %>%
  filter(sig == "sig", direction == "HFID") %>%
  arrange(logFC) %>%
  distinct(Taxa, .keep_all = TRUE) %>%
  slice_head(n = 10) %>%
  pull(Taxa)

# -----------------------------
# assign plotting groups
# -----------------------------
res_kegg <- res_kegg %>%
  mutate(
    Taxa_group = case_when(
      sig == "sig" & Taxa %in% dss_taxa ~ paste0("DSS_", Taxa),
      sig == "sig" & Taxa %in% hfid_taxa ~ paste0("HFID_", Taxa),
      TRUE ~ "NS"
    )
  )

# -----------------------------
# color palettes
# -----------------------------
dss_colors <- setNames(
  c("#ffffb2","#fed976","#feb24c","#fd8d3c","#fc4e2a",
    "#e31a1c","#bd0026","#f768a1","#fbb4b9","#f768a1")[seq_along(dss_taxa)],
  paste0("DSS_", dss_taxa)
)

hfid_colors <- setNames(
  c("#c7e9c0","#a1d99b","#74c476","#41ab5d","#238b45",
    "#66c2a4","#41b6c4","#2b8cbe","#253494","#54278f")[seq_along(hfid_taxa)],
  paste0("HFID_", hfid_taxa)
)

palette_all <- c(dss_colors, hfid_colors, NS = "grey85")

# -----------------------------
# counts
# -----------------------------
dss_count <- sum(res_kegg$direction == "DSS_HFID" & res_kegg$sig == "sig", na.rm = TRUE)
hfid_count <- sum(res_kegg$direction == "HFID" & res_kegg$sig == "sig", na.rm = TRUE)

x_min <- min(res_kegg$logFC, na.rm = TRUE)
x_max <- max(res_kegg$logFC, na.rm = TRUE)
y_max <- max(res_kegg$negLog10P, na.rm = TRUE)

# -----------------------------
# generate the volcano plot
# -----------------------------
# -----------------------------
# desired legend order
# -----------------------------
legend_order <- c(
  names(hfid_colors),   # blues/greens first
  "NS",
  names(dss_colors)     # reds/oranges second
)

# -----------------------------
# volcano plot
# -----------------------------
p_volcano <- ggplot(res_kegg, aes(x = logFC, y = negLog10P)) +
  
  # non-significant
  geom_point(
    data = subset(res_kegg, sig == "ns"),
    color = "grey85",
    size = 1.5
  ) +
  
  # significant
  geom_point(
    data = subset(res_kegg, sig == "sig"),
    aes(color = Taxa_group),
    size = 4,
    alpha = 0.85
  ) +
  
  scale_color_manual(
    values = palette_all,
    
    # FORCE LEGEND ORDER
    breaks = legend_order,
    
    name = "Taxa",
    
    labels = function(x) {
      sapply(x, function(y) {
        if (y == "NS") return("NS")
        
        taxa_name <- sub("^(DSS_|HFID_)", "", y)
        
        bquote(italic(.(taxa_name)))
      })
    },
    
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE
    )
  ) +
  
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  annotate(
    "text",
    x = x_max,
    y = y_max,
    label = paste0("DSS_HFID: ", dss_count),
    hjust = 1,
    vjust = 1,
    fontface = "bold",
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_min,
    y = y_max,
    label = paste0("HFID: ", hfid_count),
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 4
  ) +
  
  theme_classic(base_size = 13) +
  
  labs(
    title = "Volcano Plot: DSS_HFID vs HFID (males)",
    x = "log2 Fold Change",
    y = "-log10 adjusted P-value"
  )

# -----------------------------
# save plot
# -----------------------------
ggsave(
  "volcano_all_amino_taxa_colored_final.png",
  p_volcano,
  width = 9,
  height = 6,
  dpi = 600
)

# -----------------------------
# exporting data and stats
# -----------------------------

library(openxlsx)
library(dplyr)

# -----------------------------
# Build export components
# -----------------------------

# limma results
limma_all <- res %>%
  mutate(
    negLog10P = -log10(adj.P.Val)
  )

# Annotated full results
limma_annotated <- res %>%
  left_join(annot, by = "Protein")

# KEGG-filtered results 
limma_kegg <- res_kegg

# Summary stats
summary_stats <- limma_kegg %>%
  summarise(
    n_total = n(),
    n_sig = sum(sig == "sig"),
    n_up = sum(direction == "Up"),
    n_down = sum(direction == "Down"),
    n_ns = sum(direction == "NS")
  )

# Taxa-level summary
taxa_summary <- limma_kegg %>%
  group_by(Taxa, direction) %>%
  summarise(
    n = n(),
    mean_logFC = mean(logFC, na.rm = TRUE),
    min_padj = min(adj.P.Val, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(direction, desc(n))

# Top hits (for reporting)
top_hits <- limma_kegg %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50)

# -----------------------------
# Write Excel workbook
# -----------------------------
wb <- createWorkbook()

addWorksheet(wb, "limma_all")
writeData(wb, "limma_all", limma_all)

addWorksheet(wb, "limma_annotated")
writeData(wb, "limma_annotated", limma_annotated)

addWorksheet(wb, "limma_kegg_filtered")
writeData(wb, "limma_kegg_filtered", limma_kegg)

addWorksheet(wb, "summary_stats")
writeData(wb, "summary_stats", summary_stats)

addWorksheet(wb, "taxa_summary")
writeData(wb, "taxa_summary", taxa_summary)

addWorksheet(wb, "top_hits")
writeData(wb, "top_hits", top_hits)

saveWorkbook(
  wb,
  file = "limma_volcano_amino_acid_dss_hfid_vs_hfid_males_STAT_OUTPUT.xlsx",
  overwrite = TRUE
)

####################################################################################################################################
####NOW MAKING BOX PLOT OF total INTENSITY OF ALL PROTEINS BELONGING TO AMINO ACID RELEVANT PATHWAYS
####################################################################################################################################
library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(ggplot2)
library(readr)
library(ggsignif)

crosstab <- read_csv("imputed_w_taxa_kegg.csv")

# -----------------------------
# 1. Filter data & calculate sample means for STATISTICS
# -----------------------------
target_kegg <- "ko00471|ko00460|ko00400|ko00410|ko00473|ko00450|ko00290|ko00340|ko00280|ko01230|ko00260|ko00330|ko00380|ko00300|ko00310|ko00250|ko00220|ko00350|ko00430|ko00270|ko00480|ko00360"

sample_cols <- colnames(crosstab)[grepl("^LI", colnames(crosstab))]

# Convert normalized log2 intensities back to linear intensities
crosstab_linear <- crosstab %>%
  mutate(across(all_of(sample_cols), ~ 2^.))

# Calculate relative abundance within each sample
crosstab_rel <- crosstab_linear %>%
  mutate(
    across(
      all_of(sample_cols),
      ~ . / sum(., na.rm = TRUE)
    )
  )

# Sum relative abundances for proteins in amino acid pathways
summed_intensities <- crosstab_rel %>%
  filter(grepl(target_kegg, KEGG_Pathway)) %>%
  select(all_of(sample_cols)) %>%
  colSums(na.rm = TRUE)

plot_df <- data.frame(
  sample_id = names(summed_intensities),
  total_intensity = summed_intensities
) %>%
  left_join(meta, by = "sample_id") %>%
  filter(!is.na(group))

plot_df$group <- factor(plot_df$group, levels = c("control", "hfid", "dss", "dss_hfid"))

# -----------------------------
# 2. Run limma stats on sample SUMS
# -----------------------------
expr_sum_matrix <- matrix(plot_df$total_intensity, nrow = 1)
colnames(expr_sum_matrix) <- plot_df$sample_id

design_sum <- model.matrix(~ 0 + group, data = plot_df)
colnames(design_sum) <- levels(plot_df$group)

fit_sum <- lmFit(expr_sum_matrix, design_sum)

contrast_matrix_sum <- makeContrasts(
  dss_hfid_vs_hfid = dss_hfid - hfid,
  dss_hfid_vs_control = dss_hfid - control,
  dss_vs_control = dss - control,
  dss_hfid_vs_dss = dss_hfid - dss,      
  levels = design_sum
)

fit2_sum <- contrasts.fit(fit_sum, contrast_matrix_sum)
fit2_sum <- eBayes(fit2_sum)

pvals <- fit2_sum$p.value[1, ]

# -----------------------------
# 6. SAVE STATS TO EXCEL
# -----------------------------

# Raw summed intensity for every sample
sample_values <- plot_df %>%
  arrange(group) %>%
  select(sample_id, group, total_intensity)

# Summary statistics
group_summary <- plot_df %>%
  group_by(group) %>%
  summarise(
    n = n(),
    Mean = mean(total_intensity),
    SD = sd(total_intensity),
    SEM = SD/sqrt(n()),
    Median = median(total_intensity),
    Min = min(total_intensity),
    Max = max(total_intensity),
    .groups = "drop"
  )

# Limma statistics
coef_names <- colnames(contrast_matrix_sum)

stats_list <- lapply(coef_names, function(coef){
  
  res <- topTable(fit2_sum,
                  coef = coef,
                  number = Inf)
  
  res$Contrast <- coef
  
  res %>%
    select(
      Contrast,
      logFC,
      AveExpr,
      t,
      P.Value,
      adj.P.Val,
      B
    )
  
})

limma_stats <- bind_rows(stats_list)

# Write workbook
write_xlsx(
  list(
    Sample_Values = sample_values,
    Group_Summary = group_summary,
    Limma_Statistics = limma_stats
  ),
  path = "aa_pathway_stats_imp.xlsx"
)

# -----------------------------
# 3. Update factor levels (No pivot_longer needed!)
# -----------------------------
plot_df$group <- factor(
  plot_df$group, 
  levels = c("control", "hfid", "dss", "dss_hfid"),
  labels = c("Control", "HFiD", "DSS", "DSS+HFiD")
)

# -----------------------------
# 1. Update the p-value formatting function to ONLY return asterisks
# -----------------------------
format_pval <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("")
}

# -----------------------------
# 2. Re-filter for significant comparisons
# -----------------------------
# (Your existing code to calculate sig_pairs and sig_annotations remains the same)
# The sapply(all_pvals[sig_indices], format_pval) will now return only the symbols.

# -----------------------------
# 3. Generate the box plot
# -----------------------------
p_box <- ggplot(plot_df, aes(x = group, y = total_intensity, fill = group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8, color = "black") +
  scale_fill_manual(values = c(
    "Control"  = "#0072B2", 
    "HFiD"     = "#E69F00", 
    "DSS"      = "#CC79A7", 
    "DSS+HFiD" = "#009E73"
  ))

# -----------------------------
# Create significance annotations
# -----------------------------

all_pvals <- fit2_sum$p.value[1, ]

comparison_names <- list(
  c("HFiD", "DSS+HFiD"),
  c("Control", "DSS+HFiD"),
  c("Control", "DSS"),
  c("DSS", "DSS+HFiD")
)

format_pval <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else ""
}

sig_idx <- which(all_pvals < 0.05)

sig_pairs <- comparison_names[sig_idx]
sig_annotations <- sapply(all_pvals[sig_idx], format_pval)

# Add brackets with asterisk-only annotations
if (length(sig_pairs) > 0) {
  p_box <- p_box + 
    geom_signif(
      comparisons = sig_pairs,
      annotations = sig_annotations, # Now contains only "***", "**", or "*"
      step_increase = 0.12,  
      margin_top = 0.05,
      textsize = 5,                  # Slightly larger text for visibility
      vjust = -0.2
    )
}

p_box <- p_box + 
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Arginine and glutamate metabolizing proteins",
    x = "Experimental Group",
    y = "Protein relative abundance"
  )

print(p_box)

ggsave(
  "arginine_glutamate_pathway_proteins_boxplot_sum_imp.png",
  p_box,
  width=5,
  height=5.5,
  dpi=600
)

####################################################################################################################################
####NOW MAKING BOX PLOT OF total INTENSITY OF ALL PROTEINS BELONGING TO fructose and galactose KEGG pathways 
####################################################################################################################################
crosstab <- read_csv("imputed_w_taxa_kegg.csv")

library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(ggplot2)
library(readr)
library(ggsignif)
library(writexl) # Required for Excel export

# -----------------------------
# 1. Filter data & calculate sample SUMS
# -----------------------------
target_kegg <- "\\b(ko00051|ko00052)\\b"

sample_cols <- colnames(crosstab)[grepl("^LI", colnames(crosstab))]

# Convert normalized log2 intensities back to linear intensities
crosstab_linear <- crosstab %>%
  mutate(across(all_of(sample_cols), ~ 2^.))

# Calculate relative abundance within each sample
crosstab_rel <- crosstab_linear %>%
  mutate(
    across(
      all_of(sample_cols),
      ~ . / sum(., na.rm = TRUE)
    )
  )

# Sum relative abundances for proteins in amino acid pathways
summed_intensities <- crosstab_rel %>%
  filter(grepl(target_kegg, KEGG_Pathway)) %>%
  select(all_of(sample_cols)) %>%
  colSums(na.rm = TRUE)

plot_df <- data.frame(
  sample_id = names(summed_intensities),
  total_intensity = summed_intensities
) %>%
  left_join(meta, by = "sample_id") %>%
  filter(!is.na(group))

plot_df$group <- factor(plot_df$group, levels = c("control", "hfid", "dss", "dss_hfid"))

# -----------------------------
# 2. Run limma stats on sample SUMS
# -----------------------------
expr_sum_matrix <- matrix(plot_df$total_intensity, nrow = 1)
colnames(expr_sum_matrix) <- plot_df$sample_id

design_sum <- model.matrix(~ 0 + group, data = plot_df)
colnames(design_sum) <- levels(plot_df$group)

fit_sum <- lmFit(expr_sum_matrix, design_sum)

contrast_matrix_sum <- makeContrasts(
  dss_hfid_vs_hfid = dss_hfid - hfid,
  dss_hfid_vs_control = dss_hfid - control,
  dss_vs_control = dss - control,
  dss_hfid_vs_dss = dss_hfid - dss,      
  levels = design_sum
)

fit2_sum <- contrasts.fit(fit_sum, contrast_matrix_sum)
fit2_sum <- eBayes(fit2_sum)

pvals <- fit2_sum$p.value[1, ]

# -----------------------------
# 6. SAVE STATS TO EXCEL
# -----------------------------

# Raw summed intensity for every sample
sample_values <- plot_df %>%
  arrange(group) %>%
  select(sample_id, group, total_intensity)

# Summary statistics
group_summary <- plot_df %>%
  group_by(group) %>%
  summarise(
    n = n(),
    Mean = mean(total_intensity),
    SD = sd(total_intensity),
    SEM = SD/sqrt(n()),
    Median = median(total_intensity),
    Min = min(total_intensity),
    Max = max(total_intensity),
    .groups = "drop"
  )

# Limma statistics
coef_names <- colnames(contrast_matrix_sum)

stats_list <- lapply(coef_names, function(coef){
  
  res <- topTable(fit2_sum,
                  coef = coef,
                  number = Inf)
  
  res$Contrast <- coef
  
  res %>%
    select(
      Contrast,
      logFC,
      AveExpr,
      t,
      P.Value,
      adj.P.Val,
      B
    )
  
})

limma_stats <- bind_rows(stats_list)

# Write workbook
write_xlsx(
  list(
    Sample_Values = sample_values,
    Group_Summary = group_summary,
    Limma_Statistics = limma_stats
  ),
  path = "51_52_pathway_stats_imp.xlsx"
)

# -----------------------------
# 3. Update factor levels (No pivot_longer needed!)
# -----------------------------
plot_df$group <- factor(
  plot_df$group, 
  levels = c("control", "hfid", "dss", "dss_hfid"),
  labels = c("Control", "HFiD", "DSS", "DSS+HFiD")
)

# -----------------------------
# 1. Update the p-value formatting function to ONLY return asterisks
# -----------------------------
format_pval <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("")
}


# -----------------------------
# 3. Generate the box plot
# -----------------------------
p_box <- ggplot(plot_df, aes(x = group, y = total_intensity, fill = group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8, color = "black") +
  scale_fill_manual(values = c(
    "Control"  = "#0072B2", 
    "HFiD"     = "#E69F00", 
    "DSS"      = "#CC79A7", 
    "DSS+HFiD" = "#009E73"
  ))

# -----------------------------
# Create significance annotations
# -----------------------------

all_pvals <- fit2_sum$p.value[1, ]

comparison_names <- list(
  c("HFiD", "DSS+HFiD"),
  c("Control", "DSS+HFiD"),
  c("Control", "DSS"),
  c("DSS", "DSS+HFiD")
)

format_pval <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else ""
}

sig_idx <- which(all_pvals < 0.05)

sig_pairs <- comparison_names[sig_idx]
sig_annotations <- sapply(all_pvals[sig_idx], format_pval)

# Add brackets with asterisk-only annotations
if (length(sig_pairs) > 0) {
  p_box <- p_box + 
    geom_signif(
      comparisons = sig_pairs,
      annotations = sig_annotations, # Now contains only "***", "**", or "*"
      step_increase = 0.12,  
      margin_top = 0.05,
      textsize = 5,                  # Slightly larger text for visibility
      vjust = -0.2
    )
}

p_box <- p_box + 
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "ko00051 and ko00052 metabolizing proteins",
    x = "Experimental Group",
    y = "Protein relative abundance"
  )

print(p_box)

ggsave(
  "51_52_pathway_proteins_boxplot_sum_imp.png",
  p_box,
  width=5,
  height=5.5,
  dpi=600
)

####################################################################################################################################
####NOW MAKING BOX PLOT OF total INTENSITY OF ALL PROTEINS BELONGING TO all carb relevant pathways 
####################################################################################################################################
library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(ggplot2)
library(readr)
library(ggsignif)

crosstab <- read_csv("imputed_w_taxa_kegg.csv")

# -----------------------------
# 1. Filter data & calculate sample means for STATISTICS
# -----------------------------
target_kegg <- "ko00630|ko00650|ko00640|ko00500|ko01210|ko00053|ko00660|ko00660|ko00562|ko00020|ko00052|ko00010|ko00051|ko00030|ko00040|ko00620"

sample_cols <- colnames(crosstab)[grepl("^LI", colnames(crosstab))]

# Convert normalized log2 intensities back to linear intensities
crosstab_linear <- crosstab %>%
  mutate(across(all_of(sample_cols), ~ 2^.))

# Calculate relative abundance within each sample
crosstab_rel <- crosstab_linear %>%
  mutate(
    across(
      all_of(sample_cols),
      ~ . / sum(., na.rm = TRUE)
    )
  )

# Sum relative abundances for proteins in amino acid pathways
summed_intensities <- crosstab_rel %>%
  filter(grepl(target_kegg, KEGG_Pathway)) %>%
  select(all_of(sample_cols)) %>%
  colSums(na.rm = TRUE)

plot_df <- data.frame(
  sample_id = names(summed_intensities),
  total_intensity = summed_intensities
) %>%
  left_join(meta, by = "sample_id") %>%
  filter(!is.na(group))

plot_df$group <- factor(plot_df$group, levels = c("control", "hfid", "dss", "dss_hfid"))

# -----------------------------
# 2. Run limma stats on sample SUMS
# -----------------------------
expr_sum_matrix <- matrix(plot_df$total_intensity, nrow = 1)
colnames(expr_sum_matrix) <- plot_df$sample_id

design_sum <- model.matrix(~ 0 + group, data = plot_df)
colnames(design_sum) <- levels(plot_df$group)

fit_sum <- lmFit(expr_sum_matrix, design_sum)

contrast_matrix_sum <- makeContrasts(
  dss_hfid_vs_hfid = dss_hfid - hfid,
  dss_hfid_vs_control = dss_hfid - control,
  dss_vs_control = dss - control,
  dss_hfid_vs_dss = dss_hfid - dss,      
  levels = design_sum
)

fit2_sum <- contrasts.fit(fit_sum, contrast_matrix_sum)
fit2_sum <- eBayes(fit2_sum)

pvals <- fit2_sum$p.value[1, ]

# -----------------------------
# 6. SAVE STATS TO EXCEL
# -----------------------------

# Raw summed intensity for every sample
sample_values <- plot_df %>%
  arrange(group) %>%
  select(sample_id, group, total_intensity)

# Summary statistics
group_summary <- plot_df %>%
  group_by(group) %>%
  summarise(
    n = n(),
    Mean = mean(total_intensity),
    SD = sd(total_intensity),
    SEM = SD/sqrt(n()),
    Median = median(total_intensity),
    Min = min(total_intensity),
    Max = max(total_intensity),
    .groups = "drop"
  )

# Limma statistics
coef_names <- colnames(contrast_matrix_sum)

stats_list <- lapply(coef_names, function(coef){
  
  res <- topTable(fit2_sum,
                  coef = coef,
                  number = Inf)
  
  res$Contrast <- coef
  
  res %>%
    select(
      Contrast,
      logFC,
      AveExpr,
      t,
      P.Value,
      adj.P.Val,
      B
    )
  
})

limma_stats <- bind_rows(stats_list)

# Write workbook
write_xlsx(
  list(
    Sample_Values = sample_values,
    Group_Summary = group_summary,
    Limma_Statistics = limma_stats
  ),
  path = "all_carb_pathway_stats_imp.xlsx"
)

# -----------------------------
# 3. Update factor levels (No pivot_longer needed!)
# -----------------------------
plot_df$group <- factor(
  plot_df$group, 
  levels = c("control", "hfid", "dss", "dss_hfid"),
  labels = c("Control", "HFiD", "DSS", "DSS+HFiD")
)

# -----------------------------
# 1. Update the p-value formatting function to ONLY return asterisks
# -----------------------------
format_pval <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("")
}

# -----------------------------
# 3. Generate the box plot
# -----------------------------
p_box <- ggplot(plot_df, aes(x = group, y = total_intensity, fill = group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8, color = "black") +
  scale_fill_manual(values = c(
    "Control"  = "#0072B2", 
    "HFiD"     = "#E69F00", 
    "DSS"      = "#CC79A7", 
    "DSS+HFiD" = "#009E73"
  ))

# -----------------------------
# Create significance annotations
# -----------------------------

all_pvals <- fit2_sum$p.value[1, ]

comparison_names <- list(
  c("HFiD", "DSS+HFiD"),
  c("Control", "DSS+HFiD"),
  c("Control", "DSS"),
  c("DSS", "DSS+HFiD")
)

format_pval <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else ""
}

sig_idx <- which(all_pvals < 0.05)

sig_pairs <- comparison_names[sig_idx]
sig_annotations <- sapply(all_pvals[sig_idx], format_pval)

# Add brackets with asterisk-only annotations
if (length(sig_pairs) > 0) {
  p_box <- p_box + 
    geom_signif(
      comparisons = sig_pairs,
      annotations = sig_annotations, # Now contains only "***", "**", or "*"
      step_increase = 0.12,  
      margin_top = 0.05,
      textsize = 5,                  # Slightly larger text for visibility
      vjust = -0.2
    )
}

p_box <- p_box + 
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "All carb pathways metabolizing proteins",
    x = "Experimental Group",
    y = "Protein relative abundance"
  )

print(p_box)

ggsave(
  "all_carb_pathway_proteins_boxplot_sum_imp.png",
  p_box,
  width=5,
  height=5.5,
  dpi=600
)

####################################################################################################################################
# DSS+HFID VS DSS TAXA COLORED VOLCANO PLOT COMPARISON (ko00471|ko00460|ko00400|ko00410|ko00473|ko00450|ko00290|ko00340|ko00280|ko01230|ko00260|ko00330|ko00380|ko00300|ko00310|ko00250|ko00220|ko00350|ko00430|ko00270|ko00480|ko00360)
####################################################################################################################################
library(dplyr)
library(tibble)
library(limma)
library(ggplot2)
library(openxlsx)

# -----------------------------
# Load data
# -----------------------------
crosstab <- read_csv("imputed_w_taxa_kegg.csv")

sample_cols <- colnames(crosstab)[grepl("^LI", colnames(crosstab))]

expr_matrix <- crosstab %>%
  select(Protein, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), ~ as.numeric(.))) %>%
  column_to_rownames("Protein") %>%
  as.matrix()

meta <- meta %>%
  filter(sample_id %in% colnames(expr_matrix)) %>%
  arrange(match(sample_id, colnames(expr_matrix)))

# -----------------------------
# Design matrix
# -----------------------------
meta$group <- factor(meta$group, levels = c("hfid", "dss_hfid", "control", "dss"))

design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)

fit <- lmFit(expr_matrix, design)

# -----------------------------
# define contrast
# DSS_HFID vs DSS
# -----------------------------
contrast_matrix <- makeContrasts(
  dss_hfid_vs_dss = dss_hfid - dss,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# -----------------------------
# limma results
# -----------------------------
res <- topTable(fit2, coef = "dss_hfid_vs_dss", number = Inf)
res$Protein <- rownames(res)

annot <- crosstab %>%
  select(Protein, KEGG_Pathway, Taxa) %>%
  distinct() %>%
  mutate(
    Taxa = gsub("\\[|\\]", "", Taxa),
    Taxa = trimws(Taxa)
  )

# -----------------------------
# merge + filter
# -----------------------------
res_kegg <- res %>%
  left_join(annot, by = "Protein") %>%
  filter(grepl(
    "ko00471|ko00460|ko00400|ko00410|ko00473|ko00450|ko00290|ko00340|ko00280|ko01230|ko00260|ko00330|ko00380|ko00300|ko00310|ko00250|ko00220|ko00350|ko00430|ko00270|ko00480|ko00360",
    KEGG_Pathway
  )) %>%
  mutate(
    logFC = as.numeric(logFC),
    adj.P.Val = as.numeric(adj.P.Val),
    negLog10P = -log10(adj.P.Val),
    
    sig = ifelse(adj.P.Val < 0.05 & abs(logFC) >= 2, "sig", "ns"),
    
    # CORE INTERPRETATION
    direction = case_when(
      logFC >=  2 ~ "DSS_HFID",  # right
      logFC <= -2 ~ "DSS",       # left (UPDATED)
      TRUE ~ "NS"
    )
  ) %>%
  mutate(
    Taxa = gsub("\\[|\\]", "", Taxa),
    Taxa = trimws(Taxa),
    Taxa = ifelse(is.na(Taxa) | Taxa == "" | Taxa == "NA", NA, Taxa)
  ) %>%
  filter(!is.na(Taxa))

# -----------------------------
# representative taxa
# -----------------------------
dss_hfid_taxa <- res_kegg %>%
  filter(sig == "sig", direction == "DSS_HFID") %>%
  arrange(desc(logFC)) %>%
  distinct(Taxa, .keep_all = TRUE) %>%
  slice_head(n = 10) %>%
  pull(Taxa)

dss_taxa <- res_kegg %>%
  filter(sig == "sig", direction == "DSS") %>%
  arrange(logFC) %>%
  distinct(Taxa, .keep_all = TRUE) %>%
  slice_head(n = 10) %>%
  pull(Taxa)

# -----------------------------
# plotting groups
# -----------------------------
res_kegg <- res_kegg %>%
  mutate(
    Taxa_group = case_when(
      sig == "sig" & Taxa %in% dss_hfid_taxa ~ paste0("DSS_HFID_", Taxa),
      sig == "sig" & Taxa %in% dss_taxa ~ paste0("DSS_", Taxa),
      TRUE ~ "NS"
    )
  )

# -----------------------------
# colors
# -----------------------------
dss_hfid_colors <- setNames(
  c("#ffffb2","#fed976","#feb24c","#fd8d3c","#fc4e2a",
    "#e31a1c","#bd0026","#f768a1","#fbb4b9","#f768a1")[seq_along(dss_hfid_taxa)],
  paste0("DSS_HFID_", dss_hfid_taxa)
)

dss_colors <- setNames(
  c("#c7e9c0","#a1d99b","#74c476","#41ab5d","#238b45",
    "#66c2a4","#41b6c4","#2b8cbe","#253494","#54278f")[seq_along(dss_taxa)],
  paste0("DSS_", dss_taxa)
)

palette_all <- c(dss_hfid_colors, dss_colors, NS = "grey85")

# -----------------------------
# counts
# -----------------------------
dss_hfid_count <- sum(res_kegg$direction == "DSS_HFID" & res_kegg$sig == "sig", na.rm = TRUE)
dss_count <- sum(res_kegg$direction == "DSS" & res_kegg$sig == "sig", na.rm = TRUE)

x_min <- min(res_kegg$logFC, na.rm = TRUE)
x_max <- max(res_kegg$logFC, na.rm = TRUE)
y_max <- max(res_kegg$negLog10P, na.rm = TRUE)

# -----------------------------
# volcano plot
# -----------------------------
# -----------------------------
# desired legend order
# DSS (blue) first
# DSS_HFID (red) second
# -----------------------------
legend_order <- c(
  names(dss_colors),
  "NS",
  names(dss_hfid_colors)
)

# -----------------------------
# volcano plot
# -----------------------------
p_volcano <- ggplot(
  res_kegg,
  aes(x = logFC, y = negLog10P)
) +
  
  # non-significant
  geom_point(
    data = subset(res_kegg, sig=="ns"),
    color = "grey85",
    size = 1.5
  ) +
  
  # significant
  geom_point(
    data = subset(res_kegg, sig=="sig"),
    aes(color = Taxa_group),
    size = 4,
    alpha = .85
  ) +
  
  scale_color_manual(
    values = palette_all,
    
    # force grouped legend
    breaks = legend_order,
    
    name = "Taxa",
    
    labels = function(x){
      
      sapply(
        x,
        function(y){
          
          if(y=="NS") return("NS")
          
          taxa_name <- y %>%
            gsub("^DSS_HFID_","",.) %>%
            gsub("^DSS_","",.)
          
          bquote(italic(.(taxa_name)))
        }
      )
      
    },
    
    guide=guide_legend(
      ncol=1,
      byrow=TRUE
    )
    
  ) +
  
  geom_vline(
    xintercept=c(-2,2),
    linetype="dashed"
  ) +
  
  geom_hline(
    yintercept=-log10(.05),
    linetype="dashed"
  ) +
  
  # RIGHT side
  annotate(
    "text",
    x=x_max,
    y=y_max,
    label=paste0(
      "DSS_HFID: ",
      dss_hfid_count
    ),
    hjust=1,
    vjust=1,
    fontface="bold",
    size=4
  ) +
  
  # LEFT side
  annotate(
    "text",
    x=x_min,
    y=y_max,
    label=paste0(
      "DSS: ",
      dss_count
    ),
    hjust=0,
    vjust=1,
    fontface="bold",
    size=4
  ) +
  
  theme_classic(base_size=13) +
  
  labs(
    title="Volcano Plot: DSS_HFID vs DSS (males)",
    x="log2 Fold Change",
    y="-log10 adjusted P-value"
  )

ggsave(
  "volcano_DSS_HFID_vs_DSS.png",
  p_volcano,
  width=10,
  height=6,
  dpi=600
)


# -----------------------------
# generate stats outputs
# -----------------------------
library(openxlsx)
library(dplyr)

# -----------------------------
# Full limma results
# -----------------------------
limma_all <- res %>%
  mutate(negLog10P = -log10(adj.P.Val))

# -----------------------------
# Annotated full results
# -----------------------------
limma_annotated <- res %>%
  left_join(annot, by = "Protein")

# -----------------------------
# KEGG-filtered (plot dataset)
# -----------------------------
limma_kegg <- res_kegg

# -----------------------------
# Summary statistics (IMPORTANT FIXED LOGIC)
# -----------------------------
summary_stats <- limma_kegg %>%
  summarise(
    n_total = n(),
    n_sig = sum(sig == "sig", na.rm = TRUE),
    n_dss_hfid = sum(direction == "DSS_HFID" & sig == "sig", na.rm = TRUE),
    n_dss = sum(direction == "DSS" & sig == "sig", na.rm = TRUE),
    n_ns = sum(sig == "ns", na.rm = TRUE)
  )

# -----------------------------
# Taxa-level summary
# -----------------------------
taxa_summary <- limma_kegg %>%
  group_by(Taxa, direction) %>%
  summarise(
    n = n(),
    mean_logFC = mean(logFC, na.rm = TRUE),
    min_padj = min(adj.P.Val, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(direction, desc(n))

# -----------------------------
# Top hits (for reporting)
# -----------------------------
top_hits <- limma_kegg %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50)

# -----------------------------
# Write Excel workbook
# -----------------------------
wb <- createWorkbook()

addWorksheet(wb, "limma_all")
writeData(wb, "limma_all", limma_all)

addWorksheet(wb, "limma_annotated")
writeData(wb, "limma_annotated", limma_annotated)

addWorksheet(wb, "limma_kegg_filtered")
writeData(wb, "limma_kegg_filtered", limma_kegg)

addWorksheet(wb, "summary_stats")
writeData(wb, "summary_stats", summary_stats)

addWorksheet(wb, "taxa_summary")
writeData(wb, "taxa_summary", taxa_summary)

addWorksheet(wb, "top_hits")
writeData(wb, "top_hits", top_hits)

saveWorkbook(
  wb,
  file = "limma_DSS_HFID_vs_DSS_males_STATS.xlsx",
  overwrite = TRUE
)


####################################################################################################################################
# VOLCANO PLOTS COMPARING CARB KOs (ko00051|ko00052)
####################################################################################################################################

####################################################################################################################################
#DSS-HFID VS HFID comparison
####################################################################################################################################

library(dplyr)
library(tibble)

crosstab <- read_csv("imputed_w_taxa_kegg.csv")


# Identify sample columns (all LI samples)
sample_cols <- colnames(crosstab)[grepl("^LI", colnames(crosstab))]

expr_matrix <- crosstab %>%
  select(Protein, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), ~ as.numeric(.))) %>%
  column_to_rownames("Protein") %>%
  as.matrix()

meta <- meta %>%
  filter(sample_id %in% colnames(expr_matrix)) %>%
  arrange(match(sample_id, colnames(expr_matrix)))

library(limma)

meta$group <- factor(meta$group, levels = c("hfid", "dss_hfid", "control", "dss"))

design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)

fit <- lmFit(expr_matrix, design)

contrast_matrix <- makeContrasts(
  dss_hfid_vs_hfid = dss_hfid - hfid,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

library(dplyr)
library(tibble)
library(limma)
library(ggplot2)

# -----------------------------
# limma result
# -----------------------------
res <- topTable(fit2, coef = "dss_hfid_vs_hfid", number = Inf)
res$Protein <- rownames(res)

# -----------------------------
# annotation join
# -----------------------------
annot <- crosstab %>%
  select(Protein, KEGG_Pathway, Taxa) %>%
  distinct() %>%
  mutate(
    Taxa = gsub("\\[|\\]", "", Taxa),
    Taxa = trimws(Taxa)
  )

# -----------------------------
# merge + filter KEGG + clean stats
# -----------------------------
res_kegg <- res %>%
  left_join(annot, by = "Protein") %>%
  filter(grepl(
    "ko00630|ko00650|ko00640|ko00500|ko01210|ko00053|ko00660|ko00660|ko00562|ko00020|ko00052|ko00010|ko00051|ko00030|ko00040|ko00620",
    KEGG_Pathway
  )) %>%
  mutate(
    logFC = as.numeric(logFC),
    adj.P.Val = as.numeric(adj.P.Val),
    negLog10P = -log10(adj.P.Val),
    
    sig = ifelse(adj.P.Val < 0.05 & abs(logFC) >= 2, "sig", "ns"),
    
    # IMPORTANT: correct biological interpretation
    direction = case_when(
      logFC >=  2 ~ "DSS_HFID",  # right side
      logFC <= -2 ~ "HFID",      # left side
      TRUE ~ "NS"
    )
  ) %>%
  mutate(
    Taxa = gsub("\\[|\\]", "", Taxa),
    Taxa = trimws(Taxa),
    Taxa = ifelse(is.na(Taxa) | Taxa == "" | Taxa == "NA", NA, Taxa)
  ) %>%
  filter(!is.na(Taxa))

# -----------------------------
# select representative taxa
# -----------------------------
dss_taxa <- res_kegg %>%
  filter(sig == "sig", direction == "DSS_HFID") %>%
  arrange(desc(logFC)) %>%
  distinct(Taxa, .keep_all = TRUE) %>%
  slice_head(n = 10) %>%
  pull(Taxa)

hfid_taxa <- res_kegg %>%
  filter(sig == "sig", direction == "HFID") %>%
  arrange(logFC) %>%
  distinct(Taxa, .keep_all = TRUE) %>%
  slice_head(n = 10) %>%
  pull(Taxa)

# -----------------------------
# assign plotting groups
# -----------------------------
res_kegg <- res_kegg %>%
  mutate(
    Taxa_group = case_when(
      sig == "sig" & Taxa %in% dss_taxa ~ paste0("DSS_", Taxa),
      sig == "sig" & Taxa %in% hfid_taxa ~ paste0("HFID_", Taxa),
      TRUE ~ "NS"
    )
  )

# -----------------------------
# color palettes
# -----------------------------
dss_colors <- setNames(
  c("#ffffb2","#fed976","#feb24c","#fd8d3c","#fc4e2a",
    "#e31a1c","#bd0026","#f768a1","#fbb4b9","#f768a1")[seq_along(dss_taxa)],
  paste0("DSS_", dss_taxa)
)

hfid_colors <- setNames(
  c("#c7e9c0","#a1d99b","#74c476","#41ab5d","#238b45",
    "#66c2a4","#41b6c4","#2b8cbe","#253494","#54278f")[seq_along(hfid_taxa)],
  paste0("HFID_", hfid_taxa)
)

palette_all <- c(dss_colors, hfid_colors, NS = "grey85")

# -----------------------------
# counts
# -----------------------------
dss_count <- sum(res_kegg$direction == "DSS_HFID" & res_kegg$sig == "sig", na.rm = TRUE)
hfid_count <- sum(res_kegg$direction == "HFID" & res_kegg$sig == "sig", na.rm = TRUE)

x_min <- min(res_kegg$logFC, na.rm = TRUE)
x_max <- max(res_kegg$logFC, na.rm = TRUE)
y_max <- max(res_kegg$negLog10P, na.rm = TRUE)

# -----------------------------
# volcano plot
# -----------------------------
# -----------------------------
# desired legend order
# -----------------------------
legend_order <- c(
  names(hfid_colors),   # blues/greens first
  "NS",
  names(dss_colors)     # reds/oranges second
)

# -----------------------------
# volcano plot
# -----------------------------
p_volcano <- ggplot(res_kegg, aes(x = logFC, y = negLog10P)) +
  
  # non-significant
  geom_point(
    data = subset(res_kegg, sig == "ns"),
    color = "grey85",
    size = 1.5
  ) +
  
  # significant
  geom_point(
    data = subset(res_kegg, sig == "sig"),
    aes(color = Taxa_group),
    size = 4,
    alpha = 0.85
  ) +
  
  scale_color_manual(
    values = palette_all,
    
    # FORCE LEGEND ORDER
    breaks = legend_order,
    
    name = "Taxa",
    
    labels = function(x) {
      sapply(x, function(y) {
        if (y == "NS") return("NS")
        
        taxa_name <- sub("^(DSS_|HFID_)", "", y)
        
        bquote(italic(.(taxa_name)))
      })
    },
    
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE
    )
  ) +
  
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  annotate(
    "text",
    x = x_max,
    y = y_max,
    label = paste0("DSS_HFID: ", dss_count),
    hjust = 1,
    vjust = 1,
    fontface = "bold",
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_min,
    y = y_max,
    label = paste0("HFID: ", hfid_count),
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 4
  ) +
  
  theme_classic(base_size = 13) +
  
  labs(
    title = "Volcano Plot: carb pathways DSS_HFID vs HFID (males)",
    x = "log2 Fold Change",
    y = "-log10 adjusted P-value"
  )

# -----------------------------
# save
# -----------------------------
ggsave(
  "volcano_all_carb_taxa_colored_final.png",
  p_volcano,
  width = 9,
  height = 6,
  dpi = 600
)

# -----------------------------
# EXPORT THE DATA
# -----------------------------
library(openxlsx)
library(dplyr)

# -----------------------------
# Build export components
# -----------------------------

# Full limma results
limma_all <- res %>%
  mutate(
    negLog10P = -log10(adj.P.Val)
  )

# Annotated full results
limma_annotated <- res %>%
  left_join(annot, by = "Protein")

# KEGG-filtered results (your main plot dataset)
limma_kegg <- res_kegg

# Summary stats
summary_stats <- limma_kegg %>%
  summarise(
    n_total = n(),
    n_sig = sum(sig == "sig"),
    n_up = sum(direction == "Up"),
    n_down = sum(direction == "Down"),
    n_ns = sum(direction == "NS")
  )

# Taxa-level summary
taxa_summary <- limma_kegg %>%
  group_by(Taxa, direction) %>%
  summarise(
    n = n(),
    mean_logFC = mean(logFC, na.rm = TRUE),
    min_padj = min(adj.P.Val, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(direction, desc(n))

# Top hits (for reporting)
top_hits <- limma_kegg %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50)

# -----------------------------
# Write Excel workbook
# -----------------------------
wb <- createWorkbook()

addWorksheet(wb, "limma_all")
writeData(wb, "limma_all", limma_all)

addWorksheet(wb, "limma_annotated")
writeData(wb, "limma_annotated", limma_annotated)

addWorksheet(wb, "limma_kegg_filtered")
writeData(wb, "limma_kegg_filtered", limma_kegg)

addWorksheet(wb, "summary_stats")
writeData(wb, "summary_stats", summary_stats)

addWorksheet(wb, "taxa_summary")
writeData(wb, "taxa_summary", taxa_summary)

addWorksheet(wb, "top_hits")
writeData(wb, "top_hits", top_hits)

saveWorkbook(
  wb,
  file = "limma_volcano_ALL_carb_dss_hfid_vs_hfid_males_STAT_OUTPUT.xlsx",
  overwrite = TRUE
)

# =========================================================
# STACKED BAR OF A SELECT PROTEIN BY PRODUCING TAXA
# =========================================================

############################################################################################################
#KEGG BASED PROTEOMICS ANALYSIS INTEGRATING PROTEINS, KEGG PATHWAYS, AND TAXA OF ORIGIN
############################################################################################################


#set working directory
setwd("/Users/kristinasorokolet/Library/CloudStorage/Box-Box/Wright Lab Operations/Student Folders/Kris/KS_2/Proteomics/GLOBAL/RANDOM_BINS/MALES_CTABS/MS9_1P_75_unclus/KEGG_new/")

metadata_file <- "metadata.csv"

# -------------------------------
# Read metadata
# -------------------------------
meta <- read_csv(metadata_file, show_col_types = FALSE)

meta <- meta %>%
  mutate(group = recode(group,
                        "group1" = "dss_hfid",
                        "group2" = "dss",
                        "group3" = "control",
                        "group4" = "hfid"))


# Read data
crosstab <- read.csv("imputed_w_taxa_kegg.csv", header = TRUE, stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(rstatix)
library(ggpubr)
library(openxlsx)

#DO FIRST TWO STEPS FOR ANY TYPE OF PLOT BELOW

# -----------------------------
# Convert log2 → linear
# -----------------------------
df_linear <- crosstab %>%
  mutate(across(starts_with("LI"), ~ 2^. - 1))

# -----------------------------
# Long format + convert to relative abundance
# -----------------------------
df_rel <- df_linear %>%
  pivot_longer(
    cols = starts_with("LI"),
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  group_by(SampleID) %>%
  mutate(RelAbundance = Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  separate_rows(KEGG_Pathway, sep = ",\\s*") %>%
  filter(!is.na(KEGG_Pathway) & KEGG_Pathway != "") %>%
  group_by(Protein, SampleID) %>%
  mutate(RelAbundance = RelAbundance / n()) %>%
  ungroup()


library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(openxlsx)
library(ggh4x)

# ---------------------------------------------------------
# FUNCTION
# ---------------------------------------------------------
plot_protein_taxa <- function(
    protein_name,
    df_rel,
    meta,
    sex_label = "males",
    output_prefix = NULL
) {
  
  if (is.null(output_prefix)) {
    output_prefix <- gsub("[^A-Za-z0-9]+", "_", protein_name)
  }
  
  # -----------------------------
  # Filter protein
  # -----------------------------
  df_filtered <- df_rel %>%
    filter(grepl(protein_name, Function, ignore.case = TRUE)) %>%
    left_join(meta, by = c("SampleID" = "sample_id"))
  
  # -----------------------------
  # Taxa-level data
  # -----------------------------
  df_group <- df_filtered %>%
    group_by(group, Taxa) %>%
    summarise(
      RelAbundance = sum(RelAbundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Taxa = ifelse(is.na(Taxa) | Taxa == "", "unknown", Taxa),
      Taxa_clean = str_replace_all(Taxa, "\\[|\\]", ""),
      Taxa_name = ifelse(
        Taxa_clean != "unknown",
        paste0("italic('", Taxa_clean, "')"),
        "'unknown'"
      )
    )
  
  # -----------------------------
  # Sample-level data
  # -----------------------------
  df_protein <- df_filtered %>%
    group_by(SampleID, group) %>%
    summarise(
      ProteinRelAbundance = sum(RelAbundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    ungroup()
  
  group_levels <- c("control", "hfid", "dss", "dss_hfid")
  
  df_group$group <- factor(df_group$group, levels = group_levels)
  df_protein$group <- factor(df_protein$group, levels = group_levels)
  
  # -----------------------------
  # Factor levels
  # -----------------------------
  group_levels <- c("control", "hfid", "dss", "dss_hfid")
  
  group_labels <- c("Control", "HFiD", "DSS", "DSS+HFiD")

  
  # -----------------------------
  # Pairwise stats
  # -----------------------------
  stats_posthoc_all <- df_protein %>%
    pairwise_wilcox_test(
      ProteinRelAbundance ~ group,
      p.adjust.method = "BH"
    )
  
  stats_posthoc_plot <- stats_posthoc_all %>%
    filter(
      paste(group1, group2, sep = "_vs_") %in% c(
        "control_vs_dss_hfid",
        "control_vs_dss",
        "hfid_vs_dss_hfid"
      )
    )
  
  # -----------------------------
  # Bracket positions
  # -----------------------------
  # -----------------------------
  # Bracket positions (FIXED)
  # -----------------------------
  y_max_all <- df_group %>%
    group_by(group) %>%
    summarise(total = sum(RelAbundance), .groups = "drop") %>%
    summarise(max_total = max(total)) %>%
    pull(max_total)
  
  step_increase <- 0.12 * y_max_all
  start_height  <- 1.15 * y_max_all
  
  stats_posthoc_plot <- stats_posthoc_plot %>%
    mutate(
      xmin = group1,
      xmax = group2,
      y.position = start_height + (row_number() * step_increase)
    )
  
  # -----------------------------
  # Colors
  # -----------------------------
  n_taxa <- length(unique(df_group$Taxa_name))
  
  palette <- colorRampPalette(c(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
    "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7"
  ))(n_taxa)
  
  names(palette) <- unique(df_group$Taxa_name)
  
  taxa_levels <- unique(df_group$Taxa_name)
  
  taxa_labels_wrapped <- setNames(
    stringr::str_wrap(taxa_levels, width = 25),
    taxa_levels
  )
  
  # -----------------------------
  # Plot
  # -----------------------------
  # -----------------------------
  # Plot
  # -----------------------------
  p <- ggplot(
    df_group,
    aes(x = group, y = RelAbundance, fill = Taxa_name)
  ) +
    geom_col(color = "black") +
    scale_x_discrete(
      labels = c(
        control = "Control",
        hfid = "HFiD",
        dss = "DSS",
        dss_hfid = "DSS+HFiD"
      )
    ) +
    scale_fill_manual(
      name = "Taxa",
      values = palette,
      labels = setNames(
        parse(text = taxa_levels),
        taxa_levels
      )
    ) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
    theme_classic(base_size = 12) +
    theme(
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.4, "cm")
    ) +
    stat_pvalue_manual(
      stats_posthoc_plot,
      inherit.aes = FALSE,
      mapping = aes(
        x = group1,
        xend = group2,
        y = y.position,
        label = p.adj.signif
      ),
      tip.length = 0.01
    ) +
    labs(
      title = paste0(protein_name, " (", sex_label, ")"),
      x = "Treatment Group",
      y = "Relative abundance (%)"
    ) + 
    # --> ADD THIS LINE <--
    ggh4x::force_panelsizes(rows = unit(4, "in"), cols = unit(2.25, "in")) 
  
  print(p)
  
  # -----------------------------
  # Save plot
  # -----------------------------
  ggsave(
    paste0(output_prefix, "_plot.png"),
    p,
    width = 7,  # You may need to increase the total width to ensure wide legends aren't clipped
    height = 5,
    dpi = 600
  )
  
  # -----------------------------
  # Summary tables
  # -----------------------------
  group_summary <- df_protein %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean = mean(ProteinRelAbundance, na.rm = TRUE),
      median = median(ProteinRelAbundance, na.rm = TRUE),
      sd = sd(ProteinRelAbundance, na.rm = TRUE),
      min = min(ProteinRelAbundance, na.rm = TRUE),
      max = max(ProteinRelAbundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  # -----------------------------
  # Export Excel
  # -----------------------------
  wb <- createWorkbook()
  
  addWorksheet(wb, "raw_stats_input")
  writeData(wb, "raw_stats_input", df_protein)
  
  addWorksheet(wb, "group_summary")
  writeData(wb, "group_summary", group_summary)
  
  addWorksheet(wb, "pairwise_all_stats")
  writeData(wb, "pairwise_all_stats", stats_posthoc_all)
  
  addWorksheet(wb, "pairwise_plot_stats")
  writeData(wb, "pairwise_plot_stats", stats_posthoc_plot)
  
  addWorksheet(wb, "taxa_plot_input")
  writeData(wb, "taxa_plot_input", df_group)
  
  saveWorkbook(
    wb,
    paste0(output_prefix, "_results.xlsx"),
    overwrite = TRUE
  )
  
  return(p)
}

#SPECIFY PROTEINS YOU WANT TO PLOT BY TAXA

proteins <- c(
  "UDP-glucose 4-epimerase GalE",
  "Septation",
  "FAD-dependent oxidoreductase",
  "LacI",
  "GDP-mannose 4,6-dehydratase",
  "L-rhamnose isomerase",
  "UDP-glucose 4-epimerase",
  "pyruvate:ferredoxin",
  "acetate kinase",
  "NADP-specific glutamate dehydrogenase",
  "UDP-glucose--hexose-1-phosphate uridylyltransferase"
)

#CALL THE FUNCTION FOR SPECIFIED PROTEINS
for (prot in proteins) {
  plot_protein_taxa(
    protein_name = prot,
    df_rel = df_rel,
    meta = meta
  )
}

############################################################################################
#STACKED KEGG PATHWAY CONTRIBUTION PLOTS (EITHER TAXA OR PROTEIN STACKED)
############################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(rstatix)
library(ggpubr)
library(openxlsx)

#this loop does both taxa and protein breakdown of a kegg pathway

kruskal_test <- rstatix::kruskal_test
pairwise_wilcox_test <- rstatix::pairwise_wilcox_test

# =========================================================
# LOOP FUNCTION:
# protein + taxa stacked bars for a KEGG pathway
# WITH FULL STAT EXPORTS
# =========================================================

plot_ko_pathway <- function(
    df_rel,
    meta,
    ko_number,
    pathway_name,
    excel_file = NULL,
    
    protein_width = 8,
    protein_height = 12,
    
    taxa_width = 8,
    taxa_height = 12,
    
    dpi = 600
) {
  
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(openxlsx)
  library(stringr)
  
  # -------------------------------------------------------
  # Clean taxa
  # -------------------------------------------------------
  df_rel <- df_rel %>%
    mutate(
      Taxa = ifelse(
        is.na(Taxa) | Taxa == "",
        "Unknown",
        Taxa
      )
    )
  
  # -------------------------------------------------------
  # Filter KO
  # -------------------------------------------------------
  df_kos <- df_rel %>%
    filter(grepl(ko_number, KEGG_Pathway, ignore.case = TRUE)) %>%
    left_join(meta, by = c("SampleID" = "sample_id"))
  
  group_levels <- c("control", "hfid", "dss", "dss_hfid")
  
  # =======================================================
  # PROTEIN CONTRIBUTION
  # =======================================================
  
  df_protein <- df_kos %>%
    group_by(SampleID, group, Function) %>%
    summarise(
      RelAbundance = sum(RelAbundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  # -----------------------------
  # pathway abundance for stats
  # -----------------------------
  df_pathway_protein <- df_protein %>%
    group_by(SampleID, group) %>%
    summarise(
      PathwayRelAbundance = sum(RelAbundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_pathway_protein$group <- factor(
    df_pathway_protein$group,
    levels = group_levels
  )
  
  # -----------------------------
  # Kruskal-Wallis
  # -----------------------------
  stats_kw_protein <- rstatix::kruskal_test(df_pathway_protein, PathwayRelAbundance ~ group)
  
  # -----------------------------
  # Pairwise Wilcoxon
  # -----------------------------
  stats_posthoc_protein <- df_pathway_protein %>%
    pairwise_wilcox_test(
      PathwayRelAbundance ~ group,
      p.adjust.method = "BH",
      detailed = TRUE
    )
  
  # -----------------------------
  # Summary stats
  # -----------------------------
  protein_group_summary <- df_pathway_protein %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean = mean(PathwayRelAbundance, na.rm = TRUE),
      median = median(PathwayRelAbundance, na.rm = TRUE),
      sd = sd(PathwayRelAbundance, na.rm = TRUE),
      sem = sd / sqrt(n),
      min = min(PathwayRelAbundance, na.rm = TRUE),
      max = max(PathwayRelAbundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  # -----------------------------
  # Plot data
  # -----------------------------
  df_protein_group <- df_protein %>%
    group_by(group, Function) %>%
    summarise(
      RelAbundance = sum(RelAbundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_protein_group$group <- factor(
    df_protein_group$group,
    levels = group_levels
  )
  
  df_protein_group <- df_protein_group %>%
    mutate(
      Function_clean = Function %>%
        gsub("^MAG:\\s*", "", .) %>%
        gsub("\\s*n=\\d+.*$", "", .) %>%
        trimws()
    ) %>%
    mutate(
      Function_clean = case_when(
        Function_clean %in% c(
          "ornithine carbamoyltransferase",
          "Ornithine carbamoyltransferase"
        ) ~ "ornithine carbamoyltransferase",
        Function_clean == "status=active" ~ "unannotated",
        TRUE ~ Function_clean
      )
    ) %>%
    group_by(group, Function_clean) %>%
    summarise(
      RelAbundance = sum(RelAbundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  # -----------------------------
  # Protein colors
  # -----------------------------
  unique_functions <- unique(df_protein_group$Function_clean)
  
  unique_functions <- c(
    setdiff(unique_functions, "unannotated"),
    "unannotated"
  )
  
  palette_protein <- colorRampPalette(c(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
    "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7"
  ))(length(unique_functions))
  
  names(palette_protein) <- unique_functions
  
  # -----------------------------
  # y positions
  # -----------------------------
  y_max_all <- df_protein_group %>%
    group_by(group) %>%
    summarise(
      max_bar = sum(RelAbundance),
      .groups = "drop"
    ) %>%
    pull(max_bar) %>%
    max()
  
  stats_posthoc_protein <- stats_posthoc_protein %>%
    arrange(p.adj) %>%
    mutate(
      y.position = y_max_all * (1.05 + 0.07 * row_number())
    )
  
  # -----------------------------
  # Protein plot
  # -----------------------------
  p_protein <- ggplot(
    df_protein_group,
    aes(x = group, y = RelAbundance, fill = Function_clean)
  ) +
    geom_col(color = "black") +
    scale_fill_manual(
      values = palette_protein,
      name = "Protein"
    ) +
    stat_pvalue_manual(
      stats_posthoc_protein,
      inherit.aes = FALSE,
      mapping = aes(
        x = group1,
        xend = group2,
        y = y.position,
        label = p.adj.signif
      ),
      tip.length = 0.01
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = paste(
        "Proteins contributing to",
        pathway_name,
        "(",
        ko_number,
        ")"
      ),
      x = "Treatment Group",
      y = "Relative abundance %"
    )
  
  ggsave(
    paste0(ko_number, "_protein_plot.png"),
    plot = p_protein,
    width = protein_width,
    height = protein_height,
    dpi = dpi
  )
  
  # =======================================================
  # TAXA CONTRIBUTION
  # =======================================================
  
  df_taxa <- df_kos %>%
    group_by(SampleID, group, Taxa) %>%
    summarise(
      RelAbundance = sum(RelAbundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  # -----------------------------
  # pathway abundance for stats
  # -----------------------------
  df_pathway_taxa <- df_taxa %>%
    group_by(SampleID, group) %>%
    summarise(
      PathwayRelAbundance = sum(RelAbundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_pathway_taxa$group <- factor(
    df_pathway_taxa$group,
    levels = group_levels
  )
  
  # -----------------------------
  # Kruskal-Wallis
  # -----------------------------
  stats_kw_taxa <- rstatix::kruskal_test(df_pathway_taxa, PathwayRelAbundance ~ group)
  
  # -----------------------------
  # Pairwise Wilcoxon
  # -----------------------------
  stats_posthoc_taxa <- df_pathway_taxa %>%
    pairwise_wilcox_test(
      PathwayRelAbundance ~ group,
      p.adjust.method = "BH",
      detailed = TRUE
    )
  
  # -----------------------------
  # Summary stats
  # -----------------------------
  taxa_group_summary <- df_pathway_taxa %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean = mean(PathwayRelAbundance, na.rm = TRUE),
      median = median(PathwayRelAbundance, na.rm = TRUE),
      sd = sd(PathwayRelAbundance, na.rm = TRUE),
      sem = sd / sqrt(n),
      min = min(PathwayRelAbundance, na.rm = TRUE),
      max = max(PathwayRelAbundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  # -----------------------------
  # Taxa plot data
  # -----------------------------
  df_taxa_group <- df_taxa %>%
    group_by(group, Taxa) %>%
    summarise(
      RelAbundance = sum(RelAbundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_taxa_group$group <- factor(
    df_taxa_group$group,
    levels = group_levels
  )
  
  # -----------------------------
  # Taxa colors
  # -----------------------------
  unique_taxa <- unique(df_taxa_group$Taxa)
  
  palette_taxa <- colorRampPalette(c(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
    "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7"
  ))(length(unique_taxa))
  
  names(palette_taxa) <- unique_taxa
  
  # -----------------------------
  # y positions
  # -----------------------------
  y_max_all_taxa <- df_taxa_group %>%
    group_by(group) %>%
    summarise(
      max_bar = sum(RelAbundance),
      .groups = "drop"
    ) %>%
    pull(max_bar) %>%
    max()
  
  stats_posthoc_taxa <- stats_posthoc_taxa %>%
    arrange(p.adj) %>%
    mutate(
      y.position = y_max_all_taxa * (1.05 + 0.07 * row_number())
    )
  
  # -----------------------------
  # Taxa plot
  # -----------------------------
  p_taxa <- ggplot(
    df_taxa_group,
    aes(x = group, y = RelAbundance, fill = Taxa)
  ) +
    geom_col(color = "black") +
    scale_fill_manual(
      values = palette_taxa,
      name = "Taxa"
    ) +
    stat_pvalue_manual(
      stats_posthoc_taxa,
      inherit.aes = FALSE,
      mapping = aes(
        x = group1,
        xend = group2,
        y = y.position,
        label = p.adj.signif
      ),
      tip.length = 0.01
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = paste(
        "Taxa contributing to",
        pathway_name,
        "(",
        ko_number,
        ")"
      ),
      x = "Treatment Group",
      y = "Relative abundance %"
    )
  
  ggsave(
    paste0(ko_number, "_taxa_plot.png"),
    plot = p_taxa,
    width = taxa_width,
    height = taxa_height,
    dpi = dpi
  )
  
  # =======================================================
  # EXPORT EXCEL
  # =======================================================
  
  if (is.null(excel_file)) {
    
    safe_ko <- gsub("[^A-Za-z0-9_]", "_", ko_number)
    
    excel_file <- paste0(
      safe_ko,
      "_stats.xlsx"
    )
  }
  
  wb <- createWorkbook()
  
  # -----------------------------
  # RAW VALUES USED FOR STATS
  # -----------------------------
  addWorksheet(wb, "protein_raw_values")
  writeData(wb, "protein_raw_values", df_pathway_protein)
  
  addWorksheet(wb, "taxa_raw_values")
  writeData(wb, "taxa_raw_values", df_pathway_taxa)
  
  # -----------------------------
  # GROUP SUMMARIES
  # -----------------------------
  addWorksheet(wb, "protein_group_summary")
  writeData(wb, "protein_group_summary", protein_group_summary)
  
  addWorksheet(wb, "taxa_group_summary")
  writeData(wb, "taxa_group_summary", taxa_group_summary)
  
  # -----------------------------
  # KRUSKAL WALLIS
  # -----------------------------
  addWorksheet(wb, "protein_kruskal")
  writeData(wb, "protein_kruskal", stats_kw_protein)
  
  addWorksheet(wb, "taxa_kruskal")
  writeData(wb, "taxa_kruskal", stats_kw_taxa)
  
  # -----------------------------
  # PAIRWISE TESTS
  # -----------------------------
  addWorksheet(wb, "protein_pairwise")
  writeData(wb, "protein_pairwise", stats_posthoc_protein)
  
  addWorksheet(wb, "taxa_pairwise")
  writeData(wb, "taxa_pairwise", stats_posthoc_taxa)
  
  # -----------------------------
  # PLOT INPUTS
  # -----------------------------
  addWorksheet(wb, "protein_plot_input")
  writeData(wb, "protein_plot_input", df_protein_group)
  
  addWorksheet(wb, "taxa_plot_input")
  writeData(wb, "taxa_plot_input", df_taxa_group)
  
  saveWorkbook(
    wb,
    file = excel_file,
    overwrite = TRUE
  )
  
  message("Excel file saved: ", normalizePath(excel_file))
  
  return(list(
    protein_plot = p_protein,
    taxa_plot = p_taxa,
    protein_stats = stats_posthoc_protein,
    taxa_stats = stats_posthoc_taxa
  ))
}


# Usage:
result <- plot_ko_pathway(
  df_rel = df_rel,
  meta = meta,
  ko_number = "ko00051",
  pathway_name = "Fructose metabolism - males",
  
  protein_width = 9,
  protein_height = 5,
  
  taxa_width = 7,
  taxa_height = 5,
  
  dpi = 600
)

# Usage:
result <- plot_ko_pathway(
  df_rel = df_rel,
  meta = meta,
  ko_number = "ko00052",
  pathway_name = "Galactose metabolism - males",
  
  protein_width = 9,
  protein_height = 5,
  
  taxa_width = 7,
  taxa_height = 5,
  
  dpi = 600
)


#RUN THIS TO ONLY DISPLAY TOP 10 TAXA/PROTEINS IN A KEGG PATHWAY AND THEN SHADE EVERYTHING ELSE GREY AS "OTHER"

## TOP 10 COLORED CONTRIBUTIONS 
library(dplyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(openxlsx)
library(ggtext)
library(stringr)
library(ggh4x)

# -----------------------------
# Helper
# -----------------------------
get_top_features <- function(df, feature_col, value_col, top_n = 10) {
  df %>%
    group_by(.data[[feature_col]]) %>%
    summarise(total = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total)) %>%
    slice_head(n = top_n) %>%
    pull(.data[[feature_col]])
}

# -----------------------------
# MAIN FUNCTION
# -----------------------------
plot_ko_pathway <- function(df_rel, meta, ko_number, pathway_name,
                            sex = "males",
                            top_n = 10,
                            excel_file = NULL) {
  
  df_rel <- df_rel %>%
    mutate(Taxa = ifelse(is.na(Taxa) | Taxa == "", "Unknown", Taxa))
  
  df_kos <- df_rel %>%
    filter(grepl(ko_number, KEGG_Pathway, ignore.case = TRUE)) %>%
    left_join(meta, by = c("SampleID" = "sample_id"))
  
  # =========================================================
  # PROTEINS
  # =========================================================
  
  df_protein <- df_kos %>%
    group_by(SampleID, group, Function) %>%
    summarise(RelAbundance = sum(RelAbundance, na.rm = TRUE), .groups="drop")
  
  df_protein_group <- df_protein %>%
    mutate(Function_clean = Function %>%
             gsub("^MAG:\\s*", "", .) %>%
             gsub("\\s*n=\\d+.*$", "", .) %>%
             gsub("\\s*RepID=.*$", "", .) %>%  # Cleaned RepID tag
             trimws()) %>%
    mutate(Function_clean = case_when(
      grepl("status=", Function_clean, ignore.case = TRUE) ~ "unannotated",
      TRUE ~ Function_clean
    )) %>%
    group_by(group, Function_clean) %>%
    summarise(RelAbundance = sum(RelAbundance, na.rm = TRUE), .groups="drop") %>%
    mutate(group = factor(group, levels = c("control","hfid","dss","dss_hfid")))
  
  top_proteins <- get_top_features(df_protein_group, "Function_clean", "RelAbundance", top_n)
  protein_levels <- c("Other", top_proteins)
  
  # Wrap long protein names into multiple lines (e.g., max 30 characters wide)
  protein_labels_wrapped <- setNames(
    stringr::str_wrap(protein_levels, width = 30),
    protein_levels
  )
  
  df_protein_group <- df_protein_group %>%
    mutate(
      Feature_group = ifelse(Function_clean %in% top_proteins,
                             Function_clean,
                             "Other")
    ) %>%
    group_by(group, Feature_group) %>%
    summarise(RelAbundance = sum(RelAbundance, na.rm = TRUE), .groups = "drop")
  
  df_protein_group$Feature_group <- factor(df_protein_group$Feature_group,
                                           levels = protein_levels)
  
  palette_protein <- c(
    Other = "grey80",
    setNames(
      colorRampPalette(c("#4E79A7","#F28E2B","#E15759","#76B7B2",
                         "#59A14F","#EDC948","#B07AA1","#FF9DA7"))(
                           length(top_proteins)
                         ),
      top_proteins
    )
  )
  
  # =========================================================
  # TAXA
  # =========================================================
  
  df_taxa <- df_kos %>%
    group_by(SampleID, group, Taxa) %>%
    summarise(RelAbundance = sum(RelAbundance, na.rm = TRUE), .groups="drop")
  
  df_taxa_group <- df_taxa %>%
    mutate(Taxa = gsub("\\[|\\]", "", Taxa)) %>%
    group_by(group, Taxa) %>%
    summarise(RelAbundance = sum(RelAbundance, na.rm = TRUE), .groups="drop") %>%
    mutate(group = factor(group, levels = c("control","hfid","dss","dss_hfid")))
  
  top_taxa <- get_top_features(df_taxa_group, "Taxa", "RelAbundance", top_n)
  taxa_levels <- c("Other", top_taxa)
  
  # ✔ FIX: correct collapse AFTER labeling
  df_taxa_group <- df_taxa_group %>%
    mutate(
      Feature_group = ifelse(Taxa %in% top_taxa,
                             Taxa,
                             "Other")
    ) %>%
    group_by(group, Feature_group) %>%
    summarise(RelAbundance = sum(RelAbundance, na.rm = TRUE), .groups = "drop")
  
  df_taxa_group$Feature_group <- factor(df_taxa_group$Feature_group,
                                        levels = taxa_levels)
  
  taxa_labels <- setNames(
    c("Other", paste0("*", top_taxa, "*")),
    c("Other", top_taxa)
  )
  
  palette_taxa <- c(
    Other = "grey80",
    setNames(
      colorRampPalette(c("#4E79A7","#F28E2B","#E15759","#76B7B2",
                         "#59A14F","#EDC948","#B07AA1","#FF9DA7"))(
                           length(top_taxa)
                         ),
      top_taxa
    )
  )
  
  # =========================================================
  # STATS
  # =========================================================
  
  df_p <- df_protein %>%
    group_by(SampleID, group) %>%
    summarise(value = sum(RelAbundance), .groups="drop")
  
  kw_protein <- df_p %>% kruskal_test(value ~ group)
  pw_protein <- df_p %>% pairwise_wilcox_test(value ~ group, p.adjust.method="BH")
  
  df_t <- df_taxa %>%
    group_by(SampleID, group) %>%
    summarise(value = sum(RelAbundance), .groups="drop")
  
  kw_taxa <- df_t %>% kruskal_test(value ~ group)
  pw_taxa <- df_t %>% pairwise_wilcox_test(value ~ group, p.adjust.method="BH")
  
  # y positions
  y_p <- max(tapply(df_protein_group$RelAbundance,
                    df_protein_group$group, sum))
  
  pw_protein <- pw_protein %>%
    arrange(p.adj) %>%
    mutate(y.position = y_p * (1.05 + 0.07 * row_number()))
  
  y_t <- max(tapply(df_taxa_group$RelAbundance,
                    df_taxa_group$group, sum))
  
  pw_taxa <- pw_taxa %>%
    arrange(p.adj) %>%
    mutate(y.position = y_t * (1.05 + 0.07 * row_number()))
  
  # =========================================================
  # PLOTS
  # =========================================================
  
  # -----------------
  # Protein Plot
  # -----------------
  p_protein <- ggplot(df_protein_group, 
                      aes(x = group, y = RelAbundance, fill = Feature_group)) +
    geom_col(color="black", position=position_stack(reverse=TRUE)) +
    scale_x_discrete(
      labels = c(
        "control" = "Control",
        "hfid" = "HFiD",
        "dss" = "DSS",
        "dss_hfid" = "DSS+HFiD"
      )
    ) +
    scale_fill_manual(values = palette_protein, 
                      breaks = rev(protein_levels), 
                      labels = protein_labels_wrapped, # Two-line wrap
                      name = "Protein") +
    stat_pvalue_manual(
      pw_protein,
      label="p.adj.signif",
      xmin="group1",
      xmax="group2",
      y.position="y.position",
      tip.length=0.01
    ) +
    labs(
      title = paste0("Proteins contributing to ", pathway_name, 
                     " (", ko_number, ") in ", sex),
      x = "Treatment Group"
    ) +
    theme_classic() +
    # Enforce strict pixel size for axes panel (Width x Height)
    ggh4x::force_panelsizes(rows = unit(5, "in"), cols = unit(3, "in"))
  
  # -----------------
  # Taxa Plot
  # -----------------
  p_taxa <- ggplot(df_taxa_group, 
                   aes(x = group, y = RelAbundance, fill = Feature_group)) +
    geom_col(color="black", position=position_stack(reverse=TRUE)) +
    
    # --> ADD THIS BLOCK <--
    scale_x_discrete(
      labels = c(
        "control" = "Control",
        "hfid" = "HFiD",
        "dss" = "DSS",
        "dss_hfid" = "DSS+HFiD"
      )
    ) +
    
    scale_fill_manual(values = palette_taxa, 
                      breaks = rev(taxa_levels), 
                      labels = taxa_labels, 
                      name = "Taxa") +
    theme_classic() +
    theme(legend.text = element_markdown()) +
    stat_pvalue_manual(
      pw_taxa,
      label="p.adj.signif",
      xmin="group1",
      xmax="group2",
      y.position="y.position",
      tip.length=0.01
    ) +
    labs(
      title = paste0("Taxa contributing to ", pathway_name, 
                     " (", ko_number, ") in ", sex),
      x = "Treatment Group" # Optional: Cleans up the x-axis title
    )
  
  
  # =========================================================
  # EXCEL EXPORT
  # =========================================================
  if (!is.null(excel_file)) {
    wb <- createWorkbook()
    addWorksheet(wb, "KW_Protein"); writeData(wb,1,kw_protein)
    addWorksheet(wb, "Wilcox_Protein"); writeData(wb,2,pw_protein)
    addWorksheet(wb, "KW_Taxa"); writeData(wb,3,kw_taxa)
    addWorksheet(wb, "Wilcox_Taxa"); writeData(wb,4,pw_taxa)
    saveWorkbook(wb, excel_file, overwrite=TRUE)
  }
  
  # SAVE (Slightly increased total ggsave width to give wrapped text extra breathing room)
  ggsave(paste0(ko_number,"_protein.png"), p_protein, width=8, height=6, dpi=600)
  ggsave(paste0(ko_number,"_taxa.png"), p_taxa, width=8, height=8, dpi=600)
  
  return(list(
    protein_plot = p_protein,
    taxa_plot = p_taxa,
    kw_protein = kw_protein,
    kw_taxa = kw_taxa,
    pw_protein = pw_protein,
    pw_taxa = pw_taxa
  ))
}

result <- plot_ko_pathway(
  df_rel = df_rel,
  meta = meta,
  ko_number = "ko00051",
  pathway_name = "Fructose metabolism",
  top_n = 10,
  excel_file = "ko00051_stats_2.xlsx"
)

result <- plot_ko_pathway(
  df_rel = df_rel,
  meta = meta,
  ko_number = "ko00052",
  pathway_name = "Galactose metabolism",
  top_n = 10,
  excel_file = "ko00052_stats_2.xlsx"
)
