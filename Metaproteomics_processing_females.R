library(dplyr)
library(purrr)  
library(readr)
library(tidyr)
library(phyloseq)
library(readxl)
library(ggplot2)
library(writexl)

########################################################################################################################
#OUTLIER TEST on raw data
########################################################################################################################
#set working directory
setwd("/Users/kristinasorokolet/Library/CloudStorage/Box-Box/Wright Lab Operations/Student Folders/Kris/KS_2/Proteomics/GLOBAL/RANDOM_BINS/FEMALES_CTABS/MS9_1p_75_unclus")

df <- read_tsv("F_w_host_ms9_1p.tsv")
meta <- read_tsv("metadata.tsv")  # must contain columns: SampleID, Group

write_xlsx(df, "F_w_host_ms9_1p.xlsx")

meta<- meta %>%
  mutate(Group = recode(group,
                        "group1" = "DSS+HFiD",
                        "group2" = "DSS",
                        "group3" = "Control",
                        "group4" = "HFiD"))

# Sum intensities per column (sample)
intensity_cols <- grep("^LI", colnames(df), value = TRUE)
total_intensity <- colSums(df[, intensity_cols], na.rm = TRUE)
print(total_intensity)

#merge with metadata
total_df <- data.frame(
  Sample = names(total_intensity),
  TotalIntensity = total_intensity
) %>%
  left_join(meta, by = c("Sample" = "sample_id"))                    

#detect outliers per group (anything more than 3SD away)
total_df <- total_df %>%
  group_by(Group) %>%
  mutate(
    mean_group = mean(TotalIntensity),
    sd_group = sd(TotalIntensity),
    Outlier = ifelse(TotalIntensity < mean_group - 2*sd_group |
                       TotalIntensity > mean_group + 2*sd_group, "Outlier", "Normal")
  )

#plot outliers
ggplot(total_df, aes(x = Sample, y = TotalIntensity, fill = Outlier)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Group, scales = "free_x") +  # separates plots per group
  scale_fill_manual(values = c("Normal" = "blue", "Outlier" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Total Intensity per Sample by Group (pre-filtering)", y = "Total Intensity")

########################################################################################################################
###filtering all crosstabs in the directory with 50% and 75% group filters
########################################################################################################################
library(tidyverse)
library(readxl)
library(writexl)

meta <- read_tsv("metadata.tsv") %>%
  mutate(Group = recode(group,
                        "group1" = "DSS+HFiD",
                        "group2" = "DSS",
                        "group3" = "Control",
                        "group4" = "HFiD"))

# -----------------------------
# FUNCTION
# -----------------------------
row_passes_threshold <- function(row, meta, threshold = 0.5) {
  groups <- unique(meta$Group)
  
  for (g in groups) {
    samples <- meta %>% filter(Group == g) %>% pull(sample_id)
    samples <- intersect(samples, names(row))
    
    if (length(samples) == 0) next
    
    values <- as.numeric(row[samples])
    if (all(is.na(values))) next
    
    frac_present <- sum(values > 0, na.rm = TRUE) / length(samples)
    
    if (!is.na(frac_present) && frac_present >= threshold) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# ------------------------------------------
# LIST CROSS-TAB FILES
# ------------------------------------------
files <- list.files(pattern = "\\.tsv$")

# ------------------------------------------
# LOOP THROUGH FILES + THRESHOLDS
# ------------------------------------------
for (file in files) {
  message("Processing: ", file)
  
  # 1. READ THE TSV FILE
  df <- read_tsv(file)
  
  # 2. SAVE ORIGINAL FILE AS EXCEL  ← **THIS IS THE NEW LINE**
  write_xlsx(df, paste0(tools::file_path_sans_ext(file), ".xlsx"))
  
  # 3. IDENTIFY INTENSITY COLUMNS
  intensity_cols <- grep("^LI", colnames(df), value = TRUE)
  df_intensity <- df[, intensity_cols]
  
  # 4. LOOP THROUGH THRESHOLDS (50% and 75%)
  for (thr in c(0.50, 0.75)) {
    thr_label <- ifelse(thr == 0.5, "50", "75")
    message("  Threshold: ", thr_label)
    
    # Apply filtering
    keep <- apply(df_intensity, 1, row_passes_threshold,
                  meta = meta, threshold = thr)
    
    df_filtered <- df[keep, ]
    
    # Generate output names
    out_tsv  <- paste0(tools::file_path_sans_ext(file),
                       "_r_filtered", thr_label, ".tsv")
    out_xlsx <- paste0(tools::file_path_sans_ext(file),
                       "_r_filtered", thr_label, ".xlsx")
    
    # Save filtered files
    write_tsv(df_filtered, out_tsv)
    write_xlsx(df_filtered, out_xlsx)
    
    message("    Saved: ", out_tsv)
    message("    Saved: ", out_xlsx)
  }
}

message("DONE.")



##########################################################################################################################################################################################################
#start processing filtered crosstabs 
##########################################################################################################################################################################################################

library(dplyr)
library(purrr)  # make sure this is loaded
library(readr)
library(tidyr)
library(phyloseq)
library(readxl)
library(ggplot2)
library(openxlsx)
library(tidyverse)

#set working directory
setwd("/Users/kristinasorokolet/Library/CloudStorage/Box-Box/Wright Lab Operations/Student Folders/Kris/KS_2/Proteomics/GLOBAL/RANDOM_BINS/FEMALES_CTABS/MS9_1p_75_unclus/")

metadata_file <- "metadata.tsv"

meta <- read_tsv(metadata_file, show_col_types = FALSE)

meta <- meta %>%
  mutate(group = recode(group,
                        "group1" = "dss_hfid",
                        "group2" = "dss",
                        "group3" = "control",
                        "group4" = "hfid"))

#####################################################################################################
# Making taxa a separate column (move from Function in brackets to a separate column)
#####################################################################################################
library(readr)
library(dplyr)
library(stringr)
library(openxlsx)

ctab <- read_tsv("KS2_GBL_F_binned_ms9_1p_r_filtered75.tsv", show_col_types = FALSE)

# define known cofactors to exclude
cofactors <- c("ATP", "ADP", "GTP", "GDP", "NAD", "NADP", "FAD", "FMN")

ctab <- ctab %>%
  mutate(
    
    # -------------------------------------------------
    # 1. Extract Tax= format (highest confidence)
    # -------------------------------------------------
    Taxa_tax = str_extract(Function, "(?<=Tax=).*?(?=\\s+TaxID=)"),
    
    # -------------------------------------------------
    # 2. Extract ALL bracket contents
    # -------------------------------------------------
    bracket_all = str_extract_all(Function, "\\[[^\\]]+\\]"),
    
    # flatten list to string per row
    Taxa_bracket = sapply(bracket_all, function(x) {
      
      if (length(x) == 0) return(NA_character_)
      
      x <- str_remove_all(x, "\\[|\\]")
      
      # remove cofactors
      x <- x[!x %in% cofactors]
      
      # keep only plausible taxa-like entries (heuristic: contains space or lowercase)
      x <- x[str_detect(x, "[a-z]")]
      
      if (length(x) == 0) return(NA_character_)
      
      # assume remaining is taxonomy
      x[1]
    }),
    
    # -------------------------------------------------
    # 3. Final Taxa (Tax= overrides bracket)
    # -------------------------------------------------
    Taxa = coalesce(Taxa_tax, Taxa_bracket),
    
    # -------------------------------------------------
    # 4. Clean Function column
    # -------------------------------------------------
    Function = str_remove_all(Function, "\\bn=\\d+\\b"),
    Function = str_remove(Function, "\\s+Tax=.*?TaxID=\\d+\\s*"),
    Function = str_remove_all(Function, "\\[[^\\]]+\\]"),
    Function = str_remove(Function, "\\s+RepID=.*$"),
    Function = str_squish(Function)
  ) %>%
  select(-Taxa_tax, -Taxa_bracket, -bracket_all)

head(ctab)


write_csv(ctab, "ctab_w_taxa.csv")
write.xlsx(ctab, "ctab_w_taxa.xlsx")

####################################################################################
#VISUALIZE RAW DATA
####################################################################################

library(tidyverse)
library(tibble)

ctab_file <- "ctab_w_taxa.csv"

df_full <- read_csv(ctab_file, show_col_types = FALSE)

# -------------------------------------------------------------------
# Define annotation columns
# -------------------------------------------------------------------
annotation_cols <- c(
  "Cluster", "Gene", "Function", "Peptide.Number",
  "Unique Peptides.", "Shared Peptides.", "Peptides",
  "Flanked Peptides", "gene", "function", "Taxa", "KEGG_Pathway", "Peptide Number", 
  "Unique Peptide(s)", "Shared Peptide(s)"
)

# -------------------------------------------------------------------
# Split annotation + abundance
# -------------------------------------------------------------------

annotation <- df_full %>%
  select(any_of(annotation_cols)) %>%
  mutate(Protein = df_full[[1]])

df <- df_full %>%
  column_to_rownames(var = names(df_full)[1]) %>%
  select(-any_of(annotation_cols))

# -------------------------------
# Check matching
# -------------------------------
stopifnot(all(meta$sample_id %in% colnames(df)))
stopifnot(all(colnames(df) %in% meta$sample_id))

# -------------------------------
# Groups and colors
# -------------------------------
groups <- c("control", "hfid", "dss", "dss_hfid")
group_colors <- c("control" = "#E69F00",
                  "hfid" = "#56B4E9",
                  "dss" = "#009E73",
                  "dss_hfid" = "#F0E442")
#----------------
#plotting raw data for visual representation of distribution
#----------------
library(tidyverse)

df_long_raw <- df %>%
  as.data.frame() %>%
  rownames_to_column("protein") %>%
  pivot_longer(
    cols = -protein,
    names_to = "sample_id",
    values_to = "intensity"
  ) %>%
  left_join(meta, by = "sample_id") %>%
  filter(!is.na(intensity))

pdf("raw_distribution_plot.pdf", width = 10, height = 8)

ggplot(df_long_raw, aes(x = intensity, group = sample_id, color = group)) +
  geom_density(alpha = 0.4) +
  scale_x_continuous(trans = "log10") +
  theme_classic() +
  labs(
    title = "Raw DDA Intensities by Sample",
    x = "Raw intensity (log10 scale)",
    y = "Density"
  )

dev.off()


#----------------------------------------------
#fraction of missing values per sample plot
#----------------------------------------------

meta <- meta %>%
  mutate(group = tolower(group))

# Make sure group is a factor matching your color vector
group_colors <- c(
  "control" = "#E69F00",
  "hfid" = "#56B4E9",
  "dss" = "#009E73",
  "dss_hfid" = "#F0E442"
)

# Arrange samples by group and convert sample_id to factor for ordering
missing_by_sample <- df %>%
  as.data.frame() %>%
  summarise(across(everything(), ~ mean(. == 0 | is.na(.)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample_id",
    values_to = "frac_missing"
  ) %>%
  left_join(meta, by = "sample_id")

missing_by_sample <- missing_by_sample %>%
  arrange(group, sample_id) %>%
  mutate(
    group = factor(group, levels = names(group_colors)),
    sample_id = factor(sample_id, levels = unique(sample_id))
  )

# Compute overall average missingness
overall_avg <- mean(missing_by_sample$frac_missing)


pdf("fraction_missing_values.pdf", width = 10, height = 6)

ggplot(missing_by_sample, aes(x = sample_id, y = frac_missing, fill = group)) +
  geom_col() +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  geom_hline(yintercept = overall_avg, linetype = "dashed") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "top"
  ) +
  labs(
    title = "Fraction of Missing Values per Sample",
    y = "Fraction missing",
    x = "Sample"
  )


dev.off()
#-----------------------------------
#MNAR vs. MAR visual representation plot
#-----------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

# -------------------------------
# 1. Create long-format data
# -------------------------------
df_long <- df %>%
  as.data.frame() %>%
  rownames_to_column("protein") %>%
  pivot_longer(
    cols = -protein,
    names_to = "sample_id",
    values_to = "intensity"
  ) %>%
  left_join(meta %>% select(sample_id, group), by = "sample_id")

# -------------------------------
# 2. Compute per-protein per-group stats
# -------------------------------
group_detection <- df_long %>%
  group_by(protein, group) %>%
  summarize(
    detection_rate = mean(intensity > 0, na.rm = TRUE),
    med_intensity = median(intensity[intensity > 0], na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------
# 3. Compute per-group median intensity (top/bottom 50%)
# -------------------------------
group_intensity_median <- group_detection %>%
  group_by(group) %>%
  summarize(
    med_intensity_threshold = median(med_intensity, na.rm = TRUE),
    .groups = "drop"
  )

# Ensure Group types match for join
group_detection$Group <- as.character(group_detection$group)
group_intensity_median$Group <- as.character(group_intensity_median$group)

# -------------------------------
# 4. Join thresholds back
# -------------------------------
group_detection <- group_detection %>%
  left_join(group_intensity_median, by = "group")

# -------------------------------
# 5. Assign missing type (MNAR / MAR / Observed)
# -------------------------------
group_detection <- group_detection %>%
  mutate(
    missing_type = case_when(
      detection_rate >= 0.999 ~ "Observed",  # fully detected in group (floating point tolerant)
      detection_rate >= 0.75 ~ "MAR",        # ≥75% present
      detection_rate >= 0.5 & detection_rate < 0.75 & med_intensity >= med_intensity_threshold ~ "MAR",
      detection_rate >= 0.5 & detection_rate < 0.75 & med_intensity < med_intensity_threshold ~ "MNAR",
      detection_rate < 0.5 ~ "MNAR"
    )
  )

# -------------------------------
# 6. Plot
# -------------------------------
pdf("Protein_Detection_MNAR_MAR_Observed.pdf", width = 10, height = 8)

group_colors <- c(
  "MNAR" = "#E69F00",
  "MAR" = "#56B4E9",
  "Observed" = "grey50"
)

ggplot(group_detection, aes(x = med_intensity, y = detection_rate, color = missing_type)) +
  geom_point(alpha = 0.6) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "black") +
  scale_x_continuous(trans = "log10") +
  scale_color_manual(values = group_colors) +
  facet_wrap(~ group) +
  theme_classic() +
  labs(
    title = "Protein Detection vs Abundance by Group (MNAR/MAR/Observed)",
    x = "Median protein intensity",
    y = "Detection rate",
    color = "Missingness type"
  )

dev.off()  # close the PDF device

#####################################################################################################
#LOG TRANSFORM PROTEIN INTENSITIES 
#####################################################################################################
# -------------------------
# Log2 transform (DDA-safe)
# -------------------------

df_log2 <- df
df_log2[df_log2 == 0] <- NA
df_log2 <- log2(df_log2)

cat("NAs after log2 conversion:", sum(is.na(df_log2)), "\n")
cat("Percentage missing:",
    round(100 * sum(is.na(df_log2)) / length(df_log2), 1), "%\n")

# -------------------------
# Get annotation columns
# -------------------------

meta_cols <- c("Cluster", "Gene", "Function", "Peptide.Number", 
               "Unique.Peptide.s.", "Shared.Peptide.s.", "Peptides", 
               "Flanked.Peptides", "gene", "function", "Taxa", "KEGG_Pathway")

annotation <- read_csv(ctab_file, show_col_types = FALSE) %>%
  select(1, any_of(meta_cols))

colnames(annotation)[1] <- "Protein"

# -------------------------
# Join annotation to log2 data
# -------------------------

df_full <- df_log2 %>%
  as.data.frame() %>%
  rownames_to_column("Protein") %>%
  left_join(annotation, by = "Protein")

# -------------------------
# Export
# -------------------------

write.xlsx(df_full, "log2_data.xlsx", rowNames = FALSE)

################################################################################################
#NORMALIZE BASED ON MEAN
################################################################################################

# Convert to numeric matrix while keeping rownames
df_log2_mat <- as.matrix(df_log2)
mode(df_log2_mat) <- "numeric"  

# Output matrix
n_df_final <- df_log2_mat

# Make sure sample IDs in metadata match column names
meta <- meta   # rename if your metadata is called something else
all(colnames(df_log2_mat) %in% meta$sample_id)  # should be TRUE

# Create a group vector aligned to columns
groups <- meta$group[match(colnames(df_log2_mat), meta$sample_id)]
groups <- factor(groups)  # convert to factor
names(groups) <- colnames(df_log2_mat)  # name vector by sample IDs

# Loop over each group
n_df_final <- df_log2_mat

for (g in levels(groups)) {
  samp_idx <- which(groups == g)
  
  if (length(samp_idx) == 0) next  # skip empty groups
  
  group_data <- df_log2_mat[, samp_idx, drop = FALSE]
  
  sample_medians <- apply(group_data, 2, median, na.rm = TRUE)
  group_ref_median <- median(sample_medians, na.rm = TRUE)
  
  n_df_final[, samp_idx] <- sweep(group_data, 2, sample_medians, FUN = "-") + group_ref_median
  
  # Check table
  df_check <- data.frame(
    sample = colnames(group_data),
    before = round(sample_medians,3),
    after  = round(apply(n_df_final[, samp_idx], 2, median, na.rm = TRUE),3)
  )
  
  cat("\nGroup:", g, "\n")
  print(df_check)
}

print(paste("Group:", g, "Group median:", group_ref_median))

n_df_final_df <- as.data.frame(n_df_final)

n_df_final_id_full <- n_df_final_df %>%
  rownames_to_column("Protein") %>%
  left_join(annotation, by = "Protein")

write.xlsx(n_df_final_id_full, "normalized_log2_data.xlsx", rowNames = FALSE)

# Check that samples in each group now share the same median
group_check <- lapply(levels(groups), function(g) {
  med <- apply(n_df_final[, groups == g], 2, median, na.rm = TRUE)
  return(med)
})
names(group_check) <- levels(groups)
group_check

groups <- recode(groups,
                 "group1" = "dss_hfid",
                 "group2" = "dss",
                 "group3" = "control",
                 "group4" = "hfid")

group_colors <- c(
  "control" = "lightblue",
  "hfid" = "lightgreen",
  "dss" = "lightcoral",
  "dss_hfid" = "gold"
)
plot_colors <- group_colors[as.character(groups)]

pdf("01_normalized_box_plot.pdf", width = 12, height = 6)
par(mfrow = c(1,2), mar = c(8, 4, 3, 2))

# Before
boxplot(df_log2_mat,
        main = "Before Normalization",
        las = 2,
        col = plot_colors,
        ylab = "Log2 Intensity",
        cex.axis = 0.7)

# After
boxplot(n_df_final,
        main = "After Sample-Median Normalization (within groups)",
        las = 2,
        col = plot_colors,
        ylab = "Log2 Intensity",
        cex.axis = 0.7)

dev.off()

# -------------------------
# --- Density plots for before/after normalization
# -------------------------
# Colors by group
group_colors <- c("control" = "lightblue",
                  "hfid" = "lightgreen",
                  "dss" = "lightcoral",
                  "dss_hfid" = "gold")


sample_colors <- group_colors[as.character(groups)]

pdf("02_normalization_density.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))

# --- BEFORE normalization ---
d_first <- density(df_log2_mat[,1], na.rm = TRUE)
plot(d_first, main = "Before Normalization",
     xlim = range(df_log2_mat, na.rm = TRUE),
     ylim = c(0, max(sapply(1:ncol(df_log2_mat), function(i) max(density(df_log2_mat[,i], na.rm=TRUE)$y)))),
     col = sample_colors[1], lwd = 2,
     xlab = "Log2 Intensity", ylab = "Density")

for (i in 2:ncol(df_log2_mat)) {
  lines(density(df_log2_mat[,i], na.rm = TRUE), col = sample_colors[i], lwd = 2)
}

# --- AFTER normalization ---
d_first <- density(n_df_final[,1], na.rm = TRUE)
plot(d_first, main = "After Sample-Median Normalization (within groups)",
     xlim = range(n_df_final, na.rm = TRUE),
     ylim = c(0, max(sapply(1:ncol(n_df_final), function(i) max(density(n_df_final[,i], na.rm=TRUE)$y)))),
     col = sample_colors[1], lwd = 2,
     xlab = "Log2 Intensity", ylab = "Density")

for (i in 2:ncol(n_df_final)) {
  lines(density(n_df_final[,i], na.rm = TRUE), col = sample_colors[i], lwd = 2)
}

dev.off()

######################################################################################
#make a heatmap of all proteins (log transformed and normalized)
######################################################################################
library(pheatmap)
library(readr)
library(dplyr)
library(pheatmap)
library(viridis)
library(readxl)

# 1. Read data
crosstab <- read_excel("normalized_log2_data.xlsx")

# 2. Keep Protein + LI samples
df <- crosstab %>%
  select(Protein, starts_with("LI")) %>%
  as.data.frame()

# 3. Set rownames
rownames(df) <- df$Protein
df$Protein <- NULL

# 4. Force numeric conversion
df[] <- lapply(df, as.numeric)

# replace Inf/NaN with NA
mat <- as.matrix(df)
mat[!is.finite(mat)] <- NA

# 5. Align metadata
meta_li <- meta[meta$sample_id %in% colnames(mat), , drop = FALSE]
meta_li <- meta_li[match(colnames(mat), meta_li$sample_id), , drop = FALSE]

# 6. Create annotation dataframe for group colors
anno <- data.frame(Group = meta_li$group)
rownames(anno) <- meta_li$sample_id

# Make Group a factor with desired order
anno$Group <- factor(anno$Group,
                     levels = c("control", "hfid", "dss", "dss_hfid"))

# Order samples by this factor
ord <- order(anno$Group)
mat_ord <- mat[, ord]
anno_ord <- anno[ord, , drop = FALSE]

# Optional: assign custom colors
group_colors <- list(
  Group = c(
    control = "lightblue",
    hfid    = "lightcoral",
    dss     = "lightgreen",
    dss_hfid = "pink"
  )
)

# Plot
png(
  filename = "log2_norm_all_proteins_pl_z.png",
  width = 3000,
  height = 4200,
  res = 300
)

pheatmap(
  mat_ord,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = anno_ord,
  annotation_colors = group_colors,
  show_rownames = FALSE,
  na_col = "black",
  border_color = NA,
  color = viridis(100, option = "plasma"),
  main = "All proteins - log2 intensities"
)

dev.off()

##############################################################################################
#PcoA plot - not imputed
##############################################################################################
library(dplyr)
library(vegan)
library(ggplot2)

# ---------------------------
# Prepare matrix
# ---------------------------

# Select sample columns (starting with LI)
sample_cols <- colnames(n_df_final_id_full)[grepl("^LI", colnames(n_df_final_id_full))]

# Make matrix with proteins as rows, samples as columns
mat <- n_df_final_id_full %>%
  select(all_of(sample_cols), Protein) %>%  # ensure Protein column is included
  column_to_rownames(var = "Protein") %>%
  as.matrix()

# Force numeric type
storage.mode(mat) <- "numeric"

# Replace NAs with 0 (Bray-Curtis requires no NAs)
mat[is.na(mat)] <- 0

# Transpose: rows = samples, columns = proteins
pcoa_data <- t(mat)

# ---------------------------
# Compute Bray-Curtis distance
# ---------------------------
dist_mat <- vegdist(pcoa_data, method = "bray")

# ---------------------------
# Run PCoA
# ---------------------------
pcoa_res <- cmdscale(dist_mat, eig = TRUE, k = 2)

# Convert to data frame
pcoa_df <- as.data.frame(pcoa_res$points)
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$sample_id <- rownames(pcoa_df)

# % variance explained
eig <- pcoa_res$eig
eig <- eig[eig > 0]  # only positive eigenvalues
var_PC1 <- round(100 * eig[1] / sum(eig), 1)
var_PC2 <- round(100 * eig[2] / sum(eig), 1)

# Join metadata
pcoa_df <- pcoa_df %>%
  left_join(meta %>% select(sample_id, group), by = "sample_id")

pcoa_df <- pcoa_df %>%
  mutate(
    group = factor(
      group,
      levels = c("control", "hfid", "dss", "dss_hfid"),
      labels = c("Control", "HFiD", "DSS", "DSS+HFiD")
    )
  )

group_colors <- c(
  "Control"   = "#0072B2",
  "HFiD"      = "#E69F00",
  "DSS"       = "#CC79A7",
  "DSS+HFiD"  = "#009E73"
)

# ---------------------------
# Plot PCoA
# ---------------------------
pdf("PCoA_BrayCurtis_log2_norm.pdf", width = 7, height = 6)

ggplot(pcoa_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  scale_color_manual(values = group_colors) +
  theme_classic() +
  labs(
    title = "PCoA (Bray-Curtis) normalized log2 proteins",
    x = paste0("PC1 (", var_PC1, "%)"),
    y = paste0("PC2 (", var_PC2, "%)")
  )

dev.off()


################################################################################################
#IMPUTE THE DATA
################################################################################################
library(tidyverse)
library(pheatmap)       
library(ggrepel)        
library(RColorBrewer)   
library(ggpubr)         
library(viridis)        
library(limma)
library(dplyr)
library(ggplot2)
library(tibble)
library(imputeLCMD)
library(vegan)
library(openxlsx)

# Make a copy to work on
imp_df <- as.matrix(n_df_final)
mode(imp_df) <- "numeric"

# -------------------------
#impute MAR
# -------------------------

for (g in levels(groups)) {
  
  samp_idx <- which(groups == g)
  
  # subset group matrix
  group_mat <- imp_df[, samp_idx, drop = FALSE]
  
  for (i in 1:nrow(group_mat)) {
    vals <- group_mat[i, ]
    n_missing <- sum(is.na(vals))
    n_total <- length(vals)
    detection_rate <- (n_total - n_missing) / n_total
    
    if (detection_rate >= 0.75) {
      # MAR: >=75% present → group mean
      group_mat[i, is.na(vals)] <- mean(vals, na.rm = TRUE)
      
    } else if (detection_rate >= 0.5 & detection_rate < 0.75) {
      # 50–75% present → split by intensity median
      present_vals <- vals[!is.na(vals)]
      if (length(present_vals) > 0) {
        intensity_median <- median(present_vals)
        # MAR: top 50% → group mean
        group_mat[i, is.na(vals) & present_vals >= intensity_median] <- mean(present_vals[present_vals >= intensity_median], na.rm = TRUE)
        # MNAR: bottom 50% → will fill later with MinProb
      }
      
    } 
    # else detection_rate <0.5 → MNAR → leave NA for now
  }
  
  # assign back to imp_df
  imp_df[, samp_idx] <- group_mat
}

# -------------------------
#Apply MinProb for remaining NAs (MNAR proteins)
# -------------------------
remaining_na <- sum(is.na(imp_df))
if (remaining_na > 0) {
  cat("Applying MinProb imputation for remaining NAs (MNAR)...\n")
  imp_df <- impute.MinProb(imp_df, q = 0.01)
}
#---------------
#plot to check
#--------------
library(dplyr)
library(tidyr)
library(ggplot2)

# Before imputation
df_long_raw <- as.data.frame(n_df_final) %>%
  rownames_to_column("protein") %>%
  pivot_longer(-protein, names_to = "sample_id", values_to = "intensity") %>%
  left_join(meta %>% select(sample_id, group), by = "sample_id") %>%
  mutate(status = "raw")

# After imputation
df_long_imp <- as.data.frame(imp_df) %>%
  rownames_to_column("protein") %>%
  pivot_longer(-protein, names_to = "sample_id", values_to = "intensity") %>%
  left_join(meta %>% select(sample_id, group), by = "sample_id") %>%
  mutate(status = "imputed")

# Combine for plotting
df_long_plot <- bind_rows(df_long_raw, df_long_imp)

ggplot(df_long_plot, aes(x = intensity, color = status, fill = status)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ group, scales = "free") +
  scale_x_continuous(trans = "log10") +
  theme_classic() +
  labs(
    title = "Raw vs Imputed Intensities by Group",
    x = "Protein intensity (log10)",
    y = "Density",
    color = "Data",
    fill = "Data"
  )

library(dplyr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplot2)

# Raw data
df_long_raw <- as.data.frame(n_df_final) %>%
  rownames_to_column("protein") %>%
  pivot_longer(-protein, names_to = "sample_id", values_to = "intensity") %>%
  left_join(meta %>% select(sample_id, group), by = "sample_id") %>%
  mutate(
    status = ifelse(is.na(intensity) | intensity == 0, "missing", "observed"),
    plot_intensity = ifelse(status == "missing", min(intensity, na.rm = TRUE)*0.1, intensity),
    dataset = "raw"
  )

# Imputed data
df_long_imp <- as.data.frame(imp_df) %>%
  rownames_to_column("protein") %>%
  pivot_longer(-protein, names_to = "sample_id", values_to = "intensity") %>%
  left_join(meta %>% select(sample_id, group), by = "sample_id") %>%
  mutate(
    status = "imputed",
    plot_intensity = intensity,
    dataset = "imputed"
  )

# Combine
df_long_plot <- bind_rows(df_long_raw, df_long_imp)

pdf("raw_vs_imputed_density_by_group.pdf", width = 10, height = 8)

ggplot(df_long_plot, aes(x = plot_intensity, fill = status, color = status)) +
  geom_density(alpha = 0.3, adjust = 1) +
  facet_wrap(~ group, scales = "free") +
  scale_x_continuous(trans = "log10") +
  scale_fill_manual(values = c("observed" = "#56B4E9", "missing" = "#E69F00", "imputed" = "#009E73")) +
  scale_color_manual(values = c("observed" = "#56B4E9", "missing" = "#E69F00", "imputed" = "#009E73")) +
  theme_classic() +
  labs(
    title = "Raw vs Imputed Intensities by Group",
    x = "Protein intensity (log10)",
    y = "Density",
    fill = "Data type",
    color = "Data type"
  )

dev.off()
# ---------------------------
# Save imputed data
# ---------------------------
# Restore rownames just in case
rownames(imp_df) <- rownames(df_log2_mat)

tsv_file  <- "imputed_data.tsv"
xlsx_file <- "imputed_data.xlsx"

imp_df_df <- as.data.frame(imp_df)

imp_df_df_id <- as.data.frame(imp_df) %>%
  rownames_to_column("Protein") %>%
  left_join(annotation, by = "Protein")

write.table(
  imp_df_df_id,
  tsv_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.xlsx(
  imp_df_df_id,
  xlsx_file,
  rowNames = FALSE
)
