library(dplyr)
library(purrr)  
library(readr)
library(tidyr)
library(phyloseq)
library(readxl)
library(ggplot2)
library(writexl)

#set working directory
setwd("/Users/kristinasorokolet/Library/CloudStorage/Box-Box/Wright Lab Operations/Student Folders/Kris/KS_2/Proteomics/GLOBAL/RANDOM_BINS/With_host/females/")

metadata_file <- "metadata.tsv"

meta <- read_tsv(metadata_file, show_col_types = FALSE)

meta <- meta %>%
  mutate(group = recode(group,
                        "group1" = "dss_hfid",
                        "group2" = "dss",
                        "group3" = "control",
                        "group4" = "hfid"))

####################################################################################
#mouse to bacterial protein ratio
####################################################################################
library(tidyverse)

ctab_file <- "F_w_host_ms9_1p_r_filtered75.csv"

df_raw <- read_csv(ctab_file, show_col_types = FALSE)

# Store first column name (protein annotation / taxonomy)
id_col <- names(df_raw)[1]

# Identify sample columns
samples_to_keep <- unique(meta$sample_id)

df <- df_raw %>%
  select(all_of(c(id_col, samples_to_keep)))

df <- df %>%
  mutate(
    Source = ifelse(grepl("MOUSE", .[[id_col]], ignore.case = TRUE),
                    "Mouse", "Bacteria")
  )

df_long <- df %>%
  pivot_longer(
    cols = all_of(samples_to_keep),
    names_to = "sample_id",
    values_to = "intensity"
  )

df_summary <- df_long %>%
  group_by(sample_id, Source) %>%
  summarise(total_intensity = sum(intensity, na.rm = TRUE), .groups = "drop")

df_rel <- df_summary %>%
  group_by(sample_id) %>%
  mutate(rel_abundance = total_intensity / sum(total_intensity)) %>%
  ungroup()

p <- ggplot(df_rel, aes(x = sample_id, y = rel_abundance, fill = Source)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Sample",
    y = "Relative Abundance",
    title = "Mouse vs Bacterial Protein Contribution"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = "mouse_vs_bacteria_relative_abundance.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)

df_rel <- df_rel %>%
  left_join(meta, by = "sample_id")

df_group_avg <- df_rel %>%
  group_by(group, Source) %>%
  summarise(
    mean_rel_abundance = mean(rel_abundance, na.rm = TRUE),
    .groups = "drop"
  )

df_group_avg <- df_group_avg %>%
  mutate(
    group = factor(group, 
                   levels = c("control", "hfid", "dss", "dss_hfid"))
  )

p_group <- ggplot(df_group_avg,
                  aes(x = group,
                      y = mean_rel_abundance,
                      fill = Source)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = NULL,
    y = "Relative abundance (%)",
    title = "Host vs microbial protein contribution (females)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    # Text hierarchy
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.title.x = element_text(size = 11),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 10),
    
    # Clean panel
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(size = 0.3, color = "grey85"),
    
    # Axis lines
    axis.line = element_line(color = "black", size = 0.4),
    
    # Legend
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10),
    
    # Spacing
    plot.margin = margin(5, 10, 5, 5)
  )

ggsave(
  "group_average_horizontal.png",
  plot = p_group,
  width = 5,
  height = 4,
  dpi = 300
)

##########################################################################################################################################################################################################
#start processing filtered crosstabs 
##########################################################################################################################################################################################################


####################################################################################
#VISUALIZE RAW DATA
####################################################################################

library(tidyverse)
library(tibble)

ctab_file <- "F_w_host_ms9_1p_r_filtered75.csv"

df_full <- read_csv(ctab_file, show_col_types = FALSE)

# -------------------------------------------------------------------
# Define annotation columns
# -------------------------------------------------------------------
annotation_cols <- c(
  "Cluster", "Gene", "Function", "Peptide.Number",
  "Unique Peptides.", "Shared Peptides.", "Peptides",
  "Flanked Peptides", "gene", "function", "KEGG_Pathway", "Peptide Number", 
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
  filename = "log2_norm_all_proteins_pl.png",
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
#PcoA plot - before imputation
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

# ---------------------------
# Define colors
# ---------------------------
group_colors <- c(
  "control" = "lightblue",
  "hfid" = "lightcoral",
  "dss" = "lightgreen",
  "dss_hfid" = "pink"
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


####################################################################################################
### KEYWORD HEATMAPS
####################################################################################################
#Immunoglobuling proteins all (not only significant)
library(dplyr)
library(pheatmap)
library(viridis)
library(stringr)
library(tibble)
library(openxlsx)

# ---------------------------
# Keywords
# ---------------------------
keywords <- c(
  "Immunoglobulin",
  "\\bIg\\b",
  "^IG[HKL]",
  "JCHAIN"
)

## ---------------------------
# Sample columns
# ---------------------------
sample_cols <- meta$sample_id

# ---------------------------
# Build numeric matrix
# ---------------------------
heat_mat_focus <- imp_df_df_id %>%
  dplyr::select(all_of(sample_cols)) %>%
  as.matrix()

rownames(heat_mat_focus) <- imp_df_df_id$Protein

# ---------------------------
# Protein annotation (SAFE: no match(), no rowname dependency)
# ---------------------------
row_annot <- imp_df_df_id %>%
  dplyr::select(Protein, Function) %>%
  mutate(
    label = ifelse(is.na(Function) | Function == "", Protein, Function)
  )

# ---------------------------
# Filter proteins of interest
# ---------------------------
keep_proteins <- row_annot %>%
  filter(
    str_detect(
      label,
      regex(paste(keywords, collapse = "|"), ignore_case = TRUE)
    )
  ) %>%
  pull(Protein)

heat_mat_focus <- heat_mat_focus[keep_proteins, , drop = FALSE]
row_annot <- row_annot %>% filter(Protein %in% keep_proteins)

# ---------------------------
# Z-score scaling by protein
# ---------------------------
heat_mat_scaled <- t(scale(t(heat_mat_focus)))

keep <- complete.cases(heat_mat_scaled)

heat_mat_scaled <- heat_mat_scaled[keep, , drop = FALSE]
row_annot <- row_annot[keep, ]

# ---------------------------
# Order samples
# ---------------------------
meta_ord <- meta %>% arrange(group)

mat_top <- heat_mat_scaled[, meta_ord$sample_id, drop = FALSE]

anno_ord <- meta_ord %>%
  dplyr::select(sample_id, group) %>%
  tibble::column_to_rownames("sample_id")

# ---------------------------
# Group colors
# ---------------------------
group_colors <- c(
  control = "#0072B2",
  hfid = "#E69F00",
  dss = "#CC79A7",
  dss_hfid = "#009E73"
)

annotation_colors <- list(
  group = group_colors
)

# ---------------------------
# EXPORT (FIXED AND STABLE)
# ---------------------------
export_df <- mat_top %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Protein") %>%
  dplyr::left_join(row_annot, by = "Protein") %>%
  dplyr::select(Protein, Function, label, everything())

sample_annot <- meta_ord %>%
  dplyr::select(sample_id, group)

write.xlsx(
  list(
    expression_matrix = export_df,
    sample_annotation = sample_annot,
    protein_annotation = row_annot
  ),
  file = "IG_heatmap_proteins_annotated.xlsx",
  overwrite = TRUE
)

# ---------------------------
# Heatmap
# ---------------------------
png(
  "immunoglobulin_proteins_heatmap_ALL.png",
  width=2800,
  height=3000,
  res=300
)

pheatmap(
  mat_top,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 0,
  treeheight_col = 8,
  annotation_col = anno_ord,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  fontsize = 14,
  fontsize_row = 8,
  na_col = "black",
  border_color = NA,
  color = viridis(100, option = "D"),
  main = "Immunoglobulin proteins (males)",
  legend = TRUE
)

dev.off()


#------------------------------------------------------------------------------------------------------------------
#globin ALL 
#------------------------------------------------------------------------------------------------------------------


library(dplyr)
library(pheatmap)
library(viridis)
library(stringr)
library(tibble)

# ---------------------------
# Keywords
# ---------------------------
keywords <- c("globin", "ferritin")

# ---------------------------
# Sample columns
# ---------------------------
sample_cols <- meta$sample_id

# ---------------------------
# Build numeric matrix
# ---------------------------
heat_mat_focus <- imp_df_df_id %>%
  dplyr::select(all_of(sample_cols)) %>%
  as.matrix()

rownames(heat_mat_focus) <- imp_df_df_id$Protein

# ---------------------------
# Protein annotation (SAFE: no match(), no rowname dependency)
# ---------------------------
row_annot <- imp_df_df_id %>%
  dplyr::select(Protein, Function) %>%
  mutate(
    label = ifelse(is.na(Function) | Function == "", Protein, Function)
  )

# ---------------------------
# Filter proteins of interest
# ---------------------------
keep_proteins <- row_annot %>%
  filter(
    str_detect(
      label,
      regex(paste(keywords, collapse = "|"), ignore_case = TRUE)
    )
  ) %>%
  pull(Protein)

heat_mat_focus <- heat_mat_focus[keep_proteins, , drop = FALSE]
row_annot <- row_annot %>% filter(Protein %in% keep_proteins)

# ---------------------------
# Z-score scaling by protein
# ---------------------------
heat_mat_scaled <- t(scale(t(heat_mat_focus)))

keep <- complete.cases(heat_mat_scaled)

heat_mat_scaled <- heat_mat_scaled[keep, , drop = FALSE]
row_annot <- row_annot[keep, ]

# ---------------------------
# Order samples
# ---------------------------
meta_ord <- meta %>% arrange(group)

mat_top <- heat_mat_scaled[, meta_ord$sample_id, drop = FALSE]

anno_ord <- meta_ord %>%
  dplyr::select(sample_id, group) %>%
  tibble::column_to_rownames("sample_id")

# ---------------------------
# Group colors
# ---------------------------
group_colors <- c(
  control = "#0072B2",
  hfid = "#E69F00",
  dss = "#CC79A7",
  dss_hfid = "#009E73"
)

annotation_colors <- list(
  group = group_colors
)

# ---------------------------
# EXPORT (FIXED AND STABLE)
# ---------------------------
export_df <- mat_top %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Protein") %>%
  dplyr::left_join(row_annot, by = "Protein") %>%
  dplyr::select(Protein, Function, label, everything())

sample_annot <- meta_ord %>%
  dplyr::select(sample_id, group)

write.xlsx(
  list(
    expression_matrix = export_df,
    sample_annotation = sample_annot,
    protein_annotation = row_annot
  ),
  file = "globin_heatmap_proteins_annotated.xlsx",
  overwrite = TRUE
)

# ---------------------------
# HEATMAP
# ---------------------------
png(
  "globin_proteins_heatmap_ALL.png",
  width = 2800,
  height = 800,
  res = 300
)

pheatmap(
  mat_top,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 0,
  treeheight_col = 8,
  annotation_col = anno_ord,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  fontsize_row = 8,
  fontsize = 14,
  na_col = "black",
  border_color = NA,
  color = viridis::viridis(100, option = "D"),
  main = "Globin proteins (males)",
  legend = TRUE
)

dev.off()
#------------------------------------------------------------------------------------------------------------------
#proteasome 
#------------------------------------------------------------------------------------------------------------------
library(dplyr)
library(pheatmap)
library(viridis)
library(stringr)
library(tibble)
library(openxlsx)

# ---------------------------
# Keyword
# ---------------------------
keywords <- c("proteasome")

# ---------------------------
# Sample columns
# ---------------------------
sample_cols <- meta$sample_id

# ---------------------------
# Build numeric matrix
# ---------------------------
heat_mat_focus <- imp_df_df_id %>%
  dplyr::select(all_of(sample_cols)) %>%
  as.matrix()

rownames(heat_mat_focus) <- imp_df_df_id$Protein

# ---------------------------
# Protein annotation (SAFE: no match(), no rowname dependency)
# ---------------------------
row_annot <- imp_df_df_id %>%
  dplyr::select(Protein, Function) %>%
  mutate(
    label = ifelse(is.na(Function) | Function == "", Protein, Function)
  )

# ---------------------------
# Filter proteins of interest
# ---------------------------
keep_proteins <- row_annot %>%
  filter(
    str_detect(
      label,
      regex(paste(keywords, collapse = "|"), ignore_case = TRUE)
    )
  ) %>%
  pull(Protein)

heat_mat_focus <- heat_mat_focus[keep_proteins, , drop = FALSE]
row_annot <- row_annot %>% filter(Protein %in% keep_proteins)

# ---------------------------
# Z-score scaling by protein
# ---------------------------
heat_mat_scaled <- t(scale(t(heat_mat_focus)))

keep <- complete.cases(heat_mat_scaled)

heat_mat_scaled <- heat_mat_scaled[keep, , drop = FALSE]
row_annot <- row_annot[keep, ]

# ---------------------------
# Order samples
# ---------------------------
meta_ord <- meta %>% arrange(group)

mat_top <- heat_mat_scaled[, meta_ord$sample_id, drop = FALSE]

anno_ord <- meta_ord %>%
  dplyr::select(sample_id, group) %>%
  tibble::column_to_rownames("sample_id")

# ---------------------------
# Group colors
# ---------------------------
group_colors <- c(
  control = "#0072B2",
  hfid = "#E69F00",
  dss = "#CC79A7",
  dss_hfid = "#009E73"
)

annotation_colors <- list(
  group = group_colors
)

# ---------------------------
# EXPORT (FIXED AND STABLE)
# ---------------------------
export_df <- mat_top %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Protein") %>%
  dplyr::left_join(row_annot, by = "Protein") %>%
  dplyr::select(Protein, Function, label, everything())

sample_annot <- meta_ord %>%
  dplyr::select(sample_id, group)

write.xlsx(
  list(
    expression_matrix = export_df,
    sample_annotation = sample_annot,
    protein_annotation = row_annot
  ),
  file = "proteasome_heatmap_proteins_annotated.xlsx",
  overwrite = TRUE
)

# ---------------------------
# HEATMAP (labels ONLY for visualization)
# ---------------------------
mat_plot <- mat_top

# ensure alignment before relabeling
mat_plot <- mat_plot[row_annot$Protein, , drop = FALSE]

rownames(mat_plot) <- row_annot$label

png(
  "proteasome_proteins_heatmap_ALL.png",
  width = 2800,
  height = 2300,
  res = 300
)

pheatmap(
  mat_plot,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 0,
  treeheight_col = 8,
  annotation_col = anno_ord,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  fontsize_row = 8,
  fontsize = 14,
  na_col = "black",
  border_color = NA,
  color = viridis::viridis(100, option = "D"),
  main = "Proteasome proteins",
  legend = TRUE
)

dev.off()

####################################################################################
#keyword box plots
####################################################################################
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(rstatix)
library(ggpubr)

crosstab <- read_tsv("imputed_data.tsv")


# Identify sample columns (all LI samples)
sample_cols <- colnames(crosstab)[grepl("^LI", colnames(crosstab))]

expr_matrix <- crosstab %>%
  select(Protein, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), ~ as.numeric(.))) %>%
  column_to_rownames("Protein") %>%
  as.matrix()

ctab_file <- "F_w_host_ms9_1p_r_filtered75.csv" 
ctab_orig <- read_csv(ctab_file, show_col_types = FALSE)

prot_function <- ctab_orig %>%
  select(Protein = 1, Function) %>%
  mutate(
    Protein = as.character(Protein),
    Function = as.character(Function)
  )


keywords <- c("ferritin", "globin", "immunoglobulin", "\\bIg\\b")

row_annot <- data.frame(Protein = rownames(expr_matrix)) %>%
  left_join(prot_function, by = "Protein") %>%
  mutate(label = ifelse(is.na(Function), Protein, Function))

row_annot_kw <- row_annot %>%
  mutate(keyword = case_when(
    str_detect(label, regex("ferritin", ignore_case = TRUE)) ~ "Ferritin",
    str_detect(label, regex("globin", ignore_case = TRUE)) ~ "Globin",
    str_detect(label, regex("immunoglobulin|\\bIg\\b", ignore_case = TRUE)) ~ "Immunoglobulin",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(keyword))

expr_sub <- expr_matrix[rownames(expr_matrix) %in% row_annot_kw$Protein, ]

expr_long <- as.data.frame(expr_sub) %>%
  tibble::rownames_to_column("Protein") %>%
  pivot_longer(
    cols = -Protein,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  left_join(meta, by = "sample_id") %>%
  left_join(row_annot_kw %>% select(Protein, keyword), by = "Protein")

stat_tests <- expr_long %>%
  dplyr::group_by(keyword) %>%
  rstatix::wilcox_test(abundance ~ group) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj")

stat_tests <- expr_long %>%
  group_by(keyword) %>%
  pairwise_wilcox_test(abundance ~ group, p.adjust.method = "BH")

p <- ggplot(expr_long, aes(x = group, y = abundance, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
  facet_wrap(~ keyword, scales = "free_y") +
  theme_classic() +
  labs(
    title = "Protein abundance by keyword",
    x = NULL,
    y = "Expression"
  )

library(dplyr)

ypos <- expr_long %>%
  group_by(keyword) %>%
  summarise(max_y = max(abundance, na.rm = TRUE)) %>%
  mutate(y.position = max_y * 1.1)

library(rstatix)

stat_tests <- expr_long %>%
  group_by(keyword) %>%
  pairwise_wilcox_test(abundance ~ group, p.adjust.method = "BH")

stat_tests <- stat_tests %>%
  left_join(ypos %>% select(keyword, y.position), by = "keyword")

p + stat_pvalue_manual(
  stat_tests,
  label = "p.adj.signif",
  hide.ns = TRUE
)

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(rstatix)

# ---------------------------
# 1. Color palette (FIXED FORMAT)
# ---------------------------
group_colors <- c(
  control  = "#0072B2",
  hfid     = "#E69F00",
  dss      = "#CC79A7",
  dss_hfid = "#009E73"
)

# ---------------------------
# 2. Protein annotation + keyword mapping
# ---------------------------
row_annot <- data.frame(Protein = rownames(expr_matrix)) %>%
  left_join(prot_function, by = "Protein") %>%
  mutate(label = ifelse(is.na(Function), Protein, Function))

row_annot_kw <- row_annot %>%
  mutate(keyword = case_when(
    str_detect(label, regex("ferritin|globin", ignore_case = TRUE)) ~ "Ferritin/Globin",
    str_detect(label, regex("immunoglobulin|\\bIg\\b", ignore_case = TRUE)) ~ "Immunoglobulin",
    str_detect(label, regex("proteasome|psm[a-z]|psmb|psma", ignore_case = TRUE)) ~ "Proteasome",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(keyword))

# ---------------------------
# 3. Subset expression matrix
# ---------------------------
expr_sub <- expr_matrix[rownames(expr_matrix) %in% row_annot_kw$Protein, ]

# ---------------------------
# 4. Long format
# ---------------------------
expr_long <- as.data.frame(expr_sub) %>%
  tibble::rownames_to_column("Protein") %>%
  pivot_longer(
    cols = -Protein,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  left_join(meta, by = "sample_id") %>%
  left_join(row_annot_kw %>% select(Protein, keyword), by = "Protein")

expr_long <- expr_long %>%
  mutate(group = factor(
    group,
    levels = c("control", "hfid", "dss", "dss_hfid")
  ))

# ---------------------------
# 5. Plot function (CORE FIX)
# ---------------------------
plot_keyword <- function(data, kw) {
  
  df <- data %>% filter(keyword == kw)
  
  stat <- df %>%
    pairwise_wilcox_test(abundance ~ group, p.adjust.method = "BH") %>%
    add_xy_position(x = "group", step.increase = 0.2)
  
  ggplot(df, aes(x = group, y = abundance, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.5) +
    scale_fill_manual(values = group_colors) +
    stat_pvalue_manual(
      stat,
      label = "p.adj.signif",
      hide.ns = TRUE
    ) +
    theme_classic() +
    labs(
      title = kw,
      x = NULL,
      y = "Log2 expression"
    ) +
    theme(
      axis.text.x = element_text(size = 14),   # group labels
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      plot.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16, face = "bold")
    )
}
# ---------------------------
# 6. Generate plots
# ---------------------------
p_fg <- plot_keyword(expr_long, "Ferritin/Globin")
p_ig <- plot_keyword(expr_long, "Immunoglobulin")
p_prot <- plot_keyword(expr_long, "Proteasome")

# ---------------------------
# 7. Display
# ---------------------------
# Ferritin + Globin combined
ggsave(
  filename = "Ferritin_Globin_boxplot.png",
  plot = p_fg,
  width = 6,
  height = 7,
  dpi = 300
)

# Immunoglobulin
ggsave(
  filename = "Immunoglobulin_boxplot.png",
  plot = p_ig,
  width = 6,
  height = 8,
  dpi = 300
)

ggsave(
  filename = "Proteasome_boxplot.png",
  plot = p_prot,
  width = 6,
  height = 7,
  dpi = 300
)


####################################################################################
#trying to understand if dss or hfid has a stronger effect
####################################################################################
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)

df <- read_csv("normalized_log2_data.csv")

df <- df %>%
  as.data.frame()

rownames(df) <- df$Protein
df$Protein <- NULL

# treat 0 as not detected
df[df == 0] <- NA

pa <- !is.na(df)   # TRUE = detected, FALSE = not detected

pa_long <- pa %>%
  as.data.frame() %>%
  rownames_to_column("Protein") %>%
  pivot_longer(-Protein,
               names_to = "sample",
               values_to = "detected") %>%
  left_join(meta %>% select(sample_id, group),
            by = c("sample" = "sample_id"))

group_presence <- pa_long %>%
  group_by(Protein, group) %>%
  summarize(
    detection_rate = mean(detected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  distinct() %>%   # critical safety step
  pivot_wider(
    names_from = group,
    values_from = detection_rate,
    values_fill = 0   # prevents list-columns
  )

group_presence <- group_presence %>%
  mutate(across(where(is.list), ~unlist(.)))

presence_thresh <- 0.5

df_cat <- group_presence %>%
  mutate(
    DSS_HFID_present = ifelse(is.na(dss_hfid), FALSE, dss_hfid >= presence_thresh),
    DSS_present      = ifelse(is.na(dss), FALSE, dss >= presence_thresh),
    HFID_present     = ifelse(is.na(hfid), FALSE, hfid >= presence_thresh),
    
    category = case_when(
      DSS_HFID_present & !DSS_present & !HFID_present ~ "Unique DSS+HFiD",
      DSS_HFID_present & DSS_present & !HFID_present ~ "Shared with DSS",
      DSS_HFID_present & HFID_present & !DSS_present ~ "Shared with HFiD",
      TRUE ~ "Other"
    )
  )

count_df <- df_cat %>%
  dplyr::mutate(category = as.character(category)) %>%
  dplyr::count(category) %>%
  dplyr::mutate(perc = n / sum(n) * 100)

count_df$category <- factor(
  count_df$category,
  levels = c("Unique DSS+HFiD",
             "Shared with DSS",
             "Shared with HFiD",
             "Other")
)
colors <- c(
  "Unique DSS+HFiD" = "#E15759",
  "Shared with DSS" = "#4E79A7",
  "Shared with HFiD" = "#F28E2B",
  "Other" = "#BAB0AC"
)

png("DSS_HFID_protein_overlap_0.5.png",
    width = 1400, height = 1600, res = 300)

ggplot(count_df, aes(x = "DSS+HFID proteome", y = perc, fill = category)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  geom_text(aes(label = paste0(round(perc, 1), "%")),
            position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(values = colors) +
  labs(
    x = NULL,
    y = "Percentage of proteins",
    fill = "Category",
    title = "Host proteome overlap (males)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

dev.off()


# what are those proteins?
library(openxlsx)
library(dplyr)
library(readr)


# ---------------------------
# 1. Build annotated protein table
# ---------------------------
protein_table <- df_cat %>%
  left_join(
    crosstab %>% dplyr::select(Protein, Function),
    by = "Protein"
  )

# ---------------------------
# 2. Split categories
# ---------------------------
unique_dss_hfid <- protein_table %>%
  filter(category == "Unique DSS+HFiD")

shared_dss <- protein_table %>%
  filter(category == "Shared with DSS")

shared_hfid <- protein_table %>%
  filter(category == "Shared with HFiD")

# ---------------------------
# 3. Export Excel
# ---------------------------
wb <- createWorkbook()

addWorksheet(wb, "Unique_DSS_HFiD")
writeData(wb, "Unique_DSS_HFiD", unique_dss_hfid)

addWorksheet(wb, "Shared_DSS")
writeData(wb, "Shared_DSS", shared_dss)

addWorksheet(wb, "Shared_HFiD")
writeData(wb, "Shared_HFiD", shared_hfid)

saveWorkbook(
  wb,
  "DSS_HFID_protein_categories_with_function_0.5.xlsx",
  overwrite = TRUE
)

# ---------------------------
# 4. Subset proteins for heatmap
# ---------------------------
proteins_cat <- unique(unique_dss_hfid$Protein)

imp_cat <- imp_df[rownames(imp_df) %in% proteins_cat, ]

imp_cat_df <- as.data.frame(imp_cat)
imp_cat_df$Protein <- rownames(imp_cat_df)

# attach Function annotation (ONLY reliable source)
imp_cat_df <- imp_cat_df %>%
  left_join(
    crosstab %>% dplyr::select(Protein, Function),
    by = "Protein"
  )

# ---------------------------
# 5. Build matrix
# ---------------------------
mat <- imp_cat_df %>%
  column_to_rownames("Protein") %>%
  dplyr::select(-Function) %>%
  as.matrix()

# ---------------------------
# 6. Z-score scaling
# ---------------------------
mat_z <- t(scale(t(mat)))

# ---------------------------
# 7. Row labels = Function (fallback to Protein)
# ---------------------------
row_labels <- imp_cat_df$Function
row_labels[is.na(row_labels)] <- imp_cat_df$Protein

rownames(mat_z) <- make.unique(row_labels)

# ---------------------------
# 8. Column ordering + annotation
# ---------------------------
meta_ord <- meta %>% arrange(group)

mat_z <- mat_z[, meta_ord$sample_id]

anno_col <- meta_ord %>%
  dplyr::select(sample_id, group) %>%
  tibble::column_to_rownames("sample_id")

# ---------------------------
# 9. Colors
# ---------------------------
group_colors <- c(
  control  = "#0072B2",
  hfid     = "#E69F00",
  dss      = "#CC79A7",
  dss_hfid = "#009E73"
)

annotation_colors <- list(
  group = group_colors
)

# ---------------------------
# 10. Heatmap
# ---------------------------
png("unique_dss_hfid_heatmap_0.5.png",
    width = 3000, height =5000, res = 300)

pheatmap(
  mat_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 0,
  treeheight_col = 8,
  
  annotation_col = anno_col,
  annotation_colors = annotation_colors,
  
  show_rownames = TRUE,
  fontsize_row = 8,
  na_col = "black",
  border_color = NA,
  color = viridis::viridis(100, option = "D"),
  
  main = "Proteins unique to DSS+HFiD (males)"
)

dev.off()

