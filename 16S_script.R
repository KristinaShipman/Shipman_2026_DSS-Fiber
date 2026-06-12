library(dplyr)
library(purrr) 
library(readr)
library(tidyr)
library(phyloseq)
library(readxl)

base_dir <- "/Users/kristinasorokolet/Library/CloudStorage/Box-Box/Wright Lab Operations/Student Folders/Kris/KS_2/16S/RUN1/emu_count_SIMIGTv2/"
setwd("/Users/kristinasorokolet/Library/CloudStorage/Box-Box/Wright Lab Operations/Student Folders/Kris/KS_2/16S/RUN1/emu_count_SIMIGTv2/")

# List all EMU output files
files <- list.files(base_dir, pattern = ".fastq_rel-abundance.tsv$", full.names = TRUE)
samples <- gsub(".fastq_rel-abundance.tsv", "", basename(files))
meta <- read.csv("/Users/kristinasorokolet/Library/CloudStorage/Box-Box/Wright Lab Operations/Student Folders/Kris/KS_2/16S/RUN1/emu_count_SIMIGTv2/KS2_metadata.csv", header = TRUE,
                 stringsAsFactors = FALSE,
                 colClasses = "character") 

##############################################################################################
# Read all files into a single long-format data frame, remove duplicates, and anything below 0.1%
##############################################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(openxlsx)

# -----------------------------
# Step 1: Split samples by Sex
# -----------------------------
male_samples   <- meta$SampleID[meta$Sex == "M"]

female_samples <- meta$SampleID[meta$Sex == "F"]

# -----------------------------
# Step 2: Process Males
# -----------------------------
# Read all files for males into a long-format data frame
long_df1_m <- map2_dfr(files, samples, function(file, samp) {
  if(samp %in% male_samples){
    read_tsv(file, col_types = cols()) %>%
      select(tax_id, lineage, 'estimated counts') %>%
      mutate(sample = samp) %>%
      rename(counts = 'estimated counts')
  } else {
    NULL
  }
})

# save raw read data
wide_df1_m <- long_df1_m %>%
  pivot_wider(
    names_from = sample,
    values_from = counts,
    values_fill = 0,
    values_fn = sum
  )

# Save as Excel
write.xlsx(
  wide_df1_m,
  file = "raw_read_data_m.xlsx",
  overwrite = TRUE
)

# Remove unclassified/unmapped
long_df2_m <- long_df1_m %>% 
  filter(!tax_id %in% c("mapped_unclassified", "unmapped"))

# Replace NA with 0
long_df2_m$counts[is.na(long_df2_m$counts)] <- 0

# Keep only one row per tax_id per sample
long_df3_m <- long_df2_m %>%
  distinct(tax_id, lineage, sample, .keep_all = TRUE)

# Convert to wide format
merged_df_m <- long_df3_m %>%
  pivot_wider(names_from = sample, values_from = counts, values_fill = 0) %>%
  arrange(tax_id) %>%
  distinct(tax_id, .keep_all = TRUE)

# Create OTU matrix
otu_mat_m <- as.matrix(merged_df_m %>% select(-tax_id, -lineage))
rownames(otu_mat_m) <- merged_df_m$tax_id

# Collapse low abundance taxa into "Other"
threshold <- 0.001
rel_abund_mat_m <- sweep(otu_mat_m, 2, colSums(otu_mat_m), FUN = "/")
low_abundance_m <- rowSums(rel_abund_mat_m >= threshold) == 0
other_row_m <- colSums(otu_mat_m[low_abundance_m, , drop = FALSE])
otu_mat_filtered_m <- otu_mat_m[!low_abundance_m, , drop = FALSE]
otu_mat_filtered_m <- rbind(otu_mat_filtered_m, Other = other_row_m)

merged_df_filtered_m <- merged_df_m[!low_abundance_m, ]
other_df_m <- data.frame(
  tax_id = "Other",
  lineage = "Other",
  t(other_row_m),
  check.names = FALSE
)
merged_df_final_m <- bind_rows(merged_df_filtered_m, other_df_m)

# Save as Excel
write.xlsx(
  merged_df_final_m,
  file = "clean_reads_m.xlsx",
  overwrite = TRUE
)

# -----------------------------
# Step 2: Process Females
# -----------------------------
long_df1_f <- map2_dfr(files, samples, function(file, samp) {
  if(samp %in% female_samples){
    read_tsv(file, col_types = cols()) %>%
      select(tax_id, lineage, 'estimated counts') %>%
      mutate(sample = samp) %>%
      rename(counts = 'estimated counts')
  } else {
    NULL
  }
})

# save females raw read data 
wide_df1_f <- long_df1_f %>%
  pivot_wider(
    names_from = sample,
    values_from = counts,
    values_fill = 0,
    values_fn = sum
  )

# Save as Excel
write.xlsx(
  wide_df1_f,
  file = "raw_read_data_f.xlsx",
  overwrite = TRUE
)

long_df2_f <- long_df1_f %>% 
  filter(!tax_id %in% c("mapped_unclassified", "unmapped"))

long_df2_f$counts[is.na(long_df2_f$counts)] <- 0

long_df3_f <- long_df2_f %>%
  distinct(tax_id, lineage, sample, .keep_all = TRUE)

merged_df_f <- long_df3_f %>%
  pivot_wider(names_from = sample, values_from = counts, values_fill = 0) %>%
  arrange(tax_id) %>%
  distinct(tax_id, .keep_all = TRUE)

otu_mat_f <- as.matrix(merged_df_f %>% select(-tax_id, -lineage))
rownames(otu_mat_f) <- merged_df_f$tax_id

rel_abund_mat_f <- sweep(otu_mat_f, 2, colSums(otu_mat_f), FUN = "/")
low_abundance_f <- rowSums(rel_abund_mat_f >= threshold) == 0
other_row_f <- colSums(otu_mat_f[low_abundance_f, , drop = FALSE])
otu_mat_filtered_f <- otu_mat_f[!low_abundance_f, , drop = FALSE]
otu_mat_filtered_f <- rbind(otu_mat_filtered_f, Other = other_row_f)

merged_df_filtered_f <- merged_df_f[!low_abundance_f, ]
other_df_f <- data.frame(
  tax_id = "Other",
  lineage = "Other",
  t(other_row_f),
  check.names = FALSE
)
merged_df_final_f <- bind_rows(merged_df_filtered_f, other_df_f)

# Save as Excel
write.xlsx(
  merged_df_final_f,
  file = "clean_reads_f.xlsx",
  overwrite = TRUE
)


##############################################################################################
#separate lineage into ranks
##############################################################################################
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# -----------------------------
# Step 0: Define ranks
# -----------------------------
ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

# -----------------------------
# Step 1: Process Males
# -----------------------------
tax_mat_m <- merged_df_final_m %>%
  select(tax_id, lineage) %>%
  mutate(
    lineage_split = map(lineage, ~ {
      pieces <- strsplit(.x, ";")[[1]]
      pieces <- pieces[pieces != ""]       # remove empty strings
      vals <- tail(pieces, 7)              # last 7 taxonomic levels
      if (length(vals) < 7) vals <- c(rep(NA, 7 - length(vals)), vals)
      rev(vals)                            # ensure Kingdom first
    })
  ) %>%
  unnest_wider(lineage_split, names_sep = "_") %>%
  rename_with(~ ranks, starts_with("lineage_split_")) %>%
  mutate(across(everything(), trimws)) %>%
  as.data.frame()

# Remove any old lowercase "other" entries
tax_mat_m <- tax_mat_m %>%
  filter(!tolower(tax_id) == "other")

# Add single "Other" row
other_row_m <- data.frame(
  tax_id = "Other",
  matrix("Other", nrow = 1, ncol = length(ranks))
)
colnames(other_row_m)[-1] <- ranks

tax_mat_m <- bind_rows(tax_mat_m, other_row_m)

# Add " sp." where appropriate 
tax_mat_m <- tax_mat_m %>%
  mutate(Species = ifelse(
    !is.na(Species) & 
      !str_detect(Species, "\\s") & 
      Species != "Other",
    paste0(Species, " sp."),
    Species
  ))

# Assign rownames
rownames(tax_mat_m) <- merged_df_final_m$tax_id
rownames(otu_mat_filtered_m) <- merged_df_final_m$tax_id

tax_mat_phy_m <- as.matrix(tax_mat_m)

# -----------------------------
# Step 2: Process Females
# -----------------------------
tax_mat_f <- merged_df_final_f %>%
  select(tax_id, lineage) %>%
  mutate(
    lineage_split = map(lineage, ~ {
      pieces <- strsplit(.x, ";")[[1]]
      pieces <- pieces[pieces != ""]       # remove empty strings
      vals <- tail(pieces, 7)              # last 7 taxonomic levels
      if (length(vals) < 7) vals <- c(rep(NA, 7 - length(vals)), vals)
      rev(vals)                            # ensure Kingdom first
    })
  ) %>%
  unnest_wider(lineage_split, names_sep = "_") %>%
  rename_with(~ ranks, starts_with("lineage_split_")) %>%
  mutate(across(everything(), trimws)) %>%
  as.data.frame()

tax_mat_f <- tax_mat_f %>%
  filter(!tolower(tax_id) == "other")

other_row_f <- data.frame(
  tax_id = "Other",
  matrix("Other", nrow = 1, ncol = length(ranks))
)
colnames(other_row_f)[-1] <- ranks

tax_mat_f <- bind_rows(tax_mat_f, other_row_f)

tax_mat_f <- tax_mat_f %>%
  mutate(Species = ifelse(
    !is.na(Species) & 
      !str_detect(Species, "\\s") & 
      Species != "Other",
    paste0(Species, " sp."),
    Species
  ))

rownames(tax_mat_f) <- merged_df_final_f$tax_id
rownames(otu_mat_filtered_f) <- merged_df_final_f$tax_id

tax_mat_phy_f <- as.matrix(tax_mat_f)

##############################################################################################
#now prep metadata for the phyloseq object and then merge metadata, otu table and taxonomy table into a phyloseq object
##############################################################################################

library(readr)
library(dplyr)
library(tibble)
library(phyloseq)

# -----------------------------
# Step 0: Prepare metadata for males
# -----------------------------
meta_m <- meta %>%
  filter(SampleID %in% male_samples) %>%   # subset to male samples
  column_to_rownames("SampleID")           # set SampleID as rownames

# Check order matches OTU matrix columns
all(rownames(meta_m) %in% colnames(otu_mat_filtered_m))  # should be TRUE

# Reorder metadata to match OTU columns
meta_m <- meta_m[colnames(otu_mat_filtered_m), , drop = FALSE]

# -----------------------------
# Step 1: Create phyloseq object for males
# -----------------------------
OTU_m <- otu_table(otu_mat_filtered_m, taxa_are_rows = TRUE)
TAX_m <- tax_table(tax_mat_phy_m)
META_m <- sample_data(meta_m)

ps_M <- phyloseq(OTU_m, TAX_m, META_m)

# Quick checks
ps_M
sample_names(ps_M)[1:5]
taxa_names(ps_M)[1:5]
sample_data(ps_M)[1:5, ]
ntaxa(ps_M)

# -----------------------------
# Step 2: Prepare metadata for females
# -----------------------------
meta_f <- meta %>%
  filter(SampleID %in% female_samples) %>% # subset to female samples
  column_to_rownames("SampleID")

# Check order matches OTU matrix columns
all(rownames(meta_f) %in% colnames(otu_mat_filtered_f))  # should be TRUE

# Reorder metadata to match OTU columns
meta_f <- meta_f[colnames(otu_mat_filtered_f), , drop = FALSE]

# -----------------------------
# Step 3: Create phyloseq object for females
# -----------------------------
OTU_f <- otu_table(otu_mat_filtered_f, taxa_are_rows = TRUE)
TAX_f <- tax_table(tax_mat_phy_f)
META_f <- sample_data(meta_f)

ps_F <- phyloseq(OTU_f, TAX_f, META_f)

# Quick checks
ps_F
sample_names(ps_F)[1:5]
taxa_names(ps_F)[1:5]
sample_data(ps_F)[1:5, ]
ntaxa(ps_F)

########################################################################################################################
#cleanup ambiguous taxa names
########################################################################################################################
library(phyloseq)
library(dplyr)

# -----------------------------
# Step 0: Define the taxa to recode as "Other"
# -----------------------------
weird_species <- c(
  "uncultured rumen bacterium",
  "unidentified rumen bacterium",
  "unidentified sp.",
  "uncultured bacterium"
)

# -----------------------------
# Step 1: Process males
# -----------------------------
# Convert taxonomy table to data frame
tax_m <- as.data.frame(tax_table(ps_M))

# Replace weird species with "Other"
tax_m$Species <- ifelse(tax_m$Species %in% weird_species, "Other", tax_m$Species)

# Update taxonomy table in phyloseq object
tax_table(ps_M) <- as.matrix(tax_m)

# Merge all OTUs/ASVs with the same Species name
ps_M_c <- tax_glom(ps_M, taxrank = "Species")

# -----------------------------
# Step 2: Process females
# -----------------------------
tax_f <- as.data.frame(tax_table(ps_F))
tax_f$Species <- ifelse(tax_f$Species %in% weird_species, "Other", tax_f$Species)
tax_table(ps_F) <- as.matrix(tax_f)
ps_F_c <- tax_glom(ps_F, taxrank = "Species")

#############################################################################
#RAREFICATION
#############################################################################
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

#############################################################################
# check read depth
#############################################################################

# females
sums_f <- sample_sums(ps_F_c)

summary(sums_f)
min(sums_f)
median(sums_f)
quantile(sums_f, probs = c(0.05, 0.10, 0.25))

hist(
  sums_f,
  breaks = 50,
  main = "Female sequencing depth",
  xlab = "Reads"
)

# males
sums_m <- sample_sums(ps_M_c)

summary(sums_m)

hist(
  sums_m,
  breaks = 50,
  main = "Male sequencing depth",
  xlab = "Reads"
)

#############################################################################
# make read depth plots 
#############################################################################

make_reads_plot <- function(ps_obj, group_var, output_name){
  
  # sequencing depth
  depth_df <- data.frame(
    SampleID = sample_names(ps_obj),
    Reads = sample_sums(ps_obj)
  )
  
  # metadata
  meta_df <- as.data.frame(sample_data(ps_obj))
  meta_df$SampleID <- rownames(meta_df)
  
  # merge
  depth_df <- left_join(depth_df, meta_df, by = "SampleID")
  
  # plot
  p <- ggplot(depth_df,
              aes_string(x = group_var,
                         y = "Reads",
                         fill = group_var)) +
    
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    
    geom_jitter(
      aes_string(color = group_var),
      width = 0.2,
      size = 2,
      alpha = 0.8
    ) +
    
    theme_bw() +
    
    labs(
      title = "Sequencing Depth per Sample",
      x = group_var,
      y = "Read Counts"
    ) +
    
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +
    
    scale_y_continuous(
      labels = comma,
      breaks = seq(0, max(depth_df$Reads), by = 20000)
    ) +
    
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  ggsave(output_name, plot = p, width = 8, height = 6, dpi = 600)
  
  return(p)
}

#############################################################################
# plot
#############################################################################

# replace "Treatment" with your metadata column if different

reads_plot_f <- make_reads_plot(
  ps_obj = ps_F_c,
  group_var = "Treatment",
  output_name = "female_reads_plot.png"
)

reads_plot_m <- make_reads_plot(
  ps_obj = ps_M_c,
  group_var = "Treatment",
  output_name = "male_reads_plot.png"
)

#############################################################################
# Observed taxa vs. reads rarefication plots 
#############################################################################

library(phyloseq)
library(ggplot2)
library(dplyr)

#############################################################################
# make a function for rarefication curves
#############################################################################

rarefaction_points <- function(x, step = 1000, reps = 5) {
  
  x <- round(x[x > 0])   # CRITICAL FIX
  max_depth <- sum(x)
  
  depths <- seq(step, max_depth, by = step)
  
  out <- lapply(depths, function(d) {
    
    # hard safety guard
    if (d > max_depth) return(NULL)
    
    richness_reps <- replicate(reps, {
      
      pool <- rep(seq_along(x), x)
      
      # extra safety check
      if (length(pool) < d) return(NA)
      
      sampled <- sample(pool, d, replace = FALSE)
      length(unique(sampled))
    })
    
    data.frame(
      Depth = d,
      Observed = mean(richness_reps, na.rm = TRUE)
    )
  })
  
  do.call(rbind, out)
}

run_rarefaction <- function(ps_obj, step = 1000) {
  
  otu_mat <- as(otu_table(ps_obj), "matrix")
  
  if (taxa_are_rows(ps_obj)) {
    otu_mat <- t(otu_mat)
  }
  
  # CRITICAL FIX: enforce integer counts
  otu_mat <- round(otu_mat)
  
  otu_df <- as.data.frame(otu_mat)
  otu_df$Sample <- rownames(otu_df)
  
  rarefaction_df <- lapply(1:nrow(otu_df), function(i) {
    
    df <- rarefaction_points(
      as.numeric(otu_df[i, -ncol(otu_df)]),
      step = step
    )
    
    df$Sample <- otu_df$Sample[i]
    df
  }) %>% bind_rows()
  
  # attach metadata
  meta <- data.frame(sample_data(ps_obj))
  meta$Sample <- rownames(meta)
  
  left_join(rarefaction_df, meta, by = "Sample")
}

#############################################################################
# running the function
#############################################################################

rar_M <- run_rarefaction(ps_M_c, step = 1000)
rar_F <- run_rarefaction(ps_F_c, step = 1000)

#############################################################################
# plotting ratification curves
#############################################################################

plot_rarefaction <- function(df, title, filename) {
  
  p <- ggplot(df,
              aes(x = Depth,
                  y = Observed,
                  color = Treatment,
                  group = Sample)) +
    
    geom_line(alpha = 0.7) +
    
    theme_bw() +
    
    labs(
      title = title,
      x = "Sequencing depth",
      y = "Observed taxa"
    )
  
  ggsave(filename, p, width = 10, height = 7, dpi = 600)
  
  return(p)
}

#############################################################################
# generate and save plots
#############################################################################

p_m <- plot_rarefaction(
  rar_M,
  "Rarefaction curves (Males)",
  "rarefaction_males.png"
)

p_f <- plot_rarefaction(
  rar_F,
  "Rarefaction curves (Females)",
  "rarefaction_females.png"
)

#############################################################################
# rarefucation plots mean by group 
#############################################################################
plot_rarefaction_mean <- function(df, title, filename) {
  
  summary_df <- df %>%
    group_by(Treatment, Depth) %>%
    summarise(MeanObserved = mean(Observed), .groups = "drop")
  
  p <- ggplot() +
    
    geom_line(
      data = df,
      aes(Depth, Observed, group = Sample, color = Treatment),
      alpha = 0.25
    ) +
    
    geom_line(
      data = summary_df,
      aes(Depth, MeanObserved, color = Treatment),
      linewidth = 1.2
    ) +
    
    theme_bw() +
    
    labs(
      title = title,
      x = "Sequencing depth",
      y = "Observed taxa"
    )
  
  ggsave(filename, p, width = 10, height = 7, dpi = 600)
  
  return(p)
}

p_M <- plot_rarefaction_mean(
  rar_M,
  title = "Male rarefaction (mean + samples)",
  filename = "rarefaction_mean_males.png"
)

p_F <- plot_rarefaction_mean(
  rar_F,
  title = "Female rarefaction (mean + samples)",
  filename = "rarefaction_mean_females.png"
)


#############################################################################
#############################################################################
#Now rarefy based on your selected cutoff
#############################################################################

library(phyloseq)

# Look at sample read counts
sample_sums(ps_M_c)  

#convert to integers
otu_mat_int <- round(otu_table(ps_M_c))
ps_M_c_int <- ps_M_c
otu_table(ps_M_c_int) <- otu_mat_int

ps_M_rar <- rarefy_even_depth(
  ps_M_c_int,                 #phyloseq object
  sample.size = 9742, # rarefaction depth
  rngseed = 123,       # for reproducibility
  replace = FALSE,     # no replacement, standard subsampling
  trimOTUs = FALSE      # remove OTUs not present after rarefaction
)

#check if rarefication worked 
colSums(otu_table(ps_M_rar))


# Look at sample read counts
sample_sums(ps_F_c)  # ps is your phyloseq object

#convert to integers
otu_mat_int <- round(otu_table(ps_F_c))
ps_F_c_int <- ps_F_c
otu_table(ps_F_c_int) <- otu_mat_int

ps_F_rar <- rarefy_even_depth(
  ps_F_c_int,                 #phyloseq object
  sample.size = 14572, # rarefaction depth
  rngseed = 123,       # for reproducibility
  replace = FALSE,     # no replacement, standard subsampling
  trimOTUs = FALSE      # remove OTUs not present after rarefaction
)

#check if rarefication worked 
colSums(otu_table(ps_F_rar))


########################################################################################################################
########################################################################################################################
########################################################################################################################
#PLOTS
########################################################################################################################
########################################################################################################################
#MAKING ABUNDANCE PLOTS
########################################################################################################################
#PHYLUM LEVEL ABUNDANCE
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# -----------------------------
# Step 0: Collapse to Phylum level
# -----------------------------
ps_F_phylum <- tax_glom(ps_F_rar, taxrank = "Phylum")
ps_M_phylum <- tax_glom(ps_M_rar, taxrank = "Phylum")

# Convert to long-format data frames
df_F <- psmelt(ps_F_phylum)
df_M <- psmelt(ps_M_phylum)

# Make sure your sample metadata includes Treatment info
# If not, pull it from phyloseq
if(!"Treatment" %in% colnames(df_F)){
  df_F$Treatment <- sample_data(ps_F_phylum)$Treatment[df_F$Sample]
}
if(!"Treatment" %in% colnames(df_M)){
  df_M$Treatment <- sample_data(ps_M_phylum)$Treatment[df_M$Sample]
}

# Add Sex column manually
df_F$Sex <- "F"
df_M$Sex <- "M"

# -----------------------------
# Step 1: Function to plot Phylum abundance
# -----------------------------
plot_phylum_with_other <- function(df, sex_label, output_file, top_n = 20, relative = TRUE) {
  
  # Filter by sex
  df_sex <- df %>% filter(Sex == sex_label)
  
  # If plotting relative abundance, convert counts to fractions per sample
  if (relative) {
    df_sex <- df_sex %>%
      group_by(Sample) %>%
      mutate(Abundance = Abundance / sum(Abundance)) %>%
      ungroup()
  }
  
  # Aggregate at Phylum level per sample
  df_sex_agg <- df_sex %>%
    group_by(Sample, Treatment, Phylum) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # Identify top N phyla across all samples
  top_phylum <- df_sex_agg %>%
    group_by(Phylum) %>%
    summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
    arrange(desc(TotalAbundance)) %>%
    slice_head(n = top_n) %>%
    pull(Phylum)
  
  # Collapse non-top phyla into "Other"
  df_sex_top <- df_sex_agg %>%
    mutate(Phylum = ifelse(Phylum %in% top_phylum, Phylum, "Other")) %>%
    group_by(Sample, Treatment, Phylum) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # Force "Other" to be last factor level
  phy_levels <- c(sort(setdiff(unique(df_sex_top$Phylum), "Other")), "Other")
  df_sex_top$Phylum <- factor(df_sex_top$Phylum, levels = phy_levels)
  df_sex_top$Treatment <- factor(
    as.character(df_sex_top$Treatment),
    levels = c("Control", "HFiD", "DSS", "DSS+HFiD")
  )
  # Generate colors
  n_phy <- length(phy_levels) - 1
  phy_colors <- RColorBrewer::brewer.pal(min(n_phy, 12), "Set3")
  if (n_phy > 12) phy_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_phy)
  phy_colors <- c(phy_colors, "gray60")
  names(phy_colors) <- phy_levels
  
  # Plot
  p <- ggplot(df_sex_top,
              aes(x = Sample, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", color = "black", size = 0.2) +
    facet_wrap(~ Treatment, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = phy_colors) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.direction = "vertical"
    ) +
    ylab(ifelse(relative, "Relative Abundance", "Read Counts")) +
    ggtitle(paste("Phylum-level Abundance -", sex_label)) +
    guides(fill = guide_legend(ncol = 1))   # force one column
  
  # Save plot
  ggsave(output_file, plot = p, width = 10, height = 4, dpi = 600)
  
  return(p)
}

# -----------------------------
# Step 2: Plot females
# -----------------------------
p_female <- plot_phylum_with_other(
  df = df_F,
  sex_label = "F",
  output_file = "phylum_females.png",
  top_n = 20,
  relative = TRUE
)

# -----------------------------
# Step 3: Plot males
# -----------------------------
p_male <- plot_phylum_with_other(
  df = df_M,
  sex_label = "M",
  output_file = "phylum_males.png",
  top_n = 20,
  relative = TRUE
)

# Display plots
print(p_female)
print(p_male)

##################################################################################################
#FAMILY level plots

library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# -----------------------------
# Step 0: Collapse to Family level
# -----------------------------
ps_F_family <- tax_glom(ps_F_rar, taxrank = "Family")
ps_M_family <- tax_glom(ps_M_rar, taxrank = "Family")

# Convert to long-format data frames
df_F <- psmelt(ps_F_family)
df_M <- psmelt(ps_M_family)

# Add Treatment if missing
if(!"Treatment" %in% colnames(df_F)){
  df_F$Treatment <- sample_data(ps_F_family)$Treatment[df_F$Sample]
}
if(!"Treatment" %in% colnames(df_M)){
  df_M$Treatment <- sample_data(ps_M_family)$Treatment[df_M$Sample]
}

# Add Sex column
df_F$Sex <- "F"
df_M$Sex <- "M"

# -----------------------------
# Step 1: Calculate top N separately for each sex
# -----------------------------
top_n <- 20

ignore_families <- c("CAG-508", "UBA1381")

# Females
top_family_F <- df_F %>%
  filter(!Family %in% ignore_families) %>%   # exclude ignored
  group_by(Family) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 20) %>%
  pull(Family)

# Males
top_family_M <- df_M %>%
  filter(!Family %in% ignore_families) %>%   # exclude ignored
  group_by(Family) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 15) %>%
  pull(Family)

# Union of top families for consistent coloring
all_top <- union(top_family_F, top_family_M)

# Assign colors to union first
n_all <- length(all_top)
colors_base <- brewer.pal(min(n_all, 12), "Set3")
if (n_all > 12) colors_base <- colorRampPalette(brewer.pal(12, "Set3"))(n_all)
names(colors_base) <- all_top

# -----------------------------
# Step 2: Function to plot Family abundance with hybrid coloring
# -----------------------------
plot_family_with_other_hybrid <- function(df, sex_label, top_family, output_file, relative = TRUE) {
  
  df_sex <- df %>% filter(Sex == sex_label)
  
  # Convert to relative abundance if desired
  if (relative) {
    df_sex <- df_sex %>%
      group_by(Sample) %>%
      mutate(Abundance = Abundance / sum(Abundance)) %>%
      ungroup()
  }
  
  # Aggregate at Family level per sample
  df_sex_agg <- df_sex %>%
    group_by(Sample, Treatment, Family) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # Collapse non-top families into "Other"
  # Inside your plotting function, after computing top N
  df_sex_top <- df_sex_agg %>%
    mutate(
      Family = ifelse(
        Family %in% top_family, Family,   # keep top N
        "Other"                           # everything else goes to Other
      ),
      # Make sure ignored families are always "Other"
      Family = ifelse(Family %in% ignore_families, "Other", Family)
    ) %>%
    group_by(Sample, Treatment, Family) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # Get factor levels: top families first, then Other
  fam_levels <- c(sort(setdiff(top_family, "Other")), "Other")
  df_sex_top$Family <- factor(df_sex_top$Family, levels = fam_levels)
  
  df_sex_top$Treatment <- factor(
    as.character(df_sex_top$Treatment),
    levels = c("Control", "HFiD", "DSS", "DSS+HFiD")
  )
  
  # Generate colors
  # Use preassigned colors for families in all_top, assign new colors to unique families
  fam_colors <- sapply(fam_levels, function(fam) {
    if(fam %in% names(colors_base)) {
      colors_base[fam]
    } else if(fam == "Other") {
      "gray60"
    } else {
      # assign a random color if somehow new
      sample(colors(), 1)
    }
  })
  names(fam_colors) <- fam_levels
  
  # Plot
  p <- ggplot(df_sex_top, aes(x = Sample, y = Abundance, fill = Family)) +
    geom_bar(stat = "identity", color = "black", size = 0.2) +
    facet_wrap(~ Treatment, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = fam_colors) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.direction = "vertical"
    ) +
    ylab(ifelse(relative, "Relative Abundance", "Read Counts")) +
    ggtitle(paste("Family-level Abundance -", sex_label)) +
    guides(fill = guide_legend(ncol = 1))
  
  ggsave(output_file, plot = p, width = 10, height = 6, dpi = 600)
  
  return(p)
}

# -----------------------------
# Step 3: Plot females
# -----------------------------
p_female <- plot_family_with_other_hybrid(df_F, "F", top_family_F, "family_females.png")

# -----------------------------
# Step 4: Plot males
# -----------------------------
p_male <- plot_family_with_other_hybrid(df_M, "M", top_family_M, "family_males.png")

# Display
print(p_female)
print(p_male)


##################################################################################################
#SPECIES PLOT
##################################################################################################
# ===========================
# Species-level plots (top 30)
# ===========================
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# -----------------------------
# Step 0: Collapse to Species level
# -----------------------------
ps_F_species <- tax_glom(ps_F_rar, taxrank = "Species")
ps_M_species <- tax_glom(ps_M_rar, taxrank = "Species")

# Convert to long-format data frames
df_F <- psmelt(ps_F_species)
df_M <- psmelt(ps_M_species)

# Add Treatment if missing
if(!"Treatment" %in% colnames(df_F)){
  df_F$Treatment <- sample_data(ps_F_species)$Treatment[df_F$Sample]
}

if(!"Treatment" %in% colnames(df_M)){
  df_M$Treatment <- sample_data(ps_M_species)$Treatment[df_M$Sample]
}

# Add Sex column
df_F$Sex <- "F"
df_M$Sex <- "M"

# -----------------------------
# Step 1: Calculate top N separately for each sex
# -----------------------------
top_n <- 30

ignore_species <- c(
  "CAG-508 sp.",
  "UBA-1381 sp."
)

# Females
top_species_F <- df_F %>%
  filter(!Species %in% ignore_species) %>%
  group_by(Species) %>%
  summarise(
    TotalAbundance = sum(Abundance),
    .groups = "drop"
  ) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = top_n) %>%
  pull(Species)

# Males
top_species_M <- df_M %>%
  filter(!Species %in% ignore_species) %>%
  group_by(Species) %>%
  summarise(
    TotalAbundance = sum(Abundance),
    .groups = "drop"
  ) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = top_n) %>%
  pull(Species)

# -----------------------------
# Union of top species
# -----------------------------
all_top <- union(
  top_species_F,
  top_species_M
)

# -----------------------------
# Original Set3 palette
# -----------------------------
n_all <- length(all_top)

colors_base <- brewer.pal(
  min(n_all, 12),
  "Set3"
)

if (n_all > 12) {
  colors_base <- colorRampPalette(
    brewer.pal(12, "Set3")
  )(n_all)
}

# -----------------------------
# RANDOMIZE COLOR ORDER
# -----------------------------
set.seed(42)

colors_base <- sample(colors_base)

# Assign names AFTER randomization
names(colors_base) <- all_top

# -----------------------------
# Step 2: Function to plot
# -----------------------------
plot_species_with_other_hybrid <- function(
    df,
    sex_label,
    top_species,
    output_file,
    relative = TRUE
) {
  
  df_sex <- df %>%
    filter(Sex == sex_label)
  
  # Relative abundance
  if (relative) {
    
    df_sex <- df_sex %>%
      group_by(Sample) %>%
      mutate(
        Abundance = Abundance / sum(Abundance)
      ) %>%
      ungroup()
  }
  
  # Aggregate
  df_sex_agg <- df_sex %>%
    group_by(
      Sample,
      Treatment,
      Species
    ) %>%
    summarise(
      Abundance = sum(Abundance),
      .groups = "drop"
    )
  
  # Collapse into Other
  df_sex_top <- df_sex_agg %>%
    mutate(
      Species = ifelse(
        Species %in% top_species,
        Species,
        "Other"
      ),
      
      Species = ifelse(
        Species %in% ignore_species,
        "Other",
        Species
      )
    ) %>%
    group_by(
      Sample,
      Treatment,
      Species
    ) %>%
    summarise(
      Abundance = sum(Abundance),
      .groups = "drop"
    )
  
  # Treatment ordering
  df_sex_top$Treatment <- factor(
    as.character(df_sex_top$Treatment),
    levels = c(
      "Control",
      "HFiD",
      "DSS",
      "DSS+HFiD"
    )
  )
  
  # Species ordering
  sp_levels <- c(
    sort(setdiff(top_species, "Other")),
    "Other"
  )
  
  df_sex_top$Species <- factor(
    df_sex_top$Species,
    levels = sp_levels
  )
  
  # -----------------------------
  # Color assignment
  # -----------------------------
  sp_colors <- sapply(
    sp_levels,
    function(sp) {
      
      if (sp == "Other") {
        
        "grey60"
        
      } else if (sp %in% names(colors_base)) {
        
        colors_base[sp]
        
      } else {
        
        "grey85"
      }
    }
  )
  
  names(sp_colors) <- sp_levels
  
  # -----------------------------
  # Plot
  # -----------------------------
  p <- ggplot(
    df_sex_top,
    aes(
      x = Sample,
      y = Abundance,
      fill = Species
    )
  ) +
    
    geom_bar(
      stat = "identity",
      color = "black",
      size = 0.2
    ) +
    
    facet_wrap(
      ~ Treatment,
      nrow = 1,
      scales = "free_x"
    ) +
    
    scale_fill_manual(values = sp_colors) +
    
    theme_bw() +
    
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      
      strip.background = element_rect(fill = "lightgray"),
      
      strip.text = element_text(
        face = "bold",
        size = 12
      ),
      
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 15   # <-- increase title size here
      ),
      
      legend.position = "right",
      
      legend.text = element_text(
        size = 13,
        face = "italic"
      ),
      
      legend.title = element_text(
        size = 15,     # <-- "Species" label size
        face = "bold"
      )
    ) +
    
    ylab(
      ifelse(
        relative,
        "Relative Abundance",
        "Read Counts"
      )
    ) +
    
    ggtitle(
      paste(
        "Species-level Abundance -",
        sex_label
      )
    ) +
    
    guides(
      fill = guide_legend(ncol = 1)
    )
  
  # Save
  ggsave(
    output_file,
    plot = p,
    width = 15,
    height = 9,
    dpi = 600
  )
  
  return(p)
}

# -----------------------------
# Step 3: Plot females
# -----------------------------
p_female <- plot_species_with_other_hybrid(
  df_F,
  "F",
  top_species_F,
  "species_females.png",
  relative = TRUE
)

# -----------------------------
# Step 4: Plot males
# -----------------------------
p_male <- plot_species_with_other_hybrid(
  df_M,
  "M",
  top_species_M,
  "species_males.png",
  relative = TRUE
)

# Display
print(p_female)
print(p_male)

##################################################################################################
#### SAVING RELATIVE ABUNDANCE OUTPUTS
##################################################################################################
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)

get_rel_abundance_table <- function(ps, rank){
  
  # 1. Relative abundance
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  
  # 2. Extract tables
  otu <- as.data.frame(otu_table(ps_rel))
  tax <- as.data.frame(tax_table(ps_rel))
  
  # ensure taxa are rows
  if (!taxa_are_rows(ps_rel)) {
    otu <- t(otu)
  }
  
  otu <- as.data.frame(otu)
  otu$TaxonID <- rownames(otu)
  
  tax$TaxonID <- rownames(tax)
  
  # 3. join taxonomy
  df <- otu %>%
    left_join(tax, by = "TaxonID")
  
  # 4. choose rank column
  rank_col <- switch(rank,
                     "Phylum" = "Phylum",
                     "Family" = "Family",
                     "Species" = "Species",
                     stop("Invalid rank"))
  
  # 5. clean taxonomy labels
  df[[rank_col]] <- ifelse(is.na(df[[rank_col]]) | df[[rank_col]] == "",
                           "Unclassified",
                           df[[rank_col]])
  
  # 6. collapse by taxonomic rank
  df_collapsed <- df %>%
    group_by(Taxon = .data[[rank_col]]) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop")
  
  # 7. return with taxa as first column (already correct)
  df_collapsed
}

# Female
phylum_F   <- get_rel_abundance_table(ps_F_rar, "Phylum")
family_F   <- get_rel_abundance_table(ps_F_rar, "Family")
species_F  <- get_rel_abundance_table(ps_F_rar, "Species")

# Male
phylum_M   <- get_rel_abundance_table(ps_M_rar, "Phylum")
family_M   <- get_rel_abundance_table(ps_M_rar, "Family")
species_M  <- get_rel_abundance_table(ps_M_rar, "Species")

library(openxlsx)

write.xlsx(phylum_F, "phylum_F_rel_abundance.xlsx")
write.xlsx(family_F, "family_F_rel_abundance.xlsx")
write.xlsx(species_F, "species_F_rel_abundance.xlsx")

write.xlsx(phylum_M, "phylum_M_rel_abundance.xlsx")
write.xlsx(family_M, "family_M_rel_abundance.xlsx")
write.xlsx(species_M, "species_M_rel_abundance.xlsx")

##############################################################################################
#ALPHA DIVERSITY
##############################################################################################

library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)

# --- Function to calculate alpha diversity ---
calculate_alpha <- function(ps) {
  otu_mat <- as(otu_table(ps), "matrix")
  if(taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  
  shannon <- vegan::diversity(otu_mat, index = "shannon")
  simpson <- vegan::diversity(otu_mat, index = "simpson")
  observed <- rowSums(otu_mat > 0)
  
  df <- data.frame(
    SampleID = rownames(otu_mat),
    Observed = observed,
    Shannon = shannon,
    Simpson = simpson
  )
  
  # Add metadata
  meta <- data.frame(sample_data(ps))
  df$Treatment <- meta$Treatment[match(df$SampleID, rownames(meta))]
  df$Sex <- meta$Sex[match(df$SampleID, rownames(meta))]
  
  # Keep valid samples
  df <- df %>% filter(!is.na(Treatment), !is.na(Sex))
  
  df$Treatment <- factor(df$Treatment, levels = c("Control", "HFiD", "DSS", "DSS+HFiD"))
  df$Sex <- factor(df$Sex, levels = c("F","M"))
  
  df_long <- df %>%
    pivot_longer(cols = c("Observed","Shannon","Simpson"),
                 names_to = "Metric",
                 values_to = "Value") %>%
    mutate(Metric = factor(Metric, levels = c("Observed","Shannon","Simpson")))
  
  return(df_long)
}

# --- Plotting function (Cell style) ---
plot_alpha_box_signif_cell <- function(alpha_long, sex_label, output_excel = NULL) {
  
  df <- alpha_long %>% filter(Sex == sex_label)
  
  comparisons <- list(
    c("Control", "DSS"),
    c("Control", "DSS+HFiD"),
    c("HFiD", "DSS+HFiD"),
    c("DSS", "DSS+HFiD")
  )
  
  # Calculate significance per metric
  sig_df <- do.call(rbind, lapply(levels(df$Metric), function(metric) {
    df_metric <- df %>% filter(Metric == metric)
    res <- lapply(comparisons, function(comp) {
      val1 <- df_metric$Value[df_metric$Treatment == comp[1]]
      val2 <- df_metric$Value[df_metric$Treatment == comp[2]]
      p <- if(length(val1) > 0 & length(val2) > 0) wilcox.test(val1, val2, exact = FALSE)$p.value else NA
      sig <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", NA)))
      data.frame(Metric = metric, Comparison = paste(comp, collapse = " vs "), p_value = p, label = sig)
    })
    do.call(rbind, res)
  }))
  
  # Keep only significant comparisons
  sig_df_plot <- sig_df %>% filter(!is.na(label))
  if(nrow(sig_df_plot) > 0){
    sig_df_plot <- sig_df_plot %>%
      group_by(Metric) %>%
      mutate(order = row_number(),
             y = max(df$Value[df$Metric == Metric], na.rm = TRUE) * 1 +
               order * 0.03 * max(df$Value[df$Metric == Metric], na.rm = TRUE)) %>%
      ungroup()
  }
  
  # Cell-style muted palette
  # Slightly muted / dimmed Cell-style palette
  cell_palette <- c(
    "Control"   = "#2C7FB8",
    "HFiD"      = "#D9A441",
    "DSS"       = "#B96A9B",
    "DSS+HFiD"  = "#2FA187"
  )
  
  p <- ggplot(df, aes(x = Treatment, y = Value, fill = Treatment)) +
    geom_boxplot(width=0.6, outlier.shape = NA, color="black", alpha=0.85, linewidth=0.4) +
    geom_jitter(width=0.15, size=1.5, alpha=0.6, color="black") +
    facet_wrap(~Metric, scales="free_y") +
    scale_fill_manual(values = cell_palette) +
    theme_classic(base_size = 12) +
    labs(
      x = NULL,
      y = "Alpha Diversity",
      title = paste("Alpha Diversity -", ifelse(sex_label=="F","Females","Males"))
    ) +
    theme(
      legend.position = "none",
      strip.text = element_text(size=14, face="bold"),
      axis.title.y = element_text(size=12, face="bold"),
      axis.text = element_text(size=10, color="black"),
      plot.title = element_text(size=14, face="bold", hjust=0.5)
    ) +
    geom_segment(data = sig_df_plot, 
                 aes(x = match(sub(" vs .*", "", Comparison), levels(df$Treatment)),
                     xend = match(sub(".* vs ", "", Comparison), levels(df$Treatment)),
                     y = y, yend = y),
                 inherit.aes = FALSE, linewidth=0.6) +
    geom_text(data = sig_df_plot,
              aes(x = (match(sub(" vs .*", "", Comparison), levels(df$Treatment)) +
                         match(sub(".* vs ", "", Comparison), levels(df$Treatment)))/2,
                  y = y,
                  label = label),
              inherit.aes = FALSE,
              vjust=-0.5,
              size=4)
  
  return(list(plot = p, pvalues = sig_df))
}

# --- Example usage ---
alpha_F <- calculate_alpha(ps_F_rar)
alpha_M <- calculate_alpha(ps_M_rar)

res_F <- plot_alpha_box_signif_cell(alpha_F, "F")
res_M <- plot_alpha_box_signif_cell(alpha_M, "M")

# --- Save high-quality plots ---
ggsave("alpha_females_cell.png", res_F$plot, width=10, height=6, dpi=600)
ggsave("alpha_males_cell.png", res_M$plot, width=10, height=6, dpi=600)

save_alpha_outputs <- function(alpha_long, sex_label, prefix){
  
  # wide table
  wide <- alpha_long %>%
    pivot_wider(names_from = Metric, values_from = Value)
  
  write.xlsx(wide, paste0(prefix, "_alpha_values.xlsx"), rowNames = FALSE)
  
  # stats (recompute same way as in plotting function)
  comparisons <- list(
    c("Control", "DSS"),
    c("Control", "DSS+HFiD"),
    c("HFiD", "DSS+HFiD"),
    c("DSS", "DSS+HFiD")
  )
  
  df <- alpha_long %>% filter(Sex == sex_label)
  
  stats <- do.call(rbind, lapply(levels(df$Metric), function(metric) {
    
    df_metric <- df %>% filter(Metric == metric)
    
    do.call(rbind, lapply(comparisons, function(comp) {
      
      g1 <- df_metric$Value[df_metric$Treatment == comp[1]]
      g2 <- df_metric$Value[df_metric$Treatment == comp[2]]
      
      # skip incomplete comparisons
      if(length(g1) < 2 | length(g2) < 2){
        return(data.frame(
          Metric = metric,
          Group1 = comp[1],
          Group2 = comp[2],
          Test = "Wilcoxon rank-sum",
          n1 = length(g1),
          n2 = length(g2),
          median1 = median(g1, na.rm = TRUE),
          median2 = median(g2, na.rm = TRUE),
          W = NA,
          p_value = NA
        ))
      }
      
      test <- wilcox.test(g1, g2, exact = FALSE)
      
      data.frame(
        Metric = metric,
        Group1 = comp[1],
        Group2 = comp[2],
        Test = "Wilcoxon rank-sum",
        n1 = length(g1),
        n2 = length(g2),
        median1 = median(g1, na.rm = TRUE),
        median2 = median(g2, na.rm = TRUE),
        W = unname(test$statistic),
        p_value = test$p.value
      )
    }))
  }))
  
  write.xlsx(stats, paste0(prefix, "_alpha_stats.xlsx"), rowNames = FALSE)
}

save_alpha_outputs(alpha_F, "F", "females")
save_alpha_outputs(alpha_M, "M", "males")

#######################################################################################################
# beta-diversity FEMALES
#######################################################################################################

# --- Metadata ---
meta <- as(sample_data(ps_F_rar), "data.frame")
meta$SampleID <- rownames(meta)

# --- Bray–Curtis distance ---
bray_dist <- phyloseq::distance(ps_F_rar, method = "bray")

# --- PCoA ---
ord_bray <- ape::pcoa(bray_dist)

# --- Make plotting data frame ---
ord_df <- as.data.frame(ord_bray$vectors[, 1:2])
colnames(ord_df) <- c("Axis.1", "Axis.2")
ord_df$SampleID <- rownames(ord_df)
ord_df <- merge(ord_df, meta, by = "SampleID")

# --- PERMANOVA ---
meta$Treatment <- factor(meta$Treatment, levels = c("Control", "HFiD", "DSS", "DSS+HFiD"))
adonis_res <- adonis2(bray_dist ~ Treatment, data = meta, permutations = 999)
pval <- adonis_res$`Pr(>F)`[1]

# --- Plot ---
eig1 <- round(ord_bray$values$Relative_eig[1] * 100, 1)
eig2 <- round(ord_bray$values$Relative_eig[2] * 100, 1)
palette <- RColorBrewer::brewer.pal(4, "Set2")

p_beta_female <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Treatment, fill = Treatment)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = Treatment, fill = Treatment), geom = "polygon",
               alpha = 0.2, color = NA, level = 0.9, type = "norm") +
  scale_color_manual(values = setNames(palette, levels(ord_df$Treatment))) +
  scale_fill_manual(values = setNames(palette, levels(ord_df$Treatment))) +
  theme_bw() +
  labs(
    title = "Beta Diversity (PCoA - Bray–Curtis) - Females",
    subtitle = paste("PERMANOVA p =", signif(pval, 3)),
    x = paste0("PCoA 1 (", eig1, "%)"),
    y = paste0("PCoA 2 (", eig2, "%)")
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

# --- Display / Save ---
p_beta_female
ggsave("beta_females.png", plot = p_beta_female, width = 8, height = 6, dpi = 600)

library(phyloseq)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ape)
library(ggpubr)

# --- Metadata ---
meta <- as(sample_data(ps_F_rar), "data.frame")
meta$SampleID <- rownames(meta)
meta$Treatment <- factor(meta$Treatment, levels = c("Control", "HFiD", "DSS", "DSS+HFiD"))

# --- Bray–Curtis distance and PCoA ---
bray_dist <- phyloseq::distance(ps_F_rar, method = "bray")
ord_bray <- ape::pcoa(bray_dist)

# --- Prepare plotting data frame ---
ord_df <- as.data.frame(ord_bray$vectors[,1:2])
colnames(ord_df) <- c("Axis.1","Axis.2")
ord_df$SampleID <- rownames(ord_df)
ord_df <- merge(ord_df, meta, by="SampleID")

# --- PERMANOVA ---
adonis_res <- adonis2(bray_dist ~ Treatment, data = meta, permutations = 999)
pval <- adonis_res$`Pr(>F)`[1]
pval_text <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", formatC(pval, format="f", digits=3)))

# --- Variance explained ---
eig1 <- round(ord_bray$values$Relative_eig[1] * 100, 1)
eig2 <- round(ord_bray$values$Relative_eig[2] * 100, 1)

# --- Cell-style muted palette ---
cell_palette <- c(
  "Control"   = "#2C7FB8",
  "HFiD"      = "#D9A441",
  "DSS"       = "#B96A9B",
  "DSS+HFiD"  = "#2FA187"
)

# --- Plot ---
p_beta_female <- ggplot(ord_df, aes(x=Axis.1, y=Axis.2, color=Treatment, fill=Treatment)) +
  geom_point(size=4, alpha=0.8, stroke=0.5, shape=21, color="black") +
  stat_ellipse(aes(group=Treatment, fill=Treatment), geom="polygon", alpha=0.2, color=NA, level=0.9, type="norm") +
  scale_color_manual(values=cell_palette) +
  scale_fill_manual(values=cell_palette) +
  labs(
    title = "Beta Diversity (PCoA - Bray–Curtis) - Females",
    subtitle = pval_text,
    x = paste0("PCoA 1 (", eig1, "%)"),
    y = paste0("PCoA 2 (", eig2, "%)")
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title = element_text(face="bold", size=14, hjust=0.5),
    plot.subtitle = element_text(face="italic", size=11, hjust=0.5),
    axis.title = element_text(face="bold", size=12),
    axis.text = element_text(size=10, color="black"),
    legend.position="right",
    legend.title = element_text(face="bold", size=11),
    legend.text = element_text(size=10)
  )

# --- Save ---
ggsave("beta_females_cell.png", p_beta_female, width=8, height=6, dpi=600)
##########################################################################################
#BETA MALES
##########################################################################################
library(phyloseq)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ape)
library(ggpubr)

# --- Metadata ---
meta <- as(sample_data(ps_M_rar), "data.frame")
meta$SampleID <- rownames(meta)
meta$Treatment <- factor(meta$Treatment, levels = c("Control", "HFiD", "DSS", "DSS+HFiD"))

# --- Bray–Curtis distance and PCoA ---
bray_dist <- phyloseq::distance(ps_M_rar, method = "bray")
ord_bray <- ape::pcoa(bray_dist)

# --- Prepare plotting data frame ---
ord_df <- as.data.frame(ord_bray$vectors[,1:2])
colnames(ord_df) <- c("Axis.1","Axis.2")
ord_df$SampleID <- rownames(ord_df)
ord_df <- merge(ord_df, meta, by="SampleID")

# --- PERMANOVA ---
adonis_res <- adonis2(bray_dist ~ Treatment, data = meta, permutations = 999)
pval <- adonis_res$`Pr(>F)`[1]
pval_text <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", formatC(pval, format="f", digits=3)))

# --- Variance explained ---
eig1 <- round(ord_bray$values$Relative_eig[1] * 100, 1)
eig2 <- round(ord_bray$values$Relative_eig[2] * 100, 1)

# --- Cell-style muted palette ---
cell_palette <- c(
  "Control"   = "#2C7FB8",
  "HFiD"      = "#D9A441",
  "DSS"       = "#B96A9B",
  "DSS+HFiD"  = "#2FA187"
)

# --- Plot ---
p_beta_female <- ggplot(ord_df, aes(x=Axis.1, y=Axis.2, color=Treatment, fill=Treatment)) +
  geom_point(size=4, alpha=0.8, stroke=0.5, shape=21, color="black") +
  stat_ellipse(aes(group=Treatment, fill=Treatment), geom="polygon", alpha=0.2, color=NA, level=0.9, type="norm") +
  scale_color_manual(values=cell_palette) +
  scale_fill_manual(values=cell_palette) +
  labs(
    title = "Beta Diversity (PCoA - Bray–Curtis) - Males",
    subtitle = pval_text,
    x = paste0("PCoA 1 (", eig1, "%)"),
    y = paste0("PCoA 2 (", eig2, "%)")
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title = element_text(face="bold", size=14, hjust=0.5),
    plot.subtitle = element_text(face="italic", size=11, hjust=0.5),
    axis.title = element_text(face="bold", size=12),
    axis.text = element_text(size=10, color="black"),
    legend.position="right",
    legend.title = element_text(face="bold", size=11),
    legend.text = element_text(size=10)
  )

# --- Save ---
ggsave("beta_males_cell.png", p_beta_female, width=8, height=6, dpi=600)

#saving data

export_beta_excel <- function(ord_df, meta, bray_dist, sex_label, filename){
  
  library(openxlsx)
  library(dplyr)
  
  # -----------------------------
  # 1. PERMANOVA (global)
  # -----------------------------
  adonis_res <- adonis2(bray_dist ~ Treatment, data = meta, permutations = 999)
  
  adonis_df <- as.data.frame(adonis_res)
  adonis_df$Term <- rownames(adonis_df)
  adonis_df <- adonis_df %>%
    dplyr::select(Term, Df, `SumOfSqs`, R2, F, `Pr(>F)`)
  
  adonis_df$Sex <- sex_label
  adonis_df$Test <- "PERMANOVA (Bray-Curtis)"
  
  # -----------------------------
  # 2. Pairwise PERMANOVA
  # -----------------------------
  pairwise_permanova <- function(dist, meta){
    
    groups <- unique(meta$Treatment)
    combs <- combn(groups, 2, simplify = FALSE)
    
    out <- lapply(combs, function(x){
      
      sub_meta <- meta %>% filter(Treatment %in% x)
      sub_dist <- as.matrix(dist)[sub_meta$SampleID, sub_meta$SampleID]
      
      res <- adonis2(sub_dist ~ Treatment, data = sub_meta, permutations = 999)
      
      data.frame(
        Group1 = x[1],
        Group2 = x[2],
        F = res$F[1],
        R2 = res$R2[1],
        p = res$`Pr(>F)`[1]
      )
    })
    
    do.call(rbind, out)
  }
  
  pairwise_df <- pairwise_permanova(bray_dist, meta)
  pairwise_df$p_adj_BH <- p.adjust(pairwise_df$p, method = "BH")
  
  pairwise_df$Sex <- sex_label
  
  # -----------------------------
  # 3. Centroids
  # -----------------------------
  centroids <- ord_df %>%
    group_by(Treatment) %>%
    summarise(
      Axis.1 = mean(Axis.1),
      Axis.2 = mean(Axis.2),
      n = n()
    )
  
  centroids$Sex <- sex_label
  
  # -----------------------------
  # 4. PCoA coordinates
  # -----------------------------
  pcoa_df <- ord_df
  pcoa_df$Sex <- sex_label
  
  # -----------------------------
  # 5. Write multi-sheet Excel
  # -----------------------------
  wb <- createWorkbook()
  
  addWorksheet(wb, "PCoA")
  writeData(wb, "PCoA", pcoa_df)
  
  addWorksheet(wb, "Metadata")
  writeData(wb, "Metadata", meta)
  
  addWorksheet(wb, "PERMANOVA")
  writeData(wb, "PERMANOVA", adonis_df)
  
  addWorksheet(wb, "Pairwise_PERMANOVA")
  writeData(wb, "Pairwise_PERMANOVA", pairwise_df)
  
  addWorksheet(wb, "Centroids")
  writeData(wb, "Centroids", centroids)
  
  saveWorkbook(wb, filename, overwrite = TRUE)
}

# FEMALES
meta_F <- as(sample_data(ps_F_rar), "data.frame")
meta_F$SampleID <- rownames(meta_F)

bray_dist_F <- phyloseq::distance(ps_F_rar, method = "bray")

ord_bray_F <- ape::pcoa(bray_dist_F)

ord_df_F <- as.data.frame(ord_bray_F$vectors[, 1:2])
colnames(ord_df_F) <- c("Axis.1", "Axis.2")
ord_df_F$SampleID <- rownames(ord_df_F)
ord_df_F <- merge(ord_df_F, meta_F, by = "SampleID")

export_beta_excel(
  ord_df = ord_df_F,
  meta = meta_F,
  bray_dist = bray_dist_F,
  sex_label = "F",
  filename = "beta_diversity_females.xlsx"
)

meta_M <- as(sample_data(ps_M_rar), "data.frame")
meta_M$SampleID <- rownames(meta_M)

bray_dist_M <- phyloseq::distance(ps_M_rar, method = "bray")

ord_bray_M <- ape::pcoa(bray_dist_M)

ord_df_M <- as.data.frame(ord_bray_M$vectors[, 1:2])
colnames(ord_df_M) <- c("Axis.1", "Axis.2")
ord_df_M$SampleID <- rownames(ord_df_M)
ord_df_M <- merge(ord_df_M, meta_M, by = "SampleID")


export_beta_excel(
  ord_df = ord_df_M,
  meta = meta_M,
  bray_dist = bray_dist_M,
  sex_label = "M",
  filename = "beta_diversity_males.xlsx"
)

#################################################################################################
# BUBBLE PLOT - males
#################################################################################################
# -----------------------------
# 1. Collapse to Species level
# -----------------------------
# Remove stray apostrophes
tax_table(ps_M_rar)[, "Species"] <-
  gsub("'", "", tax_table(ps_M_rar)[, "Species"])

# Standardize anything starting with "Streptococcus sp."
tax_table(ps_M_rar)[, "Species"] <-
  gsub("^Streptococcus sp\\..*", 
       "Streptococcus sp.", 
       tax_table(ps_M_rar)[, "Species"])

ps_species <- tax_glom(ps_M_rar, taxrank = "Species")
df_species <- psmelt(ps_species)

# Make Treatment a factor in the desired order
df_species$Treatment <- factor(df_species$Treatment,
                               levels = c("Control", "HFiD", "DSS", "DSS+HFiD"))
# -----------------------------
# 2. Compute mean abundance and log2FC vs Control
# -----------------------------
df_bubble <- df_species %>%
  group_by(Species, Treatment, Sex) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Species, Sex) %>%
  mutate(control_abundance = mean_abundance[Treatment == "Control"],
         log2FC = log2((mean_abundance + 1e-6) / (control_abundance + 1e-6))) %>%
  ungroup()

# -----------------------------
# 3. Function to compute p-values and FDR vs Control
# -----------------------------
compute_diff <- function(df, sex, group1, group2 = "Control") {
  df_sub <- df %>% filter(Sex == sex, Treatment %in% c(group1, group2))
  
  res <- df_sub %>%
    group_by(Species) %>%
    summarise(
      pval = tryCatch(
        wilcox.test(Abundance[Treatment == group1],
                    Abundance[Treatment == group2])$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      padj = p.adjust(pval, method = "fdr"),
      Sex = sex,
      Treatment = group1
    )
  
  return(res)
}

groups_to_compare <- c("HFiD", "DSS", "DSS+HFiD")
res_list <- lapply(groups_to_compare, function(g) {
  rbind(
    compute_diff(df_species, sex = "M", group1 = g),
    compute_diff(df_species, sex = "M", group1 = g)
  )
})
res_all <- bind_rows(res_list)


exclude_species <- c(
  "Massiliimalia timonensis", 
  "Faecousia intestinalis",
  "Hominilimicola fabiformis",
  "Merdimmobilis hominis",
  "Lientehia dongpingensis",
  "Suilimivivens aceti",
  "Eisenbergiella massiliensis",
  "Romboutsia sp. G12",
  "Romboutsia timonensis",
  "Streptococcus salivarius",
  "Massiliimalia massiliensis",
  "Pelethomonas sp.",
  "Laedolimicola sp",
  "Gallimonas sp.",
  "Suilimivivens sp.",
  "Laedolimicola sp.",
  "Claveliimonas sp.",
  "Pullilachnospira sp.",
  "Brotaphodocola sp.",
  "Clostridium lamae",
  "Caproiciproducens sp.",
  "Hungatella hathewayi",
  "Enterocloster alcoholdehydrogenati",
  "Lacrimispora saccharolytica",
  "Lacrimispora brassicae",
  "Clostridium sp. B905-1",
  "Streptococcus oralis",
  "Streptocuccus salivarius",
  "Streptococcus sp.",
  "Other",
  "Fimisoma sp."
)

# For males
sig_species_m <- df_species %>%
  filter(Sex == "M") %>%
  group_by(Species) %>%
  summarise(any_sig = any(res_all$padj[res_all$Species == Species & res_all$Sex == "M"] <= 0.05, na.rm = TRUE)) %>%
  filter(any_sig) %>%
  pull(Species)

sig_species_m <- setdiff(sig_species_m, exclude_species)

df_sample_plot <- df_species %>%
  filter(Sex == "M", Species %in% sig_species_m) %>%   # keep controls for calculation
  left_join(res_all %>% filter(Sex == "M") %>% select(Species, Treatment, padj),
            by = c("Species", "Treatment")) %>%
  group_by(Species) %>%
  mutate(
    control_abundance = mean(Abundance[Treatment == "Control"], na.rm = TRUE),
    log2FC = log2((Abundance + 1e-6) / (control_abundance + 1e-6))
  ) %>%
  ungroup() %>%
  # create significance size
  mutate(
    sig_size = case_when(
      is.na(padj) | padj > 0.05 ~ "ns",
      padj <= 0.05 & padj > 0.01 ~ "≤ 0.05",
      padj <= 0.01 ~ "≤ 0.01"
    ),
    sig_size = factor(sig_size, levels = c("ns", "≤ 0.05", "≤ 0.01"))
  ) %>%
  # remove controls for plotting
  filter(Treatment != "Control")

# Ensure treatment order FIRST
df_sample_plot$Treatment <- factor(
  df_sample_plot$Treatment,
  levels = c("HFiD", "DSS", "DSS+HFiD")
)

# Compute ordering only from species present
species_order <- df_sample_plot %>%
  filter(Treatment == "DSS+HFiD") %>%
  group_by(Species) %>%
  summarise(mean_log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_log2FC) %>%   # highest → lowest
  pull(Species)

# Drop old levels before applying new ones
df_sample_plot$Species <- as.character(df_sample_plot$Species)

df_sample_plot$Species <- factor(
  df_sample_plot$Species,
  levels = species_order
)

max_fc <- max(abs(df_sample_plot$log2FC), na.rm = TRUE)
fc_lim <- quantile(abs(df_sample_plot$log2FC), probs = 0.95, na.rm = TRUE)

gg_sig_samples <- ggplot(
  df_sample_plot,
  aes(x = Sample, y = Species,
      color = log2FC, size = sig_size)
) +
  geom_point(alpha = 0.95) +
  scale_color_gradientn(
    colors = c(
      "#053061", "#2166ac", "#4393c3", "#92c5de", "#f7f7f7",
      "#f4a582", "#d6604d", "#b2182b", "#67001f"
    ),
    values = scales::rescale(
      c(-fc_lim, -fc_lim*0.6, -fc_lim*0.3, -fc_lim*0.1, 0,
        fc_lim*0.1,  fc_lim*0.3,  fc_lim*0.6,  fc_lim)
    ),
    limits = c(-fc_lim, fc_lim),
    oob = scales::squish
  ) +
  scale_size_manual(values = c("ns" = 3, "≤ 0.05" = 5, "≤ 0.01" = 7)) +
  facet_grid(~Treatment, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(
    title = "Significantly altered species abundances - males",
    x = "Sample ID (within Treatment)",
    y = "Species",
    color = "log2FC vs Control",
    size = "Adjusted p-value"
  ) +
  theme(
    axis.text.y = element_text(face = "italic"),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
    panel.grid.major.x = element_blank()
  )

gg_sig_samples

ggsave(
  filename = "significant_species_males.png",  # output file name
  plot = gg_sig_samples,                        # the plot object
  width = 9,                                   # width in inches
  height = 12,                                   # height in inches
  dpi = 300                                     # resolution
)


#saving data
library(openxlsx)
library(dplyr)

# -----------------------------
# Create workbook
# -----------------------------
wb <- createWorkbook()

# =========================================================
# 1. Differential abundance statistics (Wilcoxon + FDR)
# =========================================================
addWorksheet(wb, "Wilcoxon_DE_Stats")

writeData(
  wb,
  "Wilcoxon_DE_Stats",
  res_all %>%
    arrange(Sex, Species, Treatment)
)

# =========================================================
# 2. Effect size table (log2 fold-change vs Control)
# =========================================================
addWorksheet(wb, "Log2FC_EffectSizes")

writeData(
  wb,
  "Log2FC_EffectSizes",
  df_bubble %>%
    arrange(Sex, Species, Treatment)
)

# =========================================================
# 3. Final plotting dataset (filtered + significance)
# =========================================================
addWorksheet(wb, "BubblePlot_Data")

writeData(
  wb,
  "BubblePlot_Data",
  df_sample_plot %>%
    arrange(Sex, Species, Treatment, Sample)
)

# =========================================================
# 4. Species-level summary statistics
# =========================================================
species_summary <- df_species %>%
  group_by(Sex, Treatment, Species) %>%
  summarise(
    n_samples = n_distinct(Sample),
    mean_abundance = mean(Abundance, na.rm = TRUE),
    sd_abundance = sd(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

addWorksheet(wb, "Species_Summary")

writeData(
  wb,
  "Species_Summary",
  species_summary
)

# =========================================================
# 5. Save workbook
# =========================================================
saveWorkbook(
  wb,
  file = "bubble_plot_stats_males.xlsx",
  overwrite = TRUE
)

#################################################################################################
# BUBBLE PLOT - females
#################################################################################################
# -----------------------------
# 1. Collapse to Species level
# -----------------------------
# Remove stray apostrophes
tax_table(ps_F_rar)[, "Species"] <-
  gsub("'", "", tax_table(ps_F_rar)[, "Species"])

# Standardize anything starting with "Streptococcus sp."
tax_table(ps_F_rar)[, "Species"] <-
  gsub("^Streptococcus sp\\..*", 
       "Streptococcus sp.", 
       tax_table(ps_F_rar)[, "Species"])

ps_species <- tax_glom(ps_F_rar, taxrank = "Species")
df_species <- psmelt(ps_species)

# Make Treatment a factor in the desired order
df_species$Treatment <- factor(df_species$Treatment,
                               levels = c("Control", "HFiD", "DSS", "DSS+HFiD"))
# -----------------------------
# 2. Compute mean abundance and log2FC vs Control
# -----------------------------
df_bubble <- df_species %>%
  group_by(Species, Treatment, Sex) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Species, Sex) %>%
  mutate(control_abundance = mean_abundance[Treatment == "Control"],
         log2FC = log2((mean_abundance + 1e-6) / (control_abundance + 1e-6))) %>%
  ungroup()

# -----------------------------
# 3. Function to compute p-values and FDR vs Control
# -----------------------------
compute_diff <- function(df, sex, group1, group2 = "Control") {
  df_sub <- df %>% filter(Sex == sex, Treatment %in% c(group1, group2))
  
  res <- df_sub %>%
    group_by(Species) %>%
    summarise(
      pval = tryCatch(
        wilcox.test(Abundance[Treatment == group1],
                    Abundance[Treatment == group2])$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      padj = p.adjust(pval, method = "fdr"),
      Sex = sex,
      Treatment = group1
    )
  
  return(res)
}

groups_to_compare <- c("HFiD", "DSS", "DSS+HFiD")
res_list <- lapply(groups_to_compare, function(g) {
  rbind(
    compute_diff(df_species, sex = "F", group1 = g),
    compute_diff(df_species, sex = "F", group1 = g)
  )
})
res_all <- bind_rows(res_list)


exclude_species <- c(
  "Massiliimalia timonensis", 
  "Faecousia intestinalis",
  "Hominilimicola fabiformis",
  "Merdimmobilis hominis",
  "Lientehia dongpingensis",
  "Suilimivivens aceti",
  "Eisenbergiella massiliensis",
  "Romboutsia sp. G12",
  "Romboutsia timonensis",
  "Streptococcus salivarius",
  "Massiliimalia massiliensis",
  "Pelethomonas sp.",
  "Laedolimicola sp",
  "Gallimonas sp.",
  "Suilimivivens sp.",
  "Laedolimicola sp.",
  "Claveliimonas sp.",
  "Pullilachnospira sp.",
  "Brotaphodocola sp.",
  "Clostridium lamae",
  "Caproiciproducens sp.",
  "Hungatella hathewayi",
  "Enterocloster alcoholdehydrogenati",
  "Lacrimispora saccharolytica",
  "Lacrimispora brassicae",
  "Clostridium sp. B905-1",
  "Streptococcus oralis",
  "Streptocuccus salivarius",
  "Streptococcus sp.",
  "Other",
  "Fimisoma sp."
)

# For females
sig_species_m <- df_species %>%
  filter(Sex == "F") %>%
  group_by(Species) %>%
  summarise(any_sig = any(res_all$padj[res_all$Species == Species & res_all$Sex == "F"] <= 0.05, na.rm = TRUE)) %>%
  filter(any_sig) %>%
  pull(Species)

sig_species_m <- setdiff(sig_species_m, exclude_species)

df_sample_plot <- df_species %>%
  filter(Sex == "F", Species %in% sig_species_m) %>%   # keep controls for calculation
  left_join(res_all %>% filter(Sex == "F") %>% select(Species, Treatment, padj),
            by = c("Species", "Treatment")) %>%
  group_by(Species) %>%
  mutate(
    control_abundance = mean(Abundance[Treatment == "Control"], na.rm = TRUE),
    log2FC = log2((Abundance + 1e-6) / (control_abundance + 1e-6))
  ) %>%
  ungroup() %>%
  # create significance size
  mutate(
    sig_size = case_when(
      is.na(padj) | padj > 0.05 ~ "ns",
      padj <= 0.05 & padj > 0.01 ~ "≤ 0.05",
      padj <= 0.01 ~ "≤ 0.01"
    ),
    sig_size = factor(sig_size, levels = c("ns", "≤ 0.05", "≤ 0.01"))
  ) %>%
  # remove controls for plotting
  filter(Treatment != "Control")

# Ensure treatment order FIRST
df_sample_plot$Treatment <- factor(
  df_sample_plot$Treatment,
  levels = c("HFiD", "DSS", "DSS+HFiD")
)

# Compute ordering only from species present
species_order <- df_sample_plot %>%
  filter(Treatment == "DSS+HFiD") %>%
  group_by(Species) %>%
  summarise(mean_log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_log2FC) %>%   # highest → lowest
  pull(Species)

# Drop old levels before applying new ones
df_sample_plot$Species <- as.character(df_sample_plot$Species)

df_sample_plot$Species <- factor(
  df_sample_plot$Species,
  levels = species_order
)

max_fc <- max(abs(df_sample_plot$log2FC), na.rm = TRUE)
fc_lim <- quantile(abs(df_sample_plot$log2FC), probs = 0.95, na.rm = TRUE)

gg_sig_samples <- ggplot(
  df_sample_plot,
  aes(x = Sample, y = Species,
      color = log2FC, size = sig_size)
) +
  geom_point(alpha = 0.95) +
  scale_color_gradientn(
    colors = c(
      "#053061", "#2166ac", "#4393c3", "#92c5de", "#f7f7f7",
      "#f4a582", "#d6604d", "#b2182b", "#67001f"
    ),
    values = scales::rescale(
      c(-fc_lim, -fc_lim*0.6, -fc_lim*0.3, -fc_lim*0.1, 0,
        fc_lim*0.1,  fc_lim*0.3,  fc_lim*0.6,  fc_lim)
    ),
    limits = c(-fc_lim, fc_lim),
    oob = scales::squish
  ) +
  scale_size_manual(values = c("ns" = 3, "≤ 0.05" = 5, "≤ 0.01" = 7)) +
  facet_grid(~Treatment, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(
    title = "Significantly altered species abundances - females",
    x = "Sample ID (within Treatment)",
    y = "Species",
    color = "log2FC vs Control",
    size = "Adjusted p-value"
  ) +
  theme(
    axis.text.y = element_text(face = "italic"),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
    panel.grid.major.x = element_blank()
  )

gg_sig_samples

ggsave(
  filename = "significant_species_females.png",  # output file name
  plot = gg_sig_samples,                        # the plot object
  width = 9,                                   # width in inches
  height = 12,                                   # height in inches
  dpi = 300                                     # resolution
)


#saving bubble plot data
library(openxlsx)
library(dplyr)

# -----------------------------
# Create workbook
# -----------------------------
wb <- createWorkbook()

# =========================================================
# 1. Differential abundance statistics (Wilcoxon + FDR)
# =========================================================
addWorksheet(wb, "Wilcoxon_DE_Stats")

writeData(
  wb,
  "Wilcoxon_DE_Stats",
  res_all %>%
    arrange(Sex, Species, Treatment)
)

# =========================================================
# 2. Effect size table (log2 fold-change vs Control)
# =========================================================
addWorksheet(wb, "Log2FC_EffectSizes")

writeData(
  wb,
  "Log2FC_EffectSizes",
  df_bubble %>%
    arrange(Sex, Species, Treatment)
)

# =========================================================
# 3. Final plotting dataset (filtered + significance)
# =========================================================
addWorksheet(wb, "BubblePlot_Data")

writeData(
  wb,
  "BubblePlot_Data",
  df_sample_plot %>%
    arrange(Sex, Species, Treatment, Sample)
)

# =========================================================
# 4. Species-level summary statistics
# =========================================================
species_summary <- df_species %>%
  group_by(Sex, Treatment, Species) %>%
  summarise(
    n_samples = n_distinct(Sample),
    mean_abundance = mean(Abundance, na.rm = TRUE),
    sd_abundance = sd(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

addWorksheet(wb, "Species_Summary")

writeData(
  wb,
  "Species_Summary",
  species_summary
)

# =========================================================
# 5. Save workbook
# =========================================================
saveWorkbook(
  wb,
  file = "bubble_plot_stats_females.xlsx",
  overwrite = TRUE
)




#################################################################################################
#individual bacteria reads bar plots
#################################################################################################

library(phyloseq)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
library(openxlsx)

# -----------------------------------------------------------------------------
# FUNCTION
# -----------------------------------------------------------------------------
analyze_species <- function(ps_obj,
                            species_name,
                            sex_label,
                            output_prefix){
  
  # ---------------------------------------------------------------------------
  # Convert to relative abundance (%)
  # ---------------------------------------------------------------------------
  ps_pct <- transform_sample_counts(
    ps_obj,
    function(x) 100 * x / sum(x)
  )
  
  # ---------------------------------------------------------------------------
  # Subset species
  # ---------------------------------------------------------------------------
  taxdf <- as.data.frame(tax_table(ps_pct))
  
  keep_taxa <- rownames(taxdf)[taxdf$Species == species_name]
  
  ps_species <- prune_taxa(keep_taxa, ps_pct)
  
  df <- psmelt(ps_species)
  
  # ---------------------------------------------------------------------------
  # Order groups
  # ---------------------------------------------------------------------------
  df$Treatment <- factor(
    df$Treatment,
    levels = c("Control", "HFiD", "DSS", "DSS+HFiD")
  )
  
  # ---------------------------------------------------------------------------
  # Overall Kruskal-Wallis
  # ---------------------------------------------------------------------------
  kw_results <- df %>%
    kruskal_test(Abundance ~ Treatment)
  
  # ---------------------------------------------------------------------------
  # Pairwise Dunn test
  # ---------------------------------------------------------------------------
  pairwise_results <- df %>%
    dunn_test(
      Abundance ~ Treatment,
      p.adjust.method = "BH"
    )
  
  # ---------------------------------------------------------------------------
  # Add significance labels
  # ---------------------------------------------------------------------------
  pairwise_results <- pairwise_results %>%
    mutate(
      significance = case_when(
        p.adj <= 0.0001 ~ "****",
        p.adj <= 0.001  ~ "***",
        p.adj <= 0.01   ~ "**",
        p.adj <= 0.05   ~ "*",
        TRUE            ~ "ns"
      )
    )
  
  # ---------------------------------------------------------------------------
  # Keep only significant comparisons
  # ---------------------------------------------------------------------------
  significant_only <- pairwise_results %>%
    filter(p.adj <= 0.05)
  
  # ---------------------------------------------------------------------------
  # Save statistics to Excel
  # ---------------------------------------------------------------------------
  wb <- createWorkbook()
  
  addWorksheet(wb, "Raw_Data")
  writeData(wb, "Raw_Data", df)
  
  addWorksheet(wb, "Kruskal_Wallis")
  writeData(wb, "Kruskal_Wallis", kw_results)
  
  addWorksheet(wb, "Pairwise_Dunn")
  writeData(wb, "Pairwise_Dunn", pairwise_results)
  
  addWorksheet(wb, "Significant_Only")
  writeData(wb, "Significant_Only", significant_only)
  
  saveWorkbook(
    wb,
    paste0(output_prefix, "_stats.xlsx"),
    overwrite = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------
  cell_palette <- c(
    "Control"   = "#2C7FB8",
    "HFiD"      = "#D9A441",
    "DSS"       = "#B96A9B",
    "DSS+HFiD"  = "#2FA187"
  )
  
  p <- ggplot(df,
              aes(x = Treatment,
                  y = Abundance,
                  fill = Treatment)) +
    
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,
      alpha = 0.85,
      color = "black",
      linewidth = 0.4
    ) +
    
    geom_jitter(
      width = 0.15,
      size = 1.3,
      alpha = 0.6,
      color = "black"
    ) +
    
    scale_fill_manual(values = cell_palette) +
    
    labs(
      x = NULL,
      y = "Relative abundance (%)",
      title = paste0(species_name, " (", sex_label, ")")
    ) +
    
    theme_classic(base_size = 12) +
    
    theme(
      legend.position = "none",
      plot.title = element_text(
        face = "italic",
        size = 12,
        hjust = 0.5
      ),
      axis.title.y = element_text(size = 11),
      axis.text = element_text(size = 10, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4)
    ) +
    
    stat_compare_means(
      method = "kruskal.test",
      label.y = max(df$Abundance) * 1.05
    )
  
  ggsave(
    paste0(output_prefix, ".pdf"),
    p,
    width = 3.6,
    height = 4,
    device = cairo_pdf
  )
  
  return(list(
    raw_data = df,
    kruskal = kw_results,
    pairwise = pairwise_results,
    significant = significant_only
  ))
}

#################################################################################################
# run function
#################################################################################################

# MALES
male_results <- analyze_species(
  ps_obj = ps_M_rar,
  species_name = "Clostridium perfringens",
  sex_label = "males",
  output_prefix = "Clostridium_perfringens_males"
)

# FEMALES
female_results <- analyze_species(
  ps_obj = ps_F_rar,
  species_name = "Clostridium perfringens",
  sex_label = "females",
  output_prefix = "Clostridium_perfringens_females"
)
