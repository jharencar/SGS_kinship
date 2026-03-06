# SGS kinship: 2019 & 2022 microsatellite data
# Streptanthus glandulosus (SGS) – kinship heatmap by Site × Year

# Load libraries
library(adegenet)
library(tidyverse)
library(stringr)

## SGS ####
# Read the CSV file (new format has Site, year, and same 7 microsat loci)
# Original format: Code, Site, Sn50..Sn1618 (no year). New: Code, Site, year, Sn50..Sn1618
data <- read.csv("SGS_2019_2022_microsats.csv")
# First column is individual ID (Code); rename in case CSV has BOM
names(data)[1] <- "PlantID"

#### genind conversion ####
# New CSV has year in column 3; genotype columns start at 4 (pairs: 4-5 Sn50, 6-7 Sn347, ... 16-17 Sn1618)
# Combine each pair of allele columns and create Site_Year population labels
geno_combined <- data %>%
  mutate(
    Pop     = paste(Site, year, sep = "_"),
    Sn50    = paste(.[[4]], .[[5]], sep = ","),
    Sn347   = paste(.[[6]], .[[7]], sep = ","),
    Sn262   = paste(.[[8]], .[[9]], sep = ","),
    Sn463   = paste(.[[10]], .[[11]], sep = ","),
    Sn558b  = paste(.[[12]], .[[13]], sep = ","),
    Sn1015  = paste(.[[14]], .[[15]], sep = ","),
    Sn1618  = paste(.[[16]], .[[17]], sep = ",")
  ) %>%
  select(PlantID, Site, year, Pop, Sn50, Sn347, Sn262, Sn463, Sn558b, Sn1015, Sn1618)

# Set rownames as individual IDs
rownames(geno_combined) <- geno_combined$PlantID

#### Missing data check per microsatellite ####
# Compute % missing per marker (0,0 or NA = missing). Check 2022 subset for Sn463; exclude if < 80% call rate.
markers <- c("Sn50", "Sn347", "Sn262", "Sn463", "Sn558b", "Sn1015", "Sn1618")
missing_pct <- geno_combined %>%
  summarise(across(all_of(markers), ~ 100 * mean(.x == "0,0" | is.na(.x)), .names = "missing_{.col}"))
missing_pct_2022 <- geno_combined %>%
  filter(year == 2022) %>%
  summarise(across(all_of(markers), ~ 100 * mean(.x == "0,0" | is.na(.x)), .names = "missing_{.col}"))
cat("Percentage missing per marker (all data):\n")
print(round(t(missing_pct), 1))
cat("\nPercentage missing per marker (2022 only):\n")
print(round(t(missing_pct_2022), 1))
sn463_call_rate_2022 <- 100 - as.numeric(missing_pct_2022$missing_Sn463)
if (sn463_call_rate_2022 < 80) {
  cat("\nSn463 has", round(sn463_call_rate_2022, 1), "% call rate in 2022 (< 80%). Excluding Sn463 from analysis.\n")
  geno_combined <- geno_combined %>% select(-Sn463)
  markers <- setdiff(markers, "Sn463")
} else {
  cat("\nSn463 has", round(sn463_call_rate_2022, 1), "% call rate in 2022. Keeping Sn463.\n")
}

# Genotype columns only (no ID/pop columns) for genind
geno_only <- geno_combined %>%
  select(all_of(markers))

# Convert to genind; populations = Site_Year
SGS.genind <- df2genind(geno_only,
                        sep = ",",
                        ploidy = 2,
                        NA.char = "NA",
                        pop = geno_combined$Pop)

#### kinship matrix ####
library(SNPRelate)

# Convert to genind for kinship: use only genotype columns (Site_Year already in geno_combined$Pop).
microsat_data <- geno_combined %>% select(all_of(markers))

# Create a matrix where each cell contains the genotype string (e.g., "134,134").
genotypes_matrix <- as.matrix(microsat_data)

# Handle the NA values for proper parsing by df2genind.
genotypes_matrix[genotypes_matrix == "0,0"] <- "NA"

# Create the genind object (pop = Site_Year for grouping in heatmap).
genind_obj <- df2genind(genotypes_matrix, sep = ",", pop = geno_combined$Pop)

# Method: Using SNPRelate – converts multi-allelic microsat data to dosage matrix,
# creates temporary GDS, then IBS kinship.

genind_to_dosage <- function(genind_obj) {
  dosage_matrix <- genind_obj@tab
  dosage_matrix <- dosage_matrix[, colSums(dosage_matrix, na.rm = TRUE) > 0]
  return(t(dosage_matrix))
}

dosage_mat <- genind_to_dosage(genind_obj)

snpgdsCreateGeno(gds.fn = "microsat.gds",
                 genmat = dosage_mat,
                 snp.chromosome = rep(1, nrow(dosage_mat)),
                 snpfirstdim = TRUE)

gds_file <- snpgdsOpen("microsat.gds")

kinship_matrix_snpgds <- snpgdsIBS(gds_file,
                                   autosome.only = FALSE,
                                   snp.id = NULL,
                                   sample.id = NULL,
                                   remove.monosnp = TRUE,
                                   maf = 0.05,
                                   missing.rate = 0.05)

kinship_mat <- kinship_matrix_snpgds$ibs
dimnames(kinship_mat) <- list(indNames(genind_obj), indNames(genind_obj))

snpgdsClose(gds_file)
file.remove("microsat.gds")

# Reorder by Site_Year then ID so samples from the same pop group together.
pop_vec <- pop(genind_obj)
ordered_ids <- order(pop_vec, rownames(kinship_mat))
kinship_mat_ordered <- kinship_mat[ordered_ids, ordered_ids]
site_year_labels <- as.character(pop_vec[ordered_ids])

# Gaps: draw a line after each row/col where the population changes.
gap_idx <- which(site_year_labels[-1] != site_year_labels[-length(site_year_labels)])

# Label each population block once, centered in the block.
block_starts <- which(c(TRUE, site_year_labels[-1] != site_year_labels[-length(site_year_labels)]))
block_ends <- c(block_starts[-1] - 1L, length(site_year_labels))
block_labels <- character(length(site_year_labels))
for (i in seq_along(block_starts)) {
  mid <- round((block_starts[i] + block_ends[i]) / 2)
  block_labels[mid] <- site_year_labels[block_starts[i]]
}

library(pheatmap)

pheatmap(kinship_mat_ordered,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_row = gap_idx,
         gaps_col = gap_idx,
         labels_row = block_labels,
         labels_col = block_labels,
         main = "SGS Kinship Matrix (Site × Year)",
         fontsize = 8,
         filename = "SGS_kinship_heatmap.pdf")
cat("Kinship heatmap saved to SGS_kinship_heatmap.pdf\n")
