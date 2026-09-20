#!/usr/bin/env Rscript
# 06_03_assign_genotype_by_pearson.R — map Souporcell clusters to known genotypes.
#
# Correlates each Souporcell cluster genotype vector (cluster_genotypes.vcf) with
# each known genotype in the reference VCF (Pearson, over shared SNPs), assigns each
# genotype to its best-correlated cluster, and rescues B73 (the reference genome,
# nearly all 0/0) by hom-ref fraction. Called by 06_01_souporcell.sh.
# Usage : Rscript 06_03_assign_genotype_by_pearson.R <reference.vcf> <cluster_genotypes.vcf> <outdir>
# Output: <outdir>/ref_clust_pearson_correlations.v2.tsv, ref_clust_pearson_correlation.v2.png,
#         Genotype_ID_key.v2.txt (read together with clusters.tsv by Fig. 4H to K)

library(tidyr)
suppressMessages(library(tidyverse))
library(dplyr)
suppressMessages(library(vcfR))
suppressMessages(library(lsa))
suppressMessages(library(ComplexHeatmap))

########## Set up paths and variables ##########
args <- commandArgs(T)

reference_vcf <- as.character(args[1])   # known-genotype VCF filtered to the pooled genotypes
cluster_vcf   <- as.character(args[2])   # Souporcell cluster_genotypes.vcf
outdir        <- as.character(args[3])   # Souporcell output directory


########## Set up functions ##########
##### Calculate DS from GP if genotypes in that format #####
calculate_DS <- function(GP_df){
  columns <- c()
  for (i in 1:ncol(GP_df)){
    columns <- c(columns, paste0(colnames(GP_df)[i],"-0"), paste0(colnames(GP_df)[i],"-1"), paste0(colnames(GP_df)[i],"-2"))
  }
  df <- GP_df
  colnames(df) <- paste0("c", colnames(df))
  colnames_orig <- colnames(df)
  for (i in 1:length(colnames_orig)){
    df <- separate(df, sep = ",", col = colnames_orig[i], into = columns[(1+(3*(i-1))):(3+(3*(i-1)))])
  }
  df <- mutate_all(df, function(x) as.numeric(as.character(x)))
  for (i in 1: ncol(GP_df)){
    GP_df[,i] <- df[,(2+((i-1)*3))] + 2* df[,(3+((i-1)*3))]
  }
  return(GP_df)
}

pearson_correlation <- function(df, ref_df, clust_df){
  for (col in colnames(df)){
    for (row in rownames(df)){
      df[row,col] <- cor(as.numeric(pull(ref_df, col)), 
                         as.numeric(pull(clust_df, row)), 
                         method = "pearson", use = "complete.obs")
    }
  }
  return(df)
}


########## Read in vcf files for each of three non-reference genotype softwares ##########
ref_geno    <- read.vcfR(reference_vcf)
cluster_geno <- read.vcfR(cluster_vcf)

########## Convert to tidy data frame ##########
####### Identify which genotype FORMAT to use #######
##### Cluster VCF #####
### Check for each of the different genotype formats ##
## DS ##
format_clust=NA
cluster_geno_tidy <- as_tibble(extract.gt(element = "DS", cluster_geno, IDtoRowNames = F))

if (!all(colSums(is.na(cluster_geno_tidy)) == nrow(cluster_geno_tidy))){
  message("Found DS genotype format in cluster vcf. Will use that metric for cluster correlation.")
  format_clust = "DS"
}

## GT ##
if (is.na(format_clust)){
  cluster_geno_tidy <- as_tibble(extract.gt(element = "GT",cluster_geno, IDtoRowNames = F))
  if (!all(colSums(is.na(cluster_geno_tidy)) == nrow(cluster_geno_tidy))){
    message("Found GT genotype format in cluster vcf. Will use that metric for cluster correlation.")
    format_clust = "GT"
    
    if (any(grepl("\\|",cluster_geno_tidy[1,]))){
      separator = "|"
      message("Detected | separator for GT genotype format in cluster vcf")
    } else if (any(grepl("/",cluster_geno_tidy[1,]))) {
      separator = "/"
      message("Detected / separator for GT genotype format in cluster vcf")
    } else {
      format_clust = NA
      message("Can't identify a separator for the GT field in cluster vcf, moving on to using GP.")
    }
    
    cluster_geno_tidy <- as_tibble(lapply(cluster_geno_tidy, function(x) {gsub(paste0("0",separator,"0"),0, x)}) %>%
                                     lapply(., function(x) {gsub(paste0("0",separator,"1"),1, x)}) %>%
                                     lapply(., function(x) {gsub(paste0("1",separator,"0"),1, x)}) %>%
                                     lapply(., function(x) {gsub(paste0("1",separator,"1"),2, x)}))
    
  }
}


## GP ##
if (is.na(format_clust)){
  cluster_geno_tidy <- as_tibble(extract.gt(element = "GP",cluster_geno, IDtoRowNames =F))
  if (!all(colSums(is.na(cluster_geno_tidy)) == nrow(cluster_geno_tidy))){
    format_clust = "GP"
    cluster_geno_tidy <- calculate_DS(cluster_geno_tidy)
    message("Found GP genotype format in cluster vcf. Will use that metric for cluster correlation.")
    
  } else {
    print("Could not identify the expected genotype format fields (DS, GT or GP) in your cluster vcf. Please check the vcf file and make sure that one of the expected genotype format fields is included or run manually with your genotype format field of choice. Quitting")
    q()
  }
}


### Reference VCF ###
### Check for each of the different genotype formats ##
## DS ##
format_ref = NA
ref_geno_tidy <- as_tibble(extract.gt(element = "DS",ref_geno, IDtoRowNames = F))
if (!all(colSums(is.na(ref_geno_tidy)) == nrow(ref_geno_tidy))){
  message("Found DS genotype format in reference vcf. Will use that metric for cluster correlation.")
  format_ref = "DS"
}

## GT ##
if (is.na(format_ref)){
  ref_geno_tidy <- as_tibble(extract.gt(element = "GT",ref_geno, IDtoRowNames = F))
  if (!all(colSums(is.na(ref_geno_tidy)) == nrow(ref_geno_tidy))){
    message("Found GT genotype format in reference vcf. Will use that metric for cluster correlation.")
    format_ref = "GT"
    
    if (any(grepl("\\|",ref_geno_tidy[1,]))){
      separator = "|"
      message("Detected | separator for GT genotype format in reference vcf")
    } else if (any(grepl("/",ref_geno_tidy[1,]))) {
      separator = "/"
      message("Detected / separator for GT genotype format in reference vcf")
    } else {
      format_ref = NA
      message("Can't identify a separator for the GT field in reference vcf, moving on to using GP.")
    }
    
    ref_geno_tidy <- as_tibble(lapply(ref_geno_tidy, function(x) {gsub(paste0("0",separator,"0"),0, x)}) %>%
                                 lapply(., function(x) {gsub(paste0("0",separator,"1"),1, x)}) %>%
                                 lapply(., function(x) {gsub(paste0("1",separator,"0"),1, x)}) %>%
                                 lapply(., function(x) {gsub(paste0("1",separator,"1"),2, x)}))
    
  }
}

## GP ##
if (is.na(format_ref)){
  ref_geno_tidy <- as_tibble(extract.gt(element = "GP",ref_geno, IDtoRowNames = F))
  if (!all(colSums(is.na(ref_geno_tidy)) == nrow(ref_geno_tidy))){
    format_clust = "GP"
    ref_geno_tidy <- calculate_DS(ref_geno_tidy)
    message("Found GP genotype format in cluster vcf. Will use that metric for cluster correlation.")
    
  } else {
    print("Could not identify the expected genotype format fields (DS, GT or GP) in your cluster vcf. Please check the vcf file and make sure that one of the expected genotype format fields is included or run manually with your genotype format field of choice. Quitting")
    q()
  }
}



### Get SNP IDs that will match between reference and cluster ###
## Account for possibility that the ref or alt might be missing
if ((all(is.na(cluster_geno@fix[,'REF'])) & all(is.na(cluster_geno@fix[,'ALT']))) | (all(is.na(ref_geno@fix[,'REF'])) & all(is.na(ref_geno@fix[,'ALT'])))){
  message("The REF and ALT categories are not provided for the reference and/or the cluster vcf. Will use just the chromosome and position to match SNPs.")
  cluster_geno_tidy$ID <- paste0(cluster_geno@fix[,'CHROM'],":", cluster_geno@fix[,'POS'])
  ref_geno_tidy$ID <- paste0(ref_geno@fix[,'CHROM'],":", ref_geno@fix[,'POS'])
} else if (all(is.na(cluster_geno@fix[,'REF'])) | all(is.na(ref_geno@fix[,'REF']))){
  message("The REF categories are not provided for the reference and/or the cluster vcf. Will use the chromosome, position and ALT to match SNPs.")
  cluster_geno_tidy$ID <- paste0(cluster_geno@fix[,'CHROM'],":", cluster_geno@fix[,'POS'],"_", cluster_geno@fix[,'REF'])
  ref_geno_tidy$ID <- paste0(ref_geno@fix[,'CHROM'],":", ref_geno@fix[,'POS'],"_", ref_geno@fix[,'REF'])
} else if (all(is.na(cluster_geno@fix[,'ALT'])) | all(is.na(ref_geno@fix[,'ALT']))){
  message("The ALT categories are not provided for the reference and/or the cluster vcf. Will use the chromosome, position and REF to match SNPs.")
  cluster_geno_tidy$ID <- paste0(cluster_geno@fix[,'CHROM'],":", cluster_geno@fix[,'POS'],"_", cluster_geno@fix[,'ALT'])
  ref_geno_tidy$ID <- paste0(ref_geno@fix[,'CHROM'],":", ref_geno@fix[,'POS'],"_", ref_geno@fix[,'ALT'])
} else {
  message("Found REF and ALT in both cluster and reference genotype vcfs. Will use chromosome, position, REF and ALT to match SNPs.")
  cluster_geno_tidy$ID <- paste0(cluster_geno@fix[,'CHROM'],":", cluster_geno@fix[,'POS'],"_", cluster_geno@fix[,'REF'],"_", cluster_geno@fix[,'ALT'])
  ref_geno_tidy$ID <- paste0(ref_geno@fix[,'CHROM'],":", ref_geno@fix[,'POS'],"_", ref_geno@fix[,'REF'],"_", ref_geno@fix[,'ALT'])
}


### Update the vcf dfs to remove SNPs with no genotyopes
cluster_geno_tidy <- cluster_geno_tidy[colSums(!is.na(cluster_geno_tidy)) > 0]
ref_geno_tidy <- ref_geno_tidy[colSums(!is.na(ref_geno_tidy)) > 0]


########## Get a unique list of SNPs that is in both the reference and cluster genotypes ##########
locations  <- inner_join(ref_geno_tidy[,"ID"],cluster_geno_tidy[,"ID"])
locations <- locations[!(locations$ID %in% locations[duplicated(locations),]$ID),]

########## Keep just the SNPs that overlap ##########
ref_geno_tidy <- left_join(locations, ref_geno_tidy)
cluster_geno_tidy <- left_join(locations, cluster_geno_tidy)

########## Correlate all the cluster genotypes with the individuals genotyped ##########
##### Make a dataframe that has the clusters as the row names and the individuals as the column names #####
pearson_correlations <- as.data.frame(matrix(nrow = (ncol(cluster_geno_tidy) -1), ncol = (ncol(ref_geno_tidy) -1)))
colnames(pearson_correlations) <- colnames(ref_geno_tidy)[2:(ncol(ref_geno_tidy))]
rownames(pearson_correlations) <- colnames(cluster_geno_tidy)[2:(ncol(cluster_geno_tidy))]

ref_geno_tidy_tem <- ref_geno_tidy
cluster_geno_tidy_tem <- cluster_geno_tidy

# if ref in clusters, all position are going to be 0, and correlation will be un calculated
ref_geno_tidy_tem[,-c(1)] <- apply(ref_geno_tidy_tem[,-c(1)], 2, as.numeric)
cluster_geno_tidy_tem[,-c(1)] <- apply(cluster_geno_tidy_tem[,-c(1)], 2, as.numeric)

ref_geno_tidy_tem[,-c(1)] <- ref_geno_tidy_tem[,-c(1)]+1
cluster_geno_tidy_tem[,-c(1)] <- cluster_geno_tidy_tem[,-c(1)]+1

# Inject pseudo-SNPs (value 3) at 2% of sites in the B73 column: B73 is the
# reference genome, so its genotype vector is constant (all 0/0) and cor() would
# be undefined.
false_snp_value <- 3
if ("B73" %in% colnames(ref_geno_tidy_tem) && nrow(ref_geno_tidy_tem) >= 100) {
  n_false <- max(1, round(nrow(ref_geno_tidy_tem) * 0.02))
  Index <- sample(seq_len(nrow(ref_geno_tidy_tem)), size = n_false, replace = FALSE)
  ref_geno_tidy_tem$B73[Index] <- false_snp_value
}

pearson_correlations <- pearson_correlation(pearson_correlations, ref_geno_tidy_tem, cluster_geno_tidy_tem)
cluster <- data.frame("Cluster" = rownames(pearson_correlations))
pearson_correlations_out <- cbind(cluster, pearson_correlations)

########## Save the correlation dataframes ##########
write_delim(pearson_correlations_out, file = paste0(outdir,"/ref_clust_pearson_correlations.v2.tsv"), delim = "\t" )


########## Create correlation figures ##########
col_fun = colorRampPalette(c("white", "red"))(101)
# Replace NAs with 0 for visualization; disable clustering if any NAs remain
pearson_plot <- as.matrix(pearson_correlations)
has_na <- any(is.na(pearson_plot))
pearson_plot[is.na(pearson_plot)] <- 0
pPearsonCorrelations <- Heatmap(pearson_plot,
                                cluster_rows = !has_na,
                                cluster_columns = !has_na,
                                col = col_fun)

########## Save the correlation figures ##########
png(filename = paste0(outdir,"/ref_clust_pearson_correlation.v2.png"), width = 500)
print(pPearsonCorrelations)
dev.off()

########## Assign individual to cluster based on highest correlating individual ##########
key <- as.data.frame(matrix(nrow = ncol(pearson_correlations), ncol = 3))
colnames(key) <- c("Genotype_ID","Cluster_ID","Correlation")
key$Genotype_ID <- colnames(pearson_correlations)
for (id in key$Genotype_ID){
  col_vals <- pearson_correlations[, id]
  # Skip genotypes where all correlations are NA
  if (all(is.na(col_vals))) {
    key$Cluster_ID[which(key$Genotype_ID == id)] <- "unassigned"
    key$Correlation[which(key$Genotype_ID == id)] <- NA
    next
  }
  best_cluster <- rownames(pearson_correlations)[which.max(col_vals)]
  best_cor <- max(col_vals, na.rm = TRUE)
  row_vals <- pearson_correlations[best_cluster, ]
  if (!all(is.na(row_vals)) &&
      best_cor == max(row_vals, na.rm = TRUE)) {
    key$Cluster_ID[which(key$Genotype_ID == id)] <- best_cluster
    key$Correlation[which(key$Genotype_ID == id)] <- best_cor
  } else {
    key$Cluster_ID[which(key$Genotype_ID == id)] <- "unassigned"
    key$Correlation[which(key$Genotype_ID == id)] <- NA
  }
}

########## B73 rescue: assign via hom-ref fraction in cluster_genotypes.vcf ##########
# B73 is the reference genome — nearly all 0/0 at variant sites — so Pearson
# correlation fails (zero-variance vector). Instead, identify B73 by looking for
# the unassigned cluster with the highest fraction of hom-ref (0/0) genotype
# calls. Only applied when B73 is expected in the pool AND still unassigned.
if ("B73" %in% key$Genotype_ID && key$Cluster_ID[key$Genotype_ID == "B73"] == "unassigned") {
  message("B73 is unassigned after Pearson — attempting hom-ref rescue...")

  # Get clusters that are not already claimed by another genotype
  assigned_clusters <- key$Cluster_ID[key$Cluster_ID != "unassigned"]
  all_clusters <- colnames(cluster_geno_tidy)[colnames(cluster_geno_tidy) != "ID"]
  unassigned_clusters <- setdiff(all_clusters, assigned_clusters)

  if (length(unassigned_clusters) > 0) {
    # Compute hom-ref fraction per unassigned cluster from the already-loaded GT data
    homref_frac <- sapply(unassigned_clusters, function(cl) {
      gt_vals <- as.numeric(cluster_geno_tidy[[cl]])
      n_valid <- sum(!is.na(gt_vals))
      if (n_valid == 0) return(0)
      sum(gt_vals == 0, na.rm = TRUE) / n_valid
    })

    best_cluster <- unassigned_clusters[which.max(homref_frac)]
    best_frac <- max(homref_frac)

    # Also check margin over second-best unassigned cluster
    if (length(homref_frac) > 1) {
      sorted_fracs <- sort(homref_frac, decreasing = TRUE)
      margin <- sorted_fracs[1] - sorted_fracs[2]
    } else {
      margin <- best_frac
    }

    message(sprintf("  Best unassigned cluster: %s (%.1f%% hom-ref, margin=%.1f pp)",
                    best_cluster, best_frac * 100, margin * 100))

    # Assign B73 to the best unassigned cluster
    key$Cluster_ID[key$Genotype_ID == "B73"] <- best_cluster
    key$Correlation[key$Genotype_ID == "B73"] <- NA  # Not a Pearson value
    message(sprintf("  -> Assigned B73 to cluster %s via hom-ref rescue", best_cluster))
  } else {
    message("  No unassigned clusters available for B73 rescue.")
  }
}

########## Remainder rescue: assign leftover genotypes to leftover clusters ##########
# After B73 rescue, if any genotypes are still unassigned and there are unclaimed
# clusters, assign them greedily by best available Pearson correlation. Handles the
# case where cor() returns NA for one genotype (e.g., zero-variance reference column).
still_unassigned <- key$Genotype_ID[key$Cluster_ID == "unassigned"]
if (length(still_unassigned) > 0) {
  assigned_clusters <- key$Cluster_ID[key$Cluster_ID != "unassigned"]
  all_clusters <- rownames(pearson_correlations)
  remaining_clusters <- setdiff(all_clusters, assigned_clusters)

  if (length(still_unassigned) == 1 && length(remaining_clusters) == 1) {
    # Only one genotype and one cluster left — assign directly
    message(sprintf("Remainder rescue: assigning %s to cluster %s (only remaining pair)",
                    still_unassigned, remaining_clusters))
    key$Cluster_ID[key$Genotype_ID == still_unassigned] <- remaining_clusters
    key$Correlation[key$Genotype_ID == still_unassigned] <- NA
  } else if (length(remaining_clusters) > 0) {
    # Greedy: for each unassigned genotype, pick the remaining cluster with best cor
    for (geno in still_unassigned) {
      cors <- pearson_correlations[remaining_clusters, geno, drop = FALSE]
      if (all(is.na(cors))) {
        message(sprintf("Remainder rescue: %s has all-NA correlations, assigning to first remaining cluster %s",
                        geno, remaining_clusters[1]))
        key$Cluster_ID[key$Genotype_ID == geno] <- remaining_clusters[1]
        key$Correlation[key$Genotype_ID == geno] <- NA
      } else {
        best <- rownames(cors)[which.max(cors[,1])]
        message(sprintf("Remainder rescue: assigning %s to cluster %s (cor=%.4f)",
                        geno, best, max(cors[,1], na.rm = TRUE)))
        key$Cluster_ID[key$Genotype_ID == geno] <- best
        key$Correlation[key$Genotype_ID == geno] <- max(cors[,1], na.rm = TRUE)
      }
      remaining_clusters <- setdiff(remaining_clusters, key$Cluster_ID[key$Genotype_ID == geno])
    }
  }
}

write_delim(key, file = paste0(outdir,"/Genotype_ID_key.v2.txt"), delim = "\t")