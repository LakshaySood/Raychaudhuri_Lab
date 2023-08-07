library(tidyverse)
library(hlabud)
library(stringr)

# What I need to do for each sample: 
# Datasheet 3 contains the true identity of each allele, and the identity of each allele 
#### that each algorithm predicted 
# So, for each sample, I have to fetch the true allele AA sequence and the predicted AA sequence
#### (or maybe the nucleotide sequences and translate them myself)
# Then, I have to compare each AA predicted 

seq_extract <- function(gene_matrix, alle) {
  if (is.na(alle)) {
    return(NA)
  }
  seq <- filter(gene_matrix$sequences, str_detect(allele, fixed(alle)))
  seq <- seq[1, ]$seq # This works because data is at 2-field resolution at minimum, and 3rd field is synonymous nucleotide variants
                      # so that there's no difference in amino acids 
  return(seq)
}

accuracy <- function(seq1, seq2) { 
  if (is.na(seq1)|is.na(seq2)) {
    return(0)
  }
  
  nchar1 <- nchar(seq1)
  nchar2 <- nchar(seq2)
  
  seq1 <- strsplit(seq1, "")[[1]]
  seq2 <- strsplit(seq2, "")[[1]]
  hamming_dist <- sum(seq1 != seq2)
  return(1 - (hamming_dist/nchar1)) 
}

setwd("/Users/lakshaysood/Desktop/Raychaudhuri_Lab/hla_benchmarking")
data <- read_csv("Prediction of HLA genotypes from single-cell transcriptome data.csv")
samples <- unique(data$sample)
gene_names <- c("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")
col_titles_1 <- paste0(gene_names, "_1", collapse = "|")
col_titles_2 <- paste0(gene_names, "_2", collapse = "|")
retained_cols <- grep(paste(col_titles_1, col_titles_2, sep = "|"), names(data))
data <- data[, c(1, 2, retained_cols)]
gene_matrices <- lapply(gene_names, function(gene) hla_alignments(gene, verbose = TRUE))
names(gene_matrices) <- gene_names

raw_results <- lapply(samples, function(sample_value) {
  sample_data <- filter(data, sample == sample_value)
  if (nrow(sample_data) == 1) {return(NA)} 
  sample_dist <- list()
  for (i in 3:ncol(sample_data)) {
    ground_truth <- sample_data[sample_data$genotyper == "Ground truth", i] %>% 
      pull()
    preds <- sample_data[2:nrow(sample_data), i] %>% pull()
    gene <- sub("_.", "", names(sample_data)[i])
    ground_truth_seq <- seq_extract(gene_matrices[[gene]], ground_truth)
    pred_seq <- lapply(preds, function(pred) seq_extract(gene_matrices[[gene]], pred))
    acc <- lapply(pred_seq, function(pred) accuracy(ground_truth_seq, pred))
    names(acc) <- sample_data[2:nrow(sample_data), "genotyper"] %>% pull()
    sample_dist[[i - 2]] <- acc 
  }
  names(sample_dist) <- names(sample_data)[3:ncol(sample_data)]
  return(sample_dist)
})

names(raw_results) <- samples

summarized_results_vec <- unlist(raw_results, recursive = TRUE)
summarized_results <- as.data.frame(summarized_results_vec) %>% 
  rownames_to_column(var = "Name") %>% 
  separate("Name", into = c("Sample", "Allele", "Software"), sep = "\\.") %>% 
  separate("Allele", into = c("Gene", "Pos"), sep = "_") %>% 
  select(c(-Sample, -Pos)) %>% 
  group_by(Gene, Software) %>% 
  summarize(mean_accuracy = mean(summarized_results_vec), .groups="keep")  %>% 
  pivot_wider(names_from = "Software", values_from = "mean_accuracy")
 
#####################################################################
# Question: what parts of the sequence does software tend to predict best? 
# Each position can be marked "accurate" or "inaccurate" 
# Then you compute average accuracy at each position 

for (seq in data$A_1) {
  print(seq)
  if (nchar(seq_extract(gene_matrices$A, seq)) != 387) {
      print(c(seq, nchar(seq_extract(gene_matrices$A, seq))))
  }
}

# Question: HLAminer seems to be making multiple predictions for certain alleles. 
