#!/usr/bin/env Rscript

args = commandArgs(TRUE)
library(dada2); packageVersion("dada2")
library(dplyr)
library(tibble)

sts <- list.files(path = args[1], pattern = "*yycH_seqtab_from_dada2.rds", full.names = TRUE)
st.all <- mergeSequenceTables(tables = sts)
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
set.seed(1)
tax <- assignTaxonomy(seqtab.nochim, args[2], multithread=TRUE, taxLevels = c("Seq_number"), minBoot=35)
write.csv2(tax, args[3])
saveRDS(tax, args[4])

save.image(args[5])

seqtab.nochim <- t(seqtab.nochim)
seqtab.nochim_df <- as.data.frame(seqtab.nochim)
seqtab.nochim_df <- rownames_to_column(seqtab.nochim_df, var = "ASV")

tax_df <- as.data.frame(tax)
tax_df <- as_tibble(rownames_to_column(tax_df, var = "ASV"))
seq_taxa_spec <- dplyr::left_join(tax_df, seqtab.nochim_df, by = 'ASV')
write.csv2(seq_taxa_spec, args[6])

