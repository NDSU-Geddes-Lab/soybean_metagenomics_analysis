# Preparation of metagenomics phyloseq objects
# Imports the Kraken/Bracken BIOM table, merges metadata, and creates bacterial and fungal phyloseq objects

# Load required packages
library(phyloseq)

# Import Kraken BIOM table
merged_metagenomes <- import_biom("kraken_table_030526.biom")


# Clean taxonomy labels
merged_metagenomes@tax_table@.Data <-substring(merged_metagenomes@tax_table@.Data, 4)

colnames(merged_metagenomes@tax_table@.Data) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")



# Merge with metadata
metadata <- read.csv("metagenomics_metadata.csv",row.names = 1)
sample_data(merged_metagenomes) <- sample_data(metadata)

# Save metagenomics phyloseq object
saveRDS(merged_metagenomes,"PS.merged_metagenomes.rds")

# Subset bacterial community
ps_metag_b <- subset_taxa(merged_metagenomes, Kingdom == "Bacteria")

ps_metag_b <- prune_taxa(taxa_sums(ps_metag_b) > 0, ps_metag_b)

ps_metag_b <- prune_samples(sample_sums(ps_metag_b) > 0,ps_metag_b)

saveRDS(ps_metag_b,"PS.metagenomics_bacteria.rds")


# Subset fungal community
ps_metag_f <- subset_taxa(merged_metagenomes,
  Kingdom == "Eukaryota" &
    Phylum %in% c(
      "Ascomycota",
      "Basidiomycota",
      "Mucoromycota",
      "Microsporidia"
    )
)

ps_metag_f <- prune_taxa(taxa_sums(ps_metag_f) > 0, ps_metag_f)

ps_metag_f <- prune_samples(sample_sums(ps_metag_f) > 0, ps_metag_f)

saveRDS(ps_metag_f,"PS.metagenomics_fungi.rds")
