source("https://raw.githubusercontent.com/ssi-dk/epidome/master/scripts/epidome_functions.R?token=AHCBOOG2CAJ574PBGRRLN726IKSPK")
setwd("/Volumes/data/MPV/projects/git.repositories/epidome/")

### Load amplicon table for classification ###
ST_amplicon_table = read.table("DB/epidome_ST_amplicon_frequencies.txt",sep = "\t")

### Load dada2 output for the two primers ###
epi01_table_full = read.table("example_data/epi01_count_table_renamed.csv",sep = ";",header=TRUE,row.names=1)
epi02_table_full = read.table("example_data/epi02_count_table_renamed.csv",sep = ";",header=TRUE,row.names=1)

### Load metadata table
metadata_table = read.table("example_data/sample_metadata.txt",header=TRUE,row.names=1)
metadata_table$patient.sample.site = paste0(as.vector(metadata_table$patient.ID),' ',as.vector(metadata_table$sample.site))
metadata_table$sample.site[metadata_table$sample.site=="elbow"] = "arm"

### Setup an object for easy handling of epidome data
epidome_object = setup_epidome_object(epi01_table_full,epi02_table_full,metadata_table = metadata_table)

### Check if number of sequences from each primer for each samples match up approximately ###
compare_primer_output(epidome_object)


### Filter lowcount samples (removes any sample that has less than X sequences from one of the two primer sets, here 500) ###
epidome_filtered_samples = filter_lowcount_samples_epidome(epidome_object,500,500)

### Combine ASVs from dada output ###
epidome_ASV_combined = combine_ASVs_epidome(epidome_filtered_samples)

epidome_object_mock = prune_by_variable_epidome(epidome_ASV_combined,"sample.type",c("Mock community"))
epidome_object_clinical = prune_by_variable_epidome(epidome_ASV_combined,"sample.type",c("Clinical"))


pt_site_tbl = table(epidome_object_clinical$metada$patient.sample.site)

include_values = names(pt_site_tbl)[which(pt_site_tbl==2)]
