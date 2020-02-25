source("https://raw.githubusercontent.com/ssi-dk/epidome/master/scripts/epidome_functions.R")
setwd("/Volumes/data/MPV/projects/git.repositories/epidome/")

### Load amplicon table for classification ###
ST_amplicon_table = read.table("DB/epidome_ST_amplicon_frequencies.txt",sep = "\t")

### Load dada2 output for the two primers ###
epi01_mock_table = read.table("example_data/190920_run2_G216_seqtab_from_dada2.csv.classified.csv",sep = ";",header=TRUE,row.names=1)
epi02_mock_table = read.table("example_data/190920_run2_yycH_seqtab_from_dada2.csv.classified.csv",sep = ";",header=TRUE,row.names=1)
epi01_clinical_table = read.table("example_data/190920_run1_G216_seqtab_from_dada2.csv.classified.csv",sep = ";",header=TRUE,row.names=1)
epi02_clinical_table = read.table("example_data/190920_run1_yycH_seqtab_from_dada2.csv.classified.csv",sep = ";",header=TRUE,row.names=1)

epi01_table = combine_ASV_tables(epi01_mock_table,epi01_clinical_table)
epi02_table = combine_ASV_tables(epi02_mock_table,epi02_clinical_table)

### Load metadata table
metadata_table = read.table("example_data/article_metadata.txt",header=TRUE,row.names=1)

### Setup an object for easy handling of epidome data
epidome_object = setup_epidome_object(epi01_table,epi02_table,metadata_table = metadata_table)
#epidome_object = setup_epidome_object(epi01_table,epi02_table)

### Check if number of sequences from each primer for each samples match up approximately ###
compare_primer_output(epidome_object)
compare_primer_output(epidome_object,"sample.type")
compare_primer_output(epidome_object,"patient.sample.site")


### Filter lowcount samples (removes any sample that has less than X sequences from one of the two primer sets, here 500) ###
epidome_filtered_samples = filter_lowcount_samples_epidome(epidome_object,500,500)

### Combine ASVs from dada output ###
epidome_ASV_combined = combine_ASVs_epidome(epidome_filtered_samples)

epidome_object_mock = prune_by_variable_epidome(epidome_ASV_combined,"sample.type",c("Mock community"))
epidome_object_clinical = prune_by_variable_epidome(epidome_ASV_combined,"sample.type",c("Clinical"))


pt_site_tbl = table(epidome_object_clinical$metada$patient.sample.site)

include_values = names(pt_site_tbl)[which(pt_site_tbl==2)]

epidome_object_clinical_pruned = prune_by_variable_epidome(epidome_object_clinical,"patient.sample.site",include_values)
epidome_clinical_pruned_norm = normalize_epidome_object(epidome_object_clinical_pruned)

sample_site = epidome_clinical_pruned_norm$metadata$sample.site
pca_pt = plot_PCA_epidome(epidome_clinical_pruned_norm,"patient.ID",RColorBrewer::brewer.pal(12,"Paired")[c(1:10,12)],plot_ellipse = F)
pca_pt = pca_pt + geom_point(size=2,aes(shape=sample_site)) + scale_shape_manual(values=c(3,1))
pca_pt + geom_point()
pca_pt + geom_point(size=2)
pca_site = plot_PCA_epidome(epidome_clinical_pruned_norm,"sample.site",c())
pca_site = pca_site + geom_point(size=2)

color_variable = "patient.ID"
m = epidome_clinical_pruned_norm$metadata
color_variable_factor = m[,which(epidome_clinical_pruned_norm$meta_variables==color_variable)]
data_combined = rbind(epidome_clinical_pruned_norm$p1_table,epidome_clinical_pruned_norm$p2_table)
pca = prcomp(t(data_combined))
plot_df = data.frame(pca$x)
color_vector = setup_colors(levels(color_variable_factor),colors)
labels = c(paste0("PC1 [",sprintf("%.1f",explvar(pca)[1]),"%]"),paste0("PC2, [",sprintf("%.1f",explvar(pca)[2]),"%]"))
pca_pt = ggplot(as.data.frame(pca$x),aes(x=PC1,y=PC2,color = color_variable_factor,shape=sample_site)) + labs(color = color_variable) + geom_point(size=2, alpha=1) + scale_colour_manual(values = color_vector) + xlab(labels[1]) + ylab(labels[2]) + theme_bw()+ scale_shape_manual(values=c(16,17))



grid.arrange(pca_site,pca_pt,ncol=2)


count_table_mock = classify_epidome(epidome_object_mock,ST_amplicon_table)
count_table_mock = count_table_mock[,order(colnames(count_table_mock))]

barplot_mock = make_barplot_epidome(count_table_mock)


count_table_clinical = classify_epidome(epidome_object_clinical_pruned,ST_amplicon_table)
count_table_clinical = count_table_clinical[,order(colnames(count_table_clinical))]
barplot_clinical = make_barplot_epidome(count_table_clinical)

grid.arrange(barplot_mock, barplot_clinical, ncol=2)



color_table = read.table("example_data/epidome_color_codes.txt",sep = "\t",header=T,comment.char = "")
color_table = rbind(color_table,data.frame(c()))

mock_STs = c(2,5,14,87,215,218,"Unclassified")
count_mock_fixed = rbind(count_table_mock[rownames(count_table_mock) %in% mock_STs,],colSums(count_table_mock[!rownames(count_table_mock) %in% mock_STs,]))
rownames(count_mock_fixed)[nrow(count_mock_fixed)] = "Other"


dd<-apply(count_mock_fixed, 2, function(x) x/sum(x)*100)
count_mock_fixed<-as.data.frame(dd)

count_mock_fixed$ST = rownames(count_mock_fixed)

color_table = read.table("example_data/epidome_color_codes.txt",sep="\t",comment.char = "",header=TRUE)

melt_df = melt(count_mock_fixed)
colnames(melt_df) = c("ST","Sample","Count")
ST_levels = c(2,5,14,87,215,218,"Other","Unclassified")
melt_df$ST = factor(melt_df$ST, levels=ST_levels)
ST = unlist(lapply(ST_levels, function(x) if (x %in% color_table$ST) { as.vector(color_table$hex.code)[which(color_table$ST==x)] } else {"Missing"}))
ST[which(ST=="Missing")] = c("#f5ed5d","#e8b099")

p = ggplot() + geom_bar(aes(y = Count, x = Sample, fill = ST), data = melt_df, stat="identity") + scale_fill_manual(values = ST) + theme_bw() + ylab(ylabel) + theme(axis.text.x = element_blank()) + xlab("")
barplot_mock_fixedcol = p


#### clinical mock with article colors
rownames(count_table_clinical)[which(rownames(count_table_clinical)=='-')] = "Novel"
mock_STs = rownames(count_table_clinical)[order(rowSums(count_table_clinical),decreasing = T)][1:12]
count_mock_fixed = rbind(count_table_clinical[rownames(count_table_clinical) %in% mock_STs,],colSums(count_table_clinical[!rownames(count_table_clinical) %in% mock_STs,]))
rownames(count_mock_fixed)[nrow(count_mock_fixed)] = "Other"


dd<-apply(count_mock_fixed, 2, function(x) x/sum(x)*100)
count_mock_fixed<-as.data.frame(dd)

count_mock_fixed$ST = rownames(count_mock_fixed)


melt_df = melt(count_mock_fixed)
colnames(melt_df) = c("ST","Sample","Count")
ST_numbers = mock_STs[which(!mock_STs %in% c("Unclassified","Novel"))]
ST_levels = c(ST_numbers[order(as.numeric(ST_numbers))],c("Novel","Unclassified","Other"))
melt_df$ST = factor(melt_df$ST, levels=ST_levels)
ST = unlist(lapply(ST_levels, function(x) if (x %in% color_table$ST) { as.vector(color_table$hex.code)[which(color_table$ST==x)] } else {"Missing"}))
ST[which(ST=="Missing")] = c("#a1984d","#ba291c","#89c981")

p = ggplot() + geom_bar(aes(y = Count, x = Sample, fill = ST), data = melt_df, stat="identity") + scale_fill_manual(values = ST) + theme_bw() + ylab(ylabel) + theme(axis.text.x = element_blank()) + xlab("")
barplot_clinical_fixedcol = p




