source("https://raw.githubusercontent.com/ssi-dk/epidome/master/scripts/epidome_functions.R?token=AHCBOOG2CAJ574PBGRRLN726IKSPK")
setwd("/Volumes/data/MPV/projects/git.repositories/epidome/")

### Load amplicon table for classification ###
ST_amplicon_table = read.table("DB/epidome_ST_amplicon_frequencies.txt",sep = "\t")

### Load dada2 output for the two primers ###
epi01_table_full = read.table("example_data/epi01_count_table.csv",sep = ";",header=TRUE,row.names=1)
epi02_table_full = read.table("example_data/epi02_count_table.csv",sep = ";",header=TRUE,row.names=1)

### Load metadata table
metadata_table = read.table("example_data/sample_metadata.txt",header=TRUE,row.names=1)
metadata_table$patient.sample.site = paste0(as.vector(metadata_table$patient.ID),' ',as.vector(metadata_table$sample.site))
metadata_table$sample.site[metadata_table$sample.site=="elbow"] = "arm"

### Setup an object for easy handling of epidome data
epidome_object = setup_epidome_object(epi01_table_full,epi02_table_full,metadata_table = metadata_table)

### Check if number of sequences from each primer for each samples match up approximately ###
compare_primer_output(epidome_object)



#### Data manipulation of various sorts ####

### Combine ASVs from dada output ###
epidome_ASV_combined = combine_ASVs_epidome(epidome_object)

### Filter lowcount samples (removes any sample that has less than X sequences from one of the two primer sets, here 500) ###
epidome_filtered_samples = filter_lowcount_samples_epidome(epidome_object,500,500)

### Filter lowcount ASVs (removes any ASV that does not appear with more than X% abundance in any sample) ###
epidome_filtered_ASVs = filter_lowcount_ASVs_epidome(epidome_object,percent_threshold = 1)

### Normalize data to percent ###
epidome_object_normalized = normalize_epidome_object(epidome_filtered_samples)

### Data stratification based on metadata variables
epidome_object_clinical = prune_by_variable_epidome(epidome_object,"sample.type",c("Clinical"))
epidome_object_mock = prune_by_variable_epidome(epidome_object,"sample.type",c("Mock community"))



#### Plots and figures ###

### Setup a data frame with sequence classification based on the two primer outputs. Combine ASVs first ###
epidome_ASV_combined = combine_ASVs_epidome(epidome_object)
count_table = classify_epidome(epidome_ASV_combined,ST_amplicon_table)
### Make barplot based on classification. Set reorder=TRUE to order samples based on Bray Curtis dissimilarity and/or set normalize=FALSE to not normalize to percent ###
make_barplot_epidome(count_table,reorder=FALSE,normalize=TRUE)

### make PCA plots. Color according to variable in metadata and (optional) indicate colors to use. Set plot_ellipse=FALSE to not plot ellipse ###
epidome_object_clinical_norm = normalize_epidome_object(epidome_object_clinical)
PCA_patient_colored = plot_PCA_epidome(epidome_object_clinical_norm,color_variable = "patient.ID",colors = c(),plot_ellipse = FALSE)
PCA_patient_colored
PCA_sample_site_colored = plot_PCA_epidome(epidome_object_clinical_norm,color_variable = "sample.site",colors = c("Red","Blue"),plot_ellipse = TRUE)
PCA_sample_site_colored








p = make_barplot_epidome(count_table, reorder = FALSE, normalize = TRUE)
p


p = plot_PCA_epidome(epidome_object_normalized,"patient.ID",c(),include_ellipse = FALSE)
p
plot_PCA_epidome(epidome_object_normalized,"sample.site",c())
#plot_PCA_epidome(epidome_object_normalized,"patient.sample.site",))


### Classify based on the two primer sets and combine into a single data frame ###
count_table = classify_epidome(epidome_filtered_samples,ST_amplicon_table)

### Make a barplot of the classified data, set reorder = TRUE to sort samples based on Bray Curtis dissimilarity, set normalize = TRUE to plot percent ###
p = make_barplot_epidome(count_table, reorder = FALSE, normalize = TRUE)
p

p = make_barplot_epidome(count_table, reorder = FALSE, normalize = FALSE)
p


epidome_object_clinical = prune_by_variable_epidome(epidome_object,"sample.type",c("Clinical"))

epidome_object_clinical_norm = normalize_epidome_object(epidome_object_clinical)


pca_patient = plot_PCA_epidome(epidome_object_clinical_norm,"patient.ID",c(),include_ellipse = FALSE)
p

pca_samplesite = plot_PCA_epidome(epidome_object_clinical_norm,"sample.site",c(),include_ellipse = TRUE)
p


count_table_clinical = classify_epidome(epidome_object_clinical,ST_amplicon_table)
p_clinical = make_barplot_epidome(count_table_clinical, reorder = FALSE, normalize = TRUE)
p


epidome_object_mock = prune_by_variable_epidome(epidome_ASV_combined,"sample.type",c("Mock community"))
count_table_mock = classify_epidome(epidome_object_mock,ST_amplicon_table)
p_mock = make_barplot_epidome(count_table_mock, reorder = FALSE, normalize = TRUE)
p


p1 = epidome_object_mock$p1_table[which(rowSums(epidome_object_mock$p1_table)>1000),]
p2 = epidome_object_mock$p2_table[which(rowSums(epidome_object_mock$p2_table)>1000),]
p1_s = epidome_object_mock$p1_seqs[which(rowSums(epidome_object_mock$p1_table)>1000)]
p2_s = epidome_object_mock$p1_seqs[which(rowSums(epidome_object_mock$p2_table)>1000)]




test = dist_comparison(p1_dist,factor(t3_norm$metadata$sample.site))
test = dist_comparison(p1_dist,factor(t3_norm$metadata$patient.ID))

p <- plot_ly(type="box",data=test_df,x=~Groups,y=~dist)
p

test_vec = as.vector(t3_norm$metadata$sample.type)

g1_vec = c()
g2_vec = c()

for (i in 1:length(test_vec)) {
  g1_vec = c(g1_vec,rep(test_vec[i],(length(test_vec)-i)))
  g2_vec = c(g2_vec,test_vec[(i+1):length(test_vec)])
}
g2_vec = g2_vec[1:length(g1_vec)]
g1_vec[1:10]
g2_vec[1:10]

length(g1_vec)
length(g2_vec)

type_comb = paste0(g1_vec,' - ',g2_vec)

p = plot_ly(type='box',x=type_comb,y=as.vector(p1_dist))

p1_dist = dist(t(t3_norm$p1_table))


p2_dist = dist(t(t3_norm$p2_table))







color_variable = "patient.ID"




d1_1 = cbind(epidome_object$p1_seqs,epidome_object$p1_table[,13])
d1_2 = cbind(epidome_object$p2_seqs,epidome_object$p2_table[,13])
d2_1 = cbind(t2$p1_seqs,t2$p1_table[,13])
d2_2 = cbind(t2$p2_seqs,t2$p2_table[,13])


plot_PCA_epidome(t3_norm,"patient.ID",c())
plot_PCA_epidome(t3_norm,"sample.site",c())

d_top10 = t3_norm$p1_table[1:10,]
d_top10_melt = melt(d_top10)


ggplot(data=d_top10_melt) + geom_bar(aes(x=X2,y=value,fill=as.character(X1)), stat="identity")







sample_names = colnames(epi01_table_full)[3:ncol(epi01_table_full)]
patient_IDs = unlist(lapply(sample_names, function(x) strsplit(x,'_')[[1]][1]))
patient_IDs[which(patient_IDs=="E")] = "Even mock"
patient_IDs[which(patient_IDs=="S")] = "Staggered mock"
patient_IDs[which(patient_IDs=="nk")] = "Negative control"
sample_site = unlist(lapply(sample_names, function(x) strsplit(x,'_')[[1]][2]))
sample_site[which(patient_IDs=="Even mock")] = "Even mock"
sample_site[which(patient_IDs=="Staggered mock")] = "Staggered mock"
sample_site[which(patient_IDs=="Negative control")] = "Negative control"


metadata_table = data.frame("sample.ID"=sample_names,"patient.ID"=patient_IDs,"sample.site"=sample_site,row.names = sample_names)

metadata_table$sample.type = "Clinical"
metadata_table$sample.type[which(metadata_table$patient.ID %in% c("Staggered mock","Even mock"))] = "Mock community"
metadata_table$sample.type[which(metadata_table$patient.ID %in% c("Negative control"))] = "Negative control"

#write.table(metadata_table,file = "example_data/sample_metadata.txt",sep = "\t")

#writeLines(sample_names,"example_data/sample_names.txt")

meta_test = read.table("example_data/sample_metadata.txt",sep = "\t")

t2 = setup_epidome_object(epi01_table_full,epi02_table_full,metadata_table = metadata_table)


epi_01_seq_names = epi01_table_full$Seq_number
epi_02_seq_names = epi02_table_full$Seq_number


epi01_table = epi01_table_full[,3:ncol(epi01_table_full)]
epi02_table = epi02_table_full[,3:ncol(epi02_table_full)]

epi01_table = epi01_table_full[which(!is.na(epi01_table_full$Seq_number)),3:ncol(epi01_table_full)]
epi02_table = epi02_table_full[which(!is.na(epi02_table_full$Seq_number)),3:ncol(epi02_table_full)]

epi01_table_full = epi01_table_full[,3:ncol(epi01_table_full)]
epi02_table_full = epi02_table_full[,3:ncol(epi02_table_full)]


which(colnames(epi01_table) %in% colnames(epi02_table))
which(!colnames(epi01_table) %in% colnames(epi02_table))
which(!colnames(epi02_table) %in% colnames(epi01_table))

### Count thresholds for inclusion ###
plot(sort(colSums(epi01_table)))
sort(colSums(epi01_table))
plot(sort(colSums(epi01_table))[1:10])
epi01_threshold = 863

plot(sort(colSums(epi02_table)))
sort(colSums(epi02_table))
plot(sort(colSums(epi02_table))[1:10])
epi02_threshold = 547

include_index = which(colSums(epi01_table>=epi01_threshold) & colSums(epi02_table>=epi02_threshold))
epi_01_pruned = epi01_table[,include_index]
epi_02_pruned = epi02_table[,include_index]

epi_01_percent = apply(epi_01_pruned,2,function(x) x/sum(x)*100)
epi_02_percent = apply(epi_02_pruned,2,function(x) x/sum(x)*100)

### PCA ###


pca = prcomp(rbind(epi_01_percent,epi_02_percent))

ggplot(data=as.data.frame(pca$x), aes(x=PC1,y=PC2,color=m_subset$patient_ID)) + scale_color_manual(values = colorRamps::primary.colors(19)[c(1,2,3,4,14,6,7,8,9,12,11,10,13,5,15,16,17,18,19)]) + geom_point(size = 3) + theme(legend.position = "none") + xlab(labels[1]) + ylab(labels[2])
ggplot(data=as.data.frame(pca$x), aes(x=PC3,y=PC4)) + geom_point(size = 1)




sort(colSums(epi02_table))

plot(colSums(epi01_table),colSums(epi02_table))
cor.test(colSums(epi01_table),colSums(epi02_table))


plot(sort(colSums(epi01_table_full)))
sort(colSums(epi01_table_full))
plot(sort(colSums(epi02_table_full)))
sort(colSums(epi02_table_full))

ggplot(mapping=aes(x=c(colSums(epi01_table_full),colSums(epi01_table)),c(colSums(epi02_table_full),colSums(epi02_table)))) + geom_point(aes(colour = c(rep("All ASVs",ncol(epi01_table)),rep("ASVs matching sequence in database",ncol(epi01_table))))) + 
  xlab(label = "epi01 sequence counts") + ylab(label = "epi02 sequence counts") + labs(color = "") + ggtitle("Number of epidome sequnces classified by dada2")
cor.test(colSums(epi01_table_full),colSums(epi02_table_full))

ST_amplicon_table$epi01 = paste0("seq",ST_amplicon_table$epi01_ASV)
ST_amplicon_table$epi02 = paste0("seq",ST_amplicon_table$epi02_ASV)

assign_ST = function(epi01_seq,epi02_seqs,ST_amplicon_table) {
  sub_ST_table = ST_amplicon_table[which(ST_amplicon_table$epi01==epi01_seq),]
  print(sub_ST_table)
}

assign_ST("seq40",epi02_table_full$Seq_number,ST_amplicon_table)



classify_g216_start <- function(yycH_IDs, g216_IDs, yycH_list, g216_list, mlst, percent_threshold) {
  return_IDs = c()
  return_STs = c()
  return_percent = c()
  missing_yycH = which(!yycH_IDs %in% g216_IDs)
  missing_g216 = which(!g216_IDs %in% yycH_IDs)
  if (length(missing_yycH)>0 | length(missing_g216)>0) {
    print('mismatching IDs')
    print(paste0('yycH:'))
    print(yycH_IDs[missing_yycH])
    print(paste0('g216:'))
    print(g216_IDs[missing_g216])
  }
  print(missing_yycH)
  print(yycH_IDs)
  IDs = yycH_IDs[-missing_yycH]
  IDs = yycH_IDs
  for (ID in IDs) {
    print(ID)
    yycH_match_vec = yycH_hits$match_list[[ID]]
    g216_match_vec = g216_hits$match_list[[ID]]
    yycH_percent_vec = yycH_match_vec/sum(yycH_match_vec)*100
    g216_percent_vec = g216_match_vec/sum(g216_match_vec)*100
    yycH_percent_vec = yycH_percent_vec[which(yycH_percent_vec>percent_threshold)]
    g216_percent_vec = g216_percent_vec[which(g216_percent_vec>percent_threshold)]
    ST_vec = c()
    ST_vec_names = c()
    for (i in 1:length(g216_percent_vec)) {
      g216_group = names(g216_percent_vec)[i]
      g216_percent = g216_percent_vec[i]
      mlst_sub = mlst[which(mlst$g216==g216_group),]
      g216_STs = as.vector(mlst_sub$ST)
      g216_STs_uniq = unique(g216_STs)
      g216_ST_groups = as.vector(mlst_sub$group_ST)
      g216_ST_groups_uniq = unique(g216_ST_groups)
      if (length(g216_STs_uniq)==1) {
        if (length(g216_ST_groups_uniq)==1) {
          if (g216_STs_uniq[1] == '-') {
            print(mlst_sub)
            ST_to_add = g216_ST_groups_uniq[1]
            print(ST_to_add)
          } else {
            ST_to_add = g216_STs_uniq[1]
          }
        } else {
          ST_to_add = g216_STs_uniq[1]
        }
        ST_vec_names = c(ST_vec_names,ST_to_add)
        ST_vec = c(ST_vec,g216_percent_vec[i])
      } else {
        yycH_groups = as.vector(mlst_sub$yycH)
        within_range_yycH_groups = names(yycH_percent_vec)[which(yycH_percent_vec/g216_percent<8 & yycH_percent_vec/g216_percent>0.125)]
        groups_to_check = yycH_groups[which(yycH_groups %in% within_range_yycH_groups)]
        #print(paste0('groups to check ',groups_to_check))
        if (length(groups_to_check)==1) {
          ST_name = as.vector(mlst_sub$ST[which(mlst_sub$yycH %in% groups_to_check)])
          ST_vec_names = c(ST_vec_names,ST_name)
          ST_vec = c(ST_vec,g216_percent)
          #print(paste0('ST_name ',ST_name))
          #print(g216_percent)
        } else if (length(groups_to_check)>1) {
          ST_vec = c(ST_vec,g216_percent)
          sorted_table = sort(table(mlst_table$group_ST[which(mlst_table$group_g216==g216_group & mlst_table$group_yycH %in% groups_to_check)]),decreasing = T)
          ST_name = names(sorted_table)[1]
          ST_vec_names = c(ST_vec_names,ST_name)
        } else if (length(groups_to_check)==0) {
          ST_vec = c(ST_vec,g216_percent)
          ST_vec_names = c(ST_vec_names,'Unclassified')
        } else {
          print('HEY HERE')
        }
      }
    }
    ST_vec = ST_vec/sum(ST_vec)*100
    return_STs = c(return_STs,ST_vec_names)
    return_percent = c(return_percent,ST_vec)
    names(ST_vec) = ST_vec_names
    return_IDs = c(return_IDs,rep(ID,length(ST_vec_names)))
  }
  print(length(return_IDs))
  print(length(return_STs))
  print(length(return_percent))
  return_df = data.frame('ID'=return_IDs,'ST'=return_STs,'Percent'=return_percent)
  return(return_df)
  
}
















get_dominant_ST_group <- function(group,mlst,group_vec) {
  sub_mlst = mlst[which(group_vec==group),]
  sub_STs = mlst$group_ST[which(group_vec==group)]
  tbl = table(sub_STs)
  ordered_table = sort(tbl,decreasing=T)
  top_ST = rownames(ordered_table)[1]
  return(top_ST)
}

get_dominant_ST <- function(group,mlst,group_vec) {
  sub_mlst = mlst[which(group_vec==group),]
  sub_STs = mlst$ST[which(group_vec==group)]
  tbl = table(sub_STs)
  ordered_table = sort(tbl,decreasing=T)
  top_ST = rownames(ordered_table)[1]
  return(top_ST)
}

uniq_yycH_groups = sort(unique(mlst_table$group_yycH))
uniq_g216_groups = sort(unique(mlst_table$group_g216))

yycH_uniq_ST_group = unlist(lapply(uniq_yycH_groups, function(x) get_dominant_ST_group(x,mlst_table,mlst_table$group_yycH)))
g216_uniq_ST_group = unlist(lapply(uniq_g216_groups, function(x) get_dominant_ST_group(x,mlst_table,mlst_table$group_g216)))

yycH_uniq_ST = unlist(lapply(uniq_yycH_groups, function(x) get_dominant_ST(x,mlst_table,mlst_table$group_yycH)))
g216_uniq_ST = unlist(lapply(uniq_g216_groups, function(x) get_dominant_ST(x,mlst_table,mlst_table$group_g216)))

g216_uniq_table = data.frame("group"=uniq_g216_groups,"ST"=g216_uniq_ST,"group_ST"=g216_uniq_ST_group)
yycH_uniq_table = data.frame("group"=uniq_yycH_groups,"ST"=yycH_uniq_ST,"group_ST"=yycH_uniq_ST_group)



uniq_comb_groups = sort(unique(mlst_table$amp_comb))
comb_uniq_ST = unlist(lapply(uniq_comb_groups, function(x) get_dominant_ST(x,mlst_table,mlst_table$amp_comb)))
comb_uniq_ST_group = unlist(lapply(uniq_comb_groups, function(x) get_dominant_ST_group(x,mlst_table,mlst_table$amp_comb)))

comb_uniq_table = data.frame("group"=uniq_comb_groups,"ST"=comb_uniq_ST,"group_ST"=comb_uniq_ST_group)

comb_uniq_table = comb_uniq_table[-grep("NA",comb_uniq_table$group),]

comb_uniq_table$g216 = unlist(lapply(as.vector(comb_uniq_table$group), function(x) strsplit(x,"_")[[1]][1]))
comb_uniq_table$yycH = unlist(lapply(as.vector(comb_uniq_table$group), function(x) strsplit(x,"_")[[1]][2]))






