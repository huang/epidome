setwd("/Volumes/data/MPV/projects/git.repositories/epidome/")


mlst_table = read.table("DB/epidome_mlst_table.txt", sep = "\t")

ST_amplicon_table = read.table("DB/epidome_ST_amplicon_frequences.txt",sep = "\t")


epi01_table_full = read.table("example_data/epi01_count_table.csv",sep = ";",header=TRUE,row.names=1)
epi02_table_full = read.table("example_data/epi02_count_table.csv",sep = ";",header=TRUE,row.names=1)

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






