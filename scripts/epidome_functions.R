
primer1_table = epi01_table_full
primer2_table = epi02_table_full


setup_epidome_object <- function(primer1_table,primer2_table,metadata_table) {
  primer1_counts = primer1_table[,3:ncol(primer1_table)]
  primer2_counts = primer2_table[,3:ncol(primer2_table)]
  primer1_all_sample_names = colnames(primer1_counts)
  primer2_all_sample_names = colnames(primer2_counts)
  samples_with_both_primers = primer1_all_sample_names[which(primer1_all_sample_names %in% primer2_all_sample_names)]
  samples_missing_primer1_data = primer2_all_sample_names[which(!primer2_all_sample_names %in% primer1_all_sample_names)]
  samples_missing_primer2_data = primer1_all_sample_names[which(!primer1_all_sample_names %in% primer2_all_sample_names)]
  primer1_seqs = primer1_table$Seq_number
  primer2_seqs = primer2_table$Seq_number
  primer1_counts = primer1_table[,which(colnames(primer1_table) %in% samples_with_both_primers)]
  primer2_counts = primer2_table[,which(colnames(primer2_table) %in% samples_with_both_primers)]
  if (!missing(metadata_table)) {
    metadata_names = rownames(metadata_table)
    samples_missing_metadata = samples_with_both_primers[which(!samples_with_both_primers %in% metadata_names)]
    samples_missing_primer_data = metadata_names[which(!metadata_names %in% samples_with_both_primers)]
    include_names = metadata_names[which(metadata_names %in% samples_with_both_primers)]
    metadata_table = metadata_table[match(metadata_names,include_names),]
    epi1_table = primer1_counts[,match(colnames(primer1_counts),include_names)]
    epi2_table = primer2_counts[,match(colnames(primer2_counts),include_names)]
    epidome_object = list('p1_seqs'=primer1_seqs,'p1_table'=epi1_table,'p2_seqs'=primer2_seqs,'p2_table'=epi2_table,'metadata'=metadata_table,'sample_names'=include_names,'meta_variables'=colnames(metadata_table))
    print(paste0("Metadata loaded with ",length(metadata_names)," samples and ",ncol(metadata_table)," variables"))
    print(paste0(length(samples_missing_metadata)," samples are found in both primer sets but not in metadata: ",paste0(samples_missing_metadata,collapse = " ")))
    print(paste0(length(samples_missing_primer_data)," samples are found in metadata but is missing from one or both primer sets: ",paste0(samples_missing_primer_data,collapse = " ")))
    print(paste0(length(include_names)," samples are found in metadata and both tables and are included in epidome object"))
    
  } else {
    epi1_table = primer1_counts[,match(colnames(primer1_counts),samples_with_both_primers)]
    epi2_table = primer1_counts[,match(colnames(primer2_counts),samples_with_both_primers)]
    epidome_object = list('p1_seqs'=primer1_seqs,'p1_table'=epi1_table,'p2_seqs'=primer2_seqs,'p2_table'=epi2_table,'sample_names'=samples_with_both_primers)
    print(paste0("No metadata loaded"))
    print(paste0(length(samples_missing_primer2_data)," samples are found in p1 table but not in p2 table: ",paste0(samples_missing_primer2_data,collapse = " ")))
    print(paste0(length(samples_missing_primer1_data)," samples are found in p2 table but not in p1 table: ",paste0(samples_missing_primer2_data,collapse = " ")))
    print(paste0(length(samples_with_both_primers)," samples are found in both tables and are included in epidome object"))
  }
  return(epidome_object)
}


test = setup_epidome_object(epi01_table_full,epi02_table_full)

epidome_object = setup_epidome_object(epi01_table_full,epi02_table_full,metadata_table = metadata_table)


compare_primer_output <- function(epidome_object) {
  p1_counts = colSums(epidome_object$p1_table)
  p2_counts = colSums(epidome_object$p2_table)
  cor = cor.test(p1_counts,p2_counts)
  ggplot(data.frame(p1_counts,p2_counts),aes(x=p1_counts,y=p2_counts)) + geom_point()
  
}


filter_lowcount_epidome = function(epidome_object,p1_threshold,p2_threshold) {
  original_sample_names = epidome_object$sample_names
  include_index = which(colSums(epidome_object$p1_table) >= p1_threshold & colSums(epidome_object$p2_table) >= p2_threshold)
  epidome_object$p1_table = epidome_object$p1_table[,include_index]
  epidome_object$p2_table = epidome_object$p2_table[,include_index]
  epidome_object$sample_names = epidome_object$sample_names[include_index]
  epidome_object$metadata = epidome_object$metadata[include_index,]
  filtered_samples = original_sample_names[-include_index]
  print(paste0(length(filtered_samples)," low count samples removed from data: ",paste0(filtered_samples,collapse = " ")))
  return(epidome_object)
}

t3 = filter_lowcount_epidome(epidome_object,500,500)


normalize_epidome_object = function(epidome_object) {
  epidome_object$p1_table = apply(epidome_object$p1_table,2,function(x) x/sum(x)*100)
  epidome_object$p2_table = apply(epidome_object$p2_table,2,function(x) x/sum(x)*100)
  return(epidome_object)
}

t3_norm = normalize_epidome_object(t3)

setup_colors = function(factor_levels,colors) {
  group_count = length(factor_levels)
  if (!group_count == length(colors)) {
    if (group_count<=9) {
      colors = RColorBrewer::brewer.pal(length(factor_levels),"Set1")
    } else if (group_count<=12) {
      colors = RColorBrewer::brewer.pal(12,name="Set3")[1:group_count]
    } else {
      colors  = grDevices::rainbow(group_count)
    }
  }
}

color_variable = "patient.ID"
plot_PCA_epidome = function(epidome_object,color_variable,colors) {
  m = epidome_object$metadata
  color_variable_factor = m[,which(epidome_object$meta_variables==color_variable)]
  data_combined = rbind(epidome_object$p1_table,epidome_object$p2_table)
  pca = prcomp(t(data_combined))
  plot_df = data.frame(pca$x)
  color_vector = setup_colors(levels(color_variable_factor),colors)
  ggplot(as.data.frame(pca$x),aes(x=PC1,y=PC2,color = color_variable_factor)) + labs(color = color_variable) + geom_point(size=1, alpha=1)+ stat_ellipse(level=0.75) + scale_colour_manual(values = color_vector)
}

plot_PCA_epidome(t3_norm,"patient.ID",c())
plot_PCA_epidome(t3_norm,"sample.site",c())




p1_dist = dist(t(t3_norm$p1_table))


p2_dist = dist(t(t3_norm$p2_table))
  


