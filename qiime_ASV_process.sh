#!/bin/bash
MANIFEST = "$1"
METADATA = "$2"
###########################################



##############################################
mkdir importsequences

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path $MANIFEST \
--output-path importsequences/paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
--i-data importsequences/paired-end-demux.qza \
--o-visualization importsequences/paired-end-demux.qzv


# denoising
#####################################
mkdir dada2_result

qiime dada2 denoise-paired \
--i-demultiplexed-seqs importsequences/paired-end-demux.qza \
--p-trunc-len-f 0 \
--p-trunc-len-r 0 \
--p-n-threads 0 \
--o-representative-sequences dada2_result/rep_seqs_dada2.qza \
--o-table dada2_result/feature_table_dada2.qza \
--o-denoising-stats dada2_result/stats_dada2.qza

# tabular secuencias representativas
qiime metadata tabulate \
--m-input-file dada2_result/rep_seqs_dada2.qza \
--o-visualization dada2_result/rep_seqs_dada2.qzv

# tabular estadísticas
qiime metadata tabulate \
 --m-input-file dada2_result/stats_dada2.qza \
 --o-visualization dada2_result/stats_dada2.qzv

# tabular feature table
qiime metadata tabulate \
--m-input-file dada2_result/feature_table_dada2.qza \
--o-visualization dada2_result/feature_table_dada2.qzv

# tabular feature table con metadatos
# NOTA: de aquí se saca la sampling depth
qiime feature-table summarize \
--i-table dada2_result/feature_table_dada2.qza \
--o-visualization dada2_result/feature_table_dada2_summarized.qzv \
--m-sample-metadata-file $METADATA

# tabular feature con identificador al mapeado
qiime feature-table tabulate-seqs \
--i-data dada2_result/rep_seqs_dada2.qza \
--o-visualization dada2_result/feature_table_tabulated_seqs.qzv

############################################
mkdir phylogeny_data

# alinear al arbol
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences dada2_result/rep_seqs_dada2.qza \
--o-alignment phylogeny_data/aligned_rep_seqs.qza \
--o-masked-alignment phylogeny_data/masked_aligned_rep_seqs.qza \
--o-tree phylogeny_data/unrooted_tree.qza \
--o-rooted-tree phylogeny_data/rooted_tree.qza

###################################################
#NOTA: SAMPLING DEPTH SE SACA DE feature_table_dada2_summarized.qzv
qiime diversity core-metrics-phylogenetic \
--i-phylogeny phylogeny_data/rooted_tree.qza \
--i-table dada2_result/feature_table_dada2.qza \
--p-sampling-depth 504 \
--m-metadata-file $METADATA \
--output-dir diversity_data

# diversidad alfa con faith
# NOTA: NO SALDRÁ SI CADA GRUPO TIENE SOLO UNA MUESTRA
qiime diversity alpha-group-significance \
--i-alpha-diversity diversity_data/faith_pd_vector.qza \
--m-metadata-file $METADATA \
--o-visualization diversity_data/faith_pd_summarized.qzv

# diversidad alfa con pielou
# NOTA: NO SALDRÁ SI CADA GRUPO TIENE SOLO UNA MUESTRA
qiime diversity alpha-group-significance \
--i-alpha-diversity diversity_data/shannon_vector.qza \
--m-metadata-file $METADATA \
--o-visualization diversity_data/shannon_vector_summarized.qzv

# diversidad beta por especies
# NOTA: NO SALDRÁ SI CADA GRUPO TIENE SOLO UNA MUESTRA
# NOTA: CAMBIAR EL METADATA COLUMN SEGUN NECESIDADES
qiime diversity beta-group-significance \
--i-distance-matrix diversity_data/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file $METADATA \
--m-metadata-column group \
--o-visualization diversity_data/unweighted_unifrac_species_significance.qzv \
--p-pairwise

# beta diversidad y pcoa segun indice de bray-curtis
mkdir beta-diversity
qiime diversity-lib bray-curtis --i-table dada2_result/feature_table_dada2.qza --o-distance-matrix beta-diversity/bray_curtis_matrix.qza --p-n-jobs "auto"


# diversidad beta por blocking_primers
# NOTA: NO SALDRÁ SI CADA GRUPO TIENE SOLO UNA MUESTRA
# NOTA: CAMBIAR EL METADATA COLUMN SEGUN NECESIDADES
qiime diversity beta-group-significance \
--i-distance-matrix diversity_data/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file $METADATA \
--m-metadata-column group \
--o-visualization diversity_data/unweighted_unifrac_species_significance.qzv \
--p-pairwise

####################################################
mkdir alpha_rarefaction
# NOTA: max_depth es el valor máximo del eje X
# posibilidad: usar el numero maximo de features encontrado en la summarized table

qiime diversity alpha-rarefaction \
--i-table dada2_result/feature_table_dada2.qza \
--i-phylogeny phylogeny_data/rooted_tree.qza \
--p-max-depth 1009 \
--m-metadata-file $METADATA \
--o-visualization alpha_rarefaction/alpha_rarefaction.qzv

#######################################################
mkdir taxonomy

# source: https://docs.qiime2.org/2021.4/data-resources/
# Descarga del modelo

# V 2021.2
# Silva
wget -O taxonomy/bayes-classifier.qza https://data.qiime2.org/2021.2/common/silva-138-99-nb-classifier.qza

# V 2021.4
# Silva
wget -O taxonomy/bayes-classifier.qza https://data.qiime2.org/2021.4/common/silva-138-99-nb-classifier.qza


#######################################################
# Entrenar el modelo (alternativa)

wget -O taxonomy/ref_seqs_naive_bayes_training.qza https://data.qiime2.org/2021.4/common/silva-138-99-seqs.qza
wget -O taxonomy/ref_tax_naive_bayes_training.qza https://data.qiime2.org/2021.4/common/silva-138-99-tax.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads taxonomy/ref_seqs_naive_bayes_training.qza \
  --i-reference-taxonomy taxonomy/ref_tax_naive_bayes_training.qza \
  --o-classifier taxonomy/bayes-classifier.qza
  
#######################################################

# Utiliza el modelo para clasificar

qiime feature-classifier classify-sklearn \
--i-classifier taxonomy/bayes-classifier.qza \
--i-reads dada2_result/rep_seqs_dada2.qza \
--o-classification taxonomy/taxonomy.qza

qiime metadata tabulate \
--m-input-file taxonomy/taxonomy.qza \
--o-visualization taxonomy/taxonomy.qzv

qiime taxa barplot \
--i-table dada2_result/feature_table_dada2.qza \
--i-taxonomy taxonomy/taxonomy.qza \
--m-metadata-file $METADATA \
--o-visualization taxonomy/taxa-bar-plots.qzv

# 
qiime feature-table filter-samples \
--i-table ../04-qiime/01-dada2/table-dada2.qza \
--m-metadata-file $METADATA \
--p-where "[block]='blocked'" \
--o-filtered-table ../04-qiime/07-differencial-abundance/montana-table.qza

qiime feature-table filter-samples \
--i-table ../04-qiime/01-dada2/table-dada2.qza \
--m-metadata-file $METADATA \
--p-where "[block]='not blocked'" \
--o-filtered-table ../04-qiime/07-differencial-abundance/montana-table.qza

