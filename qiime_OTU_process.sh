#!/bin/bash
MANIFEST = "$1"
METADATA = "$2"




##############################################

mkdir importsequences

# importar secuencias
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path $MANIFEST \
--output-path importsequences/paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2

# visualizar secuencias
qiime demux summarize \
--i-data importsequences/paired-end-demux.qza \
--o-visualization importsequences/paired-end.demux.qzv

#####################################
mkdir dada2_result_OTUs

# denoising
qiime dada2 denoise-paired \
--i-demultiplexed-seqs importsequences/paired-end-demux.qza \
--p-trunc-len-f 0 \
--p-trunc-len-r 0 \
--p-n-threads 0 \ #max possible
--o-representative-sequences dada2_result_OTUs/rep_seqs_dada2.qza \
--o-table dada2_result_OTUs/feature_table_dada2.qza \
--o-denoising-stats dada2_result_OTUs/stats_dada2.qza

# se clusterizan las secuencias en un 95%
qiime vsearch cluster-features-de-novo \
--i-table dada2_result_OTUs/feature_table_dada2.qza \
--i-sequences dada2_result_OTUs/rep_seqs_dada2.qza \
--p-perc-identity 0.95 \
--o-clustered-table dada2_result_OTUs/feature_table_95.qza \
--o-clustered-sequences dada2_result_OTUs/rep_seqs_95.qza 

# tabular secuencias representativas
qiime metadata tabulate \
--m-input-file dada2_result_OTUs/rep_seqs_95.qza \
--o-visualization dada2_result_OTUs/rep_seqs_95.qzv

# tabular estadísticas
qiime metadata tabulate 
 --m-input-file dada2_result_OTUs/stats-dada2.qza 
 --o-visualization dada2_result_OTUs/stats-dada2.qzv

# tabular feature table
qiime metadata tabulate \
--m-input-file dada2_result_OTUs/feature_table_95.qza \
--o-visualization dada2_result_OTUs/feature_table_95.qzv

# tabular feature table con metadatos
qiime feature-table summarize \
--i-table dada2_result_OTUs/feature_table_95.qza \
--o-visualization dada2_result_OTUs/feature_table_95_summarized.qzv \
--m-sample-metadata-file $METADATA

# tabular feature con identificador al mapeado
qiime feature-table tabulate-seqs \
--i-data dada2_result_OTUs/feature_table_95.qza \
--o-visualization dada2_result_OTUs/feature_table_tabulated_seqs_95.qzv

############################################
mkdir phylogeny_data_OTUs

# alinear al arbol
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences dada2_result_OTUs/rep_seqs_95.qza \
--o-alignment phylogeny_data_OTUs/aligned_rep_seqs.qza \
--o-masked-alignment phylogeny_data_OTUs/masked_aligned_rep_seqs.qza \
--o-tree phylogeny_data_OTUs/unrooted_tree_95.qza \
--o-rooted-tree phylogeny_data_OTUs/rooted_tree_95.qza

###################################################
mkdir diversity_data_OTUs

qiime diversity core-metrics-phylogenetic \
--i-phylogeny phylogeny_data_OTUs/rooted_tree_95.qza \
--i-table dada2_result_OTUs/feature_table_dada2.qza \
--p-sampling-depth [] \ #DEFINIR SAMPLE DEPTH
--m-metadata-file $METADATA \
--output-dir diversity_data_OTUs

# diversidad alfa con faith
qiime diversity alpha-group-significance \
--i-alpha-diversity diversity_data_OTUs/faith_pd_vector.qza \
--m-metadata-file $METADATA \
--o-visualization diversity_data_OTUs/faith_pd_summarized.qzv

# diversidad alfa con pielou
qiime diversity alpha-group-significan 
qiime diversity alpha-group-significance \
--i-alpha-diversity diversity_data_OTUs/shannon_vector.qza \
--m-metadata-file $METADATA \
--o-visualization diversity_data_OTUs/shannon_vector_summarized.qzv

# diversidad beta por blocking_primers
qiime diversity beta-group-significance \
--i-distance-matrix diversity_data_OTUs/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file $METADATA \
--m-metadata-column block \
--o-visualization diversity_data_OTUs/unweighted_unifrac_species_significance.qzv \
--p-pairwise

####################################################
mkdir alpha_rarefaction_OTUS

qiime diversity alpha-rarefaction \
--i-table dada2_result_OTUs/featuretable_dada2.qza \
--i-phylogeny phylogeny_data_OTUs/rooted_tree_95.qza \
--p-max-depth [] \
--m-metadata-file $METADATA\
--p-iterations [] \
--o-visualization alpha_rarefaction/alpha_rarefaction.qzv

#######################################################
mkdir taxonomy_otus

wget -O silva_rep_seqs.qza https://data.qiime2.org/2021.4/common/silva-138-99-seqs.qza
get -O silva_tax.qza https://data.qiime2.org/2021.4/common/silva-138-99-tax.qza

qiime feature-classifier classify-consensus-vsearch \
--i-query dada2_result_OTUs/rep_seqs_95.qza \
--i-reference-reads silva_rep_seqs.qza \
--i-reference-taxonomy silva_tax.qza \
--p-perc-identity 0.95 \
--o-classification taxonomy_otus/taxonomy.qza

qiime metadata tabulate \
--m-input-file taxonomy_otus/taxonomy.qza \
--o-visualization taxonomy_otus/taxonomy.qzv

qiime taxa barplot \
--i-table dada2_result_OTUs/featuretable_dada2.qza \
--i-taxonomy taxonomy_otus/taxonomy.qza \
--m-metadata-file $METADATA \
--o-visualization taxonomy_otus/taxa-bar-plots.qzv




#########################################MODIFICAR AQUI##################
qiime feature-table filter-samples \
--i-table ../04-qiime/01-dada2/table-dada2.qza \
--m-metadata-file /home/scratch/irene/project_Natalia16S-bloq/ANALYSIS/samples-metadata.tsv \
--p-where "[block]='blocked'" \
--o-filtered-table ../04-qiime/07-differencial-abundance/montana-table.qza

qiime feature-table filter-samples \
--i-table ../04-qiime/01-dada2/table-dada2.qza \
--m-metadata-file /home/scratch/irene/project_Natalia16S-bloq/ANALYSIS/samples-metadata.tsv \
--p-where "[block]='not blocked'" \
--o-filtered-table ../04-qiime/07-differencial-abundance/montana-table.qza

