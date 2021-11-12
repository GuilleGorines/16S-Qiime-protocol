Welcome to this 16S analysis protocol using Qiime2, version 2021.2
##Environment settings



First step, before anything else, is setting a correct, easy-to-follow structure inside the work directory: 

```
mkdir ANALYSIS BIN REFERENCES RESULTS RAW TRIMMED
```
`BIN` will contain all scripts that will be used for this work. Right now its empty, so filling it with the appropriate files will be needed.

```
wget -O BIN/sample_catalog.py https://raw.githubusercontent.com/GuilleGorines/16S-Qiime-protocol/main/Sample_catalog.py
wget -O BIN/fastqc_directories.sh  https://raw.githubusercontent.com/GuilleGorines/16S-Qiime-protocol/main/fastqc_directories.sh
wget -O BIN/process_table.py https://raw.githubusercontent.com/GuilleGorines/16S-Qiime-protocol/main/process_table.py
```
Once that is done, place your reads inside the `5-RAW` directory, and execute


First step for qiime2 is to import the sequence that will be used. A visual file (this is, with  _.qzv_ extension) will be created as well.

```
mkdir importsequences

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path $MANIFEST \
--output-path importsequences/paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
--i-data importsequences/paired-end-demux.qza \
--o-visualization importsequences/paired-end-demux.qzv

```

```
mkdir dada2_result

qiime dada2 denoise-paired \
--i-demultiplexed-seqs importsequences/paired-end-demux.qza \
--p-trunc-len-f 0 \
--p-trunc-len-r 0 \
--p-n-threads 0 \
--o-representative-sequences dada2_result/rep_seqs_dada2.qza \
--o-table dada2_result/feature_table_dada2.qza \
--o-denoising-stats dada2_result/stats_dada2.qza

qiime metadata tabulate \
--m-input-file dada2_result/rep_seqs_dada2.qza \
--o-visualization dada2_result/rep_seqs_dada2.qzv

qiime metadata tabulate \
 --m-input-file dada2_result/stats_dada2.qza \
 --o-visualization dada2_result/stats_dada2.qzv

qiime metadata tabulate \
--m-input-file dada2_result/feature_table_dada2.qza \
--o-visualization dada2_result/feature_table_dada2.qzv

qiime feature-table summarize \
--i-table dada2_result/feature_table_dada2.qza \
--o-visualization dada2_result/feature_table_dada2_summarized.qzv \
--m-sample-metadata-file $METADATA

qiime feature-table tabulate-seqs \
--i-data dada2_result/rep_seqs_dada2.qza \
--o-visualization dada2_result/feature_table_tabulated_seqs.qzv
```

Next step is identifying the sequences. This task will be performed by using the official Naive Bayes Classifier from the Qiime2 website. 
```
mkdir taxonomy
wget -O ../REFERENCES/bayes-classifier.qza https://data.qiime2.org/2021.2/common/silva-138-99-nb-classifier.qza

qiime feature-classifier classify-sklearn \
--i-classifier ../REFERENCES/bayes-classifier.qza \
--i-reads dada2_result/rep_seqs_dada2.qza \
--o-classification taxonomy/taxonomy.qza

qiime metadata tabulate \
--m-input-file taxonomy/taxonomy.qza \
--o-visualization taxonomy/taxonomy.qzv
```

With the identified sequences, making decisions becomes easier. For example, as mitochondria and chloroplasts from the host might be present as they all present the 16S gene. Here, we will perform the analysis both with and without the mitochondria and chloroplast.

## Without mitochondria and chloroplast removal

```
mkdir analysis_with_mitochondria_chloroplast
mkdir analysis_with_mitochondria_chloroplast/identification

qiime taxa barplot \
--i-table dada2_result/feature_table_dada2.qza \
--i-taxonomy taxonomy/taxonomy.qza \
--m-metadata-file $METADATA \
--o-visualization analysis_with_mitochondria_chloroplast/identification/taxa-bar-plots.qzv
```

In order to obtain the (absolute) number of occurrences of each feature, a little juggling is required (collapsing the feature table to the seventh (species) level exporting the feature data to a biom table, converting the biom table, and parsing it with an _ad hoc_ Python script)

```
qiime taxa collapse \
  --i-table dada2_result/feature_table_dada2.qza \
  --i-taxonomy taxonomy/taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table analysis_with_mitochondria_chloroplast/identification/table-no-mitochondria-no-chloroplast_spp_lvl.qza

qiime tools export \
  --input-path analysis_with_mitochondria_chloroplast/identification/table-no-mitochondria-no-chloroplast_spp_lvl.qza \
  --output-path analysis_with_mitochondria_chloroplast/identification

biom convert \
  --input-fp analysis_with_mitochondria_chloroplast/identification/feature-table.biom \
  --output-fp analysis_with_mitochondria_chloroplast/identification/transposed_table_no_mitochondria_no_cloroplast.tsv \
  --to-tsv

python process_table.py \
   analysis_with_mitochondria_chloroplast/identification/transposed_table_no_mitochondria_no_cloroplast.tsv \
   analysis_with_mitochondria_chloroplast/identification/absolute_numbers_taxonomy_no_mitochondria_no_cloroplast.tsv 

```
Once this is done, we can make both of the trees for the sequences (rooted and unrooted)

```
mkdir analysis_with_mitochondria_chloroplast/phylogeny_data

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences dada2_result/rep_seqs_dada2.qza \
--o-alignment analysis_with_mitochondria_chloroplast/phylogeny_data/aligned_rep_seqs.qza \
--o-masked-alignment analysis_with_mitochondria_chloroplast/phylogeny_data/masked_aligned_rep_seqs.qza \
--o-tree analysis_with_mitochondria_chloroplast/phylogeny_data/unrooted_tree.qza \
--o-rooted-tree analysis_with_mitochondria_chloroplast/phylogeny_data/rooted_tree.qza
```

The next part of the analysis, the diversity data, is difficult to automate, as you need the sampling depth from `dada2_result/feature_table_dada2_summarized.qzv`.

```
mkdir analysis_with_mitochondria_chloroplast/diversity_data

qiime diversity core-metrics-phylogenetic \
--i-phylogeny phylogeny_data/rooted_tree.qza \
--i-table dada2_result/feature_table_dada2.qza \
--p-sampling-depth $SAMPLING_DEPTH \
--m-metadata-file $METADATA \
--output-dir diversity_data

qiime diversity alpha-group-significance \
--i-alpha-diversity diversity_data/faith_pd_vector.qza \
--m-metadata-file $METADATA \
--o-visualization diversity_data/faith_pd_summarized.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity diversity_data/shannon_vector.qza \
--m-metadata-file $METADATA \
--o-visualization diversity_data/shannon_vector_summarized.qzv

qiime diversity beta-group-significance \
--i-distance-matrix diversity_data/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file $METADATA \
--m-metadata-column group \
--o-visualization diversity_data/unweighted_unifrac_species_significance.qzv \
--p-pairwise

```


```
mkdir alpha_rarefaction
# NOTA: max_depth es el valor máximo del eje X
# posibilidad: usar el numero maximo de features encontrado en la summarized table

qiime diversity alpha-rarefaction \
--i-table dada2_result/feature_table_dada2.qza \
--i-phylogeny phylogeny_data/rooted_tree.qza \
--p-max-depth 1009 \
--m-metadata-file $METADATA \
--o-visualization alpha_rarefaction/alpha_rarefaction.qzv


```



## Removing mitochondria and chloroplast

```
mkdir analysis_with_mitochondria_chloroplast
mkdir analysis_with_mitochondria_chloroplast/identification




```


qiime taxa filter-table \
  --i-table feature_table_dada2.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

qiime metadata tabulate \
  --m-input-file table-no-mitochondria-no-chloroplast.qza \
  --o-visualization table-no-mitochondria-no-chloroplast.qzv
  
```
qiime taxa collapse \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table table-no-mitochondria-no-chloroplast_spp_lvl.qza

qiime tools export \
  --input-path table-no-mitochondria-no-chloroplast_spp_lvl.qza \
  --output-path analysis_with_mitochondria_chloroplast/identification/identification.biom

biom convert \
--input-fp analysis_with_mitochondria_chloroplast/identification/feature-table.biom \
--output-fp biom.tsv \
--to-tsv

python process_table.py biom.tsv output.tsv 
```
