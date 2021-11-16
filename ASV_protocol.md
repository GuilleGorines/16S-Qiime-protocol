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
Once that is done, place your reads inside the `RAW` directory, and execute

First step for Qiime2 is to import the sequence that will be used. A visual file (this is, with  _.qzv_ extension and displayable through `qiime2 tools view` or in https://view.qiime2.org/) will be created as well.

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
Once imported, Qiime2 can do quality control, and cluster the sequences by similarity to create a feature table, a stats file, and a representative sequences list. In this procedure, reads have already been through a quality control procedure, so this process will just group the sequences. 

Please note that this command is prepared for paired-end data, if single-end data was to be used, then `dada2 denoise-single` should be used instead.
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

Once this is done, we can make both of the trees (rooted and unrooted) for the representative sequences, after creating the directory for it.

```
mkdir phylogeny_data

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences dada2_result/rep_seqs_dada2.qza \
--o-alignment phylogeny_data/aligned_rep_seqs.qza \
--o-masked-alignment phylogeny_data/masked_aligned_rep_seqs.qza \
--o-tree phylogeny_data/unrooted_tree.qza \
--o-rooted-tree phylogeny_data/rooted_tree.qza
```

With the identified sequences, making decisions becomes easier. For instance, mitochondria and chloroplasts from the host might be present, as they all present the 16S gene. Their presence, however, might interfere with 

 Here, we will perform the analysis both with and without the mitochondria and chloroplast.

## Mitochondria and chloroplast NOT removed

First step is creating the directories where the files will be stored.

```
mkdir analysis_with_mit_chl
mkdir analysis_with_mit_chl/identification
```

After that, we generate the visual object containing the samples composition in a bar chart. 

```
qiime taxa barplot \
--i-table dada2_result/feature_table_dada2.qza \
--i-taxonomy taxonomy/taxonomy.qza \
--m-metadata-file $METADATA \
--o-visualization analysis_with_mit_chl/identification/taxa-bar-plots.qzv
```

In order to obtain the (absolute) number of occurrences of each feature, a little juggling is required (collapsing the feature table to the seventh (species) level exporting the feature data to a biom table, converting the biom table, and parsing it with an _ad hoc_ Python script)

```
qiime taxa collapse \
  --i-table dada2_result/feature_table_dada2.qza \
  --i-taxonomy taxonomy/taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table analysis_with_mit_chl/identification/table-with-mitochondria-chloroplast_spp_lvl.qza

qiime tools export \
  --input-path analysis_with_mit_chl/identification/table-with-mitochondria-chloroplast_spp_lvl.qza \
  --output-path analysis_with_mit_chl/identification

biom convert \
  --input-fp analysis_with_mit_chl/identification/feature-table.biom \
  --output-fp analysis_with_mit_chl/identification/transposed_table_with_mitochondria_cloroplast.tsv \
  --to-tsv

python process_table.py \
   analysis_with_mit_chl/identification/transposed_table_with_mitochondria_cloroplast.tsv \
   analysis_with_mit_chl/identification/absolute_numbers_taxonomy_with_mitochondria_cloroplast.tsv 

```

The next part of the analysis, the diversity data, is difficult to automate, as you need the sampling depth from `dada2_result/feature_table_dada2_summarized.qzv`. In this case, we will save it in the `SAMPLING_DEPTH` variable.

```
mkdir analysis_with_mit_chl/diversity_data

qiime diversity core-metrics-phylogenetic \
--i-phylogeny phylogeny_data/rooted_tree.qza \
--i-table dada2_result/feature_table_dada2.qza \
--p-sampling-depth $SAMPLING_DEPTH \
--m-metadata-file $METADATA \
--output-dir analysis_with_mit_chl/diversity_data

qiime diversity alpha-group-significance \
--i-alpha-diversity analysis_with_mit_chl/diversity_data/faith_pd_vector.qza \
--m-metadata-file $METADATA \
--o-visualization analysis_with_mit_chl/diversity_data/faith_pd_summarized.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity analysis_with_mit_chl/diversity_data/shannon_vector.qza \
--m-metadata-file $METADATA \
--o-visualization analysis_with_mit_chl/diversity_data/shannon_vector_summarized.qzv
```

Beta-significance coefficients migth be extracted as well between groups, but that depends on the `GROUP` variable, that must be the name of a column header in the metadata file. 
```
qiime diversity beta-group-significance \
--i-distance-matrix analysis_with_mit_chl/diversity_data/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file $METADATA \
--m-metadata-column $GROUP \
--o-visualization analysis_with_mit_chl/diversity_data/unweighted_unifrac_species_significance.qzv \
--p-pairwise
```

Lastly, the alpha rarefaction can be extracted. However, the max depth parameter, in the variable `MAX_DEPTH`, must be extracted from the `summarized_table.qzv`obtained previously.

```
mkdir analysis_with_mit_chl/alpha_rarefaction

qiime diversity alpha-rarefaction \
--i-table dada2_result/feature_table_dada2.qza \
--i-phylogeny phylogeny_data/rooted_tree.qza \
--p-max-depth $MAX_DEPTH \
--m-metadata-file $METADATA \
--o-visualization analysis_with_mit_chl/alpha_rarefaction/alpha_rarefaction.qzv
```

## Removing mitochondria and chloroplast

Once identified, mitochondria and chloroplast can be extracted, so they dont hinder posterior analysis. 

First step is creating the necessary directories

```
mkdir analysis_no_mit_chl
mkdir analysis_no_mit_chl/feature_table_no_mit_chl
mkdir analysis_no_mit_chl/identification
```

Once created, we generate a new feature table containing all taxids except for mitochondria and chloroplast, and the visualization files for it.

```
qiime taxa filter-table \
  --i-table dada2_result/feature_table_dada2.qza \
  --i-taxonomy taxonomy/taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table analysis_no_mit_chl/feature_table_no_mit_chl/feature_table_no_mitochondria_no_chloroplast.qza

qiime metadata tabulate \
  --m-input-file analysis_no_mit_chl/feature_table_no_mit_chl/feature_table_no_mitochondria_no_chloroplast.qza \
  --o-visualization analysis_no_mit_chl/feature_table_no_mit_chl/feature_table_no_mitochondria_no_chloroplast.qzv

qiime feature-table summarize \
--i-table analysis_no_mit_chl/feature_table_no_mit_chl/feature_table_no_mitochondria_no_chloroplast.qza \
--o-visualization analysis_no_mit_chl/feature_table_no_mit_chl/feature_table_no_mitochondria_no_chloroplast_summarized.qzv \
--m-sample-metadata-file $METADATA
```

As we did when we had the mitochondria and chloroplast, we follow the same few steps:

- Generate the bar chart with the sample composition (with the newly created, filtered feature table)

```
qiime taxa barplot \
--i-table analysis_no_mit_chl/feature_table_no_mit_chl/feature_table_no_mitochondria_no_chloroplast.qza \
--i-taxonomy taxonomy/taxonomy.qza \
--m-metadata-file $METADATA \
--o-visualization analysis_no_mit_chl/identification/taxa-bar-plots.qzv
```
- Export the absolute number of features 

```
qiime taxa collapse \
  --i-table analysis_no_mit_chl/feature_table_no_mit_chl/feature_table_no_mitochondria_no_chloroplast.qza \
  --i-taxonomy taxonomy/taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table analysis_no_mit_chl/identification/table-no-mitochondria-no-chloroplast_spp_lvl.qza

qiime tools export \
  --input-path analysis_no_mit_chl/identification/table-no-mitochondria-no-chloroplast_spp_lvl.qza \
  --output-path analysis_no_mit_chl/identification

biom convert \
--input-fp analysis_no_mit_chl/identification/feature-table.biom \
--output-fp analysis_no_mit_chl/identification/transposed_table_no_mitochondria_no_cloroplast.tsv \
--to-tsv

python process_table.py \
   analysis_no_mit_chl/identification/transposed_table_no_mitochondria_no_cloroplast.tsv \
   analysis_no_mit_chl/identification/absolute_numbers_taxonomy_no_mitochondria_no_cloroplast.tsv 

```
- Calculation of diversity data (remember to check the `SAMPLING_DEPTH`, this time from `analysis_no_mit_chl/feature_table_no_mit_chl/feature_table_no_mitochondria_no_chloroplast_summarized.qzv` instead)

```
mkdir analysis_no_mit_chl/diversity_data

qiime diversity core-metrics-phylogenetic \
--i-phylogeny phylogeny_data/rooted_tree.qza \
--i-table dada2_result/feature_table_dada2.qza \
--p-sampling-depth $SAMPLING_DEPTH \
--m-metadata-file $METADATA \
--output-dir analysis_no_mit_chl/diversity_data

qiime diversity alpha-group-significance \
--i-alpha-diversity analysis_no_mit_chl/diversity_data/faith_pd_vector.qza \
--m-metadata-file $METADATA \
--o-visualization analysis_no_mit_chl/diversity_data/faith_pd_summarized.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity analysis_not_mit_chl/diversity_data/shannon_vector.qza \
--m-metadata-file $METADATA \
--o-visualization analysis_not_mit_chl/diversity_data/shannon_vector_summarized.qzv
```

- Calculation of beta-diversity data (remember to put a valid metadata column inside of the `GROUP`variable)

```
qiime diversity beta-group-significance \
--i-distance-matrix analysis_not_mit_chl/diversity_data/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file $METADATA \
--m-metadata-column $GROUP \
--o-visualization analysis_mpt_mit_chl/diversity_data/unweighted_unifrac_species_significance.qzv \
--p-pairwise
```

- Alpha rarefaction calculation (this time, the variable `MAX_DEPTH`, must be extracted from `analysis_no_mit_chl/feature_table_no_mit_chl/feature_table_no_mitochondria_no_chloroplast_summarized.qzv`)


```
mkdir analysis_no_mit_chl/alpha_rarefaction

qiime diversity alpha-rarefaction \
--i-table analysis_no_mit_chl/feature_table_no_mit_chl/feature_table_no_mitochondria_no_chloroplast.qza \
--i-phylogeny phylogeny_data/rooted_tree.qza \
--p-max-depth $MAX_DEPTH \
--m-metadata-file $METADATA \
--o-visualization analysis_no_mit_chl/alpha_rarefaction/alpha_rarefaction.qzv
```