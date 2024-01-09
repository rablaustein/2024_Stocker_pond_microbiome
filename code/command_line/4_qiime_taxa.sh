module load qiime 2-2021.4

cd <path_to_directory_containing_the_qiime_out_and_reference_dbs_directories>

### Taxonomic analysis using the Silva 97 classifier

## taxonomic analysis
qiime feature-classifier classify-sklearn \
  --i-classifier reference_dbs/silva-97-classifier.qza \
  --i-reads qiime_out/qiime_files/dada-rep-seqs.qza \
  --o-classification qiime_out/qiime_files/taxonomy-silva-97.qza

## taxonomy print
qiime metadata tabulate \
  --m-input-file qiime_out/qiime_files/taxonomy-silva-97.qza \
  --o-visualization qiime_out/qiime_viz/taxonomy-silva-97.qzv
  
## taxonomy barplot
qiime taxa barplot \
  --i-table qiime_out/qiime_files/dada-table.qza \
  --i-taxonomy qiime_out/qiime_files/taxonomy-silva-97.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization qiime_out/qiime_viz/taxa-bar-plots-silva-97.qzv