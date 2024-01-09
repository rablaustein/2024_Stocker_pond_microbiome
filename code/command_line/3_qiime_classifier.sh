module load qiime 2-2021.4

cd <path_to_directory> # have the reference_dbs directory present and containing 'SILVA_132_QIIME_release'

### Build classifier DB for 'SILVA 97'

## import silva db reference files (download latest release from https://www.arb-silva.de/download/archive/qiime/)
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path reference_dbs/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna \
  --output-path reference_dbs/silva-97.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path reference_dbs/SILVA_132_QIIME_release/taxonomy/16S_only/97/taxonomy_all_levels.txt \
  --output-path reference_dbs/silva-97-reference-taxonomy.qza

## classifier extract sequences
qiime feature-classifier extract-reads \
  --i-sequences reference_dbs/silva-97.qza \
  --p-f-primer CCTAYGGGDBGCWGCAG --p-r-primer GACTACNVGGGTMTCTAATCC --p-trunc-len 465 \
  --o-reads reference_dbs/silva-97-reference-seqs.qza

## train classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads reference_dbs/silva-97-reference-seqs.qza \
  --i-reference-taxonomy reference_dbs/silva-97-reference-taxonomy.qza \
  --o-classifier reference_dbs/silva-97-classifier.qza