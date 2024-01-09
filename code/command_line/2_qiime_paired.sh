module load qiime 2-2021.4

cd <path_to_seq_reads_directory>

## set new directories
mkdir qiime_out
mkdir qiime_out/qiime_files
mkdir qiime_out/qiime_viz

## import sequencing reads as qiime file (input is fastq files in the 'cutadapt_out' directory)
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path cutadapt_out \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path qiime_out/qiime_files/demux-paired-end.qza

## move to qiime_out directory
cd qiime_out

## vizualize data quality
qiime demux summarize --i-data qiime_files/demux-paired-end.qza --o-visualization qiime_viz/demux-paired-end_viz.qzv

## dada de-noise paired-end
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime_files/demux-paired-end.qza \
  --p-trunc-q 20 --p-trunc-len-f 300 --p-trunc-len-r 225 \
  --p-n-threads 18 \
  --o-table qiime_files/dada-table.qza \
  --o-representative-sequences qiime_files/dada-rep-seqs.qza \
  --o-denoising-stats qiime_files/dada-denoising-stats.qza

## dada summarize feature table
qiime feature-table summarize \
  --i-table qiime_files/dada-table.qza \
  --o-visualization qiime_viz/dada-table.qzv \
  --m-sample-metadata-file ../sample-metadata.tsv

## dada tabulate feature table
qiime feature-table tabulate-seqs \
  --i-data qiime_files/dada-rep-seqs.qza \
  --o-visualization qiime_viz/dada-rep-seqs.qzv
