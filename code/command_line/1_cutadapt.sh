module load cutadapt 4.0

cd <path_to_seq_reads_directory> 

## Cutadapt: remove residual adapters, trim at phred score 20, remove reads less than 100bp  

mkdir cutadapt_out

cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG --trim-n -m 100 -q 20 -o cutadapt_out/$Read_1_cut -p cutadapt_out/$Read_2_cut reads_raw/$Read_1_raw reads_raw/$Read_2_raw