cd <path_to_halla_directory> 

## HAllA version # 0.8.20 (https://github.com/biobakery/halla)

# phylum
halla -x inputs/halla_phylum_table.txt -y inputs/halla_metadata.txt -o outputs/halla_out_phylum_q-0.05

# genus
halla -x inputs/halla_genus_table.txt -y inputs/halla_metadata.txt -o outputs/halla_out_genus_q-0.05

# genus (alpha = 0.01)
halla -x inputs/halla_genus_table.txt -y inputs/halla_metadata.txt -o outputs/halla_out_genus_q-0.01 --fdr_alpha 0.01

