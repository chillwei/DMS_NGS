# DMS_NGS
NGS analysis and data visualization of DMS library

## Merge reads of pair ended data 
use 'Fastp' tool for input and output.
see link: https://github.com/OpenGene/fastp

fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz

## Further data filtering based on v-find
