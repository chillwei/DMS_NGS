# DMS_NGS
NGS analysis and data visualization of DMS library

## Merge reads of pair ended data 
use 'Fastp' tool for input and output.
see link: https://github.com/OpenGene/fastp

    fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz

## Further data filtering based on v-find
run find-variant_QW.py in terminal

change variables in find-variant_QW.py based on the adatpor sequence for NGS and the merged fastq file name.

line 113 - 115

    if __name__ == "__main__":
        adapters = Adapters("GGCGCGGTGTTAAAT", "CATCATCACCATCACCAT" )
        fq_1 = "smurfp_ngs_merge.fq.gz"
