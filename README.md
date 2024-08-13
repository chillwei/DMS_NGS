# DMS_NGS
NGS analysis and data visualization of DMS library

## Merge reads of pair ended data 
use 'Fastp' tool for input and output.
see link: https://github.com/OpenGene/fastp

    fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz

## Further data filtering based on v-find
run find-variant_QW.py in terminal

change variables in find-variant_QW.py 
line 113 - 115 :  adatpor sequence for NGS and input file name. 
line 119: output as a .csv file, change the filename as needed


    if __name__ == "__main__":
        adapters = Adapters("GGCGCGGTGTTAAAT", "CATCATCACCATCACCAT" )
        fq_1 = "smurfp_ngs_merge.fq.gz"
        #fq_2 = "30-982598147/17L_R2_clean.fq.gz"
        
        start = time.perf_counter()
        find_variants(adapters, fq_1, save_path = "filtered_smurfp.csv")
