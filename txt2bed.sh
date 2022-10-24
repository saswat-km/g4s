# The file takes all the individual chr files in .txt 
# format and converts them to .bed format and merges 
# them and finally sorts them using sortBed of Bedtools

#!/bin/bash
for i in {1..22}
do
        awk -v r=$i '{if ($1=="DATA:") print "chr" r "\t" $2 "\t" $2+$4 "\t" $5 "\t" $3 }' output.chr$i.txt > chr$i.bed
done
cat chr*.bed > chr_all.bed
sortBed -i chr_all.bed > chr_all_sorted_hg38_g4_quadron.bed
