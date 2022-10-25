# This takes a BigWig file and a .bed 
# file, and extracts the phastcon scores
# for the regions between the regions 
# defined between the .bed file

import pyBigWig
import numpy as np
import pandas as pd
from time import time

bw = pyBigWig.open("hg38.phastCons30way.bw") #customizable

df = pd.read_csv('../input_files/Stochastic_nonCORE_origins_hg38.bed', sep='\t', comment='t', header=None) #reading the template file #customizable
header = ['chrom', 'chromStart', 'chromEnd', 'annotate'] #'quadron', 'sign', 'length', 'phastcon'] #customizable
df.columns = header[:len(df.columns)]
phastcon = np.empty([len(df)-1,4],dtype="<U9")

f = open("stochastic_non-ro_error_file.txt", "w") #customizable
s = time()
for i in range(0,len(df)-1):
    try:
        region1 = np.array(df.loc[[i]][['chrom','chromStart','chromEnd']])[0]
        region2 = np.array(df.loc[[i+1]][['chrom','chromStart','chromEnd']])[0]
        if region1[0]==region2[0]:
                region_phastcon = bw.stats(region1[0],region1[2],region2[1],exact=True)
                phastcon[i] = np.array([region1[0],region1[2],region2[1],np.round(region_phastcon,7)[0]]) #phastcon[i] = np.round(region_phastcon,7)
        else:
                f.write("Pass message| Position: %s %i %s %i \n" %(region1[0],region1[2],region2[0],region2[1])) #customizable
                pass
    except Exception as e:
        f.write("Error message: %s| Position: %s %i %i \n" %(e,region1[0],region1[2],region2[1])) #customizable
        pass

f.close()

ef = pd.DataFrame(phastcon)
header = ['chrom', 'chromStart', 'chromEnd', 'phastcon']
ef.columns = header[:len(ef.columns)]
ef.to_csv("stochastic_non-ro_chr_all_phastcon.bed", sep="\t", header=False, index=False) #customizable
e = time()
print("Time:", round((e-s)/60,3), "mins \n")
