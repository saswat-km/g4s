import pyBigWig
import numpy as np
import pandas as pd
from time import time

bw = pyBigWig.open("hg38.phastCons30way.bw")

df = pd.read_csv('chr_all_non_G+C-_stochastic_overlap.bed', sep='\t', comment='t', header=None)
header = ['chrom', 'chromStart', 'chromEnd', 'quadron', 'sign', 'length']
df.columns = header[:len(df.columns)]
phastcon = np.zeros(len(df))

f = open("chr_all_non_G+C-_stochastic_overlap_error_file.txt", "w")
s = time()
for i in range(0,len(df)):
    try:
        region = np.array(df.loc[[i]][['chrom','chromStart','chromEnd']])[0]
        region_phastcon = bw.stats(region[0],region[1],region[2],exact=True)
        phastcon[i] = np.round(region_phastcon,7)
    except Exception as e:
        f.write("Error message: %s| Position: %s %i %i \n" %(e,region[0],region[1],region[2]))
        pass

f.close()
df['phastcon'] = pd.DataFrame(np.transpose([phastcon]))
df.to_csv("chr_all_non_G+C-_stochastic_overlap_phastcon.bed", sep="\t", header=False, index=False)
e = time()
print("Time:", round((e-s)/60,3), "mins \n")
