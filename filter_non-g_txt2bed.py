# Takes all chromosome file data with quadron scores 
# as input and extracts the non-G (+ strand) and non-C  
# (- strand) and converts into a single .bed file format. 

import numpy as np
import pandas as pd

def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

def consecutive_greater2(data, stepsize=1):
    a = np.split(data, np.where(np.diff(data) != stepsize)[0]+1)
    k = []
    for i in a:
        if len(i) > 2:
            for j in i:
                k.append(j)
    return(consecutive(k))

def connect(data):
    data_new = []
    for i in range(0,len(data)-1):
        start = data[i][-1]+1
        end = data[i+1][0]-1
        series = np.arange(start,end+1,1)
        data_new.append(series)
    return(data_new)

f = open("chr_all_G/chr_all_non_G+C-.bed","w")

for chr in range(1,23):
    data = open('input_files/output.chr%s.txt' %chr,'r')
    data = data.readlines()
    begin_prev = 0

    if chr > 9: #adjust for matthias error
        adjust = 6
    else:
        adjust = 5

    for i in range(0,len(data)):
        crit = data[i].split()
        if crit[0] == "DATA:":
            start = int(crit[1])
            seq = crit[5]
            count = []
            val = 0
            for i in seq:
                val += 1
                if i == "G" and crit[2] == "+":
                    count.append(val)
                if i == "C" and crit[2] == "-":
                    count.append(val)
            cons = consecutive_greater2(count)
            inv_cons = connect(cons)
            for j in inv_cons:
                length = len(j)
                var = start - 1 - adjust
                if length > 0:
                    begin = j[0] + var
                    end = j[-1] + var
                    if begin > begin_prev:
                        f.write("chr%s\t%i\t%i\t%.2f\t%s\t%i\n" %(chr,begin,end+1,float(crit[4]),crit[2],length))
                        begin_prev = begin
    print("Completed Chr: %i" %chr)

f.close()

#df = pd.read_csv('chr_all_G/chr_all_G+C-.bed', sep='\t', comment='t', header=None)
#header = ['chrom', 'chromStart', 'chromEnd', 'quadron', 'sign', 'length']
#df.columns = header[:len(df.columns)]
