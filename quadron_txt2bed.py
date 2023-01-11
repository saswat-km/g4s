# This file converts the output file generated 
# from the Quadron analysis into a .bed format
# file for further genomic analysis


data = open('mGor_X_output.txt', 'r')
data = data.readlines()

f = open("mGor_chrX_g4_quadron.bed","w") 

adjust = 0 #adjust for quadron error
for i in range(0,len(data)):
    crit = data[i].split()
    if crit[0] == "DATA:" and crit[4] != "NA":
        start = int(crit[1])
        length = int(crit[3])
        f.write("chrX\t%i\t%i\t%.2f\t%s\t%i\n" %(start-adjust, length+start-adjust, float(crit[4]), crit[2], length))
