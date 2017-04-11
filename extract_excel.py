#gets sequence records from large fasta file based on Excel sheet list

infile = "proteins.fasta" #collection of all proteins in the genome
out = "Kseqs1.fasta"
out2 = "Kseqs1_shortids.fasta"

import xlrd #Excel handling module
from itertools import takewhile
workbook = xlrd.open_workbook('K_genes.xlsx')
worksheet = workbook.sheet_by_name('Sheet1')
#inter = [worksheet.row_values(i)[3] for i in range(worksheet.nrows)]
#this was an alternative way I tried
column_generator = (worksheet.row_values(i)[0] for i in range(worksheet.nrows))
ids = list(takewhile(str, column_generator))
#this simply writes the strings in column 1 to the list
idno = []
idstrings = []
for val in ids:
    #bits = val.split("_")
    #idstrings.append(val[2:]) #works = ID numbers of secr. prot only
    idstrings.append(str(int(val))) #make sure only numbers (integers) are used
   
from Bio import SeqIO
counter = 0
controllist = []
protseqs = SeqIO.parse(open(infile),"fasta")
with open(out,"w") as f:
    for rec in protseqs:
        if rec.id[11:16] in idstrings or rec.id[11:17] in idstrings:
            SeqIO.write(rec, f, 'fasta')
            counter += 1
            #controllist.append(rec.id)
            
fail = []          
countercontrol = 0
protseqs2 = SeqIO.parse(open(out),"fasta")
with open(out2,'w') as g:
    for rec in protseqs2:
        bits = rec.id.split("|")
        rec.id = bits[2]
        rec.description = bits[2]
        SeqIO.write(rec,g,'fasta')
        countercontrol += 1

print "fail:", countercontrol


        
