"""setting:
got a fasta file of >12000 sequences from JGI database. their IDs are in this form:
'jgi|GenSpec1|67409|estExt_fgenesh1_kg.C_1_t10001'
this is very large to work with in other programs for my taste so I want to cut
it down to the protein id which is the first number. NB this info is also stored
in 'name' and 'description' so I am keeping it in the record.
I only need to do this for <60 sequences, which obviously I want to write to a
new fasta file.
so I have this list in a txt file, however the numbers have 2 letters in front of them
for clarity but I don't need this here (yes sure I don't need python to deal
with that but it feels better)"""


infile = "proteins.fasta" #collection of all proteins in the genome
out = "cand_march16.fasta" #desired final file (duh)

idstrings = []
with open("cand_ids.txt", "r") as c: #list of proteins of interest
    for ids in c:
        idstrings.append(str(ids[2:].strip())) #now it's a python list
                                            #and first 2 characters removed
   
from Bio import SeqIO
counter = 0
controllist = []
protseqs = SeqIO.parse(open(infile),"fasta")
with open(out,"w") as f:
    for rec in protseqs:
        if rec.id[11:16] in idstrings or rec.id[11:17] in idstrings:#1
            bits = rec.id.split("|")#2
            rec.id = bits[2]
            rec.description = bits[2]#3
            SeqIO.write(rec, f, 'fasta')
            counter += 1
            controllist.append(rec.id)

#1 searching for the protein id which is 5 or 6 digits (I guess there's a more elegant way)
#2 super convenient way of dissecting the long id
#3 renaming each record while copying (rather than renaming everything)
#...which I have later done elsewhere
#changing description as well seems to be essential as it will be displayed otherwise
#rec.name still contains the full info!
#print controllist

        
