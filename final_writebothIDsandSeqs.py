#don't think I can parse through two fasta files at the same time
#plan was to add the aa seq on the same line

from Bio import SeqIO

out1 = "TESTDs_CDS_proteinshortID_ALL.fasta"
out2 = "TESTDs_proteinID_CDS.txt"
seqs = SeqIO.parse(open('Ds_CDS_proteinID_ALL.fasta'), 'fasta')
with open(out1,"w") as f:
    with open(out2,"w") as g:
        for rec in seqs:
            bits = rec.id.split("|")
            rec.id = "Ds"+bits[2] #Ds is species
            desc = rec.description.split(" ")
            rec.name = desc[1]
            rec.description = desc[1]
            SeqIO.write(rec,f,'fasta')
            g.write(rec.id)
            g.write('; ')
            g.write(rec.description)
            g.write('; ')
            g.write(str(rec.seq))
            g.write('\n') #amateur writing


        
