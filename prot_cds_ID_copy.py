#problem: dna (CDS) and protein records do not have the same ID
#need to get a fasta file with PROTEIN IDS BUT NT SEQUENCES for search

#approach: translate CDS and match with the protein records
#then I can find the transcript ID in the CDS file and change it

from Bio import SeqIO

out1 = "CDS_proteinID_ALL.fasta" 
cds = SeqIO.parse(open("CDS_WS.fasta"),"fasta")
prot = SeqIO.parse(open("proteins_WS.fasta"),"fasta")

cdsdict = {}
cdsdicttrans = {}
for rec in cds:
    bits = rec.id.split("|") 
    rec.id = bits[2] #because JGI ids are freakishly long. this is the id no.
    #rec.description = bits[2]
    cdsdict[rec.id] = [rec.seq]
    if len(rec.seq) % 3 == 0: #is the seq length a multiple of 3 (codon)?
        aaseq = rec.seq.translate()
        cdsdicttrans[rec.id] = [str(aaseq), rec.seq] #seqs with ids in dict.
print len(cdsdict)
print len(cdsdicttrans)

#cdsprotid = {}
counter = 0
with open(out1,"w") as f:
    for cand in prot:
        for key, value in cdsdicttrans.items():
            if value[0] == str(cand.seq): #remember value is list not string
                #print key
                counter += 1
                #bitsp = cand.id.split("|")
                #cand.id = bitsp[2]
                #cand.description = bits[2]
                cand.seq = value[1]
                cand.name = key #to keep the transcript id info for next step
                cand.description = key
                SeqIO.write(cand, f, 'fasta')
                #works \m/

print counter
