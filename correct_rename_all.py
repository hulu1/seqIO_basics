#this renames fasta records named 1,2,3 etc. to the id specified in the filename and adds the country code for each number
#can probably be made more quickly, this took ~40 min for 4000 files
#note: I swapped GRE and GUA1 in these procedures (i.e. correct)

from Bio import SeqIO
import glob #filename pattern matching module

counter = 0
filenamelist = []
for filename in glob.glob('rawfasta/*.txt.fasta'):
    #print filename
    counter += 1
    bits = filename.split(".")
    goi = bits[0].split("/")[1] #dummy way of getting only the actual filename
    #goi will be a str but making sure below
    filenamelist.append(goi)
print counter, "filenames processed"

goicounter = 0
def stage2 ( goi ):
    "The function makes a copy of the nt fasta file and a corresponding aa file while it changes the seq ids from 1,2,3 etc. to 1_ALP, 2_BHU etc. and adds the NZE10 seq."
    seq_in = "rawfasta/"+str(goi)+".txt.fasta" 
    seq_out = "dNdSrenamed/"+"Ds"+str(goi)+".fasta"
    seq_out2 = "dNdSrenamed/"+"Ds"+str(goi)+"aa.fasta"

    nz10 = SeqIO.parse(open("Ds_CDS_proteinshortID_ALL.fasta"),"fasta")
    for recs in nz10:
        foo = recs.id[2:]
        if goi == foo:
            nt_nz10 = recs
            nt_nz10.id = recs.id+"_NZE10"
            nt_nz10.name = recs.id+"_NZE10"
            fubar = recs.description.split(" ")
            nt_nz10.description = fubar[1]
            print goi, recs
            #description was id/name+original description which is the transcript id
##        else: #gene id might not be in catalog for some reason
##            nt_nz10 = str(goi)+"empty"
##            print goi, "NZE10 ref sequence not found"
    raw = SeqIO.parse(open(seq_in),'fasta')
    with open(seq_out,'w') as f:
##        if type(nt_nz10) is str:
##            print goi, "NZE10 ref sequence not found"
##        else:
        SeqIO.write(nt_nz10, f, 'fasta')
        for rec in raw:
            if rec.id == "1":
                rec.id = "10_GUA1"
            elif rec.id == "2":
                rec.id = "11_GRE1"
            elif rec.id == "3":
                rec.id = "12_GUA2"
            elif rec.id == "4":
                rec.id = "13_NZE2"
            elif rec.id == "5":
                rec.id = "14_NZE8"
            elif rec.id == "6":
                rec.id = "15_RUS1"
            elif rec.id == "7":
                rec.id = "16_SAF4"
            elif rec.id == "8":
                rec.id = "17_SLV1"
            elif rec.id == "9":
                rec.id = "18_USA12"
            elif rec.id == "10":
                rec.id = "01_ALP3" 
            elif rec.id == "11":
                rec.id = "02_AUS4"
            elif rec.id == "12":
                rec.id = "03_BHU1"
            elif rec.id == "13":
                rec.id = "04_CAN3"
            elif rec.id == "14":
                rec.id = "05_CHI17"
            elif rec.id == "15":
                rec.id = "06_COLN"
            elif rec.id == "16":
                rec.id = "07_COLS"
            elif rec.id == "17":
                rec.id = "08_DEN1"
            elif rec.id == "18":
                rec.id = "09_ECU13"
            rec.name = rec.id
            rec.description = "Ds"+str(goi) #keep the ID, might be needed
            SeqIO.write(rec, f, 'fasta')

        nz10aa = SeqIO.parse(open("proteins_shortid.fasta"),'fasta')
        for recsa in nz10aa:
            if goi in recsa.id:
                aa_nz10 = recsa
            else:
                pass
    ntseq = SeqIO.parse(open(seq_out),'fasta')
    with open(seq_out2,'w') as g:
        for rec in ntseq:
            aarec = rec
            aarec.seq = aarec.seq.translate()
            SeqIO.write(aarec, g, 'fasta')
            if len(str(rec.seq)) % 3 != 0: #just to be sure
                print goi, "CDS not multiple of three"
    #check nz10 translation
    fu = SeqIO.parse(open(seq_out2),'fasta')
    for rec0 in fu:
        if "_NZ10" in rec0.id:
            if rec0.seq == aa_nz10.seq:
                pass
            else:
                print "warning, translated NZ10 sequence does not match protein sequence"

    #print str(goi), "stage 2 processed"
    global goicounter
    goicounter +=1
    return;


for gois in filenamelist[:10]:
    stage2( gois )
##    if goicounter % 100 == 0: 
    print goicounter, "files of", counter, "actually copied"

print "final count:", goicounter
