#my Python approach to find BsaI binding sites for GoldenGate cloning.
#BsaI is a type IIs restriction enzyme
#it cuts 1 bp after the the recognition site, creating a 4 bp sticky overhang
#(and removing the site)
#the OH and their positions in the seq will be returned
#this could be easily changed to other REs

from Bio import SeqIO
    
nt_cand = SeqIO.parse(open("genesofinterest_nt.fasta"),"fasta")
fw = "GGTCTC"
rev = "GAGACC" #recognition sites
ohlist = ['AATG','CAAG','GCTT'] #could make this a dict when extending
                                # to indicate incompatibility
#list of overhangs that must be avoided (because other modules have them)
all_list = []
with open("cand_bsa-oh_enum.txt","w") as g:
    for recs in nt_cand:
        #wacken = [] #misplaced
        s = recs.seq
        t = s.complement() #actually the complement is not necessary!!
        for i, j in enumerate(s):
            if s[i:i +len(fw)] == fw:
                wacken = []
                g.write(recs.id)
                g.write("\tfw\t")
                g.write(str(i+1))
                g.write("\t")
                pos1 = i+7
                pos2 = pos1+1
                pos3 = pos1+2
                pos4 = pos1+3
                wacken.append([s[pos1],s[pos2],s[pos3],s[pos4],"/", \
                               t[pos1],t[pos2],t[pos3],t[pos4]])
                fu = wacken[0]
                #bar = fu[0]+fu[1]+fu[2]+fu[3]+fu[4]+fu[5]+fu[6]+fu[7]+fu[8] #crappy
                bar = ''.join(fu) #elegant
                g.write(bar[:4])
                all_list.append(bar[:4])
                all_list.append(bar[5:10])
                if bar[:4] in ohlist:#or bar[5:10] in ohlist:
                    g.write("\tcheck")
##                if bar[:4] in all_list: #this doesn't make sense it's always true
##                    g.write("\tcheck internal")
                g.write("\n")
            
            if s[i:i +len(rev)] == rev:
                wacken = []
                g.write(recs.id)
                g.write("\trev\t")
                g.write(str(i+1))
                g.write("\t")
                pos1 = i-5
                pos2 = pos1+1
                pos3 = pos1+2
                pos4 = pos1+3
                wacken.append([s[pos1],s[pos2],s[pos3],s[pos4],"/", \
                               t[pos1],t[pos2],t[pos3],t[pos4]])
                fu = wacken[0]
                bar = ''.join(fu)
                g.write(bar[:4])
                all_list.append(bar[:4])
                all_list.append(bar[5:10])
                if bar[:4] in ohlist:#"or bar[5:10] in ohlist:
                    g.write("\tcheck")
##                if bar[:4] in all_list:
##                    g.write("\tcheck internal")
                g.write("\n")

            #this works so far. watch out index is 0 based
            #need "self" check (OH which are not in modules but in gene i.e. other sites)
