import sys, random, FASTAhandle
# Creates and outputs a dictionary for nucleic acid triplets. The dictionary is what is used.
class codonfreq:
    
    def __init__(self, codonfreq):
        self.myDic = {}
        self.myRawFile = open(codonfreq).readlines()
        self.popDic()
        self.AADict = {"I":["ATT","ATC","ATA"],"L":["CTT","CTC","CTA","CTG","TTA","TTG"],"V":["GTT","GTC","GTA","GTG"],"F":["TTT","TTC"],"M":["ATG"],"C":["TGT","TGC"],"A":["GCT","GCC","GCA","GCG"],"G":["GGT","GGC","GGA","GGG"],
                  "P":["CCT","CCC","CCA","CCG"], "T":["ACT","ACC","ACA","ACG"],"S":["TCT","TCC","TCA","TCA","TCG","AGT","AGC"],"Y":["TAT","TAC"],"W":["TGG"],"Q":["CAA","CAG"],"N":["AAT","AAC"],"H":["CAT","CAC"],"E":["GAA","GAG"],
                  "D":["GAT","GAC"],"K":["AAA","AAG"],"R":["CGT","CGC","CGA","CGG","AGA","AGG"]}
    def popDic(self):
        for line in self.myRawFile:
            freq = line[3:]
            freq = freq[:len(freq)-1]
            self.myDic[line[0:3]] = freq

        
