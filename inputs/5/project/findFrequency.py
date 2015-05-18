import sys, FASTAhandle
# Takes in a protein.fasta and outputs a two dimensional array.
# The first element of the array is the single letter abriviation for the most frequent
# amino acid. The second element of the array is the frequency it occurs.
class proteinFreq:
    def __init__(self, proteinFileName):
        self.FASTAobject = FASTAhandle.FASTAhandle(proteinFileName)
        self.numOfGenes = 0
        for gene in self.FASTAobject.myGenes:
            self.numOfGenes += 1
        self.lenOfGenes = len(self.FASTAobject.myGenes[0].mySequence)
        self.makeOutputArray()
    def makeOutputArray(self):
        freq=[[0,0] for i in range(0,self.lenOfGenes)]
        for x in range(0,self.lenOfGenes):
            self.L=[]
            for gene in self.FASTAobject.myGenes:
                self.L.append(gene.mySequence[x])
            freq[x][0]=self.most_common()
            freq[x][1]=self.L.count(self.most_common())/len(self.L)
            self.freqTable = freq
    def most_common(self):
        AAfreq={}
        for AminoAcid in self.L:
            if AAfreq.__contains__(AminoAcid) == False:
                AAfreq[AminoAcid] = 1
            else:
                AAfreq[AminoAcid] += 1
        highest = ""
        highestFreq = 0
        for AA in AAfreq:
            if AAfreq[AA]>highestFreq:
                highest = AA
                highestFreq = AAfreq[AA]           
        return highest
