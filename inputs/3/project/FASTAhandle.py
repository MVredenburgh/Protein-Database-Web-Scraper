# FASTA File handeling.
class FASTAhandle:
    def __init__(self, myFileName):
        self.myGenes = []
        infile = open(myFileName)
        self.myRawFile = infile.readlines()
        self.parseGene()
    def parseGene(self):
        # Creates an array of Gene objects called myGenes.
        n=0
        HeaderList = []
        for line in self.myRawFile:
            if line[0][:1]=='>':
                HeaderList.append(n)
            n+=1
        HeaderList.append(len(self.myRawFile)+1)
        n=0
        while n < len(HeaderList)-1:
            self.myGenes.append(Gene(self.myRawFile[HeaderList[n]:HeaderList[n+1]]))
            n+=1

class Gene:
    def __init__(self, FASTAchunk):
        # FASTAchunk is the header line plus the genetic infomation.
        self.mySequence = ""
        self.header = FASTAchunk[0]
        self.makeSequence(FASTAchunk)
        self.addLength()
        self.reverse()
        
    def makeSequence(self, FASTAchunk):
        # Makes a long string out of the FASTAarray. Also removes all instances of "\n".
        FASTAgene=FASTAchunk[1:]
        for line in FASTAgene:
            line = line.replace("\n","")
            self.mySequence += line
    def addLength(self):
        myLength = len(self.mySequence)
        self.header = self.header.replace("\n","")
        self.header += " %d bp" % myLength
    def reverse(self):
        # Reverses the sequence. 
        self.reverseSequence = self.mySequence[::-1]
    def makeReverseComplement(self):
        # Creates a DNA and RNA dictionary
        DNAcomplement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'R': 'Y', 'W':'W', 'Y': 'R', 'S':'S', 'B':'V', 'V':'B', 'H':'D', 'X':'N', 'N':'N', '.':'.'}
        RNAcomplement = {'A': 'U', 'T':'A', 'U': 'A', 'C': 'G', 'G': 'C', 'R': 'Y', 'W':'W', 'Y': 'R', 'S':'S', 'B':'V', 'V':'B', 'H':'D', 'X':'N', 'N':'N', '.':'.'}
        try:
            for BasePair in self.reverseSequence:
                self.reverseComplement += DNAcomplement[BasePair]
        except KeyError:
            print ("Only FASTA files are supported. ")
