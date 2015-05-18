import sys
# Michael Vredenburgh Homework 3
# The Reverse Complement Generator

class InputData:
    fileName = ""
    typeDNAorRNA = ""
    myRawFile = [""]
    myGenes = []

    def __init__(self, myFileName, myTypeDNAorRNA):
        self.fileName = myFileName
        self.typeDNAorRNA = myTypeDNAorRNA
        self.myRawFile = InputData.extractSequence(myFileName)
        self.parseGene()
        self.toString()

    def extractSequence(myFileName):
        # Parses the FASTA file, and creates a string for the sequence.
        infile = open(myFileName)  # Creates a variable to store the opened file to reference later
        myRawFile= infile.readlines()
        return myRawFile
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
            self.myGenes.append(Gene(self.myRawFile[HeaderList[n]:HeaderList[n+1]-1],self.typeDNAorRNA))
            n+=1
        

    def toString(self):
        # prints to the command line the reverse complement FASTA file.
        for gene in self.myGenes:
            print (gene.header)
            while len(gene.reverseComplement)>0:
                if len(gene.reverseComplement)>=80:
                    print (gene.reverseComplement[0:80])
                    gene.reverseComplement = gene.reverseComplement[80:]
                else:
                    print (gene.reverseComplement)
                    break
class Gene:
    header = ""
    mySequence = ""
    reverseSequence = ""
    reverseComplement = ""
    typeDNAorRNA = ""
    def __init__(self, FASTAchunk, typeDNAorRNA):
        # FASTAchunk is the header line plus the genetic infomation.
        self.typeDNAorRNA = typeDNAorRNA
        self.header = FASTAchunk[0]
        self.makeSequence(FASTAchunk)
        self.addLength()
        self.reverse()
        self.makeReverseComplement()
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
            if self.typeDNAorRNA == "RNA":
                for BasePair in self.reverseSequence:
                    self.reverseComplement += RNAcomplement[BasePair]
            else:
                for BasePair in self.reverseSequence:
                    self.reverseComplement += DNAcomplement[BasePair]
        except KeyError:
            print ("Only DNA FASTA files are supported. You likely used a protein FASTA file. DNA coding files will still be printed.")


# Does error checking on the command line input and runs the program.
if len(sys.argv)>3:
    raise "The maximum number of inputs is two. Please enter the file path and an optional type declaration."
if len(sys.argv) == 1:
    raise "You must add a file path do your FASTA file."
if len(sys.argv) == 2:
    try:
        myInput = InputData(sys.argv[1],"")
    except FileNotFoundError:
        print ("The file path is not valid, please try again.")
if len(sys.argv) == 3:
    try:
        myInput = InputData(sys.argv[1],sys.argv[2])
    except FileNotFoundError:
        print ("The file path is not valid, please try again. Or, the second input can only be DNA or RNA.")  



