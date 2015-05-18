import sites, codonfreq, effector, sys, findFrequency, random
# Michael Vredenburgh, Project 2
# The driver runs all other associated files and contains the biologically relevant code.

class sequenceGenerator:
    def __init__(self):
        # Runs all axillary files and then executes main loop.
        self.mySequence = ""
        self.mysites = sites.sites("sites.fasta").mySitesList
        self.myInputFreq = codonfreq.codonfreq("codonfreq.txt")
        self.myCalcFreq = findFrequency.proteinFreq("protein.fasta")
        self.myEffector = effector.effector("effector.fasta").myEffector
        self.params = int(open("params.txt").readlines()[0].replace("\n",""))
        self.makeFASTA()
    def makeFASTA(self):
        # Outputs the file output.fasta in the fasta format.
        FASTAstring =""
        for i in range(0,self.params):
            FASTAstring += "> sequence %d\n" %(i+1)
            self.genSequence()
            self.sitesCheck(10)
            # Add check for restriction sites
            while len(self.mySequence)> 80:
                FASTAstring += self.mySequence[0:80] +"\n"
                self.mySequence = self.mySequence[80:]
            FASTAstring += self.mySequence + "\n"
            
        outputFile = open("output.fasta","w")
        outputFile.write(FASTAstring)
        outputFile.close()
    def sitesCheck(self, n):
        # If a forbid site is detected, new sequences will be recursively generated.
        #If a sequence without the forbidden sites cannot be generated after 10 tries it gives up and allows the process to proceed.
        termCod = ["TAG","TAA","TGA","aggaggt","agga","aggt","ggagg","caagg","caag","acaa","gt%s"%(self.myEffector[0:2]), "t%s"%(self.myEffector[0:3])]
        if n == 0:
            return
        for site in self.mysites:
            if self.mySequence.count(site.lower()) > 0:
                self.genSequence()
                return self.sitesCheck(n-1)
            if self.mySequence.count(sequenceGenerator.reverseComplement(site.lower())) > 0:
                self.genSequence()
                return self.sitesCheck(n-1)
        for site in termCod:
            if self.mySequence.count(site.lower()) > 0:
                self.genSequence()
                return self.sitesCheck(n-1)
            if self.mySequence.count(sequenceGenerator.reverseComplement(site.lower())) > 0:
                self.genSequence()
                return self.sitesCheck(n-1)
        return
    def genSequence(self):
        # This method resets the current mySequence and runs the appropriate methods (in the correct order) to generate a sequence.
        self.mySequence = ""
        self.genPromo()
        self.genGene()
        self.genTerm()
    def genPromo(self):
        # This method generates the promoter. If the effector is less than 21 base pairs long the Shine-Dalgarno sequence is put into a tight loop
        # to render it inaccessible. If the effector is too long to fit into the gap near the SD region the a hairpin that begins before the SD
        # region is used that partially uses the SD region as a part of the complement.
        promoter1 = "ttgaca"
        promoter2 = "tataat"
        AA = ["a","t","g","c"]
        startCodon = ["ATG", "GTG", "TTG"]
        self.mySequence = ""
        self.mySequence += sequenceGenerator.reverseComplement(self.myEffector)
        self.mySequence += "ag"
        for x in range(0,6):
            self.mySequence += random.choice(AA)
        self.mySequence += promoter1
        self.mySequence += 'acctcct'
        self.mySequence += self.myEffector[0:4]
        for x in range(0,21-7-6):
            self.mySequence += random.choice(AA)
        self.mySequence += promoter2
        for x in range(0,4):
            self.mySequence += random.choice(AA)
        self.mySequence = self.mySequence.lower()
        self.mySequence += random.choice(startCodon)
            
        self.mySequence = self.mySequence.replace("u","t")
        
    def genGene(self):
        # This method adds the gene to the sequence. A two dimensional array of the most conserved amino acid and its expression probability
        # is used to determine which amino acid is to be expressed using the probability to determine the frequency that amino acid will be
        # expressed in my sequences. If the most frequent amino acid is not selected (randomly based on its p value) then a random mutation
        # (a random amino acid is selected) is introduced.
        AAList = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", 
                "M", "F", "P", "S", "T", "W", "Y", "V"]
        for x in self.myCalcFreq.freqTable[1:]:
            # x[0]=AA, x[1]=p
            if x[1] + random.random() > 1:
                self.mySequence += self.genCodon(x[0])
            else:
                self.mySequence += self.genCodon(random.choice(AAList))
    def genCodon(self, AA):
        # This method generates codons based on the amino acid single letter passed in and the relative frequency of expression based on the codonfreq.txt file.
        AADict = self.myInputFreq.AADict
        codonList = AADict[AA]
        totalFreq = 0
        codonFreq = [[0,0] for x in range(0,len(codonList))]
        n=0
        for codon in codonList:
            totalFreq += float(self.myInputFreq.myDic[codon.lower()])
            codonFreq[n][0] = codon
            n+=1
        n=0
        for codon in codonList:
            codonFreq[n][1]=float(self.myInputFreq.myDic[codon.lower()])/totalFreq
            n+=1
        return self.genCodonHelper(codonFreq).lower()
    def genCodonHelper(self, codonFreq):
        # This recursive helper function helps insure the choices are random and based on the codonfreq.txt relative frequencies.
        rand = random.choice(codonFreq)
        if rand[1] + random.random() > 1:
            return rand[0]
        else:
            return self.genCodonHelper(codonFreq)
    def genTerm(self):
        # This method generates the termination sequence, created semi-random GC rich loops with trailing "a", such that the
        # translated mRNA would have trailing "U".
        termCod = ["TAG","TAA","TGA"]
        AA = ["a","t","g","c"]
        GC = ["g","c"]
        self.mySequence += random.choice(termCod)
        self.mySequence += random.choice(AA)+random.choice(AA)+random.choice(AA)
        GCHairpin = ""
        for x in range(0,5):
            GCHairpin += random.choice(GC)
        for x in range(0,4):
            GCHairpin += random.choice(AA)
        GCHairpin += sequenceGenerator.reverseComplement(GCHairpin[0:5])
        self.mySequence += GCHairpin+random.choice(AA)+random.choice(AA)
        self.mySequence += "aaaaaaaa"
        
        
    def reverseComplement(sequence):
        # This method is used to create DNA reverse complements from DNA and RNA sequences.
        sequence = sequence[::-1]
        DNAcomplement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'R': 'Y',
                         'W':'W', 'Y': 'R', 'S':'S', 'B':'V', 'V':'B', 'H':'D',
                         'X':'N', 'N':'N', '.':'.','U':'A'}
        reverseComplement = ""
        for BP in sequence:
            reverseComplement += DNAcomplement[BP.upper()]
        return reverseComplement.lower()
            
        
