import FASTAhandle
# Creates a list of sites that cannot occur in the final sequence.
class sites:
    def __init__(self, myFileName):
        self.FASTAobject = FASTAhandle.FASTAhandle(myFileName)
        self.makeList()
    def makeList(self):
        self.mySitesList = []
        for gene in self.FASTAobject.myGenes:
            self.mySitesList.append(gene.mySequence)
