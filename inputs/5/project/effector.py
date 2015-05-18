import FASTAhandle
# Creates an object to hold the imported effector
class effector:
    def __init__(self, myFileName):
        self.FASTAobject = FASTAhandle.FASTAhandle(myFileName)
        self.myEffector = self.FASTAobject.myGenes[0].mySequence
