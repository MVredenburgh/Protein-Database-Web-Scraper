import sys
from urllib.request import urlopen, urlretrieve
import urllib

# BioE 131 Project 1, Michael Vredenburgh


class myCards:
    # This object holds all the infomation for a card.
    myCards = None
    def __init__(self):
        # Object constructor: Creates HTML files with the cards ready to print.
        self.myCards = []
    def makeCards(self):
        ANGenerator = webParser()
        while len(self.myCards) < 130:
            try:
                self.myCards.append(card(ANGenerator.makeNewAN()))
            except:
                pass
    def makeAllPages(self):
        for card in self.myCards:
            try:
                if len(card.name) < 26:
                    myCards.makeWebPage(card.length, card.weight, card.stiochiometry, card.name, card.organism, card.picture)
            except:
                pass
    def makeWebPage(length, weight, stiochiometry, name, organism, picture):
        # Creates the .html and .css files from which to print from.
        html_str = """
<!DOCTYPE HTML PUBLIC>
<html lang="en">
<head>
    <title><!-- Vredenburgh Bio131 Cards --></title>
    <style>
        h1 {font-size: 40%%}
        p {line-height: 5%%}
    </style>
</head>
<body>
    <h1>%s</h1>
    <img src="%s" alt="4KAA" height="150" width="150">
    <h1><p>Organism: %s</p>
    <p>Length: %s</p>
    <p>Weight: %s</p>
    <p>Stoichiometry: %s></p></h1>
    <body background="E:\\Dropbox\\BioE131\\CardBack.jpg" bgcolor="orange"  no-repeat>
</body>
</html>    
        """ %(name,picture,organism,length,weight,stiochiometry)
        html_file = open("webpages\%s.html" % name,"w")
        html_file.write(html_str)
        html_file.close()

class card:
    length = None
    weight = None
    stiochiometry = None
    name = None
    organism = None
    picture = None
    # This page object holds all the card objects
    def __init__(self, accessionNumber):
        # Card object that stores all related data together.
        webFetcher = webParser()
        myPage = webpage(accessionNumber)
        self.length = webFetcher.getlength(myPage)
        self.stiochiometry = webFetcher.getStioch(myPage)
        self.name = webFetcher.getName(myPage)
        self.organism = webFetcher.getOrganism(myPage)
        self.picture = webFetcher.getPic(accessionNumber)
        self.weight = webFetcher.getWeight(myPage)

class webParser:
    # Parser function takes infomation from the biological database www.rcsb.org/pdb.
    # The database is accessed through accession numbers (Eg. "4KYC").
    numMade= None
    letters = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    c1=0
    c2=0
    c3=0
    c4=2
    def __init__(self):
        self.numMade = 0
        
    def getlength(self, myPage):
        # Find "Length:" then the first "<td" to the first ">" then take the string untill the first "<".
        marker = myPage.htmlArray.index(b'class="se_key">Length:</span></td>')
        containsLength = myPage.htmlArray[marker +2]
        length = [int(s) for s in containsLength[17:20].split() if s.isdigit()]
        return str(length[0])
    
    def getWeight(self, myPage):
        marker = myPage.htmlArray.index(b'Weight:</span>')
        containsWeight = myPage.htmlArray[marker +4]
        weight = containsWeight.decode("utf-8")
        return weight
    
    def getStioch(self, myPage):
        # Find "Stoichiometry:" take next string untill "<".
        # Known issues: Not all webpages brought up have stioch info, May cut off larger strings.
        marker = myPage.htmlArray.index(b'symmetry<br/>Stoichiometry:')
        containsStoich = myPage.htmlArray[marker +1]
        return containsStoich[3:10].decode("utf-8")
    
    def getName(self, myPage):
        # Find keyword "Molecule:"
        marker = myPage.htmlArray.index(b'class="se_key">Molecule:</span></td>')
        start = marker +4
        n=start
        myName = []
        name = ""
        while myPage.htmlArray[n]!= b'</td>':
            myName.append(myPage.htmlArray[n])
            n+=1

        for word in myName:
            name += word.decode("utf-8")+" "
        return name
    
    def getPic(self, accessionNumber):
        # Near the text "Biological Assemble". /pdb/images/4kzc_bio_r_500.jpg
        url = "http://www.rcsb.org/pdb/images/" + accessionNumber + "_asym_r_250.jpg"
        image = urlretrieve(url, "E:\\Dropbox\\BioE131\\images\\%s.jpg" % accessionNumber)
        return "E:\\Dropbox\\BioE131\\images\\%s.jpg" % accessionNumber
    
    def getOrganism(self, myPage):
        marker = myPage.htmlArray.index(b'Organism')
        start = marker+9
        n=start
        myOrganism = []
        organism = ""
        while myPage.htmlArray[n]!= b'</a>':
            myOrganism.append(myPage.htmlArray[n])
            n+=1
        for word in myOrganism:
            organism += word.decode("utf-8")+" "
        # chop off excess from front and back
        organism = organism[18:]
        organism = organism[:len(organism)-20]
        return organism
    
    def makeNewAN(self):
        self.c1+=1
        if self.c1 >= 26:
            self.c1 += -26
            self.c2 += 1
        if self.c2 >= 26:
            self.c2 += -26
            self.c3 += 1
        if self.c3 >= 26:
            self.c3 += -26
            self.c4 += 1
        return "%s%s%s%s" %(self.c4, self.letters[self.c3], self.letters[self.c2], self.letters[self.c1])

class webpage:
    html = None
    htmlArray = None
    def __init__(self, accessionNumber):
        downloadedData = urllib.request.urlopen('http://www.rcsb.org/pdb/explore/explore.do?structureId=%s' % accessionNumber)
        self.html = downloadedData.read()
        self.htmlArray = self.html.split()

run = myCards()
run.makeCards()
run.makeAllPages()
