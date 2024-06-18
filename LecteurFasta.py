'''
Created on 18 mars 2020

@author: felix
'''

class LecteurFasta(object):
    '''
    Classe qui définit un objet permettant de lire un fichier .fasta
    '''


    def __init__(self, cheminFichier):
        '''
        Constructor
        '''
        self.cheminFichier=cheminFichier
        
    
    def lit_seq(self):  
    # Lit un fichier Fasta et affecte les 2 premières séquences à 2 variables
        seqFile=open(self.cheminFichier,"r")
        sNuc = ""
        ligne = seqFile.readline()
        if ligne == None:
            return
        ligne = seqFile.readline()
        sNuc += ligne.strip()
        while seqFile.readline(1) != ">" and ligne != "":
            ligne = seqFile.readline()
            sNuc += ligne.strip()
        self.sequence1 = sNuc
        sNuc=""
        seqFile.readline()
        while ligne.startswith(">") == False and ligne != "":
            ligne = seqFile.readline()
            sNuc += ligne.strip()
        self.sequence2 = sNuc
        seqFile.close()

    
    