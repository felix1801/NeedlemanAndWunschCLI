'''
Created on 18 mars 2020

@author: felix
'''
import numpy as np

class MatriceTraceback(object):
    '''
        Classe qui défini le type d'une matrice de traceback

Pour modifier les valeurs de point des matchs, missmatchs et gaps, rendez-vous à la ligne 17.
Remplacer les valeurs de chaque score par les valeurs souhaitées.
    '''

                                                        #Paramètres de point en fonction d'un match, missmatch ou gap à changer ici
    def __init__(self, seq1="Pas_De", seq2="Sequence", scoreMatch = 2, scoreMissmatchIntra = 1, scoreMissmatchExtra = -1, scoreGapOuverture = -10, scoreGapExtensif = -1):
        '''
        Constructor
        '''
        self.scoreMatch = scoreMatch
        self.scoreMissmatchIntra = scoreMissmatchIntra
        self.scoreMissmatchExtra = scoreMissmatchExtra
        self.scoreGapOuverture = scoreGapOuverture
        self.scoreGapExtensif = scoreGapExtensif
        
        self.tailleSeq1 = len(seq1)
        self.tailleSeq2 = len(seq2)
        self.seq1 = seq1
        self.seq2 = seq2
        self.matrice = np.array([[None]*(self.tailleSeq2+1)]*(self.tailleSeq1+1))       # Créer la matrice de bonne dimension

    def initMTB(self):  
    # Initialise la matrice des Traceback (MTB)
        for i in range(self.tailleSeq2 + 1):
            self.matrice[0][i] = "-"  # Affecte à chaque case de la seconde ligne le symbole "-" censé représenter une flèche qui va vers la gauche
        for i in range(self.tailleSeq1 + 1):      
            self.matrice[i][0] = "|"    # Affecte à chaque case de la seconde ligne le symbole "|" censé représenter une flèche qui va vers le haut
        self.matrice[0][0] = "Done"  # Affecte à la case de coordonnées 0;0 la valeur "Done"
            
    def aligne(self):  
    # Effectue l'alignement
    
        # initialise i et j pour démarrer dans la case en bas à droite
        i = self.matrice.shape[0] - 1  
        j = self.matrice.shape[1] - 1

        #Créer les variables qui vont stocker les 3 lignes de l'alignement
        self.alignementSeq1 = ""
        self.alignementSeq2 = ""
        self.alignementQuali = ""
        
        #Parcours la matrice et rempli les variables d'alignement
        while self.matrice[i][j] != "Done":  # Parcours la matrice des traceback et s'arrête quand on atteind la case "Done"
            
            if self.matrice[i][j] == '*' or self.matrice[i][j] == "3" or self.matrice[i][j] == "DG" or self.matrice[i][j] == "DH":  
            #si la diagonale est une possibilité de chemin -> 
                self.alignementSeq1 += self.seq1[i-1]  # affecte à l'alignement de la seq 1 la lettre d'indice i de la séquence initiale
                self.alignementSeq2 += self.seq2[j-1]  # affecte à l'alignement de la seq 2 la lettre d'indice j de la séquence initiale
                if self.seq1[i-1].lower() == self.seq2[j-1].lower():  # test si c'est un match ou un missmatch
                    self.alignementQuali += "|"
                else:
                    self.alignementQuali += ":"
                    
                # se déplace dans la case de diagonale Nord-Ouest
                i -= 1  
                j -= 1
            
            if self.matrice[i][j] == '|' or self.matrice[i][j] == "GH":
            #test si le haut est une possibilité (en excluant le haut quand elle est dans le même bloc que la diagonale)
                self.alignementSeq1 += self.seq1[i-1]
                self.alignementSeq2 += '-'
                self.alignementQuali += " "
                
                # se déplace vers la case du haut
                i -= 1  
            
            elif self.matrice[i][j] == "-":
            #test si gauche est la seule possibilité
                self.alignementSeq1 += "-"
                self.alignementSeq2 += self.seq2[j-1]
                self.alignementQuali += " "
                
                # se déplace vers la case de gauche
                j -= 1  
        
        #inverse le sens des séquences alignées pour récupérer le sens initial
        self.alignementSeq1 = self.alignementSeq1[::-1]
        self.alignementQuali = self.alignementQuali[::-1]
        self.alignementSeq2 = self.alignementSeq2[::-1]
        
    def countMatch(self): 
    # compte le nombre de matchs
        nbMatch = 0
        for i in range(len(self.alignementQuali)):
            if self.alignementQuali[i] == "|":
                nbMatch += 1
            self.nbMatch=nbMatch
        
    def countMissmatch(self):  
        # compte le nombre de missmatch
        nbMissmatch = 0
        nbMissmatchIntra = 0
        nbMissmatchExtra = 0
        for i in range(len(self.alignementQuali)):
            if self.alignementQuali[i] == ":":
                nbMissmatch += 1
                if (self.alignementSeq2[i].lower() == 'a' and self.alignementSeq1[i].lower() == 'g') or (self.alignementSeq2[i].lower() == 'g' and self.alignementSeq1[i].lower() == 'a') or (self.alignementSeq2[i].lower() == 'c' and self.alignementSeq1[i].lower() == 't') or (self.alignementSeq2[i].lower() == 't' and self.alignementSeq1[i].lower() == 'c'):
                    nbMissmatchIntra += 1
                else :
                    nbMissmatchExtra += 1
        self.nbMissmatch = nbMissmatch   
        self.nbMissmatchIntra = nbMissmatchIntra
        self.nbMissmatchExtra = nbMissmatchExtra

    def countGap(self):  
        # compte le nombre de gaps
        nbGap = 0
        for i in range(len(self.alignementQuali)):
            if self.alignementQuali[i] == " ":
                nbGap += 1
        self.nbGap = nbGap
    
    def countScore(self):
        #Récuère le score de l'alignement
        score=0
        for i in range(len(self.alignementQuali)):
            if self.alignementQuali[i] == "|":
                score += self.scoreMatch
            elif self.alignementQuali[i] == "-":
                if i == 0:
                    score += self.scoreGapOuverture
                elif self.alignementQuali[i-1] == "-":
                    score += self.scoreGapExtensif
                else:
                    score += self.scoreGapOuverture
            elif self.alignementQuali[i] == ":":
                if (self.alignementSeq2[i].lower() == 'a' and self.alignementSeq1[i].lower() == 'g') or (self.alignementSeq2[i].lower() == 'g' and self.alignementSeq1[i].lower() == 'a') or (self.alignementSeq2[i].lower() == 'c' and self.alignementSeq1[i].lower() == 't') or (self.alignementSeq2[i].lower() == 't' and self.alignementSeq1[i].lower() == 'c'):  # Dans le cas d'un mismatch, test si c'est un mismatch pyrimidine/pyrimidine ou purine/purine
                    score += self.scoreMissmatchIntra
                else:
                    score += self.scoreMissmatchExtra
        self.scoreAlignement = score
        
    def getCount(self):
    #compresse les 4 fonctions précédentes en une
        self.countMatch()
        self.countMissmatch() 
        self.countGap()
        self.countScore()
        
        
        
           