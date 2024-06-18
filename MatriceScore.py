'''
Created on 18 mars 2020

@author: felix
'''
import numpy as np

class MatriceScore(object):
    '''
        Classe qui défini le type d'une matrice des scores
    
Pour modifier les valeurs de point des matchs, missmatchs et gaps, rendez-vous à la ligne 17.
Remplacer les valeurs de chaque score par les valeurs souhaitées.
    '''

                                                        #Paramètres de point en fonction d'un match, missmatch ou gap à changer ici
    def __init__(self, seq1="Pas_De", seq2="Sequence", scoreMatch = 2, scoreMissmatchIntra = 1, scoreMissmatchExtra = -1, scoreGapOuverture = -10, scoreGapExtensif = -1):
        '''
        Constructor
        '''
        self.scoreMatch = scoreMatch                    # Affecte le score à ajouter pour un match (facilement modifiable, 2 par défaut)
        self.scoreMissmatchIntra = scoreMissmatchIntra  # Affecte le score à ajouter pour un mismatch purine/purine ou pyrimidine/pyrimidine, 1 par défaut
        self.scoreMissmatchExtra = scoreMissmatchExtra  # Affecte le score à ajouter pour un mismacth purine/pyrimidine, -1 par défaut
        self.scoreGapOuverture = scoreGapOuverture      #-10 par défaut
        self.scoreGapExtensif = scoreGapExtensif       #-1 par défaut
        
        self.tailleSeq1 = len(seq1)
        self.tailleSeq2 = len(seq2)
        self.seq1 = seq1
        self.seq2 = seq2
        self.matrice=np.array([[None]*(self.tailleSeq2+1)]*(self.tailleSeq1+1))       # Créer la matrice de bonne dimension
        
    def initMS(self):  
    # Initialise la matrice de score (MS)
        
        # Affecte à chaque case de la seconde ligne une valeur de -10 à -x (x étant la taille de la séquence 2+10)
        for i in range(self.tailleSeq2 + 1):
            self.matrice[0][i] = -i - 10  
            
        # Affecte à chaque case de la seconde colonne une valeur de -10 à y (y étant la taille de la séquence 1+10)
        for i in range(self.tailleSeq1 + 1):
            self.matrice[i][0] = -i - 10     
                 
        self.matrice[0][0] = 0    
        
    def getMatch(self, i, j):
    #Récupère le score diagonal
        score = 0
        if self.seq1[i].lower() == self.seq2[j].lower():  # Test si c'est un match
            score = self.scoreMatch + self.matrice[i][j]
            
        elif (self.seq2[j].lower() == 'a' and self.seq1[i].lower() == 'g') or (self.seq2[j].lower() == 'g' and self.seq1[i].lower() == 'a') or (self.seq2[j].lower() == 'c' and self.seq1[i].lower() == 't') or (self.seq2[j].lower() == 't' and self.seq1[i].lower() == 'c'):  
        # Dans le cas d'un mismatch, test si c'est un mismatch pyrimidine/pyrimidine ou purine/purine    
            score = self.scoreMissmatchIntra + self.matrice[i][j]
        else:  
        # Si on arrive ici, c'est que c'est un mismatch purine/pyrimidine
            score = self.scoreMissmatchExtra + self.matrice[i][j]
        self.scoreDiago = score
        
    def getGapUp(self, MTB, i, j):  # Vérifie si il s'agit d'un Gap d'ouverture ou extensif puis donne le score du GapUp
    #Récupère le score gauche
        score = 0
        if i == 0:  # Si on est sur la première ligne, on en conclu que c'est forcement un gap d'ouverture
            score = self.matrice[1 + i - 1][1 + j] + self.scoreGapOuverture
        else:  # si on est pas sur la première ligne
            if MTB[1 + i - 1][1 + j] == '|':  # Si la case au dessus était déjà un gap, c'est un gap extensif
                score = self.matrice[1 + i - 1][1 + j] + self.scoreGapExtensif
            else:  # sinon c'est un gap d'ouverture
                score = self.matrice[1 + i - 1][1 + j] + self.scoreGapOuverture
        self.scoreHaut = score

    def getGapLeft(self, MTB, i, j):  # idem que isGapUp mais pour le Gap de Gauche
    #Récupère le score haut
        score = 0
        if j == 0:
            score = self.matrice[1 + i][j] + self.scoreGapOuverture
        else:
            if MTB[1 + i][j] == '-':
                score = self.matrice[1 + i][j] + self.scoreGapExtensif
            else:
                score = self.matrice[1 + i][j] + self.scoreGapOuverture
        self.scoreGauche = score
        
    def bestScore(self, MTB, i, j):  
    # Donne le meilleur score entre le score diago, gapUp et GapLeft
        self.getMatch(i, j)  # Affecte la valeur du scoreDiago
        self.getGapUp(MTB, i , j)  # Affecte la valeur du scoreGapUp
        self.getGapLeft(MTB, i , j)  # Affecte la valeur du scoreGapLeft

        if self.scoreDiago > self.scoreHaut and self.scoreDiago > self.scoreGauche:  # Test l'origine du meilleur score : ici, le score Diago est le plus élevé
            self.origineScoreMax = "*"
        elif self.scoreGauche > self.scoreDiago and self.scoreGauche > self.scoreHaut:  # scoreGauche le plus élevé
            self.origineScoreMax = "-"
        elif self.scoreHaut > self.scoreDiago and self.scoreHaut > self.scoreGauche:  # scoreHaut le plus élevé
            self.origineScoreMax = "|"
        elif self.scoreDiago == self.scoreGauche and self.scoreDiago > self.scoreHaut:  # scoreDiago et scoreGauche sont les plus elevés
            self.origineScoreMax = "DG"
        elif self.scoreDiago == self.scoreHaut and self.scoreDiago > self.scoreGauche:  # scoreDiago et scoreHaut sont les plus élevés
            self.origineScoreMax = "DH"
        elif self.scoreGauche == self.scoreHaut and self.scoreGauche > self.scoreDiago:  # scoreHaut et scoreGauche sont les plus élevés
            self.origineScoreMax = "GH"
        elif self.scoreDiago == self.scoreGauche and self.scoreDiago == self.scoreHaut:  # les 3 scores sont égaux
            self.origineScoreMax = "3"
            
        # Récupère la valeur du plus haut score   
        self.scoreMax = max(self.scoreDiago, self.scoreHaut, self.scoreGauche)  

        