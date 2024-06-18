#!/usr/bin/python3
# encoding: utf-8
'''
MortasFelixNwGUI -- shortdesc

MortasFelixNwGUI is a description

It defines classes_and_methods

@author:     user_name

@copyright:  2020 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''
from LecteurFasta import LecteurFasta
from MatriceScore import MatriceScore
from MatriceTraceback import MatriceTraceback
from appJar import gui


# the title of the button will be received as a parameter
def press():
    fic1 = app.getEntry("f1")
    fic2 = app.getEntry("f2")
    scoreMatch = 2
    scoreMissmatchIntra = 1
    scoreMissmatchExtra = -1
    scoreGapOuverture = -10
    scoreGapExtensif = -1
    
    #Lecture du fichier de données
    lecteur = LecteurFasta(fic1, fic2) #Création de l'objet "lecteur" de type "MortasFelixLecteurFasta" prenant en charge les fichiers .fasta
    lecteur.lit_seq1()   #Lecture du fichier1 .fasta
    lecteur.lit_seq2()   #Lecture du fichier1 .fasta
    print("Séquence 1: " + lecteur.sequence1)
    print("Séquence 2: " + lecteur.sequence2)
    
    #Création des matrices
    ms=MatriceScore(lecteur.sequence1, lecteur.sequence2, scoreMatch, scoreMissmatchIntra, scoreMissmatchExtra, scoreGapOuverture, scoreGapExtensif)        #Construction de la matrice des scores
    ms.initMS()              #Initialisation de la matrice en fonction de la taille des 2 séquences
                                                            #, pointMatch, pointMissmatchIntra, pointMissmatchExtra, pointGapOuverture, pointGapExtensif
    mtb=MatriceTraceback(lecteur.sequence1, lecteur.sequence2, scoreMatch, scoreMissmatchIntra, scoreMissmatchExtra, scoreGapOuverture, scoreGapExtensif)   #Construction de la matrice des traceback
    mtb.initMTB()            #Initialisation de la matrice des traceback en fonction de la taille des 2 séquences
    
    #Remplissage des 2 matrices
    for i in range(ms.tailleSeq1):
        for j in range(ms.tailleSeq2):
            ms.bestScore(mtb.matrice, i, j) #Calcul du score maximum entre celui du gapGauche, gapHaut et diagonal dans la case d'indice [i][j]
            ms.matrice[1+i][1+j] = ms.scoreMax  #Affectation du scoreMax à la case d'indice [i][j] dans la matrice des scores
            mtb.matrice[1+i][1+j] = ms.origineScoreMax #Affectation de l'origine du score max à la case d'indice [i][j] dans la matrice des traceback
    print("\nMatrice des scores :")
    print(ms.matrice)
    app.clearMessage("scores")
    app.setMessage("scores", ms.matrice)
    print("\nMatrice des traceback :")
    print(mtb.matrice)
    app.clearMessage("traceback")
    app.setMessage("traceback", mtb.matrice)

    #Alignement des séquences
    mtb.aligne()  # réalise l'alignement et le stock dans 3 variables
    print("\nAlignement :\n"+mtb.alignementSeq1)
    print(mtb.alignementQuali)
    print(mtb.alignementSeq2)
    app.clearMessage("synthese")
    alignement = "\nAlignement :\n"+mtb.alignementSeq1 + "\n" + mtb.alignementQuali  + "\n" + mtb.alignementSeq1
    app.setMessage("synthese", alignement)
    #Comptage du score et du nombre d'apparition de chaque modalité
    mtb.getCount()  #Effectuer les 4 comptages (matchs, missmatchs, gaps et score de l'alignement
    print("\nNombre de match : "+str(mtb.nbMatch))
    print("Nombre de missmatch : "+str(mtb.nbMissmatch))
    print("dont\tpurine/purine ou pyrimidine/pyrimidine : "+str(mtb.nbMissmatchIntra))
    print("\tpurine/pyrimidine : "+str(mtb.nbMissmatchExtra))
    print("Nombre de gap : "+str(mtb.nbGap))
    print("Score : "+str(mtb.scoreAlignement)+"/"+str(2*len(mtb.alignementSeq1)))


# create a GUI variable called app
app = gui()

# add & configure widgets - widgets get a name, to help referencing them later
app.addLabel("title", "Application NW")
app.setLabelBg("title", "yellow")
app.startFrame("LEFT", row=0, column=0)
app.setBg("green")
app.setSticky("NEW")
app.setStretch("COLUMN")
app.addLabel("Fichier1")
app.addFileEntry("f1")
app.addLabel("Fichier2")
app.addFileEntry("f2")
app.addLabel("Synthèse :")
app.addEmptyMessage("synthese")
app.stopFrame()

app.startFrame("RIGHT", row=0, column=1)
app.addLabel("Matrice des scores :")
app.addEmptyMessage("scores")

app.addLabel("Matrice des traceback :")
app.addEmptyMessage("traceback")
app.stopFrame()

app.startFrame("DOWN", row=1, column=0)
app.addButton("Go !", press)
app.stopFrame()
# start the GUI
app.go()

