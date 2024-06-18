
#Created on 18 mars 2020

#@author: felix

from LecteurFasta import LecteurFasta
from MatriceScore import MatriceScore
from MatriceTraceback import MatriceTraceback

#Pour modifier les valeurs de point des matchs, missmatchs et gaps, rendez-vous à la ligne 26 et 29.
#Pour chaque ligne, copiez le contenu du commentaire de la ligne du dessus (sans le #)
# et le coller juste avant la fermeture de parenthèse de la ligne où l'on se trouvait.
#Remplacer les textes 'pointXxxx' par les valeurs souhaitées.

if __name__ == '__main__':
    
    #Lecture du fichier de données
    lecteur=LecteurFasta(cheminFichier="fichier.fasta") #Création de l'objet "lecteur" de type "LecteurFasta"
    lecteur.lit_seq()   #Lecture du fichier .fasta
    print("Séquence 1: "+lecteur.sequence1)
    print("Séquence 2: "+lecteur.sequence2)
    
    #Création des matrices
    ms=MatriceScore(lecteur.sequence1, lecteur.sequence2)        #Construction de la matrice des scores
    ms.initMS()              #Initialisation de la matrice en fonction de la taille des 2 séquences
                                                            #, pointMatch, pointMissmatchIntra, pointMissmatchExtra, pointGapOuverture, pointGapExtensif
    mtb=MatriceTraceback(lecteur.sequence1, lecteur.sequence2)   #Construction de la matrice des traceback
    mtb.initMTB()            #Initialisation de la matrice des traceback en fonction de la taille des 2 séquences
    
    #Remplissage des 2 matrices
    for i in range(ms.tailleSeq1):
        for j in range(ms.tailleSeq2):
            ms.bestScore(mtb.matrice, i, j) #Calcul du score maximum entre celui du gapGauche, gapHaut et diagonal dans la case d'indice [i][j]
            ms.matrice[1+i][1+j] = ms.scoreMax  #Affectation du scoreMax à la case d'indice [i][j] dans la matrice des scores
            mtb.matrice[1+i][1+j] = ms.origineScoreMax #Affectation de l'origine du score max à la case d'indice [i][j] dans la matrice des traceback
    print("\nMatrice des scores :")
    print(ms.matrice)
    print("\nMatrice des traceback :")
    print(mtb.matrice)

    #Alignement des séquences
    mtb.aligne()  # réalise l'alignement et le stock dans 3 variables
    print("\nAlignement :\n"+mtb.alignementSeq1)
    print(mtb.alignementQuali)
    print(mtb.alignementSeq2)
    
    #Comptage du score et du nombre d'apparition de chaque modalité
    mtb.getCount()  #Effectuer les 4 comptages (matchs, missmatchs, gaps et score de l'alignement
    print("\nNombre de match : "+str(mtb.nbMatch))
    print("Nombre de missmatch : "+str(mtb.nbMissmatch))
    print("dont\tpurine/purine ou pyrimidine/pyrimidine : "+str(mtb.nbMissmatchIntra))
    print("\tpurine/pyrimidine : "+str(mtb.nbMissmatchExtra))
    print("Nombre de gap : "+str(mtb.nbGap))
    print("Score : "+str(mtb.scoreAlignement)+"/"+str(2*len(mtb.alignementSeq1)))
    
    
    
    
