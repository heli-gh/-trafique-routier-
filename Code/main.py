#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
import affichage
import matplotlib.pyplot as plt
import algo
import data

import time



"""
    Etude de modèles pour le trafic routier
"""

##----------------------------------
# But du programme :
# 1. Etude du modele d'ordre 1.
# 2. Analyse convergence.
# 3. Etude du modele d'ordre 2
butPrgm = 1
print("Que souhaitez vous faire ?")
print(" 1. Modèle trafic ordre 1.")
print(" 2. Analyse de la méthode numérique.")
print(" 3. Modèle trafic ordre 2.")
print(" 4.Multi-pas last car position")
butPrgm = input("Entrer le choix : ")

print("")
# Choix de la methode :
print("Choisir la méthode :")
print(" 1. Euler")
print(" 2. RK")
print(" 3. Multi-Pas")

choixMeth = input("Entrer votre choix : ")

NbVoitures = 1
if(butPrgm == "1" or butPrgm == "3" ):
    print("")
    print("Choisir le nombre de voiture :")
    NbVoitures = eval(input("Entre votre choix : "))
##----------------------------------


##----------------------------------
# Partie I : Modele ordre 1
if (butPrgm=="1"):
    
    espaceVt = 20.0
    x0 = np.arange(0,(NbVoitures)*espaceVt,espaceVt)
    y0 = x0 * 1.0
    
    ##----------------------------------
    # Paramètres 
    # Temps simulation
    T = 100.0
    # Discretisation
    dt = 0.01
    Times = np.arange(0.0,T,dt)
    ##----------------------------------
    
    ##----------------------------------
    # Calcul solution approchée :
    SOL = []
    start = time.time()
    if (choixMeth == "1"):
        SOL = algo.Euler(dt,Times,y0,data.fmodelO1)
    if (choixMeth == "2"):
        SOL = algo.RK(dt,Times,y0,data.fmodelO1)
    if (choixMeth == "3"):
        SOL = algo.MultiPas(dt,Times,y0,data.fmodelO1)
    print("Temps calculs : ",time.time() - start,"s") 
    ##----------------------------------
    
    # Recuperation des positions :
    xx = [e[0:NbVoitures] for e in SOL]    
    # Calculs des vitesses :
    vv = []
    for iteT,t in enumerate(Times):
        vv.append([])
        for i,x in enumerate(SOL[iteT]):
            if(i < len(SOL[iteT])-1):
                vv[iteT].append(data.V(x,SOL[iteT][i+1]-x ) )
            else:
                vv[iteT].append(data.V(x,1e3) )
                
    ##----------------------------------
    # Affichage des résultats 
    # Positions
    affichage.afficher(Times,xx)
    # Vitesses
    affichage.afficher(Times,vv)
    xx = np.array(xx)
    #affichage.afficherDensite(Times,xx)
    affichage.afficherTrafic(Times,xx,vv)
    ##----------------------------------
##----------------------------------


##----------------------------------
# Partie II :
if (butPrgm=="2"):
    list_dt = [0.1,0.05,0.00025,0.0125,0.00625]
    list_err = []
    
    T = 3.0
    for dt in list_dt:
        Times = np.arange(0.0,T,dt)

        x0 = np.array([0.0])
        SOL = []
        
        ##----------------------------------
        # Calcul solution approchée :
        start = time.time()
        if (choixMeth == "1"):
            SOL = algo.Euler(dt,Times,x0,data.fmodelO1)
        if (choixMeth == "2"):
            SOL = algo.RK(dt,Times,x0,data.fmodelO1)
        if (choixMeth == "3"):
            SOL = algo.MultiPas(dt,Times,x0,data.fmodelO1)
        print("Temps calculs : ",time.time() - start,"s")
        ##----------------------------------

        ##----------------------------------
        # Calcul de la solution exacte pour la fonction Vmax(x) = x+1
        # Modifier ICI :
        SOL_e = np.exp(Times)-1
        ##----------------------------------

        ##----------------------------------
        # Analyse erreur :
        maxerr = max([np.sqrt(np.dot(e-SOL_e[i], e-SOL_e[i])) for i,e in enumerate(SOL)])

        print("dt = ",dt," et err =",maxerr)
        list_err.append(maxerr)
        ##----------------------------------

    ## Tracer courbe erreur :
    poly = np.polyfit(np.log10(list_dt)*(-1),np.log10(list_err),1)
    print("")
    print("Droite approchee : y=",poly[0],"x - ",np.abs(poly[1]))
    affichage.afficherConv([1.0/np.array(list_dt)], [list_err])
##----------------------------------



##----------------------------------
# Partie III : Modele ordre 2
if (butPrgm=="3"):
    # Espace inter-voitures
    espaceVt = 20.0
    # Initialisation position voitures:
    x0 = np.arange(0,(NbVoitures)*espaceVt,espaceVt)
    # Vitesse initiale voiture :
    v0 = x0*0.0+90.0
    # Vecteur contenant les positions et les vitesse :
    y0 = np.concatenate((x0,v0))
    
    ##----------------------------------
    # Paramètres 
    # Temps simulation
    T = 100.0
    # Discretisation
    dt = 0.01
    Times = np.arange(0.0,T,dt)
    ##----------------------------------
    
    ##----------------------------------
    # Calcul solution approchée :
    SOL = []
    start = time.time()
    if (choixMeth == "1"):
        SOL = algo.Euler(dt,Times,y0,data.fmodelO2)
    if (choixMeth == "2"):
        SOL = algo.RK(dt,Times,y0,data.fmodelO2)
    if (choixMeth == "3"):
        SOL = algo.MultiPas(dt,Times,y0,data.fmodelO2)
    print("Temps calculs : ",time.time() - start,"s") 



    ##----------------------------------        
    xx = [e[0:NbVoitures] for e in SOL]    


    vv = [e[NbVoitures:] for e in SOL]
    ##----------------------------------
    # Affichage des résultats 
    # Positions
    affichage.afficher(Times,xx)
    # Vitesses
    affichage.afficher(Times,vv)
    xx = np.array(xx)
    # Densité
    affichage.afficherDensite(Times,xx)
    

    # Animation trafic :
    affichage.afficherTrafic(Times,xx,vv)
        ##------
##----------------------------------
    ##----------------------------------
        ##----------------------------------
    # pour calculer le distance de dernier voiture
if(butPrgm =="4"):
    DV=[]
    cars=[10,20,30,40,50,60,70,80,90,100]
    espaceVt = 20.0    
    # Temps simulation
    T = 100.0
    # Discretisation
    dt = 0.01
    Times = np.arange(0.0,T,dt)
    Times = np.arange(0.0,T,dt)


    for k in cars :

        x0 = np.arange(0,k*espaceVt,espaceVt)
        # Vitesse initiale voiture :
        v0 = x0*0.0+90.0
        # Vecteur contenant les positions et les vitesse :
        y0 = np.concatenate((x0,v0))

        SOL = []
        
        ##----------------------------------
        # Calcul solution approchée :
        
        if (choixMeth == "1"):
            SOL = algo.Euler(dt,Times,y0,data.fmodelO2)
        if (choixMeth == "2"):
            SOL = algo.RK(dt,Times,y0,data.fmodelO2)
        if (choixMeth == "3"):
            SOL = algo.MultiPas(dt,Times,y0,data.fmodelO2)
        xx = [e[0:NbVoitures] for e in SOL]  
        c=xx[len(SOL)-1]

        DV.append(c[0])
        
    print(DV)
    plt.plot(cars,DV)
    plt.show() 
print("FIN \n")
########################################