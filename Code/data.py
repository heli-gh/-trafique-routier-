#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
import copy


##----------------------------------
## Fonction vitesse max autorisee à la position x
def Vmax(x):
    # Decommenter pour analyse convergence :
    # return x+1
    # Ne pas oublier de recommenter si etude des modeles
    
    # Pour modele ordre 1 ou 2 :
    
    Vm = 90.0
    if (x <= 2000): 
        Vm = 90.0
    if (x <= 4000 and x >= 2000): 
        Vm = 90#75#90.0- 1/50*(x-2000)  
    if(x > 4000 and x < 5000):
        Vm = 50.0
    if(x > 5000):
        Vm = 110.0
    return Vm
    
##----------------------------------    
    
    
##----------------------------------    
## Vitesse d'un automobiliste: (modele ordre 1)
def V(x,dist):
    VM = Vmax(x)
    Dmin = 4.0
    
    if(dist <= Dmin):
        # Completer ICI :
        return 0.0
        # return ...
        
    if(dist <= 3*Dmin):
        # Completer ICI :
        return VM*((dist-Dmin)**2)*((5*Dmin -dist )**2)/(4*Dmin**2)**2
        # return ...
        
    if(dist >= 3*Dmin):
        # Completer ICI :
        return VM
        # return ...
##----------------------------------    
 
##----------------------------------
# Vitesses automobilistes :
def fmodelO1(t,X):
    res = X *0.0
    for j,x in enumerate(X):#why there is a x
        if(j < len(X)-1):
            # Completer ICI :
            res[j]=V(X[j],X[j+1]-X[j])
            #res[j] = ...
        if(j == len(X)-1):
            res[j] = V(X[j],1e3) # 1e3 <=> tres grand
        
    return res
##----------------------------------
        
        
##----------------------------------            
## Acceleration d'un automobiliste:
def Acc(x,dist,vit,vitRel):
    VM = Vmax(x)
    Dmin = 20.0
    if(vitRel >= 0 and vit < VM and dist >= Dmin):
        # Completer ICI :
        return 12*(VM- vit)/VM
        # return ...                
    if(vitRel <= 0 and dist <= 2*Dmin and vit >= 0):
        # Completer ICI :
        return -1*(2*Dmin-dist)*vit/(dist**2)
        #return ...
    if(dist < Dmin and vit >= 0):
        # Completer ICI :
        return -(Dmin-dist)*vit/(dist**2)
        #return ...
    if(vitRel >= 0 and vit >= VM):
        # Completer ICI :
        return -12*(vit-VM)/VM
        #return ...
    return 0.0
##----------------------------------    

##----------------------------------
# Vitesses automobilistes :
def fmodelO2(t,X):
    res = X *0.0
    nbVt = int(len(X)/2)
    res[0:nbVt] = X[nbVt:]
    for j in range(nbVt):
        if(j < nbVt-1):
            res[j+nbVt] = Acc(X[j],X[j+1]-X[j],X[j+nbVt],X[j+1+nbVt]-X[j+nbVt])
        if(j == nbVt-1):
            res[j+nbVt] = Acc(X[j],1e3,X[j+nbVt],1e3) #Acctete(t,X[j+nbVt]) #Acc(20.0,X[j],0.0)
        
    return res
##----------------------------------







