#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
import copy


##----------------------------------
#  Methode de Euler explicite :
def Euler(dt,Times,y0,f):
    sol = [y0];
    p0 = y0*1
    p1 = y0*0
    for t in Times:
        # Completer ICI :i
        p1=p0 + dt*f(t,p0)
        #p1 = ...
        
        sol.append(p1)
        p0 = p1*1.0
    sol.pop()
    return sol
##----------------------------------    


##----------------------------------
#  Methode de RK :
def RK(dt,Times,y0,f):
    sol = [y0];
    p0 = y0*1
    p1 = y0*0
    # Modifier ICI :
    a = 1/3
    b = 2/3
    c1 = 1/4
    c3 = 3/4
    
    for t in Times:
        # Completer ICI :
        #...
        k1=f(t,p0)
        k2=f(t+a*dt,p0+dt*a*k1)
        k3=f(t+b*dt,p0+dt*b*k2)
        p1=p0+dt*(c1*k1+c3*k3)
        sol.append(p1)
        p0 = p1*1.0
    sol.pop()
    return sol
##----------------------------------
    
##----------------------------------
#  Methode de Multi-pas :
    

def MultiPas(dt,Times,y0,f):
    sol = [y0];
    
    # Initialisation :
    # Completer ICI :
    p0=y0*1
    # p1=p0+dt*f(Times[0],p0)  #Euler
    s1=RK(dt,np.array([0,dt]),y0,f)
    p1=s1[1]

    # p2=p1+dt*f(Times[1],p1)  #Euler
    s2=RK(dt,np.array([dt,2*dt]),p1,f)
    p2 = s2[1]
    sol.append(p1)
    sol.append(p2)
    p3=y0*0
    
    for t in range(2,len(Times)):
        # Completer ICI :
        p3=p2+(dt/12)*(23*f(Times[t],p2)-16*f(Times[t-1],p1)+5*f(Times[t-2],p0))
        sol.append(p3)
        p0=p1*1.0
        p1=p2*1.0
        p2=p3*1.0
        
    sol.pop()
    return sol