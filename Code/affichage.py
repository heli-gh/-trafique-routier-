#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import data

from mpl_toolkits import mplot3d

from pylab import *        


def afficher(list_x,list_y):
    plt.plot(list_x,list_y)
    
    plt.show()
    
def afficherConv(list_x,list_y):
    for j,absic in enumerate(list_x):
        plt.plot(absic,list_y[j], 'x-',linewidth=2.0)
    
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(["log10(err) vs log10(1/dt)"])
    plt.show() 


def afficherDensite(times,xx):
    
    posRel = xx[0:len(times),1:]-xx[0:len(times),0:-1]
    absc = np.arange(0,len(posRel[0,:]),1)
    
    densite = 1.0/posRel
    
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
    times, absc = np.meshgrid(times, absc)
    ax.plot_surface(times,absc,densite.T,cmap='jet')
    plt.show()
    
        
def afficherTrafic(times,list_traj,vv):
    ite = 0
    
    theta = np.arange(0.0,2*np.pi,0.01)
    rb = 0.5
    xb = rb * np.cos(theta)
    yb = rb * np.sin(theta)
    sq = np.array([-2,-2])
    
    plt.ion()
    #plt.hold('on')
    dt = times[1]-times[0]
    T = times[-1]
    k = int(0.10*T*0.1/dt)
    if (k==0):
        k=1
          
    iteEnrg = 0
    distMax = list_traj[-1,-1]
    for ite,t in enumerate(times):
        if (ite % k == 0):
            plt.clf()
            
            plt.plot(1000+xb*700,2.5+yb*1,color=(1,0,0))
            plt.text(840,2.35,str(int(data.Vmax(1000)) ),fontsize=10)
            
            plt.plot([2000,2000],[-5,5],color=(0,0,0))
            
            plt.plot(3000+xb*700,2.5+yb*1,color=(1,0,0))
            plt.text(2840,2.35,str(int(data.Vmax(3000))))
            
            plt.plot([4000,4000],[-5,5],color=(0,0,0))
            
            plt.plot(4500+xb*700,2.5+yb*1,color=(1,0,0))
            plt.text(4340,2.35,str(int(data.Vmax(4500))))
            
            plt.plot([5000,5000],[-5,5],color=(0,0,0))
            
            plt.plot(7000+xb*700,2.5+yb*1,color=(1,0,0))
            plt.text(6740,2.35,str(int(data.Vmax(7000))))
            
            
            
            for j,tr in enumerate(list_traj[ite]):
                if(j<len(list_traj[ite])-1):
                    dist = list_traj[ite][j+1]-tr
                    ratio = dist/20.0
                    plt.plot(tr+xb,yb,color=(max(0.0,1.0-min(1.0,ratio)),0,0))
                    if(j == 0):
                        plt.text(tr-400.0,1.0,str(np.around(vv[ite][j])))
                else:
                    plt.plot(tr+xb,yb,color=(0,0,0))
                    plt.text(tr+0.5,-1.0,str(np.around(vv[ite][j])))
            
            plt.ylim([-5,5])
            plt.xlim([0,10000])                    
            plt.show()
            # Sauvegarde de l'image pour creer ensuite une animation gif.
            #plt.savefig('Anim/Img'+str('%03d' %(iteEnrg+1))+'.png')
            iteEnrg+=1
            plt.pause(1e-10)
            #ite = ite+k
    
    plt.ioff()