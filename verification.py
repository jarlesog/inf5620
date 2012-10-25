#-*- coding:utf-8 -*-

import os
import sys
import subprocess
import time
import shutil
import numpy as np
from scitools.std import *



#M should b a multipel of 24T



def getError(filename,t,xv,yv,Z, Nx, Ny,Lx,Ly):
        Z = np.loadtxt(filename);
        Z = np.transpose(Z);
        #exact = cos(pi*xv/Lx)*cos(pi*yv/Ly)*exp(-4*t)
        #exact = (1./3*xv-Lx/2.)*xv*xv*(1./3*yv-Ly/2.)*yv*yv*(0.7*t+0.2)
        a = 4.3; c = 2.8;
        exact = cos(pi*xv/Lx)*cos(pi*yv/Ly)*(a*t + c)
        
        error = abs(Z-exact)
        error_avr = sum(error)/((Nx+1)*(Ny+1))
        #print error
        max_error = error.max()
        #return max_error
        return max_error;
def empiricalTest():
        T = 2.32; Lx = 1; Ly = 1;
        errorlist = []
        hlist = [0.2,0.1, 0.05, 0.01,0.005]
        c_max = exp(0.3*Lx*Ly);
        for h in hlist:
                print "h = %g" %h
                #error = []
                Nx = int(Lx/h); Ny = int(Ly/h);
                dx = h; dy = h;
                C_t = 1/(sqrt(1/(dx*dx) + 1/(dy*dy))*h*c_max)*0.95

                M = int(T/(C_t *h))+1

                x = np.linspace(0,Lx,Nx+1)
                y = np.linspace(0,Ly,Ny+1)
                xv, yv = ndgrid(x,y)
                Z = np.zeros((Ny+1,Nx+1))

                subprocess.call(["./waveSquar2DNeuman.x", str(Nx), str(Ny), str(M), str(T), str(Lx), str(Ly)])
                #Gives a list og the files in the current directory
                orgdir = os.getcwd();
                listOfFile = os.listdir(orgdir);
                listOfFile.sort()

                fname = "last.dat"
                errorlist.append(getError(fname, T, xv, yv, Z, Nx, Ny, Lx, Ly)/(h*h))
                os.remove(os.path.join(orgdir, fname));
        plot(errorlist,hardcopy = "error.png")



def handcalctest():
        al = 4
        u0 = [[1,0,-1],[0,0,0],[-1,0,1]]
        tem = (al**2 + 2*pi**2)/16. - al/4.
        u1 = [[tem,0,-tem],[0,0,0],[-tem,0,tem]]
        tem = (1+exp(-al/4))*(al**2 + 2*pi**2)/16. - al/4.-1
        u2 = [[tem,0,-tem],[0,0,0],[-tem,0,tem]]
        Nx = 2; Ny = 2; M = 2; T = 2/4.;Lx = 1; Ly = 1;
        subprocess.call(["./waveSquar2DNeuman.x", str(Nx), str(Ny), str(M), str(T), str(Lx), str(Ly)])
        printmat(u0)
        print "-"*20
        printmat(u1) 
        print "-"*20
        printmat(u2) 
        print "-"*20

def printmat(a):
        for i in a:
                print i

#handcalctest();
empiricalTest();




"""
                #Remove the files that are not datafiles
                for filename in listOfFile:
                #Remove all files that are not data files.
	                if filename[:19+len(str(Nx))+len(str(Nx))] != 'wave_squar_2D_Nx' + str(Nx) + '_Ny' + str(Ny):
                                #print 'ignorering:', filename
                                pass
                                #listOfFile.remove(filename);
                        else:
                                #print 'Useing datafile:', filename;
                                ##Potensial feil telling
                                t = float(filename[(23+len(str(Nx))+len(str(Ny))+len(str(M))):-4]);
                                #Make plot from datafile
                                #error.append(getError(filename, t, xv, yv, Z, Nx, Ny, Lx, Ly)/(h*h))
        


                                os.remove(os.path.join(orgdir, filename));
"""
                #errorlist.append(error);
                #error.append(getError(filename, t, xv, yv, Z, Nx, Ny, Lx, Ly)/(h*h))
"""
        l = []
        for i in range(len(errorlist)):
                error = errorlist[i];
                l.append(error[-1]);
"""


        

"""
#Gives a list og the files in the current directory
orgdir = os.getcwd();
listOfFile = os.listdir(orgdir);


#Remove the files that are not datafiles
for filename in listOfFile:
	#Remove all files that are not data files.
	if filename[:19+len(str(Nx))+len(str(Nx))] != 'wave_squar_2D_Nx' + str(Nx) + '_Ny' + str(Ny):
		#print 'ignorering:', filename
		pass
		#listOfFile.remove(filename);
	else:
		#print 'Useing datafile:', filename;
		##Potensial feil telling
		t = float(filename[(23+len(str(Nx))+len(str(Ny))+len(str(M))):-4]);
		#Make plot from datafile
		verificate(filename, t);


		os.remove(os.path.join(orgdir, filename));


a = sorted(errordict.keys())
for i in a:
        print "error t = %4.2f: %g" %(i, errordict[i])
"""
