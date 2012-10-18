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
        exact = (1/3.*xv - Lx/2.)*xv*xv*(1/3.*yv-Ly/2.)*yv*yv*(0.7*t + 0.2)
        
        error = abs(Z-exact)
        #print error
        max_error = error.max()
        return max_error

        """
        surf(xv,yv,error,
         xlabel = 'x',
         ylabel = 'y',
         shading = 'flat',
	 #clevels=15,
	 #clabels='on',
	 colorbar='on',
	 #caxis = [0,1/36.],
	 #view = [-1,1],
         title  = 'Wave equation t = %4.2f' % t,
         axis = [0,Lx,0,Ly,-1,1],
         show = False,
	 #rstride=4, #extra
	 #cstride=4, #extra
         hardcopy= "error_t%4.2f.png"%t)
         """
def empiricalTest():
        T = 1.; Lx = 1; Ly = 1;
        errorlist = []
        hlist = [0.2,0.1, 0.05, 0.01, 0.005]
        for h in hlist:
                print "h = %g" %h
                error = []
                Nx = int(1./h); Ny = int(1./h);
                dx = h; dy = h;
                C_t = 1/(sqrt(1/(dx*dx) + 1/(dy*dy))*h)*0.95

                M = int(T/(C_t *h)) + 1

                x = np.linspace(0,Lx,Nx+1)
                y = np.linspace(0,Ly,Ny+1)
                xv, yv = ndgrid(x,y)
                Z = np.zeros((Ny+1,Nx+1))

                subprocess.call(["./waveSquar2DNeuman.x", str(Nx), str(Ny), str(M), str(T), str(Lx), str(Ly)])
                #Gives a list og the files in the current directory
                orgdir = os.getcwd();
                listOfFile = os.listdir(orgdir);
                listOfFile.sort()


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
                                error.append(getError(filename, t, xv, yv, Z, Nx, Ny, Lx, Ly)/(h*h))
        


                                os.remove(os.path.join(orgdir, filename));

                errorlist.append(error);

        l = []
        for i in range(len(errorlist)):
                error = errorlist[i];
                l.append(error[-1]);

        plot(l,hardcopy = "error.png")


empiricalTest();
        

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
