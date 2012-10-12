#-*- coding:utf-8 -*-

import os
import sys
import subprocess
import time
import shutil
import numpy as np
from scitools.all import *



#M should b a multipel of 24T
Nx = 15; Ny = 15; M = 100; T = 1; c = 1; Lx = 1; Ly = 1; 

dx = Lx/float(Nx); dy = Ly/float(Ny); dt = T/float(M);

t_0 = time.time()
subprocess.call(["./waveSquar2DNeuman.x", str(Nx), str(Ny), str(M), str(T), str(Lx), str(Ly)])
cpuTime = time.time()-t_0
print 'C++ program done!!! cpu time: ', cpuTime;

#Makeing the grid
x = np.linspace(0,Lx,Nx+1)
y = np.linspace(0,Ly,Ny+1)
xv, yv = ndgrid(x,y)
Z = np.zeros((Ny+1,Nx+1))

def openAndPlotFile(filename, t):
    #Load the datafil into a matrix
    Z = np.loadtxt(filename);
    Z = np.transpose(Z);


    #plotting
    plotfilename = 'plotWave_Nx%g_Ny%g_t%4.2f.png' % (Nx,Ny, t)
    surf(xv,yv,Z,
         xlabel = 'x',
         ylabel = 'y',
         shading = 'flat',
	 #clevels=15,
	 #clabels='on',
	 colorbar=[-1,1],
	 #view = [-1,1],
         title  = 'Wave equation t = %4.2f' % t,
         axis = [0,Lx,0,Ly,-1,1],
         show = False,
	 #rstride=4, #extra
	 #cstride=4, #extra
         hardcopy= plotfilename)
    return plotfilename


#Gives a list og the files in the current directory
orgdir = os.getcwd();
listOfFile = os.listdir(orgdir);



#Make a subdirector named plots
plotdir = os.path.join(orgdir, "plots");
try:
    #Remove old plots
    shutil.rmtree(plotdir)
    print "Old plot dirctory is removed"
except:
    pass
    #Continue
os.mkdir(plotdir);


#Remove the files that are not datafiles
for filename in listOfFile:
	#Remove all files that are not data files.
	if filename[:19+len(str(Nx))+len(str(Nx))] != 'wave_squar_2D_Nx' + str(Nx) + '_Ny' + str(Ny):
		print 'ignorering:', filename
		pass
		#listOfFile.remove(filename);
	else:
		print 'Useing datafile:', filename;
		##Potensial feil telling
		t = float(filename[(23+len(str(Nx))+len(str(Ny))+len(str(M))):-4]);
		#Make plot from datafile
		plotfilename = openAndPlotFile(filename, t);
		#copy plot to subdirctory
		shutil.copy(os.path.join(orgdir, plotfilename) ,plotdir);
		#Delete plot in current directory
		os.remove(os.path.join(orgdir, plotfilename));
		#Delete the datafiles
		##IKKE bruk denne kommandoen fÃžr du er sikker pÃ¥ at if testen er korrekt
		os.remove(os.path.join(orgdir, filename));

#change directory to where the plots are
os.chdir(plotdir)
#make a movie

##print "start time"
##t_temp = time.time()
##while time.time()-t_temp< 5:
##    pass
##print "end time!"

movie('plotWave_*.png', fps = 12, quiet = True)
#
print 'program time: ', time.time()-t_0;
