
#-*- codlin:utf-8 -*-

import os
import subprocess
import time
import shutil
import numpy as np
from scitools.all import *


#M should b a multipel of 24T
Nx = 80; Ny = 80; M = 288; T = 4; c = 1; Lx = 1; Ly = 1; 

t_0 = time.time()
subprocess.call(["waveSquar2DNeuman.x", str(Nx), str(Ny), 
        str(M), str(T), str(Lx), str(Ly)])
cpuTime = time.time()-t_0
print 'cpu time: ', cpuTime;

#Makeing the grid
x = np.linspace(0,Lx,Nx)
y = np.linspace(0,Ly,Ny)
Z = np.zeros((Ny,Nx))

def openAndPlotFile(filename, t):
    f =open(filename, 'r')
    #Filling Z with values from the file
    i = 0;
    for line in f:
        values = line.split()
        for k in xrange(len(values)):
            Z[i][k] = float(values[k]);
        i += 1
    f.close()
    #plotting
    plotfilename = 'plotWave_Nx%g_Ny%g_t%4.2f.png' % (Nx,Ny, t)
    surf(x,y,Z,
         xlabel = 'x',
         ylabel = 'y',
         shading = 'flat',
         title  = 'Wave equation t = %4.2f' % t,
         axis = [0,1,0,1,-0.5,0.5],
         show = False,
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
    #Denne if testen virker mest sannynelig ikke.
    if filename[:16+len(str(Nx))] != 'wave_squar_2D_Nx' + str(N):
        print 'Have removed the file:', filename, ' from the list!'
        listOfFile.remove(filename);
    else:
        t = float(filename[(19+len(str(N))+len(str(M))):-4]);
        #Make plot from datafile
        plotfilename = openAndPlotFile(filename, t);
        #copy plot to subdirctory
        shutil.copy(os.path.join(orgdir, plotfilename) ,plotdir);
        #Delete plot in current directory
        os.remove(os.path.join(orgdir, plotfilename));
        #Delete the datafiles
        #IKKE bruk denne kommandoen før du er sikker på at if testen er korrekt
        #os.remove(os.path.join(orgdir, filename));

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

#pydoc scitools.easyviz.movie
