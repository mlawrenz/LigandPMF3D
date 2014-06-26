from msmbuilder import Trajectory
import multiprocessing
import optparse
import pylab
import numpy
import glob
import os
from numpy import linalg
from scipy import interpolate
from scipy import integrate
import scipy.spatial

# Helper Functions
def convert(GDfree, max):
    frames=numpy.where(GDfree==numpy.inf)
    for (i,j,k) in zip(frames[0], frames[1], frames[2]):
        GDfree[i,j,k]=max
    return GDfree

# Class for 3D PMF
class PMF3D:
    def __init__(self, pops=None, xaxis=None, yaxis=None, zaxis=None, grid=None):
        self.pops=pops
        self.dx =abs(xaxis[0]-xaxis[1])
        self.dy =abs(yaxis[0]-yaxis[1])
        if zaxis!=None:
            self.dz =abs(zaxis[0]-zaxis[1])
        self.xaxis = numpy.array(xaxis)
        self.yaxis = numpy.array(yaxis)
        if zaxis!=None:
            self.zaxis = numpy.array(zaxis)
        if self.dx==3:
            self.tol=1.5
        elif self.dx==2:
            self.tol=1.0
        elif self.dx==1:
            self.tol=0.5
        # for quadrature integration
        self.a=xaxis[0]
        self.b=yaxis[-1]
        self.gfun=lambda x: yaxis[0]
        self.hfun=lambda x: yaxis[-1]
        if zaxis!=None:
            self.qfun=lambda x, y: zaxis[0]
            self.rfun=lambda x, y: zaxis[-1]
        if grid!=None:
            self.grid=grid

 
    def microstate_count_grid(self, allcoor):
        spacetrack=dict()
        for i in sorted(allcoor.keys()):
            if i not in spacetrack.keys():
                spacetrack[i]=[]
            for j in range(0, len(allcoor[i])):
                x_location=numpy.where((self.xaxis>=round(allcoor[i][j][0]))&(self.xaxis<(round(allcoor[i][j][0])+self.dx)))[0]
                y_location=numpy.where((self.yaxis>=round(allcoor[i][j][1]))&(self.yaxis<(round(allcoor[i][j][1])+self.dy)))[0]
                z_location=numpy.where((self.zaxis>=round(allcoor[i][j][2]))&(self.zaxis<(round(allcoor[i][j][2])+self.dz)))[0]
                if len(x_location)==0 or len(y_location)==0 or len(z_location)==0:
                    print "no bin available for allcoor[i][j]"
                    import pdb
                    pdb.set_trace()
                if [x_location, y_location, z_location] not in spacetrack[i]:
                    spacetrack[i].append([x_location[0], y_location[0], z_location[0]])
        return spacetrack

    def new_manual_allcoor_grid(self, spacetrack, allcoor, values, type='pops'):
        if type=='pops':
            newgrid=numpy.zeros((len(self.xaxis), len(self.yaxis), len(self.zaxis)))
        elif type=='free':
            # pad with max free energy value (pop=0)
            newgrid=max(values)*numpy.ones((len(self.xaxis), len(self.yaxis), len(self.zaxis)))
        for (n, i) in enumerate(sorted(allcoor.keys())):
            cn=len(spacetrack[i])
            for j in range(0, len(allcoor[i])):
                x_location=numpy.where((self.xaxis>=round(allcoor[i][j][0]))&(self.xaxis<(round(allcoor[i][j][0])+self.dx)))[0]
                y_location=numpy.where((self.yaxis>=round(allcoor[i][j][1]))&(self.yaxis<(round(allcoor[i][j][1])+self.dy)))[0]
                z_location=numpy.where((self.zaxis>=round(allcoor[i][j][2]))&(self.zaxis<(round(allcoor[i][j][2])+self.dz)))[0]
                if len(x_location)==0 or len(y_location)==0 or len(z_location)==0:
                    print "no bin available for allcoor[i][j]"
                    import pdb
                    pdb.set_trace()
                newgrid[x_location, y_location, z_location]+=(1.0/cn)*values[i]
        return newgrid

 

    def write_dx(self, GD, dir):
        gridspace=1
        newfile=open('%s/pmf3d_d%s.dx' % (dir, gridspace), 'w')
        newfile.write('# Data calculated MSM State Probabilities\n')
        newfile.write('object 1 class gridpositions counts %s %s %s\n' % (GD.shape[0], GD.shape[1], GD.shape[2]))
        newfile.write('origin %s %s %s\n' % (self.xaxis[0], self.yaxis[0], self.zaxis[0]))
        newfile.write('delta %s 0 0\n' % self.dx)
        newfile.write('delta 0 %s 0\n' % self.dy)
        newfile.write('delta 0 0 %s\n' % self.dz)
        newfile.write('object 2 class gridconnections counts %s %s %s\n' % (GD.shape[0], GD.shape[1], GD.shape[2]))
        newfile.write('object 3 class array type double rank 0 items %s data follows\n' % (GD.shape[0]*GD.shape[1]*GD.shape[2]))
        intergrid=numpy.zeros((GD.shape[0], GD.shape[1], GD.shape[2]))
        count=0
        for i in range(0, GD.shape[0]):
            for j in range(0, GD.shape[1]):
                for k in range(0, GD.shape[2]):
                    if count==2:
                        if GD[i][j][k]==0:
                            newfile.write('%s\n' % int(GD[i][j][k]))
                        else:
                            newfile.write('%s\n' % GD[i][j][k])
                        count=0
                    else:
                        if GD[i][j][k]==0:
                            newfile.write('%s\t' % int(GD[i][j][k]))
                        else:
                            newfile.write('%s\t' % GD[i][j][k])
                        count+=1
        newfile.write('\nobject "ligand free energy" class field')
        newfile.close()

    def pmfvolume(self, GD):
        sum=0
        oldx=self.xaxis[0]
        oldval=0
        for i in range(0, len(self.xaxis)):
            oldy=self.yaxis[0]
            for j in range(0, len(self.yaxis)):
                oldz=self.zaxis[0]
                for k in range(0, len(self.zaxis)):
                    if abs(self.zaxis[k]-oldz) !=0:
                        if GD[i,j, k] != 0:
                            sum+=GD[i,j, k]*abs(self.xaxis[i]-oldx)*abs(self.yaxis[j]-oldy)*abs(self.zaxis[k]-oldz)
                    oldz=self.zaxis[k]
                oldy=self.yaxis[j]
            oldx=self.xaxis[i]
        return float(sum)
 
