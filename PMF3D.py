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
def eval_distance(states, state_distances, map, cutoff):
    # assumes mapped states and state_distances
    mapped_cutoff_states=[]
    gen_cutoff_states=[]
    for state in states:
        if state_distances[state] < cutoff:
            location=numpy.where(map==state)[0]
            mapped_cutoff_states.append(state)
            gen_cutoff_states.append(location)
    mapped_cutoff_states=numpy.array([int(i) for i in mapped_cutoff_states])
    gen_cutoff_states=numpy.array([int(i) for i in gen_cutoff_states])
    return gen_cutoff_states, mapped_cutoff_states


def parallel_get_distance(input):
    (ligand_coors, protein_coors)=input
    minval=100000
    for ligand_coor in ligand_coors:
        for coor in protein_coors:
            val=numpy.sqrt(((ligand_coor[0]-coor[0])**2+(ligand_coor[1]-coor[1])**2+(ligand_coor[2]-coor[2])**2))
            if val < minval:
                minval=val
    return minval

def get_distance(ligand_input, protein_input, output):
    ligand_handle=open(ligand_input)
    protein_handle=open(protein_input)
    ligand_coors=ligand_handle.readlines()
    protein_coors=protein_handle.readlines()
    states=[]
    state_distances=open(output, 'w')
    for line in ligand_coors:
        n=int(line.split()[0])
        n=n-1
        print "on frame %s" % n
        #if (n % 100)==0:
        #    print "on frame %s" % n
        x=float(line.split()[1])
        y=float(line.split()[2])
        z=float(line.split()[3])
        vals=[]
        for pline in protein_coors:
            m=int(pline.split()[0])
            m=m-1
            if m==n:
                a=float(pline.split()[1])
                b=float(pline.split()[2])
                c=float(pline.split()[3])
                val=numpy.sqrt(((x-a)**2+(y-b)**2+(z-c)**2))
                vals.append(val)
            else:
                continue
        state_distances.write('%s\t%s\n' % (n, min(vals)))

def get_ref_distance(input, output, x_ref, y_ref, z_ref):
    fhandle=open(input)
    states=[]
    state_distances=open(output, 'w')
    for line in fhandle.readlines():
        n=int(line.split()[0])
        n=n-1
        if (n % 100)==0:
            print "on frame %s" % n
        x=float(line.split()[1])
        y=float(line.split()[2])
        z=float(line.split()[3])
        for (a, b, c) in zip(x_ref, y_ref, z_ref):
            val=numpy.sqrt(((x-a)**2+(y-b)**2+(z-c)**2))
            state_distances.write('%s\t%s\n' % (n, val))


def convert(GDfree, max):
    frames=numpy.where(GDfree==numpy.inf)
    for (i,j,k) in zip(frames[0], frames[1], frames[2]):
        GDfree[i,j,k]=max
    return GDfree

def get_coors(fhandle, map):
    allcoors=dict()
    for line in fhandle.readlines():
        n=int(line.split()[0])
        n=n-1
        if map[n]!=-1:
            state=map[n]
            if state not in allcoors.keys():
                allcoors[state]=[]
            x=float(line.split()[1])
            y=float(line.split()[2])
            z=float(line.split()[3])
            allcoors[state].append((x,y,z))
        else:
            continue
    return allcoors

def meshgrid3d(*arrs):
    arrs = tuple(reversed(arrs))  #edit
    lens = map(len, arrs)
    dim = len(arrs)

    sz = 1
    for s in lens:
        sz*=s

    ans = []
    for i, arr in enumerate(arrs):
        slc = [1]*dim
        slc[i] = lens[i]
        arr2 = numpy.asarray(arr).reshape(slc)
        for j, sz in enumerate(lens):
            if j!=i:
                arr2 = arr2.repeat(sz, axis=j)
        ans.append(arr2)

    return tuple(ans)

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

 
    def manual_grid(self, points, values, type='pops'):
        if type=='pops':
            newgrid=numpy.zeros((len(self.xaxis), len(self.yaxis), len(self.zaxis)))
        elif type=='free':
            # pad with max free energy value (pop=0)
            newgrid=max(values)*numpy.ones((len(self.xaxis), len(self.yaxis), len(self.zaxis)))
        for (n, coord) in enumerate(points):
            x_location=numpy.where((self.xaxis>coord[0]-self.tol)&(self.xaxis<coord[0]+self.tol))[0]
            y_location=numpy.where((self.yaxis>coord[1]-self.tol)&(self.yaxis<coord[1]+self.tol))[0]
            z_location=numpy.where((self.zaxis>coord[2]-self.tol)&(self.zaxis<coord[2]+self.tol))[0]
            newgrid[x_location, y_location, z_location]=values[n]
        return newgrid

    def hist_microstates(self, allcoor, values):
        histtrack=dict()
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
                key=(x_location[0], y_location[0], z_location[0])
                if key not in spacetrack[i]:
                    spacetrack[i].append(key)
                else:
                    pass
                if key not in histtrack.keys():
                    histtrack[key]=[]
                    histtrack[key].append(i)
                else:
                    histtrack[key].append(i)
        return histtrack, spacetrack

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

 

    def inter_microstate_count_grid(self, allcoors, pops):
        spacetrack=dict()
        for state in allcoors.keys():
            size=len(allcoors[state])
            xyz=numpy.zeros((3, size))
            for i, x in enumerate(allcoors[state]):
                xyz[0][i]=x[0]
                xyz[1][i]=x[1]
                xyz[2][i]=x[2]
            spacetrack[state]=0
            statex=range(int(round(min(xyz[0]))), int(round(max(xyz[0]))+self.dx), self.dx)
            statey=range(int(round(min(xyz[1]))), int(round(max(xyz[1]))+self.dy), self.dy)
            statez=range(int(round(min(xyz[2]))), int(round(max(xyz[2]))+self.dz), self.dz)
            #check out, order 
            mesh=meshgrid3d(statez, statey, statex)
            points=allcoors[state]
            values=[pops[state]]*len(points)
            GD=interpolate.griddata(points, values, mesh, method='linear', fill_value=0.0)
            oldx=statex[0]
            for i in range(0, len(statex)):
                oldy=statey[0]
                for j in range(0, len(statey)):
                    oldz=statez[0]
                    for k in range(0, len(statez)):
                        if (abs(statex[i]-oldx)*abs(statey[j]-oldy)*abs(statez[k]-oldz)) !=0:
                            if GD[i,j, k] != 0:
                                spacetrack[state]+=1
                        oldz=statez[k]
                    oldy=statey[j]
                oldx=statex[i]
        return spacetrack
    
    def manual2d(self, GD, axes='xy'):
        return 

    def interpolated_grid(self, points, values, type='pops'):
        mesh=meshgrid3d(self.xaxis, self.yaxis, self.zaxis)
        if type=='pops':
            GD=interpolate.griddata(points, values, mesh, method='linear', fill_value=0.0)
        elif type=='free':
            GD=interpolate.griddata(points, values, mesh, method='linear', fill_value=max(values))
        return GD

    def new_inter_microstate_count_grid(self, allcoor, pops):
        spacetrack=dict()
        mesh=meshgrid3d(self.xaxis, self.yaxis, self.zaxis)
        for (n, i) in enumerate(sorted(allcoor.keys())):
            size=len(allcoor[i])
            GD=interpolate.griddata(allcoor[i], [pops[i]]*size, mesh, method='linear', fill_value=0.0) 
            spacetrack[i]=[]
            for j in range(0, len(allcoor[i])):
                x_location=numpy.where((self.xaxis>=allcoor[i][j][0])&(self.xaxis<allcoor[i][j][0]+self.dx))[0]
                y_location=numpy.where((self.yaxis>=allcoor[i][j][1])&(self.yaxis<allcoor[i][j][1]+self.dy))[0]
                z_location=numpy.where((self.zaxis>=allcoor[i][j][2])&(self.zaxis<allcoor[i][j][2]+self.dz))[0]
                if len(x_location)==0 or len(y_location)==0 or len(z_location)==0:
                    print "no bin available for allcoor[i][j]"
                    import pdb
                    pdb.set_trace()
                if [x_location, y_location, z_location] not in spacetrack[i]:
                    spacetrack[i].append([x_location[0], y_location[0], z_location[0]])
        return spacetrack

    def interpolated_grid(self, points, values, type='pops'):
        mesh=meshgrid3d(self.xaxis, self.yaxis, self.zaxis)
        if type=='pops':
            GD=interpolate.griddata(points, values, mesh, method='linear', fill_value=0.0) 
        elif type=='free':
            GD=interpolate.griddata(points, values, mesh, method='linear', fill_value=max(values)) 
        return GD

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

    def func(self, z,y,x):
        i=numpy.where((self.xaxis>=x)&(self.xaxis<(x+self.dx)))[0]
        j=numpy.where((self.yaxis>=y)&(self.yaxis<(y+self.dy)))[0]
        k=numpy.where((self.zaxis>=z)&(self.zaxis<(z+self.dz)))[0]
        return float(self.grid[i,j,k].transpose())

    def pmfvolume2(self):
        (int, err)=scipy.integrate.tplquad(self.func, self.xaxis[0], self.xaxis[-1], self.gfun, self.hfun, self.qfun, self.rfun)
        return int, err

    def pmfvolume3(self, GD):
        sum=0
        oldval=0
        volume=0
        for k in range(0, len(self.zaxis)):
            for j in range(0, len(self.yaxis)):
                for i in range(0, len(self.xaxis)):
                    volume+=GD[i,j, k]*self.dz*self.dx*self.dy
        return float(volume)

    def pmfvolume_new(self, GD):
        start=0
        sum=0
        track=0
        for i in range(0, len(self.xaxis)):
            for j in range(0, len(self.yaxis)):
                for k in range(0, len(self.zaxis)):
                    if GD[i,j, k] != 0:
                        sum+=GD[i,j, k]
        return float(sum) 


    def pmfvolume(self, GD):
        sum=0
        oldx=self.xaxis[0]
        oldval=0
        for i in range(0, len(self.xaxis)):
            oldy=self.yaxis[0]
            for j in range(0, len(self.yaxis)):
                oldz=self.zaxis[0]
                for k in range(0, len(self.zaxis)):
                    #if (abs(self.xaxis[i]-oldx)*abs(self.yaxis[j]-oldy)*abs(self.zaxis[k]-oldz)) !=0:
                    if abs(self.zaxis[k]-oldz) !=0:
                        if GD[i,j, k] != 0:
                            sum+=GD[i,j, k]*abs(self.xaxis[i]-oldx)*abs(self.yaxis[j]-oldy)*abs(self.zaxis[k]-oldz)
                    oldz=self.zaxis[k]
                oldy=self.yaxis[j]
            oldx=self.xaxis[i]
        return float(sum)
 
    def manual2d(self, GD, axes='xy'):
        if axes=='xy':
            iaxis=self.xaxis
            jaxis=self.yaxis
        elif axes=='yz':
            iaxis=self.yaxis
            jaxis=self.zaxis
        elif axes=='xz':
            iaxis=self.xaxis
            jaxis=self.zaxis
        sum=0
        oldi=0
        oldi=0
        oldval=0
        for i in range(0, len(iaxis)):
            for j in range(0, len(jaxis)):
                sum+=(GD[i,j])*abs(iaxis[i]-oldi)*abs(iaxis[j]-oldj)
                oldj=jaxis[j]
                #oldval=grid[i,j]
            oldi=iaxis[i]
        return float(sum/2)
    
