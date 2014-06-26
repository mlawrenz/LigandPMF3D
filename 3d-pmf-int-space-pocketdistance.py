from msmbuilder import Trajectory
import PMF3D
import optparse
import pylab
from numpy import * 
import glob
import os
import sys
import pickle
from scipy import interpolate, integrate

def get_minmax(ligcoors, map):
    # build grid of min/max ligand coords from Gens.vmd_ligcoords.dat
    mapped_ligcoors=dict()
    n=0
    for i in xrange(ligcoors.shape[0]):
        if map[i]!=-1:
            mapped_ligcoors[n]=ligcoors[i]
            n+=1
    xmin=100000
    xmax=0
    ymin=100000
    ymax=0
    zmin=100000
    zmax=0
    dx=1
    dy=1
    dz=1
    convert=10 #convert nm to angstrom
    coordinates=['x', 'y', 'z']
    mins=[xmin, ymin, zmin]
    maxes=[xmax, ymax, zmax]
    for i in sorted(mapped_ligcoors.keys()):
        for j in range(0, len(mapped_ligcoors[i])):
            for (n, k) in enumerate(coordinates):
                if mapped_ligcoors[i][j][n]<=mins[n]:
                    mins[n]=mapped_ligcoors[i][j][n]
                elif mapped_ligcoors[i][j][n]>=maxes[n]:
                    maxes[n]=mapped_ligcoors[i][j][n]
    lengths=dict()
    for (n, i) in enumerate(coordinates):
        lengths[n]=int(((round(maxes[n])+1)-round(mins[n]))/1.0)
    box_volume=lengths[0]*convert*dx*lengths[1]*convert*dy*lengths[2]*convert*dz
    print "actual box volume %s" % box_volume
    total=max(lengths.values())
    ranges=dict()
    for (n, i) in enumerate(coordinates):
        pad=-(lengths[n]-total)
        ranges[n]=range(int(round(mins[n])), int(round(maxes[n]+1+(pad*1.0+1))), 1)
    return mapped_ligcoors, ranges[0], ranges[1], ranges[2], box_volume

def estimate_free_energy(space, spacetrack, mapped_com_distances, mapped_ligcoors, volume='all'):
    frees=[]
    corrs=[]
    axis=[]
    volumes=[]
    cutoffs=arange(0, 40, 1)
    if volume=='all':
        for cutoff in cutoffs:
            bound_frames=where(mapped_com_distances < cutoff)[0]
            if len(bound_frames)==0:
                print "no bound states less than reference com distance %s" % cutoff
                continue
            new_points={ key: mapped_ligcoors[key] for key in bound_frames}
            new_pops={ key: space.pops[key] for key in bound_frames}
            GD=space.new_manual_allcoor_grid(spacetrack, new_points, new_pops, type='pops')
            boundspace=space.pmfvolume(GD)

            new_pops={ key: 1.0 for key in bound_frames}
            GD=space.new_manual_allcoor_grid(spacetrack, new_points, new_pops, type='pops')
            boundvolume=space.pmfvolume(GD)
            volumes.append(boundvolume)
            print "count bound volume ", boundvolume

            # unbound frames are above COM cutoff
            unbound_frames=array([int(x) for x in mapped_states if x not in bound_frames])
            new_points={ key: mapped_ligcoors[key] for key in unbound_frames}
            new_pops={ key: space.pops[key] for key in unbound_frames}
            GD=space.new_manual_allcoor_grid(spacetrack, new_points, new_pops, type='pops')
            unboundspace=space.pmfvolume(GD)

            # for counting states
            new_pops={ key: 1.0 for key in unbound_frames}
            GD=space.new_manual_allcoor_grid(spacetrack, new_points, new_pops, type='pops')
            unboundvolume=space.pmfvolume(GD)
            #print "unbound volume ", unboundvolume

            # free energy from ratio
            depth=-0.6*log(boundspace/unboundspace)
            frees.append(depth)
            axis.append(cutoff)
            print "corrected integrated dG ratio at cutoff %s is %s" % (cutoff, depth+correction)
        k=len(space.pops)
        savetxt('%s/frees.dat' % (modeldir), [x+correction for x in frees])
        savetxt('%s/volumes.dat' % (modeldir), volumes)
        savetxt('%s/axis.dat' % (modeldir), axis)
    else:
        k=len(pops)
        standard_frees=loadtxt('%s/frees.dat' % (modeldir))
        bound_volumes=loadtxt('%s/volumes.dat' % (modeldir))
        axis=loadtxt('%s/axis.dat' % modeldir)
        cutoff=float(volume)
        index=where(axis==cutoff)[0]
        print "standard free energy is %s" % (standard_frees[index])
        print "(population-weighted) bound volume is %s A^3" % bound_volumes[index]

def main(modeldir, genfile,volume):
    dir=os.path.dirname(genfile)
    filename=genfile.split(dir)[1].split('.lh5')[0]
    gens=Trajectory.load_from_lhdf(genfile)
    print "computing PMF from model in %s" % modeldir
    map=loadtxt('%s/Mapping.dat' % modeldir)
    pops=loadtxt('%s/Populations.dat' % modeldir)
    indices=dict()
    keys=['prot', 'ligand']
    protfile='../test/AtomIndices-prot.dat'
    ligandfile='../test/AtomIndices-ligand.dat'
    indices['prot']=loadtxt(protfile, dtype=int, ndmin=1)
    indices['ligand']=loadtxt(ligandfile, dtype=int, ndmin=1)

    mapped_ligcoors, x_range, y_range, z_range, box_volume=get_minmax(gens['XYZList'][:, indices['ligand']], map)
    correction=-0.6*log(box_volume/1600.0)
    print "correction %s" % correction
    savetxt('%s/standard_correction.dat' % modeldir, array([correction,]))

    # get prot-lig distances
    if not os.path.exists('%s.com.dat' % genfile.split('.lh5')[0]):
        print "get Gen protein ligand COM distances per state!"
    else:
        print "loading COM per state info"
        com_distances=loadtxt('%s.com.dat' % genfile.split('.lh5')[0], ndmin=1)

    # Generate Mapped 3-D PMF
    frames=where(map!=-1)[0]
    mapped_states=map[frames]
    mapped_com_distances=com_distances[frames]
    space=PMF3D.PMF3D(pops, x_range, y_range, z_range)
    spacetrack=space.microstate_count_grid(mapped_ligcoors)
    new_pops={ key: pops[n] for (n, key) in enumerate(mapped_ligcoors.keys())}
    GD=space.new_manual_allcoor_grid(spacetrack, mapped_ligcoors, new_pops, type='pops')
    free=array([-0.6*log(i) for i in pops])
    subtract=min(free)
    free=array([k-subtract for k in free])
    GD=space.new_manual_allcoor_grid(spacetrack, mapped_ligcoors, new_pops, type='pops')
    GDfree=-0.6*log(GD)
    GDfree=PMF3D.convert(GDfree, max(free))
    GDfree=GDfree-min(GDfree.flatten())
    space.write_dx(GDfree, modeldir)
    #estimate_free_energy(space, spacetrack, mapped_com_distances, mapped_ligcoors, volume='all')


def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-m', '--modeldir', dest='modeldir',
                      help='directory with MSM')
    parser.add_option('-g', '--genfile', dest='genfile',
                      help='input gens file lh5')
    parser.add_option('-v', '--volume', dest='volume',
                          help='volume cutoff')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(modeldir=options.modeldir, genfile=options.genfile, volume=options.volume)

