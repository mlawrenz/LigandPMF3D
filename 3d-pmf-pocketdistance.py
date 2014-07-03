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

# Helper Functions
def eval_distance(mapped_state_distances, cutoff):
    # assumes mapped states and state_distances
    mapped_cutoff_states=[]
    for state in xrange(len(mapped_state_distances)):
        # if mindist to protein > cutoff, add to unbound frames 
        if mapped_state_distances[state] > cutoff:
            mapped_cutoff_states.append(state)
    mapped_cutoff_states=array([int(i) for i in mapped_cutoff_states])
    return mapped_cutoff_states

def get_ligand_minmax(ligcoors, map):
    # build grid of min/max ligand coords from Gens.vmd_ligcoords.dat
    mapped_ligcoors=dict()
    n=0
    for i in xrange(ligcoors.shape[0]):
        if map[i]!=-1:
            mapped_ligcoors[n]=[k*10 for k in ligcoors[i]]
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
    box_volume=lengths[0]*dx*lengths[1]*lengths[2]*dz
    print "actual box volume %s" % box_volume
    total=max(lengths.values())
    ranges=dict()
    for (n, i) in enumerate(coordinates):
        pad=-(lengths[n]-total)
        ranges[n]=range(int(round(mins[n])), int(round(maxes[n]+1+(pad*1.0+1))), 1)
    return mapped_ligcoors, ranges[0], ranges[1], ranges[2], box_volume


def estimate_com_free_energy(modeldir, space, spacetrack, mapped_states, mapped_com_distances, mapped_ligcoors, correction):
    frees=[]
    volumes=[]
    corrs=[]
    axis=[]
    #cutoffs for Protein-Ligand COM distances
    cutoffs=arange(0, 20, 1)
    if os.path.exists('%s/frees.dat' % (modeldir)):
        print "free energies already computed for %s" % modeldir
        sys.exit()
    for cutoff in cutoffs:
        bound_frames=where(mapped_com_distances < cutoff)[0]
        if len(bound_frames)==0:
            print "no bound states less than reference com distance %s" % cutoff
            continue
        print "on cutoff %s" % cutoff
        new_points={ key: mapped_ligcoors[key] for key in bound_frames}
        new_pops={ key: space.pops[key] for key in bound_frames}
        GD=space.new_manual_allcoor_grid(spacetrack, new_points, new_pops, type='pops')
        boundspace=GD.sum()

        new_pops={ key: 1.0 for key in bound_frames}
        GD=space.new_manual_allcoor_grid(spacetrack, new_points, new_pops, type='pops')
        boundvolume=GD.sum()
        volumes.append(boundspace*boundvolume)
        print "bound volume ", boundspace*boundvolume

        # unbound frames are above COM cutoff
        unbound_frames=array([int(x) for x in mapped_states if x not in bound_frames])
        new_points={ key: mapped_ligcoors[key] for key in unbound_frames}
        new_pops={ key: space.pops[key] for key in unbound_frames}
        GD=space.new_manual_allcoor_grid(spacetrack, new_points, new_pops, type='pops')
        unboundspace=GD.sum()

        # for counting states
        new_pops={ key: 1.0 for key in unbound_frames}
        GD=space.new_manual_allcoor_grid(spacetrack, new_points, new_pops, type='pops')
        unboundvolume=GD.sum()

        # free energy from ratio
        depth=-0.6*log((boundspace*boundvolume)/(unboundspace*boundvolume))
        frees.append(depth)
        axis.append(cutoff)
        print "corrected free energy ratio at cutoff %s is %s" % (cutoff, depth+correction)
    k=len(space.pops)
    savetxt('%s/frees.dat' % (modeldir), [x+correction for x in frees])
    savetxt('%s/boundvolumes.dat' % (modeldir), volumes)
    savetxt('%s/axis.dat' % (modeldir), axis)

def main(modeldir, genfile, ligandfile, writefree=False):
    dir=os.path.dirname(genfile)
    filename=genfile.split(dir)[1].split('.lh5')[0]
    gens=Trajectory.load_from_lhdf(genfile)
    print "computing PMF from model in %s" % modeldir
    map=loadtxt('%s/Mapping.dat' % modeldir)
    pops=loadtxt('%s/Populations.dat' % modeldir)
    ligandind=loadtxt(ligandfile, dtype=int, ndmin=1)

    mapped_ligcoors, x_range, y_range, z_range, box_volume=get_ligand_minmax(gens['XYZList'][:, ligandind], map)
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
    if writefree==True:
        estimate_com_free_energy(modeldir, space, spacetrack, mapped_states, mapped_state_distances, mapped_ligcoors, correction)



def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-m', '--modeldir', dest='modeldir',
                      help='directory with MSM')
    parser.add_option('-g', '--genfile', dest='genfile',
                      help='input gens file lh5')
    parser.add_option('-l', '--ligandfile', dest='ligandfile',
                      help='ligand index file')
    parser.add_option('-f', action="store_true", dest="writefree", help='perform free energy calc with P-L COM distances')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    if options.writefree==True:
        main(modeldir=options.modeldir, genfile=options.genfile, ligandfile=options.ligandfile, writefree=True)
    else:
        main(modeldir=options.modeldir, genfile=options.genfile, ligandfile=options.ligandfile)
