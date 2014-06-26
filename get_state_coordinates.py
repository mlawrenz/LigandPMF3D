#!/bin/python
import numpy
from msmbuilder import Trajectory
import pickle
import optparse


def coors_to_dict(genfile, ligandfile, protfile):
    gens=Trajectory.load_from_lhdf(genfile)
    indices=dict()
    keys=['prot', 'ligand']
    indices['prot']=numpy.loadtxt(protfile, dtype=int, ndmin=1)
    indices['ligand']=numpy.loadtxt(ligandfile, dtype=int, ndmin=1)

    data=dict()
    for key in keys:
        data[key]=dict()
        for state in xrange(gens['XYZList'].shape[0]):
            if state not in data[key].keys():
                print "on state %s" % state
                data[key][state]=[]
            for coor in gens['XYZList'][state, indices[key]]:
                data[key][state].append(coor)

    ofile=open('%s.coors.pickle' % genfile.split('.lh5')[0], 'w')
    pickle.dump(data, ofile)
    ofile.close()
    return data

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-g', '--genfile', dest='genfile',
                      help='genfile')
    parser.add_option('-l', '--ligandfile', dest='ligandfile',
                      help='ligand indices')
    parser.add_option('-p', '--protfile', dest='protfile',
                      help='protein indices')
    (options, args) = parser.parse_args()
    return (options, args)



if __name__ == "__main__":
    (options, args) = parse_commandline()
    coors_to_dict(genfile=options.genfile, ligandfile=options.ligandfile, protfile=options.protfile)

