#!/usr/bin/env python
#-*- coding: utf-8 -*-

import time as Time
import sys 
import os
import MDAnalysis
from math import cos, sin, sqrt
import math
import numpy
from numpy import matrix
from numpy import dot 

import sys
sys.path.append('../Coor')
sys.path.append('../mymath')
sys.path.append('../')
import  DNA_matrix
import  atomlib
 
def Get_RMSD_fromTRJ(param_file="", traj_file="", skip=1, dt=1, begin=0, end=-1, output_name=""):
    print "=================================================="
    print "CALCULATE THE THE RMSD OF A SERIES OF TRAJECTORIES"
    print "=================================================="
    print "  ...init..."
    START_TIME=Time.time()

    trajName_list=list() 
    for name in traj_file.split(','):
        #print name.strip()    
        trajName_list.append(name.strip())
    #print trajName_list 
   
    u_list=list()
    for name in trajName_list:
        print "  ...Reading trajectory from %s %s\t" %(param_file, name)
        u=MDAnalysis.Universe(param_file, name) 
        u_list.append(u)
    
    numFrames=0
    trajIndex=0
    for u in u_list:
        trajIndex += 1
        print "  ...%d frames in %s" %(u.trajectory.numframes, trajName_list[trajIndex-1])
        numFrames += u.trajectory.numframes

    print "  ...The total frame is %d\t" %(numFrames)


    '''set output file name'''
    if output_name == "" :
        output_name = "rmsdtraj.dat"
    else :
        pass

    try:
        fp = open(output_name , 'w')
    except:
        print "Error: Error happens in creating %s: " %(output_name)

    atomSELECT="(resname DG5 or resname DG or resname DG3 or resname DC5 or resname DC or resname DC or resname DC3 or \
                resname DA5 or resname DA or resname DA3 or resname DT5 or resname DT or resname DT or resname DT3) \
                and (not name H*)" 
    RMSD=0
    frameIDs=0
    begin +=1
    for u in u_list:
        aaSelect = u.selectAtoms(atomSELECT)
        for ts in u.trajectory:
            frameIDs += 1
            if frameIDs == begin :
                reference_coor = aaSelect.coordinates()
                RMSD=0.0
            else :
                experim_coor = aaSelect.coordinates()       
                fit_coor = DNA_matrix.Get_fit_coor(reference_coor, experim_coor)
                RMSD = DNA_matrix.Get_RMSD(fit_coor, experim_coor)

            if (frameIDs%50) == 0 :
                print  "  ...Frame Nos. %10d  RMSD=%10.5f" %(frameIDs, RMSD)
            else :
                pass

            fp.write("%10d    %10.5f\n" %(frameIDs, RMSD))

    print "  ...end...."
    print "=================================================="

if __name__=="__main__":
    Get_RMSD_fromTRJ("20-at-nuc-neutralized.prmtop", "equip1.dcd")
