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
from numpy import array 
import shutil

import sys
sys.path.append('../Coor')
sys.path.append('../mymath')
sys.path.append('../')
import  DNA_matrix
import  atomlib
 
def Trj_PBC_ION_DNA_RDF(param_file="", traj_file="", skip=1, dt=1, begin=0, end=-1, ionName='Na+',output_name="ion-dna.txt", radiusRDF=35.0, radiusDelta=0.35):

    print "=================================================="
    print "READ THE STRAND INFORMATION"
    print "=================================================="

    f = open("dnastrand.in", 'r')
    nl=0
    for line in f.readlines():
        line=line.strip('\n')
        nl=nl+1

        if nl==1:
            strand1=line.split(":")[1].split()
        elif nl==2:
            strand2=line.split(":")[1].split()
        else:
            pass
    f.close()

    print "  >>>>DNA strand>>>>"
    print "  5'-->3':  ", strand1
    print "  3'<--5':  ", strand2
    print "=================================================="


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
        output_name = "ion-dna.txt"
    else :
        pass

    try:
        fp = open(output_name , 'w')
    except:
        print "Error: Error happens in creating %s: " %(output_name)
    fp.close()

#    atomSELECT="(resname DG5 or resname DG or resname DG3 or resname DC5 or resname DC or resname DC or resname DC3 or \
#                resname DA5 or resname DA or resname DA3 or resname DT5 or resname DT or resname DT or resname DT3) \
#               and (not name H*)" 
#    RMSD=0
#    frameIDs=0
#    begin +=1
#    for u in u_list:
#        aaSelect = u.selectAtoms(atomSELECT)
#        for ts in u.trajectory:
#            frameIDs += 1
#            if frameIDs == begin :
#                reference_coor = aaSelect.coordinates()
#                RMSD=0.0
#            else :
#                experim_coor = aaSelect.coordinates()       
#                fit_coor = DNA_matrix.Get_fit_coor(reference_coor, experim_coor)
#                RMSD = DNA_matrix.Get_RMSD(fit_coor, experim_coor)
#
#            if (frameIDs%50) == 0 :
#                print  "  ...Frame Nos. %10d  RMSD=%10.5f" %(frameIDs, RMSD)
#            else :
#                pass
#
#            fp.write("%10d    %10.5f\n" %(frameIDs, RMSD))
#
#    print "  ...end...."
    print "=================================================="


    
    print "  ...init..."
    #parametes for RDF
    #radiusRDF=35.0
    #radiusDelta=0.35
    radiusN=int(radiusRDF/radiusDelta)
    arrayRDF=numpy.zeros(radiusN, numpy.int)

    dnaLength=0.0

    ionNumber=0
    boxVolumn=0.0

    #Atoms on bases used for fitting
    BASE_AG_LIST=["N9","C8","N7","C5","C6","N1","C2", "N3","C4"]
    BASE_CTU_LIST=["N1","C2","N3","C4","C5","C6"]

    if len(strand1) != len(strand2):
        print "  ...>>>Error: Different bases number in strand1 and strand2!!"
        exit()
    else:
        pass

    #renew strands:only the central 1/2 part of DNA is used for calcualtions
    resStart=len(strand1)/4
    resEnd=len(strand1)*3/4

    renewstrand1=list([strand1[0]])+strand1[resStart:resEnd]+list([strand1[-1]])
    renewstrand2=list([strand2[0]])+strand2[resStart:resEnd]+list([strand2[-1]])
    print "  ...strands renewed:"
    print "  ...strand 1:", renewstrand1
    print "  ...strand 2:", renewstrand2
    print "  ...RUNNING TRAJECTORY"

    frameIDs=0
    for u in u_list:
        for ts in u.trajectory[begin:end]:
            frameIDs +=1

            # if (frameIDs%50) == 1 :
                # print  "  ...Frame Nos. %10d  " %(frameIDs)
            # else :
                # pass
            print  "  ...Frame Nos. %10d  " %(frameIDs)

            #bo information
            boxVolumn+=ts.volume

            Rsum=numpy.zeros((1,3))
            end5center=numpy.zeros((3))
            end3center=numpy.zeros((3))
            dnacenter=numpy.zeros((3))


            ##For strand 1
            #for res in strand1:
            # renewstrand1=list([strand1[0]])+strand1[resStart:resEnd]+list([strand1[-1]])
            #print renewstrand1
            for res in renewstrand1:

                experim_coor=list()
                #print res, len(experim_coor)
                
                atomSELECT="resid" + ' ' + res
                
                residue=u.selectAtoms(atomSELECT)
                #print res, residue.resids(), residue.resnames()
                if residue.resnames()[0] in ['DA5', 'DA3', 'DA']:
                    atomSelectList=BASE_AG_LIST
                    base_name='A'
                elif residue.resnames()[0] in ['DG5', 'DG3', 'DG']:
                    atomSelectList=BASE_AG_LIST
                    base_name='G'
                elif residue.resnames()[0] in ['DT5', 'DT3', 'DT']:
                    atomSelectList=BASE_CTU_LIST
                    base_name='T'
                elif residue.resnames()[0] in ['DC5', 'DC3', 'DC']:
                    atomSelectList=BASE_CTU_LIST
                    base_name='C'
                else:
                    print ">>>Error: Do not understand the Base Name!!"
              

                for atomName in atomSelectList:
                    atomSELECT="resid" + ' ' + res + " and " + "name " + atomName
                    #print atomSELECT
                    atomu=u.selectAtoms(atomSELECT)
                    # print type(atomu.coordinates())
                    # print atomu.coordinates()[0]
                    # sys.exit()
                    experim_coor.append( atomu.coordinates()[0])
                    #print type(experim_coor)
                    tmp_coor=array(experim_coor)
                # print type(tmp_coor), tmp_coor.shape, type(base_name)
                # print tmp_coor

                R=DNA_matrix.Get_base_rotate_matrix(tmp_coor, base_name)

                if res == strand1[0]:
                    end5center+=R[1]
                elif res == strand1[len(strand1)-1]:
                    end3center+=R[1]
                else:
                    dnacenter+=R[1]


                Rv3=R[0][:,2]
                Rsum=Rsum+Rv3

            #For Strand 2
            #for res in strand2:
            # renewstrand2=list([strand2[0]])+strand2[resStart:resEnd]+list([strand2[-1]])
            for res in renewstrand2:

                experim_coor=list()
                #print res, len(experim_coor)
                
                atomSELECT="resid" + ' ' + res
                
                residue=u.selectAtoms(atomSELECT)
                #print res, residue.resids(), residue.resnames()
                if residue.resnames()[0] in ['DA5', 'DA3', 'DA']:
                    atomSelectList=BASE_AG_LIST
                    base_name='A'
                elif residue.resnames()[0] in ['DG5', 'DG3', 'DG']:
                    atomSelectList=BASE_AG_LIST
                    base_name='G'
                elif residue.resnames()[0] in ['DT5', 'DT3', 'DT']:
                    atomSelectList=BASE_CTU_LIST
                    base_name='T'
                elif residue.resnames()[0] in ['DC5', 'DC3', 'DC']:
                    atomSelectList=BASE_CTU_LIST
                    base_name='C'
                else:
                    print ">>>Error: Do not understand the Base Name!!"
              

                for atomName in atomSelectList:
                    atomSELECT="resid" + ' ' + res + " and " + "name " + atomName
                    #print atomSELECT
                    atomu=u.selectAtoms(atomSELECT)
                    # print type(atomu.coordinates())
                    # print atomu.coordinates()[0]
                    # sys.exit()
                    experim_coor.append( atomu.coordinates()[0])
                    #print type(experim_coor)
                    tmp_coor=array(experim_coor)
                # print type(tmp_coor), tmp_coor.shape, type(base_name)
                # print tmp_coor

                R=DNA_matrix.Get_base_rotate_matrix(tmp_coor, base_name)

                if res == strand2[0]:
                    end5center+=R[1]
                elif res == strand2[len(strand2)-1]:
                    end3center+=R[1]
                else:
                    dnacenter+=R[1]

                Rv3=R[0][:,2]*(-1.0)
                Rsum=Rsum+Rv3
                
            #the DNA axis of the current snapshot
            # print Rsum
            dnaaxis=numpy.array(Rsum[0])
            # print dnaaxis
            dnaaxis=numpy.dot(1/numpy.sqrt(numpy.dot(dnaaxis, dnaaxis)), dnaaxis)
            # print dnaaxis

            #the End/CENTER of DNA of the current snapshot
            # print end5center, end3center
            end5center=end5center/2.0
            end3center=end3center/2.0
            dnacenter=dnacenter/(2*(resEnd - resStart))
            dnaLength+=numpy.sqrt(numpy.dot(end5center - end3center, end5center - end3center))

            # if (frameIDs%50) == 1 :
                # print "  ...>>>>DNA 5' END CENTER: ", end5center
                # print "  ...>>>>DNA 3' END CENTER: ", end3center
                # print "  ...>>>>DNA CENTER: ", dnacenter
                # print "  ...>>>>DNA AXIS: ", dnaaxis
            # else :
                # pass
            print "  ...>>>>DNA 5' END CENTER: ", end5center
            print "  ...>>>>DNA 3' END CENTER: ", end3center
            print "  ...>>>>DNA CENTER: ", dnacenter
            print "  ...>>>>DNA AXIS: ", dnaaxis


            ##part for ions
            ionSELECT="name "+ionName
            ionsu=u.selectAtoms(ionSELECT)
            ionsCoor=ionsu.coordinates()

            if frameIDs==1:
                ionNumber=len(ionsu.atoms)

            boxsize=list([ts.dimensions[0], ts.dimensions[1], ts.dimensions[2]])
            boxsize=numpy.array(boxsize)
            for i in range(-1, 2, 1):
                for j in range(-1, 2, 1):
                    for k in range(-1, 2, 1):
                        tmp_ionsCoor=ionsCoor+(i, j, k)*boxsize
                        for ionIndex in range(len(tmp_ionsCoor)):
                            tmp_ioncoor=tmp_ionsCoor[ionIndex]
                            #print type(tmp_ioncoor), type(dnacenter)
                            vector1=tmp_ioncoor - end5center
                            vector2=tmp_ioncoor - end3center
                            # print vector1, vector2
                            # print type(vector1), type(vector2)
                            angleIon5end=DNA_matrix.Get_two_vector_angle(vector1,dnaaxis)
                            angleIon3end=DNA_matrix.Get_two_vector_angle(vector2,dnaaxis)
                            d_ion2dnaaxis=DNA_matrix.Get_distance_point2line(tmp_ioncoor, dnacenter, dnaaxis)
                            # print angleIon5end, angleIon3end, d_ion2dnaaxis
                            if angleIon5end <=90.0 and angleIon3end >= 90.0 and d_ion2dnaaxis <= radiusRDF:
                                tmpN=int(d_ion2dnaaxis/radiusDelta)
                                arrayRDF[tmpN]+=1

            #print arrayRDF
            print "  ...Time Lasted:  %f seconds.\t" %(Time.time()-START_TIME)

            #update rdf file
            if (frameIDs%100) == 1:
                tmp_dnaLength=dnaLength
                tmp_boxVolumn=boxVolumn
                tmp_averageIonDensity=ionNumber/tmp_boxVolumn

                #file bakeup
                if os.path.isfile(output_name):
                    shutil.copyfile(output_name, output_name+".bak")
                f = open(output_name, 'w')
                for m in range(0, radiusN, 1):
                    radiusR0=m*radiusDelta
                    radiusR1=(m+1)*radiusDelta
                    areaCircle=math.pi*(radiusR1*radiusR1 - radiusR0*radiusR0)
                    f.write("%5.5f \t %5.5f \n" %(m*radiusDelta+radiusDelta/2.0, arrayRDF[m]/(dnaLength *frameIDs*areaCircle*tmp_averageIonDensity)))
                f.close()


    print "=================================================="
    print "RESULTS OUTPUT"
    print "=================================================="


    dnaLength=dnaLength/(frameIDs*1.0)
    boxVolumn=boxVolumn/(frameIDs*1.0)
    averageIonDensity=ionNumber/(boxVolumn)
    f = open(output_name, 'w')
    for m in range(0, radiusN, 1):
        print m*radiusDelta+radiusDelta/2.0, arrayRDF[m]/(frameIDs*1.0)

        radiusR0=m*radiusDelta
        radiusR1=(m+1)*radiusDelta
        areaCircle=math.pi*(radiusR1*radiusR1 - radiusR0*radiusR0)

        f.write("%5.5f \t %5.5f \n" %(m*radiusDelta+radiusDelta/2.0, arrayRDF[m]/(dnaLength *frameIDs*areaCircle*averageIonDensity)))
    
    f.close()
    print "=================================================="
    print "END"
    print "=================================================="

if __name__=="__main__":
    Trj_PBC_ION_DNA_RDF("20-at-nuc-neutralized.prmtop", "../../../equip6.dcd", end=100)
