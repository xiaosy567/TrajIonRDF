#1/usr/bin/env python
#-*- coding:utf-8 -*-
'''
Created at 2012.02.23
@change:\n
    - 2012.02.23.
        - Copied those functions from parallel_analysis.py to make it simplify.
'''

import numpy
import math
import sys
import os
from numpy import matrix, linalg, dot
sys.path.append('./Coor')
import atomlib
import unit_atom
sys.path.append('./mymath')
import least_squares_fitting

def Get_base_rotate_matrix(experim_coor, base_name):
    '''
    Reading in the coor_list data and base name,  return a 3*3 rotation matrix.\n
    The algorithm come form the 3DNA. you can find it from the 3DNA manual. The\
            details are in page 22-23.\n
    note: the unit is A.
    @param experim_coor: The coordinate for the 9 atoms for A G nucleic acid base or \
            6 atoms for C T U nucleic acid base.\n
    @param base_name: Only A,T,C,G,U are allowd.
    @type  experim_coor: numpy.array
    @type  base_name:    char
    @rtype:  list
    @return:  rotation matrix in list[0] and origin coordinate in list[1]
    '''

    N = 0
    if 'A' in base_name :
        standard_coor = atomlib.BASE_A_array 
        N = 9
    elif 'T' in base_name :
        standard_coor = atomlib.BASE_T_array 
        N = 6
    elif 'C' in base_name:
        standard_coor = atomlib.BASE_C_array 
        N = 6
    elif 'G' in base_name :
        standard_coor = atomlib.BASE_G_array 
        N = 9
    elif 'U' in base_name :
        standard_coor = atomlib.BASE_U_array 
        N = 6
    else:
        print "The residue name %s not found in my standard lib." %base_name
        sys.exit()
    
    result=least_squares_fitting.Fitting(standard_coor,experim_coor)
    return result


def Get_base_RMSD(experim_coor,base_name, origin, orient):
    N = 0 
    if 'A' in base_name :
        standard_coor = atomlib.BASE_A_array 
        N = 9 
    elif 'T' in base_name :
        standard_coor = atomlib.BASE_T_array 
        N = 6 
    elif 'C' in base_name:
        standard_coor = atomlib.BASE_C_array 
        N = 6 
    elif 'G' in base_name :
        standard_coor = atomlib.BASE_G_array 
        N = 9 
    elif 'U' in base_name :
        standard_coor = atomlib.BASE_U_array 
        N = 6 
    else:
        pass
            

    Fit_coor=numpy.matrix(standard_coor) * numpy.matrix(orient).T + origin

#    print Fit_coor
#    print experim_coor
    Diff_coor=experim_coor - Fit_coor

    ssum=0.0
    for i in range(N):
        for j in range(3):
            ssum = ssum + Diff_coor[i,j]*Diff_coor[i,j]

    RMSD=math.sqrt(ssum/N)

    return RMSD


###above is copied from ZhuHong

def Get_RMSD(standard_coor, experim_coor):
    '''
     Function: used to calculate the RMSD
     Formate: Get_RMSD(standard_coor, experim_coor)
     Input: standard_coor (i.e., the reference coordinates)
            experim_coor  
     Return: RMSD 
     @type  standard_coor: numpy.array 
     @type  experim_coor: numpy.array
     ''' 
    
    Diff_coor=experim_coor - standard_coor
    #print Diff_coor
    
    ssum=0.0
    N=experim_coor.shape[0]
    for i in range(N):
        for j in range(3):
            ssum = ssum + Diff_coor[i,j]*Diff_coor[i,j]

    RMSD=math.sqrt(ssum/N)
    return RMSD     


def Get_fit_coor(standard_coor, experim_coor):
    '''
     This function is used to fit the strandard frame to the experimental frame.
     Formate: Get_fit_coor(standard_coor, experim_coor)
     Input: standard_coor >> reference coordinates
            experim_coor  >> coordiantes to be fitted
     Return fit_coor      >> fitted coordinates
     @type  standard_coor: numpy.array
     @type  experim_coor: numpy.array
     @type  fit_coor: numpy.array
     '''

    rotationMatrix, origin=least_squares_fitting.Fitting(standard_coor,experim_coor) 
    #print rotationMatrix
    #print origin
 
    fit_coor = numpy.matrix(standard_coor) * numpy.matrix(rotationMatrix).T + origin 
    #print fit_coor
    return fit_coor

def Get_two_vector_angle(vector1, vector2):
    '''
      This function is used to calculate the angle of two vectors.
      Format: Get_two_vector_angle(vector1, vector2)
      Input:  vector1
              vector2
      Return: angle (0 ~ 180 degree)
      @type   vector: numpy.array
      @type   vector: numpy.array
      @type   angle: float
    '''
    nvector1=numpy.dot(1/numpy.sqrt(numpy.dot(vector1, vector1)), vector1)
    nvector2=numpy.dot(1/numpy.sqrt(numpy.dot(vector2, vector2)), vector2)
    cosvalue=numpy.dot(nvector1, nvector2.transpose())
    angle=numpy.arccos(cosvalue)*180.0/math.pi
    return angle
    
def Get_distance_point2line(Spoint, Lpoint, Lvector):
    '''
      This function is used to calculate the distance between a point [x, y, z]
      and a line (defined by a point [x0, y0, z0] and a vector[v1, v2, v3])
            x-x0   y-y0   z-z0
      Line: -----=-----= -----
             v1     v2    v3
      Input: Spoint
             Lpoint
             Lvector
      Return distance
      @type  Spoint: numpy.array
             Lpoint: numpy.array
             Lvector: numpy.array
             Distrance: float
    '''
    tmp1=numpy.subtract(Lpoint, Spoint)
    tmp2=numpy.cross(tmp1, Lvector)
    distance=numpy.sqrt(numpy.dot(tmp2, tmp2))/(numpy.sqrt(numpy.dot(Lvector, Lvector)))
    return distance

    


#if __name__=="__main__":
#   vector1=numpy.array([ 0, 0,  1])
#   print vector1
#   vector2=numpy.array([ 0,  1,  0.5])
#   print vector2
#   print Get_two_vector_angle(vector1, vector2)
#   
#    p1=numpy.array([0, 0, 1])
#    p2=numpy.array([0, 0, 0])
#    vec=numpy.array([1, 1, 1])
#    print Get_distance_point2line(p1, p2, vec)
