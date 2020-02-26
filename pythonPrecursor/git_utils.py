# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 17:47:52 2015

@author: Christian Grønbæk

Pre-cursor to GISA C-code (github.com/ceegeeCode/GISA). This module contains the basic writhe computations (for line segments of a polygonal
curve) and some other utils. 

If you use this code 

0) it'll be "as is" -- I will have no time for support. The real alternative is to use the GISA C-code (github repo: see above)

1) To cite it: "C.Grønbæk., T.Hamelryck, P.Røgen: "GISA: Using Gauss integrals to identify rare conformations in protein structures", TO APPEAR. 

Once the paper is published, cite it with ref to that publication.  

"""

import math
import numpy as np

#from numba import int32, uint32, float32, float64
from numbapro import vectorize, cuda


import timeit
#import scipy.stats as stats

from Bio import PDB

parser = PDB.PDBParser()

targetDef = 'cpu'

#ignore warnings
import warnings
warnings.filterwarnings("ignore")

ErrorFile = r"C:/Users/Kristian/BioInformatics/projects/knots/code/errors_git.txt"

#for timing
from timeit import *

import time

##########################################################
## General utils
##########################################################

def w(segment1, segment2):
    '''To easily set the fct for computing the writhe.'''
    return w_pdb_2(segment1, segment2)
#    return w_pdb_1(segment1, segment2)


def w_pdb_1(segment1, segment2):
    '''Compute writhe contribution from two input line segments. Segments
    are pairs of vectors as returned by Biopython.'''
#    t1 = time.time()
    s1 = segment1
    s2 = segment2
    #form unit vectors at the segment extremes:
    v13 = s2[0] - s1[0]
    v23 = s2[0] - s1[1]
    v24 = s2[1] - s1[1]
    v14 = s2[1] - s1[0]
    v = [v13,v23,v24,v14,v13,v23]
    e = []
#    l13 = np.sqrt(PDB.Vector.__mul__(v13,v13))
#    l23 = np.sqrt(PDB.Vector.__mul__(v23,v23))
#    l24 = np.sqrt(PDB.Vector.__mul__(v24,v24))
#    l14 = np.sqrt(PDB.Vector.__mul__(v14,v14))
    l13 = np.sqrt(v13*v13)
    l23 = np.sqrt(v23*v23)
    l24 = np.sqrt(v24*v24)
    l14 = np.sqrt(v14*v14)
    ls = [l13,l23,l24,l14]
    for l in ls:
        if l == 0.0:
            return 0
    e13 = v13/l13
    e23 = v23/l23
    e24 = v24/l24
    e14 = v14/l14
    e = [e13,e23,e24,e14,e13,e23]
##        n_vect = np.linalg.norm(vect)
##        if n_vect == 0: #if one of the vectors is zero the two segments lie in a plane
##            return 0 
#        else:
#            n_vect = np.linalg.norm(vect)
#            e.append(vect/n_vect)
#    t2 = time.time()
    #compute the angles
    s = 0
    for i in range(1,len(e)-1):
        a = e[i-1]
        b =e[i]
        c = e[i+1]
        theta =  PDB.calc_angle(a,b,c)
#        print "angle is %f" % theta 
        #APPARENTLY THIS ANGLE IS ALWAYS POSITIVE AND BETWEEN O AND PI; PROBABLY THE ANGLE BETWEEN SEGMENT AB, AC (WITHOUT SIGN)
        s = s +theta
#        print "s is %f " % s
#    t3 = time.time
    w = np.sign(s)*2*np.pi - s
#    t4 = time.time()
#    print "t2-t1:", 10*(t2-t1)
#    print "t3-t2:", 10*(t3-t2)
#    print "t4-t3:",10*(t4-t3)
    return w  
    
    
def w_pdb_2(segment1, segment2):
    '''Compute writhe contribution from two input line segments. Segments
    are pairs of vectors as returned by Biopython.'''
#    t1 = time.time()
    s1 = segment1
    s2 = segment2
    #form unit vectors at the segment extremes:
    v13 = s2[0] - s1[0]
    v14 = s2[1] - s1[0]
    v23 = s2[0] - s1[1]
    v24 = s2[1] - s1[1]
    v = [v13,v23,v24,v14,v13,v23]
    e = []
    l13 = np.sqrt(PDB.Vector.__mul__(v13,v13))
    l23 = np.sqrt(PDB.Vector.__mul__(v23,v23))
    l24 = np.sqrt(PDB.Vector.__mul__(v24,v24))
    l14 = np.sqrt(PDB.Vector.__mul__(v14,v14))
#    l13 = np.sqrt(v13*v13)
#    l23 = np.sqrt(v23*v23)
#    l24 = np.sqrt(v24*v24)
#    l14 = np.sqrt(v14*v14)
    ls = [l13,l23,l24,l14]
    for l in ls:
        if l == 0.0:
            return 0
    e13 = v13/l13
    e23 = v23/l23
    e24 = v24/l24
    e14 = v14/l14
    e = [e13,e23,e24,e14,e13,e23]
##        n_vect = np.linalg.norm(vect)
##        if n_vect == 0: #if one of the vectors is zero the two segments lie in a plane
##            return 0 
#        else:
#            n_vect = np.linalg.norm(vect)
#            e.append(vect/n_vect)
    #compute the angles
    #compute the angles
#    t2 = time.time()
    s = 0
    for i in range(1,len(e)-1):
        a = e[i-1]
        b =e[i]
        c = e[i+1]
        aDotb = PDB.Vector.__mul__(a,b) #a*b #PDB.Vector.__mul__(a,b)
        aDotc = PDB.Vector.__mul__(a,c) #a*c #PDB.Vector.__mul__(a,c)
        bDotc = PDB.Vector.__mul__(b,c) #b*c #PDB.Vector.__mul__(b,c)
        aCrossb = PDB.Vector.__pow__(a,b)
#        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
#        cos = float(aDotb*bDotc - aDotc)/denom
#        sin = float(np.dot(aCrossb,c))/denom
        #can disregard the denominator as they are equal and positive:
        cos = float(aDotb*bDotc - aDotc)
        sin = float(PDB.Vector.__mul__(aCrossb,c))
#        theta = np.arctan2(sin,cos)
        if cos > 0:
            theta = np.arctan(sin/cos)
#            theta = sin/cos
        elif cos < 0:
            if sin > 0:
                theta = np.arctan(sin/cos) + np.pi
#             theta = sin/cos
            else:
                theta = np.arctan(sin/cos) - np.pi
        else: #iF cos ==0
            if sin > 0:
                theta = np.pi/2
            else:
                theta = -np.pi/2
#        print theta 
        #CHECK THE ANGLE COMP; MISSES SOME SIGNS
#        print "angle is %f" % theta 
        s = s +theta
#    t3 = time.time()
    w = np.sign(s)*2*np.pi - s
#    t4 = time.time()
#    print "t2-t1:", 10*(t2-t1)
#    print "t3-t2:", 10*(t3-t2)
#    print "t4-t3:",10*(t4-t3)
    return w         


def getPolygonalChain(PDBfilename, outNumpy_b = 0, optionOut_b = False, outputFile = ''):
    '''Derives the polygonal chain representation of the protein 
    from its PDB file. Output: array of line segments between the
    C-alpha atoms and ordered by these.'''
    CaChain = {}
    polyCaChain = {}
    structure = parser.get_structure(PDBfilename,PDBfilename) #file id = PDBfilename
    #we trust that the top100 files all have crystallograhic content for modelId = 0:
    model = structure[0]
    #first loop though the chains; for each residue we store the CA info: 
    cntCA = 0
    if outNumpy_b == 1:
        for chain in model:
            print "Reading chain: %s" % chain
            if chain.id.isspace():
                chain.id = '>'
                print "Chain id is blank and gets subst with: %s" % chain.id
            for residue in chain:
                if PDB.is_aa(residue):
                    if residue.id[0] == ' ': #residue.id[0] is the hetero-flag; if not blank the residue will contain hetero-atoms (HETATM in pdb format)
                        #we only add chain as key if there is an aa in the chain:
                        if not(CaChain.has_key(chain.id)):
                            print "Recording a new chain id in the structure: %s" % chain.id
                            CaChain[chain.id] = []
                            polyCaChain[chain.id] = []
                            cntCA = 0
        #                CA = residue['CA']
                        vCA = np.array(residue['CA'].get_vector()) 
                        CaChain[chain.id].append(vCA)
                        if cntCA > 0:
        #                    v = vCA - vPrev
                            polyCaChain[chain.id].append([vPrev, vCA]) #[vCA,vPrev]??
        #                    print v
                        vPrev = vCA
                        cntCA +=1
        #                print cntCA
        #                print vCA
    else: 
        for chain in model:
            print "Reading chain: %s" % chain
            if chain.id.isspace():
                chain.id = '>'
                print "Chain id is blank and gets subst with: %s" % chain.id
            for residue in chain:
                if PDB.is_aa(residue):
                    if residue.id[0] == ' ': #residue.id[0] is the hetero-flag; if not blank the residue will contain hetero-atoms (HETATM in pdb format)
                        #we only add chain as key if there is an aa in the chain:
                        if not(CaChain.has_key(chain.id)):
                            print "Recording a new chain id in the structure: %s" % chain.id
                            CaChain[chain.id] = []
                            polyCaChain[chain.id] = []
                            cntCA = 0
        #                CA = residue['CA']
                        vCA = residue['CA'].get_vector() 
                        CaChain[chain.id].append(vCA)
                        if cntCA > 0:
        #                    v = vCA - vPrev
                            polyCaChain[chain.id].append([vPrev, vCA]) #[vCA,vPrev]??
        #                    print v
                        vPrev = vCA
                        cntCA +=1
        #                print cntCA
        #                print vCA
    if optionOut_b:
        with open(outputFile, 'w') as of:
            for k in polyCaChain.keys():
                of.write('file: ' + PDBfilename + ';' + 'Chain: ' + k)
                for v in polyCaChain[k]:
                    x0, y0, z0 = v[0]
    #                print x0, y0, z0
                    x1, y1, z1 = v[1]
                    s = str(x0) + ';' + str(y0) + ';' + str(z0) + ';' + str(x1) + ';' + str(y1) + ';' + str(z1) + '\n'
                    of.write(s)
        
    return CaChain, polyCaChain 
    
#######################################################
## Utils for perturbations
#######################################################
    
def perturbPolygonalChain(CaChain, polyCaChain,
                 residueNrRange = range(0,5), 
                 perturbation = [(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3)]):
    '''Perturbs the input chain: perturbs the residues with numbers in the 
    provided range, residueNrRange, added to the number of the first residue
    in the chain; the perturbation amounts is given in 3d coords (x,y,z) in units
    of Å. 
    Obs: only residues in the chain having with indices in the residueNrRange will be
    perturbed; in particular, if no such indices exist the chain will not be perturbed.
    
    Input: as output from getPolygonalChain.
    Output: perturbed chain of CAs, corresponding polygonal chain.'''
    
    perChain = {}
    perPolyCaChain = {}
    segNrRange = {}
    for chain in CaChain.keys():
        resNr = 0
        perNr = 0
        perChain[chain] = []
        for c in CaChain[chain]:
            if resNr in residueNrRange:
                perC = PDB.Vector.__add__(c, PDB.Vector(perturbation[perNr]))
                perChain[chain].append(perC)
                perNr += 1
            else:
                perChain[chain].append(c)
            resNr +=1
        #derive the corresponding polygonal chain:
        perPolyCaChain[chain] = []
        CaCnt = 0
        for c in perChain[chain]:
            Ca = c 
            if CaCnt > 0:
                perPolyCaChain[chain].append([prevCa, Ca])
            prevCa = Ca
            CaCnt += 1
        #comes in handy later to have readily available the polychain segment numbers for
        #which the segment is perturbed:
        segNr = 0
        segNrRange[chain] = []
        print chain
        for s in polyCaChain[chain]:
            if perPolyCaChain[chain][segNr] <> s:
                segNrRange[chain].append(segNr)
                #print segNr
            segNr +=1
    return perChain, perPolyCaChain, segNrRange
    
    
#################################################################################################################
## Utils for vectorized versions
#################################################################################################################

       
@vectorize(['float64(float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64)'], 
            target = 'cpu')
def wVect_1(s101, s102, s103, s111, s112, s113, s201, s202, s203, s211, s212, s213):
    '''Compute writhe contribution from two input line segments. Segments
    are pairs of vectors as returned by Biopython.'''
#    t1 = time.time()
    w = 0
    s1 = np.array([[s101, s102, s103],[s111, s112, s113]])
    s2 = np.array([[s201, s202, s203],[s211, s212, s213]])
#    form unit vectors at the segment extremes:
    v13 = s2[0] - s1[0]
    v23 = s2[0] - s1[1]
    v24 = s2[1] - s1[1]
    v14 = s2[1] - s1[0]
    e = []
    l13 = np.sqrt(np.dot(v13,v13))
    l23 = np.sqrt(np.dot(v23,v23))
    l24 = np.sqrt(np.dot(v24,v24))
    l14 = np.sqrt(np.dot(v14,v14))
#    l13 = np.linalg.norm(v13)
#    l23 = np.linalg.norm(v23)
#    l24 = np.linalg.norm(v24)
#    l14 = np.linalg.norm(v14)
    ls = [l13,l23,l24,l14]
    for l in ls:
        if l == 0.0:
            return 0
    e13 = v13/l13
    e23 = v23/l23
    e24 = v24/l24
    e14 = v14/l14
    e = [e13,e23,e24,e14,e13,e23]
##        n_vect = np.linalg.norm(vect)
##        if n_vect == 0: #if one of the vectors is zero the two segments lie in a plane
##            return 0 
#        else:
#            n_vect = np.linalg.norm(vect)
#            e.append(vect/n_vect)
#    t2 = time.time()
    #compute the angles
    s = 0
    for i in range(1,len(e)-1):
        theta = 0 #reset
        a = e[i-1]
        b =e[i]
        c = e[i+1]
#        theta = PDB.calc_angle(a,b,c)
#        print "angle is %f" % theta 
        aDotb = np.dot(a,b) #PDB.Vector.__mul__(a,b)
#        print "aDotb: %f " % aDotb
        aDotc = np.dot(a,c) #PDB.Vector.__mul__(a,c)
        bDotc = np.dot(b,c) #PDB.Vector.__mul__(b,c)
        aCrossb = np.cross(a,b)
##        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
##        cos = float(aDotb*bDotc - aDotc)/denom
##        sin = float(np.dot(aCrossb,c))/denom
#        #can disregard the denominator as they are equal:
        cos = float(aDotb*bDotc - aDotc)
        sin = float(np.dot(aCrossb,c))
##        theta = np.arctan2(sin,cos)
        if cos > 0:
            theta = np.arctan(sin/cos)
#            theta = sin/cos
        elif cos < 0:
            if sin > 0:
                theta = np.arctan(sin/cos) + np.pi
#             theta = sin/cos
            else:
                theta = np.arctan(sin/cos) - np.pi
        else: #iF cos ==0
            if sin > 0:
                theta = np.pi/2
            else:
                theta = -np.pi/2
#        print theta 
        #CHECK THE ANGLE COMP; MISSES SOME SIGNS
        s = s +theta
#        print "s is %f " % s
#    t3 = time.time()
    w = np.sign(s)*2*np.pi - s
#    t4 = time.time()
#    print "t2-t1:", 10*(t2-t1)
#    print "t3-t2:", 10*(t3-t2)
#    print "t4-t3:",10*(t4-t3)
    return w 

    
@vectorize(['float64(float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64)'], 
            target = targetDef)
def wVect_2(s101, s102, s103, s111, s112, s113, s201, s202, s203, s211, s212, s213):
    '''Compute writhe contribution from two input line segments. Segments
    are pairs of vectors as returned by Biopython.'''
#    t1 = time.time()
    w = 0
#    form unit vectors at the segment extremes:
    #v13:
    v131 = s201 - s101
    v132 = s202 - s102
    v133 = s203 - s103
    #v23:
    v231 = s201 - s111
    v232 = s202 - s112
    v233 = s203 - s113
    #v24:
    v241 = s211 - s111
    v242 = s212 - s112
    v243 = s213 - s113
    #v14:
    v141 = s211 - s101
    v142 = s212 - s102
    v143 = s213 - s103
##    v = [v13,v23,v24,v14,v13,v23]
#    e = []
    l13 = math.sqrt(v131*v131 + v132*v132 + v133*v133)
    l23 = math.sqrt(v231*v231 + v232*v232 + v233*v233)
    l24 = math.sqrt(v241*v241 + v242*v242 + v243*v243)
    l14 = math.sqrt(v141*v141 + v142*v142 + v143*v143)
    ls = (l13,l23,l24,l14)
    for l in ls:
        if l < 1e-16:
            return 0
    e131 = v131/l13
    e132 = v132/l13
    e133 = v133/l13
#    e13 = (e131,e132,e133)
    e231 = v231/l23
    e232 = v232/l23
    e233 = v233/l23
#    e23 = (e231,e232,e233)
    e241 = v241/l24
    e242 = v242/l24
    e243 = v243/l24
#    e24 = (e241,e242,e243)
    e141 = v141/l14
    e142 = v142/l14
    e143 = v143/l14
#    e14 = (e141,e142,e143)
#    e = (e13,e23,e24,e14,e13,e23)
    #compute the angles
    s = 0
#    for i in range(1,len(e)-1):
#        theta = 0 #reset
#        a = e[i-1]
#        b =e[i]
#        c = e[i+1]
##        theta = PDB.calc_angle(a,b,c)
##        print "angle is %f" % theta 
#        aDotb = np.dot(a,b) #PDB.Vector.__mul__(a,b)
#        aDotc = np.dot(a,c) #PDB.Vector.__mul__(a,c)
#        bDotc = np.dot(b,c) #PDB.Vector.__mul__(b,c)
#        aCrossb = np.cross(a,b)
    #We have to do this explicity, ie without referrring to itereations and vectors:
    #A = e13, B= e23, C = e24:
    x = e131
    y = e132
    z = e133
    a = e231
    b = e232
    c = e233
    u = e241
    v = e242
    w = e243
    A_Dot_B = x*a + y*b + z*c
    A_Dot_C = x*u + y*v + z*w
    B_Dot_C = a*u + b*v + c*w
    A_Cross_B_1, A_Cross_B_2, A_Cross_B_3 = y*c - z*b, z*a - x*c, x*b - y*a 
###        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
###        cos = float(aDotb*bDotc - aDotc)/denom
###        sin = float(np.dot(aCrossb,c))/denom
##        #can disregard the denominator as they are equal:
    cos = float(A_Dot_B*B_Dot_C - A_Dot_C)
    sin = float(A_Cross_B_1*u + A_Cross_B_2*v + A_Cross_B_3*w)
###        theta = np.arctan2(sin,cos)
    if cos > 0:
        theta = math.atan(sin/cos)
#            theta = sin/cos
    elif cos < 0:
        if sin > 0:
            theta = math.atan(sin/cos) + math.pi
#             theta = sin/cos
        else:
            theta = math.atan(sin/cos) - math.pi
    else: #iF cos ==0
        if sin > 0:
            theta = math.pi/2
        else:
            theta = -math.pi/2
#        print theta 
    s = s +theta
#        print "s is %f " % s
    #A = e23, B= e24, C = e14:
    x = e231
    y = e232
    z = e233
    a = e241
    b = e242
    c = e243
    u = e141
    v = e142
    w = e143
    A_Dot_B = x*a + y*b + z*c
    A_Dot_C = x*u + y*v + z*w
    B_Dot_C = a*u + b*v + c*w
    A_Cross_B_1, A_Cross_B_2, A_Cross_B_3 = y*c - z*b, z*a - x*c, x*b - y*a 
###        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
###        cos = float(aDotb*bDotc - aDotc)/denom
###        sin = float(np.dot(aCrossb,c))/denom
##        #can disregard the denominator as they are equal:
    cos = float(A_Dot_B*B_Dot_C - A_Dot_C)
    sin = float(A_Cross_B_1*u + A_Cross_B_2*v + A_Cross_B_3*w)
###        theta = np.arctan2(sin,cos)
    if cos > 0:
        theta = math.atan(sin/cos)
#            theta = sin/cos
    elif cos < 0:
        if sin > 0:
            theta = math.atan(sin/cos) + math.pi
#             theta = sin/cos
        else:
            theta = math.atan(sin/cos) - math.pi
    else: #iF cos ==0
        if sin > 0:
            theta = math.pi/2
        else:
            theta = -math.pi/2
#        print theta 
    s = s +theta
#        print "s is %f " % s
    #A = e24, B= e14, C = e13:
    x = e241
    y = e242
    z = e243
    a = e141
    b = e142
    c = e143
    u = e131
    v = e132
    w = e133
    A_Dot_B = x*a + y*b + z*c
    A_Dot_C = x*u + y*v + z*w
    B_Dot_C = a*u + b*v + c*w
    A_Cross_B_1, A_Cross_B_2, A_Cross_B_3 = y*c - z*b, z*a - x*c, x*b - y*a 
###        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
###        cos = float(aDotb*bDotc - aDotc)/denom
###        sin = float(np.dot(aCrossb,c))/denom
##        #can disregard the denominator as they are equal:
    cos = float(A_Dot_B*B_Dot_C - A_Dot_C)
    sin = float(A_Cross_B_1*u + A_Cross_B_2*v + A_Cross_B_3*w)
###        theta = np.arctan2(sin,cos)
    if cos > 0:
        theta = math.atan(sin/cos)
#            theta = sin/cos
    elif cos < 0:
        if sin > 0:
            theta = math.atan(sin/cos) + math.pi
#             theta = sin/cos
        else:
            theta = math.atan(sin/cos) - math.pi
    else: #iF cos ==0
        if sin > 0:
            theta = math.pi/2
        else:
            theta = -math.pi/2
#        print theta 
    s = s +theta
#        print "s is %f " % s
    #A = e14, B= e13, C = e23:
    x = e141
    y = e142
    z = e143
    a = e131
    b = e132
    c = e133
    u = e231
    v = e232
    w = e233
    A_Dot_B = x*a + y*b + z*c
    A_Dot_C = x*u + y*v + z*w
    B_Dot_C = a*u + b*v + c*w
    A_Cross_B_1, A_Cross_B_2, A_Cross_B_3 = y*c - z*b, z*a - x*c, x*b - y*a 
###        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
###        cos = float(aDotb*bDotc - aDotc)/denom
###        sin = float(np.dot(aCrossb,c))/denom
##        #can disregard the denominator as they are equal:
    cos = float(A_Dot_B*B_Dot_C - A_Dot_C)
    sin = float(A_Cross_B_1*u + A_Cross_B_2*v + A_Cross_B_3*w)
###        theta = np.arctan2(sin,cos)
    if cos > 0:
        theta = math.atan(sin/cos)
#            theta = sin/cos
    elif cos < 0:
        if sin > 0:
            theta = math.atan(sin/cos) + math.pi
#             theta = sin/cos
        else:
            theta = math.atan(sin/cos) - math.pi
    else: #iF cos ==0
        if sin > 0:
            theta = math.pi/2
        else:
            theta = -math.pi/2
#        print theta 
    s = s +theta
#        print "s is %f " % s
#    t3 = time.time()
    #compute the sign: (could not find an available sign-fct accepted when target = 'gpu')
    if abs(s) - s > s:
        sign_s = -1
    elif abs(s) - s < s:
        sign_s = 1
    else: #sign(0) == 0
        sign_s = 0
    w = sign_s*2*math.pi - s
#    t4 = time.time()
#    print "t2-t1:", 10*(t2-t1)
#    print "t3-t2:", 10*(t3-t2)
#    print "t4-t3:",10*(t4-t3)
    return w 

@vectorize(['float64(float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64)'], 
            target = 'cpu')
def wVect_2_np(s101, s102, s103, s111, s112, s113, s201, s202, s203, s211, s212, s213):
    '''Compute writhe contribution from two input line segments. Segments
    are pairs of vectors as returned by Biopython.'''
#    t1 = time.time()
    w = 0
#    form unit vectors at the segment extremes:
    #v13:
    v131 = s201 - s101
    v132 = s202 - s102
    v133 = s203 - s103
    #v23:
    v231 = s201 - s111
    v232 = s202 - s112
    v233 = s203 - s113
    #v24
    v241 = s211 - s111
    v242 = s212 - s112
    v243 = s213 - s113
    #v14:
    v141 = s211 - s101
    v142 = s212 - s102
    v143 = s213 - s103
##    v = [v13,v23,v24,v14,v13,v23]
#    e = []
    l13 = np.sqrt(v131*v131 + v132*v132 + v133*v133)
    l23 = np.sqrt(v231*v231 + v232*v232 + v233*v233)
    l24 = np.sqrt(v241*v241 + v242*v242 + v243*v243)
    l14 = np.sqrt(v141*v141 + v142*v142 + v143*v143)
    ls = (l13,l23,l24,l14)
    for l in ls:
        if l == 0.0:
            return 0
    e131 = v131/l13
    e132 = v132/l13
    e133 = v133/l13
#    e13 = (e131,e132,e133)
    e231 = v231/l23
    e232 = v232/l23
    e233 = v233/l23
#    e23 = (e231,e232,e233)
    e241 = v241/l24
    e242 = v242/l24
    e243 = v243/l24
#    e24 = (e241,e242,e243)
    e141 = v141/l14
    e142 = v142/l14
    e143 = v143/l14
#    e14 = (e141,e142,e143)
#    e = (e13,e23,e24,e14,e13,e23)
    #compute the angles
    s = 0
#    for i in range(1,len(e)-1):
#        theta = 0 #reset
#        a = e[i-1]
#        b =e[i]
#        c = e[i+1]
##        theta = PDB.calc_angle(a,b,c)
##        print "angle is %f" % theta 
#        aDotb = np.dot(a,b) #PDB.Vector.__mul__(a,b)
#        aDotc = np.dot(a,c) #PDB.Vector.__mul__(a,c)
#        bDotc = np.dot(b,c) #PDB.Vector.__mul__(b,c)
#        aCrossb = np.cross(a,b)
    #We have to do this explicity, ie without referrring to itereations and vectors:
    #A = e13, B= e23, C = e24:
    x = e131
    y = e132
    z = e133
    a = e231
    b = e232
    c = e233
    u = e241
    v = e242
    w = e243
    A_Dot_B = x*a + y*b + z*c
    A_Dot_C = x*u + y*v + z*w
    B_Dot_C = a*u + b*v + c*w
    A_Cross_B_1, A_Cross_B_2, A_Cross_B_3 = y*c - z*b, z*a - x*c, x*b - y*a 
###        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
###        cos = float(aDotb*bDotc - aDotc)/denom
###        sin = float(np.dot(aCrossb,c))/denom
##        #can disregard the denominator as they are equal:
    cos = float(A_Dot_B*B_Dot_C - A_Dot_C)
    sin = float(A_Cross_B_1*u + A_Cross_B_2*v + A_Cross_B_3*w)
###        theta = np.arctan2(sin,cos)
    if cos > 0:
        theta = np.arctan(sin/cos)
#            theta = sin/cos
    elif cos < 0:
        if sin > 0:
            theta = np.arctan(sin/cos) + np.pi
#             theta = sin/cos
        else:
            theta = np.arctan(sin/cos) - np.pi
    else: #iF cos ==0
        if sin > 0:
            theta = np.pi/2
        else:
            theta = -np.pi/2
#        print theta 
    s = s +theta
#        print "s is %f " % s
    #A = e23, B= e24, C = e14:
    x = e231
    y = e232
    z = e233
    a = e241
    b = e242
    c = e243
    u = e141
    v = e142
    w = e143
    A_Dot_B = x*a + y*b + z*c
    A_Dot_C = x*u + y*v + z*w
    B_Dot_C = a*u + b*v + c*w
    A_Cross_B_1, A_Cross_B_2, A_Cross_B_3 = y*c - z*b, z*a - x*c, x*b - y*a 
###        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
###        cos = float(aDotb*bDotc - aDotc)/denom
###        sin = float(np.dot(aCrossb,c))/denom
##        #can disregard the denominator as they are equal:
    cos = float(A_Dot_B*B_Dot_C - A_Dot_C)
    sin = float(A_Cross_B_1*u + A_Cross_B_2*v + A_Cross_B_3*w)
###        theta = np.arctan2(sin,cos)
    if cos > 0:
        theta = np.arctan(sin/cos)
#            theta = sin/cos
    elif cos < 0:
        if sin > 0:
            theta = np.arctan(sin/cos) + np.pi
#             theta = sin/cos
        else:
            theta = np.arctan(sin/cos) - np.pi
    else: #iF cos ==0
        if sin > 0:
            theta = np.pi/2
        else:
            theta = -np.pi/2
#        print theta 
    s = s +theta
#        print "s is %f " % s
    #A = e24, B= e14, C = e13:
    x = e241
    y = e242
    z = e243
    a = e141
    b = e142
    c = e143
    u = e131
    v = e132
    w = e133
    A_Dot_B = x*a + y*b + z*c
    A_Dot_C = x*u + y*v + z*w
    B_Dot_C = a*u + b*v + c*w
    A_Cross_B_1, A_Cross_B_2, A_Cross_B_3 = y*c - z*b, z*a - x*c, x*b - y*a 
###        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
###        cos = float(aDotb*bDotc - aDotc)/denom
###        sin = float(np.dot(aCrossb,c))/denom
##        #can disregard the denominator as they are equal:
    cos = float(A_Dot_B*B_Dot_C - A_Dot_C)
    sin = float(A_Cross_B_1*u + A_Cross_B_2*v + A_Cross_B_3*w)
###        theta = np.arctan2(sin,cos)
    if cos > 0:
        theta = np.arctan(sin/cos)
#            theta = sin/cos
    elif cos < 0:
        if sin > 0:
            theta = np.arctan(sin/cos) + np.pi
#             theta = sin/cos
        else:
            theta = np.arctan(sin/cos) - np.pi
    else: #iF cos ==0
        if sin > 0:
            theta = np.pi/2
        else:
            theta = -np.pi/2
#        print theta 
    s = s +theta
#        print "s is %f " % s
    #A = e14, B= e13, C = e23:
    x = e141
    y = e142
    z = e143
    a = e131
    b = e132
    c = e133
    u = e231
    v = e232
    w = e233
    A_Dot_B = x*a + y*b + z*c
    A_Dot_C = x*u + y*v + z*w
    B_Dot_C = a*u + b*v + c*w
    A_Cross_B_1, A_Cross_B_2, A_Cross_B_3 = y*c - z*b, z*a - x*c, x*b - y*a 
###        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
###        cos = float(aDotb*bDotc - aDotc)/denom
###        sin = float(np.dot(aCrossb,c))/denom
##        #can disregard the denominator as they are equal:
    cos = float(A_Dot_B*B_Dot_C - A_Dot_C)
    sin = float(A_Cross_B_1*u + A_Cross_B_2*v + A_Cross_B_3*w)
###        theta = np.arctan2(sin,cos)
    if cos > 0:
        theta = np.arctan(sin/cos)
#            theta = sin/cos
    elif cos < 0:
        if sin > 0:
            theta = np.arctan(sin/cos) + np.pi
#             theta = sin/cos
        else:
            theta = np.arctan(sin/cos) - np.pi
    else: #iF cos ==0
        if sin > 0:
            theta = np.pi/2
        else:
            theta = -np.pi/2
#        print theta 
    s = s +theta
#        print "s is %f " % s
#    t3 = time.time()
    w = np.sign(s)*2*np.pi - s
#    t4 = time.time()
#    print "t2-t1:", 10*(t2-t1)
#    print "t3-t2:", 10*(t3-t2)
#    print "t4-t3:",10*(t4-t3)
    return w 