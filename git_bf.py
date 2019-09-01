# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 21:43:12 2015

@author: Christian Grønbæk
"""

import math
import numpy as np

import timeit

from git_utils import *


from Bio import PDB

parser = PDB.PDBParser()

#ignore warnings
import warnings
warnings.filterwarnings("ignore")

ErrorFile = r"C:/Users/Kristian/BioInformatics/projects/knots/code/errors_git.txt"

#for timing
from timeit import *

import time

from git_utils import *

#####################################################################
# "brute force" various versions 
#####################################################################

def I12_bf(pChain, m = 0, n = 0, print_b = 1):
    '''Computes the measure I1234 by brute force that is, directly off
    the defining formula.'''
    I = 0 #init
    L = len(pChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    cnt = 0
    wVal_ij = 0
    for i in range(m, n+1):
        for j in range(i+1, n+1):
            cnt +=1
            if cnt%100000 == 0:
                print cnt
            wVal_ij = w(pChain[i],pChain[j])
            I += wVal_ij
    if print_b > 0.5:
        print "I12: %f" % I
    return I
    
def I12_bf_dict(PDBfilename):
    out = {}
    CaChain, pChain = getPolygonalChain(PDBfilename)

    for chain in pChain.keys():
        L = len(pChain[chain])
        if not(out.has_key(chain)):
            out[chain] = {}
        for m in range(L):
            for n in range(m,L):
                out[chain][(m,n)] = I12_bf(pChain[chain], m=m, n=n, print_b = 0)
    return out 
    
def I1234_bf(pChain):
    '''Computes the measure I1234 by brute force that is, directly off
    the defining formula.'''
    I1234 = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    cnt = 0
    for i in range(L):
        for j in range(i+1, L):
            for k in range(j+1, L):
#                wVal_ij = w(pChain[i],pChain[j])
                for l in range(k+1, L):
                    cnt +=1
                    if cnt%100000 == 0:
                        print cnt
                    wVal_ij = w(pChain[i],pChain[j])
                    wVal_kl = w(pChain[k],pChain[l])
                    I1234 += wVal_ij*wVal_kl
    print I1234
    return I1234    

def I1234_bf2(pChain):
    '''Computes the measure I1234 by brute force that is, directly off
    the defining formula.'''
    I1234 = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    cnt = 0
    for i in range(L):
        for j in range(i+1, L):
            wVal_ij = w(pChain[i],pChain[j])
            for k in range(j+1, L):
#                wVal_ij = w(pChain[i],pChain[j])
                for l in range(k+1, L):
                    cnt +=1
                    if cnt%100000 == 0:
                        print cnt
                    wVal_kl = w(pChain[k],pChain[l])
                    I1234 += wVal_ij*wVal_kl
    return I1234             


def I1234_bf3(pChain):
    '''Computes the measure I1234 by brute that is, directly off
    the defining formula.'''
    I1234 = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    cnt = 0
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    for i in range(L):
        for j in range(i+1, L):
            wVal_ij = wDict[(i,j)]
            for k in range(j+1, L):
#                wVal_ij = w(pChain[i],pChain[j])
                for l in range(k+1, L):
                    cnt +=1
#                    if cnt%100000 == 0:
#                        print cnt
                    I1234 += wVal_ij*wDict[(k,l)]
    return I1234   

def I1234_full_bf3(pChain, m=0, n=0, print_b = 1):
    '''Computes the measure I1234_full@(m,n) by brute that is, directly off
    the defining formula.'''
    I1234 = 0 #init
    L = len(pChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    wDict = {}
    cnt = 0
    for i in range(m,n+1):
        for j in range(i+1, n+1):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    for i in range(m,n+1):
        for j in range(i+1, n+1):
            wVal_ij = wDict[(i,j)]
            for k in range(j+1, n+1):
#                wVal_ij = w(pChain[i],pChain[j])
                for l in range(k+1, n+1):
                    cnt +=1
#                    if cnt%100000 == 0:
#                        print cnt
                    I1234 += wVal_ij*wDict[(k,l)]
    if print_b > 0.5:
        print I1234
    return I1234   

def I1234_full_bf3_dict(PDBfilename):
    out = {}
    CaChain, pChain = getPolygonalChain(PDBfilename)

    for chain in pChain.keys():
        L = len(pChain[chain])
        if not(out.has_key(chain)):
            out[chain] = {}
        for m in range(L):
            for n in range(m,L):
                out[chain][(m,n)] = I1234_full_bf3(pChain[chain], m=m, n=n, print_b = 0)
    return out 

def I1234_bf4(pChain):
    '''Computes the measure I1234 by almost brute force that is, almost directly off
    the defining formula.'''
    I1234 = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    s1 = {}
    for j in range(L):
        s1[j] =0
        for i in range(0,j):
            s1[j] += wDict[(i,j)]
    s2 = {}
    for j in range(L):
        s2[j] =0
        for k in range(j+1, L):
#                wVal_ij = w(pChain[i],pChain[j])
            for l in range(k+1, L):
                s2[j] += wDict[(k,l)]
    for j in range(L):
        I1234 += s1[j]*s2[j]
    return I1234
    
def I1234_full_bf4(pChain, m = 0, n = 0, print_b = 1):
    '''Computes the measure I1234_full@(m,n) by brute force that is, directly off
    the defining formula.'''
    I1234 = 0 #init
    pChain = getPolygonalChain(PDBfilename)[1]
    L = len(pChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    wDict = {}
    for i in range(m, n+1):
        for j in range(i+1, n):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    s1 = {}
    for j in range(m,n+1):
        s1[j] =0
        for i in range(m,j):
            s1[j] += wDict[(i,j)]
    s2 = {}
    for j in range(m,n+1):
        s2[j] =0
        for k in range(j+1, n+1):
#                wVal_ij = w(pChain[i],pChain[j])
            for l in range(k+1, n+1):
                s2[j] += wDict[(k,l)]
    for j in range(m,n+1):
        I1234 += s1[j]*s2[j]
    return I1234    
        

def I1324_bf(pChain):
    '''Computes the measure I1324 by brute force that is, directly off
    the defining formula.'''
    I1324 = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    cnt = 0
    for i in range(L):
        for k in range(i+1, L):
            for j in range(k+1, L):
#                wVal_ij = w(pChain[i],pChain[j])
                for l in range(j+1, L):
                    cnt +=1
                    if cnt%100000 == 0:
                        print cnt
                    wVal_ij = w(pChain[i],pChain[j])
                    wVal_kl = w(pChain[k],pChain[l])
                    I1324 += wVal_ij*wVal_kl
    print I1324
    return I1324
            
def I1324_bf2(pChain):
    '''Computes the measure I1324 by brute force that is, directly off
    the defining formula.'''
    I1324 = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    cnt = 0
    for i in range(L):
        for k in range(i+1, L):
            for j in range(k+1, L):
                wVal_ij = w(pChain[i],pChain[j])
                for l in range(j+1, L):
                    cnt +=1
                    if cnt%100000 == 0:
                        print cnt
                    wVal_kl = w(pChain[k],pChain[l])
                    I1324 += wVal_ij*wVal_kl
    return I1324


def I1324_bf3(pChain):
    '''Computes the measure I1324 by brute force that is, directly off
    the defining formula.'''
    I1324 = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for k in range(i+1, L):
            for j in range(k+1, L):
                wVal_ij = wDict[(i,j)]
                for l in range(j+1, L):
                    cnt +=1
#                    if cnt%100000 == 0:
#                        print cnt
                    I1324 += wVal_ij*wDict[(k,l)]
    return I1324
    
def I1324_full_bf3(pChain, m=0, n =0, print_b = 1):
    '''Computes the measure I1324_full@(m,n) by brute force that is, directly off
    the defining formula.'''
    I1324 = 0 #init
    L = len(pChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    wDict = {}
    for i in range(m,n+1):
        for j in range(i+1, n+1):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(m,n+1):
        for k in range(i+1, n+1):
            for j in range(k+1, n+1):
                wVal_ij = wDict[(i,j)]
                for l in range(j+1, n+1):
                    cnt +=1
                    if cnt%100000 == 0:
                        print cnt
                    I1324 += wVal_ij*wDict[(k,l)]
    if print_b > 0.5:
        print I1324
    return I1324

def I1324_full_bf_dict(PDBfilename):
    out = {}
    CaChain, pChain = getPolygonalChain(PDBfilename)

    for chain in pChain.keys():
        L = len(pChain[chain])
        if not(out.has_key(chain)):
            out[chain] = {}
        for m in range(L):
            for n in range(m,L):
                out[chain][(m,n)] = I1324_full_bf3(pChain[chain], m=m, n=n, print_b = 0)
    return out            

def I1324_full2_bf(pChain, i= 0, m = 0, n = 0, print_b = 1):
    '''Computes the measure I1324_full2(i; m , n) by brute force that is, directly off
    the defining formula. Default values: i=0, m=0, n=0.'''
    I = 0 #init
    L = len(pChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    wDict = {}
    for k in range(L):
        for l in range(k+1, L):
            wDict[(k,l)] = w(pChain[k],pChain[l])
    cnt = 0
    for b in range(m,n+1):
#        print "b-range"
#        print range(m,n+1)
#        print "b: %d" % b
        for d in range(b+1, n+1):
#            print "d-range"
#            print range(b+1,n+1)
#            print "d: %d" % d
            wVal_bd = wDict[(b,d)]
#            print "wVal_bd: %f" % wVal_bd
            for a in range(i, b):
#                print "a-range"
#                print range(i,b)
#                print "a: %d" % a
                for c in range(b+1, d):
#                    print "c-range"
#                    print range(b+1,d)
#                    print "c: %d" % c
                    cnt +=1
                    if cnt%10000000 == 0:
                        print cnt
                    I += wVal_bd*wDict[(a,c)]
    if print_b > 0.5:
        print "I: %f" % I
    return I    

def I1324_full2_bf_dict(PDBfilename):
    out = {}
    CaChain, pChain = getPolygonalChain(PDBfilename)

    for chain in pChain.keys():
        L = len(pChain[chain])
        if not(out.has_key(chain)):
            out[chain] = {}
        for i in range(L):
            for m in range(i,L): #obs: diagonal included on purpose
                out[chain][(i,m)] = I1324_full2_bf(pChain[chain], i=i, m=m, n=L-1, print_b = 0)
    return out 

def I1423_bf(pChain):
    '''Computes the measure I1423 by brute that is, directly off
    the defining formula.'''
    I1423 = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    cnt = 0
    for i in range(L):
        for k in range(i+1, L):
            for l in range(k+1, L):
#                wVal_ij = w(pChain[i],pChain[j])
                for j in range(l+1, L):
                    cnt +=1
                    if cnt%100000 == 0:
                        print cnt
                    wVal_kl = w(pChain[k],pChain[l])
                    wVal_ij = w(pChain[i],pChain[j])
                    I1423 += wVal_ij*wVal_kl
    print I1423
    return I1423  
    
    
def I1423_bf2(pChain):
    '''Computes the measure I1423 by brute force that is, directly off
    the defining formula.'''
    I1423 = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    cnt = 0
    for i in range(L):
        for k in range(i+1, L):
            for l in range(k+1, L):
                wVal_kl = w(pChain[k],pChain[l])
                for j in range(l+1, L):
                    cnt +=1
                    if cnt%100000 == 0:
                        print cnt
                    wVal_ij = w(pChain[i],pChain[j])
                    I1423 += wVal_ij*wVal_kl
    return I1423  


def I1423_bf3(pChain):
    '''Computes the measure I1423 by brute force that is, directly off
    the defining formula.'''
    I1423 = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for k in range(i+1, L):
            for l in range(k+1, L):
                wVal_kl = wDict([k,l])
                for j in range(l+1, L):
                    cnt +=1
                    if cnt%100000 == 0:
                        print cnt 
                    I1423 += wDict[(i,j)]*wVal_kl
    return I1423

def I1423_full_bf3(pChain, m=0, n=0, print_b = 1):
    '''Computes the measure I1423@(m,n) by brute force that is, directly off
    the defining formula.'''
    I1423 = 0 #init
    L = len(pChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    wDict = {}
    for i in range(m,n+1):
        for j in range(i+1, n+1):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(m, n+1):
        for k in range(i+1, n+1):
            for l in range(k+1, n+1):
                wVal_kl = wDict[(k,l)]
                for j in range(l+1, n+1):
                    cnt +=1
                    if print_b > 0.5:
                        if cnt%100000 == 0:
                            print cnt 
                    I1423 += wDict[(i,j)]*wVal_kl
    return I1423

def I1423_full_bf3_dict(PDBfilename):
    out = {}
    CaChain, pChain = getPolygonalChain(PDBfilename)

    for chain in pChain.keys():
        L = len(pChain[chain])
        if not(out.has_key(chain)):
            out[chain] = {}
        for m in range(L):
            for n in range(m+1,L):
                out[chain][(m,n)] = I1423_full_bf3(pChain[chain], m=m, n=n, print_b = 0)
    return out 
    
    
def I1423_bf4(pChain):
    '''Computes the measure I1423 by brute force that is, directly off
    the defining formula.'''
    I1423 = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    #to avoid a key-check below we set the value of w id to zero 
    #just next to the boundary
    for m in range(L):
        wDict[(m,L)] = 0
        wDict[(-1,m)] = 0
        wDict[(m+1,m)] = 0
    wDict[(-1,L)] = 0
    s1= {}
    for k in range(L):
        s1[(k,k)] = 0
        for j in range(k+1, L):
            s1[(k,j)] = s1[(k,j-1)] + wDict[(k,j)]
    s2 = {}
    for j in range(L):
        s2[(j-1,j)] = 0 #writhe is zero at and next to the diagonal
        for i in range(0,j-1)[::-1]:    
#            print (i,j)
#            print s1[(i,j)]
            s2[(i,j)] = s2[(i+1,j)] + s1[(i,j)]
            I1423 += wDict[(i-1,j+1)]*s2[(i,j)]
    return I1423
    
    
def I1423_full0_bf3(pChain,i=0, m=0, n=0, print_b = 1):
    '''Computes the measure I1423(i,m;n) by brute force that is, directly off
    the defining formula. Defualt values: i=0, m=0, n=0.'''
    I1423 = 0 #init
    L = len(pChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    wDict = {}
    for k in range(i,n+1):
        for l in range(k+1, n+1):
            wDict[(k,l)] = w(pChain[k],pChain[l])
    cnt = 0
    for b in range(i, m+1):
        for c in range(b+1, m+1):
            wVal_bc = wDict[(b,c)]
            for a in range(i, b):
                for d in range(c+1, n+1):
                    cnt +=1
                    if print_b > 0.5:
                        if cnt%100000 == 0:
                            print cnt 
                    I1423 += wDict[(a,d)]*wVal_bc
    return I1423

def I1423_full0_bf3_dict(PDBfilename):
    '''Generates dictionary of values I1423(i,m;n=L-1).'''
    out = {}
    CaChain, pChain = getPolygonalChain(PDBfilename)

    for chain in pChain.keys():
        L = len(pChain[chain])
        if not(out.has_key(chain)):
            out[chain] = {}
        for i in range(L):
            for m in range(i+1,L):
                out[chain][(i,m)] = I1423_full0_bf3(pChain[chain], i=i, m=m, n=L-1, print_b = 0)
    return out     
    
    
def I1423_full2_bf(pChain, i= 0, m = 0, n = 0, print_b = 1):
    '''Computes the measure I1423_full2(i; m , n) by brute force that is, directly off
    the defining formula. Default values: i=0, m=0, n=0.'''
    I = 0 #init
    L = len(pChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    wDict = {}
    for k in range(L):
        for l in range( k+1, L):
            wDict[(k,l)] = w(pChain[k],pChain[l])
    cnt = 0
    for b in range(m,n+1):
        for c in range(b+1, n+1):
            wVal_bc = wDict[(b,c)]
#            print "wVal_bc: %f" % wVal_bd
            for a in range(i, b):
                for d in range(c+1, n+1):
                    cnt +=1
                    if print_b > 0.5:
                        if cnt%100 == 0:
                            print cnt
                    I += wVal_bc*wDict[(a,d)]
    if print_b > 0.5:
        print "I: %f" % I
    return I   
 
def I1423_full2_bf_dict(PDBfilename):
    '''Generates dictionary of values I1423(i;m,n=L-1).'''
    out = {}
    CaChain, pChain = getPolygonalChain(PDBfilename)

    for chain in pChain.keys():
        L = len(pChain[chain])
        if not(out.has_key(chain)):
            out[chain] = {}
        for i in range(L):
            for m in range(i,L): #obs: diagonal included on purpose
                out[chain][(i,m)] = I1423_full2_bf(pChain[chain], i=i, m=m, n=L-1, print_b =0)
    return out
    

def I123456_bf(pChain):
    '''Computes the measure I123456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L): 
        for j in range(i+1, L):
            wVal_ij = wDict[(i,j)]
            for a in range(j+1, L):
                for b in range(a+1, L):
                    wVal_ab = wDict[(a,b)]
                    for c in range(b+1, L):
                        for d in range(c+1, L):
                            I += wVal_ij*wVal_ab*wDict[(c,d)]
                            cnt +=1
#                            if cnt%1000000 == 0:
#                                print cnt
    print "number of computations: %d" % cnt
    print "I: %f" % I
    return I


def I123546_bf(pChain):
    '''Computes the measure I123456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for j in range(i+1, L):
            wVal_ij = wDict[(i,j)]
            for a in range(j+1, L):
                for b in range(a+1, L):
                    for c in range(b+1, L):
                        wVal_ac = wDict[(a,c)]
                        for d in range(c+1, L):
                            I += wVal_ij*wVal_ac*wDict[(b,d)]
                            cnt +=1
#                            if cnt%1000000 == 0:
#                                print cnt
    print "I: %f" % I
    return I


def I123645_bf(pChain):
    '''Computes the measure I123456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for j in range(i+1, L):
            wVal_ij = wDict[(i,j)]
            for a in range(j+1, L):
                for b in range(a+1, L):
                    for c in range(b+1, L):
                        wVal_bc = wDict[(b,c)]
                        for d in range(c+1, L):
                            I += wVal_ij*wVal_bc*wDict[(a,d)]
                            cnt +=1
#                            if cnt%1000000 == 0:
#                                print cnt
    print "I: %f" % I
    return I

def I132456_bf(pChain):
    '''Computes the measure I132456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for a in range(i+1, L):
            for j in range(a+1, L):
                wVal_ij = wDict[(i,j)]
                for b in range(j+1, L):
                    wVal_ab = wDict[(a,b)]
                    for c in range(b+1, L):
                        for d in range(c+1, L):
                            I += wVal_ij*wVal_ab*wDict[(c,d)]
                            cnt +=1
#                            if cnt%1000000 == 0:
#                                print cnt
    print "I: %f" % I 
    return I
    
    
def I132546_bf(pChain, m=0, n=0, print_b = 1):
    '''Computes the measure I132456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    wDict = {}
    if n==0:
        n = L-1
    for i in range(m,L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(m,n+1):
        for a in range(i+1, L):
            for j in range(a+1, n+1):
                wVal_ij = wDict[(i,j)]
                for b in range(j+1, L):
                    for c in range(b+1, L):
                        wVal_ac = wDict[(a,c)]
                        for d in range(c+1, L):
                            I += wVal_ij*wVal_ac*wDict[(b,d)]
                            cnt +=1
#                            if cnt%1000000 == 0:
#                                print cnt
    if print_b > 0.5:
        print "I(%d,%d): %f" % (m,n,I)
    return I

def I132546_bf_dict(PDBfilename):
    out = {}
    CaChain, pChain = getPolygonalChain(PDBfilename)

    for chain in pChain.keys():
        L = len(pChain[chain])
        if not(out.has_key(chain)):
            out[chain] = {}
        for m in range(L):
            for n in range(m+1,L):
                out[chain][(m,n)] = I132546_bf(pChain[chain], m=m, n=n, print_b = 0)
    return out 

def I132645_bf(pChain):
    '''Computes the measure I132456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for a in range(i+1, L):
            for j in range(a+1, L):
                wVal_ij = wDict[(i,j)]
                for b in range(j+1, L):
                    for c in range(b+1, L):
                        wVal_bc = wDict[(b,c)]
                        for d in range(c+1, L):
                            I += wVal_ij*wVal_bc*wDict[(a,d)]
                            cnt +=1
#                            if cnt%1000000 == 0:
#                                print cnt
    print "I: %f" % I
    return I       

################################################################3    
#Experiment with use of jit:
from numba import autojit, jit, float64, int64, void
    
def I132645(pChain):
    I = 0 #init
    I = I132645_jit(pChain)
    print "I: %f" % I
    return I
    
#Leaving the followng line in or out makes aa huge difference in the time
#consumed. Try it eg on 1ENH (git.PDBfilename1enh)
#@jit(float64(float64[:,:,:]))
def I132645_jit(pChain):
    L = len(pChain)
    I = 0
    wDict = np.zeros((L+1,L+1))
    for i in range(L):
        for j in range(i+1, L):
            wDict[i,j] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for a in range(i+1, L):
            for j in range(a+1, L):
                wVal_ij = wDict[i,j]
                for b in range(j+1, L):
                    for c in range(b+1, L):
                        wVal_bc = wDict[b,c]
                        for d in range(c+1, L):
                            I += wVal_ij*wVal_bc*wDict[a,d]
                            cnt +=1
#                            if cnt%1000000 == 0:
#                                print cnt
#    print "I: %f" % I
    return I       
    
#Here only a double loop:
def I12(pChain):
    '''Computes the measure I1234 by brute force that is, directly off
    the defining formula.'''
    I = 0 #init
    I = I12_jit(pChain)
#    print "I: %f" % I
    return I

#@jit(float64(float64[:,:,:]))    
def I12_jit(pChain):    
    L = len(pChain)
    I = 0    
    cnt = 0
    wVal_ij = 0
    for i in range(L):
        for j in range(i+1, L):
            cnt +=1
#            if cnt%100000 == 0:
#                print cnt
            wVal_ij = w(pChain[i],pChain[j])
            I += wVal_ij
#    print "I12: %f" % I
    return I

### Experiment done
####################################
    
def I142356_bf(pChain):
    '''Computes the measure I132456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for a in range(i+1, L):
            for b in range(a+1, L):
                wVal_ab = wDict[(a,b)]
                for j in range(b+1, L):
                    wVal_ij = wDict[(i,j)]
                    for c in range(j+1, L):
                        for d in range(c+1, L):
                            I += wVal_ij*wVal_ab*wDict[(c,d)]
                            cnt +=1
#                            if cnt%1000000 == 0:
#                                print cnt
    print "I: %f" % I
    return I    
    
    
def I142536_bf(pChain):
    '''Computes the measure I132456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for a in range(i+1, L):
            for b in range(a+1, L):
                for j in range(b+1, L):
                    wVal_ij = wDict[(i,j)]
                    for c in range(j+1, L):
                        wVal_ac = wDict[(a,c)]
                        for d in range(c+1, L):
                            I += wVal_ij*wVal_ac*wDict[(b,d)]
                            cnt +=1
#                            if cnt%1000000 == 0:
#                                print cnt
    print "I: %f" % I
    return I    

def I142635_bf(pChain):
    '''Computes the measure I132456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for a in range(i+1, L):
            for b in range(a+1, L):
                for j in range(b+1, L):
                    wVal_ij = wDict[(i,j)]
                    for c in range(j+1, L):
                        wVal_bc = wDict[(b,c)]
                        for d in range(c+1, L):
                            I += wVal_ij*wVal_bc*wDict[(a,d)]
                            cnt +=1
#                            if cnt%1000000 == 0:
#                                print cnt
    print "I: %f" % I
    return I

def I152346_bf(pChain):
    '''Computes the measure I132456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for a in range(i+1, L):
            for b in range(a+1, L):
                wVal_ab = wDict[(a,b)]
                for c in range(b+1, L):
                    for j in range(c+1, L):
                        wVal_ij = wDict[(i,j)]
                        for d in range(j+1, L):
                            I += wVal_ij*wVal_ab*wDict[(c,d)]
                            cnt +=1
#                            if cnt%1000000 == 0:
#                                print cnt
    print "I: %f" % I
    return I 


def I152436_bf(pChain):
    '''Computes the measure I132456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for a in range(i+1, L):
            for b in range(a+1, L):
                for c in range(b+1, L):
                    wVal_ac = wDict[(a,c)]
                    for j in range(c+1, L):
                        wVal_ij = wDict[(i,j)]
                        for d in range(j+1, L):
                            I += wVal_ij*wVal_ac*wDict[(b,d)]
                            cnt +=1
#                            if cnt%1000000 == 0:
#                                print cnt
    print "I: %f" % I
    return I 


def I152634_bf(pChain):
    '''Computes the measure I132456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for a in range(i+1, L):
            for b in range(a+1, L):
                for c in range(b+1, L):
                    wVal_bc = wDict[(b,c)]
                    for j in range(c+1, L):
                        wVal_ij = wDict[(i,j)]
                        for d in range(j+1, L):
                            I += wVal_ij*wVal_bc*wDict[(a,d)]
                            cnt +=1
#                            if cnt%1000000 == 0:
#                                print cnt
    print "I: %f" % I
    return I 



def I162345_bf(pChain):
    '''Computes the measure I132456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for a in range(i+1, L):
            for b in range(a+1, L):
                wVal_ab = wDict[(a,b)]
                for c in range(b+1, L):
                    for d in range(c+1, L):
                        wVal_cd = wDict[(c,d)]
                        for j in range(d+1, L):
                            I += wVal_ab*wVal_cd*wDict[(i,j)]
                            cnt +=1
    #                            if cnt%1000000 == 0:
    #                                print cnt
    print "I: %f" % I
    return I
    
def I162435_bf(pChain):
    '''Computes the measure I132456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for a in range(i+1, L):
            for b in range(a+1, L):
                for c in range(b+1, L):
                    wVal_ac = wDict[(a,c)]
                    for d in range(c+1, L):
                        wVal_bd = wDict[(b,d)]
                        for j in range(d+1, L):
                            I += wVal_ac*wVal_bd*wDict[(i,j)]
                            cnt +=1
    #                            if cnt%1000000 == 0:
    #                                print cnt
    print "I: %f" % I
    return I
    
    
def I162534_bf(pChain):
    '''Computes the measure I132456 by brute force that is, directly off
    Ithe defining formula.'''
    I = 0 #init
    L = len(pChain)
    print "Length of chain is: %d" % L
    wDict = {}
    for i in range(L):
        for j in range(i+1, L):
            wDict[(i,j)] = w(pChain[i],pChain[j])
    cnt = 0
    for i in range(L):
        for a in range(i+1, L):
            for b in range(a+1, L):
                for c in range(b+1, L):
                    wVal_bc = wDict[(b,c)]
                    for d in range(c+1, L):
                        wVal_ad = wDict[(a,d)]
                        for j in range(d+1, L):
                            I += wVal_bc*wVal_ad*wDict[(i,j)]
                            cnt +=1
    #                            if cnt%1000000 == 0:
    #                                print cnt
    print "I: %f" % I
    return I