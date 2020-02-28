# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 13:01:03 2014

@author: Christian Grønbæk

Pre-cursor to GISA C-code (github.com/ceegeeCode/GISA). This module contains, among much else, the computation of the Gauss Integrals by means
of a recursion (derived for the purpose); the very same recursion is used in the GISA C-code.

If you use this code 

0) it'll be "as is" -- I will have no time for support. The alternative is to use the GISA C-code (github repo: see above)

1) Cite it by referring

    a) the Github repository ceegeeCode/GISA
    and
    b) the paper "C.Grønbæk., T.Hamelryck, P.Røgen: "GISA: Using Gauss integrals to identify rare conformations in protein structures", TO APPEAR. 

Once the paper is published, cite it with ref to that publication.  

"""

'''Usage:

import gitP as git

### Compute measures:
#Using dictionary based version:
outDict_dict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1bio, degree =6, inNumpy_b =0, print_b = 0)

#Using array based version:
outDict, CaChain, polyCaChain = git.I(PDBfilename = git.PDBfilename1bio, degree =6, inNumpy_b =0, print_b = 0)

### Compute measures on perturbed chain:
#Using dictionary based version:
wDict  = {}
for chain in outDict_dict.keys():
    wDict[chain] = outDict_dict[chain]['w'] 
perOutDict_dict, perCaChain, perPolyCaChain, segNrRange = git.Iperturbed_dict(CaChain = CaChain, polyCaChain =polyCaChain, wDict = wDict, residueNrRange = range(20,25), degree = 6)

#Using array based version:
wDict  = {}
for chain in outDict.keys():
    wDict[chain] = outDict[chain][0] 
    
perOutDict, perCaChain, perPolyCaChain, segNrRange = git.Iperturbed(CaChain = CaChain, polyCaChain =polyCaChain, wDict = wDict, residueNrRange = range(20,25), degree = 6)

### Compute measures using vectorized version:
#dictionary based (vect_b =1!):
outDict_dict, CaChain, polyCaChain = git.I_dict(degree =6, inNumpy_b =0, vect_b =1)

#array based:
outDict, CaChain, polyCaChain = git.I(degree =6, inNumpy_b =0, vect_b = 1)

### Compute measures on perturbed chain using vectorized versions:
#Using dictionary based version:
outDict_dict, CaChain, polyCaChain = git.I_dict(degree =6, inNumpy_b =0, vect_b =1)
wDict  = {}
for chain in outDict_dict.keys():
    wDict[chain] = outDict_dict[chain]['w'] 
perOutDict_dict, perCaChain, perPolyCaChain, segNrRange = git.Iperturbed_dict(CaChain = CaChain, polyCaChain =polyCaChain, wDict = wDict, residueNrRange = range(20,25), degree = 6, vect_b =1)

#Using array based version:
outDict, CaChain, polyCaChain = git.I(degree =6, inNumpy_b =0, vect_b = 1)
wDict  = {}
for chain in outDict.keys():
    wDict[chain] = outDict[chain][0] 
perOutDict, perCaChain, perPolyCaChain, segNrRange = git.Iperturbed(CaChain = CaChain, polyCaChain =polyCaChain, wDict = wDict, residueNrRange = range(20,25), degree = 6, vect_b = 1)


### Running multi perturbation functions

#Dictionary based:

outDict_dict, CaChain, polyCaChain = git.I_dict(degree =0, inNumpy_b =0, vect_b =1 )
 
wDict  = {}
for chain in outDict_dict.keys():
    wDict[chain] = outDict_dict[chain]['w'] 
    
tot_outDict, perturbedChains = git.Iperturbed_vect_multi_dict(CaChain = CaChain,  
               polyCaChain = polyCaChain, 
               wDict = wDict,
               nrPerturbations = 10,
               residueNrRange = range(0,5),
               perturbationLength = 5, 
                degree = 4,
                full_b = 0,
                print_b = 0)

#Array based:

outDict, CaChain, polyCaChain = git.I(degree =4, inNumpy_b =0, vect_b = 1)

wDict  = {}
for chain in outDict.keys():
    wDict[chain] = outDict[chain][0] 
    
tot_outDict, perturbedChains = git.Iperturbed_vect_multi(CaChain = CaChain,  
               polyCaChain = polyCaChain, 
               wDict = wDict,
               nrPerturbations = 10,
               residueNrRange = range(0,5),
               perturbationLength = 5, 
                degree = 4,
                full_b = 0,
                print_b = 0)

#####################################################################################

THESE PARTS HAVE NOT BEEN REDONE YET, AFTER THE CHANGE OF LOADING THE STRUCTURE THE 
MULTI-SUB-CHAIN WAY

############################# Timing tests:
 
import timeit

### Comparisons with "brute force" computations:
t_I1234_bf3 = timeit.Timer('git.I1234_bf3(git.PDBfilenameDef)','import git')
t_I1234_bf3.timeit(number=1)

t_I1234_bf4 = timeit.Timer('git.I1234_bf4(git.PDBfilenameDef)','import git')
t_I1234_bf4.timeit(number=1)

t_I1324_bf3 = timeit.Timer('git.I1324_bf3(git.PDBfilenameDef)','import git')
t_I1324_bf3.timeit(number=1)

t_I1423_bf4 = timeit.Timer('git.I1423_bf4(git.PDBfilenameDef)','import git')
t_I1423_bf4.timeit(number=1)

#Clearly much slower than using the I and I_dict functions below.

### How much time does the aggregations consume?:

t_I = timeit.Timer('git.I(full_b = 0, degree =0, print_b = 0)','import git')
t_I.timeit(number=1)

t_I = timeit.Timer('git.I(full_b = 0, degree =2, print_b = 0)','import git')
t_I.timeit(number=1)

t_I = timeit.Timer('git.I(full_b = 0, print_b=0, degree =4,)','import git')
t_I.timeit(number=1)

t_I = timeit.Timer('git.I(full_b = 0, print_b=0, degree =6, vect_b = 0)','import git')
t_I.timeit(number=1)

#writhe comp, vect'ed
t_I = timeit.Timer('git.I(full_b = 0, print_b=0, degree =0, vect_b =1)','import git')
t_I.timeit(number=1)


#############################################################################
#jit expermient on speeding up for loops:
############################################################################
#Time gained by using jit to get for loops run in C:
#Obs: the times reported stem from runs done wo power cord plugged in (ie on
#battery usage along w some energy saver). If power cord is plugged in times
#reduce by a factor ~ 4.
#first we run a degree 6 measure on 1ENH wo invoking jit:
t_I = timeit.Timer('git_bf.I132645(PDBfilename = git.PDBfilename1enh)','import git, git_bf')
t_I.timeit(number=1)
#Result:
I: 11.055975 (value)
52.595803961710075 (time)
#Now we run the same again but w jit invoked (takes uncommenting the jit-line right 
#above the fct I132645_jit in git_bf and reloading the git_bf-module)
#Here's the result:
I: 11.055975
1.9628975062756524
Ie: a factor about 30 faster!

#Next let's see if there's anything gained in a 2d-loop:
#first wo jit:
t_I = timeit.Timer('git_bf.I12(PDBfilename = git.PDBfilename1enh)','import git, git_bf')
t_I.timeit(number=10)
#Result:
I: 30.136107
12.70625758181643

#Now with jit:
t_I = timeit.Timer('git_bf.I12(PDBfilename = git.PDBfilename1enh)','import git, git_bf')
t_I.timeit(number=10)
#Result:
I: 30.136107
12.809605561979197
#Nothing is gained! On 1ads (lgth ~300) then:
t_I = timeit.Timer('git_bf.I12_bf(PDBfilename = git.PDBfilename1ads, m=0, n=300)','import git, git_bf')
t_I.timeit(number=1)
#Results wo jit:
I: 119.840333
44.34628349371678
#Result w jit:
I: 119.840333
44.602786831985895
#Nothing gained here either.


#jit experiment done
############################################################################

#Which one is faster: the dictionary based version or the one based on np-arrays?

t_I = timeit.Timer('git.I(full_b = 1, degree =4)','import git')
t_I.timeit(number=1)

t_I = timeit.Timer('git.I_dict(full_b = 1, degree =4)','import git')
t_I.timeit(number=1)


##########################################################################
# Other usage:

#vectorizing:
t_wAll = timeit.Timer('git.wAll(print_b=0)','import git')
t_wAll.timeit(number=10)

t_I_vect = timeit.Timer('git.I(full_b = 0, print_b=0, degree =4, vect_b = 1)','import git')
t_I_vect.timeit(number=1)

#On lgth 5 pertubation:
t = timeit.Timer('git.timeTest(residueNrRange =range(25,30),perturbation = [(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3)],full_b = 0, print_b = 0, degree = 4)','import git')
t.timeit(number=10)  

t = timeit.Timer('git.timeTest(residueNrRange =range(25,30),perturbation = [(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3)], full_b = 0, print_b = 0, degree = 4, vect_b =1)','import git')
t.timeit(number=10)

#On lgth 10 pertubation:
t = timeit.Timer('git.timeTest(residueNrRange =range(20,30),perturbation = [(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3),(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3)],full_b = 0, print_b = 0, degree = 4)','import git')
t.timeit(number=10)  

t = timeit.Timer('git.timeTest(residueNrRange =range(20,30),perturbation = [(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3),(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3)], full_b = 0, print_b = 0, degree = 4, vect_b =1)','import git')
t.timeit(number=10)


git.timeTest_vect_multi(nrPerturbations = 1, residueNrRange = range(50,60), perturbationLength = 10, degree =4, full_b = 0, print_b =0)

lgth 5 pertubation:
t = timeit.Timer('git.timeTest_vect_multi(nrPerturbations = 10, residueNrRange = range(50,60), perturbationLength = 10, degree =4, full_b = 0, print_b =0)','import git')
t.timeit(number=1)

lgth 10 pertubation:
t = timeit.Timer('git.timeTest_vect_multi(nrPerturbations = 100, residueNrRange = range(20,30), perturbationLength = 10, degree =4, full_b = 0, print_b =0)','import git')
t.timeit(number=1)


#much other/longer chains:
#lght ~100:
PDBfilenameDef = PDBfilename2mhr
#lgth ~300:
PDBfilenameDef = PDBfilename1ads

#For these change the PDBfilenameDef accordingly below and reload the module or
#specify the desired default, PDBfilenameDef, explicitly in the function call

#running the perturbation computation in steps or using the function:
#get chains
C, PC = git.getPolygonalChain(git.PDBtop100Folder + '1adsH' )

#compute the measures
outDict, CaChain, pChain =git.I(degree=4)

#compute perturbed chain:
residueNrRange = range(10,15)
perturbation = [(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3)]
perC, perPC, segNrRange = git.perturbPolygonalChain(CaChain = C, 
                                                    polyCaChain = PC,
                     residueNrRange = residueNrRange, 
                     perturbation = perturbation)  
 
#compute measures on perturbed chain directly:
outDict, perCaChain, perPolyCaChain = git.I(CaChain = perC, polyCaChain = perPC , print_b = 1, degree =2)

#compute the measures on perturbed chain using the specific function: 
outDict_per, perCaChain, perPolyCaChain, segNrRange = git.Iperturbed(CaChain = C, polyCaChain = PC, wDict = outDict[0], residueNrRange =residueNrRange , perturbation = perturbation, degree = 2, print_b = 1)


*************************************************************************************************
*************************************************************************************************
END OF: THESE PARTS HAVE NOT BEEN REDONE YET, AFTER THE CHANGE OF LOADING THE STRUCTURE THE 
MULTI-SUB-CHAIN WAY
*************************************************************************************************
*************************************************************************************************


####################################################################################
## Checking results
####################################################################################

#Do array based version of I give the same output as the dictionary based one?
#degree = 2
git.compareIversions(PDBfilename = git.PDBfilename1not, degree = 2, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1enh, degree = 2, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1ads, degree = 2, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1bio, degree = 2, tolerance = 1e-9)


#degree = 4:
git.compareIversions(PDBfilename = git.PDBfilename1not, degree = 4, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1enh, degree = 4, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1ads, degree = 4, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1bio, degree = 4, tolerance = 1e-9)


#degree = 6:
git.compareIversions(PDBfilename = git.PDBfilename1not, degree = 6, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1enh, degree = 6, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1ads, degree = 6, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1bio, degree = 6, tolerance = 1e-9)

#Results: all ok.

#Checking results (dictionary based version) vs brute force versions:

#Re length of chain: in higher degrees the brute force versions will be too time consuming for longer
#chains; short chains: with 1not the length is L = 12; for 1mctl length is L = 27

#set default PDB-filename: e.g. change the definition in the code below and reload. 
#then compute and check results vs brute force computations ('bf'):
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilenameDef, degree =6, inNumpy_b =0 )

#to use the easily the "old code" for comparison of the output dictionaries we "transpose" the output
#from I_dict (and I) so that outer key will be the invariant name and inner key the chain id:
outDictTr = git.transposeDict(outDict)
#we use compareDictsII for the comparisons, since we have both an outer and an inner key (ie results are
#dicts of dicts)

#I12:
I12 = git_bf.I12_bf_dict(PDBfilename = git.PDBfilenameDef)
#compare
git.compareDictsII(I12, outDictTr['I12'], tolerance = 1e-6)
#result: ([[]], {'X': {}})

#for a two chain structure:
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1bio, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)
#bf:
I12 = git_bf.I12_bf_dict(PDBfilename = git.PDBfilename1bio)
#compare
git.compareDictsII(I12, outDictTr['I12'], tolerance = 1e-10)
([[], []], {'A': {}, 'B': {}})
#However, this runs too slow to live with; computing I12 by bf at all vertices for the 1bio structure is
#too much; so  we create an artificial structure consisting of two short chains (git.PDBfilenameArt)
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1art, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)
#bf:
I12 = git_bf.I12_bf_dict(PDBfilename = git.PDBfilename1art)
#compare
git.compareDictsII(I12, outDictTr['I12'], tolerance = 1e-10)
([[], []], {'A': {}, 'B': {}})

#single chain case:
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1not, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)
#bf:
I12 = git_bf.I12_bf_dict(PDBfilename = git.PDBfilename1not)
#compare
git.compareDictsII(I12, outDictTr['I12'], tolerance = 1e-10)
#result: ([[]], {'X': {}})

#Obs: the compareDicts function returns differences in keys of the two dictionaries
#if the (absolute value of the) value at the given key is larger than the set tolerance.
#Here the outDictTr['I12'][some chain id] dictionary actually has more keys than the I12[some chain id] 
#dictionary; no errors results since all values at the "extra" keys are below the tolerance (they are
#zero, in fact, and stem from an initialization in the I_dict function). The same comment
#applies to the cases below. 

#If wanted one can check the congruence of the results at some vertices (but if for
#ascertaining that non-zeroes do exist it is more obvious to simply look at the
#whole dictionary):
#some vertex:
#first get a chain id:
chainId = outDict.keys()[0]
outDict[chainId]['I12'][(5,7)]
I12 = git_bf.I12_bf_dict(PDBfilename = git.PDBfilename1art)
I12[chainId][(5,7)]
#whole chain result:
L = 9
outDict[chainId]['I12'][(0,L-1)]
I12[chainId][(0,L-1)]

#degree 4 measures:

#I1234_full:
I1234_full = git.I1234_full_bf3_dict(PDBfilename = git.PDBfilename1art)
#compare
git.compareDictsII(I1234_full, outDictTr['I1234_full'], tolerance = 1e-6)
#result: ([[], []], {'A': {}, 'B': {}})

#single chain case:
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1not, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)
I1234_full = git.I1234_full_bf3_dict(PDBfilename = git.PDBfilename1not)
#compare
git.compareDictsII(I1234_full, outDictTr['I1234_full'], tolerance = 1e-6)
#result: ([[]], {'X': {}})

#I1234, value:
CaChain, polyCaChain =  git_utils.getPolygonalChain(PDBfilename = git.PDBfilename1art)
chainId = CaChain.keys()[1]
I1234 = git_bf.I1234_bf(pChain = polyCaChain[chainId])
#compare
L = 9
outDict[chainId]['I1234'][(0,L-1)]
#result: fine

#I1324, value:
CaChain, polyCaChain =  git_utils.getPolygonalChain(PDBfilename = git.PDBfilename1art)
chainId = CaChain.keys()[1]
#compare
I1324 = git_bf.I1324_bf(pChain = polyCaChain[chainId])
L = 9
outDict[chainId]['I1324'][(0,L-1)]
#result: fine

#I1324_full:
I1324_full = git_bf.I1324_full_bf_dict(PDBfilename = git.PDBfilename1art)
#compare
git.compareDictsII(I1324_full, outDictTr['I1324_full'], tolerance = 1e-6)
#result: ([[], []], {'A': {}, 'B': {}})

#I1324_full2:
I1324_full2 = git_bf.I1324_full2_bf_dict(PDBfilename = git.PDBfilename1art)
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1art, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)
#compare
git.compareDictsII(I1324_full2, outDictTr['I1324_full2'], tolerance = 1e-6)
#result: ([[], []], {'A': {}, 'B': {}})
#same result for PDBfilename = git.PDBfilename1not (only 1 chain though)

#Check that I1324_full[(j,L-1)] = I1324_full2[(j,j+1)]:
CaChain, polyCaChain =  git_utils.getPolygonalChain(PDBfilename = git.PDBfilename1art)

chainId = CaChain.keys()[1]
L = len(polyCaChain[chainId])
x = [I1324_full[chainId][(j,L-1)] - I1324_full2[chainId][(j,j+1)] for j in range(L-1)]
x
#Result: OK; all values in x are (essentially) zero (e.g. below 1e-12 in absolute size); holds for both chains

##I1324_full and I1324_full2, single chain case:

outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1not, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)
I1324_full = git_bf.I1324_full_bf_dict(PDBfilename = git.PDBfilename1not)
#compare
git.compareDictsII(I1324_full, outDictTr['I1324_full'], tolerance = 1e-6)
#result: ([[]], {'X': {}})

#I1324_full2:
I1324_full2 = git_bf.I1324_full2_bf_dict(PDBfilename = git.PDBfilename1not)
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1not, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)
#compare
git.compareDictsII(I1324_full2, outDictTr['I1324_full2'], tolerance = 1e-6)
#result: ([[]], {'X': {}})

#Check that I1324_full[(j,L-1)] = I1324_full2[(j,j+1)]:
CaChain, polyCaChain =  git_utils.getPolygonalChain(PDBfilename = git.PDBfilename1not)

chainId = CaChain.keys()[0]
L = len(polyCaChain[chainId])
x = [I1324_full[chainId][(j,L-1)] - I1324_full2[chainId][(j,j+1)] for j in range(L-1)]
x
#Result: OK; all values in x are (essentially) zero (e.g. below 1e-12 in absolute size)



#I1423, value:
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1art, degree =6, inNumpy_b =0 )
chainId = CaChain.keys()[0]
L = len(polyCaChain[chainId])
#compare
I1423 = git_bf.I1423_bf(pChain = polyCaChain[chainId])
#compare
outDict[chainId]['I1423'][(0,L-1)]
#result: fine

#I1423, all values:
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1art, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)

I1423_full = git_bf.I1423_full_bf3_dict(PDBfilename = git.PDBfilename1art)
#compare
git.compareDictsII(I1423_full, outDictTr['I1423'], tolerance = 1e-6)
#result: ([[], []], {'A': {}, 'B': {}})

#same, on single chain structure:
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1not, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)

I1423_full = git_bf.I1423_full_bf3_dict(PDBfilename = git.PDBfilename1not)
#compare
git.compareDictsII(I1423_full, outDictTr['I1423'], tolerance = 1e-6)
#result: ([[]], {'X': {}})


#I1423_full0:
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1art, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)

I1423_full0 = git_bf.I1423_full0_bf3_dict(PDBfilename = git.PDBfilename1art)
#compare
git.compareDictsII(I1423_full0, outDictTr['I1423_full0'], tolerance = 1e-6)
#result: ([[], []], {'A': {}, 'B': {}})

#same, on single chain structure:
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1not, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)

I1423_full0 = git_bf.I1423_full0_bf3_dict(PDBfilename = git.PDBfilename1not)
#compare
git.compareDictsII(I1423_full0, outDictTr['I1423_full0'], tolerance = 1e-6)
#result: ([[]], {'X': {}})


#I1423_full2:
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1art, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)

I1423_full2 = git_bf.I1423_full2_bf_dict(PDBfilename = git.PDBfilename1art)
#compare
git.compareDictsII(I1423_full2, outDictTr['I1423_full2'], tolerance = 1e-6)
#result: ([[], []], {'A': {}, 'B': {}})

#same, on single chain structure:
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1not, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)

I1423_full2 = git_bf.I1423_full2_bf_dict(PDBfilename = git.PDBfilename1not)
#compare
git.compareDictsII(I1423_full2, outDictTr['I1423_full2'], tolerance = 1e-6)
#result: ([[]], {'X': {}})


#degree 6 measures:

#single chain case:
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1not, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)

for chainId in CaChain.keys():
    
    print "Chain id: %s" % chainId

    L = len(polyCaChain[chainId])
    
    pChain = polyCaChain[chainId]
    #12*'s:
    I = git_bf.I123456_bf(pChain)
    outDict[chainId]['I123456'][(0,L-1)] #indices: (0,lgth of chain -1); change acc'ly when changing pdb-file
    
    I = git_bf.I123546_bf(pChain)
    outDict[chainId]['I123546'][(0,L-1)]
    
    I = git.I123645_bf(pChain)
    outDict[chainId]['I123645'][(0,L-1)]
    
    #Results: all ok.
    
    #13*'s:
    I = git.I132456_bf(pChain)
    outDict[chainId]['I132456'][(0,L-1)]
    
    I = git.I132546_bf(pChain)
    outDict[chainId]['I132546'][(0,L-1)]
    
    I = git.I132645_bf(pChain)
    outDict[chainId]['I132645'][(0,L-1)]
    
    #Results: all ok.
    
    #14*'s:
    I = git.I142356_bf(pChain)
    outDict[chainId]['I142356'][(0,L-1)]
    
    I = git.I142536_bf(pChain)
    outDict[chainId]['I142536'][(0,L-1)]
    
    I = git.I142635_bf(pChain)
    outDict[chainId]['I142635'][(0,L-1)]
    
    #Results: all ok.
    
    #15*'s:
    I = git.I152346_bf(pChain)
    outDict[chainId]['I152346'][(0,L-1)]
    
    I = git.I152436_bf(pChain)
    outDict[chainId]['I152436'][(0,L-1)]
    
    I = git.I152634_bf(pChain)
    outDict[chainId]['I152634'][(0,L-1)]
    
    #Results: all ok.
    
    #16*'s:
    I = git.I162345_bf(pChain)
    outDict[chainId]['I162345'][(0,L-1)]
    
    I = git.I162435_bf(pChain)
    outDict[chainId]['I162435'][(0,L-1)]
    
    I = git.I162534_bf(pChain)
    outDict[chainId]['I162534'][(0,L-1)]
    

#Results: all ok.

#for I132546 there's also a dictionary bf-version:
I132546 = git_bf.I132546_bf_dict(PDBfilename = git.PDBfilename1not)
#compare
git.compareDictsII(I132546, outDictTr['I132546'], tolerance = 1e-6)
#result: ([[]], {'X': {}})



##two-chains case:
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1art, degree =6, inNumpy_b =0 )
outDictTr = git.transposeDict(outDict)

#The run loop above for single chain case:

for chainId in CaChain.keys():
    
    print "Chain id: %s" % chainId

    L = len(polyCaChain[chainId])
    
    pChain = polyCaChain[chainId]
    #12*'s:
    I = git_bf.I123456_bf(pChain)
    outDict[chainId]['I123456'][(0,L-1)] #indices: (0,lgth of chain -1); change acc'ly when changing pdb-file
    
    I = git_bf.I123546_bf(pChain)
    outDict[chainId]['I123546'][(0,L-1)]
    
    I = git.I123645_bf(pChain)
    outDict[chainId]['I123645'][(0,L-1)]
    
    #Results: all ok.
    
    #13*'s:
    I = git.I132456_bf(pChain)
    outDict[chainId]['I132456'][(0,L-1)]
    
    I = git.I132546_bf(pChain)
    outDict[chainId]['I132546'][(0,L-1)]
    
    I = git.I132645_bf(pChain)
    outDict[chainId]['I132645'][(0,L-1)]
    
    #Results: all ok.
    
    #14*'s:
    I = git.I142356_bf(pChain)
    outDict[chainId]['I142356'][(0,L-1)]
    
    I = git.I142536_bf(pChain)
    outDict[chainId]['I142536'][(0,L-1)]
    
    I = git.I142635_bf(pChain)
    outDict[chainId]['I142635'][(0,L-1)]
    
    #Results: all ok.
    
    #15*'s:
    I = git.I152346_bf(pChain)
    outDict[chainId]['I152346'][(0,L-1)]
    
    I = git.I152436_bf(pChain)
    outDict[chainId]['I152436'][(0,L-1)]
    
    I = git.I152634_bf(pChain)
    outDict[chainId]['I152634'][(0,L-1)]
    
    #Results: all ok.
    
    #16*'s:
    I = git.I162345_bf(pChain)
    outDict[chainId]['I162345'][(0,L-1)]
    
    I = git.I162435_bf(pChain)
    outDict[chainId]['I162435'][(0,L-1)]
    
    I = git.I162534_bf(pChain)
    outDict[chainId]['I162534'][(0,L-1)]

#Results: all ok.

#for I132546 there's also a dictionary bf-version:
I132546 = git_bf.I132546_bf_dict(PDBfilename = git.PDBfilename1art)
#compare
git.compareDictsII(I132546, outDictTr['I132546'], tolerance = 1e-6)
#result: ([[], []], {'A': {}, 'B': {}})

#Results all in all: fine for all measures, both single and two-chain cases.

#### Checks of perturbation parts

#Do array based version of Ipertubed give the same output as the dictionary based one?
#degree = 2
git.compareIversions(PDBfilename = git.PDBfilename1not, degree = 2, pert_b =1, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1enh, degree = 2, pert_b =1, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1ads, degree = 2, pert_b =1, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1bio, degree = 2, pert_b =1, tolerance = 1e-9)
 
#degree = 4
git.compareIversions(PDBfilename = git.PDBfilename1not, degree = 4, pert_b =1, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1enh, degree = 4, pert_b =1, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1ads, degree = 4, pert_b =1, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1bio, degree = 4, pert_b =1, tolerance = 1e-9)

#degree = 6
git.compareIversions(PDBfilename = git.PDBfilename1not, degree = 6, pert_b =1, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1enh, degree = 6, pert_b =1, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1ads, degree = 6, pert_b =1, tolerance = 1e-9)
git.compareIversions(PDBfilename = git.PDBfilename1bio, degree = 6, pert_b =1, tolerance = 1e-9)

#Results: all ok.


#### Check results of vectorized versions vs non-vectorized:
#Dictionary based:
#degree = 2
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1art, degree =2, inNumpy_b =0 )
outDict_vect, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1art, degree =2, inNumpy_b =0, vect_b = 1 )

#identical output?
outDict == outDict_vect
#Result is "False" but differences are tiny: we compare dictionaries of dictionaries 
for chainId in CaChain.keys():
    print "Chain id: %s" % chainId
    git.compareDictsII(outDict[chainId], outDict_vect[chainId], tolerance = 1e-6)
#Result: for each chain this: ([[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []], {'I142356': {}, 'I162345': {}, 'I152436': {}, 'I132546': {}, 'I132456': {}, 'I123546': {}, 'I1324_full2': {}, 'I1423_full0': {}, 'I1423_full2': {}, 'I123645': {}, 'I123456': {}, 'I162534': {}, 'I1324_full': {}, 'I1324': {}, 'I12': {}, 'I142536': {}, 'I152346': {}, 'I142635': {}, 'I152634': {}, 'I1423': {}, 'I1234': {}, 'I162435': {}, 'w': {}, 'I1234_full': {}, 'I132645': {}})

#degree = 4
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1art, degree =4, inNumpy_b =0 )
outDict_vect, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1art, degree =4, inNumpy_b =0, vect_b = 1 )

#identical output?
outDict == outDict_vect
#Result is "False" but differences are tiny: we compare dictionaries of dictionaries 
for chainId in CaChain.keys():
    print "Chain id: %s" % chainId
    git.compareDictsII(outDict[chainId], outDict_vect[chainId], tolerance = 1e-6)
#Result: for each chain this: ([[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []], {'I142356': {}, 'I162345': {}, 'I152436': {}, 'I132546': {}, 'I132456': {}, 'I123546': {}, 'I1324_full2': {}, 'I1423_full0': {}, 'I1423_full2': {}, 'I123645': {}, 'I123456': {}, 'I162534': {}, 'I1324_full': {}, 'I1324': {}, 'I12': {}, 'I142536': {}, 'I152346': {}, 'I142635': {}, 'I152634': {}, 'I1423': {}, 'I1234': {}, 'I162435': {}, 'w': {}, 'I1234_full': {}, 'I132645': {}})

degree = 6:
outDict, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1art, degree =6, inNumpy_b =0 )
outDict_vect, CaChain, polyCaChain = git.I_dict(PDBfilename = git.PDBfilename1art, degree =6, inNumpy_b =0, vect_b = 1 )

#identical output?
outDict == outDict_vect
#Result is "False" but differences are tiny: we compare dictionaries of dictionaries 
for chainId in CaChain.keys():
    print "Chain id: %s" % chainId
    git.compareDictsII(outDict[chainId], outDict_vect[chainId], tolerance = 1e-6)
#Result: for each chain this: ([[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []], {'I142356': {}, 'I162345': {}, 'I152436': {}, 'I132546': {}, 'I132456': {}, 'I123546': {}, 'I1324_full2': {}, 'I1423_full0': {}, 'I1423_full2': {}, 'I123645': {}, 'I123456': {}, 'I162534': {}, 'I1324_full': {}, 'I1324': {}, 'I12': {}, 'I142536': {}, 'I152346': {}, 'I142635': {}, 'I152634': {}, 'I1423': {}, 'I1234': {}, 'I162435': {}, 'w': {}, 'I1234_full': {}, 'I132645': {}})

### This holds for single chain case too (1not was cnonsidered). All results ok

#Array based, same exercise:
#degree = 2
outDict, CaChain, polyCaChain = git.I(PDBfilename = git.PDBfilename1art, degree =2, inNumpy_b =0 )
outDict_vect, CaChain, polyCaChain = git.I(PDBfilename = git.PDBfilename1art, degree =2, inNumpy_b =0, vect_b = 1 )

#identical output?
outDict == outDict_vect
#Cannot be done, array comparison ...; we compare arrays of arrays: 
for chainId in CaChain.keys():
    print "Chain id: %s" % chainId
    git.compareArrays(outDict[chainId], outDict_vect[chainId], tolerance = 1e-6)
#Result: for both chains: 
length a1: 25 a2: 25
[]

#degree = 4
outDict, CaChain, polyCaChain = git.I(PDBfilename = git.PDBfilename1art, degree =4, inNumpy_b =0 )
outDict_vect, CaChain, polyCaChain = git.I(PDBfilename = git.PDBfilename1art, degree =4, inNumpy_b =0, vect_b = 1 )

#identical output?
outDict == outDict_vect
#Cannot be done, array comparison ...; we compare arrays of arrays: 
for chainId in CaChain.keys():
    print "Chain id: %s" % chainId
    git.compareArrays(outDict[chainId], outDict_vect[chainId], tolerance = 1e-6)
#Result: for both chains: 
length a1: 25 a2: 25
[]


#degree = 6:

outDict, CaChain, polyCaChain = git.I(PDBfilename = git.PDBfilename1art, degree =6, inNumpy_b =0 )
outDict_vect, CaChain, polyCaChain = git.I(PDBfilename = git.PDBfilename1art, degree =6, inNumpy_b =0, vect_b = 1 )

#identical output?
outDict == outDict_vect
#Cannot be done, array comparison ...; we compare arrays of arrays: 
for chainId in CaChain.keys():
    print "Chain id: %s" % chainId
    git.compareArrays(outDict[chainId], outDict_vect[chainId], tolerance = 1e-6)
#Result: for both chains: 
length a1: 25 a2: 25
[]


### This holds for single chain case too (1not was cnonsidered). All results ok


######################################################################

i=0
r = ''
for chain in model:
     for residue in chain:
         if PDB.is_aa(residue):
             x = residue.id[0]
             i+=1
             if x == ' ':
                 print residue, x, i
                 #r = residue

for chain in model:
     for residue in chain:
         if PDB.is_aa(residue):
             print residue['CA']

'''

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats #not used in the first part (the Pyhton impl of the GITs via recursions)

#from numba import int32, uint32, float32, float64
from numba import autojit, jit, int64,float64, double, void
from numbapro import vectorize, cuda

from pymol import cmd

#needed in C-comparison section
import re

#import scipy.stats as stats

#ignore warnings
import warnings
warnings.filterwarnings("ignore")

ErrorFile = r"C:/Users/Kristian/BioInformatics/projects/knots/code/errors_git.txt"

#for timing
#from timeit import *
import time
import timeit

#import git modules:
from git_utils import *
from git_bf import *
#from git_test import *

from Bio import PDB

parser = PDB.PDBParser()

PDBtop100Folder = r"C:/Users/Christian/BioInformatics/projects/knots/code/top100H/"
PDBtop100Folder = "C://Users//Christian//BioInformatics//projects//knots//code//top100H//"

PDBfilename1enh = r"C:/Users/Christian/BioInformatics/projects/knots/code/1ENH/1ENH.pdb"
PDBfilename1art = r"C:/Users/Christian/BioInformatics/projects/knots/code/1art/1art.txt"

PDBfilename1ads = PDBtop100Folder + '1adsH'
PDBfilename1ifc = PDBtop100Folder + '1ifcH'
PDBfilename1knb = PDBtop100Folder + '1knbH'
PDBfilename2mhr = PDBtop100Folder + '2mhrH'
PDBfilename1tta = PDBtop100Folder + '1ttaAH'
PDBfilename2mcm = PDBtop100Folder + '2mcmH'
PDBfilename1ptx = PDBtop100Folder + '1ptxH'
PDBfilename1mct = PDBtop100Folder + '1mctIH'
PDBfilename1etm = PDBtop100Folder + '1etmH'
PDBfilename1edm = PDBtop100Folder + '1edmBH'
PDBfilename1not = PDBtop100Folder + '1notH'
PDBfilename1bio = PDBtop100Folder + 'bio1rpoH'
PDBfilename1ctj = PDBtop100Folder + '1ctjH'
PDBfilename1smd = PDBtop100Folder + '1smdH'

PDBfilenameDef = PDBfilename1art

#targetDef = 'cpu'

#@guvectorize(['void(float64[:,:],float64[:,:])'], '(m,n),(m,n)-> ()',
#           target='cpu')

#fastaggr = numba.jit(double[:,:,:](double[:,:,:,:], double)(aggr)

def I(PDBfilename = PDBfilenameDef, degree = 4,  full_b = 1, print_b =1, aggrTime_b = 0, CaChain = '', polyCaChain = '', inNumpy_b = 0, vect_b = 0):
    '''Numpy array based version of the function I_dict, which is based on dictionaries. 
    Computes for the provided PDB-file, PDBfilenameDef, the measures of the desired 
    degree, which can be 2, 4 (default) or 6. The chain of C-alpha atoms and its 
    accompaning chain of line-segments between each adjacent pair of C-alphas 
    can be provided directly, CaChain resp. polyCaChain, in which case the PDBfilename 
    can be left void ('').
    Input (other):
    CaChain: dictionary mapping each for each sub-chain of the structure (given by its 
    chain identifier) to the list (or numpy array) of 3d-coordinates of the alpha C atoms.
    polyCaChain: as CaChain but where the list consists of the 3d-coordiantes of pairs
    of neighbouring alpha C atoms.
        Obs: we use CaChain and polyCaChain both for these dictionaries (the case in this fct) 
        and, below in other fct's, for the values in these dictionaries; that is, CaChain in the 
        present fct means the dictionary, while polyCaChain in the aggr fct means the chain of 
        3d-coordinates of the alpha C atoms, ie polyCaChain[chain id] in the notation of the present 
        fct.
    inNumpy_b: boolean (1/0) to control whether the input CaChain and 
    polyCaChain are numpy arrays or not.
    full_b: only active with degree = 4; boolean (1/0) controlling whether the 
    degree 4 measures are computed on the full simplex or only the measures needed 
    for computing the degree four measures on the whole chain.
    vect_b: boolean (1/0) controlling if computation should use vectorized
    writhe-term computation (1) or not (0).
    aggrTime_b: return time spend on computation/aggregation (1) or not (0)
    Returns: a dictionary of dictionaries: outer key maps each chain-id of the 
    structure (each chain-id represents a sub-structure of the total structure) 
    to a dictionary mapping the name of each desired meaure (inner key) to its values
    across the simplex (e.g. 'I123456' is mapped to a 2d-array in which at index (i,j) 
    the value of the measure is found for the position (i,j) in the simplex.
    '''
    #for holding final output:
    outDict = {}
    #load from file if chains not supplied:
    if (len(CaChain) < 0.5 and len(polyCaChain) < 0.5):
        CaChain, pChain = getPolygonalChain(PDBfilename = PDBfilename, outNumpy_b = inNumpy_b)
    else:
        CaChain, pChain = CaChain, polyCaChain
#    wDict, I12, outDict  = aggr(pChain, degree, full_b = full_b)
    #call the computation
    aggrTime = {}
    if aggrTime_b != 1:
        for chain in CaChain.keys():
            if vect_b > 0.5:
                outDict[chain] = aggr_vect(CaChain[chain], pChain[chain], degree, full_b, print_b = print_b, aggrTime_b = aggrTime_b)
            else:
                outDict[chain] = aggr(pChain[chain], degree, full_b, print_b = print_b, aggrTime_b = aggrTime_b)
        return outDict, CaChain, pChain
    elif aggrTime_b == 1:
        for chain in CaChain.keys():
            if vect_b > 0.5:
                outDict[chain], aggrTime[chain] = aggr_vect(CaChain[chain], pChain[chain], degree, full_b, print_b = print_b, aggrTime_b = aggrTime_b)
            else:
                outDict[chain], aggrTime[chain] = aggr(pChain[chain], degree, full_b, print_b = print_b, aggrTime_b = aggrTime_b)
        return outDict, CaChain, pChain, aggrTime
        

def I_naked(PDBfilename = PDBfilenameDef, 
            degree = 4,  
            full_b = 1, 
            print_b =1, 
            aggrTime_b = 0, 
            CaChain = '', 
            polyCaChain = '', 
            inNumpy_b = 0,
            wDict = 0,
            I12 = 0, #value of I12 on full chain
            I12Dict = 0,
            I1234Dict = 0,
            I1234Dict_full = 0,
            I1234Dict_full_aid = 0,
            I1234Dict_full2 = 0,
            I1234Dict_full2_aid = 0,
            I1423Dict = 0, 
            I1423Dict_full0 = 0, 
            I1423Dict_full2 = 0,
            I1423Dict_full2_aid = 0,
            I1324Dict = 0,
            I1324Dict_full = 0,
            I1324Dict_full_aid = 0,
            I1324Dict_full2 = 0,
            I1324Dict_full2_aid = 0,
            #(12)*:
            I123456Dict = 0, 
            I123645Dict = 0,
            I123546Dict = 0,
            #(13)*:
            I132456Dict = 0,
            I132546Dict = 0,
            I132645Dict = 0,
            #(14)*
            I142356Dict = 0,
            I142536Dict = 0,
            I142635Dict = 0,
            #(15)*
            I152346Dict = 0,
            I152436Dict = 0,
            I152634Dict = 0,
            #(16)*
            I162345Dict = 0,
            I162435Dict = 0,
            I162534Dict = 0):
    '''As the function I, but takes as input also a set of existing arrays for holding the invariants' 
    values. To be used in a "global memory allocation" setting, i.e. in which these arrays are
    allocated/initialized up front. Runs only in vectorized version.
    '''
    #for holding final output:
    outDict = {}
    #load from file if chains not supplied:
    if (len(CaChain) < 0.5 and len(polyCaChain) < 0.5):
        CaChain, pChain = getPolygonalChain(PDBfilename = PDBfilename, outNumpy_b = inNumpy_b)
    else:
        CaChain, pChain = CaChain, polyCaChain
#    wDict, I12, outDict  = aggr(pChain, degree, full_b = full_b)
    #call the computation
    aggrTime = 0
    if aggrTime_b != 1:
        for chain in CaChain.keys():
            outDict[chain] = aggr_vect_naked(CaChain[chain], 
                                      pChain[chain], 
                                      degree, 
                                      full_b, 
                                      print_b = print_b, 
                                      aggrTime_b = aggrTime_b, 
                                      wDict = wDict,
                                      I12 = I12, #value of I12 on full chain
                                      I12Dict = I12Dict,
                                      I1234Dict = I1234Dict,
                                      I1234Dict_full = I1234Dict_full,
                                      I1234Dict_full_aid = I1234Dict_full_aid,
                                      I1234Dict_full2 = I1234Dict_full2,
                                      I1234Dict_full2_aid = I1234Dict_full2_aid,
                                      I1423Dict = I1423Dict, 
                                      I1423Dict_full0 = I1423Dict_full0, 
                                      I1423Dict_full2 = I1423Dict_full2,
                                      I1423Dict_full2_aid = I1423Dict_full2_aid,
                                      I1324Dict = I1324Dict,
                                      I1324Dict_full = I1324Dict_full,
                                      I1324Dict_full_aid = I1324Dict_full_aid,
                                      I1324Dict_full2 = I1324Dict_full2,
                                      I1324Dict_full2_aid = I1324Dict_full2_aid,
                                      #(12)*:
                                      I123456Dict = I123456Dict, 
                                      I123645Dict = I123645Dict,
                                      I123546Dict = I123546Dict,
                                      #(13)*:
                                      I132456Dict = I132456Dict,
                                      I132546Dict = I132546Dict,
                                      I132645Dict = I132645Dict,
                                      #(14)*
                                      I142356Dict = I142356Dict,
                                      I142536Dict = I142536Dict,
                                      I142635Dict = I142635Dict,
                                      #(15)*
                                      I152346Dict = I152346Dict,
                                      I152436Dict = I152436Dict,
                                      I152634Dict = I152634Dict,
                                      #(16)*
                                      I162345Dict = I162345Dict,
                                      I162435Dict = I162435Dict,
                                      I162534Dict = I162534Dict)
        return outDict, CaChain, pChain
    elif aggrTime_b == 1:
        aggrTime = {}
        for chain in CaChain.keys():
            outDict[chain], aggrTime[chain] = aggr_vect_naked(CaChain[chain], 
                                      pChain[chain], 
                                      degree, 
                                      full_b, 
                                      print_b = print_b, 
                                      aggrTime_b = aggrTime_b, 
                                      wDict = wDict,
                                      I12 = I12, #value of I12 on full chain
                                      I12Dict = I12Dict,
                                      I1234Dict = I1234Dict,
                                      I1234Dict_full = I1234Dict_full,
                                      I1234Dict_full_aid = I1234Dict_full_aid,
                                      I1234Dict_full2 = I1234Dict_full2,
                                      I1234Dict_full2_aid = I1234Dict_full2_aid,
                                      I1423Dict = I1423Dict, 
                                      I1423Dict_full0 = I1423Dict_full0, 
                                      I1423Dict_full2 = I1423Dict_full2,
                                      I1423Dict_full2_aid = I1423Dict_full2_aid,
                                      I1324Dict = I1324Dict,
                                      I1324Dict_full = I1324Dict_full,
                                      I1324Dict_full_aid = I1324Dict_full_aid,
                                      I1324Dict_full2 = I1324Dict_full2,
                                      I1324Dict_full2_aid = I1324Dict_full2_aid,
                                      #(12)*:
                                      I123456Dict = I123456Dict, 
                                      I123645Dict = I123645Dict,
                                      I123546Dict = I123546Dict,
                                      #(13)*:
                                      I132456Dict = I132456Dict,
                                      I132546Dict = I132546Dict,
                                      I132645Dict = I132645Dict,
                                      #(14)*
                                      I142356Dict = I142356Dict,
                                      I142536Dict = I142536Dict,
                                      I142635Dict = I142635Dict,
                                      #(15)*
                                      I152346Dict = I152346Dict,
                                      I152436Dict = I152436Dict,
                                      I152634Dict = I152634Dict,
                                      #(16)*
                                      I162345Dict = I162345Dict,
                                      I162435Dict = I162435Dict,
                                      I162534Dict = I162534Dict)
        return outDict, CaChain, pChain, aggrTime

#@jit(float64[:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:](float64[:,:,:], int64, int64, int64), target = 'cpu')
#@autojit
def aggr(pChain, degree, full_b, print_b, aggrTime_b):
    '''The computation  called in the function I. Input parameters not
    dealt with there: pChain -- a polyCaChain; outDict: both in/output; the 
    dictionary returned in the I function (see there for more).'''
    #Get length of chain etc:
    L = len(pChain)
#    print_b = 0
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    lChain = range(L)
    if print_b > 0.5:
        print lChain
    #Initialize array parameters for storing the results (the 'Dict' suffix is 
    #used for easy translation from dictionary-based version):
    wDict = np.zeros((L+1,L+1))
    I12 = 0 #value of I12 on full chain
    I12Dict = np.zeros((L+1,L+1))
    I1234Dict = np.zeros((L+1,L+1))
    I1234Dict_full = np.zeros((L+1,L+1))
    I1234Dict_full_aid = np.zeros((L+1,L+1))
    I1234Dict_full2 = np.zeros((L+1,L+1))
    I1234Dict_full2_aid = np.zeros((L+1,L+1))
    I1423Dict = np.zeros((L+1,L+1)) 
    I1423Dict_full0 = np.zeros((L+1,L+1)) 
    I1423Dict_full2 = np.zeros((L+1,L+1))
    I1423Dict_full2_aid = np.zeros((L+1,L+1))
    I1324Dict = np.zeros((L+1,L+1))
    I1324Dict_full = np.zeros((L+1,L+1))
    I1324Dict_full_aid = np.zeros((L+1,L+1))
    I1324Dict_full2 = np.zeros((L+1,L+1))
    I1324Dict_full2_aid = np.zeros((L+1,L+1))
    #(12)*:
    I123456Dict = np.zeros((L+1,L+1)) 
    I123645Dict = np.zeros((L+1,L+1))
    I123546Dict = np.zeros((L+1,L+1))
    #(13)*:
    I132456Dict = np.zeros((L+1,L+1))
    I132546Dict = np.zeros((L+1,L+1))
    I132645Dict = np.zeros((L+1,L+1))
    #(14)*
    I142356Dict = np.zeros((L+1,L+1))
    I142536Dict = np.zeros((L+1,L+1))
    I142635Dict = np.zeros((L+1,L+1))
    #(15)*
    I152346Dict = np.zeros((L+1,L+1))
    I152436Dict = np.zeros((L+1,L+1))
    I152634Dict = np.zeros((L+1,L+1))
    #(16)*
    I162345Dict = np.zeros((L+1,L+1))
    I162435Dict = np.zeros((L+1,L+1))
    I162534Dict = np.zeros((L+1,L+1))

    start = time.time()    
    
    #Main part: computing the measures using recursion formulas
    #one block for each choice of degree; only one block gets executed:   
    if degree == 0:
        if print_b > 0.5:
            print "Computes w-values only."
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = w(pChain[i],pChain[j])
                    wDict[i,j] += wVal
    elif degree == 2:
        if print_b > 0.5:
            print "Computes degree 2 measures only."
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = w(pChain[i],pChain[j])
                    wDict[i,j] += wVal
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
    elif degree == 4:
        if print_b > 0.5:
            print "Computes degree 2 and 4 measures only."
        for i in lChain[::-1]:  
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = w(pChain[i],pChain[j])
                    wDict[i,j] += wVal
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
                    #I1234: I1234(i,j;N)
                    I1234Dict[i,j] += I12Dict[j+1,L-1]*wVal + I1234Dict[i,j-1] + I1234Dict[i+1,j] - I1234Dict[i+1,j-1]
                    #I1423: I1423(i,j)
                    I1423Dict[i,j] += I12Dict[i+1,j-1]*wVal +  I1423Dict[i,j-1] + I1423Dict[i+1,j] - I1423Dict[i+1,j-1]
                    #I1324: I1324(i,j;N)
                    I1324Dict[i,j] += (I12Dict[i+1,L-1] - I12Dict[i+1,j] -I12Dict[j,L-1] )*wVal + I1324Dict[i,j-1] + I1324Dict[i+1,j] - I1324Dict[i+1,j-1]
                    if full_b == 1:
                        for b in range(i+1,j-1): #obs: the wDict values needed have been computed
                            wVal1 = wDict[i,b]
                            I1234Dict_full_aid[i,j] +=  I12Dict[b+1,j]*wVal1
                            I1324Dict_full_aid[i,j] += I12Dict[b,j]*wVal1
#                            I1234Dict_full[(i,j)] += (I12Dict[(b+1,j)] - I12Dict[(b+1,j-1)])*wDict[(i,b)]
#                            I1324Dict_full[(i,j)] += (-I12Dict[(b,j)] + I12Dict[(b,j-1)])*wDict[(i,b)]
                        #I1234Dict_full add up terms and recursion: I1234(i,j)    
                        I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                        I1234Dict_full[i,j] +=  I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                        #I1324Dict_full add up terms recursion: I1324(i,j)
                        I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                        I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]                                                   
    elif degree == 6:
        if print_b > 0.5:
            print "Computes degree 2, 4 and 6 measures." 
        for i in lChain[::-1]:  
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = w(pChain[i],pChain[j])
                    wDict[i,j] += wVal
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
                    #degree 4 measures:
                    #I1234: I1234(i,j;N)
                    I1234Dict[i,j] += I12Dict[j+1,L-1]*wVal + I1234Dict[i,j-1] + I1234Dict[i+1,j] - I1234Dict[i+1,j-1]
                    #I1423: I1423(i,j)                    
                    I1423Dict[i,j] += I12Dict[i+1,j-1]*wVal +  I1423Dict[i,j-1] + I1423Dict[i+1,j] - I1423Dict[i+1,j-1]
                    #I1324: I1324(i,j;N)
                    I1324Dict[i,j] += (I12Dict[i+1,L-1] - I12Dict[i+1,j] -I12Dict[j,L-1] )*wVal + I1324Dict[i,j-1] + I1324Dict[i+1,j] - I1324Dict[i+1,j-1]
                    for b in range(i+1,j-1): #obs: the wDict values needed have been computed
                        wVal1 = wDict[i,b]
                        I1234Dict_full_aid[i,j] +=  I12Dict[b+1,j]*wVal1
                        I1324Dict_full_aid[i,j] += I12Dict[b,j]*wVal1
                    #I1234Dict_full add up terms and recursion: I1234(i,j))    
                    I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                    I1234Dict_full[i,j] +=  I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                    #I1324Dict_full add up terms recursion: I1324(i,j)
                    I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                    I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]                            
                    #I1423_full: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet
                    I1423Dict_full0[i+1,j] += (I12Dict[i+2,j-1] - I12Dict[i+2,j])*(I12Dict[i+2,L-1] - I12Dict[i+2,j] - I12Dict[i+1,L-1] + I12Dict[i+1,j])
                    I1423Dict_full0[i+1,j] += I1423Dict_full0[i+2,j] + I1423Dict_full0[i+1,j-1] - I1423Dict_full0[i+2,j-1]
                    #to compute certain degree 6 measures we use two auxillary degree 4 
                    #measures; the recursion demands to sum j in the - direction (ie 
                    #from above and down):
                    k = L-1 + i -j+1
                    for c in range(k+1,L): #Obs: the wDict values needed have already been computed (since k,c >i)
                        wVal2 = wDict[k,c]
                        I1324Dict_full2_aid[i+1,k] += wVal2*(I12Dict[i+1,c-1] - I12Dict[i+1,k] - I12Dict[k,c-1])
                        I1423Dict_full2_aid[i+1,k] += wVal2*(I12Dict[i+1,L-1] - I12Dict[i+1,c] - I12Dict[k,L-1] + I12Dict[k,c]) 
                    #I1324Dict_full2, recursion: I1324(i;j,N)
                    I1324Dict_full2[i+1,k] += I1324Dict_full2_aid[i+1,k] + I1324Dict_full2[i+1,k+1]
                    #I1423Dict_full2, recursion: I1423(i;j,N)
                    I1423Dict_full2[i+1,k] += I1423Dict_full2_aid[i+1,k] + I1423Dict_full2[i+1,k+1] 
                    #degree 6 measures:
                    #(12)*:
                    I123456Dict[i,j] += I1234Dict[j+1,L-1]*wVal + I123456Dict[i+1,j] + I123456Dict[i,j-1] - I123456Dict[i+1,j-1] 
                    I123645Dict[i,j] += I1423Dict[j+1,L-1]*wVal + I123645Dict[i+1,j] + I123645Dict[i,j-1] - I123645Dict[i+1,j-1] 
                    I123546Dict[i,j] += I1324Dict[j+1,L-1]*wVal+ I123546Dict[i+1,j] + I123546Dict[i,j-1] - I123546Dict[i+1,j-1] 
                    #(13)*:
                    I132456Dict[i,j] += (I1234Dict[i+1,L-1] - I1234Dict[i+1,j] -I1234Dict[j,L-1])*wVal + I132456Dict[i+1,j] + I132456Dict[i,j-1] - I132456Dict[i+1,j-1]
                    I132546Dict[i+1,j] += (I1324Dict_full2[i+2,j+1] - I1324Dict_full[j,L-1])*wDict[i+1,j] + I132546Dict[i+2,j] + I132546Dict[i+1,j-1] - I132546Dict[i+2,j-1]
                    I132645Dict[i+1,j] += (I1423Dict_full2[i+2,j+1] - I1423Dict[j,L-1])*wDict[i+1,j] + I132645Dict[i+2,j] + I132645Dict[i+1,j-1] - I132645Dict[i+2,j-1]
                    #(14)*
                    #allowing direct computation:
                    I142356Dict[i,j] += I12Dict[i+1,j-1]*I12Dict[j+1,L-1]*wVal + I142356Dict[i+1,j] + I142356Dict[i,j-1] - I142356Dict[i+1,j-1]
                    #back to the more complicated parts, now I142536:
                    #three parts recognized in decompsition as already known:
                    I142536Dict[i+1,j] += wDict[i+1,j]*(I1324Dict_full[i+2,L-1] - I1324Dict[i+2,j] - I1324Dict_full2[i+2,j])
                    #recursion part:
                    I142536Dict[i+1,j] += I142536Dict[i+2,j] + I142536Dict[i+1,j-1] - I142536Dict[i+2,j-1] 
                    #Next complicated part, I142635:
                    #three parts recognized in decompsition as already known:
                    I142635Dict[i+1,j] += wDict[i+1,j]*(I1423Dict[i+2,L-1] - I1423Dict_full2[i+2,j] - I1423Dict_full0[i+2,j])
                    #recursion part:
                    I142635Dict[i+1,j] += I142635Dict[i+2,j] + I142635Dict[i+1,j-1] - I142635Dict[i+2,j-1]
                    #(15)*
                    I152346Dict[i,j] += wVal*(I1234Dict[i+1,j-1] - I1234Dict_full[i+1,j] - I12Dict[j,L-1]*I12Dict[i+1,j-1])
                    I152346Dict[i,j] += I152346Dict[i+1,j] + I152346Dict[i,j-1] - I152346Dict[i+1,j-1]
                    I152436Dict[i,j] += wVal*(I1324Dict[i+1,j-1] - I1324Dict_full[i+1,j])
                    I152436Dict[i,j] += I152436Dict[i+1,j] + I152436Dict[i,j-1] - I152436Dict[i+1,j-1]
                    I152634Dict[i,j] += wVal*(I1423Dict_full0[i+1,j-1] - I1423Dict[i+1,j])
                    I152634Dict[i,j] += I152634Dict[i+1,j] + I152634Dict[i,j-1] - I152634Dict[i+1,j-1]
                    #(16)*
                    I162345Dict[i,j] += wVal*I1234Dict_full[i+1,j-1] #ER FULL OK?
                    I162345Dict[i,j] += I162345Dict[i+1,j] + I162345Dict[i,j-1] - I162345Dict[i+1,j-1]
                    I162435Dict[i,j] += wVal*I1324Dict_full[i+1,j-1] #ER FULL OK?
                    I162435Dict[i,j] += I162435Dict[i+1,j] + I162435Dict[i,j-1] - I162435Dict[i+1,j-1]
                    I162534Dict[i,j] += wVal*I1423Dict[i+1,j-1]
                    I162534Dict[i,j] += I162534Dict[i+1,j] + I162534Dict[i,j-1] - I162534Dict[i+1,j-1]                    
        #the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:
        for j in lChain[::-1]:
            i = 0
            for c in range(j+1,L):
                wVal = wDict[j,c]
#                        I1234Dict_full2_aid[(i,j)] +=  I12Dict[(i,c-1)]*wDict[(j,c)]
                I1324Dict_full2_aid[i,j] += wVal*(I12Dict[i,c-1] - I12Dict[i,j] - I12Dict[j,c-1])
                I1423Dict_full2_aid[i,j] += wVal*(I12Dict[i,L-1] - I12Dict[i,c] - I12Dict[j,L-1] + I12Dict[j,c]) 
           #recursion:
            I1324Dict_full2[i,j] += I1324Dict_full2_aid[i,j] + I1324Dict_full2[i,j+1]
            I1423Dict_full2[i,j] += I1423Dict_full2_aid[i,j] + I1423Dict_full2[i,j+1] 
        #Missing terms for some of the degree 6 measures and I1423_full0:
        for j in lChain[1:]:
            i=0
            wVal = wDict[i,j]
            I1423Dict_full0[i,j] += (I12Dict[i+1,j-1] - I12Dict[i+1,j])*(I12Dict[i+1,L-1] - I12Dict[i+1,j] - I12Dict[i,L-1] + I12Dict[i,j])
            I1423Dict_full0[i,j] += I1423Dict_full0[i+1,j] + I1423Dict_full0[i,j-1] - I1423Dict_full0[i+1,j-1]
            #
            I132546Dict[i,j] += (I1324Dict_full2[i+1,j+1] - I1324Dict_full[j,L-1])*wVal + I132546Dict[i+1,j] + I132546Dict[i,j-1] - I132546Dict[i+1,j-1]
            I132645Dict[i,j] += (I1423Dict_full2[i+1,j+1] - I1423Dict[j,L-1])*wVal + I132645Dict[i+1,j] + I132645Dict[i,j-1] - I132645Dict[i+1,j-1]
            #
            I142536Dict[i,j] += wVal*(I1324Dict_full[i+1,L-1] - I1324Dict[i+1,j] - I1324Dict_full2[i+1,j])
            I142536Dict[i,j] += I142536Dict[i+1,j] + I142536Dict[i,j-1] - I142536Dict[i+1,j-1] 
            #
            I142635Dict[i,j] += wVal*(I1423Dict[i+1,L-1] - I1423Dict_full2[i+1,j] - I1423Dict_full0[i+1,j])
            I142635Dict[i,j] += I142635Dict[i+1,j] + I142635Dict[i,j-1] - I142635Dict[i+1,j-1]    
    else: 
        print "Degree must be 2, 4 or 6."
#    outDict = np.zeros((L+1,L+1))
#    outDict = np.array((wDict, I12Dict, I1234Dict, I1324Dict, I1423Dict)) #, I1234Dict, I1324Dict, I1423Dict)
    outDict = np.array((wDict, I12Dict, 
                        I1234Dict, I1234Dict_full, 
                        I1324Dict, I1324Dict_full, I1324Dict_full2,
                        I1423Dict, I1423Dict_full0,I1423Dict_full2,
                        I123456Dict, I123546Dict, I123645Dict,
                        I132456Dict, I132546Dict, I132645Dict,
                        I142356Dict, I142536Dict, I142635Dict,
                        I152346Dict, I152436Dict, I152634Dict,
                        I162345Dict, I162435Dict, I162534Dict)) #, I1234Dict, I1324Dict, I1423Dict)
#    return wDict #, I12, outDict

    end = time.time()

    if aggrTime_b != 1:
        return outDict
    elif aggrTime_b ==1:
        return outDict, end - start
        


def I_dict(PDBfilename = PDBfilenameDef, degree = 4,  full_b = 1, print_b =1, CaChain = '', polyCaChain = '', inNumpy_b = 0, vect_b = 0):
    '''For a version of this function based on numpy arrays see I; this version is
    based on dictionaries, i.e. all indexing is done using keys (e.g. the simplex is
    represented by a key for each vertex (i,j)). 
    I_dict computes for the provided PDB-file, PDBfilenameDef, the measures of the 
    desired degree, which can be 2, 4 (default) or 6. The chain of C-alpha atoms and its 
    accompaning chain of line-segments between each adjacent pair of C-alphas 
    can be provided directly, CaChain resp. polyCaChain, in which case the PDBfilename 
    can be left void (''). 
    Input (other):
    CaChain: dictionary mapping each for each sub-chain of the structure (given by its 
    chain identifier) to the list (or numpy array) of 3d-coordinates of the alpha C atoms.
    polyCaChain: as CaChain but where the list consists of the 3d-coordiantes of pairs
    of neighbouring alpha C atoms. 
        Obs: we use CaChain and polyCaChain both for these dictionaries (the case in this fct) 
        and, below in other fct's, for the values in these dictionaries; that is, CaChain in the 
        present fct means the dictionary, while polyCaChain in the aggr_dict fct means the chain of 
        3d-coordinates of the alpha C atoms, ie polyCaChain[chain id] in the notation of the present 
        fct.
    inNumpy_b: boolean (1/0) to control whether the input CaChain and 
    polyCaChain are numpy arrays or not.
    full_b: only active with degree = 4; boolean (1/0) controlling whether the 
    degree 4 measures are computed on the full simplex or only the measures needed 
    for computing the degree four measures on the whole chain.
    vect_b: boolean (1/0) controlling if computation should use vectorized
    writhe-term computation (1) or not (0).
    Returns: a dictionary mapping the name of each desired meaure to its values
    across the simplex (e.g. 'I123456' is mapped to a 2d-array in which at index (i,j) 
    the value of the measure is found for the position (i,j) in the simplex.'''
    #for holding final output:
    outDict = {}
    #load from file if chains not supplied:
    if (len(CaChain) < 0.5 and len(polyCaChain) < 0.5):
        CaChain, pChain = getPolygonalChain(PDBfilename = PDBfilename, outNumpy_b = inNumpy_b)
    else:
        CaChain, pChain = CaChain, polyCaChain
    #call computation:
    for chain in CaChain.keys():    
        if vect_b > 0.5:
            outDict[chain] = aggr_dict_vect(CaChain[chain], pChain[chain], degree=degree, full_b=full_b, print_b = print_b)
        else:
            outDict[chain] = aggr_dict(pChain[chain], degree=degree, full_b=full_b, print_b = print_b)
    return outDict, CaChain, pChain #wDict, I12, outDict, CaChain, pChain


def aggr_dict(pChain, degree, full_b, print_b):
    '''The computation  called in the function I_dict. Input parameters not
    dealt with there: pChain -- a polyCaChain; outDict: both in/output; the 
    dictionary returned in the I_dict function (see there for more).'''
    #Get length of chain etc:
    L = len(pChain)
#    print_b = 0
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    lChain = range(L)
    if print_b > 0.5:
        print lChain
    #init. dictionaries for storing results:
    wDict = {}
    I12 = 0
    I12Dict = {}
    I1234Dict = {}
    I1234Dict_full = {}
    I1234Dict_full_aid = {}
    I1234Dict_full2 = {}
    I1234Dict_full2_aid = {}
    I1423Dict = {} 
    I1423Dict_full0 = {} 
    I1423Dict_full2 = {}
    I1423Dict_full2_aid = {}
    I1324Dict = {}
    I1324Dict_full = {}
    I1324Dict_full_aid = {}
    I1324Dict_full2 = {}
    I1324Dict_full2_aid = {}
    #(12)*:
    I123456Dict = {} 
    I123645Dict = {}
    I123546Dict = {}
    #(13)*:
    I132456Dict = {}
    I132546Dict = {}
    I132645Dict = {}
    #(14)*
    I142356Dict = {}
    I142536Dict = {}
    I142635Dict = {}
    #(15)*
    I152346Dict = {}
    I152436Dict = {}
    I152634Dict = {}
    #(16)*
    I162345Dict = {}
    I162435Dict = {}
    I162534Dict = {}
    #init:
    for i in range(L+1):
        for j in range(L+1):
            wDict[(i,j)] = 0
            I12Dict[(i,j)] = 0
#            I12Dict[(-1,j)] = 0
            if degree > 2: 
                I1234Dict[(i,j)] = 0
                I1234Dict_full[(i,j)] = 0
                I1324Dict[(i,j)] = 0  
                I1423Dict[(i,j)] = 0
#                if full_b ==1:
                I1234Dict_full[(i,j)] = 0
                I1234Dict_full_aid[(i,j)] = 0
                I1234Dict_full2[(i,j)] = 0
                I1234Dict_full2_aid[(i,j)] = 0
                I1324Dict_full[(i,j)] = 0
                I1324Dict_full_aid[(i,j)] = 0
                I1324Dict_full2[(i,j)] = 0
#                I1324Dict_full2[(-1,j)] = 0
                I1324Dict_full2_aid[(i,j)] = 0
#                I1324Dict_full2_aid[(-1,j)] = 0
                I1423Dict_full0[(i,j)] = 0
                I1423Dict_full2[(i,j)] = 0
                I1423Dict_full2_aid[(i,j)] = 0
            if degree > 4:
                I123456Dict[(i,j)] = 0 
                I123645Dict[(i,j)] = 0
                I123546Dict[(i,j)] = 0
                #(13)*:
                I132456Dict[(i,j)] = 0
                I132546Dict[(i,j)] = 0
                I132645Dict[(i,j)] = 0
                #(14)*
                I142356Dict[(i,j)] = 0
                I142536Dict[(i,j)] = 0
                I142635Dict[(i,j)] = 0
                #(15)*
                I152346Dict[(i,j)] = 0
                I152436Dict[(i,j)] = 0
                I152634Dict[(i,j)] = 0
                #(16)*
                I162345Dict[(i,j)] = 0
                I162435Dict[(i,j)] = 0
                I162534Dict[(i,j)] = 0
    #Main part: computing the measures using recursion formulas
    #one block for each choice of degree; only one block gets executed:   
    if degree == 0:
        if print_b > 0.5:
            print "Computes degree 2 measures only."
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = w(pChain[i],pChain[j])
                    wDict[(i,j)] += wVal
    elif degree == 2:
        if print_b > 0.5:
            print "Computes degree 2 measures only."
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = w(pChain[i],pChain[j])
                    wDict[(i,j)] += wVal
                    I12 += wVal
                    I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
    elif degree == 4:
        if print_b > 0.5:
            print "Computes degree 2 and 4 measures only."
        for i in lChain[::-1]:  
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = w(pChain[i],pChain[j])
                    wDict[(i,j)] += wVal
                    I12 += wVal
                    I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
                    #I1234: I1234(i,j;N)
                    I1234Dict[(i,j)] += I12Dict[(j+1,L-1)]*wVal + I1234Dict[(i,j-1)] + I1234Dict[(i+1,j)] - I1234Dict[(i+1,j-1)]
                    #I1423: I1423(i,j)
                    I1423Dict[(i,j)] += I12Dict[(i+1,j-1)]*wVal +  I1423Dict[(i,j-1)] + I1423Dict[(i+1,j)] - I1423Dict[(i+1,j-1)]
                    #I1324: I1324(i,j;N)
                    I1324Dict[(i,j)] += (I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] -I12Dict[(j,L-1)] )*wVal + I1324Dict[(i,j-1)] + I1324Dict[(i+1,j)] - I1324Dict[(i+1,j-1)]
                    if full_b == 1: 
                        for b in range(i+1,j-1):#Obs: the wDict values needed have been computed:
                            wVal1 = wDict[(i,b)]
                            I1234Dict_full_aid[(i,j)] +=  I12Dict[(b+1,j)]*wVal1
                            I1324Dict_full_aid[(i,j)] += I12Dict[(b,j)]*wVal1
                        #I1234Dict_full add up terms and recursion: I1234(i,j)    
                        I1234Dict_full[(i,j)] += I1234Dict_full_aid[(i,j)] -I1234Dict_full_aid[(i,j-1)]
                        I1234Dict_full[(i,j)] +=  I1234Dict_full[(i,j-1)] + I1234Dict_full[(i+1,j)] - I1234Dict_full[(i+1,j-1)]
                        #I1324Dict_full add up terms recursion: I1324(i,j)
                        I1324Dict_full[(i,j)] += I1324Dict_full_aid[(i,j-1)] - I1324Dict_full_aid[(i,j)]   
                        I1324Dict_full[(i,j)] += (I12Dict[(i,j-1)] - I12Dict[(i+1,j-1)])*(I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]) + I1324Dict_full[(i,j-1)] + I1324Dict_full[(i+1,j)] - I1324Dict_full[(i+1,j-1)]                                                   
    elif degree == 6:
        if print_b > 0.5:
            print "Computes degree 2, 4 and 6 measures." 
        for i in lChain[::-1]:  
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = w(pChain[i],pChain[j])
                    wDict[(i,j)] += wVal
                    I12 += wVal
                    I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
                    #degree 4 measures:
                    #I1234: I1234(i,j;N)
                    I1234Dict[(i,j)] += I12Dict[(j+1,L-1)]*wVal + I1234Dict[(i,j-1)] + I1234Dict[(i+1,j)] - I1234Dict[(i+1,j-1)]
                    #I1423: I1423(i,j)                    
                    I1423Dict[(i,j)] += I12Dict[(i+1,j-1)]*wVal +  I1423Dict[(i,j-1)] + I1423Dict[(i+1,j)] - I1423Dict[(i+1,j-1)]
                    #I1324: I1324(i,j;N)
                    I1324Dict[(i,j)] += (I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] -I12Dict[(j,L-1)] )*wVal + I1324Dict[(i,j-1)] + I1324Dict[(i+1,j)] - I1324Dict[(i+1,j-1)]
                    for b in range(i+1,j-1): #Obs: the wDict values needed have been computed:
                            wVal1 = wDict[(i,b)]
                            I1234Dict_full_aid[(i,j)] +=  I12Dict[(b+1,j)]*wVal1
                            I1324Dict_full_aid[(i,j)] += I12Dict[(b,j)]*wVal1
                    #I1234Dict_full add up terms and recursion: I1234(i,j) 
                    I1234Dict_full[(i,j)] += I1234Dict_full_aid[(i,j)] -I1234Dict_full_aid[(i,j-1)]
                    I1234Dict_full[(i,j)] +=  I1234Dict_full[(i,j-1)] + I1234Dict_full[(i+1,j)] - I1234Dict_full[(i+1,j-1)]
                    #I1324Dict_full add up terms recursion: I1324(i,j)
                    I1324Dict_full[(i,j)] += I1324Dict_full_aid[(i,j-1)] - I1324Dict_full_aid[(i,j)]   
                    I1324Dict_full[(i,j)] += (I12Dict[(i,j-1)] - I12Dict[(i+1,j-1)])*(I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]) + I1324Dict_full[(i,j-1)] + I1324Dict_full[(i+1,j)] - I1324Dict_full[(i+1,j-1)]                            
                    #I1423_full: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet
                    I1423Dict_full0[(i+1,j)] += (I12Dict[(i+2,j-1)] - I12Dict[(i+2,j)])*(I12Dict[(i+2,L-1)] - I12Dict[(i+2,j)] - I12Dict[(i+1,L-1)] + I12Dict[(i+1,j)])
                    I1423Dict_full0[(i+1,j)] += I1423Dict_full0[(i+2,j)] + I1423Dict_full0[(i+1,j-1)] - I1423Dict_full0[(i+2,j-1)]
                    #to compute certain degree 6 measures we use two auxillary degree 4 measures; the recursion
                    #demands to sum j in the - direction (ie from above and down):
                    k = L-1 + i -j+1
                    for c in range(k+1,L): #Obs: the wDict values needed have already been computed (since k,c >i)
                        wVal2 = wDict[(k,c)]
                        I1324Dict_full2_aid[(i+1,k)] += wVal2*(I12Dict[(i+1,c-1)] - I12Dict[(i+1,k)] - I12Dict[(k,c-1)])
                        I1423Dict_full2_aid[(i+1,k)] += wVal2*(I12Dict[(i+1,L-1)] - I12Dict[(i+1,c)] - I12Dict[(k,L-1)] + I12Dict[(k,c)]) 
                    #I1324Dict_full2, recursion: I1324(i;j,N)
                    I1324Dict_full2[(i+1,k)] += I1324Dict_full2_aid[(i+1,k)] + I1324Dict_full2[(i+1,k+1)]
                    #I1423Dict_full2, recursion: I1423(i;j,N)
                    I1423Dict_full2[(i+1,k)] += I1423Dict_full2_aid[(i+1,k)] + I1423Dict_full2[(i+1,k+1)] 
                    #degree 6 measures:
                    #(12)*:
                    I123456Dict[(i,j)] += I1234Dict[(j+1,L-1)]*wVal + I123456Dict[(i+1,j)] + I123456Dict[(i,j-1)] - I123456Dict[(i+1,j-1)] 
                    I123645Dict[(i,j )] += I1423Dict[(j+1,L-1)]*wVal + I123645Dict[(i+1,j)] + I123645Dict[(i,j-1)] - I123645Dict[(i+1,j-1)] 
                    I123546Dict[(i,j)] += I1324Dict[(j+1,L-1)]*wVal+ I123546Dict[(i+1,j)] + I123546Dict[(i,j-1)] - I123546Dict[(i+1,j-1)] 
                    #(13)*:
                    I132456Dict[(i,j)] += (I1234Dict[(i+1,L-1)] - I1234Dict[(i+1,j)] -I1234Dict[(j,L-1)])*wVal + I132456Dict[(i+1,j)] + I132456Dict[(i,j-1)] - I132456Dict[(i+1,j-1)]
                    I132546Dict[(i+1,j)] += (I1324Dict_full2[(i+2,j+1)] - I1324Dict_full[(j,L-1)])*wDict[(i+1,j)] + I132546Dict[(i+2,j)] + I132546Dict[(i+1,j-1)] - I132546Dict[(i+2,j-1)]
                    I132645Dict[(i+1,j)] += (I1423Dict_full2[(i+2,j+1)] - I1423Dict[(j,L-1)])*wDict[(i+1,j)] + I132645Dict[(i+2,j)] + I132645Dict[(i+1,j-1)] - I132645Dict[(i+2,j-1)]
                    #(14)*
                    #allowing direct computation:
                    I142356Dict[(i,j)] += I12Dict[(i+1,j-1)]*I12Dict[(j+1,L-1)]*wVal + I142356Dict[(i+1,j)] + I142356Dict[(i,j-1)] - I142356Dict[(i+1,j-1)]
                    #back to the more complicated parts, now I142536:
                    #three parts recognized in decompsition as already known:
                    I142536Dict[(i+1,j)] += wDict[(i+1,j)]*(I1324Dict_full[(i+2,L-1)] - I1324Dict[(i+2,j)] - I1324Dict_full2[(i+2,j)])
                    #recursion part:
                    I142536Dict[(i+1,j)] += I142536Dict[(i+2,j)] + I142536Dict[(i+1,j-1)] - I142536Dict[(i+2,j-1)] 
                    #Next complicated part, I142635:
                    #three parts recognized in decompsition as already known:
                    I142635Dict[(i+1,j)] += wDict[(i+1,j)]*(I1423Dict[(i+2,L-1)] - I1423Dict_full2[(i+2,j)] - I1423Dict_full0[(i+2,j)])
                    #recursion part:
                    I142635Dict[(i+1,j)] += I142635Dict[(i+2,j)] + I142635Dict[(i+1,j-1)] - I142635Dict[(i+2,j-1)]
                    #(15)*
                    I152346Dict[(i,j)] += wVal*(I1234Dict[(i+1,j-1)] - I1234Dict_full[(i+1,j)] - I12Dict[(j,L-1)]*I12Dict[(i+1,j-1)])
                    I152346Dict[(i,j)] += I152346Dict[(i+1,j)] + I152346Dict[(i,j-1)] - I152346Dict[(i+1,j-1)]
                    I152436Dict[(i,j)] += wVal*(I1324Dict[(i+1,j-1)] - I1324Dict_full[(i+1,j)])
                    I152436Dict[(i,j)] += I152436Dict[(i+1,j)] + I152436Dict[(i,j-1)] - I152436Dict[(i+1,j-1)]
                    I152634Dict[(i,j)] += wVal*(I1423Dict_full0[(i+1,j-1)] - I1423Dict[(i+1,j)])
                    I152634Dict[(i,j)] += I152634Dict[(i+1,j)] + I152634Dict[(i,j-1)] - I152634Dict[(i+1,j-1)]
                    #(16)*
                    I162345Dict[(i,j)] += wVal*I1234Dict_full[(i+1,j-1)] #ER FULL OK?
                    I162345Dict[(i,j)] += I162345Dict[(i+1,j)] + I162345Dict[(i,j-1)] - I162345Dict[(i+1,j-1)]
                    I162435Dict[(i,j)] += wVal*I1324Dict_full[(i+1,j-1)] #ER FULL OK?
                    I162435Dict[(i,j)] += I162435Dict[(i+1,j)] + I162435Dict[(i,j-1)] - I162435Dict[(i+1,j-1)]
                    I162534Dict[(i,j)] += wVal*I1423Dict[(i+1,j-1)]
                    I162534Dict[(i,j)] += I162534Dict[(i+1,j)] + I162534Dict[(i,j-1)] - I162534Dict[(i+1,j-1)]                    
        #the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:
        for j in lChain[::-1]:
            i = 0
            for c in range(j+1,L):
                wVal = wDict[(j,c)]
#                        I1234Dict_full2_aid[(i,j)] +=  I12Dict[(i,c-1)]*wDict[(j,c)]
                I1324Dict_full2_aid[(i,j)] += wVal*(I12Dict[(i,c-1)] - I12Dict[(i,j)] - I12Dict[(j,c-1)])
                I1423Dict_full2_aid[(i,j)] += wVal*(I12Dict[(i,L-1)] - I12Dict[(i,c)] - I12Dict[(j,L-1)] + I12Dict[(j,c)]) 
           #recursion:
            I1324Dict_full2[(i,j)] += I1324Dict_full2_aid[(i,j)] + I1324Dict_full2[(i,j+1)]
            I1423Dict_full2[(i,j)] += I1423Dict_full2_aid[(i,j)] + I1423Dict_full2[(i,j+1)] 
        #Missing terms for some of the degree 6 measures and I1423_full0:
        for j in lChain[1:]:
            i=0
            wVal = wDict[(i,j)]
            I1423Dict_full0[(i,j)] += (I12Dict[(i+1,j-1)] - I12Dict[(i+1,j)])*(I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] - I12Dict[(i,L-1)] + I12Dict[(i,j)])
            I1423Dict_full0[(i,j)] += I1423Dict_full0[(i+1,j)] + I1423Dict_full0[(i,j-1)] - I1423Dict_full0[(i+1,j-1)]
            #
            I132546Dict[(i,j)] += (I1324Dict_full2[(i+1,j+1)] - I1324Dict_full[(j,L-1)])*wVal + I132546Dict[(i+1,j)] + I132546Dict[(i,j-1)] - I132546Dict[(i+1,j-1)]
            I132645Dict[(i,j)] += (I1423Dict_full2[(i+1,j+1)] - I1423Dict[(j,L-1)])*wVal + I132645Dict[(i+1,j)] + I132645Dict[(i,j-1)] - I132645Dict[(i+1,j-1)]
            #
            I142536Dict[(i,j)] += wVal*(I1324Dict_full[(i+1,L-1)] - I1324Dict[(i+1,j)] - I1324Dict_full2[(i+1,j)])
            I142536Dict[(i,j)] += I142536Dict[(i+1,j)] + I142536Dict[(i,j-1)] - I142536Dict[(i+1,j-1)] 
            #
            I142635Dict[(i,j)] += wVal*(I1423Dict[(i+1,L-1)] - I1423Dict_full2[(i+1,j)] - I1423Dict_full0[(i+1,j)])
            I142635Dict[(i,j)] += I142635Dict[(i+1,j)] + I142635Dict[(i,j-1)] - I142635Dict[(i+1,j-1)]    
    else: 
        print "Degree must be 2, 4 or 6."
    outDict = {}
    outDict['w'] = wDict
    outDict['I12'] = I12Dict
    outDict['I1234'] = I1234Dict
    outDict['I1234_full'] = I1234Dict_full
    outDict['I1324'] = I1324Dict
    outDict['I1324_full'] = I1324Dict_full
    outDict['I1324_full2'] = I1324Dict_full2    
    outDict['I1423'] = I1423Dict
    outDict['I1423_full0'] = I1423Dict_full0
    outDict['I1423_full2'] = I1423Dict_full2    

    outDict['I123456'] = I123456Dict
    outDict['I123546'] = I123546Dict
    outDict['I123645'] = I123645Dict

    outDict['I132456'] = I132456Dict
    outDict['I132546'] = I132546Dict
    outDict['I132645'] = I132645Dict

    outDict['I142356'] = I142356Dict
    outDict['I142536'] = I142536Dict
    outDict['I142635'] = I142635Dict
    
    outDict['I152346'] = I152346Dict
    outDict['I152436'] = I152436Dict
    outDict['I152634'] = I152634Dict

    outDict['I162345'] = I162345Dict
    outDict['I162435'] = I162435Dict
    outDict['I162534'] = I162534Dict
    return outDict #wDict, I12, outDict

#Utility for checking that the output of the array-based function, I, matches
#that of the dictionary based version, I_dict:
def compareIversions(PDBfilename = PDBfilenameDef, degree = 4,  full_b = 1, pert_b = 0, residueNrRange = range(10,15), tolerance = 1e-6):
    '''Utility for checking that the output of the array-based function, I, matches
    that of the dictionary based version, I_dict. With pert_b =1 the check is made
    for the Iperturbed functions rather than the I-functions.'''
    if pert_b > 0.5: #comparison of pertubed versions 
        #first compute measures using the two versions:
        outDict, CaChain, polyCaChain  = I(PDBfilename = PDBfilename, degree = degree,  full_b = full_b)
        outDict_dict, CaChain_dict, polyCaChain_dict  = I_dict(PDBfilename = PDBfilename, degree = degree,  full_b = full_b)
        #plug in results to Iperturbed functions and get output:
        #First get the w-input:
        #Using array based version:
        wDict_arr  = {}
        for chain in outDict.keys():
            wDict_arr[chain] = outDict[chain][0]  

        outI = Iperturbed(CaChain = CaChain, polyCaChain = polyCaChain, wDict = wDict_arr, residueNrRange = residueNrRange, degree = degree)[0]
        
        #Using the dict based:
        wDict  = {}
        for chain in outDict_dict.keys():
            wDict[chain] = outDict_dict[chain]['w'] 
               
        outIdict = Iperturbed_dict(CaChain = CaChain_dict, polyCaChain =polyCaChain_dict, wDict = wDict, residueNrRange = residueNrRange, degree = degree)[0]       

    else: #un-perturbed
        outI = I(PDBfilename = PDBfilename, degree = degree,  full_b = full_b)[0]
        outIdict = I_dict(PDBfilename = PDBfilename, degree = degree,  full_b = full_b)[0]
    #measure names in order of output from I:
    Inames = ['w', 'I12', 
                        'I1234', 'I1234_full', 
                        'I1324', 'I1324_full', 'I1324_full2',
                        'I1423', 'I1423_full0','I1423_full2',
                        'I123456', 'I123546', 'I123645',
                        'I132456', 'I132546', 'I132645',
                        'I142356', 'I142536', 'I142635',
                        'I152346', 'I152436', 'I152634',
                        'I162345', 'I162435', 'I162534']
    print "Array version output contains %d measures. Names list has length %d" % (len(outI), len(Inames))
    print "Dictionary based output contains %d measures." % len(outIdict.keys())                    
    
    #check first that both outputs cover the same sub-chains of the structure:
    for chain in outI.keys():
        if not(outIdict.has_key(chain)):
            print "Error: dict based version misses this chain: %s" % chain
        else:
            print "All chains in array based version are in dict based too"
    for chain in outIdict.keys():
        if not(outI.has_key(chain)):
            print "Error: array based version misses this chain: %s" % chain    
        else:
            print "All chains in dict based version are in array based too"
    for chain in outI.keys():
        print "Now considering chain with id: %s" % chain
        for k in range(len(Inames)):
            kName = Inames[k] #key to look up
    #        print "I'm looking at measure: %s" % kName
            cnt = 0
            for ij in outIdict[chain][kName]:
                i = ij[0]
                j = ij[1]
                d = outI[chain][k][i,j] - outIdict[chain][kName][ij]
                if abs(d) > tolerance:
                    cnt +=1 
                    print "Warning: difference of %f in ouput at %s in measure %s" % (d, ij, kName)
            print "Number of differences larger than tolerance for measure %s:%d" % (kName, cnt)
        
    

#######################################################################################
### Functions for computations on perturbed chains
#######################################################################################

def Iperturbed(CaChain, 
                polyCaChain, 
                wDict,
                residueNrRange = range(0,5),
                perturbation = [(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3)], 
                degree = 2,
                vect_b = 0,
                full_b = 1,
                print_b = 1):
    '''Array based version of Iperturbed_dict.
    Computes the change in the measures of the given polychain due to 
    the perturbation. Input: 
    CaChain: C-alpha chain (dictionary of these as in the I and I_dict fct's)
    polyCaChain: chain of segment corresponding to CaChain (dictionary of these as in the I and I_dict fct's)
    wDict: array of w-values as output from I (dictionary of these too)
        Obs: we use CaChain, polyCaChain and wDict both for these dictionaries (the case in this fct) and, below in other fct's, for the values in these
        dictionaries; that is, CaChain in the present fct means the dictionary, while perCaChain in the aggr_pert fct 
        means the chain of 3d-coordinates of the alpha C atoms, ie perCaChain[chain id] in the notation of the present 
        fct.
    residueNrRange: indices of the residues in the CaChain (C-alphas) that should be perturbed
    perturbation: the perturbation given in 3d-coordinates. Must be of lenght equal to residueNrRange
    Other: as I_dict.
    Returns: same as I, but for the perturbed chain.'''
    outDict = {}
    #compute the polygonal chain of the perturbed chain:
    perCaChain, perPolyCaChain, segNrRange = perturbPolygonalChain(CaChain = CaChain, 
                                                    polyCaChain = polyCaChain,
                     residueNrRange = residueNrRange, 
                     perturbation = perturbation)
    #call the computation
    for chain in CaChain.keys():
        if vect_b > 0.5:
            outDict[chain] = aggr_pert_vect(wDict[chain], perCaChain[chain], perPolyCaChain[chain], segNrRange[chain], degree, full_b, print_b = print_b)
        else:
            outDict[chain] = aggr_pert(wDict[chain], perCaChain[chain], perPolyCaChain[chain], segNrRange[chain], degree, full_b, print_b = print_b)
    return outDict, perCaChain, perPolyCaChain, segNrRange

#@jit(float64[:,:,:,:,:,:,:,:,:,:](float64[:,:,:], int64, int64))
def aggr_pert(wDict, perCaChain, perPolyCaChain, segNrRange, degree, full_b, print_b):
    '''The computation  called in the function I_perturbed, non-vectorized. 
    Input: passed from I_perturbed (for more see there).'''
    if print_b > 0.5:
        print segNrRange
    #Get length of chain:
    L = len(perPolyCaChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    lChain = range(L)
    #Initialize array parameters for storing the results (the 'Dict' suffix is 
    #used for easy translation from dictionary-based version):
    wPerDict = np.zeros((L+1,L+1))
    I12 = 0 #value of I12 on full chain
    I12Dict = np.zeros((L+1,L+1))
    I1234Dict = np.zeros((L+1,L+1))
    I1234Dict_full = np.zeros((L+1,L+1))
    I1234Dict_full_aid = np.zeros((L+1,L+1))
    I1234Dict_full2 = np.zeros((L+1,L+1))
    I1234Dict_full2_aid = np.zeros((L+1,L+1))
    I1423Dict = np.zeros((L+1,L+1)) 
    I1423Dict_full0 = np.zeros((L+1,L+1)) 
    I1423Dict_full2 = np.zeros((L+1,L+1))
    I1423Dict_full2_aid = np.zeros((L+1,L+1))
    I1324Dict = np.zeros((L+1,L+1))
    I1324Dict_full = np.zeros((L+1,L+1))
    I1324Dict_full_aid = np.zeros((L+1,L+1))
    I1324Dict_full2 = np.zeros((L+1,L+1))
    I1324Dict_full2_aid = np.zeros((L+1,L+1))
    #(12)*:
    I123456Dict = np.zeros((L+1,L+1)) 
    I123645Dict = np.zeros((L+1,L+1))
    I123546Dict = np.zeros((L+1,L+1))
    #(13)*:
    I132456Dict = np.zeros((L+1,L+1))
    I132546Dict = np.zeros((L+1,L+1))
    I132645Dict = np.zeros((L+1,L+1))
    #(14)*
    I142356Dict = np.zeros((L+1,L+1))
    I142536Dict = np.zeros((L+1,L+1))
    I142635Dict = np.zeros((L+1,L+1))
    #(15)*
    I152346Dict = np.zeros((L+1,L+1))
    I152436Dict = np.zeros((L+1,L+1))
    I152634Dict = np.zeros((L+1,L+1))
    #(16)*
    I162345Dict = np.zeros((L+1,L+1))
    I162435Dict = np.zeros((L+1,L+1))
    I162534Dict = np.zeros((L+1,L+1))
    #Main part: computing the measures using recursion formulas
    #one block for each choice of degree; only one block gets executed:    
    if degree == 0:
        if print_b > 0.5:
            print "Computes w-values only."
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    if ((i in segNrRange) or (j in segNrRange)):
    #                print "in change range %d , %d " % (i,j)
                        wVal = w(perPolyCaChain[i],perPolyCaChain[j])
                        wPerDict[i,j] += wVal
                    else:
                        wVal = wDict[i,j]
                        wPerDict[i,j] += wVal
    if degree == 2:
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset
                if (j > i and j <=L-1):
                    if ((i in segNrRange) or (j in segNrRange)):
    #                print "in change range %d , %d " % (i,j)
                        wVal = w(perPolyCaChain[i],perPolyCaChain[j])
                        wPerDict[i,j] += wVal
                    else:
                        wVal = wDict[i,j]
                        wPerDict[i,j] += wVal
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
    elif degree == 4:
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset
                if (j > i and j <=L-1):
                    if ((i in segNrRange) or (j in segNrRange)):
    #                print "in change range %d , %d " % (i,j)
                        wVal = w(perPolyCaChain[i],perPolyCaChain[j])
                        wPerDict[i,j] += wVal
                    else:
                        wVal = wDict[i,j]
                        wPerDict[i,j] += wVal
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
                    #I1234: I1234(i,j;N)
                    I1234Dict[i,j] += I12Dict[j+1,L-1]*wVal + I1234Dict[i,j-1] + I1234Dict[i+1,j] - I1234Dict[i+1,j-1]
                    #I1423: I1423(i,j)
                    I1423Dict[i,j] += I12Dict[i+1,j-1]*wVal +  I1423Dict[i,j-1] + I1423Dict[i+1,j] - I1423Dict[i+1,j-1]
                    #I1324: I1324(i,j;N)
                    I1324Dict[i,j] += (I12Dict[i+1,L-1] - I12Dict[i+1,j] -I12Dict[j,L-1] )*wVal + I1324Dict[i,j-1] + I1324Dict[i+1,j] - I1324Dict[i+1,j-1]
                    if full_b == 1:
                        for b in range(i+1,j-1): #obs: the wPerDict values needed have been computed
                            wVal1 = wPerDict[i,b]
                            I1234Dict_full_aid[i,j] += I12Dict[b+1,j]*wVal1
                            I1324Dict_full_aid[i,j] += I12Dict[b,j]*wVal1
#                            I1234Dict_full[(i,j)] += (I12Dict[(b+1,j)] - I12Dict[(b+1,j-1)])*wDict[(i,b)]
#                            I1324Dict_full[(i,j)] += (-I12Dict[(b,j)] + I12Dict[(b,j-1)])*wDict[(i,b)]
                        #I1234Dict_full add up terms and recursion: I1234(i,j)    
                        I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                        I1234Dict_full[i,j] += I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                        #I1324Dict_full add up terms recursion: I1324(i,j)
                        I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                        I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]                                                   
    elif degree == 6:
        if print_b > 0.5:
            print "Computes degree 2, 4 and 6 measures." 
        for i in lChain[::-1]:  
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    if ((i in segNrRange) or (j in segNrRange)):
    #                print "in change range %d , %d " % (i,j)
                        wVal = w(perPolyCaChain[i],perPolyCaChain[j])
                        wPerDict[i,j] += wVal
                    else:
                        wVal = wDict[i,j]
                        wPerDict[i,j] += wVal
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
                    #degree 4 measures:
                    #I1234: I1234(i,j;N)
                    I1234Dict[i,j] += I12Dict[j+1,L-1]*wVal + I1234Dict[i,j-1] + I1234Dict[i+1,j] - I1234Dict[i+1,j-1]
                    #I1423: I1423(i,j)                    
                    I1423Dict[i,j] += I12Dict[i+1,j-1]*wVal +  I1423Dict[i,j-1] + I1423Dict[i+1,j] - I1423Dict[i+1,j-1]
                    #I1324: I1324(i,j;N)
                    I1324Dict[i,j] += (I12Dict[i+1,L-1] - I12Dict[i+1,j] -I12Dict[j,L-1] )*wVal + I1324Dict[i,j-1] + I1324Dict[i+1,j] - I1324Dict[i+1,j-1]
                    for b in range(i+1,j-1): #obs: the wPerDict values needed have been computed
                        wVal1 = wPerDict[i,b]
                        I1234Dict_full_aid[i,j] += I12Dict[b+1,j]*wVal1
                        I1324Dict_full_aid[i,j] += I12Dict[b,j]*wVal1
                    #I1234Dict_full add up terms and recursion: I1234(i,j))    
                    I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                    I1234Dict_full[i,j] +=  I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                    #I1324Dict_full add up terms recursion: I1324(i,j)
                    I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                    I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]                            
                    #I1423_full: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet
                    I1423Dict_full0[i+1,j] += (I12Dict[i+2,j-1] - I12Dict[i+2,j])*(I12Dict[i+2,L-1] - I12Dict[i+2,j] - I12Dict[i+1,L-1] + I12Dict[i+1,j])
                    I1423Dict_full0[i+1,j] += I1423Dict_full0[i+2,j] + I1423Dict_full0[i+1,j-1] - I1423Dict_full0[i+2,j-1]
                    #to compute certain degree 6 measures we use two auxillary degree 4 
                    #measures; the recursion demands to sum j in the - direction (ie 
                    #from above and down):
                    k = L-1 + i -j+1
                    for c in range(k+1,L): #Obs: the wDict values needed have already been computed (since k,c >i)
                        wVal2 = wPerDict[k,c]
                        I1324Dict_full2_aid[i+1,k] += wVal2*(I12Dict[i+1,c-1] - I12Dict[i+1,k] - I12Dict[k,c-1])
                        I1423Dict_full2_aid[i+1,k] += wVal2*(I12Dict[i+1,L-1] - I12Dict[i+1,c] - I12Dict[k,L-1] + I12Dict[k,c]) 
                    #I1324Dict_full2, recursion: I1324(i;j,N)
                    I1324Dict_full2[i+1,k] += I1324Dict_full2_aid[i+1,k] + I1324Dict_full2[i+1,k+1]
                    #I1423Dict_full2, recursion: I1423(i;j,N)
                    I1423Dict_full2[i+1,k] += I1423Dict_full2_aid[i+1,k] + I1423Dict_full2[i+1,k+1] 
                    #degree 6 measures:
                    #(12)*:
                    I123456Dict[i,j] += I1234Dict[j+1,L-1]*wVal + I123456Dict[i+1,j] + I123456Dict[i,j-1] - I123456Dict[i+1,j-1] 
                    I123645Dict[i,j] += I1423Dict[j+1,L-1]*wVal + I123645Dict[i+1,j] + I123645Dict[i,j-1] - I123645Dict[i+1,j-1] 
                    I123546Dict[i,j] += I1324Dict[j+1,L-1]*wVal+ I123546Dict[i+1,j] + I123546Dict[i,j-1] - I123546Dict[i+1,j-1] 
                    #(13)*:
                    I132456Dict[i,j] += (I1234Dict[i+1,L-1] - I1234Dict[i+1,j] -I1234Dict[j,L-1])*wVal + I132456Dict[i+1,j] + I132456Dict[i,j-1] - I132456Dict[i+1,j-1]
                    I132546Dict[i+1,j] += (I1324Dict_full2[i+2,j+1] - I1324Dict_full[j,L-1])*wPerDict[i+1,j] + I132546Dict[i+2,j] + I132546Dict[i+1,j-1] - I132546Dict[i+2,j-1]
                    I132645Dict[i+1,j] += (I1423Dict_full2[i+2,j+1] - I1423Dict[j,L-1])*wPerDict[i+1,j] + I132645Dict[i+2,j] + I132645Dict[i+1,j-1] - I132645Dict[i+2,j-1]
                    #(14)*
                    #allowing direct computation:
                    I142356Dict[i,j] += I12Dict[i+1,j-1]*I12Dict[j+1,L-1]*wVal + I142356Dict[i+1,j] + I142356Dict[i,j-1] - I142356Dict[i+1,j-1]
                    #back to the more complicated parts, now I142536:
                    #three parts recognized in decompsition as already known:
                    I142536Dict[i+1,j] += wPerDict[i+1,j]*(I1324Dict_full[i+2,L-1] - I1324Dict[i+2,j] - I1324Dict_full2[i+2,j])
                    #recursion part:
                    I142536Dict[i+1,j] += I142536Dict[i+2,j] + I142536Dict[i+1,j-1] - I142536Dict[i+2,j-1] 
                    #Next complicated part, I142635:
                    #three parts recognized in decompsition as already known:
                    I142635Dict[i+1,j] += wPerDict[i+1,j]*(I1423Dict[i+2,L-1] - I1423Dict_full2[i+2,j] - I1423Dict_full0[i+2,j])
                    #recursion part:
                    I142635Dict[i+1,j] += I142635Dict[i+2,j] + I142635Dict[i+1,j-1] - I142635Dict[i+2,j-1]
                    #(15)*
                    I152346Dict[i,j] += wVal*(I1234Dict[i+1,j-1] - I1234Dict_full[i+1,j] - I12Dict[j,L-1]*I12Dict[i+1,j-1])
                    I152346Dict[i,j] += I152346Dict[i+1,j] + I152346Dict[i,j-1] - I152346Dict[i+1,j-1]
                    I152436Dict[i,j] += wVal*(I1324Dict[i+1,j-1] - I1324Dict_full[i+1,j])
                    I152436Dict[i,j] += I152436Dict[i+1,j] + I152436Dict[i,j-1] - I152436Dict[i+1,j-1]
                    I152634Dict[i,j] += wVal*(I1423Dict_full0[i+1,j-1] - I1423Dict[i+1,j])
                    I152634Dict[i,j] += I152634Dict[i+1,j] + I152634Dict[i,j-1] - I152634Dict[i+1,j-1]
                    #(16)*
                    I162345Dict[i,j] += wVal*I1234Dict_full[i+1,j-1] #ER FULL OK?
                    I162345Dict[i,j] += I162345Dict[i+1,j] + I162345Dict[i,j-1] - I162345Dict[i+1,j-1]
                    I162435Dict[i,j] += wVal*I1324Dict_full[i+1,j-1] #ER FULL OK?
                    I162435Dict[i,j] += I162435Dict[i+1,j] + I162435Dict[i,j-1] - I162435Dict[i+1,j-1]
                    I162534Dict[i,j] += wVal*I1423Dict[i+1,j-1]
                    I162534Dict[i,j] += I162534Dict[i+1,j] + I162534Dict[i,j-1] - I162534Dict[i+1,j-1]                    
        #the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:
        for j in lChain[::-1]:
            i = 0
            for c in range(j+1,L): #Obs: the 2d-loop above has sweeped the whole simplex so the wDictPer has been computed already 
                wVal = wPerDict[j,c]
#                        I1234Dict_full2_aid[(i,j)] +=  I12Dict[(i,c-1)]*wDict[(j,c)]
                I1324Dict_full2_aid[i,j] += wVal*(I12Dict[i,c-1] - I12Dict[i,j] - I12Dict[j,c-1])
                I1423Dict_full2_aid[i,j] += wVal*(I12Dict[i,L-1] - I12Dict[i,c] - I12Dict[j,L-1] + I12Dict[j,c]) 
           #recursion:
            I1324Dict_full2[i,j] += I1324Dict_full2_aid[i,j] + I1324Dict_full2[i,j+1]
            I1423Dict_full2[i,j] += I1423Dict_full2_aid[i,j] + I1423Dict_full2[i,j+1] 
        #Missing terms for some of the degree 6 measures and I1423_full0:
        for j in lChain[1:]:
            i=0
            wVal = wPerDict[i,j] #Obs: the 2d-loop above has sweeped the whole simplex so the wDictPer has been computed already 
            I1423Dict_full0[i,j] += (I12Dict[i+1,j-1] - I12Dict[i+1,j])*(I12Dict[i+1,L-1] - I12Dict[i+1,j] - I12Dict[i,L-1] + I12Dict[i,j])
            I1423Dict_full0[i,j] += I1423Dict_full0[i+1,j] + I1423Dict_full0[i,j-1] - I1423Dict_full0[i+1,j-1]
            #
            I132546Dict[i,j] += (I1324Dict_full2[i+1,j+1] - I1324Dict_full[j,L-1])*wVal + I132546Dict[i+1,j] + I132546Dict[i,j-1] - I132546Dict[i+1,j-1]
            I132645Dict[i,j] += (I1423Dict_full2[i+1,j+1] - I1423Dict[j,L-1])*wVal + I132645Dict[i+1,j] + I132645Dict[i,j-1] - I132645Dict[i+1,j-1]
            #
            I142536Dict[i,j] += wVal*(I1324Dict_full[i+1,L-1] - I1324Dict[i+1,j] - I1324Dict_full2[i+1,j])
            I142536Dict[i,j] += I142536Dict[i+1,j] + I142536Dict[i,j-1] - I142536Dict[i+1,j-1] 
            #
            I142635Dict[i,j] += wVal*(I1423Dict[i+1,L-1] - I1423Dict_full2[i+1,j] - I1423Dict_full0[i+1,j])
            I142635Dict[i,j] += I142635Dict[i+1,j] + I142635Dict[i,j-1] - I142635Dict[i+1,j-1]    
    else: 
        print "Degree must be 2, 4 or 6."                         
    outDict = np.array((wPerDict, I12Dict, 
                        I1234Dict, I1234Dict_full, 
                        I1324Dict, I1324Dict_full, I1324Dict_full2,
                        I1423Dict, I1423Dict_full0,I1423Dict_full2,
                        I123456Dict, I123546Dict, I123645Dict,
                        I132456Dict, I132546Dict, I132645Dict,
                        I142356Dict, I142536Dict, I142635Dict,
                        I152346Dict, I152436Dict, I152634Dict,
                        I162345Dict, I162435Dict, I162534Dict)) 
    return outDict #, perCaChain, perPolyCaChain #, wPerDict, I12Dict, I1234Dict, I1423Dict, I1324Dict,I1234Dict_full, I1324Dict_full


def Iperturbed_dict(CaChain, 
                polyCaChain, 
                wDict,
                residueNrRange = range(0,5),
                perturbation = [(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3)], 
                degree = 2,
                full_b = 1,
                vect_b = 0,
                print_b = 1):
    '''Computes the change in the measures of the given polychain due to 
    the perturbation. 
    Input: 
    CaChain: C-alpha chain (dictionary of these; outer key is chain id for sub-chain of structure; one C-alpha chain per chain id, ie chain id maps to C-alpha chain)
    polyCaChain: chain of segment corresponding to CaChain (dictionary of these; outer key is chain id for sub-chain of structure; one poly-C-alpha chain per chain id)
    wDict: dictionary of w-values as output from I_dict (dictionary of these; outer key is chain id for sub-chain of structure; one wDict per chain id)
        Obs: we use CaChain, polyCaChain and wDict both for these dictionaries (the case in this fct) and, below in other fct's, for the values in these
        dictionaries; that is, CaChain in the present fct means the dictionary, while perCaChain in the aggr_pert_dict fct 
        means the chain of 3d-coordinates of the alpha C atoms, ie perCaChain[chain id] in the notation of the present 
        fct.
    residueNrRange: indices of the residues in the CaChain (C-alphas) that should be perturbed
    perturbation: the perturbation given in 3d-coordinates. Must be of lenght equal to residueNrRange
    Other: as I_dict.

    Returns: same as I_dict, but for the perturbed chain.'''
    outDict = {}
    #compute the polygonal chain of the perturbed chain:
    perCaChain, perPolyCaChain, segNrRange = perturbPolygonalChain(CaChain = CaChain, 
                                                    polyCaChain = polyCaChain,
                                                    residueNrRange = residueNrRange, 
                                                    perturbation = perturbation)
                     
    #call computation:
    for chain in perCaChain.keys():     
        if vect_b > 0.5:
            outDict[chain] = aggr_pert_vect_dict(wDict[chain], perCaChain[chain], perPolyCaChain[chain], segNrRange = segNrRange[chain], degree=degree, full_b=full_b, print_b = print_b)
        else:
            outDict[chain] = aggr_pert_dict(wDict[chain], perCaChain[chain], perPolyCaChain[chain], segNrRange = segNrRange[chain], degree=degree, full_b=full_b, print_b = print_b)
    return outDict, perCaChain, perPolyCaChain, segNrRange #wDict, I12, outDict, CaChain, pChain


def aggr_pert_dict(wDict, perCaChain, perPolyCaChain, segNrRange, degree, full_b, print_b):                     
    '''Computation part for I_perturbed, non-vectorized.'''                 
    if print_b > 0.5:
        print segNrRange
    #Get length of chain:
    L = len(perPolyCaChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    lChain = range(L)
    #init. dictionaries for storing results:
    wPerDict = {}
    I12 = 0
    I12Dict = {}
    I1234Dict = {}
    I1234Dict_full = {}
    I1234Dict_full_aid = {}
    I1234Dict_full2 = {}
    I1234Dict_full2_aid = {}
    I1423Dict = {} 
    I1423Dict_full0 = {} 
    I1423Dict_full2 = {}
    I1423Dict_full2_aid = {}
    I1324Dict = {}
    I1324Dict_full = {}
    I1324Dict_full_aid = {}
    I1324Dict_full2 = {}
    I1324Dict_full2_aid = {}
    #(12)*:
    I123456Dict = {} 
    I123645Dict = {}
    I123546Dict = {}
    #(13)*:
    I132456Dict = {}
    I132546Dict = {}
    I132645Dict = {}
    #(14)*
    I142356Dict = {}
    I142536Dict = {}
    I142635Dict = {}
    #(15)*
    I152346Dict = {}
    I152436Dict = {}
    I152634Dict = {}
    #(16)*
    I162345Dict = {}
    I162435Dict = {}
    I162534Dict = {}
    #init:
    for i in range(L+1):
        for j in range(L+1):
            wPerDict[(i,j)] = 0
            I12Dict[(i,j)] = 0
#            I12Dict[(-1,j)] = 0
            if degree > 2: 
                I1234Dict[(i,j)] = 0
                I1234Dict_full[(i,j)] = 0
                I1324Dict[(i,j)] = 0  
                I1423Dict[(i,j)] = 0
#                if full_b ==1:
                I1234Dict_full[(i,j)] = 0
                I1234Dict_full_aid[(i,j)] = 0
                I1234Dict_full2[(i,j)] = 0
                I1234Dict_full2_aid[(i,j)] = 0
                I1324Dict_full[(i,j)] = 0
                I1324Dict_full_aid[(i,j)] = 0
                I1324Dict_full2[(i,j)] = 0
#                I1324Dict_full2[(-1,j)] = 0
                I1324Dict_full2_aid[(i,j)] = 0
#                I1324Dict_full2_aid[(-1,j)] = 0
                I1423Dict_full0[(i,j)] = 0
                I1423Dict_full2[(i,j)] = 0
                I1423Dict_full2_aid[(i,j)] = 0
            if degree > 4:
                I123456Dict[(i,j)] = 0 
                I123645Dict[(i,j)] = 0
                I123546Dict[(i,j)] = 0
                #(13)*:
                I132456Dict[(i,j)] = 0
                I132546Dict[(i,j)] = 0
                I132645Dict[(i,j)] = 0
                #(14)*
                I142356Dict[(i,j)] = 0
                I142536Dict[(i,j)] = 0
                I142635Dict[(i,j)] = 0
                #(15)*
                I152346Dict[(i,j)] = 0
                I152436Dict[(i,j)] = 0
                I152634Dict[(i,j)] = 0
                #(16)*
                I162345Dict[(i,j)] = 0
                I162435Dict[(i,j)] = 0
                I162534Dict[(i,j)] = 0
    #Main part: computing the measures using recursion formulas
    #one block for each choice of degree; only one block gets executed:   
    if degree == 2:
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset
                if (j > i and j <=L-1):
                    if ((i in segNrRange) or (j in segNrRange)):
    #                print "in change range %d , %d " % (i,j)
                        wVal = w(perPolyCaChain[i],perPolyCaChain[j])
                        wPerDict[(i,j)] += wVal
                    else:
                        wVal = wDict[(i,j)]
                        wPerDict[(i,j)] += wVal
                    I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
    elif degree == 4:
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset
                if (j > i and j <=L-1):
                    if ((i in segNrRange) or (j in segNrRange)):
    #                print "in change range %d , %d " % (i,j)
                        wVal = w(perPolyCaChain[i],perPolyCaChain[j])
                        wPerDict[(i,j)] += wVal
                    else:
                        wVal = wDict[(i,j)]
                        wPerDict[(i,j)] += wVal
                    I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
                    I1234Dict[(i,j)] += I12Dict[(j+1,L-1)]*wVal + I1234Dict[(i,j-1)] + I1234Dict[(i+1,j)] - I1234Dict[(i+1,j-1)]
                    I1423Dict[(i,j)] += I12Dict[(i+1,j-1)]*wVal +  I1423Dict[(i,j-1)] + I1423Dict[(i+1,j)] - I1423Dict[(i+1,j-1)]
                    I1324Dict[(i,j)] += (I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] -I12Dict[(j,L-1)] )*wVal + I1324Dict[(i,j-1)] + I1324Dict[(i+1,j)] - I1324Dict[(i+1,j-1)]
                    if full_b == 1:
                        for b in range(i+1,j-1): #obs: the wPerDict values needed have been computed
                            wVal1 = wPerDict[(i,b)]
                            I1234Dict_full_aid[(i,j)] += I12Dict[(b+1,j)]*wVal1
                            I1324Dict_full_aid[(i,j)] += I12Dict[(b,j)]*wVal1
#                            I1234Dict_full[(i,j)] += (I12Dict[(b+1,j)] - I12Dict[(b+1,j-1)])*wDict[(i,b)]
#                            I1324Dict_full[(i,j)] += (-I12Dict[(b,j)] + I12Dict[(b,j-1)])*wDict[(i,b)]
                        #I1234Dict_full add up terms and recursion: I1234(i,j))    
                        I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                        I1234Dict_full[i,j] +=  I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                        #I1324Dict_full add up terms recursion: I1324(i,j)
                        I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                        I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]            
    elif degree == 6:
        if print_b > 0.5:
            print "Computes degree 2, 4 and 6 measures." 
        for i in lChain[::-1]:  
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    if ((i in segNrRange) or (j in segNrRange)):
    #                print "in change range %d , %d " % (i,j)
                        wVal = w(perPolyCaChain[i],perPolyCaChain[j])
                        wPerDict[(i,j)] += wVal
                    else:
                        wVal = wDict[(i,j)]
                        wPerDict[(i,j)] += wVal
                    I12 += wVal
                    I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
                    #degree 4 measures:
                    #I1234: I1234(i,j;N)
                    I1234Dict[(i,j)] += I12Dict[(j+1,L-1)]*wVal + I1234Dict[(i,j-1)] + I1234Dict[(i+1,j)] - I1234Dict[(i+1,j-1)]
                    #I1423: I1423(i,j)                    
                    I1423Dict[(i,j)] += I12Dict[(i+1,j-1)]*wVal +  I1423Dict[(i,j-1)] + I1423Dict[(i+1,j)] - I1423Dict[(i+1,j-1)]
                    #I1324: I1324(i,j;N)
                    I1324Dict[(i,j)] += (I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] -I12Dict[(j,L-1)] )*wVal + I1324Dict[(i,j-1)] + I1324Dict[(i+1,j)] - I1324Dict[(i+1,j-1)]
                    for b in range(i+1,j-1): #obs: the wPerDict values needed have been computed
                        wVal1 = wPerDict[(i,b)]
                        I1234Dict_full_aid[(i,j)] +=  I12Dict[(b+1,j)]*wVal1
                        I1324Dict_full_aid[(i,j)] += I12Dict[(b,j)]*wVal1
                    #I1234Dict_full add up terms and recursion: I1234(i,j) 
                    I1234Dict_full[(i,j)] += I1234Dict_full_aid[(i,j)] -I1234Dict_full_aid[(i,j-1)]
                    I1234Dict_full[(i,j)] +=  I1234Dict_full[(i,j-1)] + I1234Dict_full[(i+1,j)] - I1234Dict_full[(i+1,j-1)]
                    #I1324Dict_full add up terms recursion: I1324(i,j)
                    I1324Dict_full[(i,j)] += I1324Dict_full_aid[(i,j-1)] - I1324Dict_full_aid[(i,j)]   
                    I1324Dict_full[(i,j)] += (I12Dict[(i,j-1)] - I12Dict[(i+1,j-1)])*(I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]) + I1324Dict_full[(i,j-1)] + I1324Dict_full[(i+1,j)] - I1324Dict_full[(i+1,j-1)]                            
                    #I1423_full: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet
                    I1423Dict_full0[(i+1,j)] += (I12Dict[(i+2,j-1)] - I12Dict[(i+2,j)])*(I12Dict[(i+2,L-1)] - I12Dict[(i+2,j)] - I12Dict[(i+1,L-1)] + I12Dict[(i+1,j)])
                    I1423Dict_full0[(i+1,j)] += I1423Dict_full0[(i+2,j)] + I1423Dict_full0[(i+1,j-1)] - I1423Dict_full0[(i+2,j-1)]
                    #to compute certain degree 6 measures we use two auxillary degree 4 measures; the recursion
                    #demands to sum j in the - direction (ie from above and down):
                    k = L-1 + i -j+1
                    for c in range(k+1,L): #obs: the wPerDict values needed have been computed (since k,c >i)
                        wVal2 = wPerDict[(k,c)]
                        I1324Dict_full2_aid[(i+1,k)] += wVal2*(I12Dict[(i+1,c-1)] - I12Dict[(i+1,k)] - I12Dict[(k,c-1)])
                        I1423Dict_full2_aid[(i+1,k)] += wVal2*(I12Dict[(i+1,L-1)] - I12Dict[(i+1,c)] - I12Dict[(k,L-1)] + I12Dict[(k,c)]) 
                    #I1324Dict_full2, recursion: I1324(i;j,N)
                    I1324Dict_full2[(i+1,k)] += I1324Dict_full2_aid[(i+1,k)] + I1324Dict_full2[(i+1,k+1)]
                    #I1423Dict_full2, recursion: I1423(i;j,N)
                    I1423Dict_full2[(i+1,k)] += I1423Dict_full2_aid[(i+1,k)] + I1423Dict_full2[(i+1,k+1)] 
                    #degree 6 measures:
                    #(12)*:
                    I123456Dict[(i,j)] += I1234Dict[(j+1,L-1)]*wVal + I123456Dict[(i+1,j)] + I123456Dict[(i,j-1)] - I123456Dict[(i+1,j-1)] 
                    I123645Dict[(i,j )] += I1423Dict[(j+1,L-1)]*wVal + I123645Dict[(i+1,j)] + I123645Dict[(i,j-1)] - I123645Dict[(i+1,j-1)] 
                    I123546Dict[(i,j)] += I1324Dict[(j+1,L-1)]*wVal+ I123546Dict[(i+1,j)] + I123546Dict[(i,j-1)] - I123546Dict[(i+1,j-1)] 
                    #(13)*:
                    I132456Dict[(i,j)] += (I1234Dict[(i+1,L-1)] - I1234Dict[(i+1,j)] -I1234Dict[(j,L-1)])*wVal + I132456Dict[(i+1,j)] + I132456Dict[(i,j-1)] - I132456Dict[(i+1,j-1)]
                    I132546Dict[(i+1,j)] += (I1324Dict_full2[(i+2,j+1)] - I1324Dict_full[(j,L-1)])*wPerDict[(i+1,j)] + I132546Dict[(i+2,j)] + I132546Dict[(i+1,j-1)] - I132546Dict[(i+2,j-1)]
                    I132645Dict[(i+1,j)] += (I1423Dict_full2[(i+2,j+1)] - I1423Dict[(j,L-1)])*wPerDict[(i+1,j)] + I132645Dict[(i+2,j)] + I132645Dict[(i+1,j-1)] - I132645Dict[(i+2,j-1)]
                    #(14)*
                    #allowing direct computation:
                    I142356Dict[(i,j)] += I12Dict[(i+1,j-1)]*I12Dict[(j+1,L-1)]*wVal + I142356Dict[(i+1,j)] + I142356Dict[(i,j-1)] - I142356Dict[(i+1,j-1)]
                    #back to the more complicated parts, now I142536:
                    #three parts recognized in decompsition as already known:
                    I142536Dict[(i+1,j)] += wPerDict[(i+1,j)]*(I1324Dict_full[(i+2,L-1)] - I1324Dict[(i+2,j)] - I1324Dict_full2[(i+2,j)])
                    #recursion part:
                    I142536Dict[(i+1,j)] += I142536Dict[(i+2,j)] + I142536Dict[(i+1,j-1)] - I142536Dict[(i+2,j-1)] 
                    #Next complicated part, I142635:
                    #three parts recognized in decompsition as already known:
                    I142635Dict[(i+1,j)] += wPerDict[(i+1,j)]*(I1423Dict[(i+2,L-1)] - I1423Dict_full2[(i+2,j)] - I1423Dict_full0[(i+2,j)])
                    #recursion part:
                    I142635Dict[(i+1,j)] += I142635Dict[(i+2,j)] + I142635Dict[(i+1,j-1)] - I142635Dict[(i+2,j-1)]
                    #(15)*
                    I152346Dict[(i,j)] += wVal*(I1234Dict[(i+1,j-1)] - I1234Dict_full[(i+1,j)] - I12Dict[(j,L-1)]*I12Dict[(i+1,j-1)])
                    I152346Dict[(i,j)] += I152346Dict[(i+1,j)] + I152346Dict[(i,j-1)] - I152346Dict[(i+1,j-1)]
                    I152436Dict[(i,j)] += wVal*(I1324Dict[(i+1,j-1)] - I1324Dict_full[(i+1,j)])
                    I152436Dict[(i,j)] += I152436Dict[(i+1,j)] + I152436Dict[(i,j-1)] - I152436Dict[(i+1,j-1)]
                    I152634Dict[(i,j)] += wVal*(I1423Dict_full0[(i+1,j-1)] - I1423Dict[(i+1,j)])
                    I152634Dict[(i,j)] += I152634Dict[(i+1,j)] + I152634Dict[(i,j-1)] - I152634Dict[(i+1,j-1)]
                    #(16)*
                    I162345Dict[(i,j)] += wVal*I1234Dict_full[(i+1,j-1)] #ER FULL OK?
                    I162345Dict[(i,j)] += I162345Dict[(i+1,j)] + I162345Dict[(i,j-1)] - I162345Dict[(i+1,j-1)]
                    I162435Dict[(i,j)] += wVal*I1324Dict_full[(i+1,j-1)] #ER FULL OK?
                    I162435Dict[(i,j)] += I162435Dict[(i+1,j)] + I162435Dict[(i,j-1)] - I162435Dict[(i+1,j-1)]
                    I162534Dict[(i,j)] += wVal*I1423Dict[(i+1,j-1)]
                    I162534Dict[(i,j)] += I162534Dict[(i+1,j)] + I162534Dict[(i,j-1)] - I162534Dict[(i+1,j-1)]                    
        #the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:
        for j in lChain[::-1]:
            i = 0
            for c in range(j+1,L):
                wVal = wPerDict[(j,c)] #Obs: the 2d-loop above has sweeped the whole simplex so the wDictPer has been computed already 
#                        I1234Dict_full2_aid[(i,j)] +=  I12Dict[(i,c-1)]*wDict[(j,c)]
                I1324Dict_full2_aid[(i,j)] += wVal*(I12Dict[(i,c-1)] - I12Dict[(i,j)] - I12Dict[(j,c-1)])
                I1423Dict_full2_aid[(i,j)] += wVal*(I12Dict[(i,L-1)] - I12Dict[(i,c)] - I12Dict[(j,L-1)] + I12Dict[(j,c)]) 
           #recursion:
            I1324Dict_full2[(i,j)] += I1324Dict_full2_aid[(i,j)] + I1324Dict_full2[(i,j+1)]
            I1423Dict_full2[(i,j)] += I1423Dict_full2_aid[(i,j)] + I1423Dict_full2[(i,j+1)] 
        #Missing terms for some of the degree 6 measures and I1423_full0:
        for j in lChain[1:]:
            i=0
            wVal = wPerDict[(i,j)] #Obs: the 2d-loop above has sweeped the whole simplex so the wDictPer has been computed already 
            I1423Dict_full0[(i,j)] += (I12Dict[(i+1,j-1)] - I12Dict[(i+1,j)])*(I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] - I12Dict[(i,L-1)] + I12Dict[(i,j)])
            I1423Dict_full0[(i,j)] += I1423Dict_full0[(i+1,j)] + I1423Dict_full0[(i,j-1)] - I1423Dict_full0[(i+1,j-1)]
            #
            I132546Dict[(i,j)] += (I1324Dict_full2[(i+1,j+1)] - I1324Dict_full[(j,L-1)])*wVal + I132546Dict[(i+1,j)] + I132546Dict[(i,j-1)] - I132546Dict[(i+1,j-1)]
            I132645Dict[(i,j)] += (I1423Dict_full2[(i+1,j+1)] - I1423Dict[(j,L-1)])*wVal + I132645Dict[(i+1,j)] + I132645Dict[(i,j-1)] - I132645Dict[(i+1,j-1)]
            #
            I142536Dict[(i,j)] += wVal*(I1324Dict_full[(i+1,L-1)] - I1324Dict[(i+1,j)] - I1324Dict_full2[(i+1,j)])
            I142536Dict[(i,j)] += I142536Dict[(i+1,j)] + I142536Dict[(i,j-1)] - I142536Dict[(i+1,j-1)] 
            #
            I142635Dict[(i,j)] += wVal*(I1423Dict[(i+1,L-1)] - I1423Dict_full2[(i+1,j)] - I1423Dict_full0[(i+1,j)])
            I142635Dict[(i,j)] += I142635Dict[(i+1,j)] + I142635Dict[(i,j-1)] - I142635Dict[(i+1,j-1)]    
    else: 
        print "Degree must be 2, 4 or 6."
    outDict = {}
    outDict['w'] = wPerDict
    outDict['I12'] = I12Dict
    outDict['I1234'] = I1234Dict
    outDict['I1234_full'] = I1234Dict_full
    outDict['I1324'] = I1324Dict
    outDict['I1324_full'] = I1324Dict_full
    outDict['I1324_full2'] = I1324Dict_full2    
    outDict['I1423'] = I1423Dict
    outDict['I1423_full0'] = I1423Dict_full0
    outDict['I1423_full2'] = I1423Dict_full2    

    outDict['I123456'] = I123456Dict
    outDict['I123546'] = I123546Dict
    outDict['I123645'] = I123645Dict

    outDict['I132456'] = I132456Dict
    outDict['I132546'] = I132546Dict
    outDict['I132645'] = I132645Dict

    outDict['I142356'] = I142356Dict
    outDict['I142536'] = I142536Dict
    outDict['I142635'] = I142635Dict
    
    outDict['I152346'] = I152346Dict
    outDict['I152436'] = I152436Dict
    outDict['I152634'] = I152634Dict

    outDict['I162345'] = I162345Dict
    outDict['I162435'] = I162435Dict
    outDict['I162534'] = I162534Dict
    return outDict #, wPerDict, I12Dict, I1234Dict, I1423Dict, I1324Dict,I1234Dict_full, I1324Dict_full
    
#######################################################################################
### Vectorized versions/gpu apt versions 
#######################################################################################

def wAll(PDBfilename = PDBfilenameDef, print_b =0, CaChain = '', polyCaChain = '', 
         inNumpy_b = 0, perturb_b = 0, segNrRange = '', dict_b = 1, global_b = 0, wDict = 0):
    '''Function that computes all writhe term (w(i,j)'s) in one go. The function 
    (kernel) that computes the w-values, w_Vect2, is prone for vectorization and can also
    be targeted for gpu-computation. The kernel can also be changed for testing
    vs. other kernels. See git_utils for the function(s).
    Inputs: as I-functions plus perturb_b: if 1 will only compute w's in the 
    specified segNrRange; dict_b: boolean (1/0) controlling if output is in 
    dictionary resp. array (for use in dictionary/array versions of 
    functions call wAll).'''
    #load from file if chains not supplied:
    if (len(CaChain) < 0.5 and len(polyCaChain) < 0.5):
        CaChain, pChain = getPolygonalChain(PDBfilename = PDBfilename, outNumpy_b = inNumpy_b)
    else:
        CaChain, pChain = CaChain, polyCaChain
    L = len(pChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    lChain = range(L)
    if print_b > 0.5:
        print lChain
    #init:
    if dict_b > 0.5: 
        wDict = {}
        for i in range(L+1):
            for j in range(L+1):
                wDict[(i,j)] = 0
    elif global_b == 0:
        wDict = np.zeros((L+1,L+1))
    #else: if global_b ==1 wDict must be passed as input to the function, with wDict a suitably sized array.
        
    #first derive the full set of all pairs of segments on which the writhe mus tbe computed:   
#    seg1array = []
#    seg2array = []
    s101array = []
    s102array = []
    s103array = []
    s111array = []
    s112array = []
    s113array = []
    s201array = []
    s202array = []
    s203array = []
    s211array = []
    s212array = []
    s213array = []
    if perturb_b == 0:
        for i in lChain:
            for j in lChain:
                if (j > i and j <=L-1):
                    s101array.append(pChain[i][0][0])
                    s102array.append(pChain[i][0][1])
                    s103array.append(pChain[i][0][2])
                    s111array.append(pChain[i][1][0])
                    s112array.append(pChain[i][1][1])
                    s113array.append(pChain[i][1][2])
                    s201array.append(pChain[j][0][0])
                    s202array.append(pChain[j][0][1])
                    s203array.append(pChain[j][0][2])
                    s211array.append(pChain[j][1][0])
                    s212array.append(pChain[j][1][1])
                    s213array.append(pChain[j][1][2])
    #                seg1array.append(pChain[i])
    #                seg2array.append(pChain[j])
    #    seg1array = np.array(seg1array)
    #    seg2array = np.array(seg2array)
        s101array = np.array(s101array)
        s102array = np.array(s102array)
        s103array = np.array(s103array)
        s111array = np.array(s111array)
        s112array = np.array(s112array)
        s113array = np.array(s113array)
        s201array = np.array(s201array)
        s202array = np.array(s202array)
        s203array = np.array(s203array)
        s211array = np.array(s211array)
        s212array = np.array(s212array)
        s213array = np.array(s213array)
        wAll =0
        wAll = wVect_2(s101array,s102array,s103array,s111array,s112array,s113array,s201array,s202array,s203array, s211array,s212array,s213array)
        #write the seq of results into a dictionary or array as desired
        cnt = 0
        if dict_b > 0.5:
            for i in lChain:
                for j in lChain:
                    if (j > i and j <=L-1):
                        wDict[(i,j)] += wAll[cnt]
                        cnt +=1
        else:
            for i in lChain:
                for j in lChain:
                    if (j > i and j <=L-1):
                        wDict[i,j] += wAll[cnt]
                        cnt +=1
    elif perturb_b == 1:
        lastSeg = segNrRange[::-1][0]
        firstSeg = segNrRange[0]
        for i in range(firstSeg):
            for j in segNrRange:
                s101array.append(pChain[i][0][0])
                s102array.append(pChain[i][0][1])
                s103array.append(pChain[i][0][2])
                s111array.append(pChain[i][1][0])
                s112array.append(pChain[i][1][1])
                s113array.append(pChain[i][1][2])
                s201array.append(pChain[j][0][0])
                s202array.append(pChain[j][0][1])
                s203array.append(pChain[j][0][2])
                s211array.append(pChain[j][1][0])
                s212array.append(pChain[j][1][1])
                s213array.append(pChain[j][1][2])
        for i in segNrRange:
            for j in range(i, L):
                s101array.append(pChain[i][0][0])
                s102array.append(pChain[i][0][1])
                s103array.append(pChain[i][0][2])
                s111array.append(pChain[i][1][0])
                s112array.append(pChain[i][1][1])
                s113array.append(pChain[i][1][2])
                s201array.append(pChain[j][0][0])
                s202array.append(pChain[j][0][1])
                s203array.append(pChain[j][0][2])
                s211array.append(pChain[j][1][0])
                s212array.append(pChain[j][1][1])
                s213array.append(pChain[j][1][2])
        s101array = np.array(s101array)
        s102array = np.array(s102array)
        s103array = np.array(s103array)
        s111array = np.array(s111array)
        s112array = np.array(s112array)
        s113array = np.array(s113array)
        s201array = np.array(s201array)
        s202array = np.array(s202array)
        s203array = np.array(s203array)
        s211array = np.array(s211array)
        s212array = np.array(s212array)
        s213array = np.array(s213array)
        wAll =0
        #compute the write terms:
        wAll = wVect_2(s101array,s102array,s103array,s111array,s112array,s113array,s201array,s202array,s203array, s211array,s212array,s213array)
        cnt = 0
        #write the seq of results into a dictionary or array as desired:
        if dict_b > 0.5:
            for i in lChain:
                for j in lChain:
                    if (j > i and ((i in segNrRange) or (j in segNrRange))):
                        wDict[(i,j)] += wAll[cnt]
                        cnt +=1
                    else:
                        wDict[(i,j)] = 0
        else:
            for i in lChain:
                for j in lChain:
                    if (j > i and ((i in segNrRange) or (j in segNrRange))):
                        wDict[i,j] += wAll[cnt]
                        cnt +=1
                    else:
                        wDict[i,j] = 0
    else:
        print "Param perturb_b must be 0 (non-perturbation use) or 1 (perturbation use)."
    return wDict


def aggr_dict_vect(CaChain, pChain, degree, full_b, print_b):
    '''As aggr_dict, but uses (can use) vectorized computation of the 
    writhe-terms.'''
    #Get length of chain etc:
    L = len(pChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    lChain = range(L)
    if print_b > 0.5:
        print lChain
    #compute the writhe terms in one go:
    wDict = wAll(CaChain = CaChain, polyCaChain = pChain) 
    #init. dictionaries for storing results:
    I12 = 0
    I12Dict = {}
    I1234Dict = {}
    I1234Dict_full = {}
    I1234Dict_full_aid = {}
    I1234Dict_full2 = {}
    I1234Dict_full2_aid = {}
    I1423Dict = {} 
    I1423Dict_full0 = {} 
    I1423Dict_full2 = {}
    I1423Dict_full2_aid = {}
    I1324Dict = {}
    I1324Dict_full = {}
    I1324Dict_full_aid = {}
    I1324Dict_full2 = {}
    I1324Dict_full2_aid = {}
    #(12)*:
    I123456Dict = {} 
    I123645Dict = {}
    I123546Dict = {}
    #(13)*:
    I132456Dict = {}
    I132546Dict = {}
    I132645Dict = {}
    #(14)*
    I142356Dict = {}
    I142536Dict = {}
    I142635Dict = {}
    #(15)*
    I152346Dict = {}
    I152436Dict = {}
    I152634Dict = {}
    #(16)*
    I162345Dict = {}
    I162435Dict = {}
    I162534Dict = {}
    #init:
    for i in range(L+1):
        for j in range(L+1):
            I12Dict[(i,j)] = 0
#            I12Dict[(-1,j)] = 0
            if degree > 2: 
                I1234Dict[(i,j)] = 0
                I1234Dict_full[(i,j)] = 0
                I1324Dict[(i,j)] = 0  
                I1423Dict[(i,j)] = 0
#                if full_b ==1:
                I1234Dict_full[(i,j)] = 0
                I1234Dict_full_aid[(i,j)] = 0
                I1234Dict_full2[(i,j)] = 0
                I1234Dict_full2_aid[(i,j)] = 0
                I1324Dict_full[(i,j)] = 0
                I1324Dict_full_aid[(i,j)] = 0
                I1324Dict_full2[(i,j)] = 0
#                I1324Dict_full2[(-1,j)] = 0
                I1324Dict_full2_aid[(i,j)] = 0
#                I1324Dict_full2_aid[(-1,j)] = 0
                I1423Dict_full0[(i,j)] = 0
                I1423Dict_full2[(i,j)] = 0
                I1423Dict_full2_aid[(i,j)] = 0
            if degree > 4:
                I123456Dict[(i,j)] = 0 
                I123645Dict[(i,j)] = 0
                I123546Dict[(i,j)] = 0
                #(13)*:
                I132456Dict[(i,j)] = 0
                I132546Dict[(i,j)] = 0
                I132645Dict[(i,j)] = 0
                #(14)*
                I142356Dict[(i,j)] = 0
                I142536Dict[(i,j)] = 0
                I142635Dict[(i,j)] = 0
                #(15)*
                I152346Dict[(i,j)] = 0
                I152436Dict[(i,j)] = 0
                I152634Dict[(i,j)] = 0
                #(16)*
                I162345Dict[(i,j)] = 0
                I162435Dict[(i,j)] = 0
                I162534Dict[(i,j)] = 0
    #Main part: computing the measures using recursion formulas
    #one block for each choice of degree; only one block gets executed:   
    if degree == 0:
        if print_b > 0.5:
            print "Computes degree 2 measures only."
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = wDict[(i,j)]
    elif degree == 2:
        if print_b > 0.5:
            print "Computes degree 2 measures only."
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = wDict[(i,j)]
                    I12 += wVal
                    I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
    elif degree == 4:
        if print_b > 0.5:
            print "Computes degree 2 and 4 measures only."
        for i in lChain[::-1]:  
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = wDict[(i,j)]
                    I12 += wVal
                    I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
                    #I1234: I1234(i,j;N)
                    I1234Dict[(i,j)] += I12Dict[(j+1,L-1)]*wVal + I1234Dict[(i,j-1)] + I1234Dict[(i+1,j)] - I1234Dict[(i+1,j-1)]
                    #I1423: I1423(i,j)
                    I1423Dict[(i,j)] += I12Dict[(i+1,j-1)]*wVal +  I1423Dict[(i,j-1)] + I1423Dict[(i+1,j)] - I1423Dict[(i+1,j-1)]
                    #I1324: I1324(i,j;N)
                    I1324Dict[(i,j)] += (I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] -I12Dict[(j,L-1)] )*wVal + I1324Dict[(i,j-1)] + I1324Dict[(i+1,j)] - I1324Dict[(i+1,j-1)]
                    if full_b == 1: 
                        for b in range(i+1,j-1):#Obs: the wDict values needed have been computed:
                            wVal1 = wDict[(i,b)]
                            I1234Dict_full_aid[(i,j)] +=  I12Dict[(b+1,j)]*wVal1
                            I1324Dict_full_aid[(i,j)] += I12Dict[(b,j)]*wVal1
                        #I1234Dict_full add up terms and recursion: I1234(i,j)    
                        I1234Dict_full[(i,j)] += I1234Dict_full_aid[(i,j)] -I1234Dict_full_aid[(i,j-1)]
                        I1234Dict_full[(i,j)] +=  I1234Dict_full[(i,j-1)] + I1234Dict_full[(i+1,j)] - I1234Dict_full[(i+1,j-1)]
                        #I1324Dict_full add up terms recursion: I1324(i,j)
                        I1324Dict_full[(i,j)] += I1324Dict_full_aid[(i,j-1)] - I1324Dict_full_aid[(i,j)]   
                        I1324Dict_full[(i,j)] += (I12Dict[(i,j-1)] - I12Dict[(i+1,j-1)])*(I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]) + I1324Dict_full[(i,j-1)] + I1324Dict_full[(i+1,j)] - I1324Dict_full[(i+1,j-1)]                                                   
    elif degree == 6:
        if print_b > 0.5:
            print "Computes degree 2, 4 and 6 measures." 
        for i in lChain[::-1]:  
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = wDict[(i,j)]
                    I12 += wVal
                    I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
                    #degree 4 measures:
                    #I1234: I1234(i,j;N)
                    I1234Dict[(i,j)] += I12Dict[(j+1,L-1)]*wVal + I1234Dict[(i,j-1)] + I1234Dict[(i+1,j)] - I1234Dict[(i+1,j-1)]
                    #I1423: I1423(i,j)                    
                    I1423Dict[(i,j)] += I12Dict[(i+1,j-1)]*wVal +  I1423Dict[(i,j-1)] + I1423Dict[(i+1,j)] - I1423Dict[(i+1,j-1)]
                    #I1324: I1324(i,j;N)
                    I1324Dict[(i,j)] += (I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] -I12Dict[(j,L-1)] )*wVal + I1324Dict[(i,j-1)] + I1324Dict[(i+1,j)] - I1324Dict[(i+1,j-1)]
                    for b in range(i+1,j-1): #Obs: the wDict values needed have been computed:
                            wVal1 = wDict[(i,b)]
                            I1234Dict_full_aid[(i,j)] +=  I12Dict[(b+1,j)]*wVal1
                            I1324Dict_full_aid[(i,j)] += I12Dict[(b,j)]*wVal1
                    #I1234Dict_full add up terms and recursion: I1234(i,j) 
                    I1234Dict_full[(i,j)] += I1234Dict_full_aid[(i,j)] -I1234Dict_full_aid[(i,j-1)]
                    I1234Dict_full[(i,j)] +=  I1234Dict_full[(i,j-1)] + I1234Dict_full[(i+1,j)] - I1234Dict_full[(i+1,j-1)]
                    #I1324Dict_full add up terms recursion: I1324(i,j)
                    I1324Dict_full[(i,j)] += I1324Dict_full_aid[(i,j-1)] - I1324Dict_full_aid[(i,j)]   
                    I1324Dict_full[(i,j)] += (I12Dict[(i,j-1)] - I12Dict[(i+1,j-1)])*(I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]) + I1324Dict_full[(i,j-1)] + I1324Dict_full[(i+1,j)] - I1324Dict_full[(i+1,j-1)]                            
                    #I1423_full: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet
                    I1423Dict_full0[(i+1,j)] += (I12Dict[(i+2,j-1)] - I12Dict[(i+2,j)])*(I12Dict[(i+2,L-1)] - I12Dict[(i+2,j)] - I12Dict[(i+1,L-1)] + I12Dict[(i+1,j)])
                    I1423Dict_full0[(i+1,j)] += I1423Dict_full0[(i+2,j)] + I1423Dict_full0[(i+1,j-1)] - I1423Dict_full0[(i+2,j-1)]
                    #to compute certain degree 6 measures we use two auxillary degree 4 measures; the recursion
                    #demands to sum j in the - direction (ie from above and down):
                    k = L-1 + i -j+1
                    for c in range(k+1,L): #Obs: the wDict values needed have already been computed (since k,c >i)
                        wVal2 = wDict[(k,c)]
                        I1324Dict_full2_aid[(i+1,k)] += wVal2*(I12Dict[(i+1,c-1)] - I12Dict[(i+1,k)] - I12Dict[(k,c-1)])
                        I1423Dict_full2_aid[(i+1,k)] += wVal2*(I12Dict[(i+1,L-1)] - I12Dict[(i+1,c)] - I12Dict[(k,L-1)] + I12Dict[(k,c)]) 
                    #I1324Dict_full2, recursion: I1324(i;j,N)
                    I1324Dict_full2[(i+1,k)] += I1324Dict_full2_aid[(i+1,k)] + I1324Dict_full2[(i+1,k+1)]
                    #I1423Dict_full2, recursion: I1423(i;j,N)
                    I1423Dict_full2[(i+1,k)] += I1423Dict_full2_aid[(i+1,k)] + I1423Dict_full2[(i+1,k+1)] 
                    #degree 6 measures:
                    #(12)*:
                    I123456Dict[(i,j)] += I1234Dict[(j+1,L-1)]*wVal + I123456Dict[(i+1,j)] + I123456Dict[(i,j-1)] - I123456Dict[(i+1,j-1)] 
                    I123645Dict[(i,j )] += I1423Dict[(j+1,L-1)]*wVal + I123645Dict[(i+1,j)] + I123645Dict[(i,j-1)] - I123645Dict[(i+1,j-1)] 
                    I123546Dict[(i,j)] += I1324Dict[(j+1,L-1)]*wVal+ I123546Dict[(i+1,j)] + I123546Dict[(i,j-1)] - I123546Dict[(i+1,j-1)] 
                    #(13)*:
                    I132456Dict[(i,j)] += (I1234Dict[(i+1,L-1)] - I1234Dict[(i+1,j)] -I1234Dict[(j,L-1)])*wVal + I132456Dict[(i+1,j)] + I132456Dict[(i,j-1)] - I132456Dict[(i+1,j-1)]
                    I132546Dict[(i+1,j)] += (I1324Dict_full2[(i+2,j+1)] - I1324Dict_full[(j,L-1)])*wDict[(i+1,j)] + I132546Dict[(i+2,j)] + I132546Dict[(i+1,j-1)] - I132546Dict[(i+2,j-1)]
                    I132645Dict[(i+1,j)] += (I1423Dict_full2[(i+2,j+1)] - I1423Dict[(j,L-1)])*wDict[(i+1,j)] + I132645Dict[(i+2,j)] + I132645Dict[(i+1,j-1)] - I132645Dict[(i+2,j-1)]
                    #(14)*
                    #allowing direct computation:
                    I142356Dict[(i,j)] += I12Dict[(i+1,j-1)]*I12Dict[(j+1,L-1)]*wVal + I142356Dict[(i+1,j)] + I142356Dict[(i,j-1)] - I142356Dict[(i+1,j-1)]
                    #back to the more complicated parts, now I142536:
                    #three parts recognized in decompsition as already known:
                    I142536Dict[(i+1,j)] += wDict[(i+1,j)]*(I1324Dict_full[(i+2,L-1)] - I1324Dict[(i+2,j)] - I1324Dict_full2[(i+2,j)])
                    #recursion part:
                    I142536Dict[(i+1,j)] += I142536Dict[(i+2,j)] + I142536Dict[(i+1,j-1)] - I142536Dict[(i+2,j-1)] 
                    #Next complicated part, I142635:
                    #three parts recognized in decompsition as already known:
                    I142635Dict[(i+1,j)] += wDict[(i+1,j)]*(I1423Dict[(i+2,L-1)] - I1423Dict_full2[(i+2,j)] - I1423Dict_full0[(i+2,j)])
                    #recursion part:
                    I142635Dict[(i+1,j)] += I142635Dict[(i+2,j)] + I142635Dict[(i+1,j-1)] - I142635Dict[(i+2,j-1)]
                    #(15)*
                    I152346Dict[(i,j)] += wVal*(I1234Dict[(i+1,j-1)] - I1234Dict_full[(i+1,j)] - I12Dict[(j,L-1)]*I12Dict[(i+1,j-1)])
                    I152346Dict[(i,j)] += I152346Dict[(i+1,j)] + I152346Dict[(i,j-1)] - I152346Dict[(i+1,j-1)]
                    I152436Dict[(i,j)] += wVal*(I1324Dict[(i+1,j-1)] - I1324Dict_full[(i+1,j)])
                    I152436Dict[(i,j)] += I152436Dict[(i+1,j)] + I152436Dict[(i,j-1)] - I152436Dict[(i+1,j-1)]
                    I152634Dict[(i,j)] += wVal*(I1423Dict_full0[(i+1,j-1)] - I1423Dict[(i+1,j)])
                    I152634Dict[(i,j)] += I152634Dict[(i+1,j)] + I152634Dict[(i,j-1)] - I152634Dict[(i+1,j-1)]
                    #(16)*
                    I162345Dict[(i,j)] += wVal*I1234Dict_full[(i+1,j-1)] #ER FULL OK?
                    I162345Dict[(i,j)] += I162345Dict[(i+1,j)] + I162345Dict[(i,j-1)] - I162345Dict[(i+1,j-1)]
                    I162435Dict[(i,j)] += wVal*I1324Dict_full[(i+1,j-1)] #ER FULL OK?
                    I162435Dict[(i,j)] += I162435Dict[(i+1,j)] + I162435Dict[(i,j-1)] - I162435Dict[(i+1,j-1)]
                    I162534Dict[(i,j)] += wVal*I1423Dict[(i+1,j-1)]
                    I162534Dict[(i,j)] += I162534Dict[(i+1,j)] + I162534Dict[(i,j-1)] - I162534Dict[(i+1,j-1)]                    
        #the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:
        for j in lChain[::-1]:
            i = 0
            for c in range(j+1,L):
                wVal = wDict[(j,c)]
#                        I1234Dict_full2_aid[(i,j)] +=  I12Dict[(i,c-1)]*wDict[(j,c)]
                I1324Dict_full2_aid[(i,j)] += wVal*(I12Dict[(i,c-1)] - I12Dict[(i,j)] - I12Dict[(j,c-1)])
                I1423Dict_full2_aid[(i,j)] += wVal*(I12Dict[(i,L-1)] - I12Dict[(i,c)] - I12Dict[(j,L-1)] + I12Dict[(j,c)]) 
           #recursion:
            I1324Dict_full2[(i,j)] += I1324Dict_full2_aid[(i,j)] + I1324Dict_full2[(i,j+1)]
            I1423Dict_full2[(i,j)] += I1423Dict_full2_aid[(i,j)] + I1423Dict_full2[(i,j+1)] 
        #Missing terms for some of the degree 6 measures and I1423_full0:
        for j in lChain[1:]:
            i=0
            wVal = wDict[(i,j)]
            I1423Dict_full0[(i,j)] += (I12Dict[(i+1,j-1)] - I12Dict[(i+1,j)])*(I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] - I12Dict[(i,L-1)] + I12Dict[(i,j)])
            I1423Dict_full0[(i,j)] += I1423Dict_full0[(i+1,j)] + I1423Dict_full0[(i,j-1)] - I1423Dict_full0[(i+1,j-1)]
            #
            I132546Dict[(i,j)] += (I1324Dict_full2[(i+1,j+1)] - I1324Dict_full[(j,L-1)])*wVal + I132546Dict[(i+1,j)] + I132546Dict[(i,j-1)] - I132546Dict[(i+1,j-1)]
            I132645Dict[(i,j)] += (I1423Dict_full2[(i+1,j+1)] - I1423Dict[(j,L-1)])*wVal + I132645Dict[(i+1,j)] + I132645Dict[(i,j-1)] - I132645Dict[(i+1,j-1)]
            #
            I142536Dict[(i,j)] += wVal*(I1324Dict_full[(i+1,L-1)] - I1324Dict[(i+1,j)] - I1324Dict_full2[(i+1,j)])
            I142536Dict[(i,j)] += I142536Dict[(i+1,j)] + I142536Dict[(i,j-1)] - I142536Dict[(i+1,j-1)] 
            #
            I142635Dict[(i,j)] += wVal*(I1423Dict[(i+1,L-1)] - I1423Dict_full2[(i+1,j)] - I1423Dict_full0[(i+1,j)])
            I142635Dict[(i,j)] += I142635Dict[(i+1,j)] + I142635Dict[(i,j-1)] - I142635Dict[(i+1,j-1)]    
    else: 
        print "Degree must be 2, 4 or 6."
    outDict = {}
    outDict['w'] = wDict
    outDict['I12'] = I12Dict
    outDict['I1234'] = I1234Dict
    outDict['I1234_full'] = I1234Dict_full
    outDict['I1324'] = I1324Dict
    outDict['I1324_full'] = I1324Dict_full
    outDict['I1324_full2'] = I1324Dict_full2    
    outDict['I1423'] = I1423Dict
    outDict['I1423_full0'] = I1423Dict_full0
    outDict['I1423_full2'] = I1423Dict_full2    

    outDict['I123456'] = I123456Dict
    outDict['I123546'] = I123546Dict
    outDict['I123645'] = I123645Dict

    outDict['I132456'] = I132456Dict
    outDict['I132546'] = I132546Dict
    outDict['I132645'] = I132645Dict

    outDict['I142356'] = I142356Dict
    outDict['I142536'] = I142536Dict
    outDict['I142635'] = I142635Dict
    
    outDict['I152346'] = I152346Dict
    outDict['I152436'] = I152436Dict
    outDict['I152634'] = I152634Dict

    outDict['I162345'] = I162345Dict
    outDict['I162435'] = I162435Dict
    outDict['I162534'] = I162534Dict
    return outDict #wDict, I12, outDict    
    


#@jit(float64[:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:](float64[:,:,:],float64[:,:,:], int64, int64, int64))
def aggr_vect(CaChain, pChain, degree, full_b, print_b, aggrTime_b):
    '''The computation  called in the function I. Input parameters not
    dealt with there: pChain -- a polyCaChain; outDict: both in/output; the 
    dictionary returned in the I function (see there for more).'''
    #Get length of chain etc:
    L = len(pChain)
#    print_b = 0
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    lChain = range(L)
    if print_b > 0.5:
        print lChain
    
    #Initialize array parameters for storing the results (the 'Dict' suffix is 
    #used for easy translation from dictionary-based version):
    I12 = 0 #value of I12 on full chain
    I12Dict = np.zeros((L+1,L+1))
    I1234Dict = np.zeros((L+1,L+1))
    I1234Dict_full = np.zeros((L+1,L+1))
    I1234Dict_full_aid = np.zeros((L+1,L+1))
    I1234Dict_full2 = np.zeros((L+1,L+1))
    I1234Dict_full2_aid = np.zeros((L+1,L+1))
    I1423Dict = np.zeros((L+1,L+1)) 
    I1423Dict_full0 = np.zeros((L+1,L+1)) 
    I1423Dict_full2 = np.zeros((L+1,L+1))
    I1423Dict_full2_aid = np.zeros((L+1,L+1))
    I1324Dict = np.zeros((L+1,L+1))
    I1324Dict_full = np.zeros((L+1,L+1))
    I1324Dict_full_aid = np.zeros((L+1,L+1))
    I1324Dict_full2 = np.zeros((L+1,L+1))
    I1324Dict_full2_aid = np.zeros((L+1,L+1))
    #(12)*:
    I123456Dict = np.zeros((L+1,L+1)) 
    I123645Dict = np.zeros((L+1,L+1))
    I123546Dict = np.zeros((L+1,L+1))
    #(13)*:
    I132456Dict = np.zeros((L+1,L+1))
    I132546Dict = np.zeros((L+1,L+1))
    I132645Dict = np.zeros((L+1,L+1))
    #(14)*
    I142356Dict = np.zeros((L+1,L+1))
    I142536Dict = np.zeros((L+1,L+1))
    I142635Dict = np.zeros((L+1,L+1))
    #(15)*
    I152346Dict = np.zeros((L+1,L+1))
    I152436Dict = np.zeros((L+1,L+1))
    I152634Dict = np.zeros((L+1,L+1))
    #(16)*
    I162345Dict = np.zeros((L+1,L+1))
    I162435Dict = np.zeros((L+1,L+1))
    I162534Dict = np.zeros((L+1,L+1))

    start = time.time()
    
    #compute the writhe terms in one go (call for array-structured output by setting dict_b = 0):
    wDict = wAll(CaChain = CaChain, polyCaChain = pChain, dict_b = 0) 
    
    #Main part: computing the measures using recursion formulas
    #one block for each choice of degree; only one block gets executed:   
    if degree == 0:
        if print_b > 0.5:
            print "Computes w-values only."
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = wDict[i,j]
    elif degree == 2:
        if print_b > 0.5:
            print "Computes degree 2 measures only."
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = wDict[i,j]
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
    elif degree == 4:
        if print_b > 0.5:
            print "Computes degree 2 and 4 measures only."
        for i in lChain[::-1]:  
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = wDict[i,j]
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
                    #I1234: I1234(i,j;N)
                    I1234Dict[i,j] += I12Dict[j+1,L-1]*wVal + I1234Dict[i,j-1] + I1234Dict[i+1,j] - I1234Dict[i+1,j-1]
                    #I1423: I1423(i,j)
                    I1423Dict[i,j] += I12Dict[i+1,j-1]*wVal +  I1423Dict[i,j-1] + I1423Dict[i+1,j] - I1423Dict[i+1,j-1]
                    #I1324: I1324(i,j;N)
                    I1324Dict[i,j] += (I12Dict[i+1,L-1] - I12Dict[i+1,j] -I12Dict[j,L-1] )*wVal + I1324Dict[i,j-1] + I1324Dict[i+1,j] - I1324Dict[i+1,j-1]
                    if full_b == 1:
                        for b in range(i+1,j-1): #obs: the wDict values needed have been computed
                            wVal1 = wDict[i,b]
                            I1234Dict_full_aid[i,j] +=  I12Dict[b+1,j]*wVal1
                            I1324Dict_full_aid[i,j] += I12Dict[b,j]*wVal1
#                            I1234Dict_full[(i,j)] += (I12Dict[(b+1,j)] - I12Dict[(b+1,j-1)])*wDict[(i,b)]
#                            I1324Dict_full[(i,j)] += (-I12Dict[(b,j)] + I12Dict[(b,j-1)])*wDict[(i,b)]
                        #I1234Dict_full add up terms and recursion: I1234(i,j)    
                        I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                        I1234Dict_full[i,j] +=  I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                        #I1324Dict_full add up terms recursion: I1324(i,j)
                        I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                        I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]                                                   
    elif degree == 6:
        if print_b > 0.5:
            print "Computes degree 2, 4 and 6 measures." 
        for i in lChain[::-1]:  
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = wDict[i,j]
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
                    #degree 4 measures:
                    #I1234: I1234(i,j;N)
                    I1234Dict[i,j] += I12Dict[j+1,L-1]*wVal + I1234Dict[i,j-1] + I1234Dict[i+1,j] - I1234Dict[i+1,j-1]
                    #I1423: I1423(i,j)                    
                    I1423Dict[i,j] += I12Dict[i+1,j-1]*wVal +  I1423Dict[i,j-1] + I1423Dict[i+1,j] - I1423Dict[i+1,j-1]
                    #I1324: I1324(i,j;N)
                    I1324Dict[i,j] += (I12Dict[i+1,L-1] - I12Dict[i+1,j] -I12Dict[j,L-1] )*wVal + I1324Dict[i,j-1] + I1324Dict[i+1,j] - I1324Dict[i+1,j-1]
                    for b in range(i+1,j-1): #obs: the wDict values needed have been computed
                        wVal1 = wDict[i,b]
                        I1234Dict_full_aid[i,j] +=  I12Dict[b+1,j]*wVal1
                        I1324Dict_full_aid[i,j] += I12Dict[b,j]*wVal1
                    #I1234Dict_full add up terms and recursion: I1234(i,j))    
                    I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                    I1234Dict_full[i,j] +=  I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                    #I1324Dict_full add up terms recursion: I1324(i,j)
                    I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                    I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]                            
                    #I1423_full: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet
                    I1423Dict_full0[i+1,j] += (I12Dict[i+2,j-1] - I12Dict[i+2,j])*(I12Dict[i+2,L-1] - I12Dict[i+2,j] - I12Dict[i+1,L-1] + I12Dict[i+1,j])
                    I1423Dict_full0[i+1,j] += I1423Dict_full0[i+2,j] + I1423Dict_full0[i+1,j-1] - I1423Dict_full0[i+2,j-1]
                    #to compute certain degree 6 measures we use two auxillary degree 4 
                    #measures; the recursion demands to sum j in the - direction (ie 
                    #from above and down):
                    k = L-1 + i -j+1
                    for c in range(k+1,L): #Obs: the wDict values needed have already been computed (since k,c >i)
                        wVal2 = wDict[k,c]
                        I1324Dict_full2_aid[i+1,k] += wVal2*(I12Dict[i+1,c-1] - I12Dict[i+1,k] - I12Dict[k,c-1])
                        I1423Dict_full2_aid[i+1,k] += wVal2*(I12Dict[i+1,L-1] - I12Dict[i+1,c] - I12Dict[k,L-1] + I12Dict[k,c]) 
                    #I1324Dict_full2, recursion: I1324(i;j,N)
                    I1324Dict_full2[i+1,k] += I1324Dict_full2_aid[i+1,k] + I1324Dict_full2[i+1,k+1]
                    #I1423Dict_full2, recursion: I1423(i;j,N)
                    I1423Dict_full2[i+1,k] += I1423Dict_full2_aid[i+1,k] + I1423Dict_full2[i+1,k+1] 
                    #degree 6 measures:
                    #(12)*:
                    I123456Dict[i,j] += I1234Dict[j+1,L-1]*wVal + I123456Dict[i+1,j] + I123456Dict[i,j-1] - I123456Dict[i+1,j-1] 
                    I123645Dict[i,j] += I1423Dict[j+1,L-1]*wVal + I123645Dict[i+1,j] + I123645Dict[i,j-1] - I123645Dict[i+1,j-1] 
                    I123546Dict[i,j] += I1324Dict[j+1,L-1]*wVal+ I123546Dict[i+1,j] + I123546Dict[i,j-1] - I123546Dict[i+1,j-1] 
                    #(13)*:
                    I132456Dict[i,j] += (I1234Dict[i+1,L-1] - I1234Dict[i+1,j] -I1234Dict[j,L-1])*wVal + I132456Dict[i+1,j] + I132456Dict[i,j-1] - I132456Dict[i+1,j-1]
                    I132546Dict[i+1,j] += (I1324Dict_full2[i+2,j+1] - I1324Dict_full[j,L-1])*wDict[i+1,j] + I132546Dict[i+2,j] + I132546Dict[i+1,j-1] - I132546Dict[i+2,j-1]
                    I132645Dict[i+1,j] += (I1423Dict_full2[i+2,j+1] - I1423Dict[j,L-1])*wDict[i+1,j] + I132645Dict[i+2,j] + I132645Dict[i+1,j-1] - I132645Dict[i+2,j-1]
                    #(14)*
                    #allowing direct computation:
                    I142356Dict[i,j] += I12Dict[i+1,j-1]*I12Dict[j+1,L-1]*wVal + I142356Dict[i+1,j] + I142356Dict[i,j-1] - I142356Dict[i+1,j-1]
                    #back to the more complicated parts, now I142536:
                    #three parts recognized in decompsition as already known:
                    I142536Dict[i+1,j] += wDict[i+1,j]*(I1324Dict_full[i+2,L-1] - I1324Dict[i+2,j] - I1324Dict_full2[i+2,j])
                    #recursion part:
                    I142536Dict[i+1,j] += I142536Dict[i+2,j] + I142536Dict[i+1,j-1] - I142536Dict[i+2,j-1] 
                    #Next complicated part, I142635:
                    #three parts recognized in decompsition as already known:
                    I142635Dict[i+1,j] += wDict[i+1,j]*(I1423Dict[i+2,L-1] - I1423Dict_full2[i+2,j] - I1423Dict_full0[i+2,j])
                    #recursion part:
                    I142635Dict[i+1,j] += I142635Dict[i+2,j] + I142635Dict[i+1,j-1] - I142635Dict[i+2,j-1]
                    #(15)*
                    I152346Dict[i,j] += wVal*(I1234Dict[i+1,j-1] - I1234Dict_full[i+1,j] - I12Dict[j,L-1]*I12Dict[i+1,j-1])
                    I152346Dict[i,j] += I152346Dict[i+1,j] + I152346Dict[i,j-1] - I152346Dict[i+1,j-1]
                    I152436Dict[i,j] += wVal*(I1324Dict[i+1,j-1] - I1324Dict_full[i+1,j])
                    I152436Dict[i,j] += I152436Dict[i+1,j] + I152436Dict[i,j-1] - I152436Dict[i+1,j-1]
                    I152634Dict[i,j] += wVal*(I1423Dict_full0[i+1,j-1] - I1423Dict[i+1,j])
                    I152634Dict[i,j] += I152634Dict[i+1,j] + I152634Dict[i,j-1] - I152634Dict[i+1,j-1]
                    #(16)*
                    I162345Dict[i,j] += wVal*I1234Dict_full[i+1,j-1] #ER FULL OK?
                    I162345Dict[i,j] += I162345Dict[i+1,j] + I162345Dict[i,j-1] - I162345Dict[i+1,j-1]
                    I162435Dict[i,j] += wVal*I1324Dict_full[i+1,j-1] #ER FULL OK?
                    I162435Dict[i,j] += I162435Dict[i+1,j] + I162435Dict[i,j-1] - I162435Dict[i+1,j-1]
                    I162534Dict[i,j] += wVal*I1423Dict[i+1,j-1]
                    I162534Dict[i,j] += I162534Dict[i+1,j] + I162534Dict[i,j-1] - I162534Dict[i+1,j-1]                    
        #the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:
        for j in lChain[::-1]:
            i = 0
            for c in range(j+1,L):
                wVal = wDict[j,c]
#                        I1234Dict_full2_aid[(i,j)] +=  I12Dict[(i,c-1)]*wDict[(j,c)]
                I1324Dict_full2_aid[i,j] += wVal*(I12Dict[i,c-1] - I12Dict[i,j] - I12Dict[j,c-1])
                I1423Dict_full2_aid[i,j] += wVal*(I12Dict[i,L-1] - I12Dict[i,c] - I12Dict[j,L-1] + I12Dict[j,c]) 
           #recursion:
            I1324Dict_full2[i,j] += I1324Dict_full2_aid[i,j] + I1324Dict_full2[i,j+1]
            I1423Dict_full2[i,j] += I1423Dict_full2_aid[i,j] + I1423Dict_full2[i,j+1] 
        #Missing terms for some of the degree 6 measures and I1423_full0:
        for j in lChain[1:]:
            i=0
            wVal = wDict[i,j]
            I1423Dict_full0[i,j] += (I12Dict[i+1,j-1] - I12Dict[i+1,j])*(I12Dict[i+1,L-1] - I12Dict[i+1,j] - I12Dict[i,L-1] + I12Dict[i,j])
            I1423Dict_full0[i,j] += I1423Dict_full0[i+1,j] + I1423Dict_full0[i,j-1] - I1423Dict_full0[i+1,j-1]
            #
            I132546Dict[i,j] += (I1324Dict_full2[i+1,j+1] - I1324Dict_full[j,L-1])*wVal + I132546Dict[i+1,j] + I132546Dict[i,j-1] - I132546Dict[i+1,j-1]
            I132645Dict[i,j] += (I1423Dict_full2[i+1,j+1] - I1423Dict[j,L-1])*wVal + I132645Dict[i+1,j] + I132645Dict[i,j-1] - I132645Dict[i+1,j-1]
            #
            I142536Dict[i,j] += wVal*(I1324Dict_full[i+1,L-1] - I1324Dict[i+1,j] - I1324Dict_full2[i+1,j])
            I142536Dict[i,j] += I142536Dict[i+1,j] + I142536Dict[i,j-1] - I142536Dict[i+1,j-1] 
            #
            I142635Dict[i,j] += wVal*(I1423Dict[i+1,L-1] - I1423Dict_full2[i+1,j] - I1423Dict_full0[i+1,j])
            I142635Dict[i,j] += I142635Dict[i+1,j] + I142635Dict[i,j-1] - I142635Dict[i+1,j-1]    
    else: 
        print "Degree must be 2, 4 or 6."
#    outDict = np.zeros((L+1,L+1))
#    outDict = np.array((wDict, I12Dict, I1234Dict, I1324Dict, I1423Dict)) #, I1234Dict, I1324Dict, I1423Dict)
    outDict = np.array((wDict, I12Dict, 
                        I1234Dict, I1234Dict_full, 
                        I1324Dict, I1324Dict_full, I1324Dict_full2,
                        I1423Dict, I1423Dict_full0,I1423Dict_full2,
                        I123456Dict, I123546Dict, I123645Dict,
                        I132456Dict, I132546Dict, I132645Dict,
                        I142356Dict, I142536Dict, I142635Dict,
                        I152346Dict, I152436Dict, I152634Dict,
                        I162345Dict, I162435Dict, I162534Dict)) #, I1234Dict, I1324Dict, I1423Dict)
#    return wDict #, I12, outDict
                        
    end = time.time()

    if aggrTime_b != 1:
        return outDict
    elif aggrTime_b ==1:
        return outDict, end - start
        

#"naked" version excluding the initialization of arrays for holdinng the invariants' values;
#aimed at use in "global memory allocation" version:
def aggr_vect_naked(CaChain, 
                    pChain, 
                    degree, 
                    full_b, 
                    print_b, 
                    aggrTime_b,
                    wDict, #values of w-terms
                    I12, #value of I12 on full chain
                    I12Dict,
                    I1234Dict,
                    I1234Dict_full,
                    I1234Dict_full_aid,
                    I1234Dict_full2,
                    I1234Dict_full2_aid,
                    I1423Dict, 
                    I1423Dict_full0, 
                    I1423Dict_full2,
                    I1423Dict_full2_aid,
                    I1324Dict,
                    I1324Dict_full,
                    I1324Dict_full_aid,
                    I1324Dict_full2,
                    I1324Dict_full2_aid,
                    #(12)*:
                    I123456Dict, 
                    I123645Dict,
                    I123546Dict,
                    #(13)*:
                    I132456Dict,
                    I132546Dict,
                    I132645Dict,
                    #(14)*
                    I142356Dict,
                    I142536Dict,
                    I142635Dict,
                    #(15)*
                    I152346Dict,
                    I152436Dict,
                    I152634Dict,
                    #(16)*
                    I162345Dict,
                    I162435Dict,
                    I162534Dict):
    '''The computation  called in the function I. Input parameters not
    dealt with there: pChain -- a polyCaChain; outDict: both in/output; the 
    dictionary returned in the I function (see there for more).'''
    #Get length of chain etc:
    L = len(pChain)
#    print_b = 0
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    lChain = range(L)
    if print_b > 0.5:
        print lChain
    
    #(re)init:
    for i in range(L+1):
        for j in range(L+1):
            wDict[i,j] = 0
            if degree > 0:
                I12Dict[i,j] = 0
#            I12Dict[(-1,j)] = 0
            if degree > 2: 
                I1234Dict[i,j] = 0
                I1324Dict[i,j] = 0  
                I1423Dict[i,j] = 0
                if (full_b ==1 and degree < 6):
                    I1234Dict_full[i,j] = 0
                    I1234Dict_full_aid[i,j] = 0
                    I1324Dict_full[i,j] = 0
                    I1324Dict_full_aid[i,j] = 0
                    
            if degree > 4:
                
                I1234Dict_full[i,j] = 0
                I1234Dict_full_aid[i,j] = 0
                I1234Dict_full2[i,j] = 0
                I1234Dict_full2_aid[i,j] = 0
                I1324Dict_full[i,j] = 0
                I1324Dict_full_aid[i,j] = 0
                I1324Dict_full2[i,j] = 0
#                I1324Dict_full2[(-1,j)] = 0
                I1324Dict_full2_aid[i,j] = 0
#                I1324Dict_full2_aid[(-1,j)] = 0
                I1423Dict_full0[i,j] = 0
                I1423Dict_full2[i,j] = 0
                I1423Dict_full2_aid[i,j] = 0                
                
                I123456Dict[i,j] = 0 
                I123645Dict[i,j] = 0
                I123546Dict[i,j] = 0
                #(13)*:
                I132456Dict[i,j] = 0
                I132546Dict[i,j] = 0
                I132645Dict[i,j] = 0
                #(14)*
                I142356Dict[i,j] = 0
                I142536Dict[i,j] = 0
                I142635Dict[i,j] = 0
                #(15)*
                I152346Dict[i,j] = 0
                I152436Dict[i,j] = 0
                I152634Dict[i,j] = 0
                #(16)*
                I162345Dict[i,j] = 0
                I162435Dict[i,j] = 0
                I162534Dict[i,j] = 0    

    #we want to measure the naked computation time so disregard the above initialization:
    start = time.time()
    
    #compute the writhe terms in one go (call for array-structured output by setting dict_b = 0):
    wDict = wAll(CaChain = CaChain, polyCaChain = pChain, dict_b = 0, global_b = 1, wDict = wDict) 
    
    #Main part: computing the measures using recursion formulas
    #one block for each choice of degree; only one block gets executed:   
    if degree == 0:
        if print_b > 0.5:
            print "Computes w-values only."
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = wDict[i,j]
    elif degree == 2:
        if print_b > 0.5:
            print "Computes degree 2 measures only."
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = wDict[i,j]
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
    elif degree == 4:
        if print_b > 0.5:
            print "Computes degree 2 and 4 measures only."
        for i in lChain[::-1]:  
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = wDict[i,j]
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
                    #I1234: I1234(i,j;N)
                    I1234Dict[i,j] += I12Dict[j+1,L-1]*wVal + I1234Dict[i,j-1] + I1234Dict[i+1,j] - I1234Dict[i+1,j-1]
                    #I1423: I1423(i,j)
                    I1423Dict[i,j] += I12Dict[i+1,j-1]*wVal +  I1423Dict[i,j-1] + I1423Dict[i+1,j] - I1423Dict[i+1,j-1]
                    #I1324: I1324(i,j;N)
                    I1324Dict[i,j] += (I12Dict[i+1,L-1] - I12Dict[i+1,j] -I12Dict[j,L-1] )*wVal + I1324Dict[i,j-1] + I1324Dict[i+1,j] - I1324Dict[i+1,j-1]
                    if full_b == 1:
                        for b in range(i+1,j-1): #obs: the wDict values needed have been computed
                            wVal1 = wDict[i,b]
                            I1234Dict_full_aid[i,j] +=  I12Dict[b+1,j]*wVal1
                            I1324Dict_full_aid[i,j] += I12Dict[b,j]*wVal1
#                            I1234Dict_full[(i,j)] += (I12Dict[(b+1,j)] - I12Dict[(b+1,j-1)])*wDict[(i,b)]
#                            I1324Dict_full[(i,j)] += (-I12Dict[(b,j)] + I12Dict[(b,j-1)])*wDict[(i,b)]
                        #I1234Dict_full add up terms and recursion: I1234(i,j)    
                        I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                        I1234Dict_full[i,j] +=  I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                        #I1324Dict_full add up terms recursion: I1324(i,j)
                        I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                        I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]                                                   
    elif degree == 6:
        if print_b > 0.5:
            print "Computes degree 2, 4 and 6 measures." 
        for i in lChain[::-1]:  
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    wVal = wDict[i,j]
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
                    #degree 4 measures:
                    #I1234: I1234(i,j;N)
                    I1234Dict[i,j] += I12Dict[j+1,L-1]*wVal + I1234Dict[i,j-1] + I1234Dict[i+1,j] - I1234Dict[i+1,j-1]
                    #I1423: I1423(i,j)                    
                    I1423Dict[i,j] += I12Dict[i+1,j-1]*wVal +  I1423Dict[i,j-1] + I1423Dict[i+1,j] - I1423Dict[i+1,j-1]
                    #I1324: I1324(i,j;N)
                    I1324Dict[i,j] += (I12Dict[i+1,L-1] - I12Dict[i+1,j] -I12Dict[j,L-1] )*wVal + I1324Dict[i,j-1] + I1324Dict[i+1,j] - I1324Dict[i+1,j-1]
                    for b in range(i+1,j-1): #obs: the wDict values needed have been computed
                        wVal1 = wDict[i,b]
                        I1234Dict_full_aid[i,j] +=  I12Dict[b+1,j]*wVal1
                        I1324Dict_full_aid[i,j] += I12Dict[b,j]*wVal1
                    #I1234Dict_full add up terms and recursion: I1234(i,j))    
                    I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                    I1234Dict_full[i,j] +=  I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                    #I1324Dict_full add up terms recursion: I1324(i,j)
                    I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                    I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]                            
                    #I1423_full: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet
                    I1423Dict_full0[i+1,j] += (I12Dict[i+2,j-1] - I12Dict[i+2,j])*(I12Dict[i+2,L-1] - I12Dict[i+2,j] - I12Dict[i+1,L-1] + I12Dict[i+1,j])
                    I1423Dict_full0[i+1,j] += I1423Dict_full0[i+2,j] + I1423Dict_full0[i+1,j-1] - I1423Dict_full0[i+2,j-1]
                    #to compute certain degree 6 measures we use two auxillary degree 4 
                    #measures; the recursion demands to sum j in the - direction (ie 
                    #from above and down):
                    k = L-1 + i -j+1
                    for c in range(k+1,L): #Obs: the wDict values needed have already been computed (since k,c >i)
                        wVal2 = wDict[k,c]
                        I1324Dict_full2_aid[i+1,k] += wVal2*(I12Dict[i+1,c-1] - I12Dict[i+1,k] - I12Dict[k,c-1])
                        I1423Dict_full2_aid[i+1,k] += wVal2*(I12Dict[i+1,L-1] - I12Dict[i+1,c] - I12Dict[k,L-1] + I12Dict[k,c]) 
                    #I1324Dict_full2, recursion: I1324(i;j,N)
                    I1324Dict_full2[i+1,k] += I1324Dict_full2_aid[i+1,k] + I1324Dict_full2[i+1,k+1]
                    #I1423Dict_full2, recursion: I1423(i;j,N)
                    I1423Dict_full2[i+1,k] += I1423Dict_full2_aid[i+1,k] + I1423Dict_full2[i+1,k+1] 
                    #degree 6 measures:
                    #(12)*:
                    I123456Dict[i,j] += I1234Dict[j+1,L-1]*wVal + I123456Dict[i+1,j] + I123456Dict[i,j-1] - I123456Dict[i+1,j-1] 
                    I123645Dict[i,j] += I1423Dict[j+1,L-1]*wVal + I123645Dict[i+1,j] + I123645Dict[i,j-1] - I123645Dict[i+1,j-1] 
                    I123546Dict[i,j] += I1324Dict[j+1,L-1]*wVal+ I123546Dict[i+1,j] + I123546Dict[i,j-1] - I123546Dict[i+1,j-1] 
                    #(13)*:
                    I132456Dict[i,j] += (I1234Dict[i+1,L-1] - I1234Dict[i+1,j] -I1234Dict[j,L-1])*wVal + I132456Dict[i+1,j] + I132456Dict[i,j-1] - I132456Dict[i+1,j-1]
                    I132546Dict[i+1,j] += (I1324Dict_full2[i+2,j+1] - I1324Dict_full[j,L-1])*wDict[i+1,j] + I132546Dict[i+2,j] + I132546Dict[i+1,j-1] - I132546Dict[i+2,j-1]
                    I132645Dict[i+1,j] += (I1423Dict_full2[i+2,j+1] - I1423Dict[j,L-1])*wDict[i+1,j] + I132645Dict[i+2,j] + I132645Dict[i+1,j-1] - I132645Dict[i+2,j-1]
                    #(14)*
                    #allowing direct computation:
                    I142356Dict[i,j] += I12Dict[i+1,j-1]*I12Dict[j+1,L-1]*wVal + I142356Dict[i+1,j] + I142356Dict[i,j-1] - I142356Dict[i+1,j-1]
                    #back to the more complicated parts, now I142536:
                    #three parts recognized in decompsition as already known:
                    I142536Dict[i+1,j] += wDict[i+1,j]*(I1324Dict_full[i+2,L-1] - I1324Dict[i+2,j] - I1324Dict_full2[i+2,j])
                    #recursion part:
                    I142536Dict[i+1,j] += I142536Dict[i+2,j] + I142536Dict[i+1,j-1] - I142536Dict[i+2,j-1] 
                    #Next complicated part, I142635:
                    #three parts recognized in decompsition as already known:
                    I142635Dict[i+1,j] += wDict[i+1,j]*(I1423Dict[i+2,L-1] - I1423Dict_full2[i+2,j] - I1423Dict_full0[i+2,j])
                    #recursion part:
                    I142635Dict[i+1,j] += I142635Dict[i+2,j] + I142635Dict[i+1,j-1] - I142635Dict[i+2,j-1]
                    #(15)*
                    I152346Dict[i,j] += wVal*(I1234Dict[i+1,j-1] - I1234Dict_full[i+1,j] - I12Dict[j,L-1]*I12Dict[i+1,j-1])
                    I152346Dict[i,j] += I152346Dict[i+1,j] + I152346Dict[i,j-1] - I152346Dict[i+1,j-1]
                    I152436Dict[i,j] += wVal*(I1324Dict[i+1,j-1] - I1324Dict_full[i+1,j])
                    I152436Dict[i,j] += I152436Dict[i+1,j] + I152436Dict[i,j-1] - I152436Dict[i+1,j-1]
                    I152634Dict[i,j] += wVal*(I1423Dict_full0[i+1,j-1] - I1423Dict[i+1,j])
                    I152634Dict[i,j] += I152634Dict[i+1,j] + I152634Dict[i,j-1] - I152634Dict[i+1,j-1]
                    #(16)*
                    I162345Dict[i,j] += wVal*I1234Dict_full[i+1,j-1] #ER FULL OK?
                    I162345Dict[i,j] += I162345Dict[i+1,j] + I162345Dict[i,j-1] - I162345Dict[i+1,j-1]
                    I162435Dict[i,j] += wVal*I1324Dict_full[i+1,j-1] #ER FULL OK?
                    I162435Dict[i,j] += I162435Dict[i+1,j] + I162435Dict[i,j-1] - I162435Dict[i+1,j-1]
                    I162534Dict[i,j] += wVal*I1423Dict[i+1,j-1]
                    I162534Dict[i,j] += I162534Dict[i+1,j] + I162534Dict[i,j-1] - I162534Dict[i+1,j-1]                    
        #the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:
        for j in lChain[::-1]:
            i = 0
            for c in range(j+1,L):
                wVal = wDict[j,c]
#                        I1234Dict_full2_aid[(i,j)] +=  I12Dict[(i,c-1)]*wDict[(j,c)]
                I1324Dict_full2_aid[i,j] += wVal*(I12Dict[i,c-1] - I12Dict[i,j] - I12Dict[j,c-1])
                I1423Dict_full2_aid[i,j] += wVal*(I12Dict[i,L-1] - I12Dict[i,c] - I12Dict[j,L-1] + I12Dict[j,c]) 
           #recursion:
            I1324Dict_full2[i,j] += I1324Dict_full2_aid[i,j] + I1324Dict_full2[i,j+1]
            I1423Dict_full2[i,j] += I1423Dict_full2_aid[i,j] + I1423Dict_full2[i,j+1] 
        #Missing terms for some of the degree 6 measures and I1423_full0:
        for j in lChain[1:]:
            i=0
            wVal = wDict[i,j]
            I1423Dict_full0[i,j] += (I12Dict[i+1,j-1] - I12Dict[i+1,j])*(I12Dict[i+1,L-1] - I12Dict[i+1,j] - I12Dict[i,L-1] + I12Dict[i,j])
            I1423Dict_full0[i,j] += I1423Dict_full0[i+1,j] + I1423Dict_full0[i,j-1] - I1423Dict_full0[i+1,j-1]
            #
            I132546Dict[i,j] += (I1324Dict_full2[i+1,j+1] - I1324Dict_full[j,L-1])*wVal + I132546Dict[i+1,j] + I132546Dict[i,j-1] - I132546Dict[i+1,j-1]
            I132645Dict[i,j] += (I1423Dict_full2[i+1,j+1] - I1423Dict[j,L-1])*wVal + I132645Dict[i+1,j] + I132645Dict[i,j-1] - I132645Dict[i+1,j-1]
            #
            I142536Dict[i,j] += wVal*(I1324Dict_full[i+1,L-1] - I1324Dict[i+1,j] - I1324Dict_full2[i+1,j])
            I142536Dict[i,j] += I142536Dict[i+1,j] + I142536Dict[i,j-1] - I142536Dict[i+1,j-1] 
            #
            I142635Dict[i,j] += wVal*(I1423Dict[i+1,L-1] - I1423Dict_full2[i+1,j] - I1423Dict_full0[i+1,j])
            I142635Dict[i,j] += I142635Dict[i+1,j] + I142635Dict[i,j-1] - I142635Dict[i+1,j-1]    
    else: 
        print "Degree must be 2, 4 or 6."
#    outDict = np.zeros((L+1,L+1))
#    outDict = np.array((wDict, I12Dict, I1234Dict, I1324Dict, I1423Dict)) #, I1234Dict, I1324Dict, I1423Dict)
    outDict = np.array((wDict, I12Dict, 
                        I1234Dict, I1234Dict_full, 
                        I1324Dict, I1324Dict_full, I1324Dict_full2,
                        I1423Dict, I1423Dict_full0,I1423Dict_full2,
                        I123456Dict, I123546Dict, I123645Dict,
                        I132456Dict, I132546Dict, I132645Dict,
                        I142356Dict, I142536Dict, I142635Dict,
                        I152346Dict, I152436Dict, I152634Dict,
                        I162345Dict, I162435Dict, I162534Dict)) #, I1234Dict, I1324Dict, I1423Dict)
#    return wDict #, I12, outDict
                        
    end = time.time()

    if aggrTime_b != 1:
        return outDict
    elif aggrTime_b ==1:
        return outDict, end - start



#######################################################
### Vectorized based functions for perturbations
#######################################################

#def Iperturbed_vect(CaChain, 
#                polyCaChain, 
#                wDict,
#                residueNrRange = range(0,5),
#                perturbation = [(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3)], 
#                degree = 2,
#                full_b = 1,
#                print_b = 1):
#    '''Vectorized version of Iperturbed.
#    Array based version of Iperturbed_dict.
#    Computes the change in the measures of the given polychain due to 
#    the perturbation. Input: 
#    CaChain: C-alpha chain
#    polyCaChain: chain of segment corresponding to CaChain
#    wDict: array of w-values as output from I
#    residueNrRange: indices of the residues in the CaChain (C-alphas) that should be perturbed
#    perturbation: the perturbation given in 3d-coordinates. Must be of lenght equal to residueNrRange
#    Other: as I_dict.
#    Returns: same as I, but for the perturbed chain.'''
#    #compute the polygonal chain of the perturbed chain:
#    perCaChain, perPolyCaChain, segNrRange = perturbPolygonalChain(CaChain = CaChain, 
#                                                    polyCaChain = polyCaChain,
#                     residueNrRange = residueNrRange, 
#                     perturbation = perturbation)



#@jit(float64[:,:,:,:,:,:,:,:,:,:](float64[:,:,:], int64, int64))
def aggr_pert_vect(wDict, perCaChain, perPolyCaChain, segNrRange, degree, full_b, print_b):
    '''As aggr_pert, but vectorized. 
    The computation  called in the function I_perturbed, vectorized. 
    Input: passed from I_perturbed (for more see there).'''
    if print_b > 0.5:
        print segNrRange
    #Get length of chain:
    L = len(perPolyCaChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    lChain = range(L)
    #call writhe computation (w array structured output; call with dict_b = 0):
    wPerDict = wAll(CaChain = perCaChain, polyCaChain = perPolyCaChain, perturb_b = 1, segNrRange = segNrRange, dict_b = 0)
    #Initialize array parameters for storing the results (the 'Dict' suffix is 
    #used for easy translation from dictionary-based version):
    I12 = 0 #value of I12 on full chain
    I12Dict = np.zeros((L+1,L+1))
    I1234Dict = np.zeros((L+1,L+1))
    I1234Dict_full = np.zeros((L+1,L+1))
    I1234Dict_full_aid = np.zeros((L+1,L+1))
    I1234Dict_full2 = np.zeros((L+1,L+1))
    I1234Dict_full2_aid = np.zeros((L+1,L+1))
    I1423Dict = np.zeros((L+1,L+1)) 
    I1423Dict_full0 = np.zeros((L+1,L+1)) 
    I1423Dict_full2 = np.zeros((L+1,L+1))
    I1423Dict_full2_aid = np.zeros((L+1,L+1))
    I1324Dict = np.zeros((L+1,L+1))
    I1324Dict_full = np.zeros((L+1,L+1))
    I1324Dict_full_aid = np.zeros((L+1,L+1))
    I1324Dict_full2 = np.zeros((L+1,L+1))
    I1324Dict_full2_aid = np.zeros((L+1,L+1))
    #(12)*:
    I123456Dict = np.zeros((L+1,L+1)) 
    I123645Dict = np.zeros((L+1,L+1))
    I123546Dict = np.zeros((L+1,L+1))
    #(13)*:
    I132456Dict = np.zeros((L+1,L+1))
    I132546Dict = np.zeros((L+1,L+1))
    I132645Dict = np.zeros((L+1,L+1))
    #(14)*
    I142356Dict = np.zeros((L+1,L+1))
    I142536Dict = np.zeros((L+1,L+1))
    I142635Dict = np.zeros((L+1,L+1))
    #(15)*
    I152346Dict = np.zeros((L+1,L+1))
    I152436Dict = np.zeros((L+1,L+1))
    I152634Dict = np.zeros((L+1,L+1))
    #(16)*
    I162345Dict = np.zeros((L+1,L+1))
    I162435Dict = np.zeros((L+1,L+1))
    I162534Dict = np.zeros((L+1,L+1))
    #Main part: computing the measures using recursion formulas
    #one block for each choice of degree; only one block gets executed:    
    if degree == 0:
        if print_b > 0.5:
            print "Computes w-values only."
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    if ((i in segNrRange) or (j in segNrRange)):
    #                print "in change range %d , %d " % (i,j)
                        wVal = wPerDict[i,j]
                    else:
                        wVal = wDict[i,j]
                        wPerDict[i,j] += wVal
    if degree == 2:
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset
                if (j > i and j <=L-1):
                    if ((i in segNrRange) or (j in segNrRange)):
    #                print "in change range %d , %d " % (i,j)
                        wVal = wPerDict[i,j]
                    else:
                        wVal = wDict[i,j]
                        wPerDict[i,j] += wVal
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
    elif degree == 4:
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset
                if (j > i and j <=L-1):
                    if ((i in segNrRange) or (j in segNrRange)):
    #                print "in change range %d , %d " % (i,j)
                        wVal = wPerDict[i,j]
                    else:
                        wVal = wDict[i,j]
                        wPerDict[i,j] += wVal
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
                    #I1234: I1234(i,j;N)
                    I1234Dict[i,j] += I12Dict[j+1,L-1]*wVal + I1234Dict[i,j-1] + I1234Dict[i+1,j] - I1234Dict[i+1,j-1]
                    #I1423: I1423(i,j)
                    I1423Dict[i,j] += I12Dict[i+1,j-1]*wVal +  I1423Dict[i,j-1] + I1423Dict[i+1,j] - I1423Dict[i+1,j-1]
                    #I1324: I1324(i,j;N)
                    I1324Dict[i,j] += (I12Dict[i+1,L-1] - I12Dict[i+1,j] -I12Dict[j,L-1] )*wVal + I1324Dict[i,j-1] + I1324Dict[i+1,j] - I1324Dict[i+1,j-1]
                    if full_b == 1:
                        for b in range(i+1,j-1): #obs: the wPerDict values needed have been computed
                            wVal1 = wPerDict[i,b]
                            I1234Dict_full_aid[i,j] += I12Dict[b+1,j]*wVal1
                            I1324Dict_full_aid[i,j] += I12Dict[b,j]*wVal1
#                            I1234Dict_full[(i,j)] += (I12Dict[(b+1,j)] - I12Dict[(b+1,j-1)])*wDict[(i,b)]
#                            I1324Dict_full[(i,j)] += (-I12Dict[(b,j)] + I12Dict[(b,j-1)])*wDict[(i,b)]
                        #I1234Dict_full add up terms and recursion: I1234(i,j)    
                        I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                        I1234Dict_full[i,j] += I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                        #I1324Dict_full add up terms recursion: I1324(i,j)
                        I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                        I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]                                                   
    elif degree == 6:
        if print_b > 0.5:
            print "Computes degree 2, 4 and 6 measures." 
        for i in lChain[::-1]:  
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    if ((i in segNrRange) or (j in segNrRange)):
    #                print "in change range %d , %d " % (i,j)
                        wVal = wPerDict[i,j]
                    else:
                        wVal = wDict[i,j]
                        wPerDict[i,j] += wVal
                    I12 += wVal
                    I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
                    #degree 4 measures:
                    #I1234: I1234(i,j;N)
                    I1234Dict[i,j] += I12Dict[j+1,L-1]*wVal + I1234Dict[i,j-1] + I1234Dict[i+1,j] - I1234Dict[i+1,j-1]
                    #I1423: I1423(i,j)                    
                    I1423Dict[i,j] += I12Dict[i+1,j-1]*wVal +  I1423Dict[i,j-1] + I1423Dict[i+1,j] - I1423Dict[i+1,j-1]
                    #I1324: I1324(i,j;N)
                    I1324Dict[i,j] += (I12Dict[i+1,L-1] - I12Dict[i+1,j] -I12Dict[j,L-1] )*wVal + I1324Dict[i,j-1] + I1324Dict[i+1,j] - I1324Dict[i+1,j-1]
                    for b in range(i+1,j-1): #obs: the wPerDict values needed have been computed
                        wVal1 = wPerDict[i,b]
                        I1234Dict_full_aid[i,j] += I12Dict[b+1,j]*wVal1
                        I1324Dict_full_aid[i,j] += I12Dict[b,j]*wVal1
                    #I1234Dict_full add up terms and recursion: I1234(i,j))    
                    I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                    I1234Dict_full[i,j] +=  I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                    #I1324Dict_full add up terms recursion: I1324(i,j)
                    I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                    I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]                            
                    #I1423_full: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet
                    I1423Dict_full0[i+1,j] += (I12Dict[i+2,j-1] - I12Dict[i+2,j])*(I12Dict[i+2,L-1] - I12Dict[i+2,j] - I12Dict[i+1,L-1] + I12Dict[i+1,j])
                    I1423Dict_full0[i+1,j] += I1423Dict_full0[i+2,j] + I1423Dict_full0[i+1,j-1] - I1423Dict_full0[i+2,j-1]
                    #to compute certain degree 6 measures we use two auxillary degree 4 
                    #measures; the recursion demands to sum j in the - direction (ie 
                    #from above and down):
                    k = L-1 + i -j+1
                    for c in range(k+1,L): #Obs: the wDict values needed have already been computed (since k,c >i)
                        wVal2 = wPerDict[k,c]
                        I1324Dict_full2_aid[i+1,k] += wVal2*(I12Dict[i+1,c-1] - I12Dict[i+1,k] - I12Dict[k,c-1])
                        I1423Dict_full2_aid[i+1,k] += wVal2*(I12Dict[i+1,L-1] - I12Dict[i+1,c] - I12Dict[k,L-1] + I12Dict[k,c]) 
                    #I1324Dict_full2, recursion: I1324(i;j,N)
                    I1324Dict_full2[i+1,k] += I1324Dict_full2_aid[i+1,k] + I1324Dict_full2[i+1,k+1]
                    #I1423Dict_full2, recursion: I1423(i;j,N)
                    I1423Dict_full2[i+1,k] += I1423Dict_full2_aid[i+1,k] + I1423Dict_full2[i+1,k+1] 
                    #degree 6 measures:
                    #(12)*:
                    I123456Dict[i,j] += I1234Dict[j+1,L-1]*wVal + I123456Dict[i+1,j] + I123456Dict[i,j-1] - I123456Dict[i+1,j-1] 
                    I123645Dict[i,j] += I1423Dict[j+1,L-1]*wVal + I123645Dict[i+1,j] + I123645Dict[i,j-1] - I123645Dict[i+1,j-1] 
                    I123546Dict[i,j] += I1324Dict[j+1,L-1]*wVal+ I123546Dict[i+1,j] + I123546Dict[i,j-1] - I123546Dict[i+1,j-1] 
                    #(13)*:
                    I132456Dict[i,j] += (I1234Dict[i+1,L-1] - I1234Dict[i+1,j] -I1234Dict[j,L-1])*wVal + I132456Dict[i+1,j] + I132456Dict[i,j-1] - I132456Dict[i+1,j-1]
                    I132546Dict[i+1,j] += (I1324Dict_full2[i+2,j+1] - I1324Dict_full[j,L-1])*wPerDict[i+1,j] + I132546Dict[i+2,j] + I132546Dict[i+1,j-1] - I132546Dict[i+2,j-1]
                    I132645Dict[i+1,j] += (I1423Dict_full2[i+2,j+1] - I1423Dict[j,L-1])*wPerDict[i+1,j] + I132645Dict[i+2,j] + I132645Dict[i+1,j-1] - I132645Dict[i+2,j-1]
                    #(14)*
                    #allowing direct computation:
                    I142356Dict[i,j] += I12Dict[i+1,j-1]*I12Dict[j+1,L-1]*wVal + I142356Dict[i+1,j] + I142356Dict[i,j-1] - I142356Dict[i+1,j-1]
                    #back to the more complicated parts, now I142536:
                    #three parts recognized in decompsition as already known:
                    I142536Dict[i+1,j] += wPerDict[i+1,j]*(I1324Dict_full[i+2,L-1] - I1324Dict[i+2,j] - I1324Dict_full2[i+2,j])
                    #recursion part:
                    I142536Dict[i+1,j] += I142536Dict[i+2,j] + I142536Dict[i+1,j-1] - I142536Dict[i+2,j-1] 
                    #Next complicated part, I142635:
                    #three parts recognized in decompsition as already known:
                    I142635Dict[i+1,j] += wPerDict[i+1,j]*(I1423Dict[i+2,L-1] - I1423Dict_full2[i+2,j] - I1423Dict_full0[i+2,j])
                    #recursion part:
                    I142635Dict[i+1,j] += I142635Dict[i+2,j] + I142635Dict[i+1,j-1] - I142635Dict[i+2,j-1]
                    #(15)*
                    I152346Dict[i,j] += wVal*(I1234Dict[i+1,j-1] - I1234Dict_full[i+1,j] - I12Dict[j,L-1]*I12Dict[i+1,j-1])
                    I152346Dict[i,j] += I152346Dict[i+1,j] + I152346Dict[i,j-1] - I152346Dict[i+1,j-1]
                    I152436Dict[i,j] += wVal*(I1324Dict[i+1,j-1] - I1324Dict_full[i+1,j])
                    I152436Dict[i,j] += I152436Dict[i+1,j] + I152436Dict[i,j-1] - I152436Dict[i+1,j-1]
                    I152634Dict[i,j] += wVal*(I1423Dict_full0[i+1,j-1] - I1423Dict[i+1,j])
                    I152634Dict[i,j] += I152634Dict[i+1,j] + I152634Dict[i,j-1] - I152634Dict[i+1,j-1]
                    #(16)*
                    I162345Dict[i,j] += wVal*I1234Dict_full[i+1,j-1] #ER FULL OK?
                    I162345Dict[i,j] += I162345Dict[i+1,j] + I162345Dict[i,j-1] - I162345Dict[i+1,j-1]
                    I162435Dict[i,j] += wVal*I1324Dict_full[i+1,j-1] #ER FULL OK?
                    I162435Dict[i,j] += I162435Dict[i+1,j] + I162435Dict[i,j-1] - I162435Dict[i+1,j-1]
                    I162534Dict[i,j] += wVal*I1423Dict[i+1,j-1]
                    I162534Dict[i,j] += I162534Dict[i+1,j] + I162534Dict[i,j-1] - I162534Dict[i+1,j-1]                    
        #the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:
        for j in lChain[::-1]:
            i = 0
            for c in range(j+1,L): #Obs: the 2d-loop above has sweeped the whole simplex so the wDictPer has been computed already 
                wVal = wPerDict[j,c]
#                        I1234Dict_full2_aid[(i,j)] +=  I12Dict[(i,c-1)]*wDict[(j,c)]
                I1324Dict_full2_aid[i,j] += wVal*(I12Dict[i,c-1] - I12Dict[i,j] - I12Dict[j,c-1])
                I1423Dict_full2_aid[i,j] += wVal*(I12Dict[i,L-1] - I12Dict[i,c] - I12Dict[j,L-1] + I12Dict[j,c]) 
           #recursion:
            I1324Dict_full2[i,j] += I1324Dict_full2_aid[i,j] + I1324Dict_full2[i,j+1]
            I1423Dict_full2[i,j] += I1423Dict_full2_aid[i,j] + I1423Dict_full2[i,j+1] 
        #Missing terms for some of the degree 6 measures and I1423_full0:
        for j in lChain[1:]:
            i=0
            wVal = wPerDict[i,j] #Obs: the 2d-loop above has sweeped the whole simplex so the wDictPer has been computed already 
            I1423Dict_full0[i,j] += (I12Dict[i+1,j-1] - I12Dict[i+1,j])*(I12Dict[i+1,L-1] - I12Dict[i+1,j] - I12Dict[i,L-1] + I12Dict[i,j])
            I1423Dict_full0[i,j] += I1423Dict_full0[i+1,j] + I1423Dict_full0[i,j-1] - I1423Dict_full0[i+1,j-1]
            #
            I132546Dict[i,j] += (I1324Dict_full2[i+1,j+1] - I1324Dict_full[j,L-1])*wVal + I132546Dict[i+1,j] + I132546Dict[i,j-1] - I132546Dict[i+1,j-1]
            I132645Dict[i,j] += (I1423Dict_full2[i+1,j+1] - I1423Dict[j,L-1])*wVal + I132645Dict[i+1,j] + I132645Dict[i,j-1] - I132645Dict[i+1,j-1]
            #
            I142536Dict[i,j] += wVal*(I1324Dict_full[i+1,L-1] - I1324Dict[i+1,j] - I1324Dict_full2[i+1,j])
            I142536Dict[i,j] += I142536Dict[i+1,j] + I142536Dict[i,j-1] - I142536Dict[i+1,j-1] 
            #
            I142635Dict[i,j] += wVal*(I1423Dict[i+1,L-1] - I1423Dict_full2[i+1,j] - I1423Dict_full0[i+1,j])
            I142635Dict[i,j] += I142635Dict[i+1,j] + I142635Dict[i,j-1] - I142635Dict[i+1,j-1]    
    else: 
        print "Degree must be 2, 4 or 6."                         
    outDict = np.array((wPerDict, I12Dict, 
                        I1234Dict, I1234Dict_full, 
                        I1324Dict, I1324Dict_full, I1324Dict_full2,
                        I1423Dict, I1423Dict_full0,I1423Dict_full2,
                        I123456Dict, I123546Dict, I123645Dict,
                        I132456Dict, I132546Dict, I132645Dict,
                        I142356Dict, I142536Dict, I142635Dict,
                        I152346Dict, I152436Dict, I152634Dict,
                        I162345Dict, I162435Dict, I162534Dict)) 
    return outDict #, perCaChain, perPolyCaChain #, wPerDict, I12Dict, I1234Dict, I1423Dict, I1324Dict,I1234Dict_full, I1324Dict_full




#def Iperturbed_dict_vect(CaChain, 
#                polyCaChain, 
#                wDict,
#                residueNrRange = range(0,5),
#                perturbation = [(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3)], 
#                degree = 2,
#                full_b = 1,
#                print_b = 1):
#    '''Vectorized version of Iperturbed_dict.
#    Computes the change in the measures of the given polychain due to 
#    the perturbation. Input: 
#    CaChain: C-alpha chain
#    polyCaChain: chain of segment corresponding to CaChain
#    wDict: dictionary of w-values as output from I_dict
#    residueNrRange: indices of the residues in the CaChain (C-alphas) that should be perturbed
#    perturbation: the perturbation given in 3d-coordinates. Must be of lenght equal to residueNrRange
#    Other: as I_dict.
#    Returns: same as I_dict, but for the perturbed chain.'''
#    #compute the polygonal chain of the perturbed chain:
#    perCaChain, perPolyCaChain, segNrRange = perturbPolygonalChain(CaChain = CaChain, 
#                                                    polyCaChain = polyCaChain,
#                     residueNrRange = residueNrRange, 
#                     perturbation = perturbation)



def aggr_pert_vect_dict(wDict, perCaChain, perPolyCaChain, segNrRange, degree, full_b, print_b):                     
    '''Computation part for I_perturbed, vectorized.'''   
    if print_b > 0.5:
        print segNrRange
    #Get length of chain:
    L = len(perPolyCaChain)
    if print_b > 0.5:
        print "Length of chain is: %d" % L
    lChain = range(L)
    #call writhe computation in one go:
    wPerDict = wAll(CaChain = perCaChain, polyCaChain = perPolyCaChain, perturb_b = 1, segNrRange = segNrRange)
    #init. dictionaries for storing results:
    I12 = 0
    I12Dict = {}
    I1234Dict = {}
    I1234Dict_full = {}
    I1234Dict_full_aid = {}
    I1234Dict_full2 = {}
    I1234Dict_full2_aid = {}
    I1423Dict = {} 
    I1423Dict_full0 = {} 
    I1423Dict_full2 = {}
    I1423Dict_full2_aid = {}
    I1324Dict = {}
    I1324Dict_full = {}
    I1324Dict_full_aid = {}
    I1324Dict_full2 = {}
    I1324Dict_full2_aid = {}
    #(12)*:
    I123456Dict = {} 
    I123645Dict = {}
    I123546Dict = {}
    #(13)*:
    I132456Dict = {}
    I132546Dict = {}
    I132645Dict = {}
    #(14)*
    I142356Dict = {}
    I142536Dict = {}
    I142635Dict = {}
    #(15)*
    I152346Dict = {}
    I152436Dict = {}
    I152634Dict = {}
    #(16)*
    I162345Dict = {}
    I162435Dict = {}
    I162534Dict = {}
    #init:
    for i in range(L+1):
        for j in range(L+1):
            I12Dict[(i,j)] = 0
#            I12Dict[(-1,j)] = 0
            if degree > 2: 
                I1234Dict[(i,j)] = 0
                I1234Dict_full[(i,j)] = 0
                I1324Dict[(i,j)] = 0  
                I1423Dict[(i,j)] = 0
#                if full_b ==1:
                I1234Dict_full[(i,j)] = 0
                I1234Dict_full_aid[(i,j)] = 0
                I1234Dict_full2[(i,j)] = 0
                I1234Dict_full2_aid[(i,j)] = 0
                I1324Dict_full[(i,j)] = 0
                I1324Dict_full_aid[(i,j)] = 0
                I1324Dict_full2[(i,j)] = 0
#                I1324Dict_full2[(-1,j)] = 0
                I1324Dict_full2_aid[(i,j)] = 0
#                I1324Dict_full2_aid[(-1,j)] = 0
                I1423Dict_full0[(i,j)] = 0
                I1423Dict_full2[(i,j)] = 0
                I1423Dict_full2_aid[(i,j)] = 0
            if degree > 4:
                I123456Dict[(i,j)] = 0 
                I123645Dict[(i,j)] = 0
                I123546Dict[(i,j)] = 0
                #(13)*:
                I132456Dict[(i,j)] = 0
                I132546Dict[(i,j)] = 0
                I132645Dict[(i,j)] = 0
                #(14)*
                I142356Dict[(i,j)] = 0
                I142536Dict[(i,j)] = 0
                I142635Dict[(i,j)] = 0
                #(15)*
                I152346Dict[(i,j)] = 0
                I152436Dict[(i,j)] = 0
                I152634Dict[(i,j)] = 0
                #(16)*
                I162345Dict[(i,j)] = 0
                I162435Dict[(i,j)] = 0
                I162534Dict[(i,j)] = 0
    #Main part: computing the measures using recursion formulas
    #one block for each choice of degree; only one block gets executed:   
    if degree == 2:
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset
                if (j > i and j <=L-1):
                    if ((i in segNrRange) or (j in segNrRange)):
    #                print "in change range %d , %d " % (i,j)
                        wVal = wPerDict[(i,j)]
                    else:
                        wVal = wDict[(i,j)]
                        wPerDict[(i,j)] += wVal
                    I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
    elif degree == 4:
        for i in lChain[::-1]:
            for j in lChain:
                wVal = 0 #reset
                if (j > i and j <=L-1):
                    if ((i in segNrRange) or (j in segNrRange)):
    #                print "in change range %d , %d " % (i,j)
                        wVal = wPerDict[(i,j)]
                    else:
                        wVal = wDict[(i,j)]
                        wPerDict[(i,j)] += wVal
                    I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
                    I1234Dict[(i,j)] += I12Dict[(j+1,L-1)]*wVal + I1234Dict[(i,j-1)] + I1234Dict[(i+1,j)] - I1234Dict[(i+1,j-1)]
                    I1423Dict[(i,j)] += I12Dict[(i+1,j-1)]*wVal +  I1423Dict[(i,j-1)] + I1423Dict[(i+1,j)] - I1423Dict[(i+1,j-1)]
                    I1324Dict[(i,j)] += (I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] -I12Dict[(j,L-1)] )*wVal + I1324Dict[(i,j-1)] + I1324Dict[(i+1,j)] - I1324Dict[(i+1,j-1)]
                    if full_b == 1:
                        for b in range(i+1,j-1): #obs: the wPerDict values needed have been computed
                            wVal1 = wPerDict[(i,b)]
                            I1234Dict_full_aid[(i,j)] += I12Dict[(b+1,j)]*wVal1
                            I1324Dict_full_aid[(i,j)] += I12Dict[(b,j)]*wVal1
#                            I1234Dict_full[(i,j)] += (I12Dict[(b+1,j)] - I12Dict[(b+1,j-1)])*wDict[(i,b)]
#                            I1324Dict_full[(i,j)] += (-I12Dict[(b,j)] + I12Dict[(b,j-1)])*wDict[(i,b)]
                        #I1234Dict_full add up terms and recursion: I1234(i,j))    
                        I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                        I1234Dict_full[i,j] +=  I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                        #I1324Dict_full add up terms recursion: I1324(i,j)
                        I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                        I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]            
    elif degree == 6:
        if print_b > 0.5:
            print "Computes degree 2, 4 and 6 measures." 
        for i in lChain[::-1]:  
            for j in lChain:
                wVal = 0 #reset for safety
                if (j > i and j <=L-1):
                    if ((i in segNrRange) or (j in segNrRange)):
    #                print "in change range %d , %d " % (i,j)
                        wVal = wPerDict[(i,j)]
                    else:
                        wVal = wDict[(i,j)]
                        wPerDict[(i,j)] += wVal
                    I12 += wVal
                    I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
                    #degree 4 measures:
                    #I1234: I1234(i,j;N)
                    I1234Dict[(i,j)] += I12Dict[(j+1,L-1)]*wVal + I1234Dict[(i,j-1)] + I1234Dict[(i+1,j)] - I1234Dict[(i+1,j-1)]
                    #I1423: I1423(i,j)                    
                    I1423Dict[(i,j)] += I12Dict[(i+1,j-1)]*wVal +  I1423Dict[(i,j-1)] + I1423Dict[(i+1,j)] - I1423Dict[(i+1,j-1)]
                    #I1324: I1324(i,j;N)
                    I1324Dict[(i,j)] += (I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] -I12Dict[(j,L-1)] )*wVal + I1324Dict[(i,j-1)] + I1324Dict[(i+1,j)] - I1324Dict[(i+1,j-1)]
                    for b in range(i+1,j-1): #obs: the wPerDict values needed have been computed
                        wVal1 = wPerDict[(i,b)]
                        I1234Dict_full_aid[(i,j)] +=  I12Dict[(b+1,j)]*wVal1
                        I1324Dict_full_aid[(i,j)] += I12Dict[(b,j)]*wVal1
                    #I1234Dict_full add up terms and recursion: I1234(i,j) 
                    I1234Dict_full[(i,j)] += I1234Dict_full_aid[(i,j)] -I1234Dict_full_aid[(i,j-1)]
                    I1234Dict_full[(i,j)] +=  I1234Dict_full[(i,j-1)] + I1234Dict_full[(i+1,j)] - I1234Dict_full[(i+1,j-1)]
                    #I1324Dict_full add up terms recursion: I1324(i,j)
                    I1324Dict_full[(i,j)] += I1324Dict_full_aid[(i,j-1)] - I1324Dict_full_aid[(i,j)]   
                    I1324Dict_full[(i,j)] += (I12Dict[(i,j-1)] - I12Dict[(i+1,j-1)])*(I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]) + I1324Dict_full[(i,j-1)] + I1324Dict_full[(i+1,j)] - I1324Dict_full[(i+1,j-1)]                            
                    #I1423_full: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet
                    I1423Dict_full0[(i+1,j)] += (I12Dict[(i+2,j-1)] - I12Dict[(i+2,j)])*(I12Dict[(i+2,L-1)] - I12Dict[(i+2,j)] - I12Dict[(i+1,L-1)] + I12Dict[(i+1,j)])
                    I1423Dict_full0[(i+1,j)] += I1423Dict_full0[(i+2,j)] + I1423Dict_full0[(i+1,j-1)] - I1423Dict_full0[(i+2,j-1)]
                    #to compute certain degree 6 measures we use two auxillary degree 4 measures; the recursion
                    #demands to sum j in the - direction (ie from above and down):
                    k = L-1 + i -j+1
                    for c in range(k+1,L): #obs: the wPerDict values needed have been computed (since k,c >i)
                        wVal2 = wPerDict[(k,c)]
                        I1324Dict_full2_aid[(i+1,k)] += wVal2*(I12Dict[(i+1,c-1)] - I12Dict[(i+1,k)] - I12Dict[(k,c-1)])
                        I1423Dict_full2_aid[(i+1,k)] += wVal2*(I12Dict[(i+1,L-1)] - I12Dict[(i+1,c)] - I12Dict[(k,L-1)] + I12Dict[(k,c)]) 
                    #I1324Dict_full2, recursion: I1324(i;j,N)
                    I1324Dict_full2[(i+1,k)] += I1324Dict_full2_aid[(i+1,k)] + I1324Dict_full2[(i+1,k+1)]
                    #I1423Dict_full2, recursion: I1423(i;j,N)
                    I1423Dict_full2[(i+1,k)] += I1423Dict_full2_aid[(i+1,k)] + I1423Dict_full2[(i+1,k+1)] 
                    #degree 6 measures:
                    #(12)*:
                    I123456Dict[(i,j)] += I1234Dict[(j+1,L-1)]*wVal + I123456Dict[(i+1,j)] + I123456Dict[(i,j-1)] - I123456Dict[(i+1,j-1)] 
                    I123645Dict[(i,j )] += I1423Dict[(j+1,L-1)]*wVal + I123645Dict[(i+1,j)] + I123645Dict[(i,j-1)] - I123645Dict[(i+1,j-1)] 
                    I123546Dict[(i,j)] += I1324Dict[(j+1,L-1)]*wVal+ I123546Dict[(i+1,j)] + I123546Dict[(i,j-1)] - I123546Dict[(i+1,j-1)] 
                    #(13)*:
                    I132456Dict[(i,j)] += (I1234Dict[(i+1,L-1)] - I1234Dict[(i+1,j)] -I1234Dict[(j,L-1)])*wVal + I132456Dict[(i+1,j)] + I132456Dict[(i,j-1)] - I132456Dict[(i+1,j-1)]
                    I132546Dict[(i+1,j)] += (I1324Dict_full2[(i+2,j+1)] - I1324Dict_full[(j,L-1)])*wPerDict[(i+1,j)] + I132546Dict[(i+2,j)] + I132546Dict[(i+1,j-1)] - I132546Dict[(i+2,j-1)]
                    I132645Dict[(i+1,j)] += (I1423Dict_full2[(i+2,j+1)] - I1423Dict[(j,L-1)])*wPerDict[(i+1,j)] + I132645Dict[(i+2,j)] + I132645Dict[(i+1,j-1)] - I132645Dict[(i+2,j-1)]
                    #(14)*
                    #allowing direct computation:
                    I142356Dict[(i,j)] += I12Dict[(i+1,j-1)]*I12Dict[(j+1,L-1)]*wVal + I142356Dict[(i+1,j)] + I142356Dict[(i,j-1)] - I142356Dict[(i+1,j-1)]
                    #back to the more complicated parts, now I142536:
                    #three parts recognized in decompsition as already known:
                    I142536Dict[(i+1,j)] += wPerDict[(i+1,j)]*(I1324Dict_full[(i+2,L-1)] - I1324Dict[(i+2,j)] - I1324Dict_full2[(i+2,j)])
                    #recursion part:
                    I142536Dict[(i+1,j)] += I142536Dict[(i+2,j)] + I142536Dict[(i+1,j-1)] - I142536Dict[(i+2,j-1)] 
                    #Next complicated part, I142635:
                    #three parts recognized in decompsition as already known:
                    I142635Dict[(i+1,j)] += wPerDict[(i+1,j)]*(I1423Dict[(i+2,L-1)] - I1423Dict_full2[(i+2,j)] - I1423Dict_full0[(i+2,j)])
                    #recursion part:
                    I142635Dict[(i+1,j)] += I142635Dict[(i+2,j)] + I142635Dict[(i+1,j-1)] - I142635Dict[(i+2,j-1)]
                    #(15)*
                    I152346Dict[(i,j)] += wVal*(I1234Dict[(i+1,j-1)] - I1234Dict_full[(i+1,j)] - I12Dict[(j,L-1)]*I12Dict[(i+1,j-1)])
                    I152346Dict[(i,j)] += I152346Dict[(i+1,j)] + I152346Dict[(i,j-1)] - I152346Dict[(i+1,j-1)]
                    I152436Dict[(i,j)] += wVal*(I1324Dict[(i+1,j-1)] - I1324Dict_full[(i+1,j)])
                    I152436Dict[(i,j)] += I152436Dict[(i+1,j)] + I152436Dict[(i,j-1)] - I152436Dict[(i+1,j-1)]
                    I152634Dict[(i,j)] += wVal*(I1423Dict_full0[(i+1,j-1)] - I1423Dict[(i+1,j)])
                    I152634Dict[(i,j)] += I152634Dict[(i+1,j)] + I152634Dict[(i,j-1)] - I152634Dict[(i+1,j-1)]
                    #(16)*
                    I162345Dict[(i,j)] += wVal*I1234Dict_full[(i+1,j-1)] #ER FULL OK?
                    I162345Dict[(i,j)] += I162345Dict[(i+1,j)] + I162345Dict[(i,j-1)] - I162345Dict[(i+1,j-1)]
                    I162435Dict[(i,j)] += wVal*I1324Dict_full[(i+1,j-1)] #ER FULL OK?
                    I162435Dict[(i,j)] += I162435Dict[(i+1,j)] + I162435Dict[(i,j-1)] - I162435Dict[(i+1,j-1)]
                    I162534Dict[(i,j)] += wVal*I1423Dict[(i+1,j-1)]
                    I162534Dict[(i,j)] += I162534Dict[(i+1,j)] + I162534Dict[(i,j-1)] - I162534Dict[(i+1,j-1)]                    
        #the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:
        for j in lChain[::-1]:
            i = 0
            for c in range(j+1,L):
                wVal = wPerDict[(j,c)] #Obs: the 2d-loop above has sweeped the whole simplex so the wDictPer has been computed already 
#                        I1234Dict_full2_aid[(i,j)] +=  I12Dict[(i,c-1)]*wDict[(j,c)]
                I1324Dict_full2_aid[(i,j)] += wVal*(I12Dict[(i,c-1)] - I12Dict[(i,j)] - I12Dict[(j,c-1)])
                I1423Dict_full2_aid[(i,j)] += wVal*(I12Dict[(i,L-1)] - I12Dict[(i,c)] - I12Dict[(j,L-1)] + I12Dict[(j,c)]) 
           #recursion:
            I1324Dict_full2[(i,j)] += I1324Dict_full2_aid[(i,j)] + I1324Dict_full2[(i,j+1)]
            I1423Dict_full2[(i,j)] += I1423Dict_full2_aid[(i,j)] + I1423Dict_full2[(i,j+1)] 
        #Missing terms for some of the degree 6 measures and I1423_full0:
        for j in lChain[1:]:
            i=0
            wVal = wPerDict[(i,j)] #Obs: the 2d-loop above has sweeped the whole simplex so the wDictPer has been computed already 
            I1423Dict_full0[(i,j)] += (I12Dict[(i+1,j-1)] - I12Dict[(i+1,j)])*(I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] - I12Dict[(i,L-1)] + I12Dict[(i,j)])
            I1423Dict_full0[(i,j)] += I1423Dict_full0[(i+1,j)] + I1423Dict_full0[(i,j-1)] - I1423Dict_full0[(i+1,j-1)]
            #
            I132546Dict[(i,j)] += (I1324Dict_full2[(i+1,j+1)] - I1324Dict_full[(j,L-1)])*wVal + I132546Dict[(i+1,j)] + I132546Dict[(i,j-1)] - I132546Dict[(i+1,j-1)]
            I132645Dict[(i,j)] += (I1423Dict_full2[(i+1,j+1)] - I1423Dict[(j,L-1)])*wVal + I132645Dict[(i+1,j)] + I132645Dict[(i,j-1)] - I132645Dict[(i+1,j-1)]
            #
            I142536Dict[(i,j)] += wVal*(I1324Dict_full[(i+1,L-1)] - I1324Dict[(i+1,j)] - I1324Dict_full2[(i+1,j)])
            I142536Dict[(i,j)] += I142536Dict[(i+1,j)] + I142536Dict[(i,j-1)] - I142536Dict[(i+1,j-1)] 
            #
            I142635Dict[(i,j)] += wVal*(I1423Dict[(i+1,L-1)] - I1423Dict_full2[(i+1,j)] - I1423Dict_full0[(i+1,j)])
            I142635Dict[(i,j)] += I142635Dict[(i+1,j)] + I142635Dict[(i,j-1)] - I142635Dict[(i+1,j-1)]    
    else: 
        print "Degree must be 2, 4 or 6."
    outDict = {}
    outDict['w'] = wPerDict
    outDict['I12'] = I12Dict
    outDict['I1234'] = I1234Dict
    outDict['I1234_full'] = I1234Dict_full
    outDict['I1324'] = I1324Dict
    outDict['I1324_full'] = I1324Dict_full
    outDict['I1324_full2'] = I1324Dict_full2    
    outDict['I1423'] = I1423Dict
    outDict['I1423_full0'] = I1423Dict_full0
    outDict['I1423_full2'] = I1423Dict_full2    

    outDict['I123456'] = I123456Dict
    outDict['I123546'] = I123546Dict
    outDict['I123645'] = I123645Dict

    outDict['I132456'] = I132456Dict
    outDict['I132546'] = I132546Dict
    outDict['I132645'] = I132645Dict

    outDict['I142356'] = I142356Dict
    outDict['I142536'] = I142536Dict
    outDict['I142635'] = I142635Dict
    
    outDict['I152346'] = I152346Dict
    outDict['I152436'] = I152436Dict
    outDict['I152634'] = I152634Dict

    outDict['I162345'] = I162345Dict
    outDict['I162435'] = I162435Dict
    outDict['I162534'] = I162534Dict
    return outDict, perCaChain, perPolyCaChain #, wPerDict, I12Dict, I1234Dict, I1423Dict, I1324Dict,I1234Dict_full, I1324Dict_full


####################################################################
#Versions for running a set of perturbation in one go (all vectorized):
####################################################################

#@jit(float64[:,:,:](float64[:,:,:], int64, int64, int64), target = 'cpu')
def wAll_multi(chainsInput, lengthPolyChain, print_b = 0, dict_b = 0):
    '''Version of wAll to cope with the simultaneous computation of 
    w-terms for a set of chains (aimed at computing in one call the 
    measures on a set of perturbations of the same chain.
    Input: 
    dict_b: boolean controlling if the code is array based (input:
    chainsInput as array and output as array) or dictionary based.
    Returns: array (of arrays) or dictionary (of dictionaries) 
    containing the desired w-values.'''
    L = lengthPolyChain
    lChain = range(lengthPolyChain)
    if print_b > 0.5:
        print lChain
    if dict_b > 0.5:
        #init:
        wDict_multi = {}
        #first derive the full set of all pairs of segments on which the writhe mus tbe computed:   
        s101array = []
        s102array = []
        s103array = []
        s111array = []
        s112array = []
        s113array = []
        s201array = []
        s202array = []
        s203array = []
        s211array = []
        s212array = []
        s213array = []
        for k in chainsInput.keys():
            #init
            wDict_multi[k] = {}
            for i in range(L+1):
                for j in range(L+1):
                    wDict_multi[k][(i,j)] = 0
            segNrRange = chainsInput[k][2]
            lastSeg = segNrRange[::-1][0]
            firstSeg = segNrRange[0]
            pChain = chainsInput[k][1]
            for i in range(firstSeg):
                for j in segNrRange:
                    s101array.append(pChain[i][0][0])
                    s102array.append(pChain[i][0][1])
                    s103array.append(pChain[i][0][2])
                    s111array.append(pChain[i][1][0])
                    s112array.append(pChain[i][1][1])
                    s113array.append(pChain[i][1][2])
                    s201array.append(pChain[j][0][0])
                    s202array.append(pChain[j][0][1])
                    s203array.append(pChain[j][0][2])
                    s211array.append(pChain[j][1][0])
                    s212array.append(pChain[j][1][1])
                    s213array.append(pChain[j][1][2])
            for i in segNrRange:
                for j in range(i, L):
                    s101array.append(pChain[i][0][0])
                    s102array.append(pChain[i][0][1])
                    s103array.append(pChain[i][0][2])
                    s111array.append(pChain[i][1][0])
                    s112array.append(pChain[i][1][1])
                    s113array.append(pChain[i][1][2])
                    s201array.append(pChain[j][0][0])
                    s202array.append(pChain[j][0][1])
                    s203array.append(pChain[j][0][2])
                    s211array.append(pChain[j][1][0])
                    s212array.append(pChain[j][1][1])
                    s213array.append(pChain[j][1][2])
        s101array = np.array(s101array)
        s102array = np.array(s102array)
        s103array = np.array(s103array)
        s111array = np.array(s111array)
        s112array = np.array(s112array)
        s113array = np.array(s113array)
        s201array = np.array(s201array)
        s202array = np.array(s202array)
        s203array = np.array(s203array)
        s211array = np.array(s211array)
        s212array = np.array(s212array)
        s213array = np.array(s213array)
        wAll =0
        #compute the write terms:
        wAll = wVect_2(s101array,s102array,s103array,s111array,s112array,s113array,s201array,s202array,s203array, s211array,s212array,s213array)
        cnt = 0
        #write the seq of results into a dictionary
        for k in chainsInput.keys():
            segNrRange = chainsInput[k][2]
            for i in lChain:
                for j in lChain:
                    if (j > i and ((i in segNrRange) or (j in segNrRange))):
                        wDict_multi[k][(i,j)] += wAll[cnt]
                        cnt +=1
                    else:
                        wDict_multi[k][(i,j)] = 0
    else: #dict_b != 1. 
    #Here follows essentially the same piece of code as for 
    #dict_b = 1, only is the output structured as an array                     
        #init:
        nrPer = len(chainsInput)
        #init
        wDict_multi = np.zeros((nrPer, L+1,L+1))
        #first derive the full set of all pairs of segments on which the writhe mus tbe computed:   
        s101array = []
        s102array = []
        s103array = []
        s111array = []
        s112array = []
        s113array = []
        s201array = []
        s202array = []
        s203array = []
        s211array = []
        s212array = []
        s213array = []
        for k in range(len(chainsInput)):
            segNrRange = chainsInput[k][2]
            lastSeg = segNrRange[::-1][0]
            firstSeg = segNrRange[0]
            pChain = chainsInput[k][1]
            for i in range(firstSeg):
                for j in segNrRange:
                    s101array.append(pChain[i][0][0])
                    s102array.append(pChain[i][0][1])
                    s103array.append(pChain[i][0][2])
                    s111array.append(pChain[i][1][0])
                    s112array.append(pChain[i][1][1])
                    s113array.append(pChain[i][1][2])
                    s201array.append(pChain[j][0][0])
                    s202array.append(pChain[j][0][1])
                    s203array.append(pChain[j][0][2])
                    s211array.append(pChain[j][1][0])
                    s212array.append(pChain[j][1][1])
                    s213array.append(pChain[j][1][2])
            for i in segNrRange:
                for j in range(i, L):
                    s101array.append(pChain[i][0][0])
                    s102array.append(pChain[i][0][1])
                    s103array.append(pChain[i][0][2])
                    s111array.append(pChain[i][1][0])
                    s112array.append(pChain[i][1][1])
                    s113array.append(pChain[i][1][2])
                    s201array.append(pChain[j][0][0])
                    s202array.append(pChain[j][0][1])
                    s203array.append(pChain[j][0][2])
                    s211array.append(pChain[j][1][0])
                    s212array.append(pChain[j][1][1])
                    s213array.append(pChain[j][1][2])
        s101array = np.array(s101array)
        s102array = np.array(s102array)
        s103array = np.array(s103array)
        s111array = np.array(s111array)
        s112array = np.array(s112array)
        s113array = np.array(s113array)
        s201array = np.array(s201array)
        s202array = np.array(s202array)
        s203array = np.array(s203array)
        s211array = np.array(s211array)
        s212array = np.array(s212array)
        s213array = np.array(s213array)
        wAll =0
        #compute the write terms:
        wAll = wVect_2(s101array,s102array,s103array,s111array,s112array,s113array,s201array,s202array,s203array, s211array,s212array,s213array)
        cnt = 0
        #write the seq of results into a dictionary
        for k in range(nrPer):
            segNrRange = chainsInput[k][2]
            for i in lChain:
                for j in lChain:
                    if (j > i and ((i in segNrRange) or (j in segNrRange))):
                        wDict_multi[k][i,j] += wAll[cnt]
                        cnt +=1
                    else:
                        wDict_multi[k][i,j] = 0                                    
    return wDict_multi

#@jit(float64[:,:,:](float64[:,:,:], float64[:,:,:], float64[:,:,:],int64, int64[:],int64, int64, int64, int64), target = 'cpu')
def Iperturbed_vect_multi(CaChain, 
               polyCaChain, 
               wDict,
               nrPerturbations = 1,
               residueNrRange = range(0,5),
               perturbationLength = 5, 
                degree = 4,
                full_b = 0,
                print_b = 0):
    '''Multi perturbation version of Iperturbed_vect.'''
    #for holding output:
    tot_outDict = {}
    #Generate a set of perturbations
    perturbedChains = {}
    perturbations = {}
    for i in range(nrPerturbations):
        #pick randomly the perturbations:
        perturbation = np.random.rand(perturbationLength,3)
        #compute the polygonal chain of the perturbed chain:
        perCaChain, perPolyCaChain, segNrRange = perturbPolygonalChain(CaChain = CaChain, 
                                                        polyCaChain = polyCaChain,
                         residueNrRange = residueNrRange, 
                         perturbation = perturbation)
        for chain in CaChain.keys():   
            if not(perturbedChains.has_key(chain)):
                perturbedChains[chain] = []
            perturbedChains[chain].append([perCaChain[chain], perPolyCaChain[chain], segNrRange[chain]])
            #the actual perturbations are not used beow, but could be useful to know, so we record them:
            if not(perturbations.has_key(chain)):
                perturbations[chain] = []
            perturbations[chain].append(perturbation)
    #compute invariants:
    for chain in CaChain.keys():
        perturbedChainsArr = np.array(perturbedChains[chain])
        #Get length of chain:
        L = len(polyCaChain[chain])
        if print_b > 0.5:
            print "Length of chain is: %d" % L
        #compute the writhe terms in one go (call array based version by setting dict_b = 0):
        wPerDict_multi = wAll_multi(chainsInput = perturbedChainsArr, lengthPolyChain = L, print_b = print_b, dict_b = 0)
        #compute the measures, running through the perturbations:
        lChain = range(L)       
        #init. array for storing results (all perturbations):
        tot_outDict[chain] = np.array([]) #np.zeros(nrPerturbations) 
        #loop through perturbations and compute measures:
        for p in range(len(perturbedChains[chain])): #p for perturbation
            segNrRange[chain] = perturbedChains[chain][p][2]
            if print_b > 0.5:
                print segNrRange[chain]
            firstSeg = segNrRange[chain][0]
            lastSeg = segNrRange[chain][::-1][0]
    #        if print_b > 0.5:
    #            print segNrRange    
            #Initialize/reset array parameters for storing the results (the 'Dict' suffix is 
            #used for easy translation from dictionary-based version):
            I12 = 0 #value of I12 on full chain
            I12Dict = np.zeros((L+1,L+1))
            I1234Dict = np.zeros((L+1,L+1))
            I1234Dict_full = np.zeros((L+1,L+1))
            I1234Dict_full_aid = np.zeros((L+1,L+1))
            I1234Dict_full2 = np.zeros((L+1,L+1))
            I1234Dict_full2_aid = np.zeros((L+1,L+1))
            I1423Dict = np.zeros((L+1,L+1)) 
            I1423Dict_full0 = np.zeros((L+1,L+1)) 
            I1423Dict_full2 = np.zeros((L+1,L+1))
            I1423Dict_full2_aid = np.zeros((L+1,L+1))
            I1324Dict = np.zeros((L+1,L+1))
            I1324Dict_full = np.zeros((L+1,L+1))
            I1324Dict_full_aid = np.zeros((L+1,L+1))
            I1324Dict_full2 = np.zeros((L+1,L+1))
            I1324Dict_full2_aid = np.zeros((L+1,L+1))
            #(12)*:
            I123456Dict = np.zeros((L+1,L+1)) 
            I123645Dict = np.zeros((L+1,L+1))
            I123546Dict = np.zeros((L+1,L+1))
            #(13)*:
            I132456Dict = np.zeros((L+1,L+1))
            I132546Dict = np.zeros((L+1,L+1))
            I132645Dict = np.zeros((L+1,L+1))
            #(14)*
            I142356Dict = np.zeros((L+1,L+1))
            I142536Dict = np.zeros((L+1,L+1))
            I142635Dict = np.zeros((L+1,L+1))
            #(15)*
            I152346Dict = np.zeros((L+1,L+1))
            I152436Dict = np.zeros((L+1,L+1))
            I152634Dict = np.zeros((L+1,L+1))
            #(16)*
            I162345Dict = np.zeros((L+1,L+1))
            I162435Dict = np.zeros((L+1,L+1))
            I162534Dict = np.zeros((L+1,L+1))
            #Main part: computing the measures using recursion formulas
            #one block for each choice of degree; only one block gets executed:    
            wPerDict = wPerDict_multi[p]
            if degree == 0:
                if print_b > 0.5:
                    print "Computes w-values only."
                for i in lChain[::-1]:
                    for j in lChain:
                        wVal = 0 #reset for safety
                        if (j > i and j <=L-1):
                            if ((i in segNrRange[chain]) or (j in segNrRange[chain])):
            #                print "in change range %d , %d " % (i,j)
                                wVal = wPerDict[i,j]
                            else:
                                wVal = wDict[chain][i,j]
                                wPerDict[i,j] += wVal
            if degree == 2:
                for i in lChain[::-1]:
                    for j in lChain:
                        wVal = 0 #reset
                        if (j > i and j <=L-1):
                            if ((i in segNrRange[chain]) or (j in segNrRange[chain])):
            #                print "in change range %d , %d " % (i,j)
                                wVal = wPerDict[i,j]
                            else:
                                wVal = wDict[chain][i,j]
                                wPerDict[i,j] += wVal
                            I12 += wVal
                            I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
            elif degree == 4:
                for i in lChain[::-1]:
                    for j in lChain:
                        wVal = 0 #reset
                        if (j > i and j <=L-1):
                            if ((i in segNrRange[chain]) or (j in segNrRange[chain])):
            #                print "in change range %d , %d " % (i,j)
                                wVal = wPerDict[i,j]
                            else:
                                wVal = wDict[chain][i,j]
                                wPerDict[i,j] += wVal
                            I12 += wVal
                            I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
                            #I1234: I1234(i,j;N)
                            I1234Dict[i,j] += I12Dict[j+1,L-1]*wVal + I1234Dict[i,j-1] + I1234Dict[i+1,j] - I1234Dict[i+1,j-1]
                            #I1423: I1423(i,j)
                            I1423Dict[i,j] += I12Dict[i+1,j-1]*wVal +  I1423Dict[i,j-1] + I1423Dict[i+1,j] - I1423Dict[i+1,j-1]
                            #I1324: I1324(i,j;N)
                            I1324Dict[i,j] += (I12Dict[i+1,L-1] - I12Dict[i+1,j] -I12Dict[j,L-1] )*wVal + I1324Dict[i,j-1] + I1324Dict[i+1,j] - I1324Dict[i+1,j-1]
                            if full_b == 1:
                                for b in range(i+1,j-1): #obs: the wPerDict values needed have been computed
                                    wVal1 = wPerDict[i,b]
                                    I1234Dict_full_aid[i,j] += I12Dict[b+1,j]*wVal1
                                    I1324Dict_full_aid[i,j] += I12Dict[b,j]*wVal1
        #                            I1234Dict_full[(i,j)] += (I12Dict[(b+1,j)] - I12Dict[(b+1,j-1)])*wDict[(i,b)]
        #                            I1324Dict_full[(i,j)] += (-I12Dict[(b,j)] + I12Dict[(b,j-1)])*wDict[(i,b)]
                                #I1234Dict_full add up terms and recursion: I1234(i,j)    
                                I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                                I1234Dict_full[i,j] += I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                                #I1324Dict_full add up terms recursion: I1324(i,j)
                                I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                                I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]                                                   
            elif degree == 6:
                if print_b > 0.5:
                    print "Computes degree 2, 4 and 6 measures." 
                for i in lChain[::-1]:  
                    for j in lChain:
                        wVal = 0 #reset for safety
                        if (j > i and j <=L-1):
                            if ((i in segNrRange[chain]) or (j in segNrRange[chain])):
            #                print "in change range %d , %d " % (i,j)
                                wVal = wPerDict[i,j]
                            else:
                                wVal = wDict[chain][i,j]
                                wPerDict[i,j] += wVal
                            I12 += wVal
                            I12Dict[i,j] += wVal + I12Dict[i,j-1] + I12Dict[i+1,j] - I12Dict[i+1,j-1]
                            #degree 4 measures:
                            #I1234: I1234(i,j;N)
                            I1234Dict[i,j] += I12Dict[j+1,L-1]*wVal + I1234Dict[i,j-1] + I1234Dict[i+1,j] - I1234Dict[i+1,j-1]
                            #I1423: I1423(i,j)                    
                            I1423Dict[i,j] += I12Dict[i+1,j-1]*wVal +  I1423Dict[i,j-1] + I1423Dict[i+1,j] - I1423Dict[i+1,j-1]
                            #I1324: I1324(i,j;N)
                            I1324Dict[i,j] += (I12Dict[i+1,L-1] - I12Dict[i+1,j] -I12Dict[j,L-1] )*wVal + I1324Dict[i,j-1] + I1324Dict[i+1,j] - I1324Dict[i+1,j-1]
                            for b in range(i+1,j-1): #obs: the wPerDict values needed have been computed
                                wVal1 = wPerDict[i,b]
                                I1234Dict_full_aid[i,j] += I12Dict[b+1,j]*wVal1
                                I1324Dict_full_aid[i,j] += I12Dict[b,j]*wVal1
                            #I1234Dict_full add up terms and recursion: I1234(i,j))    
                            I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                            I1234Dict_full[i,j] +=  I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                            #I1324Dict_full add up terms recursion: I1324(i,j)
                            I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                            I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]                            
                            #I1423_full: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet
                            I1423Dict_full0[i+1,j] += (I12Dict[i+2,j-1] - I12Dict[i+2,j])*(I12Dict[i+2,L-1] - I12Dict[i+2,j] - I12Dict[i+1,L-1] + I12Dict[i+1,j])
                            I1423Dict_full0[i+1,j] += I1423Dict_full0[i+2,j] + I1423Dict_full0[i+1,j-1] - I1423Dict_full0[i+2,j-1]
                            #to compute certain degree 6 measures we use two auxillary degree 4 
                            #measures; the recursion demands to sum j in the - direction (ie 
                            #from above and down):
                            k = L-1 + i -j+1
                            for c in range(k+1,L): #Obs: the wDict values needed have already been computed (since k,c >i)
                                wVal2 = wPerDict[k,c]
                                I1324Dict_full2_aid[i+1,k] += wVal2*(I12Dict[i+1,c-1] - I12Dict[i+1,k] - I12Dict[k,c-1])
                                I1423Dict_full2_aid[i+1,k] += wVal2*(I12Dict[i+1,L-1] - I12Dict[i+1,c] - I12Dict[k,L-1] + I12Dict[k,c]) 
                            #I1324Dict_full2, recursion: I1324(i;j,N)
                            I1324Dict_full2[i+1,k] += I1324Dict_full2_aid[i+1,k] + I1324Dict_full2[i+1,k+1]
                            #I1423Dict_full2, recursion: I1423(i;j,N)
                            I1423Dict_full2[i+1,k] += I1423Dict_full2_aid[i+1,k] + I1423Dict_full2[i+1,k+1] 
                            #degree 6 measures:
                            #(12)*:
                            I123456Dict[i,j] += I1234Dict[j+1,L-1]*wVal + I123456Dict[i+1,j] + I123456Dict[i,j-1] - I123456Dict[i+1,j-1] 
                            I123645Dict[i,j] += I1423Dict[j+1,L-1]*wVal + I123645Dict[i+1,j] + I123645Dict[i,j-1] - I123645Dict[i+1,j-1] 
                            I123546Dict[i,j] += I1324Dict[j+1,L-1]*wVal+ I123546Dict[i+1,j] + I123546Dict[i,j-1] - I123546Dict[i+1,j-1] 
                            #(13)*:
                            I132456Dict[i,j] += (I1234Dict[i+1,L-1] - I1234Dict[i+1,j] -I1234Dict[j,L-1])*wVal + I132456Dict[i+1,j] + I132456Dict[i,j-1] - I132456Dict[i+1,j-1]
                            I132546Dict[i+1,j] += (I1324Dict_full2[i+2,j+1] - I1324Dict_full[j,L-1])*wPerDict[i+1,j] + I132546Dict[i+2,j] + I132546Dict[i+1,j-1] - I132546Dict[i+2,j-1]
                            I132645Dict[i+1,j] += (I1423Dict_full2[i+2,j+1] - I1423Dict[j,L-1])*wPerDict[i+1,j] + I132645Dict[i+2,j] + I132645Dict[i+1,j-1] - I132645Dict[i+2,j-1]
                            #(14)*
                            #allowing direct computation:
                            I142356Dict[i,j] += I12Dict[i+1,j-1]*I12Dict[j+1,L-1]*wVal + I142356Dict[i+1,j] + I142356Dict[i,j-1] - I142356Dict[i+1,j-1]
                            #back to the more complicated parts, now I142536:
                            #three parts recognized in decompsition as already known:
                            I142536Dict[i+1,j] += wPerDict[i+1,j]*(I1324Dict_full[i+2,L-1] - I1324Dict[i+2,j] - I1324Dict_full2[i+2,j])
                            #recursion part:
                            I142536Dict[i+1,j] += I142536Dict[i+2,j] + I142536Dict[i+1,j-1] - I142536Dict[i+2,j-1] 
                            #Next complicated part, I142635:
                            #three parts recognized in decompsition as already known:
                            I142635Dict[i+1,j] += wPerDict[i+1,j]*(I1423Dict[i+2,L-1] - I1423Dict_full2[i+2,j] - I1423Dict_full0[i+2,j])
                            #recursion part:
                            I142635Dict[i+1,j] += I142635Dict[i+2,j] + I142635Dict[i+1,j-1] - I142635Dict[i+2,j-1]
                            #(15)*
                            I152346Dict[i,j] += wVal*(I1234Dict[i+1,j-1] - I1234Dict_full[i+1,j] - I12Dict[j,L-1]*I12Dict[i+1,j-1])
                            I152346Dict[i,j] += I152346Dict[i+1,j] + I152346Dict[i,j-1] - I152346Dict[i+1,j-1]
                            I152436Dict[i,j] += wVal*(I1324Dict[i+1,j-1] - I1324Dict_full[i+1,j])
                            I152436Dict[i,j] += I152436Dict[i+1,j] + I152436Dict[i,j-1] - I152436Dict[i+1,j-1]
                            I152634Dict[i,j] += wVal*(I1423Dict_full0[i+1,j-1] - I1423Dict[i+1,j])
                            I152634Dict[i,j] += I152634Dict[i+1,j] + I152634Dict[i,j-1] - I152634Dict[i+1,j-1]
                            #(16)*
                            I162345Dict[i,j] += wVal*I1234Dict_full[i+1,j-1] #ER FULL OK?
                            I162345Dict[i,j] += I162345Dict[i+1,j] + I162345Dict[i,j-1] - I162345Dict[i+1,j-1]
                            I162435Dict[i,j] += wVal*I1324Dict_full[i+1,j-1] #ER FULL OK?
                            I162435Dict[i,j] += I162435Dict[i+1,j] + I162435Dict[i,j-1] - I162435Dict[i+1,j-1]
                            I162534Dict[i,j] += wVal*I1423Dict[i+1,j-1]
                            I162534Dict[i,j] += I162534Dict[i+1,j] + I162534Dict[i,j-1] - I162534Dict[i+1,j-1]                    
                #the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:
                for j in lChain[::-1]:
                    i = 0
                    for c in range(j+1,L): #Obs: the 2d-loop above has sweeped the whole simplex so the wDictPer has been computed already 
                        wVal = wPerDict[j,c]
        #                        I1234Dict_full2_aid[(i,j)] +=  I12Dict[(i,c-1)]*wDict[(j,c)]
                        I1324Dict_full2_aid[i,j] += wVal*(I12Dict[i,c-1] - I12Dict[i,j] - I12Dict[j,c-1])
                        I1423Dict_full2_aid[i,j] += wVal*(I12Dict[i,L-1] - I12Dict[i,c] - I12Dict[j,L-1] + I12Dict[j,c]) 
                   #recursion:
                    I1324Dict_full2[i,j] += I1324Dict_full2_aid[i,j] + I1324Dict_full2[i,j+1]
                    I1423Dict_full2[i,j] += I1423Dict_full2_aid[i,j] + I1423Dict_full2[i,j+1] 
                #Missing terms for some of the degree 6 measures and I1423_full0:
                for j in lChain[1:]:
                    i=0
                    wVal = wPerDict[i,j] #Obs: the 2d-loop above has sweeped the whole simplex so the wDictPer has been computed already 
                    I1423Dict_full0[i,j] += (I12Dict[i+1,j-1] - I12Dict[i+1,j])*(I12Dict[i+1,L-1] - I12Dict[i+1,j] - I12Dict[i,L-1] + I12Dict[i,j])
                    I1423Dict_full0[i,j] += I1423Dict_full0[i+1,j] + I1423Dict_full0[i,j-1] - I1423Dict_full0[i+1,j-1]
                    #
                    I132546Dict[i,j] += (I1324Dict_full2[i+1,j+1] - I1324Dict_full[j,L-1])*wVal + I132546Dict[i+1,j] + I132546Dict[i,j-1] - I132546Dict[i+1,j-1]
                    I132645Dict[i,j] += (I1423Dict_full2[i+1,j+1] - I1423Dict[j,L-1])*wVal + I132645Dict[i+1,j] + I132645Dict[i,j-1] - I132645Dict[i+1,j-1]
                    #
                    I142536Dict[i,j] += wVal*(I1324Dict_full[i+1,L-1] - I1324Dict[i+1,j] - I1324Dict_full2[i+1,j])
                    I142536Dict[i,j] += I142536Dict[i+1,j] + I142536Dict[i,j-1] - I142536Dict[i+1,j-1] 
                    #
                    I142635Dict[i,j] += wVal*(I1423Dict[i+1,L-1] - I1423Dict_full2[i+1,j] - I1423Dict_full0[i+1,j])
                    I142635Dict[i,j] += I142635Dict[i+1,j] + I142635Dict[i,j-1] - I142635Dict[i+1,j-1]    
            else: 
                print "Degree must be 2, 4 or 6."
            #reset
            outDict = np.array([])
            #fetch results into array:                        
            outDict = np.array((wPerDict, I12Dict, 
                                I1234Dict, I1234Dict_full, 
                                I1324Dict, I1324Dict_full, I1324Dict_full2,
                                I1423Dict, I1423Dict_full0,I1423Dict_full2,
                                I123456Dict, I123546Dict, I123645Dict,
                                I132456Dict, I132546Dict, I132645Dict,
                                I142356Dict, I142536Dict, I142635Dict,
                                I152346Dict, I152436Dict, I152634Dict,
                                I162345Dict, I162435Dict, I162534Dict)) 
            #store the results:
            if p < 1:
                tot_outDict[chain] = np.array([outDict])
            else:
                tot_outDict[chain] = np.append(tot_outDict[chain], np.array([outDict]), axis =0)
    return tot_outDict, perturbedChains



def Iperturbed_vect_multi_dict(CaChain, 
               polyCaChain, 
               wDict,
               nrPerturbations = 1,
               residueNrRange = range(0,5),
               perturbationLength = 5, 
                degree = 4,
                full_b = 0,
                print_b = 0):
    '''Multi perturbation version of Iperturbed_vect_dict.'''
    #for holding output:
    tot_outDict = {}
    #Generate a set of perturbations
    perturbedChains = {}
    perturbations = {}
    for i in range(nrPerturbations):
        #pick randomly the perturbations:
        perturbation = np.random.rand(perturbationLength,3) 
        #compute the polygonal chain of the perturbed chain:
        perCaChain, perPolyCaChain, segNrRange = perturbPolygonalChain(CaChain = CaChain, 
                                                        polyCaChain = polyCaChain,
                                                        residueNrRange = residueNrRange, 
                                                        perturbation = perturbation)
        for chain in CaChain.keys():
            if not(perturbedChains.has_key(chain)):
                perturbedChains[chain] = {}
            perturbedChains[chain][i] = [perCaChain[chain], perPolyCaChain[chain], segNrRange[chain]]
            #the actual perturbations are not used beow, but could be useful to know, so we record them:
            if not(perturbations.has_key(chain)):
                perturbations[chain] = []
            perturbations[chain].append(perturbation)
 
    #compute invariants:
    perturbedChainList = []
    for chain in CaChain.keys():
        perturbedChainsList = perturbedChains[chain]
        if print_b > 0.5:
            print segNrRange[chain]
        #Get length of chain:
        L = len(polyCaChain[chain])
        if print_b > 0.5:
            print "Length of chain is: %d" % L            
    

        #compute the writhe terms in one go (call dictionary based version by setting dict_b = 1):
        wPerDict_multi = wAll_multi(chainsInput = perturbedChainsList, lengthPolyChain = L, print_b = print_b, dict_b = 1)
        #compute the measures, running through the perturbations:
        lChain = range(L)       
        #init. dictionary for storing results (all perturbations):
        tot_outDict[chain] = {}     
        for p in perturbedChains[chain].keys(): #p for perturbation
            segNrRange[chain] = perturbedChains[chain][p][2]
            firstSeg = segNrRange[chain][0]
            lastSeg = segNrRange[chain][::-1][0]
            if print_b > 0.5:
                print segNrRange[chain]
            #reset:
            #init. dictionaries for storing results for the current perturbation:
            I12 = 0
            I12Dict = {}
            I1234Dict = {}
            I1234Dict_full = {}
            I1234Dict_full_aid = {}
            I1234Dict_full2 = {}
            I1234Dict_full2_aid = {}
            I1423Dict = {} 
            I1423Dict_full0 = {} 
            I1423Dict_full2 = {}
            I1423Dict_full2_aid = {}
            I1324Dict = {}
            I1324Dict_full = {}
            I1324Dict_full_aid = {}
            I1324Dict_full2 = {}
            I1324Dict_full2_aid = {}
            #(12)*:
            I123456Dict = {} 
            I123645Dict = {}
            I123546Dict = {}
            #(13)*:
            I132456Dict = {}
            I132546Dict = {}
            I132645Dict = {}
            #(14)*
            I142356Dict = {}
            I142536Dict = {}
            I142635Dict = {}
            #(15)*
            I152346Dict = {}
            I152436Dict = {}
            I152634Dict = {}
            #(16)*
            I162345Dict = {}
            I162435Dict = {}
            I162534Dict = {}    
            for i in range(L+1):
                for j in range(L+1):
                    I12Dict[(i,j)] = 0
        #            I12Dict[(-1,j)] = 0
                    if degree > 2: 
                        I1234Dict[(i,j)] = 0
                        I1234Dict_full[(i,j)] = 0
                        I1324Dict[(i,j)] = 0  
                        I1423Dict[(i,j)] = 0
        #                if full_b ==1:
                        I1234Dict_full[(i,j)] = 0
                        I1234Dict_full_aid[(i,j)] = 0
                        I1234Dict_full2[(i,j)] = 0
                        I1234Dict_full2_aid[(i,j)] = 0
                        I1324Dict_full[(i,j)] = 0
                        I1324Dict_full_aid[(i,j)] = 0
                        I1324Dict_full2[(i,j)] = 0
        #                I1324Dict_full2[(-1,j)] = 0
                        I1324Dict_full2_aid[(i,j)] = 0
        #                I1324Dict_full2_aid[(-1,j)] = 0
                        I1423Dict_full0[(i,j)] = 0
                        I1423Dict_full2[(i,j)] = 0
                        I1423Dict_full2_aid[(i,j)] = 0
                    if degree > 4:
                        I123456Dict[(i,j)] = 0 
                        I123645Dict[(i,j)] = 0
                        I123546Dict[(i,j)] = 0
                        #(13)*:
                        I132456Dict[(i,j)] = 0
                        I132546Dict[(i,j)] = 0
                        I132645Dict[(i,j)] = 0
                        #(14)*
                        I142356Dict[(i,j)] = 0
                        I142536Dict[(i,j)] = 0
                        I142635Dict[(i,j)] = 0
                        #(15)*
                        I152346Dict[(i,j)] = 0
                        I152436Dict[(i,j)] = 0
                        I152634Dict[(i,j)] = 0
                        #(16)*
                        I162345Dict[(i,j)] = 0
                        I162435Dict[(i,j)] = 0
                        I162534Dict[(i,j)] = 0
            #Main part: computing the measures using recursion formulas
            #one block for each choice of degree; only one block gets executed:   
            wPerDict = wPerDict_multi[p]   
            if degree == 2:
                for i in lChain[::-1]:
                    for j in lChain:
                        wVal = 0 #reset
                        if (j > i and j <=L-1):
                            if ((i in segNrRange[chain]) or (j in segNrRange[chain])):
            #                print "in change range %d , %d " % (i,j)
                                wVal = wPerDict[(i,j)]
                            else:
                                wVal = wDict[chain][(i,j)]
                                wPerDict[(i,j)] += wVal
                            I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
            elif degree == 4:
                for i in lChain[::-1]:
                    for j in lChain:
                        wVal = 0 #reset
                        if (j > i and j <=L-1):
                            if ((i in segNrRange[chain]) or (j in segNrRange[chain])):
    #                            print "in change range %d , %d " % (i,j)
                                wVal = wPerDict[(i,j)]
    #                            print "wPerDict: %f " % wVal
    #                            print "wDict: %f " % wDict[(i,j)]
                            else:
                                wVal = wDict[chain][(i,j)]
                                wPerDict[(i,j)] += wVal
                            I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
                            I1234Dict[(i,j)] += I12Dict[(j+1,L-1)]*wVal + I1234Dict[(i,j-1)] + I1234Dict[(i+1,j)] - I1234Dict[(i+1,j-1)]
                            I1423Dict[(i,j)] += I12Dict[(i+1,j-1)]*wVal +  I1423Dict[(i,j-1)] + I1423Dict[(i+1,j)] - I1423Dict[(i+1,j-1)]
                            I1324Dict[(i,j)] += (I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] -I12Dict[(j,L-1)] )*wVal + I1324Dict[(i,j-1)] + I1324Dict[(i+1,j)] - I1324Dict[(i+1,j-1)]
                            if full_b == 1:
                                for b in range(i+1,j-1): #obs: the wPerDict values needed have been computed
                                    wVal1 = wPerDict[(i,b)]
                                    I1234Dict_full_aid[(i,j)] += I12Dict[(b+1,j)]*wVal1
                                    I1324Dict_full_aid[(i,j)] += I12Dict[(b,j)]*wVal1
        #                            I1234Dict_full[(i,j)] += (I12Dict[(b+1,j)] - I12Dict[(b+1,j-1)])*wDict[(i,b)]
        #                            I1324Dict_full[(i,j)] += (-I12Dict[(b,j)] + I12Dict[(b,j-1)])*wDict[(i,b)]
                                #I1234Dict_full add up terms and recursion: I1234(i,j))    
                                I1234Dict_full[i,j] += I1234Dict_full_aid[i,j] -I1234Dict_full_aid[i,j-1]
                                I1234Dict_full[i,j] +=  I1234Dict_full[i,j-1] + I1234Dict_full[i+1,j] - I1234Dict_full[i+1,j-1]
                                #I1324Dict_full add up terms recursion: I1324(i,j)
                                I1324Dict_full[i,j] += I1324Dict_full_aid[i,j-1] - I1324Dict_full_aid[i,j]   
                                I1324Dict_full[i,j] += (I12Dict[i,j-1] - I12Dict[i+1,j-1])*(I12Dict[i+1,j] - I12Dict[i+1,j-1]) + I1324Dict_full[i,j-1] + I1324Dict_full[i+1,j] - I1324Dict_full[i+1,j-1]            
            elif degree == 6:
                if print_b > 0.5:
                    print "Computes degree 2, 4 and 6 measures." 
                for i in lChain[::-1]:  
                    for j in lChain:
                        wVal = 0 #reset for safety
                        if (j > i and j <=L-1):
                            if ((i in segNrRange[chain]) or (j in segNrRange[chain])):
            #                print "in change range %d , %d " % (i,j)
                                wVal = wPerDict[(i,j)]
                            else:
                                wVal = wDict[chain][(i,j)]
                                wPerDict[(i,j)] += wVal
                            I12 += wVal
                            I12Dict[(i,j)] += wVal + I12Dict[(i,j-1)] + I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]
                            #degree 4 measures:
                            #I1234: I1234(i,j;N)
                            I1234Dict[(i,j)] += I12Dict[(j+1,L-1)]*wVal + I1234Dict[(i,j-1)] + I1234Dict[(i+1,j)] - I1234Dict[(i+1,j-1)]
                            #I1423: I1423(i,j)                    
                            I1423Dict[(i,j)] += I12Dict[(i+1,j-1)]*wVal +  I1423Dict[(i,j-1)] + I1423Dict[(i+1,j)] - I1423Dict[(i+1,j-1)]
                            #I1324: I1324(i,j;N)
                            I1324Dict[(i,j)] += (I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] -I12Dict[(j,L-1)] )*wVal + I1324Dict[(i,j-1)] + I1324Dict[(i+1,j)] - I1324Dict[(i+1,j-1)]
                            for b in range(i+1,j-1): #obs: the wPerDict values needed have been computed
                                wVal1 = wPerDict[(i,b)]
                                I1234Dict_full_aid[(i,j)] +=  I12Dict[(b+1,j)]*wVal1
                                I1324Dict_full_aid[(i,j)] += I12Dict[(b,j)]*wVal1
                            #I1234Dict_full add up terms and recursion: I1234(i,j) 
                            I1234Dict_full[(i,j)] += I1234Dict_full_aid[(i,j)] -I1234Dict_full_aid[(i,j-1)]
                            I1234Dict_full[(i,j)] +=  I1234Dict_full[(i,j-1)] + I1234Dict_full[(i+1,j)] - I1234Dict_full[(i+1,j-1)]
                            #I1324Dict_full add up terms recursion: I1324(i,j)
                            I1324Dict_full[(i,j)] += I1324Dict_full_aid[(i,j-1)] - I1324Dict_full_aid[(i,j)]   
                            I1324Dict_full[(i,j)] += (I12Dict[(i,j-1)] - I12Dict[(i+1,j-1)])*(I12Dict[(i+1,j)] - I12Dict[(i+1,j-1)]) + I1324Dict_full[(i,j-1)] + I1324Dict_full[(i+1,j)] - I1324Dict_full[(i+1,j-1)]                            
                            #I1423_full: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet
                            I1423Dict_full0[(i+1,j)] += (I12Dict[(i+2,j-1)] - I12Dict[(i+2,j)])*(I12Dict[(i+2,L-1)] - I12Dict[(i+2,j)] - I12Dict[(i+1,L-1)] + I12Dict[(i+1,j)])
                            I1423Dict_full0[(i+1,j)] += I1423Dict_full0[(i+2,j)] + I1423Dict_full0[(i+1,j-1)] - I1423Dict_full0[(i+2,j-1)]
                            #to compute certain degree 6 measures we use two auxillary degree 4 measures; the recursion
                            #demands to sum j in the - direction (ie from above and down):
                            k = L-1 + i -j+1
                            for c in range(k+1,L): #obs: the wPerDict values needed have been computed (since k,c >i)
                                wVal2 = wPerDict[(k,c)]
                                I1324Dict_full2_aid[(i+1,k)] += wVal2*(I12Dict[(i+1,c-1)] - I12Dict[(i+1,k)] - I12Dict[(k,c-1)])
                                I1423Dict_full2_aid[(i+1,k)] += wVal2*(I12Dict[(i+1,L-1)] - I12Dict[(i+1,c)] - I12Dict[(k,L-1)] + I12Dict[(k,c)]) 
                            #I1324Dict_full2, recursion: I1324(i;j,N)
                            I1324Dict_full2[(i+1,k)] += I1324Dict_full2_aid[(i+1,k)] + I1324Dict_full2[(i+1,k+1)]
                            #I1423Dict_full2, recursion: I1423(i;j,N)
                            I1423Dict_full2[(i+1,k)] += I1423Dict_full2_aid[(i+1,k)] + I1423Dict_full2[(i+1,k+1)] 
                            #degree 6 measures:
                            #(12)*:
                            I123456Dict[(i,j)] += I1234Dict[(j+1,L-1)]*wVal + I123456Dict[(i+1,j)] + I123456Dict[(i,j-1)] - I123456Dict[(i+1,j-1)] 
                            I123645Dict[(i,j )] += I1423Dict[(j+1,L-1)]*wVal + I123645Dict[(i+1,j)] + I123645Dict[(i,j-1)] - I123645Dict[(i+1,j-1)] 
                            I123546Dict[(i,j)] += I1324Dict[(j+1,L-1)]*wVal+ I123546Dict[(i+1,j)] + I123546Dict[(i,j-1)] - I123546Dict[(i+1,j-1)] 
                            #(13)*:
                            I132456Dict[(i,j)] += (I1234Dict[(i+1,L-1)] - I1234Dict[(i+1,j)] -I1234Dict[(j,L-1)])*wVal + I132456Dict[(i+1,j)] + I132456Dict[(i,j-1)] - I132456Dict[(i+1,j-1)]
                            I132546Dict[(i+1,j)] += (I1324Dict_full2[(i+2,j+1)] - I1324Dict_full[(j,L-1)])*wPerDict[(i+1,j)] + I132546Dict[(i+2,j)] + I132546Dict[(i+1,j-1)] - I132546Dict[(i+2,j-1)]
                            I132645Dict[(i+1,j)] += (I1423Dict_full2[(i+2,j+1)] - I1423Dict[(j,L-1)])*wPerDict[(i+1,j)] + I132645Dict[(i+2,j)] + I132645Dict[(i+1,j-1)] - I132645Dict[(i+2,j-1)]
                            #(14)*
                            #allowing direct computation:
                            I142356Dict[(i,j)] += I12Dict[(i+1,j-1)]*I12Dict[(j+1,L-1)]*wVal + I142356Dict[(i+1,j)] + I142356Dict[(i,j-1)] - I142356Dict[(i+1,j-1)]
                            #back to the more complicated parts, now I142536:
                            #three parts recognized in decompsition as already known:
                            I142536Dict[(i+1,j)] += wPerDict[(i+1,j)]*(I1324Dict_full[(i+2,L-1)] - I1324Dict[(i+2,j)] - I1324Dict_full2[(i+2,j)])
                            #recursion part:
                            I142536Dict[(i+1,j)] += I142536Dict[(i+2,j)] + I142536Dict[(i+1,j-1)] - I142536Dict[(i+2,j-1)] 
                            #Next complicated part, I142635:
                            #three parts recognized in decompsition as already known:
                            I142635Dict[(i+1,j)] += wPerDict[(i+1,j)]*(I1423Dict[(i+2,L-1)] - I1423Dict_full2[(i+2,j)] - I1423Dict_full0[(i+2,j)])
                            #recursion part:
                            I142635Dict[(i+1,j)] += I142635Dict[(i+2,j)] + I142635Dict[(i+1,j-1)] - I142635Dict[(i+2,j-1)]
                            #(15)*
                            I152346Dict[(i,j)] += wVal*(I1234Dict[(i+1,j-1)] - I1234Dict_full[(i+1,j)] - I12Dict[(j,L-1)]*I12Dict[(i+1,j-1)])
                            I152346Dict[(i,j)] += I152346Dict[(i+1,j)] + I152346Dict[(i,j-1)] - I152346Dict[(i+1,j-1)]
                            I152436Dict[(i,j)] += wVal*(I1324Dict[(i+1,j-1)] - I1324Dict_full[(i+1,j)])
                            I152436Dict[(i,j)] += I152436Dict[(i+1,j)] + I152436Dict[(i,j-1)] - I152436Dict[(i+1,j-1)]
                            I152634Dict[(i,j)] += wVal*(I1423Dict_full0[(i+1,j-1)] - I1423Dict[(i+1,j)])
                            I152634Dict[(i,j)] += I152634Dict[(i+1,j)] + I152634Dict[(i,j-1)] - I152634Dict[(i+1,j-1)]
                            #(16)*
                            I162345Dict[(i,j)] += wVal*I1234Dict_full[(i+1,j-1)] #ER FULL OK?
                            I162345Dict[(i,j)] += I162345Dict[(i+1,j)] + I162345Dict[(i,j-1)] - I162345Dict[(i+1,j-1)]
                            I162435Dict[(i,j)] += wVal*I1324Dict_full[(i+1,j-1)] #ER FULL OK?
                            I162435Dict[(i,j)] += I162435Dict[(i+1,j)] + I162435Dict[(i,j-1)] - I162435Dict[(i+1,j-1)]
                            I162534Dict[(i,j)] += wVal*I1423Dict[(i+1,j-1)]
                            I162534Dict[(i,j)] += I162534Dict[(i+1,j)] + I162534Dict[(i,j-1)] - I162534Dict[(i+1,j-1)]                    
                #the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:
                for j in lChain[::-1]:
                    i = 0
                    for c in range(j+1,L):
                        wVal = wPerDict[(j,c)] #Obs: the 2d-loop above has sweeped the whole simplex so the wDictPer has been computed already 
        #                        I1234Dict_full2_aid[(i,j)] +=  I12Dict[(i,c-1)]*wDict[(j,c)]
                        I1324Dict_full2_aid[(i,j)] += wVal*(I12Dict[(i,c-1)] - I12Dict[(i,j)] - I12Dict[(j,c-1)])
                        I1423Dict_full2_aid[(i,j)] += wVal*(I12Dict[(i,L-1)] - I12Dict[(i,c)] - I12Dict[(j,L-1)] + I12Dict[(j,c)]) 
                   #recursion:
                    I1324Dict_full2[(i,j)] += I1324Dict_full2_aid[(i,j)] + I1324Dict_full2[(i,j+1)]
                    I1423Dict_full2[(i,j)] += I1423Dict_full2_aid[(i,j)] + I1423Dict_full2[(i,j+1)] 
                #Missing terms for some of the degree 6 measures and I1423_full0:
                for j in lChain[1:]:
                    i=0
                    wVal = wPerDict[(i,j)] #Obs: the 2d-loop above has sweeped the whole simplex so the wDictPer has been computed already 
                    I1423Dict_full0[(i,j)] += (I12Dict[(i+1,j-1)] - I12Dict[(i+1,j)])*(I12Dict[(i+1,L-1)] - I12Dict[(i+1,j)] - I12Dict[(i,L-1)] + I12Dict[(i,j)])
                    I1423Dict_full0[(i,j)] += I1423Dict_full0[(i+1,j)] + I1423Dict_full0[(i,j-1)] - I1423Dict_full0[(i+1,j-1)]
                    #
                    I132546Dict[(i,j)] += (I1324Dict_full2[(i+1,j+1)] - I1324Dict_full[(j,L-1)])*wVal + I132546Dict[(i+1,j)] + I132546Dict[(i,j-1)] - I132546Dict[(i+1,j-1)]
                    I132645Dict[(i,j)] += (I1423Dict_full2[(i+1,j+1)] - I1423Dict[(j,L-1)])*wVal + I132645Dict[(i+1,j)] + I132645Dict[(i,j-1)] - I132645Dict[(i+1,j-1)]
                    #
                    I142536Dict[(i,j)] += wVal*(I1324Dict_full[(i+1,L-1)] - I1324Dict[(i+1,j)] - I1324Dict_full2[(i+1,j)])
                    I142536Dict[(i,j)] += I142536Dict[(i+1,j)] + I142536Dict[(i,j-1)] - I142536Dict[(i+1,j-1)] 
                    #
                    I142635Dict[(i,j)] += wVal*(I1423Dict[(i+1,L-1)] - I1423Dict_full2[(i+1,j)] - I1423Dict_full0[(i+1,j)])
                    I142635Dict[(i,j)] += I142635Dict[(i+1,j)] + I142635Dict[(i,j-1)] - I142635Dict[(i+1,j-1)]    
            else: 
                print "Degree must be 2, 4 or 6."
            #reset:    
            outDict = {}
            #fetch results into dictionary
            outDict['w'] = wPerDict
            outDict['I12'] = I12Dict
            outDict['I1234'] = I1234Dict
            outDict['I1234_full'] = I1234Dict_full
            outDict['I1324'] = I1324Dict
            outDict['I1324_full'] = I1324Dict_full
            outDict['I1324_full2'] = I1324Dict_full2    
            outDict['I1423'] = I1423Dict
            outDict['I1423_full0'] = I1423Dict_full0
            outDict['I1423_full2'] = I1423Dict_full2    
        
            outDict['I123456'] = I123456Dict
            outDict['I123546'] = I123546Dict
            outDict['I123645'] = I123645Dict
        
            outDict['I132456'] = I132456Dict
            outDict['I132546'] = I132546Dict
            outDict['I132645'] = I132645Dict
        
            outDict['I142356'] = I142356Dict
            outDict['I142536'] = I142536Dict
            outDict['I142635'] = I142635Dict
            
            outDict['I152346'] = I152346Dict
            outDict['I152436'] = I152436Dict
            outDict['I152634'] = I152634Dict
        
            outDict['I162345'] = I162345Dict
            outDict['I162435'] = I162435Dict
            outDict['I162534'] = I162534Dict
            #store the results:
            tot_outDict[chain][p] = {}
            tot_outDict[chain][p] = outDict
    return tot_outDict, perturbedChains


#########################################################################################################                    
## Utils for time testing:
#########################################################################################################

  

def timeTest(PDBfilename = PDBfilenameDef, residueNrRange =range(20,25), perturbation = [(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3)], degree =4, full_b = 1, print_b = 0, vect_b = 0, dict_b = 1):
    if dict_b > 0:        
        start1 = time.time()
        #Compute measures:
        outDict, CaChain, polyCaChain = I_dict(PDBfilename = PDBfilename,degree = degree, full_b = full_b, print_b = print_b, vect_b = vect_b)
    #    I1234_bf(PDBfilename =PDBfilename)
        
        start2 = time.time()  
        #Compute measures on perturbed chain:
        outDict, perCaChain, perPolyCaChain, segNrRange = Iperturbed_dict(CaChain = CaChain, polyCaChain =polyCaChain, wDict = outDict['w'], residueNrRange = residueNrRange, perturbation = perturbation , degree = degree, full_b = full_b, print_b = print_b, vect_b = vect_b)
        
        d = time.time()-start2
        print 'I comp took', start2-start1, 'seconds.'
        print 'Iper comp took', d, 'seconds.'
        
    else:
        start1 = time.time()
        #Compute measures:
        outDict, CaChain, polyCaChain = I(PDBfilename = PDBfilename,degree = degree, full_b = full_b, print_b = print_b, vect_b = vect_b)
    #    I1234_bf(PDBfilename =PDBfilename)
        
        start2 = time.time()  
        #Compute measures on perturbed chain:
        outDict, perCaChain, perPolyCaChain, segNrRange = Iperturbed(CaChain = CaChain, polyCaChain =polyCaChain, wDict = outDict[0], residueNrRange = residueNrRange, perturbation = perturbation , degree = degree, full_b = full_b, print_b = print_b, vect_b = vect_b)
        
        d = time.time()-start2
        print 'I comp took', start2-start1, 'seconds.'
        print 'Iper comp took', d, 'seconds.'


def timeTest_vect_multi(PDBfilename = PDBfilenameDef, nrPerturbations = 1, residueNrRange = range(20,25), perturbationLength = 5, degree =4, full_b = 1, print_b =0, vect_b = 1, dict_b = 0):
    if dict_b > 0:
        start1 = time.time()
        #Compute measures:
        outDict, CaChain, polyCaChain = I_dict(PDBfilename = PDBfilename,degree = degree, full_b = full_b, print_b = print_b, vect_b = vect_b)
    #    I1234_bf(PDBfilename =PDBfilename)
        
        start2 = time.time()  
        #Compute measures on perturbed chain:
        tot_outDict, perturbedChains = Iperturbed_vect_multi_dict(CaChain = CaChain, polyCaChain =polyCaChain, wDict = outDict['w'], nrPerturbations = nrPerturbations, residueNrRange = residueNrRange, perturbationLength = perturbationLength, degree = degree, full_b = 0, print_b = 0)
        
        d = time.time()-start2
        print 'I comp took', start2-start1, 'seconds.'
        print 'Iper comp took', d, 'seconds.'
    else:   
        start1 = time.time()
    #Compute measures:
        outDict, CaChain, polyCaChain = I(PDBfilename = PDBfilename,degree = degree, full_b = full_b, print_b = print_b, vect_b = vect_b)
    #    I1234_bf(PDBfilename =PDBfilename)
        
        start2 = time.time()  
        #Compute measures on perturbed chain:
        tot_outDict, perturbedChains = Iperturbed_vect_multi(CaChain = CaChain, polyCaChain =polyCaChain, wDict = outDict[0], nrPerturbations = nrPerturbations, residueNrRange = residueNrRange, perturbationLength = perturbationLength, degree = degree, full_b = 0, print_b = 0)
        d = time.time()-start2
        print 'I comp took', start2-start1, 'seconds.'
        print 'Iper comp took', d, 'seconds.'


def transposeDict(inDict):
    ''' for interchaning inner and outer keys in a given dictionary, inDict.'''
    out = {}
    for k1 in inDict.keys():
        for k2 in inDict[k1].keys():
            if not(out.has_key(k2)):
                out[k2] = {}
            out[k2][k1] = inDict[k1][k2]
    return out

def compareDicts(d1,d2, tolerance = 1e-6, rel_b = 0):
    '''Compares two dictionaries.'''
    k_diff = []
    diff = {}
    #first check that the keys are id:
    for k in d1.keys():
        if not(d2.has_key(k)):
#            print "Key: %s is in d1 but not in d2!" % k
            if abs(d1[k]) > tolerance:
                print "Error: dictionaries do not have identical set of keys; d1-value at key %s is above tolerance" % str(k)
                return k_diff, diff
            else:
                pass #print "Warning: dictionaries do not have identical set of keys."
    for k in d2.keys():
        if not(d1.has_key(k)):
#            print "Key: %s is in d2 but not in d1!" % k
            if abs(d2[k]) > tolerance:
                print "Error: dictionaries do not have identical set of keys; d2-value at key %s is above tolerance" % str(k)
                return k_diff, diff
            else:                
                pass #print "Warning: dictionaries do not have identical set of keys."

    #Then check the values:
    for k in d1.keys():
        if rel_b ==0:
            if abs(d2[k] - d1[k]) > tolerance:
                print "At %s d1 is: %f and d2 is %f " % (k, d1[k],d2[k])
                k_diff.append(k)
                diff[k] = d2[k] - d1[k]
        if rel_b ==1:
            if abs(d1[k]) > 1:
                if abs(float(d2[k] - d1[k])/d1[k]) > tolerance:
                    print "At %s d1 is: %f and d2 is %f " % (k, d1[k],d2[k])
                    k_diff.append(k) 
                    diff[k] = d2[k] - d1[k]
            else:
                if abs(d2[k] - d1[k]) > tolerance:
                    print "At %s d1 is: %f and d2 is %f " % (k, d1[k],d2[k])
                    k_diff.append(k)
                    diff[k] = d2[k] - d1[k]
    k_diff.sort
    return k_diff, diff
    
def compareDictsII(d1,d2, tolerance = 1e-6, rel_b = 0):
    '''Compares two dictionaries of dictionaries.'''
    k_diff = []
    diff = {}
    #first check that the keys are id:
    for k in d1.keys():
        if not(d2.has_key(k)):
            print "Key: %s is in d1 but not in d2!" % k
    for k in d2.keys():
        if not(d1.has_key(k)):
            print "Key: %s is in d2 but not in d1!" % k
    #Then check the values:
    for k in d1.keys():
        k_d, d = compareDicts(d1[k],d2[k], tolerance = tolerance, rel_b = rel_b)
        k_diff.append(k_d)
        diff[k] = d
    k_diff.sort
    return k_diff, diff    
    
def compareArrays(a1, a2, tolerance = 1e-6):
    diff = []
    print "length a1: %d a2: %d" % (len(a1), len(a2))
    for i in range(len(a1)):
        for j in range(len(a1[i])):
            for k in range(len(a1[i][j])):
                d = a1[i][j][k] - a2[i][j][k]
                if abs(d) > tolerance:
                    print "diff above tolerance at i,j,k: %d ,%d, %d" % (i,j,k)
                    diff.append(d)
    return diff                
 
########################################################################
#### plot functions for results
##########################################################################
 


#PDBlist9 = ['1notH','1mctIH', '1bpiH', '1notH','1mctIH', '1ptxH', '1notH','1mctIH', '1ptxH'] 

PDBlist9 = ['1notH','1mctIH', '1bpiH', '2mcmH', '1ttaAH', '1osaH', '1dadH',  '2cbaH', '2er7H']
#PDBlist9 = ['2mcmH', '1ttaAH', '1osaH']
#PDBlist9 = ['1dadH']
#PDBlist9 = ['2cbaH']
#PDBlist9 = ['2er7H']

def getPDBlist(path = PDBtop100Folder, fileName = r"C:\Users\Christian\Bioinformatics\projects\knots\code\top100H_source\top100names.txt") :
    '''Returns a list of the PDB filenames in the specied path/filename.
    Length of the structure is computed and a list containing (length, name)
    is returned.'''
    pdbList = []
    avg = 0
    with open(fileName, 'r') as files:
        i=0
        for f in files:
            PDBfilename = path + f.strip( )
#            print PDBfilename
            try:
                CaChain, pChain = getPolygonalChain(PDBfilename, outNumpy_b = 0)
                for chainId in CaChain.keys():
                    L = len(CaChain[chainId])
                    print (f.strip(), L)
                    pdbList.append([L,f.strip()])
                    avg = avg + L
            except ValueError:
                continue
            except KeyError:
                continue
            i+=1
        print "Number of structure read: %d" % i
        avg = float(avg)/i
        print "Average length: %f " % avg
    pdbList.sort()
    return pdbList

#top100: Average length: 187.381443 

def I_time_degree_plot(PDBnames = PDBlist9):
    '''Time for measure computation across all names in input list by 
    max-degree of measures computed; non/vectorized array/dict versions.
    Input: list of PDB-names of length 9 (else modify code).''' 
    i=0
    for name in PDBnames:
        f = PDBtop100Folder + name
        CaChain, pChain = getPolygonalChain(PDBfilename = f, outNumpy_b = 0)
        L = len(CaChain)
        i +=1
        for d in [0,1]:
            for v in [0,1]:
                timeList = []
                for degree in [0,2,4,6]:
                    if d == 0:
                        statement = 'git.I(PDBfilename =' + 'r"' + f + '", full_b = 0, print_b=0, degree =' + str(degree) + ', vect_b = ' + str(v) +')'
                    else:
                        statement = 'git.I_dict(PDBfilename =' + 'r"' + f + '", full_b = 0, print_b=0, degree =' + str(degree) + ', vect_b = ' + str(v) +')'
                    t_I = timeit.Timer(statement,'import git')
                    timeList.append(t_I.timeit(number=1))
                plt.subplot(3,3,i)
                plt.plot([0,2,4,6], timeList, label = 'dict_b =' + str(d)+ ', vect_b ='+ str(v))
#        plt.tight_layout
        if i > 6:
            plt.xlabel('degree', fontsize = 'small')
        plt.ylabel('time (s)', fontsize = 'small')
        plt.legend(loc = "upper left", fontsize = 'x-small')
        plt.title(name + ', L =' + str(L), fontsize = 'medium')

        

def Ipert_time_degree_plot(PDBnames = PDBlist9, degreeList = [2,4,6]):
    '''Time for measure computation of single perturbation across names in 
    input list by max-degree of measures computed; non/vectorized array/dict 
    versions.
    Input: list of PDB-names of length 9 (else modify code).''' 
    i=0
    for name in PDBnames:
        f = PDBtop100Folder + name
        i +=1
        for d in [0,1]:
            #we need also the w-values so we call that computation (vectorized): 
            if d < 1:
                outDict_org, CaChain, polyCaChain = I(PDBfilename = f, degree =0, inNumpy_b =0 , print_b = 0, vect_b = 1)
            else:
                outDict_org, CaChain, polyCaChain = I_dict(PDBfilename = f, degree =0, inNumpy_b =0 , print_b = 0, vect_b = 1)                
            L = len(CaChain)
            #we will perturb the chains so, as they have different lengths, we must use different ranges.
            pRange = range(int(L/2)-2, int(L/2) + 3) #a range of lgth 5 around middle of chain (all have length > 10) 
            for v in [0,1]:
                timeList = []
                for degree in degreeList:
                    if d < 1:
                        start = time.time()  
                        #Compute measures on perturbed chain:
                        outDict, perCaChain, perPolyCaChain, segNrRange = Iperturbed(CaChain = CaChain, polyCaChain =polyCaChain, wDict = outDict_org[0], residueNrRange = pRange, degree = degree, full_b = 0, print_b = 0, vect_b = v)       
                        d_t = time.time()-start
#                        statement = 'git.Iperturbed(CaChain = ' + str(CaChain) + ', polyCaChain =' + str(polyCaChain) + ', wDict = '+ str(outDict[0]) +', residueNrRange =' + str(pRange) +', degree =' + str(degree) + ', vect_b = ' + str(v) +')'
                    else:
                        start = time.time()  
                        outDict, perCaChain, perPolyCaChain, segNrRange = Iperturbed_dict(CaChain = CaChain, polyCaChain =polyCaChain, wDict = outDict_org['w'], residueNrRange = pRange, degree = degree, full_b = 0, print_b = 0, vect_b = v)       
                        d_t = time.time()-start
#                        statement = 'git.Iperturbed_dict(CaChain = ' + str(CaChain) + ', polyCaChain =' + str(polyCaChain) + ', wDict = ' + str(outDict['w']) +', residueNrRange =' + pRange + ', degree =' + str(degree) + ', vect_b = ' + str(v) +')'
#                    t_I = timeit.Timer(statement,'import git')
#                    timeList.append(t_I.timeit(number=1))
                    timeList.append(d_t)
                plt.subplot(3,3,i)
                plt.plot(degreeList, timeList, label = 'dict_b =' + str(d)+ ', vect_b ='+ str(v))
#        plt.tight_layout
        if i > 6:
            plt.xlabel('degree', fontsize = 'small')
        plt.ylabel('time (s)', fontsize = 'small')
        plt.legend(loc = "upper left", fontsize = 'x-small')
        plt.title(name + ', L =' + str(L), fontsize = 'medium')
    

def I_time_length_plot(PDBlist, degree = 4):
    '''Time for computation of measures up to specified degree across all names 
    in input list; non/vectorized array/dict versions.
    Input: as output from getPDBlist.''' 
    for d in [0,1]:
        timeList = []
        lengthList = []
        for l_item in PDBlist:
            f = PDBtop100Folder + l_item[1]
            L = l_item[0]
            if d ==0:
                statement = 'git.I(PDBfilename =' + 'r"' + f + '", full_b = 0, print_b=0, degree = '+ str(degree) +' ,vect_b = 1)'
            else:
                statement = 'git.I_dict(PDBfilename =' + 'r"' + f + '", full_b = 0, print_b=0, degree ='+ str(degree) +',vect_b = 1)'
            t_I = timeit.Timer(statement,'import gitP as git')
            timeList.append(t_I.timeit(number=1))
            lengthList.append(L)
        plt.plot(lengthList, timeList, label = 'dict_b ='+ str(d))
        plt.xlabel('length', fontsize = 'small')
        plt.ylabel('time (s)', fontsize = 'small')
        plt.legend(loc = "upper left", fontsize = 'x-small')
    plt.title("Time consumption by length", fontsize = 'medium')
        
#Usage: 
#pdbList = git.getPDBlist()
#I_time_length_Plot(pdbList)   
        
def Ipert_time_length_plot(PDBlist, degree = 4, perLength = 5):
    '''Time for computation of measure up to specified degree for single 
    perturbation across all names in input list; non/vectorized array/dict 
    versions.
    Input: as output from getPDBlist.'''      
    for d in [0,1]: #to control whether to use array or dict version
        for v in [0,1]: #to control whether to use vectorized version (v = 1) or not
            timeList = []
            lengthList = []
            for l_item in PDBlist:
                f = PDBtop100Folder + l_item[1]
                L = l_item[0]
                #we will perturb the chains so, as they have different lengths, we must use different ranges.
                pRange = range(int(L/2)-2, int(L/2) + 3) #a range of lgth 5 around middle of chain (all have length > 10) 
                #we need also the w-values so we call that computation (vectorized): 
                if d < 1:
                    outDict_org, CaChain, polyCaChain = I(PDBfilename = f, degree =0, inNumpy_b =0 , print_b = 0, vect_b = 1)
                    start = time.time()  
                    #Compute measures on perturbed chain:
                    outDict, perCaChain, perPolyCaChain, segNrRange = Iperturbed(CaChain = CaChain, polyCaChain =polyCaChain, wDict = outDict_org[0], residueNrRange = pRange, degree = degree, full_b = 0, print_b = 0, vect_b = v)       
                    d_t = time.time()-start
                else:
                    outDict_org, CaChain, polyCaChain = I_dict(PDBfilename = f, degree =0, inNumpy_b =0 , print_b = 0, vect_b = 1)                  
                    start = time.time()  
                    outDict, perCaChain, perPolyCaChain, segNrRange = Iperturbed_dict(CaChain = CaChain, polyCaChain =polyCaChain, wDict = outDict_org['w'], residueNrRange = pRange, degree = degree, full_b = 0, print_b = 0, vect_b = v)       
                    d_t = time.time()-start
                timeList.append(d_t)
                lengthList.append(L)
            plt.plot(lengthList, timeList, label = 'dict_b ='+ str(d) +', vect_b ='+ str(v))
            plt.xlabel('length', fontsize = 'small')
            plt.ylabel('time (s)', fontsize = 'small')
            plt.legend(loc = "upper left", fontsize = 'x-small')
            plt.title("Time consumption by length", fontsize = 'medium')
        
#Usage: 
#pdbList = git.getPDBlist()
#Ipert_time_length_Plot(pdbList)          

def Ipert_time_length_perLength_plot(PDBlist, degree = 4, perLengths = [5,10,20]):
    '''Time for single perturbation by length of perturbation and across all names
    in input list. Only vectorized array version.
    Input: as output from getPDBlist.'''    
    for perL in perLengths:
            #generate perturbation sizes:
            perturbation = []
            for i in range(perL):
                perturbation.append((0.1, 0.2, 0.3))
            timeList = []
            lengthList = []
            for l_item in PDBlist:
                f = PDBtop100Folder + l_item[1]
                L = l_item[0]
                #we will perturb the chains so, as they have different lengths, we must use different ranges.
                pRange = range(int(L/2)- int(perL/2), int(L/2) + perL - int(perL/2)) #a range of lgth perL around middle of chain (all have length > 10) 
#                print "perL: %d " % perL
#                print "pRange: %s" % str(pRange) #we need also the w-values so we call that computation (vectorized): 
                outDict_org, CaChain, polyCaChain = I(PDBfilename = f, degree =0, inNumpy_b =0 , print_b = 0, vect_b = 1)
                start = time.time()  
                #Compute measures on perturbed chain:
                outDict, perCaChain, perPolyCaChain, segNrRange = Iperturbed(CaChain = CaChain, polyCaChain =polyCaChain, wDict = outDict_org[0], residueNrRange = pRange, perturbation = perturbation, degree = degree, full_b = 0, print_b = 0, vect_b = 1)       
                d_t = time.time()-start
                timeList.append(d_t)
                lengthList.append(L)
            plt.plot(lengthList, timeList, label = 'perL ='+ str(perL))
            plt.xlabel('length', fontsize = 'small')
            plt.ylabel('time (s)', fontsize = 'small')
            plt.legend(loc = "upper left", fontsize = 'x-small')
            plt.title("Time consumption by length", fontsize = 'medium')

def Ipert_multi_time_nrPerms_plot(PDBnames = PDBlist9, perL = 10, nrPerturbations = [10,100,200]):
    '''Time/perturbation by nr perturbations for names in list; array or dict 
    versions.'''
    i=0
    for name in PDBnames:
        f = PDBtop100Folder + name
        i +=1
        for d in [0,1]:
            #we need also the w-values so we call that computation (vectorized): 
            if d < 1:
                outDict_org, CaChain, polyCaChain = I(PDBfilename = f, degree =0, inNumpy_b =0 , print_b = 0, vect_b = 1)
            else:
                outDict_org, CaChain, polyCaChain = I_dict(PDBfilename = f, degree =0, inNumpy_b =0 , print_b = 0, vect_b = 1)                
            L = len(CaChain)
            #we will perturb the chains so, as they have different lengths, we must use different ranges.
            pRange = range(int(L/2)- int(perL/2), int(L/2) + perL - int(perL/2)) #a range of lgth perL around middle of chain (all have length > 10) 
            timeList = []
            for nrPer in nrPerturbations:
                if d < 1:
                    start = time.time()  
                    #Compute measures on perturbed chain:
                    tot_outDict, perturbedChains = Iperturbed_vect_multi(nrPerturbations = nrPer, CaChain = CaChain, polyCaChain =polyCaChain, wDict = outDict_org[0], perturbationLength =  perL, residueNrRange = pRange, degree = 4, full_b = 0, print_b = 0)       
                    d_t = time.time()-start
                    #clear:
                    tot_outDict=[]
                    perturbedChains=[]
#                        statement = 'git.Iperturbed(CaChain = ' + str(CaChain) + ', polyCaChain =' + str(polyCaChain) + ', wDict = '+ str(outDict[0]) +', residueNrRange =' + str(pRange) +', degree =' + str(degree) + ', vect_b = ' + str(v) +')'
                else:
                    start = time.time()  
                    tot_outDict, perturbedChains = Iperturbed_vect_multi_dict(nrPerturbations = nrPer, CaChain = CaChain, polyCaChain =polyCaChain, wDict = outDict_org['w'], perturbationLength =  perL, residueNrRange = pRange, degree = 4, full_b = 0, print_b = 0)       
                    d_t = time.time()-start
                    #clear
                    tot_outDict = {}
                    perturbedChains = {}
#                    statement = 'git.Iperturbed_dict(CaChain = ' + str(CaChain) + ', polyCaChain =' + str(polyCaChain) + ', wDict = ' + str(outDict['w']) +', residueNrRange =' + pRange + ', degree =' + str(degree) + ', vect_b = ' + str(v) +')'
#                    t_I = timeit.Timer(statement,'import git')
#                    timeList.append(t_I.timeit(number=1))
                timeList.append(float(d_t)/nrPer)
            plt.subplot(3,3,i)
            plt.plot(nrPerturbations, timeList, label = 'dict_b =' + str(d))
#        plt.tight_layout
        if i > 6:
            plt.xlabel('number of perturbations', fontsize = 'small')
        plt.ylabel('time/perturbation (s)', fontsize = 'small')
        plt.legend(loc = "upper left", fontsize = 'x-small')
        plt.title(name + ', L =' + str(L), fontsize = 'medium')
        

def Ipert_multi_vs_repeat_time_nrPerms_plot(PDBnames = PDBlist9, perL = 10, nrPerturbations = [10,100,200]):
    '''Compare time/perturbation by nr perturbations for names in list: repeat 
    of single Iperturbation-fct vs. multi Iperturbation fct. Only vectorized 
    array version.'''
    version = ['Repeat', 'Multi']
    #generate perturbation sizes:
    perturbation = []
    for i in range(perL):
        perturbation.append((0.1, 0.2, 0.3))
    i=0
    for name in PDBnames:
        f = PDBtop100Folder + name
        i +=1
        for v in range(len(version)): #here controls whether to repeat Ipert of use the multi-pert fct
            #we need also the w-values so we call that computation (vectorized): 
            outDict_org, CaChain, polyCaChain = I(PDBfilename = f, degree =0, inNumpy_b =0 , print_b = 0, vect_b = 1)
            L = len(CaChain)
            #we will perturb the chains so, as they have different lengths, we must use different ranges.
            pRange = range(int(L/2)- int(perL/2), int(L/2) + perL - int(perL/2)) #a range of lgth perL around middle of chain (all have length > 10) 
            timeList = []
            for nrPer in nrPerturbations:
                if v < 1: #we repeat call of Ipert nrPer times
                    start = time.time() 
                    for j in range(nrPer):
                        #Compute measures on perturbed chain
                        outDict, perCaChain, perPolyCaChain, segNrRange = Iperturbed(CaChain = CaChain, polyCaChain =polyCaChain, wDict = outDict_org[0], residueNrRange = pRange, perturbation = perturbation, degree = 4, full_b = 0, print_b = 0, vect_b = 1)       
                    d_t = time.time()-start
                    #clear
#                    outDict = []
#                    perCaChain = []
#                    perPolyCaChain = []
#                        statement = 'git.Iperturbed(CaChain = ' + str(CaChain) + ', polyCaChain =' + str(polyCaChain) + ', wDict = '+ str(outDict[0]) +', residueNrRange =' + str(pRange) +', degree =' + str(degree) + ', vect_b = ' + str(v) +')'
                else:
                    start = time.time()  
                    tot_outDict, perturbedChains = Iperturbed_vect_multi(nrPerturbations = nrPer, CaChain = CaChain, polyCaChain =polyCaChain, wDict = outDict_org[0], perturbationLength =  perL, residueNrRange = pRange, degree = 4, full_b = 0, print_b = 0)       
                    d_t = time.time()-start
                    #clear
                    tot_outDict = {}
                    perturbedChains = {}
#                    statement = 'git.Iperturbed_dict(CaChain = ' + str(CaChain) + ', polyCaChain =' + str(polyCaChain) + ', wDict = ' + str(outDict['w']) +', residueNrRange =' + pRange + ', degree =' + str(degree) + ', vect_b = ' + str(v) +')'
#                    t_I = timeit.Timer(statement,'import git')
#                    timeList.append(t_I.timeit(number=1))
                timeList.append(float(d_t)/nrPer)
            plt.subplot(3,3,i+8)
            plt.plot(nrPerturbations, timeList, label = 'version =' + version[v])
#        plt.tight_layout
        if i > 6:
            plt.xlabel('number of perturbations', fontsize = 'small')
        plt.ylabel('time/perturbation (s)', fontsize = 'small')
        plt.legend(loc = "upper left", fontsize = 'x-small')
        plt.title(name + ', L =' + str(L), fontsize = 'medium')
        
######################################################################################
## For comparison vs C-code
######################################################################################

'''Usage:

#import the module
import gitP as git

#####################################
# comparison of results, C vs Python
#####################################

root_db = r'C:\Users\Christian\Bioinformatics\papers\geometric_anomalies_in_folds\c_code\results_C_1st\double_precision'
#or, on the very final version of the C-code for the links-paper (give the id same results; the base algo wasn't changed between
these versions, only the search part):
root_db = r'C:\Users\Christian\Bioinformatics\papers\geometric_anomalies_in_folds\c_code\results_C_1st_final\double_precision'


#Comparison of results, C vs Python (P), for all vertices:
#Get P-results for order = 3 (degree = 6) -- we only consider a few chains (should suffice):
PDBlist, oP3 = git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//SMLX//', fileName = r"C:\Users\Christian\Bioinformatics\masters_study_pgm\projects\knots\code\top100H_source\SMLXNames.txt", degree = 6, all_b =1 )
#Fetch the same for C (global mem alloc and dbl precision):
oC3All_glo = git.readResultsCcode(root_db + '\All_Ivalues_order_3_incl_abs_global_mem_alloc_SMLX.txt', all_b = 1)
#Compare:
git.compareResultsAll_PvsC_Dicts(oP3,oC3All_glo)

#Get P-results for order = 2 (degree = 4) -- we only consider a few chains (should suffice):
PDBlist, oP2 = git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//SMLX//', fileName = r"C:\Users\Christian\Bioinformatics\masters_study_pgm\projects\knots\code\top100H_source\SMLXNames.txt", degree = 4, all_b =1 )
#Fetch the same for C (global mem alloc and dbl precision):
oC2All_glo = git.readResultsCcode(root_db + '\All_Ivalues_order_2_incl_abs_global_mem_alloc_SMLX.txt', all_b = 1)
#Compare:
git.compareResultsAll_PvsC_Dicts(oP2,oC2All_glo)

#Get P-results for order = 1 (degree = 2) -- we only consider a few chains (should suffice):
PDBlist, oP1 = git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//SMLX//', fileName = r"C:\Users\Christian\Bioinformatics\masters_study_pgm\projects\knots\code\top100H_source\SMLXNames.txt", degree = 2, all_b =1 )
#Fetch the same for C (global mem alloc and dbl precision):
oC1All_glo = git.readResultsCcode(root_db + '\All_Ivalues_order_1_incl_abs_global_mem_alloc_SMLX.txt', all_b = 1)
#Compare:
git.compareResultsAll_PvsC_Dicts(oP1,oC1All_glo)

#Get P-results for order = 0 (degree = 0) -- we only consider a few chains (should suffice):
PDBlist, oP0 = git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//SMLX//', fileName = r"C:\Users\Christian\Bioinformatics\masters_study_pgm\projects\knots\code\top100H_source\SMLXNames.txt", degree = 0, all_b =1 )
#Fetch the same for C (global mem alloc and dbl precision):
oC0All_glo = git.readResultsCcode(root_db + '\All_Ivalues_order_0_incl_abs_global_mem_alloc_SMLX.txt', all_b = 1)
#Compare:
git.compareResultsAll_PvsC_Dicts(oP0,oC0All_glo)


#Results: all diffs are < 0.1% and most are much lower. Here for order 3 (SMLX set):

"
>>> git.compareResultsAll_PvsC_Dicts(oP3,oC3All_glo)
PDB: C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//SMLX//bio1rpoH
Length:61
Chain id: A
No of invariants considered: 20
Max rel diff (pct) final corner: 0.013466 attained for invariant: w
Max meanDiff-over-mean ratio is (pct): 0.000388, attained at: w
Max stdDiff-over-mean ratio is (pct): 0.000259, attained at: w
PDB: C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//SMLX//bio1rpoH
Length:61
Chain id: B
No of invariants considered: 20
Max rel diff (pct) final corner: 0.015626 attained for invariant: w
Max meanDiff-over-mean ratio is (pct): 0.000398, attained at: w
Max stdDiff-over-mean ratio is (pct): 0.000263, attained at: w
PDB: C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//SMLX//1edmBH
Length:39
Chain id: B
No of invariants considered: 20
Max rel diff (pct) final corner: 0.060324 attained for invariant: w
Max meanDiff-over-mean ratio is (pct): 0.000449, attained at: w
Max stdDiff-over-mean ratio is (pct): 0.000313, attained at: w
PDB: C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//SMLX//1benABH
Length:21
Chain id: A
No of invariants considered: 20
Max rel diff (pct) final corner: 0.000975 attained for invariant: w
Max meanDiff-over-mean ratio is (pct): 0.000234, attained at: I152634
Max stdDiff-over-mean ratio is (pct): 0.000289, attained at: I162534
PDB: C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//SMLX//1benABH
Length:30
Chain id: B
No of invariants considered: 20
Max rel diff (pct) final corner: 0.001760 attained for invariant: w
Max meanDiff-over-mean ratio is (pct): 0.000339, attained at: w
Max stdDiff-over-mean ratio is (pct): 0.000254, attained at: w
PDB: C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//SMLX//1adsH
Length:315
Chain id: >
No of invariants considered: 20
Max rel diff (pct) final corner: 0.005365 attained for invariant: w
Max meanDiff-over-mean ratio is (pct): 0.001278, attained at: w
Max stdDiff-over-mean ratio is (pct): 0.000771, attained at: w
PDB: C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//SMLX//2mcmH
Length:112
Chain id: >
No of invariants considered: 20
Max rel diff (pct) final corner: 0.012209 attained for invariant: w
Max meanDiff-over-mean ratio is (pct): 0.000710, attained at: w
Max stdDiff-over-mean ratio is (pct): 0.000438, attained at: w
No of structures considered: 7
>>> 
"


##Comparison of results, C vs Python (P), invariant-values (final corner); top100 set

#get Python results:
PDBlist, oP3_100 = git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//top100H//', fileName = r"C:\Users\Christian\Bioinformatics\masters_study_pgm\projects\knots\code\top100H_source\top100names.txt", degree = 6, all_b =0)
#fetch C-results (this time load without perturbation nr):
oC3_glo = git.readResultsCcode(root_db + '\Ivalues_order_3_incl_abs_global_mem_alloc_top100.txt', all_b = 0)
#Compare
git.compareResultsFinal_PvsC_Dicts(oP3_100, oC3_glo)

#Same, order 2:
PDBlist, oP2_100 = git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//top100H//', fileName = r"C:\Users\Christian\Bioinformatics\masters_study_pgm\projects\knots\code\top100H_source\top100names.txt", degree = 4, all_b =0)
#fetch C-results (this time load without perturbation nr):
oC2_glo = git.readResultsCcode(root_db + '\Ivalues_order_2_incl_abs_global_mem_alloc_top100.txt', all_b = 0)
#Compare
git.compareResultsFinal_PvsC_Dicts(oP2_100, oC2_glo)


#Same, order 1:
PDBlist, oP1_100 = git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//top100H//', fileName = r"C:\Users\Christian\Bioinformatics\masters_study_pgm\projects\knots\code\top100H_source\top100names.txt", degree = 2, all_b =0)
#fetch C-results (this time load without perturbation nr):
oC1_glo = git.readResultsCcode(root_db + '\Ivalues_order_1_incl_abs_global_mem_alloc_top100.txt', all_b = 0)
#Compare
git.compareResultsFinal_PvsC_Dicts(oP1_100, oC1_glo)


Res: for the vast majority diff's are below 0.001 pct (and a notch lower too) in order 3 (and in other orders too); 
the size of the diff's decrease by an order of 10 or so for every step down in order. In a handful of cases the diff 
is higher than that (order 3 fig's but others too), but only two very significant cases remain (and shared in all
orders): 1ctj and 1ifc; for these two high diff's even at high (so-so) abs values of the measures (at which the max 
rel diff is attained) are found (in fact, in order 3 the values were around 382%/0.7 and 12%/3.7 quoted as 
rel-diff-pct/abs-value); in other cases the underlying value is quite low as is the diff. (The former case of 1smd 
for which the P-code loaded 496 C-alpha's whereas the C-code loaded 495 is now corrected; it was due to an HETATOM 
being loaded in the P-code.) The cases for which the C-dict does not contain results for the given chain are all 
due to the length of the chain being lower than 6 (in the C code we exclude such cases from the computation).  


######################################################################################################
#FROM HERE AND DOWN TO LINK-SEARCHING NOT DONE# (since these checks are relevan for parts not used in the link-paper)
######################################################################################################




#Aside: comparison of local/global mem alloc in plainC and cuda-C versions (even though the cuda version does not use cuda; there could be
differences since the cuda version is a .cu file and compiles in c++, while the plainC code is id to the cuda version but with the 
gpu stuff removed; this must compile in C though!). (Comparison for perturbation code in section on that right below):
order = 3
oCAll_glo = git.readResultsCcode(root_db + '\All_Ivalues_order_' + str(order) + '_incl_abs_global_mem_alloc.txt', all_b = 1, pert_b =1)
oCAll_glo_cudaC = git.readResultsCcode(root_db + '\All_Ivalues_order_' + str(order) + '_incl_abs_global_mem_alloc_cudaC.txt', all_b = 1, pert_b =1)
oCAll_loc = git.readResultsCcode(root_db + '\All_Ivalues_order_' + str(order) + '_incl_abs_local_mem_alloc.txt', all_b = 1, pert_b =1)
oCAll_loc_cudaC = git.readResultsCcode(root_db + '\All_Ivalues_order_' + str(order) + '_incl_abs_local_mem_alloc_cudaC.txt', all_b = 1, pert_b =1)

#Compare
git.compareCresultsAll(oCAll_glo, oCAll_glo_cudaC)
#results. all zero at first 6 dec places
git.compareCresultsAll(oCAll_loc, oCAll_loc_cudaC)
#results. all zero at first 6 dec places (also in double precision)

 
#############################################
#Statements for plotting time consumption:
#############################################
 
#first get the P results (here order 3); if you want to use the global allocation version 
#set global_b = 1 else set it to 0:
global_b = 1
PDBlist, oP3 = git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//projects//knots//code//top100H//', fileName = r"C:\Users\Christian\Bioinformatics\projects\knots\code\top100H_source\top100Names.txt", degree = 6, all_b =0, global_b = global_b, repeats = 20)
#takes long time to run, so dump the results:
import cPickle as pickle
pickle.dump(oP3, open( "oP3_top100_global_1.p", "wb" ) )
to get P-results in order 2 instead (here with 20 repeats):
PDBlist, oP2 = git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//projects//knots//code//top100H//', fileName = r"C:\Users\Christian\Bioinformatics\projects\knots\code\top100H_source\top100Names.txt", degree = 4, all_b =0, global_b = global_b, repeats = 20) 
#dump:
pickle.dump(oP2, open( "oP2_top100_global.p", "wb" ) )
  
#Then load the C results:
order = 2
capRatio = 500
floorCtime = 0.001
type = 'global' #global is most fairly compared to P-results computed with global_b = 1 
prec = 'single_precision'
path = root + '\\' 
path = path + '\Ivalues_order_'+ str(order) + '_excl_abs_'+type+'_mem_alloc_top100_20repeats.txt'
oC = git.readResultsCcode(path, all_b = 0)

#Run the plot fct -- here for the cpuTime, ie including time for loading the structure:
#but first load results dumped previously:
oP2_loc= pickle.load( open( "oP2_top100_local.p", "rb" ) )
oP3_loc = pickle.load( open( "oP3_top100_local.p", "rb" ) )
oP2_glo = pickle.load( open( "oP2_top100_global.p", "rb" ) )
oP3_glo = pickle.load( open( "oP3_to p100_global.p", "rb" ) )
oP2_20 = pickle.load( open( "oP2_top100_20repeats_1.p", "rb" ) )
oP3_20 = pickle.load( open( "oP3_top100_20repeats_1.p", "rb" ) )
and get PDB-list:
PDBlist, oP0= git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//projects//knots//code//top100H//', fileName = r"C:\Users\Christian\Bioinformatics\projects\knots\code\top100H_source\top100Names.txt", degree = 0, all_b =0, global_b = 0)
#Now plot:
as_subplots_b = 1
figSize = (14,7)
git.plot_PvsC_time_length(PDBlist,oP2_20,oC, timeType= 'cpuTime', titleTimePlot = "Time cons. (total) by length, globals, order 2 (P) order 3 (C)", titlePtoCratio = "Ratio of time consumption (total) by length, Python over C", capRatio = capRatio, floorCtime = floorCtime, as_subplots_b = as_subplots_b, figSize = figSize)

#... while for the computation time only:
git.plot_PvsC_time_length(PDBlist,oP2_glo,oC, timeType= 'aggrTime', titleTimePlot = "Time cons. (comp) by length, globals, order 2 (P) order 3 (C)", titlePtoCratio = "Ratio of time consumption (comp) by length, Python over C", capRatio = capRatio, floorCtime = floorCtime,  as_subplots_b = as_subplots_b, figSize = figSize)
 
#Similarly for the other orders and for using the local memory allocation version of the C-code.


#To use the timeit-function instead we can run a full order 2 or order 3 
#run and sum up the consumption. This can be carried out by a function 
#above, slightly modified (to avoid plotting out):
#first get the PDBlist:
PDBlist, oP0= git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//projects//knots//code//top100H//', fileName = r"C:\Users\Christian\Bioinformatics\projects\knots\code\top100H_source\top100Names.txt", degree = 0, all_b =0, global_b = 0)
#the call the function:
time = git.I_time_length_plot_mod(PDBlist, degree = 4)
print time

## Compare C vs C

#time consumption by order:
#for single:
path = r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C_5th\single_precision' 
#for double:
path = r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C_5th\double_precision' 
#run:
git.plotCTimeByOrders(pathCresultFiles = path, orderList = [0,1,2,3], PDBlist = PDBlist, fileId =  r'_incl_abs_global_mem_alloc_top100.txt', timeType = 'aggrTime', title = 'C: Time consumption (comp) by length at varying order (Dbl)')
git.plotCTimeByOrders(pathCresultFiles = path, orderList = [0,1,2], PDBlist = PDBlist, fileId =  r'_incl_abs_global_mem_alloc_top100.txt', timeType = 'aggrTime', title = 'C: Time consumption (comp) by length at varying order (Dbl)')


## Aside: Compare P local vs P global:
#given that we have gotten the results for the global version above we need only the local version:
global_b = 0
#order 3:
PDBlist, oP3_loc = git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//projects//knots//code//top100H//', fileName = r"C:\Users\Christian\Bioinformatics\projects\knots\code\top100H_source\top100Names.txt", degree = 6, all_b =0, global_b = global_b)
#to get P-results in other order 2 instead:
PDBlist, oP2_loc = git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//projects//knots//code//top100H//', fileName = r"C:\Users\Christian\Bioinformatics\projects\knots\code\top100H_source\top100Names.txt", degree = 4, all_b =0, global_b = global_b)
#... or load all from dumped versions:
import cPickle as pickle
oP2_loc= pickle.load( open( "oP2_top100_local.p", "rb" ) )
oP3_loc = pickle.load( open( "oP3_top100_local.p", "rb" ) )
oP2_glo = pickle.load( open( "oP2_top100_global.p", "rb" ) )
oP3_glo = pickle.load( open( "oP3_top100_global.p", "rb" ) )

#Abuse the function by calling it on two P-dicts:
as_subplots_b = 1
figSize = (14,7)
#total time:
git.plot_PvsC_time_length(PDBlist,oP3_glo,oP3_loc, timeType= 'aggrTime', titleTimePlot = "Time cons. (comp) by length, P-global vs P-local, order 3", titlePtoCratio = "Ratio of time cons. (comp) by length, P-global over P-local", 
                          labelC = 'Python local', labelP = 'Python global', labelPoverCratio = 'Python global over Python local', capRatio = capRatio, floorCtime = floorCtime, as_subplots_b = as_subplots_b, figSize = figSize)

#comp time:
git.plot_PvsC_time_length(PDBlist,oP2_glo,oP2_loc, timeType= 'aggrTime', titleTimePlot = "Time cons. (comp) by length, P-global vs P-local, order 2", titlePtoCratio = "Ratio of time cons. (total) by length, P-global over P-local", 
                          labelC = 'Python local', labelP = 'Python global', labelPoverCratio = 'Python global over Python local', capRatio = capRatio, floorCtime = floorCtime, as_subplots_b = as_subplots_b, figSize = figSize)


######################################
#Comparison of C-code results:
######################################

##Perturbations

#Checks:
#1. local_/global_alloc vs local/global_pert, non-gpu, unpert. All orders, SML set.

#2. local/global_pert, compare all in all modes; 0-pert; order 0; top100 set;  

#3. local/global_pert, compare all in all modes; 5 random perts; All; SML set
#4. For info: rel. diff's; top100; order 3; 0-pert and fixed lgth 10 pert.
#5. Connecting the 0-pert to the unpert: compare these 
- a) for local-gpu, non-gpu or for global-gpu, non-gpu. (by #1 and #2, local and global should be id)
- b) for local-gpu, gpu or for global-gpu, gpu or gpu2 
But rendered obsolete by #2 and first point (#5a): e.g. local-gpu, gpu unpert vs 0-pert will be 
close if local-gpu, gpu on 0-pert is close to local-gpu, non-gpu 0-pert (ie.#2)
and, by #5a, local-gpu, non-gpu 0-pert is close to local-gpu, non-gpu unpert.  


root = r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C_5th\single_precision'
root_pert = r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C_5th\Perturbed\single_precision'

#1. global-gpu vs main_global_alloc on unperturbed chain:
order = 3
oCAll_glo = git.readResultsCcode(root + '\All_Ivalues_order_' + str(order) + '_incl_abs_global_mem_alloc.txt', all_b = 1, pert_b = 1)
oCAll_unpert_glo_gpu = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
oCAll_unpert_loc_gpu = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_local_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
#Compare
git.compareCresultsAll(oCAll_glo, oCAll_unpert_glo_gpu)
git.compareCresultsAll(oCAll_glo, oCAll_unpert_loc_gpu)
#Results: all are id. The same holds in all orders.


#2. local_gpu and global_gpu versions give similar results in 0-perturbation case:
#load:
oC0All_0pert_loc_top100 = git.readResultsCcode(root_pert + r'\top100\Pert_All_Ivalues_order_0_incl_abs_local_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_loc_gpu_top100 = git.readResultsCcode(root_pert + r'\top100\Pert_All_Ivalues_order_0_incl_abs_local_gpu_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_top100 = git.readResultsCcode(root_pert + r'\top100\Pert_All_Ivalues_order_0_incl_abs_global_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_gpu_top100 = git.readResultsCcode(root_pert + r'\top100\Pert_All_Ivalues_order_0_incl_abs_global_gpu_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_gpu2_top100 = git.readResultsCcode(root_pert + r'\top100\Pert_All_Ivalues_order_0_incl_abs_global_gpu2_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
#compare
git.compareCresultsAll(oC0All_0pert_loc_top100, oC0All_0pert_glo_top100, plot_b =1, title = 'local-gpu vs global-gpu, non-gpu, 0-pert, order 0 ,top100', labelDTS = 'max MeanDiffToMeanValue', labelNTS = 'max NoiseToMeanValue')
git.compareCresultsAll(oC0All_0pert_loc_top100, oC0All_0pert_loc_gpu_top100, plot_b =1, title = 'local-gpu vs local-gpu gpu, 0-pert, order 0, top100')
git.compareCresultsAll(oC0All_0pert_glo_top100, oC0All_0pert_glo_gpu_top100, plot_b =1, title = 'global-gpu vs global-gpu gpu, 0-pert, order 0, top100')
git.compareCresultsAll(oC0All_0pert_glo_top100, oC0All_0pert_glo_gpu2_top100, plot_b =1, title = 'global-gpu vs global-gpu gpu2, 0-pert, order 0, top100')
git.compareCresultsAll(oC0All_0pert_glo_gpu_top100, oC0All_0pert_glo_gpu2_top100, plot_b =1, title = 'global-gpu gpu vs global-gpu gpu2, 0-pert, order 0, top100')
#Results: see the plots!

#Bechmark comparison vs double precision: see "double precision" below.

#3. local_gpu and global_gpu versions give similar results on unperturbed and 5 random perts:
#unperturbed
#local non-gpu vs gpu:
#load results. "All" type (compare results over the whole simplex):
order = 0
oCAll_unpert_loc_SML = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_local_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
oCAll_unpert_loc_gpu_SML = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_local_gpu_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
#run comparison
git.compareCresultsAll(oCAll_unpert_loc_SML, oCAll_unpert_loc_gpu_SML)
#Results: diffs bounded by 0.1/0.35 (1ads)

#global non-gpu vs local and vs gpu vs gpu2:
order = 0
oCAll_unpert_glo_SML = git.readResultsCcode(root_pert + '\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
oCAll_unpert_glo_gpu_SML = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_gpu_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
oCAll_unpert_glo_gpu2_SML = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_gpu2_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
#run comparison
#local vs global (order 0 suffices):
git.compareCresultsAll(oCAll_unpert_loc_SML,oCAll_unpert_glo_SML)
git.compareCresultsAll(oCAll_unpert_loc_gpu_SML,oCAll_unpert_glo_gpu_SML)
#Results: all 0 (6 dec places)
#global vs global
git.compareCresultsAll(oCAll_unpert_glo_SML,oCAll_unpert_glo_gpu_SML)
git.compareCresultsAll(oCAll_unpert_glo_gpu_SML,oCAll_unpert_glo_gpu2_SML)
#Results: glo non-gpu vs glo gpu id to ditto for loc. 
#glo-gpu vs glo-gpu2 differ much less than that (bounded by 0.0008/0.021 pct, 1ads) 

#5 perts:
#local non-gpu vs gpu:
#load results. "All" type (compare results over the whole simplex):
order = 3
oCAll_5perts_loc_SML = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_local_lgth10_5_perts_rnd_SML.txt', all_b = 1, pert_b =1)
oCAll_5perts_loc_gpu_SML = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_local_gpu_lgth10_5_perts_rnd_SML.txt', all_b = 1, pert_b =1)
#run comparison
git.compareCresultsAll(oCAll_5perts_loc_SML, oCAll_5perts_loc_gpu_SML)


#global non-gpu vs local and vs gpu vs gpu2:
order = 0
oCAll_5perts_glo_SML = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_lgth10_5_perts_rnd_SML.txt', all_b = 1, pert_b =1)
oCAll_5perts_glo_gpu_SML = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_gpu_lgth10_5_perts_rnd_SML.txt', all_b = 1, pert_b =1)
oCAll_5perts_glo_gpu2_SML = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_gpu2_lgth10_5_perts_rnd_SML.txt', all_b = 1, pert_b =1)
#run comparison
#local vs global (order 0 suffices):
git.compareCresultsAll(oCAll_5perts_loc_SML,oCAll_5perts_glo_SML)
git.compareCresultsAll(oCAll_5perts_loc_gpu_SML,oCAll_5perts_glo_gpu_SML)
#Results: 0 diffs (6 dec places)
#global vs global
git.compareCresultsAll(oCAll_5perts_glo_SML,oCAll_5perts_glo_gpu_SML)
git.compareCresultsAll(oCAll_5perts_glo_gpu_SML,oCAll_5perts_glo_gpu2_SML)
#gpu vs gpu2: bounded by 0.009/0.1 pct (1ads but 1edm quite alike); varying quite a lot

#4. Compare results at the "final corner": length 0-perturbations at random positions, top100 set (final values 
#and w at all vertices); we compare also the diff's at order 3 to to the diff's seen in order 0:
#local, w-terms (order 0)
oC0All_0pert_loc_top100 = git.readResultsCcode(root_pert + r'\top100\Pert_All_Ivalues_order_0_incl_abs_local_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_loc_gpu_top100 = git.readResultsCcode(root_pert + r'\top100\Pert_All_Ivalues_order_0_incl_abs_local_gpu_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
#final values
oC3_0pert_loc_top100 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_local_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
oC3_0pert_loc_gpu_top100 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_local_gpu_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
#Compare:
#order 0 diff's (this is the same as #2 local vs local_gpu):
git.compareCresultsAll(oC0All_0pert_loc_top100, oC0All_0pert_loc_gpu_top100, plot_b = 1, title = 'local-pert vs local-pert gpu, 0-pert, order 0, top100')
#all orders, final values:
sortedFileList = git.compareCresultsFinal(oC3_0pert_loc_top100, oC3_0pert_loc_gpu_top100, title = 'local-pert vs local-pert gpu, 0-pert, rel diffs in final values, order 3, top100')

#global
#w-terms (order 0)
oC0All_0pert_glo_top100 = git.readResultsCcode(root_pert + r'\top100\Pert_All_Ivalues_order_0_incl_abs_global_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_gpu_top100 = git.readResultsCcode(root_pert + r'\top100\Pert_All_Ivalues_order_0_incl_abs_global_gpu_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_gpu2_top100 = git.readResultsCcode(root_pert + r'\top100\Pert_All_Ivalues_order_0_incl_abs_global_gpu2_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
#final values
oC3_0pert_glo_top100 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu_top100 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu2_top100 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu2_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)

#Compare:
#order 0 diff's: id to #2 on the global-gpu versions, e.g.
git.compareCresultsAll(oC0All_0pert_glo_top100, oC0All_0pert_glo_gpu_top100, plot_b = 1, title = 'global-pert vs gobal-pert gpu, 0-pert, order 0, top100')
#all orders, final values:
sortedFileList = git.compareCresultsFinal(oC3_0pert_glo_top100, oC3_0pert_glo_gpu_top100, title = 'global-pert vs global-pert gpu, 0-pert, rel diffs in final values, order 3, top100')
sortedFileList = git.compareCresultsFinal(oC3_0pert_glo_top100, oC3_0pert_glo_gpu2_top100, title = 'global-pert vs global-pert gpu2, 0-pert, rel diffs in final values, order 3, top100')

#on top8000:
#local
oC3_0pert_loc_top8000 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_local_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_loc_gpu_top8000 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_local_gpu_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
#all orders, final values:
sortedFileList = git.compareCresultsFinal(oC3_0pert_loc_top8000, oC3_0pert_loc_gpu_top8000, title = 'local-pert vs local-pert gpu, 0-pert, rel diffs in final values, order 3, top8000')

#global
oC3_0pert_glo_top8000 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu_top8000 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu2_top8000 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu2_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
#all orders, final values:
sortedFileList = git.compareCresultsFinal(oC3_0pert_glo_top8000, oC3_0pert_glo_gpu_top8000, title = 'global-pert vs global-pert gpu, 0-pert, rel diffs in final values, order 3, top8000')
sortedFileList = git.compareCresultsFinal(oC3_0pert_glo_top8000, oC3_0pert_glo_gpu2_top8000, title = 'global-pert vs global-pert gpu2, 0-pert, rel diffs in final values, order 3, top8000')

#Showing instead the differences per invariant:
#top100
means, stdDevs = git.compareCresultsFinal_2(oC3_0pert_loc_top100, oC3_0pert_loc_gpu_top100, titleDistr = 'local-gpu, non-gpu, top100', title = 'local-gpu, non-gpu vs gpu, top100',titleRatios = 'local-gpu, non-gpu vs gpu, top100') 
medians, means, stdDevs, lP, uP = git.compareCresultsFinal_3(oC3_0pert_loc_top100, oC3_0pert_loc_gpu_top100, titleDistr = 'local-gpu, non-gpu vs gpu, top100', title = 'Statistics of relative diffs, local-gpu, non-gpu vs gpu, top100')
#top8000; better since more representative:
#final values
oC3_0pert_loc_top8000 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_local_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_loc_gpu_top8000 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_local_gpu_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
#comparisons:
means, stdDevs = git.compareCresultsFinal_2(oC3_0pert_loc_top8000, oC3_0pert_loc_gpu_top8000, titleDistr = 'local-gpu, non-gpu, top8000',title = 'local-gpu, non-gpu vs gpu, top8000',titleRatios = 'local-gpu, non-gpu vs gpu, top8000') 
medians, means, stdDevs, lP, uP = git.compareCresultsFinal_3(oC3_0pert_loc_top8000, oC3_0pert_loc_gpu_top8000, titleDistr = 'local-gpu, non-gpu vs gpu, top8000',title = 'Stats of rel diffs, local-gpu, non-gpu vs gpu, top8000 (Sgl)', figSize = (10,8))
#global version:
oC3_0pert_glo_top8000 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu_top8000 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu2_top8000 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu2_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
#comp's, non-gpu vs gpu2 (should suffice):
means, stdDevs = git.compareCresultsFinal_2(oC3_0pert_glo_top8000, oC3_0pert_glo_gpu2_top8000, titleDistr = 'global-gpu, non-gpu, top8000', title = 'global-gpu, non-gpu vs gpu2, top8000',titleRatios = 'global-gpu, non-gpu vs gpu2, top8000') 
medians, means, stdDevs, lP, uP = git.compareCresultsFinal_3(oC3_0pert_glo_top8000, oC3_0pert_glo_gpu2_top8000, titleDistr = 'global-gpu, non-gpu vs gpu2, top8000', title = 'Stats of rel diffs, global-gpu, non-gpu vs gpu2, top8000 (Sgl)', figSize = (10,8))

#plots of normalized invariants (acc to Røgen 2005) and diff's:
#local:
git.compareNormalizedInv(oC3_0pert_loc_top8000 , oC3_0pert_loc_gpu_top8000, figSize = (14,10), titleDiffs = 'Differences vs values (Sgl)')
#global:
git.compareNormalizedInv(oC3_0pert_glo_top8000 , oC3_0pert_glo_gpu_top8000, figSize = (14,10), titleDiffs = 'Differences vs values (Sgl)')
git.compareNormalizedInv(oC3_0pert_glo_top8000 , oC3_0pert_glo_gpu2_top8000, figSize = (14,10), titleDiffs = 'Differences vs values (Sgl)') 


#4 on run using 4th version (sign change in w-fct):
#Compare results at the "final corner": length 0-perturbations at random positions, top100 set (final values 
#and w at all vertices); we compare also the diff's at order 3 to to the diff's seen in order 0:
#local, w-terms (order 0)
oC0All_0pert_loc_top100 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\top100\Pert_All_Ivalues_order_0_incl_abs_local_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_loc_gpu_top100 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\top100\Pert_All_Ivalues_order_0_incl_abs_local_gpu_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
#final values
oC3_0pert_loc_top100 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_local_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
oC3_0pert_loc_gpu_top100 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_local_gpu_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
#Compare:
#order 0 diff's (this is the same as #2 local vs local_gpu):
git.compareCresultsAll(oC0All_0pert_loc_top100, oC0All_0pert_loc_gpu_top100, plot_b = 1, title = 'local-pert vs local-pert gpu, 0-pert, order 0, top100')
#all orders, final values:
sortedFileList = git.compareCresultsFinal(oC3_0pert_loc_top100, oC3_0pert_loc_gpu_top100, title = 'local-pert vs local-pert gpu, 0-pert, rel diffs in final values, order 3, top100')

#global
#w-terms (order 0)
oC0All_0pert_glo_top100 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\top100\Pert_All_Ivalues_order_0_incl_abs_global_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_gpu_top100 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\top100\Pert_All_Ivalues_order_0_incl_abs_global_gpu_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_gpu2_top100 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\top100\Pert_All_Ivalues_order_0_incl_abs_global_gpu2_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
#final values
oC3_0pert_glo_top100 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu_top100 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu2_top100 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu2_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)

#Compare:
#order 0 diff's: id to #2 on the global-gpu versions, e.g.
git.compareCresultsAll(oC0All_0pert_glo_top100, oC0All_0pert_glo_gpu_top100, plot_b = 1, title = 'global-pert vs gobal-pert gpu, 0-pert, order 0, top100')
#all orders, final values:
sortedFileList = git.compareCresultsFinal(oC3_0pert_glo_top100, oC3_0pert_glo_gpu_top100, title = 'global-pert vs global-pert gpu, 0-pert, rel diffs in final values, order 3, top100')
sortedFileList = git.compareCresultsFinal(oC3_0pert_glo_top100, oC3_0pert_glo_gpu2_top100, title = 'global-pert vs global-pert gpu2, 0-pert, rel diffs in final values, order 3, top100')

#on top8000:
oC3_0pert_glo_top8000 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu_top8000 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu2_top8000 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu2_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
#all orders, final values:
sortedFileList = git.compareCresultsFinal(oC3_0pert_glo_top8000, oC3_0pert_glo_gpu_top8000, title = 'global-pert vs global-pert gpu, 0-pert, rel diffs in final values, order 3, top8000')
sortedFileList = git.compareCresultsFinal(oC3_0pert_glo_top8000, oC3_0pert_glo_gpu2_top8000, title = 'global-pert vs global-pert gpu2, 0-pert, rel diffs in final values, order 3, top8000')


#Showing instead the differences per invariant:
#top100
means, stdDevs = git.compareCresultsFinal_2(oC3_0pert_loc_top100, oC3_0pert_loc_gpu_top100, titleDistr = 'local-gpu, non-gpu, top100', title = 'local-gpu, non-gpu vs gpu, top100',titleRatios = 'local-gpu, non-gpu vs gpu, top100') 
medians, means, stdDevs, lP, uP = git.compareCresultsFinal_3(oC3_0pert_loc_top100, oC3_0pert_loc_gpu_top100, titleDistr = 'local-gpu, non-gpu vs gpu, top100', title = 'Statistics of relative diffs, local-gpu, non-gpu vs gpu, top100')
#top8000; better since more representative:
#final values
oC3_0pert_loc_top8000 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_local_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_loc_gpu_top8000 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_local_gpu_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
#comparisons:
means, stdDevs = git.compareCresultsFinal_2(oC3_0pert_loc_top8000, oC3_0pert_loc_gpu_top8000, titleDistr = 'local-gpu, non-gpu, top8000',title = 'local-gpu, non-gpu vs gpu, top8000',titleRatios = 'local-gpu, non-gpu vs gpu, top8000') 
medians, means, stdDevs, lP, uP = git.compareCresultsFinal_3(oC3_0pert_loc_top8000, oC3_0pert_loc_gpu_top8000, titleDistr = 'local-gpu, non-gpu vs gpu, top8000',title = 'Stats of rel diffs, local-gpu, non-gpu vs gpu, top8000 (Sgl)')
#global version:
oC3_0pert_glo_top8000 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu_top8000 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu2_top8000 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu2_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
#comp's, non-gpu vs gpu2 (should suffice):
means, stdDevs = git.compareCresultsFinal_2(oC3_0pert_glo_top8000, oC3_0pert_glo_gpu2_top8000, titleDistr = 'global-gpu, non-gpu, top8000', title = 'global-gpu, non-gpu vs gpu2, top8000',titleRatios = 'global-gpu, non-gpu vs gpu2, top8000') 
medians, means, stdDevs, lP, uP = git.compareCresultsFinal_3(oC3_0pert_glo_top8000, oC3_0pert_glo_gpu_top8000, titleDistr = 'global-gpu, non-gpu vs gpu, top8000', title = 'Stats of rel diffs, global-gpu, non-gpu vs gpu, top8000 (Sgl)')

#plots of normalized invariants (acc to Røgen 2005) and diff's:
#local:
git.compareNormalizedInv(oC3_0pert_loc_top8000 , oC3_0pert_loc_gpu_top8000, figSize = (14,10), titleDiffs = 'Differences vs values (Sgl)')
#global:
git.compareNormalizedInv(oC3_0pert_glo_top8000 , oC3_0pert_glo_gpu_top8000, figSize = (14,10), titleDiffs = 'Differences vs values (Sgl)')
git.compareNormalizedInv(oC3_0pert_glo_top8000 , oC3_0pert_glo_gpu2_top8000, figSize = (14,10), titleDiffs = 'Differences vs values (Sgl)') 


#5. By the reasoning above we only make 5a: compare for local-gpu and global-gpu in non-gpu mode, the 0-pert and the unpert:
order = 0
oCAll_unpert_loc_gpu = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_local_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
oCAll_unpert_glo_gpu = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
oCAll_0pert_loc_gpu = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_local_lgth10_0pert_rnd_SML.txt', all_b = 1, pert_b =1)
oCAll_0pert_glo_gpu = git.readResultsCcode(root_pert + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_lgth10_0pert_rnd_SML.txt', all_b = 1, pert_b =1)
#Compare
git.compareCresultsAll(oCAll_0pert_glo_gpu, oCAll_unpert_glo_gpu)
git.compareCresultsAll(oCAll_0pert_loc_gpu, oCAll_unpert_loc_gpu)
#Results: all diffs are zero at first 6 dec places
 


## Double precision:

root_db = r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C_5th\double_precision'
root_pert_db = r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C_5th\Perturbed\double_precision'

#1. no need to redo (#1 in single precision and #2 below should do).

#2. local_gpu and global_gpu versions give similar results in unperturbed case:
#load:
oC0All_0pert_loc_top100_db = git.readResultsCcode(root_pert_db + r'\top100\Pert_All_Ivalues_order_0_incl_abs_local_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_loc_gpu_top100_db = git.readResultsCcode(root_pert_db + r'\top100\Pert_All_Ivalues_order_0_incl_abs_local_gpu_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_top100_db = git.readResultsCcode(root_pert_db + r'\top100\Pert_All_Ivalues_order_0_incl_abs_global_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_gpu_top100_db = git.readResultsCcode(root_pert_db + r'\top100\Pert_All_Ivalues_order_0_incl_abs_global_gpu_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_gpu2_top100_db = git.readResultsCcode(root_pert_db + r'\top100\Pert_All_Ivalues_order_0_incl_abs_global_gpu2_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
#compare
git.compareCresultsAll(oC0All_0pert_loc_top100_db, oC0All_0pert_glo_top100_db, plot_b =1, title = 'local-gpu vs global-gpu, non-gpu, 0-pert, order 0, top100 (Dbl)')
git.compareCresultsAll(oC0All_0pert_loc_top100_db, oC0All_0pert_loc_gpu_top100_db, plot_b =1, title = 'local-gpu vs local-gpu gpu, 0-pert, order 0, top100 (Dbl)')
git.compareCresultsAll(oC0All_0pert_glo_top100_db, oC0All_0pert_glo_gpu_top100_db, plot_b =1, title = 'global-gpu vs global-gpu gpu, 0-pert, order 0, top100 (Dbl)')
git.compareCresultsAll(oC0All_0pert_glo_top100_db, oC0All_0pert_glo_gpu2_top100_db, plot_b =1, title = 'global-gpu vs global-gpu gpu2, 0-pert, order 0, top100 (Dbl)')
git.compareCresultsAll(oC0All_0pert_glo_gpu_top100_db, oC0All_0pert_glo_gpu2_top100_db, plot_b =1, title = 'global-gpu gpu vs global-gpu gpu2, 0-pert, order 0, top100 (Dbl)')
#Results: see the plots!


#benchmark comparison: vs main_global_alloc, order 0:
#load
oC0All_glo_top100_db = git.readResultsCcode(root_db + r'\All_Ivalues_order_0_incl_abs_global_mem_alloc_top100.txt', all_b = 1, pert_b = 1)
oC0All_0pert_glo_gpu_top100 = git.readResultsCcode(root_pert + r'\top100\Pert_All_Ivalues_order_0_incl_abs_global_gpu_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
#compare:
git.compareCresultsAll(oC0All_glo_top100_db, oC0All_0pert_glo_gpu_top100, plot_b =1, title = 'global-alloc vs global-pert gpu, 0-pert, order 0, top100', labelDTS = 'max meanDiffToValue (Sgl)', labelNTS = 'max StdDiffToValue (Sgl)')
git.compareCresultsAll(oC0All_glo_top100_db, oC0All_0pert_glo_gpu_top100_db, plot_b =1, title = 'global-alloc vs global-pert gpu, 0-pert, order 0, top100', labelDTS = 'max meanDiffToValue (Dbl)', labelNTS = 'max StdDiffToValue (Dbl)')


#3. local_gpu and global_gpu versions give similar results on unperturbed and 5 random perts:
#unperturbed
#local non-gpu vs gpu:
#load results. "All" type (compare results over the whole simplex):
order = 0
oCAll_unpert_loc_SML = git.readResultsCcode(root_pert_db + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_local_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
oCAll_unpert_loc_gpu_SML = git.readResultsCcode(root_pert_db + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_local_gpu_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
#run comparison
git.compareCresultsAll(oCAll_unpert_loc_SML, oCAll_unpert_loc_gpu_SML)
#Results: diffs bounded by 0.077/0.28 (1ads)

#global non-gpu vs local and vs gpu vs gpu2:
order = 0
oCAll_unpert_glo_SML = git.readResultsCcode(root_pert_db + '\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
oCAll_unpert_glo_gpu_SML = git.readResultsCcode(root_pert_db + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_gpu_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
oCAll_unpert_glo_gpu2_SML = git.readResultsCcode(root_pert_db + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_gpu2_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
#run comparison
#local vs global (order 0 suffices):
git.compareCresultsAll(oCAll_unpert_loc_SML,oCAll_unpert_glo_SML)
git.compareCresultsAll(oCAll_unpert_loc_gpu_SML,oCAll_unpert_glo_gpu_SML)
#Results: all 0 (6 dec places)
#global vs global
git.compareCresultsAll(oCAll_unpert_glo_SML,oCAll_unpert_glo_gpu_SML)
git.compareCresultsAll(oCAll_unpert_glo_gpu_SML,oCAll_unpert_glo_gpu2_SML)
#Results: glo non-gpu vs glo gpu id to ditto for loc. 
#glo-gpu vs glo-gpu2 differ much less than that (bounded by 0.00077/0.022 pct, 1ads) 

#5 perts:
#local non-gpu vs gpu:
#load results. "All" type (compare results over the whole simplex):
order = 0
oCAll_5perts_loc_SML_db = git.readResultsCcode(root_pert_db + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_local_lgth10_5_perts_rnd_SML.txt', all_b = 1, pert_b =1)
oCAll_5perts_loc_gpu_SML_db = git.readResultsCcode(root_pert_db + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_local_gpu_lgth10_5_perts_rnd_SML.txt', all_b = 1, pert_b =1)
#run comparison
git.compareCresultsAll(oCAll_5perts_loc_SML_db, oCAll_5perts_loc_gpu_SML_db)
#results: bounded by 0.07/0.25 pct (1ads) quite as single prec.

#global non-gpu vs local and vs gpu vs gpu2: 
order = 0
oCAll_5perts_glo_SML_db = git.readResultsCcode(root_pert_db + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_lgth10_5_perts_rnd_SML.txt', all_b = 1, pert_b =1)
oCAll_5perts_glo_gpu_SML_db = git.readResultsCcode(root_pert_db + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_gpu_lgth10_5_perts_rnd_SML.txt', all_b = 1, pert_b =1)
oCAll_5perts_glo_gpu2_SML_db = git.readResultsCcode(root_pert_db + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_gpu2_lgth10_5_perts_rnd_SML.txt', all_b = 1, pert_b =1)
#run comparison
#local vs global (order 0 suffices):
git.compareCresultsAll(oCAll_5perts_loc_SML_db,oCAll_5perts_glo_SML_db)
git.compareCresultsAll(oCAll_5perts_loc_gpu_SML_db,oCAll_5perts_glo_gpu_SML_db)
#global vs global
git.compareCresultsAll(oCAll_5perts_glo_SML_db,oCAll_5perts_glo_gpu_SML_db)
git.compareCresultsAll(oCAll_5perts_glo_gpu_SML_db,oCAll_5perts_glo_gpu2_SML_db)
#resutls: in the gpu to gpu2 comp the diffs are bounded by approx 0.01 pct/0.05 pct now att at 1edm. The 
#1ads though gives almost the same: 0.006/0.07pct.

#4. Compare results at the "final corner": length 0-perturbations at random positions, top100 set (final values 
#and w at all vertices); we compare also the diff's at order 3 to the diff's seen in order 0:
#local, w-terms (order 0)
oC0All_0pert_loc_top100_db = git.readResultsCcode(root_pert_db + r'\top100\Pert_All_Ivalues_order_0_incl_abs_local_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_loc_gpu_top100_db = git.readResultsCcode(root_pert_db + r'\top100\Pert_All_Ivalues_order_0_incl_abs_local_gpu_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
#final values
#single prec:
oC3_0pert_loc_top100 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_local_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
oC3_0pert_loc_gpu_top100 = git.readResultsCcode(root_pert + r'\for_info\Pert_Ivalues_order_3_incl_abs_local_gpu_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
#double prec:
oC3_0pert_loc_top100_db = git.readResultsCcode(root_pert_db + r'\for_info\Pert_Ivalues_order_3_incl_abs_local_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
oC3_0pert_loc_gpu_top100_db = git.readResultsCcode(root_pert_db + r'\for_info\Pert_Ivalues_order_3_incl_abs_local_gpu_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
#Compare:
#order 0 diff's: this is the same as #2 local vs local_gpu, e.g
git.compareCresultsAll(oC0All_0pert_loc_top100_db, oC0All_0pert_loc_gpu_top100_db, plot_b = 1, title = 'local-pert vs local-pert gpu, 0-pert, order 0, top100')
#all orders, final values (run both to get the comparison in one plot; use the output from the first run as input to the second so as to get the same order in both runs ):
sortedFileList = git.compareCresultsFinal(oC3_0pert_loc_top100, oC3_0pert_loc_gpu_top100, labelRelDiff = 'max rel diff (Sgl)', labelVal = 'abs val (Sgl)', use_log_b = 1)
sortedFileList = git.compareCresultsFinal(oC3_0pert_loc_top100_db, oC3_0pert_loc_gpu_top100_db, sortedFileList = sortedFileList, title = 'local-pert vs local-pert gpu, 0-pert, rel diffs in final values, order 3, top100', labelRelDiff = 'max rel diff (Dbl)', labelVal = 'abs val (Dbl)', use_log_b = 1)
#for the log10-version: set use_log_b = 1.

#global
#w-terms (order 0)
oC0All_0pert_glo_top100 = git.readResultsCcode(root_pert_db + r'\top100\Pert_All_Ivalues_order_0_incl_abs_global_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_gpu_top100 = git.readResultsCcode(root_pert_db + r'\top100\Pert_All_Ivalues_order_0_incl_abs_global_gpu_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_gpu2_top100 = git.readResultsCcode(root_pert_db + r'\top100\Pert_All_Ivalues_order_0_incl_abs_global_gpu2_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
#final values
oC3_0pert_glo_top100_db = git.readResultsCcode(root_pert_db + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu_top100_db = git.readResultsCcode(root_pert_db + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu2_top100_db = git.readResultsCcode(root_pert_db + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu2_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
#Compare:
#order 0 diff's: id to #2 on the global-gpu versions, e.g.
git.compareCresultsAll(oC0All_0pert_glo_top100, oC0All_0pert_glo_gpu_top100, plot_b = 1, title = 'global-pert vs gobal-pert gpu, 0-pert, order 0, top100')
#all orders, final values; comparison to single prec run two-and-two:
#this is actually id to the local vs local_gpu case:
sortedFileList = git.compareCresultsFinal(oC3_0pert_glo_top100, oC3_0pert_glo_gpu_top100, labelRelDiff = 'max rel diff (Sgl)', labelVal = 'abs val (Sgl)', use_log_b = 1)
sortedFileList = git.compareCresultsFinal(oC3_0pert_glo_top100_db, oC3_0pert_glo_gpu_top100_db, sortedFileList = sortedFileList, title = 'global-pert vs global-pert gpu, 0-pert, rel diffs in final values, order 3, top100', labelRelDiff = 'max rel diff (Dbl)', labelVal = 'abs val (Dbl)', use_log_b = 1)
#for the log10-version: set use_log_b = 1.

#Showing instead the differences per invariant:
#top100
means, stdDevs = git.compareCresultsFinal_2(oC3_0pert_loc_top100_db, oC3_0pert_loc_gpu_top100_db, titleDistr = 'local-gpu, non-gpu, top100', title = 'local-gpu, non-gpu vs gpu, top100',titleRatios = 'local-gpu, non-gpu vs gpu, top100') 
medians, means, stdDevs, lP, uP = git.compareCresultsFinal_3(oC3_0pert_loc_top100_db, oC3_0pert_loc_gpu_top100_db, titleDistr = 'local-gpu, non-gpu vs gpu, top100', title = 'Statistics of relative diffs, local-gpu, non-gpu vs gpu, top100')
#top8000; better since more representative:
#final values
oC3_0pert_loc_top8000_db = git.readResultsCcode(root_pert_db + r'\for_info\Pert_Ivalues_order_3_incl_abs_local_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_loc_gpu_top8000_db = git.readResultsCcode(root_pert_db + r'\for_info\Pert_Ivalues_order_3_incl_abs_local_gpu_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
#comparisons:
means, stdDevs = git.compareCresultsFinal_2(oC3_0pert_loc_top8000_db, oC3_0pert_loc_gpu_top8000_db, titleDistr = 'local-gpu, non-gpu, top8000',title = 'local-gpu, non-gpu vs gpu, top8000',titleRatios = 'local-gpu, non-gpu vs gpu, top8000') 
medians, means, stdDevs, lP, uP = git.compareCresultsFinal_3(oC3_0pert_loc_top8000_db, oC3_0pert_loc_gpu_top8000_db, titleDistr = 'local-gpu, non-gpu vs gpu, top8000',title = 'Stats of rel diffs, local-gpu, non-gpu vs gpu, top8000 (Dbl)', figSize = (10,8))
#global version:
oC3_0pert_glo_top8000_db = git.readResultsCcode(root_pert_db + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu_top8000_db = git.readResultsCcode(root_pert_db + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu2_top8000_db = git.readResultsCcode(root_pert_db + r'\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu2_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
#comp's, non-gpu vs gpu2 (should suffice):
means, stdDevs = git.compareCresultsFinal_2(oC3_0pert_glo_top8000_db, oC3_0pert_glo_gpu2_top8000_db, titleDistr = 'global-gpu, non-gpu, top8000', title = 'global-gpu, non-gpu vs gpu2, top8000',titleRatios = 'global-gpu, non-gpu vs gpu2, top8000') 
medians, means, stdDevs, lP, uP = git.compareCresultsFinal_3(oC3_0pert_glo_top8000_db, oC3_0pert_glo_gpu2_top8000_db, titleDistr = 'global-gpu, non-gpu vs gpu2, top8000', title = 'Stats of rel diffs, global-gpu, non-gpu vs gpu2, top8000 (Dbl)', figSize = (10,8))

#plots of normalized invariants (acc to Røgen 2005) and diff's:
#local:
git.compareNormalizedInv(oC3_0pert_loc_top8000_db , oC3_0pert_loc_gpu_top8000_db, figSize = (14,10), titleDiffs = 'Differences vs values (Dbl)')
#global:
git.compareNormalizedInv(oC3_0pert_glo_top8000_db , oC3_0pert_glo_gpu_top8000_db, figSize = (14,10), titleDiffs = 'Differences vs values (Dbl)')
git.compareNormalizedInv(oC3_0pert_glo_top8000_db , oC3_0pert_glo_gpu2_top8000_db, figSize = (14,10), titleDiffs = 'Differences vs values (Dbl)') 



#run using 4th version:
#and w at all vertices); we compare also the diff's at order 3 to the diff's seen in order 0:
#local, w-terms (order 0)
oC0All_0pert_loc_top10  0_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\top100\Pert_All_Ivalues_order_0_incl_abs_local_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_loc_gpu_top100_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\top100\Pert_All_Ivalues_order_0_incl_abs_local_gpu_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
#final values
#single prec:
oC3_0pert_loc_top100 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_local_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
oC3_0pert_loc_gpu_top100 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\for_info\Pert_Ivalues_order_3_incl_abs_local_gpu_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
#double prec:
oC3_0pert_loc_top100_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\for_info\Pert_Ivalues_order_3_incl_abs_local_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
oC3_0pert_loc_gpu_top100_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\for_info\Pert_Ivalues_order_3_incl_abs_local_gpu_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
#Compare:
#order 0 diff's: this is the same as #2 local vs local_gpu, e.g
git.compareCresultsAll(oC0All_0pert_loc_top100_db, oC0All_0pert_loc_gpu_top100_db, plot_b = 1, title = 'local-pert vs local-pert gpu, 0-pert, order 0, top100')
#all orders, final values (run both to get the comparison in one plot; use the output from the first run as input to the second so as to get the same order in both runs ):
sortedFileList = git.compareCresultsFinal(oC3_0pert_loc_top100, oC3_0pert_loc_gpu_top100, labelRelDiff = 'max rel diff (Sgl)', labelVal = 'abs val (Sgl)', use_log_b = 1)
sortedFileList = git.compareCresultsFinal(oC3_0pert_loc_top100_db, oC3_0pert_loc_gpu_top100_db, sortedFileList = sortedFileList, title = 'local-pert vs local-pert gpu, 0-pert, rel diffs in final values, order 3, top100', labelRelDiff = 'max rel diff (Dbl)', labelVal = 'abs val (Dbl)', use_log_b = 1)
#for the log10-version: set use_log_b = 1.

#global
#w-terms (order 0)
oC0All_0pert_glo_top100_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\top100\Pert_All_Ivalues_order_0_incl_abs_global_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_gpu_top100_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\top100\Pert_All_Ivalues_order_0_incl_abs_global_gpu_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_gpu2_top100_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\top100\Pert_All_Ivalues_order_0_incl_abs_global_gpu2_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
#final values
oC3_0pert_glo_top100_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu_top100_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu2_top100_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu2_lgth10_0pert_rnd_top100.txt', all_b = 0, pert_b =1)
#Compare:
#order 0 diff's: id to #2 on the global-gpu versions, e.g.
git.compareCresultsAll(oC0All_0pert_glo_top100_db, oC0All_0pert_glo_gpu_top100_db, plot_b = 1, title = 'global-pert vs gobal-pert gpu, 0-pert, order 0, top100')
#all orders, final values; comparison to single prec run two-and-two:
#this is actually id to the local vs local_gpu case:
sortedFileList = git.compareCresultsFinal(oC3_0pert_glo_top100, oC3_0pert_glo_gpu_top100, labelRelDiff = 'max rel diff (Sgl)', labelVal = 'abs val (Sgl)', use_log_b = 1)
sortedFileList = git.compareCresultsFinal(oC3_0pert_glo_top100_db, oC3_0pert_glo_gpu_top100_db, sortedFileList = sortedFileList, title = 'global-pert vs global-pert gpu, 0-pert, rel diffs in final values, order 3, top100', labelRelDiff = 'max rel diff (Dbl)', labelVal = 'abs val (Dbl)', use_log_b = 1)
#for the log10-version: set use_log_b = 1.

#.. while this is new:
sortedFileList = git.compareCresultsFinal(oC3_0pert_glo_top100, oC3_0pert_glo_gpu2_top100, labelRelDiff = 'max rel diff (Sgl)', labelVal = 'abs val (Sgl)', use_log_b = 0)
sortedFileList = git.compareCresultsFinal(oC3_0pert_glo_top100_db, oC3_0pert_glo_gpu2_top100_db, sortedFileList = sortedFileList, title = 'global-pert vs global-pert gpu2, 0-pert, rel diffs in final values, order 3, top100', labelRelDiff = 'max rel diff (Dbl)', labelVal = 'abs val (Dbl)', use_log_b = 0)


#Showing instead the differences per invariant:
#top100
means, stdDevs = git.compareCresultsFinal_2(oC3_0pert_loc_top100_db, oC3_0pert_loc_gpu_top100_db, titleDistr = 'local-gpu, non-gpu, top100', title = 'local-gpu, non-gpu vs gpu, top100',titleRatios = 'local-gpu, non-gpu vs gpu, top100') 
medians, means, stdDevs, lP, uP = git.compareCresultsFinal_3(oC3_0pert_loc_top100_db, oC3_0pert_loc_gpu_top100_db, titleDistr = 'local-gpu, non-gpu vs gpu, top100', title = 'Statistics of relative diffs, local-gpu, non-gpu vs gpu, top100')
#top8000; better since more representative:
#final values
oC3_0pert_loc_top8000_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\for_info\Pert_Ivalues_order_3_incl_abs_local_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_loc_gpu_top8000_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\for_info\Pert_Ivalues_order_3_incl_abs_local_gpu_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
#comparisons:
means, stdDevs = git.compareCresultsFinal_2(oC3_0pert_loc_top8000_db, oC3_0pert_loc_gpu_top8000_db, titleDistr = 'local-gpu, non-gpu, top8000',title = 'local-gpu, non-gpu vs gpu, top8000',titleRatios = 'local-gpu, non-gpu vs gpu, top8000') 
medians, means, stdDevs, lP, uP = git.compareCresultsFinal_3(oC3_0pert_loc_top8000_db, oC3_0pert_loc_gpu_top8000_db, titleDistr = 'local-gpu, non-gpu vs gpu, top8000',title = 'Stats of rel diffs, local-gpu, non-gpu vs gpu, top8000 (Dbl)')
#global version:
oC3_0pert_glo_top8000_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu_top8000_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
oC3_0pert_glo_gpu2_top8000_db = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\double_precision\for_info\Pert_Ivalues_order_3_incl_abs_global_gpu2_lgth10_0pert_rnd_top8000.txt', all_b = 0, pert_b =1)
#comp's, non-gpu vs gpu2 (should suffice):
means, stdDevs = git.compareCresultsFinal_2(oC3_0pert_glo_top8000_db, oC3_0pert_glo_gpu2_top8000_db, titleDistr = 'global-gpu, non-gpu, top8000', title = 'global-gpu, non-gpu vs gpu2, top8000',titleRatios = 'global-gpu, non-gpu vs gpu2, top8000') 
medians, means, stdDevs, lP, uP = git.compareCresultsFinal_3(oC3_0pert_glo_top8000_db, oC3_0pert_glo_gpu2_top8000_db, titleDistr = 'global-gpu, non-gpu vs gpu2, top8000', title = 'Stats of rel diffs, global-gpu, non-gpu vs gpu2, top8000 (Dbl)', figSize = (10,8))

#plots of normalized invariants (acc to Røgen 2005) and diff's:
#local:
git.compareNormalizedInv(oC3_0pert_loc_top8000_db , oC3_0pert_loc_gpu_top8000_db, figSize = (14,10), titleDistr = 'Distribution of normalized invariants (Dbl)', titleDiffs = 'Differences vs values (Dbl)')
#global:
git.compareNormalizedInv(oC3_0pert_glo_top8000_db , oC3_0pert_glo_gpu_top8000_db, figSize = (14,10), titleDistr = 'Distribution of normalized invariants (Dbl)', titleDiffs = 'Differences vs values (Dbl)')
git.compareNormalizedInv(oC3_0pert_glo_top8000_db , oC3_0pert_glo_gpu2_top8000_db, figSize = (14,10), titleDistr = 'Distribution of normalized invariants (Dbl)', titleDiffs = 'Differences vs values (Dbl)') 



#5. By the reasoning above we only make 5a: compare for local-gpu and global-gpu in non-gpu mode, the 0-pert and the unpert:
order = 0
oCAll_unpert_loc_gpu_db = git.readResultsCcode(root_pert_db + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_local_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
oCAll_unpert_glo_gpu_db = git.readResultsCcode(root_pert_db + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
oCAll_0pert_loc_gpu_db = git.readResultsCcode(root_pert_db + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_local_lgth10_0pert_rnd_SML.txt', all_b = 1, pert_b =1)
oCAll_0pert_glo_gpu_db = git.readResultsCcode(root_pert_db + r'\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_lgth10_0pert_rnd_SML.txt', all_b = 1, pert_b =1)
#Compare
git.compareCresultsAll(oCAll_0pert_glo_gpu_db, oCAll_unpert_glo_gpu_db)
git.compareCresultsAll(oCAll_0pert_loc_gpu_db, oCAll_unpert_loc_gpu_db)
#Results:
all diffs are zero at first 6 dec places


#Aside (run in 3rd version and no diffs were seen; in single prec also run in 5th with same conclusion): compare the 
#plain C version with the cuda-based (in non-gpu mode of course), #1,2,3.
#(Comparison of base code in section on that right above): 
#1
order = 0
oCAll_unpert_plainC = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_lgth0_unpert_SML_plainC.txt', all_b = 1, pert_b =1)
oCAll_unpert_glo_gpu = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_lgth0_unpert_SML.txt', all_b = 1, pert_b =1)
#Compare
git.compareCresultsAll(oCAll_unpert_plainC, oCAll_unpert_glo_gpu)
#Results: all zero at first 6 dec places

#2
oC0All_0pert_top100_plainC = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\top100\Pert_All_Ivalues_order_0_incl_abs_global_lgth10_0pert_rnd_top100_plainC.txt', all_b = 1, pert_b =1)
oC0All_0pert_glo_top100 = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\top100\Pert_All_Ivalues_order_0_incl_abs_global_lgth10_0pert_rnd_top100.txt', all_b = 1, pert_b =1)
#compare
git.compareCresultsAll(oC0All_0pert_top100_plainC, oC0All_0pert_glo_top100, plot_b =1, title = 'plainC vs global-gpu non-gpu, 0-pert, order 0, top100')
#Results: all diffs zero


#3
#load results. "All" type (compare results over the whole simplex):
order = 3
oCAll_5perts_SML_plainC = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_lgth10_5_perts_rnd_SML_plainC.txt', all_b = 1, pert_b =1)
oCAll_5perts_global_SML_cuda = git.readResultsCcode(r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\Perturbed\single_precision\SML\Pert_All_Ivalues_order_' + str(order) + '_incl_abs_global_lgth10_5_perts_rnd_SML.txt', all_b = 1, pert_b =1)
#run comparison
git.compareCresultsAll(oCAll_5perts_SML_plainC, oCAll_5perts_global_SML_cuda)
#Result: all zero to six dec places


###################################
## time consumption, perturbations
###################################

#change root if necessary:
prec = '\single_precision'
root = r'E:\Thesis\c_code\results_C_3rd' + prec
root_pert = r'E:\Thesis\c_code\results_C_3rd\Perturbed' + prec

#can also run with
root = r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C_5th' + prec
root_pert = r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C_5th\Perturbed' + prec
#but deprecated due to loading times (near full disc)


#Plots showing time consumption by varying pert lengths:
#we need the PDBlist, so run e.g
PDBlist, oP0 = git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//projects//knots//code//top100H//', fileName = r"C:\Users\Christian\Bioinformatics\projects\knots\code\top100H_source\top100Names.txt", degree = 0, all_b =0)
#then 
path = root_pert + r'\time'
fileId1 = r'_incl_abs_global_'
lengthList = ['lgth5','lgth10', 'lgth20'] #'1': base case, i.e. global_mem_alloc or local_mem_alloc
fileId2 =  r'_100_perts_rnd_top100.txt'
order = 2
#Now run the plot fct: 
cdict = git.plotCTimeByPertLength(pathCresultFiles = path, lengthList = lengthList, PDBlist = PDBlist, root = root, order = order, fileId1 =  fileId1, fileId2 =  fileId2, timeType = 'cpuTime', title = "Time consumption (C) by varying pert lengths, order 2")

#Plots comparing the different methods:

#1. non-perturbed, global vs local mem allocation:
order = 3
set = 'top100'
method2 = 'global_mem_alloc'
path2 = root
fileId2 = r'\Ivalues_order_' + str(order) + '_excl_abs_' + method2 + '_' + set + '.txt'
fileName2 = path2 + fileId2 


method3 = 'local_mem_alloc'
path3 = root
fileId3 = r'\Ivalues_order_' + str(order) + '_excl_abs_' + method3 + '_' + set + '.txt'
fileName3 = path3 + fileId3 

pathList = [fileName3, fileName2]
methodList = [method3, method2]

#run plot fct (this gives left sub plot of two; for the aggrTime use fig = fig and subPlotNr = 122):
fig = git.plotCTimeByMethod(pathCresultFilesList = pathList, fileId1 = fileId2, fileId2 =  fileId3,  methodList = methodList, PDBlist = PDBlist, order = 2, nrOfPerts =100, timeType = 'aggrTime', title = "Comp. time (C) by method, order 2", fig = fig, subPlotNr = 122)

#for doing the plots incl "I administration", do the same, but on files from later run:
order = 3
set = 'top100'
method2 = 'global_mem_alloc'
path2 = r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\single_precision'
fileId2 = r'\Ivalues_order_' + str(order) + '_incl_abs_' + method2 + '_' + set + '_incl_collecting_I.txt'
fileName2 = path2 + fileId2 

method3 = 'local_mem_alloc'
path3 = r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C\single_precision'
fileId3 = r'\Ivalues_order_' + str(order) + '_incl_abs_' + method3 + '_' + set + '_incl_freeing_and_collecting_I.txt'
fileName3 = path3 + fileId3 

pathList = [fileName3, fileName2]
methodList = [method3, method2]

#run plot fct (this gives left sub plot of two; for the aggrTime use fig = fig and subPlotNr = 122):
fig = git.plotCTimeByMethod(pathCresultFilesList = pathList, fileId1 = fileId2, fileId2 =  fileId3,  methodList = methodList, PDBlist = PDBlist, order = 2, nrOfPerts =100, timeType = 'cpuTime', title = "Total time (C) by method, incl. I-admin, order 2", fig = '', subPlotNr = 111)



#2. non-perturbed vs perturbed (global mem alloc vs gpu2 or non-gpu), the latter for varying numbers of perturbations:

root = r'E:\Thesis\c_code\results_C_3rd'
root_pert = r'E:\Thesis\c_code\results_C_4th\Perturbed'

#can also run with
root = r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C_5th'
root_pert = r'C:\Users\Christian\Bioinformatics\Thesis\c_code\results_C_5th\Perturbed'
#but deprecated due to loading times (near full disc)

prec = '\single_precision'
root = root + prec
root_pert = root_pert + prec

#we need the PDBlist, so run e.g
PDBlist, oP0 = git.getResultsPcode(path = 'C://Users//Christian//BioInformatics//projects//knots//code//top100H//', fileName = r"C:\Users\Christian\Bioinformatics\projects\knots\code\top100H_source\top100Names.txt", degree = 0, all_b =0)

order = 2
set = 'top100'
method1 = 'global_mem_alloc'
label1 = 'main_global_alloc'
path1 = root
fileId1 = r'\Ivalues_order_' + str(order) + '_incl_abs_' + method1 + '_' + set + '.txt'
fileName1 = path1 + fileId1 

nrPerts3 = 1
method3 = 'global'
label3 = method3 + '-gpu, non-gpu ' + str(nrPerts3)
path3 = root_pert +  r'\time'
fileId3 = r'\PertTime_' + str(order) + '_incl_abs_' + method3 + '_lgth10_' + str(nrPerts3) + '_perts_rnd_' + set + '.txt'
fileName3 = path3 + fileId3  

nrPerts4 = 10
method4 = 'global'
label4 = method4 + '-gpu, non-gpu ' + str(nrPerts4)
path4 = root_pert +  r'\time'
fileId4 = r'\PertTime_' + str(order) + '_incl_abs_' + method4 + '_lgth10_' + str(nrPerts4) + '_perts_rnd_' + set + '.txt'
fileName4 = path4 + fileId4  

nrPerts5 = 25
method5 = 'global'
label5 = method5 + '-gpu, non-gpu ' + str(nrPerts5)
path5 = root_pert + r'\time'
fileId5 = r'\PertTime_' + str(order) + '_incl_abs_' + method5 + '_lgth10_' + str(nrPerts5) + '_perts_rnd_' + set + '.txt'
fileName5 = path5 + fileId5  

nrPerts6 = 100
method6 = 'global'
label6 = method6 + '-gpu, non-gpu ' + str(nrPerts6)
path6 = root_pert + r'\time'
fileId6 = r'\PertTime_' + str(order) + '_incl_abs_' + method6 + '_lgth10_' + str(nrPerts6) + '_perts_rnd_' + set + '.txt'
fileName6 = path6 + fileId6  

pathList = [fileName1, fileName3, fileName4, fileName5, fileName6]
methodList = [method1, method3, method4, method5, method6]
labelList = [label1, label3, label4, label5, label6]

#avoiding the 100 perts case
pathList = [fileName1, fileName3, fileName4, fileName5]
methodList = [method1, method3, method4, method5]
labelList = [label1, label3, label4, label5]

#run plot fct:
fig = git.plotCTimeByMethod_2(pathCresultFilesList = pathList,  PDBlist = PDBlist, 
methodList = methodList, labelList = labelList, timeType = 'cpuTime', 
title = "Total time (C) by method, 1/10/25/100 perts, order 2, top100", titleRatio = "Time ratios (C, Dbl), perturbed (global-gpu) over base case (global), order 2", fig = '', perPertPlot_b =1)

#To run the same in double precision just set prec = 'double precision' and go.


#3. perturbed vs perturbed (non-gpu, gpu2):

root = r'E:\Thesis\c_code\results_C_4th\single_precision'
root_pert = r'E:\Thesis\c_code\results_C_4th\Perturbed\single_precision'

#non-gpu vs gpu2 for perturbations:
order = 2
set = 'top100'
method1 = 'global'
nrPerts = 10
path1 = root_pert + r'\time'
fileId1 = r'\PertTime_' + str(order) + '_incl_abs_' + method1 + '_lgth10_' + str(nrPerts) + '_perts_rnd_' + set + '.txt'
fileName1 = path1 + fileId1 

method4 = 'global_gpu2'
path4 = root_pert + r'\time'
fileId4 = r'\PertTime_' + str(order) + '_incl_abs_' + method4 + '_lgth10_' + str(nrPerts) + '_perts_rnd_' + set + '.txt'
fileName4 = path4 + fileId4  

pathList = [fileName1, fileName4]
methodList = [method1, method4]

#run plot fct:
fig = git.plotCTimeByMethod(pathCresultFilesList = pathList, fileId1 = fileId1, fileId2 =  fileId2,  methodList = methodList, PDBlist = PDBlist, order = 2, nrOfPerts =100, timeType = 'cpuTime', title = "Total time (C) by method, order 2, " + str(nrPerts) + " perts, top100", fig = '')


######################################################################################################
#END OF: "FROM HERE TO LINK-SEARCHING NOT DONE"
######################################################################################################



###################################
## Search for links and "pokes"
###################################

root = r'C:\Users\Christian\Bioinformatics\papers\geometric_anomalies_in_folds\c_code\results_C_1st\double_precision'


#read in results from a run on the top100 set:
fileName = r'\ClosedLoopChars_global_mem_alloc_top100.txt'
fileName = r'\ClosedLoopChars_computeGI_top100_pokeLength10.txt'
dict_top100, links_top100, pokes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + fileName, normalize_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
pokeLength = 10
titlePokesI = 'Distribution of writhe values for potential pokes of length ' + str(pokeLength) 
titlePokesII = ' writhe values for potential pokes of length ' + str(pokeLength) 
git.histoClosedLoops(links_top100, pokes_top100, binsLinks=100, binsPokes =100, color = 'blue', linksPokesDistrInOnePlot_b =1, linksDistrPlotsInSubPlots_b =1, pokesDistrPlotsInSubPlots_b =1,cutOffLinks = 11, cutOffPokes = 6, linksShowTopNr = 100, pokesShowTopNr = 100, plotLinkExs_b = 1, plotPokeExs_b = 0, plotPokesTop10_b = 1,chainsFile = chainsFile, intervalLinks = [3,4], intervalPokes = [5,5.02], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', titlePokesI = titlePokesI, titlePokesII = titlePokesII)

#or:
fig = git.histoClosedLoops(links_top100, pokes_top100, binsLinks=100, binsPokes =100, linksDistrPlotsInSubPlots_b =2, pokesDistrPlotsInSubPlots_b =1, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 12.4, cutOffPokes = 21.0, plotPokeExs_b = 0, plotLinkExs_b = 0, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-12.9,-12.7], intervalPokes = [-5,-4], plotExsInInt_b =0, titlePokesI = titlePokesI, titlePokesII = titlePokesII)

#results of counting "extreme" link values:
There were 1 links having a writhe value belo w the cutoff: -11.000000
There were 2 links having a writhe value above the cutoff: 11.000000
There were 1 links having a writhe value in the interval: [-4, -3]
There were 85 pokes having a writhe value above the cutoff: 13.000000
There were 3 pokes having a writhe value in the interval: [5, 5.02]
#This will generate 3d-plots (at 3 angles) for each of the "extreme" link and poke examples. 
#If you want to look at an example separately (here the top scoring potential link) 
#first load the chains:
chainsDict = git.readResultsCcode_chains(chainsFile)
#Then look up an example and fetch the chain; e.g.  
links_top100.sort()
val, pdb_filename, chainId, seg1, seg2 = links_top100[0] 
polyCaChain = chainsDict[pdb_filename]['>'][1]
#and call the plot; we first generate the title of the plot:
import re
pattern = re.compile('([\S]+)([^\w]+)([\w]+)')
m = pattern.match(pdb_filename)
title = m.group(3)
git.plot_3D_2segmentsHighLight(polyCaChain, segment1 = seg1, segment2 = seg2, title = title)
#for particular poke examples: use the poke_top100 list instead.

git.plot_3D_Quiver_2segmentsHighLight(polyCaChain, segment1 = seg1, segment2 = seg2, title = title)


#top8000 set:
fileName = r'\ClosedLoopChars_computeGI_top8000_pokeLength10.txt'
dict_top8000_cl, links_top8000_cl, pokes_top8000_cl = git.readResultsCcode_closedLoops(inputFileName = root + fileName)
pokeLength = 10
chainsFile = root + r'\chains_closed_loops_top8000.txt'
titlePokesI = 'Distribution of writhe values for potential pokes of length ' + str(pokeLength) 
titlePokesII = ' writhe values for potential pokes of length ' + str(pokeLength) 
fig = git.histoClosedLoops(links_top8000_cl, pokes_top8000_cl, binsLinks=100, binsPokes =100, linksDistrPlotsInSubPlots_b =2, pokesDistrPlotsInSubPlots_b =1, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 12.4, cutOffPokes = 21.0, plotPokeExs_b = 0, plotLinkExs_b = 0, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-12.9,-12.7], intervalPokes = [-5,-4], plotExsInInt_b =0, titlePokesI = titlePokesI, titlePokesII = titlePokesII)
#results of counting "extreme" link values, poke lgth 7:
There were 4 links having a writhe value below the cutoff: -12.000000
There were 7 links having a writhe value above the cutoff: 12.000000
There were 14 links having a writhe value in the interval: [-1e-06, 1e-06]
There were 12 pokes having a writhe value above the cutoff: 20.000000
There were 21 pokes having a writhe value in the interval: [5, 5.001]

#With default cut-offs instead:
There were 31 links having a writhe value below the cutoff: -11.000000
There were 22 links having a writhe value above the cutoff: 11.000000
There were 80 links having a writhe value in the interval: [-1e-06, 1e-06]
There were 461 pokes having a writhe value above the cutoff: 13.000000
There were 26 pokes having a writhe value in the interval: [5, 5.001]


#distr for links in top100 and top8000 in one plot:
fig = git.histoClosedLoops(links_top100, pokes_top100, binsLinks=100, binsPokes =100, titleDistr = 'top100', linksDistrPlotsInSubPlots_b =2, pokesDistrPlotsInSubPlots_b =1, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 12.4, cutOffPokes = 21.0, plotPokeExs_b = 0, plotLinkExs_b = 0, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-12.9,-12.7], intervalPokes = [-5,-4], plotExsInInt_b =0, titlePokesI = titlePokesI, titlePokesII = titlePokesII)
fig = git.histoClosedLoops(links_top8000_cl, pokes_top8000_cl, binsLinks=100, binsPokes =100, inFig = fig, titleDistr = 'top8000', linksDistrPlotsInSubPlots_b =2, pokesDistrPlotsInSubPlots_b =1, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 12.4, cutOffPokes = 21.0, plotPokeExs_b = 0, plotLinkExs_b = 0, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-12.9,-12.7], intervalPokes = [-5,-4], plotExsInInt_b =0, titlePokesI = titlePokesI, titlePokesII = titlePokesII)



##
## More general search: consider writhe results from all pairs of sub-chains
## of a set length (e.g 30). Philosophy: those pairs having a high (abs) 
## mutual writhe have some interesting geometry (e.g. could be links)
##

#read in results from a run on the top100 set:
dict_top100, extremeWrithes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth30.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
git.histoSubChainPairs(extremeWrithes_top100, bins=100, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth = 30)

dict_top8000, extremeWrithes_top8000 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth30.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
git.histoSubChainPairs(extremeWrithes_top8000, inFig = fig, bins=100, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 0, plotExs_b = 0, chainsFile = chainsFile, interval = [-12.9,-12.7], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+).pdb', subChainLgth = 30)
 
#Similarly for the case of sub-chains of length 15:

dict_top100, extremeWrithes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth15.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
git.histoSubChainPairs(extremeWrithes_top100, bins=100, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth = 30)
 

dict_top8000, extemeWrithes_top8000 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth15.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
git.histoSubChainPairs(extemeWrithes_top8000, bins=100, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 0, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+).pdb', subChainLgth = 15)
 

 
#read in results from a run on the top100 set:
dict_top100, extemeWrithes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_incl_abs_global_mem_alloc_top100_I1234_full.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
git.histoSubChainPairs(extemeWrithes_top100, bins=100, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth = 15)
 

##################################
# Generating plots for links-paper
##################################

import gitP as git

root = r'C:\Users\Christian\Bioinformatics\papers\geometric_anomalies_in_folds\c_code\results_C_3rd\double_precision'
root = r'C:\Users\Christian\Bioinformatics\papers\geometric_anomalies_in_folds\c_code\results_C_4th\double_precision'
root = r'D:\Bioinformatics\papers\geometric_anomalies_in_folds\c_code\results_C_5th\double_precision'



#top100 set:

#restricted search:
pokeLength = 10
loopLength = 30
#fileName = root + r'\ClosedLoopChars_global_mem_alloc_top100_pokeLgth' + str(pokeLength) + '.txt'
fileName = root + r'\ClosedLoopChars_computeGI_top100_loops' + str(loopLength) + '_pokeLength' + str(pokeLength) + '.txt'
dict_top100_cl, links_top100_cl, pokesAll_top100_cl, pokesMin_top100_cl, pokesMax_top100_cl = git.readResultsCcode_closedLoops(inputFileName = fileName, normalize_b = 1, minMaxOut_b =1)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
titlePokesI = 'Distribution of writhe values for potential pokes of length ' + str(pokeLength) 
titlePokesII = ' writhe values for potential pokes of length ' + str(pokeLength) 
pokes_top100_cl = pokesAll_top100_cl, pokesMin_top100_cl, pokesMax_top100_cl
fig = git.histoClosedLoops_minMax(links_top100_cl, pokes_top100_cl, pokesMinMax_b =1, color = 'lightgrey', binsLinks=100, binsPokes =100, linksDistrPlotsInSubPlots_b =2, pokesDistrPlotsInSubPlots_b =1, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 0.85, cutOffPokes = 0.8, plotPokeExs_b = 0, plotLinkExs_b = 1, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-1.05,-0.95], intervalPokes = [-0.5,-0.4], plotExsInInt_b =0, titleLinks = 'Distribution of writhe values', titlePokesI = titlePokesI, titlePokesII = titlePokesII, stringPattern = '([\S]+)([^\w]+)([\w]+)', lowerPercentile = 0.001, upperPercentile = 0.999)

#unrestricted search, lgth 15:
subChainLength = 15
#fileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth15.txt'
fileName= root + r'\SubChainPairChars_minmax_computeGI_top100_subChainLength' + str(subChainLength) + '.txt'
dict_top100, extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max = git.readResultsCcode_closedLoops(inputFileName = fileName, linksAndPokes_b = 0, minMaxOut_b = 1, normalize_b = 1)
###no min max:
###dict_top100, extremeWrithes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth15.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
#convert structure and then plot
extremeWrithes_top100 = extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max
git.histoSubChainPairs_minMax(extremeWrithes_top100, bins=100, color = 'lightgrey', addSubplotToFirstPlot_b = 1, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  15, lowerPercentile = 0.05, upperPercentile = 0.95)

#unrestricted search, lgth 30:
subChainLength = 30
#fileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth15.txt'
fileName= root + r'\SubChainPairChars_minmax_computeGI_top100_subChainLength' + str(subChainLength) + '.txt'
dict_top100, extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max = git.readResultsCcode_closedLoops(inputFileName = fileName, linksAndPokes_b = 0, minMaxOut_b = 1, normalize_b = 1)
###no min max:
###dict_top100, extremeWrithes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth30.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
#convert structure and then plot
extremeWrithes_top100 = extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max
git.histoSubChainPairs_minMax(extremeWrithes_top100, bins=100, color = 'lightgrey', addSubplotToFirstPlot_b = 1, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  30, lowerPercentile = 0.05, upperPercentile = 0.95)



#top8000 set:

#restricted search:
pokeLength = 10
loopLength = 30
fileName = root + r'\ClosedLoopChars_minmax_computeGI_top8000_loops' + str(loopLength) + '_pokeLength' + str(pokeLength) + '.txt'
dict_top8000_cl, links_top8000_cl, pokesAll_top8000_cl,  pokesMin_top8000_cl, pokesMax_top8000_cl = git.readResultsCcode_closedLoops(inputFileName = fileName, normalize_b = 1, minMaxOut_b = 1)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
titlePokesI = 'Distribution of writhe values for potential pokes of length ' + str(pokeLength) 
titlePokesII = ' writhe values for potential pokes of length ' + str(pokeLength) 
pokes_top8000_cl = pokesAll_top8000_cl, pokesMin_top8000_cl, pokesMax_top8000_cl
fig = git.histoClosedLoops_minMax(links_top8000_cl, pokes_top8000_cl, pokesMinMax_b =1, color = 'grey', supTitle_b = 0, binsLinks=100, binsPokes =100, linksDistrPlotsInSubPlots_b =2, pokesDistrPlotsInSubPlots_b =1, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 0.95, cutOffPokes = 0.9, plotPokeExs_b = 0, plotLinkExs_b = 0, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-1.1,-.95], intervalPokes = [-0.75,-.74], plotExsInInt_b =0, titlePokesI = titlePokesI, titlePokesII = titlePokesII, lowerPercentile = 0.00001, upperPercentile = 0.99999)

#unrestricted search, lgth 15:
subChainLength = 15
#fileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth15.txt'
fileName= root + r'\SubChainPairChars_minmax_computeGI_top8000_subChainLength' + str(subChainLength) + '.txt'
dict_top8000, extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max  = git.readResultsCcode_closedLoops(inputFileName = fileName, linksAndPokes_b = 0, minMaxOut_b = 1, normalize_b = 1)
###no min max:
###dict_top8000, extremeWrithes_top8000 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth15.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
#convert structure and then plot
extremeWrithes_top8000 = extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max
git.histoSubChainPairs_minMax(extremeWrithes_top8000, bins=100, logCount_b = 0, color = 'lightgrey', addSubplotToFirstPlot_b = 1, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  15, lowerPercentile = 0.01, upperPercentile = 0.99)

#unrestricted search, lgth 30:
subChainLength = 30
#fileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth30.txt'
fileName= root + r'\SubChainPairChars_minmax_computeGI_top8000_subChainLength' + str(subChainLength) + '.txt'
dict_top8000, extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max  = git.readResultsCcode_closedLoops(inputFileName = fileName, linksAndPokes_b = 0, minMaxOut_b = 1, normalize_b = 1)
###no min max:
###dict_top8000, extremeWrithes_top8000 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth30.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
#convert structure and then plot
extremeWrithes_top8000 = extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max
git.histoSubChainPairs_minMax(extremeWrithes_top8000, bins=100, logCount_b = 0, color = 'grey', colorMin = 'lightgrey', colorMax = 'grey', alphaMin = 1.0, alphaMax = 0.8, supTitle_b = 0, titleDistr = 'Distribution of mutual writhe values', addSubplotToFirstPlot_b = 1, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 0, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  30, lowerPercentile = 0.005, upperPercentile = 0.995)





###########
#twin plot of pymol images, 1bpi and 2cpl:
##########

root = r'C:\Users\Christian\Bioinformatics\papers\geometric_anomalies_in_folds\plots\for_main_note'

import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
#from mpl_toolkits.axes_grid1 import AxesGrid
#
#grid = AxesGrid(fig, 121,  # similar to subplot(141)
#                    nrows_ncols=(1, 2),
#                    share_all=True,
#                    axes_pad=0.05,
#                    label_mode="0",
#                    )

image_file = cbook.get_sample_data(root + r'\1bpi_unrestricted_search_pymol_resized.png')
image = plt.imread(image_file)

fig = plt.figure(figsize = (8,2))
#plt.margins(xmargin =1, ymargin = 1, tight = True)
ax = fig.add_subplot(121)

ax.imshow(image)

#grid[0].imshow(image)

plt.axis('off')  # clear x- and y-axes

image_file = cbook.get_sample_data(root + r'\2cpl_unrestricted_search_pymol_2ndview_resized.png')
image = plt.imread(image_file)

ax = fig.add_subplot(122)
ax.imshow(image)

#grid[1].imshow(image)

plt.axis('off')  # clear x- and y-axes
#plt.tight_layout
plt.subplots_adjust(top= 0.95, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.01)
plt.show()


###########
# Plots for Supplementary Data, links-paper
###########

import gitP as git

#root: as for plots for the paper
#root = r'C:\Users\Christian\Bioinformatics\papers\geometric_anomalies_in_folds\c_code\results_C_3rd\double_precision'
#first runs on: root = r'C:\Users\Christian\Bioinformatics\papers\geometric_anomalies_in_folds\c_code\results_C_1st\double_precision'
root = r'C:\Users\Christian\Bioinformatics\papers\geometric_anomalies_in_folds\c_code\results_C_4th\double_precision'

#restricted search:

#top100:

#Distr's of pokes and links for top100 set in one plot:
#read in results from a run on the top100 set:
#restricted search:
pokeLength = 10
loopLength = 30
#fileName = root + r'\ClosedLoopChars_global_mem_alloc_top100_pokeLgth' + str(pokeLength) + '.txt'
fileName = root + r'\ClosedLoopChars_computeGI_top100_loops' + str(loopLength) + '_pokeLength' + str(pokeLength) + '.txt'
dict_top100_cl, links_top100_cl, pokesAll_top100_cl, pokesMin_top100_cl, pokesMax_top100_cl = git.readResultsCcode_closedLoops(inputFileName = fileName, normalize_b = 1, minMaxOut_b =1)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
titleLinks = 'Distribution of writhe values for potential links and pokes of length ' + str(pokeLength) 
titlePokesI = 'Distribution of writhe values for potential pokes of length ' + str(pokeLength) 
titlePokesII = ' writhe values for potential pokes of length ' + str(pokeLength) 
pokes_top100_cl = pokesAll_top100_cl, pokesMin_top100_cl, pokesMax_top100_cl
fig = git.histoClosedLoops_minMax(links_top100_cl, pokes_top100_cl, pokesMinMax_b =1, colorMin = 'lightblue', color = 'blue', colorMax = 'blue', binsLinks=100, binsPokes =100, linksPokesDistrInOnePlot_b =1, linksDistrPlotsInSubPlots_b =0, pokesDistrPlotsInSubPlots_b =0, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 0.85, cutOffPokes = 0.8, plotPokeExs_b = 0, plotLinkExs_b = 1, plotLinksTop10_b = 1, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-1.05,-0.95], intervalPokes = [-0.5,-0.4], plotExsInInt_b =0, titleLinks = titleLinks, titlePokesI = titlePokesI, titlePokesII = titlePokesII, stringPattern = '([\S]+)([^\w]+)([\w]+)', lowerPercentile = 0.001, upperPercentile = 0.999)

#to look at the lower writhe value examples:
git.histoClosedLoops_minMax(links_top100_cl, pokes_top100_cl, binsLinks=100, binsPokes =100, colorMin = 'lightblue', color = 'blue', colorMax = 'blue',  elevAzimList = [[15,-45],[160,-45],[-160,-130]], linksPokesDistrInOnePlot_b =1, linksDistrPlotsInSubPlots_b =1, pokesDistrPlotsInSubPlots_b =1,cutOffLinks = 0.9, cutOffPokes = 0.85, linksShowTopNr = 100, pokesShowTopNr = 100, plotLinkExs_b = 0, plotPokeExs_b = 0, plotPokesTop10_b = 0,chainsFile = chainsFile, intervalLinks = [0,0.00001], intervalPokes = [0.25,0.251], plotExsInInt_b =1, stringPattern = '([\S]+)([^\w]+)([\w]+)', titlePokesI = titlePokesI, titlePokesII = titlePokesII)

#med pos links:
git.histoClosedLoops_minMax(links_top100_cl, pokes_top100_cl, binsLinks=100, binsPokes =100, colorMin = 'lightblue', color = 'blue', colorMax = 'blue',  elevAzimList = [[15,-45],[160,-45],[-160,-130]], linksPokesDistrInOnePlot_b =1, linksDistrPlotsInSubPlots_b =1, pokesDistrPlotsInSubPlots_b =1,cutOffLinks = 0.9, cutOffPokes = 0.85, linksShowTopNr = 100, pokesShowTopNr = 100, plotLinkExs_b = 0, plotPokeExs_b = 0, plotPokesTop10_b = 0,chainsFile = chainsFile, intervalLinks = [0.5,0.55], intervalPokes = [0.25,0.251], plotExsInInt_b =1, stringPattern = '([\S]+)([^\w]+)([\w]+)', titlePokesI = titlePokesI, titlePokesII = titlePokesII)
#med neg links:
git.histoClosedLoops_minMax(links_top100_cl, pokes_top100_cl, binsLinks=100, binsPokes =100, colorMin = 'lightblue', color = 'blue', colorMax = 'blue',  elevAzimList = [[15,-45],[160,-45],[-160,-130]], linksPokesDistrInOnePlot_b =1, linksDistrPlotsInSubPlots_b =1, pokesDistrPlotsInSubPlots_b =1,cutOffLinks = 0.9, cutOffPokes = 0.85, linksShowTopNr = 100, pokesShowTopNr = 100, plotLinkExs_b = 0, plotPokeExs_b = 0, plotPokesTop10_b = 0,chainsFile = chainsFile, intervalLinks = [-0.25, -.24], intervalPokes = [0.25,0.251], plotExsInInt_b =1, stringPattern = '([\S]+)([^\w]+)([\w]+)', titlePokesI = titlePokesI, titlePokesII = titlePokesII)

#Similarly for poke lengths 5 and 7:

pokeLength = 5
loopLength = 30
#fileName = root + r'\ClosedLoopChars_global_mem_alloc_top100_pokeLgth' + str(pokeLength) + '.txt'
fileName = root + r'\ClosedLoopChars_computeGI_top100_loops' + str(loopLength) + '_pokeLength' + str(pokeLength) + '.txt'
dict_top100_cl, links_top100_cl, pokesAll_top100_cl, pokesMin_top100_cl, pokesMax_top100_cl = git.readResultsCcode_closedLoops(inputFileName = fileName, normalize_b = 1, minMaxOut_b =1)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
titleLinks = 'Distribution of writhe values for potential links and pokes of length ' + str(pokeLength) 
titlePokesI = 'Distribution of writhe values for potential pokes of length ' + str(pokeLength) 
titlePokesII = ' writhe values for potential pokes of length ' + str(pokeLength) 
pokes_top100_cl = pokesAll_top100_cl, pokesMin_top100_cl, pokesMax_top100_cl
fig = git.histoClosedLoops_minMax(links_top100_cl, pokes_top100_cl, pokesMinMax_b =1, colorMin = 'lightblue', color = 'blue', colorMax = 'blue', binsLinks=100, binsPokes =100, linksPokesDistrInOnePlot_b =1, linksDistrPlotsInSubPlots_b =0, pokesDistrPlotsInSubPlots_b =0, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 0.85, cutOffPokes = 0.8, plotPokeExs_b = 0, plotLinkExs_b = 1, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-1.05,-0.95], intervalPokes = [-0.5,-0.4], plotExsInInt_b =0, titleLinks = titleLinks, titlePokesI = titlePokesI, titlePokesII = titlePokesII, stringPattern = '([\S]+)([^\w]+)([\w]+)', lowerPercentile = 0.001, upperPercentile = 0.999)


pokeLength = 7
loopLength = 30
#fileName = root + r'\ClosedLoopChars_global_mem_alloc_top100_pokeLgth' + str(pokeLength) + '.txt'
fileName = root + r'\ClosedLoopChars_computeGI_top100_loops' + str(loopLength) + '_pokeLength' + str(pokeLength) + '.txt'
dict_top100_cl, links_top100_cl, pokesAll_top100_cl, pokesMin_top100_cl, pokesMax_top100_cl = git.readResultsCcode_closedLoops(inputFileName = fileName, normalize_b = 1, minMaxOut_b =1)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
titleLinks = 'Distribution of writhe values for potential links and pokes of length ' + str(pokeLength) 
titlePokesI = 'Distribution of writhe values for potential pokes of length ' + str(pokeLength) 
titlePokesII = ' writhe values for potential pokes of length ' + str(pokeLength) 
pokes_top100_cl = pokesAll_top100_cl, pokesMin_top100_cl, pokesMax_top100_cl
fig = git.histoClosedLoops_minMax(links_top100_cl, pokes_top100_cl, pokesMinMax_b =1, colorMin = 'lightblue', color = 'blue', colorMax = 'blue', binsLinks=100, binsPokes =100, linksPokesDistrInOnePlot_b =1, linksDistrPlotsInSubPlots_b =0, pokesDistrPlotsInSubPlots_b =0, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 0.85, cutOffPokes = 0.8, plotPokeExs_b = 0, plotLinkExs_b = 1, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-1.05,-0.95], intervalPokes = [-0.5,-0.4], plotExsInInt_b =0, titleLinks = titleLinks, titlePokesI = titlePokesI, titlePokesII = titlePokesII, stringPattern = '([\S]+)([^\w]+)([\w]+)', lowerPercentile = 0.001, upperPercentile = 0.999)



#top8000:

#Distr's of pokes and links for top8000 set in one plot:
#read in results from a run on the top8000 set:
#restricted search:
pokeLength = 10
loopLength = 30
fileName = root + r'\ClosedLoopChars_minmax_computeGI_top8000_loops' + str(loopLength) + '_pokeLength' + str(pokeLength) + '.txt'
dict_top8000_cl, links_top8000_cl, pokesAll_top8000_cl,  pokesMin_top8000_cl, pokesMax_top8000_cl = git.readResultsCcode_closedLoops(inputFileName = fileName, normalize_b = 1, minMaxOut_b = 1)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
titleLinks = 'Distribution of writhe values for potential links and pokes of length ' + str(pokeLength) 
titlePokesI = 'Distribution of writhe values for potential pokes of length ' + str(pokeLength) 
titlePokesII = ' writhe values for potential pokes of length ' + str(pokeLength) 
pokes_top8000_cl = pokesAll_top8000_cl, pokesMin_top8000_cl, pokesMax_top8000_cl
fig = git.histoClosedLoops_minMax(links_top8000_cl, pokes_top8000_cl, pokesMinMax_b =1, pokesDistrLog_b = 0, colorMin = 'lightblue', color = 'blue', colorMax = 'blue', supTitle_b = 1, binsLinks=100, binsPokes =100, linksPokesDistrInOnePlot_b =1, linksDistrPlotsInSubPlots_b =0, pokesDistrPlotsInSubPlots_b =0, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 0.95, cutOffPokes = 0.9, plotPokeExs_b = 0, plotLinkExs_b = 0, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-1.1,-.95], intervalPokes = [-0.75,-.74], plotExsInInt_b =0, titleLinks = titleLinks, titlePokesI = titlePokesI, titlePokesII = titlePokesII, lowerPercentile = 0.00001, upperPercentile = 0.99999)

#Separate plot of links and plots:
titleLinks = 'Distribution of writhe values for potential links'
titlePokesI = 'Distribution of writhe values for potential pokes of length ' + str(pokeLength) 
titlePokesII = ' writhe values for potential pokes of length ' + str(pokeLength) 
fig = git.histoClosedLoops_minMax(links_top8000_cl, pokes_top8000_cl, pokesMinMax_b =1, pokesDistrLog_b = 1, colorMin = 'lightblue', color = 'blue', colorMax = 'blue', supTitle_b = 1, binsLinks= 50, binsPokes = 50, linksPokesDistrInOnePlot_b =0, linksDistrPlotsInSubPlots_b =1, pokesDistrPlotsInSubPlots_b =1, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 0.95, cutOffPokes = 0.9, plotPokeExs_b = 0, plotLinkExs_b = 0, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-1.1,-.95], intervalPokes = [-0.75,-.74], plotExsInInt_b =0, titleLinks = titleLinks, titlePokesI = titlePokesI, titlePokesII = titlePokesII, lowerPercentile = 0.00001, upperPercentile = 0.99999)


#Similarly for poke lengths 5 and 7:

pokeLength = 5
loopLength = 30
fileName = root + r'\ClosedLoopChars_minmax_computeGI_top8000_loops' + str(loopLength) + '_pokeLength' + str(pokeLength) + '.txt'
dict_top8000_cl, links_top8000_cl, pokesAll_top8000_cl,  pokesMin_top8000_cl, pokesMax_top8000_cl = git.readResultsCcode_closedLoops(inputFileName = fileName, normalize_b = 1, minMaxOut_b = 1)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
titleLinks = 'Distribution of writhe values for potential links and pokes of length ' + str(pokeLength) 
titlePokesI = 'Distribution of writhe values for potential pokes of length ' + str(pokeLength) 
titlePokesII = ' writhe values for potential pokes of length ' + str(pokeLength) 
pokes_top8000_cl = pokesAll_top8000_cl, pokesMin_top8000_cl, pokesMax_top8000_cl
fig = git.histoClosedLoops_minMax(links_top8000_cl, pokes_top8000_cl, pokesMinMax_b =1, pokesDistrLog_b = 1, colorMin = 'lightblue', color = 'blue', colorMax = 'blue', supTitle_b = 1, binsLinks=100, binsPokes =100, linksPokesDistrInOnePlot_b =1, linksDistrPlotsInSubPlots_b =0, pokesDistrPlotsInSubPlots_b =0, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 0.95, cutOffPokes = 0.9, plotPokeExs_b = 0, plotLinkExs_b = 0, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-1.1,-.95], intervalPokes = [-0.75,-.74], plotExsInInt_b =0, titleLinks = titleLinks, titlePokesI = titlePokesI, titlePokesII = titlePokesII, lowerPercentile = 0.00001, upperPercentile = 0.99999)

#Separate plot of links and plots:
titleLinks = 'Distribution of writhe values for potential links'
fig = git.histoClosedLoops_minMax(links_top8000_cl, pokes_top8000_cl, pokesMinMax_b =1, pokesDistrLog_b = 1, colorMin = 'lightblue', color = 'blue', colorMax = 'blue', supTitle_b = 1, binsLinks=50, binsPokes = 50, linksPokesDistrInOnePlot_b =0, linksDistrPlotsInSubPlots_b =1, pokesDistrPlotsInSubPlots_b =1, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 0.95, cutOffPokes = 0.9, plotPokeExs_b = 0, plotLinkExs_b = 0, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-1.1,-.95], intervalPokes = [-0.75,-.74], plotExsInInt_b =0, titleLinks = titleLinks, titlePokesI = titlePokesI, titlePokesII = titlePokesII, lowerPercentile = 0.00001, upperPercentile = 0.99999)


pokeLength = 7
loopLength = 30
fileName = root + r'\ClosedLoopChars_minmax_computeGI_top8000_loops' + str(loopLength) + '_pokeLength' + str(pokeLength) + '.txt'
dict_top8000_cl, links_top8000_cl, pokesAll_top8000_cl,  pokesMin_top8000_cl, pokesMax_top8000_cl = git.readResultsCcode_closedLoops(inputFileName = fileName, normalize_b = 1, minMaxOut_b = 1)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
titleLinks = 'Distribution of writhe values for potential links and pokes of length ' + str(pokeLength) 
titlePokesI = 'Distribution of writhe values for potential pokes of length ' + str(pokeLength) 
titlePokesII = ' writhe values for potential pokes of length ' + str(pokeLength) 
pokes_top8000_cl = pokesAll_top8000_cl, pokesMin_top8000_cl, pokesMax_top8000_cl
fig = git.histoClosedLoops_minMax(links_top8000_cl, pokes_top8000_cl, pokesMinMax_b =1, pokesDistrLog_b = 1, colorMin = 'lightblue', color = 'blue', colorMax = 'blue', supTitle_b = 1, binsLinks=100, binsPokes =100, linksPokesDistrInOnePlot_b =1, linksDistrPlotsInSubPlots_b =0, pokesDistrPlotsInSubPlots_b =0, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 0.95, cutOffPokes = 0.9, plotPokeExs_b = 0, plotLinkExs_b = 0, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-1.1,-.95], intervalPokes = [-0.75,-.74], plotExsInInt_b =0, titleLinks = titleLinks, titlePokesI = titlePokesI, titlePokesII = titlePokesII, lowerPercentile = 0.00001, upperPercentile = 0.99999)

#Separate plot of links and plots:
titleLinks = 'Distribution of writhe values for potential links'
fig = git.histoClosedLoops_minMax(links_top8000_cl, pokes_top8000_cl, pokesMinMax_b =1, pokesDistrLog_b = 1, colorMin = 'lightblue', color = 'blue', colorMax = 'blue', supTitle_b = 1, binsLinks=50, binsPokes = 50, linksPokesDistrInOnePlot_b =0, linksDistrPlotsInSubPlots_b =1, pokesDistrPlotsInSubPlots_b =1, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 0.95, cutOffPokes = 0.9, plotPokeExs_b = 0, plotLinkExs_b = 0, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-1.1,-.95], intervalPokes = [-0.75,-.74], plotExsInInt_b =0, titleLinks = titleLinks, titlePokesI = titlePokesI, titlePokesII = titlePokesII, lowerPercentile = 0.00001, upperPercentile = 0.99999)



#For unrestricted search: 

#top100:
#distr of writhe for lgth 15 and 30 in one plot ... plus "interesting examples":
#unrestricted search, lgth 15:
lgth = 15
#fileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_' + str(lgth) + '.txt'
fileName = root + r'\SubChainPairChars_minmax_computeGI_top100_subChainLength' + str(lgth) + '.txt'
dict_top100, extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max = git.readResultsCcode_closedLoops(inputFileName = fileName, linksAndPokes_b = 0, minMaxOut_b = 1, normalize_b = 1)
###no min max:
###dict_top100, extremeWrithes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth15.txt', linksAndPokes_b = 0, normalize_b = 1)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
#convert structure and then plot
extremeWrithes_top100 = extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max
git.histoSubChainPairs_minMax(extremeWrithes_top100, bins=100, color = 'blue', colorMin = 'lightblue', colorMax = 'blue', addSubplotToFirstPlot_b = 1, addSubPlotNr = 121, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-8.45,-8.4], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  15, lowerPercentile = 0.05, upperPercentile = 0.95)

#two projections in examples:
elevAzimList = [[15,-45],[160,-45]]
git.histoSubChainPairs_minMax(extremeWrithes_top100, bins=100, elevAzimList = elevAzimList, color = 'blue', colorMin = 'orange', colorMax = 'blue', addSubplotToFirstPlot_b = 1, addSubPlotNr = 121, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-8.45,-8.4], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  15, lowerPercentile = 0.05, upperPercentile = 0.95)


# lgth 30:
lgth = 30
#fileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_' + str(lgth) + '.txt'
fileName = root + r'\SubChainPairChars_minmax_computeGI_top100_subChainLength' + str(lgth) + '.txt'
dict_top100, extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max = git.readResultsCcode_closedLoops(inputFileName = fileName, linksAndPokes_b = 0, minMaxOut_b = 1, normalize_b = 1)
###no min max:
###dict_top100, extremeWrithes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth15.txt', linksAndPokes_b = 0, normalize_b = 1)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
#convert structure and then plot
extremeWrithes_top100 = extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max
git.histoSubChainPairs_minMax(extremeWrithes_top100, bins=100, color = 'blue', colorMin = 'lightblue', colorMax = 'blue', addSubplotToFirstPlot_b = 1, addSubPlotNr = 122, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-11,-10], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  30, lowerPercentile = 0.05, upperPercentile = 0.95)

#two projections in examples:
elevAzimList = [[15,-45],[160,-45]]
git.histoSubChainPairs_minMax(extremeWrithes_top100, bins=100, elevAzimList = elevAzimList, color = 'blue', colorMin = 'orange', colorMax = 'blue', addSubplotToFirstPlot_b = 1, addSubPlotNr = 121, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-8.45,-8.4], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  30, lowerPercentile = 0.05, upperPercentile = 0.95)


#No-min/max version (used for 1st version of plots/thesis):
#first runs on (?): 
root = r'C:\Users\Christian\Bioinformatics\papers\geometric_anomalies_in_folds\c_code\results_C_1st\double_precision'
#distr of writhe for lgth 15 and 30 in one plot ... plus "interesting examples":
#unrestricted search, lgth 15:
dict_top100, extremeWrithes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth15.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
git.histoSubChainPairs(extremeWrithes_top100, bins=100, color = 'blue', elevAzimList = [[15,-45],[160,-45],[-160,-130]], addSubplotToFirstPlot_b = 1, addSubPlotNr = 121, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  15)


#unrestricted search, lgth 30:
dict_top100, extremeWrithes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth30.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
git.histoSubChainPairs(extremeWrithes_top100, bins=100, color = 'blue',elevAzimList = [[15,-45],[160,-45],[-160,-130]], addSubplotToFirstPlot_b = 1, addSubPlotNr = 122, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth = 30)




#top8000:

#distr of writhe for lgth 15 and 30 in one plot ... plus "interesting examples":
#unrestricted search, lgth 15:
dict_top8000, extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth15.txt', linksAndPokes_b = 0, minMaxOut_b = 1, normalize_b = 1)
###no min max:
###dict_top100, extremeWrithes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth15.txt', linksAndPokes_b = 0, normalize_b = 1)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
#convert structure and then plot
extremeWrithes_top8000 = extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max
git.histoSubChainPairs_minMax(extremeWrithes_top8000, bins=100, color = 'blue', colorMin = 'orange', colorMax = 'blue', addSubplotToFirstPlot_b = 1, addSubPlotNr = 122, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-11,-10], plotExsInInt_b =0, subChainLgth =  15, lowerPercentile = 0.001, upperPercentile = 0.999)

#two projections in examples:
elevAzimList = [[15,-45],[160,-45]]
git.histoSubChainPairs_minMax(extremeWrithes_top8000, bins=100, elevAzimList = elevAzimList, color = 'blue', colorMin = 'orange', colorMax = 'blue', addSubplotToFirstPlot_b = 1, addSubPlotNr = 122, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-11,-10], plotExsInInt_b =0, subChainLgth =  15, lowerPercentile = 0.05, upperPercentile = 0.95)


# lgth 30:
dict_top8000, extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth30.txt', linksAndPokes_b = 0, minMaxOut_b = 1)
###no min max:
###dict_top100, extremeWrithes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth15.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
#convert structure and then plot
extremeWrithes_top8000 = extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max
git.histoSubChainPairs_minMax(extremeWrithes_top8000, bins=100, color = 'blue', colorMin = 'orange', colorMax = 'blue', addSubplotToFirstPlot_b = 1, addSubPlotNr = 122, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 0, plotExs_b = 0, chainsFile = chainsFile, interval = [-11,-10], plotExsInInt_b =0, subChainLgth =  30, lowerPercentile = 0.001, upperPercentile = 0.999)

#two projections in examples:
elevAzimList = [[15,-45],[160,-45]]
git.histoSubChainPairs_minMax(extremeWrithes_top8000, bins=100, elevAzimList = elevAzimList, color = 'blue', colorMin = 'orange', colorMax = 'blue', addSubplotToFirstPlot_b = 1, addSubPlotNr = 122, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-11,-10], plotExsInInt_b =0, subChainLgth =  30, lowerPercentile = 0.05, upperPercentile = 0.95)


#No-min/max version (used for 1st version of plots/thesis):
#distr of writhe for lgth 15 and 30 in one plot ... plus "interesting examples":
#unrestricted search, lgth 15:
dict_top8000, extremeWrithes_top8000 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth15.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
git.histoSubChainPairs(extremeWrithes_top8000, bins=100, color = 'blue', elevAzimList = [[15,-45],[160,-45],[-160,-130]], addSubplotToFirstPlot_b = 1, addSubPlotNr = 121, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  15)


#unrestricted search, lgth 30:
dict_top8000, extremeWrithes_top8000 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth30.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
git.histoSubChainPairs(extremeWrithes_top8000, bins=100, color = 'blue',elevAzimList = [[15,-45],[160,-45],[-160,-130]], addSubplotToFirstPlot_b = 1, addSubPlotNr = 122, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth = 30)


******************************************************
* Try out using higher order inv's instead of writhe
******************************************************

import gitP as git

root = r'C:\Users\Christian\Bioinformatics\papers\geometric_anomalies_in_folds\c_code\results_C_3rd\double_precision'

#top100 set:

invName = 'I1234'
invName = 'I1234_full'
invName = 'I1324'
invName = 'I1324_full'
invName = 'I1423'
invName = 'I152436'

(obs: see also "#Call to multi-invariant version" below)

#top100
#unrestricted search, lgth 15:
dict_top100, extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth15_' + invName + '.txt', linksAndPokes_b = 0, minMaxOut_b = 1)
###no min max: 
###dict_top100, extremeWrithes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth15_' + invName + '.txt', linksAndPokes_b = 0)

#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
extremeWrithes_top100 = extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max
git.histoSubChainPairs_minMax(extremeWrithes_top100, bins=100, color = 'lightgrey', addSubplotToFirstPlot_b = 0, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  15)



#unrestricted search, lgth 30:
dict_top100, extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth30_' + invName + '.txt', linksAndPokes_b = 0, minMaxOut_b = 1)
###no min max: 
###dict_top100, extremeWrithes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth30_' + invName + '.txt', linksAndPokes_b = 0)

#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
extremeWrithes_top100 = extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max
git.histoSubChainPairs_minMax(extremeWrithes_top100, bins=100, color = 'lightgrey', addSubplotToFirstPlot_b = 0, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  15)



#top8000
#unrestricted search, lgth 15:
dict_top8000, extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max  = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth15_' + invName + '.txt', linksAndPokes_b = 0, minMaxOut_b = 1)
###no min max:
###dict_top8000, extremeWrithes_top8000 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth15_' + invName + '.txt', linksAndPokes_b = 0)
 
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
#convert structure and then plot
extremeWrithes_top8000 = extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max
git.histoSubChainPairs_minMax(extremeWrithes_top8000, bins=100, color = 'lightgrey', addSubplotToFirstPlot_b = 1, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  15)

#unrestricted search, lgth 30:
dict_top8000, extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max  = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth30_' + invName + '.txt', linksAndPokes_b = 0, minMaxOut_b = 1)
###no min max:
###dict_top8000, extremeWrithes_top8000 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth30_' + invName + '.txt', linksAndPokes_b = 0)
 
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
#convert structure and then plot
extremeWrithes_top8000 = extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max
git.histoSubChainPairs_minMax(extremeWrithes_top8000, bins=100, color = 'lightgrey', addSubplotToFirstPlot_b = 1, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  15)

###############
#Call to multi-invariant version: 
###############

#top100

#unrestricted search, lgth 15:
order = "order2"
invList = ('I1234', 'I1324', 'I1423')
dict_top100, extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max  = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_' + order + '_minmax_global_mem_alloc_top100_lgth15.txt', linksAndPokes_b = 0, minMaxOut_b = 1, multiInv_b = 1)
chainsFile = root + r'\chains_closed_loops_top100.txt'
#convert structure and then plot
extremeWrithes_top100 = extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max
git.multiInv_histoSubChainPairs_minMax(invList = invList, subChainPairValues = extremeWrithes_top100, rowsFig = 1, colsFig = 3, logCount_b = 0, bins=100, colorMin = 'lightblue', colorMax = 'blue', distrPlotsInSubPlots_b =0, few_xticks_b =1, plotTop10_b = 0, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  15)

#unrestricted search, lgth 30:
order = "order2"
invList = ('I1234', 'I1324', 'I1423')
dict_top100, extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max  = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_' + order + '_minmax_global_mem_alloc_top100_lgth30.txt', linksAndPokes_b = 0, minMaxOut_b = 1, multiInv_b = 1)
chainsFile = root + r'\chains_closed_loops_top100.txt'
#convert structure and then plot
extremeWrithes_top100 = extremeWrithes_top100_All, extremeWrithes_top100_Min, extremeWrithes_top100_Max
git.multiInv_histoSubChainPairs_minMax(invList = invList, subChainPairValues = extremeWrithes_top100, rowsFig = 1, colsFig = 3, logCount_b = 0, bins=100, colorMin = 'lightblue', colorMax = 'blue', distrPlotsInSubPlots_b =0, few_xticks_b =1, plotTop10_b = 0, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth =  30)

#top8000

#unrestricted search, lgth 15:
order = "order2"
invList = ('I1234', 'I1324', 'I1423')
#order = "order3"
dict_top8000, extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max  = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_' + order + '_minmax_global_mem_alloc_top8000_lgth15.txt', linksAndPokes_b = 0, minMaxOut_b = 1, multiInv_b = 1)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
#convert structure and then plot
extremeWrithes_top8000 = extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max
git.multiInv_histoSubChainPairs_minMax(invList = invList, subChainPairValues = extremeWrithes_top8000, rowsFig = 1, colsFig = 3, logCount_b = 1, bins=100, colorMin = 'lightblue', colorMax = 'blue', distrPlotsInSubPlots_b =0, few_xticks_b = 1, plotTop10_b = 0, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, subChainLgth =  15)


#unrestricted search, lgth 30:
order = "order2"
invList = ('I1234', 'I1324', 'I1423')
dict_top8000, extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max  = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_' + order + '_minmax_global_mem_alloc_top8000_lgth30.txt', linksAndPokes_b = 0, minMaxOut_b = 1, multiInv_b = 1)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
#convert structure and then plot
extremeWrithes_top8000 = extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max
git.multiInv_histoSubChainPairs_minMax(invList = invList, subChainPairValues = extremeWrithes_top8000, rowsFig = 1, colsFig = 3, logCount_b = 1, bins=100, colorMin = 'lightblue', colorMax = 'blue', distrPlotsInSubPlots_b =0, lowerPercentile = 0.001, upperPercentile = 0.999, few_xticks_b =1, plotTop4_b = 1,  plotTop10_b = 0, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0,  subChainLgth =  30)

order = "order3"
invList = ''
invList = ('I1234', 'I1324', 'I1423', 'I132645','I152634', 'I162534')
invList = ( 'I132645','I152634', 'I162534')
invList = ('I1234', 'I1324', 'I1423',	
           #/*order 3*/
			"I123456",
			"I123546",
			"I123645",
			#/*13*/
			"I132456",
			"I132546",
			"I132645",
			#/*14*/
			"I142356",
			"I142536",
			"I142635",
			#/*15*/
			"I152346",
			"I152436",
			"I152634",
			#/*16*/
			"I162345",
			"I162435",
			"I162534")
invList = ('I132645', 'I123456')
dict_top8000, extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max  = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_' + order + '_minmax_global_mem_alloc_top8000_lgth30.txt', linksAndPokes_b = 0, minMaxOut_b = 1, multiInv_b = 1)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
#convert structure and then plot
extremeWrithes_top8000 = extremeWrithes_top8000_All, extremeWrithes_top8000_Min, extremeWrithes_top8000_Max

git.multiInv_histoSubChainPairs_minMax(invList = invList ,subChainPairValues = extremeWrithes_top8000, rowsFig = 2, colsFig = 3, logCount_b = 1, bins=100, color = 'lightgrey', colorMin = 'lightgrey', colorMax = 'lightgrey', alphaMin = 0.3, alphaMax = 0.6, distrPlotsInSubPlots_b =0, lowerPercentile = 0.001, upperPercentile = 0.999, few_xticks_b = 1, plotTop4_b = 0, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)\.([\w]+)', subChainLgth =  30)



###########################################################
## Comparisons vs Peter's code
###########################################################

import gitP as git

root = r'C:\Users\Christian\Bioinformatics\papers\geometric_anomalies_in_folds\c_code\results_C_3rd\double_precision'

#top100

#global, final values
cg_order2 = git.readResultsCcode(root + r'\Ivalues_order_2_global_mem_alloc_top100.txt', all_b = 0, pert_b =1)
pr_order2 = git.readResultsCcode(root + r'\Ivalues_order_2_peters_top100.txt', all_b = 0, pert_b =1)
#Compare:
sortedFileList = git.compareCresultsFinal(cg_order2, pr_order2, title = 'pr vs cg, order 2, top100', labelRelDiff = 'max rel diff (Db)', labelVal = 'abs val (Dbl)', use_log_b = 1)
#for the log10-version: set use_log_b = 1.

#.. not-log:
sortedFileList = git.compareCresultsFinal(cg_order2, pr_order2, sortedFileList = sortedFileList, title = 'pr vs cg, top100', labelRelDiff = 'max rel diff (Dbl)', labelVal = 'abs val (Dbl)', use_log_b = 0)

#order 3
cg_order3 = git.readResultsCcode(root + r'\Ivalues_order_3_global_mem_alloc_top100.txt', all_b = 0, pert_b =1)
pr_order3 = git.readResultsCcode(root + r'\Ivalues_order_3_peters_top100.txt', all_b = 0, pert_b =1)
#Compare:
sortedFileList = git.compareCresultsFinal(cg_order3, pr_order3, title = 'pr vs cg, order 3, top100', labelRelDiff = 'max rel diff (Dbl)', labelVal = 'abs val (Dbl)', use_log_b = 1)
#for the log10-version: set use_log_b = 1.

#.. while this is new:
sortedFileList = git.compareCresultsFinal(cg_order3, pr_order3, sortedFileList = sortedFileList, title = 'pr vs cg, order 3, top100', labelRelDiff = 'max rel diff (Dbl)', labelVal = 'abs val (Dbl)', use_log_b = 0)



#Showing instead the differences per invariant:
#top100
#order 2:
means, stdDevs = git.compareCresultsFinal_2(cg_order2, pr_order2, titleDistr = 'pr vs cg, order 2, top100', title = 'top100',titleRatios = 'top100') 
#to run medians etc use order 3:
#order 3:
means, stdDevs = git.compareCresultsFinal_2(cg_order3, pr_order3, titleDistr = 'pr vs cg, order 3, top100', title = 'order 3, top100',titleRatios = 'order 3, top100') 
medians, means, stdDevs, lP, uP = git.compareCresultsFinal_3(cg_order3, pr_order3, titleDistr = 'pr vs cg, order 3, top100', title = 'Statistics of relative diffs, pr vs cg, order 3, top100')


#hist2d/scatterplot comparison:
git.compareNormalizedInv(pr_order3, cg_order3, titleDistr = 'Distribution of normalized invariants', fontSizeDistr = "large", titleDiffs = 'Differences vs values', fontSizeDiffs = "large", figSize = (14,10))


#top8000

###########################################################
## Generate simplex illustrations
###########################################################

import matplotlib.pyplot as plt
import numpy as np 

## basic simplex with recursion indications:

#the large simplex
#line segments:
i0 = 0
iN = 100
simplex1x = [i0 for i in range(i0,iN+1)]
simplex1y = range(i0,iN+1)
simplex2x = range(i0,iN+1)
simplex2y = [iN for i in range(i0,iN + 1)]
simplex3x = range(i0,iN + 1)
simplex3y = range(i0,iN + 1)
simplexx = simplex1x[:]
simplexx.extend(simplex2x)
simplexx.extend(simplex3x)
simplexy = simplex1y[:]
simplexy.extend(simplex2y)
simplexy.extend(simplex3y)

#plot simplex
fig = plt.figure()
ax = fig.add_subplot(111)
a = 0.75
plt.plot(simplexx, simplexy, 'g-', alpha = a)

#add smaller simplices to illustrate basis for recursion formulas:
#1
i0 = 50
iN = 90
simplex1x = [i0 for i in range(i0,iN+1)]
simplex1y = range(i0,iN+1)
simplex2x = range(i0,iN+1)
simplex2y = [iN for i in range(i0,iN + 1)]
simplex3x = range(i0,iN + 1)
simplex3y = range(i0,iN + 1)
simplexx = simplex1x[:]
simplexx.extend(simplex2x)
simplexx.extend(simplex3x)
simplexy = simplex1y[:]
simplexy.extend(simplex2y)
simplexy.extend(simplex3y)

plt.plot(simplexx, simplexy, 'b-', alpha = a)

plt.text(53,85, '1')

#2
i0 = 50
iN = 80
simplex1x = [i0 for i in range(i0,iN+1)]
simplex1y = range(i0,iN+1)
simplex2x = range(i0,iN+1)
simplex2y = [iN for i in range(i0,iN + 1)]
simplex3x = range(i0,iN + 1)
simplex3y = range(i0,iN + 1)
simplexx = simplex1x[:]
simplexx.extend(simplex2x)
simplexx.extend(simplex3x)
simplexy = simplex1y[:]
simplexy.extend(simplex2y)
simplexy.extend(simplex3y)

plt.plot(simplexx, simplexy, 'b-', alpha = a)

plt.text(53,75, '2')

#3
i0 = 60
iN = 90
simplex1x = [i0 for i in range(i0,iN+1)]
simplex1y = range(i0,iN+1)
simplex2x = range(i0,iN+1)
simplex2y = [iN for i in range(i0,iN + 1)]
simplex3x = range(i0,iN + 1)
simplex3y = range(i0,iN + 1)
simplexx = simplex1x[:]
simplexx.extend(simplex2x)
simplexx.extend(simplex3x)
simplexy = simplex1y[:]
simplexy.extend(simplex2y)
simplexy.extend(simplex3y)

plt.plot(simplexx, simplexy, 'b-', alpha = a)

plt.text(63,85, '3')

plt.text(63,75, '4')

#add arrows to indicate recursion directions:
ax.arrow(30, 70, -9, 0, head_width=1 , head_length=1, fc='k', ec='k')
ax.arrow(20, 70, 0, 9, head_width=1 , head_length=1, fc='k', ec='k')

#make ticks (this will also extend the fig so all edges of the simplex can be seen):
arange = np.arange(-10, 111, step = 10)
xticks = plt.xticks(arange)
yticks = plt.yticks(arange)


plt.title('Sample simplex, recursion indications')
plt.show()

##Simplex w perturbation band:

#the large simplex
#line segments:
i0 = 0
iN = 100
simplex1x = [i0 for i in range(i0,iN+1)]
simplex1y = range(i0,iN+1)
simplex2x = range(i0,iN+1)
simplex2y = [iN for i in range(i0,iN + 1)]
simplex3x = range(i0,iN + 1)
simplex3y = range(i0,iN + 1)
simplexx = simplex1x[:]
simplexx.extend(simplex2x)
simplexx.extend(simplex3x)
simplexy = simplex1y[:]
simplexy.extend(simplex2y)
simplexy.extend(simplex3y)

#plot simplex
fig = plt.figure()
ax = fig.add_subplot(111)
a = 0.75
plt.plot(simplexx, simplexy, 'g-', alpha = a)

#upper boundary of perturbation band
line1x = [0 for i in range(40,51)]
line1y = [i for i in range(40,51)]
line2x = [i for i in range(0,41)]
line2y = [50 for i in range(0,41)]
line3x = [40 for i in range(50,101)]
line3y = [i for i in range(50,101)]
line4x = [i for i in range(40,51)]
line4y = [100 for i in range(40,51)]
upperx = line1x[:]
upperx.extend(line2x)
upperx.extend(line3x)
upperx.extend(line4x)
uppery = line1y[:]
uppery.extend(line2y)
uppery.extend(line3y)
uppery.extend(line4y)

#lower boundary
line5x = [i for i in range(0,41)]
line5y = [40 for i in range(0,41)]
line6x = [i for i in range(40,51)]
line6y = [i for i in range(40,51)]
line7x = [50 for i in range(50,101)]
line7y = [i for i in range(50,101)]
lowerx = line5x[:]
lowerx.extend(line5x)
lowerx.extend(line6x)
lowerx.extend(line7x)
lowery = line5y[:]
lowery.extend(line5y)
lowery.extend(line6y)
lowery.extend(line7y)

#plot simplex
fig = plt.figure()
ax = fig.add_subplot(111)
a = 0.75
color = 'orange'
plt.plot(simplexx, simplexy, 'g-', alpha = a)
#plot boundaries
plt.plot(lowerx, lowery, color, alpha = a)
plt.plot(upperx, uppery, color, alpha = a)

#fill the band
plt.fill_between(line5x, line2y, line5y, color = color, alpha = a)
plt.fill_between(line6x, line4y, line6y, color = color, alpha = a)

#make ticks (this will also extend the fig so all edges of the simplex can be seen):
arange = np.arange(-10, 111, step = 10)
xticks = plt.xticks(arange)
yticks = plt.yticks(arange)

plt.title('Sample simplex with highlighted perturbation band')
plt.show()

####################
### Simplex illu for Section meeting talk
####################

import matplotlib.pyplot as plt
import numpy as np 

## basic simplex with recursion indications:

#the large simplex
#line segments:
i0 = 0
iN = 100
simplex1x = [i0 for i in range(i0,iN+1)]
simplex1y = range(i0,iN+1)
simplex2x = range(i0,iN+1)
simplex2y = [iN for i in range(i0,iN + 1)]
simplex3x = range(i0,iN + 1)
simplex3y = range(i0,iN + 1)
simplexx = simplex1x[:]
simplexx.extend(simplex2x)
simplexx.extend(simplex3x)
simplexy = simplex1y[:]
simplexy.extend(simplex2y)
simplexy.extend(simplex3y)

#make ticks (this will also extend the fig so all edges of the simplex can be seen):
arange = np.arange(-10, 111, step = 10)
xticks = plt.xticks(arange)
yticks = plt.yticks(arange)

#plot simplex
fig = plt.figure()
ax = fig.add_subplot(111)
a = 0.75
plt.plot(simplexx, simplexy, 'g-', alpha = a)

plt.show()

#by meshgrid
i0 = -1
iN = 25
x = np.arange(i0+1,iN+1,1)
y = np.arange(i0+1,iN+1,1)
#x,y = np.meshgrid(x,y)
for i in range(iN+1):
    for j in range(iN+1):
        if i < j: #x[i]<x[j]:
            plt.scatter(i,j) #plt.scatter(x[i],x[j])
plt.xticks(np.arange(i0+6,iN+1,step= 5))
plt.yticks(np.arange(i0+6,iN+1,step= 5))
plt.show()

##########################################################
## Plots for presentation KUpostDocMay2016
##########################################################

root = r'C:\Users\Christian\Bioinformatics\papers\geometric_anomalies_in_folds\c_code\results_C_1st\double_precision'

#top100 set:
dict_top100, links_top100, pokes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\ClosedLoopChars_global_mem_alloc_top100.txt')
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
pokeLength = 10
titlePokesI = 'Distribution of writhe values for potential pokes of length ' + str(pokeLength) 
titlePokesII = ' writhe values for potential pokes of length ' + str(pokeLength) 

#top8000 set:
dict_top8000_cl, links_top8000_cl, pokes_top8000_cl = git.readResultsCcode_closedLoops(inputFileName = root + r'\ClosedLoopChars_global_mem_alloc_top8000_pokeLgth10.txt')
pokeLength = 10
chainsFile = root + r'\chains_closed_loops_top8000.txt'
titlePokesI = 'Distribution of writhe values for potential pokes of length ' + str(pokeLength) 
titlePokesII = ' writhe values for potential pokes of length ' + str(pokeLength) 

#distr for links in top100 and top8000 in one plot:
fig = git.histoClosedLoops(links_top100, pokes_top100, binsLinks=100, binsPokes =100, titleDistr = 'top100', linksDistrPlotsInSubPlots_b =2, pokesDistrPlotsInSubPlots_b =1, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 11, cutOffPokes = 21.0, plotPokeExs_b = 0, plotLinkExs_b = 1, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-12.9,-12.7], intervalPokes = [-5,-4], plotExsInInt_b =0, titlePokesI = titlePokesI, titlePokesII = titlePokesII, stringPattern = '([\S]+)([^\w]+)([\w]+)')
fig = git.histoClosedLoops(links_top8000_cl, pokes_top8000_cl, binsLinks=100, binsPokes =100, inFig = fig, titleDistr = 'top8000', linksDistrPlotsInSubPlots_b =2, pokesDistrPlotsInSubPlots_b =1, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 12.4, cutOffPokes = 21.0, plotPokeExs_b = 0, plotLinkExs_b = 0, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, intervalLinks = [-12.9,-12.7], intervalPokes = [-5,-4], plotExsInInt_b =0, titlePokesI = titlePokesI, titlePokesII = titlePokesII)

#examples from just one angle:
fig = git.histoClosedLoops(links_top100, pokes_top100, binsLinks=100, binsPokes =100, titleDistr = 'top100', linksDistrPlotsInSubPlots_b =2, pokesDistrPlotsInSubPlots_b =1, linksShowTopNr = 100, pokesShowTopNr = 100, cutOffLinks = 11, cutOffPokes = 21.0, plotPokeExs_b = 0, plotLinkExs_b = 1, plotLinksTop10_b = 0, plotPokesTop10_b = 0, chainsFile = chainsFile, elevAzimList = [[0,0]], intervalLinks = [-12.9,-12.7], intervalPokes = [-5,-4], plotExsInInt_b =0, titlePokesI = titlePokesI, titlePokesII = titlePokesII, stringPattern = '([\S]+)([^\w]+)([\w]+)')

#unrestricted search, lgth 30:

#top100:
dict_top100, extremeWrithes_top100 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top100_lgth30.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top100.txt'
git.histoSubChainPairs(extremeWrithes_top100, bins=100, color = 'blue', addSubplotToFirstPlot_b = 0, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth = 30)


#top8000:
dict_top8000, extremeWrithes_top8000 = git.readResultsCcode_closedLoops(inputFileName = root + r'\SubChainPairChars_minmax_global_mem_alloc_top8000_lgth30.txt', linksAndPokes_b = 0)
#sort the results and plot distributions (parts of) and link/poke examples (outside cutOff's)
chainsFile = root + r'\chains_closed_loops_top8000.txt'
git.histoSubChainPairs(extremeWrithes_top8000, color = 'blue', logCount_b = 0, bins=100, addSubplotToFirstPlot_b = 0, distrPlotsInSubPlots_b =0, showTopNr = 10, plotTop10_b = 1, plotExs_b = 0, chainsFile = chainsFile, interval = [-4,-3], plotExsInInt_b =0, stringPattern = '([\S]+)([^\w]+)([\w]+)', subChainLgth = 30)



'''


import csv

def readResultsCcode(inputFileName, all_b = 0, pert_b = 0, pertTime_b = 0):
    '''Read in results obtained by running the C-version.
    Output: dictionary mapping each pdb-file (name given as path to file) to
    dictionary containing the results: this dictionary maps each measure name
    (e.g. I12) to list: length of chain, degree of computation (= 2*order), 
    cpu-time for computation, measure value.
    Input:
    all_b: read in from and "all values file" (whole simplex results) or only
    from "final corner".
    pert_b: read in from file with non-perturbation results (0) or from a 
    file containing perturbation results. (in non-pert files the perturbation
    number will be 1, viz the unperterubed case).'''

    readOrder = ["PDB_file", 
                 "chain id",
                 "chain nr",
		    "str length",
			"order",
           "perturbation number",
			"cpuTime",
           "aggrTime",
			#/*order 1*/
			"I12", 
			"Ia12",
			#/*order 2:*/
			"I1234", 
			"I1324", 
			"I1423",
			#/*abs value versions:*/
			#/*12*/
			"Ia1234",
			"I12a34",
			"Ia12a34",
			#/*13*/
			"Ia1324",
			"I13a24",
			"Ia13a24",
			#/*14*/
			"Ia1423",
			"I14a23",
			"Ia14a23",
			#/*order 3*/
			"I123456",
			"I123546",
			"I123645",
			#/*13*/
			"I132456",
			"I132546",
			"I132645",
			#/*14*/
			"I142356",
			"I142536",
			"I142635",
			#/*15*/
			"I152346",
			"I152436",
			"I152634",
			#/*16*/
			"I162345",
			"I162435",
			"I162534"]

    outDict = {}

    if (all_b != 1 and pert_b !=1) :
        with open(inputFileName, 'r') as CresultsTxt:
            Cresults = csv.reader(CresultsTxt, delimiter = ';')
            for row in Cresults:
                if not(outDict.has_key(row[0])): #row[0] is the file name/path to file
                    outDict[row[0]] = {}
                if not(outDict[row[0]].has_key(row[1])): #row[1] is the chain id
                    outDict[row[0]][row[1]] = {}
                    for i in range(len(readOrder))[3:]:
                        if i!=5: #skip the pertubation number
                            #print readOrder[i]
                            name = readOrder[i]
                            outDict[row[0]][row[1]][name] = row[i]
                        i +=1
                        
                        
    if (all_b != 1 and pert_b == 1):
        with open(inputFileName, 'r') as CresultsTxt:
            Cresults = csv.reader(CresultsTxt, delimiter = ';')
            rowCnt = 0
            for row in Cresults: 
#                print row
                if rowCnt >0:#skip first "header" row
                    if not(outDict.has_key(row[0])): #row[0] is the file name/path to file
                        outDict[row[0]] = {}
                    if not(outDict[row[0]].has_key(row[1])): #row[1] is the chain id
                        outDict[row[0]][row[1]] = {}
                        if not(outDict[row[0]][row[1]].has_key(row[5])): #row[5] is the perturbation number
                            outDict[row[0]][row[1]][row[5]] = {}
                            if pertTime_b ==0:
                                #read length of structure in separately:
                                outDict[row[0]][row[1]][row[5]]["str length"] = int(row[3])
                                for i in range(len(readOrder))[8:]: #skip info in the first col's
                                    #print "row[%d]:%s" % (i, str(row[i]))
                                    name = readOrder[i]
                                    outDict[row[0]][row[1]][row[5]][name] = row[i]
                                    i +=1
                            elif pertTime_b ==1:
                                for i in range(len(readOrder))[:8]: #for pert time there's only info in first 6 col's
                                    #print "row[%d]:%s" % (i, str(row[i]))
                                    name = readOrder[i]
#                                    print name
                                    outDict[row[0]][row[1]][row[5]][name] = row[i]
                                    i +=1
                    else: #to cover the case of repeated runs of the code
                        if pertTime_b ==1:
                            outDict[row[0]][row[1]][row[5]]['cpuTime'] = float(outDict[row[0]][row[1]][row[5]]['cpuTime'])+ float(row[6])
                            outDict[row[0]][row[1]][row[5]]['aggrTime'] = float(outDict[row[0]][row[1]][row[5]]['aggrTime']) + float(row[7])
                        
                rowCnt += 1
                
                
    elif (all_b == 1 and pert_b != 1):
        strName = ''
        with open(inputFileName, 'r') as CresultsTxt:
            Cresults = csv.reader(CresultsTxt, delimiter = ';')
            for row in Cresults:
#                print row
                if row[0] == 'Header':
                    strName = row[1]
#                    print strName
                    chainId = row[2]
                    if not(outDict.has_key(strName)):
                        outDict[strName] = {}
                    if not(outDict[strName].has_key(chainId)):
                        outDict[strName][chainId] = {}
                elif row[0] == 'Values':
                    readOrder = []
                    for k in range(len(row))[3:]:
                        readOrder.append(row[k])
                        outDict[strName][chainId][row[k]] = {}
#                    print readOrder
                elif row[0] == 'Line type':
                    continue 
                else:
                    #read results:
                    i = int(row[0])
                    j = int(row[1])
                    for k in range(len(readOrder)):
                        outDict[strName][chainId][readOrder[k]][(i,j)] = float(row[k+2])
                        
                        
    elif (all_b == 1 and pert_b == 1):
        strName = ''
        pertNr = -1
        with open(inputFileName, 'r') as CresultsTxt:
            Cresults = csv.reader(CresultsTxt, delimiter = ';')
            for row in Cresults:
#                print row
                if row[0] == 'Header':
                    strName = row[1]
                    chainId = row[2]
                    pertNr = row[4]
#                    print strName
                    if not(outDict.has_key(strName)):
                        outDict[strName] = {}
                    if not(outDict[strName].has_key(chainId)):
                        outDict[strName][chainId] = {}
                    if not(outDict[strName][chainId].has_key(pertNr)):
                        outDict[strName][chainId][pertNr] = {}
                elif row[0] == 'Values':
                    readOrder = []
                    for k in range(len(row))[3:]:
                        readOrder.append(row[k])
                        if not(outDict[strName][chainId][pertNr].has_key(row[k])):
                            outDict[strName][chainId][pertNr][row[k]] = {}
#                    print readOrder
                elif row[0] == 'Line type':
                    continue 
                else: 
                    #read results:
                    i = int(row[0])
                    j = int(row[1])
                    for k in range(len(readOrder)):
                        outDict[strName][chainId][pertNr][readOrder[k]][(i,j)] = float(row[k+2])
        
    return outDict
            
def getResultsPcode(path = PDBtop100Folder, fileName = r"C:\Users\Christian\Bioinformatics\projects\knots\code\top100H_source\top100names.txt", 
                    degree = 4, all_b = 0, global_b = 0, full_b = 0, repeats = 1):
    '''Computation of measures up to specified degree across all names 
    in input list; non/vectorized array/dict versions.
    Input: 
    path and filename: as in getPDBlist.
    all_b: return results across entire simplex (1) or just the final measure, i.e.
    the value in the "final" corner of the simplex (0).
    global_b: whether to go for using the "global allocation mode" (1) in which arrays for
    storing the I-results are allocated and initialized up front sizing for the longest
    chain. Else (0) this is done per chain.
    full_b: extended order 2 measures (1) or not (0); only in use with global_b =1.
    repeats: number of repeats of the computations to run; only for all_b = 0
    Output: dictionary as output from readResultsCcode''' 
    
    readOrder = ["w", "I12", "I1234", "I1234_full", "I1324", "I1324_full", 
    "I1324_full2", "I1423", "I1423_full0", "I1423_full2", "I123456", "I123546", 
    "I123645", "I132456", "I132546", "I132645", "I142356", "I142536", "I142635", 
    "I152346", "I152436", "I152634", "I162345", "I162435", "I162534"]
    
    degree0 = ["w"]
    degree2 =  ["w", "I12"]
    degree4 = ["w", "I12", "I1234", "I1324", "I1423"]
    degree4_full = ["w", "I12", "I1234", "I1234_full", "I1324", "I1324_full", 
    "I1324_full2", "I1423", "I1423_full0", "I1423_full2"]
    degree6 = readOrder
        
    
    outDict = {}
    
    #get list of PDB file names:
    PDBlist = getPDBlist(path = path, fileName = fileName)
    
    if global_b == 0:
        for l_item in PDBlist:
            print l_item[1]
            f = path + l_item[1]
            print f
            
            cpuTime = {} #for holding the cpuTime per sub-chain
            if all_b != 1:
                for i in range(repeats):
                    print "I'm at repeat %d" % i
                    #compute measures using the array based, vectorized version: 
                    start = time.time()
                    Idict, CaChain, pChain, aggrTime = I(PDBfilename =  f, full_b = 0, print_b=0, degree = degree ,vect_b = 1, aggrTime_b = 1)
                    end  = time.time()
                    cpuTimeTotal = end - start #the time consumption for handling (incl loading) the whole structure, ie all its sub-chains
                    #We want to divide the cpuTimeTotal on the sub-chains; it seems reasonable to do this pro-rata on the length 
                    #for the laoding and mem allocation part. So we first subtract the total aggr time and then divide out:
                    aggrTimeTotal= 0
                    strL = 0                    
                    for chainId in CaChain.keys():
                        aggrTimeTotal += aggrTime[chainId]
                        strL += len(CaChain[chainId])
                    for chainId in CaChain.keys():
                        cpuTime[chainId] = float(cpuTimeTotal - aggrTimeTotal)*float(len(CaChain[chainId]))/strL + aggrTime[chainId]
                    #record results:
                    for chainId in CaChain.keys():
                        L = len(CaChain[chainId])
                        if not(outDict.has_key(f)): #f = path to pdb file
                            outDict[f] = {}
                        if not(outDict[f].has_key(chainId)): 
                            outDict[f][chainId] = {}
                        
                            outDict[f][chainId]['str length'] = L
                            outDict[f][chainId]['order'] = int(float(degree)/2)
                            outDict[f][chainId]['cpuTime'] = float(cpuTime[chainId])/repeats
                            outDict[f][chainId]['aggrTime'] = float(aggrTime[chainId])/repeats
                            k = 1
                            for k in range(len(readOrder)):
                                name = readOrder[k]
                                outDict[f][chainId][name] = Idict[chainId][k][0,L-2]
                        elif outDict[f].has_key(chainId): 
                            outDict[f][chainId]['cpuTime'] += float(cpuTime[chainId])/repeats
                            outDict[f][chainId]['aggrTime'] += float(aggrTime[chainId])/repeats
            elif all_b ==1:
                #compute measures using the dictionary based, vectorized version: 
                Idict, CaChain, pChain = I_dict(PDBfilename =  f, full_b = 0, print_b=0, degree = degree ,vect_b = 1)
                #record results:
                for chainId in CaChain.keys():
                    if not(outDict.has_key(f)): #f = path to pdb file
                        outDict[f] = {}
                    if not(outDict[f].has_key(chainId)):
                            outDict[f][chainId] = {}
                            outDict[f][chainId] = Idict[chainId]
        #                    k = 1
        #                    for k in range(len(readOrder)):
        #                        name = readOrder[k]
        #                        outDict[f][name] = Idict[name]
                        
        
    elif global_b ==1:    
        
        #very first, assign "defualt values":
        wDict = 0
        I12 = 0 #value of I12 on full chain
        I12Dict = 0
        I1234Dict = 0
        I1234Dict_full = 0
        I1234Dict_full_aid = 0
        I1234Dict_full2 = 0
        I1234Dict_full2_aid = 0
        I1423Dict = 0 
        I1423Dict_full0 = 0 
        I1423Dict_full2 = 0
        I1423Dict_full2_aid = 0
        I1324Dict = 0
        I1324Dict_full = 0
        I1324Dict_full_aid = 0
        I1324Dict_full2 = 0
        I1324Dict_full2_aid = 0
        #(12)*:
        I123456Dict = 0 
        I123645Dict = 0
        I123546Dict = 0
        #(13)*:
        I132456Dict = 0
        I132546Dict = 0
        I132645Dict = 0
        #(14)*
        I142356Dict = 0
        I142536Dict = 0
        I142635Dict = 0
        #(15)*
        I152346Dict = 0
        I152436Dict = 0
        I152634Dict = 0
        #(16)*
        I162345Dict = 0
        I162435Dict = 0
        I162534Dict = 0
            
        #first find the length of the longest chain and allocate for this size:
        maxLength = 0
        for l_item in PDBlist:
            L = l_item[0]
            if L > maxLength:
                maxLength= L
        L = maxLength
        #allocate and init:
        if degree >=0:
            wDict = np.zeros((L+1,L+1))
        if degree >= 2:
            I12 = 0 #value of I12 on full chain
            I12Dict = np.zeros((L+1,L+1))
        if degree >= 4:
    
            I1234Dict = np.zeros((L+1,L+1))
            I1324Dict = np.zeros((L+1,L+1))
            I1423Dict = np.zeros((L+1,L+1)) 

            if (full_b ==1 and degree <6):
                I1234Dict_full = np.zeros((L+1,L+1))    
                I1234Dict_full_aid = np.zeros((L+1,L+1))
                I1324Dict_full = np.zeros((L+1,L+1))
                I1324Dict_full_aid = np.zeros((L+1,L+1))

        if degree >= 6:

            I1234Dict_full = np.zeros((L+1,L+1))
            I1234Dict_full_aid = np.zeros((L+1,L+1))
            I1234Dict_full2 = np.zeros((L+1,L+1))
            I1234Dict_full2_aid = np.zeros((L+1,L+1))
            I1423Dict_full0 = np.zeros((L+1,L+1)) 
            I1423Dict_full2 = np.zeros((L+1,L+1))
            I1423Dict_full2_aid = np.zeros((L+1,L+1))
            I1324Dict_full = np.zeros((L+1,L+1))
            I1324Dict_full_aid = np.zeros((L+1,L+1))
            I1324Dict_full2 = np.zeros((L+1,L+1))
            I1324Dict_full2_aid = np.zeros((L+1,L+1))       

            #(12)*:
            I123456Dict = np.zeros((L+1,L+1)) 
            I123645Dict = np.zeros((L+1,L+1))
            I123546Dict = np.zeros((L+1,L+1))
            #(13)*:
            I132456Dict = np.zeros((L+1,L+1))
            I132546Dict = np.zeros((L+1,L+1))
            I132645Dict = np.zeros((L+1,L+1))
            #(14)*
            I142356Dict = np.zeros((L+1,L+1))
            I142536Dict = np.zeros((L+1,L+1))
            I142635Dict = np.zeros((L+1,L+1))
            #(15)*
            I152346Dict = np.zeros((L+1,L+1))
            I152436Dict = np.zeros((L+1,L+1))
            I152634Dict = np.zeros((L+1,L+1))
            #(16)*
            I162345Dict = np.zeros((L+1,L+1))
            I162435Dict = np.zeros((L+1,L+1))
            I162534Dict = np.zeros((L+1,L+1))
                      
        L = 0 #reset for safety
        
        for l_item in PDBlist:
            for i in range(repeats):
                print "I'm at repeat %d" % i
                print l_item[1]
                f = path + l_item[1]
                print f    
                cpuTime = {} #for holding the cpuTime per sub-chain
    
                if all_b != 1:
                    #compute measures using the array based, vectorized version: 
                    start = time.time()
                    Idict, CaChain, pChain, aggrTime =I_naked(PDBfilename =  f, 
                                                              full_b = 0, 
                                                              print_b=0, 
                                                              degree = degree,
                                                              aggrTime_b = 1, 
                                                              wDict = wDict,
                                                              I12 = I12, #value of I12 on full chain
                                                              I12Dict = I12Dict,
                                                              I1234Dict = I1234Dict,
                                                              I1234Dict_full = I1234Dict_full,
                                                              I1234Dict_full_aid = I1234Dict_full_aid,
                                                              I1234Dict_full2 = I1234Dict_full2,
                                                              I1234Dict_full2_aid = I1234Dict_full2_aid,
                                                              I1423Dict = I1423Dict, 
                                                              I1423Dict_full0 = I1423Dict_full0, 
                                                              I1423Dict_full2 = I1423Dict_full2,
                                                              I1423Dict_full2_aid = I1423Dict_full2_aid,
                                                              I1324Dict = I1324Dict,
                                                              I1324Dict_full = I1324Dict_full,
                                                              I1324Dict_full_aid = I1324Dict_full_aid,
                                                              I1324Dict_full2 = I1324Dict_full2,
                                                              I1324Dict_full2_aid = I1324Dict_full2_aid,
                                                              #(12)*:
                                                              I123456Dict = I123456Dict, 
                                                              I123645Dict = I123645Dict,
                                                              I123546Dict = I123546Dict,
                                                              #(13)*:
                                                              I132456Dict = I132456Dict,
                                                              I132546Dict = I132546Dict,
                                                              I132645Dict = I132645Dict,
                                                              #(14)*
                                                              I142356Dict = I142356Dict,
                                                              I142536Dict = I142536Dict,
                                                              I142635Dict = I142635Dict,
                                                              #(15)*
                                                              I152346Dict = I152346Dict,
                                                              I152436Dict = I152436Dict,
                                                              I152634Dict = I152634Dict,
                                                              #(16)*
                                                              I162345Dict = I162345Dict,
                                                              I162435Dict = I162435Dict,
                                                              I162534Dict = I162534Dict)
                    end = time.time()
                    cpuTimeTotal = end - start #the time consumption for handling (incl loading) the whole structure, ie all its sub-chains
                    #We want to divide the cpuTimeTotal on the sub-chains; it seems reasonable to do this pro-rata on the length 
                    #for the laoding and mem allocation part. So we first subtract the total aggr time and then divide out:
                    aggrTimeTotal= 0
                    strL = 0
                    for chainId in CaChain.keys():
                        aggrTimeTotal += aggrTime[chainId]
                        strL += len(CaChain[chainId])
                    for chainId in CaChain.keys():
                        cpuTime[chainId] = float(cpuTimeTotal - aggrTimeTotal)*float(len(CaChain[chainId]))/strL + aggrTime[chainId]
                    #record results:
                    for chainId in CaChain.keys():
                        L = len(CaChain[chainId])
                        if not(outDict.has_key(f)): #f = path to pdb file
                            outDict[f] = {}
                        if not(outDict[f].has_key(chainId)): 
                            outDict[f][chainId] = {}                                                                               
                                        
                            outDict[f][chainId]['str length'] = L
                            outDict[f][chainId]['order'] = int(float(degree)/2)
                            outDict[f][chainId]['cpuTime'] = float(cpuTime[chainId])/repeats
                            outDict[f][chainId]['aggrTime'] = float(aggrTime[chainId])/repeats
                            k = 1
    #                        return Idict[k]
                            for k in range(len(readOrder)):
                                name = readOrder[k]
                                if degree == 2:
                                    if degree2.count(name) >0:
                                        outDict[f][chainId][name] = Idict[chainId][k][0,L-2]
                                elif (degree == 4 and full_b ==0):
                                    if degree4.count(name) >0:
                                        outDict[f][chainId][name] = Idict[chainId][k][0,L-2]
                                elif (degree == 4 and full_b ==1):
                                    if degree4_full.count(name) >0:
                                        outDict[f][chainId][name] = Idict[chainId][k][0,L-2]
                                else: # degree ==6: degree6 = readOrder so all is covered
                                    outDict[f][chainId][name] = Idict[chainId][k][0,L-2]
                        elif outDict[f].has_key(chainId):
                            outDict[f][chainId]['cpuTime'] += float(cpuTime[chainId])/repeats
                            outDict[f][chainId]['aggrTime'] += float(aggrTime[chainId])/repeats
                else:
                    print "Code only runs with all_b =0"

    return PDBlist, outDict


def compareResultsAll_PvsC_Dicts(Pdict, Cdict):
    '''For comparing the results across the whole simplex, C vs Python (P). Computes, among
    others, for each structure: 
    1) maximum abs-relative difference in value at the "final corner" of the simplex.
    This is measured as the abs value of the difference between the C-value and the 
    P-value divided by the abs value of the P-value. The max is taken over all invariants
    computed and all structures included. 
    2) the following two ratios:
        a) the maximum "diff-to-signal": maximum over all invariants of the 
        mean-of-absolute-differences divided by the mean-absolute-value (both means 
        taken over the vertices in the simplex) abs-relative difference of cumulative 
        abs values across the simplex.
        b) the maximum "noise-to-signal": maximum over all invariants of the 
        standard deviation-of-absolute-differences divided by the mean-absolute-value 
        (both means taken over the vertices in the simplex) abs-relative difference of 
        cumulative abs values across the simplex.
    
        In these the difference is (C-value minus P-value) and the value is taken to be 
        the P-value. While a) reports an average relative difference, it is informative
        to include the variation of the differences to, which is why b) is included.
        For both a) and b) the invariant at which the max occurs is reported. 
        
    Input: both dicts must map structure name to invariant name to vertex in simplex
    to value.'''
    cntStr = 0
    L = 0
    for strKey in Pdict.keys():
        for chainId in Pdict[strKey].keys():
            L = np.sqrt(len(Pdict[strKey][chainId]['w'].keys()))
            L = int(L)
            if Cdict.has_key(strKey):
                print "PDB: %s" % strKey
                print "Length:%d" % L
                cntStr += 1
                
                if Cdict[strKey].has_key(chainId):
                    print "Chain id: %s" % chainId
                    cntI = 0
        #             maxRelDiffAll = 0
        #            maxRelDiffAllAtInv = ''
        
                    maxRelDiffFinal = 0
                    maxFinalAtInv = ''
                    
                    maxDiffToSignal = 0                    
                    invAtMaxDTS = ''
                    maxNoiseToSignal = 0
                    invAtMaxNTS  = ''
                    for Iname in Pdict[strKey][chainId].keys():
                        if Cdict[strKey][chainId].has_key(Iname):
                            #print "***************Invariant: %s" % Iname
                            #comparison at final corner:
                            relDiffFinal = 0
                            relDiffFinal = 100*abs((float(Pdict[strKey][chainId][Iname][(0,L-2)]) - float(Cdict[strKey][chainId][Iname][(0,L-2)]))/Pdict[strKey][chainId][Iname][(0,L-2)])
                            #print "Rel diff final (pct): %f" % relDiffFinal
                            if relDiffFinal > maxRelDiffFinal:
                                maxRelDiffFinal = relDiffFinal
                                maxFinalAtInv = Iname
                            #comparison over whole simplex
                            cntI +=1
                            vals = []
                            diffs = []                    
                            
                            ijCnt = 0                                        
                            relDiff = 0
                            maxRelDiffInv = 0
                            absPCum = 0
                            absDiffCum = 0
                            relDiffCum = 0
                            minAbsCum = 1e-40
                            absDiffAtMinAbsCum  = 0
                            for ijKey in Pdict[strKey][chainId][Iname].keys():
                                i,j = ijKey
                                if j > i:
            #                        relDiff = abs(float(Pdict[strKey][Iname][ijKey] - Cdict[strKey][Iname][ijKey])/Pdict[strKey][Iname][ijKey])
            #                        if relDiff > maxRelDiffInv:
            #                            maxRelDiffInv = relDiff
                                    val = float(Pdict[strKey][chainId][Iname][ijKey])
                                    absVal = abs(val)
                                    diff = float(Pdict[strKey][chainId][Iname][ijKey]) - float(Cdict[strKey][chainId][Iname][ijKey])
                                    absDiff = abs(diff)
                                    absDiffCum += absDiff 
                                    #for additional info
                                    if (absVal < minAbsCum and absVal >1e-40):
                                        minAbsCum = absVal
                                        absDiffAtMinAbsCum = absDiff
            #                        if (abs(Pdict[strKey][Iname][ijKey] - Cdict[strKey][Iname][ijKey]) > maxAbsDiff):
            #                            maxAbsDiff = abs(Pdict[st rKey][Iname][ijKey] - Cdict[strKey][Iname][ijKey]) 
                                    ijCnt +=1
                                    #collect the values for computing mean and std dev:
                                    vals.append(absVal)
                                    diffs.append(absDiff)
                                
        #                    print "Cumulativefloat(absDiff)
        #                    print "For invariant: %s the max rel diff across the simplex is: %f" % maxRelDiffInv
        #                    if maxRelDiffInv > maxRelDiffAll:
        #                        maxRelDiffAll = maxRelDiffInv
        #                        maxRelDiffAllAtInv = Iname
                            vals = np.array(vals)
                            diffs = np.array(diffs)
                            mean = np.average(vals)
                            #print "mean absVal:%f" % mean
                            meanDiff = np.average(diffs)
                            #print "meanDiff:%f" % meanDiff
                            stdDiff = np.std(diffs)
                            diffToSignal = 100*float(meanDiff)/mean
                            noiseToSignal = 100*float(stdDiff)/mean
                            if diffToSignal > maxDiffToSignal:
                                maxDiffToSignal = diffToSignal
                                invAtMaxDTS = Iname
                            if noiseToSignal > maxNoiseToSignal:
                                maxNoiseToSignal = noiseToSignal
                                invAtMaxNTS = Iname
        
                            #print "absDiffCum avg:%f" % meanDiff2
        
        #                    print "No of vertices for this invariant: %d" % ijCnt
        #                    print "Cumulative abs value: %f" % absPCum
        #                    print "Cumulative abs difference: %f" % absDiffCum
        #                    print "Relative cumulative diff (pct) for %s is: %f" % (Iname, relDiffCum)
        #                    print "Min abs value (P): %.50f" % minAbsCum
        #                    print "Abs diff at min abs: %.50f" % absDiffAtMinAbsCum
        #                    print "Max abs diff value: %.32f" % maxAbsDiff
                    print "No of invariants considered: %d" % cntI
                    print "Max rel diff (pct) final corner: %f attained for invariant: %s" % (maxRelDiffFinal, maxFinalAtInv)                    
        #            print "Max rel diff all: %f attained for invariant: %s" % (maxRelDiffAll, maxRelDiffAllAtInv)                    
                    print "Max meanDiff-over-mean ratio is (pct): %f, attained at: %s" % (maxDiffToSignal, invAtMaxDTS)                    
                    print "Max stdDiff-over-mean ratio is (pct): %f, attained at: %s" % (maxNoiseToSignal, invAtMaxNTS)
                else:
                    print "C[%s]-dict does not have chain key: %s" % (strKey, chainId) 
            else:
                print "C-dict does not have structure key: %s" % strKey
    print "No of structures considered: %d" % cntStr
                
                
def compareResultsFinal_PvsC_Dicts(Pdict, Cdict):
    '''For comparing the values of the invariants of the included structures, C vs Python (P). 
    (ie these valeus are the results at the final corner of the simplex). Computes, among
    others, for each structure: 
    1) maximum abs-relative difference in value at the "final corner" of the simplex.
    This is measured as the abs value of the difference between the C-value and the 
    P-value divided by the abs value of the P-value. The max is taken over all invariants
    and structures included. 
        
    Input: both dicts must map structure name to invariant name to value of invariant.'''
    
    readOrder = ["w", "I12", "I1234", "I1234_full", "I1324", "I1324_full", 
    "I1324_full2", "I1423", "I1423_full0", "I1423_full2", "I123456", "I123546", 
    "I123645", "I132456", "I132546", "I132645", "I142356", "I142536", "I142635", 
    "I152346", "I152436", "I152634", "I162345", "I162435", "I162534"]    
    
    cntStr = 0
    L = 0
    relDiffs = []
    lgths = []

    for strKey in Pdict.keys():
        for chainId in Pdict[strKey].keys():
            L = Pdict[strKey][chainId]['str length']
            #print L
            L = int(L)
            if Cdict.has_key(strKey):
                if Cdict[strKey].has_key(chainId):
                    print "PDB: %s chain id: %s" % (strKey, chainId)
                    print "Length P:%d  C:%s" % (L, Cdict[strKey][chainId]['str length'])
                    cntStr += 1
                    cntI = 0
        #             maxRelDiffAll = 0
        #            maxRelDiffAllAtInv = ''
        
                    maxRelDiffFinal = 0
                    maxFinalAtInv = ''
                    absValAtMax = 0
                    
                    for Iname in Pdict[strKey][chainId].keys():
                        if readOrder.count(Iname) > 0:
                            if Cdict[strKey][chainId].has_key(Iname):
                                #print "***************Invariant: %s" % Iname
                                #comparison at final corner:
                                relDiffFinal = 0
                                val = float(Pdict[strKey][chainId][Iname])
                                absVal = abs(val)
                                diff = float(Pdict[strKey][chainId][Iname]) - float(Cdict[strKey][chainId][Iname])
                                absDiff = abs(diff)
                                relDiffFinal = 100*float(absDiff)/(absVal + 1e-10)
                                #print "Rel diff final (pct): %f" % relDiffFinal
                                if relDiffFinal > maxRelDiffFinal:
                                    maxRelDiffFinal = relDiffFinal
                                    maxFinalAtInv = Iname
                                    absValAtMax = absVal
                                cntI +=1               
                    relDiffs.append([absValAtMax, maxRelDiffFinal, L])
        #
                   # print "No of invariants considered: %d" % cntI
                    print "Max rel diff (pct) final corner: %f attained for invariant: %s; absVal of: %f" % (maxRelDiffFinal, maxFinalAtInv, absValAtMax)                    
        #            print "Max rel diff all: %f attained for invariant: %s" % (maxRelDiffAll, maxRelDiffAllAtInv)                    
                else:
                    print "C[%s]-dict does not have chain key: %s" % (strKey, chainId)
                    print "Length of chain:%d" % L
            else:
                print "C-dict does not have structure key: %s" % strKey
    relDiffs.sort()
    maxRelDiffs = [relDiffs[i][1] for i in range(len(relDiffs))]
    absValAtMaxDiffs = [relDiffs[i][0] for i in range(len(relDiffs))]
    lgths = [relDiffs[i][2] for i in range(len(relDiffs))]
    logAbsValAtMaxDiffs = [np.log10(relDiffs[i][0]) for i in range(len(relDiffs))] 
    plt.figure()
    plt.xlabel('index', fontsize = 'small')
    plt.ylabel('value', fontsize = 'small')
    plt.plot(maxRelDiffs, label = 'max rel diff (pct)')
    plt.plot(absValAtMaxDiffs, label = 'abs value')
    plt.plot(lgths, label = 'chain length')
    plt.legend(loc = 'upper left')  
    plt.title('Max relative diff, value and length')
    plt.figure()
    plt.xlabel('log10(abs value)', fontsize = 'small')
    plt.ylabel('pct', fontsize = 'small')
    plt.plot(logAbsValAtMaxDiffs, maxRelDiffs, label = 'max rel diff (pct)')
    plt.title('Max relative diffs by log-val')
    plt.legend()  

    print "No of chains considered: %d" % cntStr
    
    
#********************************************************
#********************************************************
#from here this chapter ("For comparison vs C-code") has not been corrected yet
#upon change of the loading function (except otherwise mentioned)
#********************************************************
#********************************************************
    
def I_time_length_plot_mod(PDBlist, degree = 4, d = 0):
    '''Time for computation of measures up to specified degree across all names 
    in input list; non/vectorized array/dict versions.
    Input: as output from getPDBlist.
    d: with d != 0 the dict-based version runs; with d = 0 the array-based.''' 
    timeList = []
    lengthList = []
    for l_item in PDBlist:
        f = PDBtop100Folder + l_item[1]
        L = l_item[0]
        if d ==0:
            statement = 'git.I(PDBfilename =' + 'r"' + f + '", full_b = 0, print_b=0, degree = '+ str(degree) +' ,vect_b = 1)'
        else:
            statement = 'git.I_dict(PDBfilename =' + 'r"' + f + '", full_b = 0, print_b=0, degree ='+ str(degree) +',vect_b = 1)'
        t_I = timeit.Timer(statement,'import gitP as git')
        timeList.append(t_I.timeit(number=1))
        lengthList.append(L)
    timeArray = np.array(timeList)
    return lengthList, timeList, np.sum(timeArray)
    

def plot_PvsC_time_length(PDBlist, PresultsDict, CresultsDict, timeType = 'cpuTime', 
                          titleTimePlot = "Time consumption by length", 
                          titlePtoCratio = "Ratio of time consumption by length, Python over C", 
                          labelC = 'C',
                          labelP = 'Python',
                          labelPoverCratio = 'Python over C',
                          floorCtime = 0.001, 
                          capRatio = 500,
                          markerP = 'go',
                          markerC ='bo',
                          markerRatio = 'ro',
                          as_subplots_b = 0, 
                          figSize = (16,10)):
    '''Makes time-vs-length plot of the inputs.
    timeType: tyoe of time to compare for; can be "cpuTime" or "aggrTime"; when "cpuTime" the cpu time used for 
    loading the chain, allocating memory and computing the invariants is recorded; 
    when "aggrTime" only the time for computing the invariants is considered.'''
    #generate the plot lists:
    L = 0
    timeC = 0
    timeP = 0
    lengthListC = []
    lengthListP = []

    Clist = []
    Plist = []
    ratioPoverC = 0
    ratioList = []
    
    for PDBrec in PDBlist:
        L = PDBrec[0]
        f = PDBtop100Folder + PDBrec[1]
        try:
            timeC = float(CresultsDict[f][timeType])
            lengthListC.append(int(L))
            Clist.append(timeC)
            if PresultsDict.has_key(f):
                timeP = float(PresultsDict[f][timeType])
                ratioPoverC = float(timeP)/max(timeC, floorCtime)
                ratioPoverC = min(ratioPoverC, capRatio)
                ratioList.append(ratioPoverC)
        except KeyError:
            print "This file: %s was not present in C-results" % f 
        try:
            timeP = PresultsDict[f][timeType]
            lengthListP.append(int(L))
            Plist.append(timeP)
        except KeyError:
            print "This file: %s was not present in P-results" % f 
    
    #plot time vs length    
    if as_subplots_b ==1:
        fig = plt.figure(figsize = figSize)
        fig.add_subplot('121')
    else:
        fig = plt.figure()
    plt.plot(lengthListC, Clist, markerC, label = labelC)
    plt.plot(lengthListP, Plist, markerP, label = labelP)
    plt.xlabel('length', fontsize = 'medium')
    plt.ylabel('time (s)', fontsize = 'medium')
    plt.legend(loc = "upper left", fontsize = 'medium')
    plt.title(titleTimePlot, fontsize = 'medium')
    
    #plot ratio of time consumption, "python divided by C":
    if as_subplots_b ==1:
        fig.add_subplot('122')
    else: #new fig
        plt.figure()
    plt.plot(lengthListC, ratioList, markerRatio, label = labelPoverCratio)
    plt.xlabel('length', fontsize = 'medium')
    plt.ylabel('ratio', fontsize = 'medium')
    plt.legend(loc = "upper right", fontsize = 'medium')
    plt.title(titlePtoCratio, fontsize = 'medium')
    
    
def plotCTimeByOrders(pathCresultFiles, orderList, PDBlist, fileId =  r'_excl_abs_global_mem_alloc_top100.txt', timeType = 'cpuTime', title = "Time consumption (C) by length"):
    '''Plots the time (as defined by timeType) for each order in the 
    orderList of the input results from using the C code (sitting in 
    pathCresultFiles). The C-results must be for the PDB-files in the PDBlist.'''
    plt.figure()
    for order in orderList:
        #load the C-results:
        CresFile = pathCresultFiles + r'\Ivalues_order_' + str(order) + fileId
#        print CresFile
        CresultsDict = readResultsCcode(inputFileName = CresFile, all_b = 0)
        
        Clist = []
        lengthListC = []
    
        #why not: for f in CresultsDict.keys():
        #L = CresultsDict[f]["str length"]??
    
        for PDBrec in PDBlist:
            L = PDBrec[0]
            f = PDBtop100Folder + PDBrec[1]
            try:
                timeC = float(CresultsDict[f][timeType])
                lengthListC.append(int(L))
                Clist.append(timeC)
            except KeyError:
                print "This file: %s was not present in C-results for order %d" % (f,order)
        #make the plot
        plt.plot(lengthListC, Clist, label = str(order))
    
    # the plot
    plt.xlabel('length', fontsize = 'small')
    plt.ylabel('time (s)', fontsize = 'small')
    plt.legend(loc = "upper left", fontsize = 'small')
    plt.title(title, fontsize = 'medium')


def plotCTimeByPertLength(pathCresultFiles, lengthList, PDBlist, root = '',order = 2, fileId1 =  r'_incl_abs_global_', fileId2 =  r'_100_perts_rnd_top100.txt', timeType = 'cpuTime', title = "Time consumption (C) by length, varying pert lengths"):
    '''Plots the time (as defined by timeType) for each length in the 
    lengthList of the input results from using the C code (sitting in 
    pathCresultFiles). The C-results must be for the PDB-files in the PDBlist.'''
    plt.figure()
    for lgth in lengthList:
        #load the C-results:
        CresFile = pathCresultFiles + r'\PertTime_' + str(order) + fileId1 + lgth + fileId2
#        print CresFile
        if lgth == '1':
            CresultsDict = readResultsCcode(root + r'\Ivalues_order_'+ str(order) + '_excl_abs_global_mem_alloc_top100.txt', all_b = 0, pert_b = 1, pertTime_b = 1)
        else:
            CresultsDict = readResultsCcode(inputFileName = CresFile, all_b = 0, pert_b = 1, pertTime_b = 1)
        
        Clist = []
        lengthListC = []
    
        #why not: for f in CresultsDict.keys():
        #L = CresultsDict[f]["str length"]??
    
        for PDBrec in PDBlist:
            L = PDBrec[0]
            f = PDBtop100Folder + PDBrec[1]
#            print f
            try:
                cnt = 0;
                for pertNr_key in CresultsDict[f].keys():
                    timeC = float(CresultsDict[f][pertNr_key][timeType])
                    if lgth == '1':
                        timeC = 100*timeC
                    cnt +=1 
                if cnt >1:
                    print "Warning: for the PDB-file %s there were results for more than one number of perturbations. Only the last was used."
                lengthListC.append(int(L))
                Clist.append(timeC)
            except KeyError:
                print "This PDB-file: %s was not present in C-results for length: %s" % (f,lgth)
        #make the plot
        if lgth == '1':
            label = 'base case'
        else:
            label = lgth
        plt.plot(lengthListC, Clist, label = label)
    
    # the plot
    plt.xlabel('length', fontsize = 'small')
    plt.ylabel('time (s)', fontsize = 'small')
    plt.legend(loc = "upper left", fontsize = 'small')
    plt.title(title, fontsize = 'medium')


def plotCTimeByMethod(pathCresultFilesList, PDBlist, fig, fileId1 =  r'_incl_abs_', fileId2 =  r'lgth10_100_perts_rnd_top100.txt',  methodList = [''], order = 2, nrOfPerts =100, timeType = 'cpuTime', title = "Time consumption (C) by method",subPlotNr = 111):
    '''Plots the time (as defined by timeType) for each file in the pathCresultFilesList
    storing the input results from using the C code. The methodList must corrspond to the 
    methods (or: versions) used; each method is given by a pair (name, pert_b) where name is
    the name of the method (e.g 'global_gpu') and pert_b is a boolean indicating whether the
    method was used with perturbations (1) or not (0). All C-results must be for the PDB-files 
    in the PDBlist.'''
    if not(fig): #get new figure if none was supplied
        fig = plt.figure()
    fig.add_subplot(subPlotNr)
    cntMethod = 0

    for method in methodList:
        cumTime = 0
        CresFile = pathCresultFilesList[cntMethod] # + r'\PertTime_' + str(order) + fileId1 + method + fileId2
        #load the C-results:
        CresultsDict = readResultsCcode(inputFileName = CresFile, all_b = 0, pert_b = 1, pertTime_b = 1)

#        if method[1] == 1: #if perturbation method
#            CresFile = pathCresultFilesList[cntMethod] # + r'\PertTime_' + str(order) + fileId1 + method + fileId2
#            #load the C-results:
#            CresultsDict = readResultsCcode(inputFileName = CresFile, all_b = 0, pert_b = 1, pertTime_b = 1)
#        elif method[1] == 0:                      
#            CresFile = pathCresultFilesList[cntMethod] # + r'\Ivalues_order_' + str(order) + fileId1 + method + fileId2
#            #load the C-results:
#            CresultsDict = readResultsCcode(inputFileName = CresFile, all_b = 0, pert_b = 1, pertTime_b = 1)
           
        Clist = []
        lengthListC = []
    
        #why not: for f in CresultsDict.keys():
        #L = CresultsDict[f]["str length"]??    
    
        for PDBrec in PDBlist:
            L = PDBrec[0]
            f = PDBtop100Folder + PDBrec[1]
#            print f
            cnt = 0;
            try:
                for pertNr_key in CresultsDict[f].keys():
                    timeC = float(CresultsDict[f][pertNr_key][timeType])
                    cnt +=1 
                if cnt >1:
                    print "Warning: for the PDB-file %s there were results for more than one number of perturbations. Only the last was used."
 
                lengthListC.append(int(L))
                Clist.append(timeC)
                cumTime += timeC
            except KeyError:
                print "This file: %s was not present in C-results for method: %s" % (f, method)
        #make the plot
        plt.plot(lengthListC, Clist, label = method)
        print "Total %s for method %s was: %f" % (timeType, method, cumTime) 
        cntMethod +=1
    
    # the plot
    plt.xlabel('length', fontsize = 'small')
    plt.ylabel('time (s)', fontsize = 'small')
    plt.legend(loc = "upper left", fontsize = 'small')
    plt.title(title, fontsize = 'medium')

    return fig

def plotCTimeByMethod_2(pathCresultFilesList, 
                      PDBlist,
                      fig,
#                      fileId1 =  r'_incl_abs_', 
#                      fileId2 =  r'lgth10_100_perts_rnd_top100.txt',  
                      methodList = [''], 
                      labelList = [''], 
#                      order = 2, 
                      timeType = 'cpuTime', 
                      title = "Time consumption (C) by method",
                      titlePerPertPlot = "Time per perturbation",
                      titleRatio = "Time ratio: pert. computations over base case",
                      subPlotNrRows = 1,
                      subPlotNrCols = 1,
                      perPertPlot_b = 0):
    '''Plots the time (as defined by timeType) for each file in the pathCresultFilesList
    storing the input results from using the C code. The methodList must corrspond to the 
    methods (or: versions) used; each method is given by a pair (name, pert_b) where name is
    the name of the method (e.g 'global-gpu') and pert_b is a boolean indicating whether the
    method was used with perturbations (1) or not (0). The labelList allows specifying plot-labels
    for each method. All C-results must be for the PDB-files in the PDBlist.'''
   
    if not(fig): #get new figure if none was supplied
        if perPertPlot_b ==1:
            fig, axes = plt.subplots(nrows=2, ncols=1)
        else:
            fig, axes = plt.subplots(nrows=subPlotNrRows, ncols=subPlotNrCols)
#    else:
#        axes = fig.add_subplot(subPlotNr)

    
#    ax = fig.add_subplot(subPlotNr)
#    
#    if perPertPlot_b ==1:
#        fig2 = plt.figure()
#        ax2 = fig2.add_subplot(subPlotNr)
        
    cntMethod = 0
    ratioInd = 0 #to capture if pert nr of 1 is represented in the data; in that case we plot the ratios
    baseCase = 1 
    
    for method in methodList:
        cumTime = 0
        CresFile = pathCresultFilesList[cntMethod] # + r'\PertTime_' + str(order) + fileId1 + method + fileId2
        #load the C-results:
        CresultsDict = readResultsCcode(inputFileName = CresFile, all_b = 0, pert_b = 1, pertTime_b = 1)

        fProbe = CresultsDict.keys()[0]
        if (ratioInd == 0 and CresultsDict[fProbe].has_key('1')):
            fig2 = plt.figure()
            axes2 = fig2.add_subplot(111)
            ratioInd = 1
            baseCase = 1 #this is the base case; we don't want to plot the ratio for it
            CresultsDict1 = CresultsDict
        else:
            baseCase = 0 #this is not the base case
            
#        if method[1] == 1: #if perturbation method
#            CresFile = pathCresultFilesList[cntMethod] # + r'\PertTime_' + str(order) + fileId1 + method + fileId2
#            #load the C-results:
#            CresultsDict = readResultsCcode(inputFileName = CresFile, all_b = 0, pert_b = 1, pertTime_b = 1)
#        elif method[1] == 0:                      
#            CresFile = pathCresultFilesList[cntMethod] # + r'\Ivalues_order_' + str(order) + fileId1 + method + fileId2
#            #load the C-results:
#            CresultsDict = readResultsCcode(inputFileName = CresFile, all_b = 0, pert_b = 1, pertTime_b = 1)
           
        Clist = []
        ClistPerPert = []
        lengthListC = []
        ratioList = []
    
        #why not: for f in CresultsDict.keys():
        #L = CresultsDict[f]["str length"]??    
    
        for PDBrec in PDBlist:
            L = PDBrec[0]
            f = PDBtop100Folder + PDBrec[1]
#            print f
            cnt = 0;
            try:
                for pertNr_key in CresultsDict[f].keys():
                    timeC = float(CresultsDict[f][pertNr_key][timeType])
                    pertNr = max(int(pertNr_key),1) 
                    timeCperPert = float(timeC)/pertNr
                    if (ratioInd == 1 and baseCase == 0):
                        ratio = pertNr*float(CresultsDict1[f]['1'][timeType])/(float(timeC) + 0.001)
                    cnt +=1 
                if cnt >1:
                    print "Warning: for the PDB-file %s there were results for more than one number of perturbations. Only the last was used."
 
                lengthListC.append(int(L))
                Clist.append(timeC)
                ClistPerPert.append(timeCperPert)
                if (ratioInd == 1 and baseCase ==0):
                    ratioList.append(ratio)
                cumTime += timeC
            except KeyError:
                print "This file: %s was not present in C-results for method: %s" % (f, method)
        #make the plot(s)
    
#        axes[0].set_ylabel('time (s)', fontsize = 'small')
#        axes[0].legend(loc = "upper left", fontsize = 'small')
        print "Total %s for %d nr of perts in method %s was: %f" % (timeType, pertNr, method, cumTime) 

        if perPertPlot_b == 1:
            if labelList:
                axes[0].plot(lengthListC, Clist, label = labelList[cntMethod])
                axes[1].plot(lengthListC, ClistPerPert, label = labelList[cntMethod])
            else:
                axes[0].plot(lengthListC, Clist, label = method)         
                axes[1].plot(lengthListC, ClistPerPert, label = method)      
#            axes[1].set_xlabel('length', fontsize = 'small')
#            axes[1].set_ylabel('time per perturbation (s)', fontsize = 'small')
#            axes[1].legend(loc = "upper left", fontsize = 'small')
        else:
            if labelList:
                axes.plot(lengthListC, Clist, label = labelList[cntMethod])
            else:
                axes.plot(lengthListC, Clist, label = method)     

        #plot ratios:
        if (ratioInd == 1 and baseCase ==0):
#            print ratioList
            axes2.plot(lengthListC, ratioList, label = labelList[cntMethod])


        cntMethod +=1
    
    # the plot

    if perPertPlot_b ==1:
        axes[0].set_ylabel('time (s)', fontsize = 'small')
        axes[0].legend(loc = "upper left", fontsize = 'small')
        axes[0].set_title(title, fontsize = 'medium')
        axes[1].set_title(titlePerPertPlot, fontsize = 'medium')
        axes[1].set_xlabel('length', fontsize = 'small')
        axes[1].set_ylabel('time per perturbation (s)', fontsize = 'small')
        axes[1].legend(loc = "upper left", fontsize = 'small')
    else:
        axes.set_xlabel('length', fontsize = 'small')
        axes.set_ylabel('time (s)', fontsize = 'small')
        axes.legend(loc = "upper left", fontsize = 'small')
        axes.set_title(title, fontsize = 'medium')
    
    if ratioInd == 1:
        axes2.set_xlabel('length', fontsize = 'small')
        axes2.set_ylabel('factor', fontsize = 'small')
        axes2.legend(loc = "upper left", fontsize = 'small')
        axes2.set_title(titleRatio, fontsize = 'medium')
#        fig2.show()
#    plt.xlabel('length', fontsize = 'small')
#    plt.ylabel('time (s)', fontsize = 'small')
#    plt.legend(loc = "upper left", fontsize = 'small')
#    plt.title(title, fontsize = 'medium')

    return fig

def compareDictsIII(d1,d2, tolerance = 1e-6, rel_b = 0):
    '''Compares two dictionaries of dictionaries of dictionaries.'''
    k_diff = []
    diff = {}
    #first check that the keys are id:
    for k in d1.keys():
        if not(d2.has_key(k)):
            print "Key: %s is in d1 but not in d2!" % k
    for k in d2.keys():
        if not(d1.has_key(k)):
            print "Key: %s is in d2 but not in d1!" % k
    #Then check the values:
    for k in d1.keys():
        k_d, d = compareDictsII(d1[k],d2[k], tolerance = tolerance, rel_b = rel_b)
        k_diff.append(k_d)
        diff[k] = d
    k_diff.sort
    return k_diff, diff    
    

def compareDictsIV(d1,d2, tolerance = 1e-6, rel_b = 0):
    '''Compares two dictionaries of dictionaries of dictionaries of dictionaries.'''
    k_diff = []
    diff = {}
    #first check that the keys are id:
    for k in d1.keys():
        if not(d2.has_key(k)):
            print "Key: %s is in d1 but not in d2!" % k
    for k in d2.keys():
        if not(d1.has_key(k)):
            print "Key: %s is in d2 but not in d1!" % k
    #Then check the values:
    for k in d1.keys():
        k_d, d = compareDictsIII(d1[k],d2[k], tolerance = tolerance, rel_b = rel_b)
        k_diff.append(k_d)
        diff[k] = d
    k_diff.sort
    return k_diff, diff      
    
def compareCresultsAll(Cdict1, Cdict2, tolerance = 1e-6, plot_b = 0, title = '', 
                       labelDTS = 'max MeanDiffToMeanValue', labelNTS = 'max NoiseToMeanValue',
                       markerDTS = 'bo', markerNTS = 'go'):
    '''For comparing results of C-code (from two runs using different 
    code blocks or with different parameter values). The results must be
    for all vertices across the simplex.'''
    cntStr = 0
    cntPert = 0
    absVal = 0
    relDiff = 0
    
    maxDiffs = []
    
    
    for strKey in Cdict1.keys():
        print strKey
        cntPert = 0
        for pertKey in Cdict1[strKey].keys():
            maxDiffToSignal = -1
            invAtMaxDTS = ''
            maxNoiseToSignal = -1
            invAtMaxNTS = ''
            for invKey in Cdict1[strKey][pertKey].keys():
#                if invKey == 'w':
                    maxRelDiff = -1
                    maxAbsDiff = -1
                    absValAtMaxDiff = 0
                    vals = []
                    diffs = []
                    for ijKey in Cdict1[strKey][pertKey][invKey].keys():
                        i,j = ijKey
                        if j>i:
                            val = Cdict1[strKey][pertKey][invKey][ijKey]
                            diff = Cdict1[strKey][pertKey][invKey][ijKey] -  Cdict2[strKey][pertKey][invKey][ijKey]
                            absVal = abs(val)
                            absDiff = abs(diff) 
                            #relDiff = absDiff/absVal
    #                        if relDiff > maxRelDiff:
    #                            maxRelDiff = relDiff
                            if absDiff > maxAbsDiff:
                                maxAbsDiff = absDiff
                                absValAtMaxDiff = absVal
                        #collect values for computing the std dev of the diffs:
                            diffs.append(absDiff)
                            vals.append(absVal)
                    vals = np.array(vals)
                    diffs = np.array(diffs)
                    mean = np.average(vals)
                    meanDiff = np.average(diffs)
                    diffToSignal = 100*float(meanDiff)/mean
                    if diffToSignal > maxDiffToSignal:
                        maxDiffToSignal = diffToSignal
                        invAtMaxDTS = invKey
                    stdDiff = np.std(diffs)
                    noiseToSignal = 100*float(stdDiff)/mean
                    if noiseToSignal > maxNoiseToSignal:
                        maxNoiseToSignal = noiseToSignal
                        invAtMaxNTS = invKey
    #                print "For inv: %s the mean value is: %f and the std dev in diffs is: %f" % (invKey, mean, std)
    #                print "For inv: %s the max abs diff is: %f attained at the absVal: %f" % (invKey, maxAbsDiff, absValAtMaxDiff)
            maxDiffs.append([maxDiffToSignal, maxNoiseToSignal])
            print "Pert no: %s. Max meanDiff-over-mean ratio is (pct): %f, attained at: %s" % (pertKey, maxDiffToSignal, invAtMaxDTS)                    
            print "Pert no: %s. Max stdDiff-over-mean ratio is (pct): %f, attained at: %s" % (pertKey, maxNoiseToSignal, invAtMaxNTS)
            cntPert +=1
        print "No of perturbations considered: %d" % cntPert
        cntStr += 1                 
    print "No of structures considered: %d" % cntStr
    
    if plot_b ==1:
        maxDTSs = []
        maxNTSs = []
        maxDiffs.sort()        
        maxDTSs = [maxDiffs[i][0] for i in range(len(maxDiffs))]
        maxNTSs = [maxDiffs[i][1] for i in range(len(maxDiffs))]
        plt.plot(maxDTSs, markerDTS, label = labelDTS)
        plt.plot(maxNTSs, markerNTS, label = labelNTS)
        plt.xlabel('index', fontsize = 'small')
        plt.ylabel('pct', fontsize = 'small')
        plt.legend(loc = 'upper left', fontsize = "medium")
        plt.title(title)


#Modified for new load fct (loads as multimers)!!:
def compareCresultsFinal(Cdict1, Cdict2, tolerance = 1e-6, title = 'Relative diffs in final values', 
                         labelRelDiff = 'max rel diff', labelVal = 'abs val', sortedFileList = [], use_log_b = 0):
    '''For comparing results of C-code (from two runs using different 
    code blocks or with different parameter values). The results must be
    for the "final corner" of the simplex, ie the value of the invariant
    on the structure.'''
    cntStr = 0
    cntPert = 0
    absVal = 0
    relDiff = 0

    maxRelDiffs = [] #for plotting
    absValAtMaxDiffs = [] #for plotting  
    relDiffs = []
    
    if sortedFileList:
        fileKeys = sortedFileList
    else:
        fileKeys = Cdict1.keys()
    
    for strKey in fileKeys: #Cdict1.keys():
    
        for chainKey in Cdict1[strKey].keys():
            cntPert = 0
            for pertKey in Cdict1[strKey][chainKey].keys():
                maxAbsDiff = -1
                absValAtMaxDiff = 0
                maxRelDiff = -1
                maxRelDiffAtInv = ''
                maxRelDiffAtAbsVal = -1
                for invKey in Cdict1[strKey][chainKey][pertKey].keys(): #Obs: str length will also be among the inv's here, but should give rise to no diff's!
    #                if invKey == 'I162345':
    #                    print "lengths C1: %d C2: %d" % (Cdict1[strKey]['str length'], Cdict1[strKey]['str length'])
                        absVal = abs(float(Cdict1[strKey][chainKey][pertKey][invKey]))
                        absDiff = abs(float(Cdict1[strKey][chainKey][pertKey][invKey]) -  float(Cdict2[strKey][chainKey][pertKey][invKey]))         
                        relDiff = 100*float(absDiff)/(absVal + 1e-10)
                        if (relDiff > maxRelDiff):
                            maxRelDiff = relDiff
                            maxRelDiffAtInv = invKey
                            maxRelDiffAtAbsVal = absVal
                        if absDiff > maxAbsDiff:
                            maxAbsDiff = absDiff
                            absValAtMaxDiff = absVal
                maxRelDiffs.append(maxRelDiff)
                absValAtMaxDiffs.append(maxRelDiffAtAbsVal)
                relDiffs.append([maxRelDiffAtAbsVal, maxRelDiff, strKey])
                if maxRelDiff > tolerance:
                    print strKey
                    print "Pert no: %s. Max rel diff (pct): %f attained for invariant: %s. AbsVal:%f" % (pertKey, maxRelDiff, maxRelDiffAtInv, maxRelDiffAtAbsVal)
                    print "Pert no: %s. Max absDiff: %f attained at absVal:%f" % (pertKey, maxAbsDiff, absValAtMaxDiff)                  
                else:
                    print "No differences above the tolerance: %f" % tolerance
                cntPert +=1
            print "No of perturbations considered: %d" % cntPert
            cntStr += 1
    if not(sortedFileList):
        relDiffs.sort()
    if use_log_b != 1:
        maxRelDiffs = [relDiffs[i][1] for i in range(len(relDiffs))]    
        absValAtMaxDiffs = [relDiffs[i][0] for i in range(len(relDiffs))]
        plt.xlabel('index', fontsize = 'small')
        plt.ylabel('pct', fontsize = 'small')
    else:
        maxRelDiffs = [np.log10(relDiffs[i][1]) - 2 for i in range(len(relDiffs))] #subtract 2 since relDiffs are in pct    
        absValAtMaxDiffs = [np.log10(relDiffs[i][0]) for i in range(len(relDiffs))]
        labelRelDiff = 'log10 '+ labelRelDiff 
        labelVal = 'log10 ' + labelVal 
        plt.xlabel('index', fontsize = 'small')
        plt.ylabel('log10-values', fontsize = 'small')
    plt.plot(maxRelDiffs, label = labelRelDiff)
    plt.plot(absValAtMaxDiffs, label = labelVal)
    plt.legend(loc = 'upper left')                 
    plt.title(title)

    sortedFileKeysOut = [relDiffs[i][2] for i in range(len(relDiffs))]
    
    print "No of structures considered: %d" % cntStr
    
    return sortedFileKeysOut    
    
#Modified for new load fct (loads as multimers)!!:
def compareCresultsFinal_2(Cdict1, 
                           Cdict2, 
                           titleDistr = 'Distribution of invariant-values', 
                           title = 'Abs diff and std dev of abs values', 
                           titleRatios = 'max diff over std dev'):
    '''For comparing results of C-code (from two runs using different 
    code blocks or with different parameter values). The results must be
    for the "final corner" of the simplex, ie the value of the invariant
    on the structure.
    The results in Cdict1 will be understood as "base case" values. 
    The function computes for each invariant included, the std dev
    of its value across all structures covered in the input. It then
    finds, for each invariant, the structure for which the absolute 
    difference between the results supplied in Cdict1 and Cdict2 is
    largest and plots across all invariants the ratio between this
    max absolute difference and the found std dev.'''
    
    #sorted list of the invariants' names:
    sortedInvariants = [
			#/*order 1*/
			"I12", 
			"Ia12",
			#/*order 2:*/
			"I1234", 
			"I1324", 
			"I1423",
			#/*abs value versions:*/
			#/*12*/
			"Ia1234",
			"I12a34",
			"Ia12a34",
			#/*13*/
			"Ia1324",
			"I13a24",
			"Ia13a24",
			#/*14*/
			"Ia1423",
			"I14a23",
			"Ia14a23",
			#/*order 3*/
			"I123456",
			"I123546",
			"I123645",
			#/*13*/
			"I132456",
			"I132546",
			"I132645",
			#/*14*/
			"I142356",
			"I142536",
			"I142635",
			#/*15*/
			"I152346",
			"I152436",
			"I152634",
			#/*16*/
			"I162345",
			"I162435",
			"I162534"]
    
    #find the std dev per invariant:
    fileKeys = Cdict1.keys()
    chainKeysEx = Cdict1[fileKeys[0]].keys()
    nrPerts =  Cdict1[fileKeys[0]][chainKeysEx[0]].keys() #this should be 0, i.e. the code should only be run on a single perturbation case
    if len(nrPerts) > 1:
        print "Warning: The nr of perts is larger than 1!"
    #the set of covered invariant is the same for all structures covered:
    invariants = Cdict1[fileKeys[0]][chainKeysEx[0]][nrPerts[0]]
    invariantsList = invariants.keys() #Obs: will include str length
    #print invariantsList
    #loop throught invariants and find the std devs; also plot the distr's:
    stdDevs = {}
    fig = plt.figure()
#    plt.suptitle('Distribution of invariant-values')
    plt.axis('off')
    plt.title(titleDistr)
    nInvs = len(invariantsList) - 1 #subtract 1 since str length will be among the loaded invariants
    nCols = 3
    nRows = (nInvs - nInvs%nCols)/nCols +  1
    #print nRows
    n = 1
    for inv in sortedInvariants:
        if invariantsList.count(inv) > 0:
            fig.add_subplot(nRows, nCols, n)
            vals = []
            for f in fileKeys:
                for c in Cdict1[f].keys():
                    for pNr in nrPerts:
                        vals.append(float(Cdict1[f][c][pNr][inv]))
            plt.hist(vals, bins = 100, label = inv)
            plt.xlabel('value', fontsize = 'xx-small')
            plt.ylabel('count', fontsize = 'xx-small')
            plt.legend(loc = 'upper right', fontsize = 'x-small')
            plt.xticks(fontsize = 'xx-small')
            plt.yticks(fontsize = 'xx-small')
#            plt.tight_layout()
            std = np.std(np.array(vals))
            stdDevs[inv] = std
            n +=1

    #now make the comparison between the two sets of results:
    maxDiffList = []
    for inv in sortedInvariants:
        if invariantsList.count(inv) > 0:
            vals = []
            maxDiff = 0
            for f in fileKeys:
                for c in Cdict1[f].keys():
                    for pNr in nrPerts:
                        diff = float(Cdict1[f][c][pNr][inv]) - float(Cdict2[f][c][pNr][inv])
                        absDiff = abs(diff)
                        if absDiff >= maxDiff:
                            maxDiff = absDiff
                            maxAtFile = f
            maxDiffList.append([inv, maxDiff, stdDevs[inv], maxAtFile])
    #plot:()
    diffs = [maxDiffList[i][1] for i in range(len(maxDiffList))]
    stds = [maxDiffList[i][2] for i in range(len(maxDiffList))]
    xtickLabels = [maxDiffList[i][0] for i in range(len(maxDiffList))]
    #print xtickLabels
    fig, ax = plt.subplots()
#    plt.xlabel('invariant', fontsize = 'small')
    plt.ylabel('value', fontsize = 'small')
    plt.plot(diffs, label = 'max abs diff')
    plt.plot(stds, label = 'std dev')
    ax.set_xticks(ticks = range(len(xtickLabels)))
    ax.set_xticklabels(xtickLabels, rotation = 'vertical')
    plt.legend()
    plt.title(title)
    plt.show()
    #plot ratio, diff/std-dev:
    fig, ax = plt.subplots()
    ratios = [float(maxDiffList[i][1])/maxDiffList[i][2] for i in range(len(maxDiffList))]
#    plt.xlabel('invariant', fontsize = 'small')
    plt.ylabel('value', fontsize = 'small')    
    plt.plot(ratios, label = 'ratio max abs diff over std dev')
    ax.set_xticks(ticks = range(len(xtickLabels)))
    ax.set_xticklabels(xtickLabels, rotation = 'vertical')
    plt.legend()
    plt.title(titleRatios)
    plt.show()    
    
    return stdDevs, maxDiffList
    
#Modified for new load fct (loads as multimers)!!:
def compareCresultsFinal_3(Cdict1, Cdict2, lowerPct = 0.01, upperPct = 0.99, 
                           titleDistr = 'Distribution of relative differences',
                           title = 'Statistics of comparison',
                           figSize = (14,10)):
    '''For comparing results of C-code (from two runs using different 
    code blocks or with different parameter values). The results must be
    for the "final corner" of the simplex, ie the value of the invariant
    on the structure.
    The results in Cdict1 will be understood as "base case" values. 
    The function computes for each invariant included, the mean, std dev
    and {lowerPct, upperPct} percentiles in the distributions of the relative 
    difference between the two supplied results and plots the distributions 
    per invariant. It also plots the mean, std dev and percentiles per
    invariant in a separate plot.'''
    
    #sorted list of the invariants' names:
    sortedInvariants = [
			#/*order 1*/
			"I12", 
			"Ia12",
			#/*order 2:*/
			"I1234", 
			"I1324", 
			"I1423",
			#/*abs value versions:*/
			#/*12*/
			"Ia1234",
			"I12a34",
			"Ia12a34",
			#/*13*/
			"Ia1324",
			"I13a24",
			"Ia13a24",
			#/*14*/
			"Ia1423",
			"I14a23",
			"Ia14a23",
			#/*order 3*/
			"I123456",
			"I123546",
			"I123645",
			#/*13*/
			"I132456",
			"I132546",
			"I132645",
			#/*14*/
			"I142356",
			"I142536",
			"I142635",
			#/*15*/
			"I152346",
			"I152436",
			"I152634",
			#/*16*/
			"I162345",
			"I162435",
			"I162534"]
    
    #find the std dev per invariant:
    fileKeys = Cdict1.keys()   
    chainKeysEx = Cdict1[fileKeys[0]].keys()
    nrPerts =  Cdict1[fileKeys[0]][chainKeysEx[0]].keys() #this should be 0, i.e. the code should only be run on a single perturbation case
    if len(nrPerts) > 1:
        print "Warning: The nr of perts is larger than 1!"
    #the set of covered invariant is the same for all structures covered:
    invariants = Cdict1[fileKeys[0]][chainKeysEx[0]][nrPerts[0]]
    invariantsList = invariants.keys() #Obs: will include str length
    #print invariantsList
    #loop throught invariants and find the desired statitstics; also plot the distr's:
    stdDevs = {}
    means = {}
    medians = {}
    upperPercentiles = {}
    lowerPercentiles = {}
    fig = plt.figure(figsize = figSize)
    plt.axis('off')
    plt.title(titleDistr)
#    nStrs = len(fileKeys)
    #Count the number of substructres covered:
    nStrs = 0
    for f in fileKeys:
        for c in Cdict1[f].keys():
            nStrs += 1
    upperPctIdx = int(np.ceil(nStrs*upperPct))
    lowerPctIdx = int(np.floor(nStrs*lowerPct))
    nInvs = len(invariantsList) -1 #subtract 1 since str length is among loaded invariants
    nCols = 3
    nRows = (nInvs - nInvs%nCols)/nCols +  1
    #print nRows
    #When computing relative diffs we wnat to avoid cases of low values, since these give
    #non-representative high rel diffs; for this we take the 1pct lower percentile of the
    #distr of abs-values for each invariant's:
    lowerBound = {}
    lowerBoundPctIdx = int(np.floor(nStrs*0.01))
    for inv in invariants:
        if inv != "str length":
            vals = []
            for f in fileKeys:
                for c in Cdict1[f].keys():
                    for pNr in nrPerts:
                        vals.append(abs(float(Cdict1[f][c][pNr][inv])))
            vals.sort()
            lowerBound[inv] = vals[lowerBoundPctIdx]
    #now compute the stats etc    
    n = 1
    xtickLabels = []
    meansList = []
    mediansList = []
    stdDevsList = []
    upperPercentilesList = []
    lowerPercentilesList = []
    for inv in sortedInvariants:
        if invariantsList.count(inv) > 0:
            fig.add_subplot(nRows, nCols, n)
            xtickLabels.append(inv)
            vals = []
            allVals = []
            nObs = 0
            for f in fileKeys:
                for c in Cdict1[f].keys():
                    for pNr in nrPerts:
                        try:
                            allVals.append((float(Cdict1[f][c][pNr][inv]) - float(Cdict2[f][c][pNr][inv]))/float(Cdict1[f][c][pNr][inv]))
                        except ZeroDivisionError:
                            print f #Cdict1[f]
                            print c #Cdict1[f][c]
                            print pNr #Cdict1[f][c][pNr]
                            print inv #Cdict1[f][c][pNr][inv]
                        if abs(float(Cdict1[f][c][pNr][inv])) > lowerBound[inv]:
                            vals.append((float(Cdict1[f][c][pNr][inv]) - float(Cdict2[f][c][pNr][inv]))/float(Cdict1[f][c][pNr][inv]))
                        nObs +=1
            plt.hist(vals, bins = 100, label = inv)
            plt.xlabel('value', fontsize = 'xx-small')
            plt.ylabel('count', fontsize = 'xx-small')
            plt.legend(loc = 'upper right', fontsize = 'x-small')
            plt.xticks(fontsize = 'xx-small')
            plt.yticks(fontsize = 'xx-small')
#            plt.tight_layout()
            avg = 100*np.mean(np.array(vals))
            means[inv] = avg
            meansList.append(avg)
            med = 100*np.median(np.array(allVals))
            medians[inv] = med
            mediansList.append(med)
            std = 100*np.std(np.array(vals))
            stdDevs[inv] = std
            stdDevsList.append(std)
            allVals.sort()
            upperPctIdx = int(np.ceil(nObs*upperPct))
            lowerPctIdx = int(np.floor(nObs*lowerPct))
            lP = 100*allVals[lowerPctIdx]
            lowerPercentiles[inv] = lP
            lowerPercentilesList.append(lP)
            uP = 100*allVals[upperPctIdx]
            upperPercentiles[inv] = uP
            upperPercentilesList.append(uP)
            n +=1

    #now plot it.
    fig, ax = plt.subplots(figsize = figSize)
#    plt.xlabel('invariant', fontsize = 'small')
    plt.ylabel('pct', fontsize = 'medium')
    plt.plot(mediansList, 'b-', label = 'median')
    plt.plot(meansList, 'r-', label = 'mean')
    plt.plot(lowerPercentilesList, 'c-', label = str(100*lowerPct) + ' pct')
    plt.plot(upperPercentilesList, 'c-', label = str(100*upperPct) + ' pct')
#    plt.plot(stdDevsList, label = 'std dev')
    ax.set_xticks(ticks = range(len(xtickLabels)))
    ax.set_xticklabels(xtickLabels, rotation = 'vertical')
    plt.legend(loc = 'upper left', fontsize = 'small')
    plt.title(title)
    plt.show()

#    #now make the comparison between the two sets of results:
#    maxDiffList = []
#    for inv in sortedInvariants:
#        if invariantsList.count(inv) > 0:
#            vals = []
#            maxDiff = 0
#            for f in fileKeys:
#                for pNr in nrPerts:
#                    diff = float(Cdict1[f][pNr][inv]) - float(Cdict2[f][pNr][inv])
#                    absDiff = abs(diff)
#                    if absDiff >= maxDiff:
#                        maxDiff = absDiff
#                        maxAtFile = f
#            maxDiffList.append([inv, maxDiff, stdDevs[inv], maxAtFile])
#    #plot:()
#    diffs = [maxDiffList[i][1] for i in range(len(maxDiffList))]
#    stds = [maxDiffList[i][2] for i in range(len(maxDiffList))]
#    xtickLabels = [maxDiffList[i][0] for i in range(len(maxDiffList))]
#    #print xtickLabels
#    fig, ax = plt.subplots()
##    plt.xlabel('invariant', fontsize = 'small')
#    plt.ylabel('value', fontsize = 'small')
#    plt.plot(diffs, label = 'max abs diff')
#    plt.plot(stds, label = 'std dev')
#    ax.set_xticks(ticks = range(len(xtickLabels)))
#    ax.set_xticklabels(xtickLabels, rotation = 'vertical')
#    plt.legend()
#    plt.title(title)
#    plt.show()
#    #plot ratio, diff/std-dev:
#    fig, ax = plt.subplots()
#    ratios = [float(maxDiffList[i][1])/maxDiffList[i][2] for i in range(len(maxDiffList))]
##    plt.xlabel('invariant', fontsize = 'small')
#    plt.ylabel('value', fontsize = 'small')    
#    plt.plot(ratios, label = 'ratio max abs diff over std dev')
#    ax.set_xticks(ticks = range(len(xtickLabels)))
#    ax.set_xticklabels(xtickLabels, rotation = 'vertical')
#    plt.legend()
#    plt.title(titleRatios)
#    plt.show()    
#    
    return medians, means, stdDevs, lowerPercentiles, upperPercentiles


#function for normalizing the invariants acc to P.Roegen's "Evaluating protein 
#structure descriptors and tuning Gauss integral based descriptor":
RnormWeights = {"I12":0.92, 
			"Ia12":1.37,
			#/*order 2:*/
			"I1234":1.89, 
			"I1324":0.9, 
			"I1423":1.32,
			#/*abs value versions:*/
			#/*12*/
			"Ia1234":2.26,
			"I12a34":2.26,
			"Ia12a34":2.68,
			#/*13*/
			"Ia1324":1.64,
			"I13a24":1.76,
			"Ia13a24":2.56,
			#/*14*/
			"Ia1423":2.43,
			"I14a23":1.8,
			"Ia14a23":2.76,
			#/*order 3*/
			"I123456":2.9,
			"I123546":1.88,
			"I123645":2.65,
			#/*13*/
			"I132456":1.86,
			"I132546":1.69,
			"I132645":2.09,
			#/*14*/
			"I142356":2.81,
			"I142536":1.24,
			"I142635":1.38,
			#/*15*/
			"I152346":2.19,
			"I152436":2.14,
			"I152634":1.8,
			#/*16*/
			"I162345":2.04,
			"I162435":1.19,
			"I162534":2.24}
   
def RnormalizeInv(invName, invVal, L):
    '''Returns normalized value of invariant.
    Input:
    invVal: invariant value
    invName: name of invariant
    L: length of structure
    '''
    #look up the scaling factor:
    a = RnormWeights[invName]
    #scale:
    return float(invVal)/np.power(L-1,a)

from matplotlib.colors import LogNorm
from pylab import *

#Modified for new load fct (loads as multimers)!!:
def compareNormalizedInv(Cdict1, 
                         Cdict2 = {}, 
                         titleDistr = 'Distribution of normalized invariants',
                         fontSizeDistr = "large",
                         titleDiffs = 'Differences vs values',
                         fontSizeDiffs = "large",
                         figSize = (14,10)):
    '''Plots distribution of normalized invariants.'''
    
    #sorted list of the invariants' names:
    sortedInvariants = [
			#/*order 1*/
			"I12", 
			"Ia12",
			#/*order 2:*/
			"I1234", 
			"I1324", 
			"I1423",
			#/*abs value versions:*/
			#/*12*/
			"Ia1234",
			"I12a34",
			"Ia12a34",
			#/*13*/
			"Ia1324",
			"I13a24",
			"Ia13a24",
			#/*14*/
			"Ia1423",
			"I14a23",
			"Ia14a23",
			#/*order 3*/
			"I123456",
			"I123546",
			"I123645",
			#/*13*/
			"I132456",
			"I132546",
			"I132645",
			#/*14*/
			"I142356",
			"I142536",
			"I142635",
			#/*15*/
			"I152346",
			"I152436",
			"I152634",
			#/*16*/
			"I162345",
			"I162435",
			"I162534"]
    
    #find the std dev per invariant:
    fileKeys = Cdict1.keys()
    chainKeysEx = Cdict1[fileKeys[0]].keys()
    nrPerts =  Cdict1[fileKeys[0]][chainKeysEx[0]].keys() #this should be 0, i.e. the code should only be run on a single perturbation case
    if len(nrPerts) > 1:
        print "Warning: The nr of perts is larger than 1!"
     #the set of covered invariant is the same for all structures covered:
    invariants = Cdict1[fileKeys[0]][chainKeysEx[0]][nrPerts[0]]
    invariantsList = invariants.keys() #Obs: will include str length
    #print invariantsList
    #loop through invariants and plot the distr's:


    fig = plt.figure(figsize = figSize)
#    plt.suptitle('Distribution of invariant-values')
    plt.axis('off')
    plt.suptitle(titleDistr, fontsize = fontSizeDistr)
#    nInvs = len(invariantsList) -1  #subtract 1 since str length is among loaded invariants
#    nCols = 3
#    nRows = (nInvs - nInvs%nCols)/nCols +  1
    #print nRows
        
    n = 1
    m = 1
    figNr = 1
    for inv in sortedInvariants:           
        if invariantsList.count(inv) > 0:
#            m = n%16
            if (n > 15 and figNr < 2): #new fig
                fig = plt.figure(figsize = figSize)
            #    plt.suptitle('Distribution of invariant-values')
                plt.axis('off')
                plt.suptitle(titleDistr, fontsize = fontSizeDistr)
                figNr +=1 
                m = 1
            fig.add_subplot(5,3,m) #nRows, nCols, n)
            vals = []
            for f in fileKeys:
                val = 0
                L = 0
                for c in Cdict1[f].keys():
                    for pNr in nrPerts:
                        L = Cdict1[f][c][pNr]["str length"]
                        val = float(Cdict1[f][c][pNr][inv])
                        val = RnormalizeInv(inv,val,L)
                        vals.append(val)
            plt.hist(vals, bins = 100, label = inv)
            if (m%15 == 13 or m%15 == 14 or m%15 == 0):
                plt.xlabel('value', fontsize = 'small')
            if m%3 == 1:
                plt.ylabel('count', fontsize = 'small')
            plt.legend(loc = 'upper right', fontsize = 'x-small')
            plt.xticks(fontsize = 'xx-small')
            plt.yticks(fontsize = 'xx-small')
#            plt.tight_layout()
            n +=1
            m +=1

    #if Cdict2 is provided plot distributions (based on Cdict1) and 
    #make 2d-histograms of (val, difference) :
    if Cdict2:
 
        fig2 = plt.figure(figsize = figSize)
        plt.suptitle(titleDiffs, fontsize  = fontSizeDiffs)
        plt.axis('off')
#        plt.title(titleDiffs)
#        nInvs = len(invariantsList) -1  #subtract 1 since str length is among loaded invariants
#        nCols = 3
#        nRows = (nInvs - nInvs%nCols)/nCols +  1
        #print nRows       
        n = 1
        m = 1
        figNr = 1
        for inv in sortedInvariants:
#            if inv == 'I12':
            if invariantsList.count(inv) > 0:
#                m = n%10
                if (n > 9 and n < 19 and figNr == 1): #new fig
                    fig2 = plt.figure(figsize = figSize)
                #    plt.suptitle('Distribution of invariant-values')
                    plt.axis('off')
                    plt.suptitle(titleDiffs, fontsize  = fontSizeDiffs)
                    figNr +=1
                    m = 1
                elif (n >= 19 and figNr == 2): #new fig
                    fig2 = plt.figure(figsize = figSize)
                #    plt.suptitle('Distribution of invariant-values')
                    plt.axis('off')
                    plt.suptitle(titleDiffs, fontsize  = fontSizeDiffs)
                    figNr +=1
                    m = 1
#                fig2.add_subplot(5, 3, m) #nRows, nCols, n)
                if figNr == 1:
                    fig2.add_subplot(3, 3, m) #nRows, nCols, n)
                elif figNr == 2:
                    fig2.add_subplot(3, 3, m) #nRows, nCols, n)
                elif figNr == 3:
                    fig2.add_subplot(4, 3, m) #nRows, nCols, n)
                vals = []
                diffs = []
                for f in fileKeys:
                    val1 = 0
                    val2 = 0
                    diff = 0
                    L = 0
                    for c in Cdict1[f].keys():
                        for pNr in nrPerts:
                            L = Cdict1[f][c][pNr]["str length"]
                            val1 = float(Cdict1[f][c][pNr][inv])
                            val1 = RnormalizeInv(inv,val1,L)
                            vals.append(val1)
                            val2 = float(Cdict2[f][c][pNr][inv])
                            val2 = RnormalizeInv(inv,val2,L)
                            diff = val2 - val1
                            diffs.append(diff)
                #make the plot
#                plt.tight_layout()
                plt.hist2d(np.array(vals), np.array(diffs), bins=100, norm=LogNorm())    
#                plt.scatter(np.array(vals), np.array(diffs))
                plt.title(inv, fontsize = 'medium')
                if figNr != 3:
                    if (m%9 == 7 or m%9 == 8 or m%9 == 0):  
                        plt.xlabel('value', fontsize = 'small')
                else:
                    if (m%12 == 10 or m%12 == 11 or m%12 == 0):  
                        plt.xlabel('value', fontsize = 'small')
                if m%3 == 1:
                    plt.ylabel('difference', fontsize = 'small')
                plt.xticks(fontsize = 'small')
                plt.yticks(fontsize = 'small')
                plt.colorbar()
                plt.show()
                n +=1
                m +=1
#    return vals, diffs

    
# Code for handling search for links and pokes:
def readResultsCcode_chains(inputFileName):
    '''Reads in the chains in terms of 3d-coordinates of the C-alpha trace from file 
    containing chain-data as it is output by the C-code (with write_chain_b = 1).
    
    Output: dictionary mapping each pdb-file name and chain id (covered in the
    input) to the (ordered) list of xyz-coordinates of the chain of C-alphas and 
    the corresponding list of segment coordinates (of C-alpha-to-C-alpha segments).
    More precisely: 
    pdb-file -> chain id -> [[list of chain coords],[list of [segment start, segment end] coords]]
    
    Input: file name (path to file) of file containing the results from the write out
    of the chains by means of the C-code (GIsquared).
    '''
    
    outDict = {}
    
    pdb_file = ''
        
    with open(inputFileName, 'r') as CresultsTxt:
        Cresults = csv.reader(CresultsTxt, delimiter = ';')
        for row in Cresults:
            if row[0] == "PDB_file":
                pdb_file = row[2]
                chainId = row[3]
                if not(outDict.has_key(pdb_file)):
                    outDict[pdb_file] = {}
                if not(outDict[pdb_file].has_key(chainId)):
                    outDict[pdb_file][chainId]  =[]
            else:
                idx = row[0] 
                x = float(row[1])
                y = float(row[2])
                z = float(row[3])
                outDict[pdb_file][chainId].append([x,y,z])

    #generate for each entry the corresponding segments list:
    for pdb_file in outDict.keys():
        for chainId in outDict[pdb_file].keys():
            segs = []
            chain = outDict[pdb_file][chainId]
            for i in range(len(outDict[pdb_file][chainId])-1):
                startPt = chain[i]
                endPt = chain[i+1]
                segs.append([startPt, endPt])
            outDict[pdb_file][chainId] = [chain,segs] 
                
    return outDict
    
#the following fct is extended slightly from thesis-version in that it is here possible
#to also fetch the min and max values in separate lists:
def readResultsCcode_closedLoops(inputFileName, linksAndPokes_b = 1, minMaxOut_b = 0, multiInv_b = 0, normalize_b = 0):
    '''Read in results obtained by running the C-code for getting info on 
    closed loop characteristics (search for links and "pokes").
    Output: dictionary mapping each pdb-file (name given as path to file) to
    dictionary containing the results: this dictionary maps a description ("link" 
    or "poke") to a closed loop (given as a pair of indices) to another closed 
    loop or sub chain (also given by a pair of indices) which is then finally 
    mapped to the writhe value (ie  the writhe of the pair (closed loop, closed
    loop/sub chain).
    Input:
    inputFileName: path to results file.
    linksAndPokes_b: if 1 the fct will try to find data for links and pokes 
    in the input file and return these; else the fct will look for sub-chain 
    results and return the extreme of these per PDB-file in input (1).
    minMaxOut_b: if set to 1 the results for subChainPairs (ie linksAndPokes_b != 1)
    will include separate lists of the min and the max values; else only all values
    are returned.
    multiInv_b (only for linksAndPokes_b != 1): if the input for subChainPairs cover
    a set of invariants (rather than justthe writhe), and when setting this parameter to 
    1, the outputted values can be provided in dictionaries, mapping the invariant 
    name to the output for that invariant (with output structured as for the writhe).
    normalize_b (default 0): if set to 1, the order 1 Gauss Integral values will be normlized
    by 4*pi for links and pokes (when linksAndPokes_b == 1) and for pairs of sub-chains
    (linksAndPokes_b != 1 and multiInv_b !=1).'''

#    readOrder = ["PDB_file", 
#          "chain id",
#		    "description",
#			"closed loop index 1",
#           "closed loop index 2",
#			"sub chain index 1",
#           "sub chain index 2",
#           "I12"]
          
    outDict = {}
        
    with open(inputFileName, 'r') as CresultsTxt:
        Cresults = csv.reader(CresultsTxt, delimiter = ';')
        for row in Cresults:

            idxs_loop =  (int(row[3]),int(row[4]))
#            print idxs_loop
            idxs_subChain = (int(row[5]),int(row[6]))
#            print idxs_subChain
            
            if not(outDict.has_key(row[0])): #key: filename?
                outDict[row[0]] = {}
                print row[0]
            if not(outDict[row[0]].has_key(row[1])): #key: chain id?
                outDict[row[0]][row[1]] = {}
            if not(outDict[row[0]][row[1]].has_key(row[2])): #key: description? (link/poke)
                outDict[row[0]][row[1]][row[2]] = {}
                
            if not(outDict[row[0]][row[1]][row[2]].has_key( idxs_loop)):
                outDict[row[0]][row[1]][row[2]][idxs_loop] = {}
            if not(outDict[row[0]][row[1]][row[2]][idxs_loop].has_key( idxs_subChain)):
                outDict[row[0]][row[1]][row[2]][idxs_loop][idxs_subChain] = row[7]

    if linksAndPokes_b == 1:
        #gather all results from potential links (ie closed loop pairs) and
        #gather all results from potential pokes (ie (closed loop, sub chain) pairs):
        linkValues = []
        pokeValues = []
        pokeValuesMin = []
        pokeValuesMax = []
    
        
        for file_key in outDict.keys():
            for chainId in outDict[file_key].keys():
                maxVal = -1e+20 #for pokes: we want to keep which potential poke has the max and which has the min writhe values
                minVal = 1e+20 #do
                maxAt_loop_key = [] #do
                maxAt_subChain_key = [] #do
                minAt_loop_key = [] #do
                minAt_subChain_key = [] #do
                for type_key in outDict[file_key][chainId].keys():
                    for loop_key in outDict[file_key][chainId][type_key].keys():
                        for subChain_key in outDict[file_key][chainId][type_key][loop_key].keys():
                            val = float(outDict[file_key][chainId][type_key][loop_key][subChain_key])
                            if normalize_b ==1:
                                val = float(val/(4*np.pi))
                            if (minMaxOut_b == 1 and type_key == 'poke'):
                                #keep the val etc if max or min:
                                if val > maxVal:
                                    maxVal = val
                                    maxAt_loop_key = loop_key
                                    maxAt_subChain_key = subChain_key
                                if val < minVal:
                                    minVal = val
                                    minAt_loop_key = loop_key
                                    minAt_subChain_key = subChain_key
                            if type_key == 'link':
                                linkValues.append([val, file_key, chainId, loop_key, subChain_key] )
                            elif type_key == 'poke':
                                pokeValues.append([val, file_key, chainId, loop_key, subChain_key])    
                    if (minMaxOut_b == 1 and type_key == 'poke'):
                        pokeValuesMax.append([maxVal, file_key, chainId, maxAt_loop_key, maxAt_subChain_key])    
                        pokeValuesMin.append([minVal, file_key, chainId, minAt_loop_key, minAt_subChain_key])         
        if minMaxOut_b == 1:
            return outDict, linkValues, pokeValues, pokeValuesMin, pokeValuesMax
        else:
            return outDict, linkValues

    else: #if linksAndPokes_b != 1: fetch the subChain results
    
        if multiInv_b !=1:
            extremeWritheValues = []
            extremeWritheValuesMin = []
            extremeWritheValuesMax = []
    
            for file_key in outDict.keys():
                for chainId in outDict[file_key].keys():
                    maxVal = -1e+20 #for sub-chain pairs: we want to fetch only the max and min writhe values
                    minVal = 1e+20 #do
                    maxAt_loop_key = [] #do
                    maxAt_subChain_key = [] #do
                    minAt_loop_key = [] #do
                    minAt_subChain_key = [] #do
                    for type_key in outDict[file_key][chainId].keys():
                        for loop_key in outDict[file_key][chainId][type_key].keys():
                            for subChain_key in outDict[file_key][chainId][type_key][loop_key].keys():
                                val = float(outDict[file_key][chainId][type_key][loop_key][subChain_key])
                                if normalize_b ==1:
                                    val = float(val/(4*np.pi))
                                if type_key == 'subChainPair':
                                    #keep the val etc if max or min:
                                    if val > maxVal:
                                        maxVal = val
                                        maxAt_loop_key = loop_key
                                        maxAt_subChain_key = subChain_key
                                    if val < minVal:
                                        minVal = val
                                        minAt_loop_key = loop_key
                                        minAt_subChain_key = subChain_key
                        if maxVal < -10:
                            print "maxVal < -10 for: %s" % file_key
                        if type_key == 'subChainPair':
                            extremeWritheValues.append([maxVal, file_key, chainId, maxAt_loop_key, maxAt_subChain_key])    
                            extremeWritheValues.append([minVal, file_key, chainId, minAt_loop_key, minAt_subChain_key])  
                            extremeWritheValuesMax.append([maxVal, file_key, chainId, maxAt_loop_key, maxAt_subChain_key])    
                            extremeWritheValuesMin.append([minVal, file_key, chainId, minAt_loop_key, minAt_subChain_key])  
            if minMaxOut_b == 1:
                return outDict, extremeWritheValues, extremeWritheValuesMin, extremeWritheValuesMax
            else:
                return outDict, extremeWritheValues
        
        if multiInv_b ==1:
            
            extremeWritheValues = {}
            extremeWritheValuesMin = {}
            extremeWritheValuesMax = {}
    
            for file_key in outDict.keys():
                for chainId in outDict[file_key].keys():
                    for type_key in outDict[file_key][chainId].keys():
                        maxVal = -1e+20 #for sub-chain pairs: we want to fetch only the max and min writhe values
                        minVal = 1e+20 #do
                        maxAt_loop_key = [] #do
                        maxAt_subChain_key = [] #do
                        minAt_loop_key = [] #do
                        minAt_subChain_key = [] #do
                        for loop_key in outDict[file_key][chainId][type_key].keys():
                            for subChain_key in outDict[file_key][chainId][type_key][loop_key].keys():
                                val = float(outDict[file_key][chainId][type_key][loop_key][subChain_key])
                                #keep the val etc if max or min:
                                if val > maxVal:
                                    maxVal = val
                                    maxAt_loop_key = loop_key
                                    maxAt_subChain_key = subChain_key
                                if val < minVal:
                                    minVal = val
                                    minAt_loop_key = loop_key
                                    minAt_subChain_key = subChain_key 
                        if not(extremeWritheValues.has_key(type_key)):
                            extremeWritheValues[type_key] = []
                            extremeWritheValuesMin[type_key] = []
                            extremeWritheValuesMax[type_key] = []
                        extremeWritheValues[type_key].append([maxVal, file_key, chainId, maxAt_loop_key, maxAt_subChain_key])    
                        extremeWritheValues[type_key].append([minVal, file_key, chainId, minAt_loop_key, minAt_subChain_key])  
                        extremeWritheValuesMax[type_key].append([maxVal, file_key, chainId, maxAt_loop_key, maxAt_subChain_key])    
                        extremeWritheValuesMin[type_key].append([minVal, file_key, chainId, minAt_loop_key, minAt_subChain_key])  
            if minMaxOut_b == 1:
                return outDict, extremeWritheValues, extremeWritheValuesMin, extremeWritheValuesMax
            else:
                return outDict, extremeWritheValues        
        

def histoClosedLoops(linkValues, 
                     pokeValues, 
                     binsLinks = 10, 
                     binsPokes = 10, 
                     pokesDistrLog_b = 0,
                     linksShowTopNr = 100, 
                     pokesShowTopNr =100, 
                     cutOffLinks = 11.0, 
                     cutOffPokes = 13.0, 
                     lowerPercentile = 0.01,
                     upperPercentile = 0.99,
                     chainsFile = '', 
                     inFig = '',
                     titleDistr = 'Restricted search',
                     linksPokesDistrInOnePlot_b = 0,
                     linksDistrPlotsInSubPlots_b = 0, 
                     pokesDistrPlotsInSubPlots_b = 0, 
                     plotLinkExs_b = 0, 
                     plotPokeExs_b = 0 , 
                     plotLinksTop10_b = 0, 
                     plotPokesTop10_b = 0,
                     titleLinks = 'Distribution of writhe values for potential links',
                     titlePokesI = 'Distribution of writhe values for potential pokes',
                     titlePokesII = ' writhe values for potential pokes',
                     color = 'blue',
                     intervalLinks = [2,3],
                     intervalPokes = [2,3],
                     plotExsInInt_b = 0, 
                     elevAzimList =  [[15,-90], [15,0], [90,0]], 
                     color1 = 'red', 
                     markerAtom1 = 'o', 
                     markerLine1 = '-', 
                     color2 = 'blue', 
                     markerAtom2 = 'o', 
                     markerLine2 = '-', 
                     stringPattern = '([\S]+)([^\w]+)([\w]+)\.([\w]+)'):
    '''Plot histograms of the distributions of the writhe numbers from potential
    links and "pokes". Makes 3d plots of cases for which the writhe is outside the 
    interval [-cutOff, cutOff], i.e. of the "extreme writhe" cases; the relevant closed 
    loops/subChains are highlighted in each case. Examples with writhe values inside a
    set interval (intervalLinks and intervalPokes) can also be chosen (use plotExsInInt_b
    =1).'''

    #plot distributions of all values:
    #links
    linkValuesAll = [linkValues[i][0] for i in range(len(linkValues)) ]

    #We want first to find the lower and upper percentiles as desired (approx'ly):
    L = len(linkValuesAll)
    print "Number of potential links: %d" % L
    #sort the values:
    linkValuesAll.sort(reverse = False)
    #find indexes "surrounding" the percentiles:
    il_lowerPct = int(np.floor(L*lowerPercentile))
    ir_lowerPct = il_lowerPct +1 # equals int(np.ceil(L*0.01)) unless L*0.01 is integral
    print "(Links) idx's left and right of lower percentile %f: %d, %d with values: %f, %f" % (lowerPercentile, il_lowerPct, ir_lowerPct, linkValuesAll[il_lowerPct], linkValuesAll[ir_lowerPct] )
    il_upperPct = int(np.floor(L*upperPercentile))
    ir_upperPct = il_upperPct +1 # equals int(np.ceil(L*0.99)) unless L*0.99 is integral
    print "(Links) idx's left and right of upper percentile %f: %d, %d with values: %f, %f" % (upperPercentile, il_upperPct, ir_upperPct, linkValuesAll[il_upperPct], linkValuesAll[ir_upperPct] )
    #let the desired percentiles be set as the weighted average over the values left and right of the percentile:
    lowerPct = (L*lowerPercentile- il_lowerPct) *linkValuesAll[ir_lowerPct] + (ir_lowerPct - L*lowerPercentile)*linkValuesAll[il_lowerPct]
    upperPct = (L*upperPercentile- il_upperPct)*linkValuesAll[ir_upperPct] + (ir_upperPct - L*upperPercentile)*linkValuesAll[il_upperPct]
    print "(Links) Lower percentile %f pct at: %f; upper percentile %f pct at: %f" % (100*lowerPercentile, lowerPct, 100*upperPercentile, upperPct)

    
    kur = stats.kurtosis(linkValuesAll)
    print "(Links) Kurtosis is: %f" % kur
    
    #a straight histo plot: plt.hist(linkValuesAll, bins = binsLinks)
    #But: we want to transform the counts to log10(1+count)
    #Can be done like this:
#    counts, bins, patches = plt.hist(linkValuesAll, bins = binsLinks)
#    plt.close()
    counts, bins = np.histogram(linkValuesAll, bins = binsLinks) 
#    print bins
    print "len bins:%d" % len(bins)
    #bin mid points
    midPts = [float(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]
    #width of bins:
    w = [float(bins[i+1] - bins[i]) for i in range(len(bins)-1)]
    #then make a bar plot:
    if linksPokesDistrInOnePlot_b == 1:
        fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1, ncols = 2)
        ax1_1 = ax1[0]
        titleDistr = 'Links'
        titleLinks = 'Distribution of writhe values for potential links and pokes'
    elif linksDistrPlotsInSubPlots_b == 1:
        fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1, ncols = 3)
        ax1_1 = ax1[1]
        titleDistr = 'All'
    elif linksDistrPlotsInSubPlots_b == 2:
        if not(inFig): 
            fig1 = plt.figure(figsize = (14,8))
            ax1_1 = fig1.add_subplot(121)
            titleDistr = titleDistr
        else:
            fig1 = plt.figure(num = 1) #fetches the first of the list of existing figures
            ax1_1 = fig1.add_subplot(122)
            titleDistr = titleDistr
#        fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1, ncols = 2)
#        ax1_1 = ax1[0]
    else:
        fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1, ncols = 1)
        ax1_1 = ax1        
        titleDistr = 'All'
        
    ax1_1.bar([bins[i] for i in range(len(bins)-1)], np.log10(1+counts), width = w, color = color, fill = True, alpha = 0.8)
    ax1_1.set_xlabel('writhe', fontsize = "medium")
    ax1_1.set_ylabel('log(1+count)', fontsize = "medium")
    ax1_1.tick_params(labelsize = "medium")
    
    titleLinks = titleLinks
    fig1.suptitle(titleLinks, fontsize = "large")
    ax1_1.set_title(titleDistr, fontsize = "large")


    linkValues.sort(reverse = True) #21-Aug-2016
    pokeValues.sort(reverse = True) #21-Aug-2016


    #plot top link values (unless links and pokes distr are made in one plot)
    if linksPokesDistrInOnePlot_b != 1:
        
        #linkValues.sort(reverse = True) 21-Aug-2016
        
        if linksDistrPlotsInSubPlots_b ==1:
            ax1_2 = ax1[2] 
        elif linksDistrPlotsInSubPlots_b != 2:
            fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1, ncols = 1)
            ax1_2 = ax1
         
        if linksDistrPlotsInSubPlots_b != 2:
            selectedLinkValuesTop = [linkValues[i][0] for i in range(linksShowTopNr) ]
            ax1_2.hist(selectedLinkValuesTop, bins = binsLinks, color = color, fill = True, alpha = 0.8)
            titleLinks = 'Upper ' + str(linksShowTopNr) #+ ' highest writhe values for potential links'
            ax1_2.set_title(titleLinks, fontsize = "large")
            ax1_2.set_xlabel('writhe', fontsize = "medium")
            ax1_2.set_ylabel('count', fontsize = "medium")
            ax1_2.tick_params(labelsize = "medium")
    
          
        #plot bottom link values
        if linksDistrPlotsInSubPlots_b ==1:
            ax1_3 = ax1[0] 
        elif linksDistrPlotsInSubPlots_b != 2:
            fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1, ncols = 1)
            ax1_3 = ax1 
    
        if linksDistrPlotsInSubPlots_b != 2:
            selectedLinkValuesBottom = [linkValues[::-1][i][0] for i in range(linksShowTopNr) ] #linkValues were sorted in reverse above
            ax1_3.hist(selectedLinkValuesBottom, bins = binsLinks, color = color, fill = True)
            titleLinks = 'Lower ' + str(linksShowTopNr) #+ ' lowest writhe values for potential links'
            ax1_3.set_title(titleLinks, fontsize = "large")
            ax1_3.set_xlabel('writhe', fontsize = "medium")
            ax1_3.set_ylabel('count', fontsize = "medium")
            ax1_3.tick_params(labelsize = "medium")

        
#    #pokes
    pokeValuesAll = [pokeValues[i][0] for i in range(len(pokeValues)) ]
    print "Min writhe for pokes:%d" % min(pokeValuesAll)
    print "Max writhe for pokes:%d" % max(pokeValuesAll)    


    #We want first to find the lower and upper percentiles as desired (approx'ly):
    L = len(pokeValuesAll)
    print "Number of potential pokes: %d" % L
    #sort the values:
    pokeValuesAll.sort(reverse = False)
    #find indexes "surrounding" the percentiles:
    il_lowerPct = int(np.floor(L*lowerPercentile))
    ir_lowerPct = il_lowerPct +1 # equals int(np.ceil(L*0.01)) unless L*0.01 is integral
    print "(Pokes) idx's left and right of lower percentile %f: %d, %d with values: %f, %f" % (lowerPercentile, il_lowerPct, ir_lowerPct, pokeValuesAll[il_lowerPct], pokeValuesAll[ir_lowerPct] )
    il_upperPct = int(np.floor(L*upperPercentile))
    ir_upperPct = il_upperPct +1 # equals int(np.ceil(L*0.99)) unless L*0.99 is integral
    print "(Pokes) idx's left and right of upper percentile %f: %d, %d with values: %f, %f" % (upperPercentile, il_upperPct, ir_upperPct, pokeValuesAll[il_upperPct], pokeValuesAll[ir_upperPct] )
    #let the desired percentiles be set as the weighted average over the values left and right of the percentile:
    lowerPct = (L*lowerPercentile- il_lowerPct) *pokeValuesAll[ir_lowerPct] + (ir_lowerPct - L*lowerPercentile)*pokeValuesAll[il_lowerPct]
    upperPct = (L*upperPercentile- il_upperPct)*pokeValuesAll[ir_upperPct] + (ir_upperPct - L*upperPercentile)*pokeValuesAll[il_upperPct]
    print "(Pokes) Lower percentile %f pct at: %f; upper percentile %f pct at: %f" % (100*lowerPercentile, lowerPct, 100*upperPercentile, upperPct)

    
 

    #straigth histo plot: plt.hist(pokeValuesAll, bins = binsPokes)
    #But: we want to transform the counts to log10(1+count)
    #Can be done like this:
#    counts, bins, patches = plt.hist(pokeValuesAll, bins = binsPokes)
#    plt.close()
#    #bin mid points
#    midPts = [float(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]
#    #width of bins:
#    w = [float(bins[i+1] - bins[i]) for i in range(len(bins)-1)]
    #then make a bar plot:
    
    if linksPokesDistrInOnePlot_b == 1:
        
        ax1_2 = ax1[1]
        
        if pokesDistrLog_b == 1:
            #as above for the links, we want to transform the counts to log10(1+count):
            counts, bins = np.histogram(pokeValuesAll, bins = binsLinks)    
            #bin mid points
            midPts = [float(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]
            #width of bins:
            w = [float(bins[i+1] - bins[i]) for i in range(len(bins)-1)]
            #then make a bar plot:        
            ax1_2.bar(midPts, np.log10(1+counts), width = w, color = color)
            ax1_2.set_ylabel('log(1+count)', fontsize = "medium")

        else: #just make std histogram
            ax1_2.hist(pokeValuesAll, bins = binsPokes, color = color)
            ax1_2.set_ylabel('count', fontsize = "medium")

        ax1_2.set_xlabel('writhe_abs', fontsize = "medium")
        ax1_2.tick_params(labelsize = "medium")
        titleDistr = 'Pokes'
        ax1_2.set_title(titleDistr, fontsize = "large")
        
        
    else: #links and pokes distr not shown in one plot 
        if pokesDistrPlotsInSubPlots_b ==1:
            fig2, ax2 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 3)
            ax2_1 = ax2[1]
        else:
            fig2, ax2 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 1)
            ax2_1 = ax2
    #    plt.bar(midPts, np.log10(1+counts), width = w)
        ax2_1.hist(pokeValuesAll, bins = binsPokes, color = color)
        ax2_1.set_title('All', fontsize = "large")
    #    plt.ylabel('log(1+count)')
        ax2_1.set_xlabel('writhe', fontsize = "medium")
        ax2_1.set_ylabel('count', fontsize = "medium")
        ax2_1.tick_params(labelsize = "medium")
        titlePokes = titlePokesI
        fig2.suptitle(titlePokes, fontsize = "large")
    
    
        #plot top poke values   
        if pokesDistrPlotsInSubPlots_b ==1:
            ax2_2 = ax2[2]
        else:
            fig2, ax2 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 1)
            ax2_2 = ax2
            
        selectedPokeValues = [pokeValues[i][0] for i in range(pokesShowTopNr) ]
        ax2_2.hist(selectedPokeValues, bins = binsPokes, color = color)
        titlePokes= 'Upper ' + str(pokesShowTopNr) #+ titlePokesII
        ax2_2.set_title(titlePokes, fontsize = "large")
        ax2_2.set_xlabel('writhe', fontsize = "medium")
        ax2_2.set_ylabel('count', fontsize = "medium")
        ax2_2.tick_params(labelsize = "medium")
        
        #plot bottom poke values   
        if pokesDistrPlotsInSubPlots_b ==1:
            ax2_3 = ax2[0]
        else:
            fig2, ax2 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 1)
            ax2_3 = ax2
            
        selectedPokeValues = [pokeValues[::-1][i][0] for i in range(pokesShowTopNr) ]
        ax2_3.hist(selectedPokeValues, bins = binsPokes, color = color)
        titlePokes= 'Lower ' + str(pokesShowTopNr) #+ titlePokesII
        ax2_3.set_title(titlePokes, fontsize = "large")
        ax2_3.set_xlabel('writhe', fontsize = "medium")
        ax2_3.set_ylabel('count', fontsize = "medium")
        ax2_3.tick_params(labelsize = "medium")

    
    #print out the top 10 potential links and pokes:
    #load segment chains if 3d-plots desired:
    if (plotLinksTop10_b == 1 or plotPokesTop10_b ==1 or plotLinkExs_b == 1 or  plotPokeExs_b == 1 or plotExsInInt_b == 1):
        chains = readResultsCcode_chains(inputFileName = chainsFile)    
        
    print "Top 5 link values:"
    print linkValues[:5]  #Obs: sorted above (btw: moved the sorting on 21-Aug2016)
    print "Top 5 lowest link values:"
    print linkValues[::-1][:5]    
    if plotLinksTop10_b == 1:
        plotClosedLoopExamples(examples = linkValues[:5], chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern) 
        plotClosedLoopExamples(examples = linkValues[::-1][:5], chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern) 
    print "Top 5 poke values:"
    print pokeValues[:5] #Obs: sorted above (btw: moved the sorting on 21-Aug2016)
    print "Top 5 lowest poke values:"
    print pokeValues[::-1][:5] #Obs: sorted above (btw: moved the sorting on 21-Aug2016)
    if plotPokesTop10_b ==1: 
        plotClosedLoopExamples(examples = pokeValues[:5], chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern) 
        plotClosedLoopExamples(examples = pokeValues[::-1][:10], chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern) 
    #number of links below/above -cutoff/cutoff and plots of these:
    linkValuesExtremePos = []
    linkValuesExtremeNeg = []
    linkValuesInterval = []
    pokeValuesExtremeNeg = []
    pokeValuesExtremePos = []
    pokeValuesInterval = []

    cntPlus = 0
    cntNeg = 0
    cntInInt = 0
    for ex in linkValues:
        val = ex[0]
        if val > cutOffLinks:
            linkValuesExtremePos.append(ex)
            cntPlus +=1
        if val < -cutOffLinks:
            linkValuesExtremeNeg.append(ex)
            cntNeg +=1
        if (val > intervalLinks[0] and val < intervalLinks[1]):
            linkValuesInterval.append(ex)
            cntInInt +=1
    print "There were %d links having a writhe value below the cutoff: %f" % (cntNeg, -cutOffLinks)
    print "There were %d links having a writhe value above the cutoff: %f" % (cntPlus, cutOffLinks)
    print "There were %d links having a writhe value in the interval: %s" % (cntInInt, str(intervalLinks))


    if plotLinkExs_b == 1:
        plotClosedLoopExamples(examples = linkValuesExtremePos, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
        plotClosedLoopExamples(examples = linkValuesExtremeNeg, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
    if plotExsInInt_b == 1:
        plotClosedLoopExamples(examples = linkValuesInterval, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
        

    cntExPokesPos = 0
    cntExPokesNeg = 0
    cntInInt = 0 #reset
    for ex in pokeValues:
        val = ex[0]
        if val > cutOffPokes:
            pokeValuesExtremePos.append(ex)
            cntExPokesPos +=1
        if val < -cutOffPokes:
            pokeValuesExtremeNeg.append(ex)
            cntExPokesNeg +=1
        if (val > intervalPokes[0] and val < intervalPokes[1]):
            pokeValuesInterval.append(ex)
            cntInInt +=1
    print "There were %d pokes having a writhe value below the cutoff: %f" % (cntExPokesNeg, -cutOffPokes)
    print "There were %d pokes having a writhe value above the cutoff: %f" % (cntExPokesPos, cutOffPokes)
    print "There were %d pokes having a writhe value in the interval: %s" % (cntInInt, str(intervalPokes))
    if plotPokeExs_b == 1:
        plotClosedLoopExamples(examples = pokeValuesExtremeNeg, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
        plotClosedLoopExamples(examples = pokeValuesExtremePos, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
    if plotExsInInt_b == 1:
        plotClosedLoopExamples(examples = pokeValuesInterval, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    

    if linksDistrPlotsInSubPlots_b ==2:
        return fig1
          
#Version of histoClosedLoops allowing to plo min/max separately for pokes
#To use this fct get use readResultsCcode_closedLoops as follows:
#linkValues, pokeValuesAll, pokeValuesMin, pokeValuesMax  = readResultsCcode_closedLoops(inputFileName, linksAndPokes_b = 1, minMaxOut_b = 0 ...[other params])
#Then set pokeValues = pokeValuesAll, pokeValuesMin, pokeValuesMax
#and fire off histoClosedLoops_minMax(linkValues, pokeValues, ...)
def histoClosedLoops_minMax(linkValues, 
                     pokeValues, 
                     pokesMinMax_b = 1, 
                     binsLinks = 10, 
                     binsPokes = 10, 
                     pokesDistrLog_b = 0,
                     linksShowTopNr = 100, 
                     pokesShowTopNr =100, 
                     cutOffLinks = 11.0, 
                     cutOffPokes = 13.0, 
                     lowerPercentile = 0.01,
                     upperPercentile = 0.99,
                     chainsFile = '', 
                     inFig = '',
                     titleDistr = 'Restricted search',
                     linksPokesDistrInOnePlot_b = 0,
                     linksDistrPlotsInSubPlots_b = 0, 
                     pokesDistrPlotsInSubPlots_b = 0, 
                     plotLinkExs_b = 0, 
                     plotPokeExs_b = 0 , 
                     plotLinksTop10_b = 0, 
                     plotPokesTop10_b = 0,
                     titleLinks = 'Distribution of writhe values for potential links',
                     titlePokesI = 'Distribution of writhe values for potential pokes',
                     titlePokesII = ' writhe values for potential pokes',
                     color = 'blue',
                     colorMin = 'lightgrey',
                     colorMax = 'grey',
                     alphaMin = 1.0,
                     alphaMax = 0.8,
                     supTitle_b = 1,
                     fontsizeLabel = "large",
                     fontsizeSubTitle = "large",
                     fontsizeTitle = "large",
                     intervalLinks = [2,3],
                     intervalPokes = [2,3],
                     plotExsInInt_b = 0, 
                     elevAzimList =  [[15,-90], [15,0], [90,0]], 
                     color1 = 'red', 
                     markerAtom1 = 'o', 
                     markerLine1 = '-', 
                     color2 = 'blue', 
                     markerAtom2 = 'o', 
                     markerLine2 = '-', 
                     stringPattern = '([\S]+)([^\w]+)([\w]+)\.([\w]+)'):
    '''Plots histograms of the distributions of the writhe numbers from potential
    links and "pokes". 
    Makes 3d plots of cases for which the writhe is outside the interval [-cutOff, cutOff], i.e. 
    of the "exectional writhe" cases; the relevant closed loops/subChains are highlighted in each case. 
    Examples having a writhe value in top10 of all (in absolute value) of examples having a writhe inside 
    a set interval (intervalLinks and intervalPokes) can also be chosen (use plotLinksTop10_b = plotPokesTop10_b 
    = 1 and plotExsInInt_b = 1, respectively).

    Input:
     linkValues and pokeValues: data for the potential links/pokes as output by readResultsCcode_closedLoops; the structure is
     
     linkValues, pokeValuesAll, pokeValuesMin, pokeValuesMax  = readResultsCcode_closedLoops(...) 

     Then let pokeValues pokeValuesAll, pokeValuesMin, pokeValuesMax (see also usage section for explicit examples).  
     
     pokesMinMax_b (1): plots the distributions of the max-values and of the min-values in one (diff colors, colorMin and colorMax, can be chosen) 
     binsLinks (10): number of bins in links-distribution plot 
     binsPokes (10): number of bins in pokes-distribution plot   
     pokesDistrLog_b (0): if set to 1 the frequency (y-axis) in the distr's will be done by log(count +1)
     linksShowTopNr (100): plot upper and lower tails (count 100) in links' distr's  
     pokesShowTopNr (100): plot upper and lower tails (count 100) in pokes' distr's  
     cutOffLinks (0.9): links with an abs value of writhe above this scalar can be shown as 3d-plots (set plotLinksEx_b = 1) 
     cutOffPokes (0.8): pokes with an abs value of writhe above this scalar can be shown as 3d-plots (set plotPokesEx_b = 1)  
     lowerPercentile (0.01): counts for this percentile will be printed on the screen (lower: self explanatory)
     upperPercentile (0.99): as lowerPercentile 
     chainsFile (''): name (path to) of file holding the data of the chains needed for 3d-plotting out the examples. The structure
     of this file must be as output of C-code; the file will be read in by readResultsCcode_chains so see there for more 
     inFig (''): a list of (open) figures; if provided the first in the list is fetched and the plot of the distribution
     of writhes for potential links is added as a subplot next to that figure. So can be used to lot two distr's side-by-side. 
     titleDistr ('Restricted search'): "Super" title of distribution plot 
     supTitle_b (1): if set to 0 no super title will be printed (see titleDistr)
     linksPokesDistrInOnePlot_b (0): if set to 1 links and pokes distrs's will be done in one plot
     linksDistrPlotsInSubPlots_b (0): if set to 1 links distrs's will be done in one plot (tails plots added on the sides)
     pokesDistrPlotsInSubPlots_b (0): if set to 1 pokes distrs's will be done in one plot (tails plots added on the sides)
     plotLinkExs_b (0): see cutOffLinks
     plotPokeExs_b (0): see cutOffPokes/cutOffLinks 
     plotLinksTop10_b (0): if set to 1, the links with in top10 wrt the abs value of writhe will be shown as 3d-plots (the 5 top positives and the 5 top negatives)
     plotPokesTop10_b (0): if set to 1, the pokes with in top10 wrt the abs value of writhe will be shown as 3d-plots (the 5 top positives and the 5 top negatives)
     intervalLinks ([2,3]): see plotExsInInt_b; the number of potential links having a writhe in this interval is shown 
     intervalPokes ([2,3]): seeplotExsInInt_b; the number of potential pokes having a writhe in this interval is shown 
     plotExsInInt_b (0): if set to 1 the potential links and pokes having writhe values in intervalLinks and in intervalPokes will be shown (3d-plots) 
     titleLinks ('Distribution of writhe values for potential links'): title on links distr plot
     titlePokesI ('Distribution of writhe values for potential pokes'): title on pokes distr plot
     titlePokesII (' writhe values for potential pokes'): last part of title on zoom-in on tails of distr of pokes
     color ('blue'): color of bars in distr plots (if colorMin, colorMax are not used)
     colorMin ('lightgrey'): see pokesMinMax_b
     colorMax ('grey'): see pokesMinMax_b
     alphaMin (1.0): transparency factor in min/max distr plots (min-part)
     alphaMax (0.8): transparency factor in min/max distr plots (max-part)
     fontsizeLabel ("large"): font size of plot labels
     fontsizeSubTitle ("large"): font size of plot sub-titles
     fontsizeTitle ("large"): font size of plot super title (see titleDistr)
     elevAzimList ([[15,-90], [15,0], [90,0]]): list of elevation and azimutal angles used in 3d-plots. The length of the list determines the number of 
     subplots in the 3d-plots. 
     color1 ('red'): color of first high-lighted chain-segment in 3d-plot (each plotted example consists in two segments, a closed loop and a closed loop/short sub-chain)  
     markerAtom1 ('o'): marker of C-alpha atoms in first high-lighted chain-segment in 3d-plot (see also color1) 
     markerLine1 ('-'): marker of C-alpha atoms not in first or second high-lighted chain-segment in 3d-plot (see also color1)
     color2 ('blue'): color of second high-lighted chain-segment in 3d-plot (see also color1) 
     markerAtom2 ('o'): marker of C-alpha atoms in second high-lighted chain-segment in 3d-plot (see also color1)  
     markerLine2 ('-'): marker of C-alpha atoms not in first or second high-lighted chain-segment in 3d-plot (see also color1)
     stringPattern ('([\S]+)([^\w]+)([\w]+)\.([\w]+)'): string pattern o path to pdb-file; used for extracting the pdb-id for use in
     3d-plots'''

    #plot distributions of all values:
    #links
    linkValuesAll = [linkValues[i][0] for i in range(len(linkValues)) ]

    #We want first to find the lower and upper percentiles as desired (approx'ly):
    L = len(linkValuesAll)
    print "Number of potential links: %d" % L
    #sort the values:
    linkValuesAll.sort(reverse = False)
    #find indexes "surrounding" the percentiles:
    il_lowerPct = int(np.floor(L*lowerPercentile))
    ir_lowerPct = il_lowerPct +1 # equals int(np.ceil(L*0.01)) unless L*0.01 is integral
    print "(Links) idx's left and right of lower percentile %f: %d, %d with values: %f, %f" % (lowerPercentile, il_lowerPct, ir_lowerPct, linkValuesAll[il_lowerPct], linkValuesAll[ir_lowerPct] )
    il_upperPct = int(np.floor(L*upperPercentile))
    ir_upperPct = il_upperPct +1 # equals int(np.ceil(L*0.99)) unless L*0.99 is integral
    print "(Links) idx's left and right of upper percentile %f: %d, %d with values: %f, %f" % (upperPercentile, il_upperPct, ir_upperPct, linkValuesAll[il_upperPct], linkValuesAll[ir_upperPct] )
    #let the desired percentiles be set as the weighted average over the values left and right of the percentile:
    lowerPct = (L*lowerPercentile- il_lowerPct)*linkValuesAll[ir_lowerPct] + (ir_lowerPct - L*lowerPercentile)*linkValuesAll[il_lowerPct]
    upperPct = (L*upperPercentile- il_upperPct)*linkValuesAll[ir_upperPct] + (ir_upperPct - L*upperPercentile)*linkValuesAll[il_upperPct]
    print "(Links) Lower percentile %f pct at: %f; upper percentile %f pct at: %f" % (100*lowerPercentile, lowerPct, 100*upperPercentile, upperPct)

    
    kur = stats.kurtosis(linkValuesAll)
    print "(Links) Kurtosis is: %f" % kur
    
    #a straight histo plot: plt.hist(linkValuesAll, bins = binsLinks)
    #But: we want to transform the counts to log10(1+count)
    #Can be done like this:
#    counts, bins, patches = plt.hist(linkValuesAll, bins = binsLinks)
#    plt.close()
    counts, bins = np.histogram(linkValuesAll, bins = binsLinks) 
#    print counts
#    print bins
#    print "len bins:%d" % len(bins)
    #bin mid points
    midPts = [float(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]
    #width of bins:
    w = [float(bins[i+1] - bins[i]) for i in range(len(bins)-1)]
    #then make a bar plot:
    if linksPokesDistrInOnePlot_b == 1:
        fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1, ncols = 2)
        ax1_1 = ax1[0]
        titleDistr = 'Links'
        titleLinks = titleLinks
    elif linksDistrPlotsInSubPlots_b == 1:
        fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1, ncols = 3)
        ax1_1 = ax1[1]
        titleDistr = 'All'
    elif linksDistrPlotsInSubPlots_b == 2:
        if not(inFig): 
            fig1 = plt.figure(figsize = (14,8))
            ax1_1 = fig1.add_subplot(121)
            titleDistr = titleDistr
        else:
            fig1 = plt.figure(num = 1) #fetches the first of the list of existing figures
            ax1_1 = fig1.add_subplot(122)
            titleDistr = titleDistr
#        fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1, ncols = 2)
#        ax1_1 = ax1[0]
    else:
        fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1, ncols = 1)
        ax1_1 = ax1        
        titleDistr = 'All'

       
    ax1_1.bar([bins[i] for i in range(len(bins)-1)], np.log10(1+counts), width = w, color = color, fill = True, alpha = 0.8)
    #add a little x-axis to the left and to the right (multiples of 0.25):
    xLimPos = -1.00
    while(xLimPos - bins[len(bins)-1] < 0):
        xLimPos += 0.25
    xLimNeg = 1.00  
    while(bins[0] - xLimNeg < 0):
        xLimNeg -= 0.25
    ax1_1.set_xlim(xLimNeg,xLimPos)
    ax1_1.set_xlabel('writhe', fontsize = fontsizeLabel)
    ax1_1.set_ylabel('log(1+count)', fontsize = fontsizeLabel)
    ax1_1.tick_params(labelsize = fontsizeLabel)
    
    titleLinks = titleLinks
    if supTitle_b ==1:
        fig1.suptitle(titleLinks, fontsize = fontsizeTitle)
    ax1_1.set_title(titleDistr, fontsize = fontsizeSubTitle)


    linkValues.sort(reverse = True) #21-Aug-2016


    #plot top and bottom link values (unless links and pokes distr are made in one plot)
    if linksPokesDistrInOnePlot_b != 1:
        
        #linkValues.sort(reverse = True) 21-Aug-2016
        
        #plot top values
        if linksDistrPlotsInSubPlots_b ==1:
            ax1_2 = ax1[2] 
        elif linksDistrPlotsInSubPlots_b != 2:
            fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1, ncols = 1)
            ax1_2 = ax1
         
        if linksDistrPlotsInSubPlots_b != 2:
            selectedLinkValuesTop = [linkValues[i][0] for i in range(linksShowTopNr) ]
            ax1_2.hist(selectedLinkValuesTop, bins = binsLinks, color = color, fill = True, alpha = 0.8)
            titleLinks = 'Upper ' + str(linksShowTopNr) #+ ' highest writhe values for potential links'
            ax1_2.set_title(titleLinks, fontsize = fontsizeSubTitle)
            #ax1_2.set_xlabel('writhe', fontsize = fontsizeLabel)
            #set every 2nd xtick:
            j= 0
            xticks = ax1_2.get_xticks()
            xticksNew = []
            while(2*j < len(xticks)):
                xticksNew.append(xticks[2*j])
                j +=1
            ax1_2.set_xticks(xticksNew)
            ax1_2.set_ylabel('count', fontsize = fontsizeLabel)
            ax1_2.tick_params(labelsize = fontsizeLabel)
    
          
        #plot bottom link values
        if linksDistrPlotsInSubPlots_b ==1:
            ax1_3 = ax1[0] 
        elif linksDistrPlotsInSubPlots_b != 2:
            fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1, ncols = 1)
            ax1_3 = ax1 
    
        if linksDistrPlotsInSubPlots_b != 2:
            selectedLinkValuesBottom = [linkValues[::-1][i][0] for i in range(linksShowTopNr) ] #linkValues were sorted in reverse above
            ax1_3.hist(selectedLinkValuesBottom, bins = binsLinks, color = color, fill = True, alpha = 0.8)
            titleLinks = 'Lower ' + str(linksShowTopNr) #+ ' lowest writhe values for potential links'
            ax1_3.set_title(titleLinks, fontsize = fontsizeSubTitle)
            #ax1_3.set_xlabel('writhe', fontsize = fontsizeLabel)
            #set every 2nd xtick:
            j= 0
            xticks = ax1_3.get_xticks()
            xticksNew = []
            while(2*j < len(xticks)):
                xticksNew.append(xticks[2*j])
                j +=1
            ax1_3.set_xticks(xticksNew)
            ax1_3.set_ylabel('count', fontsize = fontsizeLabel)
            ax1_3.tick_params(labelsize = fontsizeLabel)

    ########################    
    #Pokes
    ########################    
        
    #Read in the data
    pokeValues, pokeValuesMin, pokeValuesMax = pokeValues    
        
    pokeValues.sort(reverse = True) #21-Aug-2016    
    pokeValuesAll = [pokeValues[i][0] for i in range(len(pokeValues))]
    print "Min writhe for pokes:%d" % min(pokeValuesAll)
    print "Max writhe for pokes:%d" % max(pokeValuesAll) 
    if pokesMinMax_b == 1:         
        pokeValuesMin = [pokeValuesMin[i][0] for i in range(len(pokeValuesMin)) ]
        pokeValuesMax = [pokeValuesMax[i][0] for i in range(len(pokeValuesMax)) ]

    #We want first to find the lower and upper percentiles as desired (approx'ly):
    L = len(pokeValuesAll)
    print "Number of potential pokes: %d" % L
    #sort the values:
    pokeValuesAll.sort(reverse = False)
    #find indexes "surrounding" the percentiles:
    il_lowerPct = int(np.floor(L*lowerPercentile))
    ir_lowerPct = il_lowerPct +1 # equals int(np.ceil(L*0.01)) unless L*0.01 is integral
    print "(Pokes) idx's left and right of lower percentile %f: %d, %d with values: %f, %f" % (lowerPercentile, il_lowerPct, ir_lowerPct, pokeValuesAll[il_lowerPct], pokeValuesAll[ir_lowerPct] )
    il_upperPct = int(np.floor(L*upperPercentile))
    ir_upperPct = il_upperPct +1 # equals int(np.ceil(L*0.99)) unless L*0.99 is integral
    print "(Pokes) idx's left and right of upper percentile %f: %d, %d with values: %f, %f" % (upperPercentile, il_upperPct, ir_upperPct, pokeValuesAll[il_upperPct], pokeValuesAll[ir_upperPct] )
    #let the desired percentiles be set as the weighted average over the values left and right of the percentile:
    lowerPct = (L*lowerPercentile- il_lowerPct) *pokeValuesAll[ir_lowerPct] + (ir_lowerPct - L*lowerPercentile)*pokeValuesAll[il_lowerPct]
    upperPct = (L*upperPercentile- il_upperPct)*pokeValuesAll[ir_upperPct] + (ir_upperPct - L*upperPercentile)*pokeValuesAll[il_upperPct]
    print "(Pokes) Lower percentile %f pct at: %f; upper percentile %f pct at: %f" % (100*lowerPercentile, lowerPct, 100*upperPercentile, upperPct)

    
 

    #straight histo plot: plt.hist(pokeValuesAll, bins = binsPokes)
    #But: we want to transform the counts to log10(1+count)
    #Can be done like this:
#    counts, bins, patches = plt.hist(pokeValuesAll, bins = binsPokes)
#    plt.close()
#    #bin mid points
#    midPts = [float(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]
#    #width of bins:
#    w = [float(bins[i+1] - bins[i]) for i in range(len(bins)-1)]
    #then make a bar plot:
    
    if linksPokesDistrInOnePlot_b == 1:
        
        ax1_2 = ax1[1]
        
        if pokesDistrLog_b == 1:
            #as above for the links, we want to transform the counts to log10(1+count):
            counts, bins = np.histogram(pokeValuesAll, bins = binsLinks)    
            #bin mid points
            midPts = [float(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]
            #width of bins:
            w = [float(bins[i+1] - bins[i]) for i in range(len(bins)-1)]
                    
            #then make a bar plot:
            if pokesMinMax_b == 1: #distinguish between the max and min values in plot
                countsMin, bins = np.histogram(pokeValuesMin, bins = bins)
                countsMax, bins = np.histogram(pokeValuesMax, bins = bins)
                        
                ax1_2.bar(midPts, np.log10(1+countsMin), width = w, color = colorMin, alpha = alphaMin)
                ax1_2.bar(midPts, np.log10(1+countsMax), width = w, color = colorMax, alpha = alphaMax)

            else: #don't distinguish min's and max's in plot
                ax1_2.bar(midPts, np.log10(1+counts), width = w, color = color)

            ax1_2.set_ylabel('log(1+count)', fontsize = fontsizeLabel)


        else: #just make std histogram
        
            if pokesMinMax_b == 1: #distinguish between the max and min values in plot
                ax1_2.hist(pokeValuesMin, bins = bins, color = colorMin, alpha = alphaMin, fill = True)
                ax1_2.hist(pokeValuesMax, bins = bins, color = colorMax, alpha = alphaMax, fill = True)
            else: #don't distinguish min's and max's in plot
                ax1_2.hist(pokeValuesAll, bins = binsPokes, color = color)
                        
            ax1_2.set_ylabel('count', fontsize = fontsizeLabel)

        ax1_2.set_xlabel('writhe', fontsize = fontsizeLabel)
        ax1_2.tick_params(labelsize = fontsizeLabel)
        titleDistr = 'Pokes'
        ax1_2.set_title(titleDistr, fontsize = fontsizeTitle)
        
        
    else: #links and pokes distr not shown in one plot 
        if pokesDistrPlotsInSubPlots_b ==1:
            fig2, ax2 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 3)
            ax2_1 = ax2[1]
            #Re y-axis label: only set for left-most plot (except when using log-count in mid plot )
        else:
            fig2, ax2 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 1)
            ax2_1 = ax2
            ax2_1.set_ylabel('count', fontsize = fontsizeLabel)

    #    plt.bar(midPts, np.log10(1+counts), width = w)
        if pokesMinMax_b == 1: #distinguish between the max and min values in plot
            ax2_1.hist(pokeValuesMin, bins = binsPokes, color = colorMin, alpha = alphaMin, fill = True)
            ax2_1.hist(pokeValuesMax, bins = binsPokes, color = colorMax, alpha = alphaMax, fill = True)    
        else: # don't distinguish
            ax2_1.hist(pokeValuesAll, bins = binsPokes, color = color)
        ax2_1.set_title('All', fontsize = fontsizeSubTitle)
#        ax2_1.set_ylabel('count', fontsize = fontsizeLabel)
    #    plt.ylabel('log(1+count)')
        ax2_1.set_xlabel('writhe', fontsize = fontsizeLabel)
        ax2_1.tick_params(labelsize = fontsizeLabel)
        titlePokes = titlePokesI
        fig2.suptitle(titlePokes, fontsize = fontsizeTitle)
    
    
        #plot top poke values   
        if pokesDistrPlotsInSubPlots_b ==1:
            ax2_2 = ax2[2]
        else:
            fig2, ax2 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 1)
            ax2_2 = ax2
            
        selectedPokeValues = [pokeValues[i][0] for i in range(pokesShowTopNr) ]
        ax2_2.hist(selectedPokeValues, bins = binsPokes, color = color)
        titlePokes= 'Upper ' + str(pokesShowTopNr) #+ titlePokesII
        ax2_2.set_title(titlePokes, fontsize = fontsizeSubTitle)
        #ax2_2.set_xlabel('writhe', fontsize = fontsizeLabel)
        #set every 2nd xtick:
        j= 0
        xticks = ax2_2.get_xticks()
        xticksNew = []
        while(2*j < len(xticks)):
            xticksNew.append(xticks[2*j])
            j +=1
        ax2_2.set_xticks(xticksNew)
#        ax2_2.set_ylabel('count', fontsize = fontsizeLabel)
        ax2_2.tick_params(labelsize = fontsizeLabel)
 

        #plot bottom poke values   
        if pokesDistrPlotsInSubPlots_b ==1:
            ax2_3 = ax2[0]
        else:
            fig2, ax2 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 1)
            ax2_3 = ax2
            
        selectedPokeValues = [pokeValues[::-1][i][0] for i in range(pokesShowTopNr) ]
        ax2_3.hist(selectedPokeValues, bins = binsPokes, color = color)
        titlePokes= 'Lower ' + str(pokesShowTopNr) #+ titlePokesII
        ax2_3.set_title(titlePokes, fontsize = fontsizeSubTitle)
        #ax2_3.set_xlabel('writhe', fontsize = fontsizeLabel)
        #set every 2nd xtick:
        j= 0
        xticks = ax2_3.get_xticks()
        xticksNew = []
        while(2*j < len(xticks)):
            xticksNew.append(xticks[2*j])
            j +=1
        ax2_3.set_xticks(xticksNew)
        ax2_3.set_ylabel('count', fontsize = fontsizeLabel)
        ax2_3.tick_params(labelsize = fontsizeLabel)


    
    #print out the top 10 potential links and pokes:
    #load segment chains if 3d-plots desired:
    if (plotLinksTop10_b == 1 or plotPokesTop10_b ==1 or plotLinkExs_b == 1 or  plotPokeExs_b == 1 or plotExsInInt_b == 1):
        chains = readResultsCcode_chains(inputFileName = chainsFile)    
        
    print "Top 5 link values:"
    print linkValues[:5]  #Obs: sorted above (btw: moved the sorting on 21-Aug2016)
    print "Top 5 lowest link values:"
    print linkValues[::-1][:5] 
    
    pymolPlotClosedExamples(outputFilePath = r'C:\Users\Christian\Bioinformatics\papers\geometric_anomalies_in_folds\c_code\results_C_4th\double_precision\pymolTest.txt' ,examples = linkValues[:5]);
    
    if plotLinksTop10_b == 1:
        plotClosedLoopExamples(examples = linkValues[:5], chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern) 
        plotClosedLoopExamples(examples = linkValues[::-1][:5], chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern) 
    print "Top 5 poke values:"
    print pokeValues[:5] #Obs: sorted above (btw: moved the sorting on 21-Aug2016)
    print "Top 5 lowest poke values:"
    print pokeValues[::-1][:5] #Obs: sorted above (btw: moved the sorting on 21-Aug2016)
    if plotPokesTop10_b ==1: 
        plotClosedLoopExamples(examples = pokeValues[:5], chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern) 
        plotClosedLoopExamples(examples = pokeValues[::-1][:10], chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern) 

    #number of links below/above -cutoff/cutoff and plots of these:
    linkValuesExtremePos = []
    linkValuesExtremeNeg = []
    linkValuesInterval = []
    pokeValuesExtremeNeg = []    
    pokeValuesExtremePos = []
    pokeValuesInterval = []

    cntPlus = 0
    cntNeg = 0
    cntInInt = 0
    for ex in linkValues:
        val = ex[0]
        if val > cutOffLinks:
            linkValuesExtremePos.append(ex)
            cntPlus +=1
        if val < -cutOffLinks:
            linkValuesExtremeNeg.append(ex)
            cntNeg +=1
        if (val > intervalLinks[0] and val < intervalLinks[1]):
            linkValuesInterval.append(ex)
            cntInInt +=1
    print "There were %d links having a writhe value below the cutoff: %f" % (cntNeg, -cutOffLinks)
    print "There were %d links having a writhe value above the cutoff: %f" % (cntPlus, cutOffLinks)
    print "There were %d links having a writhe value in the interval: %s" % (cntInInt, str(intervalLinks))


    if plotLinkExs_b == 1:
        plotClosedLoopExamples(examples = linkValuesExtremePos, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
        plotClosedLoopExamples(examples = linkValuesExtremeNeg, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
    if plotExsInInt_b == 1:
        plotClosedLoopExamples(examples = linkValuesInterval, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
        

    cntExPokesPos = 0
    cntExPokesNeg = 0
    cntInInt = 0 #reset
    for ex in pokeValues:
        val = ex[0]
        if val > cutOffPokes:
            pokeValuesExtremePos.append(ex)
            cntExPokesPos +=1
        if val < -cutOffPokes:
            pokeValuesExtremeNeg.append(ex)
            cntExPokesNeg +=1
        if (val > intervalPokes[0] and val < intervalPokes[1]):
            pokeValuesInterval.append(ex)
            cntInInt +=1
    print "There were %d pokes having a writhe value below the cutoff: %f" % (cntExPokesNeg, -cutOffPokes)
    print "There were %d pokes having a writhe value above the cutoff: %f" % (cntExPokesPos, cutOffPokes)
    print "There were %d pokes having a writhe value in the interval: %s" % (cntInInt, str(intervalPokes))
    if plotPokeExs_b == 1:
        plotClosedLoopExamples(examples = pokeValuesExtremeNeg, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
        plotClosedLoopExamples(examples = pokeValuesExtremePos, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
    if plotExsInInt_b == 1:
        plotClosedLoopExamples(examples = pokeValuesInterval, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    


    plt.show()

    
    if linksDistrPlotsInSubPlots_b ==2:
        return fig1
    

 
def histoSubChainPairs(subChainPairValues,
                       logCount_b = 0,
                       color = 'lightgrey',
                       titleDistr = 'Distr. of extreme writhe values, unrestricted search',
                       fontsizeLabel = "large",
                       fontsizeSubTitle = "large",
                       fontsizeTitle = "large",
                       subChainLgth = 30,
                       bins = 10, 
                       showTopNr = 10, 
                       cutOff = 11.0,  
                       chainsFile = '', 
                       addSubplotToFirstPlot_b = 0,
                       addSubPlotNr = 122,
                       distrPlotsInSubPlots_b = 0, 
                       plotExs_b = 0, 
                       plotTop10_b = 0, 
                       interval = [2,3], 
                       plotExsInInt_b = 0, 
                       elevAzimList =  [[15,-90], [15,0], [90,0]], 
                       color1 = 'red', 
                       markerAtom1 = 'o', 
                       markerLine1 = '-', 
                       color2 = 'blue', 
                       markerAtom2 = 'o', 
                       markerLine2 = '-', 
                       stringPattern = '([\S]+)([^\w]+)([\w]+)\.([\w]+)'):
    '''Plot histograms of the distributions of the writhe numbers among mutual 
    writhes for pairs of sub-chains. Makes 3d plots of cases for which the 
    writhe is outside the interval [-cutOff, cutOff], i.e. of the "extreme 
    extreme writhe" cases; the relevant pairs of subChains are highlighted 
    in each case. Examples with writhe values inside a set interval can also 
    be chosen (use plotExsInInt_b =1).
    
    addSubplotToFirstPlot_b == 1: will add the distr of writhe values as  a subplot 
    (122) to an existing plot. To be used for adding this distr to e.g. a plot of 
    the distr of the writhe values from a restricted search (using closedLoops ...).'''


    #plot distributions of all values:
    extremeValuesAll = [subChainPairValues[i][0] for i in range(len(subChainPairValues)) ]

    #We want first to find the lower and upper percentiles as desired (approx'ly):
    L = len(extremeValuesAll)
    print "Number of sub-chain pairs: %d" % L
    #sort the values:
    extremeValuesAll.sort(reverse = False)
    #find indexes "surrounding" the percentiles:
    il_lowerPct = int(np.floor(L*lowerPercentile))
    ir_lowerPct = il_lowerPct +1 # equals int(np.ceil(L*0.01)) unless L*0.01 is integral
    print "idx's left and right of lower percentile %f: %d, %d with values: %f, %f" % (lowerPercentile, il_lowerPct, ir_lowerPct, extremeValuesAll[il_lowerPct], extremeValuesAll[ir_lowerPct] )
    il_upperPct = int(np.floor(L*upperPercentile))
    ir_upperPct = il_upperPct +1 # equals int(np.ceil(L*0.99)) unless L*0.99 is integral
    print "idx's left and right of upper percentile %f: %d, %d with values: %f, %f" % (upperPercentile, il_upperPct, ir_upperPct, extremeValuesAll[il_upperPct], extremeValuesAll[ir_upperPct] )
    #let the desired percentiles be set as the weighted average over the values left and right of the percentile:
    lowerPct = (L*lowerPercentile- il_lowerPct) *extremeValuesAll[ir_lowerPct] + (ir_lowerPct - L*lowerPercentile)*extremeValuesAll[il_lowerPct]
    upperPct = (L*upperPercentile- il_upperPct)*extremeValuesAll[ir_upperPct] + (ir_upperPct - L*upperPercentile)*extremeValuesAll[il_upperPct]
    print "Lower percentile %f pct at: %f; upper percentile %f pct at: %f" % (100*lowerPercentile, lowerPct, 100*upperPercentile, upperPct)
    
    kur = stats.kurtosis(extremeValuesAll)
    print "Kurtosis is: %f" % kur

    if addSubplotToFirstPlot_b != 1:
#        title = 'Distr. of extreme writhe values for pairs of sub-chains of lgth '+ str(subChainLgth) 
        title = 'Sub-chains of lgth '+ str(subChainLgth) 
        if distrPlotsInSubPlots_b ==1:
            fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 3)
            ax_1 = ax1[1]
        else:
            fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 1)
            ax_1 = ax1
    else: #if addSubplotToFirstPlot_b == 1
        distrPlotsInSubPlots_b = 2
        fig1 = plt.figure(num = 1) #fetches the first of the list of existing figures
        ax_1 = fig1.add_subplot(addSubPlotNr)
#        ax_1 = ax1[1]
        title = 'Sub-chains of lgth '+ str(subChainLgth) 
        
#    if distrPlotsInSubPlots_b ==1:
#        fig.add_subplot(132)
#        title = 'Distr. of extreme writhe values for pairs of sub-chains of lgth '+ str(subChainLgth) 
#
#    if distrPlotsInSubPlots_b ==2:
#        fig.add_subplot(122)
#        title = 'Unrestricted search'
        
#    ax_1 = ax_1.plot(extremeValuesAll)
    if logCount_b != 1:
    
        ax_1.hist(extremeValuesAll, bins = bins, color = color)
        ax_1.set_ylabel('count', fontsize = fontsizeLabel)
        
    elif logCount_b == 1:
        #we want to transform the counts to log10(1+count)
        #Can be done like this:
        counts, bins = np.histogram(extremeValuesAll, bins = bins)
#        plt.close()
        #bin mid points
        midPts = [float(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]
        #width of bins:
        w = [float(bins[i+1] - bins[i]) for i in range(len(bins)-1)]

        ax_1.bar(midPts, np.log10(1+counts), width = w, color = color)
        ax_1.set_ylabel('log(1+count)', fontsize = fontsizeLabel)
        
    ax_1.set_title(title, fontsize = fontsizeSubTitle)
    plt.suptitle(titleDistr, fontsize = fontsizeTile)
    ax_1.set_xlabel('writhe', fontsize = fontsizeLabel)
    ax_1.tick_params(labelsize = fontsizeLabel)
    if addSubplotToFirstPlot_b != 0:
        plt.figure(num = 1, figureClass = fig1) #will add the plot to the existing (fig1)
    plt.show()
    
    #sort before carrying on:
    subChainPairValues.sort(reverse = True)
    
    #plot top values
    if showTopNr > 0:
        if distrPlotsInSubPlots_b ==1:
            ax_2 = ax1[2]
        else:
            fig2, ax2 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 1)
            ax_2 = ax2
     
        selectedValuesTop = [subChainPairValues[i][0] for i in range(showTopNr) ]
        ax_2.hist(selectedValuesTop, bins = bins, color = color)
        title = 'Top ' + str(showTopNr) + ' highest writhe values for pairs of sub-chains'
        ax_2.set_title(title, fontsize =fontsizeSubTitle)
        #plot bottom link values
        if distrPlotsInSubPlots_b ==1:
            ax_3 = ax1[0]
        else:
            fig3, ax3 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 1)
            ax_3 = ax3
        selectedValuesBottom = [subChainPairValues[::-1][i][0] for i in range(showTopNr) ]
        ax_3.hist(selectedValuesBottom, bins = bins, color = color)
        title = 'Top ' + str(showTopNr) + ' lowest writhe values for pairs of sub-chains'
        ax_3.set_title(title, fontsize =fontsizeSubTitle)
        
#    plt.hist(extremeValuesAll, bins = bins, color = color)
#    plt.title(title)
#    plt.show()
#    #plot top values
#    subChainPairValues.sort(reverse = True)
#    if distrPlotsInSubPlots_b ==1:
#        fig.add_subplot(133)
#    else:
#        plt.figure()
#    selectedValuesTop = [subChainPairValues[i][0] for i in range(showTopNr) ]
#    plt.hist(selectedValuesTop, bins = bins, color = color)
#    title = 'Top ' + str(showTopNr) + ' highest writhe values for pairs of sub-chains'
#    plt.title(title)
#    plt.show()
#    #plot bottom link values
#    if distrPlotsInSubPlots_b ==1:
#        fig.add_subplot(131)
#    else:
#        plt.figure()
#    selectedValuesBottom = [subChainPairValues[::-1][i][0] for i in range(showTopNr) ]
#    plt.hist(selectedValuesBottom, bins = bins, color = color)
#    titleLinks = 'Top ' + str(showTopNr) + ' lowest writhe values for pairs of sub-chains'
#    plt.title(title)
#    plt.show()
    
    #print out the top 10 cases:
    #load segment chains if 3d-plots desired:
    if (plotTop10_b == 1 or plotExs_b == 1 or plotExsInInt_b == 1):
        chains = readResultsCcode_chains(inputFileName = chainsFile)    
    
    print "Top 10 values (sub-chains):"
    print subChainPairValues[:10]  
    print "Top 10 lowest values (sub-chains):"
    print subChainPairValues[::-1][:10]    
    if plotTop10_b == 1:
        plotClosedLoopExamples(examples = subChainPairValues[:5], chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern) 
        plotClosedLoopExamples(examples = subChainPairValues[::-1][:5], chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern) 
  
    #number of cases below/above -cutoff/cutoff and plots of these:
    extremePos = []
    extremeNeg = []
    valuesInterval = []

    cntPlus = 0
    cntNeg = 0
    cntInInt = 0
    for ex in subChainPairValues:
        val = ex[0]
        if val > cutOff:
            extremePos.append(ex)
            cntPlus +=1
        if val < -cutOff:
            extremeNeg.append(ex)
            cntNeg +=1
        if (val > interval[0] and val < interval[1]):
            valuesInterval.append(ex)
            cntInInt +=1
    print "There were %d cases having a writhe value below the cutoff: %f" % (cntNeg, -cutOff)
    print "There were %d cases having a writhe value above the cutoff: %f" % (cntPlus, cutOff)
    print "There were %d cases having a writhe value in the interval: %s" % (cntInInt, str(interval))

    if plotExs_b == 1:
        plotClosedLoopExamples(examples = extremePos, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
        plotClosedLoopExamples(examples = extremeNeg, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
    if plotExsInInt_b == 1:
        plotClosedLoopExamples(examples = valuesInterval, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
        
    plt.show()
    
 
         
#version of histoSubChainPairs that distinguishes in the plot the distr's of min/max values: 
def histoSubChainPairs_minMax(subChainPairValues,
                       figIn = '', 
                       axisIn = '',
                       yLabel_b = 1,
                       logCount_b = 0,
                       color = 'grey',
                       colorMin = 'lightgrey',
                       colorMax = 'grey',
                       fontsizeLabel = "large",
                       fontsizeSubTitle = "large",
                       fontsizeTitle = "large",
                       alphaMin = 1.0,
                       alphaMax = 0.8,
                       invName = 'writhe',
                       supTitle_b = 1,
                       titleDistr = 'Distr. of extreme writhe values, unrestricted search',
                       subChainLgth = 30,
                       lowerPercentile = 0.01,
                       upperPercentile = 0.99,
                       bins = 10, 
                       showTopNr = 10,
                       cutOff = 0.9,  
                       chainsFile = '', 
                       addSubplotToFirstPlot_b = 0,
                       addSubPlotNr = 122,
                       distrPlotsInSubPlots_b = 0,
                       few_xticks_b = 0,
                       ticksSize = 'large',
                       plotExs_b = 0, 
                       plotTop10_b = 0, 
                       plotTop4_b = 0,
                       interval = [2,3], 
                       plotExsInInt_b = 0, 
                       elevAzimList =  [[15,-90], [15,0], [90,0]], 
                       color1 = 'red', 
                       markerAtom1 = 'o', 
                       markerLine1 = '-', 
                       color2 = 'blue', 
                       markerAtom2 = 'o', 
                       markerLine2 = '-', 
                       stringPattern = '([\S]+)([^\w]+)([\w]+)\.([\w]+)'):
    '''Plot histograms of the distributions of the writhe numbers (or: invariant values) 
    for pairs of sub-chains. Makes 3d plots of cases for top4 or top10 writhe value cases
    or cases for which the writhe is outside the interval [-cutOff, cutOff], i.e. 
    of the "exceptional extreme writhe" cases; the relevant pairs of subChains are 
    highlighted in each case. Examples with writhe values inside a set interval 
    can also be chosen (use plotExsInInt_b =1).
    
    Input:
       subChainPairValues: the primary input data. The structure must be as from the
       readResultsCcode_closedLoops with linksAndPokes_b = 0: first read in the data
       
       dict_X, extremeWrithes_X_All, extremeWrithes_X_Min, extremeWrithes_X_Max  = readResultsCcode_closedLoops(.., linksAndPokes_b = 0, ...)

       then convert to the structure expected in this funtion:
       subChainPairValues_X = extremeWrithes_X_All, extremeWrithes_X_Min, extremeWrithes_X_Max
       
       For explicit examples see the Usage section of this module.

       logCount_b (0): if set to 1 the frequency (y-axis) in the distr's will be done by log(count +1)

       lowerPercentile (0.01): see histoClosedLoops_minMax 
       upperPercentile (0.99): see histoClosedLoops_minMax


       figIn ('') and axisIn (''): parameters allowing to input an existing figure and an axis 
       of it (allows to pass an axis for a subplot of the figure, figIn). See also figIn in histoClosedLoops_minMax 
       
 
       addSubplotToFirstPlot_b (0): will add the distr of writhe values as a subplot (as 122 or 
       use addSubPlotNr) to an existing plot. To be used for adding this distr to e.g. a plot of 
       the distr of the writhe values from a restricted search (ie C-outpu for closedLoops and plotted 
       with histoClosedLoops_minMax)
       addSubPlotNr (122): the sub-plot number mentioned used with addSubplotToFirstPlot_b = 1

       bins (10): numbe rof bins used in distribution plot
       showTopNr (10): if this integer is set > 0, the lower and upper tails counting showTopNr examples will be shown in plots flanking the full distribution, if distrPlotsInSubPlots_b = 1, and else as separate plots 
       distrPlotsInSubPlots_b (0): see showTopNr

       few_xticks_b (0): set only three x-ticks: at first bin, at 0 and at last bin
       yLabel_b (1): if set to 1 the y-label is shown, else not.
       ticksSize ('large'): font size of tick labels
              
       titleDistr ('Distr. of extreme writhe values, unrestricted search'): see histoClosedLoops_minMax 
       supTitle_b (1): see histoClosedLoops_minMax 
       subChainLgth (30): will be included in the titles of the distributions (not in all cases)
       color ('grey'): see histoClosedLoops_minMax 
       colorMin ('lightgrey'): see histoClosedLoops_minMax 
       colorMax ('grey'): see histoClosedLoops_minMax 
       fontsizeLabel ("large"): see histoClosedLoops_minMax 
       fontsizeSubTitle ("large"): see histoClosedLoops_minMax 
       fontsizeTitle ("large"): see histoClosedLoops_minMax 
       alphaMin (1.0): see histoClosedLoops_minMax 
       alphaMax (0.8): see histoClosedLoops_minMax 
       invName ('writhe'): to allow indicating the invariant for which the data were obtained
  
       cutOff (0.9): cases having a writhe below/above -cutoff/cutoff will be shown in 3d-plots        
       chainsFile (''): see histoClosedLoops_minMax
       plotExs_b (0): see histoClosedLoops_minMax 
       plotTop10_b (0): like plotLinksTop10_b in see histoClosedLoops_minMax 
       plotTop4_b (0): like plotTop10_b, but only top4 examples are shown
       interval ([2:3]): see intervalLinks in histoClosedLoops_minMax
       plotExsInInt_b (0): see histoClosedLoops_minMax
       elevAzimList ( [[15,-90], [15,0], [90,0]]): see histoClosedLoops_minMax
       color1 ('red'): see histoClosedLoops_minMax
       markerAtom1 ('o'): see histoClosedLoops_minMax
       markerLine1 ('-'): see histoClosedLoops_minMax
       color2 ('blue'): see histoClosedLoops_minMax
       markerAtom2 ('o'): see histoClosedLoops_minMax
       markerLine2 ('-'): see histoClosedLoops_minMax
       stringPattern ('([\S]+)([^\w]+)([\w]+)\.([\w]+)'): see histoClosedLoops_minMax 
    '''

    #read input in proper structure:
    subChainPairValuesAll, subChainPairValuesMin, subChainPairValuesMax = subChainPairValues

    #plot distributions of all values:
    extremeValuesAll = [subChainPairValuesAll[i][0] for i in range(len(subChainPairValuesAll)) ]
    extremeValuesMin = [subChainPairValuesMin[i][0] for i in range(len(subChainPairValuesMin)) ]        
    extremeValuesMax = [subChainPairValuesMax[i][0] for i in range(len(subChainPairValuesMax)) ]
    
    #We want first to find the lower and upper percentiles as desired (approx'ly):
    L = len(extremeValuesAll)
    print "Number of sub-chain pairs: %d" % L
    #sort the values:
    extremeValuesAll.sort(reverse = False)
    #find indexes "surrounding" the percentiles:
    il_lowerPct = int(np.floor(L*lowerPercentile))
    ir_lowerPct = il_lowerPct +1 # equals int(np.ceil(L*0.01)) unless L*0.01 is integral
    print "idx's left and right of lower percentile %f: %d, %d with values: %f, %f" % (lowerPercentile, il_lowerPct, ir_lowerPct, extremeValuesAll[il_lowerPct], extremeValuesAll[ir_lowerPct] )
    il_upperPct = int(np.floor(L*upperPercentile))
    ir_upperPct = il_upperPct +1 # equals int(np.ceil(L*0.99)) unless L*0.99 is integral
    print "idx's left and right of upper percentile %f: %d, %d with values: %f, %f" % (upperPercentile, il_upperPct, ir_upperPct, extremeValuesAll[il_upperPct], extremeValuesAll[ir_upperPct] )
    #let the desired percentiles be set as the weighted average over the values left and right of the percentile:
    lowerPct = (L*lowerPercentile- il_lowerPct) *extremeValuesAll[ir_lowerPct] + (ir_lowerPct - L*lowerPercentile)*extremeValuesAll[il_lowerPct]
    upperPct = (L*upperPercentile- il_upperPct)*extremeValuesAll[ir_upperPct] + (ir_upperPct - L*upperPercentile)*extremeValuesAll[il_upperPct]
    print "Lower percentile %f pct at: %f; upper percentile %f pct at: %f" % (100*lowerPercentile, lowerPct, 100*upperPercentile, upperPct)
    
    kurAll = stats.kurtosis(extremeValuesAll)
    print "All: Kurtosis is: %f" % kurAll
    
    kurMin = stats.kurtosis(extremeValuesMin)
    print "Min: Kurtosis is: %f" % kurMin
    
    kurMax = stats.kurtosis(extremeValuesMax)
    print "Max: Kurtosis is: %f" % kurMax
    
    
    if addSubplotToFirstPlot_b != 1:
#        title = 'Distr. of extreme writhe values for pairs of sub-chains of lgth '+ str(subChainLgth) 
        if distrPlotsInSubPlots_b ==1:
            title = 'Sub-chains of lgth '+ str(subChainLgth) 
            if not(figIn):
                fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 3)   
                ax_1 = ax1[1]
            else:
                print("figIn should not be input when having distrPlotsInSubPlots_b =1")
            if supTitle_b == 1:
                ax_1.set_title(title, fontsize = fontsizeTitle)
        else:
            if figIn:
                title = invName
                fig1 = figIn
                ax_1 = axisIn
            else:
                fig1, ax1 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 1)
                ax_1 = ax1
    else: #if addSubplotToFirstPlot_b == 1
        distrPlotsInSubPlots_b = 2
        fig1 = plt.figure(num = 1) #fetches the first of the list of existing figures
        ax_1 = fig1.add_subplot(addSubPlotNr)
#        ax_1 = ax1[1]
        title = 'Sub-chains of lgth '+ str(subChainLgth) #Unrestricted search
        
#    if distrPlotsInSubPlots_b ==1:
#        fig.add_subplot(132)
#        title = 'Distr. of extreme writhe values for pairs of sub-chains of lgth '+ str(subChainLgth) 
#
#    if distrPlotsInSubPlots_b ==2:
#        fig.add_subplot(122)
#        title = 'Unrestricted search'
        
#    ax_1 = ax_1.plot(extremeValuesAll)
    if logCount_b != 1:
        
        #we want to plot the distr of the min vals and the distr of the max vals in one plot; 
        #so must use common bins for the two; we use bins derived for all data:
        #Can be done like this:
        counts, bins = np.histogram(extremeValuesAll, bins = bins)
       
        ax_1.hist(extremeValuesMin, bins = bins, color = colorMin, alpha = alphaMin, fill = True, label = invName)
        ax_1.hist(extremeValuesMax, bins = bins, color = colorMax, alpha = alphaMax, fill = True, label = invName)

#        ax_1.legend(fontsize = "xx-small" , markerscale = 0.5)
        ax_1.set_title(title, fontsize = fontsizeTitle)

        if yLabel_b == 1:
            ax_1.set_ylabel('count', fontsize = fontsizeLabel)
                
        if few_xticks_b == 1:
            ax_1.set_xticks((bins[0], 0, bins[len(bins)-1]))
            
    elif logCount_b == 1:
        #we want to transform the counts to log10(1+count)
        #Can be done like this:
        counts, bins = np.histogram(extremeValuesAll, bins = bins)
#        plt.close()
        #bin mid points
        midPts = [float(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]
        #width of bins:
        w = [float(bins[i+1] - bins[i]) for i in range(len(bins)-1)]

        countsMin, bins = np.histogram(extremeValuesMin, bins = bins)
        countsMax, bins = np.histogram(extremeValuesMax, bins = bins)
                        
        ax_1.bar(midPts, np.log10(1+countsMin), width = w, color = colorMin, alpha = alphaMin, label = invName)
        ax_1.bar(midPts, np.log10(1+countsMax), width = w, color = colorMax, alpha = alphaMax, label = invName)
        
#        ax_1.legend(fontsize = "xx-small" , markerscale = 0.5)
        ax_1.set_title(title, fontsize = fontsizeTitle)
       
        if yLabel_b == 1:
            ax_1.set_ylabel('log(1+count)', fontsize = fontsizeLabel)

        if few_xticks_b == 1:
            ax_1.set_xticks((bins[0], 0, bins[len(bins)-1]))

        
#    ax_1.set_title(title, fontsize = "large")
    if not(figIn):
        if supTitle_b == 1:
            plt.suptitle(titleDistr, fontsize = fontsizeTitle)
        ax_1.set_xlabel(invName, fontsize = fontsizeLabel)
    #add a little x-axis to the left and to the right (multiples of 0.25):
    xLimPos = -1.00
    while(xLimPos - bins[len(bins)-1] < 0):
        xLimPos += 0.25
    xLimNeg = 1.00  
    while(bins[0] - xLimNeg < 0):
        xLimNeg -= 0.25
        ax_1.set_xlim(xLimNeg,xLimPos)
    ax_1.tick_params(labelsize = ticksSize)
    if addSubplotToFirstPlot_b != 0:
        plt.figure(num = 1, figureClass = fig1) #will add the plot to the existing (fig1)
#    plt.tight_layout
    plt.show()
    
    #sort before carrying on:
    subChainPairValuesAll.sort(reverse = True)
    
    #plot top values
    if showTopNr > 0:
        if distrPlotsInSubPlots_b ==1:
            ax_2 = ax1[2]
        else:
            fig2, ax2 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 1)
            ax_2 = ax2
     
        selectedValuesTop = [subChainPairValuesAll[i][0] for i in range(showTopNr) ]
        ax_2.hist(selectedValuesTop, bins = bins, color = color)
        title = 'Top ' + str(showTopNr) + ' highest ' + invName + ' values for pairs of sub-chains'
        ax_2.set_title(title)
        #plot bottom link values
        if distrPlotsInSubPlots_b ==1:
            ax_3 = ax1[0]
        else:
            fig3, ax3 = plt.subplots(figsize = (14,8), nrows = 1 , ncols = 1)
            ax_3 = ax3
        selectedValuesBottom = [subChainPairValuesAll[::-1][i][0] for i in range(showTopNr) ]
        ax_3.hist(selectedValuesBottom, bins = bins, color = color)
        title = 'Top ' + str(showTopNr) + ' lowest ' + invName + ' values for pairs of sub-chains'
        ax_3.set_title(title)
        
#    plt.hist(extremeValuesAll, bins = bins, color = color)
#    plt.title(title)
#    plt.show()
#    #plot top values
#    subChainPairValues.sort(reverse = True)
#    if distrPlotsInSubPlots_b ==1:
#        fig.add_subplot(133)
#    else:
#        plt.figure()
#    selectedValuesTop = [subChainPairValues[i][0] for i in range(showTopNr) ]
#    plt.hist(selectedValuesTop, bins = bins, color = color)
#    title = 'Top ' + str(showTopNr) + ' highest writhe values for pairs of sub-chains'
#    plt.title(title)
#    plt.show()
#    #plot bottom link values
#    if distrPlotsInSubPlots_b ==1:
#        fig.add_subplot(131)
#    else:
#        plt.figure()
#    selectedValuesBottom = [subChainPairValues[::-1][i][0] for i in range(showTopNr) ]
#    plt.hist(selectedValuesBottom, bins = bins, color = color)
#    titleLinks = 'Top ' + str(showTopNr) + ' lowest writhe values for pairs of sub-chains'
#    plt.title(title)
#    plt.show()
    
    #print out the top 10 cases:
    #load segment chains if 3d-plots desired:
    if (plotTop4_b == 1 or plotTop10_b == 1 or plotExs_b == 1 or plotExsInInt_b == 1):
        chains = readResultsCcode_chains(inputFileName = chainsFile)    
    
    if (plotTop4_b == 1 or plotTop10_b ==1): 
        if plotTop4_b == 1:
            m = 2
        if plotTop10_b == 1:
            m = 5
        print "Top %d values:" % m
        print subChainPairValuesAll[:m]   
        print "Top %d lowest link values:" % m
        print subChainPairValuesAll[::-1][:m] 
        #plot the examples:
        plotClosedLoopExamples(examples = subChainPairValuesAll[:m], invName = invName, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern) 
        plotClosedLoopExamples(examples = subChainPairValuesAll[::-1][:m], invName = invName, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern) 
  
    #number of cases below/above -cutoff/cutoff and plots of these:
    extremePos = []
    extremeNeg = []
    valuesInterval = []

    cntPlus = 0
    cntNeg = 0
    cntInInt = 0
    for ex in subChainPairValuesAll:
        val = ex[0]
        if val > cutOff:
            extremePos.append(ex)
            cntPlus +=1
        if val < -cutOff:
            extremeNeg.append(ex)
            cntNeg +=1
        if (val > interval[0] and val < interval[1]):
            valuesInterval.append(ex)
            cntInInt +=1
    print "There were %d cases having a %s value below the cutoff: %f" % (cntNeg, invName, -cutOff)
    print "There were %d cases having a %s value above the cutoff: %f" % (cntPlus, invName, cutOff)
    print "There were %d cases having a %s value in the interval: %s" % (cntInInt, invName, str(interval))

    if plotExs_b == 1:
        plotClosedLoopExamples(examples = extremePos, invName = invName, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
        plotClosedLoopExamples(examples = extremeNeg, invName = invName, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
    if plotExsInInt_b == 1:
        plotClosedLoopExamples(examples = valuesInterval, invName = invName, chainsDict = chains, elevAzimList = elevAzimList, color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1, color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, stringPattern = stringPattern)    
        
    plt.show()

#multi-call to histoSubChainPairs_minMax, for use when input contains sub-chain pair
#values for several invariants:
def multiInv_histoSubChainPairs_minMax(subChainPairValues,
                                       rowsFig = 1,
                                       colsFig = 1,
                                       invList = '',
                                       yLabelLeft_b = 1,
                                       logCount_b = 0,
                                       color = 'lightgrey',
                                       colorMin = 'lightgrey',
                                       colorMax = 'grey',
                                       alphaMin = 1.0,
                                       alphaMax = 0.8,
                                       titleDistr = '',
                                       lowerPercentile = 0.01,
                                       upperPercentile = 0.99,
                                       subChainLgth = 30,
                                       bins = 10, 
                                       cutOff = 11.0,  
                                       chainsFile = '', 
                                       addSubPlotNr = 122,
                                       distrPlotsInSubPlots_b = 0,
                                       few_xticks_b = 0,
                                       ticksSize = 'small',
                                       plotExs_b = 0, 
                                       plotTop10_b = 0,
                                       plotTop4_b = 0,
                                       interval = [2,3], 
                                       plotExsInInt_b = 0, 
                                       elevAzimList =  [[15,-90], [15,0], [90,0]], 
                                       color1 = 'red', 
                                       markerAtom1 = 'o', 
                                       markerLine1 = '-', 
                                       color2 = 'blue', 
                                       markerAtom2 = 'o', 
                                       markerLine2 = '-', 
                                       stringPattern = '([\S]+)([^\w]+)([\w]+)\.([\w]+)'):
    '''Input, subChainPairValues, must have exactly the structure as 
    output from readResultsCcode_closedLoops with multiInv_b set to 1.
    
    rowsFig, colsFig: parameters allowing to get the plots placed as subplots in an overall
    plot as subinput an existing figure (subplot structure given by rowsFig rows and colsFig cols).
    yLabelLeft_b: show only y-labels at left-most plots (1).'''
    
    #read input in proper structure:
    subChainPairValuesAll, subChainPairValuesMin, subChainPairValuesMax = subChainPairValues
    
    #if provided, only consider the invariants in the invList and else all:
    if invList:
        keysConsidered = invList
    else:
        keysConsidered = subChainPairValuesAll.keys()

    #plotting is done filling a set of figures with nRows*nCols sufficient to 
    #cover all the invariants in the input invList. The sub-plots are filled up
    #according to a fixed ordering of the invariants:

    #sorted list of the invariants' names:
    sortedInvariants = [
			#/*order 1*/
			"I12", 
			"Ia12",
			#/*order 2:*/
			"I1234", 
			"I1324", 
			"I1423",
			#/*abs value versions:*/
			#/*12*/
			"Ia1234",
			"I12a34",
			"Ia12a34",
			#/*13*/
			"Ia1324",
			"I13a24",
			"Ia13a24",
			#/*14*/
			"Ia1423",
			"I14a23",
			"Ia14a23",
			#/*order 3*/
			"I123456",
			"I123546",
			"I123645",
			#/*13*/
			"I132456",
			"I132546",
			"I132645",
			#/*14*/
			"I142356",
			"I142536",
			"I142635",
			#/*15*/
			"I152346",
			"I152436",
			"I152634",
			#/*16*/
			"I162345",
			"I162435",
			"I162534"]
   
    print keysConsidered

    n = 0 #counts the number of sub-plots in current figure
    for invKey in sortedInvariants:   
        if keysConsidered.count(invKey) > 0:

            #(re)convert data and call the plot fct:
            subChainPairValues_invKey = subChainPairValuesAll[invKey], subChainPairValuesMin[invKey], subChainPairValuesMax[invKey]

            if (n >= rowsFig*colsFig or n < 1): #new fig
                #create outer figure:    
                bigFig = plt.figure() #, axes = plt.subplots(nrows = rowsFig, ncols = colsFig)
                #gen. an appropriate title:
                titleDistr = 'Distr. of extreme inv. values, unrestricted search' #, sub-chain lgth = ' +  str(subChainLgth)
                plt.suptitle(titleDistr, fontsize = 'large')
                n = 0
                subPltCnt = 0 
                subPltRowCnt = 0 #subplot row count
                subPltColCnt = 0 #subplot col count

#            if (m%15 == 13 or m%15 == 14 or m%15 == 0):
#                plt.xlabel('value', fontsize = 'small')
#            if m%3 == 1:
#                plt.ylabel('count', fontsize = 'small')
#            plt.legend(loc = 'upper right', fontsize = 'x-small')
#            plt.xticks(fontsize = 'xx-small')
#            plt.yticks(fontsize = 'xx-small')
#            plt.tight_layout()
    
    

    
            if rowsFig == 1:
                axisIdx = subPltColCnt
            elif colsFig == 1:
                axisIdx = subPltRowCnt
            else:
                axisIdx = subPltRowCnt, subPltColCnt
    
            print "Axis index: %d,%d" % (subPltRowCnt, subPltColCnt)
            
            n +=1  
            axis = bigFig.add_subplot(rowsFig,colsFig,n) 
    
            #only want y-labels on the left-most plots:
            if subPltColCnt == 0:
                yLabelLeft_b = 1
            else:
                yLabelLeft_b = 0
    
            histoSubChainPairs_minMax(subChainPairValues = subChainPairValues_invKey,
                       figIn = bigFig, 
                       axisIn = axis, #axes[axisIdx],
                       yLabel_b = yLabelLeft_b,
                       logCount_b = logCount_b,
                       color = color,
                       colorMin = colorMin,
                       colorMax = colorMax,
                       alphaMin = alphaMin,
                       alphaMax = alphaMax,
                       invName = invKey,
                       titleDistr = '', #titleDistr,
                       lowerPercentile = lowerPercentile,
                       upperPercentile = upperPercentile, 
                       subChainLgth = subChainLgth,
                       bins = bins, 
                       showTopNr = 0, 
                       cutOff = cutOff,  
                       chainsFile = chainsFile, 
                       addSubplotToFirstPlot_b = 0,
                       addSubPlotNr = addSubPlotNr,
                       distrPlotsInSubPlots_b = distrPlotsInSubPlots_b, 
                       few_xticks_b = few_xticks_b,
                       ticksSize = ticksSize, 
                       plotExs_b = plotExs_b, 
                       plotTop10_b = plotTop10_b, 
                       plotTop4_b = plotTop4_b, 
                       interval = interval, 
                       plotExsInInt_b = plotExsInInt_b, 
                       elevAzimList =  elevAzimList, 
                       color1 = color1, 
                       markerAtom1 = markerAtom1, 
                       markerLine1 = markerLine1, 
                       color2 = color2, 
                       markerAtom2 = markerAtom2, 
                       markerLine2 = markerLine2, 
                       stringPattern = stringPattern) 
                      
            subPltCnt += 1                   
            subPltColCnt = (subPltCnt)%colsFig
            if subPltColCnt == 0:
                subPltRowCnt += 1




#Original multi-call to histoSubChainPairs_minMax, for use when input contains sub-chain pair
#values for several invariants:
def orig_multiInv_histoSubChainPairs_minMax(subChainPairValues,
                                       rowsFig = 1,
                                       colsFig = 1,
                                       invList = '',
                                       yLabelLeft_b = 1,
                                       logCount_b = 0,
                                       color = 'lightgrey',
                                       colorMin = 'lightgrey',
                                       colorMax = 'grey',
                                       titleDistr = '',
                                       lowerPercentile = 0.01,
                                       upperPercentile = 0.99,
                                       subChainLgth = 30,
                                       bins = 10, 
                                       cutOff = 11.0,  
                                       chainsFile = '', 
                                       addSubPlotNr = 122,
                                       distrPlotsInSubPlots_b = 0,
                                       few_xticks_b = 0,
                                       plotExs_b = 0, 
                                       plotTop10_b = 0,
                                       plotTop4_b = 0,
                                       interval = [2,3], 
                                       plotExsInInt_b = 0, 
                                       elevAzimList =  [[15,-90], [15,0], [90,0]], 
                                       color1 = 'red', 
                                       markerAtom1 = 'o', 
                                       markerLine1 = '-', 
                                       color2 = 'blue', 
                                       markerAtom2 = 'o', 
                                       markerLine2 = '-', 
                                       stringPattern = '([\S]+)([^\w]+)([\w]+)\.([\w]+)'):
    '''Input, subChainPairValues, must have exactly the structure as 
    output from readResultsCcode_closedLoops with multiInv_b set to 1.
    
    rowsFig, colsFig: parameters allowing to get the plots placed as subplots in an overall
    plot as subinput an existing figure (subplot structure given by rowsFig rows and colsFig cols).
    yLabelLeft_b: show only y-labels at left-most plots (1).'''
    
    #read input in proper structure:
    subChainPairValuesAll, subChainPairValuesMin, subChainPairValuesMax = subChainPairValues
    
    #create outer figure:    
    bigFig, axes = plt.subplots(nrows = rowsFig, ncols = colsFig)
    #gen. an appropriate title:
    titleDistr = 'Distr. of extreme inv. values, unrestricted search, sub-chain lgth = ' +  str(subChainLgth)
    plt.suptitle(titleDistr, fontsize = 'large')

    subPltCnt = 0 
    subPltRowCnt = 0 #subplot row count
    subPltColCnt = 0 #subplot col count
    
    #if provided, only consider the invariants in the invList and else all:
    if invList:
        keysConsidered = invList
    else:
        keysConsidered = subChainPairValuesAll.keys()
    
    #loop over the keys:
    for invKey in keysConsidered:
        
        #(re)convert data and call the plot fct:
        subChainPairValues_invKey = subChainPairValuesAll[invKey], subChainPairValuesMin[invKey], subChainPairValuesMax[invKey]
    
#        #gen. an appropriate title:
#        titleDistr = invKey

        if rowsFig == 1:
            axisIdx = subPltColCnt
        elif colsFig == 1:
            axisIdx = subPltRowCnt
        else:
            axisIdx = subPltRowCnt, subPltColCnt

        print "Axis index: %d,%d" % (subPltRowCnt, subPltColCnt)

        #only want y-labels on the left-most plots:
        if subPltColCnt == 0:
            yLabelLeft_b = 1
        else:
            yLabelLeft_b = 0

        histoSubChainPairs_minMax(subChainPairValues = subChainPairValues_invKey,
                   figIn = bigFig, 
                   axisIn = axes[axisIdx],
                   yLabel_b = yLabelLeft_b,
                   logCount_b = logCount_b,
                   color = color,
                   colorMin = colorMin,
                   colorMax = colorMax,
                   invName = invKey,
                   titleDistr = '', #titleDistr,
                   lowerPercentile = lowerPercentile,
                   upperPercentile = upperPercentile, 
                   subChainLgth = subChainLgth,
                   bins = bins, 
                   showTopNr = 0, 
                   cutOff = cutOff,  
                   chainsFile = chainsFile, 
                   addSubplotToFirstPlot_b = 0,
                   addSubPlotNr = addSubPlotNr,
                   distrPlotsInSubPlots_b = distrPlotsInSubPlots_b, 
                   few_xticks_b = few_xticks_b,
                   plotExs_b = plotExs_b, 
                   plotTop10_b = plotTop10_b, 
                   plotTop4_b = plotTop4_b, 
                   interval = interval, 
                   plotExsInInt_b = plotExsInInt_b, 
                   elevAzimList =  elevAzimList, 
                   color1 = color1, 
                   markerAtom1 = markerAtom1, 
                   markerLine1 = markerLine1, 
                   color2 = color2, 
                   markerAtom2 = markerAtom2, 
                   markerLine2 = markerLine2, 
                   stringPattern = stringPattern) 
                  
        subPltCnt += 1                   
        subPltColCnt = (subPltCnt)%colsFig
        if subPltColCnt == 0:
            subPltRowCnt += 1    


#utilities to create tube around piecewise linear curve (plot of chain):
#3x3 matrix repr rotation of angle theta around z-axis:
def T_theta(theta):
    '''returns: 3x3 matrix repr rotation of angle theta around z-axis:'''
    
    T0 = np.array([cos(theta), -sin(theta), 0])
    T1 = np.array([sin(theta), cos(theta), 0])
    T2 = np.array([0,0, 1])

    T = np.matrix([T0,T1,T2])    
    
    return T

def T_thetaCosSin(cosTheta, sinTheta):
    '''returns: 3x3 matrix repr rotation of angle theta around z-axis:'''
    
    T0 = np.array([cosTheta, -sinTheta, 0])
    T1 = np.array([sinTheta, cosTheta, 0])
    T2 = np.array([0,0, 1])

    T = np.matrix([T0,T1,T2])    
    
    return T
      

#3x3 matrix repr rotation of angle phi around x-axis:
def T_phi(phi):
    '''returns: 3x3 matrix repr rotation of angle phi around x-axis:'''
    T0 = np.array([1,0, 0])    
    T1 = np.array([0, cos(phi), -sin(phi)])
    T2 = np.array([0, sin(phi), cos(phi)])
    
    T = np.matrix([T0,T1,T2])
    
    #print T
    
    return T
    

def T_phiCosSin(cosPhi,sinPhi):
    '''returns: 3x3 matrix repr rotation of angle phi around x-axis:'''
    T0 = np.array([1,0, 0])    
    T1 = np.array([0, cosPhi, -sinPhi])
    T2 = np.array([0, sinPhi, cosPhi])
    
    T = np.matrix([T0,T1,T2])    
    
    return T  
    
    
def linTransform(T,x):
    return dot(T,x)
    
#to generate n coordinates of circle with radius r around center c and
#tilted by rotataion of angles (theta, phi) around z resp x-axis:
def tiltedCircle(c, r, theta, phi, n):
    '''Returns n coordinates of circle with radius r around center c and
    tilted by rotataion of angles (theta, phi) around z resp x-axis.
    Input:
    c: center must be np.array of shape 3'''
    
    out = np.zeros(shape = (n+1,3))
    
    #generate n pts on unit circle in x-y-plane
    for i in range(n+1): #first and last points are id: this is so as to get a full circle when using plot
        #pt on unit circle at angle 2pi/n*i:
        out[i] = np.array([cos(2*np.pi/n*i),sin(2*np.pi/n*i),0])
        #scale radius to r
        out[i] = r*out[i]
        #tilt it
        out[i] = np.dot(T_phi(phi),out[i])
        out[i] = np.dot(T_theta(theta),out[i])
        #translate to center c:
        out[i] = out[i] + c
    
    return out
    
    
def tiltedCircleCosSin(c, r, cosTheta, sinTheta, cosPhi, sinPhi, n):
    '''Returns n coordinates of circle with radius r around center c and
    tilted by rotataion of angles (theta, phi) around z resp x-axis.
    Input:
    c: center must be np.array of shape 3
    cosTheta, sinTheta: instead of angle Theta this version uses the pair
    cos(theta), sin(theta) instead.'''
    
    out = np.zeros(shape = (n+1,3))
    
    #generate n pts on unit circle in x-y-plane
    for i in range(n+1): #first and last points are id: this is so as to get a full circle when using plot
        #pt on unit circle at angle 2pi/n*i:
        out[i] = np.array([cos(2*np.pi/n*i),sin(2*np.pi/n*i),0])
        #scale radius to r
        out[i] = r*out[i]
        #tilt it
        out[i] = np.dot(T_phiCosSin(cosPhi, sinPhi),out[i])
        out[i] = np.dot(T_thetaCosSin(cosTheta, sinTheta),out[i])       
        #translate to center c:
        out[i] = out[i] + c
    
    return out
    
        
#For generating small Pymol scripts for a set of examples:
def pymolClosedExamples(outputFilePath,
                           examples, 
                           color = 'lightorange',
                           color1 = 'red', 
                           color2 = 'blue', 
                           background = 'white'):
    '''Generate a small Pymol script to generate in Pymol 
    the example with the two sub strings highlighted in the set colors'''

    with open(outputFilePath, 'w') as outfile:   
    
        for ex in examples:
            val = ex[0]
            pdb_filename = ex[1]
            chainId = ex[2]
    #        print pdb_filename
            segment1 = ex[3]
            segment2 = ex[4]
        
            loadPDB = 'PyMOL> load ' + pdb_filename + '\n'
            cmd.load(pdb_filename)
            
#            outfile.write(script)
        


#the follwing plot function is "stolen" from the LinuxPython exam (Oct 2013) though slightly
#modified for our purpose: to generate 3d-pictures of a structure with the ability to color
#two given segments so as to highlight them. The targeted usage is to color a closed loop
#and a sub chain/closed loop, forming a potential link or poke example:

from mpl_toolkits.mplot3d import Axes3D #needed, see matplotlib.org

def plot_3D_2segmentsHighLight(polyCaChain, 
                               segment1, 
                               segment2, 
                               color1 = 'red', 
                               markerAtom1 = 'o', 
                               markerLine1 = '-', 
                               markerAtomStartPt1 = 'o',
                               startPtSize1 = 75,
                               color2 = 'blue', 
                               markerAtom2 = 'o', 
                               markerLine2 = '-', 
                               markerAtomStartPt2 = 'o',
                               startPtSize2 = 75,
                               elevAzimList = [[0, 36]], 
                               title = '', 
                               outFile = '', 
                               optionOut = 'screen',
                               tubePlot_b = 0,
                               tubeRadius = 0.75,
                               nrOfCircles = 10,
                               nrOfPtsOnCircles = 10):
    '''Creates one ore more 3D plot for the supplied pdb file/protein and colors the two given segments 
    in the set colors. Called by plotClosedLoopExamples and, in turn, by the functions histoSubChainPairs_minMax
    and histoClosedLoops_minMax.
    
    Input:
       polyCaChain: list of [segment start, segment end]; as output by readResultsCcode_chains (see there for more) 
       segment1: [start index, end index] of 1st segment 
       segment2: [start index, end index] of 2nd segment 
       color1 ('red'): high-light color for 1st segment
       markerAtom1 ('o'): marker of C-alphas in 1st segment
       markerLine1 ('-'): marker of C-alphas not in 1st or 2nd segment
       markerAtomStartPt1 ('o'): marker of first C-alpha in 1st segment
       startPtSize1 (75): size of marker of first C-alpha in 1st segment
       color2 ('blue'): as for 1st segment (color1)
       markerAtom2 ('o'): as for 1st segment
       markerLine2 ('-'): as for 1st segment
       markerAtomStartPt2 ('o'): as for 1st segment
       startPtSize2 (75): as for 1st segment
       elevAzimList ([[0, 36]]): list of angles for which to plot; the length of the list should be 1, 2 or 3.
       There will be one subplot for each angle
       title (''): (first part of) title on plot 
       outFile (''): if not blank ('') plot will be written out to this pdf
       optionOut ('screen'): if set to 'pdf' plot will be written out [outFile].pdf
       tubePlot_b (0): for plotting the C-alpha trace (curve) as a tube 
       tubeRadius (0.75): radius of the tube
       nrOfCircles (10): skeleton for tube (nr of circles around each C-alpha-to-C-alpha segment)
       nrOfPtsOnCircles (10): skeleton for tube (approx of circle around the curve)
       '''          
        
#    try:
#        #CaChain, polyCaChain = getPolygonalChain(PDBfilename = pdb_filename, outNumpy_b = 0)
#        CaChain, polyCaChain = readResultsCcode_chains(inputFileName = pdb_filename)
#    except KeyError:
#        print 'Load of this PDB-file: %s failed' % pdb_filename
#        return    
    
    fig = plt.figure()
    L = len(elevAzimList)
    cnt = 1
    i_e = 0
    for i_e in range(L): #elevAzim in elevAzimList:
        elevAzim = elevAzimList[i_e]
#        print "Subplot nr:%d " % 100 + L*10  + cnt
        ax = fig.add_subplot(100 + L*10  + cnt ,projection ='3d') #adds 3rd axis
    
        seg0 = []
        seg1 = []
        seg2 = []
        startPtSeg1 = []
        startPtSeg2 = []
        for i in range(len(polyCaChain)):
            pt0 = polyCaChain[i][0]
            x1 = [pt0[0], pt0[1],pt0[2]]
            pt1 = polyCaChain[i][1]        
            x2 = [pt1[0], pt1[1],pt1[2]]
            seg0.append(x1)
            if (segment1[0]<= i and i <= segment1[1]): #use color1 if the C-alpha (residue) is in segment1
    #            print "at i:%d color is : %s" % (i, color1)
                seg1.append(x1)
                #we have appended the starting pt of each segment; for the last segment append alos the end point:
                if i == segment1[1]:
                    seg1.append(x2)
    #            ax.scatter(x1, 'o', c=color1)
                #we want to indicate the starting point of the segment, so store it:
                if i == segment1[0]:
                    startPtSeg1 = x1
            #Obs: we use an "if" next rather than an "elif": only diff is that it can happen that 
            #segment1[1] = segment2[0] in which case the elif will imply that there'll be no
            #startPtSeg2! This does not happen when using an if instead:
            if (segment2[0]<= i and i <= segment2[1]): #use color2 if the C-alpha (residue) is in segment2
    #            ax.scatter(x1, 'o', c=color2)
    #            print "at i:%d color is : %s" % (i,color2)
                seg2.append(x1)
                #we have appended the starting pt of each segment; for the last segment append alos the end point:
                if i == segment2[1]:
                    seg2.append(x2)
                #we want to indicate the starting point of the segment, so store it:
                if i == segment2[0]:
                    startPtSeg2 = x1
        segs = [seg0,seg1,seg2]
        for i in range(3):
            thisPart = np.array(segs[i])   
            x = thisPart[:,0]
            y = thisPart[:,1]
            z = thisPart[:,2]
            #first color all black and then overwrite in segment1 & 2:
            if i == 0:
                col = "black"
                ax.plot(x,y,z, '-', c=col, linewidth = 2)    
    #            ax.scatter(x,y,z, ".", c = col)  
            elif i ==1:
                mrkLine = markerLine1
                mrkAtom = markerAtom1
                col = color1
                ax.scatter(x,y,z, mrkAtom, c = col)
                #indicate the start point:
                x0 = startPtSeg1[0]
                x1 = startPtSeg1[1]
                x2 = startPtSeg1[2]
#                print "start point 1: %f, %f, %f" % (x0,x1,x2)
                ax.scatter(x0, x1, x2, markerAtomStartPt1, c = col, s = startPtSize1) 
                #plot the set of line segments, or, if desired, "tubing them" 
                if tubePlot_b != 1:            
                    ax.plot(x,y,z, mrkLine, c=col, linewidth = 10, alpha = 0.5)    
                else: #generate a tube around the segments
                    for i in range(len(seg1)-1):
                        #generate a series of nrOfCircles circles placed around each segment, in planes perpendicular to the direction of the segment
                        #compute angles for tilting the circles to be placed around the line segments:
                        v = np.array(seg1[i+1]) - np.array(seg1[i]) #vector along line segment
                        theta = np.arctan2(v[1], v[0])
                        phi = np.arctan2(v[2], np.sqrt(v[0]*v[0] + v[1]*v[1]) )
#                        X = np.zeros(shape = (nrOfCircles,nrOfPtsOnCircles +1))
#                        Y = np.zeros(shape = (nrOfCircles,nrOfPtsOnCircles +1))
#                        Z = np.zeros(shape = (nrOfCircles,nrOfPtsOnCircles +1))

                        for j in range(nrOfCircles):
                            c = float(j)/nrOfCircles*np.array(seg1[i+1]) + (1- float(j)/nrOfCircles)*np.array(seg1[i])
                            #pts = tiltedCircleCosSin(np.array(c), 1, cosTheta, sinTheta, cosPhi, sinPhi, 10)
                            pts = tiltedCircle(np.array(c), tubeRadius, theta - np.pi*0.5, phi - np.pi*0.5 , nrOfPtsOnCircles)
                            x_pts = pts.T[0]
                            y_pts = pts.T[1]
                            z_pts = pts.T[2]
                            #ax.scatter(x_pts, y_pts, z_pts, "*", c = col)
                            ax.plot(x_pts, y_pts, z_pts, "-", c = col, alpha = 0.5)
                            #keep record of the first and last circle for next step
                            if j == 0:
                                x_pts_start = pts.T[0]
                                y_pts_start = pts.T[1]
                                z_pts_start = pts.T[2]
                            if j == nrOfCircles -1:
                                x_pts_end = pts.T[0]
                                y_pts_end = pts.T[1]
                                z_pts_end = pts.T[2]
#                            X[j] = x_pts
#                            Y[j] = y_pts
#                            Z[j] = z_pts
#                        surf = ax.plot_surface(X.T, Y.T,Z.T,  color = col, facecolors = col, alpha = 0.5)
                        #generate a series of straight line segments connecting the circles:
                        for j in range(nrOfPtsOnCircles-1):
                            x = [x_pts_start[j], x_pts_end[j]]
                            y = [y_pts_start[j], y_pts_end[j]]
                            z = [z_pts_start[j], z_pts_end[j]]
                            #ax.plot(x,y,z, mrkLine, c=col, linewidth = 2, alpha = 0.5)  
                            ax.plot_trisurf([x_pts_start[j], x_pts_end[j], x_pts_start[j+1], x_pts_end[j+1]], [y_pts_start[j], y_pts_end[j], y_pts_start[j+1], y_pts_end[j+1]], [z_pts_start[j], z_pts_end[j], z_pts_start[j+1], z_pts_end[j+1]] , color = 'green', alpha = 0.5)
                            #surf = ax.plot_surface([x_pts_start[j], x_pts_end[j], x_pts_start[j+1], x_pts_end[j+1]], [y_pts_start[j], y_pts_end[j], y_pts_start[j+1], y_pts_end[j+1]], [z_pts_start[j], z_pts_end[j], z_pts_start[j+1], z_pts_end[j+1]] , color = col, facecolors = col, alpha = 0.5)

                
            elif i == 2:
                mrkLine = markerLine2
                mrkAtom = markerAtom2
                col = color2
                ax.scatter(x,y,z, mrkAtom, c = col)            
                #indicate the start point:
                x0 = startPtSeg2[0]
                x1 = startPtSeg2[1]
                x2 = startPtSeg2[2]
#                print "start point 2: %f, %f, %f" % (x0,x1,x2)
                ax.scatter(x0, x1, x2, markerAtomStartPt2, c = col, s = startPtSize2) 
                #plot the set of line segments, or, if desired, "tubing them" 
                if tubePlot_b != 1:            
                    ax.plot(x,y,z, mrkLine, c=col, linewidth = 10, alpha = 0.5)
                else: #generate a tube around the segments
                    for i in range(len(seg2)-1):
                        #generate a series of nrOfCircles circles placed around each segment, in planes perpendicular to the direction of the segment
                        #compute angles for tilting the circles to be placed around the line segments:
                        v = np.array(seg2[i+1]) - np.array(seg2[i]) #vector along line segment
                        theta = np.arctan2(v[1], v[0])
                        phi = np.arctan2(v[2], np.sqrt(v[0]*v[0] + v[1]*v[1]) )
                        for j in range(nrOfCircles):
                            c = float(j)/nrOfCircles*np.array(seg2[i+1]) + (1- float(j)/nrOfCircles)*np.array(seg2[i])
                            #pts = tiltedCircleCosSin(np.array(c), 1, cosTheta, sinTheta, cosPhi, sinPhi, 10)
                            pts = tiltedCircle(np.array(c), tubeRadius, theta - np.pi*0.5, phi - np.pi*0.5 , nrOfPtsOnCircles)
                            x_pts = pts.T[0]
                            y_pts = pts.T[1]
                            z_pts = pts.T[2]
                            #ax.scatter(x_pts, y_pts, z_pts, "*", c = col)
                            ax.plot(x_pts, y_pts, z_pts, "-", c = col, alpha = 0.3)
                            #keep record of the first and last circle for next step
                            if j == 0:
                                x_pts_start = pts.T[0]
                                y_pts_start = pts.T[1]
                                z_pts_start = pts.T[2]
                            if j == nrOfCircles -1:
                                x_pts_end = pts.T[0]
                                y_pts_end = pts.T[1]
                                z_pts_end = pts.T[2]       
                        #generate a series of straight line segments connecting the circles:
                        for j in range(nrOfPtsOnCircles):
                            x = [x_pts_start[j], x_pts_end[j]]
                            y = [y_pts_start[j], y_pts_end[j]]
                            z = [z_pts_start[j], z_pts_end[j]]
                            #ax.plot(x,y,z, mrkLine, c=col, linewidth = 2, alpha = 0.3)  
                            ax.plot_trisurf([x_pts_start[j], x_pts_end[j], x_pts_start[j+1], x_pts_end[j+1]], [y_pts_start[j], y_pts_end[j], y_pts_start[j+1], y_pts_end[j+1]], [z_pts_start[j], z_pts_end[j], z_pts_start[j+1], z_pts_end[j+1]] , color = col, alpha = 0.5)

                
#        print elevAzim
        elev = elevAzim[0]
        azim = elevAzim[1]          
        ax.view_init(elev=elev,azim =azim)
            
        cnt +=1
        
#    if not(title):
#        #find the PDB name (avoiding the file path)
#        pattern = re.compile('([\S]+)([^\w]+)([\w]+).([\S]+)')
#        m = pattern.match(pdb_filename)
#        title = m.group(3)
    plt.suptitle(title + ', '+ str(segment1) + ', '+ str(segment2), fontsize = "x-large")
    
    
    if optionOut == 'pdf':
        plt.savefig(outFile)
        
    return ax
    


def plotClosedLoopExamples(examples, 
                           invName = 'writhe',
                           chainsDict = {}, 
                           elevAzimList = [[0,36]], 
                           color1 = 'black', 
                           color2 = 'grey',
                           markerAtom1 = 'o', 
                           markerLine1 = '-', 
                           markerAtom2 = 'o', 
                           markerLine2 = '-', 
                           stringPattern = '([\S]+)([^\w]+)([\w]+)\.([\w]+)'):
    '''Plot (3d) a set of closed loop examples such as those returned when loading the C-results
    using readResultsCcode_closedLoops (e.g. look at a subset such as linkValues[:10]); calls 
    the fct plot_3D_2segmentsHighLight to plot each example. In each example the segments are
    high-lighted by distinct colors/markers.
    
    The function is itself called by the functions by the functions histoClosedLoops_minMax and 
    histoSubChainPairs_minMax.    
    
    Input:
    
       examples: as output from readResultsCcode_closedLoops (or part of it); see histoClosedLoops_minMax and histoSubChainPairs_minMax.
       chainsDict: as output from readResultsCcode_chains; see histoClosedLoops_minMax and histoSubChainPairs_minMax.
    
       invName ('writhe'): see histoSubChainPairs_minMax
       elevAzimList ([[0, 36]]): see plot_3D_2segmentsHighLight
       color1 ('red'): see plot_3D_2segmentsHighLight
       color2 ('blue'): see plot_3D_2segmentsHighLight
       markerAtom1 ('o'): see plot_3D_2segmentsHighLight    
       markerLine1 ('-'): see plot_3D_2segmentsHighLight    
       markerAtom2 ('o'): see plot_3D_2segmentsHighLight    
       markerLine2 ('-'): see plot_3D_2segmentsHighLight    
       stringPattern ('([\S]+)([^\w]+)([\w]+)\.([\w]+)'): see histoClosedLoops_minMax
    '''
    for ex in examples:
        val = ex[0]
        pdb_filename = ex[1]
        chainId = ex[2]
#        print pdb_filename
        segment1 = ex[3]
        segment2 = ex[4]
        #gen plot title
        #first find the PDB name (avoiding the file path)
        pattern = re.compile(stringPattern)
#        print pdb_filename
        m = pattern.match(pdb_filename)
        title = m.group(3)
        #now add the writhe value too:
        title = title + ' chain: ' + chainId + ', ' + invName + ': ' + str(val)
        #fetch the segments chain:
        try: 
            polyCaChain = chainsDict[pdb_filename][chainId][1]
        except KeyError:
            print 'Load of this PDB-file: %s chain: %s failed' % (pdb_filename, chainId)
            continue        
#        print "Now at structure: %s chain: %s" %  (pdb_filename, chainId)
        #call the plot fct:
        ax = plot_3D_2segmentsHighLight(polyCaChain, segment1, segment2, 
                                   color1 = color1, markerAtom1 = markerAtom1, markerLine1 = markerLine1,
                                   color2 = color2, markerAtom2 = markerAtom2, markerLine2 = markerLine2, 
                                   elevAzimList = elevAzimList, title = title, outFile = '', optionOut = 'screen')

#        return ax




def plot_3D_Quiver_2segmentsHighLight(polyCaChain, 
                               segment1, 
                               segment2, 
                               color1 = 'red', 
                               markerAtom1 = 'o', 
                               markerLine1 = '-', 
                               color2 = 'blue', 
                               markerAtom2 = 'o', 
                               markerLine2 = '-', 
                               elevAzimList = [[0, 36]], 
                               title = '', 
                               outFile = '', 
                               optionOut = 'screen'):
    '''Creates one ore more 3D plot for the supplied pdb file/protein and colors the two given segments 
    in the set colors. 
    Input:
    the elevAzimList sets the angles for which to plot; the length of the list should be 1, 2 or 3.'''     
        
#    try:
#        #CaChain, polyCaChain = getPolygonalChain(PDBfilename = pdb_filename, outNumpy_b = 0)
#        CaChain, polyCaChain = readResultsCcode_chains(inputFileName = pdb_filename)
#    except KeyError:
#        print 'Load of this PDB-file: %s failed' % pdb_filename
#        return    
    
    fig = plt.figure()
    L = len(elevAzimList)
    cnt = 0
    for elevAzim in elevAzimList:
#        print "Subplot nr:%d " % 100 + L*10  + cnt
        ax = fig.add_subplot(100 + L*10  + cnt ,projection ='3d') #adds 3rd axis
    
        seg0 = []
        seg1 = []
        seg2 = []
        
        minLgth = 1e16
        for i in range(len(polyCaChain)):
            pt0 = polyCaChain[i][0]
            x1 = [pt0[0], pt0[1],pt0[2]]
            pt1 = polyCaChain[i][1]        
            x2 = [pt1[0], pt1[1],pt1[2]]
            seg0.append(x1)
            v = np.array(pt1) - np.array(pt0)
            l = np.linalg.norm(v)
            if l < minLgth:
                minLgth = l
            if (segment1[0]<= i and i <= segment1[1]): #use color1 if the C-alpha (residue) is in segment1
    #            print "at i:%d color is : %s" % (i, color1)
                seg1.append(x1)
                #we have appended the starting pt of each segment; for the last segment append alos the end point:
                if i == segment1[1]:
                    seg1.append(x2)
    #            ax.scatter(x1, 'o', c=color1)
            elif (segment2[0]<= i and i <= segment2[1]): #use color2 if the C-alpha (residue) is in segment2
    #            ax.scatter(x1, 'o', c=color2)
    #            print "at i:%d color is : %s" % (i,color2)
                seg2.append(x1)
                #we have appended the starting pt of each segment; for the last segment append alos the end point:
                if i == segment2[1]:
                    seg2.append(x2)
        print "min lgth:%d" % minLgth
        segs = [seg0,seg1,seg2]
        for i in range(3):
            thisPart = np.array(segs[i])   
            x = thisPart[:,0]
            y = thisPart[:,1]
            z = thisPart[:,2]
            u = x[1:] - x[:(len(x)-1)]
            v = y[1:] - y[:(len(y)-1)]
            w = z[1:] - z[:(len(z)-1)]
            x = x[1:]
            y = y[1:]
            z = z[1:]
            #first color all black and then overwrite in segment1 & 2:
            if i == 0:
                col = "black"
                ax.plot(x,y,z, '-', c=col, linewidth = 0.5)    
    #            ax.scatter(x,y,z, ".", c = col)  
            elif i ==1:
                mrkLine = markerLine1
                mrkAtom = markerAtom1
                col = color1
                cols = [(0,0,1) for i in range(len(x))]
                ax.quiver(x,y,z,u,v,w, colors =cols, length = minLgth, linewidth = 2)  
#                ax.plot(x,y,z, mrkLine, c=col, linewidth = 2)    
#                ax.scatter(x,y,z, mrkAtom, c = col)  
            elif i == 2:
                mrkLine = markerLine2
                mrkAtom = markerAtom2
                col = color2
                cols = [(1,0,0) for i in range(len(x))]
                ax.quiver(x,y,z,u,v,w, colors =cols, length = minLgth, linewidth = 2) #, c=col, linewidth = 2) 
#                ax.plot(x,y,z, mrkLine, c=col, linewidth = 4, alpha = 1)    
#                ax.scatter(x,y,z, mrkAtom, c = col)            
                
#        print elevAzim
        elev = elevAzim[0]
        azim = elevAzim[1]          
        ax.view_init(elev=elev,azim =azim)
            
        cnt +=1
        
#    if not(title):
#        #find the PDB name (avoiding the file path)
#        pattern = re.compile('([\S]+)([^\w]+)([\w]+).([\S]+)')
#        m = pattern.match(pdb_filename)
#        title = m.group(3)
    plt.title(title + ', '+ str(segment1) + ', '+ str(segment2), fontsize = "x-large")
    
    
    if optionOut == 'pdf':
        plt.savefig(outFile)


            
def PCA(inputFileName, incl_abs_b =0):
    '''1. Computes and compares the variance in the supplied invariant values at order 2
    and order 3.
    2. Performs a principal component analysis of the values at order 3.
    Input: path to input file as for the fct readResultsCcode.
    Output: ...TBD
    '''
    
    #get the values of the invariants:
    resDict = readResultsCcode(inputFileName = inputFileName, all_b = 0, pert_b = 0, pertTime_b = 0)


    order1 = ["I12", "Ia12"]
    order2 = ["I1234", 
			"I1324", 
			"I1423",
			#/*abs value versions:*/
			#/*12*/
			"Ia1234",
			"I12a34",
			"Ia12a34",
			#/*13*/
			"Ia1324",
			"I13a24",
			"Ia13a24",
			#/*14*/
			"Ia1423",
			"I14a23",
			"Ia14a23"]
    order3 = ["I123456",
			"I123546",
			"I123645",
			#/*13*/
			"I132456",
			"I132546",
			"I132645",
			#/*14*/
			"I142356",
			"I142536",
			"I142635",
			#/*15*/
			"I152346",
			"I152436",
			"I152634",
			#/*16*/
			"I162345",
			"I162435",
			"I162534"] 
    
    if incl_abs_b ==1:

        readOrder_3 = ["str length",
                       "I12", 
    			"Ia12",
    			#/*order 2:*/
    			"I1234", 
    			"I1324", 
    			"I1423",
    			#/*abs value versions:*/
    			#/*12*/
    			"Ia1234",
    			"I12a34",
    			"Ia12a34",
    			#/*13*/
    			"Ia1324",
    			"I13a24",
    			"Ia13a24",
    			#/*14*/
    			"Ia1423",
    			"I14a23",
    			"Ia14a23",
    			#/*order 3*/
    			"I123456",
    			"I123546",
    			"I123645",
    			#/*13*/
    			"I132456",
    			"I132546",
    			"I132645",
    			#/*14*/
    			"I142356",
    			"I142536",
    			"I142635",
    			#/*15*/
    			"I152346",
    			"I152436",
    			"I152634",
    			#/*16*/
    			"I162345",
    			"I162435",
    			"I162534"]
    
    
        readOrder_2 = ["str length",
                "I12", 
    			"Ia12",
    			#/*order 2:*/
    			"I1234", 
    			"I1324", 
    			"I1423",
    			#/*abs value versions:*/
    			#/*12*/
    			"Ia1234",
    			"I12a34",
    			"Ia12a34",
    			#/*13*/
    			"Ia1324",
    			"I13a24",
    			"Ia13a24",
    			#/*14*/
    			"Ia1423",
    			"I14a23",
    			"Ia14a23"]
       
    else:
        readOrder_3 = ["str length",
                   "I12", 
			#/*order 2:*/
			"I1234", 
			"I1324", 
			"I1423",
			#/*order 3*/
			"I123456",
			"I123546",
			"I123645",
			#/*13*/
			"I132456",
			"I132546",
			"I132645",
			#/*14*/
			"I142356",
			"I142536",
			"I142635",
			#/*15*/
			"I152346",
			"I152436",
			"I152634",
			#/*16*/
			"I162345",
			"I162435",
			"I162534"]


        readOrder_2 = ["str length", 
    			"I12",
    			#/*order 2:*/
    			"I1234", 
    			"I1324", 
    			"I1423"]



    #Compute the covariance matrix at order 3 and at order 2:
    N2 = len(readOrder_2) -1       
    N3 = len(readOrder_3) -1
    print "Number of order 2 invariants: %d order 3: %d" % (N2, N3)
    #create matrix in which the rows are the set of values for one of the invariants:
    Imatrix_2 = np.zeros(shape = (N2,8000))
    Imatrix_3 = np.zeros(shape = (N3,8000))
    cntFile = 0
    for fileKey in resDict.keys():
        cntName = 0
        for name in readOrder_3: #resDict[fileKey].keys():
            if name == "str length":
                #print resDict[fileKey][name]
                try:
                    L = int(resDict[fileKey][name])
                except ValueError:
                    continue
            if (name != "cpuTime" and name != "aggrTime" and name != "order" and name != "str length"):
                vol = 1
                if order1.count(name) > 0:
                    vol = L-1
                elif order2.count(name) > 0:
                    vol = 0.5*((L-1)*(L-1) - L)
                else: #order 3
                    vol = 0.5*((L-1)*(L-1)*(L-1) - 2*L*L)
#                    for i in range(L-1):
#                        vol += 0.5*((L-1)*(L-1) - L)
#                vol = 1
                try:
                    Imatrix_3[cntName, cntFile] = float(resDict[fileKey][name])/vol
#                    if cntFile == 0:
#                        print fileKey, name, Imatrix[cntName, cntFile]                        
                except ValueError:
                    continue
                cntName +=1 
        cntName = 0
        for name in readOrder_2: #resDict[fileKey].keys():
            if name == "str length":
                try:
                    L = int(resDict[fileKey][name])
                except ValueError:
                    continue
            if (name != "cpuTime" and name != "aggrTime" and name != "order" and name != "str length"):
                vol = 1
                if order1.count(name) > 0:
                    vol = L-1
                elif order2.count(name) > 0:
                    vol = 0.5*((L-1)*(L-1) - L)
                else: #order 3
                    vol = 0.5*((L-1)*(L-1)*(L-1) - 2*L*L)
#                    for i in range(L-1):
#                        vol += 0.5*((L-1)*(L-1) - L)
#                vol = 1
                try:
                    Imatrix_2[cntName, cntFile] = float(resDict[fileKey][name])/vol
#                    if cntFile == 0:
#                        print fileKey, name, Imatrix[cntName, cntFile]                        
                except ValueError:
                    continue
                cntName +=1 
#        print cntName
        cntFile +=1
    print "Data loaded for number of files: %d" % cntFile

    #covariance matrix, order 3:
    corr = np.corrcoef(Imatrix_3) 
    cov = np.cov(Imatrix_3)
    
    var3 = np.trace(cov) 
    print "Variance in order 3 values: %f" % var3
    
    #PCA
    eigVals, eigVects = np.linalg.eig(corr)
    fig1 = plt.figure()
    plt.plot(eigVals, label = "Eigenvalues, order 3")
    fig2 = plt.figure()
    for i in range(N3)[:14]:
        plt.plot(eigVects[i], label = i)
    plt.legend()
    
    fig3 = plt.figure()
    for i in range(N3)[14:]:
        plt.plot(eigVects[i], label = i)
    plt.legend()

    v3 = []
    v = 0
    eigVals_cov, eigVects_cov = np.linalg.eig(cov)
    for i in range(N3):
        v +=  eigVals_cov[i]
        v3.append(float(v)/var3)
    fig4 = plt.figure()
    plt.plot(v3, label = "order 3, variance explanation degree")
    plt.legend()

    #covariance matrix, order 2:
    corr = np.corrcoef(Imatrix_2) 
    cov = np.cov(Imatrix_2)
    
    #print corr
    
    var2 = np.trace(cov) 
    print "Variance in order 3 values: %f" % var2
    
    #PCA
    eigVals, eigVects = np.linalg.eig(corr)
    fig5 = plt.figure()
    plt.plot(eigVals, label = "Eigenvalues, order 2")
    fig6 = plt.figure()
    for i in range(N2):
        plt.plot(eigVects[i], label = i)
    plt.legend()

    v2 = []
    v = 0
    eigVals_cov, eigVects_cov = np.linalg.eig(cov)
    for i in range(N2):
        v +=  eigVals_cov[i]
        v2.append(float(v)/var2)
    fig7 = plt.figure()
    plt.plot(v2, label = "order 2, variance explanation degree")
    plt.legend()

    ratio = float(var2)/var3
    print "Variance ratio, order 2 over order 3:%f" % ratio

    return Imatrix_2, Imatrix_3


def scatterTwoDicts(d1, d2):
    x = []
    y = []
    for f1 in d1.keys():
        for s1 in d1[f1].keys():
            for sc1 in d1[f1][s1].keys():
                for sc2 in d1[f1][s1][sc1].keys():
                    x.append(float(d1[f1][s1][sc1][sc2]))
                    y.append(float(d2[f1][s1][sc1][sc2]))
    plt.scatter(x,y)
    
    
    
###############################################################
##    For enabling execution from Pymol 
###############################################################

cmd.extend('readResultsCcode_closedLoops', readResultsCcode_closedLoops)
cmd.extend('histoClosedLoops_minMax',histoClosedLoops_minMax)    