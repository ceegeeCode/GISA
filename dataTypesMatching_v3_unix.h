/* **********************************************************************************************************************************
*
*
* ***************************************  Header for: Data types util for GISA_v* and SAonGISA_v* *******************************************************************************
*
* See main code (GISA_v*) for do and how to cite.
*
* Present version: version _v3, > 1st Jan, 2019.
*
* Author: Christian Grønbæk
*
**********************************************************************************************/



#ifndef DATATYPESMATCHING_V3_UNIX
#define DATATYPESMATCHING_V3_UNIX

//#pragma once 


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
//#include "GISA_v9_unix.h"
#include "mathUtils.h"


/**************************************************
******************* Structs ***********************
***************************************************/

/*Struct holding a list of floats serving the purpose of being an array of floats e.g. tolerances/epsilons alogn with its length*/
typedef struct floatArrayWLength{

	int length; /* to hold the length of the epsRange array*/
	float *floatArray;

} floatArrayWLength;

/*Struct for holding a word*/
typedef struct word{

	char *word;

} word;

typedef struct window{

	int windowNr; /*to hold the window nr of a window*/

	int *segIndices; /*to hold the segment indices of the corr segment*/

} window;


typedef struct windowPair{

	int windowNr_1;
	/*segment indices (bounds of each window considered)*/
	int *segIndices_1;

	int windowNr_2;
	/*segment indices (bounds of each window considered)*/
	int *segIndices_2;
	
	//int *windowNr; /*to hold the window nr's for two windows*/

	//int **segIndices; /*to hold the segment indices of the two corr segments*/

} windowPair;


/*Struct for holding the score and pValue for a query*/
typedef struct queryRawScore{

	char *fileName;

	char *structureName; /*name of structure in PDB-jargon, particular to application, so e.g. domain name in CATH*/
	char *classId; /*classification id in e.g. CATH or SCOP*/
	
	char *chainId;
	int  chainNr;
	int  chainLen;

	int nrOfWindows; 

	//int nrOfWindowPairs; 

	char *windowInfo; /*to hold a list of window nr's with corr segment indices, or, similarly, or window pairs*/ 

	double score;

	double pValue; 

	double optionalValue; //to hold a writhe value or other

} queryRawScore;

/*Struct for storing, for a query, the output of a match on a window*/ 
typedef struct matchCandidate{

	char *fileName;

	char *structureName; /*name of structure in PDB-jargon, particular to application, so e.g. domain name in CATH*/
	char *classId; /*classification id in e.g. CATH or SCOP*/

	int fileNr; /*internal file number for ease of look-ups*/ 
	
	char *chainId;
	int  chainNr;
	int  chainLen;

	int nrOfWindows; 

	struct window window;

	//int windowsNr;
	///*segment indices (bounds of each window considered)*/
	//int *segIndices;

	double distMin;

	double score;

} matchCandidate;


/*Struct for storing, for a query, the output of a match on a window pair*/ 
typedef struct matchCandidate_windowPair{

	char *fileName;

	char *structureName; /*name of structure in PDB-jargon, particular to application, so e.g. domain name in CATH*/
	char *classId; /*classification id in e.g. CATH or SCOP*/

	int fileNr; /*internal file number for ease of look-ups*/ 
	
	char *chainId;
	int  chainNr;
	int  chainLen;

	int nrOfWindows; 

	struct windowPair windowPair;

	//int windowsNr_1;
	/*segment indices (bounds of each window considered)*/
	//int *segIndices_1;

	//int windowsNr_2;
	/*segment indices (bounds of each window considered)*/
	//int *segIndices_2;

	double distMin;

} matchCandidate_windowPair;

/*Struct for storing final scores per query*/ 
typedef struct matchScore{

	/*char *queryFileName;
	
	char *queryChainId;*/

	char *matchFileName;

	char *matchStructureName; /*name of structure in PDB-jargon, particular to application, so e.g. domain name in CATH*/
	char *matchClassId; /*classification id in e.g. CATH or SCOP*/
	
	char *matchChainId;

	int  matchChainLen;

	int matchNrOfWindows;

	int *places; //for keeping the number of 1st, 2nd, ... placed window/window pair matches of the given match structure (matchFileName)

	int *queryWindowNr; //for keeping the window number of each matched query-window.
	int *windowNr; //for keeping the window number (in the matching structure) matching the matched query-window.

	float queryCoverage; //a measure of how large a fraction of the query is covered by the match 

	double countScore;/*to contain a score based on places array (right above) and normalized e.g. by length of match candidate (usually)*/

	double aggrWinScore; /*a score obtained by summing up the scores for individual windows*/ 

	double score; /*to contain a score based on the distance of the invariant values and normalized e.g. by length of match candidate (usually)*/

} matchScore;

/*struct for holding array of invariant values on windows plus a few descriptives 
for each structure in a data base, containing I-values for use in structural alignment
(here an excerpt of invariants up to order 2)*/
typedef struct db_I_windows_order2_ptr{

	char *fileName;

	char *structureName; /*name of structure in PDB-jargon, particular to application, so e.g. domain name in CATH*/
	char *classId; /*classification id in e.g. CATH or SCOP*/

	char *chainId;
	int  chainNr;
	int  chainLen;

	int  nrOfWindows;

	struct window *window;
	
	//int *windowsNr;
	/*segment indices (bounds of each window considered)*/
	//int **segIndices;

	/*pointers for order 1 to 3 measures:*/
	double *I12;
	double *Ia12;
	
	double *I1423;

	/* for full computation of the order 2 measures across the simplex:*/
    double *I1234_full;
    double *I1324_full;

	double *Ia1234_full;
	double *I12a34_full;
	double *Ia12a34_full;

	double *Ia1324_full;
	double *I13a24_full;
	double *Ia13a24_full;

	double *Ia14a23;
	double *Ia1423;
	double *I14a23;

} db_I_windows_order2_ptr;

/*struct for holding flattened version of the records in e.g. db_I_windows_order2_ptr, with I-values
binned; db_I_windows_order2_ptr will then consist of a ptr to the binned data. To be used for fast 
look-up of I-values vector*/
typedef struct binned_I_windows{

	//int recKey; /*to hold an id of the binned record; the recKey will run from 0 and up to the total nr of windows in the underlying I_windows_ptr*/

	int fileNr; /*id ref to the fileName, chainId in the underlying I_windows_ptr*/ 

	int windowNr; /*windowNr in the underlying I_windows_ptr; the pair (fileNr, windowNr) identifies record in  the underlying I_windows_ptr*/ 
	
	int *binnedIvector; /*vector of integers resulting from binning of the present I_windows_ptr record*/ 

} binned_I_windows;

/*As binned_I_windows but for pairs of windows*/
typedef struct binned_I_windowPairs{

	//int recKey; /*to hold an id of the binned record; the recKey will run from 0 and up to the total nr of windows in the underlying I_windows_ptr*/

	int fileNr; /*id ref to the fileName, chainId in the underlying I_windows_ptr*/ 

	int windowNr_1; /*windowNr in the underlying I_windows_ptr; the pair (fileNr, windowNr) identifies record in  the underlying I_windows_ptr*/ 

	int windowNr_2; /*windowNr in the underlying I_windows_ptr; the pair (fileNr, windowNr) identifies record in  the underlying I_windows_ptr*/ 

	int *binnedIvector; /*vector of integers resulting from binning of the present I_windows_ptr record*/ 

} binned_I_windowPairs;


/*For storing the counts of "word occurences"*/
typedef struct binned_I_counts{

	int *binnedIvector; /*vector of integers resulting from binning of the present I_windows_ptr record*/ 

	int count;

} binned_I_counts;


/*as db_I_windows_order2_ptr, but for window pairs*/
typedef struct db_I_windowPairs_order2_ptr{

	char *fileName;

	char *structureName; /*name of structure in PDB-jargon, particular to application, so e.g. domain name in CATH*/
	char *classId; /*classification id in e.g. CATH or SCOP*/

	
	char *chainId;
	int  chainNr;
	int  chainLen;

	int  nrOfWindows;
	int  nrOfWindowPairs;

	struct windowPair **windowPair;

	//int *windowsNr_1;
	///*segment indices (bounds of each window considered)*/
	//int **segIndices_1;

	//int *windowsNr_2;
	///*segment indices (bounds of each window considered)*/
	//int **segIndices_2;

	/*pointers for order 1 to 3 measures:*/
	double **I12;
	double **Ia12;
	
	double **I1423;

	///* for full computation of the order 2 measures across the simplex:*/
	//double * I1234_full_aid;
	//double * I1324_full_aid;
    double **I1234_full;
    double **I1324_full;

	double **Ia1234_full;
	double **I12a34_full;
	double **Ia12a34_full;

	double **Ia1324_full;
	double **I13a24_full;
	double **Ia13a24_full;

	double **Ia1423;
	double **I14a23;
	double **Ia14a23;

} db_I_windowPairs_order2_ptr;

/*for holding I-values on windows for a set of records (aimed at computing mean and std dev for each invariant)*/
typedef struct I_windows_order2_raw_ptr{

	/*pointers for order 1 to 2 measures:*/
	double *I12;
	double *Ia12;
	
	double *I1423;

	/* for full computation of the order 2 measures across the simplex:*/
    double *I1234_full;
    double *I1324_full;

	double *Ia1234_full;
	double *I12a34_full;
	double *Ia12a34_full;

	double *Ia1324_full;
	double *I13a24_full;
	double *Ia13a24_full;

	double *Ia14a23;
	double *Ia1423;
	double *I14a23;

} I_windows_order2_raw_ptr;

/*for holding I-values on window pairs for a set of records (aimed at computing mean and std dev for each invariant)*/
typedef struct I_windowPairs_order2_raw_ptr{

	/*pointers for ord	er 1 to 2 measures:*/
	double **I12;
	double **Ia12;
	
	double **I1423;

	///* for full computation of the order 2 measures across the simplex:*/
	//double * I1234_full_aid;
	//double * I1324_full_aid;
    double **I1234_full;
    double **I1324_full;

	double **Ia1234_full;
	double **I12a34_full;
	double **Ia12a34_full;

	double **Ia1324_full;
	double **I13a24_full;
	double **Ia13a24_full;

	double **Ia14a23;
	double **Ia1423;
	double **I14a23;



} I_windowPairs_order2_raw_ptr;

typedef struct I_windows_order2_meanStddev{

	/*pointers for order 1 to 3 measures:*/
	double mean_I12;
	double mean_Ia12;
	
	double mean_I1423;

	/* for full computation of the order 2 measures across the simplex:*/
    double mean_I1234_full;
    double mean_I1324_full;

	double mean_Ia12a34_full;
    double mean_Ia13a24_full;
	double mean_Ia14a23;

	double mean_Ia1234_full;
    double mean_Ia1324_full;
	double mean_Ia1423;

	double mean_I12a34_full;
    double mean_I13a24_full;
	double mean_I14a23;


	/*pointers for order 1 to 3 measures:*/
	double stddev_I12;
	double stddev_Ia12;
	
	double stddev_I1423;

	///* for full computation of the order 2 measures across the simplex:*/
    double stddev_I1234_full;
    double stddev_I1324_full;

	double stddev_Ia12a34_full;
    double stddev_Ia13a24_full;
	double stddev_Ia14a23;

	double stddev_Ia1234_full;
    double stddev_Ia1324_full;
	double stddev_Ia1423;

	double stddev_I12a34_full;
    double stddev_I13a24_full;
	double stddev_I14a23;

} I_windows_order2_meanStddev;

typedef struct I_windowPairs_order2_meanStddev{

	/*pointers for order 1 to 3 measures:*/
	double mean_I12;
	double mean_Ia12;
	
	double mean_I1423;

	///* for full computation of the order 2 measures across the simplex:*/
	//double * I1234_full_aid;
	//double * I1324_full_aid;
    double mean_I1234_full;
    double mean_I1324_full;

	double mean_Ia12a34_full;
    double mean_Ia13a24_full;
	double mean_Ia14a23;

	double mean_Ia1234_full;
    double mean_Ia1324_full;
	double mean_Ia1423;

	double mean_I12a34_full;
    double mean_I13a24_full;
	double mean_I14a23;

	/*pointers for order 1 to 3 measures:*/
	double stddev_I12;
	double stddev_Ia12;
	
	double stddev_I1423;

	///* for full computation of the order 2 measures across the simplex:*/
    double stddev_I1234_full;
    double stddev_I1324_full;

	double stddev_Ia12a34_full;
    double stddev_Ia13a24_full;
	double stddev_Ia14a23;

	double stddev_Ia1234_full;
    double stddev_Ia1324_full;
	double stddev_Ia1423;

	double stddev_I12a34_full;
    double stddev_I13a24_full;
	double stddev_I14a23;

} I_windowPairs_order2_meanStddev;

typedef struct matchRange{

	int start;
	int end;
	int nrOfMismatches;

} matchRange;


#endif