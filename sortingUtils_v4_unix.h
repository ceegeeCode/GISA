/* **********************************************************************************************************************************
*
*
* ***************************************  Header for: Sorting utils for GISA_v* and SAonGISA_v* *******************************************************************************
*
* See main code (GISA_v*) for do and how to cite.
*

Author: Christian Grønbæk
Version: v1, 3rd July 2017; v2 Nov 2017; _v3 Jan 2019

This version: _v4, > 1st Jan, 2019.

***********************************************************************************************/

#ifndef SORTINGUTILS_V4_UNIX
#define SORTINGUTILS_V4_UNIX


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include "GISA_v10_unix.h"
#include "dataTypesMatching_v3_unix.h"
#include "mathUtils.h"


/**************************************************
************ Prototype functions ******************
***************************************************/


void bubbleSort(struct matchCandidate *, int);

void bubbleSort2(struct matchCandidate *, int );

void bubbleSort2_windowPairs(struct matchCandidate_windowPair *, int);


/*scores bubble sort, left to rigth*/
void bubbleSortScoresLR(struct matchScore *, int );

void bubbleSortScoresRL(struct matchScore *, int );

/*DOESN'T WORK! 
modified version of the lexicographic ordering allowing
a pre-set total nr of mismatches; THE INTENTION WAS ONLY TO ALLOW a mismatch 
of 1 per letter, but that results in an "order" which isn't transitive!
The same holds for the present version (which just allows some nr of mismatches) 
--- it's not a proper order.
*/
int lexicographicMod(int *, int *, int *, int , int );


void bubbleSortLexicoBinnedDB_LR(struct binned_I_windows *, int , int , int );

/*Sort binned DB results by lexicographic ordering;
nrOfDBrecords: the number of all windows in the DB (for which there is a set of I-values), summed up over
all structures (file, chain) that is*/
void bubbleSortLexicoBinnedDB_RL(struct binned_I_windows *, int , int , int );

void bubbleSort_binned_I_counts(struct binned_I_counts *, int , int );

/*heap sort match candidates (on the dist-value "distMin")*/
void siftDownMatchCandidates(struct matchCandidate *, int, int);

void heapifyMatchCandidates(struct matchCandidate *, int );

void heapSortMatchCandidates(struct matchCandidate *, int );



/*heap sort for binned DB*/
void siftDownBinnedDB(struct binned_I_windows *, int , int , int );

void heapifyBinnedDB(struct binned_I_windows *, int , int );

void heapSortBinnedDB(struct binned_I_windows *, int , int );

/*heap sort for binned DB pairs*/ 
void siftDownBinnedDBPairs(struct binned_I_windowPairs *, int , int , int );

void heapifyBinnedDBPairs(struct binned_I_windowPairs *, int , int );

void heapSortBinnedDBPairs(struct binned_I_windowPairs *, int , int );



/*heap sort rarity scores*/
void siftDownRarityScores(struct Irarity_windows *, int , int);

void heapifyRarityScores(struct Irarity_windows *, int);

void heapSortRarityScores(struct Irarity_windows *, int);


/*heap sort rarity scores, for pairs*/
void siftDownRarityScoresPairs(struct Irarity_windowPairs *, int , int);

void heapifyRarityScoresPairs(struct Irarity_windowPairs *, int);

void heapSortRarityScoresPairs(struct Irarity_windowPairs *, int);


/*heap sort array of raw scores (over set of queries)*/
void siftDownRawScores(struct queryRawScore *ptr_queryRawScore, int start, int end);

void heapifyRawScores(struct queryRawScore *ptr_queryRawScore, int nrOfStructures);

void heapSortRawScores(struct queryRawScore *ptr_queryRawScore, int nrOfStructures);


/*heap sort for match set*/
void siftDownMatchSet(struct binned_I_windows *, int , int , int );

void heapifyMatchSet(struct binned_I_windows *, int , int );

void heapSortMatchSet(struct binned_I_windows *, int , int );

/*heap sort for match scores set*/
void siftDownMatchScores(struct matchScore *, int , int , int , int );

void heapifyMatchScores(struct matchScore *, int , int , int );

void heapSortMatchScores(struct matchScore *, int , int , int );



int bisectionalSearchMatchScores(char *, struct matchScore *, int , int );

/*Search through the (lexicographically pre-sorted) DB of binned I-vectors for a given
query (an I-window vector repr as binned I-vector)
returns: first and last index in sorted list (ptr_binned_I_windows_DB) matching the query*/
int bisectionalSearch(int *, struct binned_I_windows , struct binned_I_windows *, int , int , int );

/*older version*/
int bisectionalSearch_1(int *, struct binned_I_windows , struct binned_I_windows *, int , int , int );

/*version of bisectional search which finds the range of hits directly, i.e. without out 
the final left- and right extensions.
29/3 2017: changed by adding "&& mid != lastRight" in first if-clause */
int bisectionalSearch2(int *, struct binned_I_windows , struct binned_I_windows *, int , int , int );

/*bisectional search for "letter" at position position of query, binned_I_windows_query, among the letters at the same position
in the range of the array ptr_binned_I_windows_DB from  ptr_range.start to ptr_range.end.
Returns: index of the hit; updates ptr_range.start and .end with the last left resp. right interval bounds of the search*/
int bisectSingle(struct binned_I_windows , int , int *, struct binned_I_windows *);

/*Before Jan 22, '18 (new version includes a trailing if-clause!)*/
int bisectSingle_old_2(struct binned_I_windows , int , int *, struct binned_I_windows *);

/*Before april 4, '17:*/
int bisectSingle_old(struct binned_I_windows , int , int *, struct binned_I_windows *);


/*for extending left or right from a hit (with index hitIndex); searches for extension left from ptr_range.start and up to 
hitIndex, and right from hitIndex and to ptr_range.end. In particular hitIndex must sit in the interval [ptr_range.start,
ptr_range.end] for this to work.
Input: leftRigth_b -- if 0 both extensions are done; if -1 only left, if 1 only right.
Transforms: updates ptr_range .start and .end to left resp right indexes found.*/
int bisectSingleExtend(int hitIndex, struct binned_I_windows , int , int *, struct binned_I_windows *, int , int );

/*Before april 4, '17:*/
int bisectSingleExtend_old(int hitIndex, struct binned_I_windows , int , int *, struct binned_I_windows *, int );


/*version of bisectional search which finds the range of hits directly, i.e. without out the final left- and right extensions*/
int bisectionalSearchSingle(int *, int , struct matchRange **, struct binned_I_windows , int , struct matchRange , struct binned_I_windows *, int , int );

/*Before Jan 22, '18; extension from the ends in case no hit was found was missing */
int bisectionalSearchSingle_old_2(int *, int , struct matchRange **, struct binned_I_windows , int , struct matchRange , struct binned_I_windows *, int , int );

/*Before april 4, '17:*/
int bisectionalSearchSingle_old(int *, int , struct matchRange **, struct binned_I_windows , int , struct matchRange , struct binned_I_windows *, int , int );


/*bisectional fct's, pair versions:*/
int bisectSingle_Pairs(struct binned_I_windowPairs , int , int *, struct binned_I_windowPairs *);

/*for extending left or right from a hit (with index hitIndex); searches for extension left from ptr_range.start and up to 
hitIndex, and right from hitIndex and to ptr_range.end. In particular hitIndex must sit in the interval [ptr_range.start,
ptr_range.end] for this to work.
Input: leftRigth_b -- if 0 both extensions are done; if -1 only left, if 1 only right.
Transforms: updates ptr_range .start and .end to left resp right indexes found.*/
int bisectSingleExtend_Pairs(int , struct binned_I_windowPairs , int , int *, struct binned_I_windowPairs *, int , int );

int bisectionalSearchSingle_Pairs(int *, int , struct matchRange **, struct binned_I_windowPairs , int , struct matchRange , struct binned_I_windowPairs *, int , int );


int wordSearch_I_binned_I_counts(int *, struct binned_I_counts *, int , int , int);

#endif