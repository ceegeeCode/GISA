/* **********************************************************************************************************************************
*
*
* ***************************************  Header for GISA_v10 *******************************************************************************
*
* See main code (GISA_v*) for do and how to cite.
*
* Author: Christian Grønbæk
*
**********************************************************************************************/


#ifndef GISA_V10_unix
#define GISA_V10_unix

//#pragma once 

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>
#include<unistd.h>
#include<sys/types.h>
#include<dirent.h>
#include "dataTypesMatching_v3_unix.h"

#define pi acos(-1.0)
#define pathLength 2048
#define fileNameLength 500
#define maxNrOfChains 100
#define maxNrOfDomains 100
#define lengthCLFformat20 7 //number of char's in CATH file List Format v 2.0
#define fileNameCharSize 1000
#define structureNameCharSize 100
#define classIdCharSize 50
#define chainIdCharSize 20
#define classifierSize 20 /*for classifier, e.g. SCOP, as array of int's: max allowed length of this array (for SCOP it's 4) */
#define stdRealSegLength 49 /* square of max length of segments in alphaC-trace; usually about 4 Aangstroem should do; we set it somewhat higher  */
#define strLengthCutOff 10 /*chains shorter than this int will not be included in computation of invariants*/

/**************************************************
******************* Structs ***********************
***************************************************/

typedef struct double2{

	double val[2];

} double2;

/*struct for holding content of directory*/
typedef struct dirContent{

	int numberOfFiles;
	char ** ptr_dirList;
	char ** ptr_fileNameList;

} dirContent;

/*Struct or keeping info on SCOP, CATH, .. class. 
For SCOP the classSize wil be the number of id's with the same class identiier and fold id-number 
(in x.1.2.3 x is the class identifier and 1 is the fold number).*/
typedef struct classInfo{

	char *classId; /*classification id in e.g. CATH or SCOP*/

	int *classIdConv; /*for converting a class id to an array of int's; e.g. a SCOP class a.1.2.3.4 to [0,1,2,3,4]*/

	int classSize; /*number of chains in the class with id classId; helpful e.g. for pinpointing classes of a single element ("singletons")*/

} classInfo;


/*struct for keeping rough info about possible multi-mer PDB-file: the number of chains in the structure
and an array indexed by chain number and with entries (chainId, chainlength), ptr_chainInStr: */

typedef struct chainInfo{

	char *structureName; /*name of structure in PDB-jargon, particular to application, so e.g. domain name in CATH*/
	
	char *classId; /*classification id in e.g. CATH or SCOP*/
	int *classIdConv; /*for converting a class id to an array of int's; e.g. a SCOP class a.1.2.3.4 to [0,1,2,3,4]*/

	char *chainId;
	//char chainId[2]; /*we shall assume that the chain id is always a single character (may be a blank)*/ 
	int chainLength;


} chainInfo;

typedef struct chainInStructure{

	int numberOfChains;
	struct chainInfo *ptr_chainInStr;

 } chainInStructure;

typedef struct cAlphaCoords{
	double x;
	double y;
	double z;

} cAlphaCoords;

/*struct for C-alpha providning residue nr and 3d-coords*/
typedef struct cAlpha{
	
	struct cAlphaCoords coords;

	int residueNr;

} cAlpha;

/* struct for holding the two points determining a chain segment*/
typedef struct segment{
	struct cAlphaCoords s1;
	struct cAlphaCoords s2;

} segment;

/*struct for holding two segments (writhe contribution works on pairs of segments):*/
typedef struct twoSegments{

	struct segment seg1;

	struct segment seg2;
} twoSegments;

/*two segments given as twelve doubles (four 3d points):*/
typedef struct twoSegmentCoords{

	/* coords of 1st point in 1st segment*/
	double s101;
	double s102;
	double s103;
	/* coords of 2nd point in 1st segment*/
	double s111;
	double s112;
	double s113;
	/* coords of 1st point in 2nd segment*/
	double s201;
	double s202;
	double s203;
	/* coords of 2nd point in 2nd segment*/
	double s211;
	double s212;
	double s213;
} twoSegmentCoords;

/*struct for holding the indices of two segments*/
typedef struct twoSegmentIndex{

	int idx1;
	int idx2;

} twoSegmentIndex;


/*struct for holding double array of invariant values plus a few descriptives of the structure*/
typedef struct I_ptr{

	char *fileName;
	
	char *structureName; /*name of structure in PDB-jargon, particular to application, so e.g. domain name in CATH*/
	char *classId; /*classification id in e.g. CATH or SCOP

	char *chainId; /*we shall assume that the chain id is always a single character (may be a blank)*/ 

	char *chainId;
	int *chainNr;
	int *chainLen;

	int *order;

	/*pointer to w-values*/
	double **wVal;
	/*pointers for order 1 to 3 measures:*/
	double **I12;
	double **Ia12;
	
	double **I1234;
	double **I1324;
	double **I1423;

	/* for full computation of the order 2 measures across the simplex:*/
	double **I1234_full_aid;
	double **I1324_full_aid;
    double **I1234_full;
    double **I1324_full;
	/* assisting the order 3 case:*/
	double **I1324_full2_aid;
	double **I1324_full2;
	double **I1423_full0;
	double **I1423_full2_aid;
	double **I1423_full2;

	/* for "absolute value version" of order 2 measures: */
	double **Ia1234;
	double **I12a34;
	double **Ia12a34;

	double **Ia1324;
	double **I13a24;
	double **Ia13a24;

	double **Ia1423;
	double **I14a23;
	double **Ia14a23;

	/* for full computation of the order 2 measures across the simplex:*/
	double **Ia12a34_full_aid;
	double **Ia1234_full_aid;
	double **I12a34_full_aid;

	double **Ia13a24_full_aid;
	double **Ia1324_full_aid;
	double **I13a24_full_aid;

	double **Ia1234_full;
	double **I12a34_full;
	double **Ia12a34_full;

	double **Ia1324_full;
	double **I13a24_full;
	double **Ia13a24_full;


	/*order 3*/
	/*12*/
	double **I123456;
	double **I123645;
	double **I123546;
	/*13*/
	double **I132456;
	double **I132546;
	double **I132645;
	/*14*/
	double **I142356;
	double **I142536;
	double **I142635;
	/*15*/
	double **I152346;
	double **I152436;
	double **I152634;
	/*16*/
	double **I162345;
	double **I162435;
	double **I162534;


} I_ptr;


/*Struct for holding the values of the invariants for one structure (one PDB file):*/
typedef struct I_values{

	char *fileName;

	char *structureName; /*name of structure in PDB-jargon, particular to application, so e.g. domain name in CATH*/
	char *classId; /*classification id in e.g. CATH or SCOP

	char *chainId; /*we shall assume that the chain id is always a single character (may be a blank)*/ 

	char *chainId;
	int chainNr;
	int chainLen;

	int order;

	int pertNo; /*index -- or number -- of perturbation*/

	double cpuTime;
	double aggrTime;

	/*doubles for order 1 to 3 measures:*/
	double I12;
	double Ia12;
	
	double I1234;
	double I1324;
	double I1423;

	/* for "absolute value version" of order 2 measures: */
	double Ia1234;
	double I12a34;
	double Ia12a34;

	double Ia1324;
	double I13a24;
	double Ia13a24;

	double Ia1423;
	double I14a23;
	double Ia14a23;

	/*order 3*/
	/*12*/
	double I123456;
	double I123645;
	double I123546;
	/*13*/
	double I132456;
	double I132546;
	double I132645;
	/*14*/
	double I142356;
	double I142536;
	double I142635;
	/*15*/
	double I152346;
	double I152436;
	double I152634;
	/*16*/
	double I162345;
	double I162435;
	double I162534;

} I_values;



/*struct for holding a set of "coordinate arrays" as given by the set of 
pairs of segment end points in a chain; there are therefore 12 coordinates 
(two segments consist have four end points, each having three coord's) */
typedef struct Segment_ptr{

	/*struct twoSegmentCoords *ptr_segmentPairCoords;
	struct twoSegmentCoords segCoords;*/
	double *ptr_segmentPairCoords_s101;
	double *ptr_segmentPairCoords_s102;
	double *ptr_segmentPairCoords_s103;

	double *ptr_segmentPairCoords_s111;
	double *ptr_segmentPairCoords_s112;
	double *ptr_segmentPairCoords_s113;

	double *ptr_segmentPairCoords_s201;
	double *ptr_segmentPairCoords_s202;
	double *ptr_segmentPairCoords_s203;

	double *ptr_segmentPairCoords_s211;
	double *ptr_segmentPairCoords_s212;
	double *ptr_segmentPairCoords_s213;
} Segment_ptr;

/*struct for holding output from closed loop search, per chain:*/
typedef struct closedLoopCharacteristic{

	char *fileName;
	char *description; /*will be either "link" or "poke"*/
	struct twoSegmentIndex closedLoop; /*the two segment indices define the start and end of the sub chain, thought to be a closed loop*/
	struct twoSegmentIndex subChain; /*the two segment indices define the start and end of the sub chain*/
	double I12; /*value of writhe of the (closedLoop, subChain) pair*/

} closedLoopCharacteristic;


/*struct for holding array of invariant values on windows plus a few descriptives of the structure*/
typedef struct I_windows_ptr{

	char *fileName;
	
	char *structureName; /*name of structure in PDB-jargon, particular to application, so e.g. domain name in CATH*/
	char *classId; /*classification id in e.g. CATH or SCOP

	char *chainId; /*we shall assume that the chain id is always a single character (may be a blank)*/ 
	char *chainId;
	int  chainNr;
	int  chainLen;

	int  order;

	int  windowLgth;
	int  stepSize;
	int  nrOfWindows;

	struct window *window;

	//int *windowNr;
	///*segment indices (bounds of each window considered)*/
	//int **segIndices;


	///*pointer to w-values*/
	//double *wVal;

	/*pointers for order 1 to 3 measures:*/
	double *I12;
	double *Ia12;
	
	double *I1234;
	double *I1324;
	double *I1423;

	///* for full computation of the order 2 measures across the simplex:*/
	//double * I1234_full_aid;
	//double * I1324_full_aid;
    double *I1234_full;
    double *I1324_full;
	///* assisting the order 3 case:*/
	//double * I1324_full2_aid;
	//double * I1324_full2;
	//double *I1423_full0;
	//double *I1423_full2_aid;
	//double *I1423_full2;

	/* for "absolute value version" of order 2 measures: */
	double *Ia1234;
	double *I12a34;
	double *Ia12a34;
	
	double *Ia1324;
	double *I13a24;
	double *Ia13a24;

	double *Ia1423;
	double *I14a23;
	double *Ia14a23;

	/* for full computation of the order 2 measures across the simplex:*/
	double *Ia1234_full;
	double *I12a34_full;
	double *Ia12a34_full;

	double *Ia1324_full;
	double *I13a24_full;
	double *Ia13a24_full;


	/*order 3*/
	/*12*/
	double *I123456;
	double *I123546;
	double *I123645;
	/*13*/
	double *I132456;
	double *I132546;
	double *I132645;
	/*14*/
	double *I142356;
	double *I142536;
	double *I142635;
	/*15*/
	double *I152346;
	double *I152436;
	double *I152634;
	/*16*/
	double *I162345;
	double *I162435;
	double *I162534;


} I_windows_ptr;

/*struct for holding array of invariant values on pairs windows, incuding mutual values, plus a few descriptives of the structure*/
typedef struct I_windowPairs_ptr{

	char *fileName;
	
	char *structureName; /*name of structure in PDB-jargon, particular to application, so e.g. domain name in CATH*/
	char *classId; /*classification id in e.g. CATH or SCOP

	char *chainId; /*we shall assume that the chain id is always a single character (may be a blank)*/ 

	char *chainId;
	int  chainNr;
	int  chainLen;

	int  order;

	int  windowLgth;
	int  stepSize;
	int  nrOfWindows;
	int  nrOfWindowPairs;

	struct windowPair **windowPair;

	//int *windowsNr_1;
	///*segment indices (bounds of each window considered)*/
	//int **segIndices_1;

	//int *windowsNr_2;
	///*segment indices (bounds of each window considered)*/
	//int **segIndices_2;

	/////*pointer to w-values*/
	////double *wVal;

	///* For 1st window:*/
	///*pointers for order 1 to 2 measures:*/
	//double *I12_1;
	//double *Ia12_1;
	//
	//double *I1234_1;
	//double *I1324_1;
	//double *I1423_1;

	/////* for full computation of the order 2 measures across the simplex:*/
	////double * I1234_full_aid_1;
	////double * I1324_full_aid_1;
 //   double *I1234_full_1;
 //   double *I1324_full_1;
	/////* assisting the order 3 case:*/
	////double * I1324_full2_aid_1;
	////double * I1324_full2_1;
	////double *I1423_full0_1;
	////double *I1423_full2_aid_1;
	////double *I1423_full2_1;

	///* for "absolute value version" of order 2 measures: */
	//double *Ia1234_1;
	//double *I12a34_1;
	//double *Ia12a34_1;
	//
	//double *Ia1324_1;
	//double *I13a24_1;
	//double *Ia13a24_1;

	//double *Ia1423_1;
	//double *I14a23_1;
	//double *Ia14a23_1;


	///* For 2nd window:*/
	///*pointers for order 1 to 2 measures:*/
	//double *I12_2;
	//double *Ia12_2;
	//
	//double *I1234_2;
	//double *I1324_2;
	//double *I1423_2;

	/////* for full computation of the order 2 measures across the simplex:*/
	////double * I1234_full_aid_2;
	////double * I1324_full_aid_2;
 //   double *I1234_full_2;
 //   double *I1324_full_2;
	/////* assisting the order 3 case:*/
	////double * I1324_full2_aid_2;
	////double * I1324_full2_2;
	////double *I1423_full0_2;
	////double *I1423_full2_aid_2;
	////double *I1423_full2_2;

	///* for "absolute value version" of order 2 measures: */
	//double *Ia1234_2;
	//double *I12a34_2;
	//double *Ia12a34_2;
	//
	//double *Ia1324_2;
	//double *I13a24_2;
	//double *Ia13a24_2;

	//double *Ia1423_2;
	//double *I14a23_2;
	//double *Ia14a23_2;


	
	/* For (window1, window2) mutual values:*/
		/* For (window1, window2) mutual values:*/
	/*These pointers are "double arrays": each index 
	is for a particular window, and the index pair will thus 
	define a pair of windows; the aim is to store values of
	mutual kind; for index1 = index2 the value will be that of
	the invariant on the window (at index1) itself.*/
	/*pointers for order 1 to 2 measures:*/
	double **I12;
	double **Ia12;
	
	double **I1234;
	double **I1324;
	double **I1423;

	///* for full computation of the order 2 measures across the simplex:*/
	//double **I1234_full_aid;
	//double **I1324_full_aid;
    double **I1234_full;
    double **I1324_full;
	///* assisting the order 3 case:*/
	//double **I1324_full2_aid;
	//double **I1324_full2;
	//double **I1423_full0;
	//double **I1423_full2_aid;
	//double **I1423_full2;

	/* for "absolute value version" of order 2 measures: */
	double **Ia1234;
	double **I12a34;
	double **Ia12a34;
	

	double **Ia1324;
	double **I13a24;
	double **Ia13a24;

	double **Ia1423;
	double **I14a23;
	double **Ia14a23;

	/* for full computation of the order 2 measures across the simplex:*/
	double **Ia1234_full;
	double **I12a34_full;
	double **Ia12a34_full;

	double **Ia1324_full;
	double **I13a24_full;
	double **Ia13a24_full;


} I_windowPairs_ptr;

typedef struct Irarity_windows{

	struct window window;
	
	//int windowsNr;
	///*segment indices (bounds of each window considered)*/
	//int *segIndices;

	double rarityScore;

} Irarity_windows;

/*struct for holding array of rarity scores on windows plus a few descriptives of the structure*/
typedef struct Irarity_windows_ptr{

	char *fileName;
	
	char *structureName; /*name of structure in PDB-jargon, particular to application, so e.g. domain name in CATH*/
	char *classId; /*classification id in e.g. CATH or SCOP

	//char *chainId; /*we shall assume that the chain id is always a single character (may be a blank)*/ 
	char *chainId;
	int  chainNr;
	int  chainLen;

	int  order;

	int  windowLgth;
	int  stepSize;
	int  nrOfWindows;

	struct window *window;

	//int *windowsNr;
	///*segment indices (bounds of each window considered)*/
	//int **segIndices;


	/*pointers for order 1 to 3 measures:*/
	double *rarityScore;

} Irarity_windows_ptr;

typedef struct Irarity_windowPairs{

	struct windowPair windowPair;

	//int windowsNr_1;

	///*segment indices (bounds of each window considered)*/
	//int *segIndices_1;

	//int windowsNr_2;

	///*segment indices (bounds of each window considered)*/
	//int *segIndices_2;

	double rarityScore;

} Irarity_windowPairs;

/*struct for holding array of rarity scores on pairs windows, incuding mutual values, plus a few descriptives of the structure*/
typedef struct Irarity_windowPairs_ptr{

	char *fileName;
	
	char *structureName; /*name of structure in PDB-jargon, particular to application, so e.g. domain name in CATH*/
	char *classId; /*classification id in e.g. CATH or SCOP

	char *chainId; /*we shall assume that the chain id is always a single character (may be a blank)*/ 

	char *chainId;
	int  chainNr;
	int  chainLen;

	int  order;

	int  windowLgth;
	int  nrOfWindows;
	int  nrOfWindowPairs;

	struct windowPair **windowPair;

	//int *windowsNr_1;
	///*segment indices (bounds of each window considered)*/
	//int **segIndices_1;

	//int *windowsNr_2;
	///*segment indices (bounds of each window considered)*/
	//int **segIndices_2;

	/*These pointers are "double arrays": each index 
	is for a particular window, and the index pair will thus 
	define a pair of windows; the aim is to store values of
	mutual kind; for index1 = index2 the value will be that of
	the invariant on the window (at index1) itself.*/
	double **rarityScore;


} Irarity_windowPairs_ptr;




/**************************************************
************ Prototype functions ******************
***************************************************/


//max and min are not part of standard C:
int min(int a, int b);

int max(int a, int b);

/*Util fct: input a filename with extension; returns (in/out params): file name wo extension
and the extension*/
void splitFileName(char *fileNameWExt, char *fileName, char *ext);


/*functions for getting the contents of a named directory*/
struct dirContent ListDirectoryContents(char *);

int ListDirectoryContents2(char *, struct dirContent *, int *);

/*functions for reading content of a PDB-file and outputting the array of C-alpha coordinates
The first, readPDBChainStructure, gets the number of chains in the file (structure) and the length of 
each if the chains and returns this info in an array {...,(chainNr, chainLength), ...}. The 
second, main_readPDB, reads in the array of C-alpha coordinates in a given file (structure) 
along with chain length and chain number.*/
struct chainInStructure readPDBChainStructure(FILE *);

int readPDBChainStructure2(FILE * , struct chainInStructure * );

int readPDBChainStructureSCOP(FILE * , struct chainInStructure * );

int convertSCOP(struct classInfo *, int );

int convertCATH(struct classInfo *, int);

int readPDBDomainStructureCATH(FILE * , struct chainInStructure *);

int readCATHDomainList(FILE * , struct chainInfo **, int  );

struct cAlpha * main_readPDB(FILE *, int, int);

int main_readPDB2(FILE * , struct cAlpha *, int , int );


/*function computing the square of the distance between two C-alphas; used for finding closed loops
in a structure; we compute the square dist rather than the dist simply to avoid a sqrt (which is
unneccessary for our purposes)*/
double distCalphas(struct cAlphaCoords, struct cAlphaCoords);

double I_mutual(double **, int , int , int , int , double );

/*for computing mean and std dev of 1-dim array*/
struct double2 meanStddev(double *, int);

/*simple fct scaling all entries in a 1-d array*/
void scaleArray(double , double **, int );

/*simple fct normalizing all entries in a 1-d array to mean zero and stddev 1*/
void normalizeArray(double , double , double **, int );

/*simple fct normalizing all entries in a 2-d "simplex" array to mean zero and stddev 1;
simplex means: the 2-d array has only entries at (i,j) having i <= j and with i ranging from
0 to inOutSimplexSideLength and j up to inOutSimplexSideLength.*/
void normalizeSimplexArray(double mean, double stddev, double ***inOutArray, int simplexSideLength);

/*function ordering two input vectors of integers, v and w, lexicographically; 
returns: 1 if v < w, -1 if v > w and 0 if v == w
nrOfEntries: length of v
*/
int lexicographic(int *, int *, int );

int lexicographicChar(char *, char *, int );

void siftDown1dArray(double *ptr_values, int start, int end, int decreasing_b);

void heapify1dArray(double *ptr_values, int nrOfRecords, int decreasing_b);

void heapSort1dArray(double *ptr_values, int nrOfRecords, int decreasing_b);

int bisectionalSearchValueL(double value, double *ptr_values, int nrOfRecords);

/*function siftDown(a, start, end) is
   (end represents the limit of how far down the heap to sift)
   root := start

   while root * 2 + 1 ≤ end do       (While the root has at least one child)
      child := root * 2 + 1           (root*2+1 points to the left child)
      (If the child has a sibling and the child's value is less than its sibling's...)
      if child + 1 ≤ end and a[child] < a[child + 1] then
         child := child + 1           (... then point to the right child instead)
      if a[root] < a[child] then     (out of max-heap order)
         swap(a[root], a[child])
         root := child                (repeat to continue sifting down the child now)
      else
         return
		*/
void siftDownDBClass(struct classInfo *, int , int , int , int );

/*function heapify(a,count) is
   (start is assigned the index in a of the last parent node)
   start := (count - 2) / 2
   
   while start ≥ 0 do
      (sift down the node at index start to the proper place
       such that all nodes below the start index are in heap
       order)
      siftDown(a, start, count-1)
      start := start - 1
   (after sifting down the root all nodes/elements are in heap order)
 */
void heapifyDBClass(struct classInfo *, int , int , int );

/*function heapSort(a, count) is
   input: an unordered array a of length count
 
   (first place a in max-heap order)
   heapify(a, count)
 
   end := count - 1
   while end > 0 do
      (swap the root(maximum value) of the heap with the
       last element of the heap)
      swap(a[end], a[0])
      (decrement the size of the heap so that the previous
       max value will stay in its proper place)
      end := end - 1
      (put the heap back in max-heap order)
      siftDown(a, 0, end)
*/
void heapSortDBClass(struct classInfo *, int , int , int );

/*Search through the lexicographically pre-sorted list of classifiers (supposed to be represented as int-arrays of length
nrOfEntries for a given query
returns: an index in sorted list matching the query (if any) else -1*/
int bisectionalSearchClass(struct classInfo , struct classInfo *, int , int );


/*for extending left or right from a hit (with index hitIndex); searches for extension left from ptr_range.start and up to 
hitIndex, and right from hitIndex and to ptr_range.end. In particular hitIndex must sit in the interval [ptr_range.start,
ptr_range.end] for this to work.
Input: leftRigth_b -- if 0 both extensions are done; if -1 only left, if 1 only right.
Transforms: updates ptr_range .start and .end to left resp right indexes found.*/
int bisectSingleExtendClass(int , struct classInfo , int , int *, struct classInfo *, int );


/*Heap sort of CATH domain list on domain name entry (if sortClass_b = 0) or on convClassId (sortClass_b =1):*/
void siftDownCATHDomain(struct chainInfo *, int , int , int , int );

void heapifyCATHDomain(struct chainInfo *, int , int , int );

void heapSortCATHDomain(struct chainInfo *, int , int , int );

/*search for a domain (query) in the list/array ptr_cathDomain of CATH domains. 
Returns: the hit index (-1 if no hit) */
int bisectionalSearchCATHDomain(char *, struct chainInfo *, int , int );



/*Functions for computing w-terms*/

/*function computing the w-value of a pair of segments, the latter input as 12 doubles*/
double w(
	/* coords of 1st point in 1st segment*/
	double,
	double,
	double,
	/* coords of 2nd point in 1st segment*/
	double,
	double,
	double,
	/* coords of 1st point in 2nd segment*/
	double,
	double,
	double,
	/* coords of 2nd point in 2nd segment*/
	double,
	double,
	double);


/*Function computing writhe contribution from two input line segments. Segments
are pairs of vectors as returned by.*/
double w_seg12(struct segment, struct segment);


/*function for computing all w-values across the simplex for a given sturcture*/
int wAll(struct segment *, int, double **);


/*Functions for carrying out the aggregations (recursions) across the simplex:*/

/*Aggregation/computation of invariants (including the abolute value versions) by recursion
through the simplex (traversing the invariants and the simplex). Takes w-values for the chain 
as input (wVal "double array"):*/
int aggr(int, int, int, double **, struct I_ptr);

/*Both omputation of w-terms and aggregation/computation of invariants excluding the abolute value versions;
the aggregation is done by recursion traversing the invariants and the simplex:*/
int aggrAndW_ExAbs(struct segment *, int, int, int, struct I_ptr);

/*Both computation of w-terms and aggregation/computation of invariants including the w-terms and the abolute 
value versions. The aggregation is done by recursion traversing the invariants and the simplex:*/
int aggrAndW(struct segment *, int , int, int, struct I_ptr );

/*Functions aimed for links/pokes searches:*/
int aggrAndW_wClosedLoops(struct segment *, int, int, int, struct I_ptr, struct twoSegmentIndex *, int, double, int *);

int examineClosedLoops(char *, struct twoSegmentIndex *, int, struct I_ptr, int, struct segment *, char *, struct cAlpha *, int, double, double);

int writheSubChainPairs(char *, struct segment *, int , struct I_ptr, int, char *, struct cAlpha *, int);

int genInvariantSubChainPairs(char *, struct segment *, int, struct I_ptr, double **, char *, int, int);

/*Functions aimed at alignment use*/
int setWindowCharacteristics(int , int *,  int , int *, int *, int *);

int getInvariantsOnWindows(struct I_windows_ptr *, int , int , int , int, int , struct I_ptr );


int getInvariantsOnWindows_1(struct I_windows_ptr *, int , int ,  int , struct I_ptr );


int writeInvariantsOnWindows(char *, struct I_windows_ptr , int , int );

//OBS: THIS ONLY WORKS/MAKES SENSE UP TO AND INCL ORDER 2:
int getInvariantsOnWindowPairs(struct I_windowPairs_ptr *, int , int , int ,  int , int, struct I_ptr );

int getInvariantsOnWindowPairs_1(struct I_windowPairs_ptr *, int , int ,  int , struct I_ptr );


//OBS: THIS ONLY WORKS/MAKES SENSE UP TO AND INCL ORDER 2:
int getInvariantsOnWindowPairsAdjWindowLength(struct I_windowPairs_ptr *, int , int ,  int , struct I_ptr );

int writeInvariantsOnwindowPairs(char *, struct I_windowPairs_ptr , int , int , int);


/*Functions for writing out the results to .txt files*/

int collectIvalues(struct I_values **, struct I_ptr, int, int, double, double);


int writeChain(char *, char *,  char *, char *, struct cAlpha *, int);

int writeIvaluesToFile(char *, struct I_values **, int, int);

int writeAllIvaluesToFile(char *, int, struct I_ptr, int);


/*Code blocks for making code in main function (below) more transparent; take care of memory 
allocation, initialization and freeing memeory:*/

int alloc_init_I_values(struct I_values ***, int, int);

int alloc_init_I_measures(struct I_ptr *, int, int, int, int, int, int, struct twoSegmentIndex ** );


int alloc_init_segment_ptr(struct Segment_ptr *, int, int);

int init_segment_ptr(struct Segment_ptr, struct segment *, int, int *);

int init_I_measures(struct I_ptr, int, int, int);


int alloc_init_I_windows_ptr(struct I_windows_ptr *, int , int , int , int , int, int, int );

int alloc_init_I_windowPairs_ptr(struct I_windowPairs_ptr *, int , int , int , int , int, int, int );


int alloc_init_ptr_Irarity_windows(struct Irarity_windows **, int , int , int , int , int, int, int );

int alloc_init_ptr_Irarity_windowPairs(struct Irarity_windowPairs **, int, int, int, int, int, int, int);

int alloc_init_Irarity_windowPairs_ptr(struct Irarity_windowPairs_ptr *, int , int , int , int , int );

int init_Irarity_windowPairs_ptr(struct Irarity_windowPairs_ptr , struct I_windowPairs_ptr *);


int alloc_init_chainInStr(struct chainInStructure *, int );

int reinit_chainInStr(struct chainInStructure , int );

/*
** functions for freeing memory allocated to pointers
*/
void freeDblArray(double **, size_t);

void free_I_measures(struct I_ptr, int, int, int);

void free_I_values(struct I_values **, int );



/*main program blocks*/
int computeGI(char[100] , char *, int, 
				int, char *, int, 
				char *, int , int , 
				int , int , int ,
				int , int, int , 
				int , int , int , 
				int , int , int , 
				int , int );


int computeGI_windows(char[100] , char *, int, int, char *, 
					  int, char *, int, int, int , 
					  int , int , int , int , int, 
					  int , int , int , int , int ,
					  int , int, int , int , int , 
					  int , int , int , int, int , int, 
					  int, char[2000], int , int );


#endif
