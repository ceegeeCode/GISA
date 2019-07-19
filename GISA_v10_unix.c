/* **********************************************************************************************************************************
*
*
* ***************************************   GISA_v10 *******************************************************************************
*
*
*
*     Present version: version _v10, > 21th June, 2019.
*     Previous version: version _v3 _v4, > 16th Jan, 2017; version _v5 > 5th Apr, 2017; version _v6, > 30th Oct, 2017; version _v7, > 17th nov, 2017; version _v8, > 17th Jan, 2018; version _v9, > 17th Sept, 2018.
*
*     Changed from v9: 
*     * Apparently an old version of the fct examineClosedLoops sat there; the right version of the fct has now been copied from GISA_v2
*     * Also from v2 the file name generation was copied (for the closed loop/inv's on sub-chain write outs)
*     New in version _v8: 
*     * Now possible to include all order 2 invariants in the matching (ie all 12 "full" order 2's, abs value versions incl'ed) 
*	  Changes from v2:
*     * structureName and classId added throughout for id of structures/domains and their classification id e.g. in SCOP or CATH.*
*
**************************************************************************************************************************************
* 
*     C-implementation of algorithm for computing Gauss integral based invariants for protein fold desciption. Code aimed 
*     for searching for "rare geometries" is included ("pokes"/"links"). Code for computing invariants on windows of a fixed length (subchains
*     of fixed length) is included; similarly for pairs of windows of a fixed length. The latter is (also) aimed for application to
*     structural alignment.
*
*     Author: Christian Grønbæk
*
*     If using this code cite this source:
*     C.Grønbæk., T.Hamelryck, P.Røgen: "GISA: Using Gauss integrals to identify rare conformations in 
*     protein structures", TO APPEAR.
*     
*     Notes:
*     !See ReadMe, code outlines and example runs on www.github.com, repository = ceeGeeCode/GISA! 
*     
*     In brief: The function to use for computing the invariants is computeGI, which uses global memory 
*     allocation (memory is allocated up front); the code allows running the computations across a 
*     set (directory) of PDB-files; memory is then allocated to cope with the longest chain in the set 
*	  (is done in chunks/batches of PPB-files, so as to keep memory consumption down by only allocating 
*     for longest chain in the chunk ). 
*     To use computeGI: run GISA_main_unix.o from command line; see examples and ReadMe. 
*
*     The code for computeGI is preceded by a description/documentation. An outline of the code 
*     can found oin the Github repo.
*     
*     The function computeGI_windows works quite as computeGI, only it derives from the computation of each GI
*     the values of that GI on every subchain of the set length of the surrounding structure. The same function
*     can be used to compute the GIs on every pair of windows (the value I[i][j] for the pair (i,j) will be either
*     the value of I-mutual on the pair; if i = j the value equals the I on the ith window). It is called from
*     a command line by executing SAonGISA_v10_uni.o in flavour makeDB (see material on github and/or use
*     the --help option).
*
*	  OBS: pointers are allocated memory with explicit cast for the sake of avoiding troubles with missing 
*     casts if code is augmented with e.g. CUDA C for gpu-computing. 

**********************************************************************************************/


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




/*struct for 3d point*/

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

int init_ptr_Irarity_windowPairs(struct Irarity_windowPairs **,  int );

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
				int , int);


int computeGI_windows(char[100] , char *, int, int, char *, 
					  int, char *, int, int, int , 
					  int , int , int , int , int, 
					  int , int , int , int , int ,
					  int , int, int , int , int , 
					  int , int , int , int, int , int,
					  int, char[2000], int , int );


/****************************************************
*****************Coding section *********************
*****************************************************/

//max and min are not part of standard C:
int min(int a, int b){
	
	if(a < b){ 
		return a;
	}
	else{
		return b;
	}
}

int max(int a, int b){

	if(a > b){ 
		return a;
	}
	else{
		return b;
	}
}

/*Util fct: input a filename with extension; returns (in in/out params): file name wo extension
and the extension*/
void splitFileName(char *fileNameWExt, char *fileName, char *ext){

	int extLen = 0, len = 0, lenFileName = 0, n = 0;
	char * fromDot ;

	char * token;
	const char dot[2] = "."; 
	
	/*first find the extension*/
    fromDot = strrchr(fileNameWExt, '.'); /*the part of fileNameWExt from the last dot*/
	printf("fromDot %s\n", fromDot);
	if(!fromDot){ 
		
		ext = "NN\0";

		/*extLen = strlen(ext);
		len = strlen(fileNameWExt) + extLen;*/

		strcpy(fileName, fileNameWExt); 

	}
	else{
		if(fromDot == fileNameWExt){/*fromDot == fileNameWExt happens if the file name has the form ".extension"*/

			printf("Fatal error: file name is blank, only the extension is present: %s\n",  fileNameWExt);
			return;
		}
		else{

			strcpy(ext, fromDot + 1);
			printf("ext %s\n", ext);

			/*extLen = strlen(ext);
			len = strlen(fileNameWExt);
			printf("len %d, ext len %d\n", len, extLen);

			n = len - extLen - 1;
			printf("n: %d", n);
			strncpy(fileName, fileNameWExt, n); */

			token = strtok(fileNameWExt, dot); //first token before dot
			strcpy(fileName, token);
			/*while( token != NULL ) { //ext wil catch the last token after the first one
				printf( " %s\n", token );
				strcpy(ext, token);

				token = strtok(NULL, dot);
			}*/

			printf("fileName %s\n", fileName);
			printf("fileNameWExt %s\n", fileNameWExt);
		
		}
	}

	lenFileName = strlen(fileName);
	fileName[lenFileName] = '\0';

}

/* Directory walk for Unix.
Parts inspired heavily by a post on http://stackoverflow.com by "lloydm".
Function for getting the contents of a named directory*/
struct dirContent ListDirectoryContents(char *sDir){
    
    DIR *dir;
    struct dirent *entry;

	int fileNumber = 0;
	int fileUnknownTypeNumber = 0;
	int i = 0;
	char * sPath; /*for making the code more readable; will point to address with value of fileName*/
	char * fileName; /*to contain the file name, without the path in front, that is.*/
	char * fileNameWExt;
	char * separator;
	char * ext;
	char ** ptr_dirList; /*pointer to contain path-filename list*/
	char ** ptr_fileNameList; /*pointer to contain filename list*/

	/*for holding output:*/
	struct dirContent output; 
	/*initialize:*/
	output.numberOfFiles = 0;
	output.ptr_dirList = NULL;
	output.ptr_fileNameList = NULL;

	sPath = (char *)calloc(pathLength, sizeof(char));
	fileName = (char *)calloc(fileNameLength, sizeof(char));
	fileNameWExt = (char *)calloc(fileNameLength, sizeof(char));
	ext = (char *)calloc(fileNameLength, sizeof(char));
	/*initialize these ptr's*/
	strcpy(fileName, "NN\0");
	strcpy(fileNameWExt, "NN\0");
	strcpy(ext, "NN\0");

	/*Check if directory exists, and, second, of it is non-empty*/
    if (!(dir = opendir(sDir))){
		printf("Sorry - could not open directory <%s>\n",sDir);
        return output;
	}
    if (!(entry = readdir(dir))){
		printf("Sorry -- the directory <%s >is apparently empty\n", sDir);
        return output;
	}

	/*First loop through the directory and count the number of files we want to consider*/
    do { /*While there's something to read, read the contents if it's a directory*/ 
        if (entry->d_type == DT_DIR) { 
            
			snprintf(sPath, sizeof(sPath)-1, "%s/%s", sDir, entry->d_name);
			printf("sPath: %s\n", sPath);
            
			if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0){
                continue;
			}

            //ListDirectoryContents(sPath); /*recurse*/
        
		}
		else{
			if (entry->d_type == DT_REG) { /*entry a regular file*/ 

				printf("regular file\n");
				fileNumber += 1;

			}
			else{
				if (entry->d_type == DT_UNKNOWN) { /*entry is file of nknown type*/ 

					printf("Obs: file %s has unknown type\n", entry->d_name);
					fileUnknownTypeNumber += 1;

				}
			}
		}

    } while (entry = readdir(dir));

	printf("Number of files of unknown type found in the directory: %d\n", fileUnknownTypeNumber);
	fileNumber += fileUnknownTypeNumber;
	printf("Total number of files found in the directory: %d\n", fileNumber);

	/*Loop through the directory again and collect the file names*/
	/*first allocate memory to our output pointer:*/
	ptr_dirList = (char **)malloc(fileNumber*pathLength*sizeof(char));
	ptr_fileNameList = (char **)malloc(fileNumber*fileNameLength*sizeof(char));

	/*init*/
	for (i= 0; i < fileNumber; i++){

			 ptr_dirList[i] = (char *)malloc(pathLength*sizeof(char));
			 ptr_fileNameList[i] = (char *)malloc(fileNameLength*sizeof(char));

			 strcpy(ptr_dirList[i], "NN\0");
			 strcpy(ptr_fileNameList[i], "NN\0");

	}


	//Reset.
	sPath = ptr_dirList[0]; /* let sPath use the first available address of ptr_dirList */
	fileName = ptr_fileNameList[0]; /* let fileName use the first available address of ptr_fileNameList */
	sprintf(sPath, "%s", sDir); //sPath now becomes the path sDir*

	printf("sPath: %s\n", sPath);

	/*Rewind to start from top:*/
	closedir(dir);
	dir = opendir(sDir);
	if (!(entry = readdir(dir))){
        return output;
	}

	fileNumber = 0;

	/*Loop through the directory again and record data in ptr*/
	do { /*While there's something to read, read the contents if it's a directory*/ 
        if (entry->d_type == DT_DIR) { 

			snprintf(sPath, sizeof(sPath)-1, "%s/%s", sDir, entry->d_name);
            
			if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0){
                continue;
			}

            //ListDirectoryContents(sPath); /*recurse*/
        
		}
		else{
			if (entry->d_type == DT_REG) { /*entry a regular file*/ 
				
				fileNumber += 1;
				sprintf(sPath, "%s/%s", sDir, entry->d_name);
				sprintf(fileNameWExt, "%s", entry->d_name);
				/*sscanf(fileNameWExt, "%s%[.]%3s",fileName, separator, ext);
				printf("fileName %s\n", fileName);
				printf("Ext %s\n", ext);*/
				splitFileName(fileNameWExt, fileName, ext);
				printf("fileName: %s has ext: %s\n", fileName, ext);

				/*record the names in output pointer*/
				sPath = ptr_dirList[fileNumber]; /*get sPath a new address: the next in line*/
				fileName = ptr_fileNameList[fileNumber]; /*get fileName a new address: the next in line*/
				//getchar();
			}
			else{
				if (entry->d_type == DT_UNKNOWN && strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0) { /*entry is file of unknown type*/ 

					//printf("Unknown file type\n");

					fileNumber += 1;
					sprintf(sPath, "%s/%s", sDir, entry->d_name);
					sprintf(fileNameWExt, "%s", entry->d_name);
					//printf("fileNameWExt %s\n", fileNameWExt);
					/*sscanf(fileNameWExt, "%s%[.]%3s",fileName, separator, ext);
					printf("fileName %s\n", fileName);
					printf("Ext %s\n", ext);*/
					splitFileName(fileNameWExt, fileName, ext);
					printf("fileName: %s ext: %s is seen as UNKNOWN type by readdir\n", fileName, ext);

					/*record the names in output pointer*/
					sPath = ptr_dirList[fileNumber]; /*get sPath a new address: the next in line*/
					fileName = ptr_fileNameList[fileNumber]; /*get fileName a new address: the next in line*/

				}
			}
		}


    } while (entry = readdir(dir));

	printf("dir walk ok\n");
	closedir(dir);

	output.numberOfFiles = fileNumber;
	output.ptr_dirList = ptr_dirList;
	output.ptr_fileNameList = ptr_fileNameList;

	return output;


}



int ListDirectoryContents2(char *sDir, struct dirContent *ptr_output, int *ptr_subDirCnt){
    
	int returnVal = 0;

	DIR *dir;
	struct dirent *entry;

	char * token;
	const char dot[2] = "."; 

	int fileNumber = 0;
	int fileUnknownTypeNumber = 0;
	int fileNumberStart;
	int fileNumberEnd;
	char * sPath; /*for making the code more readable; will point to address with value of fileName*/
	char * fileName; /*to contain the file name, without the path in front, that is.*/	
	char * fileNameWExt;
	char * separator;
	char * ext;
	char ** ptr_dirList; /*pointer to contain filename list*/
	char ** ptr_fileNameList; /*pointer to contain filename list*/
	int i = 0;
	int l = 0;

	/*initialize if no directories have been scanned yet:*/
	/*if(*ptr_subDirCnt == 0){
		ptr_output -> numberOfFiles = 0;
		ptr_output -> ptr_dirList = NULL;
		ptr_output -> ptr_fileNameList = NULL;
	}*/

	fileNumberStart = 0; //ptr_output -> numberOfFiles;
	fileNumberEnd = fileNumberStart;

	sPath = (char *)malloc(pathLength*sizeof(char));
	fileName = (char *)malloc(fileNameLength*sizeof(char));
	fileNameWExt = (char *)malloc(fileNameLength*sizeof(char));
	ext = (char *)malloc(fileNameLength*sizeof(char));
	/*initialize these ptr's*/
	strcpy(sPath, "NN\0");
	strcpy(fileName, "NN\0");
	strcpy(fileNameWExt, "NN\0");
	strcpy(ext, "NN\0");

	/*Check if directory can be found*/
	if (!(dir = opendir(sDir))){
		printf("Sorry - could not open directory <%s>\n",sDir);
        return 1;
	}
    if (!(entry = readdir(dir))){
		printf("Sorry -- the directory <%s >is apparently empty\n", sDir);
        return 1;
	} 


	sprintf(sPath, "%s//*.*", sDir); //first directory "is" sDir//*.* 


	/*The code here uses a recursion, which results in moving through a given tree of sub-directories.
	Upon execution the 1st do-while will either scan through the files if the directory holds nothing but 
	files, or, if there is a sub-directory, the control will enter another execution of the function itself
	for that sub-directory (ie a recursion). If that sub-directory contains only files it will scan through these
	or, if this sub-dir also contains a sub-dir that will call another execution of the function it self. And so
	on. The recursive calls will stop when entering a sub-sub-... -dir only containing files and, when done, the 
	control will resume in the "outer do-while" one step out. Here it will then seek any other sub-dirs or scan
	files if no other sub-dirs are remain to be scanned. Then done control will then move another step put (or back) to 
	that do-while ... and so on. In the final step the control returns to the do-while where the recursion was 
	initiated, moving on to scan any yet un-scanned sub-dirs or, finally, scanning any files in the dir itself
	(if no files are present in a (sub)dir, control will resume where it left, so move one step up; ultimately
	the outer if-clause "if(strcmp(fdFile.cFileName, ".") != 0 && strcmp(fdFile.cFileName, "..") != 0)" will terminate
	the exectuion).
	The code has two do-whiles: the first contains the recursion and control will move through the dir-tree as explained;
	in the 2nd the a ptr is (re)allocated for and fed with the ids (names) of the files in the current directory.*/ 
		
	/*First loop through the directory and count the number of files we want to consider*/
    do { /*While there's something to read, read the contents if it's a directory*/ 
        if (entry->d_type == DT_DIR) { 
            
			snprintf(sPath, sizeof(sPath)-1, "%s/%s", sDir, entry->d_name);
            
			if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0){
                continue;
			}

            printf("Count files in sub-directory: %s\n", sPath);
                ListDirectoryContents2(sPath, ptr_output, ptr_subDirCnt); /*recursion*/
        
		}
		else{
			if (entry->d_type == DT_REG) { /*entry a regular file*/ 

				fileNumber += 1;

			}
			else{
				if (entry->d_type == DT_UNKNOWN) { /*entry is file of unknown type*/ 

					printf("Obs: file type %s has unknown type\n", entry->d_name);
					fileUnknownTypeNumber += 1;

				}
			}
		}

    } while (entry = readdir(dir));


	printf("Number of files of unknown type found in this sub-directory: %d\n", fileUnknownTypeNumber);
	fileNumber += fileUnknownTypeNumber;
	printf("Total number of files found in this sub-directory: %d\n", fileNumber);


	if(fileNumber>0){

		fileNumberEnd += fileNumber;

		printf("SubdirCnt:%d\n",*ptr_subDirCnt);

		/*loop through the directory again and collect the file names*/
		/*first allocate memory to our output pointer:*/
		if(*ptr_subDirCnt ==0){
			ptr_output -> ptr_dirList = (char **)malloc(fileNumber*pathLength*sizeof(char));
			ptr_output -> ptr_fileNameList = (char **)malloc(fileNumber*fileNameLength*sizeof(char));
		}
		else{

			ptr_output -> ptr_dirList = (char **)realloc(ptr_output -> ptr_dirList, fileNumberEnd*pathLength*sizeof(char));
			ptr_output -> ptr_fileNameList = (char **)realloc(ptr_output -> ptr_fileNameList, fileNumberEnd*fileNameLength*sizeof(char));

		}
		/*init*/
		ptr_output -> numberOfFiles = fileNumberEnd;
		printf("I start at: %d and end at:%d\n",fileNumberStart, fileNumberEnd);
		for (i= fileNumberStart; i < fileNumberEnd; i++){

			/*ptr_output -> ptr_dirList[i] = (char *)calloc(pathLength, sizeof(char));
			ptr_output -> ptr_fileNameList[i] = (char *)calloc(fileNameLength, sizeof(char));*/

			ptr_output -> ptr_dirList[i] = (char *)malloc(pathLength*sizeof(char));
			ptr_output -> ptr_fileNameList[i] = (char *)malloc(fileNameLength*sizeof(char));


			strcpy(ptr_output -> ptr_dirList[i], "NN\0");
			strcpy(ptr_output -> ptr_fileNameList[i], "NN\0");

		}

		*ptr_subDirCnt +=1;

		printf("SubdirCnt:%d\n",*ptr_subDirCnt);
		
		/*Loop through the directory again and record data in ptr*/

		//Reset.
		sPath = ptr_output -> ptr_dirList[fileNumberStart]; /* let sPath use the first available address of ptr_dirList */
		fileName = ptr_output -> ptr_fileNameList[fileNumberStart]; /* let fileName use the first available address of ptr_fileNameList */
		sprintf(sPath, "%s//*.*", sDir); //sPath now becomes the directory sDir//*.*, where sDir is the directory lead to in the do-while above via the recursion!

		printf("sPath: %s\n", sPath);

		/*Rewind to start from top:*/
		closedir(dir);
		dir = opendir(sDir);
		if (!(entry = readdir(dir))){
			return 1;
		}

		fileNumber = fileNumberStart;

		/*Loop through the directory again and record data in ptr*/
		do { /*While there's something to read, read the contents if it's a directory*/ 
			if (entry->d_type == DT_DIR) { 

				snprintf(sPath, sizeof(sPath)-1, "%s/%s", sDir, entry->d_name);
            
				if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0){
					continue;
				}
        
			}
			else{
				if (entry->d_type == DT_REG) { /*entry a regular file*/ 
				
					fileNumber += 1;
					sprintf(sPath, "%s/%s", sDir, entry->d_name);
					sprintf(fileNameWExt, "%s", entry->d_name);
					splitFileName(fileNameWExt, fileName, ext);
					printf("fileName: %s ext: %s\n", fileName, ext);

					/*record the names in output pointer*/
					sPath = ptr_output -> ptr_dirList[fileNumber]; /*get sPath a new address: the next in line*/
					fileName = ptr_output -> ptr_fileNameList[fileNumber]; /*get fileName a new address: the next in line*/

				}
				else{
					if (entry->d_type == DT_UNKNOWN && strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0) { /*entry is file of nknown type*/ 

					//printf("Unknown file type\n");

					fileNumber += 1;
					sprintf(sPath, "%s/%s", sDir, entry->d_name);
					sprintf(fileNameWExt, "%s", entry->d_name);
					//printf("fileNameWExt %s\n", fileNameWExt);
					splitFileName(fileNameWExt, fileName, ext);
					printf("fileName: %s ext: %s is seen as UNKNOWN type by readdir\n", fileName, ext);

					/*record the names in output pointer*/
					sPath = ptr_output -> ptr_dirList[fileNumber]; /*get sPath a new address: the next in line*/
					fileName = ptr_output -> ptr_fileNameList[fileNumber]; /*get fileName a new address: the next in line*/
					
					}
				}
			}

		} while (entry = readdir(dir));

		closedir(dir);

	} /* end "if fileNumber > 0"

	/*if(fileNumber != fileNumberTotal){
		printf("Warning: less files recorded than found in directory");
	}*/

	/*
	for (i=0;i<10;i++){
		printf("check :%s ",ptr_dirList[i]);
	}
	*/


	return returnVal;

}


/*functions for reading content of a PDB-file and outputting the array of C-alpha coordinates
The first, readPDBChainStructure, gets the number of chains in the file (structure) and the length of 
each of the chains and returns this info in an array {...,(chainNr, chainLength), ...}. The 
second, main_readPDB, reads in the array of C-alpha coordinates in a given file (structure) 
along with chain length and chain number.*/
struct chainInStructure readPDBChainStructure(FILE * ptr_file){

		char buf[1000];
		char *currentChainId; /*to contain the id of the current chain when looping through the lines of the pdb file*/
		char chainIdInit[2] = ">"; /*some "unlikely" value for init of the chain identifier*/
		int  nrOfChains = 0;
		int  chainNr =0;
		int  chainIdIsBlank = 0; /*in some PDB file the chain id is blank; when this happens we set chainIdIsBlank = 1*/

		int chainLen = 0;
		int resNr = 0;
		int resNrPrior = -999;

		struct chainInfo *ptr_chainInStr; /*to contain the info -- number, length, id -- of the sub-chains making up the structure*/
		struct chainInStructure chainInStr;

		int i =0;

		currentChainId = (char *) malloc(sizeof(char)*2);

		/*We need to find the number of chains in the structure before we can allocate memory*/
		strcpy(currentChainId, chainIdInit);/*init to unlikely value*/
		while (fgets(buf,1000, ptr_file)!=NULL){
			if (strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && buf[21] != *currentChainId){ 
				chainNr += 1;
				/* if chain indicator (at position 21 in buf) is blank we record this:*/
				if (buf[21] == ' '){
					chainIdIsBlank = 1;
					exit;
				}
				//printf("buf at 21:%c\n", buf[21]);
				sscanf(buf+21, "%c", currentChainId);
				//printf("Chain Id:%c\n", *chainId);
			}
		}

		nrOfChains = chainNr;
		chainInStr.numberOfChains = nrOfChains;
		//printf("Nr of chains:%d\n", nrOfChains);

		/*Get the (chainNr, chainLength) structure of the file:*/

		/*allocate memory for chain structure pointer and initialize*/
		ptr_chainInStr = (struct chainInfo *) calloc (nrOfChains, sizeof(struct chainInfo));
		for(i = 0; i< nrOfChains; i++){
			ptr_chainInStr[i].chainId = (char *) malloc(sizeof(char)*2);
			strcpy(ptr_chainInStr[i].chainId, chainIdInit); /*init to unlikely value*/
			//printf("Chain Id initial value:%c\n", *ptr_chainInStr[i].chainId);
		}

		rewind(ptr_file);

		strcpy(currentChainId, chainIdInit); /*reset to unlikely value*/
		chainNr = -1;
		while (fgets(buf,1000, ptr_file)!=NULL){
			/*find the residue number*/
			sscanf(buf +23, "%d", &resNr);
			if (chainIdIsBlank == 1 && strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior || buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
				chainIdIsBlank = 2; /*just means that the id is blank, but has now been handled: loop never enters this if-clause again*/
				chainNr +=1;
				resNrPrior = resNr;
				chainLen =1; /*reset*/
				//printf("Chain id in file is blank, so the recorded value is id to initial value:%c\n", *ptr_chainInStr[chainNr].chainId);
				ptr_chainInStr[chainNr].chainLength = chainLen; /*for the event that the chain only has 1 residue! ...*/
			}
			else if (chainIdIsBlank == 0 && buf[21] != *currentChainId && strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior || buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
				sscanf(buf+21, "%c", currentChainId);
				//printf("Current chain Id in loop:%c\n", *currentChainId);
				chainNr +=1;
				resNrPrior = resNr;
				chainLen =1; /*reset*/
				sscanf(buf+21, "%c", ptr_chainInStr[chainNr].chainId);
				//printf("Chain id recorded:%c\n", *ptr_chainInStr[chainNr].chainId);
				ptr_chainInStr[chainNr].chainLength = chainLen; /*for the event that the chain only has 1 residue! ...*/
			}
			else if ((chainIdIsBlank == 2 || buf[21] == *currentChainId) && strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior || buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
				resNrPrior = resNr;
				chainLen +=1;
				ptr_chainInStr[chainNr].chainLength = chainLen;
			}
		}

		chainInStr.ptr_chainInStr = ptr_chainInStr;

		/*for (i=0;i<= nrOfChains-1; i++){
			printf("chain len of chain %d of id %c is:%d\n", i, *chainInStr.ptr_chainInStr[i].chainId, chainInStr.ptr_chainInStr[i].chainLength);
		}*/

    	return chainInStr;
	}

int readPDBChainStructure2(FILE * ptr_file, struct chainInStructure *ptr_chainInStr ){

		int returnVal =0;

		char buf[1000];
		char *currentStructureName;
		char *structureNameInit = "NN";
		char *currentClassId;
		char *classIdInit = "NN";
		char currentChainId[chainIdCharSize]; /*to contain the id of the current chain when looping through the lines of the pdb file*/
		char chainIdInit[2] = ">"; /*some "unlikely" value for init of the chain identifier*/
		int  nrOfChains = 0;
		int  chainNr =0;
		int  chainIdIsBlank = 0; /*in some PDB file the chain id is blank; when this happens we set chainIdIsBlank = 1*/

		int chainLen = 0;
		int resNr = 0;
		int resNrPrior = -999;

		struct chainInfo *ptr_chainInStr_temp; /*to contain the info -- number, length, id -- of the sub-chains making up the structure*/
		//struct chainInStructure chainInStr;

		int i,j =0;

		int isModel_b = 0; 

		//currentChainId = (char *) malloc(sizeof(char)*2);

		/*We need to find the number of chains in the structure before we can (assess whether to re)allocate more memory*/
		strcpy(currentChainId, chainIdInit);/*init to unlikely value*/
		while (fgets(buf,1000, ptr_file)!=NULL){

			/*if the structure is really a set of models we handle it separately*/
			if (strstr(buf, "MODEL") == buf){
				chainNr += 1;
				isModel_b = 1;
			}
			/*and else we scan through the file and try capturing the nr of chains; if the chain id is left 
			blank is the file, only one chain is allowed*/
			else if(isModel_b == 0){

				if ((strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf) && buf[13] == 'C' && buf[14] == 'A' && buf[21] != *currentChainId){ 
					chainNr += 1;
					/* if chain indicator (at position 21 in buf) is blank we record this:*/
					if (buf[21] == ' '){
						chainIdIsBlank = 1;
						exit;
					}
					//printf("buf at 21:%c\n", buf[21]);
					sscanf(buf+21, "%c", currentChainId);
					//printf("ChainNr: %d Chain Id:%c\n", chainNr, *currentChainId);
				}

			}
		}

				
		nrOfChains = chainNr;
		ptr_chainInStr -> numberOfChains = nrOfChains;
		//printf("Nr of chains:%d\n", nrOfChains);
		//getchar();
		ptr_chainInStr_temp = ptr_chainInStr -> ptr_chainInStr;

		/*if the nr of chains is larger than what was globally allocated for (to *ptr_ptr_chainInStr)
		we reallocate some more, and initialize too:*/
		if(nrOfChains > maxNrOfChains){
			/*reallocate memory to chain structure pointer and initialize*/
			ptr_chainInStr_temp = (struct chainInfo *) realloc (ptr_chainInStr ->ptr_chainInStr, nrOfChains*sizeof(struct chainInfo));
			for(i = 0; i< nrOfChains; i++){
				ptr_chainInStr_temp[i].chainId = (char *) malloc(sizeof(char)*chainIdCharSize);
				ptr_chainInStr_temp[i].structureName = (char *) malloc(sizeof(char)*structureNameCharSize);
				ptr_chainInStr_temp[i].classId = (char *) malloc(sizeof(char)*classIdCharSize);
				ptr_chainInStr_temp[i].classIdConv = (int *) malloc(classifierSize*sizeof(int));
				/*init*/
				strcpy(ptr_chainInStr_temp[i].chainId, chainIdInit); /*init to unlikely value*/
				strcpy(ptr_chainInStr_temp[i].structureName, structureNameInit);
				strcpy(ptr_chainInStr_temp[i].classId, classIdInit);
				//printf("Chain Id initial value:%c\n", *ptr_chainInStr[i].chainId);
				
				for(j=0; j< classifierSize; j++){ptr_chainInStr_temp[i].classIdConv[j] = -1;} //to be updated later, if at all


			}
		}

		rewind(ptr_file);

		/*Get the (chainNr, chainLength) structure of the file:*/
		/*Again we handle the not-/model cases separately:*/
		if(isModel_b == 0){

			strcpy(currentChainId, chainIdInit); /*reset to unlikely value*/
			chainNr = -1;
			while (fgets(buf,1000, ptr_file)!=NULL){
				/*find the residue number*/
				sscanf(buf +23, "%d", &resNr);
				if (chainIdIsBlank == 1 && (strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf)  && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior|| buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
					//printf("(saaledes stadig ganske hemmelig), at start of chain: %d resNr: %d prior:%d\n", chainNr, resNr, resNrPrior);
					//getchar();
					chainIdIsBlank = 2; /*just means that the id is blank, but has now been handled: loop never enters this if-clause again*/
					chainNr +=1;
					resNrPrior = resNr;
					chainLen =1; /*reset*/
					//printf("Chain id in file is blank, so the recorded value is id to initial value:%c\n", *ptr_chainInStr[chainNr].chainId);
					ptr_chainInStr_temp[chainNr].chainLength = chainLen; 
					strcpy(ptr_chainInStr_temp[chainNr].chainId, currentChainId);
				}
				else if (chainIdIsBlank == 0 && buf[21] != *currentChainId && (strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf)  && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior|| buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
					//printf("At start of chain:%d , (langmodigt venter ...), resNr: %d prior:%d\n", chainNr, resNr, resNrPrior);
					//getchar();
					sscanf(buf+21, "%c", currentChainId);
					//printf("Current chain Id in loop:%c\n", *currentChainId);
					chainNr +=1;
					resNrPrior = resNr;
					chainLen =1; /*reset*/
					//sscanf(buf+21, "%c", ptr_chainInStr_temp[chainNr].chainId);
					//printf("Chain id recorded:%c\n", *ptr_chainInStr[chainNr].chainId);
					ptr_chainInStr_temp[chainNr].chainLength = chainLen; 
					strcpy(ptr_chainInStr_temp[chainNr].chainId, currentChainId);
				}
				//else if ((chainIdIsBlank == 2 || buf[21] == *currentChainId) && strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior || buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
				//	resNrPrior = resNr;
				//	chainLen +=1;
				//	ptr_chainInStr_temp[chainNr].chainLength = chainLen;
				//}
				/*NOT ... CHANGED TO "resNr >= resNrPrior+1" ON 14TH JUNE 2019; HOLES ARE NOW HANDLED BY REMOVING CHAINS CONTAINING TOO LONG SEGMENTS: The "resNr == resNrPrior+1" serves to avoid holes in structures*/
				else if ((chainIdIsBlank == 2 || buf[21] == *currentChainId) && (strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf)  && buf[13] == 'C' && buf[14] == 'A' && (resNr == resNrPrior +1|| buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
					//printf("Inside chain: %d, (Sommeren mild ...), resNr: %d prior:%d\n", chainNr, resNr, resNrPrior);
					//getchar();
					resNrPrior = resNr;
					chainLen +=1;
					//
					ptr_chainInStr_temp[chainNr].chainLength = chainLen;
					//printf("In chain: Current chain Id in loop:%c ResNr: %d chain Len\n", *currentChainId, resNr, chainLen);
				}
				/*if the chain has a hole and the chain recorded so far is below a global cut-off we reset the chain to allow for the case where the records to follow consitute a proper chain:*/  
				else if ((chainIdIsBlank == 2 || buf[21] == *currentChainId) && chainLen < strLengthCutOff && (strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf)  && buf[13] == 'C' && buf[14] == 'A' &&  (resNr != resNrPrior +1 && buf[26] != 'A' && buf[26] != 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
					//printf("Inside chain: %d, (saaledes stadig ganske hemmelig), resNr: %d prior:%d .. resetting chain due to hole\n", chainNr, resNr, resNrPrior);
					//getchar();
					resNrPrior = resNr;
					chainLen =1;
					//
					ptr_chainInStr_temp[chainNr].chainLength = chainLen;
					//printf("In chain: Current chain Id in loop:%c ResNr: %d chain Len\n", *currentChainId, resNr, chainLen);
				}
				/*record the characteristics found when arriving at a TER or an END; chainLen> 0 is there to filter away cases of "empty" domains*/
				else if (chainLen > 0 && strstr(buf, "TER") == buf){
					//printf("The buf, (Langmodigt ...):%s chain nr:%d chainId:%s\n", buf, chainNr, currentChainId);
					//strcpy(ptr_chainInStr_temp[chainNr].chainId, currentChainId);
					//printf("Chain id recorded:%c\n", *ptr_chainInStr_temp[chainNr].chainId);
					ptr_chainInStr_temp[chainNr].chainLength = chainLen;
					resNrPrior = -999; /*reset*/
				}
			}
		}
		else if(isModel_b == 1){ 

			chainNr = -1;
			while (fgets(buf,1000, ptr_file)!=NULL){
				
				/*we use the running chainNr as chainId*/
				if (strstr(buf, "MODEL") == buf){
					chainNr += 1;
					chainLen = 0;/*reset*/
					snprintf(currentChainId,10,"%d",chainNr);
					//printf("Langmodigt venter bysvalen .. model case; chainNr: %d, currentChainId:  %s\n", chainNr, currentChainId);
					//getchar();
				}
				else{

					/*find the residue number*/
					sscanf(buf +23, "%d", &resNr);
					/*record the first residue nr in the model/chain:*/
					if(chainLen==0){
						resNrPrior = resNr -1; /*technical init*/
						//printf("(Sommeren mild oedsles ..) .. model case; init, the resNr is %d and the prior is %d\n", resNr, resNrPrior);
						//getchar();
					}

					/*For MODELs we soften the "resNr == resNrPrior+1" (which serves to avoid holes in structures) to just resNr != resNrPrior*/
					/*We also allow HETATOMs in the chains*/
					if ((strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf)  && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior || buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
						//printf("(saaledes stadig ganske hemmelig) .. model case, inside chain; the resNr is %d and the prior is %d\n", resNr, resNrPrior);
						//getchar();
						resNrPrior = resNr;
						chainLen +=1; 
					}
					/*record the characteristics when arriving at a TER*/
					if (strstr(buf, "TER") == buf){
						//printf("Langmodigt venter bysvalen .. model case, at TER; chainNr: %d, currentChainId:  %s\n", chainNr, currentChainId);
						//getchar();
						strcpy(ptr_chainInStr_temp[chainNr].chainId, currentChainId);
						//printf("Chain id recorded:%c\n", *ptr_chainInStr[chainNr].chainId);
						ptr_chainInStr_temp[chainNr].chainLength = chainLen; 
						resNrPrior = -999; //reset
					}
			
				}

			} 

		}
			



		ptr_chainInStr -> ptr_chainInStr = ptr_chainInStr_temp;

		/* for (i=0;i<= nrOfChains-1; i++){
			printf("chain len of chain %d of id %c is:%d\n", i, ptr_chainInStr_temp[i].chainId, ptr_chainInStr_temp[i].chainLength);
		}*/

    	return returnVal;
	}

int readPDBChainStructureSCOP(FILE * ptr_file, struct chainInStructure *ptr_chainInStr ){

		int returnVal =0;

		char buf[1000];
		char *currentStructureName;
		char *structureNameInit = "NN";
		char *currentClassId;
		char *classIdInit = "NN";
		char currentChainId[2]; /*to contain the id of the current chain when looping through the lines of the pdb file*/
		char chainIdInit[2] = ">"; /*some "unlikely" value for init of the chain identifier*/
		int  nrOfChains = 0;
		int  chainNr =0;
		int  chainIdIsBlank = 0; /*in some PDB file the chain id is blank; when this happens we set chainIdIsBlank = 1*/

		int chainLen = 0;
		int resNr = 0;
		int resNrPrior = -999;

		struct chainInfo *ptr_chainInStr_temp; /*to contain the info -- number, length, id -- of the sub-chains making up the structure*/
		//struct chainInStructure chainInStr;

		int i,j =0;

		//currentChainId = (char *) malloc(sizeof(char)*2);

		/*We need to find the number of chains in the structure before we can (assess whether to re)allocate more memory*/
		strcpy(currentChainId, chainIdInit);/*init to unlikely value*/
		while (fgets(buf,1000, ptr_file)!=NULL){

			if (strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && buf[21] != *currentChainId){ 
				chainNr += 1;
				/* if chain indicator (at position 21 in buf) is blank we record this:*/
				if (buf[21] == ' '){
					chainIdIsBlank = 1;
					exit;
				}
				//printf("buf at 21:%c\n", buf[21]);
				sscanf(buf+21, "%c", currentChainId);
				//printf("Chain Id:%c\n", *chainId);
			}
		}

		nrOfChains = chainNr;
		ptr_chainInStr -> numberOfChains = nrOfChains;
		//printf("Nr of chains:%d\n", nrOfChains);
		ptr_chainInStr_temp = ptr_chainInStr -> ptr_chainInStr;

		/*Get the (chainNr, chainLength) structure of the file:*/

		/*if the nr of chains is larger than what was globally allocated for (to *ptr_ptr_chainInStr)
		we reallocate some more, and initialize too:*/
		if(nrOfChains > maxNrOfChains){
			/*reallocate memory to chain structure pointer and initialize*/
			ptr_chainInStr_temp = (struct chainInfo *) realloc (ptr_chainInStr ->ptr_chainInStr, nrOfChains*sizeof(struct chainInfo));
			for(i = 0; i< nrOfChains; i++){
				ptr_chainInStr_temp[i].chainId = (char *) malloc(sizeof(char)*chainIdCharSize);
				ptr_chainInStr_temp[i].structureName = (char *) malloc(sizeof(char)*structureNameCharSize);
				ptr_chainInStr_temp[i].classId = (char *) malloc(sizeof(char)*classIdCharSize);
				ptr_chainInStr_temp[i].classIdConv = (int *) malloc(classifierSize*sizeof(int));
				/*init*/
				strcpy(ptr_chainInStr_temp[i].chainId, chainIdInit); /*init to unlikely value*/
				strcpy(ptr_chainInStr_temp[i].structureName, structureNameInit);
				strcpy(ptr_chainInStr_temp[i].classId, classIdInit);
				//printf("Chain Id initial value:%c\n", *ptr_chainInStr[i].chainId);
				
				for(j=0; j< classifierSize; j++){ptr_chainInStr_temp[i].classIdConv[j] = -1;} //to be updated later, if at all


			}
		}

		rewind(ptr_file);

		strcpy(currentChainId, chainIdInit); /*reset to unlikely value*/
		chainNr = -1;
		while (fgets(buf,1000, ptr_file)!=NULL){
			/*find the residue number (if any)*/
			sscanf(buf +23, "%d", &resNr);

			if (strstr(buf, "HEADER") == buf){
				chainNr +=1;
				sscanf(buf+30, "%s", ptr_chainInStr_temp[chainNr].structureName);
			}
			else if (strstr(buf, "REMARK") == buf && strstr(buf, "sccs") ){ //&& strstr(buf + 25, "sccs") == buf + 25
				sscanf(buf+ 30, "%s.%s", ptr_chainInStr_temp[chainNr].classId);
			}
			else if (chainIdIsBlank == 1 && strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior || buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
				chainIdIsBlank = 2; /*just means that the id is blank, but has now been handled: loop never enters this if-clause again*/
				//chainNr +=1;
				resNrPrior = resNr;
				chainLen =1; /*reset*/
				//printf("Chain id in file is blank, so the recorded value is id to initial value:%c\n", *ptr_chainInStr[chainNr].chainId);
				ptr_chainInStr_temp[chainNr].chainLength = chainLen; /*for the event that the chain only has 1 residue! ...*/
			}
			else if (chainIdIsBlank == 0 && buf[21] != *currentChainId && strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior || buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
				sscanf(buf+21, "%c", currentChainId);
				//printf("Current chain Id in loop:%c\n", *currentChainId);
				//chainNr +=1;
				resNrPrior = resNr;
				chainLen =1; /*reset*/
				sscanf(buf+21, "%c", ptr_chainInStr_temp[chainNr].chainId);
				//printf("Chain id recorded:%c\n", *ptr_chainInStr[chainNr].chainId);
				ptr_chainInStr_temp[chainNr].chainLength = chainLen; /*for the event that the chain only has 1 residue! ...*/
			}
			else if ((chainIdIsBlank == 2 || buf[21] == *currentChainId) && strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior || buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
				resNrPrior = resNr;
				chainLen +=1;
				ptr_chainInStr_temp[chainNr].chainLength = chainLen;
			}
		}

		ptr_chainInStr -> ptr_chainInStr = ptr_chainInStr_temp;

		/*for (i=0;i<= nrOfChains-1; i++){
			printf("chain len of chain %d of id %c is:%d\n", i, *chainInStr.ptr_chainInStr[i].chainId, chainInStr.ptr_chainInStr[i].chainLength);
		}*/

    	return returnVal;
	}

/*Currently id to readPDBChainStructure2:*/
int readPDBDomainStructureCATH(FILE * ptr_file, struct chainInStructure *ptr_chainInStr ){

		int returnVal =0;

		char buf[1000];
		char *currentStructureName;
		char *structureNameInit = "NN";
		char *currentClassId;
		char *classIdInit = "NN";
		char currentChainId[chainIdCharSize]; /*to contain the id of the current chain when looping through the lines of the pdb file*/
		char chainIdInit[2] = ">"; /*some "unlikely" value for init of the chain identifier*/
		int  nrOfChains = 0;
		int  chainNr =0;
		int  chainIdIsBlank = 0; /*in some PDB file the chain id is blank; when this happens we set chainIdIsBlank = 1*/

		int chainLen = 0;
		int resNr = 0;
		int resNrPrior = -999;

		struct chainInfo *ptr_chainInStr_temp; /*to contain the info -- number, length, id -- of the sub-chains making up the structure*/
		//struct chainInStructure chainInStr;

		int i,j =0;

		int isModel_b = 0; 

		//currentChainId = (char *) malloc(sizeof(char)*2);

		/*We need to find the number of chains in the structure before we can (assess whether to re)allocate more memory*/
		strcpy(currentChainId, chainIdInit);/*init to unlikely value*/
		while (fgets(buf,1000, ptr_file)!=NULL){

			/*if the structure is really a set of models we handle it separately*/
			if (strstr(buf, "MODEL") == buf){
				chainNr += 1;
				isModel_b = 1;
			}
			/*and else we scan through the file and try capturing the nr of chains; if the chain id is left 
			blank is the file, only one chain is allowed*/
			else if(isModel_b == 0){

				if (strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && buf[21] != *currentChainId){ 
					chainNr += 1;
					/* if chain indicator (at position 21 in buf) is blank we record this:*/
					if (buf[21] == ' '){
						chainIdIsBlank = 1;
						exit;
					}
					//printf("buf at 21:%c\n", buf[21]);
					sscanf(buf+21, "%c", currentChainId);
					//printf("Chain Id:%c\n", *chainId);
				}
			}
		}

				
		nrOfChains = chainNr; /*if the structure is a set of models, this will be the number of models included*/
		ptr_chainInStr -> numberOfChains = nrOfChains;
		//printf("Nr of chains:%d\n", nrOfChains);
		ptr_chainInStr_temp = ptr_chainInStr -> ptr_chainInStr;

		/*if the nr of chains is larger than what was globally allocated for (to *ptr_ptr_chainInStr)
		we reallocate some more, and initialize too:*/
		if(nrOfChains > maxNrOfChains){

			/*reallocate memory to chain structure pointer and initialize*/
			ptr_chainInStr_temp = (struct chainInfo *) realloc (ptr_chainInStr ->ptr_chainInStr, nrOfChains*sizeof(struct chainInfo));
			for(i = 0; i< nrOfChains; i++){
				ptr_chainInStr_temp[i].chainId = (char *) malloc(sizeof(char)*chainIdCharSize);
				ptr_chainInStr_temp[i].structureName = (char *) malloc(sizeof(char)*structureNameCharSize);
				ptr_chainInStr_temp[i].classId = (char *) malloc(sizeof(char)*classIdCharSize);
				ptr_chainInStr_temp[i].classIdConv = (int *) malloc(classifierSize*sizeof(int));
				/*init*/
				strcpy(ptr_chainInStr_temp[i].chainId, chainIdInit); /*init to unlikely value*/
				strcpy(ptr_chainInStr_temp[i].structureName, structureNameInit);
				strcpy(ptr_chainInStr_temp[i].classId, classIdInit);
				//printf("Chain Id initial value:%c\n", *ptr_chainInStr[i].chainId);
				ptr_chainInStr_temp[i].chainLength = 0;

				
				for(j=0; j< classifierSize; j++){ptr_chainInStr_temp[i].classIdConv[j] = -1;} //to be updated later, if at all


			}
		}

		rewind(ptr_file);

		/*Get the (chainNr, chainLength) structure of the file:*/
		/*Again we handle the not-/model cases separately:*/
		if(isModel_b == 0){

			strcpy(currentChainId, chainIdInit); /*reset to unlikely value*/
			chainNr = -1;
			while (fgets(buf,1000, ptr_file)!=NULL){
				/*find the residue number*/
				sscanf(buf +23, "%d", &resNr);

				if (chainIdIsBlank == 1 && (strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf) && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior || buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
					chainIdIsBlank = 2; /*just means that the id is blank, but has now been handled: loop never enters this if-clause again*/
					chainNr +=1;
					resNrPrior = resNr;
					chainLen = 1; /*reset*/
					//printf("Chain id in file is blank, so the recorded value is id to initial value:%c resNr:%d\n", *ptr_chainInStr_temp[chainNr].chainId, resNr);
					ptr_chainInStr_temp[chainNr].chainLength = chainLen; 
					strcpy(ptr_chainInStr_temp[chainNr].chainId, currentChainId);
				}
				else if (chainIdIsBlank == 0 && buf[21] != *currentChainId && (strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf) && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior|| buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
					sscanf(buf+21, "%c", currentChainId);
					//printf("At start of chain: Current chain Id in loop:%c ResNr: %d \n", *currentChainId, resNr, resNrPrior);
					//getchar();
					chainNr +=1;
					resNrPrior = resNr;
					chainLen =1; /*reset*/
					//sscanf(buf+21, "%c", ptr_chainInStr_temp[chainNr].chainId);
					//printf("Chain id recorded:%c\n", *ptr_chainInStr[chainNr].chainId);
					ptr_chainInStr_temp[chainNr].chainLength = chainLen;
					strcpy(ptr_chainInStr_temp[chainNr].chainId, currentChainId);
				}
				/*The "resNr == resNrPrior+1" serves to avoid holes in structures*/
				else if ((chainIdIsBlank == 2 || buf[21] == *currentChainId) && (strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf) && buf[13] == 'C' && buf[14] == 'A' && (resNr == resNrPrior +1|| buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
					//printf("Inside chain: Current chain Id in loop:%c ResNr: %d Prior resNr: %d\n", *currentChainId, resNr, resNrPrior);
					//getchar();
					resNrPrior = resNr;
					chainLen +=1;
					ptr_chainInStr_temp[chainNr].chainLength = chainLen;
					//printf("In chain: Current chain Id in loop:%c ResNr: %d chain Len\n", *currentChainId, resNr, chainLen);
				}
				/*if the chain has a hole and the chain recorded so far is below a global cut-off we reset the chain to allow for the case where the records to follow consitute a proper chain:*/  
				else if ((chainIdIsBlank == 2 || buf[21] == *currentChainId) && chainLen < strLengthCutOff && (strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf)  && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior +1 && buf[26] != 'A' && buf[26] != 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
					//printf("Inside chain: %d, Sommeren mild ... , resNr: %d prior:%d .. resetting chain due to hole\n", chainNr, resNr, resNrPrior);
					//getchar();
					resNrPrior = resNr;
					chainLen =1;
					//
					ptr_chainInStr_temp[chainNr].chainLength = chainLen;
					//printf("In chain: Current chain Id in loop:%c ResNr: %d chain Len\n", *currentChainId, resNr, chainLen);
				}
				/*record the characteristics found when arriving at a TER or an END; chainLen> 0 is there to filter away cases of "empty" domains*/
				else if (chainLen > 0 && (strstr(buf, "TER") == buf||strstr(buf, "END") == buf||strstr(buf, "END OF SEGMENT OF DOMAIN") == buf)){
					//printf("The buf, (saaledes stadig ganske hemmelig):%s chain nr:%d chainId:%s\n", buf, chainNr, currentChainId);
					//strcpy(ptr_chainInStr_temp[chainNr].chainId, currentChainId);
					//printf("Chain id recorded:%c\n", *ptr_chainInStr_temp[chainNr].chainId);
					ptr_chainInStr_temp[chainNr].chainLength = chainLen;
					resNrPrior = -999; /*reset*/
				}
			

			
			} 

		}/*end isModel_b ==0*/
		else if(isModel_b == 1){ 

			chainNr = -1;
			while (fgets(buf,1000, ptr_file)!=NULL){
				
				/*we use the running chainNr as chainId*/
				if (strstr(buf, "MODEL") == buf){
					chainNr += 1;
					chainLen = 0;/*reset*/
					snprintf(currentChainId,10,"%d",chainNr);
				}
				else{

					/*find the residue number*/
					sscanf(buf +23, "%d", &resNr);
					/*record the first residue nr in the model/chain:*/
					if(chainLen==0){
						resNrPrior = resNr -1; /*technical init*/
						//printf("Langmodigt venter bysvalen .. model case; init, the resNr is %d and the prior is %d\n", resNr, resNrPrior);
					}

					/*For MODELs we soften the "resNr == resNrPrior+1" (which serves to avoid holes in structures) to just resNr != resNrPrior*/
					/*We also allow HETATOMs in the chains*/
					if ((strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf)  && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior || buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
						resNrPrior = resNr;
						chainLen +=1; 
					}
					/*record the characteristics when arriving at a TER*/
					if (strstr(buf, "TER") == buf){
						strcpy(ptr_chainInStr_temp[chainNr].chainId, currentChainId);
						//printf("Chain id recorded:%c\n", *ptr_chainInStr[chainNr].chainId);
						ptr_chainInStr_temp[chainNr].chainLength = chainLen; 
						resNrPrior = -999; //reset
					}
			
				}

			} 

		}

		ptr_chainInStr -> ptr_chainInStr = ptr_chainInStr_temp;

		/*for (i=0;i<= nrOfChains-1; i++){
			printf("chain len of chain %d of id %c is:%d\n", i, *chainInStr.ptr_chainInStr[i].chainId, chainInStr.ptr_chainInStr[i].chainLength);
		}*/

		//getchar();

    	return returnVal;
	}

int readPDBDomainStructureCATH_old(FILE * ptr_file, struct chainInStructure *ptr_chainInStr ){

		int returnVal =0;

		char buf[1000];
		//char *currentStructureName;
		//char *structureNameInit = "NN";
		char *currentClassId;
		char *classIdInit = "NN";
		char currentChainId[2]; /*to contain the id of the current chain when looping through the lines of the pdb file*/
		char chainIdInit[2] = ">"; /*some "unlikely" value for init of the chain identifier*/
		int  nrOfChains = 0;
		int  chainNr =0; //static
		int  chainIdIsBlank = 0; /*in some PDB file the chain id is blank; when this happens we set chainIdIsBlank = 1*/

		int chainLen = 0;
		int resNr = 0;
		int resNrPrior = -999;

		struct chainInfo *ptr_chainInStr_temp; /*to contain the info -- number, length, id -- of the sub-chains making up the structure*/
		//struct chainInStructure chainInStr;

		int i,j =0;

		//currentChainId = (char *) malloc(sizeof(char)*2);

		/*Since CATH domains are pieces of a single PDB chain, we know the nr of chains (=1) and, also,  
		need not allocate more memory*/
				
		nrOfChains = 1;
		ptr_chainInStr -> numberOfChains = nrOfChains;
		//printf("Nr of chains:%d\n", nrOfChains);
		ptr_chainInStr_temp = ptr_chainInStr -> ptr_chainInStr;

		/*Get the (chainNr, chainLength) structure of the file:*/

		strcpy(currentChainId, chainIdInit); /*reset to unlikely value*/
		while (fgets(buf,1000, ptr_file)!=NULL){
			/*find the residue number*/
			sscanf(buf +23, "%d", &resNr);
			//printf("resNr: %d\n", resNr);
			if (strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior || buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
				resNrPrior = resNr;
				/* if chain indicator (at position 21 in buf) isn't blank we record this as chain id (else it becomes the 
				default, ie chainIdInit:*/
				if (buf[21] != ' '){
					sscanf(buf+21, "%c", currentChainId);
				}
				chainLen +=1;
				//printf("Chain len: %d\n", chainLen);
			}
		}
		/*record the values (last read in)*/
		strcpy(ptr_chainInStr_temp[chainNr].chainId, currentChainId);
		ptr_chainInStr_temp[chainNr].chainLength = chainLen;

		ptr_chainInStr -> ptr_chainInStr = ptr_chainInStr_temp;

		/*for (i=0;i<= nrOfChains-1; i++){
			printf("chain len of chain %d of id %c is:%d\n", i, *chainInStr.ptr_chainInStr[i].chainId, chainInStr.ptr_chainInStr[i].chainLength);
		}*/

    	return returnVal;
	}

int readPDBHeaderSCOP(FILE * ptr_file, struct chainInStructure *ptr_chainInStr ){

		int returnVal =0;

		char buf[1000];
		char *currentStructureName;
		char *structureNameInit = "NN";
		char *currentClassId;
		char *classIdInit = "NN";
		char currentChainId[2]; /*to contain the id of the current chain when looping through the lines of the pdb file*/
		char chainIdInit[2] = ">"; /*some "unlikely" value for init of the chain identifier*/
		int  nrOfChains = 0;
		int  chainNr =0;
		int  chainIdIsBlank = 0; /*in some PDB file the chain id is blank; when this happens we set chainIdIsBlank = 1*/

		int chainLen = 0;
		int resNr = 0;
		int resNrPrior = -999;

		struct chainInfo *ptr_chainInStr_temp; /*to contain the info -- number, length, id -- of the sub-chains making up the structure*/
		//struct chainInStructure chainInStr;

		int i =0;

		//currentChainId = (char *) malloc(sizeof(char)*2);

		/*We need to find the number of chains in the structure before we can (assess whether to re)allocate more memory*/
		strcpy(currentChainId, chainIdInit);/*init to unlikely value*/
		while (fgets(buf,1000, ptr_file)!=NULL){

			if (strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && buf[21] != *currentChainId){ 
				chainNr += 1;
				/* if chain indicator (at position 21 in buf) is blank we record this:*/
				if (buf[21] == ' '){
					chainIdIsBlank = 1;
					exit;
				}
				//printf("buf at 21:%c\n", buf[21]);
				sscanf(buf+21, "%c", currentChainId);
				//printf("Chain Id:%c\n", *chainId);
			}
		}

		nrOfChains = chainNr;
		ptr_chainInStr -> numberOfChains = nrOfChains;
		//printf("Nr of chains:%d\n", nrOfChains);
		ptr_chainInStr_temp = ptr_chainInStr -> ptr_chainInStr;

		/*Get the (chainNr, chainLength) structure of the file:*/

		/*if the nr of chains is larger than what was globally allocated for (to *ptr_ptr_chainInStr)
		we reallocate some more, and initialize too:*/
		if(nrOfChains > maxNrOfChains){
			/*reallocate memory to chain structure pointer and initialize*/
			ptr_chainInStr_temp = (struct chainInfo *) realloc (ptr_chainInStr ->ptr_chainInStr, nrOfChains*sizeof(struct chainInfo));
			for(i = 0; i< nrOfChains; i++){
				ptr_chainInStr_temp[i].chainId = (char *) malloc(sizeof(char)*chainIdCharSize);
				ptr_chainInStr_temp[i].structureName = (char *) malloc(sizeof(char)*structureNameCharSize);
				ptr_chainInStr_temp[i].classId = (char *) malloc(sizeof(char)*classIdCharSize);
				/*init*/
				strcpy(ptr_chainInStr_temp[i].chainId, chainIdInit); /*init to unlikely value*/
				strcpy(ptr_chainInStr_temp[i].structureName, structureNameInit);
				strcpy(ptr_chainInStr_temp[i].classId, classIdInit);
				//printf("Chain Id initial value:%c\n", *ptr_chainInStr[i].chainId);
			}
		}

		rewind(ptr_file);

		strcpy(currentChainId, chainIdInit); /*reset to unlikely value*/
		chainNr = -1;
		while (fgets(buf,1000, ptr_file)!=NULL){
			/*find the residue number (if any)*/
			sscanf(buf +23, "%d", &resNr);

			if (strstr(buf, "HEADER") == buf){
				chainNr +=1;
				sscanf(buf+30, "%s", ptr_chainInStr_temp[chainNr].structureName);
			}
			else if (strstr(buf, "REMARK") == buf && strstr(buf, "sccs") ){ //&& strstr(buf + 25, "sccs") == buf + 25
				sscanf(buf+ 30, "%s", ptr_chainInStr_temp[chainNr].classId);
			}
			else if (chainIdIsBlank == 1 && strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior || buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
				chainIdIsBlank = 2; /*just means that the id is blank, but has now been handled: loop never enters this if-clause again*/
				//chainNr +=1;
				resNrPrior = resNr;
				chainLen =1; /*reset*/
				//printf("Chain id in file is blank, so the recorded value is id to initial value:%c\n", *ptr_chainInStr[chainNr].chainId);
				ptr_chainInStr_temp[chainNr].chainLength = chainLen; /*for the event that the chain only has 1 residue! ...*/
			}
			else if (chainIdIsBlank == 0 && buf[21] != *currentChainId && strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior || buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
				sscanf(buf+21, "%c", currentChainId);
				//printf("Current chain Id in loop:%c\n", *currentChainId);
				//chainNr +=1;
				resNrPrior = resNr;
				chainLen =1; /*reset*/
				sscanf(buf+21, "%c", ptr_chainInStr_temp[chainNr].chainId);
				//printf("Chain id recorded:%c\n", *ptr_chainInStr[chainNr].chainId);
				ptr_chainInStr_temp[chainNr].chainLength = chainLen; /*for the event that the chain only has 1 residue! ...*/
			}
			else if ((chainIdIsBlank == 2 || buf[21] == *currentChainId) && strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior || buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
				resNrPrior = resNr;
				chainLen +=1;
				ptr_chainInStr_temp[chainNr].chainLength = chainLen;
			}
		}

		ptr_chainInStr -> ptr_chainInStr = ptr_chainInStr_temp;

		/*for (i=0;i<= nrOfChains-1; i++){
			printf("chain len of chain %d of id %c is:%d\n", i, *chainInStr.ptr_chainInStr[i].chainId, chainInStr.ptr_chainInStr[i].chainLength);
		}*/

    	return returnVal;
	}

int convertSCOP(struct classInfo *ptr_classInfo, int structureCnt){

	char classId[2]  = "N";
	int  classInt = -1;
	int fold = -1;
	int superFamily = -1;
	int family = -1;

	int i;

	for(i = 0; i < structureCnt; i++){

		//printf("class at %d: %s\n", i, ptr_classInfo[i].classId);
		sscanf(ptr_classInfo[i].classId, "%c.%d.%d.%d", &classId, &fold, &superFamily, &family);
		/*printf("class read in:%c\n", classId[0]);
		printf("fold read in:%d\n", fold);*/

		if(strcmp(classId,"a") == 0){ classInt = 0;}
		else if (strcmp(classId,"b") == 0){ classInt = 1;}
		else if (strcmp(classId,"c") == 0){ classInt = 2;}
		else if (strcmp(classId,"d") == 0){ classInt = 3;}
		else if (strcmp(classId,"e") == 0){ classInt = 4;}
		else if (strcmp(classId,"f") == 0){ classInt = 5;}
		else if (strcmp(classId,"g") == 0){ classInt = 6;}
		else if (strcmp(classId,"h") == 0){ classInt = 7;}
		else if (strcmp(classId,"i") == 0){ classInt = 8;}
		else if (strcmp(classId,"j") == 0){ classInt = 9;}
		else if (strcmp(classId,"k") == 0){ classInt = 10;}
		else{
			printf("Warning: class:%c unrecognized\n", classId);
		}

		ptr_classInfo[i].classIdConv[0] = classInt;
		ptr_classInfo[i].classIdConv[1] = fold;
		ptr_classInfo[i].classIdConv[2] = superFamily; 
		ptr_classInfo[i].classIdConv[3] = family; 

	}

}

int convertCATH(struct classInfo *ptr_classInfo, int structureCnt){

	int  classInt = -1;
	int  architecture = -1;
	int  topology = -1;
	//int  homology = -1;

	int i;

	for(i = 0; i < structureCnt; i++){

		//printf("class at %d: %s\n", i, ptr_classInfo[i].classId);
		sscanf(ptr_classInfo[i].classId, "%d_%d_%d", &classInt, &architecture, &topology);
		/*printf("class read in:%c\n", classId[0]);
		printf("fold read in:%d\n", fold);*/

		ptr_classInfo[i].classIdConv[0] = classInt;
		ptr_classInfo[i].classIdConv[1] = architecture;
		ptr_classInfo[i].classIdConv[2] = topology; 

	}

}

/*for reading and converting a particular list domain names and CAT-values, like the list of 
Kolodny et al. at http://csb.stanford.edu/~rachel/comparison/subset_list_web; the domain names
in this list is as of CATH_2.4 and concist in 6-char's; since that version the format has
changed to 7 char': the two last give the chain nr, so that chain nr 1 is 01, etc. In the
old format only the last (6th) char held the chain nr; we here convert the old format to the
new by inserting a 0 between the 5th and 6th char.
Returns: nr of domains in the list.
Updates: ptr_ptr_cathDomain */
int readCATHDomainList(FILE * ptr_file, struct chainInfo **ptr_ptr_cathDomain, int convToCLF20_b ){ //CLF: CATH list file format 

		int returnVal =0;

		char line[1000];
		int recCnt = 0;
		char *linePart;
		char *line_temp;
		char * tempInit = "N";
		int idxInLine = -1;
		const char *underScore = "_";

		int currentExt = maxNrOfDomains;
		int extSize = 1000; //nr of domains allocated for in each memory extension (reallocation)
		char *structureNameInit = "NN";
		char *classIdInit = "NN";
		char chainIdInit[2] = ">"; /*some "unlikely" value for init of the chain identifier*/

		struct chainInfo *ptr_cathDomain_temp; /*to contain the info -- number, length, id -- of the sub-chains making up the structure*/
		char structureName_temp[lengthCLFformat20 +2] = "\0";
		//char *structureName_temp;

		int i, j =0;

		ptr_cathDomain_temp = *ptr_ptr_cathDomain;

		line_temp = (char *) malloc(200*sizeof(char));


		while (fgets(line,1000, ptr_file)!=NULL){

			/*reallocate mem if necessary*/
			if(recCnt >= currentExt){

				ptr_cathDomain_temp = (struct chainInfo *) realloc (ptr_cathDomain_temp, (currentExt +extSize)*sizeof(struct chainInfo));
				
				for(i = currentExt; i< currentExt + extSize; i++){
					ptr_cathDomain_temp[i].chainId = (char *) malloc(sizeof(char)*chainIdCharSize);
					ptr_cathDomain_temp[i].structureName = (char *) malloc(sizeof(char)*structureNameCharSize);
					ptr_cathDomain_temp[i].classId = (char *) malloc(sizeof(char)*classIdCharSize);
					ptr_cathDomain_temp[i].classIdConv = (int *) malloc(classifierSize*sizeof(int));
					/*init*/
					strcpy(ptr_cathDomain_temp[i].chainId, chainIdInit); /*init to unlikely value*/
					strcpy(ptr_cathDomain_temp[i].structureName, structureNameInit);
					strcpy(ptr_cathDomain_temp[i].classId, classIdInit);
					//printf("Chain Id initial value:%c\n", *ptr_chainInStr[i].chainId);
				
					for(j=0; j< classifierSize; j++){ptr_cathDomain_temp[i].classIdConv[j] = -1;} //to be updated later, if at all
	
				}

				currentExt = recCnt + extSize;

			}

			/*first lines of file may contain start by "#" or a blank (line ahs then length 1) and we skip these)*/
			if(strstr(line, "#") == line || strlen(line) == 1){ continue;}

			//printf("line: %s\n", line);
			//printf("line lgth:%d\n",strlen(line));

			line_temp = strtok(line, "\n"); //effect of this is to take part of line which is not the newline, "\n"

			linePart = strtok(line, " "); //" " is used as the delimiter throughout
			idxInLine = 0;
			//find the nr of windows for this chain/structure:
			while(linePart){

				/*printf("linepart: %s %d", linePart, strcspn(linePart, "\n"));
				getchar();*/

				/*the domain name:*/
				if (idxInLine == 0){
					strcpy(ptr_cathDomain_temp[recCnt].structureName,linePart);
					//printf("str name:%s\n", ptr_cathDomain_temp[recCnt].structureName);
				}

				if (idxInLine == 1){
					strcpy(ptr_cathDomain_temp[recCnt].classId,linePart);
				    ptr_cathDomain_temp[recCnt].classIdConv[0] = atoi(linePart);
				}
			
				if (idxInLine == 2){
					strcat(ptr_cathDomain_temp[recCnt].classId,underScore);
					strcat(ptr_cathDomain_temp[recCnt].classId, linePart);
					ptr_cathDomain_temp[recCnt].classIdConv[1] = atoi(linePart);
				}

				if (idxInLine == 3){
					strcat(ptr_cathDomain_temp[recCnt].classId,underScore);
					strcat(ptr_cathDomain_temp[recCnt].classId, linePart);
					ptr_cathDomain_temp[recCnt].classIdConv[2] = atoi(linePart);
				}

				if (idxInLine == 10){
					ptr_cathDomain_temp[recCnt].chainLength = atoi(linePart);
				}

				linePart = strtok(NULL, " "); //this gets the next token of the string (which is line)
				idxInLine += 1;
			}
			recCnt +=1;
		
		}

		/*if desired convert to CLF 2.0 from previous format:*/
		if(convToCLF20_b ==1){

			//structureName_temp = (char *)malloc(20*sizeof(char));

			for(i=0; i < recCnt; i++){

				structureName_temp[0] = '\0'; //reset
				//printf("structureName_temp:%s\n",structureName_temp);
				strncpy(structureName_temp, ptr_cathDomain_temp[i].structureName,5); //copy first 5 char's
				//scanf(structureName_temp, "%s", ptr_cathDomain_temp[i].structureName);
				structureName_temp[5] = '0'; //insert 0
				structureName_temp[6] = ptr_cathDomain_temp[i].structureName[5] ; //copy last char
				//structureName_temp[7] = '\0'; //add null termination
				//if(i <10){printf("structureName_temp:%s\n",structureName_temp);}
				strcpy(ptr_cathDomain_temp[i].structureName, structureName_temp); //replace 

			}

		}

		/*record results "transformation  style"*/
		*ptr_ptr_cathDomain = ptr_cathDomain_temp;

		return recCnt; 

}

struct cAlpha * main_readPDB(FILE * ptr_file, int chainNumber, int chainLength){
    
		char buf[1000];
		char *chainId;
		char chainIdInit[3] = "Ø"; /*some "unlikely" value for init of the chain identifier*/
		int chainNr = -1;
		int  n =0;
		int resNr = 0;
		int resNrPrior = -999;
		int firstResNr = -999;

		struct cAlpha *ptr_chain = NULL;

		chainId = (char *) malloc(sizeof(char)*3);

		/*allocate memory for array of coordinates*/
		//ptr_chain = (struct cAlpha *) calloc (chainLength, sizeof(struct cAlpha));
		ptr_chain = (struct cAlpha *) malloc (chainLength*sizeof(struct cAlpha));

		/*printf("chainNumber:%d\n",chainNumber);
		printf("chainLength:%d\n",chainLength);*/

		/*populate the array*/
		strcpy(chainId, chainIdInit); /*reset to unlikely value*/
		while (fgets(buf,1000, ptr_file)!=NULL && chainNr <= chainNumber){

			/*get the chainNr:*/
			if (buf[21] != *chainId && strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A'){ 
				/* if chain indicator (at position 21 in buf) is blank, chainNumber must be 0 (there can only
				be 1 chain in case the chain identifier is left blank in the PDB file):*/
				if (buf[21] == ' '){
					if (chainNumber != 0){ 
						printf("Error: chain identifier is blank but chainNumber <> 0\n");
						exit;
					}
				}
				sscanf(buf+21, "%c", chainId);
				//printf("Current chain Id:%c\n", *chainId);
				chainNr +=1;
			}
			/*read in the coordinates if the chainNr is right:*/
			if (chainNr == chainNumber){
				//printf("Chain nr:%d\n", chainNr);
				/*find the residue number*/
				sscanf(buf + 23, "%d\n", &resNr);
				/*record the first residue nr:*/
				if(n==0){
					firstResNr = resNr;
					resNrPrior = resNr -1; /*technical init*/
				}
				//printf("Res nr:%d\n", resNr);
				/*The "resNr == resNrPrior+1" serves to avoid holes in structures*/
				if (strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && (resNr == resNrPrior+1|| buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
					sscanf(buf+30,"%lf %lf %lf",&ptr_chain[n].coords.x, &ptr_chain[n].coords.y, &ptr_chain[n].coords.z);
					ptr_chain[n].residueNr = resNr;
					/*if (n > 100 && n < 160){ 
					printf("n: %d Res nr:%d x y z:%lf %lf %lf\n", n, resNr, ptr_chain[n].x, ptr_chain[n].y, ptr_chain[n].z);
					}*/
					//printf("Chain nr:%d\n", chainNr);
					//printf(buf+30,"%lf %lf %lf\n");
					//printf("x: %lf y:%lf z:%lf \n", ptr_chain[n].x, ptr_chain[n].y, ptr_chain[n].z );
					//printf(buf +9, "%d\n");
					resNrPrior = resNr;
				
					n +=1;
				}
				/*terminate the read in when TER is reached (for the first time):*/
				if (strstr(buf, "TER") == buf){break;}
			}
		}

		//printf("Populated chain Length: %d\n", n);


		/*printf("Chain: \n");
		for ( n=0; n<10; n++ ){
			printf("x: %lf y:%lf z:%lf \n", ptr_chain[n].x, ptr_chain[n].y, ptr_chain[n].z );
		}*/

    	return ptr_chain;

}



int main_readPDB2(FILE * ptr_file, struct cAlpha *ptr_chain, int chainNumber, int chainLength){
    
		int returnVal = 0;

		char buf[1000];
		//char *chainId;
		char chainId[3];
		char chainIdInit[3] = "Ø"; /*some "unlikely" value for init of the chain identifier*/
		int isModel_b = 0;
		int chainNr = -1;
		int  n =0;
		int resNr = 0;
		int resNrPrior = -999;
		int firstResNr = -999;

		//struct cAlpha *ptr_chain = NULL;

		//chainId = (char *) malloc(sizeof(char)*3);
		 
		/*allocate memory for array of coordinates*/
		//ptr_chain = (struct cAlpha *) calloc (chainLength, sizeof(struct cAlpha));
		//ptr_chain = (struct cAlpha *) malloc (chainLength*sizeof(struct cAlpha));

		/*printf("chainNumber:%d\n",chainNumber);
		printf("chainLength:%d\n",chainLength);*/

		/*populate the array*/
		strcpy(chainId, chainIdInit); /*reset to unlikely value*/
		while (fgets(buf,1000, ptr_file)!=NULL && chainNr <= chainNumber){

			/*if the structure is really a set of models we handle it separately*/
			if (strstr(buf, "MODEL") == buf){
				isModel_b = 1;
				chainNr += 1;
				n = 0;
			}
			else{
				/*get the chainNr:*/
				if (isModel_b == 0 && buf[21] != *chainId && (strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf)  && buf[13] == 'C' && buf[14] == 'A'){ 
					/* if chain indicator (at position 21 in buf) is blank, chainNumber must be 0 (there can only
					be 1 chain in case the chain identifier is left blank in the PDB file):*/
					if (buf[21] == ' '){
						if (chainNumber != 0){ 
							printf("Error: chain identifier is blank but chainNumber <> 0\n");
							exit;
						}
					}
					sscanf(buf+21, "%c", chainId);
					//printf("Current chain Id:%c\n", *chainId);
					chainNr +=1;
					n = 0;
				}
			}
			/*read in the coordinates if the chainNr is right:*/
			if (isModel_b == 0){
				if (chainNr == chainNumber){
					/*find the residue number*/
					sscanf(buf + 23, "%d\n", &resNr);
					/*record the first residue nr:*/
					if(n==0){
						resNrPrior = resNr -1; /*technical init*/
					}
					/*The "resNr == resNrPrior+1" serves to avoid holes in structures*/
					if ((strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf)  && buf[13] == 'C' && buf[14] == 'A' && (resNr == resNrPrior+1|| buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
				
						sscanf(buf+30,"%lf %lf %lf",&ptr_chain[n].coords.x, &ptr_chain[n].coords.y, &ptr_chain[n].coords.z);
						ptr_chain[n].residueNr = resNr;
						/*if (n > 100 && n < 160){ 
						printf("n: %d Res nr:%d x y z:%lf %lf %lf\n", n, resNr, ptr_chain[n].x, ptr_chain[n].y, ptr_chain[n].z);
						}*/
						//printf("Chain nr:%d\n", chainNr);
						//printf(buf+30,"%lf %lf %lf\n");
						//printf("x: %lf y:%lf z:%lf \n", ptr_chain[n].x, ptr_chain[n].y, ptr_chain[n].z );
						//printf(buf +9, "%d\n");
						resNrPrior = resNr;
				
						n +=1;
					}
					/*if the chain has a hole and the chain recorded so far is below a global cut-off we reset the chain to allow for the case where the records to follow consitute a proper chain:*/  
					else if ((strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf)  && n < strLengthCutOff  && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior +1)){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
						//printf("Inside chain: %d, (langmodigt venter ...), resNr: %d prior:%d\n", chainNr, resNr, resNrPrior);
						//getchar();
						/*We do not reset the chain coord's already read in, but just the length obtained so far (n)
						We rely on having read in the right chain length, when fetching the chain info:*/
						n = 0;
						sscanf(buf+30,"%lf %lf %lf",&ptr_chain[n].coords.x, &ptr_chain[n].coords.y, &ptr_chain[n].coords.z);
						ptr_chain[n].residueNr = resNr;
						resNrPrior = resNr;
						n += 1;
						
					}
					/*terminate the read in when TER is reached (for the first time):*/
					if (strstr(buf, "TER") == buf||strstr(buf, "END") == buf){break;}
				}
			}
			else{
				if (isModel_b == 1){
					if (chainNr == chainNumber){
						/*find the residue number*/
						sscanf(buf + 23, "%d\n", &resNr);
						/*record the first residue nr:*/
						if(n==0){
							resNrPrior = resNr -1; /*technical init*/
						}
						/*For MODELs we soften the "resNr == resNrPrior+1" (which serves to avoid holes in structures) to just resNr != resNrPrior*/
						/*We also allow HETATMs in the chains*/
						if ((strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf) && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior|| buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
				
							sscanf(buf+30,"%lf %lf %lf",&ptr_chain[n].coords.x, &ptr_chain[n].coords.y, &ptr_chain[n].coords.z);

							ptr_chain[n].residueNr = resNr;

							/*if (n > 100 && n < 160){ 
							printf("n: %d Res nr:%d x y z:%lf %lf %lf\n", n, resNr, ptr_chain[n].x, ptr_chain[n].y, ptr_chain[n].z);
							}*/
							//printf("Chain nr:%d\n", chainNr);
							//printf(buf+30,"%lf %lf %lf\n");
							//printf("x: %lf y:%lf z:%lf \n", ptr_chain[n].x, ptr_chain[n].y, ptr_chain[n].z );
							//printf(buf +9, "%d\n");
							resNrPrior = resNr;
				
							n +=1;
						}
						///*if the chain has a hole and the chain recorded so far is below a global cut-off we reset the chain to allow for the case where the records to follow consitute a proper chain:*/  
						//else if ((strstr(buf, "ATOM") == buf||strstr(buf, "HETATM") == buf)  && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior +1)){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
						//	printf("Inside chain: %d, Sommeren mild oedsles ..., resNr: %d prior:%d ... resetting due to hole\n", chainNr, resNr, resNrPrior);
						//	getchar();
						//	/*We do not reset the chain coord's already read in, but just the length obtained so far (n)
						//	We rely on having read in the right chain length, when fetching the chain info:*/
						//	n = 0;
						//	sscanf(buf+30,"%lf %lf %lf",&ptr_chain[n].x, &ptr_chain[n].y, &ptr_chain[n].z);
						//	resNrPrior = resNr;
						//	n += 1;
						//
						//}
						/*terminate the read in when TER is reached (for the first time):*/
						if (strstr(buf, "TER") == buf){break;}
					}
				}
			}

		}

		if(n != chainLength){
			printf("... Populated chain Length: %d while chainLength: %d\n", n, chainLength);
			//getchar();
		}


		/*printf("Chain: \n");
		for ( n=0; n<10; n++ ){
			printf("x: %lf y:%lf z:%lf \n", ptr_chain[n].x, ptr_chain[n].y, ptr_chain[n].z );
		}*/

    	return returnVal;

	}

/*function computing the square of the distance between two Calphas; used for finding closed loops
in a structure; we compute the square dist rather than the dist simply to avoid a sqrt (which is
unneccessary for our purposes)*/
double distCalphas(struct cAlphaCoords cA1, struct cAlphaCoords cA2){

	double d = 0;

	d = (cA1.x - cA2.x)*(cA1.x - cA2.x) + (cA1.y - cA2.y)*(cA1.y - cA2.y) + (cA1.z - cA2.z)*(cA1.z - cA2.z);

	return d;

}

double I_mutual(double **I, int i1, int i2, int j1, int j2, double stdDev){

	double mutVal  = 0.0;
	//double x,y,z,w;
	//x = I[i1][j2];


	//Could use 
	//mutVal = I[i1][j2] - I[i1][j1-1] - I[i2+1][j2] + I[i2+1][j1-1];
	//but since we're not allocating to e.g. handling i1 = j1 = 0, we use the following if-clause:
	if(i1==j1 && i2 ==j2){
		//mutVal = I[i1][i2];
		mutVal = 0;
	}
	else{
		mutVal = I[i1][j2] - I[i1][j1-1] - I[i2+1][j2] + I[i2+1][j1-1];
	}


	mutVal = (double) mutVal/stdDev;



    return mutVal;

};

/*for computing mean and std dev over 1-dim array*/
struct double2 meanStddev(double *ptr_values, int N){

	int    n  = 0;
	double mean = 0;
	double var = 0;
	double stddev = 0;

	double2 out;

	/*Compute mean and std deviation in two loops (can be done in one but
	may cause numm problems)*/

	for(n = 0; n < N; n++){

		mean += ptr_values[n];

	}

	mean = (double) mean/N;

	for(n = 0; n < N; n++){

		var += pow((ptr_values[n] - mean),2);

	}

	stddev = (double) sqrt(var/N);
	//printf("std dev: %lf\n", stddev);

	out.val[0] = mean;
	out.val[1] = stddev;

	return out;

}


/*simple fct scaling all entries in a 1-d array*/
void scaleArray(double scaleFactor, double **inOutArray, int inOurArrayLength){

	int i = 0;
	double fraction = (double) 1/scaleFactor;

	//printf("Frac:%lf\n", fraction);

	for(i = 0; i < inOurArrayLength; i++){

		/*if(i==0){
			printf("before:%lf\n", (*inOutArray)[i]);
		}*/
		(*inOutArray)[i] =  fraction*(*inOutArray)[i];
		/*if(i==0){
			printf("after:%lf\n", (*inOutArray)[i]);
		}*/
	}

}


/*simple fct normalizing all entries in a 1-d array to mean zero and stddev 1*/
void normalizeArray(double mean, double stddev, double **inOutArray, int inOutArrayLength){

	int i = 0;
	double fraction = (double) 1/stddev;

	//printf("Frac:%lf\n", fraction);

	for(i = 0; i < inOutArrayLength; i++){

		/*if(i==0){
			printf("before:%lf\n", (*inOutArray)[i]);
		}*/
		(*inOutArray)[i] =  fraction*((*inOutArray)[i] - mean);
		/*if(i==0){
			printf("after:%lf\n", (*inOutArray)[i]);
		}*/
	}

}

/*simple fct normalizing all entries in a 2-d "simplex" array to mean zero and stddev 1;
simplex means: the 2-d array has only entries at (i,j) having i <= j and with i ranging from
0 to inOutSimplexSideLength and j up to inOutSimplexSideLength.*/
void normalizeSimplexArray(double mean, double stddev, double ***inOutArray, int simplexSideLength){

	int i = 0, j = 0;
	double fraction = (double) 1/stddev;

	//printf("Frac:%lf\n", fraction);

	for(i = 0; i < simplexSideLength; i++){

		for(j = i; j < simplexSideLength; j++){

		/*if(i==0){
			printf("before:%lf\n", (*inOutArray)[i]);
		}*/

		(*inOutArray)[i][j] =  fraction*((*inOutArray)[i][j] - mean);
		
		/*if(i==0){
			printf("after:%lf\n", (*inOutArray)[i]);
		}*/
		}
	}

}

/*function ordering two input characters, v and w, lexicographically; 
returns: 1 if v < w, -1 if v > w and 0 if v == w
nrOfEntries: length of v
*/
int lexicographicChar(char *v, char *w, int nrOfEntries){

	int returnVal = 0;

	if(nrOfEntries >= 1){

		//printf("v0: %d w0: %d\n", v[0], w[0]);

		if(v[0] < w[0]){returnVal = 1;}

		else if(v[0] > w[0]){returnVal = -1;}

		else{//recurse if v[0] == w[0]

			//printf("recurse");
	
			returnVal = lexicographicChar(v+1,w+1, nrOfEntries -1);

		}

	}

	return returnVal;
}

/*function ordering two input vectors of integers, v and w, lexicographically; 
returns: 1 if v < w, -1 if v > w and 0 if v == w
nrOfEntries: length of v
*/
int lexicographic(int *v, int *w, int nrOfEntries){

	int returnVal = 0;

	if(nrOfEntries >= 1){

		//printf("v0: %d w0: %d\n", v[0], w[0]);

		if(v[0] < w[0]){returnVal = 1;}

		else if(v[0] > w[0]){returnVal = -1;}

		else{//recurse if v[0] == w[0]

			//printf("recurse");
	
			returnVal = lexicographic(v+1,w+1, nrOfEntries -1);

		}

	}

	return returnVal;
}

/*heap sort 1-d array of floats:*/
void siftDown1dArray(double *ptr_values, int start, int end, int decreasing_b){

	int root = start;
	int child;

	double temp;


	if(decreasing_b == 0){

		while(2*root + 1 <= end){ //while root has a child

			child = 2*root + 1; //"left" child
		
			//printf("in sift down ... child:%d\n", child);

			/*if right child exists and its value is larger than left child's: move to right child instead*/
			if(child + 1 <= end && ptr_values[child] < ptr_values[child+1]){
				child = child +1;
			}
			/*swap root and child if value at child > value at root*/
			if(ptr_values[root] < ptr_values[child]){
			
				temp = ptr_values[root];
				ptr_values[root] =  ptr_values[child];
				ptr_values[child] = temp;

				root = child;
			}
			else{
				return;
			}

		} 

	}


	if(decreasing_b == 1){

		while(2*root + 1 <= end){ //while root has a child

			child = 2*root + 1; //"left" child
		
			//printf("in sift down ... child:%d\n", child);

			/*if right child exists and its value is larger than left child's: move to right child instead*/
			if(child + 1 <= end && ptr_values[child] > ptr_values[child+1]){
				child = child +1;
			}
			/*swap root and child if value at child > value at root*/
			if(ptr_values[root] > ptr_values[child]){
			
				temp = ptr_values[root];
				ptr_values[root] =  ptr_values[child];
				ptr_values[child] = temp;

				root = child;
			}
			else{
				return;
			}

		} 

	}

}

void heapify1dArray(double *ptr_values, int nrOfRecords, int decreasing_b){

	int start;

	/*start at last possible root:*/
	if(nrOfRecords%2 == 0){
		
		start = (int) (nrOfRecords - 2)/2;
	}
	else{
		start = (int) (nrOfRecords - 1)/2;
	}

	while(start >= 0){

		//printf("in heapify .. start:%d\n", start);

		/*sift down the element at index start until all nodes
		below that index are in heap order*/
		siftDown1dArray(ptr_values, start, nrOfRecords - 1, decreasing_b);

		start = start - 1;
		/*when root has been sifted down all elts are in heap order*/ 
	}
}

void heapSort1dArray(double *ptr_values, int nrOfRecords, int decreasing_b){
	
	int end;

	double temp;

	/*first pt whole array in max-heap order*/
	heapify1dArray(ptr_values, nrOfRecords, decreasing_b);

	end = nrOfRecords -1;
	while(end > 0){

		//printf("end:%d\n" , end);

		/*swap the present root and the last value of the heap
		(the root is then placed at the end of the heap/array)*/
		temp = ptr_values[end];
		ptr_values[end] =  ptr_values[0];
		ptr_values[0] = temp;

		/*cut the array size by one so as to leave the final elt (which just came from the root)
		out of next step*/
		end = end -1;

		/*sift down the new root so as to heapify the array terminating at the (new) end*/
		siftDown1dArray(ptr_values, 0, end, decreasing_b);

	}

}


/*Bisectional search in a 1-d array, ptr_values, for the interval [i,i+1] for which the value sits
in [ptr_values[i], ptr_values[i+1]]. Returns: i.
nrOfRecords = length of the 1-d array, ptr_values.*/ 
int bisectionalSearchValueL(double value, double *ptr_values, int nrOfRecords){


	int left = 0;
	int right = nrOfRecords - 1 ; //highest index in the array

	int mid = -1;

	double tol = 0.000001;
	
	if(value >= ptr_values[right]){
		//printf("... here's the largest value: %lf and the value at the highest index: %lf\n", value, ptr_values[right]);
		//getchar();
		return right;
	}

	while(right - left > 1){//while(right -left > 1){ //while(mid  != lastRight && left <= right)

		mid = (int) ceil( (double) (right + left)/2 );

		//printf("left: %d rigth: %d mid:%d value left: %lf mid: %lf right: %lf\n", left, right, mid, ptr_values[left], ptr_values[mid], ptr_values[right]);

		if( value < ptr_values[mid]){

			right = mid;

		}
		else{ // if( value >= ptr_values[mid] ){

			left = mid;
		}
	
	}

	if(ptr_values[right] -value < tol){

		return right;

	}
	else{

		return left;

	}



}


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
void siftDownDBClass(struct classInfo *ptr_classInfo, int start, int end, int nrOfEntries, int use_classIdConv_b){

	int nrOfMismatches = 0; //(int *)malloc(2*sizeof(int));

	int root = start;
	int child;

	struct classInfo temp;

	/*if we use the classifier as it is, ie without converting it to an array of int's, we use lex ordering on characters:*/ 
	if(use_classIdConv_b == 0){
		while(2*root + 1 <= end){ //while root has a child

			child = 2*root + 1; //"left" child
		
			//printf("in sift down ... child:%d\n", child);

			/*if right child exists and it svalue is larger than left child's: move to right child instead*/
			nrOfMismatches = 0; //reset
			if(child + 1 <= end && lexicographicChar(ptr_classInfo[child].classId, ptr_classInfo[child +1].classId, nrOfEntries) == 1){
				child = child +1;
				nrOfMismatches = 0; //reset
			}
			nrOfMismatches = 0; //reset
			/*swap root and child if value at child > value at root*/
			if(lexicographicChar(ptr_classInfo[root].classId, ptr_classInfo[child].classId, nrOfEntries) == 1){
			
				temp = ptr_classInfo[root];
				ptr_classInfo[root] =  ptr_classInfo[child];
				ptr_classInfo[child] = temp;

				root = child;
				nrOfMismatches = 0; //reset
			}
			else{
				nrOfMismatches = 0; //reset
				return;
			}

		}
	}

	else{ /*if use converted classifiers, ie classifiers as arrays of int's*/

		while(2*root + 1 <= end){ //while root has a child

			child = 2*root + 1; //"left" child
		
			//printf("in sift down ... child:%d\n", child);

			/*if right child exists and it svalue is larger than left child's: move to right child instead*/
			nrOfMismatches = 0; //reset
			if(child + 1 <= end && lexicographic(ptr_classInfo[child].classIdConv, ptr_classInfo[child +1].classIdConv, nrOfEntries) == 1){
				child = child +1;
				nrOfMismatches = 0; //reset
			}
			nrOfMismatches = 0; //reset
			/*swap root and child if value at child > value at root*/
			if(lexicographic(ptr_classInfo[root].classIdConv, ptr_classInfo[child].classIdConv, nrOfEntries) == 1){
			
				temp = ptr_classInfo[root];
				ptr_classInfo[root] =  ptr_classInfo[child];
				ptr_classInfo[child] = temp;

				root = child;
				nrOfMismatches = 0; //reset
			}
			else{
				nrOfMismatches = 0; //reset
				return;
			}

		}

	}

}

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
void heapifyDBClass(struct classInfo *ptr_classInfo, int nrOfDBrecords, int nrOfEntries, int use_classIdConv_b){

	int start;

	/*start at last possible root:*/
	if(nrOfDBrecords%2 == 0){
		
		start = (int) (nrOfDBrecords - 2)/2;
	}
	else{
		start = (int) (nrOfDBrecords - 1)/2;
	}

	while(start >= 0){

		//printf("in heapify .. start:%d\n", start);

		/*sift down the element at index start until all nodes
		below that index are in heap order*/
		siftDownDBClass(ptr_classInfo, start, nrOfDBrecords - 1, nrOfEntries, use_classIdConv_b);

		start = start - 1;
		/*when root has been sifted down all elts are in heap order*/ 
	}
}

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
void heapSortDBClass(struct classInfo *ptr_classInfo, int nrOfDBrecords, int nrOfEntries, int use_classIdConv_b){
	
	int end;

	struct classInfo temp;

	/*first pt whole array in max-heap order*/
	heapifyDBClass(ptr_classInfo, nrOfDBrecords, nrOfEntries, use_classIdConv_b);

	end = nrOfDBrecords -1;
	while(end > 0){

		//printf("end:%d\n" , end);

		/*swap the present root and the last value of the heap
		(the root is then placed at the end of the heap/array)*/
		temp = ptr_classInfo[end];
		ptr_classInfo[end] =  ptr_classInfo[0];
		ptr_classInfo[0] = temp;

		/*cut the array size by one so as to leave the final elt (which just came from the root)
		out of next step*/
		end = end -1;

		/*sift down the new root so as to heapify the array terminating at the (new) end*/
		siftDownDBClass(ptr_classInfo, 0, end, nrOfEntries, use_classIdConv_b);

	}

}

/*Search through the lexicographically pre-sorted list of classifiers (supposed to be represented as int-arrays of length
nrOfEntries for a given query
returns: an index in sorted list matching the query (if any) else -1*/
int bisectionalSearchClass(struct classInfo classInfo_query, struct classInfo *ptr_classInfo, int nrOfStructures, int nrOfEntries){

	int returnVal = 0;


	int left = 0;
	int right = nrOfStructures -1;
	int lastLeft = left;
	int lastRight = right;

	int mid = 0;
	int midNew = -1;

	int lex = -2;

	int hit; 

	int m, n;

	while(right >= left && mid != lastRight){

		lastLeft = left;
		lastRight = right;

		mid = (int) ceil( (double) (right + left)/2 );

		//printf("mid:%d\n", mid);
		lex = lexicographic(classInfo_query.classIdConv, ptr_classInfo[mid].classIdConv, nrOfEntries);

		if( lex == 1){

			right = mid;

		}
		else if( lex == -1){

			left = mid;

		}
		else{ //lex == 0 so we have a hit "already" 
			 hit = mid;
			 //printf("break\n");
			 break;
		}

		//midNew = (int) ceil( (double) (right + left)/2 );

	}

	if(lex == 1){
		if(lexicographic(classInfo_query.classIdConv, ptr_classInfo[left].classIdConv, nrOfEntries) == 0){hit = left;}
	}
	else if(lex == -1){
		if(lexicographic(classInfo_query.classIdConv, ptr_classInfo[right].classIdConv, nrOfEntries) == 0){hit = right;}
	}

	return hit;

}


/*for extending left or right from a hit (with index hitIndex); searches for extension left from ptr_range.start and up to 
hitIndex, and right from hitIndex and to ptr_range.end. In particular hitIndex must sit in the interval [ptr_range.start,
ptr_range.end] for this to work.
Input: leftRigth_b -- if 0 both extensions are done; if -1 only left, if 1 only right.
Transforms: updates ptr_range .start and .end to left resp right indexes found.*/
int bisectSingleExtendClass(int hitIndex, struct classInfo classInfo, int position, int *ptr_range, struct classInfo *ptr_classInfo, int leftRight_b){

	int returnVal = 0;

	int left;
	int lastLeft;
	int right;
	int lastRight;
	int mid= -1;
	int hit = -1;

	int hitLR[2];

	hitLR[0] = hitIndex;
	hitLR[1] = hitIndex;
	
	if(hitIndex < 0){
		hitLR[0] = -1;
		hitLR[1] = -2;
	}
	else{
		/*left-ward extension:*/
		if(leftRight_b <= 0){

			left = ptr_range[0];
			right = hitIndex;
			lastLeft = left;
			mid = right;
			while(lastLeft < mid) {

				lastLeft = left;
				mid = (int) floor( (double) (right + left)/2 );

				//printf("mid:%d\n", mid);
		
				if( classInfo.classIdConv[position] == ptr_classInfo[mid].classIdConv[position]){

					right = mid;

				}
				else{

					left = mid;

				}

				hit = right;

				//printf("extend right, n:%d\n", n);
				hitLR[0] = hit;

			}

		}


		if(leftRight_b >= 0){
			/*right-ward extension:*/
			left = hitIndex;
			right = ptr_range[1];
			lastRight = right;
			mid = -1;
			while(lastRight - mid > 0 ) {

				lastRight =  right;
				mid = (int) ceil( (double) (right + left)/2 );

				//printf("mid:%d\n", mid);
		
				if( classInfo.classIdConv[position] == ptr_classInfo[mid].classIdConv[position]){

					left = mid;

				}
				else{

					right = mid;

				}

				hit = left;

				//printf("extend right, n:%d\n", n);
				hitLR[1] = hit;

			}
		}
	}
	
	ptr_range[0] = hitLR[0];
	ptr_range[1] = hitLR[1];

	return returnVal;

}


/*Heap sort of CATH domain list on domain name entry (if sortClass_b = 0) or on convClassId (sortClass_b =1):*/
void siftDownCATHDomain(struct chainInfo *ptr_cathDomain, int start, int end, int nrOfEntries, int sortClass_b){

	
	int root = start;
	int child;

	struct chainInfo temp;

	/*if sortClass_b = 0, we sort on the domain name using lex ordering on characters:*/ 
	if(sortClass_b == 0){
		while(2*root + 1 <= end){ //while root has a child

			child = 2*root + 1; //"left" child
		
			//printf("in sift down ... child:%d\n", child);

			/*if right child exists and its value is larger than left child's: move to right child instead*/
			if(child + 1 <= end && lexicographicChar(ptr_cathDomain[child].structureName, ptr_cathDomain[child+1].structureName, nrOfEntries) == 1){
				child = child +1;
			}
			/*swap root and child if value at child > value at root*/
			if(lexicographicChar(ptr_cathDomain[root].structureName, ptr_cathDomain[child].structureName, nrOfEntries) == 1){
			
				temp = ptr_cathDomain[root];
				ptr_cathDomain[root] =  ptr_cathDomain[child];
				ptr_cathDomain[child] = temp;

				root = child;
			}
			else{
				return;
			}

		}
	}

	else{ /*if use "converted classifiers" (for CATH just the array [C, A, T]), ie classifiers as arrays of int's*/

		while(2*root + 1 <= end){ //while root has a child

			child = 2*root + 1; //"left" child
		
			//printf("in sift down ... child:%d\n", child);

			/*if right child exists and it svalue is larger than left child's: move to right child instead*/
			if(child + 1 <= end && lexicographic(ptr_cathDomain[child].classIdConv, ptr_cathDomain[child +1].classIdConv, nrOfEntries) == 1){
				child = child +1;
			}
			/*swap root and child if value at child > value at root*/
			if(lexicographic(ptr_cathDomain[root].classIdConv, ptr_cathDomain[child].classIdConv, nrOfEntries) == 1){
			
				temp = ptr_cathDomain[root];
				ptr_cathDomain[root] =  ptr_cathDomain[child];
				ptr_cathDomain[child] = temp;

				root = child;
			}
			else{
				return;
			}

		}

	}

}

void heapifyCATHDomain(struct chainInfo *ptr_cathDomain, int nrOfDBrecords, int nrOfEntries, int sortClass_b){

	int start;

	/*start at last possible root:*/
	if(nrOfDBrecords%2 == 0){
		
		start = (int) (nrOfDBrecords - 2)/2;
	}
	else{
		start = (int) (nrOfDBrecords - 1)/2;
	}

	while(start >= 0){

		//printf("in heapify .. start:%d\n", start);

		/*sift down the element at index start until all nodes
		below that index are in heap order*/
		siftDownCATHDomain(ptr_cathDomain, start, nrOfDBrecords - 1, nrOfEntries, sortClass_b);

		start = start - 1;
		/*when root has been sifted down all elts are in heap order*/ 
	}
}

void heapSortCATHDomain(struct chainInfo *ptr_cathDomain, int nrOfDBrecords, int nrOfEntries, int sortClass_b){
	
	int end;

	struct chainInfo temp;

	/*first pt whole array in max-heap order*/
	heapifyCATHDomain(ptr_cathDomain, nrOfDBrecords, nrOfEntries, sortClass_b);

	end = nrOfDBrecords -1;
	while(end > 0){

		//printf("end:%d\n" , end);

		/*swap the present root and the last value of the heap
		(the root is then placed at the end of the heap/array)*/
		temp = ptr_cathDomain[end];
		ptr_cathDomain[end] =  ptr_cathDomain[0];
		ptr_cathDomain[0] = temp;

		/*cut the array size by one so as to leave the final elt (which just came from the root)
		out of next step*/
		end = end -1;

		/*sift down the new root so as to heapify the array terminating at the (new) end*/
		siftDownCATHDomain(ptr_cathDomain, 0, end, nrOfEntries, sortClass_b);

	}

}

/*search for a domain (query) in the list/array ptr_cathDomain of CATH domains. 
Returns: the hit index (-1 if no hit) */
int bisectionalSearchCATHDomain(char *cathDomain_query, struct chainInfo *ptr_cathDomain, int nrOfStructures, int nrOfLetters){

	int returnVal = 0;

	int left = 0;
	int right = nrOfStructures -1;
	int lastLeft = left;
	int lastRight = right;

	int mid = -1;

	int lex = -2;

	int hit = -1; 

	int m, n;

	while(right >= left && mid != lastRight){

		lastLeft = left;
		lastRight = right;

		mid = (int) ceil( (double) (right + left)/2 );

		//printf("mid:%d\n", mid);
		lex = lexicographicChar(cathDomain_query, ptr_cathDomain[mid].structureName, nrOfLetters);

		if( lex == 1){

			right = mid;

		}
		else if( lex == -1){

			left = mid;

		}
		else{ //lex == 0 so we have a hit "already" 
			 hit = mid;
			 //printf("break\n");
			 break;
		}

		//midNew = (int) ceil( (double) (right + left)/2 );

	}

	if(lex == 1){
		if(lexicographicChar(cathDomain_query, ptr_cathDomain[left].structureName, nrOfLetters) == 0){hit = left;}
	}
	else if(lex == -1){
		if(lexicographicChar(cathDomain_query, ptr_cathDomain[right].structureName, nrOfLetters) == 0){hit = right;}
	}

	return hit;

}


/*Functions for computing w-terms*/

/*function computing the w-value of a pair of segments, the latter input as 12 doubles*/
double w( 
	/* coords of 1st point in 1st segment*/
	double s101,
	double s102,
	double s103,
	/* coords of 2nd point in 1st segment*/
	double s111,
	double s112,
	double s113,
	/* coords of 1st point in 2nd segment*/
	double s201,
	double s202,
	double s203,
	/* coords of 2nd point in 2nd segment*/
	double s211,
	double s212,
	double s213
	)
{ 

	/*form unit vectors at the segment extremes:*/
	/*v13:*/
    double v131 = s201 - s101;
    double v132 = s202 - s102;
    double v133 = s203 - s103;
    /*v23:*/
    double v231 = s201 - s111;
    double v232 = s202 - s112;
    double v233 = s203 - s113;
    /*v24:*/
    double v241 = s211 - s111;
    double v242 = s212 - s112;
    double v243 = s213 - s113;
    /*v14:*/
    double v141 = s211 - s101;
    double v142 = s212 - s102;
    double v143 = s213 - s103;

    /*   v = [v13,v23,v24,v14,v13,v23]
    e = []*/
    
	double l13 = sqrt(v131*v131 + v132*v132 + v133*v133);
    double l23 = sqrt(v231*v231 + v232*v232 + v233*v233);
    double l24 = sqrt(v241*v241 + v242*v242 + v243*v243);
    double l14 = sqrt(v141*v141 + v142*v142 + v143*v143);
   
	//double ls[4] = {l13,l23,l24,l14};

	double e131;
    double e132;
    double e133;
	/*e13 = (e131,e132,e133)*/
    double e231;
    double e232;
    double e233;
	/*e23 = (e231,e232,e233)*/
    double e241;
    double e242;
    double e243;
	/*e24 = (e241,e242,e243)*/
    double e141;
    double e142;
    double e143;
	/*e14 = (e141,e142,e143)
	e = (e13,e23,e24,e14,e13,e23)*/
    /*compute the angles*/
    double s = 0;
    double signS = 0;

	/*what we want to do, here in python style:
#    for i in range(1,len(e)-1):
#        theta = 0 #reset
#        a = e[i-1]
#        b =e[i]
#        c = e[i+1]
##        theta = PDB.calc_angle(a,b,c)
##        print "angle is %lf" % theta 
#        aDotb = np.dot(a,b) #PDB.Vector.__mul__(a,b)
#        aDotc = np.dot(a,c) #PDB.Vector.__mul__(a,c)
#        bDotc = np.dot(b,c) #PDB.Vector.__mul__(b,c)
#        aCrossb = np.cross(a,b)*/
    /*We have to do this explicity, ie without referrring to itereations and vectors:*/
    /*A = e13, B= e23, C = e24:*/
    double x;
    double y;
    double z;
    double a;
    double b;
    double c;
    double u;
    double v;
    double w;
    double A_Dot_B;
    double A_Dot_C;
    double B_Dot_C;
    double A_Cross_B_1;
	double A_Cross_B_2;
	double A_Cross_B_3;
/* python style:
###        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
###        cos = double(aDotb*bDotc - aDotc)/denom
###        sin = double(np.dot(aCrossb,c))/denom
##        #can disregard the denominator as they are equal:*/
    double cos;
    double sin;
	double theta;

	l13 = v131*v131 + v132*v132 + v133*v133;
    l23 = v231*v231 + v232*v232 + v233*v233;
    l24 = v241*v241 + v242*v242 + v243*v243;
    l14 = v141*v141 + v142*v142 + v143*v143;

	/* check if segments are adjacent, i.e. that vectors v has length "zero":*/
	if (l13 < 1e-6){
			return 0;
	}
	if (l23 < 1e-6){
			return 0;
	}
	if (l24 < 1e-6){
			return 0;
	}
	if (l14 < 1e-6){
			return 0;
	}

	l13 = sqrt(l13);
    l23 = sqrt(l23);
    l24 = sqrt(l24);
    l14 = sqrt(l14);

	e131 = v131/l13;
    e132 = v132/l13;
    e133 = v133/l13;
	/*e13 = (e131,e132,e133)*/
    e231 = v231/l23;
    e232 = v232/l23;
    e233 = v233/l23;
	/*e23 = (e231,e232,e233)*/
    e241 = v241/l24;
    e242 = v242/l24;
    e243 = v243/l24;
	/*e24 = (e241,e242,e243)*/
    e141 = v141/l14;
    e142 = v142/l14;
    e143 = v143/l14;

	x = e131;
    y = e132;
    z = e133;
    a = e231;
    b = e232;
    c = e233;
    u = e241;
    v = e242;
    w = e243;
    A_Dot_B = x*a + y*b + z*c;
    A_Dot_C = x*u + y*v + z*w;
    B_Dot_C = a*u + b*v + c*w;
    A_Cross_B_1 = y*c - z*b;
	A_Cross_B_2 = z*a - x*c;
	A_Cross_B_3 = x*b - y*a;
/* python style:
###        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
###        cos = double(aDotb*bDotc - aDotc)/denom
###        sin = double(np.dot(aCrossb,c))/denom
##        #can disregard the denominator as they are equal:*/
    cos = A_Dot_B*B_Dot_C - A_Dot_C;
    sin = A_Cross_B_1*u + A_Cross_B_2*v + A_Cross_B_3*w;

	 /*for angle in [-pi ,pi]:*/
	theta = atan2(sin,cos);

    /*printf("%lf %lf %lf", cos, sin, theta);*/

    s = s + theta;

	/*A = e23, B= e24, C = e14:*/
    x = e231;
    y = e232;
    z = e233;
    a = e241;
    b = e242;
    c = e243;
    u = e141;
    v = e142;
    w = e143;
	A_Dot_B = x*a + y*b + z*c;
    A_Dot_C = x*u + y*v + z*w;
    B_Dot_C = a*u + b*v + c*w;
    A_Cross_B_1 = y*c - z*b;
	A_Cross_B_2 = z*a - x*c;
	A_Cross_B_3 = x*b - y*a;
/* python style:
###        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
###        cos = double(aDotb*bDotc - aDotc)/denom
###        sin = double(np.dot(aCrossb,c))/denom
##        #can disregard the denominator as they are equal:*/
    cos = A_Dot_B*B_Dot_C - A_Dot_C;
    sin = A_Cross_B_1*u + A_Cross_B_2*v + A_Cross_B_3*w;

	/*for angle in [-pi ,pi]:*/
	theta = atan2(sin,cos);

    /*printf("%lf %lf %lf", cos, sin, theta);*/

	s += theta;

    /*A = e24, B= e14, C = e13:*/
    x = e241;
    y = e242;
    z = e243;
    a = e141;
    b = e142;
    c = e143;
    u = e131;
    v = e132;
    w = e133;
    A_Dot_B = x*a + y*b + z*c;
    A_Dot_C = x*u + y*v + z*w;
    B_Dot_C = a*u + b*v + c*w;
    A_Cross_B_1 = y*c - z*b;
	A_Cross_B_2 = z*a - x*c;
	A_Cross_B_3 = x*b - y*a;
/* python style:
###        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
###        cos = double(aDotb*bDotc - aDotc)/denom
###        sin = double(np.dot(aCrossb,c))/denom
##        #can disregard the denominator as they are equal:*/
    cos = A_Dot_B*B_Dot_C - A_Dot_C;
    sin = A_Cross_B_1*u + A_Cross_B_2*v + A_Cross_B_3*w;

	/*for angle in [-pi ,pi]:*/
	theta = atan2(sin,cos);

	/*printf("%lf %lf %lf", cos, sin, theta);*/

	s += theta;

	/*A = e14, B= e13, C = e23:*/
    x = e141;
    y = e142;
    z = e143;
    a = e131;
    b = e132;
    c = e133;
    u = e231;
    v = e232;
    w = e233;
    A_Dot_B = x*a + y*b + z*c;
    A_Dot_C = x*u + y*v + z*w;
    B_Dot_C = a*u + b*v + c*w;
    A_Cross_B_1 = y*c - z*b;
	A_Cross_B_2 = z*a - x*c;
	A_Cross_B_3 = x*b - y*a;
/* python style:
###        denom = np.sqrt((1 - aDotb*aDotb)*(1- bDotc*bDotc))
###        cos = double(aDotb*bDotc - aDotc)/denom
###        sin = double(np.dot(aCrossb,c))/denom
##        #can disregard the denominator as they are equal:*/
    cos = A_Dot_B*B_Dot_C - A_Dot_C;
    sin = A_Cross_B_1*u + A_Cross_B_2*v + A_Cross_B_3*w;

	/*for angle in [-pi ,pi]:*/
	theta = atan2(sin,cos);

   /*printf("%lf %lf %lf", cos, sin, theta);*/
    
	s += theta;

	/*compute the sign of s:*/
	signS = (s>0) - (s<0);

	return signS*2*pi - s;
	
}


double w_seg12(struct segment seg1, struct segment seg2)
{ 
	/*struct segment seg1 = twoSeg.seg1;
	struct segment seg2 = twoSeg.seg2;*/

	/* coords of 1st point in 1st segment*/
	double s101 = seg1.s1.x;
	double s102 = seg1.s1.y;
	double s103 = seg1.s1.z;

	/* coords of 2nd point in 1st segment*/
	double s111 = seg1.s2.x;
	double s112 = seg1.s2.y;
	double s113 = seg1.s2.z;

	/* coords of 1st point in 2nd segment*/
	double s201 = seg2.s1.x;
	double s202 = seg2.s1.y;
	double s203 = seg2.s1.z;

	/* coords of 2nd point in 2nd segment*/
	double s211 = seg2.s2.x;
	double s212 = seg2.s2.y;
	double s213 = seg2.s2.z;

	double l1, l2;
    
	/*If one of the two segments is actually/practically the null-vector we force the fct to return 0 (as it should):*/    
	l1 = (seg1.s2.z - seg1.s1.z)*(seg1.s2.z - seg1.s1.z) + (seg1.s2.y - seg1.s1.y)*(seg1.s2.y - seg1.s1.y) + (seg1.s2.x - seg1.s1.x)*(seg1.s2.x - seg1.s1.x);

	if(l1 < 1e-6){
		return 0.0;
	}
	else{
		l2 = (seg2.s2.z - seg2.s1.z)*(seg2.s2.z - seg2.s1.z) + (seg2.s2.y - seg2.s1.y)*(seg2.s2.y - seg2.s1.y) + (seg2.s2.x - seg2.s1.x)*(seg2.s2.x - seg2.s1.x);
		if(l2 < 1e-6){
			return 0.0;
		}
		else{/*call the computation*/
			return w(s101, s102, s103, s111, s112, s113, s201, s202, s203, s211, s212, s213);
		}
	}
	
	}

/*function for computing all w-values across the simplex for a given sturcture. 
This presupposes that adequate memory is 
allocated to the "double array" wVal*/
int wAll(struct segment *ptr_segment, int chainLen, double **wVal){
	
	int i = 0;
	int j = 0;
	int L = 0;
	
	L = chainLen - 1;

	for ( i=0; i<=L-1; i++ ){
			for ( j=i; j<=L-1; j++ ){ /* We include i = j for the purpose of checking that w[diagonal] =0. Only in order 0 */
					wVal[i][j] = w_seg12(ptr_segment[i], ptr_segment[j]);
			}
	}

	return 1;
}

/*Aggregation/computation of invariants (including the abolute value versions) by recursion
through the simplex (traversing the invariants and the simplex). Takes w-values for the chain 
as input (wVal "double array"):*/
int aggr(int chainLen, int order, int full_b, double **wVal, struct I_ptr I_measures ){

	int L = 0;
	int i = 0;
	int j = 0;
	int b = 0;
	int k = 0;
	int c = 0;
	
	/* the parts of the output structure to hold the output:*/
	double wValue;
	double abs_wValue; 
	double wVal1;
	double abs_wVal1;
	double wVal2;
	double abs_wVal2;

	/*next a series of pointers, all copies of the input/output I_measures.* pointers;
	they are merely there for convenience -- readability:*/

	double **wValCopy = I_measures.wVal;
	/*pointers for order 1 to 3 measures:*/
	double **I12 = I_measures.I12;
	double **Ia12 = I_measures.Ia12;

	double **I1234 = I_measures.I1234;
	double **I1324 = I_measures.I1324;
	double **I1423  = I_measures.I1423;

	/* for "absolute value version" of order 2 measures: */
	double **Ia1234 = I_measures.Ia1234;
	double **I12a34 = I_measures.I12a34;
	double **Ia12a34 = I_measures.Ia12a34;

	double **Ia1324 = I_measures.Ia1324;
	double **I13a24 = I_measures.I13a24;
	double **Ia13a24 = I_measures.Ia13a24;

	double **Ia1423  = I_measures.Ia1423;
	double **I14a23  = I_measures.I14a23;
	double **Ia14a23  = I_measures.Ia14a23;


	/* for full computation of the order 2 measures across the simplex:*/
	double ** I1234_full_aid = I_measures.I1234_full_aid;
	double ** I1324_full_aid = I_measures.I1324_full_aid;
    double ** I1234_full = I_measures.I1234_full;
    double ** I1324_full = I_measures.I1324_full;
	/* for full computation of the "absolute value version" order 2 measures across the simplex:*/
	double ** Ia12a34_full_aid = I_measures.Ia12a34_full_aid;
	double ** Ia1234_full_aid = I_measures.Ia1234_full_aid;
	double ** I12a34_full_aid = I_measures.I12a34_full_aid;

    double ** Ia12a34_full = I_measures.Ia12a34_full;
	double ** Ia1234_full = I_measures.Ia1234_full;
	double ** I12a34_full = I_measures.I12a34_full;

	double ** Ia13a24_full_aid = I_measures.Ia13a24_full_aid;
	double ** Ia1324_full_aid = I_measures.Ia1324_full_aid;
	double ** I13a24_full_aid = I_measures.I13a24_full_aid;

    double ** Ia13a24_full = I_measures.Ia13a24_full;
	double ** Ia1324_full = I_measures.Ia1324_full;
	double ** I13a24_full = I_measures.I13a24_full;


	/* assisting the order 3 case:*/
	double ** I1324_full2_aid = I_measures.I1324_full2_aid;
	double ** I1324_full2 = I_measures.I1324_full2;
	double **I1423_full0 = I_measures.I1423_full0;
	double **I1423_full2_aid = I_measures.I1423_full2_aid;
	double **I1423_full2 = I_measures.I1423_full2;

	/*order 3*/
	/*12*/
	double **I123456 = I_measures.I123456;
	double **I123645 = I_measures.I123645;
	double **I123546 = I_measures.I123546;
	/*13*/
	double **I132456 = I_measures.I132456;
	double **I132546 = I_measures.I132546;
	double **I132645 = I_measures.I132645;
	/*14*/
	double **I142356 = I_measures.I142356;
	double **I142536 = I_measures.I142536;
	double **I142635 = I_measures.I142635;
	/*15*/
	double **I152346 = I_measures.I152346;
	double **I152436 = I_measures.I152436;
	double **I152634 = I_measures.I152634;
	/*16*/
	double **I162345 = I_measures.I162345;
	double **I162435 = I_measures.I162435;
	double **I162534 = I_measures.I162534;

	L = chainLen -1;

	I_measures.chainLen = &chainLen; 
	I_measures.order = &order;

	if (order == 0){

		for ( i=0; i<=L-1; i++ ){
			for ( j=i+1; j<=L-1; j++ ){ 
					wValCopy[i][j] = wVal[i][j];
			}
		}
	}

	if (order == 1){

		for ( i=L-1; i>= 0; i-- ){
			for ( j=i+1; j<=L-1; j++ ){
				wValue = wVal[i][j];
				wValCopy[i][j] = wValue;
				I12[i][j] = wValue + I12[i][j-1] + I12[i+1][j] - I12[i+1][j-1];
				Ia12[i][j] = fabsf(wValue) + Ia12[i][j-1] + Ia12[i+1][j] - Ia12[i+1][j-1];

			}
		}
	}


	if (order == 2){

		for ( i=L-1; i>= 0; i-- ){
			for ( j=i+1; j<=L-1; j++ ){
					wValue = wVal[i][j];
					wValCopy[i][j] = wValue;
					abs_wValue = fabsf(wValue);
					I12[i][j] = wValue + I12[i][j-1] + I12[i+1][j] - I12[i+1][j-1];
					Ia12[i][j] = abs_wValue + Ia12[i][j-1] + Ia12[i+1][j] - Ia12[i+1][j-1];

					/*I1234: I1234(i,j;N)*/
                    I1234[i][j] = I12[j+1][L-1]*wValue + I1234[i][j-1] + I1234[i+1][j] - I1234[i+1][j-1];
                    /*I1423: I1423(i,j)*/
					I1423[i][j] = I12[i+1][j-1]*wValue +  I1423[i][j-1] + I1423[i+1][j] - I1423[i+1][j-1];
                    /*I1324: I1324(i,j;N)*/
                    I1324[i][j] = (I12[i+1][L-1] - I12[i+1][j] -I12[j][L-1] )*wValue + I1324[i][j-1] + I1324[i+1][j] - I1324[i+1][j-1];

					/*absolute value versions:*/
                    Ia1234[i][j] = I12[j+1][L-1]*abs_wValue + Ia1234[i][j-1] + Ia1234[i+1][j] - Ia1234[i+1][j-1];
                    I12a34[i][j] = Ia12[j+1][L-1]*wValue + I12a34[i][j-1] + I12a34[i+1][j] - I12a34[i+1][j-1];
                    Ia12a34[i][j] = Ia12[j+1][L-1]*abs_wValue + Ia12a34[i][j-1] + Ia12a34[i+1][j] - Ia12a34[i+1][j-1];

					Ia1324[i][j] = (I12[i+1][L-1] - I12[i+1][j] -I12[j][L-1] )*abs_wValue + Ia1324[i][j-1] + Ia1324[i+1][j] - Ia1324[i+1][j-1];
					I13a24[i][j] = (Ia12[i+1][L-1] - Ia12[i+1][j] -Ia12[j][L-1] )*wValue + I13a24[i][j-1] + I13a24[i+1][j] - I13a24[i+1][j-1];
					Ia13a24[i][j] = (Ia12[i+1][L-1] - Ia12[i+1][j] -Ia12[j][L-1] )*abs_wValue + Ia13a24[i][j-1] + Ia13a24[i+1][j] - Ia13a24[i+1][j-1];

					Ia1423[i][j] = I12[i+1][j-1]*abs_wValue +  Ia1423[i][j-1] + Ia1423[i+1][j] - Ia1423[i+1][j-1];
					I14a23[i][j] = Ia12[i+1][j-1]*wValue +  I14a23[i][j-1] + I14a23[i+1][j] - I14a23[i+1][j-1];
					Ia14a23[i][j] = Ia12[i+1][j-1]*abs_wValue +  Ia14a23[i][j-1] + Ia14a23[i+1][j] - Ia14a23[i+1][j-1];


					if (full_b == 1){
						b = i+1;
                        for (b = i+1; b <= j-2; b++){ /*#obs: the w-values needed have been computed*/
						 	wVal1 = wVal[i][b];
							abs_wVal1 = fabsf(wVal1);
							I1234_full_aid[i][j] += I12[b+1][j]*wVal1;
                            I1324_full_aid[i][j] += I12[b][j]*wVal1;

							Ia12a34_full_aid[i][j] += Ia12[b+1][j]*abs_wVal1;
							Ia1234_full_aid[i][j] += Ia12[b+1][j]*wVal1;
							I12a34_full_aid[i][j] += I12[b+1][j]*abs_wVal1;

							Ia13a24_full_aid[i][j] += I12[b][j]*abs_wVal1;
							Ia1324_full_aid[i][j] += Ia12[b][j]*wVal1;
							I13a24_full_aid[i][j] += I12[b][j]*abs_wVal1;


						}
                        /*I1234_full add up terms and recursion: I1234(i,j)*/
                        I1234_full[i][j] = I1234_full_aid[i][j] - I1234_full_aid[i][j-1];
                        I1234_full[i][j] += I1234_full[i][j-1] + I1234_full[i+1][j] - I1234_full[i+1][j-1];

                        /*I1324_full add up terms recursion: I1324(i,j)*/
                        I1324_full[i][j] = I1324_full_aid[i][j-1] - I1324_full_aid[i][j];
                        I1324_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + I1324_full[i][j-1] + I1324_full[i+1][j] - I1324_full[i+1][j-1];
						
						/*abs value versions:*/
						Ia12a34_full[i][j] = Ia12a34_full_aid[i][j] - Ia12a34_full_aid[i][j-1];
                        Ia12a34_full[i][j] += Ia12a34_full[i][j-1] + Ia12a34_full[i+1][j] - Ia12a34_full[i+1][j-1];

						Ia1234_full[i][j] = Ia1234_full_aid[i][j] - Ia1234_full_aid[i][j-1];
                        Ia1234_full[i][j] += Ia1234_full[i][j-1] + Ia1234_full[i+1][j] - Ia1234_full[i+1][j-1];

						I12a34_full[i][j] = I12a34_full_aid[i][j] -I12a34_full_aid[i][j-1];
                        I12a34_full[i][j] += I12a34_full[i][j-1] + I12a34_full[i+1][j] - I12a34_full[i+1][j-1];

						Ia13a24_full[i][j] = Ia13a24_full_aid[i][j-1] - Ia13a24_full_aid[i][j];
                        Ia13a24_full[i][j] += (Ia12[i][j-1] - Ia12[i+1][j-1])*(Ia12[i+1][j] - Ia12[i+1][j-1]) + Ia13a24_full[i][j-1] + Ia13a24_full[i+1][j] - Ia13a24_full[i+1][j-1];

						Ia1324_full[i][j] = Ia1324_full_aid[i][j-1] - Ia1324_full_aid[i][j];
                        Ia1324_full[i][j] += (Ia12[i][j-1] - Ia12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + Ia1324_full[i][j-1] + Ia1324_full[i+1][j] - Ia1324_full[i+1][j-1];

						I13a24_full[i][j] = I13a24_full_aid[i][j-1] - I13a24_full_aid[i][j];
                        I13a24_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(Ia12[i+1][j] - Ia12[i+1][j-1]) + I13a24_full[i][j-1] + I13a24_full[i+1][j] - I13a24_full[i+1][j-1];

					}
			}
		}
	}




	if (order == 3){

		for ( i=L-1; i>= 0; i-- ){
			for ( j=i+1; j<=L-1; j++ ){
					wValue = wVal[i][j];
					/*sscanf(wValue,"%lf",&wVal[cnt]); */
					wValCopy[i][j] = wValue;
					abs_wValue = fabsf(wValue);
					I12[i][j] = wValue + I12[i][j-1] + I12[i+1][j] - I12[i+1][j-1];
					Ia12[i][j] = abs_wValue + Ia12[i][j-1] + Ia12[i+1][j] - Ia12[i+1][j-1];

					/*I1234: I1234(i,j;N)*/
                    I1234[i][j] = I12[j+1][L-1]*wValue + I1234[i][j-1] + I1234[i+1][j] - I1234[i+1][j-1];
                    /*I1423: I1423(i,j)*/
					I1423[i][j] = I12[i+1][j-1]*wValue +  I1423[i][j-1] + I1423[i+1][j] - I1423[i+1][j-1];
                    /*I1324: I1324(i,j;N)*/
                    I1324[i][j] = (I12[i+1][L-1] - I12[i+1][j] -I12[j][L-1] )*wValue + I1324[i][j-1] + I1324[i+1][j] - I1324[i+1][j-1];
					
					/*absolute value versions:*/
                    Ia1234[i][j] = I12[j+1][L-1]*abs_wValue + Ia1234[i][j-1] + Ia1234[i+1][j] - Ia1234[i+1][j-1];
                    I12a34[i][j] = Ia12[j+1][L-1]*wValue + I12a34[i][j-1] + I12a34[i+1][j] - I12a34[i+1][j-1];
                    Ia12a34[i][j] = Ia12[j+1][L-1]*abs_wValue + Ia12a34[i][j-1] + Ia12a34[i+1][j] - Ia12a34[i+1][j-1];

					Ia1324[i][j] = (I12[i+1][L-1] - I12[i+1][j] -I12[j][L-1] )*abs_wValue + Ia1324[i][j-1] + Ia1324[i+1][j] - Ia1324[i+1][j-1];
					I13a24[i][j] = (Ia12[i+1][L-1] - Ia12[i+1][j] -Ia12[j][L-1] )*wValue + I13a24[i][j-1] + I13a24[i+1][j] - I13a24[i+1][j-1];
					Ia13a24[i][j] = (Ia12[i+1][L-1] - Ia12[i+1][j] -Ia12[j][L-1] )*abs_wValue + Ia13a24[i][j-1] + Ia13a24[i+1][j] - Ia13a24[i+1][j-1];

					Ia1423[i][j] = I12[i+1][j-1]*abs_wValue +  Ia1423[i][j-1] + Ia1423[i+1][j] - Ia1423[i+1][j-1];
					I14a23[i][j] = Ia12[i+1][j-1]*wValue +  I14a23[i][j-1] + I14a23[i+1][j] - I14a23[i+1][j-1];
					Ia14a23[i][j] = Ia12[i+1][j-1]*abs_wValue +  Ia14a23[i][j-1] + Ia14a23[i+1][j] - Ia14a23[i+1][j-1];


					b = i+1;
                    for (b = i+1; b <= j-2; b++){ /*#obs: the w-values needed have been computed*/
						wVal1 = wVal[i][b];
						abs_wVal1 = fabsf(wVal1);
						I1234_full_aid[i][j] += I12[b+1][j]*wVal1;
                        I1324_full_aid[i][j] += I12[b][j]*wVal1;

						Ia12a34_full_aid[i][j] += Ia12[b+1][j]*abs_wVal1;
						Ia1234_full_aid[i][j] += Ia12[b+1][j]*wVal1;
						I12a34_full_aid[i][j] += I12[b+1][j]*abs_wVal1;

						Ia13a24_full_aid[i][j] += I12[b][j]*abs_wVal1;
						Ia1324_full_aid[i][j] += Ia12[b][j]*wVal1;
						I13a24_full_aid[i][j] += I12[b][j]*abs_wVal1;


					}
                    /*I1234_full add up terms and recursion: I1234(i,j)*/
                    I1234_full[i][j] = I1234_full_aid[i][j] -I1234_full_aid[i][j-1];
                    I1234_full[i][j] += I1234_full[i][j-1] + I1234_full[i+1][j] - I1234_full[i+1][j-1];
                    /*I1324_full add up terms recursion: I1324(i,j)*/
                    I1324_full[i][j] = I1324_full_aid[i][j-1] - I1324_full_aid[i][j];
                    I1324_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + I1324_full[i][j-1] + I1324_full[i+1][j] - I1324_full[i+1][j-1];
					
					/*abs value versions:*/
					Ia12a34_full[i][j] = Ia12a34_full_aid[i][j] - Ia12a34_full_aid[i][j-1];
                    Ia12a34_full[i][j] += Ia12a34_full[i][j-1] + Ia12a34_full[i+1][j] - Ia12a34_full[i+1][j-1];

					Ia1234_full[i][j] = Ia1234_full_aid[i][j] - Ia1234_full_aid[i][j-1];
                    Ia1234_full[i][j] += Ia1234_full[i][j-1] + Ia1234_full[i+1][j] - Ia1234_full[i+1][j-1];

					I12a34_full[i][j] = I12a34_full_aid[i][j] -I12a34_full_aid[i][j-1];
                    I12a34_full[i][j] += I12a34_full[i][j-1] + I12a34_full[i+1][j] - I12a34_full[i+1][j-1];

					Ia13a24_full[i][j] = Ia13a24_full_aid[i][j-1] - Ia13a24_full_aid[i][j];
                    Ia13a24_full[i][j] += (Ia12[i][j-1] - Ia12[i+1][j-1])*(Ia12[i+1][j] - Ia12[i+1][j-1]) + Ia13a24_full[i][j-1] + Ia13a24_full[i+1][j] - Ia13a24_full[i+1][j-1];

					Ia1324_full[i][j] = Ia1324_full_aid[i][j-1] - Ia1324_full_aid[i][j];
                    Ia1324_full[i][j] += (Ia12[i][j-1] - Ia12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + Ia1324_full[i][j-1] + Ia1324_full[i+1][j] - Ia1324_full[i+1][j-1];

					I13a24_full[i][j] = I13a24_full_aid[i][j-1] - I13a24_full_aid[i][j];
                    I13a24_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(Ia12[i+1][j] - Ia12[i+1][j-1]) + I13a24_full[i][j-1] + I13a24_full[i+1][j] - I13a24_full[i+1][j-1];

					
					/*I1423_full0: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet*/
                    I1423_full0[i+1][j] = (I12[i+2][j-1] - I12[i+2][j])*(I12[i+2][L-1] - I12[i+2][j] - I12[i+1][L-1] + I12[i+1][j]);
                    I1423_full0[i+1][j] += I1423_full0[i+2][j] + I1423_full0[i+1][j-1] - I1423_full0[i+2][j-1];
                    /*to compute certain degree 6 measures we use two auxiliary order 2
                    measures; the recursion demands to sum j in the - direction (ie 
                    from above and down):*/
                    k = L-1 + i -j+1;
                    for (c =k+1;c <= L-1; c++){ /*Obs: the w values needed have already been computed (since k,c >i)*/
                        wVal2 = wVal[k][c];
                        I1324_full2_aid[i+1][k] += wVal2*(I12[i+1][c-1] - I12[i+1][k] - I12[k][c-1]);
                        I1423_full2_aid[i+1][k] += wVal2*(I12[i+1][L-1] - I12[i+1][c] - I12[k][L-1] + I12[k][c]); 
					}
					/*I1324_full2, recursion: I1324(i;j,N)*/
                    I1324_full2[i+1][k] = I1324_full2_aid[i+1][k] + I1324_full2[i+1][k+1];
                    /*I1423_full2, recursion: I1423(i;j,N)*/
					I1423_full2[i+1][k] = I1423_full2_aid[i+1][k] + I1423_full2[i+1][k+1];

					/*order 3:*/
                    /*(12)*/
                    I123456[i][j] = I1234[j+1][L-1]*wValue + I123456[i+1][j] + I123456[i][j-1] - I123456[i+1][j-1];
                    I123645[i][j] = I1423[j+1][L-1]*wValue + I123645[i+1][j] + I123645[i][j-1] - I123645[i+1][j-1];
                    I123546[i][j] = I1324[j+1][L-1]*wValue+ I123546[i+1][j] + I123546[i][j-1] - I123546[i+1][j-1];
                    /*(13)*/
                    I132456[i][j] = (I1234[i+1][L-1] - I1234[i+1][j] -I1234[j][L-1])*wValue + I132456[i+1][j] + I132456[i][j-1] - I132456[i+1][j-1];
					I132546[i+1][j] = (I1324_full2[i+2][j+1] - I1324_full[j][L-1])*wVal[i+1][j] + I132546[i+2][j] + I132546[i+1][j-1] - I132546[i+2][j-1];
					I132645[i+1][j] = (I1423_full2[i+2][j+1] - I1423[j][L-1])*wVal[i+1][j] + I132645[i+2][j] + I132645[i+1][j-1] - I132645[i+2][j-1];
                    /*(14)*/
                    /*allowing direct computation:*/
                    I142356[i][j] = I12[i+1][j-1]*I12[j+1][L-1]*wValue + I142356[i+1][j] + I142356[i][j-1] - I142356[i+1][j-1];
                    /*back to the more complicated parts, now I142536:
                    three parts recognized in decompsition as already known:*/
                    I142536[i+1][j] = wVal[i+1][j]*(I1324_full[i+2][L-1] - I1324[i+2][j] - I1324_full2[i+2][j]);
                    /*recursion part:*/
                    I142536[i+1][j] += I142536[i+2][j] + I142536[i+1][j-1] - I142536[i+2][j-1];
                    /*Next complicated part, I142635:*/
                    /*three parts recognized in decompsition as already known:*/
					I142635[i+1][j] = wVal[i+1][j]*(I1423[i+2][L-1] - I1423_full2[i+2][j] - I1423_full0[i+2][j]);
                    /*recursion part:*/
                    I142635[i+1][j] += I142635[i+2][j] + I142635[i+1][j-1] - I142635[i+2][j-1];
                    /*(15)*/
                    I152346[i][j] = wValue*(I1234[i+1][j-1] - I1234_full[i+1][j] - I12[j][L-1]*I12[i+1][j-1]);
                    I152346[i][j] += I152346[i+1][j] + I152346[i][j-1] - I152346[i+1][j-1];
                    I152436[i][j] = wValue*(I1324[i+1][j-1] - I1324_full[i+1][j]);
                    I152436[i][j] += I152436[i+1][j] + I152436[i][j-1] - I152436[i+1][j-1];
                    I152634[i][j] = wValue*(I1423_full0[i+1][j-1] - I1423[i+1][j]);
                    I152634[i][j] += I152634[i+1][j] + I152634[i][j-1] - I152634[i+1][j-1];
                    /*(16)*/
                    I162345[i][j] = wValue*I1234_full[i+1][j-1];
                    I162345[i][j] += I162345[i+1][j] + I162345[i][j-1] - I162345[i+1][j-1];
                    I162435[i][j] = wValue*I1324_full[i+1][j-1];
                    I162435[i][j] += I162435[i+1][j] + I162435[i][j-1] - I162435[i+1][j-1];
                    I162534[i][j] = wValue*I1423[i+1][j-1];
                    I162534[i][j] += I162534[i+1][j] + I162534[i][j-1] - I162534[i+1][j-1];
			}
		}
		/*the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:*/
		j=L-1;
		for (j=L-1; j>= 0; j--){
			i = 0;
			c = j+1;
			for (c = j+1; c <= L-1; c++){
				wValue = wVal[j][c];
				I1324_full2_aid[i][j] += wValue*(I12[i][c-1] - I12[i][j] - I12[j][c-1]);
				I1423_full2_aid[i][j] += wValue*(I12[i][L-1] - I12[i][c] - I12[j][L-1] + I12[j][c]);
			}
			/*recursion:*/
			I1324_full2[i][j] = I1324_full2_aid[i][j] + I1324_full2[i][j+1];
			I1423_full2[i][j] = I1423_full2_aid[i][j] + I1423_full2[i][j+1]; 
		}


		/*#Missing terms for some of the degree 6 measures and I1423_full0:*/
		j= 1;
		for (j=1; j<= L-1; j++ ){
			i=0;
			wValue = wVal[i][j];
			I1423_full0[i][j] = (I12[i+1][j-1] - I12[i+1][j])*(I12[i+1][L-1] - I12[i+1][j] - I12[i][L-1] + I12[i][j]);
			I1423_full0[i][j] += I1423_full0[i+1][j] + I1423_full0[i][j-1] - I1423_full0[i+1][j-1];
			/*#*/
			I132546[i][j] = (I1324_full2[i+1][j+1] - I1324_full[j][L-1])*wValue + I132546[i+1][j] + I132546[i][j-1] - I132546[i+1][j-1];
			I132645[i][j] = (I1423_full2[i+1][j+1] - I1423[j][L-1])*wValue + I132645[i+1][j] + I132645[i][j-1] - I132645[i+1][j-1];
			/*#*/
			I142536[i][j] = wValue*(I1324_full[i+1][L-1] - I1324[i+1][j] - I1324_full2[i+1][j]);
			I142536[i][j] += I142536[i+1][j] + I142536[i][j-1] - I142536[i+1][j-1]; 
			/*#*/
			I142635[i][j] = wValue*(I1423[i+1][L-1] - I1423_full2[i+1][j] - I1423_full0[i+1][j]);
			I142635[i][j] += I142635[i+1][j] + I142635[i][j-1] - I142635[i+1][j-1];

		}
	}

	if (order !=0 && order != 1 && order != 2 && order != 3){ 
		printf("Order must be 0,1,2 or 3");
	}

	return 1; 
}

/*Both omputation of w-terms and aggregation/computation of invariants excluding the abolute value versions;
the aggregation is done by recursion traversing the invariants and the simplex:*/
int aggrAndW_ExAbs(struct segment *ptr_segment, int chainLen, int order, int full_b, struct I_ptr I_measures ){

	int L = 0;
	int i = 0;
	int j = 0;
	int b = 0;
	int k = 0;
	int c = 0;
	
	/* the parts of the output structure to hold the output:*/
	double wValue;
	double wVal1;
	double wVal2;

	/*next a series of pointers, all copies of the input/output I_measures.* pointers;
	they are merely there for convenience -- readability:*/

	double **wVal = I_measures.wVal;
	/*pointers for order 1 to 3 measures:*/
	double **I12 = I_measures.I12;
	
	double **I1234 = I_measures.I1234;
	double **I1324 = I_measures.I1324;
	double **I1423  = I_measures.I1423;

	/* for full computation of the order 2 measures across the simplex:*/
	double ** I1234_full_aid = I_measures.I1234_full_aid;
	double ** I1324_full_aid = I_measures.I1324_full_aid;
    double ** I1234_full = I_measures.I1234_full;
    double ** I1324_full = I_measures.I1324_full;
	/* assisting the order 3 case:*/
	double ** I1324_full2_aid = I_measures.I1324_full2_aid;
	double ** I1324_full2 = I_measures.I1324_full2;
	double **I1423_full0 = I_measures.I1423_full0;
	double **I1423_full2_aid = I_measures.I1423_full2_aid;
	double **I1423_full2 = I_measures.I1423_full2;

	/* for "absolute value version" of order 2 measures: */
	double **Ia1234 = I_measures.Ia1234;
	double **I12a34 = I_measures.I12a34;
	double **Ia1324 = I_measures.Ia1324;
	
	double **I13a24 = I_measures.I13a24;
	double **Ia1423 = I_measures.Ia1423;
	double **I14a23 = I_measures.I14a23;

	/*order 3*/
	/*12*/
	double **I123456 = I_measures.I123456;
	double **I123645 = I_measures.I123645;
	double **I123546 = I_measures.I123546;
	/*13*/
	double **I132456 = I_measures.I132456;
	double **I132546 = I_measures.I132546;
	double **I132645 = I_measures.I132645;
	/*14*/
	double **I142356 = I_measures.I142356;
	double **I142536 = I_measures.I142536;
	double **I142635 = I_measures.I142635;
	/*15*/
	double **I152346 = I_measures.I152346;
	double **I152436 = I_measures.I152436;
	double **I152634 = I_measures.I152634;
	/*16*/
	double **I162345 = I_measures.I162345;
	double **I162435 = I_measures.I162435;
	double **I162534 = I_measures.I162534;

	L = chainLen -1;
	
	I_measures.chainLen = &chainLen; 
	I_measures.order = &order;

	if (order == 0){

		for ( i=0; i<=L-1; i++ ){
			for ( j=i; j<=L-1; j++ ){ /* We include i = j for the purpose of checking that w[diagonal] =0. Only in order 0 */
					wValue = w_seg12(ptr_segment[i], ptr_segment[j]);
					wVal[i][j] = wValue;
			}
		}
	}

	if (order == 1){

		for ( i=L-1; i>= 0; i-- ){
			for ( j=i+1; j<=L-1; j++ ){
				wValue = w_seg12(ptr_segment[i], ptr_segment[j]);
				wVal[i][j] = wValue;
				I12[i][j] = wValue + I12[i][j-1] + I12[i+1][j] - I12[i+1][j-1];
			}
		}
	}


	if (order == 2){


		for ( i=L-1; i>= 0; i-- ){

			for ( j=i+1; j<=L-1; j++ ){

					wValue = w_seg12(ptr_segment[i], ptr_segment[j]);
					wVal[i][j] = wValue;
					I12[i][j] = wValue + I12[i][j-1] + I12[i+1][j] - I12[i+1][j-1];

					/*I1234: I1234(i,j;N)*/
                    I1234[i][j] = I12[j+1][L-1]*wValue + I1234[i][j-1] + I1234[i+1][j] - I1234[i+1][j-1];
                    /*I1423: I1423(i,j)*/
					I1423[i][j] = I12[i+1][j-1]*wValue +  I1423[i][j-1] + I1423[i+1][j] - I1423[i+1][j-1];
                    /*I1324: I1324(i,j;N)*/
                    I1324[i][j] = (I12[i+1][L-1] - I12[i+1][j] -I12[j][L-1] )*wValue + I1324[i][j-1] + I1324[i+1][j] - I1324[i+1][j-1];

					if (full_b == 1){
						b = i+1;
                        for (b = i+1; b <= j-2; b++){ /*#obs: the w-values needed have been computed*/
						 	wVal1 = wVal[i][b];
							I1234_full_aid[i][j] += I12[b+1][j]*wVal1;
                            I1324_full_aid[i][j] += I12[b][j]*wVal1;
						}
                        /*I1234_full add up terms and recursion: I1234(i,j)*/
                        I1234_full[i][j] = I1234_full_aid[i][j] -I1234_full_aid[i][j-1];
                        I1234_full[i][j] += I1234_full[i][j-1] + I1234_full[i+1][j] - I1234_full[i+1][j-1];
                        /*I1324_full add up terms recursion: I1324(i,j)*/
                        I1324_full[i][j] = I1324_full_aid[i][j-1] - I1324_full_aid[i][j];
                        I1324_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + I1324_full[i][j-1] + I1324_full[i+1][j] - I1324_full[i+1][j-1];
					}
			}
		}

	}




	if (order == 3){

		for ( i=L-1; i>= 0; i-- ){
			for ( j=i+1; j<=L-1; j++ ){
					wValue = w_seg12(ptr_segment[i], ptr_segment[j]);
					wVal[i][j] = wValue;
					I12[i][j] = wValue + I12[i][j-1] + I12[i+1][j] - I12[i+1][j-1];

					/*I1234: I1234(i,j;N)*/
                    I1234[i][j] = I12[j+1][L-1]*wValue + I1234[i][j-1] + I1234[i+1][j] - I1234[i+1][j-1];
                    /*I1423: I1423(i,j)*/
					I1423[i][j] = I12[i+1][j-1]*wValue +  I1423[i][j-1] + I1423[i+1][j] - I1423[i+1][j-1];
                    /*I1324: I1324(i,j;N)*/
                    I1324[i][j] = (I12[i+1][L-1] - I12[i+1][j] -I12[j][L-1] )*wValue + I1324[i][j-1] + I1324[i+1][j] - I1324[i+1][j-1];
					
					b = i+1;
                    for (b = i+1; b <= j-2; b++){ /*#obs: the w-values needed have been computed*/
						wVal1 = wVal[i][b];
						I1234_full_aid[i][j] += I12[b+1][j]*wVal1;
                        I1324_full_aid[i][j] += I12[b][j]*wVal1;
					}
                    /*I1234_full add up terms and recursion: I1234(i,j)*/
                    I1234_full[i][j] = I1234_full_aid[i][j] -I1234_full_aid[i][j-1];
                    I1234_full[i][j] += I1234_full[i][j-1] + I1234_full[i+1][j] - I1234_full[i+1][j-1];
                    /*I1324_full add up terms recursion: I1324(i,j)*/
                    I1324_full[i][j] = I1324_full_aid[i][j-1] - I1324_full_aid[i][j];
                    I1324_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + I1324_full[i][j-1] + I1324_full[i+1][j] - I1324_full[i+1][j-1];
					/*I1423_full: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet*/
                    I1423_full0[i+1][j] = (I12[i+2][j-1] - I12[i+2][j])*(I12[i+2][L-1] - I12[i+2][j] - I12[i+1][L-1] + I12[i+1][j]);
                    I1423_full0[i+1][j] += I1423_full0[i+2][j] + I1423_full0[i+1][j-1] - I1423_full0[i+2][j-1];
                    /*to compute certain degree 6 measures we use two auxillary degree 4 
                    measures; the recursion demands to sum j in the - direction (ie 
                    from above and down):*/
                    k = L-1 + i -j+1;
                    for (c =k+1;c <= L-1; c++){ /*Obs: the w values needed have already been computed (since k,c >i)*/
                        wVal2 = wVal[k][c];
                        I1324_full2_aid[i+1][k] += wVal2*(I12[i+1][c-1] - I12[i+1][k] - I12[k][c-1]);
                        I1423_full2_aid[i+1][k] += wVal2*(I12[i+1][L-1] - I12[i+1][c] - I12[k][L-1] + I12[k][c]); 
					}
					/*I1324_full2, recursion: I1324(i;j,N)*/
                    I1324_full2[i+1][k] = I1324_full2_aid[i+1][k] + I1324_full2[i+1][k+1];
                    /*I1423_full2, recursion: I1423(i;j,N)*/
					I1423_full2[i+1][k] = I1423_full2_aid[i+1][k] + I1423_full2[i+1][k+1];

					/*order 3:*/
                    /*(12)*/
                    I123456[i][j] = I1234[j+1][L-1]*wValue + I123456[i+1][j] + I123456[i][j-1] - I123456[i+1][j-1];
                    I123645[i][j] = I1423[j+1][L-1]*wValue + I123645[i+1][j] + I123645[i][j-1] - I123645[i+1][j-1];
                    I123546[i][j] = I1324[j+1][L-1]*wValue+ I123546[i+1][j] + I123546[i][j-1] - I123546[i+1][j-1];
                    /*(13)*/
                    I132456[i][j] = (I1234[i+1][L-1] - I1234[i+1][j] -I1234[j][L-1])*wValue + I132456[i+1][j] + I132456[i][j-1] - I132456[i+1][j-1];
					I132546[i+1][j] = (I1324_full2[i+2][j+1] - I1324_full[j][L-1])*wVal[i+1][j] + I132546[i+2][j] + I132546[i+1][j-1] - I132546[i+2][j-1];
					I132645[i+1][j] = (I1423_full2[i+2][j+1] - I1423[j][L-1])*wVal[i+1][j] + I132645[i+2][j] + I132645[i+1][j-1] - I132645[i+2][j-1];
                    /*(14)*/
                    /*allowing direct computation:*/
                    I142356[i][j] = I12[i+1][j-1]*I12[j+1][L-1]*wValue + I142356[i+1][j] + I142356[i][j-1] - I142356[i+1][j-1];
                    /*back to the more complicated parts, now I142536:
                    three parts recognized in decompsition as already known:*/
                    I142536[i+1][j] = wVal[i+1][j]*(I1324_full[i+2][L-1] - I1324[i+2][j] - I1324_full2[i+2][j]);
                    /*recursion part:*/
                    I142536[i+1][j] += I142536[i+2][j] + I142536[i+1][j-1] - I142536[i+2][j-1];
                    /*Next complicated part, I142635:*/
                    /*three parts recognized in decompsition as already known:*/
					I142635[i+1][j] = wVal[i+1][j]*(I1423[i+2][L-1] - I1423_full2[i+2][j] - I1423_full0[i+2][j]);
                    /*recursion part:*/
                    I142635[i+1][j] += I142635[i+2][j] + I142635[i+1][j-1] - I142635[i+2][j-1];
                    /*(15)*/
                    I152346[i][j] = wValue*(I1234[i+1][j-1] - I1234_full[i+1][j] - I12[j][L-1]*I12[i+1][j-1]);
                    I152346[i][j] += I152346[i+1][j] + I152346[i][j-1] - I152346[i+1][j-1];
                    I152436[i][j] = wValue*(I1324[i+1][j-1] - I1324_full[i+1][j]);
                    I152436[i][j] += I152436[i+1][j] + I152436[i][j-1] - I152436[i+1][j-1];
                    I152634[i][j] = wValue*(I1423_full0[i+1][j-1] - I1423[i+1][j]);
                    I152634[i][j] += I152634[i+1][j] + I152634[i][j-1] - I152634[i+1][j-1];
                    /*(16)*/
                    I162345[i][j] = wValue*I1234_full[i+1][j-1];
                    I162345[i][j] += I162345[i+1][j] + I162345[i][j-1] - I162345[i+1][j-1];
                    I162435[i][j] = wValue*I1324_full[i+1][j-1];
                    I162435[i][j] += I162435[i+1][j] + I162435[i][j-1] - I162435[i+1][j-1];
                    I162534[i][j] = wValue*I1423[i+1][j-1];
                    I162534[i][j] += I162534[i+1][j] + I162534[i][j-1] - I162534[i+1][j-1];
			}
		}
		/*the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:*/
		j=L-1;
		for (j=L-1; j>= 0; j--){
			i = 0;
			c = j+1;
			for (c = j+1; c <= L-1; c++){
				wValue = wVal[j][c];
				I1324_full2_aid[i][j] += wValue*(I12[i][c-1] - I12[i][j] - I12[j][c-1]);
				I1423_full2_aid[i][j] += wValue*(I12[i][L-1] - I12[i][c] - I12[j][L-1] + I12[j][c]);
			}
			/*recursion:*/
			I1324_full2[i][j] = I1324_full2_aid[i][j] + I1324_full2[i][j+1];
			I1423_full2[i][j] = I1423_full2_aid[i][j] + I1423_full2[i][j+1]; 
		}


		/*#Missing terms for some of the degree 6 measures and I1423_full0:*/
		j= 1;
		for (j=1; j<= L-1; j++ ){
			i=0;
			wValue = wVal[i][j];
			I1423_full0[i][j] = (I12[i+1][j-1] - I12[i+1][j])*(I12[i+1][L-1] - I12[i+1][j] - I12[i][L-1] + I12[i][j]);
			I1423_full0[i][j] += I1423_full0[i+1][j] + I1423_full0[i][j-1] - I1423_full0[i+1][j-1];
			/*#*/
			I132546[i][j] = (I1324_full2[i+1][j+1] - I1324_full[j][L-1])*wValue + I132546[i+1][j] + I132546[i][j-1] - I132546[i+1][j-1];
			I132645[i][j] = (I1423_full2[i+1][j+1] - I1423[j][L-1])*wValue + I132645[i+1][j] + I132645[i][j-1] - I132645[i+1][j-1];
			/*#*/
			I142536[i][j] = wValue*(I1324_full[i+1][L-1] - I1324[i+1][j] - I1324_full2[i+1][j]);
			I142536[i][j] += I142536[i+1][j] + I142536[i][j-1] - I142536[i+1][j-1]; 
			/*#*/
			I142635[i][j] = wValue*(I1423[i+1][L-1] - I1423_full2[i+1][j] - I1423_full0[i+1][j]);
			I142635[i][j] += I142635[i+1][j] + I142635[i][j-1] - I142635[i+1][j-1];

		}
	}

	if (order !=0 && order != 1 && order != 2 && order != 3){ 
		printf("Order must be 0,1,2 or 3");
	}

	return 1; 
}

/*Both computation of w-terms and aggregation/computation of invariants including the w-terms and the abolute 
value versions. The aggregation is done by recursion traversing the invariants and the simplex:*/
int aggrAndW(struct segment *ptr_segment, int chainLen, int order, int full_b, struct I_ptr I_measures ){

	int L = 0;
	int i = 0;
	int j = 0;
	int b = 0;
	int k = 0;
	int c = 0;
	
	/* the parts of the output structure to hold the output:*/
	double wValue;
	double abs_wValue; 
	double wVal1;
	double abs_wVal1; 
	double wVal2;

	/*next a series of pointers, all copies of the input/output I_measures.* pointers;
	they are merely there for convenience -- readaility:*/

	double **wVal = I_measures.wVal;
	/*pointers for order 1 to 3 measures:*/
	double **I12 = I_measures.I12;
	double **Ia12 = I_measures.Ia12;

	double **I1234 = I_measures.I1234;
	double **I1324 = I_measures.I1324;
	double **I1423  = I_measures.I1423;

	/* for "absolute value version" of order 2 measures: */
	double **Ia1234 = I_measures.Ia1234;
	double **I12a34 = I_measures.I12a34;
	double **Ia12a34 = I_measures.Ia12a34;

	double **Ia1324 = I_measures.Ia1324;
	double **I13a24 = I_measures.I13a24;
	double **Ia13a24 = I_measures.Ia13a24;

	double **Ia1423  = I_measures.Ia1423;
	double **I14a23  = I_measures.I14a23;
	double **Ia14a23  = I_measures.Ia14a23;


	/* for full computation of the order 2 measures across the simplex:*/
	double ** I1234_full_aid = I_measures.I1234_full_aid;
	double ** I1324_full_aid = I_measures.I1324_full_aid;
    double ** I1234_full = I_measures.I1234_full;
    double ** I1324_full = I_measures.I1324_full;
	/* for full computation of the "absolute value version" order 2 measures across the simplex:*/
	double ** Ia12a34_full_aid = I_measures.Ia12a34_full_aid;
	double ** Ia1234_full_aid = I_measures.Ia1234_full_aid;
	double ** I12a34_full_aid = I_measures.I12a34_full_aid;

    double ** Ia12a34_full = I_measures.Ia12a34_full;
	double ** Ia1234_full = I_measures.Ia1234_full;
	double ** I12a34_full = I_measures.I12a34_full;

	double ** Ia13a24_full_aid = I_measures.Ia13a24_full_aid;
	double ** Ia1324_full_aid = I_measures.Ia1324_full_aid;
	double ** I13a24_full_aid = I_measures.I13a24_full_aid;

    double ** Ia13a24_full = I_measures.Ia13a24_full;
	double ** Ia1324_full = I_measures.Ia1324_full;
	double ** I13a24_full = I_measures.I13a24_full;

	/* assisting the order 3 case:*/
	double ** I1324_full2_aid = I_measures.I1324_full2_aid;
	double ** I1324_full2 = I_measures.I1324_full2;
	double **I1423_full0 = I_measures.I1423_full0;
	double **I1423_full2_aid = I_measures.I1423_full2_aid;
	double **I1423_full2 = I_measures.I1423_full2;

	/*order 3*/
	/*12*/
	double **I123456 = I_measures.I123456;
	double **I123645 = I_measures.I123645;
	double **I123546 = I_measures.I123546;
	/*13*/
	double **I132456 = I_measures.I132456;
	double **I132546 = I_measures.I132546;
	double **I132645 = I_measures.I132645;
	/*14*/
	double **I142356 = I_measures.I142356;
	double **I142536 = I_measures.I142536;
	double **I142635 = I_measures.I142635;
	/*15*/
	double **I152346 = I_measures.I152346;
	double **I152436 = I_measures.I152436;
	double **I152634 = I_measures.I152634;
	/*16*/
	double **I162345 = I_measures.I162345;
	double **I162435 = I_measures.I162435;
	double **I162534 = I_measures.I162534;

	L = chainLen -1;

	I_measures.chainLen = &chainLen; 
	I_measures.order = &order;

	if (order == 0){

		for ( i=0; i<=L-1; i++ ){
			for ( j=i; j<=L-1; j++ ){ /* We include i = j for the purpose of checking that w[diagonal] =0. Only in order 0 */
					wValue = w_seg12(ptr_segment[i], ptr_segment[j]);
					/*sscanf(wValue,"%lf",&wVal[cnt]); */
					wVal[i][j] = wValue;
			}
		}
	}

	if (order == 1){

		for ( i=L-1; i>= 0; i-- ){
			for ( j=i+1; j<=L-1; j++ ){
				wValue = w_seg12(ptr_segment[i], ptr_segment[j]);
				/*sscanf(wValue,"%lf",&wVal[cnt]); */
				wVal[i][j] = wValue;
				I12[i][j] = wValue + I12[i][j-1] + I12[i+1][j] - I12[i+1][j-1];
				Ia12[i][j] = fabsf(wValue) + Ia12[i][j-1] + Ia12[i+1][j] - Ia12[i+1][j-1];

			}
		}
	}


	if (order == 2){

		for ( i=L-1; i>= 0; i-- ){
			for ( j=i+1; j<=L-1; j++ ){
					wValue = w_seg12(ptr_segment[i], ptr_segment[j]);
					/*sscanf(wValue,"%lf",&wVal[cnt]); */
					wVal[i][j] = wValue;
					abs_wValue = fabsf(wValue);
					I12[i][j] = wValue + I12[i][j-1] + I12[i+1][j] - I12[i+1][j-1];
					Ia12[i][j] = abs_wValue + Ia12[i][j-1] + Ia12[i+1][j] - Ia12[i+1][j-1];

					/*I1234: I1234(i,j;N)*/
                    I1234[i][j] = I12[j+1][L-1]*wValue + I1234[i][j-1] + I1234[i+1][j] - I1234[i+1][j-1];
                    /*I1423: I1423(i,j)*/
					I1423[i][j] = I12[i+1][j-1]*wValue +  I1423[i][j-1] + I1423[i+1][j] - I1423[i+1][j-1];
                    /*I1324: I1324(i,j;N)*/
                    I1324[i][j] = (I12[i+1][L-1] - I12[i+1][j] -I12[j][L-1] )*wValue + I1324[i][j-1] + I1324[i+1][j] - I1324[i+1][j-1];

					/*absolute value versions:*/
                    Ia1234[i][j] = I12[j+1][L-1]*abs_wValue + Ia1234[i][j-1] + Ia1234[i+1][j] - Ia1234[i+1][j-1];
                    I12a34[i][j] = Ia12[j+1][L-1]*wValue + I12a34[i][j-1] + I12a34[i+1][j] - I12a34[i+1][j-1];
                    Ia12a34[i][j] = Ia12[j+1][L-1]*abs_wValue + Ia12a34[i][j-1] + Ia12a34[i+1][j] - Ia12a34[i+1][j-1];

					Ia1324[i][j] = (I12[i+1][L-1] - I12[i+1][j] -I12[j][L-1] )*abs_wValue + Ia1324[i][j-1] + Ia1324[i+1][j] - Ia1324[i+1][j-1];
					I13a24[i][j] = (Ia12[i+1][L-1] - Ia12[i+1][j] -Ia12[j][L-1] )*wValue + I13a24[i][j-1] + I13a24[i+1][j] - I13a24[i+1][j-1];
					Ia13a24[i][j] = (Ia12[i+1][L-1] - Ia12[i+1][j] -Ia12[j][L-1] )*abs_wValue + Ia13a24[i][j-1] + Ia13a24[i+1][j] - Ia13a24[i+1][j-1];

					Ia1423[i][j] = I12[i+1][j-1]*abs_wValue +  Ia1423[i][j-1] + Ia1423[i+1][j] - Ia1423[i+1][j-1];
					I14a23[i][j] = Ia12[i+1][j-1]*wValue +  I14a23[i][j-1] + I14a23[i+1][j] - I14a23[i+1][j-1];
					Ia14a23[i][j] = Ia12[i+1][j-1]*abs_wValue +  Ia14a23[i][j-1] + Ia14a23[i+1][j] - Ia14a23[i+1][j-1];


					if (full_b == 1){
						b = i+1;
                        for (b = i+1; b <= j-2; b++){ /*#obs: the w-values needed have been computed*/
						 	wVal1 = wVal[i][b];
							abs_wVal1 = fabsf(wVal1);
							I1234_full_aid[i][j] += I12[b+1][j]*wVal1;
                            I1324_full_aid[i][j] += I12[b][j]*wVal1;

							Ia12a34_full_aid[i][j] += Ia12[b+1][j]*abs_wVal1;
							Ia1234_full_aid[i][j] += Ia12[b+1][j]*wVal1;
							I12a34_full_aid[i][j] += I12[b+1][j]*abs_wVal1;

							Ia13a24_full_aid[i][j] += I12[b][j]*abs_wVal1;
							Ia1324_full_aid[i][j] += Ia12[b][j]*wVal1;
							I13a24_full_aid[i][j] += I12[b][j]*abs_wVal1;

						}
                        /*I1234_full add up terms and recursion: I1234(i,j)*/
                        I1234_full[i][j] = I1234_full_aid[i][j] -I1234_full_aid[i][j-1];
                        I1234_full[i][j] += I1234_full[i][j-1] + I1234_full[i+1][j] - I1234_full[i+1][j-1];
                        /*I1324_full add up terms recursion: I1324(i,j)*/
                        I1324_full[i][j] = I1324_full_aid[i][j-1] - I1324_full_aid[i][j];
                        I1324_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + I1324_full[i][j-1] + I1324_full[i+1][j] - I1324_full[i+1][j-1];

						/*abs value versions:*/
						Ia12a34_full[i][j] = Ia12a34_full_aid[i][j] - Ia12a34_full_aid[i][j-1];
						Ia12a34_full[i][j] += Ia12a34_full[i][j-1] + Ia12a34_full[i+1][j] - Ia12a34_full[i+1][j-1];

						Ia1234_full[i][j] = Ia1234_full_aid[i][j] - Ia1234_full_aid[i][j-1];
						Ia1234_full[i][j] += Ia1234_full[i][j-1] + Ia1234_full[i+1][j] - Ia1234_full[i+1][j-1];

						I12a34_full[i][j] = I12a34_full_aid[i][j] -I12a34_full_aid[i][j-1];
						I12a34_full[i][j] += I12a34_full[i][j-1] + I12a34_full[i+1][j] - I12a34_full[i+1][j-1];

						Ia13a24_full[i][j] = Ia13a24_full_aid[i][j-1] - Ia13a24_full_aid[i][j];
						Ia13a24_full[i][j] += (Ia12[i][j-1] - Ia12[i+1][j-1])*(Ia12[i+1][j] - Ia12[i+1][j-1]) + Ia13a24_full[i][j-1] + Ia13a24_full[i+1][j] - Ia13a24_full[i+1][j-1];

						Ia1324_full[i][j] = Ia1324_full_aid[i][j-1] - Ia1324_full_aid[i][j];
						Ia1324_full[i][j] += (Ia12[i][j-1] - Ia12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + Ia1324_full[i][j-1] + Ia1324_full[i+1][j] - Ia1324_full[i+1][j-1];

						I13a24_full[i][j] = I13a24_full_aid[i][j-1] - I13a24_full_aid[i][j];
						I13a24_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(Ia12[i+1][j] - Ia12[i+1][j-1]) + I13a24_full[i][j-1] + I13a24_full[i+1][j] - I13a24_full[i+1][j-1];


					}
			}
		}
	}




	if (order == 3){

		for ( i=L-1; i>= 0; i-- ){
			for ( j=i+1; j<=L-1; j++ ){
					wValue = w_seg12(ptr_segment[i], ptr_segment[j]);
					/*sscanf(wValue,"%lf",&wVal[cnt]); */
					wVal[i][j] = wValue;
					abs_wValue = fabsf(wValue);
					I12[i][j] = wValue + I12[i][j-1] + I12[i+1][j] - I12[i+1][j-1];
					Ia12[i][j] = abs_wValue + Ia12[i][j-1] + Ia12[i+1][j] - Ia12[i+1][j-1];

					/*I1234: I1234(i,j;N)*/
                    I1234[i][j] = I12[j+1][L-1]*wValue + I1234[i][j-1] + I1234[i+1][j] - I1234[i+1][j-1];
                    /*I1423: I1423(i,j)*/
					I1423[i][j] = I12[i+1][j-1]*wValue +  I1423[i][j-1] + I1423[i+1][j] - I1423[i+1][j-1];
                    /*I1324: I1324(i,j;N)*/
                    I1324[i][j] = (I12[i+1][L-1] - I12[i+1][j] -I12[j][L-1] )*wValue + I1324[i][j-1] + I1324[i+1][j] - I1324[i+1][j-1];
					
					/*absolute value versions:*/
                    Ia1234[i][j] = I12[j+1][L-1]*abs_wValue + Ia1234[i][j-1] + Ia1234[i+1][j] - Ia1234[i+1][j-1];
                    I12a34[i][j] = Ia12[j+1][L-1]*wValue + I12a34[i][j-1] + I12a34[i+1][j] - I12a34[i+1][j-1];
                    Ia12a34[i][j] = Ia12[j+1][L-1]*abs_wValue + Ia12a34[i][j-1] + Ia12a34[i+1][j] - Ia12a34[i+1][j-1];

					Ia1324[i][j] = (I12[i+1][L-1] - I12[i+1][j] -I12[j][L-1] )*abs_wValue + Ia1324[i][j-1] + Ia1324[i+1][j] - Ia1324[i+1][j-1];
					I13a24[i][j] = (Ia12[i+1][L-1] - Ia12[i+1][j] -Ia12[j][L-1] )*wValue + I13a24[i][j-1] + I13a24[i+1][j] - I13a24[i+1][j-1];
					Ia13a24[i][j] = (Ia12[i+1][L-1] - Ia12[i+1][j] -Ia12[j][L-1] )*abs_wValue + Ia13a24[i][j-1] + Ia13a24[i+1][j] - Ia13a24[i+1][j-1];

					Ia1423[i][j] = I12[i+1][j-1]*abs_wValue +  Ia1423[i][j-1] + Ia1423[i+1][j] - Ia1423[i+1][j-1];
					I14a23[i][j] = Ia12[i+1][j-1]*wValue +  I14a23[i][j-1] + I14a23[i+1][j] - I14a23[i+1][j-1];
					Ia14a23[i][j] = Ia12[i+1][j-1]*abs_wValue +  Ia14a23[i][j-1] + Ia14a23[i+1][j] - Ia14a23[i+1][j-1];


					b = i+1;
                    for (b = i+1; b <= j-2; b++){ /*#obs: the w-values needed have been computed*/
						wVal1 = wVal[i][b];
						abs_wVal1 = fabsf(wVal1);
						I1234_full_aid[i][j] += I12[b+1][j]*wVal1;
                        I1324_full_aid[i][j] += I12[b][j]*wVal1;

						Ia12a34_full_aid[i][j] += Ia12[b+1][j]*abs_wVal1;
						Ia1234_full_aid[i][j] += Ia12[b+1][j]*wVal1;
						I12a34_full_aid[i][j] += I12[b+1][j]*abs_wVal1;

						Ia13a24_full_aid[i][j] += I12[b][j]*abs_wVal1;
						Ia1324_full_aid[i][j] += Ia12[b][j]*wVal1;
						I13a24_full_aid[i][j] += I12[b][j]*abs_wVal1;

					}
                    /*I1234_full add up terms and recursion: I1234(i,j)*/
                    I1234_full[i][j] = I1234_full_aid[i][j] -I1234_full_aid[i][j-1];
                    I1234_full[i][j] += I1234_full[i][j-1] + I1234_full[i+1][j] - I1234_full[i+1][j-1];
                    /*I1324_full add up terms recursion: I1324(i,j)*/
                    I1324_full[i][j] = I1324_full_aid[i][j-1] - I1324_full_aid[i][j];
                    I1324_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + I1324_full[i][j-1] + I1324_full[i+1][j] - I1324_full[i+1][j-1];
					
					/*abs value versions:*/
					Ia12a34_full[i][j] = Ia12a34_full_aid[i][j] - Ia12a34_full_aid[i][j-1];
                    Ia12a34_full[i][j] += Ia12a34_full[i][j-1] + Ia12a34_full[i+1][j] - Ia12a34_full[i+1][j-1];

					Ia1234_full[i][j] = Ia1234_full_aid[i][j] - Ia1234_full_aid[i][j-1];
                    Ia1234_full[i][j] += Ia1234_full[i][j-1] + Ia1234_full[i+1][j] - Ia1234_full[i+1][j-1];

					I12a34_full[i][j] = I12a34_full_aid[i][j] -I12a34_full_aid[i][j-1];
                    I12a34_full[i][j] += I12a34_full[i][j-1] + I12a34_full[i+1][j] - I12a34_full[i+1][j-1];

					Ia13a24_full[i][j] = Ia13a24_full_aid[i][j-1] - Ia13a24_full_aid[i][j];
                    Ia13a24_full[i][j] += (Ia12[i][j-1] - Ia12[i+1][j-1])*(Ia12[i+1][j] - Ia12[i+1][j-1]) + Ia13a24_full[i][j-1] + Ia13a24_full[i+1][j] - Ia13a24_full[i+1][j-1];

					Ia1324_full[i][j] = Ia1324_full_aid[i][j-1] - Ia1324_full_aid[i][j];
                    Ia1324_full[i][j] += (Ia12[i][j-1] - Ia12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + Ia1324_full[i][j-1] + Ia1324_full[i+1][j] - Ia1324_full[i+1][j-1];

					I13a24_full[i][j] = I13a24_full_aid[i][j-1] - I13a24_full_aid[i][j];
                    I13a24_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(Ia12[i+1][j] - Ia12[i+1][j-1]) + I13a24_full[i][j-1] + I13a24_full[i+1][j] - I13a24_full[i+1][j-1];

					
					/*I1423_full: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet*/
                    I1423_full0[i+1][j] = (I12[i+2][j-1] - I12[i+2][j])*(I12[i+2][L-1] - I12[i+2][j] - I12[i+1][L-1] + I12[i+1][j]);
                    I1423_full0[i+1][j] += I1423_full0[i+2][j] + I1423_full0[i+1][j-1] - I1423_full0[i+2][j-1];
                    /*to compute certain degree 6 measures we use two auxillary degree 4 
                    measures; the recursion demands to sum j in the - direction (ie 
                    from above and down):*/
                    k = L-1 + i -j+1;
                    for (c =k+1;c <= L-1; c++){ /*Obs: the w values needed have already been computed (since k,c >i)*/
                        wVal2 = wVal[k][c];
                        I1324_full2_aid[i+1][k] += wVal2*(I12[i+1][c-1] - I12[i+1][k] - I12[k][c-1]);
                        I1423_full2_aid[i+1][k] += wVal2*(I12[i+1][L-1] - I12[i+1][c] - I12[k][L-1] + I12[k][c]); 
					}
					/*I1324_full2, recursion: I1324(i;j,N)*/
                    I1324_full2[i+1][k] = I1324_full2_aid[i+1][k] + I1324_full2[i+1][k+1];
                    /*I1423_full2, recursion: I1423(i;j,N)*/
					I1423_full2[i+1][k] = I1423_full2_aid[i+1][k] + I1423_full2[i+1][k+1];

					/*order 3:*/
                    /*(12)*/
                    I123456[i][j] = I1234[j+1][L-1]*wValue + I123456[i+1][j] + I123456[i][j-1] - I123456[i+1][j-1];
                    I123645[i][j] = I1423[j+1][L-1]*wValue + I123645[i+1][j] + I123645[i][j-1] - I123645[i+1][j-1];
                    I123546[i][j] = I1324[j+1][L-1]*wValue+ I123546[i+1][j] + I123546[i][j-1] - I123546[i+1][j-1];
                    /*(13)*/
                    I132456[i][j] = (I1234[i+1][L-1] - I1234[i+1][j] -I1234[j][L-1])*wValue + I132456[i+1][j] + I132456[i][j-1] - I132456[i+1][j-1];
					I132546[i+1][j] = (I1324_full2[i+2][j+1] - I1324_full[j][L-1])*wVal[i+1][j] + I132546[i+2][j] + I132546[i+1][j-1] - I132546[i+2][j-1];
					I132645[i+1][j] = (I1423_full2[i+2][j+1] - I1423[j][L-1])*wVal[i+1][j] + I132645[i+2][j] + I132645[i+1][j-1] - I132645[i+2][j-1];
                    /*(14)*/
                    /*allowing direct computation:*/
                    I142356[i][j] = I12[i+1][j-1]*I12[j+1][L-1]*wValue + I142356[i+1][j] + I142356[i][j-1] - I142356[i+1][j-1];
                    /*back to the more complicated parts, now I142536:
                    three parts recognized in decompsition as already known:*/
                    I142536[i+1][j] = wVal[i+1][j]*(I1324_full[i+2][L-1] - I1324[i+2][j] - I1324_full2[i+2][j]);
                    /*recursion part:*/
                    I142536[i+1][j] += I142536[i+2][j] + I142536[i+1][j-1] - I142536[i+2][j-1];
                    /*Next complicated part, I142635:*/
                    /*three parts recognized in decompsition as already known:*/
					I142635[i+1][j] = wVal[i+1][j]*(I1423[i+2][L-1] - I1423_full2[i+2][j] - I1423_full0[i+2][j]);
                    /*recursion part:*/
                    I142635[i+1][j] += I142635[i+2][j] + I142635[i+1][j-1] - I142635[i+2][j-1];
                    /*(15)*/
                    I152346[i][j] = wValue*(I1234[i+1][j-1] - I1234_full[i+1][j] - I12[j][L-1]*I12[i+1][j-1]);
                    I152346[i][j] += I152346[i+1][j] + I152346[i][j-1] - I152346[i+1][j-1];
                    I152436[i][j] = wValue*(I1324[i+1][j-1] - I1324_full[i+1][j]);
                    I152436[i][j] += I152436[i+1][j] + I152436[i][j-1] - I152436[i+1][j-1];
                    I152634[i][j] = wValue*(I1423_full0[i+1][j-1] - I1423[i+1][j]);
                    I152634[i][j] += I152634[i+1][j] + I152634[i][j-1] - I152634[i+1][j-1];
                    /*(16)*/
                    I162345[i][j] = wValue*I1234_full[i+1][j-1];
                    I162345[i][j] += I162345[i+1][j] + I162345[i][j-1] - I162345[i+1][j-1];
                    I162435[i][j] = wValue*I1324_full[i+1][j-1];
                    I162435[i][j] += I162435[i+1][j] + I162435[i][j-1] - I162435[i+1][j-1];
                    I162534[i][j] = wValue*I1423[i+1][j-1];
                    I162534[i][j] += I162534[i+1][j] + I162534[i][j-1] - I162534[i+1][j-1];
			}
		}
		/*the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:*/
		j=L-1;
		for (j=L-1; j>= 0; j--){
			i = 0;
			c = j+1;
			for (c = j+1; c <= L-1; c++){
				wValue = wVal[j][c];
				I1324_full2_aid[i][j] += wValue*(I12[i][c-1] - I12[i][j] - I12[j][c-1]);
				I1423_full2_aid[i][j] += wValue*(I12[i][L-1] - I12[i][c] - I12[j][L-1] + I12[j][c]);
			}
			/*recursion:*/
			I1324_full2[i][j] = I1324_full2_aid[i][j] + I1324_full2[i][j+1];
			I1423_full2[i][j] = I1423_full2_aid[i][j] + I1423_full2[i][j+1]; 
		}


		/*#Missing terms for some of the degree 6 measures and I1423_full0:*/
		j = 1;
		for (j=1; j<= L-1; j++ ){
			i=0;
			wValue = wVal[i][j];
			I1423_full0[i][j] = (I12[i+1][j-1] - I12[i+1][j])*(I12[i+1][L-1] - I12[i+1][j] - I12[i][L-1] + I12[i][j]);
			I1423_full0[i][j] += I1423_full0[i+1][j] + I1423_full0[i][j-1] - I1423_full0[i+1][j-1];
			/*#*/
			I132546[i][j] = (I1324_full2[i+1][j+1] - I1324_full[j][L-1])*wValue + I132546[i+1][j] + I132546[i][j-1] - I132546[i+1][j-1];
			I132645[i][j] = (I1423_full2[i+1][j+1] - I1423[j][L-1])*wValue + I132645[i+1][j] + I132645[i][j-1] - I132645[i+1][j-1];
			/*#*/
			I142536[i][j] = wValue*(I1324_full[i+1][L-1] - I1324[i+1][j] - I1324_full2[i+1][j]);
			I142536[i][j] += I142536[i+1][j] + I142536[i][j-1] - I142536[i+1][j-1]; 
			/*#*/
			I142635[i][j] = wValue*(I1423[i+1][L-1] - I1423_full2[i+1][j] - I1423_full0[i+1][j]);
			I142635[i][j] += I142635[i+1][j] + I142635[i][j-1] - I142635[i+1][j-1];

		}
	}

	if (order !=0 && order != 1 && order != 2 && order != 3){ 
		printf("Order must be 0,1,2 or 3");
	}

	return 1; 
}


/* As aggrAndW but with a search for closed loops added: inside the traversal of the simplex (taking part of the 
recursion for computing the invariants) a 1d-loop is inserted traversing indices above the current "i" to search 
for the first place where a closed-loop distance criterium is fulfilled (if at any index). We neglect loops in 
which a segment has "too long" length ie which consists in an artificial peptide bond introduced by the way we 
load structures (as one long chain disregarding if it consists of several sub-chains (this can be avoided if 
changing the load procedure). A standard length of sqare root of 20 is used for this purpose. Also we neglect 
loops which are "too short": the minLoopLength allows varying this.*/
int aggrAndW_wClosedLoops(struct segment *ptr_segment, int chainLen, int order, int full_b, struct I_ptr I_measures, struct twoSegmentIndex *ptr_closedLoopInd, int closedLoopLength, double closedLoopDist, int * ptr_nrOfClosedLoops){

	int L = 0;
	int i = 0;
	int j = 0;
	int b = 0;
	int k = 0;
	int c = 0;
	int minLoopLength = 7;
	int loopTooLongInd = 0;
	//double stdRealSegLength = 20; /*close to square of 4.5*/
	int n = 0;
	int last_j = 0;
	
	int cnt = 0; /*for counting and indexing closed loops*/ 

	/* the parts of the output structure to hold the output:*/
	double wValue;
	double abs_wValue; 
	double wVal1;
	double abs_wVal1;
	double wVal2;

	/*next a series of pointers, all copies of the input/output I_measures.* pointers;
	they are merely there for convenience -- readability:*/

	double **wVal = I_measures.wVal;
	/*pointers for order 1 to 3 measures:*/
	double **I12 = I_measures.I12;
	double **Ia12 = I_measures.Ia12;

	double **I1234 = I_measures.I1234;
	double **I1324 = I_measures.I1324;
	double **I1423  = I_measures.I1423;

	/* for "absolute value version" of order 2 measures: */
	double **Ia1234 = I_measures.Ia1234;
	double **I12a34 = I_measures.I12a34;
	double **Ia12a34 = I_measures.Ia12a34;

	double **Ia1324 = I_measures.Ia1324;
	double **I13a24 = I_measures.I13a24;
	double **Ia13a24 = I_measures.Ia13a24;

	double **Ia1423  = I_measures.Ia1423;
	double **I14a23  = I_measures.I14a23;
	double **Ia14a23  = I_measures.Ia14a23;


	/* for full computation of the order 2 measures across the simplex:*/
	double ** I1234_full_aid = I_measures.I1234_full_aid;
	double ** I1324_full_aid = I_measures.I1324_full_aid;
    double ** I1234_full = I_measures.I1234_full;
    double ** I1324_full = I_measures.I1324_full;

	/* for full computation of the "absolute value version" order 2 measures across the simplex:*/
	double ** Ia12a34_full_aid = I_measures.Ia12a34_full_aid;
	double ** Ia1234_full_aid = I_measures.Ia1234_full_aid;
	double ** I12a34_full_aid = I_measures.I12a34_full_aid;

    double ** Ia12a34_full = I_measures.Ia12a34_full;
	double ** Ia1234_full = I_measures.Ia1234_full;
	double ** I12a34_full = I_measures.I12a34_full;

	double ** Ia13a24_full_aid = I_measures.Ia13a24_full_aid;
	double ** Ia1324_full_aid = I_measures.Ia1324_full_aid;
	double ** I13a24_full_aid = I_measures.I13a24_full_aid;

    double ** Ia13a24_full = I_measures.Ia13a24_full;
	double ** Ia1324_full = I_measures.Ia1324_full;
	double ** I13a24_full = I_measures.I13a24_full;

	/* assisting the order 3 case:*/
	double ** I1324_full2_aid = I_measures.I1324_full2_aid;
	double ** I1324_full2 = I_measures.I1324_full2;
	double **I1423_full0 = I_measures.I1423_full0;
	double **I1423_full2_aid = I_measures.I1423_full2_aid;
	double **I1423_full2 = I_measures.I1423_full2;

	/*order 3*/
	/*12*/
	double **I123456 = I_measures.I123456;
	double **I123645 = I_measures.I123645;
	double **I123546 = I_measures.I123546;
	/*13*/
	double **I132456 = I_measures.I132456;
	double **I132546 = I_measures.I132546;
	double **I132645 = I_measures.I132645;
	/*14*/
	double **I142356 = I_measures.I142356;
	double **I142536 = I_measures.I142536;
	double **I142635 = I_measures.I142635;
	/*15*/
	double **I152346 = I_measures.I152346;
	double **I152436 = I_measures.I152436;
	double **I152634 = I_measures.I152634;
	/*16*/
	double **I162345 = I_measures.I162345;
	double **I162435 = I_measures.I162435;
	double **I162534 = I_measures.I162534;

	L = chainLen -1;

	I_measures.chainLen = &chainLen; 
	I_measures.order = &order;

	if (order == 0){

		printf("Closed loops facility only works with order >= 1. The function returns without effect.");
		
	}

	last_j = L; 

	if (order == 1){

		for ( i=L-1; i>= 0; i-- ){
			
			/*find possible closed loops starting at first coord of segment[i]
			We discard loops of minLoopLength by letting n run from i+minLoopLength. 
			Also, we want to find the "longest irreducible loop". Since we reverse 
			in the i-loop, we try to extend a loop by adding to its "i head": if we 
			for i have found a loop terminating at n, we will not search further at 
			i (because of the last_j = n) and in the next step (at i-1) we will ultimately 
			check if (i-1,n-1) is also closed (but this will only happen if we have not hit 
			a loop (i-1,n') for some n' lower than n*/
			for (n=i+minLoopLength; n <= min(i+ closedLoopLength,last_j-1); n++){ /*should it be last_j rather than last_j -1?*/
				/*if segment n is "too long" we skip the rest of the search at this i:*/
				if (distCalphas(ptr_segment[n].s1, ptr_segment[n].s2) > stdRealSegLength){
					break;
				} 
				/* if(i>=10 && i < 20){
					printf("i:%d n:%d last_j:%d\n", i,n, last_j);
					printf("At i:%d n:%d dist is: %lf\n", i,n, distCalphas(ptr_segment[i].s1, ptr_segment[n].s2));
					getchar();
				}*/
				if (distCalphas(ptr_segment[i].s1, ptr_segment[n].s2) < closedLoopDist){
					ptr_closedLoopInd[cnt].idx1 = i;
					ptr_closedLoopInd[cnt].idx2 = n;
					/* printf("Seg i: %d 1st pt is x y z : %lf %lf %lf\n", i, ptr_segment[i].s1.x,ptr_segment[i].s1.y, ptr_segment[i].s1.z );
					printf("Seg n: %d 2nd pt is x y z : %lf %lf %lf\n", n, ptr_segment[n].s2.x,ptr_segment[n].s2.y, ptr_segment[n].s2.z );
					printf("At i:%d n:%d dist is: %lf\n", i,n, distCalphas(ptr_segment[i].s1, ptr_segment[n].s2));*/
					last_j = n; /*for getting only the irreducible closed loops. Obs: we reverse in the i-loop!*/
					cnt +=1;
				}
			}
			
			for ( j=i+1; j<=L-1; j++ ){
				wValue = w_seg12(ptr_segment[i], ptr_segment[j]);
				wVal[i][j] = wValue;
				I12[i][j] = wValue + I12[i][j-1] + I12[i+1][j] - I12[i+1][j-1];
				Ia12[i][j] = fabsf(wValue) + Ia12[i][j-1] + Ia12[i+1][j] - Ia12[i+1][j-1];

			}
		}
	}


	if (order == 2){

		for ( i=L-1; i>= 0; i-- ){
			
			/*find possible closed loops starting at first coord of segment[i]
			We discard "length three loops" by letting n run from i+4. Also, we
			want to find the "longest "irreducible loop". Note that
			since we reverse in the i-loop, so we try to extend a loop by
			adding to its "i head": if we for i have found a loop terminating
			at n, we will not search further at i (becasue of the last_j = n)
			and in the next step we will ultimately check if (i-1,n) is also
			closed (but this will only happen if we have not hit a loop 
			(i-1,n') for some n' lower than n*/
			for (n=i+minLoopLength; n <= min(i+ closedLoopLength,last_j-1); n++){
				/*if segment n is "too long" we skip the rest of the search at this i:*/
				if (distCalphas(ptr_segment[n].s1, ptr_segment[n].s2) > stdRealSegLength){
					break;
				}
				//printf("i:%d n:%d last_j:%d\n", i,n, last_j);
				if (distCalphas(ptr_segment[i].s1, ptr_segment[n].s2) < closedLoopDist){
					ptr_closedLoopInd[cnt].idx1 = i;
					ptr_closedLoopInd[cnt].idx2 = n;
					last_j = n; /*for getting only the irreducible closed loops. Obs: we reverse in the i-loop!*/
					cnt +=1;
				}
			}
			
			for ( j=i+1; j<=L-1; j++ ){
					wValue = w_seg12(ptr_segment[i], ptr_segment[j]);
					wVal[i][j] = wValue;
					abs_wValue = fabsf(wValue);
					I12[i][j] = wValue + I12[i][j-1] + I12[i+1][j] - I12[i+1][j-1];
					Ia12[i][j] = abs_wValue + Ia12[i][j-1] + Ia12[i+1][j] - Ia12[i+1][j-1];

					/*I1234: I1234(i,j;N)*/
                    I1234[i][j] = I12[j+1][L-1]*wValue + I1234[i][j-1] + I1234[i+1][j] - I1234[i+1][j-1];
                    /*I1423: I1423(i,j)*/
					I1423[i][j] = I12[i+1][j-1]*wValue +  I1423[i][j-1] + I1423[i+1][j] - I1423[i+1][j-1];
                    /*I1324: I1324(i,j;N)*/
                    I1324[i][j] = (I12[i+1][L-1] - I12[i+1][j] -I12[j][L-1] )*wValue + I1324[i][j-1] + I1324[i+1][j] - I1324[i+1][j-1];

					/*absolute value versions:*/
                    Ia1234[i][j] = I12[j+1][L-1]*abs_wValue + Ia1234[i][j-1] + Ia1234[i+1][j] - Ia1234[i+1][j-1];
                    I12a34[i][j] = Ia12[j+1][L-1]*wValue + I12a34[i][j-1] + I12a34[i+1][j] - I12a34[i+1][j-1];
                    Ia12a34[i][j] = Ia12[j+1][L-1]*abs_wValue + Ia12a34[i][j-1] + Ia12a34[i+1][j] - Ia12a34[i+1][j-1];

					Ia1324[i][j] = (I12[i+1][L-1] - I12[i+1][j] -I12[j][L-1] )*abs_wValue + Ia1324[i][j-1] + Ia1324[i+1][j] - Ia1324[i+1][j-1];
					I13a24[i][j] = (Ia12[i+1][L-1] - Ia12[i+1][j] -Ia12[j][L-1] )*wValue + I13a24[i][j-1] + I13a24[i+1][j] - I13a24[i+1][j-1];
					Ia13a24[i][j] = (Ia12[i+1][L-1] - Ia12[i+1][j] -Ia12[j][L-1] )*abs_wValue + Ia13a24[i][j-1] + Ia13a24[i+1][j] - Ia13a24[i+1][j-1];

					Ia1423[i][j] = I12[i+1][j-1]*abs_wValue +  Ia1423[i][j-1] + Ia1423[i+1][j] - Ia1423[i+1][j-1];
					I14a23[i][j] = Ia12[i+1][j-1]*wValue +  I14a23[i][j-1] + I14a23[i+1][j] - I14a23[i+1][j-1];
					Ia14a23[i][j] = Ia12[i+1][j-1]*abs_wValue +  Ia14a23[i][j-1] + Ia14a23[i+1][j] - Ia14a23[i+1][j-1];


					if (full_b == 1){
						b = i+1;
                        for (b = i+1; b <= j-2; b++){ /*#obs: the w-values needed have been computed*/
						 	wVal1 = wVal[i][b];
							abs_wVal1 = fabsf(wVal1);
							I1234_full_aid[i][j] += I12[b+1][j]*wVal1;
                            I1324_full_aid[i][j] += I12[b][j]*wVal1;

							Ia12a34_full_aid[i][j] += Ia12[b+1][j]*abs_wVal1;
							Ia1234_full_aid[i][j] += Ia12[b+1][j]*wVal1;
							I12a34_full_aid[i][j] += I12[b+1][j]*abs_wVal1;

							Ia13a24_full_aid[i][j] += I12[b][j]*abs_wVal1;
							Ia1324_full_aid[i][j] += Ia12[b][j]*wVal1;
							I13a24_full_aid[i][j] += I12[b][j]*abs_wVal1;



						}
                        /*I1234_full add up terms and recursion: I1234(i,j)*/
                        I1234_full[i][j] = I1234_full_aid[i][j] -I1234_full_aid[i][j-1];
                        I1234_full[i][j] += I1234_full[i][j-1] + I1234_full[i+1][j] - I1234_full[i+1][j-1];
                        /*I1324_full add up terms recursion: I1324(i,j)*/
                        I1324_full[i][j] = I1324_full_aid[i][j-1] - I1324_full_aid[i][j];
                        I1324_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + I1324_full[i][j-1] + I1324_full[i+1][j] - I1324_full[i+1][j-1];

						/*abs value versions:*/
						Ia12a34_full[i][j] = Ia12a34_full_aid[i][j] - Ia12a34_full_aid[i][j-1];
						Ia12a34_full[i][j] += Ia12a34_full[i][j-1] + Ia12a34_full[i+1][j] - Ia12a34_full[i+1][j-1];

						Ia1234_full[i][j] = Ia1234_full_aid[i][j] - Ia1234_full_aid[i][j-1];
						Ia1234_full[i][j] += Ia1234_full[i][j-1] + Ia1234_full[i+1][j] - Ia1234_full[i+1][j-1];

						I12a34_full[i][j] = I12a34_full_aid[i][j] -I12a34_full_aid[i][j-1];
						I12a34_full[i][j] += I12a34_full[i][j-1] + I12a34_full[i+1][j] - I12a34_full[i+1][j-1];

						Ia13a24_full[i][j] = Ia13a24_full_aid[i][j-1] - Ia13a24_full_aid[i][j];
						Ia13a24_full[i][j] += (Ia12[i][j-1] - Ia12[i+1][j-1])*(Ia12[i+1][j] - Ia12[i+1][j-1]) + Ia13a24_full[i][j-1] + Ia13a24_full[i+1][j] - Ia13a24_full[i+1][j-1];

						Ia1324_full[i][j] = Ia1324_full_aid[i][j-1] - Ia1324_full_aid[i][j];
						Ia1324_full[i][j] += (Ia12[i][j-1] - Ia12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + Ia1324_full[i][j-1] + Ia1324_full[i+1][j] - Ia1324_full[i+1][j-1];

						I13a24_full[i][j] = I13a24_full_aid[i][j-1] - I13a24_full_aid[i][j];
						I13a24_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(Ia12[i+1][j] - Ia12[i+1][j-1]) + I13a24_full[i][j-1] + I13a24_full[i+1][j] - I13a24_full[i+1][j-1];

					}
			}
		}
	}


	if (order == 3){

		for ( i=L-1; i>= 0; i-- ){
			
			/*find possible closed loops starting at first coord of segment[i]
			We discard "length three loops" by letting n run from i+4. Also, we
			want to find the "longest "irreducible loop". Note that
			since we reverse in the i-loop, so we try to extend a loop by
			adding to its "i head": if we for i have found a loop terminating
			at n, we will not search further at i (becasue of the last_j = n)
			and in the next step we will ultimately check if (i-1,n) is also
			closed (but this will only happen if we have not hit a loop 
			(i-1,n') for some n' lower than n*/
			for (n=i+minLoopLength; n <= min(i+ closedLoopLength,last_j-1); n++){
				/*if segment n is "too long" we skip the rest of the search at this i:*/
				if (distCalphas(ptr_segment[n].s1, ptr_segment[n].s2) > stdRealSegLength){
					break;
				}
				//printf("i:%d n:%d last_j:%d\n", i,n, last_j);
				if (distCalphas(ptr_segment[i].s1, ptr_segment[n].s2) < closedLoopDist){
					ptr_closedLoopInd[cnt].idx1 = i;
					ptr_closedLoopInd[cnt].idx2 = n;
					last_j = n; /*for getting only the irreducible closed loops. Obs: we reverse in the i-loop!*/
					cnt +=1;
				}
			}

			for ( j=i+1; j<=L-1; j++ ){
					wValue = w_seg12(ptr_segment[i], ptr_segment[j]);
					wVal[i][j] = wValue;
					abs_wValue = fabsf(wValue);
					I12[i][j] = wValue + I12[i][j-1] + I12[i+1][j] - I12[i+1][j-1];
					Ia12[i][j] = abs_wValue + Ia12[i][j-1] + Ia12[i+1][j] - Ia12[i+1][j-1];

					/*I1234: I1234(i,j;N)*/
                    I1234[i][j] = I12[j+1][L-1]*wValue + I1234[i][j-1] + I1234[i+1][j] - I1234[i+1][j-1];
                    /*I1423: I1423(i,j)*/
					I1423[i][j] = I12[i+1][j-1]*wValue +  I1423[i][j-1] + I1423[i+1][j] - I1423[i+1][j-1];
                    /*I1324: I1324(i,j;N)*/
                    I1324[i][j] = (I12[i+1][L-1] - I12[i+1][j] -I12[j][L-1] )*wValue + I1324[i][j-1] + I1324[i+1][j] - I1324[i+1][j-1];
					
					/*absolute value versions:*/
                    Ia1234[i][j] = I12[j+1][L-1]*abs_wValue + Ia1234[i][j-1] + Ia1234[i+1][j] - Ia1234[i+1][j-1];
                    I12a34[i][j] = Ia12[j+1][L-1]*wValue + I12a34[i][j-1] + I12a34[i+1][j] - I12a34[i+1][j-1];
                    Ia12a34[i][j] = Ia12[j+1][L-1]*abs_wValue + Ia12a34[i][j-1] + Ia12a34[i+1][j] - Ia12a34[i+1][j-1];

					Ia1324[i][j] = (I12[i+1][L-1] - I12[i+1][j] -I12[j][L-1] )*abs_wValue + Ia1324[i][j-1] + Ia1324[i+1][j] - Ia1324[i+1][j-1];
					I13a24[i][j] = (Ia12[i+1][L-1] - Ia12[i+1][j] -Ia12[j][L-1] )*wValue + I13a24[i][j-1] + I13a24[i+1][j] - I13a24[i+1][j-1];
					Ia13a24[i][j] = (Ia12[i+1][L-1] - Ia12[i+1][j] -Ia12[j][L-1] )*abs_wValue + Ia13a24[i][j-1] + Ia13a24[i+1][j] - Ia13a24[i+1][j-1];

					Ia1423[i][j] = I12[i+1][j-1]*abs_wValue +  Ia1423[i][j-1] + Ia1423[i+1][j] - Ia1423[i+1][j-1];
					I14a23[i][j] = Ia12[i+1][j-1]*wValue +  I14a23[i][j-1] + I14a23[i+1][j] - I14a23[i+1][j-1];
					Ia14a23[i][j] = Ia12[i+1][j-1]*abs_wValue +  Ia14a23[i][j-1] + Ia14a23[i+1][j] - Ia14a23[i+1][j-1];


					b = i+1;
                    for (b = i+1; b <= j-2; b++){ /*#obs: the w-values needed have been computed*/
						wVal1 = wVal[i][b];
						abs_wVal1 = fabsf(wVal1);
						I1234_full_aid[i][j] += I12[b+1][j]*wVal1;
                        I1324_full_aid[i][j] += I12[b][j]*wVal1;

						Ia12a34_full_aid[i][j] += Ia12[b+1][j]*abs_wVal1;
						Ia1234_full_aid[i][j] += Ia12[b+1][j]*wVal1;
						I12a34_full_aid[i][j] += I12[b+1][j]*abs_wVal1;

						Ia13a24_full_aid[i][j] += I12[b][j]*abs_wVal1;
						Ia1324_full_aid[i][j] += Ia12[b][j]*wVal1;
						I13a24_full_aid[i][j] += I12[b][j]*abs_wVal1;

					}
                    /*I1234_full add up terms and recursion: I1234(i,j)*/
                    I1234_full[i][j] = I1234_full_aid[i][j] -I1234_full_aid[i][j-1];
                    I1234_full[i][j] += I1234_full[i][j-1] + I1234_full[i+1][j] - I1234_full[i+1][j-1];
                    /*I1324_full add up terms recursion: I1324(i,j)*/
                    I1324_full[i][j] = I1324_full_aid[i][j-1] - I1324_full_aid[i][j];
                    I1324_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + I1324_full[i][j-1] + I1324_full[i+1][j] - I1324_full[i+1][j-1];
					/*abs value versions:*/
					Ia12a34_full[i][j] = Ia12a34_full_aid[i][j] - Ia12a34_full_aid[i][j-1];
					Ia12a34_full[i][j] += Ia12a34_full[i][j-1] + Ia12a34_full[i+1][j] - Ia12a34_full[i+1][j-1];

					Ia1234_full[i][j] = Ia1234_full_aid[i][j] - Ia1234_full_aid[i][j-1];
					Ia1234_full[i][j] += Ia1234_full[i][j-1] + Ia1234_full[i+1][j] - Ia1234_full[i+1][j-1];

					I12a34_full[i][j] = I12a34_full_aid[i][j] -I12a34_full_aid[i][j-1];
					I12a34_full[i][j] += I12a34_full[i][j-1] + I12a34_full[i+1][j] - I12a34_full[i+1][j-1];

					Ia13a24_full[i][j] = Ia13a24_full_aid[i][j-1] - Ia13a24_full_aid[i][j];
					Ia13a24_full[i][j] += (Ia12[i][j-1] - Ia12[i+1][j-1])*(Ia12[i+1][j] - Ia12[i+1][j-1]) + Ia13a24_full[i][j-1] + Ia13a24_full[i+1][j] - Ia13a24_full[i+1][j-1];

					Ia1324_full[i][j] = Ia1324_full_aid[i][j-1] - Ia1324_full_aid[i][j];
					Ia1324_full[i][j] += (Ia12[i][j-1] - Ia12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + Ia1324_full[i][j-1] + Ia1324_full[i+1][j] - Ia1324_full[i+1][j-1];

					I13a24_full[i][j] = I13a24_full_aid[i][j-1] - I13a24_full_aid[i][j];
					I13a24_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(Ia12[i+1][j] - Ia12[i+1][j-1]) + I13a24_full[i][j-1] + I13a24_full[i+1][j] - I13a24_full[i+1][j-1];

					
					/*I1423_full: I1423(i,j;N). Has to be lagged in i since we do not know I2[(i,L-1)] yet*/
                    I1423_full0[i+1][j] = (I12[i+2][j-1] - I12[i+2][j])*(I12[i+2][L-1] - I12[i+2][j] - I12[i+1][L-1] + I12[i+1][j]);
                    I1423_full0[i+1][j] += I1423_full0[i+2][j] + I1423_full0[i+1][j-1] - I1423_full0[i+2][j-1];
                    /*to compute certain degree 6 measures we use two auxillary degree 4 
                    measures; the recursion demands to sum j in the - direction (ie 
                    from above and down):*/
                    k = L-1 + i -j+1;
                    for (c =k+1;c <= L-1; c++){ /*Obs: the w values needed have already been computed (since k,c >i)*/
                        wVal2 = wVal[k][c];
                        I1324_full2_aid[i+1][k] += wVal2*(I12[i+1][c-1] - I12[i+1][k] - I12[k][c-1]);
                        I1423_full2_aid[i+1][k] += wVal2*(I12[i+1][L-1] - I12[i+1][c] - I12[k][L-1] + I12[k][c]); 
					}
					/*I1324_full2, recursion: I1324(i;j,N)*/
                    I1324_full2[i+1][k] = I1324_full2_aid[i+1][k] + I1324_full2[i+1][k+1];
                    /*I1423_full2, recursion: I1423(i;j,N)*/
					I1423_full2[i+1][k] = I1423_full2_aid[i+1][k] + I1423_full2[i+1][k+1];

					/*order 3:*/
                    /*(12)*/
                    I123456[i][j] = I1234[j+1][L-1]*wValue + I123456[i+1][j] + I123456[i][j-1] - I123456[i+1][j-1];
                    I123645[i][j] = I1423[j+1][L-1]*wValue + I123645[i+1][j] + I123645[i][j-1] - I123645[i+1][j-1];
                    I123546[i][j] = I1324[j+1][L-1]*wValue+ I123546[i+1][j] + I123546[i][j-1] - I123546[i+1][j-1];
                    /*(13)*/
                    I132456[i][j] = (I1234[i+1][L-1] - I1234[i+1][j] -I1234[j][L-1])*wValue + I132456[i+1][j] + I132456[i][j-1] - I132456[i+1][j-1];
					I132546[i+1][j] = (I1324_full2[i+2][j+1] - I1324_full[j][L-1])*wVal[i+1][j] + I132546[i+2][j] + I132546[i+1][j-1] - I132546[i+2][j-1];
					I132645[i+1][j] = (I1423_full2[i+2][j+1] - I1423[j][L-1])*wVal[i+1][j] + I132645[i+2][j] + I132645[i+1][j-1] - I132645[i+2][j-1];
                    /*(14)*/
                    /*allowing direct computation:*/
                    I142356[i][j] = I12[i+1][j-1]*I12[j+1][L-1]*wValue + I142356[i+1][j] + I142356[i][j-1] - I142356[i+1][j-1];
                    /*back to the more complicated parts, now I142536:
                    three parts recognized in decompsition as already known:*/
                    I142536[i+1][j] = wVal[i+1][j]*(I1324_full[i+2][L-1] - I1324[i+2][j] - I1324_full2[i+2][j]);
                    /*recursion part:*/
                    I142536[i+1][j] += I142536[i+2][j] + I142536[i+1][j-1] - I142536[i+2][j-1];
                    /*Next complicated part, I142635:*/
                    /*three parts recognized in decompsition as already known:*/
					I142635[i+1][j] = wVal[i+1][j]*(I1423[i+2][L-1] - I1423_full2[i+2][j] - I1423_full0[i+2][j]);
                    /*recursion part:*/
                    I142635[i+1][j] += I142635[i+2][j] + I142635[i+1][j-1] - I142635[i+2][j-1];
                    /*(15)*/
                    I152346[i][j] = wValue*(I1234[i+1][j-1] - I1234_full[i+1][j] - I12[j][L-1]*I12[i+1][j-1]);
                    I152346[i][j] += I152346[i+1][j] + I152346[i][j-1] - I152346[i+1][j-1];
                    I152436[i][j] = wValue*(I1324[i+1][j-1] - I1324_full[i+1][j]);
                    I152436[i][j] += I152436[i+1][j] + I152436[i][j-1] - I152436[i+1][j-1];
                    I152634[i][j] = wValue*(I1423_full0[i+1][j-1] - I1423[i+1][j]);
                    I152634[i][j] += I152634[i+1][j] + I152634[i][j-1] - I152634[i+1][j-1];
                    /*(16)*/
                    I162345[i][j] = wValue*I1234_full[i+1][j-1];
                    I162345[i][j] += I162345[i+1][j] + I162345[i][j-1] - I162345[i+1][j-1];
                    I162435[i][j] = wValue*I1324_full[i+1][j-1];
                    I162435[i][j] += I162435[i+1][j] + I162435[i][j-1] - I162435[i+1][j-1];
                    I162534[i][j] = wValue*I1423[i+1][j-1];
                    I162534[i][j] += I162534[i+1][j] + I162534[i][j-1] - I162534[i+1][j-1];
			}
		}
		/*the missing i==0 boundary for I1324_full2, I1423_full0 and I1423_full2:*/
		j=L-1;
		for (j=L-1; j>= 0; j--){
			i = 0;
			c = j+1;
			for (c = j+1; c <= L-1; c++){
				wValue = wVal[j][c];
				I1324_full2_aid[i][j] += wValue*(I12[i][c-1] - I12[i][j] - I12[j][c-1]);
				I1423_full2_aid[i][j] += wValue*(I12[i][L-1] - I12[i][c] - I12[j][L-1] + I12[j][c]);
			}
			/*recursion:*/
			I1324_full2[i][j] = I1324_full2_aid[i][j] + I1324_full2[i][j+1];
			I1423_full2[i][j] = I1423_full2_aid[i][j] + I1423_full2[i][j+1]; 
		}


		/*#Missing terms for some of the degree 6 measures and I1423_full0:*/
		j= 1;
		for (j=1; j<= L-1; j++ ){
			i=0;
			wValue = wVal[i][j];
			I1423_full0[i][j] = (I12[i+1][j-1] - I12[i+1][j])*(I12[i+1][L-1] - I12[i+1][j] - I12[i][L-1] + I12[i][j]);
			I1423_full0[i][j] += I1423_full0[i+1][j] + I1423_full0[i][j-1] - I1423_full0[i+1][j-1];
			/*#*/
			I132546[i][j] = (I1324_full2[i+1][j+1] - I1324_full[j][L-1])*wValue + I132546[i+1][j] + I132546[i][j-1] - I132546[i+1][j-1];
			I132645[i][j] = (I1423_full2[i+1][j+1] - I1423[j][L-1])*wValue + I132645[i+1][j] + I132645[i][j-1] - I132645[i+1][j-1];
			/*#*/
			I142536[i][j] = wValue*(I1324_full[i+1][L-1] - I1324[i+1][j] - I1324_full2[i+1][j]);
			I142536[i][j] += I142536[i+1][j] + I142536[i][j-1] - I142536[i+1][j-1]; 
			/*#*/
			I142635[i][j] = wValue*(I1423[i+1][L-1] - I1423_full2[i+1][j] - I1423_full0[i+1][j]);
			I142635[i][j] += I142635[i+1][j] + I142635[i][j-1] - I142635[i+1][j-1];

		}
	}


	*ptr_nrOfClosedLoops = cnt;

	//printf("number of closed loops found in load-PDB: %d" , *ptr_nrOfClosedLoops);


	if (order !=0 && order != 1 && order != 2 && order != 3){ 
		printf("Order must be 0,1,2 or 3");
	}

	return 1; 
}

/*Fct for examing the set of closed loops found when running aggrAndW_wClosedLoops. Two tasks are
carried out: 1) in a 2d-loop over the closed loops the mutual writhe of each pair is computed and 
written to the file defined by ptr_fileNameOut (each pair getting "type" set to "link") and 2)
in a 1d-loop over all closed loops we search upstream and downstream from the loop for "pokes",
that is short sub-chains (of length pokeLength) for which the pair (closed loop, sub-chain) has
hight writhe. These two share the "outer" for-loop over closed loops. 
Only links/pokes having a writhe which in absolute value exceeds the set threshold, thresholdLinks/thresholdPokes,
are written file.
As for the closed loops (see aggrAndW_wClosedLoops) we neglect pokes in which a segment has "too 
long" length ie which consists in an artificial peptide bond introduced by the way we load 
structures (as one long chain disregarding if it consists of several sub-chains (this can be 
avoided if changing the load procedure)*/
//OBS: This was copied from GISA_v2:
int examineClosedLoops(char *ptr_fileNameOut, struct twoSegmentIndex *ptr_closedLoopInd, int nrOfClosedLoops, struct I_ptr I_measures, int chainLength, struct segment *ptr_segment, char *ptr_fileNameOutChains, struct cAlpha *ptr_chain, int pokeLength, double thresholdLinks, double thresholdPokes){

	int returnVal = 0;

	int cntClosedLoop1 = 0;
	int cntClosedLoop2 = 0;
	int i_1 = 0;
	int j_1 = 0;
	int i_2 = 0;
	int j_2 = 0;
	int i_2_prime = 0;
	int j_1_prime = 0;
	int n = 0;
	int writeChain_b = 0; 

	int L = 0;

	double I_closedLoop;
	double realSegLength_first;
	double realSegLength;
	double realSegLength_last;
	//double stdRealSegLength = 20; /*close to square of 4.5: to avoid sub-chains with segments that are "strangely long" */
	int pokeTooLongInd_first = 0;
	int pokeTooLongInd = 0;
	int pokeTooLongInd_last = 0;
	double I_subChain;
	int cntSubChain = 0;
	double I_subChain_max = -1e+20; /*just some large negative number*/ 
	int maxAt_n = -1;
	double I_subChain_min = 1e+20; /*just some large positive number*/ 
	int minAt_n = -1;

	FILE *ptr_fileOut;

	/*results will be written to this file. Obs: we use a-mode so as NOT to clear the file -- we 
	are writing multiple times*/
	ptr_fileOut = fopen(ptr_fileNameOut, "a");


	L = chainLength -1; //number of segments in chain


	/*loop through the set of closed loops found*/
	for (cntClosedLoop1 = 0; cntClosedLoop1 <= nrOfClosedLoops-1; cntClosedLoop1 ++){

		/*start and end index of the loop, loop1:*/
		i_1 = ptr_closedLoopInd[cntClosedLoop1].idx1;
		j_1 = ptr_closedLoopInd[cntClosedLoop1].idx2;
		

		//printf("cntClosedLoop1 %d indxs %d,%d\n", cntClosedLoop1, i_1, j_1);

		/*loop through the closed loops again to find the writhe contributions
		of each pair of closed loops in the structure; we only look for pairs
		of disjoint loops in this version (cf. the if clause j_1 < i_2):*/

		for (cntClosedLoop2 = 0; cntClosedLoop2 <= nrOfClosedLoops-1; cntClosedLoop2 ++){

			/*start and end index of the loop, loop2:*/
			i_2 = ptr_closedLoopInd[cntClosedLoop2].idx1;
			j_2 = ptr_closedLoopInd[cntClosedLoop2].idx2;

			//if (i_1 < i_2){
			if (j_1 < i_2){
				
				//printf("idx 2: %d,%d ", i_2, j_2);

				/* this is actually superflous with the present if-clause surrounding it, but
				allows softening up the clause:*/
				i_2_prime = max(i_2,j_1); /*if i_2 < j_1 this must eval to j_1; else to i_2*/
				j_1_prime = min(j_1,j_2); /*if j_2 < j_1 this must eval to j_2; else to j_1*/

				/*compute writhe of the (closed loop1,closed loop2) pair*/
				I_closedLoop = I_measures.I12[i_1][j_2] - I_measures.I12[i_1][i_2-1] - I_measures.I12[j_1_prime+1][j_2] + I_measures.I12[j_1+1][i_2_prime-1]; //omit the "self terms": + I_measures.I12[i_1][j_1] + I_measures.I12[i_2][j_2];
					

				/*write results to file if writhe above threshold:*/
				if(fabs(I_closedLoop) > thresholdLinks){
					fprintf(ptr_fileOut, 
					"%s;%s;%s;%d;%d;%d;%d;%lf\n",
					I_measures.fileName,
					I_measures.chainId,
					"link",
					i_1,
					j_1,
					i_2,
					j_2,
					I_closedLoop
					);

					/*when pair of closed loops exists (as here) we want to write out the chain too:*/
					writeChain_b +=1;
				}

			}
		}

		/*search for "pokes": we look at a short sub chain of length pokeLength (disjoint 
		from closed loop 1) and compute the writhe of the pair (closed loop 1, sub chain); 
		we keep the sub chains giving the largest (+/-) writhe .
		First part of sub chain "upstream of the closed loop":*/
		/*reset*/
		cntSubChain = 0;
		pokeTooLongInd_first = 0;
		pokeTooLongInd = 0;
		pokeTooLongInd_last = 0;
		I_subChain_max = -1e+20;
		I_subChain_min = 1e+20;
		/*we want to avoid pokes in which one segment connects one chain in the structure to 
		another; we do this by a requirement on the length of each segment of a potential poke:*/
		realSegLength_first = distCalphas(ptr_segment[j_1+1].s1,ptr_segment[j_1+1].s2); /*real length (squared) of first segment in potential poke*/
		if (realSegLength_first > stdRealSegLength){
			pokeTooLongInd_first = 1;
			pokeTooLongInd = 1;
		}
		/* compute an indicator showing whether one or more of the segments within a potential 
		poke is longer than the std length of a segment (~square of about 4-5 �ngstr�m). Obs:
		segments are indexed sarting at 0 so last segment has idx L - 1*/ 
		for (n = j_1+2; n <= min(j_1 + pokeLength, L-1); n++){
			realSegLength = distCalphas(ptr_segment[n].s1,ptr_segment[n].s2);
			if (realSegLength > stdRealSegLength){
				pokeTooLongInd +=1;
			}
		}


		//printf("Chain L: %d", chainLength);
		for (n = j_1+1; n+pokeLength <= L -1; n++){
			if (pokeTooLongInd < 1){
				cntSubChain +=1; 
				/*compute writhe of the (closed loop1,subChain) pair*/
				//I_subChain = I_measures.Ia12[i_1][n+pokeLength] - I_measures.Ia12[i_1][n-1] - I_measures.Ia12[j_1+1][n+pokeLength] + I_measures.Ia12[j_1+1][n-1]; //omit the "self terms": +  I_measures.Ia12[i_1][j_1] +  I_measures.Ia12[n][n+5];
				I_subChain = I_measures.I12[i_1][n+pokeLength] - I_measures.I12[i_1][n-1] - I_measures.I12[j_1+1][n+pokeLength] + I_measures.I12[j_1+1][n-1]; //omit the "self terms": +  I_measures.Ia12[i_1][j_1] +  I_measures.Ia12[n][n+5];
				if (I_subChain > I_subChain_max){

					I_subChain_max = I_subChain;
					maxAt_n = n;
					//if (I_subChain < 0){printf("Warning! At n: %d I max:%lf\n",n, I_subChain);}

				}
				if (I_subChain < I_subChain_min){

					I_subChain_min = I_subChain;
					minAt_n = n;
					//if (I_subChain < 0){printf("Warning! At n: %d I max:%lf\n",n, I_subChain);}

				}

			}
			/*recalc the pokeTooLongInd indicator*/
			realSegLength_last = distCalphas(ptr_segment[n+pokeLength+1].s1,ptr_segment[n+pokeLength+1].s2);
			if (realSegLength_last > stdRealSegLength){
				pokeTooLongInd_last = 1;
			}
			pokeTooLongInd = pokeTooLongInd_last + pokeTooLongInd - pokeTooLongInd_first;
			/*new value of "first" to be used in next loop-step:*/
			realSegLength_first = distCalphas(ptr_segment[n+1].s1,ptr_segment[n+1].s2); /*real length (squared) of first segment in potential poke*/
			if (realSegLength_first > stdRealSegLength){
				pokeTooLongInd_first = 1;
			}


		}


		/*Second part: sub chain is downstream of the closed loop:*/
		/*reset. Obs: cntSubChain and I_subChain_max must not be reset */
		pokeTooLongInd_first = 0;
		pokeTooLongInd = 0;
		pokeTooLongInd_last = 0;
		/*first compute the pokeTooLongInd indicator:*/
		realSegLength_first = distCalphas(ptr_segment[0].s1,ptr_segment[0].s2); /*real length (squared) of first segment in potential poke*/
		if (realSegLength_first > stdRealSegLength){
			pokeTooLongInd_first = 1;
			pokeTooLongInd = 1;
		}
		/* and complete the computation:*/ 
		for (n = 1; n <= min(pokeLength, L -1); n++){
			realSegLength = distCalphas(ptr_segment[n].s1,ptr_segment[n].s2);
			if (realSegLength > stdRealSegLength){
				pokeTooLongInd +=1;
			}
		}
		/*loop thorugh the potential pokes*/
		for (n = 0; n + pokeLength <= i_1; n++){
			if (pokeTooLongInd < 1){
				cntSubChain +=1;
				/*compute writhe of the (closed loop1,subChain) pair*/
				//I_subChain = I_measures.Ia12[n][j_1] - I_measures.Ia12[n][i_1-1] - I_measures.Ia12[n+1 +pokeLength][j_1] + I_measures.Ia12[n+1 + pokeLength][i_1-1]; //omit the "self terms": + I_measures.Ia12[i_1][j_1] + I_measures.Ia12[n][n+5];
				I_subChain = I_measures.I12[n][j_1] - I_measures.I12[n][i_1-1] - I_measures.I12[n+1 +pokeLength][j_1] + I_measures.I12[n+1 + pokeLength][i_1-1]; //omit the "self terms": + I_measures.Ia12[i_1][j_1] + I_measures.Ia12[n][n+5];
				
				if (I_subChain > I_subChain_max){

					I_subChain_max = I_subChain;
					maxAt_n = n;
				}
				if (I_subChain < I_subChain_min){

					I_subChain_min = I_subChain;
					minAt_n = n;

				}
			}
			/*recalc the pokeTooLongInd indicator*/
			realSegLength_last = distCalphas(ptr_segment[n+pokeLength+1].s1,ptr_segment[n+pokeLength+1].s2);
			if (realSegLength_last > stdRealSegLength){
				pokeTooLongInd_last = 1;
			}
			pokeTooLongInd = pokeTooLongInd_last + pokeTooLongInd - pokeTooLongInd_first;
			/*new value of "first" to be used in next loop-step. Obs: n+1 is less than chainLen - 1 since inside loop:*/
			realSegLength_first = distCalphas(ptr_segment[n+1].s1,ptr_segment[n+1].s2); /*real length (squared) of first segment in potential poke*/
			if (realSegLength_first > stdRealSegLength){
				pokeTooLongInd_first = 1;
			}
		}

		/*write results to file - if any:*/
		if (cntSubChain >0 && fabs(I_subChain_max) > thresholdPokes){
			fprintf(ptr_fileOut, 
			"%s;%s;%s;%d;%d;%d;%d;%lf\n",
			I_measures.fileName,
			I_measures.chainId,
			"poke",
			i_1,
			j_1,
			maxAt_n,
			maxAt_n + pokeLength,
			I_subChain_max
			);

			/*when pair of closed loops and poke exists (as here) we want to write out the chain too:*/
			writeChain_b +=1;
		}
		if (cntSubChain >0 && fabs(I_subChain_min) > thresholdPokes){
			fprintf(ptr_fileOut, 
			"%s;%s;%s;%d;%d;%d;%d;%lf\n",
			I_measures.fileName,
			I_measures.chainId,
			"poke",
			i_1,
			j_1,
			minAt_n,
			minAt_n + pokeLength,
			I_subChain_min
			);

			/*when pair of closed loops and poke exists (as here) we want to write out the chain too:*/
			writeChain_b +=1;
		}

	}

	/*write out the chain if relevant:*/
	if (writeChain_b > 0){

		returnVal = writeChain(ptr_fileNameOutChains, I_measures.fileName, I_measures.structureName, I_measures.chainId, ptr_chain, *I_measures.chainLen);

	}

	fclose(ptr_fileOut);

	return returnVal;
}

/*Simple fct which computes the mutual writhe of all (disjoint) pairs of sub-chains; the sub-chains 
considered all have the same fixed length, subChainLength (typically between 15 and 30). Writes for 
each pair the mutual writhe value to the file defined by ptr_fileNameOut.*/
int writheSubChainPairs(char *ptr_fileNameOut, struct segment *ptr_segment, int chainLength, struct I_ptr I_measures, int subChainLength, char *ptr_fileNameOutChains, struct cAlpha *ptr_chain, int writeAll_b){

	int returnVal = 0;

	int L = 0;
	int i = 0;
	int j = 0;
	//int i_1 = 0;
	//int j_1 = 0;
	//int i_2 = 0;
	//int j_2 = 0;

    int iMin_start = 0;
	int iMin_end = subChainLength;
	int jMin_start = 0;
	int jMin_end = subChainLength;
	double Imin = 1e16; /*some gigantic number ...*/
	int iMax_start = 0;
	int iMax_end = subChainLength;
	int jMax_start = 0;
	int jMax_end = subChainLength;
	double Imax = -1e16; /*some gigantic negative number ...*/

	double I_subChainPair = 0;

	FILE *ptr_fileOut;

	/*results will be written to this file. Obs: we use a-mode so as NOT to clear the file -- we 
	are writing multiple times*/
	ptr_fileOut = fopen(ptr_fileNameOut, "a");


	L = chainLength -1; //number of segments in chain

	for(i = 0; i + subChainLength < L; i++){

		/*let us do with the case of disjoint sub-chains, so let the loop start at i + subChainLength;
		can of course easily be softened if desired:*/
		for(j = i+subChainLength; j + subChainLength < L; j++){ /*consider only j's above i to avoid double counting and j's above i+subChainLength to only get disjoint cases*/  

			
			/*compute mutual writhe of the pair sub-chain (i,i+subChainLen) and (j,j+subChainLen)*/
			I_subChainPair = I_measures.I12[i][j+subChainLength] - I_measures.I12[i][j-1] - I_measures.I12[i+subChainLength+1][j+subChainLength] + I_measures.I12[i+subChainLength+1][j-1]; //omit the "self terms" like: + I_measures.I12[i][j-1] + I_measures.I12[i+subChainLen+1][j+subChainLen+1];

			//printf("(saaledes stadig ganske hemmelig) ... %d %d %lf\n", i,j, I_subChainPair );

			if(writeAll_b == 1){
				/*write result to file:*/
				fprintf(ptr_fileOut, 
				"%s;%s;%s;%d;%d;%d;%d;%lf\n",
				I_measures.fileName,
				I_measures.chainId,
				"subChainPair",
				i,
				i+subChainLength,
				j,
				j+subChainLength,
				I_subChainPair
				);
			}
			else{
				if(I_subChainPair > Imax){
					iMax_start = i;
				    iMax_end = i+subChainLength,
					jMax_start = j;
					jMax_end = j+subChainLength;
					Imax = I_subChainPair;
				}
				if(I_subChainPair < Imin){
					iMin_start = i;
				    iMin_end = i+subChainLength,
					jMin_start = j;
					jMin_end = j+subChainLength;
					Imin = I_subChainPair;
				}
			}

		}
	
	}

	if(writeAll_b == 0){

		/*write Imin result to file:*/
		fprintf(ptr_fileOut, 
		"%s;%s;%s;%d;%d;%d;%d;%lf\n",
		I_measures.fileName,
		I_measures.chainId,
		"subChainPair",
		iMin_start,
		iMin_end,
		jMin_start,
		jMin_end,
		Imin
		);

		/*write Imax result to file:*/
		fprintf(ptr_fileOut, 
		"%s;%s;%s;%d;%d;%d;%d;%lf\n",
		I_measures.fileName,
		I_measures.chainId,
		"subChainPair",
		iMax_start,
		iMax_end,
		jMax_start,
		jMax_end,
		Imax
		);
	}

	/*write out the chain:*/
	returnVal = writeChain(ptr_fileNameOutChains, I_measures.fileName, I_measures.structureName, I_measures.chainId, ptr_chain, *I_measures.chainLen);

	fclose(ptr_fileOut);

	return returnVal;
};

/*Version of writheSubChainPairs that does the same but for a generic invariant.
Additional input: a specific measure, e.g. I_measure.I1324, as "I_measure" and its name, I_measure_name, which in the example would be the string "I1324"; 
if the specified measure is I12 the fct works as writheSubChainPairs except that it does not write out the chains.*/
int genInvariantSubChainPairs(char *ptr_fileNameOut, struct segment *ptr_segment, int chainLength, struct I_ptr I_measures, double **I_measure, char *I_measure_name, int subChainLength, int writeAll_b){

	int returnVal = 0;

	int L = 0;
	int i = 0;
	int j = 0;
	//int i_1 = 0;
	//int j_1 = 0;
	//int i_2 = 0;
	//int j_2 = 0;

    int iMin_start = 0;
	int iMin_end = subChainLength;
	int jMin_start = 0;
	int jMin_end = subChainLength;
	double Imin = 1e16; /*some gigantic number ...*/
	int iMax_start = 0;
	int iMax_end = subChainLength;
	int jMax_start = 0;
	int jMax_end = subChainLength;
	double Imax = -1e16; /*some gigantic negative number ...*/

	double I_subChainPair = 0;

	FILE *ptr_fileOut;

	/*results will be written to this file. Obs: we use a-mode so as NOT to clear the file -- we 
	are writing multiple times*/
	ptr_fileOut = fopen(ptr_fileNameOut, "a");


	L = chainLength -1; //number of segments in chain

	for(i = 0; i + subChainLength < L; i++){

		/*let us do with the case of disjoint sub-chains, so let the loop start at i + subChainLength;
		can of course easily be softened if desired:*/
		for(j = i+subChainLength; j + subChainLength < L; j++){ /*consider only j's above i to avoid double counting and j's above i+subChainLength to only get disjoint cases*/  

			
			/*compute mutual writhe of the pair sub-chain (i,i+subChainLen) and (j,j+subChainLen)*/
			I_subChainPair = I_measure[i][j+subChainLength] - I_measure[i][j-1] - I_measure[i+subChainLength+1][j+subChainLength] + I_measure[i+subChainLength+1][j-1]; //omit the "self terms" like: + I_measures.I12[i][j-1] + I_measures.I12[i+subChainLen+1][j+subChainLen+1];

			if(writeAll_b == 1){
				/*write result to file:*/
				fprintf(ptr_fileOut, 
				"%s;%s;%s;%d;%d;%d;%d;%lf\n",
				I_measures.fileName,
				I_measures.chainId,
				I_measure_name, /*we here plug in the invariant's name rather than a desciption as in writheSubChainPairs; this allows using the same load-fct in subsequent plotting as for data generated by writheSubChainPairs */
				i,
				i+subChainLength,
				j,
				j+subChainLength,
				I_subChainPair
				);
			}
			else{
				if(I_subChainPair > Imax){
					iMax_start = i;
				    iMax_end = i+subChainLength,
					jMax_start = j;
					jMax_end = j+subChainLength;
					Imax = I_subChainPair;
				}
				if(I_subChainPair < Imin){
					iMin_start = i;
				    iMin_end = i+subChainLength,
					jMin_start = j;
					jMin_end = j+subChainLength;
					Imin = I_subChainPair;
				}
			}

		}
	
	}

	if(writeAll_b == 0){

		/*write Imin result to file:*/
		fprintf(ptr_fileOut, 
		"%s;%s;%s;%d;%d;%d;%d;%lf\n",
		I_measures.fileName,
		I_measures.chainId,
		I_measure_name, /*we here plug in the invariant's name rather than a desciption as in writheSubChainPairs; this allows using the same load-fct in subsequent plotting as for data generated by writheSubChainPairs */
		iMin_start,
		iMin_end,
		jMin_start,
		jMin_end,
		Imin
		);

		/*write Imax result to file:*/
		fprintf(ptr_fileOut, 
		"%s;%s;%s;%d;%d;%d;%d;%lf\n",
		I_measures.fileName,
		I_measures.chainId,
		I_measure_name, /*we here plug in the invariant's name rather than a desciption as in writheSubChainPairs; this allows using the same load-fct in subsequent plotting as for data generated by writheSubChainPairs */
		iMax_start,
		iMax_end,
		jMax_start,
		jMax_end,
		Imax
		);
	}


	fclose(ptr_fileOut);

	returnVal = 0;

	return returnVal;
};


/*Functions aimed at alignment use*/

/*Aux function for the fcts getting invariant values on single windows or pairs of windows. The
purpose is to allow getting a good set of windows, given a set window length and a given chain length.
One may either modify the set chain length, by increasing or decreasing it, so that the given chain length
is well covered; or one may use the set window length and place the largest number of (half) windows
in the chain, with an offset from the start of the chain to the beginning of the first window (default).
The input "type" allows switching between different versions:
type = 1 : this will return an offset and the largest number of overlapping windows that can be placed in the 
given chain (overlaps by half the set window length). The offset is chosen so that the gap at the end of the 
chain will have a value smaller then the offset but as close to it as possible
type = 2: find a good halfwindowlength: the shortest length such that ceil of chainLength/pre-set 
halfwindowlength number of windows of this length will cover the chainLength (ie be longer
than the chainLength)*/
int setWindowCharacteristics(int type, int *ptr_windowLength,  int chainLength,  int *ptr_stepSize, int *ptr_offset, int *ptr_numberOfWindows){

	int returnVal = 0;

	int l; /*halfWindowLength;*/
	int windowLength = *ptr_windowLength;
	int windowLengthNew = *ptr_windowLength;
	int stepSize = *ptr_stepSize;

	int maxNrOfSteps = 0;
	int maxNrOfHalfWindows = 0;
	int gap = 0;
	int offset = 0;
	int numberOfWindows = *ptr_numberOfWindows;

	double x = 0, x12 = 0, x0 = 0, x1 = 0;

	int print_final_b = 0;

	/*OBS: thoughout this fct we replace chainLen by chainLen - 1: it's not the chain of residues we are considering but the
	chain of segments*/


	if(type == 0){

		/*Max number of half-windows in chain:*/
		numberOfWindows = (int) floor((double) (chainLength - 1 -  windowLength)/stepSize) + 1; //+1: the subtracted windowLength corr to the first window; the nr of steps fitting in the remaining chainLength - windowLength is the fist summand 

		/*Find a good offset:*/
		gap = chainLength - 1 - (windowLength + (numberOfWindows-1)*stepSize);
		offset = (int) floor((double) gap*0.5);

		//printf("number of wins %d gap: %d\n", numberOfWindows, gap);
		
	}

	if(type == 1){

		l = floor(windowLength*0.5);

		stepSize = l;

		/*Max number of half-windows in chain:*/
		maxNrOfHalfWindows = floor((double) chainLength - 1/l);

		/*Find a good offset:*/
		gap = chainLength - 1 - maxNrOfHalfWindows*l;
		offset = floor((double) gap*0.5);
		
		//No matter whether maxNrOfHalfWindows is odd or even:
		numberOfWindows = maxNrOfHalfWindows - 1;

	}

	if(type == 2){

		l = floor(windowLength*0.5);

		if(print_final_b ==1){

			printf("Initial window length:%d\n", windowLength);
			printf("Initial half-window length:%d\n", l);
		}

		x = ((double)chainLength - 1/windowLength); /*intial nr of windows*/
		x0 = floor(x);
		x1 = ceil(x);
		x12 = (x0 + x1)*0.5;

		if(print_final_b ==1){

			printf("x: %f floor:%f ceil:%f x12:%f\n" , x, x0,x1,x12);

		}

		windowLengthNew = windowLength;

		if (chainLength - 1 >= windowLengthNew){ /*else the floor, x0, is 0 and it becomes possible to enter an inf loop; and: there's no need for recalc of l in this case*/

			if(chainLength - 1 > 2*l*x12){
				while(chainLength - 1 > 2*l*x0 ){ 
					l = l+1;
				}

			}

			if(chainLength - 1 <= 2*l*x12){
				while(chainLength - 1 <= 2*l*x1){
					l = l-1;
				}
				l = l+1;
			}

			/*(re)set the windowslength*/
			windowLengthNew = 2*l;

			stepSize = l;
			
			/*Derive the corresponding number of windows:*/
			/*Max number of half-windows in chain:*/
			maxNrOfHalfWindows = floor((double) chainLength - 1/l);
			numberOfWindows = maxNrOfHalfWindows - 1;

		}
	
	}

	*ptr_windowLength = windowLengthNew;
	*ptr_offset = offset;
	*ptr_stepSize = stepSize;
	*ptr_numberOfWindows = numberOfWindows;
	 
	if(print_final_b ==1){

		printf("setWindowCharacteristics output: windowLengthNew %d, offset: %d, numberOfWindows: %d \n", windowLengthNew, offset, numberOfWindows);

	}

	return returnVal;

}

int getInvariantsOnWindows(struct I_windows_ptr *ptr_I_windows, int order, int type, int windowLength,  int stepSize, int chainLength, struct I_ptr I_measures){
					
	/*Missing: normalize the invariants? Maybe not: rather normalize to std dev 1 over a large set of proteins*/

	/*OBS: if order = 2, this function can only be run with full_b = 1 !!*/ 

	int returnVal = 0;

	struct I_windows_ptr I_windows;

	int windowLengthNew;
	int offset = 0;
	int stepSizeNew = 0;
	int numberOfWindows = 0;
	
	int halfWindowLength = 0;
	int m,i1, i2;
	
	int print_final_b = 0;

	double stdDev[37] = {1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

	//{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	//window lgth 10: {2.8794183730387206, 3.1373627520852181, 251.87019317644356, 2.7544917488659966, 0.47301920743168391, 3.8497904145625688, 2.1784048338891364, 400.20921877115393, 2037.6031948047, 3257.3206209926161, 4.7203358510219697, 16.518798056770912, 29.216804503605964, 1.801520668458215, 1.0331604384066702, 2.0787892397180379, 21952.093983224724, 206.74864790956718, 929.94712645316554, 231.7257303569755, 3.1180518751356461, 42.852396883777075, 33.438321327026095, 0.5596441088540719, 0.62681652219491957, 0.4348036340119264, 0.27232630330540153, 0.29442378583176942, 0.099114767962417319, 0.33554987618933851, 0.088768225856408936};
	//window lgth 16: {4.4641099154513864, 4.920892915993436, 398.74398063425411, 4.0894212008222288, 1.5076672087020262, 12.874152081421339, 3.5089490627166833, 737.45779331567678, 3228.6523449016881, 6018.7337513209841, 8.7033642421426336, 32.786392205572056, 64.696119779696602, 5.9926969652357576, 3.8498805893342247, 7.3886713042331769, 35367.198751878968, 330.76504427734363, 1475.7093035028699, 354.66080480341117, 5.0016344992252355, 69.738986611449434, 106.03625057380583, 1.343638488723234, 1.6022416124506742, 1.6714059972318025, 0.79802966567923139, 1.4092649366256791, 0.85681309514609094, 0.90785801037173508, 0.55810075421482297};
	//window lgth 20: {5.3901399770479355, 5.861571637452208, 488.48519950899049, 4.8291856160488846, 2.6354948407460408, 20.652139701104851, 4.1766169732877003, 977.20939163810806, 3951.0438587733843, 7992.8542996233682, 11.665622320959576, 46.347680950488332, 97.59755169311471, 10.278999427155021, 7.028785350832754, 13.080755471847077, 43765.756516314708, 407.56627475774997, 1811.0060595203281, 427.83802805200924, 6.2285917892647085, 86.909696052493516, 185.65893178650711, 2.122085995464611, 2.6302599979070207, 3.2703387581019276, 1.3219550074824022, 2.963513771457376, 2.4721148444560801, 1.5228104995423042, 1.2113966085273393};


	I_windows = *ptr_I_windows;

	if(print_final_b ==1){

		printf("Initial window length:%d\n", windowLength);
	}


	windowLengthNew = windowLength;
	stepSizeNew = stepSize;

	/*get an appropriate window length, using the type provided:*/

	returnVal = setWindowCharacteristics(type, &windowLengthNew,  chainLength, &stepSizeNew, &offset, &numberOfWindows);

	halfWindowLength = floor((double) windowLengthNew/2);

	//printf("numberOfWindows in getInvariantsOnWindows: %d\n", numberOfWindows);

	I_windows.fileName = I_measures.fileName;
	I_windows.structureName = I_measures.structureName;
	I_windows.classId = I_measures.classId;
	I_windows.chainId = I_measures.chainId;
	I_windows.chainNr = *I_measures.chainNr;
	I_windows.chainLen = *I_measures.chainLen;
	I_windows.order = *I_measures.order;
	I_windows.windowLgth = windowLengthNew;
	I_windows.stepSize = stepSizeNew;
	
	if(print_final_b ==1){

		printf("chain len: %d; window lgth: %d\n", chainLength, I_windows.windowLgth);

	}

	/*fetch the invariant values*/
	for (m = 0; m < numberOfWindows ; m++){

		I_windows.window[m].windowNr = m;

		//printf("window nr:%d\n",I_windows.windowsNr[m]);

		i1= (int)(m*stepSizeNew + offset);
		i2 = (int)(min(chainLength-1, i1 + windowLength));
		if(i1 + windowLength > chainLength-1){
			printf("Warning: window length %d and nr of windows %d inapproriate for chain length %d!",windowLength, numberOfWindows, chainLength );
		}
		if(print_final_b ==1){
			printf("i1: %d, i2:%d\n", i1, i2);
		}

		I_windows.window[m].segIndices[0] = i1;
		I_windows.window[m].segIndices[1] = i2;
		
		if (order >= 1){
			I_windows.I12[m] = I_measures.I12[i1][i2]/stdDev[0];
			I_windows.Ia12[m] = I_measures.Ia12[i1][i2]/stdDev[1];
		}
		/*order 2*/
		if (order >= 2){

			I_windows.I1234[m] = I_measures.I1234[i1][i2]/stdDev[2];
			I_windows.I1324[m] = I_measures.I1324[i1][i2]/stdDev[3];
			I_windows.I1423[m] = I_measures.I1423[i1][i2]/stdDev[4];

			I_windows.I1234_full[m] = I_measures.I1234_full[i1][i2]/stdDev[5];
			I_windows.I1324_full[m] = I_measures.I1324_full[i1][i2]/stdDev[6];

			I_windows.Ia12a34_full[m] = I_measures.Ia12a34_full[i1][i2]/stdDev[7];
			I_windows.Ia13a24_full[m] = I_measures.Ia13a24_full[i1][i2]/stdDev[8];

			I_windows.Ia1234_full[m] = I_measures.Ia1234_full[i1][i2]/stdDev[9];
			I_windows.Ia1324_full[m] = I_measures.Ia1324_full[i1][i2]/stdDev[10];

			I_windows.I12a34_full[m] = I_measures.I12a34_full[i1][i2]/stdDev[11];
			I_windows.I13a24_full[m] = I_measures.I13a24_full[i1][i2]/stdDev[12];

			//*abs values versions*/
			//*12*/
			I_windows.Ia1234[m] = I_measures.Ia1234[i1][i2]/stdDev[13];
			I_windows.I12a34[m] = I_measures.I12a34[i1][i2]/stdDev[14];
			I_windows.Ia12a34[m] = I_measures.Ia12a34[i1][i2]/stdDev[15];
			//*13*/
			I_windows.Ia1324[m] = I_measures.Ia1324[i1][i2]/stdDev[16];
			I_windows.I13a24[m] = I_measures.I13a24[i1][i2]/stdDev[17];
			I_windows.Ia13a24[m] = I_measures.Ia13a24[i1][i2]/stdDev[18];
			//*14*/
			I_windows.Ia1423[m] = I_measures.Ia1423[i1][i2]/stdDev[19];
			I_windows.I14a23[m] = I_measures.I14a23[i1][i2]/stdDev[20];
			I_windows.Ia14a23[m] = I_measures.Ia14a23[i1][i2]/stdDev[21];


		}
		//*order 3*/
		if (order >= 3){
			I_windows.I123456[m] = I_measures.I123456[i1][i2]/stdDev[22];
			I_windows.I123546[m] = I_measures.I123546[i1][i2]/stdDev[23];
			I_windows.I123645[m] = I_measures.I123645[i1][i2]/stdDev[24];

			I_windows.I132456[m] = I_measures.I132456[i1][i2]/stdDev[25];
			I_windows.I132546[m] = I_measures.I132546[i1][i2]/stdDev[26];
			I_windows.I132645[m] = I_measures.I132645[i1][i2]/stdDev[27];

			I_windows.I142356[m] = I_measures.I142356[i1][i2]/stdDev[28];
			I_windows.I142536[m] = I_measures.I142536[i1][i2]/stdDev[29];
			I_windows.I142635[m] = I_measures.I142635[i1][i2]/stdDev[30];

			I_windows.I152346[m] = I_measures.I152346[i1][i2]/stdDev[31];
			I_windows.I152436[m] = I_measures.I152436[i1][i2]/stdDev[32];
			I_windows.I152634[m] = I_measures.I152634[i1][i2]/stdDev[33];

			I_windows.I162345[m] = I_measures.I162345[i1][i2]/stdDev[34];
			I_windows.I162435[m] = I_measures.I162435[i1][i2]/stdDev[35];
			I_windows.I162534[m] = I_measures.I162534[i1][i2]/stdDev[36];
		}

	}

	if(m!=numberOfWindows){
		printf("Warning: in getInvariantsOnWindows the populated nr of windows was %d while the derived numberOfWindows was %d\n", m, numberOfWindows);
	}

	I_windows.nrOfWindows = m;


	*ptr_I_windows = I_windows;

	if(print_final_b ==1){
		printf("In getInvariantsOnWindows: final numberOfWindows %d\n", ptr_I_windows -> nrOfWindows);
	}

	return 1;

}

int writeInvariantsOnWindows(char *ptr_fileNameOut, struct I_windows_ptr I_windows, int chainLength, int order){

	FILE *ptr_fileOut;

	int m = 0, i1 = 0, i2 = 0;

	int nrOfWindows;

	int print_final_b = 0;

	///*MISSING: NORMALIZE THE INVS?? Maybe not: rather normalize to std dev 1 over a large set of proteins*/

	//int l; /*halfWindowLength;*/
	//int windowLengthNew;
	////int nrWindows = 0; 

	//double x = 0, x12 = 0, x0 = 0, x1 = 0;

	//int m = 0;  
	//int i1 = 0, i2 = 0; 

	//int print_final_b = 0;

	//FILE *ptr_fileOut;


	///*find a good halfwindowlength: the shortest length such that ceil of chainLength/pre-set 
	//halfwindowlength number of windows of this length will cover the chainLength (ie be longer
	//than the chainLength)*/

	//l = floor(windowLength*0.5);

	//if(print_final_b ==1){

	//	printf("Initial window length:%d\n", windowLength);
	//	printf("Initial half-window length:%d\n", l);
	//}

	//x = ((double)chainLength/windowLength); /*intial nr of windows*/
	//x0 = floor(x);
	//x1 = ceil(x);
	//x12 = (x0 + x1)*0.5;

	//if(print_final_b ==1){

	//	printf("x: %f floor:%f ceil:%f x12:%f\n" , x, x0,x1,x12);

	//}

	//windowLengthNew = windowLength;

	//if (chainLength >= windowLengthNew){ /*else the floor, x0, is 0 and it becomes possible to enter an inf loop; and: there's no need for recalc of l in this case*/

	//	if(chainLength > 2*l*x12){
	//		while(chainLength > 2*l*x0 ){ 
	//			l = l+1;
	//		}

	//	}

	//	if(chainLength <= 2*l*x12){
	//		while(chainLength <= 2*l*x1){
	//			l = l-1;
	//		}
	//		l = l+1;
	//	}

	//	/*(re)set the windowslength*/
	//	windowLengthNew = 2*l;
	//	/*x = ((double)chainLength/windowLength);
	//	nrWindows = 2*ceil(x) - 1;*/

	//}


	//nrOfWindows = (int) floor((double) (chainLength-1)/(0.5*I_windows.windowLgth));



	//printf("nr of Ws:%d\n", nrOfWindows);

	/*results will be written to this file. Obs: we use a-mode since we 
	are writing multiple times*/
	ptr_fileOut = fopen(ptr_fileNameOut, "a");

	/* write header line first:*/
	fprintf(ptr_fileOut, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;\n", 
			"Line type",
			"PDBfile",
			"structureId",
			"classId",
			"chainId",
			"chainNr",
		    "chainLength",
			"order",
			"nrOfWindows",
			"windowLength",
			"stepSize");
	/* write the header content:*/
	fprintf(ptr_fileOut, "%s;%s;%s;%s;%s;%d;%d;%d;%d;%d;%d;\n", 
			"Header",
			I_windows.fileName,
			I_windows.structureName,
			I_windows.classId,
			I_windows.chainId,
			I_windows.chainNr,
			I_windows.chainLen,
			I_windows.order,
			I_windows.nrOfWindows,
			I_windows.windowLgth,
			I_windows.stepSize
			);


	if (order == 1){
		/* write header line first:*/
		fprintf(ptr_fileOut, "%s;%s;%s;%s;%s;%s;\n", 
				"Values",
				"windowNr",
				"i1",
				"i2",
				/*order 1*/
				"I12", 
				"Ia12");

		if(print_final_b ==1){

			printf("chain len: %d; window lgth: %d\n", chainLength,I_windows.windowLgth );

		}

		/*write out the invariant values*/
		for (m = 0; m < I_windows.nrOfWindows ; m++){

			i1 = I_windows.window[m].segIndices[0];
			i2 = I_windows.window[m].segIndices[1];
			
			if(print_final_b ==1){
				printf("i1: %d, i2:%d\n", i1, i2);
			}

			fprintf(ptr_fileOut, 
			"%d;%d;%d;%lf;%lf\n",
			I_windows.window[m].windowNr,
			i1,
			i2,
			I_windows.I12[m],
			I_windows.Ia12[m]
			)
			;

		}
	}

	
	if (order == 2){

		/* write header line first:*/
		fprintf(ptr_fileOut, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;\n", 
				"Values",
				"windowNr",
				"i1",
				"i2",
				/*order 1*/
				"I12", 
				"Ia12",
				/*order 2 full:*/
				"I1234_full", 
				"I1324_full", 
				/*abs value versions order 2 full:*/
				"Ia12a34_full", 
				"Ia1234_full", 
				"I12a34_full", 
				"Ia13a24_full", 
				"Ia1324_full", 
				"I13a24_full", 
				/*order 2:*/
				"I1234", 
				"I1324", 
				"I1423",
				/*abs value versions:*/
				/*12*/
				"Ia1234",
				"I12a34",
				"Ia12a34",
				/*13*/
				"Ia1324",
				"I13a24",
				"Ia13a24",
				/*14*/
				"Ia1423",
				"I14a23",
				"Ia14a23");

		if(print_final_b ==1){

			printf("chain len: %d; window lgth: %d\n", chainLength,I_windows.windowLgth );

		}

		/*write out the invariant values*/
		for (m = 0; m < I_windows.nrOfWindows ; m++){

			i1 = I_windows.window[m].segIndices[0];
			i2 = I_windows.window[m].segIndices[1];
			if(print_final_b ==1){
				printf("i1: %d, i2:%d\n", i1, i2);
			}

			fprintf(ptr_fileOut, 
			"%d;%d;%d;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;\n",
			I_windows.window[m].windowNr,
			i1,
			i2,
			I_windows.I12[m],
			I_windows.Ia12[m],
			/*order 2, full*/
			I_windows.I1234_full[m],
			I_windows.I1324_full[m],
			/*abs values versions order 2 full*/
			I_windows.Ia12a34_full[m],
			I_windows.Ia1234_full[m],
			I_windows.I12a34_full[m],
			I_windows.Ia13a24_full[m],
			I_windows.Ia1324_full[m],
			I_windows.I13a24_full[m],
			/*order 2*/
			I_windows.I1234[m],
			I_windows.I1324[m],
			I_windows.I1423[m],
			/*abs values versions*/
			/*12*/
			I_windows.Ia1234[m],
			I_windows.I12a34[m],
			I_windows.Ia12a34[m],
			/*13*/
			I_windows.Ia1324[m],
			I_windows.I13a24[m],
			I_windows.Ia13a24[m],
			/*14*/
			I_windows.Ia1423[m],
			I_windows.I14a23[m],
			I_windows.Ia14a23[m]
			)
			;
		}
	}


	if (order == 3){

		/* write header line first:*/
		fprintf(ptr_fileOut,"%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;\n", 
				"Values",
				"windowNr",
				"i1",
				"i2",
				/*order 1*/
				"I12", 
				"Ia12",
				/*order 2 full:*/
				"I1234_full", 
				"I1324_full", 
				/*abs value versions order 2 full:*/
				"Ia12a34_full", 
				"Ia1234_full", 
				"I12a34_full", 
				"Ia13a24_full", 
				"Ia1324_full", 
				"I13a24_full", 
				/*order 2:*/
				"I1234", 
				"I1324", 
				"I1423",
				/*abs value versions:*/
				/*12*/
				"Ia1234",
				"I12a34",
				"Ia12a34",
				/*13*/
				"Ia1324",
				"I13a24",
				"Ia13a24",
				/*14*/
				"Ia1423",
				"I14a23",
				"Ia14a23",
				/*order 3*/
				"I123456",
				"I123546",
				"I123645",
				/*13*/
				"I132456",
				"I132546",
				"I132645",
				/*14*/
				"I142356",
				"I142536",
				"I142635",
				/*15*/
				"I152346",
				"I152436",
				"I152634",
				/*16*/
				"I162345",
				"I162435",
				"I162534");
	

		if(print_final_b ==1){

			printf("chain len: %d; window lgth: %d\n", chainLength,I_windows.windowLgth );

		}

		/*write out the invariant values*/
		for (m = 0; m < I_windows.nrOfWindows ; m++){

			i1 = I_windows.window[m].segIndices[0];
			i2 = I_windows.window[m].segIndices[1];
			if(print_final_b ==1){
				printf("i1: %d, i2:%d\n", i1, i2);
			}

			fprintf(ptr_fileOut, 
			"%d;%d;%d;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;\n",
			I_windows.window[m].windowNr,
			i1,
			i2,
			I_windows.I12[m],
			I_windows.Ia12[m],
			/*order 2, full*/
			I_windows.I1234_full[m],
			I_windows.I1324_full[m],
			/*abs values versions order 2 full*/
			I_windows.Ia12a34_full[m],
			I_windows.Ia1234_full[m],
			I_windows.I12a34_full[m],
			I_windows.Ia13a24_full[m],
			I_windows.Ia1324_full[m],
			I_windows.I13a24_full[m],
			/*order 2*/
			I_windows.I1234[m],
			I_windows.I1324[m],
			I_windows.I1423[m],
			/*abs values versions*/
			/*12*/
			I_windows.Ia1234[m],
			I_windows.I12a34[m],
			I_windows.Ia12a34[m],
			/*13*/
			I_windows.Ia1324[m],
			I_windows.I13a24[m],
			I_windows.Ia13a24[m],
			/*14*/
			I_windows.Ia1423[m],
			I_windows.I14a23[m],
			I_windows.Ia14a23[m],
			/*order 3*/
			I_windows.I123456[m],
			I_windows.I123546[m],
			I_windows.I123645[m],

			I_windows.I132456[m],
			I_windows.I132546[m],
			I_windows.I132645[m],

			I_windows.I142356[m],
			I_windows.I142536[m],
			I_windows.I142635[m],

			I_windows.I152346[m],
			I_windows.I152436[m],
			I_windows.I152634[m],

			I_windows.I162345[m],
			I_windows.I162435[m],
			I_windows.I162534[m]
			)
			;

		}
	}

	fclose(ptr_fileOut);

	return 1;

}

//OBS: THIS ONLY WORKS/MAKES SENSE UP TO AND INCL ORDER 2:
int getInvariantsOnWindowPairs(struct I_windowPairs_ptr *ptr_I_windowPairs, int order, int type, int windowLength,  int stepSize,  int chainLength, struct I_ptr I_measures){
					
	/*Missing: normalize the invariants? Maybe not: rather normalize to std dev 1 over a large set of proteins*/
	
	int returnVal = 0;

	struct I_windowPairs_ptr I_windowPairs = *ptr_I_windowPairs;

	int windowLengthNew;
	int stepSizeNew = 0;
	int offset = 0;
	int numberOfWindows = 0;
	
	int halfWindowLength = 0;
	int m, i1, i2, n, j1, j2;

	int cntWinPairs = 0;

	int print_final_b = 0;


	double stdDev[22] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	//window lgth 10: {2.8794183730387206, 3.1373627520852181, 251.87019317644356, 2.7544917488659966, 0.47301920743168391, 3.8497904145625688, 2.1784048338891364, 400.20921877115393, 2037.6031948047, 3257.3206209926161, 4.7203358510219697, 16.518798056770912, 29.216804503605964, 1.801520668458215, 1.0331604384066702, 2.0787892397180379, 21952.093983224724, 206.74864790956718, 929.94712645316554, 231.7257303569755, 3.1180518751356461, 42.852396883777075, 33.438321327026095, 0.5596441088540719, 0.62681652219491957, 0.4348036340119264, 0.27232630330540153, 0.29442378583176942, 0.099114767962417319, 0.33554987618933851, 0.088768225856408936};
	//window lgth 16: {2.1113146875549202, 4.5080284706900366, 54.938331791844057, 2.8371922664641245,  37.451849322379303}
		//window lgth 20: {5.3901399770479355, 5.861571637452208, 488.48519950899049, 4.8291856160488846, 2.6354948407460408, 20.652139701104851, 4.1766169732877003, 977.20939163810806, 3951.0438587733843, 7992.8542996233682, 11.665622320959576, 46.347680950488332, 97.59755169311471, 10.278999427155021, 7.028785350832754, 13.080755471847077, 43765.756516314708, 407.56627475774997, 1811.0060595203281, 427.83802805200924, 6.2285917892647085, 86.909696052493516, 185.65893178650711, 2.122085995464611, 2.6302599979070207, 3.2703387581019276, 1.3219550074824022, 2.963513771457376, 2.4721148444560801, 1.5228104995423042, 1.2113966085273393};

	
	windowLengthNew = windowLength;
	stepSizeNew = stepSize;

	/*get an appropriate window length, using the type provided:*/
	returnVal = setWindowCharacteristics(type, &windowLengthNew,  chainLength, &stepSizeNew, &offset, &numberOfWindows);
	//printf("in getInvariantsOnWindowPairs; numberOfWindows %d\n", numberOfWindows);
	halfWindowLength = floor((double) windowLengthNew/2);


	if(print_final_b ==1){

		printf("Window length:%d\n", windowLength);
	}

	strcpy(I_windowPairs.fileName, I_measures.fileName);
	strcpy(I_windowPairs.structureName, I_measures.structureName);
	strcpy(I_windowPairs.classId , I_measures.classId);
	strcpy(I_windowPairs.chainId , I_measures.chainId);

	/* I_windowPairs.fileName = I_measures.fileName;
	I_windowPairs.structureName = I_measures.structureName;
	I_windowPairs.classId = I_measures.classId;
	I_windowPairs.chainId = I_measures.chainId;*/

	//printf("Sommeren mild oedlses ...\n");

	I_windowPairs.chainNr = *I_measures.chainNr;
	I_windowPairs.chainLen = *I_measures.chainLen;
	I_windowPairs.order = *I_measures.order;
	I_windowPairs.windowLgth = windowLength;
	
	if(print_final_b ==1){

		printf("chain len: %d; window lgth: %d\n", chainLength, I_windowPairs.windowLgth);

	}

	cntWinPairs = 0;
	/*fetch the invariant values*/
	for (m = 0; m < numberOfWindows; m++){

		//I_windowPairs.windowPair.windowsNr_1[m] = m;

		//printf("window nr:%d\n",*I_windowPairs.windowsNr);
		i1= (int)(m*stepSizeNew + offset);
		i2 = (int)(min(chainLength-1, i1 + windowLength));
		if(i1 + windowLength > chainLength-1){
			printf("Warning: window length %d and nr of windows %d inapproriate for chain length %d!",windowLength, numberOfWindows, chainLength );
		}
		if(print_final_b ==1){
			printf("i1: %d, i2:%d\n", i1, i2);
		}

		//if (m==0) printf("i1: %d, i2:%d\n", i1, i2);

		//I_windowPairs.windowPair.segIndices_1[m][0] = i1;
		//I_windowPairs.windowPair.segIndices_1[m][1] = i2;
				
		for (n = m + 1 ; n < numberOfWindows; n++){

			I_windowPairs.windowPair[m][n].windowNr_1 = m;
			I_windowPairs.windowPair[m][n].segIndices_1[0] = i1;
			I_windowPairs.windowPair[m][n].segIndices_1[1] = i2;

			//I_windowPairs.windowPair.windowsNr_2[n] = n;

			j1= (int)(n*stepSizeNew + offset);
			j2 = (int)(min(chainLength-1, j1 + windowLength));
			if(print_final_b ==1){
				printf("j1: %d, j2:%d\n", j1, j2);
			}

			//if (n == 0) printf("j1: %d, j2:%d\n", j1, j2);
				
			//I_windowPairs.segIndices_2[n][0] = j1;
			//I_windowPairs.segIndices_2[n][1] = j2;

			I_windowPairs.windowPair[m][n].windowNr_2 = n;
			I_windowPairs.windowPair[m][n].segIndices_2[0] = j1;
			I_windowPairs.windowPair[m][n].segIndices_2[1] = j2;

			if (order >= 1){
				I_windowPairs.I12[m][n] = I_mutual(I_measures.I12, i1, i2, j1, j2, stdDev[0]);
				I_windowPairs.Ia12[m][n] = I_mutual(I_measures.Ia12, i1, i2, j1, j2, stdDev[1]); 
			}

			if (order >= 2){

				/*order 2*/
				I_windowPairs.I1234[m][n] = I_mutual(I_measures.I1234, i1, i2, j1, j2, stdDev[2]);
				I_windowPairs.I1324[m][n] = I_mutual(I_measures.I1324, i1, i2, j1, j2, stdDev[3]);
				I_windowPairs.I1423[m][n] = I_mutual(I_measures.I1423, i1, i2, j1, j2, stdDev[4]);

				I_windowPairs.I1234_full[m][n] = I_mutual(I_measures.I1234_full, i1, i2, j1, j2, stdDev[5]);
				I_windowPairs.I1324_full[m][n] = I_mutual(I_measures.I1324_full, i1, i2, j1, j2, stdDev[6]);

				I_windowPairs.Ia12a34_full[m][n] = I_mutual(I_measures.Ia12a34_full, i1, i2, j1, j2, stdDev[7]);
				I_windowPairs.Ia13a24_full[m][n] = I_mutual(I_measures.Ia13a24_full, i1, i2, j1, j2, stdDev[8]);

				I_windowPairs.Ia1234_full[m][n] = I_mutual(I_measures.Ia1234_full, i1, i2, j1, j2, stdDev[9]);
				I_windowPairs.Ia1324_full[m][n] = I_mutual(I_measures.Ia1324_full, i1, i2, j1, j2, stdDev[10]);

				I_windowPairs.I12a34_full[m][n] = I_mutual(I_measures.I12a34_full, i1, i2, j1, j2, stdDev[11]);
				I_windowPairs.I13a24_full[m][n] = I_mutual(I_measures.I13a24_full, i1, i2, j1, j2, stdDev[12]);

				//*abs values versions*/
				//*12*/
				I_windowPairs.Ia1234[m][n] = I_mutual(I_measures.Ia1234, i1, i2, j1, j2, stdDev[13]);
				I_windowPairs.I12a34[m][n] = I_mutual(I_measures.I12a34, i1, i2, j1, j2, stdDev[14]);
				I_windowPairs.Ia12a34[m][n] = I_mutual(I_measures.Ia12a34, i1, i2, j1, j2, stdDev[15]);
				//*13*/
				I_windowPairs.Ia1324[m][n] = I_mutual(I_measures.Ia1324, i1, i2, j1, j2, stdDev[16]);
				I_windowPairs.I13a24[m][n] = I_mutual(I_measures.I13a24, i1, i2, j1, j2, stdDev[17]);
				I_windowPairs.Ia13a24[m][n] = I_mutual(I_measures.Ia13a24, i1, i2, j1, j2, stdDev[18]);
				//*14*/
				I_windowPairs.Ia1423[m][n] = I_mutual(I_measures.Ia1423, i1, i2, j1, j2, stdDev[19]);
				I_windowPairs.I14a23[m][n] = I_mutual(I_measures.I14a23, i1, i2, j1, j2, stdDev[20]);
				I_windowPairs.Ia14a23[m][n] = I_mutual(I_measures.Ia14a23, i1, i2, j1, j2, stdDev[21]);
			
			}


			cntWinPairs += 1;
	
		}
	}

	if(m!=numberOfWindows){
		printf("Warning: in getInvariantsOnWindows the populated nr of windows was %d while the derived numberOfWindows was %d\n", m, numberOfWindows);
	}
	I_windowPairs.nrOfWindows = m;
	I_windowPairs.nrOfWindowPairs = cntWinPairs;

	*ptr_I_windowPairs = I_windowPairs;

	return 1;

}


int writeInvariantsOnwindowPairs(char *ptr_fileNameOut, struct I_windowPairs_ptr I_windowPairs, int chainLength, int order, int writeOnlyDisjointPairs_b){

	FILE *ptr_fileOut;

	int m = 0, n= 0, i1 = 0, i2 = 0, j1 = 0, j2 = 0;

	int nrOfWindows;

	int print_final_b = 0;

	/*results will be written to this file. Obs: we use a-mode since we 
	are writing multiple times*/
	ptr_fileOut = fopen(ptr_fileNameOut, "a");

	/* write header line first:*/
	fprintf(ptr_fileOut, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;\n", 
			"Line type",
			"PDBfile",
			"structureId",
			"classId",
			"chainId",
			"chainNr",
		    "chainLength",
			"order",
			"nrOfWindows",
			"nrOfWindowPairs",
			"windowLength");
	/* write the header content:*/
	fprintf(ptr_fileOut, "%s;%s;%s;%s;%s;%d;%d;%d;%d;%d;%d;\n", 
			"Header",
			I_windowPairs.fileName,
			I_windowPairs.structureName,
			I_windowPairs.classId,
			I_windowPairs.chainId,
			I_windowPairs.chainNr,
			I_windowPairs.chainLen,
			I_windowPairs.order,
			I_windowPairs.nrOfWindows,
			I_windowPairs.nrOfWindowPairs,
			I_windowPairs.windowLgth
			);


	if (order == 1){
		/* write header line first:*/
		fprintf(ptr_fileOut, "%s;%s;%s;%s;%s;%s;%s;%s;%s;\n", 
				"Values",
				"windowNr1",
				"i1",
				"i2",
				"windowNr2",
				"j1",
				"j2",
				/*order 1*/
				"I12", 
				"Ia12");

		if(print_final_b ==1){

			printf("chain len: %d; window lgth: %d\n", chainLength,I_windowPairs.windowLgth );

		}

		/*write out the invariant values*/
		for (m = 0; m < I_windowPairs.nrOfWindows ; m++){

			/*i1 = I_windowPairs.segIndices_1[m][0];
			i2 = I_windowPairs.segIndices_1[m][1];
			
			if(print_final_b ==1){
				printf("i1: %d, i2:%d\n", i1, i2);
			}*/


			for (n = m + 1; n < I_windowPairs.nrOfWindows ; n++){

				i1 = I_windowPairs.windowPair[m][n].segIndices_1[0];
				i2 = I_windowPairs.windowPair[m][n].segIndices_1[1];
			
				if(print_final_b ==1){
					printf("i1: %d, i2:%d\n", i1, i2);
				}


				/*j1 = I_windowPairs.segIndices_2[n][0];
				j2 = I_windowPairs.segIndices_2[n][1];*/

				j1 = I_windowPairs.windowPair[m][n].segIndices_2[0];
				j2 = I_windowPairs.windowPair[m][n].segIndices_2[1];

				if(writeOnlyDisjointPairs_b == 1){

					if(i2 > j1){ 

						continue;
					}
				}

				fprintf(ptr_fileOut, 
				"%d;%d;%d;%d;%d;%d;%lf;%lf;\n",
				I_windowPairs.windowPair[m][n].windowNr_1,
				i1,
				i2,
				I_windowPairs.windowPair[m][n].windowNr_2,
				j1, 
				j2,
				I_windowPairs.I12[m][n],
				I_windowPairs.Ia12[m][n]
				)
				;

			}
		}
	}

	
	if (order == 2){

		/* write header line first:*/
		fprintf(ptr_fileOut, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;\n", 
				"Values",
				"windowNr1",
				"i1",
				"i2",
				"windowNr2",
				"j1",
				"j2",
				/*order 1*/
				"I12", 
				"Ia12",
				/*order 2 full:*/
				"I1234_full", 
				"I1324_full", 
				/*abs value versions order 2 full:*/
				"Ia12a34_full", 
				"Ia1234_full", 
				"I12a34_full", 
				"Ia13a24_full", 
				"Ia1324_full", 
				"I13a24_full", 
				/*order 2:*/
				"I1234", 
				"I1324", 
				"I1423", 
				/*abs value versions:*/
				/*12*/
				"Ia1234",
				"I12a34",
				"Ia12a34",
				/*13*/
				"Ia1324",
				"I13a24",
				"Ia13a24",
				/*14*/
				"Ia1423",
				"I14a23",
				"Ia14a23");

		if(print_final_b ==1){

			printf("chain len: %d; window lgth: %d\n", chainLength,I_windowPairs.windowLgth );

		}

		/*write out the invariant values*/
		for (m = 0; m < I_windowPairs.nrOfWindows ; m++){

			/*i1 = I_windowPairs.segIndices_1[m][0];
			i2 = I_windowPairs.segIndices_1[m][1];*/
			
			if(print_final_b ==1){
				printf("i1: %d, i2:%d\n", i1, i2);
			}

			for (n = m +1 ; n < I_windowPairs.nrOfWindows ; n++){

				i1 = I_windowPairs.windowPair[m][n].segIndices_1[0];
				i2 = I_windowPairs.windowPair[m][n].segIndices_1[1];

				/*j1 = I_windowPairs.segIndices_2[n][0];
				j2 = I_windowPairs.segIndices_2[n][1];*/

				j1 = I_windowPairs.windowPair[m][n].segIndices_2[0];
				j2 = I_windowPairs.windowPair[m][n].segIndices_2[1];

				if(writeOnlyDisjointPairs_b == 1){

					if(i2 > j1){ 

						continue;
					}
				}

				fprintf(ptr_fileOut, 
				"%d;%d;%d;%d;%d;%d;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;\n",
				I_windowPairs.windowPair[m][n].windowNr_1,
				i1,
				i2,
				I_windowPairs.windowPair[m][n].windowNr_2,
				j1, 
				j2,
				/*order 1*/
				I_windowPairs.I12[m][n],
				I_windowPairs.Ia12[m][n],
				/*order 2, full*/
				I_windowPairs.I1234_full[m][n],
				I_windowPairs.I1324_full[m][n],
				/*abs values versions order 2 full*/
				I_windowPairs.Ia12a34_full[m][n],
				I_windowPairs.Ia1234_full[m][n],
				I_windowPairs.I12a34_full[m][n],
				I_windowPairs.Ia13a24_full[m][n],
				I_windowPairs.Ia1324_full[m][n],
				I_windowPairs.I13a24_full[m][n],
				/*order 2*/
				I_windowPairs.I1234[m][n],
				I_windowPairs.I1324[m][n],
				I_windowPairs.I1423[m][n],
				/*abs values versions*/
				/*12*/
				I_windowPairs.Ia1234[m][n],
				I_windowPairs.I12a34[m][n],
				I_windowPairs.Ia12a34[m][n],
				/*13*/
				I_windowPairs.Ia1324[m][n],
				I_windowPairs.I13a24[m][n],
				I_windowPairs.Ia13a24[m][n],
				/*14*/
				I_windowPairs.Ia1423[m][n],
				I_windowPairs.I14a23[m][n],
				I_windowPairs.Ia14a23[m][n]
				)
				;

			}
		}
	}

	fclose(ptr_fileOut);

	return 1;

}


/*Functions for writing out to files:*/

/*convenience fct's for writing out values:*/
int collectIvalues(struct I_values **ptr_I_values, struct I_ptr I_measures, int subStructureNr, int perturbationNumber, double chainTime, double compTime){

	int order = *I_measures.order;
	int L = *I_measures.chainLen - 1;

	ptr_I_values[subStructureNr][perturbationNumber].fileName = I_measures.fileName;
	/*ptr_I_values[subStructureNr][perturbationNumber].structureName = I_measures.structureName;
	ptr_I_values[subStructureNr][perturbationNumber].classId = I_measures.classId;*/
	ptr_I_values[subStructureNr][perturbationNumber].chainId = I_measures.chainId;
	ptr_I_values[subStructureNr][perturbationNumber].chainNr = *I_measures.chainNr;
	ptr_I_values[subStructureNr][perturbationNumber].chainLen = *I_measures.chainLen;
	ptr_I_values[subStructureNr][perturbationNumber].order = order;
	ptr_I_values[subStructureNr][perturbationNumber].pertNo = perturbationNumber;
	ptr_I_values[subStructureNr][perturbationNumber].cpuTime = chainTime;
	ptr_I_values[subStructureNr][perturbationNumber].aggrTime = compTime;

	if (order >=1){
		ptr_I_values[subStructureNr][perturbationNumber].I12 = I_measures.I12[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].Ia12 = I_measures.Ia12[0][L-1];
	}

	if (order >= 2){
			
		ptr_I_values[subStructureNr][perturbationNumber].I1234 =  I_measures.I1234[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].Ia1234 =  I_measures.Ia1234[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].I12a34 = I_measures.I12a34[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].Ia12a34 = I_measures.Ia12a34[0][L-1];

		ptr_I_values[subStructureNr][perturbationNumber].I1324 =  I_measures.I1324[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].Ia1324 = I_measures.Ia1324[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].I13a24 = I_measures.I13a24[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].Ia13a24 = I_measures.Ia13a24[0][L-1];

		ptr_I_values[subStructureNr][perturbationNumber].I1423 =  I_measures.I1423[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].Ia1423 = I_measures.Ia1423[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].I14a23 = I_measures.I14a23[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].Ia14a23 = I_measures.Ia14a23[0][L-1];

	}

	if (order >=3){

		ptr_I_values[subStructureNr][perturbationNumber].I123456 =  I_measures.I123456[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].I123546 =  I_measures.I123546[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].I123645 =  I_measures.I123645[0][L-1];

		ptr_I_values[subStructureNr][perturbationNumber].I132456 =  I_measures.I132456[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].I132546 =  I_measures.I132546[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].I132645 =  I_measures.I132645[0][L-1];

		ptr_I_values[subStructureNr][perturbationNumber].I142356 =  I_measures.I142356[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].I142536 =  I_measures.I142536[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].I142635 =  I_measures.I142635[0][L-1];

		ptr_I_values[subStructureNr][perturbationNumber].I152346 =  I_measures.I152346[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].I152436 =  I_measures.I152436[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].I152634 =  I_measures.I152634[0][L-1];

		ptr_I_values[subStructureNr][perturbationNumber].I162345 =  I_measures.I162345[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].I162435 =  I_measures.I162435[0][L-1];
		ptr_I_values[subStructureNr][perturbationNumber].I162534 =  I_measures.I162534[0][L-1];
	}

	return 1;
}


/*functions for writing out the results to .txt files*/

int writeChain(char *ptr_fileNameOut, char *fileName, char *structureName, char *chainId, struct cAlpha * ptr_chain, int chainLen){

	FILE *ptr_fileOut;

	int i;

	/*results will be written to this file. We use a-mode since we are calling this writing-fct
	repeatedly (in main):*/
	ptr_fileOut = fopen(ptr_fileNameOut, "a");

	fprintf(ptr_fileOut, "%s;%s;%s;%s;%s;%s\n",
		"PDBfile",
		"structureId",
		"chainId",
		fileName,
		structureName,
		chainId
		);

	//Col names:
	fprintf(ptr_fileOut, "%s;%s;%s;%s;%s;\n",
		"Index",
		"x",
		"y",
		"z",
		"residueNr"
		);

	for (i=0; i < chainLen; i++){
		fprintf(ptr_fileOut,"%d;%lf;%lf;%lf;%d\n",
			i,
			ptr_chain[i].coords.x,
			ptr_chain[i].coords.y,
			ptr_chain[i].coords.z,
			ptr_chain[i].residueNr);


	}


	fclose(ptr_fileOut);

	return 1;

}

int writeIvaluesToFile(char *ptr_fileNameOut, struct I_values **ptr_I_values, int numberOfSubStructures, int numberOfPerturbations){
	
	/*char *fileNameOut, */

	int strCnt = 0; /*structure count*/
	int pertCnt = 0; /*perturbation count*/

	FILE *ptr_fileOut;

	/*results will be written to this file. Obs: we use w-mode so as to clear the file -- we 
	are writing only once*/
	ptr_fileOut = fopen(ptr_fileNameOut, "w");

	/* write header line first:*/
	fprintf(ptr_fileOut, "%s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s\n", 
			"PDB_file", 
			"chainId",
			"chainNr",
		    "chainLength",
			"order",
			"perturbationNr",
			"cpuTime",
			"aggrTime",
			/*order 1*/
			"I12", 
			"Ia12",
			/*order 2:*/
			"I1234", 
			"I1324", 
			"I1423",
			/*abs value versions:*/
			/*12*/
			"Ia1234",
			"I12a34",
			"Ia12a34",
			/*13*/
			"Ia1324",
			"I13a24",
			"Ia13a24",
			/*14*/
			"Ia1423",
			"I14a23",
			"Ia14a23",
			/*order 3*/
			"I123456",
			"I123546",
			"I123645",
			/*13*/
			"I132456",
			"I132546",
			"I132645",
			/*14*/
			"I142356",
			"I142536",
			"I142635",
			/*15*/
			"I152346",
			"I152436",
			"I152634",
			/*16*/
			"I162345",
			"I162435",
			"I162534");


	for (strCnt = 0; strCnt < numberOfSubStructures; strCnt ++){
	
		for (pertCnt = 0; pertCnt < numberOfPerturbations; pertCnt ++){

			//printf("In write-to file: %s; %d; %d; %lf; %lf\n", ptr_I_values[strCnt][pertCnt].fileName, ptr_I_values[strCnt][pertCnt].chainLen, ptr_I_values[strCnt][pertCnt].order, ptr_I_values[strCnt][pertCnt].cpuTime,ptr_I_values[strCnt][pertCnt].I12 );
			//printf("In write-to... : %s\n", ptr_I_values[strCnt][pertCnt].fileName);

			fprintf(ptr_fileOut, 
				"%s;%s;%d;%d;%d;%d;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf\n",
				ptr_I_values[strCnt][pertCnt].fileName,
				ptr_I_values[strCnt][pertCnt].chainId,
				ptr_I_values[strCnt][pertCnt].chainNr,
				ptr_I_values[strCnt][pertCnt].chainLen,
				ptr_I_values[strCnt][pertCnt].order,
				ptr_I_values[strCnt][pertCnt].pertNo,
				ptr_I_values[strCnt][pertCnt].cpuTime,
				ptr_I_values[strCnt][pertCnt].aggrTime,
				/*order 1*/
				ptr_I_values[strCnt][pertCnt].I12,
				ptr_I_values[strCnt][pertCnt].Ia12,
				/*order 2*/
				ptr_I_values[strCnt][pertCnt].I1234,
				ptr_I_values[strCnt][pertCnt].I1324,
				ptr_I_values[strCnt][pertCnt].I1423,
				/*abs 12*/
				ptr_I_values[strCnt][pertCnt].Ia1234,
				ptr_I_values[strCnt][pertCnt].I12a34,
				ptr_I_values[strCnt][pertCnt].Ia12a34,
				/*abs 13*/
				ptr_I_values[strCnt][pertCnt].Ia1324,
				ptr_I_values[strCnt][pertCnt].I13a24,
				ptr_I_values[strCnt][pertCnt].Ia13a24,
				/*abs 14*/
				ptr_I_values[strCnt][pertCnt].Ia1423,
				ptr_I_values[strCnt][pertCnt].I14a23,
				ptr_I_values[strCnt][pertCnt].Ia14a23,
				/*order 3*/
				/*12*/
				ptr_I_values[strCnt][pertCnt].I123456,
				ptr_I_values[strCnt][pertCnt].I123546,
				ptr_I_values[strCnt][pertCnt].I123645,
				/*13*/
				ptr_I_values[strCnt][pertCnt].I132456,
				ptr_I_values[strCnt][pertCnt].I132546,
				ptr_I_values[strCnt][pertCnt].I132645,
				/*14*/
				ptr_I_values[strCnt][pertCnt].I142356,
				ptr_I_values[strCnt][pertCnt].I142536,
				ptr_I_values[strCnt][pertCnt].I142635,
				/*15*/
				ptr_I_values[strCnt][pertCnt].I152346,
				ptr_I_values[strCnt][pertCnt].I152436,
				ptr_I_values[strCnt][pertCnt].I152634,
				/*16*/
				ptr_I_values[strCnt][pertCnt].I162345,
				ptr_I_values[strCnt][pertCnt].I162435,
				ptr_I_values[strCnt][pertCnt].I162534
				);
		}
	}

	//if (order == 1){
	//	fprintf(ptr_fileOut, "%s;%lf;%lf\n",ptr_dirList[fileNr], I_measures.I12[0][L-1]);
	//}

	//if (order == 2){
	//	fprintf(ptr_fileOut, "%s;%lf;%lf;%lf;%lf\n",ptr_dirList[fileNr], I_measures.I12[0][L-1], I_measures.I1234[0][L-1],I_measures.I1324[0][L-1], I_measures.I1423[0][L-1]);

	//			
	//}

	//if (order == 3){
	//	fprintf(ptr_fileOut, "%s;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf\n",
	//		ptr_dirList[fileNr], 
	//		I_measures.I12[0][L-1], 
	//		I_measures.I1234[0][L-1],
	//		I_measures.I1324[0][L-1], 
	//		I_measures.I1423[0][L-1],
	//				
	//		I_measures.I123456[0][L-1],
	//		I_measures.I123645[0][L-1],
	//		I_measures.I123546[0][L-1],
	//		/*13*/
	//		I_measures.I132456[0][L-1],
	//		I_measures.I132546[0][L-1],
	//		I_measures.I132645[0][L-1],
	//		/*14*/
	//		I_measures.I142356[0][L-1],
	//		I_measures.I142536[0][L-1],
	//		I_measures.I142635[0][L-1],
	//		/*15*/
	//		I_measures.I152346[0][L-1],
	//		I_measures.I152436[0][L-1],
	//		I_measures.I152634[0][L-1],
	//		/*16*/
	//		I_measures.I162345[0][L-1],
	//		I_measures.I162435[0][L-1],
	//		I_measures.I162534[0][L-1]);
	//}

	fclose(ptr_fileOut);

	return 1;
}

int writeAllIvaluesToFile(char *ptr_fileNameOut, int full_b, struct I_ptr I_measures, int perturbationNumber){
	 
	int i = 0;
	int j = 0;

	int L = 0;
	int order = *I_measures.order;

	/*char *fileNameOut, */

	FILE *ptr_fileOut;

	/*results will be written to this file. We use a-mode since we are calling this writing-fct
	repeatedly (in main)*/
	ptr_fileOut = fopen(ptr_fileNameOut, "a");

	/* write header line first:*/
	fprintf(ptr_fileOut, "%s;%s;%s;%s;%s;%s;%s\n", 
			"Line type",
			"PDBfile", 
			"chainId",
			"chainNr",
		    "chainLength",
			"order",
			"perturbationNr");
	/* write the header content:*/
	fprintf(ptr_fileOut, "%s;%s;%s;%d;%d;%d;%d\n", 
			"Header",
			I_measures.fileName,
			I_measures.chainId,
			*I_measures.chainNr,
			*I_measures.chainLen,
			*I_measures.order,
			perturbationNumber);

	if (order == 0){
		/* write header for results lines:*/
		fprintf(ptr_fileOut, "%s;%s;%s;%s\n", 
				"Values",
				"i",
				"j",
				/*order 0*/
				"w");
		 
		L = *I_measures.chainLen - 1;

		for (i = 0; i <= L; i++){
			for (j = 0; j <= L; j++){

			fprintf(ptr_fileOut, 
				"%d;%d;%lf\n",
				i,
				j,
				/*order 0*/
				I_measures.wVal[i][j]
				);
			}
		}
	}

	if (order == 1){
		/* write header for results lines:*/
		fprintf(ptr_fileOut, "%s;%s;%s;%s;%s;%s\n", 
				"Values",
				"i",
				"j",
				/*order 0*/
				"w",
				/*order 1*/
				"I12", 
				"Ia12");
		 
		L = *I_measures.chainLen - 1;

		for (i = 0; i <= L; i++){
			for (j = 0; j <= L; j++){

			fprintf(ptr_fileOut, 
				"%d;%d;%lf;%lf;%lf\n",
				i,
				j,
				/*order 0*/
				I_measures.wVal[i][j],
				/*order 1*/
				I_measures.I12[i][j],
				I_measures.Ia12[i][j]
				);
			}
		}
	}

	if (order == 2){
		/* write header for results lines:*/
		fprintf(ptr_fileOut, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s\n", 
				"Values",
				"i",
				"j",
				/*order 0*/
				"w",
				/*order 1*/
				"I12", 
				"Ia12",
				/*order 2:*/
				"I1234", 
				"I1324", 
				"I1423",
				/*abs value versions:*/
				/*12*/
				"Ia1234",
				"I12a34",
				"Ia12a34",
				/*13*/
				"Ia1324",
				"I13a24",
				"Ia13a24",
				/*14*/
				"Ia1423",
				"I14a23",
				"Ia14a23");
		 
		L = *I_measures.chainLen - 1;

		for (i = 0; i <= L; i++){
			for (j = 0; j <= L; j++){

			fprintf(ptr_fileOut, 
				"%d;%d;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf\n",
				i,
				j,
				/*order 0*/
				I_measures.wVal[i][j],
				/*order 1*/
				I_measures.I12[i][j],
				I_measures.Ia12[i][j],
				/*order 2*/
				I_measures.I1234[i][j],
				I_measures.I1324[i][j],
				I_measures.I1423[i][j],
				/*absolute value versions*/
				/*12*/
				I_measures.Ia1234[i][j],
				I_measures.I12a34[i][j],
				I_measures.Ia12a34[i][j],
				/*13*/
				I_measures.Ia1324[i][j],
				I_measures.I13a24[i][j],
				I_measures.Ia13a24[i][j],
				/*14*/
				I_measures.Ia1423[i][j],
				I_measures.I14a23[i][j],
				I_measures.Ia14a23[i][j]
				);
			}
		}
	}

	if (order == 3){
		/* write header for results lines:*/
		fprintf(ptr_fileOut, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s\n", 
				"Values",
				"i",
				"j",
				/*order 0*/
				"w",
				/*order 1*/
				"I12", 
				"Ia12",
				/*order 2:*/
				"I1234", 
				"I1324", 
				"I1423",
				/*abs value versions:*/
				/*12*/
				"Ia1234",
				"I12a34",
				"Ia12a34",
				/*13*/
				"Ia1324",
				"I13a24",
				"Ia13a24",
				/*14*/
				"Ia1423",
				"I14a23",
				"Ia14a23",
				/*order 3*/
				"I123456",
				"I123546",
				"I123645",
				/*13*/
				"I132456",
				"I132546",
				"I132645",
				/*14*/
				"I142356",
				"I142536",
				"I142635",
				/*15*/
				"I152346",
				"I152436",
				"I152634",
				/*16*/
				"I162345",
				"I162435",
				"I162534");
		 
		L = *I_measures.chainLen - 1;

		for (i = 0; i <= L; i++){
			for (j = 0; j <= L; j++){

			fprintf(ptr_fileOut, 
				"%d;%d;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf\n",
				i,
				j,
				/*order 0*/
				I_measures.wVal[i][j],
				/*order 1*/
				I_measures.I12[i][j],
				I_measures.Ia12[i][j],
				/*order 2*/
				I_measures.I1234[i][j],
				I_measures.I1324[i][j],
				I_measures.I1423[i][j],
				/*absolute value versions*/
				I_measures.Ia1234[i][j],
				I_measures.I12a34[i][j],
				I_measures.Ia12a34[i][j],
				/*13*/
				I_measures.Ia1324[i][j],
				I_measures.I13a24[i][j],
				I_measures.Ia13a24[i][j],
				/*14*/
				I_measures.Ia1423[i][j],
				I_measures.I14a23[i][j],
				I_measures.Ia14a23[i][j],
				/*order 3*/
				/*12*/
				I_measures.I123456[i][j],
				I_measures.I123546[i][j],
				I_measures.I123645[i][j],
				/*13*/
				I_measures.I132456[i][j],
				I_measures.I132546[i][j],
				I_measures.I132645[i][j],
				/*14*/
				I_measures.I142356[i][j],
				I_measures.I142536[i][j],
				I_measures.I142635[i][j],
				/*15*/
				I_measures.I152346[i][j],
				I_measures.I152436[i][j],
				I_measures.I152634[i][j],
				/*16*/
				I_measures.I162345[i][j],
				I_measures.I162435[i][j],
				I_measures.I162534[i][j]
				);
			}
		}
	}

	fclose(ptr_fileOut);

	return 1;
}

/*Code blocks for making code in main function (below) more transparent; takes care of memory 
allocation, initialization and freeing memeory:*/
 
int alloc_init_I_values(struct I_values ***ptr_ptr_I_values, int numberOfSubStructures, int numberOfPerturbations){

	int i = 0;
	int j = 0;
	struct I_values **ptr_I_values;

	/*allocate memory to pointer ptr_I_values and initialize:*/
	ptr_I_values = (struct I_values **) malloc (numberOfSubStructures*numberOfPerturbations*sizeof(struct I_values *));
	for (i=0; i< numberOfSubStructures; i++){
		ptr_I_values[i] = (struct I_values *) malloc (numberOfPerturbations*sizeof(struct I_values));
		for (j=0;j<numberOfPerturbations;j++){
			ptr_I_values[i][j].fileName = "NN";
			ptr_I_values[i][j].chainId = "NN";
			ptr_I_values[i][j].chainNr = 0;
			ptr_I_values[i][j].chainLen = 0;
			ptr_I_values[i][j].order = -1;
			ptr_I_values[i][j].pertNo = 0;
			ptr_I_values[i][j].cpuTime = -1.0;
			ptr_I_values[i][j].aggrTime = -1.0;


			ptr_I_values[i][j].I12 = 0.0;
			ptr_I_values[i][j].Ia12 = 0.0;
		
			ptr_I_values[i][j].I1234 = 0.0;
			ptr_I_values[i][j].Ia1234 = 0.0;
			ptr_I_values[i][j].I12a34 = 0.0;
			ptr_I_values[i][j].Ia12a34 = 0.0;

			ptr_I_values[i][j].I1324 = 0.0;
			ptr_I_values[i][j].Ia1324 = 0.0;
			ptr_I_values[i][j].I13a24 = 0.0;
			ptr_I_values[i][j].Ia13a24 = 0.0;

			ptr_I_values[i][j].I1423 = 0.0;
			ptr_I_values[i][j].Ia1423 = 0.0;
			ptr_I_values[i][j].I14a23 = 0.0;
			ptr_I_values[i][j].Ia14a23 = 0.0;

			ptr_I_values[i][j].I123456 = 0.0;
			ptr_I_values[i][j].I123546 = 0.0;
			ptr_I_values[i][j].I123645 = 0.0;

			ptr_I_values[i][j].I132456 = 0.0;
			ptr_I_values[i][j].I132546 = 0.0;
			ptr_I_values[i][j].I132645 = 0.0;

			ptr_I_values[i][j].I142356 = 0.0;
			ptr_I_values[i][j].I142536 = 0.0;
			ptr_I_values[i][j].I142635 = 0.0;

			ptr_I_values[i][j].I152346 = 0.0;
			ptr_I_values[i][j].I152436 = 0.0;
			ptr_I_values[i][j].I152634 = 0.0;

			ptr_I_values[i][j].I162345 = 0.0;
			ptr_I_values[i][j].I162435 = 0.0;
			ptr_I_values[i][j].I162534 = 0.0;
		}
	}

	
	*ptr_ptr_I_values = ptr_I_values;

	return 1;
}

/*Allocation of memory to/init of pointers in struct of double ptr's (dbl arrays) for containing 
invariants' values across the simplex. Could also be written as fct returning an allocated and
init'ed ptr, but here we use the "transform" style:*/ 
int alloc_init_I_measures(struct I_ptr *ptr_I_measures, int order, int full_b, int chainNr, int chainLen, int closed_loops_b, int alreadyAllocatedToChainLen, struct twoSegmentIndex ** ptr_ptr_closedLoopInd){

	int i = 0;
	int j = 0;
	int L = chainLen - 1;
	int cnt = 0; /*only used if closed_loops_b == 1*/

	struct I_ptr I_measures;

	/*if closed_loops_b ==1 this ptr is in use:*/
	struct twoSegmentIndex *ptr_closedLoopInd;

	if(alreadyAllocatedToChainLen == 0){

		I_measures.fileName = (char *) malloc (sizeof(char)*fileNameCharSize);
		I_measures.structureName = (char *) malloc (sizeof(char)*structureNameCharSize);
		I_measures.classId = (char *) malloc (sizeof(char)*classIdCharSize);
		I_measures.chainId = (char *) malloc (sizeof(char)*chainIdCharSize);

		/*init*/
		strcpy(I_measures.fileName,"NN");
		strcpy(I_measures.structureName, "NN");
		strcpy(I_measures.classId,"NN");
		strcpy(I_measures.chainId, "NN");
		

		//I_measures.chainNr = (int *) malloc(sizeof(int));
		I_measures.chainNr = &chainNr;
		//I_measures.chainLen = (int *) malloc(sizeof(int));
		I_measures.chainLen = &chainLen;
		//I_measures.order = (int *) malloc(sizeof(int));
		I_measures.order = &order;

	}
	
	 
	if (order == 0){

		/*allocate memory for double array of w-values */
		if(alreadyAllocatedToChainLen == 0){

			I_measures.wVal = (double **) malloc ((L+1)*sizeof(double *));
			for (i=0; i <= L; i++){
				I_measures.wVal[i] =  (double *)malloc((L+1) * sizeof(double));
				/*initialize to zero everywhere:*/
				for (j=0; j <= L; j++){
				I_measures.wVal[i][j] = 0.0;
				cnt +=1;
				}
			}
		}
		//else reallocate some more memory if necc:
		else if(alreadyAllocatedToChainLen < chainLen){ 

			ptr_I_measures -> wVal  = realloc(ptr_I_measures->wVal, (L+1)*(L+1)*sizeof(double));

			for (i=0; i <= L; i++){

				ptr_I_measures -> wVal[i]  = realloc(ptr_I_measures->wVal[i], (L+1)*sizeof(double));
				/* Initialize to zero everywhere: */
				for (j=0; j <= L; j++){
					ptr_I_measures -> wVal[i][j] = 0.0;
					cnt +=1;
				}
			}
		}

	}
	


	if (order == 1){

		/*allocate memory for double arrays of w and I12 values (order <= 1 measures)*/
		if(alreadyAllocatedToChainLen == 0){
		
			I_measures.wVal = (double **) malloc ((L+1)*(L+1)*sizeof(double ));
			I_measures.I12 = (double **) malloc ((L+1)*(L+1)*sizeof(double));	
			I_measures.Ia12 = (double **) malloc ((L+1)*(L+1)*sizeof(double ));

			for (i=0; i <= L; i++){
				I_measures.wVal[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I12[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.Ia12[i] = (double *)malloc((L+1)*sizeof(double));
				/* Initialize to zero everywhere: */
				for (j=0; j <= L; j++){
					I_measures.wVal[i][j] = 0.0;
					I_measures.I12[i][j] = 0.0;
					I_measures.Ia12[i][j] = 0.0;
					cnt +=1;

				}
			}

		}
		//else reallocate some more memory if necc:
		else if(alreadyAllocatedToChainLen < chainLen){ 

			ptr_I_measures -> wVal  = realloc(ptr_I_measures->wVal, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I12  = realloc(ptr_I_measures->I12, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> Ia12  = realloc(ptr_I_measures->Ia12, (L+1)*(L+1)*sizeof(double));
			//printf("(saaledes stadig ganske hemmelig)?\n");
			for (i=0; i <= L; i++){

				if(i < alreadyAllocatedToChainLen){	
					ptr_I_measures -> wVal[i]  = realloc(ptr_I_measures->wVal[i], (L+1)*sizeof(double));
					ptr_I_measures -> I12[i]  = realloc(ptr_I_measures->I12[i], (L+1)*sizeof(double));
					ptr_I_measures -> Ia12[i]  = realloc(ptr_I_measures->Ia12[i], (L+1)*sizeof(double));
				}
				else if(i >= alreadyAllocatedToChainLen){	
					ptr_I_measures -> wVal[i]  = (double *) malloc((L+1)*sizeof(double));
					ptr_I_measures -> I12[i]  = (double *) malloc((L+1)*sizeof(double));
					ptr_I_measures -> Ia12[i]  = (double *) malloc((L+1)*sizeof(double));
				}

				/*
				ptr_I_measures -> wVal[i]  = realloc(ptr_I_measures->wVal[i], (L+1)*sizeof(double));
				ptr_I_measures -> I12[i]  = realloc(ptr_I_measures->I12[i], (L+1)*sizeof(double));
				ptr_I_measures -> Ia12[i]  = realloc(ptr_I_measures->Ia12[i], (L+1)*sizeof(double));
				*/

				/* Initialize to zero everywhere: */
				for (j=0; j <= L; j++){
					ptr_I_measures -> wVal[i][j] = 0.0;
					ptr_I_measures -> I12[i][j] = 0.0;
					ptr_I_measures -> Ia12[i][j] = 0.0;
					cnt +=1;
				}
			}

		}
		
	}


	if (order == 2){

		/*allocate memory for double arrays of for measures of order <= 2: */
		if(alreadyAllocatedToChainLen == 0){
			
			I_measures.wVal = (double **) malloc ((L+1)*sizeof(double *));

			I_measures.I12 = (double **) malloc ((L+1)*sizeof(double *));

			I_measures.I1234 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I1324 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I1423 = (double **) malloc ((L+1)*sizeof(double *));

			/*absolute value versions*/
			I_measures.Ia12 = (double **) malloc ((L+1)*sizeof(double *));
				
			I_measures.Ia1234 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I12a34 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.Ia12a34 = (double **) malloc ((L+1)*sizeof(double *));
					
			I_measures.Ia1324 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I13a24 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.Ia13a24 = (double **) malloc ((L+1)*sizeof(double *));
					
			I_measures.Ia1423 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I14a23 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.Ia14a23 = (double **) malloc ((L+1)*sizeof(double *));

			/*for "full" order 2*/
			if (full_b == 1){

				I_measures.I1234_full_aid = (double **) malloc ((L+1)*sizeof(double *));
				I_measures.I1324_full_aid = (double **) malloc ((L+1)*sizeof(double *));
				I_measures.I1234_full = (double **) malloc ((L+1)*sizeof(double *));
				I_measures.I1324_full = (double **) malloc ((L+1)*sizeof(double *));

				/*absolute value versions*/
				I_measures.Ia12a34_full_aid = (double **) malloc ((L+1)*sizeof(double *));
				I_measures.Ia13a24_full_aid = (double **) malloc ((L+1)*sizeof(double *));
				I_measures.Ia12a34_full = (double **) malloc ((L+1)*sizeof(double *));
				I_measures.Ia13a24_full = (double **) malloc ((L+1)*sizeof(double *));

				I_measures.Ia1234_full_aid = (double **) malloc ((L+1)*sizeof(double *));
				I_measures.Ia1324_full_aid = (double **) malloc ((L+1)*sizeof(double *));
				I_measures.Ia1234_full = (double **) malloc ((L+1)*sizeof(double *));
				I_measures.Ia1324_full = (double **) malloc ((L+1)*sizeof(double *));

				I_measures.I12a34_full_aid = (double **) malloc ((L+1)*sizeof(double *));
				I_measures.I13a24_full_aid = (double **) malloc ((L+1)*sizeof(double *));
				I_measures.I12a34_full = (double **) malloc ((L+1)*sizeof(double *));
				I_measures.I13a24_full = (double **) malloc ((L+1)*sizeof(double *));

			}


			for (i=0; i <= L; i++){
				
				I_measures.wVal[i] = (double *)malloc((L+1) * sizeof(double));
				
				I_measures.I12[i] = (double *)malloc((L+1) * sizeof(double));

				I_measures.I1234[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I1324[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I1423[i] = (double *)malloc((L+1) * sizeof(double));

				/*absolute value versions */
				I_measures.Ia12[i] = (double *) malloc ((L+1)*sizeof(double));
				
				I_measures.Ia1234[i] = (double *) malloc ((L+1)*sizeof(double));
				I_measures.I12a34[i] = (double *) malloc ((L+1)*sizeof(double));
				I_measures.Ia12a34[i] = (double *) malloc ((L+1)*sizeof(double));
					
				I_measures.Ia1324[i] = (double *) malloc ((L+1)*sizeof(double));
				I_measures.I13a24[i] = (double *) malloc ((L+1)*sizeof(double));
				I_measures.Ia13a24[i] = (double *) malloc ((L+1)*sizeof(double));
					
				I_measures.Ia1423[i] = (double *) malloc ((L+1)*sizeof(double));
				I_measures.I14a23[i] = (double *) malloc ((L+1)*sizeof(double));
				I_measures.Ia14a23[i] = (double *) malloc ((L+1)*sizeof(double));

				/*for "full" order 2*/
				if (full_b == 1){

					I_measures.I1234_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
					I_measures.I1324_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
					I_measures.I1234_full[i] = (double *)malloc((L+1) * sizeof(double));
					I_measures.I1324_full[i] = (double *)malloc((L+1) * sizeof(double));

					/*absolute value versions*/
					I_measures.Ia12a34_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
					I_measures.Ia13a24_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
					I_measures.Ia12a34_full[i] = (double *)malloc((L+1) * sizeof(double));
					I_measures.Ia13a24_full[i] = (double *)malloc((L+1) * sizeof(double));

					I_measures.Ia1234_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
					I_measures.Ia1324_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
					I_measures.Ia1234_full[i] = (double *)malloc((L+1) * sizeof(double));
					I_measures.Ia1324_full[i] = (double *)malloc((L+1) * sizeof(double));

					I_measures.I12a34_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
					I_measures.I13a24_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
					I_measures.I12a34_full[i] = (double *)malloc((L+1) * sizeof(double));
					I_measures.I13a24_full[i] = (double *)malloc((L+1) * sizeof(double));

				}

				/* Initialize to zero everywhere: */
				for (j=0; j <= L; j++){

					I_measures.wVal[i][j] = 0.0;
					
					I_measures.I12[i][j] = 0.0;
					I_measures.Ia12[i][j] = 0.0;

					I_measures.I1234[i][j] = 0.0;
					I_measures.I1324[i][j] = 0.0;
					I_measures.I1423[i][j] = 0.0;

					I_measures.Ia1234[i][j] = 0.0;
					I_measures.I12a34[i][j] = 0.0;
					I_measures.Ia12a34[i][j] = 0.0;

					I_measures.Ia1324[i][j] = 0.0;
					I_measures.I13a24[i][j] = 0.0;
					I_measures.Ia13a24[i][j] = 0.0;
							
					I_measures.Ia1423[i][j] = 0.0;
					I_measures.I14a23[i][j] = 0.0;
					I_measures.Ia14a23[i][j] = 0.0;

					if (full_b == 1){
						I_measures.I1234_full_aid[i][j] = 0.0;
						I_measures.I1324_full_aid[i][j] = 0.0;
						I_measures.I1234_full[i][j] = 0.0;
						I_measures.I1324_full[i][j] = 0.0;

						/*absolute value versions*/
						I_measures.Ia12a34_full_aid[i][j] = 0.0;
						I_measures.Ia13a24_full_aid[i][j] = 0.0;
						I_measures.Ia12a34_full[i][j] = 0.0;
						I_measures.Ia13a24_full[i][j] = 0.0;

						I_measures.Ia1234_full_aid[i][j] = 0.0;
						I_measures.Ia1324_full_aid[i][j] = 0.0;
						I_measures.Ia1234_full[i][j] = 0.0;
						I_measures.Ia1324_full[i][j] = 0.0;

						I_measures.I12a34_full_aid[i][j] = 0.0;
						I_measures.I13a24_full_aid[i][j] = 0.0;
						I_measures.I12a34_full[i][j] = 0.0;
						I_measures.I13a24_full[i][j] = 0.0;
					}
					cnt +=1;
				}
			}
		}
		//else reallocate some more memory if necc:
		else if(alreadyAllocatedToChainLen < chainLen){ 

			ptr_I_measures -> wVal  = realloc(ptr_I_measures->wVal, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I12  = realloc(ptr_I_measures->I12, (L+1)*(L+1)*sizeof(double));
		
			ptr_I_measures -> I1234 = realloc(ptr_I_measures->I1234, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I1324 = realloc(ptr_I_measures->I1324, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I1423 = realloc(ptr_I_measures->I1423, (L+1)*(L+1)*sizeof(double));


			/*absolute value versions*/
			ptr_I_measures -> Ia12  = realloc(ptr_I_measures->Ia12, (L+1)*(L+1)*sizeof(double));
			
			ptr_I_measures -> Ia1234 = realloc(ptr_I_measures->Ia1234, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I12a34 = realloc(ptr_I_measures->I12a34, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> Ia12a34 = realloc(ptr_I_measures->Ia12a34, (L+1)*(L+1)*sizeof(double));
					
			ptr_I_measures -> Ia1324 = realloc(ptr_I_measures->Ia1324, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I13a24 = realloc(ptr_I_measures->I13a24, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> Ia13a24 = realloc(ptr_I_measures->Ia13a24, (L+1)*(L+1)*sizeof(double));
					
			ptr_I_measures -> Ia1423 = realloc(ptr_I_measures->Ia1423, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I14a23 = realloc(ptr_I_measures->I14a23, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> Ia14a23 = realloc(ptr_I_measures->Ia14a23, (L+1)*(L+1)*sizeof(double));

			/*for "full" order 2*/
			if (full_b == 1){

				ptr_I_measures -> I1234_full_aid = realloc(ptr_I_measures->I1234_full_aid, (L+1)*(L+1)*sizeof(double));
				ptr_I_measures -> I1324_full_aid = realloc(ptr_I_measures->I1324_full_aid, (L+1)*(L+1)*sizeof(double));
				ptr_I_measures -> I1234_full = realloc(ptr_I_measures->I1234_full, (L+1)*(L+1)*sizeof(double));
				ptr_I_measures -> I1324_full = realloc(ptr_I_measures->I1324_full, (L+1)*(L+1)*sizeof(double));

				/*absolute value versions*/
				ptr_I_measures -> Ia12a34_full_aid = realloc(ptr_I_measures->Ia12a34_full_aid, (L+1)*(L+1)*sizeof(double));
				ptr_I_measures -> Ia13a24_full_aid = realloc(ptr_I_measures->Ia13a24_full_aid, (L+1)*(L+1)*sizeof(double));
				ptr_I_measures -> Ia12a34_full = realloc(ptr_I_measures->Ia12a34_full, (L+1)*(L+1)*sizeof(double));
				ptr_I_measures -> Ia13a24_full = realloc(ptr_I_measures->Ia13a24_full, (L+1)*(L+1)*sizeof(double));

				ptr_I_measures -> Ia1234_full_aid = realloc(ptr_I_measures->Ia1234_full_aid, (L+1)*(L+1)*sizeof(double));
				ptr_I_measures -> Ia1324_full_aid = realloc(ptr_I_measures->Ia1324_full_aid, (L+1)*(L+1)*sizeof(double));
				ptr_I_measures -> Ia1234_full = realloc(ptr_I_measures->Ia1234_full, (L+1)*(L+1)*sizeof(double));
				ptr_I_measures -> Ia1324_full = realloc(ptr_I_measures->Ia1324_full, (L+1)*(L+1)*sizeof(double));

				ptr_I_measures -> I12a34_full_aid = realloc(ptr_I_measures->I12a34_full_aid, (L+1)*(L+1)*sizeof(double));
				ptr_I_measures -> I13a24_full_aid = realloc(ptr_I_measures->I13a24_full_aid, (L+1)*(L+1)*sizeof(double));
				ptr_I_measures -> I12a34_full = realloc(ptr_I_measures->I12a34_full, (L+1)*(L+1)*sizeof(double));
				ptr_I_measures -> I13a24_full = realloc(ptr_I_measures->I13a24_full, (L+1)*(L+1)*sizeof(double));
			
			}
			for (i=0; i <= L; i++){

				if(i < alreadyAllocatedToChainLen){

					ptr_I_measures -> wVal[i]  = realloc(ptr_I_measures->wVal[i], (L+1)*sizeof(double));
					ptr_I_measures -> I12[i]  = realloc(ptr_I_measures->I12[i], (L+1)*sizeof(double));
				
					ptr_I_measures -> I1234[i] = realloc(ptr_I_measures->I1234[i], (L+1)*sizeof(double));
					ptr_I_measures -> I1324[i] = realloc(ptr_I_measures->I1324[i], (L+1)*sizeof(double));
					ptr_I_measures -> I1423[i] = realloc(ptr_I_measures->I1423[i], (L+1)*sizeof(double));


					/*absolute value versions*/
					ptr_I_measures -> Ia12[i]  = realloc(ptr_I_measures->Ia12[i], (L+1)*sizeof(double));
					
					ptr_I_measures -> Ia1234[i] = realloc(ptr_I_measures->Ia1234[i], (L+1)*sizeof(double));
					ptr_I_measures -> I12a34[i] = realloc(ptr_I_measures->I12a34[i], (L+1)*sizeof(double));
					ptr_I_measures -> Ia12a34[i] = realloc(ptr_I_measures->Ia12a34[i], (L+1)*sizeof(double));
							
					ptr_I_measures -> Ia1324[i] = realloc(ptr_I_measures->Ia1324[i], (L+1)*sizeof(double));
					ptr_I_measures -> I13a24[i] = realloc(ptr_I_measures->I13a24[i], (L+1)*sizeof(double));
					ptr_I_measures -> Ia13a24[i] = realloc(ptr_I_measures->Ia13a24[i], (L+1)*sizeof(double));
							
					ptr_I_measures -> Ia1423[i] = realloc(ptr_I_measures->Ia1423[i], (L+1)*sizeof(double));
					ptr_I_measures -> I14a23[i] = realloc(ptr_I_measures->I14a23[i], (L+1)*sizeof(double));
					ptr_I_measures -> Ia14a23[i] = realloc(ptr_I_measures->Ia14a23[i], (L+1)*sizeof(double));


					/*for "full" order 2*/
					if (full_b == 1){

						ptr_I_measures -> I1234_full_aid[i] = realloc(ptr_I_measures->I1234_full_aid[i], (L+1)*sizeof(double));
						ptr_I_measures -> I1324_full_aid[i] = realloc(ptr_I_measures->I1324_full_aid[i], (L+1)*sizeof(double));
						ptr_I_measures -> I1234_full[i] = realloc(ptr_I_measures->I1234_full[i], (L+1)*sizeof(double));
						ptr_I_measures -> I1324_full[i] = realloc(ptr_I_measures->I1324_full[i], (L+1)*sizeof(double));

						/*absolute value versions*/
						ptr_I_measures -> Ia12a34_full_aid[i] = realloc(ptr_I_measures->Ia12a34_full_aid[i], (L+1)*sizeof(double));
						ptr_I_measures -> Ia13a24_full_aid[i] = realloc(ptr_I_measures->Ia13a24_full_aid[i], (L+1)*sizeof(double));
						ptr_I_measures -> Ia12a34_full[i] = realloc(ptr_I_measures->Ia12a34_full[i], (L+1)*sizeof(double));
						ptr_I_measures -> Ia13a24_full[i] = realloc(ptr_I_measures->Ia13a24_full[i], (L+1)*sizeof(double));

						ptr_I_measures -> Ia1234_full_aid[i] = realloc(ptr_I_measures->Ia1234_full_aid[i], (L+1)*sizeof(double));
						ptr_I_measures -> Ia1324_full_aid[i] = realloc(ptr_I_measures->Ia1324_full_aid[i], (L+1)*sizeof(double));
						ptr_I_measures -> Ia1234_full[i] = realloc(ptr_I_measures->Ia1234_full[i], (L+1)*sizeof(double));
						ptr_I_measures -> Ia1324_full[i] = realloc(ptr_I_measures->Ia1324_full[i], (L+1)*sizeof(double));

						ptr_I_measures -> I12a34_full_aid[i] = realloc(ptr_I_measures->I12a34_full_aid[i], (L+1)*sizeof(double));
						ptr_I_measures -> I13a24_full_aid[i] = realloc(ptr_I_measures->I13a24_full_aid[i], (L+1)*sizeof(double));
						ptr_I_measures -> I12a34_full[i] = realloc(ptr_I_measures->I12a34_full[i], (L+1)*sizeof(double));
						ptr_I_measures -> I13a24_full[i] = realloc(ptr_I_measures->I13a24_full[i], (L+1)*sizeof(double));

					}
				}
				else if(i >= alreadyAllocatedToChainLen){

					ptr_I_measures -> wVal[i]  = (double *) malloc((L+1)*sizeof(double));
					ptr_I_measures -> I12[i]  = (double *) malloc((L+1)*sizeof(double));
				
					ptr_I_measures -> I1234[i] = (double *) malloc((L+1)*sizeof(double));
					ptr_I_measures -> I1324[i] = (double *) malloc((L+1)*sizeof(double));
					ptr_I_measures -> I1423[i] = (double *) malloc((L+1)*sizeof(double));


					/*absolute value versions*/
					ptr_I_measures -> Ia12[i]  = (double *) malloc((L+1)*sizeof(double));
					
					ptr_I_measures -> Ia1234[i] = (double *) malloc((L+1)*sizeof(double));
					ptr_I_measures -> I12a34[i] = (double *) malloc((L+1)*sizeof(double));
					ptr_I_measures -> Ia12a34[i] = (double *) malloc((L+1)*sizeof(double));
							
					ptr_I_measures -> Ia1324[i] = (double *) malloc((L+1)*sizeof(double));
					ptr_I_measures -> I13a24[i] = (double *) malloc((L+1)*sizeof(double));
					ptr_I_measures -> Ia13a24[i] = (double *) malloc((L+1)*sizeof(double));
							
					ptr_I_measures -> Ia1423[i] = (double *) malloc((L+1)*sizeof(double));
					ptr_I_measures -> I14a23[i] = (double *) malloc((L+1)*sizeof(double));
					ptr_I_measures -> Ia14a23[i] = (double *) malloc((L+1)*sizeof(double));


					/*for "full" order 2*/
					if (full_b == 1){

						ptr_I_measures -> I1234_full_aid[i] = (double *) malloc((L+1)*sizeof(double));
						ptr_I_measures -> I1324_full_aid[i] = (double *) malloc((L+1)*sizeof(double));
						ptr_I_measures -> I1234_full[i] = (double *) malloc((L+1)*sizeof(double));
						ptr_I_measures -> I1324_full[i] = (double *) malloc((L+1)*sizeof(double));

						/*absolute value versions*/
						ptr_I_measures -> Ia12a34_full_aid[i] = (double *) malloc((L+1)*sizeof(double));
						ptr_I_measures -> Ia13a24_full_aid[i] = (double *) malloc((L+1)*sizeof(double));
						ptr_I_measures -> Ia12a34_full[i] = (double *) malloc((L+1)*sizeof(double));
						ptr_I_measures -> Ia13a24_full[i] = (double *) malloc((L+1)*sizeof(double));

						ptr_I_measures -> Ia1234_full_aid[i] = (double *) malloc((L+1)*sizeof(double));
						ptr_I_measures -> Ia1324_full_aid[i] = (double *) malloc((L+1)*sizeof(double));
						ptr_I_measures -> Ia1234_full[i] = (double *) malloc((L+1)*sizeof(double));
						ptr_I_measures -> Ia1324_full[i] = (double *) malloc((L+1)*sizeof(double));

						ptr_I_measures -> I12a34_full_aid[i] = (double *) malloc((L+1)*sizeof(double));
						ptr_I_measures -> I13a24_full_aid[i] = (double *) malloc((L+1)*sizeof(double));
						ptr_I_measures -> I12a34_full[i] = (double *) malloc((L+1)*sizeof(double));
						ptr_I_measures -> I13a24_full[i] = (double *) malloc((L+1)*sizeof(double));

					}


				
				}


				/* Initialize to zero everywhere: */
				for (j=0; j <= L; j++){
					
					ptr_I_measures -> wVal[i][j] = 0.0;
					ptr_I_measures -> I12[i][j] = 0.0;

					ptr_I_measures -> I1234[i][j] = 0.0;
					ptr_I_measures -> I1324[i][j] = 0.0;
					ptr_I_measures -> I1423[i][j] = 0.0;

					/*absolute value versions*/
					ptr_I_measures -> Ia12[i][j] = 0.0;

					ptr_I_measures -> Ia1234[i][j] = 0.0;
					ptr_I_measures -> I12a34[i][j] = 0.0;
					ptr_I_measures -> Ia12a34[i][j] = 0.0;

					ptr_I_measures -> Ia1324[i][j] = 0.0;
					ptr_I_measures -> I13a24[i][j] = 0.0;
					ptr_I_measures -> Ia13a24[i][j] = 0.0;

					ptr_I_measures -> Ia1423[i][j] = 0.0;
					ptr_I_measures -> I14a23[i][j] = 0.0;
					ptr_I_measures -> Ia14a23[i][j] = 0.0;

					/*for "full" order 2*/
					if (full_b == 1){
					
						ptr_I_measures -> I1234_full_aid[i][j] = 0.0;
						ptr_I_measures -> I1324_full_aid[i][j] = 0.0;
						ptr_I_measures -> I1234_full[i][j] = 0.0;
						ptr_I_measures -> I1324_full[i][j] = 0.0;

						/*absolute value versions*/
						ptr_I_measures -> Ia12a34_full_aid[i][j] = 0.0;
						ptr_I_measures -> Ia13a24_full_aid[i][j] = 0.0;
						ptr_I_measures -> Ia12a34_full[i][j] = 0.0;
						ptr_I_measures -> Ia13a24_full[i][j] = 0.0;

						ptr_I_measures -> Ia1234_full_aid[i][j] = 0.0;
						ptr_I_measures -> Ia1324_full_aid[i][j] = 0.0;
						ptr_I_measures -> Ia1234_full[i][j] = 0.0;
						ptr_I_measures -> Ia1324_full[i][j] = 0.0;

						ptr_I_measures -> I12a34_full_aid[i][j] = 0.0;
						ptr_I_measures -> I13a24_full_aid[i][j] = 0.0;
						ptr_I_measures -> I12a34_full[i][j] = 0.0;
						ptr_I_measures -> I13a24_full[i][j] = 0.0;
					
					}

					cnt +=1;
				}
			}

		}

		//printf("Langmodigt venter bysvalen stadig stædig længes ... hertil ok ..\n");

	}


	if (order == 3){

		/*allocate memory for double arrays of for measures of order <= 3: */
		if(alreadyAllocatedToChainLen == 0){
			
			I_measures.wVal = (double **) malloc ((L+1)*sizeof(double *));

			/*order 1*/
			I_measures.I12 = (double **) malloc ((L+1)*sizeof(double *));

			/*order 2*/
			I_measures.I1234 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I1324 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I1423 = (double **) malloc ((L+1)*sizeof(double *));

			I_measures.I1234_full_aid = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I1324_full_aid = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I1234_full = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I1324_full = (double **) malloc ((L+1)*sizeof(double *));

			/*assisting*/
			I_measures.I1324_full2_aid = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I1324_full2 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I1423_full0 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I1423_full2_aid = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I1423_full2 = (double **) malloc ((L+1)*sizeof(double *));

			/*absolute value versions*/
			I_measures.Ia12 = (double **) malloc ((L+1)*sizeof(double *));
				
			I_measures.Ia1234 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I12a34 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.Ia12a34 = (double **) malloc ((L+1)*sizeof(double *));
					
			I_measures.Ia1324 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I13a24 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.Ia13a24 = (double **) malloc ((L+1)*sizeof(double *));
					
			I_measures.Ia1423 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I14a23 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.Ia14a23 = (double **) malloc ((L+1)*sizeof(double *));

			/*absolute value versions of "full" order 2's*/
			I_measures.Ia12a34_full_aid = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.Ia13a24_full_aid = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.Ia12a34_full = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.Ia13a24_full = (double **) malloc ((L+1)*sizeof(double *));

			I_measures.Ia1234_full_aid = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.Ia1324_full_aid = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.Ia1234_full = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.Ia1324_full = (double **) malloc ((L+1)*sizeof(double *));

			I_measures.I12a34_full_aid = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I13a24_full_aid = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I12a34_full = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I13a24_full = (double **) malloc ((L+1)*sizeof(double *));

			/*order 3*/
			/*12*/
			I_measures.I123456 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I123645 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I123546 = (double **) malloc ((L+1)*sizeof(double *));
			/*13*/
			I_measures.I132456 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I132546 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I132645 = (double **) malloc ((L+1)*sizeof(double *));
			/*14*/
			I_measures.I142356 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I142536 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I142635 = (double **) malloc ((L+1)*sizeof(double *));
			/*15*/
			I_measures.I152346 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I152436 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I152634 = (double **) malloc ((L+1)*sizeof(double *));
			/*16*/
			I_measures.I162345 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I162435 = (double **) malloc ((L+1)*sizeof(double *));
			I_measures.I162534 = (double **) malloc ((L+1)*sizeof(double *));


			for (i=0; i <= L; i++){
				
				I_measures.wVal[i] = (double *)malloc((L+1) * sizeof(double));
				
				I_measures.I12[i] = (double *)malloc((L+1) * sizeof(double));

				I_measures.I1234[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I1324[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I1423[i] = (double *)malloc((L+1) * sizeof(double));

				I_measures.I1234_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I1324_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I1234_full[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I1324_full[i] = (double *)malloc((L+1) * sizeof(double));

				/*assisting*/
				I_measures.I1324_full2_aid[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I1324_full2[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I1423_full0[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I1423_full2_aid[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I1423_full2[i] = (double *)malloc((L+1) * sizeof(double));

				/*absolute value versions*/
				I_measures.Ia12[i] = (double *) malloc ((L+1)*sizeof(double));
				
				I_measures.Ia1234[i] = (double *) malloc ((L+1)*sizeof(double));
				I_measures.I12a34[i] = (double *) malloc ((L+1)*sizeof(double));
				I_measures.Ia12a34[i] = (double *) malloc ((L+1)*sizeof(double));
					
				I_measures.Ia1324[i] = (double *) malloc ((L+1)*sizeof(double));
				I_measures.I13a24[i] = (double *) malloc ((L+1)*sizeof(double));
				I_measures.Ia13a24[i] = (double *) malloc ((L+1)*sizeof(double));
					
				I_measures.Ia1423[i] = (double *) malloc ((L+1)*sizeof(double));
				I_measures.I14a23[i] = (double *) malloc ((L+1)*sizeof(double));
				I_measures.Ia14a23[i] = (double *) malloc ((L+1)*sizeof(double));

				/*absolute value versions of "full" order 2's*/
				I_measures.Ia12a34_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.Ia13a24_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.Ia12a34_full[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.Ia13a24_full[i] = (double *)malloc((L+1) * sizeof(double));

				I_measures.Ia1234_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.Ia1324_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.Ia1234_full[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.Ia1324_full[i] = (double *)malloc((L+1) * sizeof(double));

				I_measures.I12a34_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I13a24_full_aid[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I12a34_full[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I13a24_full[i] = (double *)malloc((L+1) * sizeof(double));

				/*order 3*/
				/*12*/
				I_measures.I123456[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I123645[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I123546[i] = (double *)malloc((L+1) * sizeof(double));
				/*13*/
				I_measures.I132456[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I132546[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I132645[i] = (double *)malloc((L+1) * sizeof(double));
				/*14*/
				I_measures.I142356[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I142536[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I142635[i] = (double *)malloc((L+1) * sizeof(double));
				/*15*/
				I_measures.I152346[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I152436[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I152634[i] = (double *)malloc((L+1) * sizeof(double));
				/*16*/
				I_measures.I162345[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I162435[i] = (double *)malloc((L+1) * sizeof(double));
				I_measures.I162534[i] = (double *)malloc((L+1) * sizeof(double));

				/* Initialize to zero everywhere: */
				for (j=0; j <= L; j++){

					I_measures.wVal[i][j] = 0.0;
					
					I_measures.I12[i][j] = 0.0;

					I_measures.I1234[i][j] = 0.0;
					I_measures.I1324[i][j] = 0.0;
					I_measures.I1423[i][j] = 0.0;
					
					/*assisting*/
					I_measures.I1324_full2_aid[i][j] = 0.0;
					I_measures.I1324_full2[i][j] = 0.0;
					I_measures.I1423_full0[i][j] = 0.0;
					I_measures.I1423_full2_aid[i][j] = 0.0;
					I_measures.I1423_full2[i][j] = 0.0;

					/*absolute value versions*/
					I_measures.Ia12[i][j] = 0.0;

					I_measures.Ia1234[i][j] = 0.0;
					I_measures.I12a34[i][j] = 0.0;
					I_measures.Ia12a34[i][j] = 0.0;

					I_measures.Ia1324[i][j] = 0.0;
					I_measures.I13a24[i][j] = 0.0;
					I_measures.Ia13a24[i][j] = 0.0;
							
					I_measures.Ia1423[i][j] = 0.0;
					I_measures.I14a23[i][j] = 0.0;
					I_measures.Ia14a23[i][j] = 0.0;


					I_measures.I1234_full_aid[i][j] = 0.0;
					I_measures.I1324_full_aid[i][j] = 0.0;
					I_measures.I1234_full[i][j] = 0.0;
					I_measures.I1324_full[i][j] = 0.0;

					/*absolute value versions of "full" order 2's*/
					I_measures.Ia12a34_full_aid[i][j] = 0.0;
					I_measures.Ia13a24_full_aid[i][j] = 0.0;
					I_measures.Ia12a34_full[i][j] = 0.0;
					I_measures.Ia13a24_full[i][j] = 0.0;

					I_measures.Ia1234_full_aid[i][j] = 0.0;
					I_measures.Ia1324_full_aid[i][j] = 0.0;
					I_measures.Ia1234_full[i][j] = 0.0;
					I_measures.Ia1324_full[i][j] = 0.0;

					I_measures.I12a34_full_aid[i][j] = 0.0;
					I_measures.I13a24_full_aid[i][j] = 0.0;
					I_measures.I12a34_full[i][j] = 0.0;
					I_measures.I13a24_full[i][j] = 0.0;


					/*order 3*/
					/*12*/
					I_measures.I123456[i][j] = 0.0;
					I_measures.I123645[i][j] = 0.0;
					I_measures.I123546[i][j] = 0.0;
					/*13*/
					I_measures.I132456[i][j] = 0.0;
					I_measures.I132546[i][j] = 0.0;
					I_measures.I132645[i][j] = 0.0;
					/*14*/
					I_measures.I142356[i][j] = 0.0;
					I_measures.I142536[i][j] = 0.0;
					I_measures.I142635[i][j] = 0.0;
					/*15*/
					I_measures.I152346[i][j] = 0.0;
					I_measures.I152436[i][j] = 0.0;
					I_measures.I152634[i][j] = 0.0;
					/*16*/
					I_measures.I162345[i][j] = 0.0;
					I_measures.I162435[i][j] = 0.0;
					I_measures.I162534[i][j] = 0.0;

					cnt +=1;

				}

			}

		}
		//else reallocate some more memory if necc:
		else if(alreadyAllocatedToChainLen < chainLen){ 

			ptr_I_measures -> wVal  = realloc(ptr_I_measures->wVal, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I12  = realloc(ptr_I_measures->I12, (L+1)*(L+1)*sizeof(double));
		
			ptr_I_measures -> I1234 = realloc(ptr_I_measures->I1234, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I1324 = realloc(ptr_I_measures->I1324, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I1423 = realloc(ptr_I_measures->I1423, (L+1)*(L+1)*sizeof(double));

			ptr_I_measures -> I1234_full_aid = realloc(ptr_I_measures->I1234_full_aid, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I1324_full_aid = realloc(ptr_I_measures->I1324_full_aid, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I1234_full = realloc(ptr_I_measures->I1234_full, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I1324_full = realloc(ptr_I_measures->I1324_full, (L+1)*(L+1)*sizeof(double));

			/*assisting*/
			ptr_I_measures -> I1324_full2_aid = realloc(ptr_I_measures->I1324_full2_aid, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I1324_full2 = realloc(ptr_I_measures->I1324_full2, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I1423_full0 = realloc(ptr_I_measures->I1423_full0, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I1423_full2_aid = realloc(ptr_I_measures->I1423_full2_aid, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I1423_full2 = realloc(ptr_I_measures->I1423_full2, (L+1)*(L+1)*sizeof(double));


			/*absolute value versions*/
			ptr_I_measures -> Ia12  = realloc(ptr_I_measures->Ia12, (L+1)*(L+1)*sizeof(double));
			
			ptr_I_measures -> Ia1234 = realloc(ptr_I_measures->Ia1234, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I12a34 = realloc(ptr_I_measures->I12a34, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> Ia12a34 = realloc(ptr_I_measures->Ia12a34, (L+1)*(L+1)*sizeof(double));
					
			ptr_I_measures -> Ia1324 = realloc(ptr_I_measures->Ia1324, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I13a24 = realloc(ptr_I_measures->I13a24, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> Ia13a24 = realloc(ptr_I_measures->Ia13a24, (L+1)*(L+1)*sizeof(double));
					
			ptr_I_measures -> Ia1423 = realloc(ptr_I_measures->Ia1423, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I14a23 = realloc(ptr_I_measures->I14a23, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> Ia14a23 = realloc(ptr_I_measures->Ia14a23, (L+1)*(L+1)*sizeof(double));

			/*absolute value versions of "full" order 2's*/
			ptr_I_measures -> Ia12a34_full_aid = realloc(ptr_I_measures->Ia12a34_full_aid, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> Ia13a24_full_aid = realloc(ptr_I_measures->Ia13a24_full_aid, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> Ia12a34_full = realloc(ptr_I_measures->Ia12a34_full, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> Ia13a24_full = realloc(ptr_I_measures->Ia13a24_full, (L+1)*(L+1)*sizeof(double));

			ptr_I_measures -> Ia1234_full_aid = realloc(ptr_I_measures->Ia1234_full_aid, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> Ia1324_full_aid = realloc(ptr_I_measures->Ia1324_full_aid, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> Ia1234_full = realloc(ptr_I_measures->Ia1234_full, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> Ia1324_full = realloc(ptr_I_measures->Ia1324_full, (L+1)*(L+1)*sizeof(double));

			ptr_I_measures -> I12a34_full_aid = realloc(ptr_I_measures->I12a34_full_aid, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I13a24_full_aid = realloc(ptr_I_measures->I13a24_full_aid, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I12a34_full = realloc(ptr_I_measures->I12a34_full, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I13a24_full = realloc(ptr_I_measures->I13a24_full, (L+1)*(L+1)*sizeof(double));

			/*order 3*/
			/*12*/
			ptr_I_measures -> I123456 = realloc(ptr_I_measures->I123456, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I123645 = realloc(ptr_I_measures->I123645, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I123546 = realloc(ptr_I_measures->I123546, (L+1)*(L+1)*sizeof(double));
			/*13*/
			ptr_I_measures -> I132456 = realloc(ptr_I_measures->I132456, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I132546 = realloc(ptr_I_measures->I132546, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I132645 = realloc(ptr_I_measures->I132645, (L+1)*(L+1)*sizeof(double));
			/*14*/
			ptr_I_measures -> I142356 = realloc(ptr_I_measures->I142356, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I142536 = realloc(ptr_I_measures->I142536, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I142635 = realloc(ptr_I_measures->I142635, (L+1)*(L+1)*sizeof(double));
			/*15*/
			ptr_I_measures -> I152346 = realloc(ptr_I_measures->I152346, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I152436 = realloc(ptr_I_measures->I152436, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I152634 = realloc(ptr_I_measures->I152634, (L+1)*(L+1)*sizeof(double));
			/*16*/
			ptr_I_measures -> I162345 = realloc(ptr_I_measures->I162345, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I162435 = realloc(ptr_I_measures->I162435, (L+1)*(L+1)*sizeof(double));
			ptr_I_measures -> I162534 = realloc(ptr_I_measures->I162534, (L+1)*(L+1)*sizeof(double));


			for (i=0; i <= L; i++){

				if(i < alreadyAllocatedToChainLen){

					ptr_I_measures -> wVal[i]  = realloc(ptr_I_measures->wVal[i], (L+1)*sizeof(double));
					ptr_I_measures -> I12[i]  = realloc(ptr_I_measures->I12[i], (L+1)*sizeof(double));
				
					ptr_I_measures -> I1234[i] = realloc(ptr_I_measures->I1234[i], (L+1)*sizeof(double));
					ptr_I_measures -> I1324[i] = realloc(ptr_I_measures->I1324[i], (L+1)*sizeof(double));
					ptr_I_measures -> I1423[i] = realloc(ptr_I_measures->I1423[i], (L+1)*sizeof(double));

					ptr_I_measures -> I1234_full_aid[i] = realloc(ptr_I_measures->I1234_full_aid[i], (L+1)*sizeof(double));
					ptr_I_measures -> I1324_full_aid[i] = realloc(ptr_I_measures->I1324_full_aid[i], (L+1)*sizeof(double));
					ptr_I_measures -> I1234_full[i] = realloc(ptr_I_measures->I1234_full[i], (L+1)*sizeof(double));
					ptr_I_measures -> I1324_full[i] = realloc(ptr_I_measures->I1324_full[i], (L+1)*sizeof(double));

					/*assisting*/
					ptr_I_measures -> I1324_full2_aid = realloc(ptr_I_measures->I1324_full2_aid, (L+1)*(L+1)*sizeof(double));
					ptr_I_measures -> I1324_full2 = realloc(ptr_I_measures->I1324_full2, (L+1)*(L+1)*sizeof(double));
					ptr_I_measures -> I1423_full0 = realloc(ptr_I_measures->I1423_full0, (L+1)*(L+1)*sizeof(double));
					ptr_I_measures -> I1423_full2_aid = realloc(ptr_I_measures->I1423_full2_aid, (L+1)*(L+1)*sizeof(double));
					ptr_I_measures -> I1423_full2 = realloc(ptr_I_measures->I1423_full2, (L+1)*(L+1)*sizeof(double));

					/*absolute value versions*/
					ptr_I_measures -> Ia12[i]  = realloc(ptr_I_measures->Ia12[i], (L+1)*sizeof(double));
					
					ptr_I_measures -> Ia1234[i] = realloc(ptr_I_measures->Ia1234[i], (L+1)*sizeof(double));
					ptr_I_measures -> I12a34[i] = realloc(ptr_I_measures->I12a34[i], (L+1)*sizeof(double));
					ptr_I_measures -> Ia12a34[i] = realloc(ptr_I_measures->Ia12a34[i], (L+1)*sizeof(double));
							
					ptr_I_measures -> Ia1324[i] = realloc(ptr_I_measures->Ia1324[i], (L+1)*sizeof(double));
					ptr_I_measures -> I13a24[i] = realloc(ptr_I_measures->I13a24[i], (L+1)*sizeof(double));
					ptr_I_measures -> Ia13a24[i] = realloc(ptr_I_measures->Ia13a24[i], (L+1)*sizeof(double));
							
					ptr_I_measures -> Ia1423[i] = realloc(ptr_I_measures->Ia1423[i], (L+1)*sizeof(double));
					ptr_I_measures -> I14a23[i] = realloc(ptr_I_measures->I14a23[i], (L+1)*sizeof(double));
					ptr_I_measures -> Ia14a23[i] = realloc(ptr_I_measures->Ia14a23[i], (L+1)*sizeof(double));

					ptr_I_measures -> Ia12a34_full_aid[i] = realloc(ptr_I_measures->Ia12a34_full_aid[i], (L+1)*sizeof(double));
					ptr_I_measures -> Ia13a24_full_aid[i] = realloc(ptr_I_measures->Ia13a24_full_aid[i], (L+1)*sizeof(double));
					ptr_I_measures -> Ia12a34_full[i] = realloc(ptr_I_measures->Ia12a34_full[i], (L+1)*sizeof(double));
					ptr_I_measures -> Ia13a24_full[i] = realloc(ptr_I_measures->Ia13a24_full[i], (L+1)*sizeof(double));

					ptr_I_measures -> Ia1234_full_aid[i] = realloc(ptr_I_measures->Ia1234_full_aid[i], (L+1)*sizeof(double));
					ptr_I_measures -> Ia1324_full_aid[i] = realloc(ptr_I_measures->Ia1324_full_aid[i], (L+1)*sizeof(double));
					ptr_I_measures -> Ia1234_full[i] = realloc(ptr_I_measures->Ia1234_full[i], (L+1)*sizeof(double));
					ptr_I_measures -> Ia1324_full[i] = realloc(ptr_I_measures->Ia1324_full[i], (L+1)*sizeof(double));

					ptr_I_measures -> I12a34_full_aid[i] = realloc(ptr_I_measures->I12a34_full_aid[i], (L+1)*sizeof(double));
					ptr_I_measures -> I13a24_full_aid[i] = realloc(ptr_I_measures->I13a24_full_aid[i], (L+1)*sizeof(double));
					ptr_I_measures -> I12a34_full[i] = realloc(ptr_I_measures->I12a34_full[i], (L+1)*sizeof(double));
					ptr_I_measures -> I13a24_full[i] = realloc(ptr_I_measures->I13a24_full[i], (L+1)*sizeof(double));

					/*order 3*/
					/*12*/
					ptr_I_measures -> I123456[i] = realloc(ptr_I_measures->I123456[i], (L+1)*sizeof(double));
					ptr_I_measures -> I123645[i] = realloc(ptr_I_measures->I123645[i], (L+1)*sizeof(double));
					ptr_I_measures -> I123546[i] = realloc(ptr_I_measures->I123546[i], (L+1)*sizeof(double));
					/*13*/
					ptr_I_measures -> I132456[i] = realloc(ptr_I_measures->I132456[i], (L+1)*sizeof(double));
					ptr_I_measures -> I132546[i] = realloc(ptr_I_measures->I132546[i], (L+1)*sizeof(double));
					ptr_I_measures -> I132645[i] = realloc(ptr_I_measures->I132645[i], (L+1)*sizeof(double));
					/*14*/
					ptr_I_measures -> I142356[i] = realloc(ptr_I_measures->I142356[i], (L+1)*sizeof(double));
					ptr_I_measures -> I142536[i] = realloc(ptr_I_measures->I142536[i], (L+1)*sizeof(double));
					ptr_I_measures -> I142635[i] = realloc(ptr_I_measures->I142635[i], (L+1)*sizeof(double));
					/*15*/
					ptr_I_measures -> I152346[i] = realloc(ptr_I_measures->I152346[i], (L+1)*sizeof(double));
					ptr_I_measures -> I152436[i] = realloc(ptr_I_measures->I152436[i], (L+1)*sizeof(double));
					ptr_I_measures -> I152634[i] = realloc(ptr_I_measures->I152634[i], (L+1)*sizeof(double));
					/*16*/
					ptr_I_measures -> I162345[i] = realloc(ptr_I_measures->I162345[i], (L+1)*sizeof(double));
					ptr_I_measures -> I162435[i] = realloc(ptr_I_measures->I162435[i], (L+1)*sizeof(double));
					ptr_I_measures -> I162534[i] = realloc(ptr_I_measures->I162534[i], (L+1)*sizeof(double));

				}
				else if(i >= alreadyAllocatedToChainLen){

					ptr_I_measures -> wVal[i]  = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I12[i]  = (double *)malloc((L+1)*sizeof(double));
				
					ptr_I_measures -> I1234[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I1324[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I1423[i] = (double *)malloc((L+1)*sizeof(double));

					ptr_I_measures -> I1234_full_aid[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I1324_full_aid[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I1234_full[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I1324_full[i] = (double *)malloc((L+1)*sizeof(double));

					/*assisting*/
					ptr_I_measures -> I1324_full2_aid[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I1324_full2[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I1423_full0[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I1423_full2_aid[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I1423_full2[i] = (double *)malloc((L+1)*sizeof(double));

					/*absolute value versions*/
					ptr_I_measures -> Ia12[i]  = (double *)malloc((L+1)*sizeof(double));
					
					ptr_I_measures -> Ia1234[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I12a34[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> Ia12a34[i] = (double *)malloc((L+1)*sizeof(double));
							
					ptr_I_measures -> Ia1324[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I13a24[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> Ia13a24[i] = (double *)malloc((L+1)*sizeof(double));
							
					ptr_I_measures -> Ia1423[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I14a23[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> Ia14a23[i] = (double *)malloc((L+1)*sizeof(double));

					ptr_I_measures -> Ia12a34_full_aid[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> Ia13a24_full_aid[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> Ia12a34_full[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> Ia13a24_full[i] = (double *)malloc((L+1)*sizeof(double));

					ptr_I_measures -> Ia1234_full_aid[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> Ia1324_full_aid[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> Ia1234_full[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> Ia1324_full[i] = (double *)malloc((L+1)*sizeof(double));

					ptr_I_measures -> I12a34_full_aid[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I13a24_full_aid[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I12a34_full[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I13a24_full[i] = (double *)malloc((L+1)*sizeof(double));

					/*order 3*/
					/*12*/
					ptr_I_measures -> I123456[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I123645[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I123546[i] = (double *)malloc((L+1)*sizeof(double));
					/*13*/
					ptr_I_measures -> I132456[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I132546[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I132645[i] = (double *)malloc((L+1)*sizeof(double));
					/*14*/
					ptr_I_measures -> I142356[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I142536[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I142635[i] = (double *)malloc((L+1)*sizeof(double));
					/*15*/
					ptr_I_measures -> I152346[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I152436[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I152634[i] = (double *)malloc((L+1)*sizeof(double));
					/*16*/
					ptr_I_measures -> I162345[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I162435[i] = (double *)malloc((L+1)*sizeof(double));
					ptr_I_measures -> I162534[i] = (double *)malloc((L+1)*sizeof(double));


				}



				/* Initialize to zero everywhere: */
				for (j=0; j <= L; j++){
					
					ptr_I_measures -> wVal[i][j] = 0.0;
					ptr_I_measures -> I12[i][j] = 0.0;

					ptr_I_measures -> I1234[i][j] = 0.0;
					ptr_I_measures -> I1324[i][j] = 0.0;
					ptr_I_measures -> I1423[i][j] = 0.0;

					ptr_I_measures -> I1234_full_aid[i][j] = 0.0;
					ptr_I_measures -> I1324_full_aid[i][j] = 0.0;
					ptr_I_measures -> I1234_full[i][j] = 0.0;
					ptr_I_measures -> I1324_full[i][j] = 0.0;

					/*assisting*/
					ptr_I_measures -> I1324_full2_aid[i][j] = 0.0;
					ptr_I_measures -> I1324_full2[i][j] = 0.0;
					ptr_I_measures -> I1423_full0[i][j] = 0.0;
					ptr_I_measures -> I1423_full2_aid[i][j] = 0.0;
					ptr_I_measures -> I1423_full2[i][j] = 0.0;


					/*absolute value versions*/
					ptr_I_measures -> Ia12[i][j] = 0.0;

					ptr_I_measures -> Ia1234[i][j] = 0.0;
					ptr_I_measures -> I12a34[i][j] = 0.0;
					ptr_I_measures -> Ia12a34[i][j] = 0.0;

					ptr_I_measures -> Ia1324[i][j] = 0.0;
					ptr_I_measures -> I13a24[i][j] = 0.0;
					ptr_I_measures -> Ia13a24[i][j] = 0.0;

					ptr_I_measures -> Ia1423[i][j] = 0.0;
					ptr_I_measures -> I14a23[i][j] = 0.0;
					ptr_I_measures -> Ia14a23[i][j] = 0.0;
				
					ptr_I_measures -> Ia12a34_full_aid[i][j] = 0.0;
					ptr_I_measures -> Ia13a24_full_aid[i][j] = 0.0;
					ptr_I_measures -> Ia12a34_full[i][j] = 0.0;
					ptr_I_measures -> Ia13a24_full[i][j] = 0.0;

					ptr_I_measures -> Ia1234_full_aid[i][j] = 0.0;
					ptr_I_measures -> Ia1324_full_aid[i][j] = 0.0;
					ptr_I_measures -> Ia1234_full[i][j] = 0.0;
					ptr_I_measures -> Ia1324_full[i][j] = 0.0;

					ptr_I_measures -> I12a34_full_aid[i][j] = 0.0;
					ptr_I_measures -> I13a24_full_aid[i][j] = 0.0;
					ptr_I_measures -> I12a34_full[i][j] = 0.0;
					ptr_I_measures -> I13a24_full[i][j] = 0.0;						


					ptr_I_measures -> I123456[i][j] = 0.0;
					ptr_I_measures -> I123645[i][j] = 0.0;
					ptr_I_measures -> I123546[i][j] = 0.0;
					/*13*/
					ptr_I_measures -> I132456[i][j] = 0.0;
					ptr_I_measures -> I132546[i][j] = 0.0;
					ptr_I_measures -> I132645[i][j] = 0.0;
					/*14*/
					ptr_I_measures -> I142356[i][j] = 0.0;
					ptr_I_measures -> I142536[i][j] = 0.0;
					ptr_I_measures -> I142635[i][j] = 0.0;
					/*15*/
					ptr_I_measures -> I152346[i][j] = 0.0;
					ptr_I_measures -> I152436[i][j] = 0.0;
					ptr_I_measures -> I152634[i][j] = 0.0;
					/*16*/
					ptr_I_measures -> I162345[i][j] = 0.0;
					ptr_I_measures -> I162435[i][j] = 0.0;
					ptr_I_measures -> I162534[i][j] = 0.0;

					cnt +=1;
				}
			}

		}

	}

	if(alreadyAllocatedToChainLen == 0){
		*ptr_I_measures = I_measures;
	}

	if (closed_loops_b ==1){

		ptr_closedLoopInd = (struct twoSegmentIndex *)calloc(cnt , sizeof(struct twoSegmentIndex));

		*ptr_ptr_closedLoopInd = ptr_closedLoopInd;

	}

	return 1;
}

/*Needed?*/
/*allocate ptr's to segment coord arrays; initialize host side arrays to zero:*/ 
int alloc_init_segment_ptr(struct Segment_ptr *ptr_Segment_ptr, int chainLen, int simplexCnt){

	int i = 0;
	int j = 0;
	int L = chainLen - 1;
	int cnt = 0;

	struct Segment_ptr Coord_arrays;

	size_t size = simplexCnt * sizeof(double);

	/*allocate*/
	Coord_arrays.ptr_segmentPairCoords_s101 = (double *) calloc (simplexCnt, sizeof(double));
	Coord_arrays.ptr_segmentPairCoords_s102 = (double *) calloc (simplexCnt, sizeof(double));
	Coord_arrays.ptr_segmentPairCoords_s103 = (double *) calloc (simplexCnt, sizeof(double));

	Coord_arrays.ptr_segmentPairCoords_s111 = (double *) calloc (simplexCnt, sizeof(double));
	Coord_arrays.ptr_segmentPairCoords_s112 = (double *) calloc (simplexCnt, sizeof(double));
	Coord_arrays.ptr_segmentPairCoords_s113 = (double *) calloc (simplexCnt, sizeof(double));

	Coord_arrays.ptr_segmentPairCoords_s201 = (double *) calloc (simplexCnt, sizeof(double));
	Coord_arrays.ptr_segmentPairCoords_s202 = (double *) calloc (simplexCnt, sizeof(double));
	Coord_arrays.ptr_segmentPairCoords_s203 = (double *) calloc (simplexCnt, sizeof(double));

	Coord_arrays.ptr_segmentPairCoords_s211 = (double *) calloc (simplexCnt, sizeof(double));
	Coord_arrays.ptr_segmentPairCoords_s212 = (double *) calloc (simplexCnt, sizeof(double));
	Coord_arrays.ptr_segmentPairCoords_s213 = (double *) calloc (simplexCnt, sizeof(double));

	/*init the arrays with zeros:*/
	for (cnt = 0; cnt < simplexCnt; cnt++){
		/*start point of seg i, x,y,z:*/
		Coord_arrays.ptr_segmentPairCoords_s101[cnt] = 0.0;
		Coord_arrays.ptr_segmentPairCoords_s102[cnt] = 0.0;
		Coord_arrays.ptr_segmentPairCoords_s103[cnt] = 0.0;
		/*end point of seg i, x,y,z:*/
		Coord_arrays.ptr_segmentPairCoords_s111[cnt] = 0.0;
		Coord_arrays.ptr_segmentPairCoords_s112[cnt] = 0.0;
		Coord_arrays.ptr_segmentPairCoords_s113[cnt] = 0.0;
		/*start point of seg j, x,y,z:*/
		Coord_arrays.ptr_segmentPairCoords_s201[cnt] = 0.0;
		Coord_arrays.ptr_segmentPairCoords_s202[cnt] = 0.0;
		Coord_arrays.ptr_segmentPairCoords_s203[cnt] = 0.0;
		/*start point of seg j, x,y,z:*/
		Coord_arrays.ptr_segmentPairCoords_s211[cnt] = 0.0;
		Coord_arrays.ptr_segmentPairCoords_s212[cnt] = 0.0;
		Coord_arrays.ptr_segmentPairCoords_s213[cnt] = 0.0;

	}

	/*transform the input: set to the allocated copy:*/
	*ptr_Segment_ptr = Coord_arrays;

	return 1;
}


/*initialize ptr's to segment coord arrays to the values given by an input segment chain:*/
int init_segment_ptr(struct Segment_ptr Coord_arrays, struct segment *ptr_segment, int chainLen, int *ptr_simplexCnt){

	int i = 0;
	int j = 0;
	int L = chainLen - 1;
	int cnt = 0;

	/*populate the array:*/
	for ( i=L-1; i>=0; i--){
		for ( j=i+1; j<=L-1; j++ ){
			/*start point of seg i, x,y,z:*/
			Coord_arrays.ptr_segmentPairCoords_s101[cnt] = ptr_segment[i].s1.x;
			Coord_arrays.ptr_segmentPairCoords_s102[cnt] = ptr_segment[i].s1.y;
			Coord_arrays.ptr_segmentPairCoords_s103[cnt] = ptr_segment[i].s1.z;
			/*end point of seg i, x,y,z:*/
			Coord_arrays.ptr_segmentPairCoords_s111[cnt] = ptr_segment[i].s2.x;
			Coord_arrays.ptr_segmentPairCoords_s112[cnt] = ptr_segment[i].s2.y;
			Coord_arrays.ptr_segmentPairCoords_s113[cnt] = ptr_segment[i].s2.z;
			/*start point of seg j, x,y,z:*/
			Coord_arrays.ptr_segmentPairCoords_s201[cnt] = ptr_segment[j].s1.x;
			Coord_arrays.ptr_segmentPairCoords_s202[cnt] = ptr_segment[j].s1.y;
			Coord_arrays.ptr_segmentPairCoords_s203[cnt] = ptr_segment[j].s1.z;
			/*start point of seg j, x,y,z:*/
			Coord_arrays.ptr_segmentPairCoords_s211[cnt] = ptr_segment[j].s2.x;
			Coord_arrays.ptr_segmentPairCoords_s212[cnt] = ptr_segment[j].s2.y;
			Coord_arrays.ptr_segmentPairCoords_s213[cnt] = ptr_segment[j].s2.z;

			
			cnt += 1;
		}
	}

	*ptr_simplexCnt = cnt; /*keep the simplex count: the number of vertices in the simplex */

	return 1;
}


int init_I_measures(struct I_ptr I_measures, int order, int full_b, int chainLen){

	int i = 0;
	int j = 0;
	int L = chainLen - 1;

	/*(RE)INITIALIZATION*/

	if (order == 0){

		/*allocate memory for double array of w-values */
		for (i=0; i <= L; i++){
			//I_measures.wVal[i] = (double *)realloc(I_measures.wVal[i], (L+1) * sizeof(double)); /* (double *)malloc((L+1) * sizeof(double));*/
			/* Initialize to zero everywhere: */
			for (j=0; j <= L; j++){
				I_measures.wVal[i][j] = 0.0;
			}
		}
	}


	if (order == 1){

		for (i=0; i <= L; i++){

			/* Initialize to zero everywhere: */
			for (j=0; j <= L; j++){
			 	I_measures.wVal[i][j] = 0.0;
				I_measures.I12[i][j] = 0.0;
				I_measures.Ia12[i][j] = 0.0;

			}
		}
	}


	if (order == 2){

		for (i=0; i <= L; i++){
			 
			/* Initialize to zero everywhere: */
			for (j=0; j <= L; j++){

			 	I_measures.wVal[i][j] = 0.0;
				
				I_measures.I12[i][j] = 0.0;
				I_measures.Ia12[i][j] = 0.0;

				I_measures.I1234[i][j] = 0.0;
				I_measures.I1324[i][j] = 0.0;
				I_measures.I1423[i][j] = 0.0;

				I_measures.Ia1234[i][j] = 0.0;
				I_measures.I12a34[i][j] = 0.0;
				I_measures.Ia12a34[i][j] = 0.0;

				I_measures.Ia1324[i][j] = 0.0;
				I_measures.I13a24[i][j] = 0.0;
				I_measures.Ia13a24[i][j] = 0.0;
						
				I_measures.Ia1423[i][j] = 0.0;
				I_measures.I14a23[i][j] = 0.0;
				I_measures.Ia14a23[i][j] = 0.0;

				if (full_b == 1){
					I_measures.I1234_full_aid[i][j] = 0.0;
					I_measures.I1324_full_aid[i][j] = 0.0;
					I_measures.I1234_full[i][j] = 0.0;
					I_measures.I1324_full[i][j] = 0.0;

					/*absolute value versions*/
					I_measures.Ia12a34_full_aid[i][j] = 0.0;
					I_measures.Ia13a24_full_aid[i][j] = 0.0;
					I_measures.Ia12a34_full[i][j] = 0.0;
					I_measures.Ia13a24_full[i][j] = 0.0;

					I_measures.Ia1234_full_aid[i][j] = 0.0;
					I_measures.Ia1324_full_aid[i][j] = 0.0;
					I_measures.Ia1234_full[i][j] = 0.0;
					I_measures.Ia1324_full[i][j] = 0.0;

					I_measures.I12a34_full_aid[i][j] = 0.0;
					I_measures.I13a24_full_aid[i][j] = 0.0;
					I_measures.I12a34_full[i][j] = 0.0;
					I_measures.I13a24_full[i][j] = 0.0;
					
				}
			}
		}
	}


	if (order == 3){

		for (i=0; i <= L; i++){

			/* Initialize to zero everywhere: */
			for (j=0; j <= L; j++){

			 	I_measures.wVal[i][j] = 0.0;
				
				I_measures.I12[i][j] = 0.0;
				I_measures.Ia12[i][j] = 0.0;

				I_measures.I1234[i][j] = 0.0;
				I_measures.I1324[i][j] = 0.0;
				I_measures.I1423[i][j] = 0.0;
				
				I_measures.Ia1234[i][j] = 0.0;
				I_measures.I12a34[i][j] = 0.0;
				I_measures.Ia12a34[i][j] = 0.0;

				I_measures.Ia1324[i][j] = 0.0;
				I_measures.I13a24[i][j] = 0.0;
				I_measures.Ia13a24[i][j] = 0.0;
						
				I_measures.Ia1423[i][j] = 0.0;
				I_measures.I14a23[i][j] = 0.0;
				I_measures.Ia14a23[i][j] = 0.0;


				I_measures.I1234_full_aid[i][j] = 0.0;
				I_measures.I1324_full_aid[i][j] = 0.0;
				I_measures.I1234_full[i][j] = 0.0;
				I_measures.I1324_full[i][j] = 0.0;

				I_measures.I1324_full2_aid[i][j] = 0.0;
				I_measures.I1324_full2[i][j] = 0.0;
				I_measures.I1423_full0[i][j] = 0.0;
				I_measures.I1423_full2_aid[i][j] = 0.0;
				I_measures.I1423_full2[i][j] = 0.0;

				/*absolute value versions of "full" order 2's*/
				I_measures.Ia12a34_full_aid[i][j] = 0.0;
				I_measures.Ia13a24_full_aid[i][j] = 0.0;
				I_measures.Ia12a34_full[i][j] = 0.0;
				I_measures.Ia13a24_full[i][j] = 0.0;

				I_measures.Ia1234_full_aid[i][j] = 0.0;
				I_measures.Ia1324_full_aid[i][j] = 0.0;
				I_measures.Ia1234_full[i][j] = 0.0;
				I_measures.Ia1324_full[i][j] = 0.0;

				I_measures.I12a34_full_aid[i][j] = 0.0;
				I_measures.I13a24_full_aid[i][j] = 0.0;
				I_measures.I12a34_full[i][j] = 0.0;
				I_measures.I13a24_full[i][j] = 0.0;

				/*order 3*/
				/*12*/
				I_measures.I123456[i][j] = 0.0;
				I_measures.I123645[i][j] = 0.0;
				I_measures.I123546[i][j] = 0.0;
				/*13*/
				I_measures.I132456[i][j] = 0.0;
				I_measures.I132546[i][j] = 0.0;
				I_measures.I132645[i][j] = 0.0;
				/*14*/
				I_measures.I142356[i][j] = 0.0;
				I_measures.I142536[i][j] = 0.0;
				I_measures.I142635[i][j] = 0.0;
				/*15*/
				I_measures.I152346[i][j] = 0.0;
				I_measures.I152436[i][j] = 0.0;
				I_measures.I152634[i][j] = 0.0;
				/*16*/
				I_measures.I162345[i][j] = 0.0;
				I_measures.I162435[i][j] = 0.0;
				I_measures.I162534[i][j] = 0.0;
			}
		}
	}	

	return 1;
}


int alloc_init_I_windows_ptr(struct I_windows_ptr *ptr_I_windows_ptr, int order, int full_b, int chainLen, int windowLength, int stepSize, int nrOfWindows, int alreadyAllocatedTo_nrOfWindows){


	struct I_windows_ptr I_windows;

	int N = nrOfWindows + 1; 
	int i;

	if(alreadyAllocatedTo_nrOfWindows == 0){ //indicates that no allocation was done yet

		I_windows.fileName = (char *) malloc (sizeof(char)*fileNameCharSize);
		strcpy(I_windows.fileName, "NN");

		I_windows.structureName = (char *) malloc (sizeof(char)*structureNameCharSize);
		strcpy(I_windows.structureName, "NN");

		I_windows.classId = (char *) malloc (sizeof(char)*classIdCharSize);
		strcpy(I_windows.classId, "NN");

		I_windows.chainId = (char *) malloc (chainIdCharSize*sizeof(char));
		strcpy(I_windows.chainId, "NN");
		
		I_windows.chainNr = 0;
		I_windows.chainLen = chainLen;
		I_windows.order = order;
		I_windows.windowLgth = windowLength;
		I_windows.stepSize = stepSize;
		I_windows.nrOfWindows = 0; /*this will store the true nr of windows for which invariants are kept*/

		I_windows.window = (struct window *) malloc((N+1)*sizeof(struct window));

		/*ini segment indices to zero everywhere:*/

		for (i=0; i <= N; i++){

			I_windows.window[i].windowNr = 0; 

			I_windows.window[i].segIndices = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

			I_windows.window[i].segIndices[0] = 0;
			I_windows.window[i].segIndices[1] = 0;

		}

		printf("nrOfWindows allocated for:%d\n", N);
		//getchar();

	}
	else if(nrOfWindows > alreadyAllocatedTo_nrOfWindows){ //indicates that reallocation must be done

		I_windows = *ptr_I_windows_ptr;

		strcpy(I_windows.fileName, "NN");
		strcpy(I_windows.structureName, "NN");
		strcpy(I_windows.classId, "NN");
		strcpy(I_windows.chainId, "NN");
		
		I_windows.chainNr = 0;
		I_windows.chainLen = chainLen;
		I_windows.order = order;
		I_windows.windowLgth = windowLength;
		I_windows.stepSize = stepSize;
		I_windows.nrOfWindows = 0; /*this will store the true nr of windows for which invariants are kept*/

		printf("nrOfWindows allocated for:%d\n", N);
		//getchar();

		/*reallocate*/
		I_windows.window = realloc(I_windows.window, (N+1)*sizeof(struct window));

		//printf("Langmodigt venter bysvalen stadig stædigt længes\n");

		/*ini segment indices to zero everywhere:*/

		for (i=0; i <= N; i++){

			//if(i< alreadyAllocatedTo_nrOfWindows){printf("Sommeren mild oedsles ... %d %d\n", i, I_windows.window[i].windowNr);}

			I_windows.window[i].windowNr = 0; 

			if(i >= alreadyAllocatedTo_nrOfWindows){
				I_windows.window[i].segIndices = (int *) malloc(2*sizeof(int));
			}
			else if(i < alreadyAllocatedTo_nrOfWindows){
				I_windows.window[i].segIndices = realloc(I_windows.window[i].segIndices, 2*sizeof(int));
			}

			//I_windows.window[i].segIndices = realloc(I_windows.window[i].segIndices, 2*sizeof(int));

			//printf("... kyssene nemmere pluds'ligt\n");

			I_windows.window[i].segIndices[0] = 0;
			I_windows.window[i].segIndices[1] = 0;

			//printf("I_windows.window[i].segIndices[0] %d, I_windows.window[i].segIndices[1] %d\n", I_windows.window[i].segIndices[0], I_windows.window[i].segIndices[1]);
			//getchar();

		}

	}
	else { //just reinit

		I_windows = *ptr_I_windows_ptr;

		strcpy(I_windows.fileName, "NN");
		strcpy(I_windows.structureName, "NN");
		strcpy(I_windows.classId, "NN");
		strcpy(I_windows.chainId, "NN");
		
		I_windows.chainNr = 0;
		I_windows.chainLen = chainLen;
		I_windows.order = order;
		I_windows.windowLgth = windowLength;
		I_windows.stepSize = stepSize;
		I_windows.nrOfWindows = 0; /*this will store the true nr of windows for which invariants are kept*/

		for (i=0; i <= N; i++){

			I_windows.window[i].windowNr = 0; 

			I_windows.window[i].segIndices[0] = 0;
			I_windows.window[i].segIndices[1] = 0;

			//printf("I_windows.window[i].segIndices[0] %d, I_windows.window[i].segIndices[1] %d\n", I_windows.window[i].segIndices[0], I_windows.window[i].segIndices[1]);
			//getchar();

		}


	}


	if (order == 1){

		/*allocate memory for double arrays of w and I12 values (order <= 1 measures)*/
		if(alreadyAllocatedTo_nrOfWindows == 0){ //indicates that no allocation was done yet

			I_windows.I12 = (double *) malloc ((N+1)*sizeof(double));
				
			I_windows.Ia12 = (double *) malloc ((N+1)*sizeof(double));
				
			for (i=0; i <= N; i++){
				
				/* Initialize to zero everywhere: */
				I_windows.I12[i] = 0.0;
				I_windows.Ia12[i] = 0.0;

			}

		}
		else if(nrOfWindows > alreadyAllocatedTo_nrOfWindows){ //indicates that reallocation must be done
	
			I_windows.I12 = realloc (I_windows.I12, (N+1)*sizeof(double));
				
			I_windows.Ia12 = realloc (I_windows.Ia12, (N+1)*sizeof(double));
				
			for (i=0; i <= N; i++){
				
				/* Initialize to zero everywhere: */
				I_windows.I12[i] = 0.0;
				I_windows.Ia12[i] = 0.0;

			}


		}
		else{ //just reinit
	

			for (i=0; i <= N; i++){
				
				/* Initialize to zero everywhere: */
				I_windows.I12[i] = 0.0;
				I_windows.Ia12[i] = 0.0;

			}

		}

	}


	if (order == 2){

		/*allocate memory for double arrays of for measures of order <= 2: */
		if(alreadyAllocatedTo_nrOfWindows == 0){ //indicates that no allocation was done yet

			I_windows.I12 = (double *) malloc ((N+1)*sizeof(double));

			I_windows.I1234 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.I1324 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.I1423 = (double *) malloc ((N+1)*sizeof(double));

			/*absolute value versions*/
			I_windows.Ia12 = (double *) malloc ((N+1)*sizeof(double));
				
			I_windows.Ia1234 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.I12a34 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.Ia12a34 = (double *) malloc ((N+1)*sizeof(double));
					
			I_windows.Ia1324 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.I13a24 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.Ia13a24 = (double *) malloc ((N+1)*sizeof(double));
					
			I_windows.Ia1423 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.I14a23 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.Ia14a23 = (double *) malloc ((N+1)*sizeof(double));

			/*for "full" order 2*/
			if (full_b == 1){
				I_windows.I1234_full = (double *) malloc ((N+1)*sizeof(double));
				I_windows.I1324_full = (double *) malloc ((N+1)*sizeof(double));

				/*abs value versions*/
				I_windows.Ia12a34_full = (double *) malloc ((N+1)*sizeof(double));
				I_windows.Ia13a24_full = (double *) malloc ((N+1)*sizeof(double));

				I_windows.Ia1234_full = (double *) malloc ((N+1)*sizeof(double));
				I_windows.Ia1324_full = (double *) malloc ((N+1)*sizeof(double));

				I_windows.I12a34_full = (double *) malloc ((N+1)*sizeof(double));
				I_windows.I13a24_full = (double *) malloc ((N+1)*sizeof(double));

			}

			/* Initialize to zero everywhere: */
			for (i=0; i <= N; i++){
					
				I_windows.I12[i] = 0.0;
				I_windows.Ia12[i] = 0.0;

				I_windows.I1234[i] = 0.0;
				I_windows.I1324[i] = 0.0;
				I_windows.I1423[i] = 0.0;

				I_windows.Ia1234[i] = 0.0;
				I_windows.I12a34[i] = 0.0;
				I_windows.Ia12a34[i] = 0.0;

				I_windows.Ia1324[i] = 0.0;
				I_windows.I13a24[i] = 0.0;
				I_windows.Ia13a24[i] = 0.0;
							
				I_windows.Ia1423[i] = 0.0;
				I_windows.I14a23[i] = 0.0;
				I_windows.Ia14a23[i] = 0.0;

				if (full_b == 1){
					I_windows.I1234_full[i] = 0.0;
					I_windows.I1324_full[i] = 0.0;

					/*abs value versions*/
					I_windows.Ia12a34_full[i] = 0.0;
					I_windows.Ia13a24_full[i] = 0.0;

					I_windows.Ia1234_full[i] = 0.0;
					I_windows.Ia1324_full[i] = 0.0;

					I_windows.I12a34_full[i] = 0.0;
					I_windows.I13a24_full[i] = 0.0;
				}

			}
		
		}
		else if(nrOfWindows > alreadyAllocatedTo_nrOfWindows){ //indicates that reallocation must be done

			I_windows.I12 = realloc(I_windows.I12 ,(N+1)*sizeof(double)); 

			I_windows.I1234 = realloc(I_windows.I1234 ,(N+1)*sizeof(double)); 
			I_windows.I1324 = realloc(I_windows.I1324 ,(N+1)*sizeof(double)); 
			I_windows.I1423 = realloc(I_windows.I1423 ,(N+1)*sizeof(double)); 

			/*absolute value versions*/
			I_windows.Ia12 = realloc(I_windows.Ia12 ,(N+1)*sizeof(double)); 
				
			I_windows.Ia1234 = realloc(I_windows.Ia1234 ,(N+1)*sizeof(double)); 
			I_windows.I12a34 = realloc(I_windows.I12a34 ,(N+1)*sizeof(double)); 
			I_windows.Ia12a34 = realloc(I_windows.Ia12a34 ,(N+1)*sizeof(double)); 
					
			I_windows.Ia1324 = realloc(I_windows.Ia1324 ,(N+1)*sizeof(double)); 
			I_windows.I13a24 = realloc(I_windows.I13a24 ,(N+1)*sizeof(double)); 
			I_windows.Ia13a24 = realloc(I_windows.Ia13a24 ,(N+1)*sizeof(double)); 
					
			I_windows.Ia1423 = realloc(I_windows.Ia1423 ,(N+1)*sizeof(double)); 
			I_windows.I14a23 = realloc(I_windows.I14a23 ,(N+1)*sizeof(double)); 
			I_windows.Ia14a23 = realloc(I_windows.Ia14a23 ,(N+1)*sizeof(double)); 

			/*for "full" order 2*/
			if (full_b == 1){
				I_windows.I1234_full = realloc(I_windows.I1234_full ,(N+1)*sizeof(double)); 
				I_windows.I1324_full = realloc(I_windows.I1324_full ,(N+1)*sizeof(double)); 

				/*abs value versions*/
				I_windows.Ia12a34_full = realloc(I_windows.Ia12a34_full ,(N+1)*sizeof(double)); 
				I_windows.Ia13a24_full = realloc(I_windows.Ia13a24_full ,(N+1)*sizeof(double)); 

				I_windows.Ia1234_full = realloc(I_windows.Ia1234_full ,(N+1)*sizeof(double)); 
				I_windows.Ia1324_full = realloc(I_windows.Ia1324_full ,(N+1)*sizeof(double)); 

				I_windows.I12a34_full = realloc(I_windows.I12a34_full ,(N+1)*sizeof(double)); 
				I_windows.I13a24_full = realloc(I_windows.I13a24_full ,(N+1)*sizeof(double)); 

			}


			/* Initialize to zero everywhere: */
			for (i=0; i <= N; i++){
					
				I_windows.I12[i] = 0.0;
				I_windows.Ia12[i] = 0.0;

				I_windows.I1234[i] = 0.0;
				I_windows.I1324[i] = 0.0;
				I_windows.I1423[i] = 0.0;

				I_windows.Ia1234[i] = 0.0;
				I_windows.I12a34[i] = 0.0;
				I_windows.Ia12a34[i] = 0.0;

				I_windows.Ia1324[i] = 0.0;
				I_windows.I13a24[i] = 0.0;
				I_windows.Ia13a24[i] = 0.0;
							
				I_windows.Ia1423[i] = 0.0;
				I_windows.I14a23[i] = 0.0;
				I_windows.Ia14a23[i] = 0.0;

				if (full_b == 1){
					I_windows.I1234_full[i] = 0.0;
					I_windows.I1324_full[i] = 0.0;

					/*abs value versions*/
					I_windows.Ia12a34_full[i] = 0.0;
					I_windows.Ia13a24_full[i] = 0.0;

					I_windows.Ia1234_full[i] = 0.0;
					I_windows.Ia1324_full[i] = 0.0;

					I_windows.I12a34_full[i] = 0.0;
					I_windows.I13a24_full[i] = 0.0;
				}

			}

		}
		else{ //just reinit

			/* Initialize to zero everywhere: */
			for (i=0; i <= N; i++){
					
				I_windows.I12[i] = 0.0;
				I_windows.Ia12[i] = 0.0;

				I_windows.I1234[i] = 0.0;
				I_windows.I1324[i] = 0.0;
				I_windows.I1423[i] = 0.0;

				I_windows.Ia1234[i] = 0.0;
				I_windows.I12a34[i] = 0.0;
				I_windows.Ia12a34[i] = 0.0;

				I_windows.Ia1324[i] = 0.0;
				I_windows.I13a24[i] = 0.0;
				I_windows.Ia13a24[i] = 0.0;
							
				I_windows.Ia1423[i] = 0.0;
				I_windows.I14a23[i] = 0.0;
				I_windows.Ia14a23[i] = 0.0;

				if (full_b == 1){
					I_windows.I1234_full[i] = 0.0;
					I_windows.I1324_full[i] = 0.0;

					/*abs value versions*/
					I_windows.Ia12a34_full[i] = 0.0;
					I_windows.Ia13a24_full[i] = 0.0;

					I_windows.Ia1234_full[i] = 0.0;
					I_windows.Ia1324_full[i] = 0.0;

					I_windows.I12a34_full[i] = 0.0;
					I_windows.I13a24_full[i] = 0.0;
				}

			}

		}

		
	}


	if (order == 3){


		/*allocate memory for double arrays of for measures of order <= 2: */
		if(alreadyAllocatedTo_nrOfWindows == 0){ //indicates that no allocation was done yet

			I_windows.I12 = (double *) malloc ((N+1)*sizeof(double));

			I_windows.I1234 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.I1324 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.I1423 = (double *) malloc ((N+1)*sizeof(double));

			/*absolute value versions*/
			I_windows.Ia12 = (double *) malloc ((N+1)*sizeof(double));
				
			I_windows.Ia1234 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.I12a34 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.Ia12a34 = (double *) malloc ((N+1)*sizeof(double));
					
			I_windows.Ia1324 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.I13a24 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.Ia13a24 = (double *) malloc ((N+1)*sizeof(double));
					
			I_windows.Ia1423 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.I14a23 = (double *) malloc ((N+1)*sizeof(double));
			I_windows.Ia14a23 = (double *) malloc ((N+1)*sizeof(double));

			/*for "full" order 2*/
			I_windows.I1234_full = (double *) malloc ((N+1)*sizeof(double));
			I_windows.I1324_full = (double *) malloc ((N+1)*sizeof(double));

			/*abs value versions*/
			I_windows.Ia12a34_full = (double *) malloc ((N+1)*sizeof(double));
			I_windows.Ia13a24_full = (double *) malloc ((N+1)*sizeof(double));

			I_windows.Ia1234_full = (double *) malloc ((N+1)*sizeof(double));
			I_windows.Ia1324_full = (double *) malloc ((N+1)*sizeof(double));

			I_windows.I12a34_full = (double *) malloc ((N+1)*sizeof(double));
			I_windows.I13a24_full = (double *) malloc ((N+1)*sizeof(double));

			/*order 3*/
			/*12*/
			I_windows.I123456 = (double *) malloc ((N+1)*sizeof(double)); 
			I_windows.I123645 = (double *) malloc ((N+1)*sizeof(double)); 
			I_windows.I123546 = (double *) malloc ((N+1)*sizeof(double)); 
			/*13*/
			I_windows.I132456 = (double *) malloc ((N+1)*sizeof(double)); 
			I_windows.I132546 = (double *) malloc ((N+1)*sizeof(double)); 
			I_windows.I132645 = (double *) malloc ((N+1)*sizeof(double)); 
			/*14*/
			I_windows.I142356 = (double *) malloc ((N+1)*sizeof(double)); 
			I_windows.I142536 = (double *) malloc ((N+1)*sizeof(double)); 
			I_windows.I142635 = (double *) malloc ((N+1)*sizeof(double)); 
			/*15*/
			I_windows.I152346 = (double *) malloc ((N+1)*sizeof(double)); 
			I_windows.I152436 = (double *) malloc ((N+1)*sizeof(double)); 
			I_windows.I152634 = (double *) malloc ((N+1)*sizeof(double)); 
			/*16*/
			I_windows.I162345 = (double *) malloc ((N+1)*sizeof(double)); 
			I_windows.I162435 = (double *) malloc ((N+1)*sizeof(double)); 
			I_windows.I162534 = (double *) malloc ((N+1)*sizeof(double)); 

			/* Initialize to zero everywhere: */
			for (i=0; i <= N; i++){
					
				I_windows.I12[i] = 0.0;
				I_windows.Ia12[i] = 0.0;

				I_windows.I1234[i] = 0.0;
				I_windows.I1324[i] = 0.0;
				I_windows.I1423[i] = 0.0;
					
				I_windows.Ia1234[i] = 0.0;
				I_windows.I12a34[i] = 0.0;
				I_windows.Ia12a34[i] = 0.0;

				I_windows.Ia1324[i] = 0.0;
				I_windows.I13a24[i] = 0.0;
				I_windows.Ia13a24[i] = 0.0;
							
				I_windows.Ia1423[i] = 0.0;
				I_windows.I14a23[i] = 0.0;
				I_windows.Ia14a23[i] = 0.0;

				I_windows.I1234_full[i] = 0.0;
				I_windows.I1324_full[i] = 0.0;

				/*abs value versions order 2 full*/
				I_windows.Ia12a34_full[i] = 0.0;
				I_windows.Ia13a24_full[i] = 0.0;

				I_windows.Ia1234_full[i] = 0.0;
				I_windows.Ia1324_full[i] = 0.0;

				I_windows.I12a34_full[i] = 0.0;
				I_windows.I13a24_full[i] = 0.0;

				//I_windows.I1324_full2[i] = 0.0;
				//I_windows.I1423_full0[i] = 0.0;
				//I_windows.I1423_full2[i] = 0.0;


				/*order 3*/
				/*12*/
				I_windows.I123456[i] = 0.0;
				I_windows.I123645[i] = 0.0;
				I_windows.I123546[i] = 0.0;
				/*13*/
				I_windows.I132456[i] = 0.0;
				I_windows.I132546[i] = 0.0;
				I_windows.I132645[i] = 0.0;
				/*14*/
				I_windows.I142356[i] = 0.0;
				I_windows.I142536[i] = 0.0;
				I_windows.I142635[i] = 0.0;
				/*15*/
				I_windows.I152346[i] = 0.0;
				I_windows.I152436[i] = 0.0;
				I_windows.I152634[i] = 0.0;
				/*16*/
				I_windows.I162345[i] = 0.0;
				I_windows.I162435[i] = 0.0;
				I_windows.I162534[i] = 0.0;


			}

		
		}
		else if(nrOfWindows > alreadyAllocatedTo_nrOfWindows){ //indicates that reallocation must be done
		
			/*order 1*/
			I_windows.I12 = realloc(I_windows.I12 ,(N+1)*sizeof(double)); 

			/*order 2*/
			I_windows.I1234 = realloc(I_windows.I1234 ,(N+1)*sizeof(double)); 
			I_windows.I1324 = realloc(I_windows.I1324 ,(N+1)*sizeof(double)); 
			I_windows.I1423 = realloc(I_windows.I1423 ,(N+1)*sizeof(double)); 

			I_windows.I1234_full = realloc(I_windows.I1234_full ,(N+1)*sizeof(double)); 
			I_windows.I1324_full = realloc(I_windows.I1324_full ,(N+1)*sizeof(double)); 

			/*abs value versions*/
			I_windows.Ia12a34_full = realloc(I_windows.Ia12a34_full ,(N+1)*sizeof(double)); 
			I_windows.Ia13a24_full = realloc(I_windows.Ia13a24_full ,(N+1)*sizeof(double)); 

			I_windows.Ia1234_full = realloc(I_windows.Ia1234_full ,(N+1)*sizeof(double)); 
			I_windows.Ia1324_full = realloc(I_windows.Ia1324_full ,(N+1)*sizeof(double)); 

			I_windows.I12a34_full = realloc(I_windows.I12a34_full ,(N+1)*sizeof(double)); 
			I_windows.I13a24_full = realloc(I_windows.I13a24_full ,(N+1)*sizeof(double)); 

			/*I_windows.I1324_full2 = realloc(I_windows. ,(N+1)*sizeof(double)); 
			I_windows.I1423_full0 = realloc(I_windows. ,(N+1)*sizeof(double)); 
			I_windows.I1423_full2 = realloc(I_windows. ,(N+1)*sizeof(double)); */

			/*absolute value versions*/
			I_windows.Ia12 = realloc(I_windows.Ia12 ,(N+1)*sizeof(double)); 
				
			I_windows.Ia1234 = realloc(I_windows.Ia1234 ,(N+1)*sizeof(double)); 
			I_windows.I12a34 = realloc(I_windows.I12a34 ,(N+1)*sizeof(double)); 
			I_windows.Ia12a34 = realloc(I_windows.Ia12a34 ,(N+1)*sizeof(double)); 
					
			I_windows.Ia1324 = realloc(I_windows.Ia1324 ,(N+1)*sizeof(double)); 
			I_windows.I13a24 = realloc(I_windows.I13a24 ,(N+1)*sizeof(double)); 
			I_windows.Ia13a24 = realloc(I_windows.Ia13a24 ,(N+1)*sizeof(double)); 
					
			I_windows.Ia1423 = realloc(I_windows.Ia1423 ,(N+1)*sizeof(double)); 
			I_windows.I14a23 = realloc(I_windows.I14a23,(N+1)*sizeof(double)); 
			I_windows.Ia14a23 = realloc(I_windows.Ia14a23 ,(N+1)*sizeof(double)); 

			/*order 3*/
			/*12*/
			I_windows.I123456 = realloc(I_windows.I123456 ,(N+1)*sizeof(double)); 
			I_windows.I123645 = realloc(I_windows.I123645 ,(N+1)*sizeof(double)); 
			I_windows.I123546 = realloc(I_windows.I123546 ,(N+1)*sizeof(double)); 
			/*13*/
			I_windows.I132456 = realloc(I_windows.I132456 ,(N+1)*sizeof(double)); 
			I_windows.I132546 = realloc(I_windows.I132546 ,(N+1)*sizeof(double)); 
			I_windows.I132645 = realloc(I_windows.I132645 ,(N+1)*sizeof(double)); 
			/*14*/
			I_windows.I142356 = realloc(I_windows.I142356 ,(N+1)*sizeof(double)); 
			I_windows.I142536 = realloc(I_windows.I142536 ,(N+1)*sizeof(double)); 
			I_windows.I142635 = realloc(I_windows.I142635 ,(N+1)*sizeof(double)); 
			/*15*/
			I_windows.I152346 = realloc(I_windows.I152346 ,(N+1)*sizeof(double)); 
			I_windows.I152436 = realloc(I_windows.I152436 ,(N+1)*sizeof(double)); 
			I_windows.I152634 = realloc(I_windows.I152634 ,(N+1)*sizeof(double)); 
			/*16*/
			I_windows.I162345 = realloc(I_windows.I162345 ,(N+1)*sizeof(double)); 
			I_windows.I162435 = realloc(I_windows.I162435 ,(N+1)*sizeof(double)); 
			I_windows.I162534 = realloc(I_windows.I162534 ,(N+1)*sizeof(double)); 


			/* Initialize to zero everywhere: */
			for (i=0; i <= N; i++){
					
				I_windows.I12[i] = 0.0;
				I_windows.Ia12[i] = 0.0;

				I_windows.I1234[i] = 0.0;
				I_windows.I1324[i] = 0.0;
				I_windows.I1423[i] = 0.0;
					
				I_windows.Ia1234[i] = 0.0;
				I_windows.I12a34[i] = 0.0;
				I_windows.Ia12a34[i] = 0.0;

				I_windows.Ia1324[i] = 0.0;
				I_windows.I13a24[i] = 0.0;
				I_windows.Ia13a24[i] = 0.0;
							
				I_windows.Ia1423[i] = 0.0;
				I_windows.I14a23[i] = 0.0;
				I_windows.Ia14a23[i] = 0.0;

				I_windows.I1234_full[i] = 0.0;
				I_windows.I1324_full[i] = 0.0;

				/*abs value versions order 2 full*/
				I_windows.Ia12a34_full[i] = 0.0;
				I_windows.Ia13a24_full[i] = 0.0;

				I_windows.Ia1234_full[i] = 0.0;
				I_windows.Ia1324_full[i] = 0.0;

				I_windows.I12a34_full[i] = 0.0;
				I_windows.I13a24_full[i] = 0.0;

				//I_windows.I1324_full2[i] = 0.0;
				//I_windows.I1423_full0[i] = 0.0;
				//I_windows.I1423_full2[i] = 0.0;


				/*order 3*/
				/*12*/
				I_windows.I123456[i] = 0.0;
				I_windows.I123645[i] = 0.0;
				I_windows.I123546[i] = 0.0;
				/*13*/
				I_windows.I132456[i] = 0.0;
				I_windows.I132546[i] = 0.0;
				I_windows.I132645[i] = 0.0;
				/*14*/
				I_windows.I142356[i] = 0.0;
				I_windows.I142536[i] = 0.0;
				I_windows.I142635[i] = 0.0;
				/*15*/
				I_windows.I152346[i] = 0.0;
				I_windows.I152436[i] = 0.0;
				I_windows.I152634[i] = 0.0;
				/*16*/
				I_windows.I162345[i] = 0.0;
				I_windows.I162435[i] = 0.0;
				I_windows.I162534[i] = 0.0;


			}
		}
		else { //just reinit
	
			/* Initialize to zero everywhere: */
			for (i=0; i <= N; i++){
					
				I_windows.I12[i] = 0.0;
				I_windows.Ia12[i] = 0.0;

				I_windows.I1234[i] = 0.0;
				I_windows.I1324[i] = 0.0;
				I_windows.I1423[i] = 0.0;
					
				I_windows.Ia1234[i] = 0.0;
				I_windows.I12a34[i] = 0.0;
				I_windows.Ia12a34[i] = 0.0;

				I_windows.Ia1324[i] = 0.0;
				I_windows.I13a24[i] = 0.0;
				I_windows.Ia13a24[i] = 0.0;
							
				I_windows.Ia1423[i] = 0.0;
				I_windows.I14a23[i] = 0.0;
				I_windows.Ia14a23[i] = 0.0;

				I_windows.I1234_full[i] = 0.0;
				I_windows.I1324_full[i] = 0.0;

				/*abs value versions order 2 full*/
				I_windows.Ia12a34_full[i] = 0.0;
				I_windows.Ia13a24_full[i] = 0.0;

				I_windows.Ia1234_full[i] = 0.0;
				I_windows.Ia1324_full[i] = 0.0;

				I_windows.I12a34_full[i] = 0.0;
				I_windows.I13a24_full[i] = 0.0;

				//I_windows.I1324_full2[i] = 0.0;
				//I_windows.I1423_full0[i] = 0.0;
				//I_windows.I1423_full2[i] = 0.0;


				/*order 3*/
				/*12*/
				I_windows.I123456[i] = 0.0;
				I_windows.I123645[i] = 0.0;
				I_windows.I123546[i] = 0.0;
				/*13*/
				I_windows.I132456[i] = 0.0;
				I_windows.I132546[i] = 0.0;
				I_windows.I132645[i] = 0.0;
				/*14*/
				I_windows.I142356[i] = 0.0;
				I_windows.I142536[i] = 0.0;
				I_windows.I142635[i] = 0.0;
				/*15*/
				I_windows.I152346[i] = 0.0;
				I_windows.I152436[i] = 0.0;
				I_windows.I152634[i] = 0.0;
				/*16*/
				I_windows.I162345[i] = 0.0;
				I_windows.I162435[i] = 0.0;
				I_windows.I162534[i] = 0.0;


			}
		}

	}

	if(alreadyAllocatedTo_nrOfWindows == 0 || nrOfWindows > alreadyAllocatedTo_nrOfWindows){
	
		printf("nrOfWindows allocated for: %d\n", i);

	}
	else{

		printf("nrOfWindows allocated to previously (no new alloc done): %d\n", alreadyAllocatedTo_nrOfWindows);

	}

	*ptr_I_windows_ptr = I_windows;

	//getchar();


	return N + 1;

}

int alloc_init_I_windowPairs_ptr(struct I_windowPairs_ptr *ptr_I_windowPairs_ptr, int order, int full_b, int chainLen, int windowLength, int stepSize, int nrOfWindows, int alreadyAllocatedTo_nrOfWindows){

	int returnVal = 0;

	struct I_windowPairs_ptr I_windowPairs;

	int N = nrOfWindows + 1; /*should suffice: the nrOfWindows will be set for serving the longest chain in the DB-set*/
	int i = 0,j = 0;

	int cntPairs = 0;

	//printf("nrOfWindowPairs allocated for:%d\n", (N+1)*(N+1));

	if(alreadyAllocatedTo_nrOfWindows == 0){ //indicates that no allocation was done yet

		/*allocate and init*/
		I_windowPairs.fileName = (char *) malloc (sizeof(char)*fileNameCharSize);
		strcpy(I_windowPairs.fileName, "NN");
		
		I_windowPairs.structureName = (char *) malloc(sizeof(char)*structureNameCharSize);
		strcpy(I_windowPairs.structureName, "NN");

		I_windowPairs.classId = (char *) malloc (sizeof(char)*classIdCharSize);
		strcpy(I_windowPairs.classId, "NN");

		I_windowPairs.chainId = (char *) malloc (sizeof(char)*chainIdCharSize);
		strcpy(I_windowPairs.chainId, "NN");

		I_windowPairs.chainNr = 0;
		I_windowPairs.chainLen = chainLen;
		I_windowPairs.order = order;
		I_windowPairs.windowLgth = windowLength;
		I_windowPairs.stepSize = stepSize;
		I_windowPairs.nrOfWindows = 0; /*this will store the true nr of windows involved (giving rise to nrOfWindowPairs) */
		I_windowPairs.nrOfWindowPairs = 0; /*this will store the true nr of window pairs for which invariants are kept*/

		I_windowPairs.windowPair = (struct windowPair **) malloc((N+1)*(N+1)*sizeof(struct windowPair));

		/*init segment indices to zero everywhere:*/
		for (i=0; i <= N; i++){

			I_windowPairs.windowPair[i] = (struct windowPair *) malloc((N+1)*sizeof(struct windowPair));

			for (j=0; j <= N; j++){

				I_windowPairs.windowPair[i][j].windowNr_1 = 0; 
				I_windowPairs.windowPair[i][j].windowNr_2 = 0;

				I_windowPairs.windowPair[i][j].segIndices_1 = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 
				I_windowPairs.windowPair[i][j].segIndices_2 = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

				I_windowPairs.windowPair[i][j].segIndices_1[0] = 0;
				I_windowPairs.windowPair[i][j].segIndices_1[1] = 0;

				I_windowPairs.windowPair[i][j].segIndices_2[0] = 0;
				I_windowPairs.windowPair[i][j].segIndices_2[1] = 0;

			}
		}

	}
	else if(nrOfWindows > alreadyAllocatedTo_nrOfWindows){ //indicates that reallocation must be done

		I_windowPairs = *ptr_I_windowPairs_ptr;
		
		/*(re)init*/
		strcpy(I_windowPairs.fileName, "NN");
		strcpy(I_windowPairs.structureName, "NN");
		strcpy(I_windowPairs.classId, "NN");
		strcpy(I_windowPairs.chainId, "NN");

		I_windowPairs.chainNr = 0;
		I_windowPairs.chainLen = chainLen;
		I_windowPairs.order = order;
		I_windowPairs.windowLgth = windowLength;
		I_windowPairs.stepSize = stepSize;
		I_windowPairs.nrOfWindows = 0; /*this will store the true nr of windows involved (giving rise to nrOfWindowPairs) */
		I_windowPairs.nrOfWindowPairs = 0; /*this will store the true nr of window pairs for which invariants are kept*/

		/*reallocate*/
		I_windowPairs.windowPair = realloc(I_windowPairs.windowPair, (N+1)*(N+1)*sizeof(struct windowPair));

		/*init segment indices to zero everywhere:*/
		for (i=0; i <= N; i++){

			if(i < alreadyAllocatedTo_nrOfWindows){
				I_windowPairs.windowPair[i] = realloc(I_windowPairs.windowPair[i], (N+1)*sizeof(struct windowPair));
			} else if(i >= alreadyAllocatedTo_nrOfWindows){

				I_windowPairs.windowPair[i] = (struct windowPair *) malloc((N+1)*sizeof(struct windowPair));

			}	


			for (j=0; j <= N; j++){

				I_windowPairs.windowPair[i][j].windowNr_1 = 0; 
				I_windowPairs.windowPair[i][j].windowNr_2 = 0;

				if(i >= alreadyAllocatedTo_nrOfWindows || j >= alreadyAllocatedTo_nrOfWindows){
					I_windowPairs.windowPair[i][j].segIndices_1 = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 
					I_windowPairs.windowPair[i][j].segIndices_2 = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 
				}

				I_windowPairs.windowPair[i][j].segIndices_1[0] = 0;
				I_windowPairs.windowPair[i][j].segIndices_1[1] = 0;

				I_windowPairs.windowPair[i][j].segIndices_2[0] = 0;
				I_windowPairs.windowPair[i][j].segIndices_2[1] = 0;

			}
		}


	}
	else{ //just reinit

		I_windowPairs = *ptr_I_windowPairs_ptr;
		
		/*(re)init*/
		strcpy(I_windowPairs.fileName, "NN");
		strcpy(I_windowPairs.structureName, "NN");
		strcpy(I_windowPairs.classId, "NN");
		strcpy(I_windowPairs.chainId, "NN");

		I_windowPairs.chainNr = 0;
		I_windowPairs.chainLen = chainLen;
		I_windowPairs.order = order;
		I_windowPairs.windowLgth = windowLength;
		I_windowPairs.stepSize = stepSize;
		I_windowPairs.nrOfWindows = 0; /*this will store the true nr of windows involved (giving rise to nrOfWindowPairs) */
		I_windowPairs.nrOfWindowPairs = 0; /*this will store the true nr of window pairs for which invariants are kept*/

		/*init segment indices to zero everywhere:*/
		for (i=0; i <= N; i++){

			for (j=0; j <= N; j++){

				I_windowPairs.windowPair[i][j].windowNr_1 = 0; 
				I_windowPairs.windowPair[i][j].windowNr_2 = 0;

				I_windowPairs.windowPair[i][j].segIndices_1[0] = 0;
				I_windowPairs.windowPair[i][j].segIndices_1[1] = 0;

				I_windowPairs.windowPair[i][j].segIndices_2[0] = 0;
				I_windowPairs.windowPair[i][j].segIndices_2[1] = 0;

			}
		}


	}




	if (order == 1){

		/*Mutuals:*/
		/*allocate memory for double arrays of w and I12 values (order <= 1 measures)*/
		if(alreadyAllocatedTo_nrOfWindows == 0){ //indicates that no allocation was done yet
			
			I_windowPairs.I12 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
				
			I_windowPairs.Ia12 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
				
			for (i=0; i <= N; i++){
				
				I_windowPairs.I12[i] = (double *)malloc((N+1) * sizeof(double));
				I_windowPairs.Ia12[i] = (double *) malloc ((N+1)*sizeof(double));
				
				/* Initialize to zero everywhere: */
				for (j=i; j<= N; j++){
					I_windowPairs.I12[i][j] = 0.0;
					I_windowPairs.Ia12[i][j] = 0.0;

					cntPairs +=1;
				}

			}

		}
		else if(nrOfWindows > alreadyAllocatedTo_nrOfWindows){ //indicates that reallocation must be done

			I_windowPairs.I12 = realloc (I_windowPairs.I12, (N+1)*(N+1)*sizeof(double ));
				
			I_windowPairs.Ia12 = realloc (I_windowPairs.Ia12, (N+1)*(N+1)*sizeof(double ));
				
			for (i=0; i <= N; i++){
				
				if(i < alreadyAllocatedTo_nrOfWindows){

					I_windowPairs.I12[i] = realloc(I_windowPairs.I12[i], (N+1)*sizeof(double));
					I_windowPairs.Ia12[i] = realloc (I_windowPairs.Ia12[i], (N+1)*sizeof(double));

				} else if(i >= alreadyAllocatedTo_nrOfWindows){

					I_windowPairs.I12[i] = (double *) malloc((N+1)*sizeof(double));
					I_windowPairs.Ia12[i] =(double *) malloc((N+1)*sizeof(double));

				}
				
				/* Initialize to zero everywhere: */
				for (j=i; j<= N; j++){

					I_windowPairs.I12[i][j] = 0.0;
					I_windowPairs.Ia12[i][j] = 0.0;

					cntPairs +=1;
				}

			}

		}
		else { //just reinit
				
			for (i=0; i <= N; i++){
				
				/* Initialize to zero everywhere: */
				for (j=i; j<= N; j++){

					I_windowPairs.I12[i][j] = 0.0;
					I_windowPairs.Ia12[i][j] = 0.0;

					cntPairs +=1;
				}

			}

		}

	}


	if (order == 2){

		/*Mutuals:*/
		/*allocate memory for double arrays of for measures of order <= 2: */
		if(alreadyAllocatedTo_nrOfWindows == 0){ //indicates that no allocation was done yet	

			I_windowPairs.I12 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));

			I_windowPairs.I1234 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
			I_windowPairs.I1324 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
			I_windowPairs.I1423 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));

			/*absolute value versions*/
			I_windowPairs.Ia12 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
				
			I_windowPairs.Ia1234 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
			I_windowPairs.I12a34 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
			I_windowPairs.Ia12a34 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
					
			I_windowPairs.Ia1324 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
			I_windowPairs.I13a24 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
			I_windowPairs.Ia13a24 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
					
			I_windowPairs.Ia1423 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
			I_windowPairs.I14a23 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
			I_windowPairs.Ia14a23 = (double **) malloc ((N+1)*(N+1)*sizeof(double ));

			/*for "full" order 2*/
			if (full_b == 1){
				I_windowPairs.I1234_full = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
				I_windowPairs.I1324_full = (double **) malloc ((N+1)*(N+1)*sizeof(double ));

				/*abs value versions*/
				I_windowPairs.Ia12a34_full = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
				I_windowPairs.Ia13a24_full = (double **) malloc ((N+1)*(N+1)*sizeof(double ));

				I_windowPairs.Ia1234_full = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
				I_windowPairs.Ia1324_full = (double **) malloc ((N+1)*(N+1)*sizeof(double ));

				I_windowPairs.I12a34_full = (double **) malloc ((N+1)*(N+1)*sizeof(double ));
				I_windowPairs.I13a24_full = (double **) malloc ((N+1)*(N+1)*sizeof(double ));

			}

			/* Initialize to zero everywhere: */
			for (i=0; i <= N; i++){

				I_windowPairs.I12[i] = (double *)malloc((N+1)*sizeof(double));
				I_windowPairs.Ia12[i] = (double *) malloc((N+1)*sizeof(double));

				I_windowPairs.I1234[i] = (double *)malloc((N+1)*sizeof(double));
				I_windowPairs.I1324[i] = (double *) malloc((N+1)*sizeof(double));
				I_windowPairs.I1423[i] = (double *)malloc((N+1)*sizeof(double));

				I_windowPairs.Ia1234[i] = (double *)malloc((N+1)*sizeof(double));
				I_windowPairs.I12a34[i] = (double *)malloc((N+1)*sizeof(double));
				I_windowPairs.Ia12a34[i] = (double *)malloc((N+1)*sizeof(double));

				I_windowPairs.Ia1324[i] = (double *) malloc((N+1)*sizeof(double));
				I_windowPairs.I13a24[i] = (double *) malloc((N+1)*sizeof(double));
				I_windowPairs.Ia13a24[i] = (double *) malloc((N+1)*sizeof(double));

				I_windowPairs.Ia1423[i] = (double *)malloc((N+1)*sizeof(double));
				I_windowPairs.I14a23[i] = (double *)malloc((N+1)*sizeof(double));
				I_windowPairs.Ia14a23[i] = (double *)malloc((N+1)*sizeof(double));

				/*for "full" order 2*/
				if (full_b == 1){
					I_windowPairs.I1234_full[i] = (double *) malloc((N+1)*sizeof(double));
					I_windowPairs.I1324_full[i] = (double *) malloc((N+1)*sizeof(double));

					/*abs value versions*/
					I_windowPairs.Ia12a34_full[i] = (double *) malloc((N+1)*sizeof(double));
					I_windowPairs.Ia13a24_full[i] = (double *) malloc((N+1)*sizeof(double));

					I_windowPairs.Ia1234_full[i] = (double *) malloc((N+1)*sizeof(double));
					I_windowPairs.Ia1324_full[i] = (double *) malloc((N+1)*sizeof(double));

					I_windowPairs.I12a34_full[i] = (double *) malloc((N+1)*sizeof(double));
					I_windowPairs.I13a24_full[i] = (double *) malloc((N+1)*sizeof(double));
				}

				for (j=0; j<= N; j++){
					
					I_windowPairs.I12[i][j] = 0.0;
					I_windowPairs.Ia12[i][j] = 0.0;

					I_windowPairs.I1234[i][j] = 0.0;
					I_windowPairs.I1324[i][j] = 0.0;
					I_windowPairs.I1423[i][j] = 0.0;

					I_windowPairs.Ia1234[i][j] = 0.0;
					I_windowPairs.I12a34[i][j] = 0.0;
					I_windowPairs.Ia12a34[i][j] = 0.0;

					I_windowPairs.Ia1324[i][j] = 0.0;
					I_windowPairs.I13a24[i][j] = 0.0;
					I_windowPairs.Ia13a24[i][j] = 0.0;
							
					I_windowPairs.Ia1423[i][j] = 0.0;
					I_windowPairs.I14a23[i][j] = 0.0;
					I_windowPairs.Ia14a23[i][j] = 0.0;

					if (full_b == 1){
						I_windowPairs.I1234_full[i][j] = 0.0;
						I_windowPairs.I1324_full[i][j] = 0.0;

						/*abs value versions*/
						I_windowPairs.Ia12a34_full[i][j] = 0.0;
						I_windowPairs.Ia13a24_full[i][j] = 0.0;

						I_windowPairs.Ia1234_full[i][j] = 0.0;
						I_windowPairs.Ia1324_full[i][j] = 0.0;

						I_windowPairs.I12a34_full[i][j] = 0.0;
						I_windowPairs.I13a24_full[i][j] = 0.0;
					}

					cntPairs +=1;

				}

			}

		}
		else if(nrOfWindows > alreadyAllocatedTo_nrOfWindows){ //indicates that reallocation must be done

			I_windowPairs.I12 = realloc (I_windowPairs.I12 ,(N+1)*(N+1)*sizeof(double ));

			I_windowPairs.I1234 = realloc (I_windowPairs.I1234 ,(N+1)*(N+1)*sizeof(double ));
			I_windowPairs.I1324 = realloc (I_windowPairs.I1324 ,(N+1)*(N+1)*sizeof(double ));
			I_windowPairs.I1423 = realloc (I_windowPairs.I1423 ,(N+1)*(N+1)*sizeof(double ));

			/*absolute value versions*/
			I_windowPairs.Ia12 = realloc (I_windowPairs.Ia12 ,(N+1)*(N+1)*sizeof(double ));
				
			I_windowPairs.Ia1234 = realloc (I_windowPairs.Ia1234 ,(N+1)*(N+1)*sizeof(double ));
			I_windowPairs.I12a34 = realloc (I_windowPairs.I12a34 ,(N+1)*(N+1)*sizeof(double ));
			I_windowPairs.Ia12a34 = realloc (I_windowPairs.Ia12a34,(N+1)*(N+1)*sizeof(double ));
					
			I_windowPairs.Ia1324 = realloc (I_windowPairs.Ia1324 ,(N+1)*(N+1)*sizeof(double ));
			I_windowPairs.I13a24 = realloc (I_windowPairs.I13a24 ,(N+1)*(N+1)*sizeof(double ));
			I_windowPairs.Ia13a24 = realloc (I_windowPairs.Ia13a24 ,(N+1)*(N+1)*sizeof(double ));
					
			I_windowPairs.Ia1423 = realloc (I_windowPairs.Ia1423 ,(N+1)*(N+1)*sizeof(double ));
			I_windowPairs.I14a23 = realloc (I_windowPairs.I14a23 ,(N+1)*(N+1)*sizeof(double ));
			I_windowPairs.Ia14a23 = realloc (I_windowPairs.Ia14a23 ,(N+1)*(N+1)*sizeof(double ));

			/*for "full" order 2*/
			if (full_b == 1){
				I_windowPairs.I1234_full = realloc (I_windowPairs.I1234_full ,(N+1)*(N+1)*sizeof(double ));
				I_windowPairs.I1324_full = realloc (I_windowPairs.I1324_full ,(N+1)*(N+1)*sizeof(double ));

				/*abs value versions*/
				I_windowPairs.Ia12a34_full = realloc (I_windowPairs.Ia12a34_full ,(N+1)*(N+1)*sizeof(double ));
				I_windowPairs.Ia13a24_full = realloc (I_windowPairs.Ia13a24_full ,(N+1)*(N+1)*sizeof(double ));

				I_windowPairs.Ia1234_full = realloc (I_windowPairs.Ia1234_full ,(N+1)*(N+1)*sizeof(double ));
				I_windowPairs.Ia1324_full = realloc (I_windowPairs.Ia1324_full ,(N+1)*(N+1)*sizeof(double ));

				I_windowPairs.I12a34_full = realloc (I_windowPairs.I12a34_full ,(N+1)*(N+1)*sizeof(double ));
				I_windowPairs.I13a24_full = realloc (I_windowPairs.I13a24_full ,(N+1)*(N+1)*sizeof(double ));

			}

			/* Initialize to zero everywhere: */
			for (i=0; i <= N; i++){

				if(i < alreadyAllocatedTo_nrOfWindows){

					I_windowPairs.I12[i] = realloc(I_windowPairs.I12[i] ,(N+1)*sizeof(double));
					I_windowPairs.Ia12[i] = realloc(I_windowPairs.Ia12[i] ,(N+1)*sizeof(double));

					I_windowPairs.I1234[i] = realloc(I_windowPairs.I1234[i] ,(N+1)*sizeof(double));
					I_windowPairs.I1324[i] = realloc(I_windowPairs.I1324[i] ,(N+1)*sizeof(double));
					I_windowPairs.I1423[i] = realloc(I_windowPairs.I1423[i] ,(N+1)*sizeof(double));

					I_windowPairs.Ia1234[i] = realloc(I_windowPairs.Ia1234[i] ,(N+1)*sizeof(double));
					I_windowPairs.I12a34[i] = realloc(I_windowPairs.I12a34[i] ,(N+1)*sizeof(double));
					I_windowPairs.Ia12a34[i] = realloc(I_windowPairs.Ia12a34[i] ,(N+1)*sizeof(double));

					I_windowPairs.Ia1324[i] = realloc(I_windowPairs.Ia1324[i] ,(N+1)*sizeof(double));
					I_windowPairs.I13a24[i] = realloc(I_windowPairs.I13a24[i] ,(N+1)*sizeof(double));
					I_windowPairs.Ia13a24[i] = realloc(I_windowPairs.Ia13a24[i] ,(N+1)*sizeof(double));

					I_windowPairs.Ia1423[i] = realloc(I_windowPairs.Ia1423[i] ,(N+1)*sizeof(double));
					I_windowPairs.I14a23[i] = realloc(I_windowPairs.I14a23[i] ,(N+1)*sizeof(double));
					I_windowPairs.Ia14a23[i] = realloc(I_windowPairs.Ia14a23[i] ,(N+1)*sizeof(double));

					/*for "full" order 2*/
					if (full_b == 1){
						I_windowPairs.I1234_full[i] = realloc(I_windowPairs.I1234_full[i] ,(N+1)*sizeof(double));
						I_windowPairs.I1324_full[i] = realloc(I_windowPairs.I1324_full[i] ,(N+1)*sizeof(double));

						/*abs value versions*/
						I_windowPairs.Ia12a34_full[i] = realloc(I_windowPairs.Ia12a34_full[i] ,(N+1)*sizeof(double));
						I_windowPairs.Ia13a24_full[i] = realloc(I_windowPairs.Ia13a24_full[i] ,(N+1)*sizeof(double));

						I_windowPairs.Ia1234_full[i] = realloc(I_windowPairs.Ia1234_full[i] ,(N+1)*sizeof(double));
						I_windowPairs.Ia1324_full[i] = realloc(I_windowPairs.Ia1324_full[i] ,(N+1)*sizeof(double));

						I_windowPairs.I12a34_full[i] = realloc(I_windowPairs.I12a34_full[i] ,(N+1)*sizeof(double));
						I_windowPairs.I13a24_full[i] = realloc(I_windowPairs.I13a24_full[i] ,(N+1)*sizeof(double));
					}

				}
				else if(i >= alreadyAllocatedTo_nrOfWindows){

					I_windowPairs.I12[i] = (double *) malloc((N+1)*sizeof(double));
					I_windowPairs.Ia12[i] = (double *) malloc((N+1)*sizeof(double));

					I_windowPairs.I1234[i] = (double *) malloc((N+1)*sizeof(double));
					I_windowPairs.I1324[i] = (double *) malloc((N+1)*sizeof(double));
					I_windowPairs.I1423[i] = (double *) malloc((N+1)*sizeof(double));

					I_windowPairs.Ia1234[i] = (double *) malloc((N+1)*sizeof(double));
					I_windowPairs.I12a34[i] = (double *) malloc((N+1)*sizeof(double));
					I_windowPairs.Ia12a34[i] = (double *) malloc((N+1)*sizeof(double));

					I_windowPairs.Ia1324[i] = (double *) malloc((N+1)*sizeof(double));
					I_windowPairs.I13a24[i] = (double *) malloc((N+1)*sizeof(double));
					I_windowPairs.Ia13a24[i] = (double *) malloc((N+1)*sizeof(double));

					I_windowPairs.Ia1423[i] = (double *) malloc((N+1)*sizeof(double));
					I_windowPairs.I14a23[i] = (double *) malloc((N+1)*sizeof(double));
					I_windowPairs.Ia14a23[i] = (double *) malloc((N+1)*sizeof(double));

					/*for "full" order 2*/
					if (full_b == 1){
						I_windowPairs.I1234_full[i] = (double *) malloc((N+1)*sizeof(double));
						I_windowPairs.I1324_full[i] = (double *) malloc((N+1)*sizeof(double));

						/*abs value versions*/
						I_windowPairs.Ia12a34_full[i] = (double *) malloc((N+1)*sizeof(double));
						I_windowPairs.Ia13a24_full[i] = (double *) malloc((N+1)*sizeof(double));

						I_windowPairs.Ia1234_full[i] = (double *) malloc((N+1)*sizeof(double));
						I_windowPairs.Ia1324_full[i] = (double *) malloc((N+1)*sizeof(double));

						I_windowPairs.I12a34_full[i] = (double *) malloc((N+1)*sizeof(double));
						I_windowPairs.I13a24_full[i] = (double *) malloc((N+1)*sizeof(double));
					}

				}


				for (j=0; j<= N; j++){
				

					I_windowPairs.I12[i][j] = 0.0;
					I_windowPairs.Ia12[i][j] = 0.0;

					I_windowPairs.I1234[i][j] = 0.0;
					I_windowPairs.I1324[i][j] = 0.0;
					I_windowPairs.I1423[i][j] = 0.0;

					I_windowPairs.Ia1234[i][j] = 0.0;
					I_windowPairs.I12a34[i][j] = 0.0;
					I_windowPairs.Ia12a34[i][j] = 0.0;

					I_windowPairs.Ia1324[i][j] = 0.0;
					I_windowPairs.I13a24[i][j] = 0.0;
					I_windowPairs.Ia13a24[i][j] = 0.0;
							
					I_windowPairs.Ia1423[i][j] = 0.0;
					I_windowPairs.I14a23[i][j] = 0.0;
					I_windowPairs.Ia14a23[i][j] = 0.0;

					if (full_b == 1){
						I_windowPairs.I1234_full[i][j] = 0.0;
						I_windowPairs.I1324_full[i][j] = 0.0;

						/*abs value versions*/
						I_windowPairs.Ia12a34_full[i][j] = 0.0;
						I_windowPairs.Ia13a24_full[i][j] = 0.0;

						I_windowPairs.Ia1234_full[i][j] = 0.0;
						I_windowPairs.Ia1324_full[i][j] = 0.0;

						I_windowPairs.I12a34_full[i][j] = 0.0;
						I_windowPairs.I13a24_full[i][j] = 0.0;
					}

					cntPairs +=1;

				}
				
			}

			//printf("(saaledes stadig ganske hemmelig) hertil osse ok .. %d\n", i);
			//getchar();
		
		}
		else { //just reinit 

			/* Initialize to zero everywhere: */
			for (i=0; i <= N; i++){

				for (j=0; j<= N; j++){
				

					I_windowPairs.I12[i][j] = 0.0;
					I_windowPairs.Ia12[i][j] = 0.0;

					I_windowPairs.I1234[i][j] = 0.0;
					I_windowPairs.I1324[i][j] = 0.0;
					I_windowPairs.I1423[i][j] = 0.0;

					I_windowPairs.Ia1234[i][j] = 0.0;
					I_windowPairs.I12a34[i][j] = 0.0;
					I_windowPairs.Ia12a34[i][j] = 0.0;

					I_windowPairs.Ia1324[i][j] = 0.0;
					I_windowPairs.I13a24[i][j] = 0.0;
					I_windowPairs.Ia13a24[i][j] = 0.0;
							
					I_windowPairs.Ia1423[i][j] = 0.0;
					I_windowPairs.I14a23[i][j] = 0.0;
					I_windowPairs.Ia14a23[i][j] = 0.0;

					if (full_b == 1){
						I_windowPairs.I1234_full[i][j] = 0.0;
						I_windowPairs.I1324_full[i][j] = 0.0;

						/*abs value versions*/
						I_windowPairs.Ia12a34_full[i][j] = 0.0;
						I_windowPairs.Ia13a24_full[i][j] = 0.0;

						I_windowPairs.Ia1234_full[i][j] = 0.0;
						I_windowPairs.Ia1324_full[i][j] = 0.0;

						I_windowPairs.I12a34_full[i][j] = 0.0;
						I_windowPairs.I13a24_full[i][j] = 0.0;
					}

					cntPairs +=1;

				}
				
			}


		}


	}


	if(alreadyAllocatedTo_nrOfWindows == 0 || nrOfWindows > alreadyAllocatedTo_nrOfWindows){

		printf("nrOfWindowPairs allocated to: %d whereof initiated (simplex/upper triangular matrix): %d\n", (N+1)*(N+1) , cntPairs);

	}
	else{

		printf("nrOfWindowPairs allocated to previously (no new alloc done): %d whereof initiated (simplex/upper triangular matrix): %d\n", returnVal, cntPairs);

	}

	*ptr_I_windowPairs_ptr = I_windowPairs; /*necc??*/

	returnVal = (N+1)*(N+1);


	return returnVal;

}



/*Fcts for the rarity method:*/ 
int alloc_init_ptr_Irarity_windows(struct Irarity_windows **ptr_ptr_Irarity_windows, int order, int full_b, int chainLen, int windowLength, int stepSize, int nrOfWindows, int alreadyAllocatedTo_nrOfWindows){


	struct Irarity_windows *ptr_Irarity_windows;

	int N = nrOfWindows + 1; 
	int i;

	if(alreadyAllocatedTo_nrOfWindows == 0){ //indicates that no allocation was done yet

		ptr_Irarity_windows = (struct Irarity_windows *) malloc((N+1)*sizeof(struct Irarity_windows)); 

		/*ini segment indices to zero everywhere:*/

		for (i=0; i <= N; i++){

			ptr_Irarity_windows[i].window.windowNr = 0; 

			ptr_Irarity_windows[i].window.segIndices = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

			ptr_Irarity_windows[i].window.segIndices[0] = 0;
			ptr_Irarity_windows[i].window.segIndices[1] = 0;

			/* Initialize to zero everywhere: */
			ptr_Irarity_windows[i].rarityScore = 0.0;

		}

		*ptr_ptr_Irarity_windows = ptr_Irarity_windows;
	}
	else if(nrOfWindows > alreadyAllocatedTo_nrOfWindows){ //indicates that reallocation must be done

		ptr_Irarity_windows = *ptr_ptr_Irarity_windows;

		ptr_Irarity_windows = realloc(ptr_Irarity_windows, (N+1)*sizeof(struct Irarity_windows)); 

		/*ini segment indices to zero everywhere:*/

		for (i=0; i <= N; i++){

			ptr_Irarity_windows[i].window.windowNr = 0; 

			ptr_Irarity_windows[i].window.segIndices = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

			ptr_Irarity_windows[i].window.segIndices[0] = 0;
			ptr_Irarity_windows[i].window.segIndices[1] = 0;

			/* Initialize to zero everywhere: */
			ptr_Irarity_windows[i].rarityScore = 0.0;

		}
	
		*ptr_ptr_Irarity_windows = ptr_Irarity_windows; //needed??

	} else{ //just re-init to zero everywhere:

		ptr_Irarity_windows = *ptr_ptr_Irarity_windows;

		for (i=0; i < N; i++){

			ptr_Irarity_windows[i].window.windowNr = 0; 

			ptr_Irarity_windows[i].window.segIndices[0] = 0;
			ptr_Irarity_windows[i].window.segIndices[1] = 0;

			ptr_Irarity_windows[i].rarityScore = 0.0;

		}

		*ptr_ptr_Irarity_windows = ptr_Irarity_windows;

	}


	return N + 1;

}


int init_ptr_Irarity_windows(struct Irarity_windows **ptr_ptr_Irarity_windows,  int nrOfWindows){


	struct Irarity_windows *ptr_Irarity_windows;

	int N = nrOfWindows + 1; 
	int i;

	//just re-init to zero everywhere:

	ptr_Irarity_windows = *ptr_ptr_Irarity_windows;

	for (i=0; i < N; i++){

		ptr_Irarity_windows[i].window.windowNr = 0; 

		ptr_Irarity_windows[i].window.segIndices[0] = 0;
		ptr_Irarity_windows[i].window.segIndices[1] = 0;

		ptr_Irarity_windows[i].rarityScore = 0.0;

	}

	*ptr_ptr_Irarity_windows = ptr_Irarity_windows;

	return N + 1;

}


int alloc_init_ptr_Irarity_windowPairs(struct Irarity_windowPairs **ptr_ptr_Irarity_windowPairs, int order, int full_b, int chainLen, int windowLength, int stepSize, int nrOfWindows, int alreadyAllocatedTo_nrOfWindows){


	struct Irarity_windowPairs *ptr_Irarity_windowPairs;

	int N = nrOfWindows; 
	int i, j;

	int cntPairs = 0;

	if(alreadyAllocatedTo_nrOfWindows == 0){ //indicates that no allocation was done yet

		//ptr_Irarity_windowPairs = (struct Irarity_windowPairs *) malloc((N*(N-1)/2 + N)*sizeof(struct Irarity_windowPairs)); 
		ptr_Irarity_windowPairs = (struct Irarity_windowPairs *) malloc(N*N*sizeof(struct Irarity_windowPairs)); 
		/*ini segment indices to zero everywhere:*/

		for (i=0; i < N; i++){

			for (j=0; j < N; j++){

				ptr_Irarity_windowPairs[cntPairs].windowPair.windowNr_1 = 0; 

				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_1 = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_1[0] = 0;
				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_1[1] = 0;

				ptr_Irarity_windowPairs[cntPairs].windowPair.windowNr_2 = 0; 

				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_2 = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_2[0] = 0;
				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_2[1] = 0;

		
				/* Initialize to zero everywhere: */
				ptr_Irarity_windowPairs[cntPairs].rarityScore = 0.0;

				cntPairs += 1;
			}
		}

		*ptr_ptr_Irarity_windowPairs = ptr_Irarity_windowPairs;

	}
	else if(nrOfWindows > alreadyAllocatedTo_nrOfWindows){ //indicates that reallocation must be done

		ptr_Irarity_windowPairs = *ptr_ptr_Irarity_windowPairs;

		//ptr_Irarity_windowPairs = realloc(ptr_Irarity_windowPairs, (N*(N-1)/2 + N)*sizeof(struct Irarity_windowPairs)); 
		ptr_Irarity_windowPairs = realloc(ptr_Irarity_windowPairs, N*N*sizeof(struct Irarity_windowPairs)); 

		/*ini segment indices to zero everywhere:*/

		for (i=0; i < N; i++){

			for (j=0; j < N; j++){

				ptr_Irarity_windowPairs[cntPairs].windowPair.windowNr_1 = 0; 

				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_1 = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_1[0] = 0;
				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_1[1] = 0;

				ptr_Irarity_windowPairs[cntPairs].windowPair.windowNr_2 = 0; 

				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_2 = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_2[0] = 0;
				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_2[1] = 0;

		
				/* Initialize to zero everywhere: */
				ptr_Irarity_windowPairs[cntPairs].rarityScore = 0.0;

				cntPairs += 1;
			}
		}

		*ptr_ptr_Irarity_windowPairs = ptr_Irarity_windowPairs; //needed??

	} else{ //just re-init to zero everywhere:

		ptr_Irarity_windowPairs = *ptr_ptr_Irarity_windowPairs;

		/*ini segment indices to zero everywhere:*/

		for (i=0; i < N; i++){

			for (j=0; j < N; j++){

				ptr_Irarity_windowPairs[cntPairs].windowPair.windowNr_1 = 0; 

				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_1[0] = 0;
				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_1[1] = 0;

				ptr_Irarity_windowPairs[cntPairs].windowPair.windowNr_2 = 0; 

				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_2[0] = 0;
				ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_2[1] = 0;

		
				/* Initialize to zero everywhere: */
				ptr_Irarity_windowPairs[cntPairs].rarityScore = 0.0;

				cntPairs += 1;
			}
		}

		*ptr_ptr_Irarity_windowPairs = ptr_Irarity_windowPairs; //needed??

	}




	return N;

}


int init_ptr_Irarity_windowPairs(struct Irarity_windowPairs **ptr_ptr_Irarity_windowPairs,  int nrOfWindows){

	struct Irarity_windowPairs *ptr_Irarity_windowPairs;

	int N = nrOfWindows; 
	int i, j;

	int cntPairs = 0;

	//just re-init to zero everywhere:

	ptr_Irarity_windowPairs = *ptr_ptr_Irarity_windowPairs;

	/*ini segment indices to zero everywhere:*/

	for (i=0; i < N; i++){

		for (j=0; j < N; j++){

			ptr_Irarity_windowPairs[cntPairs].windowPair.windowNr_1 = 0; 

			ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_1[0] = 0;
			ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_1[1] = 0;

			ptr_Irarity_windowPairs[cntPairs].windowPair.windowNr_2 = 0; 

			ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_2[0] = 0;
			ptr_Irarity_windowPairs[cntPairs].windowPair.segIndices_2[1] = 0;

	
			/* Initialize to zero everywhere: */
			ptr_Irarity_windowPairs[cntPairs].rarityScore = 0.0;

			cntPairs += 1;
		}
	}

	*ptr_ptr_Irarity_windowPairs = ptr_Irarity_windowPairs; //needed??


	return N;

}

int alloc_init_Irarity_windowPairs_ptr(struct Irarity_windowPairs_ptr *ptr_Irarity_windowPairs_ptr, int order, int full_b, int chainLen, int windowLength, int nrOfWindows){


	struct Irarity_windowPairs_ptr Irarity_windowPairs;

	int N = nrOfWindows + 1; /*should suffice: the nrOfWindows will be set for serving the longest chain in the DB-set*/
	int i,j;

	Irarity_windowPairs.fileName = (char *) malloc (sizeof(char));
	Irarity_windowPairs.fileName = "NN";
	
	Irarity_windowPairs.structureName = (char *) malloc (sizeof(char)*1000);
	strcpy(Irarity_windowPairs.structureName, "NN");

	Irarity_windowPairs.classId = (char *) malloc (sizeof(char)*classIdCharSize);
	strcpy(Irarity_windowPairs.classId, "NN");

	Irarity_windowPairs.chainId = (char *) malloc (sizeof(char));
	Irarity_windowPairs.chainId = "NN";

	Irarity_windowPairs.chainNr = 0;
	Irarity_windowPairs.chainLen = chainLen;
	Irarity_windowPairs.order = order;
	Irarity_windowPairs.windowLgth = windowLength;
	Irarity_windowPairs.nrOfWindows = 0; /*this will store the true nr of windows involved (giving rise to nrOfWindowPairs) */
	Irarity_windowPairs.nrOfWindowPairs = 0; /*this will store the true nr of window pairs for which invariants are kept*/

	printf("nrOfWindowPairs allocated for rarity scoring:%d\n", (N+1)*(N+1));

	Irarity_windowPairs.windowPair = (struct windowPair **) malloc( (N+1)*(N+1)*sizeof(windowPair));

	for (i=0; i <= N; i++){

		for (j=0; j <= N; j++){

			Irarity_windowPairs.windowPair[i][j].windowNr_1 = 0;
			Irarity_windowPairs.windowPair[i][j].windowNr_2 = 0;

			Irarity_windowPairs.windowPair[i][j].segIndices_1[0] = 0;
			Irarity_windowPairs.windowPair[i][j].segIndices_1[1] = 0;

			Irarity_windowPairs.windowPair[i][j].segIndices_2[0] = 0;
			Irarity_windowPairs.windowPair[i][j].segIndices_2[1] = 0;

		}

	}

	///*Window 1:*/
	//Irarity_windowPairs.windowsNr_1 = (int *) malloc((N+1)*sizeof(int)); 
	//Irarity_windowPairs.segIndices_1 = (int **) malloc((N+1)*sizeof(int*)); 

	///*ini segment indices to zero everywhere:*/

	//for (i=0; i <= N; i++){

	//	Irarity_windowPairs.windowsNr_1[i] = 0; 

	//	Irarity_windowPairs.segIndices_1[i] = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

	//	Irarity_windowPairs.segIndices_1[i][0] = 0;
	//	Irarity_windowPairs.segIndices_1[i][1] = 0;

	//}


	///*Window 2:*/
	//Irarity_windowPairs.windowsNr_2 = (int *) malloc((N+1)*sizeof(int)); 
	//Irarity_windowPairs.segIndices_2 = (int **) malloc((N+1)*sizeof(int *)); 

	///*ini segment indices to zero everywhere:*/

	//for (i=0; i <= N; i++){

	//	Irarity_windowPairs.windowsNr_2[i] = 0; 

	//	Irarity_windowPairs.segIndices_2[i] = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

	//	Irarity_windowPairs.segIndices_2[i][0] = 0;
	//	Irarity_windowPairs.segIndices_2[i][1] = 0;

	//}


	Irarity_windowPairs.rarityScore = (double **) malloc ((N+1)*sizeof(double*));
					
	for (i=0; i <= N; i++){

		Irarity_windowPairs.rarityScore[i] = (double *)malloc((N+1) * sizeof(double));
			
		/* Initialize to zero everywhere: */
		for (j=0; j<= N; j++){
			Irarity_windowPairs.rarityScore[i][j] = 0.0;
		}


	}

	*ptr_Irarity_windowPairs_ptr = Irarity_windowPairs;

	return (N+1)*(N+1);

	}

int init_Irarity_windowPairs_ptr(struct Irarity_windowPairs_ptr Irarity_windowPairs, struct I_windowPairs_ptr *ptr_I_windowPairs){

	int k, l = 0;

	strcpy(Irarity_windowPairs.fileName, ptr_I_windowPairs -> fileName);
	strcpy(Irarity_windowPairs.structureName, ptr_I_windowPairs -> structureName);
	strcpy(Irarity_windowPairs.classId, ptr_I_windowPairs -> classId);
	strcpy(Irarity_windowPairs.chainId, ptr_I_windowPairs -> chainId);
	Irarity_windowPairs.chainNr = ptr_I_windowPairs -> chainNr;
	Irarity_windowPairs.chainLen = ptr_I_windowPairs -> chainLen;
	Irarity_windowPairs.order = ptr_I_windowPairs -> order;
	Irarity_windowPairs.nrOfWindows = ptr_I_windowPairs -> nrOfWindows;
	Irarity_windowPairs.windowLgth = ptr_I_windowPairs -> windowLgth;
	//Irarity_windowPairs.stepSize = ptr_I_windowPairs -> stepSize;

	for(k=0; k< Irarity_windowPairs.nrOfWindows;k++){
		
		/*Irarity_windowPairs.windowsNr_1[k] = ptr_I_windowPairs -> windowsNr_1[k];
		Irarity_windowPairs.segIndices_1[k][0] = ptr_I_windowPairs -> segIndices_1[k][0];
		Irarity_windowPairs.segIndices_1[k][1] = ptr_I_windowPairs -> segIndices_1[k][1];


		Irarity_windowPairs.windowsNr_2[k] = ptr_I_windowPairs -> windowsNr_2[k];
		Irarity_windowPairs.segIndices_2[k][0] = ptr_I_windowPairs -> segIndices_2[k][0];
		Irarity_windowPairs.segIndices_2[k][1] = ptr_I_windowPairs -> segIndices_2[k][1];*/

		for(l=0; l< Irarity_windowPairs.nrOfWindows;l++){

			Irarity_windowPairs.windowPair[k][l].windowNr_1 = ptr_I_windowPairs -> windowPair[k][l].windowNr_1;
			Irarity_windowPairs.windowPair[k][l].windowNr_2 = ptr_I_windowPairs -> windowPair[k][l].windowNr_2;

			Irarity_windowPairs.windowPair[k][l].segIndices_1[0] = ptr_I_windowPairs -> windowPair[k][l].segIndices_1[0];	
			Irarity_windowPairs.windowPair[k][l].segIndices_1[1] = ptr_I_windowPairs -> windowPair[k][l].segIndices_1[1];

			Irarity_windowPairs.windowPair[k][l].segIndices_2[0] = ptr_I_windowPairs -> windowPair[k][l].segIndices_2[0];	
			Irarity_windowPairs.windowPair[k][l].segIndices_2[1] = ptr_I_windowPairs -> windowPair[k][l].segIndices_2[1];

			Irarity_windowPairs.rarityScore[k][l] = 0.0;
		}

	}

}


int alloc_init_chainInStr(struct chainInStructure *ptr_chainInStr, int nrOfChains){

	struct chainInStructure chainInStr;

	char *structureNameInit = "NN";
	char *classIdInit = "NN";
	char chainIdInit[2] = ">";

	int i, j = 0;


	chainInStr.numberOfChains = nrOfChains;
	chainInStr.ptr_chainInStr = (struct chainInfo *) malloc (nrOfChains*sizeof(struct chainInfo));
	/*init*/
	for(i = 0; i< nrOfChains; i++){
		chainInStr.ptr_chainInStr[i].chainId = (char *) malloc(sizeof(char)*chainIdCharSize);
		chainInStr.ptr_chainInStr[i].structureName = (char *) malloc(sizeof(char)*structureNameCharSize);
		chainInStr.ptr_chainInStr[i].classId = (char *) malloc(sizeof(char)*classIdCharSize);
		chainInStr.ptr_chainInStr[i].classIdConv = (int *) malloc(classifierSize*sizeof(int));

		strcpy(chainInStr.ptr_chainInStr[i].chainId, chainIdInit); /*init to unlikely value*/
		strcpy(chainInStr.ptr_chainInStr[i].structureName, structureNameInit);
		strcpy(chainInStr.ptr_chainInStr[i].classId, classIdInit);
		chainInStr.ptr_chainInStr[i].chainLength = 0;
		/*printf("Chain Id initial value:%c\n", *chainInStr.ptr_chainInStr[i].chainId);
		printf("Class Id initial value:%c\n", *chainInStr.ptr_chainInStr[i].classId);*/
		for(j=0; j< classifierSize; j++){chainInStr.ptr_chainInStr[i].classIdConv[j] = -1;}

	}

	*ptr_chainInStr = chainInStr;

	return 0;

}

int reinit_chainInStr(struct chainInStructure chainInStr, int nrOfChains){

	char *structureNameInit = "NN";
	char *classIdInit = "NN";
	char chainIdInit[2] = ">";

	int i, j = 0;

	chainInStr.numberOfChains = nrOfChains;

	/*init*/
	for(i = 0; i< nrOfChains; i++){

		strcpy(chainInStr.ptr_chainInStr[i].chainId, chainIdInit); /*init to unlikely value*/
		strcpy(chainInStr.ptr_chainInStr[i].structureName, structureNameInit);
		strcpy(chainInStr.ptr_chainInStr[i].classId, classIdInit);
		chainInStr.ptr_chainInStr[i].chainLength = 0;
		/*printf("Chain Id initial value:%c\n", *chainInStr.ptr_chainInStr[i].chainId);
		printf("Class Id initial value:%c\n", *chainInStr.ptr_chainInStr[i].classId);*/
		for(j=0; j< classifierSize; j++){chainInStr.ptr_chainInStr[i].classIdConv[j] = -1;}

	}


	return 0;

}


/*
** functions for freeing memory allocated to pointers
*/

void freeDblArray(double **dblArray, size_t size){
    
	size_t i = 0;

	for (i = 0; i < size; i++){
        free(dblArray[i]);
	}
    free(dblArray);
}

void free_I_measures(struct I_ptr I_measures, int order, int full_b, int chainLen){

	int L = chainLen - 1;
	int i = 0;

	int j = 0;


	if (order == 0){
		for (i = 0; i < L+1; i++){
			free(I_measures.wVal[i]);
		}
		free(I_measures.wVal);
	}

	if (order == 1){

		for (i = 0; i < L + 1 ; i++){
			//printf("freeing I12... at %d\n", i);
			free(I_measures.wVal[i]);
			free(I_measures.I12[i]);
			free(I_measures.Ia12[i]);
		}
		free(I_measures.wVal);
		free(I_measures.I12);
		free(I_measures.Ia12);
		

	}

	if (order == 2){

		for (i = 0; i < L+1; i++){

			free(I_measures.wVal[i]);
			 
			free(I_measures.I12[i]);

		  	free(I_measures.I1234[i]);
			free(I_measures.I1324[i]);
			free(I_measures.I1423[i]);

			/*absolute value versions */
			free(I_measures.Ia12[i]);
			
			free(I_measures.Ia1234[i]);
			free(I_measures.I12a34[i]);
			free(I_measures.Ia12a34[i]);
				
			free(I_measures.Ia1324[i]);
			free(I_measures.I13a24[i]);
			free(I_measures.Ia13a24[i]);
				
			free(I_measures.Ia1423[i]);
			free(I_measures.I14a23[i]);
			free(I_measures.Ia14a23[i]);

			/*for "full" order 2*/
			if (full_b == 1){
				free(I_measures.I1234_full_aid[i]);
				free(I_measures.I1324_full_aid[i]);
				free(I_measures.I1234_full[i]);
				free(I_measures.I1324_full[i]);

				/*absolute value versions*/
				free(I_measures.Ia1234_full_aid[i]);
				free(I_measures.Ia1324_full_aid[i]);
				free(I_measures.Ia1234_full[i]);
				free(I_measures.Ia1324_full[i]);

				free(I_measures.I12a34_full_aid[i]);
				free(I_measures.I13a24_full_aid[i]);
				free(I_measures.I12a34_full[i]);
				free(I_measures.I13a24_full[i]);

				free(I_measures.Ia12a34_full_aid[i]);
				free(I_measures.Ia13a24_full_aid[i]);
				free(I_measures.Ia12a34_full[i]);
				free(I_measures.Ia13a24_full[i]);


			}

		}
		free(I_measures.wVal);
			 
		free(I_measures.I12);

		free(I_measures.I1234);
		free(I_measures.I1324);
		free(I_measures.I1423);

		/*absolute value versions */
		free(I_measures.Ia12);
			
		free(I_measures.Ia1234);
		free(I_measures.I12a34);
		free(I_measures.Ia12a34);
				
		free(I_measures.Ia1324);
		free(I_measures.I13a24);
		free(I_measures.Ia13a24);
				
		free(I_measures.Ia1423);
		free(I_measures.I14a23);
		free(I_measures.Ia14a23);

		/*for "full" order 2*/
		if (full_b == 1){

			free(I_measures.I1234_full_aid);
			free(I_measures.I1324_full_aid);
			free(I_measures.I1234_full);
			free(I_measures.I1324_full);

			/*absolute value versions*/
			free(I_measures.Ia1234_full_aid);
			free(I_measures.Ia1324_full_aid);
			free(I_measures.Ia1234_full);
			free(I_measures.Ia1324_full);

			free(I_measures.I12a34_full_aid);
			free(I_measures.I13a24_full_aid);
			free(I_measures.I12a34_full);
			free(I_measures.I13a24_full);

			free(I_measures.Ia12a34_full_aid);
			free(I_measures.Ia13a24_full_aid);
			free(I_measures.Ia12a34_full);
			free(I_measures.Ia13a24_full);
			

		}

	}

	if (order == 3){

		for (i = 0; i < L+1; i++){

			free(I_measures.wVal[i]);
			 
			free(I_measures.I12[i]);

		  	free(I_measures.I1234[i]);
			free(I_measures.I1324[i]);
			free(I_measures.I1423[i]);

			/*absolute value versions */
			free(I_measures.Ia12[i]);
			
			free(I_measures.Ia1234[i]);
			free(I_measures.I12a34[i]);
			free(I_measures.Ia12a34[i]);
				
			free(I_measures.Ia1324[i]);
			free(I_measures.I13a24[i]);
			free(I_measures.Ia13a24[i]);
				
			free(I_measures.Ia1423[i]);
			free(I_measures.I14a23[i]);
			free(I_measures.Ia14a23[i]);

			/*for "full" order 2*/
			free(I_measures.I1234_full_aid[i]);
			free(I_measures.I1324_full_aid[i]);
			free(I_measures.I1234_full[i]);
			free(I_measures.I1324_full[i]);

			/*... and the absolute value versions*/
			free(I_measures.Ia1234_full_aid[i]);
			free(I_measures.Ia1324_full_aid[i]);
			free(I_measures.Ia1234_full[i]);
			free(I_measures.Ia1324_full[i]);

			free(I_measures.I12a34_full_aid[i]);
			free(I_measures.I13a24_full_aid[i]);
			free(I_measures.I12a34_full[i]);
			free(I_measures.I13a24_full[i]);

			free(I_measures.Ia12a34_full_aid[i]);
			free(I_measures.Ia13a24_full_aid[i]);
			free(I_measures.Ia12a34_full[i]);
			free(I_measures.Ia13a24_full[i]);

			/*order 3*/
			/*12*/
			free(I_measures.I123456[i]);
			free(I_measures.I123645[i]);
			free(I_measures.I123546[i]);
			/*13*/
			free(I_measures.I132456[i]);
			free(I_measures.I132546[i]);
			free(I_measures.I132645[i]);
			/*14*/
			free(I_measures.I142356[i]);
			free(I_measures.I142536[i]);
			free(I_measures.I142635[i]);
			/*15*/
			free(I_measures.I152346[i]);
			free(I_measures.I152436[i]);
			free(I_measures.I152634[i]);
			/*16*/
			free(I_measures.I162345[i]);
			free(I_measures.I162435[i]);
			free(I_measures.I162534[i]);

		}

		free(I_measures.wVal);
			 
		free(I_measures.I12);

		free(I_measures.I1234);
		free(I_measures.I1324);
		free(I_measures.I1423);

		/*absolute value versions */
		free(I_measures.Ia12);
			
		free(I_measures.Ia1234);
		free(I_measures.I12a34);
		free(I_measures.Ia12a34);
				
		free(I_measures.Ia1324);
		free(I_measures.I13a24);
		free(I_measures.Ia13a24);
				
		free(I_measures.Ia1423);
		free(I_measures.I14a23);
		free(I_measures.Ia14a23);

		/*for "full" order 2*/
		free(I_measures.I1234_full_aid);
		free(I_measures.I1324_full_aid);
		free(I_measures.I1234_full);
		free(I_measures.I1324_full);

		/*... and the absolute value versions*/
		free(I_measures.Ia1234_full_aid);
		free(I_measures.Ia1324_full_aid);
		free(I_measures.Ia1234_full);
		free(I_measures.Ia1324_full);

		free(I_measures.I12a34_full_aid);
		free(I_measures.I13a24_full_aid);
		free(I_measures.I12a34_full);
		free(I_measures.I13a24_full);

		free(I_measures.Ia12a34_full_aid);
		free(I_measures.Ia13a24_full_aid);
		free(I_measures.Ia12a34_full);
		free(I_measures.Ia13a24_full);

		/*order 3*/
		/*12*/
		free(I_measures.I123456);
		free(I_measures.I123645);
		free(I_measures.I123546);
		/*13*/
		free(I_measures.I132456);
		free(I_measures.I132546);
		free(I_measures.I132645);
		/*14*/
		free(I_measures.I142356);
		free(I_measures.I142536);
		free(I_measures.I142635);
		/*15*/
		free(I_measures.I152346);
		free(I_measures.I152436);
		free(I_measures.I152634);
		/*16*/
		free(I_measures.I162345);
		free(I_measures.I162435);
		free(I_measures.I162534);

	}

	

	free(I_measures.fileName);
	free(I_measures.structureName);
	free(I_measures.chainId);
	free(I_measures.classId);

	printf("done freeing I-mea ...\n");

}

void free_I_values(struct I_values **ptr_I_values, int numberOfFiles){

	int i = 0;
	int j = 0;

	for (i=0; i< numberOfFiles; i++){
		free(ptr_I_values[i]);
	}
	free(ptr_I_values);
}


int main_GISA(){

	int returnVal = 0;

	int maxChainLength = 1500; /*chains longer than this will be discarded; can be set acc to RAM availability*/

	/*compute invariants up to and including this order:*/
	int order = 2;
	/*inlcude absolute value versions of order one and two invariants:*/
	int incl_abs_b = 1;
	/* compute full order 2 measures across the simplex*/
	int full_b = 1;

	/*Detect and record closed loops; if chosen I12 values corr 
	to the loop will be stored for alter use (finding "pokes").
	Can only be used with order >= 1 (since I12 is needed):*/
	int closed_loops_b = 0;
	int closedLoopLength = 30;
	int pokeLength = 10;
	double closedLoopDist = 7*7; /*square of distance in nr of Ångstrøm -- we use square of dist*/
	/*For computing and recording the mutual invariant values (e.g writhe) of pairs of sub-chains (1) or not (0); aimed
	at more generally searching for particular geometric shapes than with the closed-loops
	search:*/
	int invValSubChainPairs_b = 0; /*if set to 1: computes the mutual writhe; if set to 2 computes the mutual invariant value for order 1 and order 2 "full" invariant (takes setting order = 3 or full_b = 1); if set to 3: all invariants (takes setting order = 3)*/
	int writeSubChainPairsAll_b = 0; /*if set to 1 the results for all sub-chain pairs will be attempted to be written out (may be lots of data!)*/
	int subChainLength = 30; /* the lenght of the sub-chains in the pairs considered*/
	

	/*split versions: w computation split away from aggregation:*/
	int split_b = 0;
	/* print various intermediates on the screen (1):*/
	int print_b = 0;
	/*print put basic stuff to the screen:*/
	int print_basic_b = 0;
	/*write out to .txt file the values of the measures on the given chain (i.e. the "top corner" value for each chain)*/
	int write_final_b = 0;
	/*write out to .txt file the values of the measures on the given chain
	and all its sub-chains. Obs: this can become very large data sets*/
	int write_all_b = 0;


	/*name of folder holding the PDB files (the input):*/
	int loadFromSubDirs_b = 0; //if PDB files sit in sub-directories at some level of the "root" directory, dirName, use loadFromSubDirs_b = 1; else 0.
	int use_scop_b = 0; //if using data from SCOP (particular data structure when loading files)
	int use_cath_b =0;//if using data from CATH (particular data structure when loading files)
	char DBName[100] = "top100";
	//char DBName[100] = "top8000";
	//char DBName[100] = "SCOPe_ASTRAL_2.06_40pct";
	//char DBName[100] = "SCOP_test";
	//char DBName[100] = "CATH_v2_4_Kolodny";

	char dirPath[1000] = "/data/tkj375/data/structural/Kinemage/top100H";

	//char CATHListPath[1000] = "/data/tkj375/data/structural/CATH/CATH-v2_4_Kolodny/subset_list_web.txt";
	char CATHListPath[1000] = "NN";

	char outputPath[1000]  = "/data/tkj375/GISA_unix/results_C_7th/double_precision";


	//char fileNameChains[1000] = "//chains_subChainPairs_top100.txt";
	char fileNameChains[1000] = "/chains_top100.txt";


	/*For getting and writing out the invariant values for a set windowlength:*/
	int get_windows_b = 0;
	int get_windowPairs_b = 0; /*at most one of get_windows_b and get_windowPairs_b should be set to 1*/
	int windowLength = 16;
	int write_windows_b = 0;
	int write_windowPairs_b = 0; /*at most one of write_windows_b and write_windowPairs_b should be set to 1; corr to get_windos_b resp. get_windowPairs_b*/
	int maxNrOfWindows = 1;


	returnVal = computeGI(DBName, dirPath, use_scop_b, 
							use_cath_b, CATHListPath, loadFromSubDirs_b, 
							outputPath, maxChainLength, order, 
							incl_abs_b, full_b, split_b,
							write_final_b, write_all_b, closed_loops_b, 
							closedLoopLength, closedLoopDist, pokeLength,invValSubChainPairs_b, 
							writeSubChainPairsAll_b, subChainLength, print_b, 
							print_basic_b);

	getchar();

	return returnVal;

}
 
/*Main function: 
Computes all invariants up to and including the desired order (e.g. order = 2) for 
all PDB-files in a directory (defined by dirPath). Results are written out as desired: 
for each PDB-file, all results (write_all_b = 1) meaning the invariants' values for all 
vertices in the simplex or only final values (write_final_b =1) meaning the invariants' values 
of the given structure (PDB-file), i.e. the value at the top left corner of the simplex. 
(Simplex: each vertex (i,j) represents a sub-chain of the structure; see also the Supplementary 
for an explanation). 
As an additional feature it is possible to envoke a particular version of the computation of the 
invariants in which closed loops are detected and recorded (ie a restricted search), along with the 
writhe value (I12) of the loop. For running this set closed_loops_b to 1. A far more general version 
of this (a unrestricted search) is also callable in which the mutual writhe (I12) of all (disjoint) 
sub-chains of a specified length is computed and recorded (: written to a file).


Input:
For the usage/params of the compiled code use the --help on GISA_main_*. 
But here is some of it (if not all, and maybe not absolutely up to date):

Input:
char DBName[100]: name of data base (e.g. 'PDB')
char *dirPath: path to (top) directory in which the PDB-files of the DB are located 
char *outputPath: path to folder to hold the output
char *fileNameChains: name of file containing the chain info written out

int order: invariants will be computed up to and including this value
int incl_abs_b: if 1, absolute value versions of invariants of order 2 and below will be included
int full_b: if 1, all order 2 measures will be computed across the simplex
int split_b: if 1 the w computation is split away from aggregation

int write_final_b: if 1, writes out to .txt file the values of the measures on the given chain. Name 
of file is generated inside the code, but path to folder is set by outputPath

int write_all_b: writes out to .txt file the values of the measures on the given chain
and all its sub-chains. Obs: this can become very large data sets. Name of file is set inside the code, 
but path to folder is set by outputPath

int closed_loops_b: if 1, closed loops will be detected and recorded (file written). 
Works only with order at least 1 (and works fastest with 1).

int closedLoopLength: max length of closed loops in number of segments; a number of residues (30)
int closedLoopDist: square of distance; used for qualifying sub-chains as closed loops: if 
end pts of a particular sub-schain are within a squared distance less than this number, the 
sub-chain is recorded as a closed loop.

int pokeLength: the length of the pokes to be searched for (number of residues).

double thresholdLinks: only cases (pairs of closed loops) having a writhe which in absolute value is 
above this threshold are written to file. (Only effective if closed_loops_b = 1.)
double thresholdPokes: only cases (potential pokes, ie a closed loops and a sub-chain) having a writhe 
which in absolute value is above this threshold are written to file.(Only effective if closed_loops_b = 1.)

int invValSubChainPairs_b: if set to 1, computes the mutual writhe; if set to 2, computes the mutual invariant value 
for order 1 and order 2 "full" invariant (takes setting order = 3 or full_b = 1); if set to 3, computes for all invariants 
(takes setting order = 3)

int writeSubChainPairsAll_b: if set to 1 the results for all sub-chain pairs will be attempted to be written out (may be lots of data!)
int subChainLength: the lenght of the sub-chains in the pairs considered		  
double thresholdWrithe: only cases (pairs of sub-chains) having a writhe which in absolute value is 
above this threshold are written to file. (Only effective if invValSubChainPairs_b = 1.)

int print_b: if 1 some intermediate values are output to the screen
int print_basic_b: if 1 some basic intermediate values are output to the screen (e.g. length of loaded chain)

Output:
The function writes results (invariant values) to .txt files as set with the inputs: write_all_b, write_final_b, 
closed_loops_b and writheSubChainPairs_b. If desired the function can write excerpts of results to the screen 
(use print_b or print_basic_b). The function also writes certain time consumption values to the screen at the 
end of the run.

Notes:
Memory sufficient for holding the invariants' values of the the chains is allocated for in blocks/chunks of structures
in the given directory (of PDB-files); for each block the longest chain is observed and corresponding memeory
is allocated; this allocation is done up front. 
The directory is looped through and, for each structure (chain or, possibly, sub-chain of the given PDB-file) 
the invariants are computed after having being reinitialized to zero everywhere relevant (i.e. set to zero in 
the simplex of size given by the current structure). For each structure (chain or, possibly, sub-chain of the 
given PDB-file) memory is allocated for the chain of segments and initialized with values as read in from the 
PDB-file. This memory is freed again at the completion of the computations/write-outs of the structure (chain or 
sub-chain of PDB-file).

Disregarding chains with "holes": chains will be read in eventhough the series of residue indices in the PDB-file is
interrupted. A chain will though be disregarded if it contains one or more "geometric holes", meaning segments 
in the alpha-C trace of length longer than the square root of globally set stdRealSegLength (should not be set 
lower than square of 3.5) 

The allocated memory (for the invariants' values and some more) is freed after the completion of the last 
structure in the directory. 
*/
int computeGI(char DBName[100], char *dirPath, int use_scop_b, 
				int use_cath_b, char *cathListPath, int loadFromSubDirs_b, 
				char *outputPath, int maxChainLength, int order, 
				int incl_abs_b, int full_b, int split_b, 
				int write_final_b, int write_all_b, int closed_loops_b, 
				int closedLoopLength, int closedLoopDist, int pokeLength,
				int invValSubChainPairs_b, int writeSubChainPairsAll_b, int subChainLength, 
				int print_b, int print_basic_b){

	int returnval = 0;

	//int strLengthCutOff = 7; /*chains shorter than this int will not be included in computation of invariants*/

	struct dirContent dirContent;
	char ** ptr_dirList;
	char ** ptr_fileNameList;
	int subDirCnt = 0;
	int numberOfFiles = 0;
	int subStructureCnt = 0; /*some structures are multi-mers; we compute the invariants for each sub-structure (each mono-mer)*/
	const int numberOfPerturbations = 1; /*only there since needed when writing out the final values (preparation for perturbations code)*/
	const int pertNo = 0; /*ditto*/

	
	/*for "dynamically" generating the file names for output:*/
	char orderStr[10] = "-1"; //value is placeholder
	char inclAbsStr[10] = "-1"; //value is placeholder

	char fileNameOut1[1000] = "/Ivalues_order_";
	char fileNameOutAll1[1000] = "/All_Ivalues_order_";
	char fileNameOutClosedLoops1[1000] = "/ClosedLoopChars"; 
	char fileNameOutSubChainPairs1[1000] = "/SubChainPairChars"; 
	char fileNameOut2[1000] = "_minmax_"; //; //; //"_incl_abs_"; //"_"; //

	//char fileNameOut3[1000] = "_mem_alloc_top100.txt"; 
	//char fileNameOutChains1[1000] = "\\chains_closed_loops_top8000.txt";
	char fileNameOutChains1[1000] =  "/chains_closed_loops_";
	char fileNameOutChains2[1000] =  "/chains_subchain_pairs_";
	char fileNameOut3[1000] = "computeGI_"; 
	char extensionName[20] = ".txt";

	char fileNameOutAll[2000] = "\0";
	char fileNameOut[2000] = "\0";
	char fileNameOutClosedLoops[2000] = "\0";
	char fileNameOutSubChainPairs[2000] = "\0";
	char fileNameOutChains[2000] = "\0";

	char *ptr_fileNameOut; /*pointer to file for holding final measure values (i.e. at final corner of simplex)*/ 
	char *ptr_fileNameOutAll; /*pointer to file for holding all measure values (i.e. across whole simplex)*/ 
	char *ptr_fileNameOutClosedLoops; /*pointer to file for holding mutual writhe value of closed loop pairs*/ 
	char *ptr_fileNameOutSubChainPairs; /*pointer to file for holding mutual writhe value of sub-chain pairs*/ 
	char *ptr_fileNameOutChains; /*pointer to file for holding chain coordinates*/

	char closedLoopLengthStr[10] = "\0"; /*for holding closedLoopLength converted to a string*/
	char pokeLengthStr[10] = "\0"; /*for holding pokeLength converted to a string*/
	char closedLoopDistStr[10] = "\0"; /*for holding closedLoopdist converted to a string*/
	char subChainLengthStr[10] = "\0"; /*for holding subChainLength converted to a string*/
	
	clock_t startChain, endChain, startComp, endComp, startComplete, endComplete;
	double timeChain = 0, compTime = 0, compTimeComplete = 0, timeComplete = 0, max_compTime = 0, max_timeChain = 0;
	//long double timeChain = 0, compTime = 0, compTimeComplete = 0, timeComplete = 0, max_compTime = 0, max_timeChain = 0;

	FILE *ptr_fileIn;

	FILE *ptr_cath_domain_list_file;
	int convToCLF20_b = 1;
	struct chainInfo *ptr_cathDomain;
	int nrOfCATHDomainsInList;
	int hitIndexDomain = -1;
	int lastHitIndexDomain = -1;

	struct chainInStructure chainInStr;
	char *currentStructureName;
	char *structureNameInit = "NN";
	char *currentClassId;
	char *classIdInit = "NN";
	char chainIdInit[2] = ">";
	int chainsSkippedTooShort = 0;
	int chainsSkippedTooLong = 0;
	int chainSkippedSegmTooLong = 0;
	int segTooLong_b = 0;
	int resNr = 0; /*residue number*/
	int resNrPrior = -999;
	int chainLen = 0; /*length of C-alpha chain*/
	int maxChainLen =0; /*max length of chain among all loaded*/
	int L = 0; /*length of segment chain = chainLen -1*/
	int simplexCnt = 0;/*size of simplex*/
	int maxSimplexCnt = 0; /*size of simplex corr to max chain length*/

	struct cAlpha *ptr_chain = NULL; /*for holding C-alpha chain */
	struct segment *ptr_segment = NULL; /* for holding chain of segments, ith C-alpha to jth C-alpha*/
	struct segment segCoords;	 

	int n = 0;
	int i = 0;
	int j = 0;
	int cnt = 0;
	int fileNr = 0;
	int chainNr = 0;

	double avgChainLength = 0;
	int nrOfChainsAbove500 = 0;
	int nrOfChainsAbove1000 = 0;

	/* for holding meausure values (i.e. the core of the final output); 
	content to be received from aggr function:*/
	struct I_ptr I_measures; 

	int returnVal = 0;

	/*pointer to hold the results (that this is a ptr of a ptr is due to that the write-out can also
	handle the perturbation case; we use pertNo = 0 here):*/
	struct I_values **ptr_I_values;

	/*for closed loops finding:*/
	struct twoSegmentIndex *ptr_closedLoopInd;
	int initValNrOfClosedLoops = 0;
	int nrOfClosedLoops = 0;
	int cntClosedLoop = 0;

	/****************************************************************************************************************/
	/*Generation of output file names:*/
	/****************************************************************************************************************/

	if (write_all_b == 1){

		sprintf(orderStr, "%d", order);
		strcat(fileNameOutAll, outputPath);
		strcat(fileNameOutAll, fileNameOutAll1);
		strcat(fileNameOutAll, orderStr);
		sprintf(inclAbsStr, "%d", incl_abs_b);
		strcat(fileNameOutAll, inclAbsStr);
		strcat(fileNameOutAll, "_");
		strcat(fileNameOutAll, fileNameOut2);
		strcat(fileNameOutAll, fileNameOut3);
		if (DBName != ""){
			strcat(fileNameOutAll, DBName);
		}
		strcat(fileNameOutAll, extensionName);
		printf("Results for all vertices will be written to: %s\n", fileNameOutAll);

		ptr_fileNameOutAll = fileNameOutAll;
		
		/*clear the file by opening it in w-mode and then closing it again*/
	    fclose(fopen(ptr_fileNameOutAll, "w"));

	}

	if (closed_loops_b == 1){

		/*generate file name for writing out closed loops/poke details:*/
		sprintf(pokeLengthStr,"%d",pokeLength);
		sprintf(closedLoopLengthStr,"%d",closedLoopLength);

		strcat(fileNameOutClosedLoops, outputPath);
		strcat(fileNameOutClosedLoops, fileNameOutClosedLoops1);

		strcat(fileNameOutClosedLoops, fileNameOut2);
		strcat(fileNameOutClosedLoops, fileNameOut3);

		if (DBName != ""){
			strcat(fileNameOutClosedLoops, DBName);
			strcat(fileNameOutClosedLoops, "_");
		}

		strcat(fileNameOutClosedLoops, "closedLoopLength_");
		sprintf(closedLoopLengthStr, "%d", closedLoopLength);
		strcat(fileNameOutClosedLoops, closedLoopLengthStr);

		strcat(fileNameOutClosedLoops, "_closedLoopDist_");
		sprintf(closedLoopDistStr, "%0.1f", sqrt(closedLoopDist));
		strcat(fileNameOutClosedLoops, closedLoopDistStr);
		strcat(fileNameOutClosedLoops, "_pokeLength");
		strcat(fileNameOutClosedLoops, pokeLengthStr);
		
		strcat(fileNameOutClosedLoops, extensionName);
		printf("Results for closed loops will be written to: %s\n", fileNameOutClosedLoops);

		ptr_fileNameOutClosedLoops = fileNameOutClosedLoops;
		
		/*clear the file by opening it in w-mode and then closing it again*/
	    fclose(fopen(ptr_fileNameOutClosedLoops, "w"));

		/*generate file name for writing out chains:*/
		strcat(fileNameOutChains, outputPath);
		strcat(fileNameOutChains, fileNameOutChains1);
		if (DBName != ""){
			strcat(fileNameOutChains, DBName);
		}
		strcat(fileNameOutChains, extensionName);		
	
		ptr_fileNameOutChains = fileNameOutChains;
		
		/*clear the file by opening it in w-mode and then closing it again*/
	    fclose(fopen(ptr_fileNameOutChains, "w"));

	}
	
	if (invValSubChainPairs_b != 0 ){

		/*generate file name for writing out writhe results for sub-chain pairs:*/
		strcat(fileNameOutSubChainPairs, outputPath);
		strcat(fileNameOutSubChainPairs, fileNameOutSubChainPairs1);
		strcat(fileNameOutSubChainPairs, fileNameOut2);
		strcat(fileNameOutSubChainPairs, fileNameOut3);

		if (DBName != ""){
			strcat(fileNameOutSubChainPairs, DBName);
			strcat(fileNameOutSubChainPairs, "_");
		}

		strcat(fileNameOutSubChainPairs, "subChainLength_");
		sprintf(subChainLengthStr, "%d", subChainLength);
		strcat(fileNameOutSubChainPairs, subChainLengthStr);

		strcat(fileNameOutSubChainPairs, extensionName);
		printf("Results for sub-chain pairs will be written to: %s\n", fileNameOutSubChainPairs);

		ptr_fileNameOutSubChainPairs = fileNameOutSubChainPairs;
		
		/*clear the file by opening it in w-mode and then closing it again*/
	    fclose(fopen(ptr_fileNameOutSubChainPairs, "w"));

		/*generate file name for writing out chains:*/
		strcat(fileNameOutChains, outputPath);
		strcat(fileNameOutChains, fileNameOutChains2);
		if (DBName != ""){
			strcat(fileNameOutChains, DBName);
		}
		strcat(fileNameOutChains, extensionName);	
		
		ptr_fileNameOutChains = fileNameOutChains;
		
		/*clear the file by opening it in w-mode and then closing it again*/
	    fclose(fopen(ptr_fileNameOutChains, "w"));

	}


	/****************************************************************************************************************/
	/*Get content of the desired directory: number of files and a list of the file names:*/
	/****************************************************************************************************************/
	
	if(loadFromSubDirs_b ==0){
		dirContent = ListDirectoryContents(dirPath);
	}
	else{

		returnVal = ListDirectoryContents2(dirPath, &dirContent, &subDirCnt);
	}

	numberOfFiles = dirContent.numberOfFiles; 
	ptr_dirList = dirContent.ptr_dirList;
	ptr_fileNameList = dirContent.ptr_fileNameList;

	printf("Number of files in directory:%d \n",numberOfFiles);
	//printf("First file:%s\n",ptr_dirList[0]);

	startComplete = clock();


	/****************************************************************************************************************/
	/*loop through the list of filenames in the directory to find chain info, eg the max chain length (for setting
	the size in the (global) mem alloc to I_measures's ptr's). First though we make a bulk allocation
	to the chain-info keeping pointer, using the pre-set maxNrOfChains. We also allocate to a 
	ptr keeping the class info for every chain/structure:*/
	/****************************************************************************************************************/

	chainInStr.numberOfChains = maxNrOfChains;
	chainInStr.ptr_chainInStr = (struct chainInfo *) malloc (maxNrOfChains*sizeof(struct chainInfo));
	/*init*/
	for(i = 0; i< maxNrOfChains; i++){
		chainInStr.ptr_chainInStr[i].chainId = (char *) malloc(sizeof(char)*2);
		chainInStr.ptr_chainInStr[i].structureName = (char *) malloc(sizeof(char)*10);
		chainInStr.ptr_chainInStr[i].classId = (char *) malloc(sizeof(char)*classIdCharSize);
		chainInStr.ptr_chainInStr[i].classIdConv = (int *) malloc(classifierSize*sizeof(int));


		strcpy(chainInStr.ptr_chainInStr[i].chainId, chainIdInit); /*init to unlikely value*/
		strcpy(chainInStr.ptr_chainInStr[i].structureName, structureNameInit);
		strcpy(chainInStr.ptr_chainInStr[i].classId, classIdInit);
		/*printf("Chain Id initial value:%c\n", *chainInStr.ptr_chainInStr[i].chainId);
		printf("Class Id initial value:%c\n", *chainInStr.ptr_chainInStr[i].classId);*/
		for(j=0; j< classifierSize; j++){chainInStr.ptr_chainInStr[i].classIdConv[j] = -1;}

	}

	/*Fork: whether we use data from CATH or don't. In CATH data domain info (name, class id, chain lenght) is
	stored in file separately from the file containing the PDB 3d-coord data. In SCOP and Kinemage data
	the PDB-files contain the info (either directly or, else, extracted by the code).*/

	/*if CATH data, read in CATH domain list. Note: we abuse the chainInStrcture struct to accomodate
	domains rather than chains -- in CATH each chain may be split in several domains*/ 
	if(use_cath_b == 1){


		/*first allocate a block of memory (if not large enough it will be extended when exectuing readCATHDomainList)*/
		ptr_cathDomain = (struct chainInfo *) malloc (maxNrOfDomains*sizeof(struct chainInfo));
		/*init*/
		for(i = 0; i< maxNrOfDomains; i++){
			ptr_cathDomain[i].chainId = (char *) malloc(sizeof(char)*2);
			ptr_cathDomain[i].structureName = (char *) malloc(sizeof(char)*10); //to hold domain name
			ptr_cathDomain[i].classId = (char *) malloc(sizeof(char)*classIdCharSize);
			ptr_cathDomain[i].classIdConv = (int *) malloc(classifierSize*sizeof(int));

			strcpy(ptr_cathDomain[i].chainId, chainIdInit); /*init to unlikely value*/
			strcpy(ptr_cathDomain[i].structureName, structureNameInit);
			strcpy(ptr_cathDomain[i].classId, classIdInit);
			/*printf("Chain Id initial value:%c\n", *ptr_cathDomain[i].chainId);
			printf("Class Id initial value:%c\n", *ptr_cathDomain[i].classId);*/
			for(j=0; j< classifierSize; j++){ptr_cathDomain[i].classIdConv[j] = -1;}
			ptr_cathDomain[i].chainLength = -1;

		}

		/*open CATH file list:*/
		ptr_cath_domain_list_file = fopen(cathListPath, "r");

		/*read in contents of CATH domain list file*/
		nrOfCATHDomainsInList = readCATHDomainList(ptr_cath_domain_list_file, &ptr_cathDomain, convToCLF20_b );

		printf("nr of cath doms: %d\n", nrOfCATHDomainsInList);

		/*lex-sort the ptr_cath_domain_list_file on the domain name entry. Below we
		need to look up each file for which we have PDB-coords in this list. In other
		words, this list serves as a "positive list" of domains, for which to compute
		the invariants.*/
		heapSortCATHDomain(ptr_cathDomain, nrOfCATHDomainsInList, lengthCLFformat20,0);

		//for(i=0;i <nrOfCATHDomainsInList;i++){printf("dom %d: %s\n", i, ptr_cathDomain[i].structureName);}

	}

	/*loop through the list of filenames in the directory to find the max chain length (for setting
	the size in the (global) mem alloc to I_measures's ptr's). (As some CATH domain lists (like the one
	from Kolodny et al's web site) do not contain the chain lengths, we also include a read-through 
	of CATH pdb-files here.)*/

	for (fileNr = 0; fileNr <numberOfFiles; fileNr ++){

		if (print_b == 1){
			printf("File name:%s fileNr:%d\n",ptr_dirList[fileNr], fileNr);
		}

		/*if using CATH data: check if the file, ptr_fileNameList[fileNr], is in the
		domain list specified (and loaded and lex-sorted above)*/
		if(use_cath_b == 1){
			hitIndexDomain = bisectionalSearchCATHDomain(ptr_fileNameList[fileNr], ptr_cathDomain, nrOfCATHDomainsInList, lengthCLFformat20);
			if(hitIndexDomain == -1){
				//printf("file not in domain list\n");
				continue;} //file name was not in the list
			else{ //Not really needed here, since we don't capture class info in this function: read into chainInStr some available values from ptr_cathDomain[hitIndexDomain]; obs: each CATH domain sits in one chain only (therefore: ptr_chainInStr[0]) 
				chainInStr.ptr_chainInStr[0].structureName = ptr_cathDomain[hitIndexDomain].structureName;
				chainInStr.ptr_chainInStr[0].classId = ptr_cathDomain[hitIndexDomain].classId;
				chainInStr.ptr_chainInStr[0].classIdConv = ptr_cathDomain[hitIndexDomain].classIdConv;
			}
		}

		//printf("hitIndexDomain:%d\n", hitIndexDomain);

		ptr_fileIn = fopen(ptr_dirList[fileNr], "r");


		if(!ptr_fileIn){
			printf("Sorry, the file: %s, fileNr: %d could not be read\n", ptr_dirList[fileNr], fileNr);
			continue;
		}


		/*get number of chains, their id's and lengths for this file:*/
		//chainInStr = readPDBChainStructure(ptr_fileIn);
		if(use_scop_b ==1){
			readPDBChainStructureSCOP(ptr_fileIn, &chainInStr );
		}
		else if(use_cath_b ==1){
			readPDBDomainStructureCATH(ptr_fileIn, &chainInStr ); //will read chain length and chain id into chainInStr
		}
		else{
			readPDBChainStructure2(ptr_fileIn, &chainInStr );
		}

		fclose(ptr_fileIn);

		//printf("structure: %s has chainInStr.numberOfChains:%d of length: %d\n",ptr_fileNameList[fileNr], chainInStr.numberOfChains, chainInStr.ptr_chainInStr[0].chainLength);

		for (i=0;i <= chainInStr.numberOfChains -1;i++){
			if(chainInStr.ptr_chainInStr[i].chainLength <= maxChainLength){ //we disregard very long chains for RAM availability reasons
				maxChainLen = max(chainInStr.ptr_chainInStr[i].chainLength,maxChainLen);
				if(chainInStr.ptr_chainInStr[i].chainLength >= strLengthCutOff){
					subStructureCnt += 1; /*count and index the sub-structures for which the invariants will be computed*/

					avgChainLength += chainInStr.ptr_chainInStr[i].chainLength;
					if(chainInStr.ptr_chainInStr[i].chainLength >500){
						nrOfChainsAbove500 += 1;
						if(chainInStr.ptr_chainInStr[i].chainLength >1000){
							nrOfChainsAbove1000 += 1;
						}
					}
				}
				else{printf("structure: %s is shorter than %d\n", ptr_fileNameList[fileNr], strLengthCutOff);}
			}
			else{printf("structure: %s is longer than %d\n", ptr_fileNameList[fileNr], maxChainLength);}
		}

	}

	L = maxChainLen -1;

	printf("max chain len:%d\n", maxChainLen);
	printf("sub structure count:%d\n", subStructureCnt);

	avgChainLength = (double) avgChainLength/subStructureCnt;

	printf("Avg chain length: %lf\n", avgChainLength);
	//printf("Nr of chains above lgth 500: %d, legth 1000: %d\n", nrOfChainsAbove500, nrOfChainsAbove1000);
	//return 0;


	/****************************************************************************************************************/
	/*MEMORY ALLOCATION*/
	/****************************************************************************************************************/

	/*We allocate to longest chain*/
	returnVal = alloc_init_I_measures(&I_measures, order, full_b, chainNr, maxChainLen, closed_loops_b, 0, &ptr_closedLoopInd);

	ptr_chain = (struct cAlpha *) malloc (maxChainLen*sizeof(struct cAlpha));
	/*init*/
	for ( n=0; n <= maxChainLen -1 ; n++ ){
		
		ptr_chain[n].coords.x = 0.0;
		ptr_chain[n].coords.y = 0.0;
		ptr_chain[n].coords.z = 0.0;

		ptr_chain[n].residueNr = -123;

	}
	ptr_segment = (struct  segment *) malloc (maxChainLen*sizeof(struct segment)); /* for the case that re-allocating of mem to ptr_segment is used instead*/
	/*init*/
	segCoords.s1.x = 0.0;
	segCoords.s1.y = 0.0;
	segCoords.s1.z = 0.0;
	segCoords.s2.x = 0.0;
	segCoords.s2.y = 0.0;
	segCoords.s2.z = 0.0;
	for ( n=0; n <= maxChainLen -1 ; n++ ){
			segCoords.s1 = ptr_chain[n].coords;
			segCoords.s2 = ptr_chain[n+1].coords;

			ptr_segment[n] = segCoords;

	}

	if (write_final_b == 1){

		/*allocate memory to pointer ptr_I_values and initialize (obs: numberOfPerturbations is preset to 1):*/
		returnVal = alloc_init_I_values(&ptr_I_values, subStructureCnt, numberOfPerturbations);

	}

	subStructureCnt =0; /*reset*/
	

	//startComplete = clock();

	lastHitIndexDomain = 0;

	/****************************************************************************************************************/
	/* MAIN LOOP */
	/****************************************************************************************************************/

	/*loop through the list of filenames in the directory and compute the desired measures all along:*/
	for (fileNr = 0; fileNr <numberOfFiles; fileNr ++){

		if (fileNr%100 == 0){printf("Now at file nr:%d\n", fileNr);}

		if (print_basic_b == 1){
			printf("Now handling file:%s\nfileNr:%d\n",ptr_dirList[fileNr], fileNr);
		}

		if(use_cath_b == 1){
			hitIndexDomain = bisectionalSearchCATHDomain(ptr_fileNameList[fileNr], ptr_cathDomain, nrOfCATHDomainsInList, lengthCLFformat20);
			/*if((hitIndexDomain - lastHitIndexDomain) >1){printf("Str:%s missing\n", ptr_cathDomain[hitIndexDomain -1]);}
			lastHitIndexDomain = hitIndexDomain;
			getchar();*/
			if(hitIndexDomain == -1){
				//printf("file not in domain list\n");
				continue;} //file name was not in the list
			else{ //read into chainInStr some available values from ptr_cathDomain[hitIndexDomain]; obs: each CATH domain sits in one chain only (therefore: ptr_chainInStr[0]) 
				chainInStr.ptr_chainInStr[0].structureName = ptr_cathDomain[hitIndexDomain].structureName;
				chainInStr.ptr_chainInStr[0].classId = ptr_cathDomain[hitIndexDomain].classId;
				chainInStr.ptr_chainInStr[0].classIdConv = ptr_cathDomain[hitIndexDomain].classIdConv;
			}
		}

		ptr_fileIn = fopen(ptr_dirList[fileNr], "r");


		if (!ptr_fileIn){
			printf("Sorry, the file could not be found\n");
			//getchar();
			return 0;
		}

		/*reinit chain info*/
		returnVal = reinit_chainInStr(chainInStr, maxNrOfChains);


		/*get number of chains and their lengths for this file:*/
		//chainInStr = readPDBChainStructure(ptr_fileIn);
		if(use_scop_b ==1){
			readPDBChainStructureSCOP(ptr_fileIn, &chainInStr );
		}
		else if(use_cath_b ==1){
			readPDBDomainStructureCATH(ptr_fileIn, &chainInStr ); //will read chain length and chain id into chainInStr
		}
		else{
			readPDBChainStructure2(ptr_fileIn, &chainInStr );
		}

		/****************************************************************************************************************/
		/* INNER LOOP */
		/****************************************************************************************************************/

		/*loop through the chains in the current structure/file ... and compute ... :*/
		chainNr = 0; /*reset*/
		for (chainNr = 0; chainNr <= chainInStr.numberOfChains -1 ; chainNr ++){

			/*printf("Now handling: Chain: %s of %s of length %d\n", chainInStr.ptr_chainInStr[chainNr].chainId, ptr_dirList[fileNr], chainInStr.ptr_chainInStr[chainNr].chainLength);
			getchar();*/


			/*We only handle structures of length at most maxChainLength*/
			if(chainInStr.ptr_chainInStr[chainNr].chainLength > maxChainLength){
				printf("Chain: %s of %s has length %d and was skipped\n", chainInStr.ptr_chainInStr[chainNr].chainId, ptr_dirList[fileNr], chainInStr.ptr_chainInStr[chainNr].chainLength);
				chainsSkippedTooLong += 1;
				continue;
			}

			startChain = clock();

			I_measures.fileName = ptr_dirList[fileNr];
			I_measures.structureName = chainInStr.ptr_chainInStr[chainNr].structureName;
			I_measures.classId = chainInStr.ptr_chainInStr[chainNr].classId;
			I_measures.chainId = chainInStr.ptr_chainInStr[chainNr].chainId;
			I_measures.chainNr = &chainNr;

			chainLen = chainInStr.ptr_chainInStr[chainNr].chainLength;

			if (print_basic_b ==1){
				printf("sub str cnt:%d\n",subStructureCnt);
				printf("Chain id:%s\n", chainInStr.ptr_chainInStr[chainNr].chainId);
				printf("Chain nr:%d\n", chainNr);
				printf("Chain Length: %d\n", chainLen);
			}

			/*if the chain length is less than strLengthCutOff we skip the computation*/
			if(chainLen < strLengthCutOff){
				//if (print_basic_b ==1){
				printf("Chain %s of %s skipped since length is < %d\n", I_measures.chainId, I_measures.fileName,  strLengthCutOff);
				//}
				chainsSkippedTooShort += 1;
				continue;
			}

			rewind(ptr_fileIn);
			returnVal = main_readPDB2(ptr_fileIn, ptr_chain, chainNr, chainLen);

			/*length of segments chain:*/
			L = chainLen - 1;
			/*allocate memory for array of segments*/
			//ptr_segment = (struct  segment *) calloc (L, sizeof(struct segment));
			//ptr_segment = (struct  segment *) malloc (L*sizeof(struct segment));
			//ptr_segment = realloc(ptr_segment,L*sizeof(struct segment));


			/*populate ptr_segment */
			for ( n=0; n <= L -1 ; n++ ){
				segCoords.s1 = ptr_chain[n].coords;
				segCoords.s2 = ptr_chain[n+1].coords;

				ptr_segment[n] = segCoords;

			}

			/*We skip the chain if it contains a segment which is longer than the globally set
			 stdRealSegLength:*/
			segTooLong_b = 0;
			for ( n=0; n <= L -1 ; n++ ){
				if (distCalphas(ptr_segment[n].s1, ptr_segment[n].s2) > stdRealSegLength){
					segTooLong_b = n+1;
					break;
				}
			} 
			if(segTooLong_b >= 1 ){
				chainSkippedSegmTooLong +=1;
				printf("Chain: %s in file: %s removed due to too long %d'th segment\n", I_measures.chainId, ptr_dirList[fileNr],  segTooLong_b-1);
				continue;
			}


			if(print_b == 1){
				printf ("Segments: ");
				for ( n=0; n<5; n++ ){
					printf ("s1, x: %lf y:%lf z:%lf \n",ptr_segment[n].s1.x, ptr_segment[n].s1.y, ptr_segment[n].s1.z);
					printf ("s2, x: %lf y:%lf z:%lf \n",ptr_segment[n].s2.x, ptr_segment[n].s2.y, ptr_segment[n].s2.z);
				}
			}

			I_measures.chainLen = &chainLen;

			/*(RE)INITIALIZATION*/
			returnVal = init_I_measures(I_measures, order, full_b, chainLen);

			I_measures.order = &order;

			startComp = clock();

			/*call omputation of the invariants*/
			if (incl_abs_b == 1 && split_b == 0 && closed_loops_b == 0){
				returnVal = aggrAndW(ptr_segment, chainLen, order, full_b, I_measures);
				}
			if (incl_abs_b == 0 && split_b == 0 && closed_loops_b == 0){
				returnVal = aggrAndW_ExAbs(ptr_segment, chainLen, order, full_b, I_measures);
				}
			if (incl_abs_b == 1 && split_b == 1 && closed_loops_b == 0){
				returnVal = wAll(ptr_segment, chainLen, I_measures.wVal);
				returnVal = aggr(chainLen, order, full_b, I_measures.wVal, I_measures);
			}


			/*If desired compute and record the mutual writhe of pairs of closed-loops and 
			pairs (closed loops, segment) at set lengths. Aim is to search for particular geometries 
			(links and "pokes")*/
			if (closed_loops_b == 1){
				returnVal = aggrAndW_wClosedLoops(ptr_segment, chainLen, order, full_b, I_measures, ptr_closedLoopInd, closedLoopLength, closedLoopDist, &nrOfClosedLoops);

				/*examine the closed loops for links and "pokes" and store/write out the results:*/
				returnVal = examineClosedLoops(ptr_fileNameOutClosedLoops, ptr_closedLoopInd, nrOfClosedLoops, I_measures, chainLen, ptr_segment, ptr_fileNameOutChains, ptr_chain, pokeLength, -1, -1);

			}

			//printf("Nr of closed loops found in this structure: %d\n", nrOfClosedLoops);


			/*If desired compute and record the mutual writhe/invariant value of pairs of sub-chains of 
			a set length. Aim is a more general search for particular geometries (e.g. links)*/
			if (invValSubChainPairs_b != 0){

				if(chainLen > 2*subChainLength + 1){

					if (invValSubChainPairs_b == 1){ /*compute only mutual writhe*/ 
					//printf("(langmodigt venter bysvalen stadig stædigt længes ... ) sei stesso qui ..?%d\n", subChainLength);
					//getchar();
						returnVal = writheSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, subChainLength, ptr_fileNameOutChains, ptr_chain, writeSubChainPairsAll_b);
					} 

					if (invValSubChainPairs_b == 2){ /*compute only mutual invariant values for order 1 and order 2; this takes setting the order parameter to 3 or the full_b parameter to 1*/ 

						if (order > 2 || full_b ==1){

							/*genInvariantSubChainPairs(char *ptr_fileNameOut, struct segment *ptr_segment, int chainLength, struct I_ptr I_measures, double **I_measure, char **I_measure_name, int subChainLength, char *ptr_fileNameOutChains, struct cAlpha *ptr_chain, int writeAll_b)*/

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I12, "I12", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I1234_full, "I1234", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I1324_full, "I1324", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I1423 , "I1423", subChainLength, writeSubChainPairsAll_b);

							/*write out the chain:*/
							returnVal = writeChain(ptr_fileNameOutChains, I_measures.fileName,I_measures.structureName, I_measures.chainId, ptr_chain, *I_measures.chainLen);

						}
						else{
							printf("Set full_b to 1 or order = 3 to do this!\n");
						}

					}

					if (invValSubChainPairs_b == 3){ /*compute mutual invariant values for all orders; takes setting order = 3*/ 

						if (order == 3){
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I12, "I12", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.Ia12, "Ia12", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I1234, "I1234_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.Ia1234, "Ia1234_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I12a34, "I12a34_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.Ia12a34, "Ia12a34_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I1234_full, "I1234", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I1324, "I1324_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.Ia1324, "Ia1324_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I13a24, "I13a24_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.Ia13a24, "Ia13a24_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I1324_full, "I1324", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I1423, "I1423_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.Ia1423, "Ia1423_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I14a23, "I14a23_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.Ia14a23, "Ia14a23_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I1423, "I1423", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I123456, "I123456", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I123546, "I123546", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I123645, "I123645", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I132456, "I132456", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I132546, "I132546", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I132645, "I132645", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I142356, "I142356", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I142536, "I142536", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I142635, "I142635", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I152346, "I152346", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I152436, "I152436", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I152634, "I152634", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I162345, "I162345", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I162435, "I162435", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I162534, "I162534", subChainLength, writeSubChainPairsAll_b);

							/*write out the chain:*/
							returnVal = writeChain(ptr_fileNameOutChains, I_measures.fileName,I_measures.structureName, I_measures.chainId, ptr_chain, *I_measures.chainLen);

						}
						else{
							printf("Set order = 3 to do this\n");
						}
					}


				}

				else{
					printf("Twice the SubChainLength in search longer than chain length; file: %s, chain id: %s\n", I_measures.fileName, I_measures.chainId);
				}

			}


			endComp = clock();

			compTime = ((double) (endComp - startComp))/CLOCKS_PER_SEC;
			//compTime = ((long double) (endComp - startComp)/QueryPerformanceFrequency);
			if(print_b == 1){
				printf("Comp time: %lf\n",compTime);
				//printf("CPU time spend for measures on this file (only computation/aggregation):%lf \n", compTime);
			}
			if (compTime > max_compTime){
				max_compTime = compTime;
			}
			compTimeComplete += compTime;

			/*write all I-values to file if wanted:*/
			if (write_all_b == 1){

				returnVal = writeAllIvaluesToFile(ptr_fileNameOutAll, full_b, I_measures, pertNo);

			}


			/*look at the w-values etc:*/
			if (print_b == 1){

				/*for ( i=0; i<=L-1; i++ ){
					for ( j=0; j<=L-1; j++ ){
						printf("w[%d][%d]:%lf ", i,j, I_measures.wVal[i][j]);		
					}
				}*/

				if (order >= 1){
					printf("I12[%d][%d]:%lf ", 0,L-1, I_measures.I12[0][L-1]);
					printf("Ia12[%d][%d]:%lf ", 0,L-1, I_measures.Ia12[0][L-1]);
					/*for ( i=0; i<=L-1; i++ ){
						for ( j=0; j<=L-1; j++ ){
							printf("I12[%d][%d]:%lf ", i,j, I_measures.I12[i][j]);
						}
					}*/
				}

				if (order >= 2){
					printf("I1234[%d][%d]:%lf ", 0,L-1, I_measures.I1234[0][L-1]);
					/*for ( i=0; i<=L-1; i++ ){
						for ( j=0; j<=L-1; j++ ){
							printf("I1234[%d][%d]:%lf ", i,j, I_measures.I1234[i][j]);
						}
					}*/
					if (full_b == 1 && order ==2){
						printf("I1234_full[%d][%d]:%lf ", 0,L-1, I_measures.I1234_full[0][L-1]);
						printf("I1324_full[%d][%d]:%lf ", 0,L-1, I_measures.I1324_full[0][L-1]);
					}
				}

				if (order >= 3){
					printf("I123456[%d][%d]:%lf ", 0,L-1, I_measures.I123456[0][L-1]);
					printf("I142356[%d][%d]:%lf ", 0,L-1, I_measures.I142356[0][L-1]);
					printf("I132645[%d][%d]:%lf ", 0,L-1, I_measures.I132645[0][L-1]);
					printf("I132546[%d][%d]:%lf ", 0,L-1, I_measures.I132546[0][L-1]);
					printf("I142536[%d][%d]:%lf ", 0,L-1, I_measures.I142536[0][L-1]);
					printf("I142635[%d][%d]:%lf ", 0,L-1, I_measures.I142635[0][L-1]);
					printf("I1423_full0[%d][%d]:%lf ", 0,L-1, I_measures.I1423_full0[0][L-1]);
					printf("I1423_full2[%d][%d]:%lf ", 0,L-1, I_measures.I1423_full2[0][L-1]);
					printf("I1324_full2[%d][%d]:%lf ", 0,L-1, I_measures.I1324_full2[0][L-1]);
					printf("I1234_full[%d][%d]:%lf ", 0,L-1, I_measures.I1234_full[0][L-1]);
					printf("I1324_full[%d][%d]:%lf ", 0,L-1, I_measures.I1324_full[0][L-1]);

				}
			}

			endChain = clock();
			timeChain = ((double) (endChain - startChain)) / CLOCKS_PER_SEC;
			//printf("CPU time spend for measures on this file incl load of chain and more:%lf \n", timeChain);
			if (timeChain > max_timeChain){
				max_timeChain =timeChain;
			}

			/*collect final results for write-out if that's desired*/
			if (write_final_b == 1){
				/*reads I-values at corner 0,L-1 into the ptr_I_values ptr at index = subStructureCnt; pertNo is preset to 0:*/
				returnVal = collectIvalues(ptr_I_values, I_measures, subStructureCnt, pertNo, timeChain, compTime);
		
			}

			subStructureCnt += 1; /*will count and index the sub-structures for which we compute the invariants*/
		
		} /*end chainNr loop*/

		
		fclose(ptr_fileIn);

	} /*end fileNr loop*/


	/*free the globally allocated memory*/
	free_I_measures(I_measures, order, full_b, maxChainLen);
	free(ptr_segment);	
	free(ptr_chain);
	free(dirContent.ptr_dirList); //this is allocated in the ListDirectoryContents2 fct 

	if (write_final_b == 1){

		/*generate file name for the output*/ 
		sprintf(orderStr, "%d", order);
		strcat(fileNameOut, outputPath);
		strcat(fileNameOut, fileNameOut1);
		strcat(fileNameOut, orderStr);
		strcat(fileNameOut, fileNameOut2);
		strcat(fileNameOut, fileNameOut3);
		strcat(fileNameOut, DBName);
		strcat(fileNameOut, extensionName);
		printf("Results for chain (only) will be written to: %s\n", fileNameOut);

		ptr_fileNameOut = fileNameOut;

		/*clear the file by opening it in w-mode and then closing it again*/
	    fclose(fopen(ptr_fileNameOut, "w"));

		/*call writing the output (obs: numberOfPerturbations is preset to 1)*/
		returnVal = writeIvaluesToFile(ptr_fileNameOut, ptr_I_values, subStructureCnt, numberOfPerturbations);

		/*free the memory used for the purpose:*/
		free(ptr_I_values);
		 
	}

	endComplete = clock();
	//printf("start clocks:%d, end clocks: %d, clocks_per sec:%d", endComplete, startComplete, CLOCKS_PER_SEC);
	timeComplete = ((double) (endComplete - startComplete))/CLOCKS_PER_SEC;
	//compTimeComplete = ((double) (compTimeComplete)) / CLOCKS_PER_SEC;
	printf("CPU time spend for all files incl loading data:%lf \n", timeComplete);
	printf("CPU time spend for computations on all files:%lf \n", compTimeComplete);
	printf("Max CPU time spend for measures across all files:%lf \n", max_timeChain);

	printf("Done number of files:%d \n", fileNr);
	printf("Done number of chains/sub-structures:%d \n", subStructureCnt);
	printf("%d too short and %d too long chains were skipped\n", chainsSkippedTooShort,chainsSkippedTooLong);
	printf("%d chains were skipped due to a too long segment (longer than sqrt of %d)\n",chainSkippedSegmTooLong, stdRealSegLength);


	//getchar(); 
	
	returnVal = 1;

	return returnVal;
}



int computeGI_windows(char DBName[100], char *dirPath, int use_scop_b, int use_cath_b, char *cathListPath, 
					  int loadFromSubDirs_b, char *outputPath, int maxChainLength, int omitSingletons_b, int order, 
					  int incl_abs_b, int full_b, int split_b, int write_final_b, int write_all_b, 
					  int closed_loops_b, int closedLoopLength, int closedLoopDist, int pokeLength, int invValSubChainPairs_b, int writeSubChainPairsAll_b,
					  int subChainLength, int windowCoveringType, int get_windows_b, int get_windowPairs_b, int windowLength,
					  int stepSize, int write_windows_b, int write_windowPairs_b, int writeOnlyDisjointPairs_b, int maxNrOfWindows, 
					  int write_chains_b, char fileNameChains[2000], int print_b, int print_basic_b){

	int returnval = 0;

	//int strLengthCutOff = 10; /*chains shorter than this int will not be included in computation of invariants*/

	struct dirContent dirContent;
	int subDirCnt = 0;
	int numberOfFiles = 0;
	int subStructureCnt = 0; /*some structures are multi-mers; we compute the invariants for each sub-structure (each mono-mer)*/
	int nrOfSingletons = 0;
	const int numberOfPerturbations = 1; /*only there since needed when writing out the final values (preparation for perturbations code)*/
	const int pertNo = 0; /*ditto*/
	char ** ptr_dirList;
	char ** ptr_fileNameList;

	/*for "dynamically" generating the file names for output:*/
	char orderStr[10] = "-1"; //value is placeholder
	char inclAbsStr[10] = "-1"; //value is placeholder
	char windowLengthStr[10] = "-1"; //value is placeholder
	char windowStepSizeStr[10] = "-1"; //value is placeholder
	char windowCovTypeStr[10] = "-1"; //value is placeholder
	char writeOnlyDisjointPairsStr[10] = "-1"; //value is placeholder


	char fileNameOut1[1000] = "//Ivalues_order_";
	char fileNameOutAll1[1000] = "//All_Ivalues_order_";
	char fileNameOutClosedLoops1[1000] = "//ClosedLoopChars"; 
	char fileNameOutSubChainPairs1[1000] = "//SubChainPairChars"; 
	char fileNameOutWindows1[1000] = "//Invariants_windowlgth_";
	char fileNameOutwindowPairs1[1000] = "//Invariants_Pairs_windowlgth_";
	char fileNameOut2[1000] = "_order_"; //; //; //"_minmax_";//"_incl_abs_"; //"_"; //

	//char fileNameOut3[1000] = "_mem_alloc_top100.txt"; 
	//char fileNameOutChains1[1000] = "\\chains_closed_loops_top8000.txt";
	char fileNameOutChains1[1000] =  "//chains_closed_loops_";
	char fileNameOutChains2[1000] =  "/chains_subchain_pairs_";
	char fileNameOut3[1000] = "computeGI_windows_"; 
	char extensionName[30] = ".txt";

	char fileNameOutAll[2000] = "\0";
	char fileNameOut[2000] = "\0";
	char fileNameOutClosedLoops[2000] = "\0";
	char fileNameOutSubChainPairs[2000] = "\0";
	char fileNameOutWindows[2000] = "\0";
	char fileNameOutwindowPairs[2000] = "\0";
	char fileNameOutChains[2000] = "\0";

	char *ptr_fileNameOut; /*pointer to file for holding final measure values (i.e. at final corner of simplex)*/ 
	char *ptr_fileNameOutAll; /*pointer to file for holding all measure values (i.e. across whole simplex)*/ 
	char *ptr_fileNameOutClosedLoops; /*pointer to file for holding mutual writhe value of closed loop pairs*/ 
	char *ptr_fileNameOutSubChainPairs; /*pointer to file for holding mutual writhe value of sub-chain pairs*/ 
	char *ptr_fileNameOutWindows; /*pointer to file for holding measures' values on sub-chains of pre-set windowlength*/
	FILE *ptr_fileOutWindows;
	char *ptr_fileNameOutwindowPairs; /*pointer to file for holding measures' values on pairs of sub-chains of pre-set windowlength*/
	FILE *ptr_fileOutwindowPairs;
	char *ptr_fileNameOutChains; /*pointer to file for holding chain coordinates*/

	char closedLoopLengthStr[10] = "\0"; /*for holding closedLoopLength converted to a string*/
	char pokeLengthStr[10] = "\0"; /*for holding pokeLength converted to a string*/
	char closedLoopDistStr[10] = "\0"; /*for holding closedLoopdist converted to a string*/
	char subChainLengthStr[10] = "\0"; /*for holding subChainLength converted to a string*/

	clock_t startChain, endChain, startComp, endComp, startComplete, endComplete;
	double timeChain = 0, compTime = 0, compTimeComplete = 0, timeComplete = 0, max_compTime = 0, max_timeChain = 0;

	FILE *ptr_fileIn;
		
	FILE *ptr_cath_domain_list_file;
	int convToCLF20_b = 0;
	struct chainInfo *ptr_cathDomain;
	int nrOfCATHDomainsInList;
	int hitIndexDomain = -1;
	int lastHitIndexDomain = -1;

	struct chainInStructure chainInStr;
	int structureNameLength = 7;
	char *currentStructureName;
	char *lastStructureName;
	char *structureNameInit = "NN";
	char *currentClassId;
	char *classIdInit = "NN";
	char chainIdInit[2] = ">";

	struct classInfo *ptr_classInfo;
	int maxNrOfStructures;
	int *ptr_range = (int *) malloc(2*sizeof(int));
	int hitIndex;
	int classSize = 0;
	int pos = 0; 


	int chainsSkippedTooShort = 0;
	int chainsSkippedTooLong = 0;
	int chainSkippedSegmTooLong = 0;
	int segTooLong_b = 0;
	int resNr = 0; /*residue number*/
	int resNrPrior = -999;
	int chainLen = 0; /*length of C-alpha chain*/
	int maxChainLen =0; /*max length of chain among all loaded*/
	int L = 0; /*length of segment chain = chainLen -1*/
	int simplexCnt = 0;/*size of simplex*/
	int maxSimplexCnt = 0; /*size of simplex corr to max chain length*/

	struct cAlpha * ptr_chain = NULL; /*for holding C-alpha chain */
	struct segment *ptr_segment = NULL; /* for holding chain of segments, ith C-alpha to jth C-alpha*/
	struct segment segCoords;	 

	int n = 0;
	int i = 0;
	int j = 0;
	int cnt = 0;
	int fileNr = 0;
	int chainNr = 0;

	/* for holding meausure values (i.e. the core of the final output); 
	content to be received from aggr function:*/
	struct I_ptr I_measures; 

	int returnVal = 0;

	/*pointer to hold the results (that this is a ptr of a ptr is due to that the write-out can also
	handle the perturbation case; we use pertNo = 0 here):*/
	struct I_values **ptr_I_values;

	/*struct for holding invariants on windows resp. window pairs:*/
	struct I_windows_ptr I_windows;
	struct I_windowPairs_ptr I_windowPairs;

	/*for closed loops finding:*/
	struct twoSegmentIndex *ptr_closedLoopInd;
	int initValNrOfClosedLoops = 0;
	int nrOfClosedLoops = 0;
	int cntClosedLoop = 0;

	FILE *fp; //dummy file ptr used when checking file names

	/*declaration ends here*/

	/*General assertions*/
	if (order == 2 && full_b == 0 && (get_windows_b == 1 || get_windowPairs_b == 1)){

		printf("When invariants on windows or window pairs are required and order is 2, full_b must be 1!\n");
		printf("To continue full_b is set to 1; press anything to continue\n");
		full_b = 1;
		getchar();
	}

	/****************************************************************************************************************/
	/*Generation of output file names:*/
	/****************************************************************************************************************/

	if (write_all_b == 1){

		sprintf(orderStr, "%d", order);
		strcat(fileNameOutAll, outputPath);
		strcat(fileNameOutAll, fileNameOutAll1);
		strcat(fileNameOutAll, orderStr);
		sprintf(inclAbsStr, "%d", incl_abs_b);
		strcat(fileNameOutAll, inclAbsStr);
		strcat(fileNameOutAll, "_");
		strcat(fileNameOutAll, fileNameOut2);
		strcat(fileNameOutAll, fileNameOut3);
		if (DBName != ""){
			strcat(fileNameOutAll, DBName);
		}
		strcat(fileNameOutAll, extensionName);
		printf("Results for all vertices will be written to: %s\n", fileNameOutAll);

		ptr_fileNameOutAll = fileNameOutAll;
		
		/*clear the file by opening it in w-mode and then closing it again*/
	    fp = fopen(ptr_fileNameOutAll, "w");
		if (fp == NULL) {
			printf("Sorry, file name %s is invalid. Could the path be wrong? A blank space?", ptr_fileNameOutAll);
			getchar();
			exit(EXIT_FAILURE);
		}
		else {
			fclose(fp);
		}

	}

	if (closed_loops_b == 1){

		/*generate file name for writing out closed loops/poke details:*/
		sprintf(pokeLengthStr,"%d",pokeLength);
		sprintf(closedLoopLengthStr,"%d",closedLoopLength);

		strcat(fileNameOutClosedLoops, outputPath);
		strcat(fileNameOutClosedLoops, fileNameOutClosedLoops1);

		strcat(fileNameOutClosedLoops, fileNameOut2);
		strcat(fileNameOutClosedLoops, fileNameOut3);

		if (DBName != ""){
			strcat(fileNameOutClosedLoops, DBName);
			strcat(fileNameOutClosedLoops, "_");
		}

		strcat(fileNameOutClosedLoops, "closedLoopLength_");
		sprintf(closedLoopLengthStr, "%d", closedLoopLength);
		strcat(fileNameOutClosedLoops, closedLoopLengthStr);

		strcat(fileNameOutClosedLoops, "_closedLoopDist_");
		sprintf(closedLoopDistStr, "%0.1f", sqrt(closedLoopDist));
		strcat(fileNameOutClosedLoops, closedLoopDistStr);
		strcat(fileNameOutClosedLoops, "_pokeLength");
		strcat(fileNameOutClosedLoops, pokeLengthStr);
		
		strcat(fileNameOutClosedLoops, extensionName);
		printf("Results for closed loops will be written to: %s\n", fileNameOutClosedLoops);

		ptr_fileNameOutClosedLoops = fileNameOutClosedLoops;
		
		/*clear the file by opening it in w-mode and then closing it again*/
	    fclose(fopen(ptr_fileNameOutClosedLoops, "w"));

		/*generate file name for writing out chains:*/
		strcat(fileNameOutChains, outputPath);
		strcat(fileNameOutChains, fileNameOutChains1);
		if (DBName != ""){
			strcat(fileNameOutChains, DBName);
		}
		strcat(fileNameOutChains, extensionName);		
	
		ptr_fileNameOutChains = fileNameOutChains;
		
		/*clear the file by opening it in w-mode and then closing it again*/
	    fclose(fopen(ptr_fileNameOutChains, "w"));

	}
	
	if (invValSubChainPairs_b != 0 ){

		/*generate file name for writing out writhe results for sub-chain pairs:*/
		strcat(fileNameOutSubChainPairs, outputPath);
		strcat(fileNameOutSubChainPairs, fileNameOutSubChainPairs1);
		strcat(fileNameOutSubChainPairs, fileNameOut2);
		strcat(fileNameOutSubChainPairs, fileNameOut3);

		if (DBName != ""){
			strcat(fileNameOutSubChainPairs, DBName);
			strcat(fileNameOutSubChainPairs, "_");
		}

		strcat(fileNameOutSubChainPairs, "subChainLength_");
		sprintf(subChainLengthStr, "%d", subChainLength);
		strcat(fileNameOutSubChainPairs, subChainLengthStr);

		strcat(fileNameOutSubChainPairs, extensionName);
		printf("Results for sub-chain pairs will be written to: %s\n", fileNameOutSubChainPairs);

		ptr_fileNameOutSubChainPairs = fileNameOutSubChainPairs;
		
		/*clear the file by opening it in w-mode and then closing it again*/
	    fclose(fopen(ptr_fileNameOutSubChainPairs, "w"));

		/*generate file name for writing out chains:*/
		strcat(fileNameOutChains, outputPath);
		strcat(fileNameOutChains, fileNameOutChains2);
		if (DBName != ""){
			strcat(fileNameOutChains, DBName);
		}
		strcat(fileNameOutChains, extensionName);	
		
		ptr_fileNameOutChains = fileNameOutChains;
		
		/*clear the file by opening it in w-mode and then closing it again*/
	    fclose(fopen(ptr_fileNameOutChains, "w"));

	}

	
	if (write_windows_b != 0){

		/*generate file name for writing out invariants' values on subchains*/
		sprintf(orderStr, "%d", order);
		strcat(fileNameOutWindows, outputPath);
		strcat(fileNameOutWindows, fileNameOutWindows1);
		sprintf(windowLengthStr, "%d", windowLength);
		strcat(fileNameOutWindows, windowLengthStr);
		strcat(fileNameOutWindows, "_");
		sprintf(windowStepSizeStr, "%d", stepSize);
		strcat(fileNameOutWindows, windowStepSizeStr);
		strcat(fileNameOutWindows, fileNameOut2);
		strcat(fileNameOutWindows, orderStr);
		strcat(fileNameOutWindows, "_");
		sprintf(inclAbsStr, "%d", incl_abs_b);
		strcat(fileNameOutWindows, inclAbsStr);
		strcat(fileNameOutWindows, "_");
		strcat(fileNameOutWindows, fileNameOut3);
		if (DBName != ""){
			strcat(fileNameOutWindows, DBName);
			strcat(fileNameOutWindows, "_");
		}
		strcat(fileNameOutWindows,"winCovType");
		sprintf(windowCovTypeStr, "%d", windowCoveringType);
		strcat(fileNameOutWindows, windowCovTypeStr);
		strcat(fileNameOutWindows, extensionName);
		printf("Results for windows will be written to: %s\n", fileNameOutWindows);

		 
		ptr_fileNameOutWindows = fileNameOutWindows;
		
		/*clear the file by opening it in w-mode and then closing it again*/
	    fp = fopen(ptr_fileNameOutWindows, "w");
		if (fp == NULL) {
			printf("Sorry, file name %s is invalid. Could the path be wrong? A blank space?", ptr_fileNameOutWindows);
			getchar();
			exit(EXIT_FAILURE);
		}
		else {
			fclose(fp);
		}

	}

	if (write_windowPairs_b != 0){

		/*generate file name for writing out invariants' values on subchains*/
		sprintf(orderStr, "%d", order);
		strcat(fileNameOutwindowPairs, outputPath);
		strcat(fileNameOutwindowPairs, fileNameOutwindowPairs1);
		sprintf(windowLengthStr, "%d", windowLength);
		strcat(fileNameOutwindowPairs, windowLengthStr);
		strcat(fileNameOutwindowPairs, "_");
		sprintf(windowStepSizeStr, "%d", stepSize);
		strcat(fileNameOutwindowPairs, windowStepSizeStr);
		strcat(fileNameOutwindowPairs, fileNameOut2);
		strcat(fileNameOutwindowPairs, orderStr);
		strcat(fileNameOutwindowPairs, "_");
		sprintf(inclAbsStr, "%d", incl_abs_b);
		strcat(fileNameOutwindowPairs, inclAbsStr);
		strcat(fileNameOutwindowPairs, "_");
		strcat(fileNameOutwindowPairs, fileNameOut3);
		if (DBName != ""){
			strcat(fileNameOutwindowPairs, DBName);
			strcat(fileNameOutwindowPairs, "_");
		}
		
		strcat(fileNameOutwindowPairs, "winCovType");
		sprintf(windowCovTypeStr, "%d", windowCoveringType);
		strcat(fileNameOutwindowPairs, windowCovTypeStr);
		strcat(fileNameOutwindowPairs, "_onlyDisjointPairs");
		sprintf(writeOnlyDisjointPairsStr, "%d", writeOnlyDisjointPairs_b);
		strcat(fileNameOutwindowPairs, writeOnlyDisjointPairsStr);
		strcat(fileNameOutwindowPairs, extensionName);
		printf("Results for window pairs will be written to: %s\n", fileNameOutwindowPairs);

		 
		ptr_fileNameOutwindowPairs = fileNameOutwindowPairs;
		
		/*clear the file by opening it in w-mode and then closing it again*/
	    fp = fopen(ptr_fileNameOutwindowPairs, "w");
		if (fp == NULL) {
			printf("Sorry, file name %s is invalid. Could the path be wrong? A blank space?", ptr_fileNameOutwindowPairs);
			getchar();
			exit(EXIT_FAILURE);
		}
		else {
			fclose(fp);
		}

	}

		
	/****************************************************************************************************************/
	/*Get content of the desired directory: number of files and a list of the file names:*/
	/****************************************************************************************************************/

	if(loadFromSubDirs_b ==0){
		dirContent = ListDirectoryContents(dirPath);
	}
	else{
		returnVal = ListDirectoryContents2(dirPath, &dirContent, &subDirCnt);
	}

	numberOfFiles = dirContent.numberOfFiles; 
	ptr_dirList = dirContent.ptr_dirList;
	ptr_fileNameList = dirContent.ptr_fileNameList;

	printf("Number of files in directory:%d \n",numberOfFiles);
	printf("First file: %s\n", ptr_fileNameList[0]);
	//getchar();

	startComplete = clock();

	/****************************************************************************************************************/
	/*loop through the list of filenames in the directory to find chain info, eg the max chain length (for setting
	the size in the (global) mem alloc to I_measures's ptr's). First though we make a bulk allocation
	to the chain-info keeping pointer, using the pre-set maxNrOfChains. We also allocate to a 
	ptr keeping the class info for every chain/structure:*/
	/****************************************************************************************************************/

	chainInStr.numberOfChains = maxNrOfChains;
	chainInStr.ptr_chainInStr = (struct chainInfo *) calloc (maxNrOfChains, sizeof(struct chainInfo));
	/*init*/

	for(i = 0; i< maxNrOfChains; i++){

		chainInStr.ptr_chainInStr[i].chainId = (char *) malloc(sizeof(char)*2);
		chainInStr.ptr_chainInStr[i].structureName = (char *) malloc(sizeof(char)*10);
		chainInStr.ptr_chainInStr[i].classId = (char *) malloc(sizeof(char)*classIdCharSize);
		chainInStr.ptr_chainInStr[i].classIdConv = (int *) malloc(classifierSize*sizeof(int));

		strcpy(chainInStr.ptr_chainInStr[i].chainId, chainIdInit); /*init to unlikely value*/
		strcpy(chainInStr.ptr_chainInStr[i].structureName, structureNameInit);
		strcpy(chainInStr.ptr_chainInStr[i].classId, classIdInit);
		//printf("Chain Id initial value: %c\n", *ptr_chainInStr[i].chainId);
		for(j=0; j< classifierSize; j++){chainInStr.ptr_chainInStr[i].classIdConv[j] = -1;}

	}

	/*Fork: whether we use data from CATH or don't. In CATH data domain info (name, class id, chain length) is
	stored in file separately from the file containing the PDB 3d-coord data. In SCOP and Kinemage data
	the PDB-files contain the info (either directly or, else, extracted by the code).*/

	/*if CATH data, read in CATH domain list. Note: we abuse the chainInStrcture struct to accomodate
	domains rather than chains -- in CATH each chain may be split in several domains*/ 
	if(use_cath_b == 1){

		if(convToCLF20_b == 1){
			structureNameLength =  lengthCLFformat20;
		}
		else{
			structureNameLength =  6;
		}

		/*first allocate a block of memory (if not large enough it will be extended when exectuing readCATHDomainList)*/
		ptr_cathDomain = (struct chainInfo *) malloc (maxNrOfDomains*sizeof(struct chainInfo));
		/*init*/
		for(i = 0; i< maxNrOfDomains; i++){
			ptr_cathDomain[i].chainId = (char *) malloc(sizeof(char)*2);
			ptr_cathDomain[i].structureName = (char *) malloc(sizeof(char)*(structureNameLength +2)); //to hold domain name
			ptr_cathDomain[i].classId = (char *) malloc(sizeof(char)*classIdCharSize);
			ptr_cathDomain[i].classIdConv = (int *) malloc(classifierSize*sizeof(int));

			strcpy(ptr_cathDomain[i].chainId, chainIdInit); /*init to unlikely value*/
			strcpy(ptr_cathDomain[i].structureName, structureNameInit);
			strcpy(ptr_cathDomain[i].classId, classIdInit);
			/*printf("Chain Id initial value:%c\n", *ptr_cathDomain[i].chainId);
			printf("Class Id initial value:%c\n", *ptr_cathDomain[i].classId);*/
			for(j=0; j< classifierSize; j++){ptr_cathDomain[i].classIdConv[j] = -1;}
			ptr_cathDomain[i].chainLength = -1;

		}

		//getchar(); 

		/*open CATH file list:*/
		ptr_cath_domain_list_file = fopen(cathListPath, "r");

		/*read in contents of CATH domain list file*/
		nrOfCATHDomainsInList = readCATHDomainList(ptr_cath_domain_list_file, &ptr_cathDomain, convToCLF20_b );
		//getchar();


		printf("nr of cath doms: %d\n", nrOfCATHDomainsInList);

		/*lex-sort the ptr_cath_domain_list_file on the domain name entry. Below we
		need to look up each file for which we have PDB-coords in this list. In other
		words, this list serves as a "positive list" of domains, for which to compute
		the invariants.*/
		if(convToCLF20_b == 1){
			heapSortCATHDomain(ptr_cathDomain, nrOfCATHDomainsInList, structureNameLength,0);
		}
		else{
			heapSortCATHDomain(ptr_cathDomain, nrOfCATHDomainsInList, 6,0);
		}

		//for(i=0;i <nrOfCATHDomainsInList;i++){printf("dom %d: %s\n", i, ptr_cathDomain[i].structureName);}
		//getchar();
	}

	/*the ptr to hold the class info; we initialize it in the loop over fileNr right below.
	We want an array allowing to be sorted directly on the classId; thereafter the size
	of each class (classSize) can be easily had:*/
	if(use_cath_b == 1){
		maxNrOfStructures = nrOfCATHDomainsInList*maxNrOfChains;
	}
	else{
		maxNrOfStructures = numberOfFiles*maxNrOfChains;
	}

	ptr_classInfo = (struct classInfo *) malloc(maxNrOfStructures*sizeof(struct classInfo));
	
	for(subStructureCnt=0; subStructureCnt< maxNrOfStructures;subStructureCnt++){

		ptr_classInfo[subStructureCnt].classId = (char *) malloc(sizeof(char)*classIdCharSize);
		ptr_classInfo[subStructureCnt].classIdConv = (int *)malloc(classifierSize*sizeof(int)); 
		//ptr_classInfo[subStructureCnt].classSize = (int)malloc(sizeof(int));

    }

	/*then loop over nr of files; aim: to find the number of structures to handle, subStructureCnt, and 
	the max chain length, maxChainLen:*/
	printf("Fetching chain info for data base structures ...\n");
	subStructureCnt = 0;
	lastStructureName = (char *) malloc(sizeof(char)*(structureNameLength + 2));
	strcpy(lastStructureName,structureNameInit);
	for (fileNr = 0; fileNr <numberOfFiles; fileNr ++){

		if (fileNr%100 == 0){printf(".. now fetched the chain info for %d structures\n", fileNr);}

		if (print_b == 1){
			printf("File name:%s fileNr:%d\n",ptr_dirList[fileNr], fileNr);
		}

		/*if using CATH data: check if the file, ptr_fileNameList[fileNr], is in the
		domain list specified (and loaded and lex-sorted above)*/
		if(use_cath_b == 1){ //if not using CATH, chainInStr.ptr_chainInStr[i].classId was read in in the pdb-file reads above
			hitIndexDomain = bisectionalSearchCATHDomain(ptr_fileNameList[fileNr], ptr_cathDomain, nrOfCATHDomainsInList, structureNameLength);
			if(hitIndexDomain == -1){
				if(print_basic_b > 0){
					printf("File %s not in domain list\n",ptr_fileNameList[fileNr]);
					//getchar();
				}
				continue;
			} //file name was not in the domain list
		}

		ptr_fileIn = fopen(ptr_dirList[fileNr], "r");

		if (!ptr_fileIn){
			printf("Sorry, the file: %s, fileNr: %d could not be read\n", ptr_dirList[fileNr], fileNr);
			//getchar();
			continue;
		}

		/*get number of chains, their id's and lengths for this file:*/
		//chainInStr = readPDBChainStructure(ptr_fileIn);
		if(use_scop_b ==1){
			readPDBChainStructureSCOP(ptr_fileIn, &chainInStr );
		}
		else if(use_cath_b ==1){
			readPDBDomainStructureCATH(ptr_fileIn, &chainInStr );//will read chain length and chain id into chainInStr
			//getchar();
		}
		else{
			readPDBChainStructure2(ptr_fileIn, &chainInStr );
		}

		fclose(ptr_fileIn); 

		for (i=0;i <= chainInStr.numberOfChains -1;i++){
			if(chainInStr.ptr_chainInStr[i].chainLength <= maxChainLength){ //we disregard very long chains for RAM availability reasons
				maxChainLen = max(chainInStr.ptr_chainInStr[i].chainLength,maxChainLen);
				if(chainInStr.ptr_chainInStr[i].chainLength >= strLengthCutOff){
				
					//printf("sub str cnt 0: %d\n", subStructureCnt);

					/*if using CATH data: we checked inside the fileNr-loop if the file, ptr_fileNameList[fileNr], is in the
					domain list specified (and loaded and lex-sorted above). We need another check to avoid loading in 
					several chains for the structureName under consideration:*/
					if(use_cath_b == 1){ //if not using CATH, chainInStr.ptr_chainInStr[i].classId was read in in the pdb-file reads above
						
						/*In case the pdb-file contains several chains: the domain is attached to a single chain (is a subset of 
						a single chain) and we want only to consider that chain. we therefore skip "wrong chains". Obs: hitIndexDomain
						was found above:*/
						if(chainInStr.numberOfChains > 1 && ptr_cathDomain[hitIndexDomain].structureName[4] != chainInStr.ptr_chainInStr[i].chainId[0]){

							if(print_basic_b > 0){
								printf("CATH domain: %s does not sit in this chain: %s\n", ptr_cathDomain[hitIndexDomain].structureName, chainInStr.ptr_chainInStr[i].chainId);
							}
							continue;
						}

						/*update ptr_classInfo with the same data:*/
						strcpy(ptr_classInfo[subStructureCnt].classId, ptr_cathDomain[hitIndexDomain].classId);
						for(j=0; j< classifierSize; j++){ptr_classInfo[subStructureCnt].classIdConv[j] = ptr_cathDomain[hitIndexDomain].classIdConv[j];}

						//printf("unsorted, converted class %d: %d %d %d\n", i, ptr_classInfo[subStructureCnt].classIdConv[0], ptr_classInfo[subStructureCnt].classIdConv[1],ptr_classInfo[subStructureCnt].classIdConv[2]);

					}
					
					/*read available/initial values into ptr_classInfo if the data source is not CATH data
					(if it is CATH we have read all values in above): */ 
					if(use_cath_b != 1){
						strcpy(ptr_classInfo[subStructureCnt].classId, chainInStr.ptr_chainInStr[i].classId);
						for(j=0; j< classifierSize; j++){ptr_classInfo[subStructureCnt].classIdConv[j] = -1;}
					}

					//printf("class id: %s\n", chainInStr.ptr_chainInStr[i].classId);

					ptr_classInfo[subStructureCnt].classSize = -1; /*to be updated to proper number below*/
					/*reallocate some memory if necessary*/
					if(subStructureCnt > maxNrOfStructures - 1){
						maxNrOfStructures += maxNrOfStructures + 100;
						ptr_classInfo = (struct classInfo *) realloc(ptr_classInfo, maxNrOfStructures);
					}

					//printf("sub str cnt 1: %d\n", subStructureCnt);
					
					subStructureCnt += 1; /*count and index the sub-structures for which the invariants will be computed*/

				}

			}
		}

	}

	printf("... done.\n");

	L = maxChainLen -1;

	printf("max chain len:%d\n", maxChainLen);
	printf("sub structure count:%d\n", subStructureCnt);

	//getchar();
	//return 0;

	/*sort the class-info array lexicographically on the classId; then get the number of structures
	in each file and update the ptr_classInfo*/
	if(use_scop_b == 1 || use_cath_b == 1){

		if(use_scop_b == 1){
			/*convert classId's to array's of int's:*/
			//for(i=0; i< 20; i++){printf("class %d: %s\n", i, ptr_classInfo[i].classId);}
			returnVal = convertSCOP(ptr_classInfo, subStructureCnt);
			//for(i=0; i< subStructureCnt; i++){printf("at %d class %s was conv to: %d %d %d %d\n", i, ptr_classInfo[i].classId, ptr_classInfo[i].classIdConv[0], ptr_classInfo[i].classIdConv[1],ptr_classInfo[i].classIdConv[2],ptr_classInfo[i].classIdConv[3]);}
			//for(i=0; i< 30; i++){printf("at %d class %s was conv to: %d %d %d %d\n", i, ptr_classInfo[i].classId, ptr_classInfo[i].classIdConv[0], ptr_classInfo[i].classIdConv[1],ptr_classInfo[i].classIdConv[2],ptr_classInfo[i].classIdConv[3]);}

			/*now sort lexicographically*/
			heapSortDBClass(ptr_classInfo, subStructureCnt, 4, 1); 
			//4: class, fold, super-family, family. Could do with first two
			//for(i=0; i< min(30,subStructureCnt); i++){printf("sorted, converted class %d: %d %d %d %d\n", i, ptr_classInfo[i].classIdConv[0], ptr_classInfo[i].classIdConv[1],ptr_classInfo[i].classIdConv[2],ptr_classInfo[i].classIdConv[3]);}
		}
		else if(use_cath_b == 1){
			/*sort lexicographically*/
			//for(i=0; i< min(30,subStructureCnt); i++){printf("unsorted, converted class %d: %d %d %d\n", i, ptr_classInfo[i].classIdConv[0], ptr_classInfo[i].classIdConv[1],ptr_classInfo[i].classIdConv[2]);}
			heapSortDBClass(ptr_classInfo, subStructureCnt, 3, 1); //3: C, A, T (H is omitted here). Could do with first two
			//for(i=0; i< min(30,subStructureCnt); i++){printf("sorted, converted class %d: %d %d %d\n", i, ptr_classInfo[i].classIdConv[0], ptr_classInfo[i].classIdConv[1],ptr_classInfo[i].classIdConv[2]);}
		}


		/*now we count the number of elts in each category [class letter.fold nr]
		We do this by starting from the first value in ptr_classInfo and extending 
		from the current value as far as possible to the right; this gives the number
		of elt's having that value; do this over by stepping one to the right; keep
		doing this until end of ptr_classInfo*/
		ptr_range[0] = 0;
		ptr_range[1] = subStructureCnt -1;
		hitIndex = 0;
		while(hitIndex < subStructureCnt){
			
			//if(hitIndex%1000 == 0){ getchar();}

			for(pos = 0; pos < 2; pos ++){ 
				//printf("before, range 0: %d 1: %d\n", ptr_range[0], ptr_range[1]);
				bisectSingleExtendClass(hitIndex, ptr_classInfo[hitIndex], pos, ptr_range, ptr_classInfo, 1);	
				//printf("range 0: %d 1: %d\n", ptr_range[0], ptr_range[1]);
			}
			/*record the result*/
			classSize = ptr_range[1] - ptr_range[0] + 1;
			//printf("classSize:%d\n",classSize);
			//printf("range 0: %d 1: %d\n", ptr_range[0], ptr_range[1]);
			for(i = ptr_range[0]; i <= ptr_range[1] ; i++){ ptr_classInfo[i].classSize = classSize;}

			hitIndex = ptr_range[1] + 1;
			ptr_range[0] = hitIndex;
			ptr_range[1] = subStructureCnt -1;

		}
	}

	/*cnt = 0;
	for(i=0; i< 100; i++){
		getchar();
		for(j=0; j < 100; j++){
			printf("class of chain %d: %s was conv to: %d %d %d %d and has size %d\n", cnt, ptr_classInfo[cnt].classId, ptr_classInfo[cnt].classIdConv[0], ptr_classInfo[cnt].classIdConv[1],ptr_classInfo[cnt].classIdConv[2],ptr_classInfo[cnt].classIdConv[3], ptr_classInfo[cnt].classSize );
			cnt += 1;
		}
	}
	cnt = 0;*/
	//return 0;

	/****************************************************************************************************************/
	/*MEMORY ALLOCATION*/
	/****************************************************************************************************************/
	
	/*We allocate to longest chain*/
	returnVal = alloc_init_I_measures(&I_measures, order, full_b, chainNr, maxChainLen, closed_loops_b, 0, &ptr_closedLoopInd);

	ptr_chain = (struct cAlpha *) malloc (maxChainLen*sizeof(struct cAlpha));
	/*init*/
	for ( n=0; n <= maxChainLen -1 ; n++ ){
		
		ptr_chain[n].coords.x = 0.0;
		ptr_chain[n].coords.y = 0.0;
		ptr_chain[n].coords.z = 0.0;

		ptr_chain[n].residueNr = -123;

	}
	ptr_segment = (struct  segment *) malloc (maxChainLen*sizeof(struct segment)); /* for the case that re-allocating of mem to ptr_segment is used instead*/
	/*init*/
	segCoords.s1.x = 0.0;
	segCoords.s1.y = 0.0;
	segCoords.s1.z = 0.0;
	segCoords.s2.x = 0.0;
	segCoords.s2.y = 0.0;
	segCoords.s2.z = 0.0;
	for ( n=0; n <= maxChainLen -1 ; n++ ){

		segCoords.s1 = ptr_chain[n].coords;
		segCoords.s2 = ptr_chain[n+1].coords;

		ptr_segment[n] = segCoords;

	}


	if (write_final_b == 1){

		/*allocate memory to pointer ptr_I_values and initialize (obs: numberOfPerturbations is preset to 1):*/
		returnVal = alloc_init_I_values(&ptr_I_values, subStructureCnt, numberOfPerturbations);

	}


	if (get_windows_b != 0){
		printf("Window lgth:%d\n", windowLength);

		maxNrOfWindows = ceil((double) 1.1*(maxChainLen -windowLength)/stepSize); //times 1.1 to be on the safe side
		/*allocate memory to struct I_windows and initialize*/
		returnVal = alloc_init_I_windows_ptr(&I_windows, order, full_b, maxChainLen, windowLength, stepSize, maxNrOfWindows, 0);

	}

	if (get_windowPairs_b != 0){

		maxNrOfWindows = ceil((double) 1.1*(maxChainLen -windowLength)/stepSize); //times 1.1 to be on the safe side
		/*allocate memory to struct I_windows and initialize*/
		returnVal = alloc_init_I_windowPairs_ptr(&I_windowPairs, order, full_b, maxChainLen, windowLength, stepSize, maxNrOfWindows, 0);

	}

	subStructureCnt =0; /*reset*/
	

	//startComplete = clock();

	/****************************************************************************************************************/
	/* MAIN LOOP */
	/****************************************************************************************************************/
	
	/*loop through the list of filenames in the directory and compute the desired measures all along:*/
	for (fileNr = 0; fileNr <numberOfFiles; fileNr ++){

		if (fileNr%100 == 0){printf("Now at file nr:%d\n", fileNr);}

		if(print_basic_b < 10){
			printf("Now handling file:%s\nfileNr:%d\n",ptr_dirList[fileNr], fileNr);
		}

		if(use_cath_b == 1){
			hitIndexDomain = bisectionalSearchCATHDomain(ptr_fileNameList[fileNr], ptr_cathDomain, nrOfCATHDomainsInList, structureNameLength);
			/*if((hitIndexDomain - lastHitIndexDomain) >1){printf("Str:%s missing\n", ptr_cathDomain[hitIndexDomain -1]);}
			lastHitIndexDomain = hitIndexDomain;
			getchar();*/
			if(hitIndexDomain == -1){
				if(print_basic_b > 0){
					printf("File %s not in domain list\n",ptr_fileNameList[fileNr]);
					//getchar();
				}
				continue;
			} //file name was not in the list
		}


		ptr_fileIn = fopen(ptr_dirList[fileNr], "r");

		I_measures.fileName = ptr_dirList[fileNr];

		if (!ptr_fileIn){
			printf("Sorry, the file: %s, fileNr: %d could not be read\n", ptr_dirList[fileNr], fileNr);
			continue;
		}

		//if (!ptr_fileIn){
		//	printf("Sorry, the file could not be found\n");
		//	getchar();
		//	continue;
		//	//return 0;
		//}

		/*reinit chain info*/
		returnVal = reinit_chainInStr(chainInStr, maxNrOfChains);

		/*get number of chains and their lengths for this file:*/
		//chainInStr = readPDBChainStructure(ptr_fileIn);
		if(use_scop_b ==1){
			readPDBChainStructureSCOP(ptr_fileIn, &chainInStr );
		}
		else if(use_cath_b ==1){
			readPDBDomainStructureCATH(ptr_fileIn, &chainInStr );//will read chain length and chain id into chainInStr
		}
		else{
			readPDBChainStructure2(ptr_fileIn, &chainInStr );
			/*assign a structure name (id to fileName) in case data are not SCOP nor CATH (in these cases the structure name is read in by read */
			for (i=0;i <= chainInStr.numberOfChains -1;i++){
				/*reset*/
				strcpy(chainInStr.ptr_chainInStr[i].structureName, structureNameInit);
				//printf("str name, 1: %s\n", chainInStr.ptr_chainInStr[i].structureName);
				strcpy(chainInStr.ptr_chainInStr[i].structureName, ptr_fileNameList[fileNr]);
				//printf("str name, 2: %s\n", chainInStr.ptr_chainInStr[i].structureName);
				strcat(chainInStr.ptr_chainInStr[i].structureName, chainInStr.ptr_chainInStr[i].chainId);
				//printf("str name, 3: %s\n", chainInStr.ptr_chainInStr[i].structureName);
				//getchar();
			}
		}

		/****************************************************************************************************************/
		/* INNER LOOP */
		/****************************************************************************************************************/

		/*loop through the chains in the current structure/file ... and compute ... :*/
		chainNr = 0; /*reset*/
		for (chainNr = 0; chainNr <= chainInStr.numberOfChains -1 ; chainNr ++){

			printf("Now handling: Chain: %s of %s of length %d\n", chainInStr.ptr_chainInStr[chainNr].chainId, ptr_dirList[fileNr], chainInStr.ptr_chainInStr[chainNr].chainLength);
			//getchar();

			/*We only handle structures of length at most maxChainLength*/
			if(chainInStr.ptr_chainInStr[chainNr].chainLength > maxChainLength){
				printf("Chain: %s of %s has length %d and was skipped\n", chainInStr.ptr_chainInStr[chainNr].chainId, ptr_dirList[fileNr], chainInStr.ptr_chainInStr[chainNr].chainLength);
				chainsSkippedTooLong += 1;
				continue;
			}


			startChain = clock();

			if(use_cath_b == 1){
				/*In case the pdb-file contains several chains: the domain is attached to a single chain (is a subset of 
				a single chain) and we want only to consider that chain. we therefore skip "wrong chains":*/
				if(chainInStr.numberOfChains > 1 && ptr_cathDomain[hitIndexDomain].structureName[4] != chainInStr.ptr_chainInStr[chainNr].chainId[0]){
					if(print_basic_b > 0){
						printf("CATH domain: %s does not sit in this chain: %s\n", ptr_cathDomain[hitIndexDomain].structureName, chainInStr.ptr_chainInStr[chainNr].chainId);
					}
					continue;
				}

				//read into chainInStr some available values from ptr_cathDomain[hitIndexDomain]; obs: hitIndexDomain
				//was found above at beginning of loop over fileNr
				//printf("file is in domain list\n");
				strcpy(chainInStr.ptr_chainInStr[chainNr].structureName, ptr_cathDomain[hitIndexDomain].structureName);
				////adjoin the chainId to he strucure name (in CATH the chain id sits in the structure name, but we do not bother
				////here to single out which chain (if there are several) and so we consider more chains should the pdb-file not
				////be confined to the "right" chain (in some case there are several id chains, and the domain then sits in one of
				////them)
				//strcat(chainInStr.ptr_chainInStr[chainNr].structureName, "_");
				//strcat(chainInStr.ptr_chainInStr[chainNr].structureName, chainInStr.ptr_chainInStr[chainNr].chainId);
				strcpy(chainInStr.ptr_chainInStr[chainNr].classId, ptr_cathDomain[hitIndexDomain].classId);
				for(i=0; i< classifierSize; i++){chainInStr.ptr_chainInStr[chainNr].classIdConv[i] =  ptr_cathDomain[hitIndexDomain].classIdConv[i];}

			}


			/*for some structures (nmr ensembles) a single chain may be represented several times in the 
			same PDB-file; we want to skip multiple copies:*/
			if(strcmp(chainInStr.ptr_chainInStr[chainNr].structureName, lastStructureName) == 0){
				continue;
			}

			I_measures.structureName = chainInStr.ptr_chainInStr[chainNr].structureName;
			I_measures.classId = chainInStr.ptr_chainInStr[chainNr].classId;
			I_measures.chainId = chainInStr.ptr_chainInStr[chainNr].chainId;
			I_measures.chainNr = &chainNr;

			chainLen = chainInStr.ptr_chainInStr[chainNr].chainLength;

			if (print_basic_b ==1){
				printf("sub str cnt:%d\n",subStructureCnt);
				printf("Chain id:%s\n", chainInStr.ptr_chainInStr[chainNr].chainId);
				printf("Chain nr:%d\n", chainNr);
				printf("Chain Length: %d\n", chainLen);
			}

			/*if the chain length is less than strLengthCutOff we skip the computation*/
			if(chainLen < strLengthCutOff){
				//if (print_basic_b ==1){
				printf("Chain %s of %s skipped since length is < %d\n", I_measures.chainId, I_measures.fileName,  strLengthCutOff);
				//}
				//getchar();
				chainsSkippedTooShort += 1;
				continue;
			}

			/*if desired skip the current chain/structure if its class only contains one element ("singleton").
			For this look up the current chain's class in the ptr_classInfo list and fetch the size of the class*/
			if(omitSingletons_b == 1 && ptr_classInfo[subStructureCnt].classSize <2){
				printf("Class of chain: %s of %s was skipped since it is contains only %d elts (singleton)\n", chainInStr.ptr_chainInStr[chainNr].chainId, ptr_dirList[fileNr], ptr_classInfo[subStructureCnt].classSize);
				nrOfSingletons += 1;
				subStructureCnt += 1; //!
				continue;
			}

			if(use_scop_b ==1){
				if(strstr(ptr_classInfo[subStructureCnt].classId, "i") == ptr_classInfo[subStructureCnt].classId || 
					strstr(ptr_classInfo[subStructureCnt].classId, "j") == ptr_classInfo[subStructureCnt].classId ||
					strstr(ptr_classInfo[subStructureCnt].classId, "k") == ptr_classInfo[subStructureCnt].classId ||
					strstr(ptr_classInfo[subStructureCnt].classId, "h") == ptr_classInfo[subStructureCnt].classId){
					printf("Class %s of chain of %s was skipped since it is h, i, j or k\n", chainInStr.ptr_chainInStr[chainNr].chainId, ptr_dirList[fileNr], ptr_classInfo[subStructureCnt].classId);
					subStructureCnt += 1; //!
					continue;
				}
			}

			I_measures.chainLen = &chainLen;

			rewind(ptr_fileIn);
			//ptr_chain = main_readPDB(ptr_fileIn, chainNr, chainLen);
			main_readPDB2(ptr_fileIn, ptr_chain, chainNr, chainLen);

			/*length of segments chain:*/
			L = chainLen - 1;
			/*allocate memory for array of segments*/
			ptr_segment = (struct  segment *) calloc (L, sizeof(struct segment));
			//ptr_segment = realloc(ptr_segment,L*sizeof(struct segment));


			/*populate ptr_segment */
			for ( n=0; n <= L -1 ; n++ ){
				segCoords.s1 = ptr_chain[n].coords;
				segCoords.s2 = ptr_chain[n+1].coords;

				ptr_segment[n] = segCoords;
			}

			/*We skip the chain if it contains a segment which is longer than the globally set
			 stdRealSegLength:*/
			segTooLong_b = 0;
			for ( n=0; n <= L -1 ; n++ ){
				if (distCalphas(ptr_segment[n].s1, ptr_segment[n].s2) > stdRealSegLength){
					segTooLong_b = n+1;
					break;
				}
			} 
			if(segTooLong_b >= 1 ){
				chainSkippedSegmTooLong +=1;
				printf("Chain: %s in structure: %s removed due to too long %d 'th segment\n", chainInStr.ptr_chainInStr[chainNr].chainId, chainInStr.ptr_chainInStr[chainNr].structureName,  segTooLong_b-1);
				continue;
			}

			if(print_b == 1){
				printf ("Segments: ");
				for ( n=0; n<5; n++ ){
					printf ("s1, x: %lf y:%lf z:%lf \n",ptr_segment[n].s1.x, ptr_segment[n].s1.y, ptr_segment[n].s1.z);
					printf ("s2, x: %lf y:%lf z:%lf \n",ptr_segment[n].s2.x, ptr_segment[n].s2.y, ptr_segment[n].s2.z);
				}
			}

			/*(RE)INITIALIZATION*/
			returnVal = init_I_measures(I_measures, order, full_b, chainLen);

			I_measures.order = &order;


			startComp = clock();


			/*call computation of the invariants*/
			if (incl_abs_b == 1 && split_b == 0 && closed_loops_b == 0){
				returnVal = aggrAndW(ptr_segment, chainLen, order, full_b, I_measures);
				}
			if (incl_abs_b == 0 && split_b == 0 && closed_loops_b == 0){
				returnVal = aggrAndW_ExAbs(ptr_segment, chainLen, order, full_b, I_measures);
				}
			if (incl_abs_b == 1 && split_b == 1 && closed_loops_b == 0){
				returnVal = wAll(ptr_segment, chainLen, I_measures.wVal);
				returnVal = aggr(chainLen, order, full_b, I_measures.wVal, I_measures);
			}


			/*if desired compute and record the mutual writhe of pairs of closed-loops and 
			pairs (closed loops, segment) at set lengths. Aim is search for particular geometries 
			(links and "pokes")*/
			if (closed_loops_b == 1){
				returnVal = aggrAndW_wClosedLoops(ptr_segment, chainLen, order, full_b, I_measures, ptr_closedLoopInd, closedLoopLength, closedLoopDist, &nrOfClosedLoops);

				/*examine the closed loops for links and "pokes" and store/write out the results:*/
				returnVal = examineClosedLoops(ptr_fileNameOutClosedLoops, ptr_closedLoopInd, nrOfClosedLoops, I_measures, chainLen, ptr_segment, ptr_fileNameOutChains, ptr_chain, pokeLength, -1,-1);

			}

			//printf("Nr of closed loops found in this structure: %d\n", nrOfClosedLoops);

			/*if desired compute and record the mutual writhe/invariant value of pairs of sub-chains of 
			a set length. Aim is a more general search for particular geometries (e.g. links)*/
			if (invValSubChainPairs_b != 0){

				if(chainLen > 2*subChainLength + 1){

					if (invValSubChainPairs_b == 1){ /*compute only mutual writhe*/ 
						returnVal = writheSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, subChainLength, ptr_fileNameOutChains, ptr_chain, writeSubChainPairsAll_b);
					} 

					if (invValSubChainPairs_b == 2){ /*compute only mutual invariant values for order 1 and order 2; this takes setting the order parameter to 3 or the full_b parameter to 1*/ 

						if (order > 2 || full_b ==1){

							/*genInvariantSubChainPairs(char *ptr_fileNameOut, struct segment *ptr_segment, int chainLength, struct I_ptr I_measures, double **I_measure, char **I_measure_name, int subChainLength, char *ptr_fileNameOutChains, struct cAlpha *ptr_chain, int writeAll_b)*/

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I12, "I12", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I1234_full, "I1234", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I1324_full, "I1324", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I1423 , "I1423", subChainLength, writeSubChainPairsAll_b);

							/*write out the chain:*/
							returnVal = writeChain(ptr_fileNameOutChains, I_measures.fileName,I_measures.structureName, I_measures.chainId, ptr_chain, *I_measures.chainLen);

						}
						else{
							printf("Set full_b to 1 or order = 3 to do this\n");
						}

					}

					if (invValSubChainPairs_b == 3){ /*compute mutual invariant values for all orders; takes setting order = 3*/ 

						if (order == 3){
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I12, "I12", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.Ia12, "Ia12", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I1234, "I1234_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.Ia1234, "Ia1234_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I12a34, "I12a34_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.Ia12a34, "Ia12a34_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I1234_full, "I1234", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I1324, "I1324_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.Ia1324, "Ia1324_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I13a24, "I13a24_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.Ia13a24, "Ia13a24_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I1324_full, "I1324", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I1423, "I1423_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.Ia1423, "Ia1423_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I14a23, "I14a23_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.Ia14a23, "Ia14a23_rel", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I1423, "I1423", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I123456, "I123456", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I123546, "I123546", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I123645, "I123645", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I132456, "I132456", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I132546, "I132546", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I132645, "I132645", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I142356, "I142356", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I142536, "I142536", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I142635, "I142635", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I152346, "I152346", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I152436, "I152436", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I152634, "I152634", subChainLength, writeSubChainPairsAll_b);

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I162345, "I162345", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I162435, "I162435", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I162534, "I162534", subChainLength, writeSubChainPairsAll_b);

							/*write out the chain:*/
							returnVal = writeChain(ptr_fileNameOutChains, I_measures.fileName,I_measures.structureName, I_measures.chainId, ptr_chain, *I_measures.chainLen);

						}
						else{
							printf("Set order = 3 to do this\n");
						}
					}


				}

				else{
					printf("Twice the SubChainLength in search longer than chain length; file: %s, chain id: %s\n", I_measures.fileName, I_measures.chainId);
				}

			}



			/*for structural alignment on windows*/
			if (get_windows_b != 0){


				returnVal = getInvariantsOnWindows(&I_windows, order, windowCoveringType, windowLength, stepSize, L, I_measures);

				if (write_chains_b != 0){
					returnVal = writeChain(fileNameChains, I_measures.fileName, I_measures.structureName, I_measures.chainId, ptr_chain, *I_measures.chainLen);
				}

				
				if (write_windows_b != 0){
					returnVal = writeInvariantsOnWindows(ptr_fileNameOutWindows, I_windows, L, order);
				}
			}


			/*for structural alignment on window pairs*/
			if (get_windowPairs_b != 0){

				returnVal = getInvariantsOnWindowPairs(&I_windowPairs, order, windowCoveringType, windowLength, stepSize, L, I_measures);

				//printf("I ...:%lf\n",I_windowPairs.I12[0][0]);

				if (write_chains_b != 0){
					returnVal = writeChain(fileNameChains, I_measures.fileName, I_measures.structureName, I_measures.chainId, ptr_chain, *I_measures.chainLen);
				}

				
				if (write_windowPairs_b != 0){

					returnVal = writeInvariantsOnwindowPairs(ptr_fileNameOutwindowPairs, I_windowPairs, L, order, writeOnlyDisjointPairs_b);
				}
			}


			endComp = clock();

			//compTime = ((long double) (endComp - startComp))/ CLOCKS_PER_SEC;
			compTime = ((double) (endComp - startComp));
			if(print_b == 1){
				printf("Comp time: %f\n",compTime);
				//printf("CPU time spend for measures on this file (only computation/aggregation):%lf \n", compTime);
			}
			if (compTime > max_compTime){
				max_compTime = compTime;
			}
			compTimeComplete += compTime;

			/*write all I-values to file if wanted:*/
			if (write_all_b == 1){

				returnVal = writeAllIvaluesToFile(ptr_fileNameOutAll, full_b, I_measures, pertNo);

			}


			/*look at the w-values etc:*/
			if (print_b == 1){

				/*for ( i=0; i<=L-1; i++ ){
					for ( j=0; j<=L-1; j++ ){
						printf("w[%d][%d]:%lf ", i,j, I_measures.wVal[i][j]);		
					}
				}*/

				if (order >= 1){
					printf("I12[%d][%d]:%lf ", 0,L-1, I_measures.I12[0][L-1]);
					printf("Ia12[%d][%d]:%lf ", 0,L-1, I_measures.Ia12[0][L-1]);
					/*for ( i=0; i<=L-1; i++ ){
						for ( j=0; j<=L-1; j++ ){
							printf("I12[%d][%d]:%lf ", i,j, I_measures.I12[i][j]);
						}
					}*/
				}

				if (order >= 2){
					printf("I1234[%d][%d]:%lf ", 0,L-1, I_measures.I1234[0][L-1]);
					/*for ( i=0; i<=L-1; i++ ){
						for ( j=0; j<=L-1; j++ ){
							printf("I1234[%d][%d]:%lf ", i,j, I_measures.I1234[i][j]);
						}
					}*/
					if (full_b == 1 && order ==2){
						printf("I1234_full[%d][%d]:%lf ", 0,L-1, I_measures.I1234_full[0][L-1]);
						printf("I1324_full[%d][%d]:%lf ", 0,L-1, I_measures.I1324_full[0][L-1]);
					}
				}

				if (order >= 3){
					printf("I123456[%d][%d]:%lf ", 0,L-1, I_measures.I123456[0][L-1]);
					printf("I142356[%d][%d]:%lf ", 0,L-1, I_measures.I142356[0][L-1]);
					printf("I132645[%d][%d]:%lf ", 0,L-1, I_measures.I132645[0][L-1]);
					printf("I132546[%d][%d]:%lf ", 0,L-1, I_measures.I132546[0][L-1]);
					printf("I142536[%d][%d]:%lf ", 0,L-1, I_measures.I142536[0][L-1]);
					printf("I142635[%d][%d]:%lf ", 0,L-1, I_measures.I142635[0][L-1]);
					printf("I1423_full0[%d][%d]:%lf ", 0,L-1, I_measures.I1423_full0[0][L-1]);
					printf("I1423_full2[%d][%d]:%lf ", 0,L-1, I_measures.I1423_full2[0][L-1]);
					printf("I1324_full2[%d][%d]:%lf ", 0,L-1, I_measures.I1324_full2[0][L-1]);
					printf("I1234_full[%d][%d]:%lf ", 0,L-1, I_measures.I1234_full[0][L-1]);
					printf("I1324_full[%d][%d]:%lf ", 0,L-1, I_measures.I1324_full[0][L-1]);

				}
			}


			/*free(ptr_chain);
			free(ptr_segment);	*/

			endChain = clock();
			//timeChain = ((double) (endChain - startChain)) / CLOCKS_PER_SEC;
			timeChain = ((double) (endChain - startChain));
			//printf("CPU time spend for measures on this file incl load of chain and more:%lf \n", timeChain);
			if (timeChain > max_timeChain){
				max_timeChain =timeChain;
			}

			/*collect final results for write-out if that's desired*/
			if (write_final_b == 1){
				/*reads I-values at corner 0,L-1 into the ptr_I_values ptr at index = subStructureCnt; pertNo is preset to 0:*/
				returnVal = collectIvalues(ptr_I_values, I_measures, subStructureCnt, pertNo, timeChain, compTime);
		
			}

			subStructureCnt += 1; /*will count and index the sub-structures for which we compute the invariants*/

			strcpy(lastStructureName, chainInStr.ptr_chainInStr[chainNr].structureName); 
		
		} /*end chainNr loop*/

		
		fclose(ptr_fileIn);

	} /*end fileNr loop*/


	/*free the globally allocated memory*/
	free_I_measures(I_measures, order, full_b, maxChainLen);
	free(ptr_segment);	
	free(ptr_chain);
	free(dirContent.ptr_dirList); //this is allocated in the ListDirectoryContents2 fct 

	if (write_final_b == 1){

		/*generate file name for the output*/ 
		sprintf(orderStr, "%d", order);
		strcat(fileNameOut, outputPath);
		strcat(fileNameOut, fileNameOut1);
		strcat(fileNameOut, orderStr);
		strcat(fileNameOut, fileNameOut2);
		strcat(fileNameOut, fileNameOut3);
		strcat(fileNameOut, DBName);
		strcat(fileNameOut, extensionName);
		printf("Results for chain (only) will be written to: %s\n", fileNameOut);

		ptr_fileNameOut = fileNameOut;

		/*clear the file by opening it in w-mode and then closing it again*/
	    fclose(fopen(ptr_fileNameOut, "w"));

		/*call writing the output (obs: numberOfPerturbations is preset to 1)*/
		returnVal = writeIvaluesToFile(ptr_fileNameOut, ptr_I_values, subStructureCnt, numberOfPerturbations);

		/*free the memory used for the purpose:*/
		free(ptr_I_values);
		 
	}

	endComplete = clock();
	//printf("start clocks:%d, end clocks: %d, clocks_per sec:%d", endComplete, startComplete, CLOCKS_PER_SEC);
	//timeComplete = ((long double) (endComplete - startComplete)) / CLOCKS_PER_SEC;
	timeComplete = ((double) (endComplete - startComplete))/ CLOCKS_PER_SEC;
	compTimeComplete = ((double) (compTimeComplete)) / CLOCKS_PER_SEC;
	max_timeChain = max_timeChain/ CLOCKS_PER_SEC;
	printf("CPU time spend for all files incl loading data:%lf \n", timeComplete);
	printf("CPU time spend for computations on all files:%lf \n", compTimeComplete);
	printf("Max CPU time spend for measures across all files:%lf \n", max_timeChain);

	printf("Done number of files:%d \n", fileNr);
	printf("Done number of chains/sub-structures:%d \n", subStructureCnt);
	printf("%d too short and %d too long chains were skipped\n", chainsSkippedTooShort, chainsSkippedTooLong);
	printf("%d chains were skipped due to a too long segment (longer than sqrt of %d)\n",chainSkippedSegmTooLong, stdRealSegLength);
	if(omitSingletons_b == 1){
		printf("Nr of singletons (skipped): %d\n", nrOfSingletons);
	}

	//getchar(); 
	
	returnVal = 1;

	return returnVal;
}
 

/***************************************************************************
* License text
*****************************************************************************/
/*
                    GNU GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.

                            Preamble

  The GNU General Public License is a free, copyleft license for
software and other kinds of works.

  The licenses for most software and other practical works are designed
to take away your freedom to share and change the works.  By contrast,
the GNU General Public License is intended to guarantee your freedom to
share and change all versions of a program--to make sure it remains free
software for all its users.  We, the Free Software Foundation, use the
GNU General Public License for most of our software; it applies also to
any other work released this way by its authors.  You can apply it to
your programs, too.

  When we speak of free software, we are referring to freedom, not
price.  Our General Public Licenses are designed to make sure that you
have the freedom to distribute copies of free software (and charge for
them if you wish), that you receive source code or can get it if you
want it, that you can change the software or use pieces of it in new
free programs, and that you know you can do these things.

  To protect your rights, we need to prevent others from denying you
these rights or asking you to surrender the rights.  Therefore, you have
certain responsibilities if you distribute copies of the software, or if
you modify it: responsibilities to respect the freedom of others.

  For example, if you distribute copies of such a program, whether
gratis or for a fee, you must pass on to the recipients the same
freedoms that you received.  You must make sure that they, too, receive
or can get the source code.  And you must show them these terms so they
know their rights.

  Developers that use the GNU GPL protect your rights with two steps:
(1) assert copyright on the software, and (2) offer you this License
giving you legal permission to copy, distribute and/or modify it.

  For the developers' and authors' protection, the GPL clearly explains
that there is no warranty for this free software.  For both users' and
authors' sake, the GPL requires that modified versions be marked as
changed, so that their problems will not be attributed erroneously to
authors of previous versions.

  Some devices are designed to deny users access to install or run
modified versions of the software inside them, although the manufacturer
can do so.  This is fundamentally incompatible with the aim of
protecting users' freedom to change the software.  The systematic
pattern of such abuse occurs in the area of products for individuals to
use, which is precisely where it is most unacceptable.  Therefore, we
have designed this version of the GPL to prohibit the practice for those
products.  If such problems arise substantially in other domains, we
stand ready to extend this provision to those domains in future versions
of the GPL, as needed to protect the freedom of users.

  Finally, every program is threatened constantly by software patents.
States should not allow patents to restrict development and use of
software on general-purpose computers, but in those that do, we wish to
avoid the special danger that patents applied to a free program could
make it effectively proprietary.  To prevent this, the GPL assures that
patents cannot be used to render the program non-free.

  The precise terms and conditions for copying, distribution and
modification follow.

                       TERMS AND CONDITIONS

  0. Definitions.

  "This License" refers to version 3 of the GNU General Public License.

  "Copyright" also means copyright-like laws that apply to other kinds of
works, such as semiconductor masks.

  "The Program" refers to any copyrightable work licensed under this
License.  Each licensee is addressed as "you".  "Licensees" and
"recipients" may be individuals or organizations.

  To "modify" a work means to copy from or adapt all or part of the work
in a fashion requiring copyright permission, other than the making of an
exact copy.  The resulting work is called a "modified version" of the
earlier work or a work "based on" the earlier work.

  A "covered work" means either the unmodified Program or a work based
on the Program.

  To "propagate" a work means to do anything with it that, without
permission, would make you directly or secondarily liable for
infringement under applicable copyright law, except executing it on a
computer or modifying a private copy.  Propagation includes copying,
distribution (with or without modification), making available to the
public, and in some countries other activities as well.

  To "convey" a work means any kind of propagation that enables other
parties to make or receive copies.  Mere interaction with a user through
a computer network, with no transfer of a copy, is not conveying.

  An interactive user interface displays "Appropriate Legal Notices"
to the extent that it includes a convenient and prominently visible
feature that (1) displays an appropriate copyright notice, and (2)
tells the user that there is no warranty for the work (except to the
extent that warranties are provided), that licensees may convey the
work under this License, and how to view a copy of this License.  If
the interface presents a list of user commands or options, such as a
menu, a prominent item in the list meets this criterion.

  1. Source Code.

  The "source code" for a work means the preferred form of the work
for making modifications to it.  "Object code" means any non-source
form of a work.

  A "Standard Interface" means an interface that either is an official
standard defined by a recognized standards body, or, in the case of
interfaces specified for a particular programming language, one that
is widely used among developers working in that language.

  The "System Libraries" of an executable work include anything, other
than the work as a whole, that (a) is included in the normal form of
packaging a Major Component, but which is not part of that Major
Component, and (b) serves only to enable use of the work with that
Major Component, or to implement a Standard Interface for which an
implementation is available to the public in source code form.  A
"Major Component", in this context, means a major essential component
(kernel, window system, and so on) of the specific operating system
(if any) on which the executable work runs, or a compiler used to
produce the work, or an object code interpreter used to run it.

  The "Corresponding Source" for a work in object code form means all
the source code needed to generate, install, and (for an executable
work) run the object code and to modify the work, including scripts to
control those activities.  However, it does not include the work's
System Libraries, or general-purpose tools or generally available free
programs which are used unmodified in performing those activities but
which are not part of the work.  For example, Corresponding Source
includes interface definition files associated with source files for
the work, and the source code for shared libraries and dynamically
linked subprograms that the work is specifically designed to require,
such as by intimate data communication or control flow between those
subprograms and other parts of the work.

  The Corresponding Source need not include anything that users
can regenerate automatically from other parts of the Corresponding
Source.

  The Corresponding Source for a work in source code form is that
same work.

  2. Basic Permissions.

  All rights granted under this License are granted for the term of
copyright on the Program, and are irrevocable provided the stated
conditions are met.  This License explicitly affirms your unlimited
permission to run the unmodified Program.  The output from running a
covered work is covered by this License only if the output, given its
content, constitutes a covered work.  This License acknowledges your
rights of fair use or other equivalent, as provided by copyright law.

  You may make, run and propagate covered works that you do not
convey, without conditions so long as your license otherwise remains
in force.  You may convey covered works to others for the sole purpose
of having them make modifications exclusively for you, or provide you
with facilities for running those works, provided that you comply with
the terms of this License in conveying all material for which you do
not control copyright.  Those thus making or running the covered works
for you must do so exclusively on your behalf, under your direction
and control, on terms that prohibit them from making any copies of
your copyrighted material outside their relationship with you.

  Conveying under any other circumstances is permitted solely under
the conditions stated below.  Sublicensing is not allowed; section 10
makes it unnecessary.

  3. Protecting Users' Legal Rights From Anti-Circumvention Law.

  No covered work shall be deemed part of an effective technological
measure under any applicable law fulfilling obligations under article
11 of the WIPO copyright treaty adopted on 20 December 1996, or
similar laws prohibiting or restricting circumvention of such
measures.

  When you convey a covered work, you waive any legal power to forbid
circumvention of technological measures to the extent such circumvention
is effected by exercising rights under this License with respect to
the covered work, and you disclaim any intention to limit operation or
modification of the work as a means of enforcing, against the work's
users, your or third parties' legal rights to forbid circumvention of
technological measures.

  4. Conveying Verbatim Copies.

  You may convey verbatim copies of the Program's source code as you
receive it, in any medium, provided that you conspicuously and
appropriately publish on each copy an appropriate copyright notice;
keep intact all notices stating that this License and any
non-permissive terms added in accord with section 7 apply to the code;
keep intact all notices of the absence of any warranty; and give all
recipients a copy of this License along with the Program.

  You may charge any price or no price for each copy that you convey,
and you may offer support or warranty protection for a fee.

  5. Conveying Modified Source Versions.

  You may convey a work based on the Program, or the modifications to
produce it from the Program, in the form of source code under the
terms of section 4, provided that you also meet all of these conditions:

    a) The work must carry prominent notices stating that you modified
    it, and giving a relevant date.

    b) The work must carry prominent notices stating that it is
    released under this License and any conditions added under section
    7.  This requirement modifies the requirement in section 4 to
    "keep intact all notices".

    c) You must license the entire work, as a whole, under this
    License to anyone who comes into possession of a copy.  This
    License will therefore apply, along with any applicable section 7
    additional terms, to the whole of the work, and all its parts,
    regardless of how they are packaged.  This License gives no
    permission to license the work in any other way, but it does not
    invalidate such permission if you have separately received it.

    d) If the work has interactive user interfaces, each must display
    Appropriate Legal Notices; however, if the Program has interactive
    interfaces that do not display Appropriate Legal Notices, your
    work need not make them do so.

  A compilation of a covered work with other separate and independent
works, which are not by their nature extensions of the covered work,
and which are not combined with it such as to form a larger program,
in or on a volume of a storage or distribution medium, is called an
"aggregate" if the compilation and its resulting copyright are not
used to limit the access or legal rights of the compilation's users
beyond what the individual works permit.  Inclusion of a covered work
in an aggregate does not cause this License to apply to the other
parts of the aggregate.

  6. Conveying Non-Source Forms.

  You may convey a covered work in object code form under the terms
of sections 4 and 5, provided that you also convey the
machine-readable Corresponding Source under the terms of this License,
in one of these ways:

    a) Convey the object code in, or embodied in, a physical product
    (including a physical distribution medium), accompanied by the
    Corresponding Source fixed on a durable physical medium
    customarily used for software interchange.

    b) Convey the object code in, or embodied in, a physical product
    (including a physical distribution medium), accompanied by a
    written offer, valid for at least three years and valid for as
    long as you offer spare parts or customer support for that product
    model, to give anyone who possesses the object code either (1) a
    copy of the Corresponding Source for all the software in the
    product that is covered by this License, on a durable physical
    medium customarily used for software interchange, for a price no
    more than your reasonable cost of physically performing this
    conveying of source, or (2) access to copy the
    Corresponding Source from a network server at no charge.

    c) Convey individual copies of the object code with a copy of the
    written offer to provide the Corresponding Source.  This
    alternative is allowed only occasionally and noncommercially, and
    only if you received the object code with such an offer, in accord
    with subsection 6b.

    d) Convey the object code by offering access from a designated
    place (gratis or for a charge), and offer equivalent access to the
    Corresponding Source in the same way through the same place at no
    further charge.  You need not require recipients to copy the
    Corresponding Source along with the object code.  If the place to
    copy the object code is a network server, the Corresponding Source
    may be on a different server (operated by you or a third party)
    that supports equivalent copying facilities, provided you maintain
    clear directions next to the object code saying where to find the
    Corresponding Source.  Regardless of what server hosts the
    Corresponding Source, you remain obligated to ensure that it is
    available for as long as needed to satisfy these requirements.

    e) Convey the object code using peer-to-peer transmission, provided
    you inform other peers where the object code and Corresponding
    Source of the work are being offered to the general public at no
    charge under subsection 6d.

  A separable portion of the object code, whose source code is excluded
from the Corresponding Source as a System Library, need not be
included in conveying the object code work.

  A "User Product" is either (1) a "consumer product", which means any
tangible personal property which is normally used for personal, family,
or household purposes, or (2) anything designed or sold for incorporation
into a dwelling.  In determining whether a product is a consumer product,
doubtful cases shall be resolved in favor of coverage.  For a particular
product received by a particular user, "normally used" refers to a
typical or common use of that class of product, regardless of the status
of the particular user or of the way in which the particular user
actually uses, or expects or is expected to use, the product.  A product
is a consumer product regardless of whether the product has substantial
commercial, industrial or non-consumer uses, unless such uses represent
the only significant mode of use of the product.

  "Installation Information" for a User Product means any methods,
procedures, authorization keys, or other information required to install
and execute modified versions of a covered work in that User Product from
a modified version of its Corresponding Source.  The information must
suffice to ensure that the continued functioning of the modified object
code is in no case prevented or interfered with solely because
modification has been made.

  If you convey an object code work under this section in, or with, or
specifically for use in, a User Product, and the conveying occurs as
part of a transaction in which the right of possession and use of the
User Product is transferred to the recipient in perpetuity or for a
fixed term (regardless of how the transaction is characterized), the
Corresponding Source conveyed under this section must be accompanied
by the Installation Information.  But this requirement does not apply
if neither you nor any third party retains the ability to install
modified object code on the User Product (for example, the work has
been installed in ROM).

  The requirement to provide Installation Information does not include a
requirement to continue to provide support service, warranty, or updates
for a work that has been modified or installed by the recipient, or for
the User Product in which it has been modified or installed.  Access to a
network may be denied when the modification itself materially and
adversely affects the operation of the network or violates the rules and
protocols for communication across the network.

  Corresponding Source conveyed, and Installation Information provided,
in accord with this section must be in a format that is publicly
documented (and with an implementation available to the public in
source code form), and must require no special password or key for
unpacking, reading or copying.

  7. Additional Terms.

  "Additional permissions" are terms that supplement the terms of this
License by making exceptions from one or more of its conditions.
Additional permissions that are applicable to the entire Program shall
be treated as though they were included in this License, to the extent
that they are valid under applicable law.  If additional permissions
apply only to part of the Program, that part may be used separately
under those permissions, but the entire Program remains governed by
this License without regard to the additional permissions.

  When you convey a copy of a covered work, you may at your option
remove any additional permissions from that copy, or from any part of
it.  (Additional permissions may be written to require their own
removal in certain cases when you modify the work.)  You may place
additional permissions on material, added by you to a covered work,
for which you have or can give appropriate copyright permission.

  Notwithstanding any other provision of this License, for material you
add to a covered work, you may (if authorized by the copyright holders of
that material) supplement the terms of this License with terms:

    a) Disclaiming warranty or limiting liability differently from the
    terms of sections 15 and 16 of this License; or

    b) Requiring preservation of specified reasonable legal notices or
    author attributions in that material or in the Appropriate Legal
    Notices displayed by works containing it; or

    c) Prohibiting misrepresentation of the origin of that material, or
    requiring that modified versions of such material be marked in
    reasonable ways as different from the original version; or

    d) Limiting the use for publicity purposes of names of licensors or
    authors of the material; or

    e) Declining to grant rights under trademark law for use of some
    trade names, trademarks, or service marks; or

    f) Requiring indemnification of licensors and authors of that
    material by anyone who conveys the material (or modified versions of
    it) with contractual assumptions of liability to the recipient, for
    any liability that these contractual assumptions directly impose on
    those licensors and authors.

  All other non-permissive additional terms are considered "further
restrictions" within the meaning of section 10.  If the Program as you
received it, or any part of it, contains a notice stating that it is
governed by this License along with a term that is a further
restriction, you may remove that term.  If a license document contains
a further restriction but permits relicensing or conveying under this
License, you may add to a covered work material governed by the terms
of that license document, provided that the further restriction does
not survive such relicensing or conveying.

  If you add terms to a covered work in accord with this section, you
must place, in the relevant source files, a statement of the
additional terms that apply to those files, or a notice indicating
where to find the applicable terms.

  Additional terms, permissive or non-permissive, may be stated in the
form of a separately written license, or stated as exceptions;
the above requirements apply either way.

  8. Termination.

  You may not propagate or modify a covered work except as expressly
provided under this License.  Any attempt otherwise to propagate or
modify it is void, and will automatically terminate your rights under
this License (including any patent licenses granted under the third
paragraph of section 11).

  However, if you cease all violation of this License, then your
license from a particular copyright holder is reinstated (a)
provisionally, unless and until the copyright holder explicitly and
finally terminates your license, and (b) permanently, if the copyright
holder fails to notify you of the violation by some reasonable means
prior to 60 days after the cessation.

  Moreover, your license from a particular copyright holder is
reinstated permanently if the copyright holder notifies you of the
violation by some reasonable means, this is the first time you have
received notice of violation of this License (for any work) from that
copyright holder, and you cure the violation prior to 30 days after
your receipt of the notice.

  Termination of your rights under this section does not terminate the
licenses of parties who have received copies or rights from you under
this License.  If your rights have been terminated and not permanently
reinstated, you do not qualify to receive new licenses for the same
material under section 10.

  9. Acceptance Not Required for Having Copies.

  You are not required to accept this License in order to receive or
run a copy of the Program.  Ancillary propagation of a covered work
occurring solely as a consequence of using peer-to-peer transmission
to receive a copy likewise does not require acceptance.  However,
nothing other than this License grants you permission to propagate or
modify any covered work.  These actions infringe copyright if you do
not accept this License.  Therefore, by modifying or propagating a
covered work, you indicate your acceptance of this License to do so.

  10. Automatic Licensing of Downstream Recipients.

  Each time you convey a covered work, the recipient automatically
receives a license from the original licensors, to run, modify and
propagate that work, subject to this License.  You are not responsible
for enforcing compliance by third parties with this License.

  An "entity transaction" is a transaction transferring control of an
organization, or substantially all assets of one, or subdividing an
organization, or merging organizations.  If propagation of a covered
work results from an entity transaction, each party to that
transaction who receives a copy of the work also receives whatever
licenses to the work the party's predecessor in interest had or could
give under the previous paragraph, plus a right to possession of the
Corresponding Source of the work from the predecessor in interest, if
the predecessor has it or can get it with reasonable efforts.

  You may not impose any further restrictions on the exercise of the
rights granted or affirmed under this License.  For example, you may
not impose a license fee, royalty, or other charge for exercise of
rights granted under this License, and you may not initiate litigation
(including a cross-claim or counterclaim in a lawsuit) alleging that
any patent claim is infringed by making, using, selling, offering for
sale, or importing the Program or any portion of it.

  11. Patents.

  A "contributor" is a copyright holder who authorizes use under this
License of the Program or a work on which the Program is based.  The
work thus licensed is called the contributor's "contributor version".

  A contributor's "essential patent claims" are all patent claims
owned or controlled by the contributor, whether already acquired or
hereafter acquired, that would be infringed by some manner, permitted
by this License, of making, using, or selling its contributor version,
but do not include claims that would be infringed only as a
consequence of further modification of the contributor version.  For
purposes of this definition, "control" includes the right to grant
patent sublicenses in a manner consistent with the requirements of
this License.

  Each contributor grants you a non-exclusive, worldwide, royalty-free
patent license under the contributor's essential patent claims, to
make, use, sell, offer for sale, import and otherwise run, modify and
propagate the contents of its contributor version.

  In the following three paragraphs, a "patent license" is any express
agreement or commitment, however denominated, not to enforce a patent
(such as an express permission to practice a patent or covenant not to
sue for patent infringement).  To "grant" such a patent license to a
party means to make such an agreement or commitment not to enforce a
patent against the party.

  If you convey a covered work, knowingly relying on a patent license,
and the Corresponding Source of the work is not available for anyone
to copy, free of charge and under the terms of this License, through a
publicly available network server or other readily accessible means,
then you must either (1) cause the Corresponding Source to be so
available, or (2) arrange to deprive yourself of the benefit of the
patent license for this particular work, or (3) arrange, in a manner
consistent with the requirements of this License, to extend the patent
license to downstream recipients.  "Knowingly relying" means you have
actual knowledge that, but for the patent license, your conveying the
covered work in a country, or your recipient's use of the covered work
in a country, would infringe one or more identifiable patents in that
country that you have reason to believe are valid.

  If, pursuant to or in connection with a single transaction or
arrangement, you convey, or propagate by procuring conveyance of, a
covered work, and grant a patent license to some of the parties
receiving the covered work authorizing them to use, propagate, modify
or convey a specific copy of the covered work, then the patent license
you grant is automatically extended to all recipients of the covered
work and works based on it.

  A patent license is "discriminatory" if it does not include within
the scope of its coverage, prohibits the exercise of, or is
conditioned on the non-exercise of one or more of the rights that are
specifically granted under this License.  You may not convey a covered
work if you are a party to an arrangement with a third party that is
in the business of distributing software, under which you make payment
to the third party based on the extent of your activity of conveying
the work, and under which the third party grants, to any of the
parties who would receive the covered work from you, a discriminatory
patent license (a) in connection with copies of the covered work
conveyed by you (or copies made from those copies), or (b) primarily
for and in connection with specific products or compilations that
contain the covered work, unless you entered into that arrangement,
or that patent license was granted, prior to 28 March 2007.

  Nothing in this License shall be construed as excluding or limiting
any implied license or other defenses to infringement that may
otherwise be available to you under applicable patent law.

  12. No Surrender of Others' Freedom.

  If conditions are imposed on you (whether by court order, agreement or
otherwise) that contradict the conditions of this License, they do not
excuse you from the conditions of this License.  If you cannot convey a
covered work so as to satisfy simultaneously your obligations under this
License and any other pertinent obligations, then as a consequence you may
not convey it at all.  For example, if you agree to terms that obligate you
to collect a royalty for further conveying from those to whom you convey
the Program, the only way you could satisfy both those terms and this
License would be to refrain entirely from conveying the Program.

  13. Use with the GNU Affero General Public License.

  Notwithstanding any other provision of this License, you have
permission to link or combine any covered work with a work licensed
under version 3 of the GNU Affero General Public License into a single
combined work, and to convey the resulting work.  The terms of this
License will continue to apply to the part which is the covered work,
but the special requirements of the GNU Affero General Public License,
section 13, concerning interaction through a network will apply to the
combination as such.

  14. Revised Versions of this License.

  The Free Software Foundation may publish revised and/or new versions of
the GNU General Public License from time to time.  Such new versions will
be similar in spirit to the present version, but may differ in detail to
address new problems or concerns.

  Each version is given a distinguishing version number.  If the
Program specifies that a certain numbered version of the GNU General
Public License "or any later version" applies to it, you have the
option of following the terms and conditions either of that numbered
version or of any later version published by the Free Software
Foundation.  If the Program does not specify a version number of the
GNU General Public License, you may choose any version ever published
by the Free Software Foundation.

  If the Program specifies that a proxy can decide which future
versions of the GNU General Public License can be used, that proxy's
public statement of acceptance of a version permanently authorizes you
to choose that version for the Program.

  Later license versions may give you additional or different
permissions.  However, no additional obligations are imposed on any
author or copyright holder as a result of your choosing to follow a
later version.

  15. Disclaimer of Warranty.

  THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY
OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM
IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF
ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

  16. Limitation of Liability.

  IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS
THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY
GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE
USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF
DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD
PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

  17. Interpretation of Sections 15 and 16.

  If the disclaimer of warranty and limitation of liability provided
above cannot be given local legal effect according to their terms,
reviewing courts shall apply local law that most closely approximates
an absolute waiver of all civil liability in connection with the
Program, unless a warranty or assumption of liability accompanies a
copy of the Program in return for a fee.

                     END OF TERMS AND CONDITIONS

            How to Apply These Terms to Your New Programs

  If you develop a new program, and you want it to be of the greatest
possible use to the public, the best way to achieve this is to make it
free software which everyone can redistribute and change under these terms.

  To do so, attach the following notices to the program.  It is safest
to attach them to the start of each source file to most effectively
state the exclusion of warranty; and each file should have at least
the "copyright" line and a pointer to where the full notice is found.

    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

Also add information on how to contact you by electronic and paper mail.

  If the program does terminal interaction, make it output a short
notice like this when it starts in an interactive mode:

    <program>  Copyright (C) <year>  <name of author>
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.

The hypothetical commands `show w' and `show c' should show the appropriate
parts of the General Public License.  Of course, your program's commands
might be different; for a GUI interface, you would use an "about box".

  You should also get your employer (if you work as a programmer) or school,
if any, to sign a "copyright disclaimer" for the program, if necessary.
For more information on this, and how to apply and follow the GNU GPL, see
<https://www.gnu.org/licenses/>.

  The GNU General Public License does not permit incorporating your program
into proprietary programs.  If your program is a subroutine library, you
may consider it more useful to permit linking proprietary applications with
the library.  If this is what you want to do, use the GNU Lesser General
Public License instead of this License.  But first, please read
<https://www.gnu.org/licenses/why-not-lgpl.html>.

*/
