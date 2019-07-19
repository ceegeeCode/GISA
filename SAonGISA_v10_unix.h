/* **********************************************************************************************************************************
*
*
* ***************************************  Header for SAonGISA_v10 *******************************************************************************
*
* See main code (SAonGISA_v*) for do and how to cite.
*
* Author: Christian Grønbæk
*
**********************************************************************************************/


#ifndef SAonGISA_V10_unix
#define SAonGISA_V10_unix

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>
#include<unistd.h>
#include<sys/types.h>
#include "GISA_v10_unix.h"
#include "sortingUtils_v4_unix.h"
#include "dataTypesMatching_v3_unix.h"
#include "mathUtils.h"

#define _SECURE_SCL 0
#define maxNrOfInvariants 14
#define maxWindowInfoSize 50


/**************************************************
************ Prototype functions ******************
***************************************************/


int * readDBresultsWindowsToPtr(struct db_I_windows_order2_ptr **, char[1000], int, int *);

int * loadDBWindowsFromDBPairs(struct db_I_windows_order2_ptr **, struct db_I_windowPairs_order2_ptr *, int, int );


struct I_windows_order2_meanStddev normalizeDBresultsWindowsPtr(struct db_I_windows_order2_ptr *, int  );

int binDBresultsWindows(struct binned_I_windows **, struct db_I_windows_order2_ptr *, int , int , int, double **, int );

int binQueryResultsWindows(struct binned_I_windows *, struct I_windows_ptr , double **,  int , double **, int );

void normalizeQueryWindowsPtr(struct I_windows_ptr , struct I_windows_order2_meanStddev, int);

int binQueryResultsWindowPairs(struct binned_I_windowPairs *, struct I_windowPairs_ptr , double ***,  int , double **, int , int );

int *readDBresultsWindowPairsToPtr(struct db_I_windowPairs_order2_ptr **, char[1000], int, int * );

int readSortDistData(FILE *, double **);

int generateBins(int , struct db_I_windows_order2_ptr *, int , int , int , double ***, int, double );

int binArray(int *, double *, int , double **, int );

double distEuklid(double *, double *, int);

double relDistEuklid(double *, double *, int);

void collect_Irarity_windows(struct Irarity_windows *, struct I_windows_ptr *, int , double );

void collect_Irarity_windowPairs(struct Irarity_windowPairs *, struct I_windowPairs_ptr *, int , int, int, double );

void collect_queryRawScore(struct queryRawScore *, int , struct I_windows_ptr , struct I_windowPairs_ptr , char *, double , double , int , double );

void simpleScore(int , int , int, int, int, int, struct matchScore *ptr_Scores, int nrOfDbRec, struct matchCandidate *ptr_matchCandidates, int topN, int lengthStructureName, float rarityFactor);

void simpleScore_1(int, struct matchScore *, int, struct matchCandidate *, int);

void distDistrScore(double *, int , int, int, int, int, int , int , struct matchScore *, int , struct matchCandidate *, int , int , int);

//void scoreFct(double *, int , int , struct matchScore *, int , int , struct matchCandidate *, int , int , int , struct floatArrayWLength , struct floatArrayWLength );
void scoreFct(double *, int , int , int, int , int, int , int, struct matchScore *, int , int , struct matchCandidate *, int , int , int , struct floatArrayWLength , struct floatArrayWLength );

int getHitStatistics(char * , int , int , int );

int getHitStatisticsCAT(char * scoresFileName, int useSCOP_b, int useCATH_b, int topN);


/*for housekeeping in matching */
double collectMatchCandidates2(struct matchCandidate *, int , int , int , double , struct I_windows_ptr );

double collectMatchCandidates(struct matchCandidate *, int , int *, int lMin, double , int , struct db_I_windows_order2_ptr *);

double collectMatchCandidates_windowPairs(struct matchCandidate_windowPair *, int , int *, int , int , double , int, struct db_I_windowPairs_order2_ptr *);

int alloc_init_matchCandidates(struct matchCandidate **, int );

int reinit_matchCandidates(struct matchCandidate *, int );

int alloc_init_matchScores(struct matchScore **, struct db_I_windows_order2_ptr *, int , int , int);

int reinit_matchScores(struct matchScore *, int , int , int );

/*for holding DB results*/
int alloc_init_DBresultsWindows(struct db_I_windows_order2_ptr **, int , int , int );

/*as alloc_init_DBresultsWindows, but for window pairs:*/
int alloc_init_DBresultsWindowPairs(struct db_I_windowPairs_order2_ptr **, int , int , int , int );

/*allocate and init ptr to binned I-windows, ie to "words of GIs"*/
int alloc_init_ptr_binned_I_windows(struct binned_I_windows **, int , int, int );

/*ptr to array over windows of ptr to binned I-windows; to hold the match set for all windows in a query*/
int alloc_init_ptr_ptr_binned_I_windows(struct binned_I_windows ***, int, int , int, int );

/*allocate and init ptr to binned I-windows, ie to "words of GIs"*/
int alloc_init_ptr_binned_I_windowPairs(struct binned_I_windowPairs **, int , int, int );




/*Main blocks*/

int main_makeDB();

int main_match();


int rawRarity0(int, char[100], char[100], 
				char[200], char[1000], char[2000], 
				char[100], char[1000], int , 
				int ,  int , int, 
				double , int , int , 
				char *, int , int, 
				char *, int , int , 
				int , int , int , 
				int ,  int , int, 
				int , int , int, 
				int );



int rawRarity1(int, char[1000], char[100], 
				char[100], char[100], char[200], 
				char[1000], char[2000], int ,
			    double, int , int , 
				int , char *, int, 
				int , char *, int, 
				int , int , int , 
				double , int , int,
				int , int , int , 
				int , int , int , 
				int , int , int , 
				int ,  int , int , 
				int, int);



int rawRarity2(int, char[1000], char[100], 
				char[1000], char[100], char[100], 
				char[100], char[200], char[1000],  
				char[1000], char[2000], 
				double, int, 
			    int, int , int , 
			   	char *, int , int, 
			   	char *, int , int,
			    int , int ,  int , 
			    double , int , int,
			   int , int ,  int , 
			   int , int , int , 
			   int , int , int,
			   int , int , int , 
			   int , int , int , 
			   int, int , int , 
			   int );


/****************************************************/

#endif