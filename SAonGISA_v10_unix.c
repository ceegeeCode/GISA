/* **********************************************************************************************************************************
*
*
* ***************************************   Strucutral Analysis based on GISA  (SAonGISA) *******************************************************************************
*
*
*
*     Present version: version _v10, > 1st Jan, 2019.
*     Previous versions: version _v3 _v4, > 19th Jan, 2017; version _v5 _v6, > 5th Apr, 2017; version _v7, > 12th Jul, 2017: version _v8, > 17th Nov, 2017, _v9, > 17th Jan, 2018.
*
*     New in version _v9: 
*     * Now possible to include all order 2 invariants in the matching (ie all 12 "full" order 2's, abs value versions incl'ed) 
*
**************************************************************************************************************************************
* 
*     C-implementation of algorithm for computating Gauss integral based invariants for protein fold desciption. Code aimed for 
*     searching for "rare geometries" is included.
*
*     Author: Christian Grønbæk
*
*     If using this code cite this source:
*     C.Grønbæk., T.Hamelryck, P.Røgen: "GISA: Using Gauss integrals to identify rare conformations in 
*     protein structures", TO APPEAR.
*     
*     Notes:
*     !See ReadMe, code outlines and example runs on www.github.com, repository = ceeGeeCode/GISA !
*
*	  OBS: pointers are allocated memory with explicit cast for the sake of avoiding troubles with missing 
*     casts if code is augmented with e.g. CUDA C for gpu-computing. 
*
**********************************************************************************************/

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
			   int, int , int , int );

/****************************************************
*****************Coding section *********************
*****************************************************/


int * readDBresultsWindowsToPtr(struct db_I_windows_order2_ptr **ptr_ptr_dbResults, char DBResultsFilePath[1000], int chunkSize, int *ptr_maxNrOfWindows ){

	int *out;

	/*pointer to hold DB results*/
	struct db_I_windows_order2_ptr *ptr_dbResults;

	/*for holding means and std devs (if desired)*/
	struct I_windows_order2_meanStddev I_windows_order2_meanStddev_out; 

	FILE *ptr_fileDB;

	int returnVal;

	/*for reading DB files:*/
	char line[1000];
	char *linePart;
	char linePartFileName[200];
	char linePartChainId[10];
	int idxInLine = -1;
	int idx_fileName = -1;
	int idx_strName = -1;
	int idx_classId = -1;
	int idx_chainId = -1;
	int idx_chainNr = -1; 
	int idx_chainLength = -1;
	int idx_nrOfWindows = -1;
	int idx_windowNr = -1; 
	int winNr = -1;
	int idx_i1 = -1;
	int idx_i2 = -1;
	int idx_I12 = -1;
	int idx_Ia12 = -1;
	int idx_I1234_full = -1;
	int idx_I1324_full = -1;
	int idx_I1423 = -1;

	int idx_Ia12a34_full = -1;
	int idx_Ia13a24_full = -1;
	int idx_Ia14a23 = -1;
	
	int idx_Ia1234_full = -1;
	int idx_Ia1324_full = -1;
	int idx_Ia1423 = -1;
	
	int idx_I12a34_full = -1;
	int idx_I13a24_full = -1;
	int idx_I14a23 = -1;

	int i_passedHeadingLines = 0;
	int i_passedValueLines = 0;
	int i_matchLineCnt = 0;

	int recNr = 0; //record number for each record in DB (ie for each loaded structure)
	int chunkRecNr; //cnt number of records within each chunk (when loading)
	int chunkNr; //cnt number of chunks
	long int currFilePos; //integer giving the current position in file (when reading)
	int nrOfWindows = 0;
	int maxNrOfWindows = 0;
	int rewSuccess = 1;
	int totalNrOfWins = 0;

	/*alloc and init the output ptr*/
	out = (int *) malloc(2*sizeof(int));
	out[0] = 0;
	out[1] = 1;

	/*Open the file containing the GIs on windows/pairs of the data base structures:*/
	ptr_fileDB = fopen(DBResultsFilePath, "r");


	if(ptr_fileDB == NULL){
		printf("Sorry -- cannot find the file: %s\n", DBResultsFilePath);
		printf("Press enter if you want to continue.\n");
		getchar();

	}

	/*Read first entry of DB file to record its structure*/
	while (fgets(line,1000, ptr_fileDB)!=NULL){

		idxInLine = 0;
		if (strstr(line, "Line type") == line){

			/*loop trough the line and record important indices:*/
			linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
			while(linePart)
			{
				if (strcmp(linePart, "PDBfile")==0){					
					idx_fileName = idxInLine;
				}

				if (strcmp(linePart, "structureId")==0){					
					idx_strName = idxInLine;
				}

				if (strcmp(linePart, "classId")==0){					
					idx_classId = idxInLine;
				}

				if (strcmp(linePart, "chainId")==0){					
					idx_chainId = idxInLine;
				}

				if (strcmp(linePart, "chainNr")==0){					
					idx_chainNr = idxInLine;
				}

				if (strcmp(linePart, "chainLength")==0){					
					idx_chainLength = idxInLine;
				}

				if (strcmp(linePart, "nrOfWindows")==0){					
					idx_nrOfWindows = idxInLine;
				}

				linePart=strtok(NULL, ";"); //this gets the next token of the string (which is line)
				idxInLine += 1;
			}
		}

		else {

			idxInLine = -1; //start at -1 since the first token in lines starting with the string "Values" is, of course, always "Values"; the vlaues line following such a line however have no placeholder for this word (ie it ois displaced by 1 rel to the "Values"-starting line) 

			if (strstr(line, "Values") == line){

				/*loop trough the line and record important indices:*/
				linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
				while(linePart)
				{

					if (strcmp(linePart, "windowNr") == 0){					
						idx_windowNr = idxInLine;
					}
					if (strcmp(linePart,"i1") == 0){					
						idx_i1 = idxInLine;
					}
					if (strcmp(linePart,"i2") == 0){					
						idx_i2 = idxInLine;
					}
					if (strcmp(linePart,"I12") == 0){					
						idx_I12 = idxInLine;
					}
					if (strcmp(linePart,"Ia12") == 0){					
						idx_Ia12 = idxInLine;
					}
					if (strcmp(linePart,"I1234_full") == 0){					
						idx_I1234_full = idxInLine;
					}
					if (strcmp(linePart,"Ia12a34_full") == 0){					
						idx_Ia12a34_full = idxInLine;
					}
					if (strcmp(linePart,"Ia1234_full") == 0){					
						idx_Ia1234_full = idxInLine;
					}
					if (strcmp(linePart,"I12a34_full") == 0){					
						idx_I12a34_full = idxInLine;
					}

					if (strcmp(linePart,"I1324_full") == 0){					
						idx_I1324_full = idxInLine;
					}
					if (strcmp(linePart,"Ia13a24_full") == 0){					
						idx_Ia13a24_full = idxInLine;
					}
					if (strcmp(linePart,"Ia1324_full") == 0){					
						idx_Ia1324_full = idxInLine;
					}
					if (strcmp(linePart,"I13a24_full") == 0){					
						idx_I13a24_full = idxInLine;
					}
					
					if (strcmp(linePart,"I1423") == 0){					
						idx_I1423 = idxInLine;
					}
					if (strcmp(linePart,"Ia14a23") == 0){					
						idx_Ia14a23 = idxInLine;
					}
					if (strcmp(linePart,"Ia1423") == 0){					
						idx_Ia1423 = idxInLine;
					}
					if (strcmp(linePart,"I14a23") == 0){					
						idx_I14a23 = idxInLine;
					}

					linePart=strtok(NULL, ";"); //this gets the next token of the string (which is line)
					idxInLine += 1;
				}

				break;
			}

		}

	}

	/*get back to first line of file*/
	rewind(ptr_fileDB);

	/*loop through DB and keep the records in the pointer declared for the purpose;
	memory is allocated chunk-wise. This is handled as follows: in an outer while-loop
	the file's records are traversed; inside this two while-loops are placed, each 
	will traverse a "current chunk" once: in the first the highest nr of windows used
	for covering the structure within a chunk is found; this nr is used to re-allocate
	memory to the ptr holding the DB values; loading these into the ptr is done in the 
	second inner while-loop. For handling each chunk the first inner while-loop traverses
	it, upon exiting this loop the ptr to the file is reset to the beginning of the chunk
	and the records are then loaded in in the next inner while-loop */

	recNr = 0; //init
	chunkRecNr = 0; //init
	chunkNr = 0; //init
	i_passedHeadingLines = 0; //reset
	while(i_passedHeadingLines == 0){ //outer while-loop
		
		currFilePos = ftell(ptr_fileDB); //this will be executed every time a new chunk is started at
		
		chunkRecNr = 0; //reset
		maxNrOfWindows = 0; //reset

		printf("chunk: %d\n", chunkNr);

		while (fgets(line,1000, ptr_fileDB)!=NULL){ //1st inner while-loop
			
			/*exit the inner while-loop for allocating memory and recording the DB values  */
			if (strstr(line, "Line type") == line && i_passedHeadingLines != 0 && chunkRecNr == chunkSize){
				
				//printf("exit pt, chunk Nr: %d\n", chunkNr);
				//getchar();
				i_passedHeadingLines = 0;
				break;
			}
			else{
				if (strstr(line, "Header") == line){

					//printf("line: %s\n", line);

					chunkRecNr += 1;
					i_passedHeadingLines = 1;

					linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
					idxInLine = 0;
					//find the nr of windows for this chain/structure:
					while(linePart){

						if (idxInLine == idx_nrOfWindows){nrOfWindows = atoi(linePart);}

						if (nrOfWindows > maxNrOfWindows){maxNrOfWindows = nrOfWindows;}

						linePart = strtok(NULL, ";"); //this gets the next token of the string (which is line)
						idxInLine += 1;
					}

					//printf("nr of winds: %d max so far: %d\n", nrOfWindows, maxNrOfWindows );
				}
				
			}

		} //1st inner while-loop 

		if(maxNrOfWindows > *ptr_maxNrOfWindows){
			*ptr_maxNrOfWindows = maxNrOfWindows;
		}

		//printf("maxNrOfWindows: %d\n", *ptr_maxNrOfWindows);

		/*Now load the I values into ptr's declared for the purpose;
		we first allocate memory:*/
		returnVal = alloc_init_DBresultsWindows(&ptr_dbResults, chunkNr, chunkSize, maxNrOfWindows);

		/*for(m=0;m<chunkSize; m++){
			for(i=0; i< nrOfWindows; i++){

				printf("at rec %d and win %d seg inds are: %d %d\n", m, i, ptr_dbResults[m].segIndices[i][0], ptr_dbResults[m].segIndices[i][1]); 

			}
		}*/

		//rewind to the position at the beginning of the current chunk:
		rewSuccess = fseek(ptr_fileDB, currFilePos, SEEK_SET);
		//carry on with the 2nd inner while-loop, and load the records into the pointer
		//In case the file contains less records than the set chunkSize or if the last piece
		//of a file has been reached, we may have chunkRecNr < chunkSize and still have reaced
		//this spot. Therefore we reest also (in case chunkRecNr =chunkSize this was done in 
		//if-clause above):
		i_passedHeadingLines = 0;
		chunkRecNr = 0; //reset
		while (fgets(line,1000, ptr_fileDB)!=NULL){ //2nd inner while-loop
						 
			//printf("line %s", line);

			/*exit the inner while-loop for storing the record*/
			if (strstr(line, "Line type") == line && i_passedHeadingLines != 0 ){

				//increments and resets
				recNr += 1;
				chunkRecNr += 1;
				//i_matchLineCnt = 0;
				i_passedHeadingLines = 0;
				//printf("Rec nr: %d\n", recNr);
				//printf("Chunk rec nr: %d\n", chunkRecNr);
				//getchar();
				if(chunkRecNr == chunkSize){
					break;
				}


			}
			else{
				if (strstr(line, "Header") == line){

					linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
					idxInLine = 0;
					while(linePart){

						if (idxInLine == idx_fileName){
							//strcpy(linePartFileName,linePart); //have to strcpy via a "fixed" string (char array)
							//printf("file %s", linePart);
							strcpy(ptr_dbResults[recNr].fileName, linePart);
						}
						if (idxInLine == idx_strName){
							//strcpy(linePartFileName,linePart); //have to strcpy via a "fixed" string (char array)
							//printf("str %s", linePart);
							strcpy(ptr_dbResults[recNr].structureName, linePart);
						}
						if (idxInLine == idx_classId){ 
							//strcpy(linePartChainId,linePart); //have to strcpy via a "fixed" string (char array)
							//printf("class %s", linePart);
							strcpy(ptr_dbResults[recNr].classId, linePart);
						}
						if (idxInLine == idx_chainId){ 
							//strcpy(linePartChainId,linePart); //have to strcpy via a "fixed" string (char array)
							//printf("chain %s", linePart);
							strcpy(ptr_dbResults[recNr].chainId, linePart);
						}
						if (idxInLine == idx_chainNr){
							ptr_dbResults[recNr].chainNr = atoi(linePart);
						}
						if (idxInLine == idx_chainLength){
							ptr_dbResults[recNr].chainLen = atoi(linePart);
						}
						if (idxInLine == idx_nrOfWindows){
							ptr_dbResults[recNr].nrOfWindows = atoi(linePart);
						}

						linePart = strtok(NULL, ";"); //this gets the next token of the string (which is line)
						idxInLine += 1;
					}
					//printf("file nr %d loaded: %s\n", recNr, ptr_dbResults[recNr].fileName);
					//printf(".....%s\n", ptr_dbResults[recNr].chainId);
					//printf(".....%d\n", ptr_dbResults[recNr].nrOfWindows);
				}
				else{

					if (strstr(line, "Values") == line){

						i_passedHeadingLines = 1;
						/*recNr += 1;
						chunkRecNr += 1;
						printf("Rec nr: %d\n", recNr);
						printf("Chunk rec nr: %d\n", chunkRecNr);*/
						i_matchLineCnt = 0;  //reset

					}
					else{

						if (i_passedHeadingLines == 1){

							linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
							idxInLine = 0;
							while(linePart){

								if (idxInLine == idx_windowNr){
									winNr = atoi(linePart);
									ptr_dbResults[recNr].window[winNr].windowNr = winNr;
									//ptr_dbResults[recNr].windowsNr[winNr] = winNr;
									//printf("w-nr .. %d\n", ptr_dbResults[recNr].windowsNr[i_matchLineCnt]); 
									totalNrOfWins +=1;
								}
								if (idxInLine == idx_i1){
									ptr_dbResults[recNr].window[winNr].segIndices[0] = atoi(linePart);
									//ptr_dbResults[recNr].segIndices[winNr][0] = atoi(linePart);
									//printf("i1 .. %d\n", ptr_dbResults[recNr].segIndices[i_matchLineCnt][0]); 
								}
								if (idxInLine == idx_i2){
									//printf("part ..:%s", linePart);
									ptr_dbResults[recNr].window[winNr].segIndices[1] = atoi(linePart);
									//ptr_dbResults[recNr].segIndices[winNr][1] = atoi(linePart);
									//printf("i2 .. %d\n", ptr_dbResults[recNr].segIndices[i_matchLineCnt][1]); 
								}
								if (idxInLine == idx_I12){ptr_dbResults[recNr].I12[winNr] = atof(linePart);}
								if (idxInLine == idx_Ia12){ptr_dbResults[recNr].Ia12[winNr] = atof(linePart);}
								if (idxInLine == idx_I1234_full){ptr_dbResults[recNr].I1234_full[winNr] = atof(linePart);}
								if (idxInLine == idx_I1324_full){ptr_dbResults[recNr].I1324_full[winNr] = atof(linePart);}
								if (idxInLine == idx_I1423){ptr_dbResults[recNr].I1423[winNr] = atof(linePart);}

								if (idxInLine == idx_Ia12a34_full){ptr_dbResults[recNr].Ia12a34_full[winNr] = atof(linePart);}
								if (idxInLine == idx_Ia13a24_full){ptr_dbResults[recNr].Ia13a24_full[winNr] = atof(linePart);}
								if (idxInLine == idx_Ia14a23){ptr_dbResults[recNr].Ia14a23[winNr] = atof(linePart);}

								if (idxInLine == idx_Ia1234_full){ptr_dbResults[recNr].Ia1234_full[winNr] = atof(linePart);}
								if (idxInLine == idx_Ia1324_full){ptr_dbResults[recNr].Ia1324_full[winNr] = atof(linePart);}
								if (idxInLine == idx_Ia1423){ptr_dbResults[recNr].Ia1423[winNr] = atof(linePart);}

								if (idxInLine == idx_I12a34_full){ptr_dbResults[recNr].I12a34_full[winNr] = atof(linePart);}
								if (idxInLine == idx_I13a24_full){ptr_dbResults[recNr].I13a24_full[winNr] = atof(linePart);}
								if (idxInLine == idx_I14a23){ptr_dbResults[recNr].I14a23[winNr] = atof(linePart);}
								
								linePart = strtok(NULL, ";"); //this gets the next token of the string (which is line)
								idxInLine += 1;
							}

							i_matchLineCnt +=1;

						}
					}
				}
			}

		} //2nd inner while-loop

		chunkNr += 1;
		//printf("chunkNr: %d\n", chunkNr);
	}

	/*for(i=0; i<13;i++){
		printf("ptr_dbResults filename at %d: %s\n", i, ptr_dbResults[i].fileName);
		printf("............. nrOfWindows at %d: %d\n", i, ptr_dbResults[i].nrOfWindows);
	}*/

	//pass the ptr_dbResults to in/out variable:
	*ptr_ptr_dbResults = ptr_dbResults;

	printf("Nr of records loaded: %d\n", recNr +1); 

	fclose(ptr_fileDB);

	out[0] =  recNr +1; //return the nr of records loaded (starts at idx 0 so we add 1) 
	out[1] = totalNrOfWins;
	

	return out;

}

int * loadDBWindowsFromDBPairs(struct db_I_windows_order2_ptr **ptr_ptr_dbResults, struct db_I_windowPairs_order2_ptr *ptr_dbResultsPairs, int nrOfDBstructures_pairsLoad, int chunkSize){

	int *out;

	int returnVal;

	/*pointer to hold DB results*/
	struct db_I_windows_order2_ptr *ptr_dbResults;

	int recNr = 0; //record number for each record in DB (ie for each loaded structure)
	int chunkRecNr = 0; //cnt number of records within each chunk (when loading)
	int chunkNr = 0; //cnt number of chunks
	int maxNrOfWindows = 0;
	int totalNrOfWins = 0;

	int m = 0; //m will be the "record nr"
	int i = 0, j = 0;

	/*alloc and init the output ptr*/
	out = (int *) malloc(2*sizeof(int));
	out[0] = 0;
	out[1] = 1;


	for(m = 0; m < nrOfDBstructures_pairsLoad; m++){

		if(m%chunkSize*chunkNr == 0){

			/*find the max nr of windows needed for the next chunk:*/
			maxNrOfWindows = 0; //init/reset
			for(j = m; j < min(m + chunkSize, nrOfDBstructures_pairsLoad); j++){

				if(maxNrOfWindows < ptr_dbResultsPairs[j].nrOfWindows){
					maxNrOfWindows = ptr_dbResultsPairs[j].nrOfWindows;
				}

				printf("Str name %s maxNrOfWindows %d\n", ptr_dbResultsPairs[j].structureName, maxNrOfWindows);

			}

			printf("maxNrOfWindows: %d\n", maxNrOfWindows);

			/*the (re)allocate:*/
			returnVal = alloc_init_DBresultsWindows(&ptr_dbResults, chunkNr, chunkSize, maxNrOfWindows);

			/*increment the chunkNr (for the next chunk); this will also count the nr of chunks allocated for:*/
			chunkNr += 1;
		}

		/*Now load the pairs data into the single windows ptr:*/ 
		//printf("In  load of pairs to singel wins: chain len from pairs: %d\n", ptr_dbResultsPairs[m].chainLen);
		strcpy(ptr_dbResults[m].fileName, ptr_dbResultsPairs[m].fileName);
		strcpy(ptr_dbResults[m].structureName, ptr_dbResultsPairs[m].structureName);
		strcpy(ptr_dbResults[m].chainId, ptr_dbResultsPairs[m].chainId);
		strcpy(ptr_dbResults[m].classId, ptr_dbResultsPairs[m].classId);
		ptr_dbResults[m].chainLen = ptr_dbResultsPairs[m].chainLen;
		ptr_dbResults[m].chainNr = ptr_dbResultsPairs[m].chainNr;
		ptr_dbResults[m].nrOfWindows =  ptr_dbResultsPairs[m].nrOfWindows;

		for(i = 0; i< ptr_dbResultsPairs[m].nrOfWindows;i++){

			ptr_dbResults[m].window[i].windowNr =  ptr_dbResultsPairs[m].windowPair[i][i].windowNr_1;
			ptr_dbResults[m].window[i].segIndices[0] = ptr_dbResultsPairs[m].windowPair[i][i].segIndices_1[0];
			ptr_dbResults[m].window[i].segIndices[1] = ptr_dbResultsPairs[m].windowPair[i][i].segIndices_1[1];
			/*ptr_dbResults[m].windowsNr[i] =  ptr_dbResultsPairs[m].windowsNr_1[i];
			ptr_dbResults[m].segIndices[i][0] = ptr_dbResultsPairs[m].segIndices_1[i][0];
			ptr_dbResults[m].segIndices[i][1] = ptr_dbResultsPairs[m].segIndices_1[i][1];*/
			ptr_dbResults[m].I12[i] = ptr_dbResultsPairs[m].I12[i][i];
			ptr_dbResults[m].Ia12[i] = ptr_dbResultsPairs[m].Ia12[i][i];
			ptr_dbResults[m].I1234_full[i] = ptr_dbResultsPairs[m].I1234_full[i][i];
			ptr_dbResults[m].I1324_full[i] = ptr_dbResultsPairs[m].I1324_full[i][i];
			ptr_dbResults[m].I1423[i] = ptr_dbResultsPairs[m].I1423[i][i];

			ptr_dbResults[m].Ia12a34_full[i] = ptr_dbResultsPairs[m].Ia12a34_full[i][i];
			ptr_dbResults[m].Ia13a24_full[i] = ptr_dbResultsPairs[m].Ia13a24_full[i][i];
			ptr_dbResults[m].Ia14a23[i] = ptr_dbResultsPairs[m].Ia14a23[i][i];

			ptr_dbResults[m].Ia1234_full[i] = ptr_dbResultsPairs[m].Ia1234_full[i][i];
			ptr_dbResults[m].Ia1324_full[i] = ptr_dbResultsPairs[m].Ia1324_full[i][i];
			ptr_dbResults[m].Ia1423[i] = ptr_dbResultsPairs[m].Ia1423[i][i];

			ptr_dbResults[m].I12a34_full[i] = ptr_dbResultsPairs[m].I12a34_full[i][i];
			ptr_dbResults[m].I13a24_full[i] = ptr_dbResultsPairs[m].I13a24_full[i][i];
			ptr_dbResults[m].I14a23[i] = ptr_dbResultsPairs[m].I14a23[i][i];
		}

		
		totalNrOfWins += ptr_dbResultsPairs[m].nrOfWindows;
		
	}


	//pass the ptr_dbResults to in/out variable:
	*ptr_ptr_dbResults = ptr_dbResults;

	printf("Nr of single window records (structures) populated from pairs results: %d\n", m ); 

	out[0] = m; //return the nr of records loaded (starts at idx 0 so we add 1) 
	out[1] = totalNrOfWins;
	
	return out;

}

struct I_windows_order2_meanStddev normalizeDBresultsWindowsPtr(struct db_I_windows_order2_ptr *ptr_dbResults, int nrOfRecords ){

	int n = 0;
	int i = 0;

	int totalNrOfWindows = nrOfRecords*100;
	
	int cnt  = 0;

	struct I_windows_order2_raw_ptr I_windows_order2_raw_ptr;
	struct I_windows_order2_meanStddev I_windows_order2_meanStddev_out;

	struct double2 meanStddevOut;

	/*allocate memory*/
	I_windows_order2_raw_ptr.I12 = (double *) malloc(totalNrOfWindows*sizeof(double));
	I_windows_order2_raw_ptr.Ia12 = (double *) malloc(totalNrOfWindows*sizeof(double));

	I_windows_order2_raw_ptr.I1234_full = (double *) malloc(totalNrOfWindows*sizeof(double));
	I_windows_order2_raw_ptr.I1324_full = (double *) malloc(totalNrOfWindows*sizeof(double));
	I_windows_order2_raw_ptr.I1423 = (double *) malloc(totalNrOfWindows*sizeof(double));

	I_windows_order2_raw_ptr.Ia12a34_full = (double *) malloc(totalNrOfWindows*sizeof(double));
	I_windows_order2_raw_ptr.Ia13a24_full = (double *) malloc(totalNrOfWindows*sizeof(double));
	I_windows_order2_raw_ptr.Ia14a23 = (double *) malloc(totalNrOfWindows*sizeof(double));

	I_windows_order2_raw_ptr.Ia1234_full = (double *) malloc(totalNrOfWindows*sizeof(double));
	I_windows_order2_raw_ptr.Ia1324_full = (double *) malloc(totalNrOfWindows*sizeof(double));
	I_windows_order2_raw_ptr.Ia1423 = (double *) malloc(totalNrOfWindows*sizeof(double));

	I_windows_order2_raw_ptr.I12a34_full = (double *) malloc(totalNrOfWindows*sizeof(double));
	I_windows_order2_raw_ptr.I13a24_full = (double *) malloc(totalNrOfWindows*sizeof(double));
	I_windows_order2_raw_ptr.I14a23 = (double *) malloc(totalNrOfWindows*sizeof(double));

	/*init: read in all values to ptr*/
	for(n=0;n<nrOfRecords;n++){

		/*reallocate memory if necessary*/
		if(cnt > totalNrOfWindows){

			totalNrOfWindows = floor(totalNrOfWindows*1.25);
			I_windows_order2_raw_ptr.I12 = (double *) realloc(I_windows_order2_raw_ptr.I12,  totalNrOfWindows*sizeof(double));
			I_windows_order2_raw_ptr.Ia12 = (double *) realloc(I_windows_order2_raw_ptr.Ia12,  totalNrOfWindows*sizeof(double));

			I_windows_order2_raw_ptr.I1234_full = (double *) realloc(I_windows_order2_raw_ptr.I1234_full,  totalNrOfWindows*sizeof(double));
			I_windows_order2_raw_ptr.I1324_full = (double *) realloc(I_windows_order2_raw_ptr.I1324_full,  totalNrOfWindows*sizeof(double));
			I_windows_order2_raw_ptr.I1423 = (double *) realloc(I_windows_order2_raw_ptr.I1423,  totalNrOfWindows*sizeof(double));

			I_windows_order2_raw_ptr.Ia12a34_full = (double *) realloc(I_windows_order2_raw_ptr.Ia12a34_full,  totalNrOfWindows*sizeof(double));
			I_windows_order2_raw_ptr.Ia13a24_full = (double *) realloc(I_windows_order2_raw_ptr.Ia13a24_full,  totalNrOfWindows*sizeof(double));
			I_windows_order2_raw_ptr.Ia14a23 = (double *) realloc(I_windows_order2_raw_ptr.Ia14a23,  totalNrOfWindows*sizeof(double));

			I_windows_order2_raw_ptr.Ia1234_full = (double *) realloc(I_windows_order2_raw_ptr.Ia1234_full,  totalNrOfWindows*sizeof(double));
			I_windows_order2_raw_ptr.Ia1324_full = (double *) realloc(I_windows_order2_raw_ptr.Ia1324_full,  totalNrOfWindows*sizeof(double));
			I_windows_order2_raw_ptr.Ia1423 = (double *) realloc(I_windows_order2_raw_ptr.Ia1423,  totalNrOfWindows*sizeof(double));

			I_windows_order2_raw_ptr.I12a34_full = (double *) realloc(I_windows_order2_raw_ptr.I12a34_full,  totalNrOfWindows*sizeof(double));
			I_windows_order2_raw_ptr.I13a24_full = (double *) realloc(I_windows_order2_raw_ptr.I13a24_full,  totalNrOfWindows*sizeof(double));
			I_windows_order2_raw_ptr.I14a23 = (double *) realloc(I_windows_order2_raw_ptr.I14a23,  totalNrOfWindows*sizeof(double));

		}

		/*??? before 27 March 2017 this was: for(i= cnt; i< ptr_dbResults[n].nrOfWindows; i++){*/
		for(i= 0; i< ptr_dbResults[n].nrOfWindows; i++){

			I_windows_order2_raw_ptr.I12[cnt] = ptr_dbResults[n].I12[i];
			I_windows_order2_raw_ptr.Ia12[cnt] = ptr_dbResults[n].Ia12[i];
			
			I_windows_order2_raw_ptr.I1234_full[cnt] = ptr_dbResults[n].I1234_full[i];
			I_windows_order2_raw_ptr.I1324_full[cnt] = ptr_dbResults[n].I1324_full[i];
			I_windows_order2_raw_ptr.I1423[cnt] = ptr_dbResults[n].I1423[i];

			I_windows_order2_raw_ptr.Ia12a34_full[cnt] = ptr_dbResults[n].Ia12a34_full[i];
			I_windows_order2_raw_ptr.Ia13a24_full[cnt] = ptr_dbResults[n].Ia13a24_full[i];
			I_windows_order2_raw_ptr.Ia14a23[cnt] = ptr_dbResults[n].Ia14a23[i];

			I_windows_order2_raw_ptr.Ia1234_full[cnt] = ptr_dbResults[n].Ia1234_full[i];
			I_windows_order2_raw_ptr.Ia1324_full[cnt] = ptr_dbResults[n].Ia1324_full[i];
			I_windows_order2_raw_ptr.Ia1423[cnt] = ptr_dbResults[n].Ia1423[i];

			I_windows_order2_raw_ptr.I12a34_full[cnt] = ptr_dbResults[n].I12a34_full[i];
			I_windows_order2_raw_ptr.I13a24_full[cnt] = ptr_dbResults[n].I13a24_full[i];
			I_windows_order2_raw_ptr.I14a23[cnt] = ptr_dbResults[n].I14a23[i];

			cnt +=1;
		}

	}

	/*compute mean and std dev's:*/
	meanStddevOut = meanStddev(I_windows_order2_raw_ptr.I12, cnt);
	I_windows_order2_meanStddev_out.mean_I12 = meanStddevOut.val[0];
	I_windows_order2_meanStddev_out.stddev_I12 = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windows_order2_raw_ptr.Ia12, cnt);
	I_windows_order2_meanStddev_out.mean_Ia12 = meanStddevOut.val[0];
	I_windows_order2_meanStddev_out.stddev_Ia12 = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windows_order2_raw_ptr.I1234_full, cnt);
	I_windows_order2_meanStddev_out.mean_I1234_full = meanStddevOut.val[0];
	I_windows_order2_meanStddev_out.stddev_I1234_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windows_order2_raw_ptr.I1324_full, cnt);
	I_windows_order2_meanStddev_out.mean_I1324_full = meanStddevOut.val[0];
	I_windows_order2_meanStddev_out.stddev_I1324_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windows_order2_raw_ptr.I1423, cnt);
	I_windows_order2_meanStddev_out.mean_I1423 = meanStddevOut.val[0];
	I_windows_order2_meanStddev_out.stddev_I1423 = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windows_order2_raw_ptr.Ia12a34_full, cnt);
	I_windows_order2_meanStddev_out.mean_Ia12a34_full = meanStddevOut.val[0];
	I_windows_order2_meanStddev_out.stddev_Ia12a34_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windows_order2_raw_ptr.Ia13a24_full, cnt);
	I_windows_order2_meanStddev_out.mean_Ia13a24_full = meanStddevOut.val[0];
	I_windows_order2_meanStddev_out.stddev_Ia13a24_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windows_order2_raw_ptr.Ia14a23, cnt);
	I_windows_order2_meanStddev_out.mean_Ia14a23 = meanStddevOut.val[0];
	I_windows_order2_meanStddev_out.stddev_Ia14a23 = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windows_order2_raw_ptr.Ia1234_full, cnt);
	I_windows_order2_meanStddev_out.mean_Ia1234_full = meanStddevOut.val[0];
	I_windows_order2_meanStddev_out.stddev_Ia1234_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windows_order2_raw_ptr.Ia1324_full, cnt);
	I_windows_order2_meanStddev_out.mean_Ia1324_full = meanStddevOut.val[0];
	I_windows_order2_meanStddev_out.stddev_Ia1324_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windows_order2_raw_ptr.Ia1423, cnt);
	I_windows_order2_meanStddev_out.mean_Ia1423 = meanStddevOut.val[0];
	I_windows_order2_meanStddev_out.stddev_Ia1423 = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windows_order2_raw_ptr.I12a34_full, cnt);
	I_windows_order2_meanStddev_out.mean_I12a34_full = meanStddevOut.val[0];
	I_windows_order2_meanStddev_out.stddev_I12a34_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windows_order2_raw_ptr.I13a24_full, cnt);
	I_windows_order2_meanStddev_out.mean_I13a24_full = meanStddevOut.val[0];
	I_windows_order2_meanStddev_out.stddev_I13a24_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windows_order2_raw_ptr.I14a23, cnt);
	I_windows_order2_meanStddev_out.mean_I14a23 = meanStddevOut.val[0];
	I_windows_order2_meanStddev_out.stddev_I14a23 = meanStddevOut.val[1];

	/*Now scale the input by the std-dev*/
	for(n=0;n<nrOfRecords;n++){

		normalizeArray(I_windows_order2_meanStddev_out.mean_I12, I_windows_order2_meanStddev_out.stddev_I12, &ptr_dbResults[n].I12, ptr_dbResults[n].nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_out.mean_Ia12, I_windows_order2_meanStddev_out.stddev_Ia12, &ptr_dbResults[n].Ia12, ptr_dbResults[n].nrOfWindows);
		
		normalizeArray(I_windows_order2_meanStddev_out.mean_I1234_full, I_windows_order2_meanStddev_out.stddev_I1234_full, &ptr_dbResults[n].I1234_full, ptr_dbResults[n].nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_out.mean_I1324_full, I_windows_order2_meanStddev_out.stddev_I1324_full, &ptr_dbResults[n].I1324_full, ptr_dbResults[n].nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_out.mean_I1423, I_windows_order2_meanStddev_out.stddev_I1423, &ptr_dbResults[n].I1423, ptr_dbResults[n].nrOfWindows);

		normalizeArray(I_windows_order2_meanStddev_out.mean_Ia12a34_full, I_windows_order2_meanStddev_out.stddev_Ia12a34_full, &ptr_dbResults[n].Ia12a34_full, ptr_dbResults[n].nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_out.mean_Ia13a24_full, I_windows_order2_meanStddev_out.stddev_Ia13a24_full, &ptr_dbResults[n].Ia13a24_full, ptr_dbResults[n].nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_out.mean_Ia14a23, I_windows_order2_meanStddev_out.stddev_Ia14a23, &ptr_dbResults[n].Ia14a23, ptr_dbResults[n].nrOfWindows);

		normalizeArray(I_windows_order2_meanStddev_out.mean_Ia1234_full, I_windows_order2_meanStddev_out.stddev_Ia1234_full, &ptr_dbResults[n].Ia1234_full, ptr_dbResults[n].nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_out.mean_Ia1324_full, I_windows_order2_meanStddev_out.stddev_Ia1324_full, &ptr_dbResults[n].Ia1324_full, ptr_dbResults[n].nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_out.mean_Ia1423, I_windows_order2_meanStddev_out.stddev_Ia1423, &ptr_dbResults[n].Ia1423, ptr_dbResults[n].nrOfWindows);

		normalizeArray(I_windows_order2_meanStddev_out.mean_I12a34_full, I_windows_order2_meanStddev_out.stddev_I12a34_full, &ptr_dbResults[n].I12a34_full, ptr_dbResults[n].nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_out.mean_I13a24_full, I_windows_order2_meanStddev_out.stddev_I13a24_full, &ptr_dbResults[n].I13a24_full, ptr_dbResults[n].nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_out.mean_I14a23, I_windows_order2_meanStddev_out.stddev_I14a23, &ptr_dbResults[n].I14a23, ptr_dbResults[n].nrOfWindows);

	}

	/*free memory*/
	free(I_windows_order2_raw_ptr.I12);
	free(I_windows_order2_raw_ptr.Ia12);

	free(I_windows_order2_raw_ptr.I1234_full);
	free(I_windows_order2_raw_ptr.I1324_full);
	free(I_windows_order2_raw_ptr.I1423);

	free(I_windows_order2_raw_ptr.Ia12a34_full);
	free(I_windows_order2_raw_ptr.Ia13a24_full);
	free(I_windows_order2_raw_ptr.Ia14a23);

	free(I_windows_order2_raw_ptr.Ia1234_full);
	free(I_windows_order2_raw_ptr.Ia1324_full);
	free(I_windows_order2_raw_ptr.Ia1423);

	free(I_windows_order2_raw_ptr.I12a34_full);
	free(I_windows_order2_raw_ptr.I13a24_full);
	free(I_windows_order2_raw_ptr.I14a23);


	return I_windows_order2_meanStddev_out;

}

void normalizeQueryWindowsPtr(struct I_windows_ptr I_windows, struct I_windows_order2_meanStddev I_windows_order2_meanStddev_in, int nrOfInvs){

	/*normalize using the input means and std devs*/
	normalizeArray(I_windows_order2_meanStddev_in.mean_I12, I_windows_order2_meanStddev_in.stddev_I12, &I_windows.I12, I_windows.nrOfWindows);
	normalizeArray(I_windows_order2_meanStddev_in.mean_Ia12, I_windows_order2_meanStddev_in.stddev_Ia12, &I_windows.Ia12, I_windows.nrOfWindows);

	if(nrOfInvs > 2){
		normalizeArray(I_windows_order2_meanStddev_in.mean_I1234_full, I_windows_order2_meanStddev_in.stddev_I1234_full, &I_windows.I1234_full, I_windows.nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_in.mean_I1324_full, I_windows_order2_meanStddev_in.stddev_I1324_full, &I_windows.I1324_full, I_windows.nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_in.mean_I1423, I_windows_order2_meanStddev_in.stddev_I1423, &I_windows.I1423, I_windows.nrOfWindows);

		normalizeArray(I_windows_order2_meanStddev_in.mean_Ia12a34_full, I_windows_order2_meanStddev_in.stddev_Ia12a34_full, &I_windows.Ia12a34_full, I_windows.nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_in.mean_Ia13a24_full, I_windows_order2_meanStddev_in.stddev_Ia13a24_full, &I_windows.Ia13a24_full, I_windows.nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_in.mean_Ia14a23, I_windows_order2_meanStddev_in.stddev_Ia14a23, &I_windows.Ia14a23, I_windows.nrOfWindows);

		normalizeArray(I_windows_order2_meanStddev_in.mean_Ia1234_full, I_windows_order2_meanStddev_in.stddev_Ia1234_full, &I_windows.Ia1234_full, I_windows.nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_in.mean_Ia1324_full, I_windows_order2_meanStddev_in.stddev_Ia1324_full, &I_windows.Ia1324_full, I_windows.nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_in.mean_Ia1423, I_windows_order2_meanStddev_in.stddev_Ia1423, &I_windows.Ia1423, I_windows.nrOfWindows);

		normalizeArray(I_windows_order2_meanStddev_in.mean_I12a34_full, I_windows_order2_meanStddev_in.stddev_I12a34_full, &I_windows.I12a34_full, I_windows.nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_in.mean_I13a24_full, I_windows_order2_meanStddev_in.stddev_I13a24_full, &I_windows.I13a24_full, I_windows.nrOfWindows);
		normalizeArray(I_windows_order2_meanStddev_in.mean_I14a23, I_windows_order2_meanStddev_in.stddev_I14a23, &I_windows.I14a23, I_windows.nrOfWindows);


	}

}

/*nrOfDBrecords: the number of all windows in the DB (for which there is a set of I-values), summed up over
all structures (file, chain) that is*/
int binDBresultsWindows(struct binned_I_windows **ptr_ptr_binned_I_windows, struct db_I_windows_order2_ptr *ptr_dbResults, int nrOfDBstructures, int nrOfDBrecords, int nrOfEntries, double **bins, int nrOfBins){

	int j;

	int n = 0;
	int i = 0;
	int cnt = 0;

	int returnVal;

	struct binned_I_windows *ptr_binned_I_windows; /*to be output*/

	int *binnedVector;
	double * Ivector;

	/*alloc and init*/
	returnVal = alloc_init_ptr_binned_I_windows(&ptr_binned_I_windows, nrOfDBrecords, 0, nrOfEntries);


	//binnedVector = (int *) malloc(nrOfEntries*sizeof(int));
	Ivector = (double *) malloc(nrOfEntries*sizeof(double));

	/*init*/
	for(i=0; i<nrOfEntries;i++){Ivector[i] = 0.0;}

	/*bin and read in*/
	for(n = 0; n < nrOfDBstructures; n++){

		for(i = 0; i < ptr_dbResults[n].nrOfWindows; i++){

			ptr_binned_I_windows[cnt].fileNr = n; //fileNr is actually the structure id (ie file and chain)
			//ptr_binned_I_windows[cnt].windowsNr = i;
			ptr_binned_I_windows[cnt].windowNr = ptr_dbResults[n].window[i].windowNr;
			//ptr_binned_I_windows[cnt].recKey = cnt;

			Ivector[0] =  ptr_dbResults[n].I12[i];
			if(nrOfEntries >= 2){
				Ivector[1] =  ptr_dbResults[n].Ia12[i];
			}
			if(nrOfEntries >= 3){		
				Ivector[2] =  ptr_dbResults[n].I1234_full[i];
			}
			if(nrOfEntries >= 4){	
				Ivector[3] =  ptr_dbResults[n].I1324_full[i];
			}
			if(nrOfEntries >= 5){	
				Ivector[4] =  ptr_dbResults[n].I1423[i];
			}
			if(nrOfEntries >= 6){	
				Ivector[5] =  ptr_dbResults[n].Ia12a34_full[i];
			}
			if(nrOfEntries >= 7){	
				Ivector[6] =  ptr_dbResults[n].Ia13a24_full[i];
			}
			if(nrOfEntries >= 8){	
				Ivector[7] =  ptr_dbResults[n].Ia14a23[i];
			}
			if(nrOfEntries >= 9){	
				Ivector[8] =  ptr_dbResults[n].Ia1234_full[i];
			}
			if(nrOfEntries >= 10){	
				Ivector[9] =  ptr_dbResults[n].Ia1324_full[i];
			}
			if(nrOfEntries >= 11){	
				Ivector[10] =  ptr_dbResults[n].Ia1423[i];
			}
			if(nrOfEntries >= 12){	
				Ivector[11] =  ptr_dbResults[n].I12a34_full[i];
			}
			if(nrOfEntries >= 13){	
				Ivector[12] =  ptr_dbResults[n].I13a24_full[i];
			}
			if(nrOfEntries >= 14){	
				Ivector[13] =  ptr_dbResults[n].I14a23[i];
			}


			/*bin the values*/

			returnVal = binArray(ptr_binned_I_windows[cnt].binnedIvector, Ivector, nrOfEntries, bins, nrOfBins);

			//ptr_binned_I_windows[cnt].binnedIvector = binnedVector;

			/*for(j = 0; j < nrOfEntries; j++){ printf("db-window: %d: %d 'th entry %lf placed in %d\n",i, j, Ivector[j], ptr_binned_I_windows[i].binnedIvector[j]);}
			getchar();*/

			cnt += 1;
		}
	}

	*ptr_ptr_binned_I_windows = ptr_binned_I_windows;

	return returnVal;

}

/*OBS: unlike the above binDBresultsWindows fct, the fct for binning the query's window results, binQueryResultsWindows, does not 
contain allocation to the ptr ptr_binned_I_windows; also te results are not read into a placeholder on the fly ("Ivector") but
rather provided as input (ptr_Ivector) */
int binQueryResultsWindows(struct binned_I_windows *ptr_binned_I_windows, struct I_windows_ptr I_windows, double **ptr_Ivector,  int nrOfEntries, double **bins, int nrOfBins){

	int j,k = 0;

	int returnVal = 0;

	for(k = 0; k < I_windows.nrOfWindows; k++){

		ptr_binned_I_windows[k].fileNr = -117;
		ptr_binned_I_windows[k].windowNr = I_windows.window[k].windowNr;

		/*bin the values*/

		returnVal = binArray(ptr_binned_I_windows[k].binnedIvector, ptr_Ivector[k], nrOfEntries, bins, nrOfBins);

		//ptr_binned_I_windows[cnt].binnedIvector = binnedVector;

		/*for(j = 0; j < nrOfEntries; j++){ printf("q-window: %d: %d 'th entry %lf placed in %d\n",k, j, ptr_Ivector[k][j], ptr_binned_I_windows[k].binnedIvector[j]);}
		getchar();*/
	}

	return returnVal;

}

/*Similar to binDBresultsWindows, pairs version*/
/*nrOfDBrecords: the number of all window pairs in the DB (for which there is a set of I-values), summed up over
all structures (file, chain) that is*/
int binDBresultsWindowPairs(struct binned_I_windowPairs **ptr_ptr_binned_I_windowPairs, struct db_I_windowPairs_order2_ptr *ptr_dbResultPairs, int nrOfDBstructures, int nrOfDBrecords, int nrOfEntries, double **bins, int nrOfBins, int onlyDisjointPairs_b){


	int j = 0;

	int n = 0;
	int i = 0;
	int cnt = 0;
	int readIn_b = 1;
	int k = 0;

	int returnVal;

	struct binned_I_windowPairs *ptr_binned_I_windowPairs; /*to be output*/

	int *binnedVector;
	double * Ivector;

	/*alloc and init*/
	returnVal = alloc_init_ptr_binned_I_windowPairs(&ptr_binned_I_windowPairs, nrOfDBrecords, 0, nrOfEntries);


	//binnedVector = (int *) malloc(nrOfEntries*sizeof(int));
	Ivector = (double *) malloc(nrOfEntries*sizeof(double));

	/*init*/
	for(i=0; i<nrOfEntries;i++){Ivector[i] = 0.0;}

	/*bin and read in*/
	for(n = 0; n < nrOfDBstructures; n++){

		for(i = 0; i < ptr_dbResultPairs[n].nrOfWindows; i++){

			for(j = i +1 ; j < ptr_dbResultPairs[n].nrOfWindows; j++){

				readIn_b = 1;
				
				//if(onlyDisjointPairs_b ==1 && (ptr_dbResultPairs[n].segIndices_1[i][1] > ptr_dbResultPairs[n].segIndices_2[j][0]||ptr_dbResultPairs[n].segIndices_1[i][0] < 0||ptr_dbResultPairs[n].segIndices_1[i][1] < 0)){ //the <0 parts here are due to init values of seg indices are negative; a pair may appear as disjoint just because the first pair has init value and the second hasn't)

				//	readIn_b = 0;
				//	
				//}

				if(onlyDisjointPairs_b ==1 && (ptr_dbResultPairs[n].windowPair[i][j].segIndices_1[1] > ptr_dbResultPairs[n].windowPair[i][j].segIndices_2[0]||ptr_dbResultPairs[n].windowPair[i][j].segIndices_1[0] < 0||ptr_dbResultPairs[n].windowPair[i][j].segIndices_1[1] < 0)){ //the <0 parts here are due to init values of seg indices are negative; a pair may appear as disjoint just because the first pair has init value and the second hasn't)

					readIn_b = 0;
					
				}

				if(readIn_b ==1){

					ptr_binned_I_windowPairs[cnt].fileNr = n; //fileNr is actually the structure id (ie file and chain)
					ptr_binned_I_windowPairs[cnt].windowNr_1 = ptr_dbResultPairs[n].windowPair[i][j].windowNr_1;
					ptr_binned_I_windowPairs[cnt].windowNr_2 = ptr_dbResultPairs[n].windowPair[i][j].windowNr_2;
					/*ptr_binned_I_windowPairs[cnt].windowsNr_1 = i;
					ptr_binned_I_windowPairs[cnt].windowsNr_2 = j;*/
					//ptr_binned_I_windowPairs[cnt].recKey = cnt;

					if(ptr_dbResultPairs[n].windowPair[i][j].segIndices_2[0] < 0){
						printf("(Langmodigt venter bysvalen ...) w1: %d (%d,%d) w2: %d (%d,%d)\n", i, ptr_dbResultPairs[n].windowPair[i][j].segIndices_1[0], ptr_dbResultPairs[n].windowPair[i][j].segIndices_1[1], j, ptr_dbResultPairs[n].windowPair[i][j].segIndices_2[0], ptr_dbResultPairs[n].windowPair[i][j].segIndices_2[1]);
						//getchar();
					}
					Ivector[0] =  ptr_dbResultPairs[n].I12[i][j];
					if(nrOfEntries >= 2){		
						Ivector[1] =  ptr_dbResultPairs[n].Ia12[i][j];
					}
					if(nrOfEntries >= 3){		
						Ivector[2] =  ptr_dbResultPairs[n].I1234_full[i][j];
					}
					if(nrOfEntries >= 4){	
						Ivector[3] =  ptr_dbResultPairs[n].I1324_full[i][j];
					}
					if(nrOfEntries >= 5){	
						Ivector[4] =  ptr_dbResultPairs[n].I1423[i][j];
					}
					if(nrOfEntries >= 6){	
						Ivector[5] =  ptr_dbResultPairs[n].Ia12a34_full[i][j];
					}
					if(nrOfEntries >= 7){	
						Ivector[6] =  ptr_dbResultPairs[n].Ia13a24_full[i][j];
					}
					if(nrOfEntries >= 8){	
						Ivector[7] =  ptr_dbResultPairs[n].Ia14a23[i][j];
					}
					if(nrOfEntries >= 9){	
						Ivector[8] =  ptr_dbResultPairs[n].Ia1234_full[i][j];
					}
					if(nrOfEntries >= 10){	
						Ivector[9] =  ptr_dbResultPairs[n].Ia1324_full[i][j];
					}
					if(nrOfEntries >= 11){	
						Ivector[10] =  ptr_dbResultPairs[n].Ia1423[i][j];
					}
					if(nrOfEntries >= 12){	
						Ivector[11] =  ptr_dbResultPairs[n].I12a34_full[i][j];
					}
					if(nrOfEntries >= 13){	
						Ivector[12] =  ptr_dbResultPairs[n].I13a24_full[i][j];
					}
					if(nrOfEntries >= 14){	
						Ivector[13] =  ptr_dbResultPairs[n].I14a23[i][j];
					}


					/*bin the values*/

					returnVal = binArray(ptr_binned_I_windowPairs[cnt].binnedIvector, Ivector, nrOfEntries, bins, nrOfBins);

					//ptr_binned_I_windowPairs[cnt].binnedIvector = binnedVector;

					/*for(k = 0; k < nrOfEntries; k++){ printf("db-window: %d: %d 'th entry %lf placed in %d\n",i, j, Ivector[k], ptr_binned_I_windowPairs[cnt].binnedIvector[k]);}
					getchar();*/

					cnt += 1;

					//printf("(Sommeren mild oedsles k.... ) but C the count is: %d\n", cnt); 

				}

			}
		}
	}

	*ptr_ptr_binned_I_windowPairs = ptr_binned_I_windowPairs;

	return returnVal;

}


/*Similar to binQueryResultsWindows, pairs version:*/
int binQueryResultsWindowPairs(struct binned_I_windowPairs *ptr_binned_I_windowPairs, struct I_windowPairs_ptr I_windowPairs, double ***ptr_Ivector,  int nrOfEntries, double **bins, int nrOfBins, int onlyDisjointPairs_b){

	int j, k, l = 0;

	int returnVal = 0;

	int cntPairs = 0;

	int readIn_b = 1;

	for(k = 0; k < I_windowPairs.nrOfWindows; k++){

		for(l = k + 1 ; l < I_windowPairs.nrOfWindows; l++){

			readIn_b = 1;

			if(onlyDisjointPairs_b ==1 && (I_windowPairs.windowPair[k][l].segIndices_1[1] > I_windowPairs.windowPair[k][l].segIndices_2[0]||I_windowPairs.windowPair[k][l].segIndices_1[0] < 0||I_windowPairs.windowPair[k][l].segIndices_1[1] < 0)){ //the <0 parts here are due to init values of seg indices are negative; a pair may appear as disjoint just because the first pair has init value and the second hasn't)

					readIn_b = 0;
					
			}

			if(readIn_b == 1){

				ptr_binned_I_windowPairs[cntPairs].fileNr = -117;
				ptr_binned_I_windowPairs[cntPairs].windowNr_1 = I_windowPairs.windowPair[k][l].windowNr_1;
				ptr_binned_I_windowPairs[cntPairs].windowNr_2 = I_windowPairs.windowPair[k][l].windowNr_2;
				/*ptr_binned_I_windowPairs[cntPairs].windowsNr_1 = I_windowPairs.windowsNr_1[k];
				ptr_binned_I_windowPairs[cntPairs].windowsNr_2 = I_windowPairs.windowsNr_2[l];*/
			
				/*bin the values*/

				returnVal = binArray(ptr_binned_I_windowPairs[cntPairs].binnedIvector, ptr_Ivector[k][l], nrOfEntries, bins, nrOfBins);

				//ptr_binned_I_windows[cnt].binnedIvector = binnedVector;
				//for(j = 0; j < nrOfEntries; j++){ printf("%d'th entry %lf placed in %d\n", j, ptr_Ivector[k][l][j], ptr_binned_I_windowPairs[cntPairs].binnedIvector[j]);}
				//getchar();

				cntPairs += 1;

			}
	
		}
	}

	return returnVal;

}



int * readDBresultsWindowPairsToPtr(struct db_I_windowPairs_order2_ptr **ptr_ptr_dbResults, char DBResultsFilePath[1000], int chunkSize, int *ptr_maxNrOfWindows){

	int *out;

	int returnVal = 0;

	/*pointer to hold DB results*/
	struct db_I_windowPairs_order2_ptr *ptr_dbResults;

	FILE *ptr_fileDB;

	/*for reading DB files:*/
	char line[1000];
	char *linePart;
	char linePartFileName[200];
	char linePartChainId[10];
	int idxInLine = 0;
	int idx_fileName = -1;
	int idx_strName = -1; //0?
	int idx_chainLength = -1;
	int idx_classId = -1;
	int idx_chainId = -1;
	int idx_nrOfWindows = -1;
	int idx_nrOfWindowPairs = -1;
	int idx_windowNr1 = -1, winNr1 = -1;
	int segIndices_1[2];
	int idx_i1 = -1;
	int idx_i2 = -2;
	int idx_windowNr2= -1, winNr2 = -1;
	int segIndices_2[2];
	int idx_j1 = -3;
	int idx_j2 = -4;
	int idx_I12 = -1;
	int idx_Ia12 = -1;
	int idx_I1234_full = -1;
	int idx_I1324_full = -1;
	int idx_I1423 = -1;
	int idx_Ia12a34_full = -1;
	int idx_Ia13a24_full = -1;
	int idx_Ia14a23 = -1;
	int idx_Ia1234_full = -1;
	int idx_Ia1324_full = -1;
	int idx_Ia1423 = -1;
	int idx_I12a34_full = -1;
	int idx_I13a24_full = -1;
	int idx_I14a23 = -1;

	int i_passedHeadingLines = 0;
	int i_passedValueLines = 0;
	int i_matchLineCnt = 0;

	int recNr = 0; //record number for each record in DB (ie for each loaded structure)
	int totalNrOfWinPairs = 0; //count of total nr of window pairs loaded
	int chunkRecNr; //cnt number of records within each chunk (when loading)
	int chunkNr; //cnt number of chunks
	long int currFilePos; //integer giving the current position in file (when reading)
	int nrOfWindows = 0;
	int maxNrOfWindows = 0;
	int nrOfWindowPairs = 0;
	int maxNrOfWindowPairs = 0;
	int rewSuccess = 1;

	/*alloc and init the output ptr*/
	out = (int *) malloc(2*sizeof(int));
	out[0] = 0;
	out[1] = 0;

	/*init values:*/
	winNr1 = -1;
	segIndices_1[0] = -1;
	segIndices_1[1] = -2;
	winNr2 = -1;
	segIndices_2[0] = -3;
	segIndices_2[1] = -4;


	/*Open the file containing the GIs on windows/pairs of the data base structures:*/
	ptr_fileDB = fopen(DBResultsFilePath, "r");

	if(ptr_fileDB == NULL){
		printf("Sorry -- cannot find the file holding the DB pairs results: %s\n", DBResultsFilePath);
		printf("Press enter if you want to continue.\n");
		getchar();

	}


	/*Read first entry of DB file to record its structure*/
	while (fgets(line,1000, ptr_fileDB)!=NULL){

		idxInLine = 0;
		if (strstr(line, "Line type") == line){

			/*loop trough the line and record important indices:*/
			linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
			while(linePart)
			{
				//printf("linePart: %s\n",linePart);

				if (strcmp(linePart, "PDBfile")==0){					
					idx_fileName = idxInLine;
				}

				if (strcmp(linePart, "structureId")==0){					
					idx_strName = idxInLine;
				}

				if (strcmp(linePart, "chainLength")==0){					
					idx_chainLength = idxInLine;
				}

				if (strcmp(linePart, "classId")==0){					
					idx_classId = idxInLine;
				}

				if (strcmp(linePart, "chainId")==0){					
					idx_chainId = idxInLine;
				}

				if (strcmp(linePart, "nrOfWindows")==0){					
					idx_nrOfWindows = idxInLine;
					/*printf("nrOfWindows idx: %d line: %s\n",idx_nrOfWindows, line);
					getchar();*/
				}

				if (strcmp(linePart, "nrOfWindowPairs")==0){					
					idx_nrOfWindowPairs = idxInLine;
					/*printf("nrOfWindowPairs idx: %d\n",idx_nrOfWindowPairs);
					getchar();*/
				}

				linePart=strtok(NULL, ";"); //this gets the next token of the string (which is line)
				idxInLine += 1;
			}
		}

		else {

			idxInLine = -1; //start at -1 since the first token in lines starting with the string "Values" is, of course, always "Values"; the values line following such a line however have no placeholder for this word (ie it ois displaced by 1 rel to the "Values"-starting line) 

			if (strstr(line, "Values") == line){

				/*loop trough the line and record important indices:*/
				linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
				while(linePart)
				{				
					if (strcmp(linePart, "windowNr1") == 0){					
						idx_windowNr1 = idxInLine;
					}
					if (strcmp(linePart,"i1") == 0){					
						idx_i1 = idxInLine;
					}
					if (strcmp(linePart,"i2") == 0){					
						idx_i2 = idxInLine;
					}
					if (strcmp(linePart, "windowNr2") == 0){					
						idx_windowNr2 = idxInLine;
					}
					if (strcmp(linePart,"j1") == 0){					
						idx_j1 = idxInLine;
					}
					if (strcmp(linePart,"j2") == 0){					
						idx_j2 = idxInLine;
					}
					if (strcmp(linePart,"I12") == 0){					
						idx_I12 = idxInLine;
						//printf("idx_I12  found:%d\n", idx_I12 );
					}
					if (strcmp(linePart,"Ia12") == 0){					
						idx_Ia12 = idxInLine;
					}
					
					if (strcmp(linePart,"I1234_full") == 0){					
						idx_I1234_full = idxInLine;
					}
					if (strcmp(linePart,"I1324_full") == 0){					
						idx_I1324_full = idxInLine;
					}
					if (strcmp(linePart,"I1423") == 0){					
						idx_I1423 = idxInLine;
					}

					if (strcmp(linePart,"Ia12a34_full") == 0){					
						idx_Ia12a34_full = idxInLine;
					}
					if (strcmp(linePart,"Ia13a24_full") == 0){					
						idx_Ia13a24_full = idxInLine;
					}
					if (strcmp(linePart,"Ia14a23") == 0){					
						idx_Ia14a23 = idxInLine;
					}

					if (strcmp(linePart,"Ia1234_full") == 0){					
						idx_Ia1234_full = idxInLine;
					}
					if (strcmp(linePart,"Ia1324_full") == 0){					
						idx_Ia1324_full = idxInLine;
					}
					if (strcmp(linePart,"Ia1423") == 0){					
						idx_Ia1423 = idxInLine;
					}

					if (strcmp(linePart,"I12a34_full") == 0){					
						idx_I12a34_full = idxInLine;
					}
					if (strcmp(linePart,"I13a24_full") == 0){					
						idx_I13a24_full = idxInLine;
					}
					if (strcmp(linePart,"I14a23") == 0){					
						idx_I14a23 = idxInLine;
					}

					linePart=strtok(NULL, ";"); //this gets the next token of the string (which is line)
					idxInLine += 1;
				}

				break;
			}

		}
	}



	/*get back to first line of file*/
	rewind(ptr_fileDB);

	/*loop through DB and keep the records in the pointer declared for the purpose;
	memory is allocated chunk-wise. This is handled as follows: in an outer while-loop
	the file's records are traversed; inside this two while-loops are placed, each 
	will traverse a "current chunk" once: in the first the highest nr of windows used
	for covering the structure within a chunk is found; this nr is used to re-allocate
	memory to the ptr holding the DB values; loading these into the ptr is done in the 
	second inner while-loop. Forhandling each chunk the first inner while-loop traverses
	it, upon exting this loop the ptr to the file is reset to the beginning of the chunk
	and the records are then loaded in in the next inner while-loop */

	recNr = 0; //init
	chunkRecNr = 0; //init
	chunkNr = 0; //init
	i_passedHeadingLines = 0; //reset
	while(i_passedHeadingLines == 0){ //outer while-loop
		
		currFilePos = ftell(ptr_fileDB); //this will be executed every time a new chunk is started at
		chunkRecNr = 0; //reset
		maxNrOfWindowPairs = 0; //reset

		while (fgets(line,1000, ptr_fileDB)!=NULL){ //1st inner while-loop

			/*exit the inner while-loop for allocating memory and recording the DB values  */
			if (strstr(line, "Line type") == line && i_passedHeadingLines != 0 && chunkRecNr == chunkSize){
				
				//printf("exit pt, chunk Nr: %d\n", chunkNr);
				//getchar();
				i_passedHeadingLines = 0;
				break;
			}
			else{
				if (strstr(line, "Header") == line){

					//printf("line: %s\n", line);

					chunkRecNr += 1;
					i_passedHeadingLines = 1;

					linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
					idxInLine = 0;
					//find the nr of windows for this chain/structure:
					while(linePart){

						if (idxInLine == idx_nrOfWindows){nrOfWindows = atoi(linePart);}

						if (nrOfWindows > maxNrOfWindows){maxNrOfWindows = nrOfWindows;}

						if (idxInLine == idx_nrOfWindowPairs){nrOfWindowPairs = atoi(linePart);}

						if (nrOfWindowPairs > maxNrOfWindowPairs){maxNrOfWindowPairs = nrOfWindowPairs;}

						linePart = strtok(NULL, ";"); //this gets the next token of the string (which is line)
						idxInLine += 1;
					}

					//printf("nr of winds: %d max so far: %d\n", nrOfWindowPairs, maxnrOfWindowPairs );
				}
				
			}

		} //1st inner while-loop 

		if(maxNrOfWindows > *ptr_maxNrOfWindows){

			*ptr_maxNrOfWindows = maxNrOfWindows; 

		}

		/*Now load the I values into ptr's declared for the purpose;
		we first allocate memory:*/
		returnVal = alloc_init_DBresultsWindowPairs(&ptr_dbResults, chunkNr, chunkSize, maxNrOfWindows, maxNrOfWindowPairs);

		/*for(m=0;m<chunkSize; m++){
			for(i=0; i< nrOfWindowPairs; i++){

				printf("at rec %d and win %d seg inds are: %d %d\n", m, i, ptr_dbResults[m].segIndices[i][0], ptr_dbResults[m].segIndices[i][1]); 

			}
		}*/

		//rewind to the position at the beginning of the current chunk:
		rewSuccess = fseek(ptr_fileDB, currFilePos, SEEK_SET);
		//carry on with the 2nd inner while-loop, and load the records into the pointer
		i_passedHeadingLines = 0;
		chunkRecNr = 0; //reset
		while (fgets(line,1000, ptr_fileDB)!=NULL){ //2nd inner while-loop
						 
			/*exit the inner while-loop for storing the record*/
			if (strstr(line, "Line type") == line && i_passedHeadingLines != 0 ){
				
				//increments and resets
				recNr += 1;
				chunkRecNr += 1;
				//i_matchLineCnt = 0;
				i_passedHeadingLines = 0;
				//printf("Rec nr: %d\n", recNr);
				//printf("Chunk rec nr: %d\n", chunkRecNr);
				//getchar();
				if(chunkRecNr == chunkSize){
					break;
				}


			}
			else{
				if (strstr(line, "Header") == line){

					linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
					idxInLine = 0;
					while(linePart){

						if (idxInLine == idx_fileName){
							//strcpy(linePartFileName,linePart); //have to strcpy via a "fixed" string (char array)
							strcpy(ptr_dbResults[recNr].fileName, linePart);
						}
						if (idxInLine == idx_strName){
							//strcpy(linePartFileName,linePart); //have to strcpy via a "fixed" string (char array)
							strcpy(ptr_dbResults[recNr].structureName, linePart);
						}
						if (idxInLine == idx_chainLength){
							ptr_dbResults[recNr].chainLen = atoi(linePart);
						}

						if (idxInLine == idx_classId){ 
							//strcpy(linePartChainId,linePart); //have to strcpy via a "fixed" string (char array)
							strcpy(ptr_dbResults[recNr].classId, linePart);
						}
						if (idxInLine == idx_chainId){ 
							//strcpy(linePartChainId,linePart); //have to strcpy via a "fixed" string (char array)
							strcpy(ptr_dbResults[recNr].chainId, linePart);
						}

						
						if (idxInLine == idx_nrOfWindows){
							ptr_dbResults[recNr].nrOfWindows = atoi(linePart);
						}

						if (idxInLine == idx_nrOfWindowPairs){
							ptr_dbResults[recNr].nrOfWindowPairs = atoi(linePart);
							//printf("found nrOfWindowPairs: %d\n",ptr_dbResults[recNr].nrOfWindowPairs);
						}

						linePart = strtok(NULL, ";"); //this gets the next token of the string (which is line)
						idxInLine += 1;
					}
					//printf("rec Nr 1 %d\n", recNr);
					/*printf(".....%s\n", ptr_dbResults[recNr].fileName);
					printf(".....%s\n", ptr_dbResults[recNr].chainId);
					printf(".....%d\n", ptr_dbResults[recNr].nrOfWindowPairs);*/
					//getchar();
				}
				else{

					if (strstr(line, "Values") == line){

						i_passedHeadingLines = 1;
						/*recNr += 1;
						chunkRecNr += 1;*/
						//printf("Rec nr: %d\n", recNr);
						//printf("Chunk rec nr: %d\n", chunkRecNr);

						i_matchLineCnt = 0;  //reset

					}
					else{

						if (i_passedHeadingLines == 1){

							//printf("rec Nr 2 %d\n", recNr);

							linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
							idxInLine = 0;
							while(linePart){

								if (idxInLine == idx_windowNr1){
									winNr1 = atoi(linePart);
									//ptr_dbResults[recNr].windowsNr_1[winNr1] = winNr1;
									//printf("w-nr .. %d\n", ptr_dbResults[recNr].windowsNr[i_matchLineCnt]); 
								}
								if (idxInLine == idx_i1){
									segIndices_1[0] = atoi(linePart);
									//ptr_dbResults[recNr].segIndices_1[winNr1][0] = atoi(linePart);
									//printf("i1 .. %d\n", ptr_dbResults[recNr].segIndices[i_matchLineCnt][0]); 
								}
								if (idxInLine == idx_i2){
									segIndices_1[1] = atoi(linePart);
									//printf("part ..:%s", linePart);
									//ptr_dbResults[recNr].segIndices_1[winNr1][1] = atoi(linePart);
	 								//printf("i2 .. %d\n", ptr_dbResults[recNr].segIndices[i_matchLineCnt][1]); 
								}

								if (idxInLine == idx_windowNr2){
									winNr2 = atoi(linePart);
									ptr_dbResults[recNr].windowPair[winNr1][winNr2].windowNr_1 = winNr1;
									ptr_dbResults[recNr].windowPair[winNr1][winNr2].segIndices_1[0] = segIndices_1[0];
									ptr_dbResults[recNr].windowPair[winNr1][winNr2].segIndices_1[1] = segIndices_1[1];
									ptr_dbResults[recNr].windowPair[winNr1][winNr2].windowNr_2 = winNr2;
									//ptr_dbResults[recNr].windowsNr_2[winNr2] = winNr2;
									//printf("w-nr .. %d\n", ptr_dbResults[recNr].windowsNr[i_matchLineCnt]); 
									totalNrOfWinPairs += 1;
								}
								if (idxInLine == idx_j1){
									//printf("winNr2: %d\n",winNr2);
									ptr_dbResults[recNr].windowPair[winNr1][winNr2].segIndices_2[0] = atoi(linePart);
									//ptr_dbResults[recNr].segIndices_2[winNr2][0] = atoi(linePart);
									//printf("i1 .. %d\n", ptr_dbResults[recNr].segIndices[i_matchLineCnt][0]); 
								}
								if (idxInLine == idx_j2){
									//printf("part ..:%s", linePart);
									ptr_dbResults[recNr].windowPair[winNr1][winNr2].segIndices_2[1] = atoi(linePart);
									//ptr_dbResults[recNr].segIndices_2[winNr2][1] = atoi(linePart);
									//printf("i2 .. %d\n", ptr_dbResults[recNr].segIndices[i_matchLineCnt][1]); 
								}



								if (idxInLine == idx_I12){ptr_dbResults[recNr].I12[winNr1][winNr2] = atof(linePart);}
								if (idxInLine == idx_Ia12){ptr_dbResults[recNr].Ia12[winNr1][winNr2] = atof(linePart);}
								
								if (idxInLine == idx_I1234_full){ptr_dbResults[recNr].I1234_full[winNr1][winNr2] = atof(linePart);}
								if (idxInLine == idx_I1324_full){ptr_dbResults[recNr].I1324_full[winNr1][winNr2] = atof(linePart);}
								if (idxInLine == idx_I1423){ptr_dbResults[recNr].I1423[winNr1][winNr2] = atof(linePart);}
								
								if (idxInLine == idx_Ia12a34_full){ptr_dbResults[recNr].Ia12a34_full[winNr1][winNr2] = atof(linePart);}
								if (idxInLine == idx_Ia13a24_full){ptr_dbResults[recNr].Ia13a24_full[winNr1][winNr2] = atof(linePart);}
								if (idxInLine == idx_Ia14a23){ptr_dbResults[recNr].Ia14a23[winNr1][winNr2] = atof(linePart);}
								
								if (idxInLine == idx_Ia1234_full){ptr_dbResults[recNr].Ia1234_full[winNr1][winNr2] = atof(linePart);}
								if (idxInLine == idx_Ia1324_full){ptr_dbResults[recNr].Ia1324_full[winNr1][winNr2] = atof(linePart);}
								if (idxInLine == idx_Ia1423){ptr_dbResults[recNr].Ia1423[winNr1][winNr2] = atof(linePart);}
								
								if (idxInLine == idx_I12a34_full){ptr_dbResults[recNr].I12a34_full[winNr1][winNr2] = atof(linePart);}
								if (idxInLine == idx_I13a24_full){ptr_dbResults[recNr].I13a24_full[winNr1][winNr2] = atof(linePart);}
								if (idxInLine == idx_I14a23){ptr_dbResults[recNr].I14a23[winNr1][winNr2] = atof(linePart);}
								
								linePart = strtok(NULL, ";"); //this gets the next token of the string (which is line)
								idxInLine += 1;
							}

							i_matchLineCnt +=1;

						}
					}
				}
			}

		} //2nd inner while-loop

		chunkNr += 1;
		//printf("chunkNr: %d\n", chunkNr);
	}

	/*for(i=0; i<13;i++){
		printf("ptr_dbResults filename at %d: %s\n", i, ptr_dbResults[i].fileName);
		printf("............. nrOfWindowPairs at %d: %d\n", i, ptr_dbResults[i].nrOfWindowPairs);
	}*/

	//pass the ptr_dbResults to in/out variable:
	*ptr_ptr_dbResults = ptr_dbResults;

	printf("Nr of records loaded: %d\n", recNr+1); 

	fclose(ptr_fileDB);

	out[0] = recNr + 1; /*the recNr is init'ed at 0, and the first rec has that index*/ 
	out[1] = totalNrOfWinPairs;
	//printf("Total nr of window pairs loaded: %d\n", out[1]); 

	return out;

}

struct I_windowPairs_order2_meanStddev normalizeDBresultsWindowPairsPtr(struct db_I_windowPairs_order2_ptr *ptr_dbResults, int nrOfStructures, int nrOfRecords){

	int n = 0;
	int i = 0, j = 0;

	int totalNrOfWindowPairs = nrOfRecords;  /* nrOfStructures*1000; /*if a structure is covered by 50 windows, there will be about 50*50/2 = 50*25 = 1250 window pairs */
	
	int cnt  = 0;

	struct I_windows_order2_raw_ptr I_windowPairs_order2_raw_ptr; /*We "abuse" the struct aimed for single windows!*/ 
	struct I_windowPairs_order2_meanStddev I_windowPairs_order2_meanStddev_out;

	struct double2 meanStddevOut;

	/*allocate memory*/
	printf("nrOfRecords: %d \n", nrOfRecords);
	totalNrOfWindowPairs = floor(nrOfRecords*1.05); /*just to be on the safe side*/
	printf("totalNrOfWindowPairs: %d \n", totalNrOfWindowPairs );
	//getchar();

	I_windowPairs_order2_raw_ptr.I12 = (double *) malloc(totalNrOfWindowPairs*sizeof(double));
	I_windowPairs_order2_raw_ptr.Ia12 = (double *) malloc(totalNrOfWindowPairs*sizeof(double));
	
	I_windowPairs_order2_raw_ptr.I1234_full = (double *) malloc(totalNrOfWindowPairs*sizeof(double));
	I_windowPairs_order2_raw_ptr.I1324_full = (double *) malloc(totalNrOfWindowPairs*sizeof(double));
	I_windowPairs_order2_raw_ptr.I1423 = (double *) malloc(totalNrOfWindowPairs*sizeof(double));

	I_windowPairs_order2_raw_ptr.Ia12a34_full = (double *) malloc(totalNrOfWindowPairs*sizeof(double));
	I_windowPairs_order2_raw_ptr.Ia13a24_full = (double *) malloc(totalNrOfWindowPairs*sizeof(double));
	I_windowPairs_order2_raw_ptr.Ia14a23 = (double *) malloc(totalNrOfWindowPairs*sizeof(double));

	I_windowPairs_order2_raw_ptr.Ia1234_full = (double *) malloc(totalNrOfWindowPairs*sizeof(double));
	I_windowPairs_order2_raw_ptr.Ia1324_full = (double *) malloc(totalNrOfWindowPairs*sizeof(double));
	I_windowPairs_order2_raw_ptr.Ia1423 = (double *) malloc(totalNrOfWindowPairs*sizeof(double));

	I_windowPairs_order2_raw_ptr.I12a34_full = (double *) malloc(totalNrOfWindowPairs*sizeof(double));
	I_windowPairs_order2_raw_ptr.I13a24_full = (double *) malloc(totalNrOfWindowPairs*sizeof(double));
	I_windowPairs_order2_raw_ptr.I14a23 = (double *) malloc(totalNrOfWindowPairs*sizeof(double));

	/*init: read in all values to ptr*/
	for(n=0;n<nrOfStructures;n++){

		//printf(".. (saaledes stadig ganske hemmelig) ... ok at %d\n", n);

		for(i= 0; i< ptr_dbResults[n].nrOfWindows; i++){
			for(j= i; j< ptr_dbResults[n].nrOfWindows; j++){

				I_windowPairs_order2_raw_ptr.I12[cnt] = ptr_dbResults[n].I12[i][j];
				I_windowPairs_order2_raw_ptr.Ia12[cnt] = ptr_dbResults[n].Ia12[i][j];
				
				I_windowPairs_order2_raw_ptr.I1234_full[cnt] = ptr_dbResults[n].I1234_full[i][j];
				I_windowPairs_order2_raw_ptr.I1324_full[cnt] = ptr_dbResults[n].I1324_full[i][j];
				I_windowPairs_order2_raw_ptr.I1423[cnt] = ptr_dbResults[n].I1423[i][j];

				I_windowPairs_order2_raw_ptr.Ia12a34_full[cnt] = ptr_dbResults[n].Ia12a34_full[i][j];
				I_windowPairs_order2_raw_ptr.Ia13a24_full[cnt] = ptr_dbResults[n].Ia13a24_full[i][j];
				I_windowPairs_order2_raw_ptr.Ia14a23[cnt] = ptr_dbResults[n].Ia14a23[i][j];

				I_windowPairs_order2_raw_ptr.Ia1234_full[cnt] = ptr_dbResults[n].Ia1234_full[i][j];
				I_windowPairs_order2_raw_ptr.Ia1324_full[cnt] = ptr_dbResults[n].Ia1324_full[i][j];
				I_windowPairs_order2_raw_ptr.Ia1423[cnt] = ptr_dbResults[n].Ia1423[i][j];

				I_windowPairs_order2_raw_ptr.I12a34_full[cnt] = ptr_dbResults[n].I12a34_full[i][j];
				I_windowPairs_order2_raw_ptr.I13a24_full[cnt] = ptr_dbResults[n].I13a24_full[i][j];
				I_windowPairs_order2_raw_ptr.I14a23[cnt] = ptr_dbResults[n].I14a23[i][j];

				cnt +=1;

				/*reallocate memory if necessary*/
				if(cnt > totalNrOfWindowPairs){

					totalNrOfWindowPairs = ((int) floor(cnt*1.25));
					printf("totalNrOfWindowPairs now: %d \n", totalNrOfWindowPairs);

					I_windowPairs_order2_raw_ptr.I12 = (double *) realloc(I_windowPairs_order2_raw_ptr.I12,  totalNrOfWindowPairs*sizeof(double));
					I_windowPairs_order2_raw_ptr.Ia12 = (double *) realloc(I_windowPairs_order2_raw_ptr.Ia12,  totalNrOfWindowPairs*sizeof(double));
					
					I_windowPairs_order2_raw_ptr.I1234_full = (double *) realloc(I_windowPairs_order2_raw_ptr.I1234_full,  totalNrOfWindowPairs*sizeof(double));
					I_windowPairs_order2_raw_ptr.I1324_full = (double *) realloc(I_windowPairs_order2_raw_ptr.I1324_full,  totalNrOfWindowPairs*sizeof(double));
					I_windowPairs_order2_raw_ptr.I1423 = (double *) realloc(I_windowPairs_order2_raw_ptr.I1423,  totalNrOfWindowPairs*sizeof(double));

					I_windowPairs_order2_raw_ptr.Ia12a34_full = (double *) realloc(I_windowPairs_order2_raw_ptr.Ia12a34_full,  totalNrOfWindowPairs*sizeof(double));
					I_windowPairs_order2_raw_ptr.Ia13a24_full = (double *) realloc(I_windowPairs_order2_raw_ptr.Ia13a24_full,  totalNrOfWindowPairs*sizeof(double));
					I_windowPairs_order2_raw_ptr.Ia14a23 = (double *) realloc(I_windowPairs_order2_raw_ptr.Ia14a23,  totalNrOfWindowPairs*sizeof(double));

					I_windowPairs_order2_raw_ptr.Ia1234_full = (double *) realloc(I_windowPairs_order2_raw_ptr.Ia1234_full,  totalNrOfWindowPairs*sizeof(double));
					I_windowPairs_order2_raw_ptr.Ia1324_full = (double *) realloc(I_windowPairs_order2_raw_ptr.Ia1324_full,  totalNrOfWindowPairs*sizeof(double));
					I_windowPairs_order2_raw_ptr.Ia1423 = (double *) realloc(I_windowPairs_order2_raw_ptr.Ia1423,  totalNrOfWindowPairs*sizeof(double));

					I_windowPairs_order2_raw_ptr.I12a34_full = (double *) realloc(I_windowPairs_order2_raw_ptr.I12a34_full,  totalNrOfWindowPairs*sizeof(double));
					I_windowPairs_order2_raw_ptr.I13a24_full = (double *) realloc(I_windowPairs_order2_raw_ptr.I13a24_full,  totalNrOfWindowPairs*sizeof(double));
					I_windowPairs_order2_raw_ptr.I14a23 = (double *) realloc(I_windowPairs_order2_raw_ptr.I14a23,  totalNrOfWindowPairs*sizeof(double));

				}

			}
		}

	}

	/*compute mean and std dev's:*/
	meanStddevOut = meanStddev(I_windowPairs_order2_raw_ptr.I12, cnt);
	I_windowPairs_order2_meanStddev_out.mean_I12 = meanStddevOut.val[0];
	I_windowPairs_order2_meanStddev_out.stddev_I12 = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windowPairs_order2_raw_ptr.Ia12, cnt);
	I_windowPairs_order2_meanStddev_out.mean_Ia12 = meanStddevOut.val[0];
	I_windowPairs_order2_meanStddev_out.stddev_Ia12 = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windowPairs_order2_raw_ptr.I1234_full, cnt);
	I_windowPairs_order2_meanStddev_out.mean_I1234_full = meanStddevOut.val[0];
	I_windowPairs_order2_meanStddev_out.stddev_I1234_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windowPairs_order2_raw_ptr.I1324_full, cnt);
	I_windowPairs_order2_meanStddev_out.mean_I1324_full = meanStddevOut.val[0];
	I_windowPairs_order2_meanStddev_out.stddev_I1324_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windowPairs_order2_raw_ptr.I1423, cnt);
	I_windowPairs_order2_meanStddev_out.mean_I1423 = meanStddevOut.val[0];
	I_windowPairs_order2_meanStddev_out.stddev_I1423 = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windowPairs_order2_raw_ptr.Ia12a34_full, cnt);
	I_windowPairs_order2_meanStddev_out.mean_Ia12a34_full = meanStddevOut.val[0];
	I_windowPairs_order2_meanStddev_out.stddev_Ia12a34_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windowPairs_order2_raw_ptr.Ia13a24_full, cnt);
	I_windowPairs_order2_meanStddev_out.mean_Ia13a24_full = meanStddevOut.val[0];
	I_windowPairs_order2_meanStddev_out.stddev_Ia13a24_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windowPairs_order2_raw_ptr.Ia14a23, cnt);
	I_windowPairs_order2_meanStddev_out.mean_Ia14a23 = meanStddevOut.val[0];
	I_windowPairs_order2_meanStddev_out.stddev_Ia14a23 = meanStddevOut.val[1];

		meanStddevOut = meanStddev(I_windowPairs_order2_raw_ptr.Ia1234_full, cnt);
	I_windowPairs_order2_meanStddev_out.mean_Ia1234_full = meanStddevOut.val[0];
	I_windowPairs_order2_meanStddev_out.stddev_Ia1234_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windowPairs_order2_raw_ptr.Ia1324_full, cnt);
	I_windowPairs_order2_meanStddev_out.mean_Ia1324_full = meanStddevOut.val[0];
	I_windowPairs_order2_meanStddev_out.stddev_Ia1324_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windowPairs_order2_raw_ptr.Ia1423, cnt);
	I_windowPairs_order2_meanStddev_out.mean_Ia1423 = meanStddevOut.val[0];
	I_windowPairs_order2_meanStddev_out.stddev_Ia1423 = meanStddevOut.val[1];

		meanStddevOut = meanStddev(I_windowPairs_order2_raw_ptr.I12a34_full, cnt);
	I_windowPairs_order2_meanStddev_out.mean_I12a34_full = meanStddevOut.val[0];
	I_windowPairs_order2_meanStddev_out.stddev_I12a34_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windowPairs_order2_raw_ptr.I13a24_full, cnt);
	I_windowPairs_order2_meanStddev_out.mean_I13a24_full = meanStddevOut.val[0];
	I_windowPairs_order2_meanStddev_out.stddev_I13a24_full = meanStddevOut.val[1];

	meanStddevOut = meanStddev(I_windowPairs_order2_raw_ptr.I14a23, cnt);
	I_windowPairs_order2_meanStddev_out.mean_I14a23 = meanStddevOut.val[0];
	I_windowPairs_order2_meanStddev_out.stddev_I14a23 = meanStddevOut.val[1];


	/*Now scale the input by the std-dev*/
	for(n=0;n<nrOfStructures;n++){

		normalizeSimplexArray(I_windowPairs_order2_meanStddev_out.mean_I12, I_windowPairs_order2_meanStddev_out.stddev_I12, &ptr_dbResults[n].I12, ptr_dbResults[n].nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_out.mean_Ia12, I_windowPairs_order2_meanStddev_out.stddev_Ia12, &ptr_dbResults[n].Ia12, ptr_dbResults[n].nrOfWindows);
		
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_out.mean_I1234_full, I_windowPairs_order2_meanStddev_out.stddev_I1234_full, &ptr_dbResults[n].I1234_full, ptr_dbResults[n].nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_out.mean_I1324_full, I_windowPairs_order2_meanStddev_out.stddev_I1324_full, &ptr_dbResults[n].I1324_full, ptr_dbResults[n].nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_out.mean_I1423, I_windowPairs_order2_meanStddev_out.stddev_I1423, &ptr_dbResults[n].I1423, ptr_dbResults[n].nrOfWindows);

		normalizeSimplexArray(I_windowPairs_order2_meanStddev_out.mean_Ia12a34_full, I_windowPairs_order2_meanStddev_out.stddev_Ia12a34_full, &ptr_dbResults[n].Ia12a34_full, ptr_dbResults[n].nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_out.mean_Ia13a24_full, I_windowPairs_order2_meanStddev_out.stddev_Ia13a24_full, &ptr_dbResults[n].Ia13a24_full, ptr_dbResults[n].nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_out.mean_Ia14a23, I_windowPairs_order2_meanStddev_out.stddev_Ia14a23, &ptr_dbResults[n].Ia14a23, ptr_dbResults[n].nrOfWindows);

		normalizeSimplexArray(I_windowPairs_order2_meanStddev_out.mean_Ia1234_full, I_windowPairs_order2_meanStddev_out.stddev_Ia1234_full, &ptr_dbResults[n].Ia1234_full, ptr_dbResults[n].nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_out.mean_Ia1324_full, I_windowPairs_order2_meanStddev_out.stddev_Ia1324_full, &ptr_dbResults[n].Ia1324_full, ptr_dbResults[n].nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_out.mean_Ia1423, I_windowPairs_order2_meanStddev_out.stddev_Ia1423, &ptr_dbResults[n].Ia1423, ptr_dbResults[n].nrOfWindows);

		normalizeSimplexArray(I_windowPairs_order2_meanStddev_out.mean_I12a34_full, I_windowPairs_order2_meanStddev_out.stddev_I12a34_full, &ptr_dbResults[n].I12a34_full, ptr_dbResults[n].nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_out.mean_I13a24_full, I_windowPairs_order2_meanStddev_out.stddev_I13a24_full, &ptr_dbResults[n].I13a24_full, ptr_dbResults[n].nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_out.mean_I14a23, I_windowPairs_order2_meanStddev_out.stddev_I14a23, &ptr_dbResults[n].I14a23, ptr_dbResults[n].nrOfWindows);

	}

	/*free memory*/
	free(I_windowPairs_order2_raw_ptr.I12);
	free(I_windowPairs_order2_raw_ptr.Ia12);

	free(I_windowPairs_order2_raw_ptr.I1234_full);
	free(I_windowPairs_order2_raw_ptr.I1324_full);
	free(I_windowPairs_order2_raw_ptr.I1423);

	free(I_windowPairs_order2_raw_ptr.Ia12a34_full);
	free(I_windowPairs_order2_raw_ptr.Ia13a24_full);
	free(I_windowPairs_order2_raw_ptr.Ia14a23);

	free(I_windowPairs_order2_raw_ptr.Ia1234_full);
	free(I_windowPairs_order2_raw_ptr.Ia1324_full);
	free(I_windowPairs_order2_raw_ptr.Ia1423);

	free(I_windowPairs_order2_raw_ptr.I12a34_full);
	free(I_windowPairs_order2_raw_ptr.I13a24_full);
	free(I_windowPairs_order2_raw_ptr.I14a23);

	return I_windowPairs_order2_meanStddev_out;

}

void normalizeQueryWindowPairsPtr(struct I_windowPairs_ptr I_windowPairs, struct I_windowPairs_order2_meanStddev I_windowPairs_order2_meanStddev_in, int nrOfInvs ){

	/*normalize using the input means and std devs*/
	normalizeSimplexArray(I_windowPairs_order2_meanStddev_in.mean_I12, I_windowPairs_order2_meanStddev_in.stddev_I12, &I_windowPairs.I12, I_windowPairs.nrOfWindows);
	normalizeSimplexArray(I_windowPairs_order2_meanStddev_in.mean_Ia12, I_windowPairs_order2_meanStddev_in.stddev_Ia12, &I_windowPairs.Ia12, I_windowPairs.nrOfWindows);
	if(nrOfInvs > 2){
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_in.mean_I1234_full, I_windowPairs_order2_meanStddev_in.stddev_I1234_full, &I_windowPairs.I1234_full, I_windowPairs.nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_in.mean_I1324_full, I_windowPairs_order2_meanStddev_in.stddev_I1324_full, &I_windowPairs.I1324_full, I_windowPairs.nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_in.mean_I1423, I_windowPairs_order2_meanStddev_in.stddev_I1423, &I_windowPairs.I1423, I_windowPairs.nrOfWindows);

		normalizeSimplexArray(I_windowPairs_order2_meanStddev_in.mean_Ia12a34_full, I_windowPairs_order2_meanStddev_in.stddev_Ia12a34_full, &I_windowPairs.Ia12a34_full, I_windowPairs.nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_in.mean_Ia13a24_full, I_windowPairs_order2_meanStddev_in.stddev_Ia13a24_full, &I_windowPairs.Ia13a24_full, I_windowPairs.nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_in.mean_Ia14a23, I_windowPairs_order2_meanStddev_in.stddev_Ia14a23, &I_windowPairs.Ia14a23, I_windowPairs.nrOfWindows);

		normalizeSimplexArray(I_windowPairs_order2_meanStddev_in.mean_Ia1234_full, I_windowPairs_order2_meanStddev_in.stddev_Ia1234_full, &I_windowPairs.Ia1234_full, I_windowPairs.nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_in.mean_Ia1324_full, I_windowPairs_order2_meanStddev_in.stddev_Ia1324_full, &I_windowPairs.Ia1324_full, I_windowPairs.nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_in.mean_Ia1423, I_windowPairs_order2_meanStddev_in.stddev_Ia1423, &I_windowPairs.Ia1423, I_windowPairs.nrOfWindows);

		normalizeSimplexArray(I_windowPairs_order2_meanStddev_in.mean_I12a34_full, I_windowPairs_order2_meanStddev_in.stddev_I12a34_full, &I_windowPairs.I12a34_full, I_windowPairs.nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_in.mean_I13a24_full, I_windowPairs_order2_meanStddev_in.stddev_I13a24_full, &I_windowPairs.I13a24_full, I_windowPairs.nrOfWindows);
		normalizeSimplexArray(I_windowPairs_order2_meanStddev_in.mean_I14a23, I_windowPairs_order2_meanStddev_in.stddev_I14a23, &I_windowPairs.I14a23, I_windowPairs.nrOfWindows);

	}
}


/*Read in dist-data (for allowing a score based on distribution of dist's obtained in some search) and sort them 
(lowest first). 
Updates the input ptr_ptr_distances (pointing to a sorted array of dist-values) and returns the length of the array. 
A raw, tabular structure is expected of the file: 
1) it must be ";" separated; 
2) the first line must provide the field names and in proper order, and one field must be named "dist"
3) all following lines containing data*/
int readSortDistData(FILE *ptr_fileDistData, double **ptr_ptr_distances){

	int i;

	double *ptr_distances_temp;

	int chunkSize = 10000;
	int cnt = 0;

	char line[1000];
	char *linePart;

	int idxInLine;
	int idx_dist = -1;

	/*Allocate memory initial memory:*/
	ptr_distances_temp = (double *) malloc(chunkSize*sizeof(double));

	///*init*/
	//for(i=0;i=chunkSize;i++){

	//	ptr_distances_temp[i] = 1e10; //just some pos number alrger than any expected distance

	//}

	/*First find the index of the "dist" field inthe first line of the file:*/
	while (fgets(line,1000, ptr_fileDistData)!=NULL && cnt == 0 ){

		idxInLine = 0;
		
		/*loop trough the line and record important indices:*/
		linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
		while(linePart)
		{
			if (strcmp(linePart, "dist")==0){					
				idx_dist = idxInLine;
			}
			linePart=strtok(NULL, ";"); //this gets the next token of the string (which is line)
			idxInLine += 1;
		}

		cnt +=1;
	}

	/*printf("idx_dist: %d\n", idx_dist);
	getchar();*/

	cnt = 0; //reset

	/*Read in:*/
	while (fgets(line,1000, ptr_fileDistData)!=NULL){
		
		/*Allocate or re-allocate memory as needed:*/
		if(cnt%(chunkSize - 1) == 0){
			
			ptr_distances_temp = (double *) realloc(ptr_distances_temp, (cnt+chunkSize)*sizeof(double));

			/*init*/
			/*for(i=cnt;i=cnt+chunkSize;i++){

				ptr_distances_temp[i] = 1e10;

			}*/
		}


		/*loop trough the line and record important indices:*/
		linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
		idxInLine = 0;
		while(linePart){
			if (idxInLine == idx_dist){
				ptr_distances_temp[cnt] = atof(linePart);
				cnt += 1;
			}
			linePart=strtok(NULL, ";"); //this gets the next token of the string (which is line)
			idxInLine += 1;
		}

	}

	//printf("cnt: %d\n", cnt);
	//for(i=0;i<50;i++){printf("dist at %d: %lf\n", i, ptr_distances_temp[i]);}
	//getchar();

	/*Sort the array:*/
	heapSort1dArray(ptr_distances_temp, cnt, 0);
	//for(i=0;i<50;i++){printf("dist at %d after sorting: %lf\n", i, ptr_distances_temp[i]);}
	//getchar();
	/*Record and return:*/
	*ptr_ptr_distances = ptr_distances_temp;

	return cnt;

}


int readSortRarityScoreData(FILE *ptr_fileRarityScoreData, double **ptr_ptr_rarityScores){

	int i;

	double *ptr_rarityScores_temp;
	//double *ptr_rarityScores_temp2;

	int chunkSize = 10000;
	int cnt = 0;

	char line[20000];
	char *linePart;

	int idxInLine;
	int idx_rarityScore = -1;

	int idx_score = -1;

	/*Allocate memory initial memory:*/
	ptr_rarityScores_temp = (double *) malloc(chunkSize*sizeof(double));

	///*init*/
	//for(i=0;i=chunkSize;i++){

	//	ptr_distances_temp[i] = 1e10; //just some pos number larger than any expected distance

	//}

	/*First find the index of the "dist" field in the first line of the file:*/
	while (fgets(line,10000, ptr_fileRarityScoreData)!=NULL && cnt == 0){

		if(cnt == 0){

			idxInLine = 0;
			
			/*loop trough the line and record important indices:*/
			linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
			while(linePart)
			{
				if (strcmp(linePart, "rarityScore")==0){					
					idx_rarityScore = idxInLine;
					//printf("The rarityScore sits at: %d\n", idx_rarityScore);
				}
				if (strcmp(linePart, "score")==0){					
					idx_score = idxInLine;
				}
				linePart=strtok(NULL, ";"); //this gets the next token of the string (which is line)
				idxInLine += 1;
			}
			
		}

		cnt +=1;
	}

	//printf("line", line);
	//printf("idx_dist: %d\n", idx_dist);
	//getchar();

	cnt = 0; //reset
	/*get back to first line of file*/
	rewind(ptr_fileRarityScoreData);
	

	//Read in:
	while (fgets(line,10000, ptr_fileRarityScoreData)!=NULL){
		
		
		/*Allocate or re-allocate memory as needed:*/
		if(cnt%(chunkSize - 1) == 0 && cnt > 0){
			
			ptr_rarityScores_temp = (double *) realloc(ptr_rarityScores_temp, (cnt+chunkSize)*sizeof(double));

		}

		//printf("line %s\n", line);
		//getchar();

		/*loop trough the line and record important indices:*/
		linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
		idxInLine = 0;
		while(linePart){
			//printf("Line part %s\n", linePart);
			if (idxInLine == idx_rarityScore){
				if(cnt >0){ //the first line is a header and will erroneously contribute a 0
					ptr_rarityScores_temp[cnt] = atof(linePart);
					//printf("(saaledes stadig ganske hemmelig), at cnt %d it read in: %lf\n", cnt, ptr_rarityScores_temp[cnt]);
					//getchar();
				}
				cnt += 1;
			}
			else if (idxInLine == idx_score){
				if(cnt >0){ //the first line is a header and will erroneously contribute a 0
					ptr_rarityScores_temp[cnt] = atof(linePart);
				}
				cnt += 1;
			}
			linePart=strtok(NULL, ";"); //this gets the next token of the string (which is line)
			idxInLine += 1;
		}

	}

	/*
	//We use a 2nd temp ptr to avoid getting a longer array than needed to hold the scores:
	ptr_rarityScores_temp2 = (double *) malloc(cnt*sizeof(double));
	//Read in the values:
	for(i=0;i<cnt;i++){

		ptr_rarityScores_temp2[i] = ptr_rarityScores_temp[i];

	}
	*/

	//let go of the mem to the temp ptr:
	//free(ptr_rarityScores_temp);

	/*printf("cnt: %d\n", cnt);
	for(i=0;i<30;i++){printf("score at %d: %lf\n", i, ptr_rarityScores_temp[i]);}
	getchar();
	*/

	/*Sort the array:*/
	heapSort1dArray(ptr_rarityScores_temp, cnt, 0);
	
	/*for(i=0;i<30;i++){
		printf("score at %d after sorting: %lf\n", i, ptr_rarityScores_temp[i]);
		//if(i%100 == 0){getchar();}
	}
	getchar();
	*/

	/*Record and return:*/
	*ptr_ptr_rarityScores = ptr_rarityScores_temp;

	return cnt; 

}

/*Fcts for binning*/

/*Generate a binning of nrOfBins bins using one of two methods
binningType = 1: for each entry (invariant) divide in equi-distant bins: divide the range from the minimum value of the 
invariant to the maximum in nrOfBins equally long bins
binningType = 2: for each entry (invariant) divide the data in "equal percentiles": divide the ordered data points in 
nrOfBins equally large portions*/
int generateBins(int binningType, struct db_I_windows_order2_ptr *ptr_dbResults, int nrOfDBstructures, int nrOfDBrecords, int nrOfEntries, double ***ptr_bins, int nrOfBins, double compensationForRounding){

	int returnVal = 0;

	double **bins;
	int i, j, k, m, n;
	int allocateForNrOfDBrecords;


	double 	*minIVal, *maxIVal;

	double *Ivector; /*convenience ptr*/

	double **Ivector_2; /*another convenience ptr; to be used in the type 2 method*/
	int cnt = 0;
	int idx = 0;

	double binSize;

	int binCount;


	/*Allocate memory to the bins pointer:*/
	bins = (double **) malloc(nrOfEntries*nrOfBins*sizeof(double));
	for(j = 0; j< nrOfEntries; j++){

		bins[j] = (double *) malloc(nrOfBins*sizeof(double));

		/*init*/
		for(k=0; k < nrOfBins; k++){

			bins[j][k] = -117.0;
		
		}


	}

	if(binningType == 0){ //Custom type; change bin definitions below when nrOf bins change!

	
		for(j = 0; j< nrOfEntries; j++){

			//bins = {-2.0,-1.0, -0.5, 0.0, 0.5, 1.0, 2.0};
			/*bins[0] = -2;
			bins[1] = -1;
			bins[2] = -0.5;
			bins[3] = 0;
			bins[4] = 0.5;
			bins[5] = 1;
			bins[6] = 2;*/


			/*bins[0] = -2;
			bins[1] = -1.5;
			bins[2] = -1;
			bins[3] = -0.5;
			bins[4] = -0.25;
			bins[5] = 0.0;
			bins[6] = 0.25;
			bins[7] = 0.5;
			bins[8] = 1.0;
			bins[9] = 1.5;
			bins[10] = 2;*/

			/*bins[0] = -2.0;
			bins[1] = -1.5;
			bins[2] = -1;
			bins[3] = -0.75;
			bins[4] = -0.5;
			bins[5] = -0.25;
			bins[6] = -0.1;
			bins[7] = 0.0;
			bins[8] = 0.1;
			bins[9] = 0.25;
			bins[10] = 0.5;
			bins[11] = 0.75;
			bins[12] = 1.0;
			bins[13] = 1.5;
			bins[14] = 2.0;*/

			/*bins[j][0] = -2;
			bins[j][1] = -1.6;
			bins[j][2] = -1.2;
			bins[j][3] = -0.8;
			bins[j][4] = -0.4;
			bins[j][5] = 0.0;
			bins[j][6] = 0.4;
			bins[j][7] = 0.8;
			bins[j][8] = 1.2;
			bins[j][9] = 1.6;
			bins[j][10] = 2;*/

			bins[j][0] = -3.0;
			bins[j][1] = -2.5;
			bins[j][2] = -2;
			bins[j][3] = -1.5;
			bins[j][4] = -1;
			bins[j][5] = -0.5;
			bins[j][6] = -0.25;
			bins[j][7] = -0.1;
			bins[j][8] =  0.0;
			bins[j][9] = 0.1;
			bins[j][10] = 0.25;
			bins[j][11] = 0.5;
			bins[j][12] = 1;
			bins[j][13] = 1.5;
			bins[j][14] = 2;
			bins[j][15] = 2.5;
			bins[j][16] = 3.0;

			//bins[0] = -3.0;
			//bins[1] = -2.75;
			//bins[2] = -2.5;
			//bins[3] = -2.25;
			//bins[4] = -2.0;
			//bins[5] = -1.75;
			//bins[6] = -1.5;
			//bins[7] = -1.25;
			//bins[8] = -1.0;
			//bins[9] = -0.75;
			//bins[10] = -0.5;
			//bins[11] = -0.25;
			//bins[12] = 0;
			//bins[13] = 0.25;
			//bins[14] = 0.5;
			//bins[15] = 0.75;
			//bins[16] = 1.0;
			//bins[17] = 1.25;
			//bins[18] = 1.5;
			//bins[19] = 1.75;
			//bins[20] = 2.0;
			//bins[21] = 2.25;
			//bins[22] = 2.5;
			//bins[23] = 2.75;
			//bins[24] = 3.0;

			/*bins[0] = -3.0;
			bins[1] = -2.5;
			bins[2] = -2.0;
			bins[3] = -1.75;
			bins[4] = -1.5;
			bins[5] = -1.25;
			bins[6] = -1.0;
			bins[7] = -0.75;
			bins[8] = -0.5;
			bins[9] = -0.25;
			bins[10] = -0.125;
			bins[11] = 0;
			bins[12] = 0.125;
			bins[13] = 0.25;
			bins[14] = 0.5;
			bins[15] = 0.75;
			bins[16] = 1.0;
			bins[17] = 1.25;
			bins[18] = 1.5;
			bins[19] = 1.75;
			bins[20] = 2.0;
			bins[21] = 2.5;
			bins[22] = 3.0;*/

		}

	}

	/*Equi-distant binning*/
	if(binningType == 1){

		/*alloc and init*/
		minIVal = (double *) malloc(nrOfEntries*sizeof(double));
		maxIVal = (double *) malloc(nrOfEntries*sizeof(double));	

		for(j = 0; j< nrOfEntries; j++){

			minIVal[j] = 1e20;
			maxIVal[j] = -1e20;

		}


		Ivector = (double *) malloc(maxNrOfInvariants*sizeof(double));

		/*init*/
		for(j=0;j<maxNrOfInvariants;j++){ 
			Ivector[j] = 0.0; 
		}


		/*loop through the data, and record for each entry (invariant) its min and max value. Then  
		derive the binning:*/
		for(n = 0; n < nrOfDBstructures; n++){

			for(i = 0; i < ptr_dbResults[n].nrOfWindows; i++){


				Ivector[0] =  ptr_dbResults[n].I12[i];
				Ivector[1] =  ptr_dbResults[n].Ia12[i];

				Ivector[2] =  ptr_dbResults[n].I1234_full[i];
				Ivector[3] =  ptr_dbResults[n].I1324_full[i];
				Ivector[4] =  ptr_dbResults[n].I1423[i];

				Ivector[5] =  ptr_dbResults[n].Ia12a34_full[i];
				Ivector[6] =  ptr_dbResults[n].Ia13a24_full[i];
				Ivector[7] =  ptr_dbResults[n].Ia14a23[i];

				Ivector[8] =  ptr_dbResults[n].Ia1234_full[i];
				Ivector[9] =  ptr_dbResults[n].Ia1324_full[i];
				Ivector[10] =  ptr_dbResults[n].Ia1423[i];

				Ivector[11] =  ptr_dbResults[n].I12a34_full[i];
				Ivector[12] =  ptr_dbResults[n].I13a24_full[i];
				Ivector[13] =  ptr_dbResults[n].I14a23[i];

				for(j=0; j < nrOfEntries;j++){

					if(Ivector[j] < minIVal[j]){
						minIVal[j] = Ivector[j];
					}
					if(Ivector[j] > maxIVal[j]){
						maxIVal[j] = Ivector[j];
					}

				}
			}
		}

		/*Derive the bins*/
		for(j=0; j < nrOfEntries;j++){

			binSize = (double) (maxIVal[j] - minIVal[j])/(nrOfBins+1);

			//printf("bin %d size: %lf\n", j, binSize);
			//getchar();

			for(k = 0; k < nrOfBins;k++){

				bins[j][k] = (k+1)*binSize + minIVal[j];

			}

		}

	}


	/*Equi-percentile binning*/
	if(binningType == 2){

		/*alloc and init*/
		Ivector_2 = (double **) malloc(maxNrOfInvariants*nrOfDBrecords*sizeof(double)); 

		for(j=0;j<maxNrOfInvariants;j++){ 

			allocateForNrOfDBrecords = (int) floor(1.05*nrOfDBrecords); //1.05* to be on the safe side
			Ivector_2[j] = (double *) malloc(allocateForNrOfDBrecords*sizeof(double)); 
					
			for(m=0; m<allocateForNrOfDBrecords;m++){
				Ivector_2[j][m] = 0.0; 
			}


		}


		/*first read the invariants' values into a pointer; then sort:*/
		for(n = 0; n < nrOfDBstructures; n++){

			for(i = 0; i < ptr_dbResults[n].nrOfWindows; i++){

				Ivector_2[0][cnt] =  ptr_dbResults[n].I12[i];
				Ivector_2[1][cnt] =  ptr_dbResults[n].Ia12[i];

				Ivector_2[2][cnt] =  ptr_dbResults[n].I1234_full[i];
				Ivector_2[3][cnt] =  ptr_dbResults[n].I1324_full[i];
				Ivector_2[4][cnt] =  ptr_dbResults[n].I1423[i];


				Ivector_2[5][cnt] =  ptr_dbResults[n].Ia12a34_full[i];
				Ivector_2[6][cnt] =  ptr_dbResults[n].Ia13a24_full[i];
				Ivector_2[7][cnt] =  ptr_dbResults[n].Ia14a23[i];

				Ivector_2[8][cnt] =  ptr_dbResults[n].Ia1234_full[i];
				Ivector_2[9][cnt] =  ptr_dbResults[n].Ia1324_full[i];
				Ivector_2[10][cnt] =  ptr_dbResults[n].Ia1423[i];

				Ivector_2[11][cnt] =  ptr_dbResults[n].I12a34_full[i];
				Ivector_2[12][cnt] =  ptr_dbResults[n].I13a24_full[i];
				Ivector_2[13][cnt] =  ptr_dbResults[n].I14a23[i];

				
				cnt += 1;
			}
		}

		/*Now sort the values for each entry and derive the bins:*/
		binCount = floor((double) nrOfDBrecords/(nrOfBins+1)); 
		for(j=0; j < nrOfEntries;j++){ 
			
			heapSort1dArray(Ivector_2[j], nrOfDBrecords, 0);

			/*loop over the nr of bins and record the bins (ie the bin divide points):*/
			for(k = 0; k < nrOfBins; k++){

				idx = (k+1)*binCount;
				bins[j][k] = Ivector_2[j][idx] + compensationForRounding; 
			}

		}
	}

	//if(binningType == 3){
	//}

	/*for(j=0; j < nrOfEntries;j++){ 
		for(k = 0; k < nrOfBins; k++){
			printf("For entry %d bin %d is: %lf\n", j, k, bins[j][k]);
			getchar();
		}
	}*/

	if(binningType == 2){
		//free(Ivector); .. not worth it
		freeDblArray(Ivector_2, maxNrOfInvariants);
	}

	*ptr_bins = bins;

	return returnVal;

}

/*Similar function for pairs of I-values; intended to be used for pairs of disjoint windows:*/
int generateBinsPairs(int binningType, struct db_I_windowPairs_order2_ptr *ptr_dbResultsPairs, int nrOfDBstructures, int nrOfDBrecords, int nrOfEntries, double ***ptr_bins, int nrOfBins, double compensationForRounding){

	int returnVal = 0;

	double **bins;
	int i, j, k, m, n;
	int allocateForNrOfDBrecords;

	double 	*minIVal, *maxIVal;

	double *Ivector; /*convenience ptr*/

	double **Ivector_2; /*another convenience ptr; to be used in the type 2 method*/
	int cnt = 0;
	int idx = 0;

	double binSize;

	int binCount;

	minIVal = (double *) malloc(nrOfEntries*sizeof(double));
	maxIVal = (double *) malloc(nrOfEntries*sizeof(double));

	Ivector = (double *) malloc(maxNrOfInvariants*sizeof(double));
	Ivector_2 = (double **) malloc(maxNrOfInvariants*sizeof(double*)); 

	/*Allocate memory to the bins pointer:*/
	bins = (double **) malloc(nrOfEntries*sizeof(double *));
	for(j = 0; j< nrOfEntries; j++){

		bins[j] = (double *) malloc(nrOfBins*sizeof(double));

		/*init*/
		minIVal[j] = 1e20;
		maxIVal[j] = -1e20;

		for(k = 0; k < nrOfBins;k++){

				bins[j][k] = -117.0;

		}

	}

	if(binningType == 0){ //Custom type; change bin definitions below when nrOf bins change!

	
		for(j = 0; j< nrOfEntries; j++){

			bins[j] = (double *) malloc(17*sizeof(double));

			//bins = {-2.0,-1.0, -0.5, 0.0, 0.5, 1.0, 2.0};
			/*bins[0] = -2;
			bins[1] = -1;
			bins[2] = -0.5;
			bins[3] = 0;
			bins[4] = 0.5;
			bins[5] = 1;
			bins[6] = 2;*/


			/*bins[0] = -2;
			bins[1] = -1.5;
			bins[2] = -1;
			bins[3] = -0.5;
			bins[4] = -0.25;
			bins[5] = 0.0;
			bins[6] = 0.25;
			bins[7] = 0.5;
			bins[8] = 1.0;
			bins[9] = 1.5;
			bins[10] = 2;*/

			/*bins[0] = -2.0;
			bins[1] = -1.5;
			bins[2] = -1;
			bins[3] = -0.75;
			bins[4] = -0.5;
			bins[5] = -0.25;
			bins[6] = -0.1;
			bins[7] = 0.0;
			bins[8] = 0.1;
			bins[9] = 0.25;
			bins[10] = 0.5;
			bins[11] = 0.75;
			bins[12] = 1.0;
			bins[13] = 1.5;
			bins[14] = 2.0;*/

			/*bins[j][0] = -2;
			bins[j][1] = -1.6;
			bins[j][2] = -1.2;
			bins[j][3] = -0.8;
			bins[j][4] = -0.4;
			bins[j][5] = 0.0;
			bins[j][6] = 0.4;
			bins[j][7] = 0.8;
			bins[j][8] = 1.2;
			bins[j][9] = 1.6;
			bins[j][10] = 2;*/

			bins[j][0] = -3.0;
			bins[j][1] = -2.5;
			bins[j][2] = -2;
			bins[j][3] = -1.5;
			bins[j][4] = -1;
			bins[j][5] = -0.5;
			bins[j][6] = -0.25;
			bins[j][7] = -0.1;
			bins[j][8] =  0.0;
			bins[j][9] = 0.1;
			bins[j][10] = 0.25;
			bins[j][11] = 0.5;
			bins[j][12] = 1;
			bins[j][13] = 1.5;
			bins[j][14] = 2;
			bins[j][15] = 2.5;
			bins[j][16] = 3.0;

			//bins[0] = -3.0;
			//bins[1] = -2.75;
			//bins[2] = -2.5;
			//bins[3] = -2.25;
			//bins[4] = -2.0;
			//bins[5] = -1.75;
			//bins[6] = -1.5;
			//bins[7] = -1.25;
			//bins[8] = -1.0;
			//bins[9] = -0.75;
			//bins[10] = -0.5;
			//bins[11] = -0.25;
			//bins[12] = 0;
			//bins[13] = 0.25;
			//bins[14] = 0.5;
			//bins[15] = 0.75;
			//bins[16] = 1.0;
			//bins[17] = 1.25;
			//bins[18] = 1.5;
			//bins[19] = 1.75;
			//bins[20] = 2.0;
			//bins[21] = 2.25;
			//bins[22] = 2.5;
			//bins[23] = 2.75;
			//bins[24] = 3.0;

			/*bins[0] = -3.0;
			bins[1] = -2.5;
			bins[2] = -2.0;
			bins[3] = -1.75;
			bins[4] = -1.5;
			bins[5] = -1.25;
			bins[6] = -1.0;
			bins[7] = -0.75;
			bins[8] = -0.5;
			bins[9] = -0.25;
			bins[10] = -0.125;
			bins[11] = 0;
			bins[12] = 0.125;
			bins[13] = 0.25;
			bins[14] = 0.5;
			bins[15] = 0.75;
			bins[16] = 1.0;
			bins[17] = 1.25;
			bins[18] = 1.5;
			bins[19] = 1.75;
			bins[20] = 2.0;
			bins[21] = 2.5;
			bins[22] = 3.0;*/

		}

	}

	/*Equi-distant binning*/
	if(binningType == 1){

		/*init*/
		for(j=0;j<maxNrOfInvariants;j++){ Ivector[j] = 0.0; }


		/*loop through the data, and record for each entry (invariant) its min and max value. Then  
		derive the binning:*/
		for(n = 0; n < nrOfDBstructures; n++){

			for(i = 0; i < ptr_dbResultsPairs[n].nrOfWindows; i++){

				for(j = i; j < ptr_dbResultsPairs[n].nrOfWindows; j++){

					Ivector[0] =  ptr_dbResultsPairs[n].I12[i][j];
					Ivector[1] =  ptr_dbResultsPairs[n].Ia12[i][j];

					Ivector[2] =  ptr_dbResultsPairs[n].I1234_full[i][j];
					Ivector[3] =  ptr_dbResultsPairs[n].I1324_full[i][j];
					Ivector[4] =  ptr_dbResultsPairs[n].I1423[i][j];

					Ivector[5] =  ptr_dbResultsPairs[n].Ia12a34_full[i][j];
					Ivector[6] =  ptr_dbResultsPairs[n].Ia13a24_full[i][j];
					Ivector[7] =  ptr_dbResultsPairs[n].Ia14a23[i][j];

					Ivector[8] =  ptr_dbResultsPairs[n].Ia1234_full[i][j];
					Ivector[9] =  ptr_dbResultsPairs[n].Ia1324_full[i][j];
					Ivector[10] =  ptr_dbResultsPairs[n].Ia1423[i][j];

					Ivector[11] =  ptr_dbResultsPairs[n].I12a34_full[i][j];
					Ivector[12] =  ptr_dbResultsPairs[n].I13a24_full[i][j];
					Ivector[13] =  ptr_dbResultsPairs[n].I14a23[i][j];

					for(k=0; k < nrOfEntries;k++){

						if(Ivector[k] < minIVal[k]){
							minIVal[k] = Ivector[k];
						}
						if(Ivector[k] > maxIVal[k]){
							maxIVal[k] = Ivector[k];
						}

					}

				}
			}
		}


		/*Derive the bins*/
		for(j=0; j < nrOfEntries;j++){

			//printf("(langmodig venter bysvalen stadig stædig længes) ...the min was: %lf and the max: %lf ...what on ..?\n", minIVal[j],maxIVal[j]);
			binSize = maxIVal[j] - minIVal[j];
			//printf(".. so the bins cover a range of: %lf\n", binSize);
			binSize /= nrOfBins-1;
			//printf(".. which is then divided into bins of size: %lf\n", binSize);

			//printf("bin %d size: %lf\n", j, binSize);
			//getchar();


			for(k = 0; k < nrOfBins;k++){

				bins[j][k] = k*binSize + minIVal[j]; /*lie this 1st bin-divide is at minIVal[j], last at maxIVal[j] */

			}

		}

	}


	/*Equi-percentile binning*/
	if(binningType == 2){

		/*init*/
		for(j=0;j<maxNrOfInvariants;j++){ 
			
				allocateForNrOfDBrecords = (int) floor(1.05*nrOfDBrecords); //1.05* to be on the safe side
				Ivector_2[j] = (double *) malloc(allocateForNrOfDBrecords*sizeof(double)); 
					
				for(m=0; m<allocateForNrOfDBrecords;m++){
					Ivector_2[j][m] = 0.0; 
				}

		}


		/*first read the invariants' values into a pointer; then sort:*/
		for(n = 0; n < nrOfDBstructures; n++){

			for(i = 0; i < ptr_dbResultsPairs[n].nrOfWindows; i++){

				for(j = i; j < ptr_dbResultsPairs[n].nrOfWindows; j++){

					Ivector_2[0][cnt] =  ptr_dbResultsPairs[n].I12[i][j];
					Ivector_2[1][cnt] =  ptr_dbResultsPairs[n].Ia12[i][j];

					Ivector_2[2][cnt] =  ptr_dbResultsPairs[n].I1234_full[i][j];
					Ivector_2[3][cnt] =  ptr_dbResultsPairs[n].I1324_full[i][j];
					Ivector_2[4][cnt] =  ptr_dbResultsPairs[n].I1423[i][j];


					Ivector_2[5][cnt] =  ptr_dbResultsPairs[n].Ia12a34_full[i][j];
					Ivector_2[6][cnt] =  ptr_dbResultsPairs[n].Ia13a24_full[i][j];
					Ivector_2[7][cnt] =  ptr_dbResultsPairs[n].Ia14a23[i][j];

					Ivector_2[8][cnt] =  ptr_dbResultsPairs[n].Ia1234_full[i][j];
					Ivector_2[9][cnt] =  ptr_dbResultsPairs[n].Ia1324_full[i][j];
					Ivector_2[10][cnt] =  ptr_dbResultsPairs[n].Ia1423[i][j];

					Ivector_2[11][cnt] =  ptr_dbResultsPairs[n].I12a34_full[i][j];
					Ivector_2[12][cnt] =  ptr_dbResultsPairs[n].I13a24_full[i][j];
					Ivector_2[13][cnt] =  ptr_dbResultsPairs[n].I14a23[i][j];

				
					cnt += 1;

				}
			}
		}
		
		//printf("Is this ok? cnt: %d nrOfDBrecords: %d\n", cnt, nrOfDBrecords);
		//getchar();

		/*Now sort the values for each entry and derive the bins:*/
		binCount = floor((double) nrOfDBrecords/(nrOfBins+1)); 
		for(j=0; j < nrOfEntries;j++){ 
			
			heapSort1dArray(Ivector_2[j], nrOfDBrecords, 0);

			/*loop over the nr of bins and record the bins (ie the bin divide points):*/
			for(k = 0; k < nrOfBins; k++){

				idx = (k+1)*binCount;
				bins[j][k] = Ivector_2[j][idx] + compensationForRounding; 
			}

		}
	}

	//if(binningType == 3){
	//}

	/*
	for(j=0; j < nrOfEntries;j++){ 
		for(k = 0; k < nrOfBins; k++){
			printf("For entry %d bin %d is: %lf\n", j, k, bins[j][k]);
			//getchar();
		}
	}
	*/


	//free(Ivector); .. not worth it
	if(binningType == 2){
		freeDblArray(Ivector_2, maxNrOfInvariants);
	}

	*ptr_bins = bins;

	return returnVal;

}


/*place each entry of an an array in a bin*/
int binArray(int *binnedVector, double *Ivector, int nrOfEntries, double **bins, int nrOfBins){

	int returnVal= 0;

	int i = 0;
	int m = 0;

	for(i = 0; i < nrOfEntries; i++){
		
		//printf("I-val at %d: %lf\n", i, Ivector[i]);

		m = 0; //reset
		binnedVector[i] = 0;
		while( m < nrOfBins && Ivector[i] > bins[i][m] ){ m +=1; }
		binnedVector[i] = m;
		//printf("%d 'th entry placed in %d\n",i, binnedVector[i]);
	}

	return returnVal;

}



/*Distance fcts:*/
double distEuklid(double *ptr_x, double *ptr_y, int length){
	
	int i = 0;
	double d = 0;

	for(i=0; i<length;i++){

		d += pow(ptr_y[i] - ptr_x[i],2);
	}

	//d = sqrt(d);

	return d;
}

double relDistEuklid(double *ptr_x,  double *ptr_y, int length){
	
	int i = 0;
	double d = 0;

	for(i=0; i<length;i++){

		d += pow((ptr_y[i] - ptr_x[i])/ptr_x[i],2);
	}

	d = sqrt(d);

	return d;
}


/*Fcts for the rarity method:*/ 
void collect_Irarity_windows(struct Irarity_windows *ptr_Irarity_windows, struct I_windows_ptr *ptr_I_windows, int windowIndex, double rarityScore){


	/*strcpy(Irarity_windows.fileName, ptr_I_windows -> fileName);
	strcpy(Irarity_windows.structureName, ptr_I_windows -> structureName);
	strcpy(Irarity_windows.classId, ptr_I_windows -> classId);
	strcpy(Irarity_windows.chainId, ptr_I_windows -> chainId);
	Irarity_windows.chainNr = ptr_I_windows -> chainNr;
	Irarity_windows.chainLen = ptr_I_windows -> chainLen;
	Irarity_windows.order = ptr_I_windows -> order;
	Irarity_windows.nrOfWindows = ptr_I_windows -> nrOfWindows;
	Irarity_windows.windowLgth = ptr_I_windows -> windowLgth;
	Irarity_windows.stepSize = ptr_I_windows -> stepSize;*/

		
	/*ptr_Irarity_windows[windowIndex].windowsNr = ptr_I_windows -> windowsNr[windowIndex];
	ptr_Irarity_windows[windowIndex].segIndices[0] = ptr_I_windows -> segIndices[windowIndex][0];
	ptr_Irarity_windows[windowIndex].segIndices[1] = ptr_I_windows -> segIndices[windowIndex][1];
*/

	ptr_Irarity_windows[windowIndex].window.windowNr = ptr_I_windows -> window[windowIndex].windowNr;
	ptr_Irarity_windows[windowIndex].window.segIndices[0] = ptr_I_windows -> window[windowIndex].segIndices[0];
	ptr_Irarity_windows[windowIndex].window.segIndices[1] = ptr_I_windows -> window[windowIndex].segIndices[1];

	ptr_Irarity_windows[windowIndex].rarityScore = rarityScore;


}

void collect_Irarity_windowPairs(struct Irarity_windowPairs *ptr_Irarity_windowPairs, struct I_windowPairs_ptr *ptr_I_windowPairs, int windowPairsIndex, int windowIndex_1, int windowIndex_2, double rarityScorePairs){

	/*ptr_Irarity_windowPairs[windowPairsIndex].windowsNr_1 = ptr_I_windowPairs -> windowsNr_1[windowIndex_1];
	ptr_Irarity_windowPairs[windowPairsIndex].segIndices_1[0] = ptr_I_windowPairs -> segIndices_1[windowIndex_1][0];
	ptr_Irarity_windowPairs[windowPairsIndex].segIndices_1[1] = ptr_I_windowPairs -> segIndices_1[windowIndex_1][1];

	ptr_Irarity_windowPairs[windowPairsIndex].windowsNr_2 = ptr_I_windowPairs -> windowsNr_2[windowIndex_2];
	ptr_Irarity_windowPairs[windowPairsIndex].segIndices_2[0] = ptr_I_windowPairs -> segIndices_2[windowIndex_2][0];
	ptr_Irarity_windowPairs[windowPairsIndex].segIndices_2[1] = ptr_I_windowPairs -> segIndices_2[windowIndex_2][1];
	*/
	ptr_Irarity_windowPairs[windowPairsIndex].windowPair.windowNr_1 = ptr_I_windowPairs -> windowPair[windowIndex_1][windowIndex_2].windowNr_1;
	ptr_Irarity_windowPairs[windowPairsIndex].windowPair.segIndices_1[0] = ptr_I_windowPairs -> windowPair[windowIndex_1][windowIndex_2].segIndices_1[0];
	ptr_Irarity_windowPairs[windowPairsIndex].windowPair.segIndices_1[1] = ptr_I_windowPairs -> windowPair[windowIndex_1][windowIndex_2].segIndices_1[1];

	ptr_Irarity_windowPairs[windowPairsIndex].windowPair.windowNr_2 = ptr_I_windowPairs -> windowPair[windowIndex_1][windowIndex_2].windowNr_2;
	ptr_Irarity_windowPairs[windowPairsIndex].windowPair.segIndices_2[0] = ptr_I_windowPairs -> windowPair[windowIndex_1][windowIndex_2].segIndices_2[0];
	ptr_Irarity_windowPairs[windowPairsIndex].windowPair.segIndices_2[1] = ptr_I_windowPairs -> windowPair[windowIndex_1][windowIndex_2].segIndices_2[1];

	ptr_Irarity_windowPairs[windowPairsIndex].rarityScore = rarityScorePairs;


}

void collect_queryRawScore(struct queryRawScore *ptr_queryRawScore, int pairs_b, struct I_windows_ptr I_windows, struct I_windowPairs_ptr I_windowPairs, char *windowInfo, double score, double pValue, int structureIndex, double optionalValue){

	if(pairs_b ==1){

		strcpy(ptr_queryRawScore[structureIndex].fileName, I_windowPairs.fileName);
		strcpy(ptr_queryRawScore[structureIndex].structureName, I_windowPairs.structureName);
		strcpy(ptr_queryRawScore[structureIndex].chainId, I_windowPairs.chainId);
		strcpy(ptr_queryRawScore[structureIndex].classId, I_windowPairs.classId);
		strcpy(ptr_queryRawScore[structureIndex].windowInfo, windowInfo);
		ptr_queryRawScore[structureIndex].nrOfWindows = I_windowPairs.nrOfWindows;
		ptr_queryRawScore[structureIndex].chainLen = I_windowPairs.chainLen;
		ptr_queryRawScore[structureIndex].chainNr = I_windowPairs.chainNr;

		ptr_queryRawScore[structureIndex].score = score;
		ptr_queryRawScore[structureIndex].pValue = pValue;

		ptr_queryRawScore[structureIndex].optionalValue = optionalValue;

	}
	else if(pairs_b == 0){
		strcpy(ptr_queryRawScore[structureIndex].fileName, I_windows.fileName);
		strcpy(ptr_queryRawScore[structureIndex].structureName, I_windows.structureName);
		strcpy(ptr_queryRawScore[structureIndex].chainId, I_windows.chainId);
		strcpy(ptr_queryRawScore[structureIndex].classId, I_windows.classId);
		strcpy(ptr_queryRawScore[structureIndex].windowInfo, windowInfo);
		ptr_queryRawScore[structureIndex].nrOfWindows = I_windows.nrOfWindows;
		ptr_queryRawScore[structureIndex].chainLen = I_windows.chainLen;
		ptr_queryRawScore[structureIndex].chainNr = I_windows.chainNr;

		ptr_queryRawScore[structureIndex].score = score;
		ptr_queryRawScore[structureIndex].pValue = pValue;

		ptr_queryRawScore[structureIndex].optionalValue = optionalValue;
	}

}

/*score topN window/window pair matches and update results:*/
void simpleScore(int queryLength, int normalizeWithQueryLength_b, int queryWindowNr, int queryNrOfWindows, int windowLength, int stepSize, struct matchScore *ptr_Scores, int nrOfDbStructures, struct matchCandidate *ptr_matchCandidates, int topN, int lengthStructureName, float rarityFactor){
	
	double scoreBaseVal = 2;
	double logBaseVal = 4;
	double offset = 0.1;

	int i, j, m, k = 0;
	int hit = 0;

	int nrOfWindows;
	int nrOfWindowsOverlap;
	double lengthRatio;
	double logLengthRatio;
	double ratio;
	int rawCountScore;
	double countScore;
	double rawDistScore;
	double distScore;
	double qLengthModdedDistScore;

	int hitWindowAlreadyUsed_b = 0;

	nrOfWindowsOverlap = ((int) floor((float) windowLength/stepSize));

	for(i=0;i < topN; i++){

		nrOfWindows = ptr_matchCandidates[i].nrOfWindows;

		rawCountScore = pow(scoreBaseVal, topN - i);

		lengthRatio = (double) queryLength/ptr_matchCandidates[i].chainLen;

		logLengthRatio = logbase(lengthRatio,logBaseVal);

		ratio = (double) 1/(1 + pow(logLengthRatio,2));

		//ratio = (double) 1/nrOfWindows;

		countScore = (double) rawCountScore*ratio;

		/*we let the inverse of the dist be the score; to cap the score we put an offset in 
		(this will also prevent the score from blowing up). The motivation is to let bad dist-hits
		get a small score, and to some extent independently of how bad it is, while promoting fine
		dist-hits. The typical situation is that a good candidate matches a query in some stretch of the 
		chain and this we want to capture*/
		distScore = (double) (1/(offset + ptr_matchCandidates[i].distMin))*ratio; //- log10(rarityFactor)

		ptr_matchCandidates[i].score = distScore;

		/*search through all candidates and find the match, if any. */
		//printf("Score fct, match cand str name: %s\n", ptr_matchCandidates[i].structureName);
		hit = bisectionalSearchMatchScores(ptr_matchCandidates[i].structureName, ptr_Scores, nrOfDbStructures, lengthStructureName); //since the structure name includes both the pdb id and the chain id, matching on it is equiv to: strcmp(ptr_Scores[m].matchFileName,ptr_matchCandidates[i].fileName) ==0  && strcmp(ptr_Scores[m].matchChainId, ptr_matchCandidates[i].chainId) == 0

		if(hit != -1){ //if hit == -1 no match was found

			if(strcmp(ptr_matchCandidates[i].structureName,ptr_Scores[hit].matchStructureName) != 0 || strcmp(ptr_matchCandidates[i].chainId, ptr_Scores[hit].matchChainId) != 0){
				printf("STOPSTOPSTOP: cand: %s match id:%s chains %s and %s \n", ptr_matchCandidates[i].structureName, ptr_Scores[hit].matchStructureName, ptr_matchCandidates[i].chainId, ptr_Scores[hit].matchChainId );
				getchar();
			}

			/*normalize by query length if desired:*/
			qLengthModdedDistScore = distScore;
			if(normalizeWithQueryLength_b == 1){ 
				qLengthModdedDistScore = distScore/queryLength;
			}
				
				
			/*A given hit (and q-window) can have several windows matching the q-window well, ie the same
				structure may occur several times in the ptr_matchCandidates-array. We only want to let the best 
				of these hits contribute to the hit's total score, however only if that hit was not already used to 
				hit another query-window:*/

			/*check if the matching window was used already:*/
			hitWindowAlreadyUsed_b = 0;
			for(k = 0; k < queryNrOfWindows; k++){
					
				/*if(ptr_Scores[hit].windowNr[k] == ptr_matchCandidates[i].windowsNr){
					hitWindowAlreadyUsed_b = 1;
				}*/
				if(ptr_Scores[hit].windowNr[k] == ptr_matchCandidates[i].window.windowNr){
					hitWindowAlreadyUsed_b = 1;
				}

			}

			if(ptr_Scores[hit].windowNr[queryWindowNr] == -1 && hitWindowAlreadyUsed_b == 0){

				ptr_Scores[hit].matchNrOfWindows = nrOfWindows;

				ptr_Scores[hit].countScore += countScore;

				ptr_Scores[hit].aggrWinScore += qLengthModdedDistScore;

				ptr_Scores[hit].places[i] +=1;

				/*compute the coverage so far:*/
				/*This is slightly tricky: for queryWindowNr >= nrOfWindowsOverlap the clause
				"queryWindowNr - j >= -nrOfWindowsOverlap +1  && j < nrOfWindowsOverlap" is eq to
				"j < nrOfWindowsOverlap" and queryWindowNr - j >= 0  is automatic, why the max at
				the end is superfluous too; for queryWindowNr < nrOfWindowsOverlap the situation is
				though different since we cannot go further back with queryWindowNr - j than to 0.
				So queryWindowNr = 0, 1, ..., nrOfWindowsOverlap -1 needs special care, and this is
				the reason for the somewhat odd clause "queryWindowNr - j >= -nrOfWindowsOverlap +1"
				and the max at the end*/

					
				/*if(queryWindowNr == 0){

					ptr_Scores[hit].queryCoverage = 1;

				}
				else{*/


				j=1;
				while(queryWindowNr - j >= -nrOfWindowsOverlap +1  && j < nrOfWindowsOverlap && ptr_Scores[hit].windowNr[max(queryWindowNr-j,0)] == -1){j +=1;} 
						
				ptr_Scores[hit].queryCoverage += j*((float) stepSize/windowLength);

				//}

				ptr_Scores[hit].score += ptr_Scores[hit].queryCoverage*qLengthModdedDistScore;

				ptr_Scores[hit].queryWindowNr[queryWindowNr] = queryWindowNr;
				ptr_Scores[hit].windowNr[queryWindowNr] = ptr_matchCandidates[i].window.windowNr;
				//ptr_Scores[hit].windowNr[queryWindowNr] = ptr_matchCandidates[i].windowNr;
			
			}

		}
		else{ 
			printf("Warning: no match with this structure name: %s\n", ptr_matchCandidates[i].structureName);
			//getchar();
		}

	}

}

void simpleScore_windowPairs(int queryLength, struct matchScore *ptr_Scores, int nrOfDbStructures, struct matchCandidate_windowPair *ptr_matchCandidates_windowPairs, int topN, int lengthStructureName){
	
	double scoreBaseVal = 2;
	double logBaseVal = 4;
	double offset = 0.1;

	int i, m = 0;
	int hit = 0;

	int nrOfWindows;
	double lengthRatio;
	double logLengthRatio;
	double ratio;
	int rawCountScore;
	double countScore;
	double rawDistScore;
	double distScore;

	int print_warning_b = 0;

	for(i=0;i < topN; i++){

		nrOfWindows = ptr_matchCandidates_windowPairs[i].nrOfWindows;

		rawCountScore = pow(scoreBaseVal, topN - i);

		lengthRatio = (double) queryLength/ptr_matchCandidates_windowPairs[i].chainLen;

		logLengthRatio = 2*logbase(lengthRatio,logBaseVal);

		ratio = (double) 1/(1 + pow(logLengthRatio,2));

		//ratio = (double) 1/nrOfWindows;

		countScore = (double) rawCountScore*ratio;

		/*we let the inverse of the dist be the score; to cap the score we put an offset in 
		(this will also prevent the score from blowing up). The motivation is to let bad dist-hits
		get a small score, and to some extent independently of how bad it is, while promoting fine
		dist-hits. The typical situation is that a good candidate matches a query in some stretch of the 
		chain and this we want to capture*/
		distScore = (double) 1/(offset + ptr_matchCandidates_windowPairs[i].distMin)*ratio;

		/*search through all candidates and find the match, if any. */
		//printf("Score fct, match cand str name: %s\n", ptr_matchCandidates[i].structureName);
		hit = bisectionalSearchMatchScores(ptr_matchCandidates_windowPairs[i].structureName, ptr_Scores, nrOfDbStructures, lengthStructureName); //since the structure name icnludes both the pdb id and the chain id, matching on it is equiv to: strcmp(ptr_Scores[m].matchFileName,ptr_matchCandidates[i].fileName) ==0  && strcmp(ptr_Scores[m].matchChainId, ptr_matchCandidates[i].chainId) == 0

		if(hit != -1){ //if hit == -1 no match was found
				
				ptr_Scores[hit].matchNrOfWindows = nrOfWindows;

				ptr_Scores[hit].countScore += countScore;

				ptr_Scores[hit].score += distScore;

				ptr_Scores[hit].places[i] +=1;
		}
		else{ 
			if(print_warning_b == 1){
				printf("Warning: no match with this structure name: %s\n", ptr_matchCandidates_windowPairs[i].structureName);
			}
		}

	}

}


/*older version not using bisectional search*/
void simpleScore_1(int queryLength, struct matchScore *ptr_scores, int nrOfDbStructures, struct matchCandidate *ptr_matchCandidates, int topN){
	
	double scoreBase = 2;
	double offset = 0.1;

	int i, m = 0;
	int hit = 0;

	int nrOfWindows;
	double lengthRatio;
	double logLengthRatio;
	double ratio;
	int rawCountScore;
	double countScore;
	double rawDistScore;
	double distScore;

	for(i=0;i < topN; i++){

		nrOfWindows = ptr_matchCandidates[i].nrOfWindows;

		rawCountScore = pow(scoreBase, topN - i);

		lengthRatio = (double) queryLength/ptr_matchCandidates[i].chainLen;

		logLengthRatio = logbase(lengthRatio,4);

		ratio = (double) 1/(1 + pow(logLengthRatio,2));

		//ratio = (double) 1/nrOfWindows;

		countScore = (double) rawCountScore*ratio;

		/*we let the inverse of the dist be the score; to cap the score we put an offset in 
		(this will also prevent the score from blowing up). The motivation is to let bad dist-hits
		get a small score, and to some extent independently of how bad it is, while promoting fine
		dist-hits. The typical situation is that a good candidate matches a query in some stretch of the 
		chain and this we want to capture*/
		distScore = (double) 1/(offset + ptr_matchCandidates[i].distMin)*ratio;



		/*run through all candidates and find the match, if any. 
		COULD BE DONE MUCH FASTER BY BISECTION IN LEXI ORDERED ptr_Scores! Takes removing the path-name from the fileNames though */
		for(m = 0; m < nrOfDbStructures; m++){
		
			//seek hit, score it and update the results ptr:
			if(strcmp(ptr_scores[m].matchFileName,ptr_matchCandidates[i].fileName) ==0  && strcmp(ptr_scores[m].matchChainId, ptr_matchCandidates[i].chainId) == 0){
				
				ptr_scores[m].matchNrOfWindows = nrOfWindows;

				ptr_scores[m].countScore += countScore;

				ptr_scores[m].score += distScore;

				ptr_scores[m].places[i] +=1;


			}

		}
	}

}

/*scoring fct's using the distribution of the dist's throughout the data base*/
void distDistrScore(double *ptr_distDistr, int sizeDistDistr, int queryLength, int normalizeWithQueryLength_b, int queryWindowNr, int queryNrOfWindows, int windowLength, int stepSize, struct matchScore *ptr_Scores, int nrOfDbStructures, struct matchCandidate *ptr_matchCandidates, int topN, int lengthStructureName, int lengthMod_b){
	
	double scoreBase = 2;
	int logBase = 4;
	double offset = 0.1;

	int i, j, k = 0;
	int hitDistValue;
	double pValue;

	int hit = 0;

	int nrOfWindows;
	int nrOfWindowsOverlap;
	double lengthRatio;
	double logLengthRatio;
	double ratio;
	int rawCountScore;
	double countScore;
	double rawDistScore;
	double distScore;
	double qLengthModdedDistScore;

	int hitWindowAlreadyUsed_b = 0;

	nrOfWindowsOverlap = ((int) floor((float) windowLength/stepSize));

	//score by means of distribution of the dist's occuring when searching the db against itself.
	for(i=0;i < topN; i++){

		nrOfWindows = ptr_matchCandidates[i].nrOfWindows;

		rawCountScore = pow(scoreBase, topN - i);

		lengthRatio = (double) queryLength/ptr_matchCandidates[i].chainLen;

		logLengthRatio = logbase(lengthRatio,logBase);

		/*Ex's of following ratio: if logBase = 2, and the query is twice (half) as long as 
		the match, the lengthRatio is 2 (1/2) and the logLengthRatio is 1 (-1). The ratio
		then becomes 1/5. If logBase = 4 the ratio is 1/5 if the query is 4 (1/4) times as long
		as the match*/
		ratio = (double) 1/(1 + pow(logLengthRatio,2));

		//ratio = (double) 1/nrOfWindows;

		/*a simple count score:*/
		countScore = (double) rawCountScore*ratio;

		/*To score the hit (properly) we look up the distMin value in the distribution, ptr_distDistr, and 
		thereby fetch its position in that (sorted array); from there the "percentile" (p-value)
		is obtained just by dividing by the size of the distribution (ie the length of the ptr_distDistr
		array). The score is then defined as -log10(p-value) divided by the query-length (so as to normalize
		for query length):*/
		hitDistValue = bisectionalSearchValueL(ptr_matchCandidates[i].distMin, ptr_distDistr, sizeDistDistr);

		pValue = (double) (hitDistValue +1)/sizeDistDistr; //hitDistValue is the left end of the interval in which distMin sits; we add 1 to avoid numm disturbance

		distScore = (double) -log10(pValue);

		ptr_matchCandidates[i].score = distScore;

		if(lengthMod_b == 1){

			distScore = distScore*ratio;

		}

		/*Search through all candidates and find the match, if any (clearly, there ought to be one). */
		hit = bisectionalSearchMatchScores(ptr_matchCandidates[i].structureName, ptr_Scores, nrOfDbStructures, lengthStructureName); //since the structure name icnludes both the pdb id and the chain id, matching on it is equiv to: strcmp(ptr_Scores[m].matchFileName,ptr_matchCandidates[i].fileName) ==0  && strcmp(ptr_Scores[m].matchChainId, ptr_matchCandidates[i].chainId) == 0

		if(hit != -1){ //if hit == -1 no match was found

			/*normalize by query length if desired:*/
			qLengthModdedDistScore = distScore;
			if(normalizeWithQueryLength_b == 1){ 
				qLengthModdedDistScore = distScore/queryLength;
			}
				
				
			/*A given hit (and q-window) can have several windows matching the q-window well, ie the same
			structure may occur several times in the ptr_matchCandidates-array. We only want to let the best 
			of these hits contribute to the hit's total score, however only if that hit was not already used to 
			hit another query-window:*/

			/*check if the matching window was used already:*/
			hitWindowAlreadyUsed_b = 0;
			for(k = 0; k < queryNrOfWindows; k++){
					
				/*if(ptr_Scores[hit].windowNr[k] == ptr_matchCandidates[i].windowsNr){
					hitWindowAlreadyUsed_b = 1;
				}*/
				if(ptr_Scores[hit].windowNr[k] == ptr_matchCandidates[i].window.windowNr){
					hitWindowAlreadyUsed_b = 1;
				}

			}

			if(ptr_Scores[hit].windowNr[queryWindowNr] == -1 && hitWindowAlreadyUsed_b == 0){

				ptr_Scores[hit].matchNrOfWindows = nrOfWindows;

				ptr_Scores[hit].countScore += countScore;

				ptr_Scores[hit].aggrWinScore += distScore/queryLength;

				ptr_Scores[hit].places[i] +=1;

				/*compute the coverage so far:*/
				/*compute the coverage so far:*/
				/*This is slightly tricky: for queryWindowNr >= nrOfWindowsOverlap the clause
				"queryWindowNr - j >= -nrOfWindowsOverlap +1  && j < nrOfWindowsOverlap" is eq to
				"j < nrOfWindowsOverlap" and queryWindowNr - j >= 0  is automatic, why the max at
				the end is superfluous too; for queryWindowNr < nrOfWindowsOverlap the situation is
				though different since we cannot go further back with queryWindowNr - j than to 0.
				So queryWindowNr = 0, 1, ..., nrOfWindowsOverlap -1 needs special care, and this is
				the reason for the somewhat odd clause "queryWindowNr - j >= -nrOfWindowsOverlap +1"
				and the max at the end*/

					
				/*if(queryWindowNr == 0){

					ptr_Scores[hit].queryCoverage = 1;

				}
				else{*/


				j=1;
				while(queryWindowNr - j >= -nrOfWindowsOverlap +1  && j < nrOfWindowsOverlap && ptr_Scores[hit].windowNr[max(queryWindowNr-j,0)] == -1){j +=1;} 
						
				ptr_Scores[hit].queryCoverage += j*((float) stepSize/windowLength);

				//}

				ptr_Scores[hit].score += ptr_Scores[hit].queryCoverage*distScore/queryLength;

				ptr_Scores[hit].queryWindowNr[queryWindowNr] = queryWindowNr;
				ptr_Scores[hit].windowNr[queryWindowNr] = ptr_matchCandidates[i].window.windowNr;
				//ptr_Scores[hit].windowNr[queryWindowNr] = ptr_matchCandidates[i].windowNr;

			}

		}
		else{ 
			printf("Warning: this structure name: %s was not found in the scores list (ptr_Score)\n", ptr_matchCandidates[i].structureName);
		}
	}


}

void distDistrScore_windowPairs(double *ptr_distDistr, int sizeDistDistr, int queryLength, struct matchScore *ptr_Scores, int nrOfDbStructures, struct matchCandidate_windowPair *ptr_matchCandidates_windowPairs, int topN, int lengthStructureName, int lengthMod_b){
	
	double scoreBase = 2;
	int logBase = 4;
	double offset = 0.1;

	int i = 0;
	int hitDistValue;
	double pValue;

	int hit = 0;

	int nrOfWindows;
	double lengthRatio;
	double logLengthRatio;
	double ratio;
	int rawCountScore;
	double countScore;
	double rawDistScore;
	double distScore;

	//score by means of 
	for(i=0;i < topN; i++){

		nrOfWindows = ptr_matchCandidates_windowPairs[i].nrOfWindows;

		rawCountScore = pow(scoreBase, topN - i);

		lengthRatio = (double) queryLength/ptr_matchCandidates_windowPairs[i].chainLen;

		logLengthRatio = logbase(lengthRatio,logBase);

		/*Ex's of following ratio: if logBase = 2, and the query is twice (half) as long as 
		the match, the lengthRatio is 2 (1/2) and the logLengthRatio is 1 (-1). The ratio
		then becomes 1/5. If logBase = 4 the ratio is 1/5 if the query is 4 (1/4) times as long
		as the match*/
		ratio = (double) 1/(1 + pow(logLengthRatio,2));

		//ratio = (double) 1/nrOfWindows;

		/*a simple count score:*/
		countScore = (double) rawCountScore*ratio;

		/*To score the hit (properly) we look up the distMin value in the distribution, ptr_distDistr, and 
		thereby fetch its position in that (sorted array); from there the "percentile" (p-value)
		is obtained just by dividing by the size of the distribution (ie the length of the ptr_distDistr
		array). The score is then defined as -log10(p-value) divided by the query-length (so as to normalize
		for query length):*/
		hitDistValue = bisectionalSearchValueL(ptr_matchCandidates_windowPairs[i].distMin, ptr_distDistr, sizeDistDistr);

		pValue = (double) (hitDistValue +1)/sizeDistDistr; //hitDistValue is the left end of the interval in which distMin sits; we add 1 to avoid numm disturbance

		distScore = (double) -log10(pValue)/queryLength;

		if(lengthMod_b == 1){

			distScore = distScore*ratio;

		}

		/*Search through all candidates and find the match, if any (clearly, there ought to be one). */
		hit = bisectionalSearchMatchScores(ptr_matchCandidates_windowPairs[i].structureName, ptr_Scores, nrOfDbStructures, lengthStructureName); //since the structure name icnludes both the pdb id and the chain id, matching on it is equiv to: strcmp(ptr_Scores[m].matchFileName,ptr_matchCandidates[i].fileName) ==0  && strcmp(ptr_Scores[m].matchChainId, ptr_matchCandidates[i].chainId) == 0

		if(hit != -1){ //if hit == -1 no match was found
				
				ptr_Scores[hit].matchNrOfWindows = nrOfWindows;

				ptr_Scores[hit].countScore += countScore;

				ptr_Scores[hit].score += distScore;

				ptr_Scores[hit].places[i] +=1;
		}
		else{ 
			printf("Warning: this structure name: %s was not found in the scores list (ptr_Score)\n", ptr_matchCandidates_windowPairs[i].structureName);
		}
	}


}

/*Scoring fct's using an approximation of the frequency of the obtained hits across the their dist's*/
//void scoreFct(double *ptr_distDistr, int sizeDistDistr, int queryLength, struct matchScore *ptr_Scores, int nrOfDbStructures, int nrOfDbRecords, struct matchCandidate *ptr_matchCandidates, int topN, int lengthStructureName, int lengthMod_b, struct floatArrayWLength epsArray, struct floatArrayWLength frqArray){
void scoreFct(double *ptr_distDistr, int sizeDistDistr, int queryLength, int normalizeWithQueryLength_b, int queryWindowNr, int queryNrOfWindows, int windowLength, int stepSize,  struct matchScore *ptr_Scores, int nrOfDbStructures, int nrOfDbRecords, struct matchCandidate *ptr_matchCandidates, int topN, int lengthStructureName, int lengthMod_b, struct floatArrayWLength epsArray, struct floatArrayWLength frqArray){
	
	double scoreBaseVal = 2;
	double logBaseVal = 4;
	double offset = 0.1;

	int i, j, k, m = 0;
	int hit = 0;

	int nrOfWindowsOverlap; 

	int nrOfWindows;
	double lengthRatio;
	double logLengthRatio;
	double ratio;
	int rawCountScore;
	double countScore;
	double distScore;
	double qLengthModdedDistScore;

	int hitDistValue;
	double pValue = 1.0;

	int hitWindowAlreadyUsed_b = 0;

	nrOfWindowsOverlap = ((int) floor((float) windowLength/stepSize));

	/*reset the frequency frqArray:*/
	for(k=0; k < frqArray.length; k++){

		frqArray.floatArray[k] = 0.0;

	}

	/*first find the actual number of occurences of hits up to each level in the epsList*/
	for(i=0; i < topN; i++){

		j = 0;
		while(j < epsArray.length && ptr_matchCandidates[i].distMin > epsArray.floatArray[j]){j += 1;}

		//build the cumulative frequencies:
		for(k=j; k<frqArray.length; k++){frqArray.floatArray[k] += 1;}

	}

		
	//score by means of distribution of the dist's occuring when searching the db against itself. And: the probability
	// of getting a hit at the obtained dist for the qiven query 
	for(i=0;i < topN; i++){

		nrOfWindows = ptr_matchCandidates[i].nrOfWindows;

		rawCountScore = pow(scoreBaseVal, topN - i);

		lengthRatio = (double) queryLength/ptr_matchCandidates[i].chainLen;

		logLengthRatio = logbase(lengthRatio,logBaseVal);

		/*Ex's of following ratio: if logBase = 2, and the query is twice (half) as long as 
		the match, the lengthRatio is 2 (1/2) and the logLengthRatio is 1 (-1). The ratio
		then becomes 1/5. If logBase = 4 the ratio is 1/5 if the query is 4 (1/4) times as long
		as the match*/
		ratio = (double) 1/(1 + pow(logLengthRatio,2));

		//ratio = (double) 1/nrOfWindows;

		/*a simple count score:*/
		countScore = (double) rawCountScore*ratio;

		/*To score the hit (properly??) we use
		
		p(q' hit to q at dist) = p(dist|q' hit to q)*p(q' hit to q) 

		To get p(dist|q' hit to q): for a range of d values (epsRange) count the frequency of hits with dist < d 
		among the hits (ie the ones send to this fct). For a given value of dist the sought prob is then
		estimated as frq(dist)/#hits.
		An estimate of p(q' hit to q) is #hits/#all fragments.

		All in all we then get p(q' hit to q at dist) = frq(dist)/#all fragments ?OK?

		
		Wrong:
		p(q' hit to q at dist < d) = p(q' hit to q | dist < d)*p(dist) 
		To get p(dist) we look up the distMin value in the distribution, ptr_distDistr, and 
		thereby fetch its position in that sorted array; from there the "percentile" (p-value)
		is obtained just by dividing by the size of the distribution (ie the length of the ptr_distDistr
		array). Secondly, p(q' hit to q | dist < d) is estimated from the frq's in the actual list of hits
		to the query, q. BUT THAT'S NO GOOD: TO GET THIS PROB YOU RATHER NEED: THE SET OF ALL HITS TO A QUERY 
		WITH DIST < D AND THEN, AMONG THESE, TO FIND THE ONES BEING A HIT TO Q. AND WE DON'T HAVE THIS INFO. 
		
		The score is then defined 
		as -log10(p-value) divided by the query-length (so as to normalize for query length thereby allowing to compare the scores for different queries). 
		*/
		
		hitDistValue = bisectionalSearchValueL(ptr_matchCandidates[i].distMin, ptr_distDistr, sizeDistDistr);

		pValue = ((double) (hitDistValue +1)/sizeDistDistr); //hitDistValue is the left end of the interval in which distMin sits; we add 1 to avoid numm disturbance

		/*printf("p value %lf\n", pValue);
		getchar();*/

		/*next find estimate of p(dist|q' hit to q)*/
		j = 0;
		while(j < epsArray.length && ptr_matchCandidates[i].distMin > epsArray.floatArray[j]){j += 1;}

		pValue = ((double) frqArray.floatArray[j]/nrOfDbRecords)*pValue;

		/*if(pValue == 0.0){

			for(k=0;k<frqArray.length;k++){
				printf("frqArray %d: %lf\n",j,frqArray.floatArray[j]);
				getchar();
			}
			
		}*/

		distScore = ((double) -log10(pValue)); 

		//if(distScore > ptr_matchCandidates[i].score){ //there might already be a score > 0 recorded, since the same structure can have several hits for a given q-window
		ptr_matchCandidates[i].score = distScore;
		//}

		if(lengthMod_b == 1){

			distScore = distScore*ratio;

		}

		/*Search through all candidates and find the (index of the) match, if any (clearly, there ought to be one). */
		hit = bisectionalSearchMatchScores(ptr_matchCandidates[i].structureName, ptr_Scores, nrOfDbStructures, lengthStructureName); //since the structure name icnludes both the pdb id and the chain id, matching on it is equiv to: strcmp(ptr_Scores[m].matchFileName,ptr_matchCandidates[i].fileName) ==0  && strcmp(ptr_Scores[m].matchChainId, ptr_matchCandidates[i].chainId) == 0

		if(hit != -1){ //if hit == -1 no match was found

			/*normalize by query length if desired:*/
			qLengthModdedDistScore = distScore;
			if(normalizeWithQueryLength_b == 1){ 
				qLengthModdedDistScore = distScore/queryLength;
			}
				
			/*A given hit (and q-window) can have several windows matching the q-window well, ie the same
			structure may occur several times in the ptr_matchCandidates-array. We only want to let the best 
			of these hits contribute to the hit's total score, however only if that hit was not already used to 
			hit another query-window:*/

			/*check if the matching window was used already:*/
			hitWindowAlreadyUsed_b = 0;
			for(k = 0; k < queryNrOfWindows; k++){
					
				/*if(ptr_Scores[hit].windowNr[k] == ptr_matchCandidates[i].windowsNr){
					hitWindowAlreadyUsed_b = 1;
				}*/
				if(ptr_Scores[hit].windowNr[k] == ptr_matchCandidates[i].window.windowNr){
					hitWindowAlreadyUsed_b = 1;
				}

			}

			if(ptr_Scores[hit].queryWindowNr[queryWindowNr] == -1 && hitWindowAlreadyUsed_b == 0){

				ptr_Scores[hit].matchNrOfWindows = nrOfWindows;

				ptr_Scores[hit].countScore += countScore;

				ptr_Scores[hit].aggrWinScore += qLengthModdedDistScore;

				ptr_Scores[hit].places[i] +=1;

				/*compute the coverage so far:*/
				/*This is slightly tricky: for queryWindowNr >= nrOfWindowsOverlap the clause
				"queryWindowNr - j >= -nrOfWindowsOverlap +1  && j < nrOfWindowsOverlap" is eq to
				"j < nrOfWindowsOverlap" and queryWindowNr - j >= 0  is automatic, why the max at
				the end is superfluous too; for queryWindowNr < nrOfWindowsOverlap the situation is
				though different since we cannot go further back with queryWindowNr - j than to 0.
				So queryWindowNr = 0, 1, ..., nrOfWindowsOverlap -1 needs special care, and this is
				the reason for the somewhat odd clause "queryWindowNr - j >= -nrOfWindowsOverlap +1"
				and the max at the end*/

					
				/*if(queryWindowNr == 0){

					ptr_Scores[hit].queryCoverage = 1;

				}
				else{*/


				j=1;
				while(queryWindowNr - j >= -nrOfWindowsOverlap +1  && j < nrOfWindowsOverlap && ptr_Scores[hit].windowNr[max(queryWindowNr-j,0)] == -1){j +=1;} 
						
				ptr_Scores[hit].queryCoverage += j*((float) stepSize/windowLength);

				//}

				ptr_Scores[hit].score += ptr_Scores[hit].queryCoverage*qLengthModdedDistScore;

				ptr_Scores[hit].queryWindowNr[queryWindowNr] = queryWindowNr;
				ptr_Scores[hit].windowNr[queryWindowNr] = ptr_matchCandidates[i].window.windowNr;
				//ptr_Scores[hit].windowNr[queryWindowNr] = ptr_matchCandidates[i].windowNr;
			}

		}
		else{ 

			if(strcmp(ptr_matchCandidates[i].structureName, "NN") != 0){
				printf("Warning: this structure name: %s was not found in the scores list (ptr_Score)\n", ptr_matchCandidates[i].structureName);
			}

		}
	}

}


int getHitStatistics(char * scoresFileName, int useSCOP_b, int useCATH_b, int topN){

	int returnVal = 1;

	FILE *ptr_fileScores;

	/*for reading DB files:*/
	int nrOfQueries = 0;
	int lineCnt;
	int lineCntThisQuery;
	char line[10000];
	char *linePart;
	int idxInLine = -1;
	int idx_dbName = -1;
	int idx_queryFileName = -1;
	int idx_qStructureId = -1;
	int idx_qClassId = -1;
	int idx_qChainId = -1;
	int idx_qChainLength = -1;
	int idx_qNrOfWindows = -1;

	int idx_dbFileName = -1;
	int idx_dbStructureId = -1;
	int idx_dbClassId = -1;
	int idx_dbChainId = -1;
	int idx_dbChainLength = -1;
	int idx_dbNrOfWindows = -1;

	int idx_rawScore = -1;
	int idx_score = -1;

	int score;
	int hiScore;
	int rawScore;
	int nrOfHits = 0;
	int tiedHit_b = 0;
	int nrOfTiedHits = 0;
	int nrOfHitsNot_ijk = 0; /*disregard the SCOP "non classes", i j and k */

	int not_self_match = 0; // to be set to 1 if the query is not itself the first (and best) match); else 0 
	int nrOfNotSelfMatch =0; //no of occurences with not_self_match = 1

	char *qStructureId;
	char *qStructureIdLast;
	char *dbStructureId;

	struct classInfo qClassInfo;
	struct classInfo dbClassInfo;

	qClassInfo.classId = (char *) malloc(classifierSize*sizeof(char));
	qClassInfo.classIdConv = (int *) malloc(classifierSize*sizeof(int));
	
	dbClassInfo.classId = (char *) malloc(classifierSize*sizeof(char));
	dbClassInfo.classIdConv = (int *) malloc(classifierSize*sizeof(int));

	/*init*/
	if(useSCOP_b ==1){

		qStructureId = (char *) malloc(10*sizeof(char));
		qStructureIdLast = (char *) malloc(10*sizeof(char));
	    dbStructureId = (char *) malloc(10*sizeof(char));
		strcpy(qStructureId,"NN");
		strcpy(qStructureIdLast,"NN");
		strcpy(dbStructureId, "NN");

		strcpy(qClassInfo.classId,"NN");
		qClassInfo.classIdConv[0] = -1;
		qClassInfo.classIdConv[1] = -1;
		qClassInfo.classIdConv[2] = -1;
		qClassInfo.classIdConv[3] = -1;

		strcpy(dbClassInfo.classId,"NN");
		dbClassInfo.classIdConv[0] = -1;
		dbClassInfo.classIdConv[1] = -1;
		dbClassInfo.classIdConv[2] = -1;
		dbClassInfo.classIdConv[3] = -1;
	}

	if(useCATH_b ==1){

		qStructureId = (char *) malloc(lengthCLFformat20*sizeof(char));
		qStructureIdLast = (char *) malloc(lengthCLFformat20*sizeof(char));
	    dbStructureId = (char *) malloc(lengthCLFformat20*sizeof(char));
		strcpy(qStructureId,"NN");
		strcpy(qStructureIdLast,"NN");
		strcpy(dbStructureId, "NN");

		strcpy(qClassInfo.classId,"NN");
		qClassInfo.classIdConv[0] = -1;
		qClassInfo.classIdConv[1] = -1;
		qClassInfo.classIdConv[2] = -1;
		qClassInfo.classIdConv[3] = -1;

		strcpy(dbClassInfo.classId,"NN");
		dbClassInfo.classIdConv[0] = -1;
		dbClassInfo.classIdConv[1] = -1;
		dbClassInfo.classIdConv[2] = -1;
		dbClassInfo.classIdConv[3] = -1;
	}

	/*Open the file containing the GIs on windows/pairs of the data base structures:*/
	ptr_fileScores = fopen(scoresFileName, "r");


	if(ptr_fileScores == NULL){
		printf("Sorry -- cannot find the file: %s\n", scoresFileName);
	}


	/*Read first line of file to record its structure*/
	lineCnt = 0;
	while (fgets(line,1000, ptr_fileScores)!=NULL  && lineCnt < 1){

		idxInLine = 0;

		/*loop trough the line and record important indices:*/
		linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
		while(linePart)
		{
			if (strcmp(linePart, "DB name")==0){					
				idx_dbName = idxInLine;
			}

			if (strcmp(linePart, "Query file")==0){					
				idx_queryFileName = idxInLine;
			}

			if (strcmp(linePart, "qStructureId")==0){					
				idx_qStructureId = idxInLine;
			}

			if (strcmp(linePart, "qClassId")==0){					
				idx_qClassId = idxInLine;
			}

			if (strcmp(linePart, "qChainId")==0){					
				idx_qChainId = idxInLine;
			}

			if (strcmp(linePart, "qChainLength")==0){					
				idx_qChainLength = idxInLine;
			}

			if (strcmp(linePart, "qNrOfWindows")==0){					
				idx_qNrOfWindows = idxInLine;
			}

			if (strcmp(linePart, "DB file")==0){					
				idx_dbFileName = idxInLine;
			}

			if (strcmp(linePart, "dbStructureId")==0){					
				idx_dbStructureId = idxInLine;
			}

			if (strcmp(linePart, "dbChainId")==0){					
				idx_dbChainId = idxInLine;
			}

			if (strcmp(linePart, "dbClassId")==0){					
				idx_dbClassId = idxInLine;
			}

			if (strcmp(linePart, "dbChainLength")==0){					
				idx_dbChainLength = idxInLine;
			}

			if (strcmp(linePart, "dbNrOfWindows")==0){					
				idx_dbNrOfWindows = idxInLine;
			}

			if (strcmp(linePart, "rawScore")==0){					
				idx_rawScore = idxInLine;
			}

			if (strcmp(linePart, "score")==0){					
				idx_score = idxInLine;
			}

			linePart=strtok(NULL, ";"); //this gets the next token of the string (which is line)
			idxInLine += 1;
			
		}

		lineCnt +=1;

	}

	/*get back to first line of file*/
	rewind(ptr_fileScores);

	/*loop through the file and record the hit statistics*/
	lineCnt = 0;
	lineCntThisQuery = 0;
	while (fgets(line,1000, ptr_fileScores)!=NULL){

		/*skip first line (is header)*/
		if(lineCnt == 0){
			lineCnt +=1; 
			continue;
		}

		//printf("lineCntThisQuery: %d\n", lineCntThisQuery);
		//printf("qStructureId: %s,  qStructureIdLast: %s\n", qStructureId, qStructureIdLast);
		//getchar();

		if(lineCntThisQuery >= 1){

			linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
			idxInLine = 0;
			while(linePart){

				if (idxInLine == idx_qStructureId){
					strcpy(qStructureId, linePart);
				}
				if (idxInLine == idx_qClassId){
					//strcpy(linePartFileName,linePart); //have to strcpy via a "fixed" string (char array)
					//printf("file %s", linePart);
					strcpy(qClassInfo.classId, linePart);
					/*printf("qClassInfo.classId:%s\n", qClassInfo.classId);
					getchar();*/
					if(useSCOP_b ==1){
						convertSCOP(&qClassInfo, 1);
					}
					if(useCATH_b ==1){
						convertCATH(&qClassInfo, 1);
					}
					/*printf("qClass: %s qClass conv: %d %d %d\n", qClassInfo.classId, qClassInfo.classIdConv[0], qClassInfo.classIdConv[1], qClassInfo.classIdConv[2]);
					getchar();*/
				}
				if (idxInLine == idx_dbStructureId){
					strcpy(dbStructureId, linePart);
				}
				if (idxInLine == idx_dbClassId){
					//strcpy(linePartFileName,linePart); //have to strcpy via a "fixed" string (char array)
					//printf("file %s", linePart);
					strcpy(dbClassInfo.classId, linePart);
					/*printf("dbClassInfo.classId:%s\n", dbClassInfo.classId);
					getchar();*/
					if(useSCOP_b ==1){
						convertSCOP(&dbClassInfo, 1);
					}
					if(useCATH_b ==1){
						convertCATH(&dbClassInfo, 1);
					}
					/*printf("dbClass: %s dbClass conv: %d %d %d\n", dbClassInfo.classId, dbClassInfo.classIdConv[0], dbClassInfo.classIdConv[1], dbClassInfo.classIdConv[2]);
					getchar();*/
				}

				if (idxInLine == idx_rawScore){
					//strcpy(linePartFileName,linePart); //have to strcpy via a "fixed" string (char array)
					//printf("file %s", linePart);
					rawScore = atoi(linePart);
				}

				linePart=strtok(NULL, ";"); //this gets the next token of the string (which is line)
				idxInLine += 1;
			}

			//printf("qClassId: %s, dbClassId: %s\n",qClassInfo.classId,dbClassInfo.classId );

			if(strcmp(qStructureId, qStructureIdLast)  == 0){
				
				if(not_self_match == 0){

					//if(useSCOP_b == 1){
				
					if(lineCntThisQuery ==1){
					
						if(qClassInfo.classIdConv[0] == dbClassInfo.classIdConv[0] && qClassInfo.classIdConv[1] == dbClassInfo.classIdConv[1] && qClassInfo.classIdConv[2] == dbClassInfo.classIdConv[2]){
						//if(qClassInfo.classIdConv[2] == dbClassInfo.classIdConv[2]){

							nrOfHits +=1;
					
						}

						hiScore = rawScore;

					}


					else{ //
			
						if(rawScore == hiScore){
				
							if(qClassInfo.classIdConv[0] == dbClassInfo.classIdConv[0] && qClassInfo.classIdConv[1] == dbClassInfo.classIdConv[1] && qClassInfo.classIdConv[2] == dbClassInfo.classIdConv[2]){
							//if(qClassInfo.classIdConv[2] == dbClassInfo.classIdConv[2]){

								tiedHit_b = 1;
							}
						}
					}
					//}
				}

				lineCntThisQuery +=1;

			}
			else{ //reset

				nrOfTiedHits += tiedHit_b;
				/*reset*/
				tiedHit_b = 0;
				lineCntThisQuery = 0;
				not_self_match = 0;
				strcpy(qStructureIdLast, qStructureId);

			}

		}

		if(lineCntThisQuery == 0){

			//printf("line cnt:%d\n",lineCnt);
			//printf("line: %s\n", line);
			linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
			idxInLine = 0;
			while(linePart){

				if (idxInLine == idx_qClassId){
					//strcpy(linePartFileName,linePart); //have to strcpy via a "fixed" string (char array)
					//printf("file %s", linePart);
					strcpy(qClassInfo.classId, linePart);
					/*printf("qClassInfo.classId:%s\n", qClassInfo.classId);
					getchar();*/
					if(useSCOP_b ==1){
						convertSCOP(&qClassInfo, 1);
					}
					if(useCATH_b ==1){
						convertCATH(&qClassInfo, 1);
					}
				}
				if (idxInLine == idx_qStructureId){
					strcpy(qStructureId, linePart);
					strcpy(qStructureIdLast, linePart);
				}
				if (idxInLine == idx_dbStructureId){
					strcpy(dbStructureId, linePart);
				}
				if (idxInLine == idx_dbClassId){
					//strcpy(linePartFileName,linePart); //have to strcpy via a "fixed" string (char array)
					//printf("file %s", linePart);
					strcpy(dbClassInfo.classId, linePart);
					if(useSCOP_b ==1){
						convertSCOP(&dbClassInfo, 1);
					}
					if(useCATH_b ==1){
						convertCATH(&dbClassInfo, 1);
					}
				}

				linePart=strtok(NULL, ";"); //this gets the next token of the string (which is line)
				idxInLine += 1;
			}
	
			//printf("qStructureId: %s,  dbStructureId: %s\n", qStructureId, dbStructureId);

			//check if query is itself first match; if not we check if the first match is a hit:
			if(strcmp(qStructureId,dbStructureId) !=0){
				printf("query: %s db:% s\n", qStructureId, dbStructureId);
				not_self_match = 1;
				nrOfNotSelfMatch +=1;
				if(qClassInfo.classIdConv[0] == dbClassInfo.classIdConv[0] && qClassInfo.classIdConv[1] == dbClassInfo.classIdConv[1] && qClassInfo.classIdConv[2] == dbClassInfo.classIdConv[2]){
				//if(qClassInfo.classIdConv[2] == dbClassInfo.classIdConv[2]){

						nrOfHits +=1;
					
				}
			}

			lineCntThisQuery +=1;

			nrOfQueries +=1;

		}

		lineCnt += 1;

	}

	printf("Nr of queries: %d, hits: %d, tied hits:%d, query itself not first match: %d \n",  nrOfQueries, nrOfHits, nrOfTiedHits, nrOfNotSelfMatch);

	getchar();

	return returnVal;

}

int getHitStatisticsCAT(char * scoresFileName, int useSCOP_b, int useCATH_b, int topN){

	int returnVal = 1;

	FILE *ptr_fileScores;

	/*for reading DB files:*/
	int nrOfQueries = 0;
	int lineCnt;
	int lineCntThisQuery;
	char line[10000];
	char *linePart;
	int idxInLine = -1;
	int idx_dbName = -1;
	int idx_queryFileName = -1;
	int idx_qStructureId = -1;
	int idx_qClassId = -1;
	int idx_qChainId = -1;
	int idx_qWindows = -1;
	int idx_qChainLength = -1;
	int idx_qNrOfWindows = -1;

	int idx_dbFileName = -1;
	int idx_dbStructureId = -1;
	int idx_dbClassId = -1;
	int idx_dbChainId = -1;
	int idx_dbWindows = -1;
	int idx_dbChainLength = -1;
	int idx_dbNrOfWindows = -1;

	int idx_rawScore = -1;
	int idx_score = -1;
	int idx_coverage = -1;

	int score;
	int hiScore;
	int rawScore;
	int tiedHit_b = 0;
	int nrOfTiedHits = 0;
	int nrOfHitsC = 0;
	int tiedHitC_b = 0;
	int nrOfHitsCA = 0;
	int nrOfHitsCAT = 0;
	int nrOfHitsNot_ijk = 0; /*disregard the SCOP "non classes", i j and k */

	int not_self_match = 0; // to be set to 1 if the query is not itself the first (and best) match); else 0 
	int nrOfNotSelfMatch =0; //no of occurences with not_self_match = 1

	char *qStructureId;
	char *qStructureIdLast;
	char *dbStructureId;

	struct classInfo qClassInfo;
	struct classInfo dbClassInfo;

	qClassInfo.classId = (char *) malloc(classifierSize*sizeof(char));
	qClassInfo.classIdConv = (int *) malloc(classifierSize*sizeof(int));
	
	dbClassInfo.classId = (char *) malloc(classifierSize*sizeof(char));
	dbClassInfo.classIdConv = (int *) malloc(classifierSize*sizeof(int));

	linePart = (char *) malloc(1000*sizeof(char));

	/*init*/
	if(useSCOP_b ==1){

		qStructureId = (char *) malloc(10*sizeof(char));
		qStructureIdLast = (char *) malloc(10*sizeof(char));
	    dbStructureId = (char *) malloc(10*sizeof(char));
		strcpy(qStructureId,"NN");
		strcpy(qStructureIdLast,"NN");
		strcpy(dbStructureId, "NN");

		strcpy(linePart, "bysvalen");

		strcpy(qClassInfo.classId,"NN");
		qClassInfo.classIdConv[0] = -1;
		qClassInfo.classIdConv[1] = -1;
		qClassInfo.classIdConv[2] = -1;
		qClassInfo.classIdConv[3] = -1;

		strcpy(dbClassInfo.classId,"NN");
		dbClassInfo.classIdConv[0] = -1;
		dbClassInfo.classIdConv[1] = -1;
		dbClassInfo.classIdConv[2] = -1;
		dbClassInfo.classIdConv[3] = -1;
	}

	if(useCATH_b ==1){

		qStructureId = (char *) malloc(lengthCLFformat20*sizeof(char));
		qStructureIdLast = (char *) malloc(lengthCLFformat20*sizeof(char));
	    dbStructureId = (char *) malloc(lengthCLFformat20*sizeof(char));
		strcpy(qStructureId,"NN");
		strcpy(qStructureIdLast,"NN");
		strcpy(dbStructureId, "NN");

		strcpy(qClassInfo.classId,"NN");
		qClassInfo.classIdConv[0] = -1;
		qClassInfo.classIdConv[1] = -1;
		qClassInfo.classIdConv[2] = -1;
		qClassInfo.classIdConv[3] = -1;

		strcpy(dbClassInfo.classId,"NN");
		dbClassInfo.classIdConv[0] = -1;
		dbClassInfo.classIdConv[1] = -1;
		dbClassInfo.classIdConv[2] = -1;
		dbClassInfo.classIdConv[3] = -1;
	}

	/*Open the file containing the GIs on windows/pairs of the data base structures:*/
	ptr_fileScores = fopen(scoresFileName, "r");


	if(ptr_fileScores == NULL){
		printf("Sorry -- cannot find the file: %s\n", scoresFileName);
	}


	/*Read first line of file to record its structure*/
	lineCnt = 0;
	while (fgets(line,10000, ptr_fileScores)!=NULL  && lineCnt < 1){


		//printf("line: %s\n", line);
	    //getchar();

		idxInLine = 0;
	

		/*loop trough the line and record important indices:*/
		linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
		while(linePart)
		{
			/*printf("linePart: %s\n", linePart);
			getchar();*/

			if (strcmp(linePart, "DB name")==0){					
				idx_dbName = idxInLine;
			}

			if (strcmp(linePart, "Query file")==0){					
				idx_queryFileName = idxInLine;
			}

			if (strcmp(linePart, "qStructureId")==0){					
				idx_qStructureId = idxInLine;
			}

			if (strcmp(linePart, "qClassId")==0){					
				idx_qClassId = idxInLine;
			}

			if (strcmp(linePart, "qChainId")==0){					
				idx_qChainId = idxInLine;
			}

			if (strcmp(linePart, "qChainLength")==0){					
				idx_qChainLength = idxInLine;
			}

			if (strcmp(linePart, "qNrOfWindows")==0){					
				idx_qNrOfWindows = idxInLine;
			}

			if (strcmp(linePart, "qWindows")==0){					
				idx_qWindows = idxInLine;
			}


			if (strcmp(linePart, "DB file")==0){					
				idx_dbFileName = idxInLine;
			}

			if (strcmp(linePart, "dbStructureId")==0){					
				idx_dbStructureId = idxInLine;
			}

			if (strcmp(linePart, "dbChainId")==0){					
				idx_dbChainId = idxInLine;
			}

			if (strcmp(linePart, "dbClassId")==0){					
				idx_dbClassId = idxInLine;
			}

			if (strcmp(linePart, "dbWindows")==0){					
				idx_dbWindows = idxInLine;
			}

			if (strcmp(linePart, "dbChainLength")==0){					
				idx_dbChainLength = idxInLine;
			}

			if (strcmp(linePart, "dbNrOfWindows")==0){					
				idx_dbNrOfWindows = idxInLine;
			}

			if (strcmp(linePart, "rawScore")==0){					
				idx_rawScore = idxInLine;
			}

			if (strcmp(linePart, "score")==0){					
				idx_score = idxInLine;
			}

			if (strcmp(linePart, "coverage")==0){					
				idx_coverage = idxInLine;
			}

			linePart=strtok(NULL, ";"); //this gets the next token of the string (which is line)
			idxInLine += 1;
			
		}

		lineCnt +=1;

	}

	/*get back to first line of file*/
	rewind(ptr_fileScores);

	/*loop through the file and record the hit statistics*/
	lineCnt = 0;
	lineCntThisQuery = 0;
	while (fgets(line,10000, ptr_fileScores)!=NULL){

		/*skip first line (is header)*/
		if(lineCnt == 0){
			lineCnt +=1; 
			continue;
		}

		//printf("lineCntThisQuery: %d\n", lineCntThisQuery);
		//printf("qStructureId: %s,  qStructureIdLast: %s\n", qStructureId, qStructureIdLast);
		//getchar();

		if(lineCntThisQuery >= 1){

			linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
			idxInLine = 0;
			while(linePart){

				/*printf("lineCntThisQuery >=1, linePart: %s\n", linePart);
				getchar();*/

				if (idxInLine == idx_qStructureId){
					strcpy(qStructureId, linePart);
				}
				if (idxInLine == idx_qClassId){
					//strcpy(linePartFileName,linePart); //have to strcpy via a "fixed" string (char array)
					//printf("file %s", linePart);
					strcpy(qClassInfo.classId, linePart);
					/*printf("qClassInfo.classId:%s\n", qClassInfo.classId);
					getchar();*/
					if(useSCOP_b ==1){
						convertSCOP(&qClassInfo, 1);
					}
					if(useCATH_b ==1){
						convertCATH(&qClassInfo, 1);
					}
					/*printf("qClass: %s qClass conv: %d %d %d\n", qClassInfo.classId, qClassInfo.classIdConv[0], qClassInfo.classIdConv[1], qClassInfo.classIdConv[2]);
					getchar();*/
				}
				if (idxInLine == idx_dbStructureId){
					strcpy(dbStructureId, linePart);
				}
				if (idxInLine == idx_dbClassId){
					//strcpy(linePartFileName,linePart); //have to strcpy via a "fixed" string (char array)
					//printf("file %s", linePart);
					strcpy(dbClassInfo.classId, linePart);
					/*printf("dbClassInfo.classId:%s\n", dbClassInfo.classId);
					getchar();*/
					if(useSCOP_b ==1){
						convertSCOP(&dbClassInfo, 1);
					}
					if(useCATH_b ==1){
						convertCATH(&dbClassInfo, 1);
					}
					/*printf("dbClass: %s dbClass conv: %d %d %d\n", dbClassInfo.classId, dbClassInfo.classIdConv[0], dbClassInfo.classIdConv[1], dbClassInfo.classIdConv[2]);
					getchar();*/
				}

				if (idxInLine == idx_rawScore){
					//strcpy(linePartFileName,linePart); //have to strcpy via a "fixed" string (char array)
					//printf("file %s", linePart);
					rawScore = atoi(linePart);
				}

				linePart=strtok(NULL, ";"); //this gets the next token of the string (which is line)
				idxInLine += 1;
			}

			//printf("qClassId: %s, dbClassId: %s\n",qClassInfo.classId,dbClassInfo.classId );

			if(strcmp(qStructureId, qStructureIdLast)  == 0){
				
				if(not_self_match == 0){

					//if(useSCOP_b == 1){
				
					if(lineCntThisQuery ==1){
					
						if(qClassInfo.classIdConv[0] == dbClassInfo.classIdConv[0]){
							
							nrOfHitsC +=1;

							if(qClassInfo.classIdConv[1] == dbClassInfo.classIdConv[1]){

								nrOfHitsCA +=1;

								if(qClassInfo.classIdConv[2] == dbClassInfo.classIdConv[2]){

									nrOfHitsCAT +=1;

								}
								
							}
						//if(qClassInfo.classIdConv[2] == dbClassInfo.classIdConv[2]){
					
						}

						hiScore = rawScore;

					}


					else{ //
			
						if(rawScore == hiScore){
				
							if(qClassInfo.classIdConv[0] == dbClassInfo.classIdConv[0]){
							
							nrOfHitsC +=1;

							if(qClassInfo.classIdConv[1] == dbClassInfo.classIdConv[1]){

								nrOfHitsCA +=1;

								if(qClassInfo.classIdConv[2] == dbClassInfo.classIdConv[2]){

									nrOfHitsCAT +=1;

								}
								
							}
						//if(qClassInfo.classIdConv[2] == dbClassInfo.classIdConv[2]){
					
						}
						}
					}
					//}
				}

				lineCntThisQuery +=1;

			}
			else{ //reset

				nrOfTiedHits += tiedHit_b;
				/*reset*/
				tiedHit_b = 0;
				lineCntThisQuery = 0;
				not_self_match = 0;
				strcpy(qStructureIdLast, qStructureId);

			}

		}

		if(lineCntThisQuery == 0){

			//printf("line cnt:%d\n",lineCnt);
			//printf("line: %s\n", line);
			linePart = strtok(line, ";"); //";" is used as the delimiter throughout 
			idxInLine = 0;
			while(linePart){

				/*printf("lineCntThisQuery = 0, idxInLine %d linePart: %s\n", idxInLine, linePart);
				getchar();*/

				if (idxInLine == idx_qClassId){
					//strcpy(linePartFileName,linePart); //have to strcpy via a "fixed" string (char array)
					//printf("file %s", linePart);
					strcpy(qClassInfo.classId, linePart);
					/*printf("qClassInfo.classId:%s\n", qClassInfo.classId);
					getchar();*/
					if(useSCOP_b ==1){
						convertSCOP(&qClassInfo, 1);
					}
					if(useCATH_b ==1){
						convertCATH(&qClassInfo, 1);
					}
				}
				if (idxInLine == idx_qStructureId){
					strcpy(qStructureId, linePart);
					strcpy(qStructureIdLast, linePart);
				}
				if (idxInLine == idx_dbStructureId){
					strcpy(dbStructureId, linePart);
				}
				if (idxInLine == idx_dbClassId){
					//strcpy(linePartFileName,linePart); //have to strcpy via a "fixed" string (char array)
					//printf("file %s", linePart);
					strcpy(dbClassInfo.classId, linePart);
					if(useSCOP_b ==1){
						convertSCOP(&dbClassInfo, 1);
					}
					if(useCATH_b ==1){
						convertCATH(&dbClassInfo, 1);
					}
				}

				linePart=strtok(NULL, ";"); //this gets the next token of the string (which is line)
				idxInLine += 1;
			}
	
			//printf("qStructureId: %s,  dbStructureId: %s\n", qStructureId, dbStructureId);

			//check if query is itself first match; if not we check if the first match is a hit:
			if(strcmp(qStructureId,dbStructureId) !=0){
				printf("query: %s db:% s\n", qStructureId, dbStructureId);
				//getchar();
				not_self_match = 1;
				nrOfNotSelfMatch +=1;
				if(qClassInfo.classIdConv[0] == dbClassInfo.classIdConv[0]){
							
					nrOfHitsC +=1;

					if(qClassInfo.classIdConv[1] == dbClassInfo.classIdConv[1]){

						nrOfHitsCA +=1;

						if(qClassInfo.classIdConv[2] == dbClassInfo.classIdConv[2]){

							nrOfHitsCAT +=1;

						}
								
					}
					
				}

			}

			lineCntThisQuery +=1;

			nrOfQueries +=1;

		}

		lineCnt += 1;

	}

	printf("Nr of queries: %d, hits C A T: %d %d %d, tied hits:%d, query itself not first match: %d \n",  nrOfQueries, nrOfHitsC, nrOfHitsCA, nrOfHitsCAT, nrOfTiedHits, nrOfNotSelfMatch);

	getchar();

	return returnVal;

}


double collectMatchCandidates2(struct matchCandidate *ptr_matchCandidate, int topN, int matchCandidateIdx, int lMin, double distMin, struct I_windows_ptr db_I_windows){

	/*pseudo code:
	1) if matchCandidateIdx < topN: collect
	2) if matchCandidateIdx = topN: collect and sort
	3) else: find the right index/position of the match in the ordered array of collected matches; 
	shift all indices below by 1; remove the last entry which are now not among the topN best hits 
	(this will have index topN+1) 
	*/

	double returnVal = 1e20; // some gigantic number

	if (matchCandidateIdx < topN){

		ptr_matchCandidate[matchCandidateIdx].fileName = db_I_windows.fileName;
		ptr_matchCandidate[matchCandidateIdx].structureName = db_I_windows.structureName;
		ptr_matchCandidate[matchCandidateIdx].classId = db_I_windows.classId;
		ptr_matchCandidate[matchCandidateIdx].chainId = db_I_windows.chainId;
		ptr_matchCandidate[matchCandidateIdx].chainNr = db_I_windows.chainNr;
		ptr_matchCandidate[matchCandidateIdx].chainLen = db_I_windows.chainLen;
		ptr_matchCandidate[matchCandidateIdx].nrOfWindows = db_I_windows.nrOfWindows;
		//ptr_matchCandidate[matchCandidateIdx].windowsNr = db_I_windows.windowsNr[lMin];
		//ptr_matchCandidate[matchCandidateIdx].segIndices = db_I_windows.segIndices[lMin];
		ptr_matchCandidate[matchCandidateIdx].window.windowNr = db_I_windows.window[lMin].windowNr;
		ptr_matchCandidate[matchCandidateIdx].window.segIndices = db_I_windows.window[lMin].segIndices;
		ptr_matchCandidate[matchCandidateIdx].distMin = distMin;

		returnVal = 1e20; // some gigantic number

	}

	return returnVal;

}

double collectMatchCandidates(struct matchCandidate *ptr_matchCandidate, int topN, int *ptr_matchCandidateIdx, int lMin, double distMin, int fileNr, struct db_I_windows_order2_ptr *dbResults){

	/*pseudo code:
	1) if *ptr_matchCandidateIdx < topN: collect
	2) if *ptr_matchCandidateIdx = topN: collect and sort
	3) else: find the right index/position of the match in the ordered array of collected matches; 
	shift all indices below by 1; remove the last entry which are now not among the topN best hits 
	(this will have index topN+1) 
	*/

	double returnVal = 1e20; // some gigantic number

	if (*ptr_matchCandidateIdx < topN-1){

		//getchar();

		//printf("In collect, db file: %s\n", dbResults -> fileName);
		//getchar();

		//ptr_matchCandidate[*ptr_matchCandidateIdx].fileName = strdup(dbResults.fileName);
		//ptr_matchCandidate[*ptr_matchCandidateIdx].fileName = dbResults -> fileName;
		strcpy(ptr_matchCandidate[*ptr_matchCandidateIdx].fileName, dbResults -> fileName);
		//ptr_matchCandidate[*ptr_matchCandidateIdx].structureName = dbResults -> structureName;
		strcpy(ptr_matchCandidate[*ptr_matchCandidateIdx].structureName, dbResults -> structureName);
		ptr_matchCandidate[*ptr_matchCandidateIdx].fileNr = fileNr;
		strcpy(ptr_matchCandidate[*ptr_matchCandidateIdx].classId, dbResults -> classId);
		strcpy(ptr_matchCandidate[*ptr_matchCandidateIdx].chainId, dbResults -> chainId);
		ptr_matchCandidate[*ptr_matchCandidateIdx].chainNr = dbResults -> chainNr;
		ptr_matchCandidate[*ptr_matchCandidateIdx].chainLen = dbResults -> chainLen;
		ptr_matchCandidate[*ptr_matchCandidateIdx].nrOfWindows = dbResults -> nrOfWindows;
		/*ptr_matchCandidate[*ptr_matchCandidateIdx].windowsNr = dbResults -> windowsNr[lMin];
		ptr_matchCandidate[*ptr_matchCandidateIdx].segIndices[0] = dbResults -> segIndices[lMin][0];
		ptr_matchCandidate[*ptr_matchCandidateIdx].segIndices[1] = dbResults -> segIndices[lMin][1];*/

		ptr_matchCandidate[*ptr_matchCandidateIdx].window.windowNr = dbResults -> window[lMin].windowNr;
		ptr_matchCandidate[*ptr_matchCandidateIdx].window.segIndices[0] = dbResults -> window[lMin].segIndices[0];
		ptr_matchCandidate[*ptr_matchCandidateIdx].window.segIndices[1] = dbResults -> window[lMin].segIndices[1];

		ptr_matchCandidate[*ptr_matchCandidateIdx].distMin = distMin;
		ptr_matchCandidate[*ptr_matchCandidateIdx].score = 0.0;
		//printf("P1: match cand file name: %s\n", dbResults -> fileName);

		//sort (we do this here only to serve the case in which there are less than topN match candidates
		//bubbleSort2(ptr_matchCandidate, topN);
		//heapSortMatchCandidates(ptr_matchCandidate, topN);

		*ptr_matchCandidateIdx += 1;
		returnVal = returnVal; // some gigantic number

		//printf("match cand file name: %s\n", dbResults -> fileName);

		
	}
	else{

		if (*ptr_matchCandidateIdx == topN-1){

			//collect
			//ptr_matchCandidate[*ptr_matchCandidateIdx].fileName = dbResults -> fileName;
			strcpy(ptr_matchCandidate[*ptr_matchCandidateIdx].fileName, dbResults -> fileName);
			//ptr_matchCandidate[*ptr_matchCandidateIdx].structureName = dbResults -> structureName;
			strcpy(ptr_matchCandidate[*ptr_matchCandidateIdx].structureName, dbResults -> structureName);
			ptr_matchCandidate[*ptr_matchCandidateIdx].fileNr = fileNr;
			strcpy(ptr_matchCandidate[*ptr_matchCandidateIdx].classId, dbResults -> classId);
			strcpy(ptr_matchCandidate[*ptr_matchCandidateIdx].chainId, dbResults -> chainId);
			ptr_matchCandidate[*ptr_matchCandidateIdx].chainNr = dbResults -> chainNr;
			ptr_matchCandidate[*ptr_matchCandidateIdx].chainLen = dbResults -> chainLen;
			ptr_matchCandidate[*ptr_matchCandidateIdx].nrOfWindows = dbResults -> nrOfWindows;
			//ptr_matchCandidate[*ptr_matchCandidateIdx].windowsNr = dbResults -> windowsNr[lMin];
			//ptr_matchCandidate[*ptr_matchCandidateIdx].segIndices[0] = dbResults -> segIndices[lMin][0];
			//ptr_matchCandidate[*ptr_matchCandidateIdx].segIndices[1] = dbResults -> segIndices[lMin][1];

			ptr_matchCandidate[*ptr_matchCandidateIdx].window.windowNr = dbResults -> window[lMin].windowNr;
			ptr_matchCandidate[*ptr_matchCandidateIdx].window.segIndices[0] = dbResults -> window[lMin].segIndices[0];
			ptr_matchCandidate[*ptr_matchCandidateIdx].window.segIndices[1] = dbResults -> window[lMin].segIndices[1];

			ptr_matchCandidate[*ptr_matchCandidateIdx].distMin = distMin;
			ptr_matchCandidate[*ptr_matchCandidateIdx].score = 0.0;

			/*printf("P1: match cand file name: %s\n", dbResults -> fileName);*/
			//sort
			bubbleSort2(ptr_matchCandidate, topN);
			//heapSortMatchCandidates(ptr_matchCandidate, topN);
			//highest distMin:
			returnVal = ptr_matchCandidate[topN-1].distMin;

			//printf("highest distMin:%f\n", returnVal);
			//getchar();
			//printf("P2: match cand file name: %s\n", dbResults -> fileName);

		}

	}

	
	return returnVal;

}

double collectMatchCandidates_windowPairs(struct matchCandidate_windowPair *ptr_matchCandidate_windowPairs, int topN, int *ptr_matchCandidateIdx, int iMin, int jMin, double distMin,  int fileNr, struct db_I_windowPairs_order2_ptr *dbResults){

	/*pseudo code:
	1) if *ptr_matchCandidateIdx < topN: collect
	2) if *ptr_matchCandidateIdx = topN: collect and sort
	3) else: find the right index/position of the match in the ordered array of collected matches; 
	shift all indices below by 1; remove the last entry which are now not among the topN best hits 
	(this will have index topN+1) 
	*/

	double returnVal = 1e20; // some gigantic number

	if (*ptr_matchCandidateIdx < topN-1){

		//getchar();

		//printf("In collect, db file: %s\n", dbResults -> fileName);
		//getchar();

		//ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].fileName = strdup(dbResults.fileName);
		//ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].fileName = dbResults -> fileName;
		strcpy(ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].fileName, dbResults -> fileName);
		//ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].structureName = dbResults -> structureName;
		strcpy(ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].structureName, dbResults -> structureName);
		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].fileNr = fileNr;
		strcpy(ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].classId, dbResults -> classId);
		strcpy(ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].chainId, dbResults -> chainId);
		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].chainNr = dbResults -> chainNr;
		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].chainLen = dbResults -> chainLen;
		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].nrOfWindows = dbResults -> nrOfWindows;

		/*ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowsNr_1 = dbResults -> windowsNr_1[iMin];
		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].segIndices_1[0] = dbResults -> segIndices_1[iMin][0];
		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].segIndices_1[1] = dbResults -> segIndices_1[iMin][1];
		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowsNr_2 = dbResults -> windowsNr_2[jMin];
		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].segIndices_2[0] = dbResults -> segIndices_2[jMin][0];
		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].segIndices_2[1] = dbResults -> segIndices_2[jMin][1];*/

		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowPair.windowNr_1 = dbResults -> windowPair[iMin][jMin].windowNr_1;
		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowPair.segIndices_1[0] = dbResults -> windowPair[iMin][jMin].segIndices_1[0];
		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowPair.segIndices_1[1] = dbResults -> windowPair[iMin][jMin].segIndices_1[1];
		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowPair.windowNr_2 = dbResults -> windowPair[iMin][jMin].windowNr_2;
		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowPair.segIndices_2[0] = dbResults -> windowPair[iMin][jMin].segIndices_2[0];
		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowPair.segIndices_2[1] = dbResults -> windowPair[iMin][jMin].segIndices_2[1];

		ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].distMin = distMin;

		//printf("P1: match cand file name: %s chain len: %d\n", dbResults -> fileName, dbResults -> chainLen);

		*ptr_matchCandidateIdx += 1;
		returnVal = 1e20; // some gigantic number

		//printf("match cand file name: %s\n", dbResults -> fileName);

		//sort (we do this here only to serve the case in which there are less than topN match candidates
		bubbleSort2_windowPairs(ptr_matchCandidate_windowPairs, topN);

	}
	else{

		if (*ptr_matchCandidateIdx == topN-1){

			//collect
			//ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].fileName = dbResults -> fileName;
			strcpy(ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].fileName, dbResults -> fileName);
			//ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].structureName = dbResults -> structureName;
			strcpy(ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].structureName, dbResults -> structureName);
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].fileNr = fileNr;
			strcpy(ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].classId, dbResults -> classId);
			strcpy(ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].chainId, dbResults -> chainId);
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].chainNr = dbResults -> chainNr;
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].chainLen = dbResults -> chainLen;
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].nrOfWindows = dbResults -> nrOfWindows;
			
			/*ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowsNr_1 = dbResults -> windowsNr_1[iMin];
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].segIndices_1[0] = dbResults -> segIndices_1[iMin][0];
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].segIndices_1[1] = dbResults -> segIndices_1[iMin][1];
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowsNr_2 = dbResults -> windowsNr_2[jMin];
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].segIndices_2[0] = dbResults -> segIndices_2[jMin][0];
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].segIndices_2[1] = dbResults -> segIndices_2[jMin][1];*/
			
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowPair.windowNr_1 = dbResults -> windowPair[iMin][jMin].windowNr_1;
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowPair.segIndices_1[0] = dbResults -> windowPair[iMin][jMin].segIndices_1[0];
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowPair.segIndices_1[1] = dbResults -> windowPair[iMin][jMin].segIndices_1[1];
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowPair.windowNr_2 = dbResults -> windowPair[iMin][jMin].windowNr_2;
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowPair.segIndices_2[0] = dbResults -> windowPair[iMin][jMin].segIndices_2[0];
			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].windowPair.segIndices_2[1] = dbResults -> windowPair[iMin][jMin].segIndices_2[1];


			ptr_matchCandidate_windowPairs[*ptr_matchCandidateIdx].distMin = distMin;

			//printf("P1: match cand file name: %s\n", dbResults -> fileName);
			//sort
			bubbleSort2_windowPairs(ptr_matchCandidate_windowPairs, topN);
			//highest distMin:
			returnVal = ptr_matchCandidate_windowPairs[topN-1].distMin;

			//printf("highest distMin:%f\n", returnVal);
			//getchar();
			//printf("P2: match cand file name: %s\n", dbResults -> fileName);

		}

	}


	return returnVal;

}


int alloc_init_matchCandidates(struct matchCandidate **ptr_ptr_matchCandidates, int topN){

	int i = 0;
	int j = 0;
	struct matchCandidate *ptr_matchCandidates;

	/*allocate memory to pointer ptr_I_values and initialize:*/
	ptr_matchCandidates = (struct matchCandidate *) malloc (topN*sizeof(struct matchCandidate));
	for (i=0; i< topN; i++){

		ptr_matchCandidates[i].fileName = (char *) malloc (sizeof(char)*1000);
		ptr_matchCandidates[i].structureName = (char *) malloc (sizeof(char)*20);
		ptr_matchCandidates[i].classId = (char *) malloc (sizeof(char)*20);
		ptr_matchCandidates[i].chainId = (char *) malloc (sizeof(char)*2);

		//ptr_matchCandidates[i].segIndices = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 
		ptr_matchCandidates[i].window.segIndices = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/

		strcpy(ptr_matchCandidates[i].fileName, "NN");
		strcpy(ptr_matchCandidates[i].structureName, "NN");
		ptr_matchCandidates[i].fileNr = -17;
		strcpy(ptr_matchCandidates[i].classId,"NN");
		strcpy(ptr_matchCandidates[i].chainId,"NN");
		ptr_matchCandidates[i].chainNr = -1171;
		ptr_matchCandidates[i].chainLen = -1171;
		ptr_matchCandidates[i].nrOfWindows = -1171;
		//ptr_matchCandidates[i].windowsNr = -1171;
		//ptr_matchCandidates[i].segIndices[0] = -1171;
		//ptr_matchCandidates[i].segIndices[1] = -1171;
		ptr_matchCandidates[i].window.windowNr = -1171;
		ptr_matchCandidates[i].window.segIndices[0] = -1171;
		ptr_matchCandidates[i].window.segIndices[1] = -1171;
		ptr_matchCandidates[i].distMin = 1e20;
		ptr_matchCandidates[i].score = 0;
	}


	*ptr_ptr_matchCandidates = ptr_matchCandidates;


	return 1;
}

int reinit_matchCandidates(struct matchCandidate *ptr_matchCandidates, int topN){

	int i = 0;
	int j = 0;

	for (i=0; i< topN; i++){

		strcpy(ptr_matchCandidates[i].fileName, "NN");
		strcpy(ptr_matchCandidates[i].structureName, "NN");
		ptr_matchCandidates[i].fileNr = -17;
		strcpy(ptr_matchCandidates[i].classId,"NN");
		strcpy(ptr_matchCandidates[i].chainId,"NN");
		ptr_matchCandidates[i].chainNr = -1172;
		ptr_matchCandidates[i].chainLen = -1172;
		ptr_matchCandidates[i].nrOfWindows = -1172;
		//ptr_matchCandidates[i].windowsNr = -1172;
		//ptr_matchCandidates[i].segIndices[0] = -1172;
		//ptr_matchCandidates[i].segIndices[1] = -1172;
		ptr_matchCandidates[i].window.windowNr = -1172;
		ptr_matchCandidates[i].window.segIndices[0] = -1172;
		ptr_matchCandidates[i].window.segIndices[1] = -1172;
		ptr_matchCandidates[i].distMin = 1e20;
		ptr_matchCandidates[i].score = 0;
	}

	return 1;
}


int alloc_init_ptr_matchCandidates(struct matchCandidate ***ptr_ptr_ptr_matchCandidates, int nrOfWindows, int topN){

	int returnVal = 0;
	
	int k = 0;
	struct matchCandidate **ptr_ptr_matchCandidates;

	/*allocate memory to pointer ptr_I_values and initialize:*/
	ptr_ptr_matchCandidates = (struct matchCandidate **) malloc (nrOfWindows*sizeof(struct matchCandidate *));
	
	for (k=0; k < nrOfWindows; k++){

		returnVal = alloc_init_matchCandidates(&ptr_ptr_matchCandidates[k], topN); 

	}


	*ptr_ptr_ptr_matchCandidates = ptr_ptr_matchCandidates;

	return returnVal;
}

int reinit_ptr_matchCandidates(struct matchCandidate **ptr_ptr_matchCandidates, int nrOfWindows, int topN){

	int returnVal = 0;
	
	int k = 0;

	for (k=0; k < topN; k++){

		reinit_matchCandidates(ptr_ptr_matchCandidates[k], topN);

	}

	return returnVal;
}



int alloc_init_matchCandidates_windowPairs(struct matchCandidate_windowPair **ptr_ptr_matchCandidates_windowPairs, int topN){

	int i = 0;
	int j = 0;
	struct matchCandidate_windowPair *ptr_matchCandidates_windowPairs;

	/*allocate memory to pointer ptr_I_values and initialize:*/
	ptr_matchCandidates_windowPairs = (struct matchCandidate_windowPair *) malloc (topN*sizeof(struct matchCandidate_windowPair));
	for (i=0; i< topN; i++){

		ptr_matchCandidates_windowPairs[i].fileName = (char *) malloc (sizeof(char)*1000);
		ptr_matchCandidates_windowPairs[i].structureName = (char *) malloc (sizeof(char)*20);
		ptr_matchCandidates_windowPairs[i].classId = (char *) malloc (sizeof(char)*20);
		ptr_matchCandidates_windowPairs[i].chainId = (char *) malloc (sizeof(char)*2);

		//ptr_matchCandidates_windowPairs[i].segIndices_1 = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 
		//ptr_matchCandidates_windowPairs[i].segIndices_2 = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

		ptr_matchCandidates_windowPairs[i].windowPair.segIndices_1 = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 
		ptr_matchCandidates_windowPairs[i].windowPair.segIndices_2 = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/

		strcpy(ptr_matchCandidates_windowPairs[i].fileName, "NN");
		strcpy(ptr_matchCandidates_windowPairs[i].structureName, "NN");
		strcpy(ptr_matchCandidates_windowPairs[i].classId,"NN");
		strcpy(ptr_matchCandidates_windowPairs[i].chainId,"NN");

		ptr_matchCandidates_windowPairs[i].chainNr = -1171;
		ptr_matchCandidates_windowPairs[i].chainLen = -1171;
		ptr_matchCandidates_windowPairs[i].nrOfWindows = -1171;

		/*ptr_matchCandidates_windowPairs[i].windowsNr_1 = -1171;
		ptr_matchCandidates_windowPairs[i].segIndices_1[0] = -1171;
		ptr_matchCandidates_windowPairs[i].segIndices_1[1] = -1171;
		ptr_matchCandidates_windowPairs[i].windowsNr_2 = -1171;
		ptr_matchCandidates_windowPairs[i].segIndices_2[0] = -1171;
		ptr_matchCandidates_windowPairs[i].segIndices_2[1] = -1171;*/

		ptr_matchCandidates_windowPairs[i].windowPair.windowNr_1 = -1171;
		ptr_matchCandidates_windowPairs[i].windowPair.segIndices_1[0] = -1171;
		ptr_matchCandidates_windowPairs[i].windowPair.segIndices_1[1] = -1171;
		ptr_matchCandidates_windowPairs[i].windowPair.windowNr_2 = -1171;
		ptr_matchCandidates_windowPairs[i].windowPair.segIndices_2[0] = -1171;
		ptr_matchCandidates_windowPairs[i].windowPair.segIndices_2[1] = -1171;

		ptr_matchCandidates_windowPairs[i].distMin = 1e20;
	}


	*ptr_ptr_matchCandidates_windowPairs = ptr_matchCandidates_windowPairs;


	return 1;
}

int reinit_matchCandidates_windowPairs(struct matchCandidate_windowPair *ptr_matchCandidates_windowPairs, int topN){

	int i = 0;
	int j = 0;

	for (i=0; i< topN; i++){

		strcpy(ptr_matchCandidates_windowPairs[i].fileName, "NN");
		strcpy(ptr_matchCandidates_windowPairs[i].structureName, "NN");
		strcpy(ptr_matchCandidates_windowPairs[i].classId,"NN");
		strcpy(ptr_matchCandidates_windowPairs[i].chainId,"NN");
		ptr_matchCandidates_windowPairs[i].chainNr = -1172;
		ptr_matchCandidates_windowPairs[i].chainLen = -1172;
		ptr_matchCandidates_windowPairs[i].nrOfWindows = -1172;

		/*ptr_matchCandidates_windowPairs[i].windowsNr_1 = -1172;
		ptr_matchCandidates_windowPairs[i].segIndices_1[0] = -1172;
		ptr_matchCandidates_windowPairs[i].segIndices_1[1] = -1172;
		ptr_matchCandidates_windowPairs[i].windowsNr_2 = -1172;
		ptr_matchCandidates_windowPairs[i].segIndices_2[0] = -1172;
		ptr_matchCandidates_windowPairs[i].segIndices_2[1] = -1172;*/

		ptr_matchCandidates_windowPairs[i].windowPair.windowNr_1 = -1172;
		ptr_matchCandidates_windowPairs[i].windowPair.segIndices_1[0] = -1172;
		ptr_matchCandidates_windowPairs[i].windowPair.segIndices_1[1] = -1172;
		ptr_matchCandidates_windowPairs[i].windowPair.windowNr_2 = -1172;
		ptr_matchCandidates_windowPairs[i].windowPair.segIndices_2[0] = -1172;
		ptr_matchCandidates_windowPairs[i].windowPair.segIndices_2[1] = -1172;

		ptr_matchCandidates_windowPairs[i].distMin = 1e20;
	}

	return 1;
}


int alloc_init_matchScores(struct matchScore **ptr_ptr_matchScores, struct db_I_windows_order2_ptr *ptr_dbResults,int nrOfDbStructures, int topN, int maxNrOfWindows){

	int i = 0;
	int j = 0;
	struct matchScore *ptr_matchScores;

	/*allocate memory to pointer ptr_I_values and initialize:*/
	ptr_matchScores = (struct matchScore *) malloc (nrOfDbStructures*sizeof(struct matchScore));
	for (i=0; i< nrOfDbStructures; i++){

		/*ptr_matchScores[i].queryFileName = (char *) malloc (sizeof(char)*1000);
		ptr_matchScores[i].queryChainId = (char *) malloc (sizeof(char));*/

		ptr_matchScores[i].matchFileName = (char *) malloc (sizeof(char)*1000);
		ptr_matchScores[i].matchStructureName = (char *) malloc (sizeof(char)*1000);
		ptr_matchScores[i].matchClassId = (char *) malloc (sizeof(char)*20);
		ptr_matchScores[i].matchChainId = (char *) malloc (sizeof(char));

		//ptr_matchScores[i].matchNrOfWindows = (int) malloc(sizeof(int));
		//ptr_matchScores[i].rawScore = (int) malloc(sizeof(int));
		//ptr_matchScores[i].score =  malloc(sizeof(double));


		ptr_matchScores[i].places = (int *) malloc(topN*sizeof(int)); /*just a 2-d int array*/ 

		ptr_matchScores[i].queryWindowNr = (int *) malloc(maxNrOfWindows*sizeof(int)); /*just a 2-d int array*/ 
		ptr_matchScores[i].windowNr = (int *) malloc(maxNrOfWindows*sizeof(int)); /*just a 2-d int array*/ 


		/*init*/
		strcpy(ptr_matchScores[i].matchFileName, ptr_dbResults[i].fileName);
		strcpy(ptr_matchScores[i].matchStructureName, ptr_dbResults[i].structureName);
		strcpy(ptr_matchScores[i].matchClassId, ptr_dbResults[i].classId);
		strcpy(ptr_matchScores[i].matchChainId, ptr_dbResults[i].chainId);

		ptr_matchScores[i].matchChainLen = ptr_dbResults[i].chainLen;
		ptr_matchScores[i].matchNrOfWindows = ptr_dbResults[i].nrOfWindows;

		ptr_matchScores[i].queryCoverage = 0;

		ptr_matchScores[i].countScore = 0.0;
		ptr_matchScores[i].aggrWinScore = 0.0;
		ptr_matchScores[i].score = 0.0;


		for(j = 0; j< topN;j++){

			ptr_matchScores[i].places[j] = 0;

		}

		for(j = 0; j< maxNrOfWindows;j++){

			ptr_matchScores[i].queryWindowNr[j] = -1;
			ptr_matchScores[i].windowNr[j] = -1;
		}

	}


	*ptr_ptr_matchScores = ptr_matchScores;


	return 1;
}

/*for re-initializing an existing ptr_matchScores*/
int reinit_matchScores(struct matchScore *ptr_matchScores, int nrOfDbStructures, int topN, int maxNrOfWindows){

	int i = 0;
	int j = 0;

	for (i=0; i< nrOfDbStructures; i++){

		/*re-init the fields containing scores*/
		ptr_matchScores[i].countScore = 0.0;
		ptr_matchScores[i].aggrWinScore = 0.0;
		ptr_matchScores[i].score = 0.0;
		ptr_matchScores[i].queryCoverage = 0;

		for(j = 0; j< topN;j++){

			ptr_matchScores[i].places[j] = 0;

		}

		for(j = 0; j< maxNrOfWindows;j++){

			ptr_matchScores[i].queryWindowNr[j] = -1;
			ptr_matchScores[i].windowNr[j] = -1;

		}

	}

	return 1;
}

/*ptr for holding the array of scores across a set of queries:*/
int alloc_init_queryRawScore(struct queryRawScore **ptr_ptr_queryRawScore, int topN, int alreadyAllocatedTo_topN, int windowInfoSize, int alreadyAllocatedTo_windowInfoSize){

	int i = 0;
	struct queryRawScore *ptr_queryRawScore;


	if(alreadyAllocatedTo_windowInfoSize == 0){ //indicates that no allocation was done yet
	
		/*allocate memory to pointer ptr_I_values and initialize:*/
		ptr_queryRawScore = (struct queryRawScore *) malloc (topN*sizeof(struct queryRawScore));

		for (i=0; i< topN; i++){

			ptr_queryRawScore[i].fileName = (char *) malloc (sizeof(char)*fileNameCharSize);
			ptr_queryRawScore[i].structureName = (char *) malloc (sizeof(char)*structureNameCharSize);
			ptr_queryRawScore[i].classId = (char *) malloc (sizeof(char)*classIdCharSize);
			ptr_queryRawScore[i].chainId = (char *) malloc (sizeof(char)*chainIdCharSize);
			ptr_queryRawScore[i].windowInfo = (char *) malloc (sizeof(char)*windowInfoSize);
		
			/*init*/
			strcpy(ptr_queryRawScore[i].fileName, "NN");
			strcpy(ptr_queryRawScore[i].structureName, "NN");
			strcpy(ptr_queryRawScore[i].classId,"NN");
			strcpy(ptr_queryRawScore[i].chainId,"NN");
			strcpy(ptr_queryRawScore[i].windowInfo,"-1");

			ptr_queryRawScore[i].chainNr = -1171;
			ptr_queryRawScore[i].chainLen = -1171;
			ptr_queryRawScore[i].nrOfWindows = -1171;

			ptr_queryRawScore[i].score = 0.0;
			ptr_queryRawScore[i].pValue = 1.0;

			ptr_queryRawScore[i].optionalValue = -9999;
		}
	}
	else if(windowInfoSize > alreadyAllocatedTo_windowInfoSize || topN > alreadyAllocatedTo_topN){ //indicates that reallocation must be done

		ptr_queryRawScore = *ptr_ptr_queryRawScore;

		//printf("(saaledes stadig  ..\n");

		printf("topN %d alreadyAllocatedTo_topN %d\n", topN, alreadyAllocatedTo_topN);
		//getchar();

		if(topN > alreadyAllocatedTo_topN){
			ptr_queryRawScore = realloc (ptr_queryRawScore, topN*sizeof(struct queryRawScore));
		}

		//printf("... ganske hemmelig)\n");

		for (i=0; i< topN; i++){

			if(i >= alreadyAllocatedTo_topN){	
				ptr_queryRawScore[i].fileName = (char *) malloc (sizeof(char)*fileNameCharSize);
				ptr_queryRawScore[i].structureName = (char *) malloc (sizeof(char)*structureNameCharSize);
				ptr_queryRawScore[i].classId = (char *) malloc (sizeof(char)*classIdCharSize);
				ptr_queryRawScore[i].chainId = (char *) malloc (sizeof(char)*chainIdCharSize);
				ptr_queryRawScore[i].windowInfo = (char *) malloc (sizeof(char)*alreadyAllocatedTo_windowInfoSize);	
			}
			
			//printf("Sommeren mild oedsles ...\n");

			if(windowInfoSize > alreadyAllocatedTo_windowInfoSize){
				ptr_queryRawScore[i].windowInfo = realloc (ptr_queryRawScore[i].windowInfo, sizeof(char)*windowInfoSize);
			}

			//printf(".. kyssene nemmere pluds'ligt\n");
		
			/*init*/
			strcpy(ptr_queryRawScore[i].fileName, "NN");
			strcpy(ptr_queryRawScore[i].structureName, "NN");
			strcpy(ptr_queryRawScore[i].classId,"NN");
			strcpy(ptr_queryRawScore[i].chainId,"NN");
			strcpy(ptr_queryRawScore[i].windowInfo,"-1");

			ptr_queryRawScore[i].chainNr = -1171;
			ptr_queryRawScore[i].chainLen = -1171;
			ptr_queryRawScore[i].nrOfWindows = -1171;

			ptr_queryRawScore[i].score = 0.0;
			ptr_queryRawScore[i].pValue = 1.0;

			ptr_queryRawScore[i].optionalValue = -9999;

		}

	} else { //just re-init

		ptr_queryRawScore = *ptr_ptr_queryRawScore;

		for (i=0; i< topN; i++){
		
			/*init*/
			strcpy(ptr_queryRawScore[i].fileName, "NN");
			strcpy(ptr_queryRawScore[i].structureName, "NN");
			strcpy(ptr_queryRawScore[i].classId,"NN");
			strcpy(ptr_queryRawScore[i].chainId,"NN");
			strcpy(ptr_queryRawScore[i].windowInfo,"-1");

			ptr_queryRawScore[i].chainNr = -1171;
			ptr_queryRawScore[i].chainLen = -1171;
			ptr_queryRawScore[i].nrOfWindows = -1171;

			ptr_queryRawScore[i].score = 0.0;
			ptr_queryRawScore[i].pValue = 1.0;

			ptr_queryRawScore[i].optionalValue = -9999;
		}

		*ptr_ptr_queryRawScore = ptr_queryRawScore;

	}

	if(alreadyAllocatedTo_windowInfoSize == 0 ||windowInfoSize > alreadyAllocatedTo_windowInfoSize || topN > alreadyAllocatedTo_topN){
		
		*ptr_ptr_queryRawScore = ptr_queryRawScore;

	}

	return 1;
}


int reinit_queryRawScore(struct queryRawScore *ptr_queryRawScore, int topN){

	int i = 0;

	for (i=0; i< topN; i++){

		strcpy(ptr_queryRawScore[i].fileName, "NN");
		strcpy(ptr_queryRawScore[i].structureName, "NN");
		strcpy(ptr_queryRawScore[i].classId,"NN");
		strcpy(ptr_queryRawScore[i].chainId,"NN");
		strcpy(ptr_queryRawScore[i].windowInfo,"-1");

		ptr_queryRawScore[i].chainNr = -1171;
		ptr_queryRawScore[i].chainLen = -1171;
		ptr_queryRawScore[i].nrOfWindows = -1171;

		ptr_queryRawScore[i].score = 0.0;
		ptr_queryRawScore[i].pValue = 1.0;
		ptr_queryRawScore[i].optionalValue = -9999;

	}

	return 0;

}


/*In this allocation procedure (for blocks of records of DB results) extension of previously allocated 
memory (to chunkNr = 0) can also be made. For reallocation chunkNr > 0.*/ 
int alloc_init_DBresultsWindows(struct db_I_windows_order2_ptr **ptr_ptr_dbResults, int chunkNr, int chunkSize, int nrOfWindows){

	int i = 0;
	int m = 0;
	int startAt;
	int endAt;

	struct db_I_windows_order2_ptr *ptr_dbResults;

	int order = 2;

	if (chunkNr == 0){
		/*allocate initial memory:*/
		ptr_dbResults = (struct db_I_windows_order2_ptr *) malloc (chunkSize*sizeof(struct db_I_windows_order2_ptr));
		//printf("I've allocated initially\n");
	}
	else{ //reallocate another chunk

		ptr_dbResults = (struct db_I_windows_order2_ptr *) realloc(*ptr_ptr_dbResults,  (chunkNr+1)*chunkSize*sizeof(struct db_I_windows_order2_ptr));
		//printf("I've re-allocated\n");
		
	}

	/*initialize*/
	startAt = chunkNr*chunkSize;
	endAt = (chunkNr + 1)*chunkSize;
	printf("I start at: %d and end at %d\n", startAt, endAt);
	for (m=startAt; m< endAt; m++){

		ptr_dbResults[m].fileName = (char *) malloc(sizeof(char)*fileNameCharSize);
		ptr_dbResults[m].structureName = (char *) malloc(sizeof(char)*structureNameCharSize);
		ptr_dbResults[m].classId = (char *) malloc (sizeof(char)*classIdCharSize);
		ptr_dbResults[m].chainId = (char *) malloc (sizeof(char)*chainIdCharSize); 

		/*ptr_dbResults[m].windowsNr = (int *) malloc((nrOfWindows+1)*sizeof(int)); 
		ptr_dbResults[m].segIndices = (int **) malloc((nrOfWindows+1)*sizeof(int *)); */

		ptr_dbResults[m].window = (struct window *) malloc((nrOfWindows+1)*sizeof(window)); 
		
		/*init*/
		strcpy(ptr_dbResults[m].fileName, "NN");
		strcpy(ptr_dbResults[m].structureName, "NN");
		strcpy(ptr_dbResults[m].classId, "NN");
		strcpy(ptr_dbResults[m].chainId, "NN");

		ptr_dbResults[m].chainNr = -1;
		ptr_dbResults[m].chainLen = -1;
		ptr_dbResults[m].nrOfWindows = nrOfWindows;

		for (i=0; i <= nrOfWindows; i++){

			//ptr_dbResults[m].windowsNr[i] = 0; 

			//ptr_dbResults[m].segIndices[i] = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

			///*init segment indices to zero everywhere:*/
			//ptr_dbResults[m].segIndices[i][0] = 0;
			//ptr_dbResults[m].segIndices[i][1] = 0;

			ptr_dbResults[m].window[i].windowNr = 0; 

			ptr_dbResults[m].window[i].segIndices = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

			/*init segment indices to zero everywhere:*/
			ptr_dbResults[m].window[i].segIndices[0] = 0;
			ptr_dbResults[m].window[i].segIndices[1] = 0;

		}

		if (order == 1){

			/*allocate memory for double arrays of w and I12 values (order <= 1 measures)*/
			ptr_dbResults[m].I12 = (double *) malloc ((nrOfWindows+1)*sizeof(double));
			
			ptr_dbResults[m].Ia12 = (double *) malloc ((nrOfWindows+1)*sizeof(double));
			
			for (i=0; i <= nrOfWindows; i++){
			
				/* Initialize to zero everywhere: */
					ptr_dbResults[m].I12[i] = 0.0;
					ptr_dbResults[m].Ia12[i] = 0.0;

			}
		}

		if (order == 2){

		/*allocate memory for double arrays for measures of order <= 2: */

			ptr_dbResults[m].I12 = (double *) malloc ((nrOfWindows+1)*sizeof(double));

			/*absolute value versions*/
			ptr_dbResults[m].Ia12 = (double *) malloc ((nrOfWindows+1)*sizeof(double));
		
			ptr_dbResults[m].I1234_full = (double *) malloc ((nrOfWindows+1)*sizeof(double));
			ptr_dbResults[m].I1324_full = (double *) malloc ((nrOfWindows+1)*sizeof(double));
			ptr_dbResults[m].I1423 = (double *) malloc ((nrOfWindows+1)*sizeof(double));

			ptr_dbResults[m].Ia12a34_full = (double *) malloc ((nrOfWindows+1)*sizeof(double));
			ptr_dbResults[m].Ia13a24_full = (double *) malloc ((nrOfWindows+1)*sizeof(double));
			ptr_dbResults[m].Ia14a23 = (double *) malloc ((nrOfWindows+1)*sizeof(double));

			ptr_dbResults[m].Ia1234_full = (double *) malloc ((nrOfWindows+1)*sizeof(double));
			ptr_dbResults[m].Ia1324_full = (double *) malloc ((nrOfWindows+1)*sizeof(double));
			ptr_dbResults[m].Ia1423 = (double *) malloc ((nrOfWindows+1)*sizeof(double));

			ptr_dbResults[m].I12a34_full = (double *) malloc ((nrOfWindows+1)*sizeof(double));
			ptr_dbResults[m].I13a24_full = (double *) malloc ((nrOfWindows+1)*sizeof(double));
			ptr_dbResults[m].I14a23 = (double *) malloc ((nrOfWindows+1)*sizeof(double));

			/* Initialize to zero everywhere: */
			for (i=0; i <= nrOfWindows; i++){
				
				ptr_dbResults[m].I12[i] = 0.0;
				ptr_dbResults[m].Ia12[i] = 0.0;

				ptr_dbResults[m].I1234_full[i] = 0.0;
				ptr_dbResults[m].I1324_full[i] = 0.0;
				ptr_dbResults[m].I1423[i] = 0.0;

				ptr_dbResults[m].Ia12a34_full[i] = 0.0;
				ptr_dbResults[m].Ia13a24_full[i] = 0.0;
				ptr_dbResults[m].Ia14a23[i] = 0.0;

				ptr_dbResults[m].Ia1234_full[i] = 0.0;
				ptr_dbResults[m].Ia1324_full[i] = 0.0;
				ptr_dbResults[m].Ia1423[i] = 0.0;

				ptr_dbResults[m].I12a34_full[i] = 0.0;
				ptr_dbResults[m].I13a24_full[i] = 0.0;
				ptr_dbResults[m].I14a23[i] = 0.0;

			}
		}

	}

	*ptr_ptr_dbResults = ptr_dbResults;

	return 1;
}

/*as alloc_init_DBresultsWindows, but for window pairs:*/
int alloc_init_DBresultsWindowPairs(struct db_I_windowPairs_order2_ptr **ptr_ptr_dbResults, int chunkNr, int chunkSize, int nrOfWindows, int nrOfWindowPairs){

	int i = 0, j = 0;
	int m = 0;
	int startAt;
	int endAt;

	struct db_I_windowPairs_order2_ptr *ptr_dbResults;

	int order = 2;

	int allocToNrOfWindowPairs =1;

	/*Would prefer to allocate to a "simplex" of size nrOfWindowPairs, but app'ly has to be the square*/
	allocToNrOfWindowPairs = (nrOfWindows+1)*(nrOfWindows+1); //Obs: we could do with allocating to a triangular shape of size N*(N+1)/2, but the allocation seems only to work for rectangulars 

	if(allocToNrOfWindowPairs < nrOfWindowPairs){

		printf("Warning: allocToNrOfWindowPairs is %d so lower than nrOfWindowPairs (%d)\n",allocToNrOfWindowPairs, nrOfWindowPairs);
	
	}

	if (chunkNr == 0){
		/*allocate initial memory:*/
		ptr_dbResults = (struct db_I_windowPairs_order2_ptr *) malloc (chunkSize*sizeof(struct db_I_windowPairs_order2_ptr));
		//printf("I've allocated initially\n");
	}
	else{ //reallocate another chunk

		ptr_dbResults = (struct db_I_windowPairs_order2_ptr *) realloc(*ptr_ptr_dbResults,  (chunkNr+1)*chunkSize*sizeof(struct db_I_windowPairs_order2_ptr));
		//printf("I've re-allocated\n");
		
	}

	/*initialize*/
	startAt = chunkNr*chunkSize;
	endAt = (chunkNr + 1)*chunkSize;
	printf("I start at: %d and end at %d\n", startAt, endAt);
	for (m=startAt; m< endAt; m++){

		ptr_dbResults[m].fileName = (char *) malloc(sizeof(char)*500);
		ptr_dbResults[m].structureName = (char *) malloc(sizeof(char)*100);
		ptr_dbResults[m].classId = (char *) malloc (sizeof(char)*20);
		ptr_dbResults[m].chainId = (char *) malloc (sizeof(char)*3);

		//ptr_dbResults[m].windowsNr_1 = (int *) malloc((nrOfWindows+1)*sizeof(int)); /*_1: refers to the 1st window in the pair*/
		//ptr_dbResults[m].windowsNr_2 = (int *) malloc((nrOfWindows+1)*sizeof(int)); /*_2: refers to the 2nd window in the pair*/
		//ptr_dbResults[m].segIndices_1 = (int **) malloc((nrOfWindows+1)*sizeof(int *)); 
		//ptr_dbResults[m].segIndices_2 = (int **) malloc((nrOfWindows+1)*sizeof(int *)); 

		ptr_dbResults[m].windowPair = (struct windowPair **) malloc((nrOfWindows+1)*(nrOfWindows+1)*sizeof(struct windowPair));
		
		/*init*/
		strcpy(ptr_dbResults[m].fileName, "NN");
		strcpy(ptr_dbResults[m].structureName, "NN");
		strcpy(ptr_dbResults[m].chainId, "NN");
		strcpy(ptr_dbResults[m].classId, "NN");

		ptr_dbResults[m].chainNr = -1;
		ptr_dbResults[m].chainLen = -1;
		ptr_dbResults[m].nrOfWindows = nrOfWindows;
		ptr_dbResults[m].nrOfWindowPairs = nrOfWindowPairs;

		for (i=0; i <= nrOfWindows; i++){

			ptr_dbResults[m].windowPair[i] = (struct windowPair *) malloc((nrOfWindows+1)*sizeof(struct windowPair));

			//ptr_dbResults[m].windowsNr_1[i] = -1; 
			//ptr_dbResults[m].windowsNr_2[i] = -1; 

			//ptr_dbResults[m].segIndices_1[i] = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 
			//ptr_dbResults[m].segIndices_2[i] = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

			///*init segment indices to unlikely values everywhere:*/
			//ptr_dbResults[m].segIndices_1[i][0] = -1;
			//ptr_dbResults[m].segIndices_1[i][1] = -2;

			//ptr_dbResults[m].segIndices_2[i][0] = -3;
			//ptr_dbResults[m].segIndices_2[i][1] = -4;

			for (j=0; j <= nrOfWindows; j++){

				ptr_dbResults[m].windowPair[i][j].windowNr_1 = -1; 
				ptr_dbResults[m].windowPair[i][j].windowNr_2 = -1; 

				ptr_dbResults[m].windowPair[i][j].segIndices_1 = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 
				ptr_dbResults[m].windowPair[i][j].segIndices_2 = (int *) malloc(2*sizeof(int)); /*just a 2-d int array*/ 

				/*init segment indices to unlikely values everywhere:*/
				ptr_dbResults[m].windowPair[i][j].segIndices_1[0] = -1;
				ptr_dbResults[m].windowPair[i][j].segIndices_1[1] = -2;

				ptr_dbResults[m].windowPair[i][j].segIndices_2[0] = -3;
				ptr_dbResults[m].windowPair[i][j].segIndices_2[1] = -4;
			}

		}


		if (order == 1){

			/*allocate memory for double arrays of w and I12 values (order <= 1 measures)*/
			ptr_dbResults[m].I12 = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));		
			ptr_dbResults[m].Ia12 = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));
			
			for (i=0; i <= nrOfWindows; i++){

				ptr_dbResults[m].I12[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));
				ptr_dbResults[m].Ia12[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));

				for (j=0; j <= nrOfWindows; j++){
			
				/* Initialize to zero everywhere: */
					ptr_dbResults[m].I12[i][j] = 0.0;
					ptr_dbResults[m].Ia12[i][j] = 0.0;
				}
			}
		}

		if (order == 2){

		/*allocate memory for double arrays for measures of order <= 2: */

			ptr_dbResults[m].I12 = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));		
			ptr_dbResults[m].Ia12 = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));
		
			ptr_dbResults[m].I1234_full = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));
			ptr_dbResults[m].I1324_full = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));
			ptr_dbResults[m].I1423 = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));

			ptr_dbResults[m].Ia12a34_full = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));
			ptr_dbResults[m].Ia13a24_full = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));
			ptr_dbResults[m].Ia14a23 = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));

			ptr_dbResults[m].Ia1234_full = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));
			ptr_dbResults[m].Ia1324_full = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));
			ptr_dbResults[m].Ia1423 = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));

			ptr_dbResults[m].I12a34_full = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));
			ptr_dbResults[m].I13a24_full = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));
			ptr_dbResults[m].I14a23 = (double **) malloc (allocToNrOfWindowPairs*sizeof(double));
			//ptr_dbResults[m].I13a24_full = (double **) malloc ((nrOfWindows + 1)*sizeof(double *));
			//ptr_dbResults[m].I12a34_full = (double **) malloc ((nrOfWindows + 1)*sizeof(double *));
			//ptr_dbResults[m].I14a23 = (double **) malloc ((nrOfWindows + 1)*sizeof(double *));

			/* Initialize to zero everywhere: */
			for (i=0; i < nrOfWindows; ++i){
			
				//printf("I'm allocating at %d of nrOfWindows: %d \n", i, nrOfWindows);
				ptr_dbResults[m].I12[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));
				ptr_dbResults[m].Ia12[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));

				ptr_dbResults[m].I1234_full[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));
				ptr_dbResults[m].I1324_full[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));
				ptr_dbResults[m].I1423[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));

				ptr_dbResults[m].Ia12a34_full[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));
				ptr_dbResults[m].Ia13a24_full[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));
				ptr_dbResults[m].Ia14a23[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));

				ptr_dbResults[m].Ia1234_full[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));
				ptr_dbResults[m].Ia1324_full[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));
				ptr_dbResults[m].Ia1423[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));

				ptr_dbResults[m].I12a34_full[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));
				ptr_dbResults[m].I13a24_full[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));
				ptr_dbResults[m].I14a23[i] = (double *) malloc ((nrOfWindows + 1)*sizeof(double));

				for (j=0; j < nrOfWindows; j++){
			
					ptr_dbResults[m].I12[i][j] = 0.0;
					ptr_dbResults[m].Ia12[i][j] = 0.0;

					ptr_dbResults[m].I1234_full[i][j] = 0.0;
					ptr_dbResults[m].I1324_full[i][j] = 0.0;
					ptr_dbResults[m].I1423[i][j] = 0.0;

					ptr_dbResults[m].Ia12a34_full[i][j] = 0.0;
					ptr_dbResults[m].Ia13a24_full[i][j] = 0.0;
					ptr_dbResults[m].Ia14a23[i][j] = 0.0;

					ptr_dbResults[m].Ia1234_full[i][j] = 0.0;
					ptr_dbResults[m].Ia1324_full[i][j] = 0.0;
					ptr_dbResults[m].Ia1423[i][j] = 0.0;

					ptr_dbResults[m].I12a34_full[i][j] = 0.0;
					ptr_dbResults[m].I13a24_full[i][j] = 0.0;
					ptr_dbResults[m].I14a23[i][j] = 0.0;

				}

			}

			//printf("file nr %d done\n", m);
			//getchar();	

		}

	}

	*ptr_ptr_dbResults = ptr_dbResults;


	return 1;
}

/*allocate and init ptr to binned I-windows, ie to "words of GIs"*/
int alloc_init_ptr_binned_I_windows(struct binned_I_windows **ptr_ptr_binned_I_windows, int nrOfRecords, int alreadyAllocatedTo_nrOfRecords, int nrOfEntries){

	int returnVal = 0;

	struct binned_I_windows *ptr_binned_I_windows;

	int n = 0;
	int i = 0;

	if(alreadyAllocatedTo_nrOfRecords == 0){ //no mem allocated yet

		/*allocate suff memory*/
		ptr_binned_I_windows = (struct binned_I_windows *) malloc(nrOfRecords*sizeof(struct binned_I_windows)); 

		for(n = 0; n < nrOfRecords; n++){

			//ptr_binned_I_windows[n].recKey = n;
			ptr_binned_I_windows[n].fileNr = -1; //init to unlikely value
			ptr_binned_I_windows[n].windowNr = -1; //init to unlikely value
			
			ptr_binned_I_windows[n].binnedIvector = (int *) malloc(nrOfEntries*sizeof(int));

			/*init*/
			for(i=0; i < nrOfEntries; i++){
				ptr_binned_I_windows[n].binnedIvector[i] = -1; //init to unlikely value
			}

		}

		*ptr_ptr_binned_I_windows = ptr_binned_I_windows;

	}
	else if(nrOfRecords > alreadyAllocatedTo_nrOfRecords){ //need to allocate some more mem

		ptr_binned_I_windows = *ptr_ptr_binned_I_windows;

		/*allocate suff memory*/
		ptr_binned_I_windows = realloc(ptr_binned_I_windows, nrOfRecords*sizeof(struct binned_I_windows)); 

		for(n = 0; n < nrOfRecords; n++){

			//ptr_binned_I_windows[n].recKey = n;
			ptr_binned_I_windows[n].fileNr = -1; //init to unlikely value
			ptr_binned_I_windows[n].windowNr = -1; //init to unlikely value
			
			if(n >= alreadyAllocatedTo_nrOfRecords){ //not allocated yet
				ptr_binned_I_windows[n].binnedIvector = (int *) malloc(nrOfEntries*sizeof(int));
			}
			else if(n < alreadyAllocatedTo_nrOfRecords){ //not allocated yet
				ptr_binned_I_windows[n].binnedIvector = realloc(ptr_binned_I_windows[n].binnedIvector, nrOfEntries*sizeof(int));
			}

			/*init*/
			for(i=0; i < nrOfEntries; i++){
				ptr_binned_I_windows[n].binnedIvector[i] = -1; //init to unlikely value
			}

		}

		*ptr_ptr_binned_I_windows = ptr_binned_I_windows;

	}

	return returnVal;

}

/*ptr to array over windows of ptr to binned I-windows; to hold the match set for all windows in a query (nrOfRecords intended to be the nr of rec's in a DB)*/
int alloc_init_ptr_ptr_binned_I_windows(struct binned_I_windows ***ptr_ptr_ptr_binned_I_windows, int nrOfRecords, int nrOfEntries, int nrOfWindows, int alreadyAllocatedTo_nrOfWindows){

	int returnVal = 0;
	
	int k = 0, n = 0, i = 0;
	struct binned_I_windows **ptr_ptr_binned_I_windows;

	/*allocate memory to pointer ptr_I_values and initialize:*/
	if(alreadyAllocatedTo_nrOfWindows == 0){ //no mem allocated yet

		ptr_ptr_binned_I_windows = (struct binned_I_windows **) malloc (nrOfWindows*sizeof(struct binned_I_windows *));
		
		for (k=0; k < nrOfWindows; k++){
			
			returnVal = alloc_init_ptr_binned_I_windows(&ptr_ptr_binned_I_windows[k], nrOfRecords, 0, nrOfEntries); 

		}

		*ptr_ptr_ptr_binned_I_windows = ptr_ptr_binned_I_windows;

	}
	else if(nrOfWindows > alreadyAllocatedTo_nrOfWindows){ //need to allocate some more mem

		ptr_ptr_binned_I_windows = *ptr_ptr_ptr_binned_I_windows;

		ptr_ptr_binned_I_windows = realloc(ptr_ptr_binned_I_windows, nrOfWindows*sizeof(struct binned_I_windows *));
		
		for (k=alreadyAllocatedTo_nrOfWindows; k < nrOfWindows; k++){ //OBS the start id of the loop!
			
			returnVal = alloc_init_ptr_binned_I_windows(&ptr_ptr_binned_I_windows[k], nrOfRecords, 0, nrOfEntries); 

		}
		/*To be on the safe side we re-init for all k's allocated to before:*/
		for (k=0; k < alreadyAllocatedTo_nrOfWindows; k++){

			for(n = 0; n < nrOfRecords; n++){

				//ptr_binned_I_windows[n].recKey = n;
				ptr_ptr_binned_I_windows[k][n].fileNr = -1; //init to unlikely value
				ptr_ptr_binned_I_windows[k][n].windowNr = -1; //init to unlikely value
				
				/*init*/
				for(i=0; i < nrOfEntries; i++){
					ptr_ptr_binned_I_windows[k][n].binnedIvector[i] = -1; //init to unlikely value
				}

			}


		}
		
		*ptr_ptr_ptr_binned_I_windows = ptr_ptr_binned_I_windows;

	}

	return returnVal;
}

/*allocate and init ptr to binned I-windows, ie to "words of GIs"*/
int alloc_init_ptr_binned_I_windowPairs(struct binned_I_windowPairs **ptr_ptr_binned_I_windowPairs, int nrOfRecords, int alreadyAllocatedTo_nrOfRecords, int nrOfEntriesForPairs){

	int returnVal = 0;

	struct binned_I_windowPairs *ptr_binned_I_windowPairs;

	int n = 0;
	int i = 0;

	if(alreadyAllocatedTo_nrOfRecords == 0){ //no mem allocated yet

		/*allocate suff memory*/
		ptr_binned_I_windowPairs = (struct binned_I_windowPairs *) malloc(nrOfRecords*sizeof(struct binned_I_windowPairs)); 
		
		/*init*/
		for(n = 0; n < nrOfRecords; n++){

			//ptr_binned_I_windows[n].recKey = n;
			ptr_binned_I_windowPairs[n].fileNr = -1; //init to unlikely value
			ptr_binned_I_windowPairs[n].windowNr_1 = -1; //init to unlikely value
			ptr_binned_I_windowPairs[n].windowNr_2 = -1; //init to unlikely value
			
			ptr_binned_I_windowPairs[n].binnedIvector = (int *) malloc(nrOfEntriesForPairs*sizeof(int));

			for(i=0; i < nrOfEntriesForPairs; i++){
				ptr_binned_I_windowPairs[n].binnedIvector[i] = -1; //init to unlikely value
			}

		}

		*ptr_ptr_binned_I_windowPairs = ptr_binned_I_windowPairs;

	}
	else if(nrOfRecords > alreadyAllocatedTo_nrOfRecords){ //need to allocate some more mem

		ptr_binned_I_windowPairs = *ptr_ptr_binned_I_windowPairs;

		ptr_binned_I_windowPairs = realloc(ptr_binned_I_windowPairs, nrOfRecords*sizeof(struct binned_I_windowPairs)); 


		/*init*/
		for(n = 0; n < nrOfRecords; n++){

			//ptr_binned_I_windows[n].recKey = n;
			ptr_binned_I_windowPairs[n].fileNr = -1; //init to unlikely value
			ptr_binned_I_windowPairs[n].windowNr_1 = -1; //init to unlikely value
			ptr_binned_I_windowPairs[n].windowNr_2 = -1; //init to unlikely value
			
			if(n >= alreadyAllocatedTo_nrOfRecords){
				ptr_binned_I_windowPairs[n].binnedIvector = (int *) malloc(nrOfEntriesForPairs*sizeof(int));
			}
			else if(n < alreadyAllocatedTo_nrOfRecords){ //not allocated yet
				ptr_binned_I_windowPairs[n].binnedIvector = realloc(ptr_binned_I_windowPairs[n].binnedIvector, nrOfEntriesForPairs*sizeof(int));
			}

			for(i=0; i < nrOfEntriesForPairs; i++){
				ptr_binned_I_windowPairs[n].binnedIvector[i] = -1; //init to unlikely value
			}

		}

		*ptr_ptr_binned_I_windowPairs = ptr_binned_I_windowPairs;

	}

	//printf("(saaledes stadig ganske ... ) ... in allocation to ptr_binned_I_windowPairs, n was: %d while nrOfRecords:%d\n", n, nrOfRecords);

	return returnVal;

}

/*Nowhere in use!: ptr to array over window pairs of ptr to binned I-windowPairs; to hold the match set for all window pairs in a query*/
int alloc_init_ptr_ptr_binned_I_windowPairs(struct binned_I_windowPairs ***ptr_ptr_ptr_binned_I_windowPairs, int nrOfRecords, int nrOfEntries, int nrOfWindowPairs){

	int returnVal = 0;
	
	int k = 0;
	struct binned_I_windowPairs **ptr_ptr_binned_I_windowPairs;

	/*allocate memory to pointer ptr_I_values and initialize:*/
	ptr_ptr_binned_I_windowPairs = (struct binned_I_windowPairs **) malloc (nrOfWindowPairs*sizeof(struct binned_I_windowPairs *));
	
	for (k=0; k < nrOfWindowPairs; k++){
		
		returnVal = alloc_init_ptr_binned_I_windowPairs(&ptr_ptr_binned_I_windowPairs[k], nrOfRecords, 0, nrOfEntries); 

	}


	*ptr_ptr_ptr_binned_I_windowPairs = ptr_ptr_binned_I_windowPairs;

	return returnVal;
}






/*Main program blocks. Change the main by changing the fct called, chosing between
the main_* main-functions below.*/

int main_makeDB(){

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
	int use_scop_b = 0; //if using data from SCOP (particular data structure when loading files)
	int use_cath_b = 1;//if using data from CATH (particular data structure when loading files)

	int loadFromSubDirs_b = 0; //if PDB files sit in sub-directories at some level of the "root" directory, dirName, use loadFromSubDirs_b = 1; else 0.

	int omitSingletons_b = 0;

	//char DBName[100] = "3chains";
	//char DBName[100] = "SCOPe_ASTRAL_2.06_40pct";


	//char DBName[100] = "top100";
	//char DBName[100] = "top8000";
	//char DBName[100] = "SCOPe_ASTRAL_2.06_40pct";
	//char DBName[100] = "SCOP_test";
	//char DBName[100] = "CATH_v2_4_Kolodny";
	//char DBName[100] = "CATH_v2_4";
	char DBName[100] = "CATH_test1";
	//char DBName[100] = "CATH_test";

	//char dirPath[1000] = "/home/tkj375/masters_study_pgm/projects/knots/code/top100H";
	//char dirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/Kinemage/top100H/";
	//char dirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/Kinemage/top8000_chains_70/";
	//char dirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH-v2_4/CATH2.4AllAtomsConnectedDomains";
	char dirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH_test1/";
	//char dirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH_test/";


	char CATHListPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH-v2_4_Kolodny/subset_list_web.txt";
	//char CATHListPath[1000] = "NN";

	//char outputPath[1000]  = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top100";
	//char outputPath[1000]  = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000";
	//char outputPath[1000]  = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH-v2_4";
	//char outputPath[1000]  = "/isdata/kroghgrp/tkj375/GISA_unix/results/v9/CATH_test";
	char outputPath[1000]  = "/isdata/kroghgrp/tkj375/GISA_unix/results/v9/CATH_test1";
	//char outputPath[1000]  = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top100";
	//char outputPath[1000]  = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000/";

	//char fileNameChains[1000] = "//chains_subChainPairs_top100.txt";
	//char fileNameChains[1000] = "/chains_top100.txt";
	//char fileNameChains[1000] = "/chains_top8000.txt";
	char fileNameChains[1000] = "/CATH-v2_4.txt";
	//char fileNameChains[1000] = "/CATH_test.txt";
	//char fileNameChains[1000] = "/CATH_test1.txt";

	/*For getting and writing out the invariant values for a set windowlength:*/
	int windowCoveringType = 0; /*which method to use for deriving a "good" window length; used in setWindowCharacteristics*/
	int get_windows_b = 1;
	int get_windowPairs_b = 0; /*at most one of get_windows_b and get_windowPairs_b should be set to 1*/
	int windowLength = 16;
	int stepSize = 2;
	int write_windows_b = 1;
	int write_windowPairs_b = 0; /*at most one of write_windows_b and write_windowPairs_b should be set to 1; corr to get_windos_b resp. get_windowPairs_b*/
	int writeOnlyDisjointPairs_b = 1;
	int maxNrOfWindows = 1;
	int write_chains_b = 0;
	char pathChainsDB[2000] = "\0";
	char fileNameOutChainsDB[1000] =  "/chains_";
	char fileNameOutExtensionChains[20] = ".txt";

	
	if (write_chains_b !=0){
		/*generate file name for writing out the data base chains:*/
		strcat(pathChainsDB , outputPath);
		strcat(pathChainsDB , fileNameOutChainsDB);
		strcat(pathChainsDB , DBName);
		strcat(pathChainsDB , fileNameOutExtensionChains);

		printf("DB chains will be written to: %s\n", pathChainsDB);
		
		/*clear the file by opening it in w-mode and then closing it again*/
		fclose(fopen(pathChainsDB , "w"));
	}

	printf("I'll be running a DB-run. Press enter to carry on."); 
	getchar();

	returnVal = computeGI_windows(DBName, dirPath, use_scop_b, use_cath_b, CATHListPath, loadFromSubDirs_b, outputPath, maxChainLength, omitSingletons_b ,
			  order, incl_abs_b, full_b, split_b, write_final_b, write_all_b, 
			  closed_loops_b, closedLoopLength, closedLoopDist, pokeLength, invValSubChainPairs_b, writeSubChainPairsAll_b, subChainLength, 
			  windowCoveringType, get_windows_b, get_windowPairs_b, windowLength, stepSize, write_windows_b, write_windowPairs_b, 
			  writeOnlyDisjointPairs_b, maxNrOfWindows, write_chains_b, pathChainsDB,
			  print_b, print_basic_b);


	return returnVal;

}

int main_rarity0(){

	int returnVal = 0;

	int maxChainLength = 1500; /*chains longer than this will be discarded; can be set acc to RAM availability*/
	
	/*compute invariants up to and including this order IN THIS FCT ONLY ORDER 1 IS USED:*/
	int order = 1;
	/*inlcude absolute value versions of order one and two invariants:*/
	int incl_abs_b = 1;
	/* compute full order 2 measures across the simplex*/
	int full_b = 1;


	/*PLACE HOLDER VALUES, NOT IN USE HERE:*/
	int closed_loops_b = 0;
	int closedLoopLength = 30;
	double closedLoopDist = 7*7; /*square of distance in nr of Ångstrøm -- we use square of dist*/
	/*For computing and recording the mutual invariant values (e.g writhe) of pairs of sub-chains (1) or not (0); aimed
	at more generally searching for particular geometric shapes than with the closed-loops
	search:*/
	int invValSubChainPairs_b = 0; /*if set to 1: computes the mutual writhe; if set to 2 computes the mutual invariant value for order 1 and order 2 "full" invariant (takes setting order = 3 or full_b = 1); if set to 3: all invariants (takes setting order = 3)*/
	int writeSubChainPairsAll_b = 0; /*if set to 1 the results for all sub-chain pairs will be attempted to be written out (may be lots of data!)*/
	int subChainLength = 30; /* the lenght of the sub-chains in the pairs considered*/


	/*THESE SHOULD ONLY BE CHANGED IF YOU REALLY WANT IT:
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
	/*Want to halt the execution at various places?*/
	int halted_b = 0; 

	/*"GENUINE" PARAMETERS, TO BE SET*/
	/*queries in this drectory*/
	//char queriesName[100] = "SCOP_test";
	//char queriesName[100] = "testQuery1";
	//char queriesName[100] = "3chains";
	//char queriesName[100] = "3links";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/Kinemage/3links";
	//char queriesName[100] = "top100";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/Kinemage/top100H";
	//char queriesDirPath[1000] = "C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//top100H";
	//char queriesName[100] = "top8000";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/Kinemage/top8000_chains_70";
	//char queriesDirPath[1000] = "C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//top8000_chains_70";
	//char queriesName[100] = "SCOPe_ASTRAL_206_40pct";
	//char queriesDirPath[1000] = "D:\\Bioinformatics\\data\\structural\\SCOP\\SCOP_206_40pct";
	//char queriesName[100] = "CATH_test3";
	//char queriesName[100] = "CATH_v2_4_Kolodny";
	//char queriesDirPath[1000] = "D://Bioinformatics//isdata/kroghgrp//structural//CATH//CATH-v2_4//CATH2_4AllAtomsConnectedDomains";
	//char queriesName[100] = "CATH_test";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH_test";
	char queriesName[100] = "CATH_test1";
	char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH_test1";
	//char queriesName[100] = "CATH-v2_4";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH-v2_4/CATH2.4AllAtomsConnectedDomains";
	
	char CATHListPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH-v2_4_Kolodny/subset_list_web.txt";
	//char CATHListPath[1000] = "NN";

	/*name of folder holding the PDB files (the input):*/
	int loadFromSubDirs_b = 0; //if PDB files sit in sub-directories at some level of the "root" directory, dirName, use loadFromSubDirs_b = 1; else 0.
	
	int chunkSize = 100; /* load in chunks of this size (nr of structures) */

	int use_scop_b = 0; /*if using query data from SCOP (particular data structure when loading files)*/
	int use_cath_b = 0; /*if using query data from CATH (particular data structure when loading files)*/

	//if using data from CATH (particular data structure when loading files)
	//char DBName[100] = "3chains";
	//char DBName[100] = "top100";
	//char DBName[100] = "top8000";
	//char DBName[100] = "SCOPe_ASTRAL_206_40pct";
	//char DBName[100] = "CATH_v2_4_Kolodny";
	//char DBName[100] = "CATH_test";
	char DBName[100] = "CATH_test1";
	//char DBName[100] = "CATH-v2_4";
	//char DBName[100] = "SCOP_test";


	/*For window results:*/
	//char DBResultsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top100/Invariants_windowlgth_16_2_computeGI_windows_top100_winCovType1.txt";

	//char DBResultsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000/Invariants_windowlgth_16_2_computeGI_windows_top8000_winCovType1.txt";


	/*For window pair results:*/
	//char DBResultsPairsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top100/Invariants_Pairs_windowlgth_16_2_computeGI_windows_top100_winCovType1.txt";
	//char DBResultsPairsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000/Invariants_Pairs_windowlgth_16_2_order_1_computeGI_windows_top8000_winCovType0_onlyDisjointPairs1.txt";
	//char DBResultsPairsFilePath[1000] = "NN";
	//char DBResultsPairsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH_test/Invariants_Pairs_windowlgth_16_2_order_2_computeGI_windows_CATH_test_winCovType0_onlyDisjointPairs1.txt";
	char DBResultsPairsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v9/CATH_test1/Invariants_Pairs_windowlgth_20_4_order_1_computeGI_windows_CATH_test1_winCovType0_onlyDisjointPairs1.txt";

	//char DBResultsPairsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH-v2_4/Invariants_Pairs_windowlgth_16_2_order_1_computeGI_windows_CATH_v2_4_winCovType0_onlyDisjointPairs1.txt";

	//char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/3links"; 
	//char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top100";
	//char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000";
	//char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH_test";
	char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v9/CATH_test1";
    //char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH-v2_4";

	
	int onlyDisjointPairs_b = 1; 

	int normalize_b = 0;  /*only in its pla ce when using more than one invariant for the scan/matching; for this function usually only the (mutual) writhe is used, so the normalization can (should) be skipped*/

	int useOnlyMaxMutVal_b = 0;
	int signed_b = 0;
	int absMutValScorePValues_b = 1; /*if set to 1 and useOnlyMaxMutVal_b = 0, provide a absMutValScoresDistr*/

	double thresholdMutual = 5.0;

	//char absMutValScoresDistrName[100] = "_CATH24MatchCATH24";
	//char absMutValScoresDistrFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH-v2_4/RarityScan0_ScoresPairs_windowslgth_16_2_order_1_CATH-v2_4_CATH-v2_4_notnorm_1invs__winCovType0_threshMut5.000000_scoreByAbsMutualWrithe.txt";
	
	char absMutValScoresDistrName[100] = "_CATHtest1MatchCATHtest1";
	char absMutValScoresDistrFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v9/CATH_test1/RarityScan0_ScoresPairs_windowslgth_16_2_order_1_CATH_test1_CATH_test1_notnorm_1invs__winCovType0_threshMut5.000000_scoreByAbsMutualWrithe.txt";

	//char absMutValScoresDistrName[100] = "_scoreDistrTop8000";
	//char absMutValScoresDistrFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000/RarityScan0_ScoresPairs_windowslgth_16_2_order_1_top8000_top8000_notnorm_1invs__winCovType0_threshMut10.000000_scoreByAbsMutualWrithe.txt";
	

	/*For getting and writing out the invariant values for a set windowlength:*/
	int matchWindowPairs_b = 1; /*match on window pairs (1) or not (0); if set to 1 also set matchWindows_b to 1*/

	int windowCoveringType = 0;
	int windowLength = 16;
	int stepSize = 2;

	int getAvgScore_b = 0;

	int write_windowPairs_b = 0.0; 
	
	int write_scoresPairs_b = 1; //set to 1 if scores should be output
	
	int writeWindowInfo_b = 0;

	int write_chains_b = 0;
	char pathChainsQueries[2000] = "\0";
	char fileNameOutChains[1000] =  "/RarityScan0_chains_";
	//char fileNameOutExtensionChains[20] = "_top100.txt";
	//char fileNameOutExtensionChains[20] = "_CATH_test.txt";
	char fileNameOutExtensionChains[20] = "_CATH_test1.txt";
	//char fileNameOutExtensionChains[20] = "_CATH-v2_4.txt";
 

	if (write_chains_b !=0){
		/*generate file name for writing out the data base chains:*/
		strcat(pathChainsQueries , outputPath);
		strcat(pathChainsQueries , fileNameOutChains);
		strcat(pathChainsQueries , queriesName);
		strcat(pathChainsQueries , fileNameOutExtensionChains);

		printf("Query chains will be written to: %s\n", pathChainsQueries);
		
		/*clear the file by opening it in w-mode and then closing it again*/
		fclose(fopen(pathChainsQueries , "w"));
	}

	printf("I'll be running a rarity scan, type 0, for pairs of windows. Press enter to carry on."); 
	getchar();

	returnVal = rawRarity0(halted_b , DBName, queriesName, 
					queriesDirPath, DBResultsPairsFilePath, pathChainsQueries,
					absMutValScoresDistrName, absMutValScoresDistrFilePath, onlyDisjointPairs_b, 
					useOnlyMaxMutVal_b, signed_b, absMutValScorePValues_b, thresholdMutual,
					use_scop_b, use_cath_b, CATHListPath, 
					loadFromSubDirs_b, normalize_b, outputPath, maxChainLength, 
					chunkSize, windowCoveringType, windowLength, 
					stepSize, getAvgScore_b, write_windowPairs_b, 
					write_scoresPairs_b, write_chains_b, writeWindowInfo_b, 
					print_b, print_basic_b);




	return returnVal;

}

int main_rarity1(){

	int returnVal = 0;

	int maxChainLength = 1500; /*chains longer than this will be discarded; can be set acc to RAM availability*/
	
	/*compute invariants up to and including this order:*/
	int order = 1;
	/*inlcude absolute value versions of order one and two invariants:*/
	int incl_abs_b = 1;
	/* compute full order 2 measures across the simplex*/
	int full_b = 1;

	/*Detect and record closed loops; if chosen I12 values corr 
	to the loop will be stored for alter use (finding "pokes").
	Can only be used with order >= 1 (since I12 is needed):*/
	int closed_loops_b = 0;
	int closedLoopLength = 30;
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
	/*Want to halt the execution at various places?*/
	int halted_b = 0; 


	/*queries in this drectory*/
	//char queriesName[100] = "SCOP_test";
	//char queriesName[100] = "testQuery1";
	//char queriesName[100] = "3chains";
	//char queriesName[100] = "3links";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/Kinemage/3links";
	//char queriesName[100] = "top100";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/Kinemage/top100H";
	//char queriesDirPath[1000] = "C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//top100H";
	//char queriesName[100] = "top8000";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/Kinemage/top8000_chains_70";
	//char queriesDirPath[1000] = "C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//top8000_chains_70";
	//char queriesName[100] = "SCOPe_ASTRAL_206_40pct";
	//char queriesDirPath[1000] = "D:\\Bioinformatics\\data\\structural\\SCOP\\SCOP_206_40pct";
	//char queriesName[100] = "CATH_test3";
	//char queriesName[100] = "CATH_v2_4_Kolodny";
	//char queriesDirPath[1000] = "D://Bioinformatics//isdata/kroghgrp//structural//CATH//CATH-v2_4//CATH2_4AllAtomsConnectedDomains";
	char queriesName[100] = "CATH_test";
	char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH_test";
	//char queriesName[100] = "CATH_test1";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH_test1";
	//char queriesName[100] = "CATH-v2_4";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH-v2_4/CATH2.4AllAtomsConnectedDomains";
	
	char CATHListPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH-v2_4_Kolodny/subset_list_web.txt";
	//char CATHListPath[1000] = "NN";

	/*name of folder holding the PDB files (the input):*/
	int loadFromSubDirs_b = 0; //if PDB files sit in sub-directories at some level of the "root" directory, dirName, use loadFromSubDirs_b = 1; else 0.
	
	int chunkSize = 100; /* load in chunks of this size (nr of structures) */

	int use_scop_b = 0; /*if using query data from SCOP (particular data structure when loading files)*/
	int use_cath_b = 1; /*if using query data from CATH (particular data structure when loading files)*/

	//if using data from CATH (particular data structure when loading files)
	//char DBName[100] = "3chains";
	//char DBName[100] = "top100";
	//char DBName[100] = "top8000";
	//char DBName[100] = "SCOPe_ASTRAL_206_40pct";
	//char DBName[100] = "CATH_v2_4_Kolodny";
	char DBName[100] = "CATH_test";
	//char DBName[100] = "CATH-v2_4";
	//char DBName[100] = "SCOP_test";


	/*For window results:*/
	//char DBResultsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top100/Invariants_windowlgth_16_2_computeGI_windows_top100_winCovType1.txt";

	//char DBResultsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000/Invariants_windowlgth_16_2_computeGI_windows_top8000_winCovType1.txt";


	/*For window pair results:*/
	//char DBResultsPairsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top100/Invariants_Pairs_windowlgth_16_2_computeGI_windows_top100_winCovType1.txt";
	//char DBResultsPairsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000/Invariants_Pairs_windowlgth_16_2_computeGI_windows_top8000_winCovType1.txt";
	//char DBResultsPairsFilePath[1000] = "NN";
	char DBResultsPairsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH_test/Invariants_Pairs_windowlgth_16_2_order_2_computeGI_windows_CATH_test_winCovType0_onlyDisjointPairs1.txt";
	//char DBResultsPairsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH-v2_4/Invariants_Pairs_windowlgth_16_2_order_2_computeGI_windows_CATH_v2_4_winCovType0_onlyDisjointPairs1.txt";

	//char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/3links"; 
	//char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top100";
	//char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000";
	char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH_test";
	//char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH_test1";
    //char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH-v2_4";
	
	/*score using just a simple, raw scoring fct or use additional input giving a distribution of distance-values (between
	GIs) allowing to score a value as "-log (p-value)":*/
	int rarityScorePValues_b = 0; /*if set to 1 provide rarityScore-distribution FileName; else set to 0 (and set write_scores_b = 1 to have the rarityScores written to file; this can then be used as the source of the rarityScores-distribution*/
	
	//char distDistrNamePairs[100] = "_top100MatchPairsTop8000";
	//char distDistrFileNamePairs[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000/MatchPairs_windowlgth_16_order_2_top100_top8000_1wins_norm_GISAv7_DB6_5mms_17bins_thrMut02.txt";
	char rarityScorePairsDistrName[100] = "_CATH24MatchCATH24";
	char rarityScorePairsDistrFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH-v2_4/RarityScan1_ScoresPairs_windowslgth_16_2_order_1_CATH-v2_4_CATH-v2_4_1wins_notnorm_1invs_0mmsPairs_20bins1_winCovType0_threshMut5.000000.txt";
	//char rarityScorePairsDistrFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH_test/RarityScan_Scores_windowslgth_16_2_order_2_CATH_test_CATH_test_1wins_norm_5invs_0mms_0mmsPairs_10bins2_winCovType0_CATH24MatchCATH24.txt";

	int onlyDisjointPairs_b = 1; 

	int normalize_b = 0;  /*only in its place when using more than one invariant for the scan/matching; for this function usually only the (mutual) writhe is used, so the normalization can (should) be skipped*/

	double thresholdMutual = 5.0;

	/*For getting and writing out the invariant values for a set windowlength:*/
	int matchWindowPairs_b = 1; /*match on window pairs (1) or not (0); if set to 1 also set matchWindows_b to 1*/
	int write_matchWindowPairs_b = 0; /*write the match results for each window pair to file or not (0) */

	int windowCoveringType = 0;
	int windowLength = 16;
	int stepSize = 2;

	int write_windowPairs_b = 0; 
	
	int write_scoresPairs_b = 1; //set to 1 if scores should be output
	
	int writeWindowInfo_b = 0;

	int write_chains_b = 0;
	char pathChainsQueries[2000] = "\0";
	char fileNameOutChains[1000] =  "/RarityScan_chains_";
	//char fileNameOutExtensionChains[20] = "_top100.txt";
	//char fileNameOutExtensionChains[20] = "_CATH_test.txt";
	char fileNameOutExtensionChains[20] = "_CATH-v2_4.txt";

	int nrOfWindowsForMatch = 1; /*just leave it at 1, has no influence apart from appearing in file names*/

	int allowedNrOfMismatchesPairs = 1;
	int nrOfEntriesForPairs = 1; /*  setting this to 1 means that only the mutual writhe is considered; >1: mutuals of higher order inv's too*/

	/*For generating a set of bins*/
	int binningType = 1; /*don't use type 2 since this makes all bins equally popupated (not what is wanted in a rarity search); don't use type 0 either, unless you want to go directly into the code and change the hardcoded values ...*/
	int nrOfBins = 20;
	double compensationForRounding = 0.0000009; //When using equi-percentile binning (binningType 2) bin interval end points will be shifted to the right by this amount; to prevent separation of e.g. self-hits


	///*generate the proper nrOfEntriesForPairs based on the settings of the parameters order and incl_abs_b:*/
	//if(order == 1){

	//	if(incl_abs_b != 1){ 

	//		nrOfEntriesForPairs = 1;
	//		printf("Search will be based on mutual writhe only (recommended)");

	//	}
	//	if(incl_abs_b == 1){ 

	//		nrOfEntriesForPairs = 2;
	//		printf("Search will be based on mutual writhe and its absolute value version (average crossing number)");

	//	}
	//}

	//if(order == 2){


	//	if(incl_abs_b != 1){ 

	//		nrOfEntriesForPairs = 4;
	//		printf("Search will be based on four invarinats: mutual writhe and mutual I1234, I1324 and I1423.");

	//	}
	//	if(incl_abs_b == 1){ 

	//		nrOfEntriesForPairs = 14;
	//		printf("Search will be based on 14 invariants: mutual writhe and mutual I1234, I1324 and I1423 and all their absolute value versions");

	//	}
	//}


	if (write_chains_b !=0){
		/*generate file name for writing out the data base chains:*/
		strcat(pathChainsQueries , outputPath);
		strcat(pathChainsQueries , fileNameOutChains);
		strcat(pathChainsQueries , queriesName);
		strcat(pathChainsQueries , fileNameOutExtensionChains);

		printf("Query chains will be written to: %s\n", pathChainsQueries);
		
		/*clear the file by opening it in w-mode and then closing it again*/
		fclose(fopen(pathChainsQueries , "w"));
	}

	printf("I'll be running a rarity scan, type 1, for pairs of windows. Press enter to carry on."); 
	getchar();

	returnVal = rawRarity1(halted_b, rarityScorePairsDistrFilePath, 
					rarityScorePairsDistrName, DBName, queriesName, 
					queriesDirPath, DBResultsPairsFilePath, 
					pathChainsQueries,
					onlyDisjointPairs_b,
					thresholdMutual,
					rarityScorePValues_b,  use_scop_b, 
					use_cath_b, CATHListPath, loadFromSubDirs_b, 
					normalize_b, outputPath, maxChainLength, 
					nrOfEntriesForPairs, binningType, nrOfBins, 
					compensationForRounding, allowedNrOfMismatchesPairs, nrOfWindowsForMatch, 
					chunkSize, order, incl_abs_b, 
					full_b, write_matchWindowPairs_b, 
					windowCoveringType, windowLength, stepSize, 
					write_windowPairs_b, write_scoresPairs_b, write_chains_b, 
					writeWindowInfo_b, print_b, print_basic_b);



	return returnVal;

}

int main_rarity2(){

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
	/*Want to halt the execution at various places?*/
	int halted_b = 0; 


	/*queries in this drectory*/
	//char queriesName[100] = "SCOP_test";
	//char queriesName[100] = "testQuery1";
	//char queriesName[100] = "3chains";
	//char queriesName[100] = "top100";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/Kinemage/top100H";
	//char queriesDirPath[1000] = "C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//top100H";
	//char queriesName[100] = "top8000";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/Kinemage/top8000_chains_70";
	//char queriesDirPath[1000] = "C://Users//Christian//BioInformatics//masters_study_pgm//projects//knots//code//top8000_chains_70";
	//char queriesName[100] = "SCOPe_ASTRAL_206_40pct";
	//char queriesDirPath[1000] = "D:\\Bioinformatics\\data\\structural\\SCOP\\SCOP_206_40pct";
	//char queriesName[100] = "CATH_test3";
	//char queriesName[100] = "CATH_v2_4_Kolodny";
	//char queriesDirPath[1000] = "D://Bioinformatics//isdata/kroghgrp//structural//CATH//CATH-v2_4//CATH2_4AllAtomsConnectedDomains";
	char queriesName[100] = "CATH_test";
	char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH_test";
	//char queriesName[100] = "CATH_test1";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH_test1";
	//char queriesName[100] = "CATH-v2_4";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH-v2_4/CATH2.4AllAtomsConnectedDomains";
	
	char CATHListPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH-v2_4_Kolodny/subset_list_web.txt";
	//char CATHListPath[1000] = "NN";

	/*name of folder holding the PDB files (the input):*/
	int loadFromSubDirs_b = 0; //if PDB files sit in sub-directories at some level of the "root" directory, dirName, use loadFromSubDirs_b = 1; else 0.
	
	int chunkSize = 100; /* load in chunks of this size (nr of structures) */

	int use_scop_b = 0; /*if using data from SCOP (particular data structure when loading files)*/
	int use_cath_b = 1;

	//if using data from CATH (particular data structure when loading files)
	//char DBName[100] = "3chains";
	//char DBName[100] = "top100";
	//char DBName[100] = "top8000";
	//char DBName[100] = "SCOPe_ASTRAL_206_40pct";
	//char DBName[100] = "CATH_v2_4_Kolodny";
	char DBName[100] = "CATH_test";
	//char DBName[100] = "CATH-v2_4";
	//char DBName[100] = "SCOP_test";


	/*For window results:*/
	//char DBResultsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top100/Invariants_windowlgth_16_2_computeGI_windows_top100_winCovType1.txt";

	//char DBResultsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000/Invariants_windowlgth_16_2_computeGI_windows_top8000_winCovType1.txt";

	//char DBResultsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH-v2_4/Invariants_windowlgth_16_2_computeGI_windows_CATH_v2_4_winCovType1.txt";
	//char DBResultsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH-v2_4/Invariants_windowlgth_32_2_computeGI_windows_CATH_v2_4_winCovType1.txt";
	//char DBResultsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH_test/Invariants_windowlgth_16_2_order_2_computeGI_windows_CATH_test_winCovType0.txt";
	char DBResultsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH_test/Invariants_windowlgth_30_4_order_2_computeGI_windows_CATH_test_winCovType0.txt";


	/*For window pair results:*/
	//char DBResultsPairsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top100/Invariants_Pairs_windowlgth_16_2_computeGI_windows_top100_winCovType1.txt";
	//char DBResultsPairsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000/Invariants_Pairs_windowlgth_16_2_computeGI_windows_top8000_winCovType1.txt";
	//char DBResultsPairsFilePath[1000] = "NN";
	//char DBResultsPairsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH_test/Invariants_Pairs_windowlgth_16_2_order_2_computeGI_windows_CATH_test_winCovType0.txt";
	char DBResultsPairsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH_test/Invariants_Pairs_windowlgth_30_4_order_2_computeGI_windows_CATH_test_winCovType0_onlyDisjointPairs0.txt";
	 
	//char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top100";
	//char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000";
	char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH_test";
	//char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH_test1";
    //char outputPath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH-v2_4";
	
	/*score using just a simple, raw scoring fct or use additional input giving a distribution of distance-values (between
	GIs) allowing to score a value as "-log (p-value)":*/
	int rarityScorePValues_b = 0; /*if set to 1 provide rarityScore-distribution FileName; else set to 0 (and set write_scores_b = 1 to have the rarityScores written to file; this can then be used as the source of the rarityScores-distribution*/
	
	//char distDistrName[100] = "CATH_2_4";
	//char distDistrFileName[1000] = "D:\\Bioinformatics\\papers\\structural_alignment_by_GIs\\c_code\\CATH\\CATH-v2_4_Kolodny\\MatchQueriesDB_windowlgth_16_order_2_CATH_v2_4_Kolodny_CATH_v2_4_Kolodny_3wins_norm_GISAv4_5mms_13bins_distDistr.txt";
	//char distDistrName[100] = "SCOPe_ASTRAL_206_40pct";
	//char distDistrName[100] = "_top100MatchTop8000";
	//char distDistrFileName[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000/MatchQueriesDB_windowlgth_16_order_2_top100_top8000_1wins_norm_GISAv7_DB6_5mms_17bins_thrMut02.txt";
	char rarityScoreDistrName[100] = "_CATH24MatchCATH24";
	char rarityScoreDistrFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH-v2_4/RarityScan_Scores_windowslgth_16_2_order_2_CATH-v2_4_CATH-v2_4_1wins_norm_5invs_0mms_0mmsPairs_10bins2_winCovType0.txt";
	//char rarityScoreDistrFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH_test/RarityScan_Scores_windowslgth_16_2_order_2_CATH_test_CATH_test_1wins_norm_5invs_0mms_0mmsPairs_10bins2_winCovType0_CATH24MatchCATH24.txt";

	//char distDistrNamePairs[100] = "_top100MatchPairsTop8000";
	//char distDistrFileNamePairs[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top8000/MatchPairs_windowlgth_16_order_2_top100_top8000_1wins_norm_GISAv7_DB6_5mms_17bins_thrMut02.txt";
	char rarityScorePairsDistrName[100] = "_CATH24MatchCATH24";
	char rarityScorePairsDistrFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH-v2_4/RarityScan_ScoresPairs_windowslgth_16_2_order_2_CATH-v2_4_CATH-v2_4_1wins_norm_5invs_0mms_0mmsPairs_10bins2_winCovType0.txt";
	//char rarityScorePairsDistrFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/CATH_test/RarityScan_Scores_windowslgth_16_2_order_2_CATH_test_CATH_test_1wins_norm_5invs_0mms_0mmsPairs_10bins2_winCovType0_CATH24MatchCATH24.txt";

	int normalize_b = 1;

	int onlyDisjointPairs_b = 1; 

	double thresholdMutual = 4.0;

	/*For getting and writing out the invariant values for a set windowlength:*/
	int matchWindows_b = 1; /*match on single windows or not (0); will write the match results for each window to file */
	int write_matchWindows_b = 0; /*write the match results for each window to file or not (0) */
	int matchWindowPairs_b = 1; /*match on window pairs (1) or not (0); if set to 1 also set matchWindows_b to 1*/
	int write_matchWindowPairs_b = 0; /*write the match results for each window pair to file or not (0) */

	int windowCoveringType = 0;
	int windowLength = 30;
	int stepSize = 4;

	int write_windows_b = 0;
	int write_windowPairs_b = 0; /*if you use DB4 or DB5, at most one of write_windows_b and write_windowPairs_b should be set to 1; corr to matchWindows_b and matchWindowPairs_b resp.*/

	int write_scores_b = 1; //set to 1 if scores should be output
	int write_scoresPairs_b = 1; //set to 1 if scores should be output

	int writeWindowInfo_b = 0;
	
	int write_chains_b = 0;
	char pathChainsQueries[2000] = "\0";
	char fileNameOutChains[1000] =  "/RarityScan2_chains_";
	//char fileNameOutExtensionChains[20] = "_top100.txt";
	//char fileNameOutExtensionChains[20] = "_CATH_test.txt";
	char fileNameOutExtensionChains[20] = "_CATH-v2_4.txt";

	int nrOfWindowsForMatch = 1; /*just leave it at 1, has no influence apart from appearing in file names*/
	int topN = 10; // topN best matches for each window/window pair in the query will be written out

	int allowedNrOfMismatches = 0; //at most the number of invariants used in a match (more has no effect)
	int nrOfEntries = 5; //equals the number of invariants used in a match, ie nrOfInvsForMatch
	int allowedNrOfMismatchesPairs = 0;
	int nrOfEntriesForPairs = 1; /*  setting this to 1 means that only the mutual writhe is considered; >1: mutuals of higher order inv's too*/

	/*For generating a set of bins*/
	int binningType = 1;
	int nrOfBins = 10;
	double compensationForRounding = 0.0000009; //When using equi-percentile binning bin interval end points will be shifted to the right by this amount; to prevent separation of e.g. self-hits


	if (write_chains_b !=0){
		/*generate file name for writing out the data base chains:*/
		strcat(pathChainsQueries , outputPath);
		strcat(pathChainsQueries , fileNameOutChains);
		strcat(pathChainsQueries , queriesName);
		strcat(pathChainsQueries , fileNameOutExtensionChains);

		printf("Query chains will be written to: %s\n", pathChainsQueries);
		
		/*clear the file by opening it in w-mode and then closing it again*/
		fclose(fopen(pathChainsQueries , "w"));
	}

	//returnVal = computeGI(DBName, dirPath, outputPath, order, incl_abs_b, full_b, split_b,write_final_b, write_all_b, 
	//		  closed_loops_b, closedLoopLength, closedLoopDist, invValSubChainPairs_b, writeSubChainPairsAll_b, subChainLength, 
	//		  print_b, print_basic_b);

	/*printf("DBResultsFilePath:%s\n", DBResultsFilePath);

	getchar();*/


    /*returnVal = matchQueriesDB3(DBName, queriesName, queriesDirPath, DBResultsFilePath, use_scop_b, normalize_b, outputPath, maxChainLength,
					topN, chunkSize,
				   order, incl_abs_b, full_b, split_b,
				   matchWindows_b, matchWindowPairs_b, pathChainsQueries, write_chains_b,
				   print_b, print_basic_b);
*/

	if (matchWindows_b == 1 && matchWindowPairs_b == 0){

		printf("I'll be running a rarity scan, type 2, in single windows mode. Press enter to carry on."); 
		getchar();
	}
	else{
		if(matchWindows_b == 1 && matchWindowPairs_b == 1){

			printf("I'll be running a rarity scan, type 2, in single and pairs windows mode. Press enter to carry on."); 
			getchar();
		}
	}

	returnVal = rawRarity2(halted_b, rarityScoreDistrFilePath, rarityScoreDistrName, rarityScorePairsDistrFilePath, 
					rarityScorePairsDistrName, DBName, queriesName, queriesDirPath, DBResultsFilePath, 
					DBResultsPairsFilePath, pathChainsQueries,
					thresholdMutual, onlyDisjointPairs_b, 
					rarityScorePValues_b, use_scop_b, 
					use_cath_b, CATHListPath, loadFromSubDirs_b, 
					normalize_b, outputPath, maxChainLength,
					nrOfEntries, nrOfEntriesForPairs, binningType, 
					nrOfBins, compensationForRounding, allowedNrOfMismatches, 
					allowedNrOfMismatchesPairs, nrOfWindowsForMatch, chunkSize, 
					order, incl_abs_b, full_b, 
					write_matchWindows_b, matchWindowPairs_b, 
					write_matchWindowPairs_b, windowCoveringType, windowLength, 
					stepSize, write_windows_b, write_windowPairs_b, 
					write_scores_b, write_scoresPairs_b, write_chains_b,
					writeWindowInfo_b, print_b, print_basic_b);



	return returnVal;

}


/*
Function for deciding if a structure (query) stands out as "rare" compared to a back ground (based on a set of PDB
files). Rareness is decided by means of mutual writhe (possibly absolute value thereof, but signed by default). 
We often use the term "data base" for the back ground. The code allows running a set of queries against the data base; 
in fact it is set up so that one may run a series of query sets against the same back ground (withou reloading the back ground). 

The usage is:
1) Create the data base: this is done by running the main function here (SAonGISA_main_*) in the makeDB flavour. This wraps a call 
to the GI_windows function of the GISA_v*_unix-code. Must be run in "window pairs mode" (and) so that a file of mutual writhe values 
for all pairs of disjoint windows for each structure in the PDB-set is created. Set the order to 1; set windowing parameters as 
desired (e.g. windowLength 20, stepSize 4); use windowType 0.
2) Now run the rarity scan with rawRarity0: there are two methods of scoring available: 
a) Score by the maximal mutual writhe: for a given query (q) we pick out the highest absolute
mutual writhe among all the positive (negative) writhe values in the structure, or, if running in"unsigned mode" we just 
keep the pair for which the absolute value of the mutual writhe is highest (so that in this unsigned mode there is one value 
,max-abs-mutual values for each q, while in the signed case we consider the highest positive and the highest negative). 
The obtained value is then held up against the distribution of similarly obtained max mutual writhe values across the structures 
in the (background) data base (signed: two distributions, one for positive writhes and one for negative; unsigned: a single distribution
of max-abs values). This provides directly a p-value (and a score = -log10(p-value)), viz the frequency of max abs mutual values in the 
data base higher than the obtained value.
b) Score by absolute mutual writhe above a set threshold: for a given query (q) we run through all pairs of windows in q; 
for each pair for which the absolute mutual writhe (amw) is above a set threshold, T (e.g. 5), we find the probability that 
a pair in the database has an absolute mutual writhe higher than the absolute mutual writhe (amw) for the given 
pair of windows (prob(abs mutual writhe > amw)). This is just a look-up in the background distribution of amw's. The final score 
for the query is then the sum of -log(p)-values in this set of pairs

score(q) = - Sum( log(p-values))

where the sum is over all pairs in q with amw > T as just explained. If wanted, the average of this score can be had (as final score
output) -- ie 

score = - Sum( log(p-values))/#pairs in query 

The reason for using the average score is that the probability of having some (rare) 3d-configuration should increase
with the the length of the structure (here quadratically as the number of pair is quadratic as a function of the length of query). 
However, using the average score has some disadvantages, eg a short structure with one rare window pair will get a higher score than 
a longer structure with exactly the same odd pair and no other particularities; on the other hand, such short structures ought to 
appear rare if "weird 3d-configurations" are distributed uniformly over all pairs in the database. In addition, the option (a) allows getting 
the score based on the most extreme pair rather than just an average consideration (and, in addition, in the unsiged mode).

In case b, to obtain a final p-value corresponding to the obtained score, it is necessary to first run rawRarity0 
with query set = data base and with absMutValScorePValues_b = 0. This generates the background distribution of scores (and of amw's); 
if we believe that the data base consists of a representative set of structures (for a given purpose), we can with reason score any query
set against this background. So when this run for a background is done, we run rawRarity0 for the desired query set now with 
absMutValScorePValues_b = 1. In case a, p-values and corresponding scores are had directly.

The scoring method a) scans for structures having one (or more) exceptional mutual writhe pairs and allows distinguishing bewteen the 
positive writhe pairs and the negative writhe pairs; the b) method scans for structures having possibly several pairs of high 
absolute mutual writhe (above the threshold) but maybe none of an exceptional level. With b) a high treshold on the mutual writhe 
should be used, e.g. 10.


Final scores are written to a file, in decreasing order ( ie with lowest p-value at the top). 
*/
int rawRarity0( int halted_b, char DBName[100], char queriesName[100], 
				char queriesDirPath[1000], char DBResultsPairsFileName[1000], char fileNameQueryChains[2000], 
				char absMutValScoresDistrName[100], char absMutValScoresDistrFilePath[1000], int onlyDisjointPairs_b, 
				int useOnlyMaxMutVal_b, int signed_b, int absMutValScorePValues_b, 
				double thresholdMutual, int use_scop_b, int use_cath_b, 
				char *cathListPath, int loadFromSubDirs_b, int normalize_b, 
				char *outputPath, int maxChainLength, int chunkSize, 
				int windowCoveringType, int windowLength, int stepSize, 
				int getAvgScore_b, int write_windowPairs_b,  int write_rarityScoresPairs_b, 
				int write_chains_b, int writeWindowInfo_b, int print_b, 
				int print_basic_b){

	int nrOfEntriesForPairs = 1; /*MUST NOT BE CHANGED! ONLY THE MUTUAL WRITHE WILL BE USED IN THIS FCT*/
	int order = 1; /*DITTO*/
	int incl_abs_b = 0; /*DITTO*/
	int full_b = 0; /*no effect, but must be defined, ie placeholder*/
	//int normalize_b = 0; /*DITTO*/

	int use_lex = 0;

	int returnVal = 0;

	int split_b = 0;

	int cntMutAboveThreshold = 0;
	//int windowLength = 16;
	//int write_windows_b = 1; /*for writing out I-values on windows for queries*/
	//int write_windowPairs_b = 0; /*at most one of write_windows_b and write_windowPairs_b should be set to 1; corr to get_windos_b resp. get_windowPairs_b*/
	/*for matching query to DB (set of structures in the set directory:*/
	int k,l, queryNrOfWindowPairs;
	int u;
	double maxAbsMutValQuery, mutVal, absMutVal;
	double maxPosMutValQuery, maxNegMutValQuery ;
	//int nrOfWindowsForMatch = 3;
	
	double ***ptr_query_pairs, ***ptr_match_pairs; 

	//int strLengthCutOff = 10; /*chains shorter than this int will not be included in computation of invariants*/

	struct dirContent dirContent;
	int subDirCnt = 0;
	int numberOfFiles = 0;
	int subStructureCnt = 0; /*some structures are multi-mers; we compute the invariants for each sub-structure (each mono-mer)*/
	int maxNrOfWindows = 0; /*max nr of windows across loaded structures*/
	int maxNrOfWindowPairs = 0; /*max nr of window pairs across loaded structures*/
	const int numberOfPerturbations = 1; /*only there since needed when writing out the final values (preparation for perturbations code)*/
	const int pertNo = 0; /*ditto*/
	char **ptr_dirList;
	char **ptr_fileNameList;

	/*for "dynamically" generating the file names for output:*/
	char orderStr[10] = "-1"; //place holder value
	char inclAbsStr[10] = "-1"; //value is placeholder
	char nrInvsForMatchStr[10] = "-1"; //place holder value
	char windowLengthStr[10] = "-1"; //place holder value
	char windowStepSizeStr[10] = "-1"; //value is placeholder
	char windowCovTypeStr[10] = "-1"; //value is placeholder
	char topNStr[10] = "-1"; //value is placeholder
	char thresholdMutualStr[10] = "-1"; //value is placeholder

	char fileNameOutScoresPairs1[1000] = "/RarityScan0_ScoresPairs_windowslgth_";
	char fileNameOutwindowPairs1[1000] = "/RarityScan0_MutualWrithe_Pairs_windowlgth_";	
	char fileNameOut2[1000] = "_order_";

	char fileNameOut3[1000] = "\0";  //"_testQuery_top8000_5invs_norm.txt"; 
	char extensionName[50] = ".txt";

	char fileNameOutwindowPairs[2000] = "\0";
	char fileNameOutScoresPairs[2000] = "\0";
	char fileNameOutScoresPairs_pos[2000] = "\0";
	char fileNameOutScoresPairs_neg[2000] = "\0";


	char *ptr_fileNameOutwindowPairs; /*pointer to file for holding measures' values on pairs of sub-chains of pre-set windowlength*/
	FILE *ptr_fileOutwindowPairs;
	char *ptr_fileNameOutScoresPairs;
	FILE *ptr_fileOutScoresPairs;
	char *ptr_fileNameOutScoresPairs_pos;
	FILE *ptr_fileOutScoresPairs_pos;
	char *ptr_fileNameOutScoresPairs_neg;
	FILE *ptr_fileOutScoresPairs_neg;


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
	int structureNameLength = 7; /*string-length of structure name in DB, e.g SCOP, CATh*/
	char *currentStructureName;
	char *structureNameInit = "NN";
	char *currentClassId;
	char *classIdInit = "NN";
	char chainIdInit[2] = ">";
	int chainsSkippedTooShort = 0;
	int chainsSkippedTooLong = 0;
	int chainSkippedSegmTooLong = 0;
	int segTooLong_b = 0;
	//int resNr = 0; /*residue number*/
	//int resNrPrior = -999;
	int chainLen = 0; /*length of C-alpha chain*/
	int maxChainLen =0; /*max length of chain among all loaded*/
	int L = 0; /*length of segment chain = chainLen -1*/
	int simplexCnt = 0;/*size of simplex*/
	int maxSimplexCnt = 0; /*size of simplex corr to max chain length*/

	struct cAlpha  *ptr_chain = NULL; /*for holding C-alpha chain */
	struct segment *ptr_segment = NULL; /* for holding chain of segments, ith C-alpha to jth C-alpha*/
	struct segment segCoords;	 

	int m = 0, n = 0, mMin, nMin;
	int i = 0;
	int j = 0;
	int cnt = 0;
	int fileNr = 0;
	int chainNr = 0;

	/* for holding meausure values (i.e. the core of the final output); 
	content to be received from aggr function:*/
	struct I_ptr I_measures;
	double ***I_measures_raw; /*triple pointer to replace the struct I_ptr use*/ 
	int alreadyAllocatedTo_chainLen  = 0;
	int alreadyAllocatedTo_maxNrOfWindows = 0;
	int alreadyAllocatedTo_maxNrOfWindowPairs = 0;
	int alreadyAllocatedTo_subStructureCnt = 0;

	/*pointer to hold the results (that this is a ptr of a ptr is due to that the write-out can also
	handle the perturbation case; we use pertNo = 0 here):*/
	struct I_values **ptr_I_values;

	/*struct for holding invariants on windows resp. window pairs for the query and for each DB entry:*/
	struct I_windowPairs_ptr I_windowPairs; //for query
	struct I_windows_ptr I_windows_dummy; //just placeholder for call to collect_queryRawScore
	struct I_windowPairs_ptr db_I_windowPairs;
	int nrOfDBstructuresLoaded = 0; 
	int nrOfDBrecordsLoaded = 0;
	int maxNrOfWindowsDB = 0;
	int maxNrOfWindowPairsDB = 0;	
	/*pointer to hold DB results*/
	struct db_I_windowPairs_order2_ptr *ptr_dbResultWindowPairs;
	struct I_windowPairs_order2_meanStddev I_windowPairs_order2_meanStddev_DB;
	double *ptr_mutualWritheDistr;
	double *ptr_absMutualWritheDistr;
	double *ptr_maxAbsMutualWritheDistr;
	double *ptr_maxPosMutualWritheDistr;
	double *ptr_maxNegMutualWritheDistr;
	double maxAbsMutVal;
	double maxPosMutVal;
	double maxNegMutVal;
	double optionalValue = -9999;
	double optionalValue_pos = -9999;
	double optionalValue_neg = 9999;
	double *ptr_absMutValScoresDistr; //sorted list of dist-values (to be had from ...)
	FILE *ptr_fileAbsMutValScoresDistr;
	int sizeOfAbsMutValScoresDistr = -1;


	/*for searching*/
	int *intArray3; //convenience array
	int db_fileNr;//for window pairs searching
	int db_wNr1;//for window pairs searching
	int db_wNr2;//for window pairs searching
	int cntMatch;
	struct Irarity_windowPairs *ptr_Irarity_windowPairs;
	struct Irarity_windowPairs *ptr_Irarity_windowPairs_pos;
	struct Irarity_windowPairs *ptr_Irarity_windowPairs_neg;


	int cntMatch_12 = 0;
	int cntPairs = 0;
	int cntPairsConsidered = 0;
	
	double score = 0.0; 
	int cntScoredQueries = 0; 
	FILE *ptr_fileAbsMutValDistr;
	int sizeOfAbsMutValDistr = -1;
	struct queryRawScore *ptr_queryRawScore;
	struct queryRawScore *ptr_queryRawScore_pos;
	struct queryRawScore *ptr_queryRawScore_neg;
	int hitAbsMutValPair;
	double pValueAbsMutVal = 1;
	int hitAbsMutValScore;
	int hitMaxAbsMutValQuery;
	int hitMaxPosMutValQuery;
	int hitMaxNegMutValQuery;
	double pValueQuery = -1;
	char * qWindowsString;
	char * qWindowsStringAid;
	char * rarityScoreString;


	/*for closed loops finding, but only place holders here:*/
	int closed_loops_b = 0; //must be kept like this!
	struct twoSegmentIndex *ptr_closedLoopInd;
	int initValNrOfClosedLoops = 0;
	int nrOfClosedLoops = 0;
	int cntClosedLoop = 0;

	FILE *fp; //dummy file ptr used when checking file names

	/*params to control the execution*/
	int keepGoing_b = 1; 
	char keepGoingStr[10] = "-1"; //place holder value
	int killroyWasHere_b = 0;
	int dbLoaded_b = 0;
	char dataTypeIndicatorStr[3] = "-1"; //place holder value
	int dataTypeIndicator = -1;
	int cntRuns = 0; //only for testing

	/*Declarations stop here*/


	/*The code for running one set of queries against the data base is wrapped in an outer while-loop, which allows the user to
	enter a new queris set, without having to reload the data base*/
	while(keepGoing_b == 1){

		/*********************************************************************************************************/
		/*Generation of output file names:*/
		/*********************************************************************************************************/

		/*we first build a string (fileNameOut3) shared by most file names, which contains various params values;
		fileNameOut3 will be concatenated into the file name strings below:*/

		if(normalize_b == 1){
			strcat(fileNameOut3, "norm_");
		}
		else{
			strcat(fileNameOut3, "notnorm_");
		}
		sprintf(nrInvsForMatchStr, "%d", nrOfEntriesForPairs);
		strcat(fileNameOut3, nrInvsForMatchStr);
		strcat(fileNameOut3, "invs_");
		strcat(fileNameOut3, "_winCovType");
		sprintf(windowCovTypeStr, "%d", windowCoveringType);
		strcat(fileNameOut3, windowCovTypeStr);
		strcat(fileNameOut3, "_threshMut");
		sprintf(thresholdMutualStr, "%0.2f", thresholdMutual);
		strcat(fileNameOut3, thresholdMutualStr);


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
			strcat(fileNameOutwindowPairs, queriesName);
			strcat(fileNameOutwindowPairs, "_");
			strcat(fileNameOutwindowPairs, DBName);
			strcat(fileNameOutwindowPairs, "_");
			strcat(fileNameOutwindowPairs, fileNameOut3);
			strcat(fileNameOutwindowPairs, extensionName);
			printf("Results for window pairs will be written to: %s\n", fileNameOutwindowPairs);

			 
			ptr_fileNameOutwindowPairs = fileNameOutwindowPairs;
			
			/*clear the file by opening it in w-mode and then closing it again*/
			fp= fopen(ptr_fileNameOutwindowPairs, "w");
			if(fp == NULL){
				printf("Sorry, file name %s is invalid. Could the path be wrong? A blank space?", ptr_fileNameOutwindowPairs);
				getchar();
				exit(EXIT_FAILURE);
			}
			else{
				fclose(fp);
			}
		

		}



		if (write_rarityScoresPairs_b !=0){

			/*generate file name for writing out match values*/
			sprintf(orderStr, "%d", order);
			strcat(fileNameOutScoresPairs, outputPath);
			strcat(fileNameOutScoresPairs, fileNameOutScoresPairs1);
			sprintf(windowLengthStr, "%d", windowLength);
			strcat(fileNameOutScoresPairs, windowLengthStr);
			strcat(fileNameOutScoresPairs, "_");
			sprintf(windowStepSizeStr, "%d", stepSize);
			strcat(fileNameOutScoresPairs, windowStepSizeStr);
			strcat(fileNameOutScoresPairs, fileNameOut2);
			strcat(fileNameOutScoresPairs, orderStr);
			strcat(fileNameOutScoresPairs, "_");
			sprintf(inclAbsStr, "%d", incl_abs_b);
			strcat(fileNameOutScoresPairs, inclAbsStr);
			strcat(fileNameOutScoresPairs, "_");
			strcat(fileNameOutScoresPairs, queriesName);
			strcat(fileNameOutScoresPairs, "_");
			strcat(fileNameOutScoresPairs, DBName);
			strcat(fileNameOutScoresPairs, "_");
			strcat(fileNameOutScoresPairs, fileNameOut3);
	
			/*when getting the rarity score p-values, include info on which distribution was used:*/
			if(useOnlyMaxMutVal_b == 0){
				strcat(fileNameOutScoresPairs, "_scoreByAbsMutualWrithe");

				/*when getting the scores' p-values, include info on which distribution was used:*/
				if(absMutValScorePValues_b == 1){
					strcat(fileNameOutScoresPairs, absMutValScoresDistrName);
				}
			}
			else if(useOnlyMaxMutVal_b == 1){

				if(signed_b == 0){
					strcat(fileNameOutScoresPairs, "_scoreByMaxAbsMutualWrithe");
				}
				else{ //when signed, we need two files for output:

					strcpy(fileNameOutScoresPairs_pos, fileNameOutScoresPairs);
					strcat(fileNameOutScoresPairs_pos, "_scoreByMaxPosMutualWrithe");
					strcpy(fileNameOutScoresPairs_neg, fileNameOutScoresPairs);
					strcat(fileNameOutScoresPairs_neg, "_scoreByMaxNegMutualWrithe");

				}
			}
			
			if(useOnlyMaxMutVal_b == 0||signed_b == 0){
				strcat(fileNameOutScoresPairs, extensionName);
				printf("Rarity scores will be written to: %s\n", fileNameOutScoresPairs);

				ptr_fileNameOutScoresPairs = fileNameOutScoresPairs;
			
				/*clear the file by opening it in w-mode and then closing it again*/
				ptr_fileOutScoresPairs = fopen(ptr_fileNameOutScoresPairs, "w");

				fprintf(ptr_fileOutScoresPairs, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;\n",
					"DB name",
					"Query file",
					"structureName",
					"chainId",
					"classId",
					"chainLen",
					"nrOfWindows",
					"windowInfo",
					"rarityScore",
					"pValue",
					"optionalvalue"
					);  		

			}
			else{
				strcat(fileNameOutScoresPairs_pos, extensionName);
				strcat(fileNameOutScoresPairs_neg, extensionName);
				printf("Rarity scores for the positive writhes will be written to: %s\n", fileNameOutScoresPairs_pos);
				printf("Rarity scores for the negative writhes will be written to: %s\n", fileNameOutScoresPairs_neg);

				ptr_fileNameOutScoresPairs_pos = fileNameOutScoresPairs_pos;
				ptr_fileNameOutScoresPairs_neg = fileNameOutScoresPairs_neg;
			
				/*clear the files by opening them in w-mode and then closing it again*/
				ptr_fileOutScoresPairs_pos = fopen(ptr_fileNameOutScoresPairs_pos, "w");
				ptr_fileOutScoresPairs_neg = fopen(ptr_fileNameOutScoresPairs_neg, "w");

				fprintf(ptr_fileOutScoresPairs_pos, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;\n",
					"DB name",
					"Query file",
					"structureName",
					"chainId",
					"classId",
					"chainLen",
					"nrOfWindows",
					"windowInfo",
					"rarityScore",
					"pValue",
					"optionalvalue"
					);  	

				fprintf(ptr_fileOutScoresPairs_neg, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;\n",
					"DB name",
					"Query file",
					"structureName",
					"chainId",
					"classId",
					"chainLen",
					"nrOfWindows",
					"windowInfo",
					"rarityScore",
					"pValue",
					"optionalvalue"
					);  		

			
			}


			//printf("If fine, press enter to continue.\n");
			//getchar();


		}

		if(halted_b == 1){
			printf("If the output paths are fine, press enter to continue.\n");
			getchar();
		}

		/*********************************************************************************************************/
		/*Get content of the desired directory: number of files and a list of the file names:*/
		/*********************************************************************************************************/

		if(loadFromSubDirs_b ==0){
			dirContent = ListDirectoryContents(queriesDirPath);
		}
		else{
			returnVal = ListDirectoryContents2(queriesDirPath, &dirContent, &subDirCnt);
		}
		numberOfFiles = dirContent.numberOfFiles; 
		ptr_dirList = dirContent.ptr_dirList;
		ptr_fileNameList = dirContent.ptr_fileNameList;

		printf("Number of query files in directory:%d \n",numberOfFiles);
		printf("1st file: %s\n",ptr_fileNameList[0] );
		
		//printf("Press enter to continue\n");
		//getchar();

		startComplete = clock();


		/****************************************************************************************************************/
		/*loop through the list of filenames in the directory to find chain info, eg the max chain length (for setting
		the size in the (global) mem alloc to I_measures's ptr's). First though we make a bulk allocation
		to the chain-info keeping pointer, using the pre-set maxNrOfChains:*/
		/****************************************************************************************************************/

		if(killroyWasHere_b == 0){
			returnVal = alloc_init_chainInStr(&chainInStr, maxNrOfChains);
		}


		/*if CATH data, read in CATH domain list. Note: we abuse the struct chainInStructure to accomodate
		domains rather than chains -- in CATH each chain may be split in several domains*/ 
		if(use_cath_b == 1){

			if(convToCLF20_b == 1){
				structureNameLength = lengthCLFformat20;
			}
			else{
				structureNameLength = 6;
			}

			/*first allocate a block of memory (if not large enough it will be extended when exectuing readCATHDomainList)*/
			ptr_cathDomain = (struct chainInfo *) malloc (maxNrOfDomains*sizeof(struct chainInfo));
			/*init*/
			for(i = 0; i< maxNrOfDomains; i++){
				ptr_cathDomain[i].chainId = (char *) malloc(sizeof(char)*2);
				ptr_cathDomain[i].structureName = (char *) malloc(sizeof(char)*(structureNameLength+2)); //to hold domain name
				ptr_cathDomain[i].classId = (char *) malloc(sizeof(char)*10);
				ptr_cathDomain[i].classIdConv = (int *) malloc(classifierSize*sizeof(int));

				strcpy(ptr_cathDomain[i].chainId, chainIdInit); /*init to unlikely value*/
				strcpy(ptr_cathDomain[i].structureName, structureNameInit);
				strcpy(ptr_cathDomain[i].classId, classIdInit);
				/*printf("Chain Id initial value:%c\n", *ptr_cathDomain[i].chainId);
				printf("Class Id initial value:%c\n", *ptr_cathDomain[i].classId);*/
				ptr_cathDomain[i].chainLength = -1;

				for(j=0; j< classifierSize; j++){ptr_cathDomain[i].classIdConv[j] = -1;}

			}

			/*open CATH file list:*/
			ptr_cath_domain_list_file = fopen(cathListPath, "r");

			if (!ptr_cath_domain_list_file){
				printf("Sorry, the file containing the positive-list of CATH domains was not found (%s was provided). Enter a valid file name and run again.\n", cathListPath);
				return 0;
			}

			/*read in contents of CATH domain list file*/
			nrOfCATHDomainsInList = readCATHDomainList(ptr_cath_domain_list_file, &ptr_cathDomain, convToCLF20_b );

			printf("nr of cath doms: %d\n", nrOfCATHDomainsInList);

			/*lex-sort the ptr_cath_domain_list_file on the domain name entry. Below we
			need to look up each file for which we have PDB-coords in this list. In other
			words, this list serves as a "positive list" of domains, for which to compute
			the invariants.*/
			heapSortCATHDomain(ptr_cathDomain, nrOfCATHDomainsInList, structureNameLength,0);

			/*for(i=0;i <nrOfCATHDomainsInList;i++){
				if(i%100 == 0){
					for(j=i;j <i+100;j++){
						printf("dom %d: %s\n", j, ptr_cathDomain[j].structureName);
					}
					getchar();
				}
			}*/
		}

		printf("Fetching chain info for data base structures ...\n");
		for(fileNr = 0; fileNr <numberOfFiles; fileNr ++){

			if (fileNr%100 == 0){printf(".. now fetched the chain info for %d structures\n", fileNr);}

			if (print_b == 1){
				printf("File name:%s path:%s fileNr:%d\n",ptr_fileNameList[fileNr], ptr_dirList[fileNr], fileNr);
			}

			/*if using CATH data: check if the file, ptr_fileNameList[fileNr], is in the
			domain list specified (and loaded and lex-sorted above)*/
			if(use_cath_b == 1){
				hitIndexDomain = bisectionalSearchCATHDomain(ptr_fileNameList[fileNr], ptr_cathDomain, nrOfCATHDomainsInList, lengthCLFformat20);
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
				continue;
			}

			/*get number of chains, their id's and lengths for this file:*/
			//chainInStr = readPDBChainStructure(ptr_fileIn);
			if(use_scop_b ==1){
				readPDBChainStructureSCOP(ptr_fileIn, &chainInStr );
			}
			else if(use_cath_b ==1){
				readPDBDomainStructureCATH(ptr_fileIn, &chainInStr );//will read chain length and chain id into chainInStr
			}
			else{
				readPDBChainStructure2(ptr_fileIn, &chainInStr );
			}


			fclose(ptr_fileIn);  


			//printf("chainInStr.numberOfChains: %d\n", chainInStr.numberOfChains);
			for (i=0;i <= chainInStr.numberOfChains -1;i++){
				//printf("chainInStr.ptr_chainInStr[i].chainLength: %d\n", chainInStr.ptr_chainInStr[i].chainLength);
				if(chainInStr.ptr_chainInStr[i].chainLength <= maxChainLength){ //we disregard very long chains for RAM availability reasons
					maxChainLen = max(chainInStr.ptr_chainInStr[i].chainLength,maxChainLen);
					if(chainInStr.ptr_chainInStr[i].chainLength >= strLengthCutOff){
										
						/*if using CATH data: we checked inside the fileNr-loop if the file, ptr_fileNameList[fileNr], is in the
						domain list specified (and loaded and lex-sorted above). We need another check to avoid loading in 
						several chains for the structureName under consideration:*/
						if(use_cath_b == 1){ //if not using CATH, chainInStr.ptr_chainInStr[i].classId was read in in the pdb-file reads above
							
							/*In case the pdb-file contains several chains: the domain is attached to a single chain (is a subset of 
							a single chain) and we want only to consider that chain. we therefore skip "wrong chains". Obs: hitIndexDomain
							was found above:*/
							if(chainInStr.numberOfChains > 1 && ptr_cathDomain[hitIndexDomain].structureName[4] != chainInStr.ptr_chainInStr[i].chainId[0]){

								if(print_basic_b > 0){
									printf("CATH domain: %s is regarded not to sit in this chain: %s\n", ptr_cathDomain[hitIndexDomain].structureName, chainInStr.ptr_chainInStr[i].chainId);
								}
								continue;
							}

						}

						subStructureCnt += 1; /*count and index the sub-structures for which the invariants will be computed*/
					
					}
				}
			}

		}
		printf("... done.\n");

		L = maxChainLen -1;

		printf("max chain len:%d\n", maxChainLen);
		printf("sub structure count:%d\n", subStructureCnt);
		if(halted_b == 1){
			getchar();
		}

		/*********************************************************************************************************/
		/*MEMORY ALLOCATION*/
		/*********************************************************************************************************/

		/*We allocate memory for ptr holding the I-values and a ptr for holding the 3d-coord's, both for the longest chain*/
		printf("alreadyAllocatedTo_chainLen %d\n", alreadyAllocatedTo_chainLen);
		returnVal = alloc_init_I_measures(&I_measures, order, full_b, chainNr, maxChainLen, closed_loops_b, alreadyAllocatedTo_chainLen, &ptr_closedLoopInd);

		if(alreadyAllocatedTo_chainLen == 0){ //no allocation was done yet

			ptr_chain = (struct cAlpha *) malloc (maxChainLen*sizeof(struct cAlpha));
			/*init*/
			for ( n=0; n <= maxChainLen -1 ; n++ ){
				
				ptr_chain[n].coords.x = 0.0;
				ptr_chain[n].coords.y = 0.0;
				ptr_chain[n].coords.z = 0.0;

				ptr_chain[n].residueNr = -123;

			}

			ptr_segment = (struct  segment *) malloc ((maxChainLen-1)*sizeof(struct segment)); 
			/*init*/
			segCoords.s1.x = 0.0;
			segCoords.s1.y = 0.0;
			segCoords.s1.z = 0.0;
			segCoords.s2.x = 0.0;
			segCoords.s2.y = 0.0;
			segCoords.s2.z = 0.0;
			for ( n=0; n <= maxChainLen -2 ; n++ ){
					segCoords.s1 = ptr_chain[n].coords;
					segCoords.s2 = ptr_chain[n+1].coords;

					ptr_segment[n] = segCoords;

			}

		}
		else if(maxChainLen > alreadyAllocatedTo_chainLen){ //reallocation necc

			ptr_chain = realloc (ptr_chain, maxChainLen*sizeof(struct cAlpha));

			/*(re)init*/
			for ( n=0; n <= maxChainLen -1 ; n++ ){
				
				ptr_chain[n].coords.x = 0.0;
				ptr_chain[n].coords.y = 0.0;
				ptr_chain[n].coords.z = 0.0;

				ptr_chain[n].residueNr = -123;
			}

			ptr_segment = realloc (ptr_segment, (maxChainLen-1)*sizeof(struct segment));

			/*init*/
			segCoords.s1.x = 0.0;
			segCoords.s1.y = 0.0;
			segCoords.s1.z = 0.0;
			segCoords.s2.x = 0.0;
			segCoords.s2.y = 0.0;
			segCoords.s2.z = 0.0;
			for ( n=0; n <= maxChainLen -2 ; n++ ){

					segCoords.s1 = ptr_chain[n].coords;
					segCoords.s2 = ptr_chain[n+1].coords;

					ptr_segment[n] = segCoords;

			}


		}
		alreadyAllocatedTo_chainLen = max(alreadyAllocatedTo_chainLen, maxChainLen);

		/****************************************************************************************************************/
		/* Allocate mem to some ptr's for query: to hold I-values on window pairs and for the scores	
		/****************************************************************************************************************/

		maxNrOfWindows =  (int) ceil((double) 1.1*(maxChainLen-windowLength)/stepSize); /*this should just suffice for the longest chain -- we may be using overlapping windows*/

		maxNrOfWindowPairs = alloc_init_I_windowPairs_ptr(&I_windowPairs, order, full_b, maxChainLen, windowLength, stepSize, maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows);

		///*allocate and initialize a ptr to hold the scores results for each query*/
		printf("subStructureCnt %d  alreadyAllocatedTo_subStructureCnt %d\n", subStructureCnt,alreadyAllocatedTo_subStructureCnt );
		printf("maxWindowInfoSize*maxNrOfWindowPairs %d maxWindowInfoSize*alreadyAllocatedTo_maxNrOfWindowPairs %d\n", maxWindowInfoSize*maxNrOfWindowPairs, maxWindowInfoSize*alreadyAllocatedTo_maxNrOfWindowPairs);

		/*for holding rarity scores to be written out:*/
		if(useOnlyMaxMutVal_b == 1 && signed_b == 1){

			alloc_init_queryRawScore(&ptr_queryRawScore_pos, subStructureCnt, alreadyAllocatedTo_subStructureCnt, maxWindowInfoSize*maxNrOfWindowPairs, maxWindowInfoSize*alreadyAllocatedTo_maxNrOfWindowPairs ); /*maxWindowInfoSize = 50 ... should suffice!*/
			alloc_init_queryRawScore(&ptr_queryRawScore_neg, subStructureCnt, alreadyAllocatedTo_subStructureCnt, maxWindowInfoSize*maxNrOfWindowPairs, maxWindowInfoSize*alreadyAllocatedTo_maxNrOfWindowPairs ); /*maxWindowInfoSize = 50 ... should suffice!*/

			returnVal = alloc_init_ptr_Irarity_windowPairs(&ptr_Irarity_windowPairs_pos, order, full_b, maxChainLen, windowLength, stepSize, maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows);
			returnVal = alloc_init_ptr_Irarity_windowPairs(&ptr_Irarity_windowPairs_neg, order, full_b, maxChainLen, windowLength, stepSize, maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows);
		}
		else{

			alloc_init_queryRawScore(&ptr_queryRawScore, subStructureCnt, alreadyAllocatedTo_subStructureCnt, maxWindowInfoSize*maxNrOfWindowPairs, maxWindowInfoSize*alreadyAllocatedTo_maxNrOfWindowPairs ); /*maxWindowInfoSize = 50 ... should suffice!*/

			returnVal = alloc_init_ptr_Irarity_windowPairs(&ptr_Irarity_windowPairs, order, full_b, maxChainLen, windowLength, stepSize, maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows);
		}

		//update:
		alreadyAllocatedTo_maxNrOfWindows = max(alreadyAllocatedTo_maxNrOfWindows, maxNrOfWindows);
		alreadyAllocatedTo_maxNrOfWindowPairs = max(alreadyAllocatedTo_maxNrOfWindowPairs, maxNrOfWindowPairs);
		alreadyAllocatedTo_subStructureCnt = max(alreadyAllocatedTo_subStructureCnt, subStructureCnt);

		/***************************************************/
		/*Load DB results for pairs*/
		/***************************************************/

		if(dbLoaded_b == 0){

			/*Load DB for the pair I-values:*/
			intArray3 = readDBresultsWindowPairsToPtr(&ptr_dbResultWindowPairs, DBResultsPairsFileName, chunkSize, &maxNrOfWindowsDB);

			nrOfDBstructuresLoaded = intArray3[0];
			nrOfDBrecordsLoaded =  intArray3[1]; //equals nr of window pairs loaded
				
			printf("nrOfDBstructuresLoaded %d nrOfDBrecordsLoaded %d\n", nrOfDBstructuresLoaded, nrOfDBrecordsLoaded);
			//getchar();


			/*normalize the DB results if desired (it usually is NOT for the use of this fct)*/
			if(normalize_b == 1){

				printf("normalizing DB result (pairs data) ..\n");


				I_windowPairs_order2_meanStddev_DB = normalizeDBresultsWindowPairsPtr(ptr_dbResultWindowPairs, nrOfDBstructuresLoaded, nrOfDBrecordsLoaded);

				printf("Pairs normalization: I12 move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_I12, I_windowPairs_order2_meanStddev_DB.stddev_I12);


				printf(".. done \n");
				if(halted_b == 1){
					getchar();
				}

			}


			/*Next read the mutual writhe values into a 1d-array and sort it.
			We also load a similar array containing the absolute values of the 
			writhe numbers, and a particular array containing for each structure
			the max absolute value of the writhe (across all its window pairs):*/
			ptr_mutualWritheDistr = (double *) malloc(nrOfDBrecordsLoaded*sizeof(double));
			ptr_absMutualWritheDistr = (double *) malloc(nrOfDBrecordsLoaded*sizeof(double));
			if(useOnlyMaxMutVal_b == 1){
				if(signed_b == 0){
					ptr_maxAbsMutualWritheDistr = (double *) malloc(nrOfDBstructuresLoaded*sizeof(double));
				}
				else{
					ptr_maxPosMutualWritheDistr = (double *) malloc(nrOfDBstructuresLoaded*sizeof(double));
					ptr_maxNegMutualWritheDistr = (double *) malloc(nrOfDBstructuresLoaded*sizeof(double));
				}
			}
			cnt = 0;
			for(n=0;n<nrOfDBstructuresLoaded; n++){

				if(useOnlyMaxMutVal_b == 1){
					if(signed_b == 0){
						maxAbsMutVal = 0.0;
						ptr_maxAbsMutualWritheDistr[n] = 0.0;
					}
					else{
						maxPosMutVal = 0.0;
						ptr_maxPosMutualWritheDistr[n] = 0.0;
						maxNegMutVal = 0.0;
						ptr_maxNegMutualWritheDistr[n] = 0.0;
					}
				}


				for(k=0; k < ptr_dbResultWindowPairs[n].nrOfWindows; k++){


					/*if(ptr_dbResultWindowPairs[n].windowNr_1[k] == -1){
							continue;
					}*/

					for(l=k+1; l < ptr_dbResultWindowPairs[n].nrOfWindows; l++){

						if(ptr_dbResultWindowPairs[n].windowPair[k][l].windowNr_1 == -1||ptr_dbResultWindowPairs[n].windowPair[k][l].windowNr_2 == -1){
							continue;
						}
						if(onlyDisjointPairs_b == 1 && ptr_dbResultWindowPairs[n].windowPair[k][l].segIndices_1[1] >  ptr_dbResultWindowPairs[n].windowPair[k][l].segIndices_2[0]){

							continue;

						}
			
						
						ptr_mutualWritheDistr[cnt] = ptr_dbResultWindowPairs[n].I12[k][l];
						ptr_absMutualWritheDistr[cnt] = fabs(ptr_dbResultWindowPairs[n].I12[k][l]);

						if(useOnlyMaxMutVal_b == 1){
							if(signed_b == 0){	
								if(ptr_absMutualWritheDistr[cnt] > maxAbsMutVal){
									maxAbsMutVal = ptr_absMutualWritheDistr[cnt];
									//printf("Langmodigt venter bysvalen stadig stædigt længes ...maxAbsMutVal:%lf\n", maxAbsMutVal);
									ptr_maxAbsMutualWritheDistr[n] = maxAbsMutVal;
								}
							}
							else{

								if(ptr_mutualWritheDistr[cnt] > 0 && ptr_absMutualWritheDistr[cnt] > maxPosMutVal){
									maxPosMutVal = ptr_absMutualWritheDistr[cnt];
									ptr_maxPosMutualWritheDistr[n] = maxPosMutVal;
								}
								else if(ptr_mutualWritheDistr[cnt] < 0 && ptr_absMutualWritheDistr[cnt] > maxNegMutVal){
									maxNegMutVal = ptr_absMutualWritheDistr[cnt];
									ptr_maxNegMutualWritheDistr[n] = maxNegMutVal; //we use the abs-value here! -- this will allow us to use the same bisectional-search in the neg case as in the pos
								}

							}
						}
						
						cnt += 1;
					}
				}
			}

			printf("Obs: mutual writhe distr loaded with %d values while nrOfDBrecordsLoaded was: %d\n", cnt, nrOfDBrecordsLoaded);
			printf("Obs: max mutual writhe distr loaded with %d values while nrOfDBstructuresLoaded was: %d\n", n, nrOfDBstructuresLoaded);
			if(cnt != nrOfDBrecordsLoaded){
				printf("Warning: mutual writhe distr loaded with %d values while nrOfDBrecordsLoaded was: %d\n", cnt, nrOfDBrecordsLoaded);
			}

			printf("sorting DB results (pairs data) ..\n");

			//for(i=0;i<10; i++){printf("ptr_absMutualWritheDistr[%d]: %lf\n", i,  ptr_absMutualWritheDistr[i]);}
			//getchar();
			if(halted_b == 1){
				if(useOnlyMaxMutVal_b == 1){
					
					printf("Before sorting: \n");

					if(signed_b == 0){	
						for(i=0;i<10; i++){printf("ptr_maxAbsMutualWritheDistr[%d]: %lf\n", i,  ptr_maxAbsMutualWritheDistr[i]);}
					}
					else{
						for(i=0;i<10; i++){printf("ptr_maxPosMutualWritheDistr[%d]: %lf\n", n-1-i,  ptr_maxPosMutualWritheDistr[n-1-i]);}
						for(i=0;i<10; i++){printf("ptr_maxNegMutualWritheDistr[%d]: %lf\n", n-1-i,  ptr_maxNegMutualWritheDistr[n-1-i]);}
					}
				
				}
			}

			/*sort in increasing order (obs: if changing to decreasing you must change the bisectional search below too!)*/
			if(useOnlyMaxMutVal_b == 0){
				heapSort1dArray(ptr_mutualWritheDistr, cnt, 0);
				heapSort1dArray(ptr_absMutualWritheDistr, cnt, 0);
			}
			else if(useOnlyMaxMutVal_b == 1){
				if(signed_b == 0){	
					heapSort1dArray(ptr_maxAbsMutualWritheDistr, n, 0);
				}
				else{
					heapSort1dArray(ptr_maxPosMutualWritheDistr, n, 0);
					heapSort1dArray(ptr_maxNegMutualWritheDistr, n, 0);
				}
			}

			printf(".. done \n");

			//for(i=0;i<10; i++){printf("ptr_absMutualWritheDistr[%d]: %lf\n", cnt-1-i,  ptr_absMutualWritheDistr[cnt-1-i]);}
			if(halted_b == 1){
				if(useOnlyMaxMutVal_b == 1){

					printf("After sorting: \n");

					if(signed_b == 0){	
						for(i=0;i<10; i++){printf("ptr_maxAbsMutualWritheDistr[%d]: %lf\n", n-1-i,  ptr_maxAbsMutualWritheDistr[n-1-i]);}
					}
					else{
						for(i=0;i<10; i++){printf("ptr_maxPosMutualWritheDistr[%d]: %lf\n", n-1-i,  ptr_maxPosMutualWritheDistr[n-1-i]);}
						for(i=0;i<10; i++){printf("ptr_maxNegMutualWritheDistr[%d]: %lf\n", n-1-i,  ptr_maxNegMutualWritheDistr[n-1-i]);}
					}
				
				}
			}

			//getchar();

			/****************************************************************************************************************/
			/* Read in rarityScore-data and generate the distribution of these (just a sort):*/
			/****************************************************************************************************************/
			/*Load a scores distribution (abs mut val based) if needed*/  
			if(absMutValScorePValues_b == 1){
				
						ptr_fileAbsMutValScoresDistr = fopen(absMutValScoresDistrFilePath, "r");

						if (!ptr_fileAbsMutValScoresDistr){
							printf("Sorry, the file containing the background scores, absMutValScoresDistrFilePath, was not found (%s was provided). Enter a valid file name and run again.\n", absMutValScoresDistrFilePath);
							return 0;
						}


						/*this will provide a sorted list (smallest to highest), ptr_distDistr, of dist-values.
						Will be used to score each window hit (as -log(p-value), where the p-value is obtained by a
						look up in this array: */
						sizeOfAbsMutValScoresDistr = readSortRarityScoreData(ptr_fileAbsMutValScoresDistr, &ptr_absMutValScoresDistr);

						printf("sizeOfAbsMutValScoresDistr %d\n", sizeOfAbsMutValScoresDistr);
						if(halted_b == 1){getchar();}

						for(i=0;i<min(500, sizeOfAbsMutValScoresDistr);i++){
							printf("abs mut val based score at %d after sorting: %lf\n", i, ptr_absMutValScoresDistr[i]);
							//if(i%100 == 0){getchar();}
						}
						if(halted_b == 1){getchar();}


				}


			/*Finally a couple of allocations to strings for read-out of scores:*/
			if(write_rarityScoresPairs_b == 1){
				qWindowsStringAid = (char *) malloc(5*sizeof(char));
				qWindowsString = (char *) malloc(20*maxNrOfWindowPairs*sizeof(char));
				rarityScoreString = (char *) malloc(10*sizeof(char));

				/*init*/
				strcpy(qWindowsStringAid, "\0");
				strcpy(qWindowsString, "\0");
				strcpy(rarityScoreString, "\0");
			}

			/*Change indicator of db-load:*/
			dbLoaded_b = 1;

		}
		
		subStructureCnt = 0; /*reset*/

		/****************************************************************************************************************/
		/* MAIN LOOP */
		/****************************************************************************************************************/
		//startComplete = clock();

		/*loop through the list of filenames in the queries' directory and compute the desired measures all along:*/
		for (fileNr = 0; fileNr <numberOfFiles; fileNr ++){

			if (fileNr%10 == 0){printf("Now at file nr: %d name: %s\n", fileNr, ptr_fileNameList[fileNr]);}

			//printf("Query: %s\n", ptr_dirList[fileNr]);
			//getchar();

			if (print_basic_b == 1){
				printf("Now handling file:%s\n fileNr:%d\n",ptr_dirList[fileNr], fileNr);
			}

			if(use_cath_b == 1){
				hitIndexDomain = bisectionalSearchCATHDomain(ptr_fileNameList[fileNr], ptr_cathDomain, nrOfCATHDomainsInList, lengthCLFformat20);
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

			strcpy(I_measures.fileName, ptr_dirList[fileNr]);

			if (!ptr_fileIn){
				printf("Sorry, the file %s could not be found\n", ptr_dirList[fileNr]);
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
				readPDBDomainStructureCATH(ptr_fileIn, &chainInStr );//will read chain length and chain id into chainInStr; other values (structureName etc) will be read in below
				//printf("Sommeren mild ... ... help .. file: %s chainInStr.numberOfChains: %d \n",I_measures.fileName, chainInStr.numberOfChains);
				//getchar();
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
				}
			}

			/****************************************************************************************************************/
			/* INNER LOOP */
			/****************************************************************************************************************/
			/*loop through the chains in the current structure/file ... and compute ... :*/

			chainNr = 0; /*reset*/
			for (chainNr = 0; chainNr <= chainInStr.numberOfChains -1 ; chainNr ++){

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
							printf("CATH domain: %s is regarded not to sit in this chain: %s\n", ptr_cathDomain[hitIndexDomain].structureName, chainInStr.ptr_chainInStr[chainNr].chainId);
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

				strcpy(I_measures.structureName, chainInStr.ptr_chainInStr[chainNr].structureName);
				strcpy(I_measures.classId, chainInStr.ptr_chainInStr[chainNr].classId);
				strcpy(I_measures.chainId, chainInStr.ptr_chainInStr[chainNr].chainId);

				I_measures.chainNr = &chainNr;

				chainLen = chainInStr.ptr_chainInStr[chainNr].chainLength;

				if (print_basic_b ==1){
					printf("sub str:%s\n",I_measures.structureName);
					printf("sub str cnt:%d\n",subStructureCnt);
					printf("Chain id:%s\n", chainInStr.ptr_chainInStr[chainNr].chainId);
					printf("Chain nr:%d\n", chainNr);
					printf("Chain Length: %d\n", chainLen);
				}

				/*if the chain length is less than strLengthCutOff we skip the computation*/
				if(chainLen < strLengthCutOff){
					//if (print_basic_b ==1){
					printf("Chain %s of %s skipped since length is < %d\n", I_measures.chainId, I_measures.fileName,  strLengthCutOff);
					//getchar();
					//}
					chainsSkippedTooShort += 1;
					continue;
				}

				
				if(use_scop_b ==1){
					if(strstr(chainInStr.ptr_chainInStr[chainNr].classId, "i") == chainInStr.ptr_chainInStr[chainNr].classId || 
						strstr(chainInStr.ptr_chainInStr[chainNr].classId, "j") == chainInStr.ptr_chainInStr[chainNr].classId ||
						strstr(chainInStr.ptr_chainInStr[chainNr].classId, "k") == chainInStr.ptr_chainInStr[chainNr].classId ||
						strstr(chainInStr.ptr_chainInStr[chainNr].classId, "h") == chainInStr.ptr_chainInStr[chainNr].classId){
						printf("Class %s of chain of %s was skipped since it is h, i, j or k\n", chainInStr.ptr_chainInStr[chainNr].chainId, ptr_dirList[fileNr], chainInStr.ptr_chainInStr[chainNr].classId);
						subStructureCnt += 1; //!
						continue;
					}
				}

				I_measures.chainLen = &chainLen;

				rewind(ptr_fileIn);
				returnVal = main_readPDB2(ptr_fileIn, ptr_chain, chainNr, chainLen);

				/*length of segments chain:*/
				L = chainLen - 1;
				/*allocate memory for array of segments*/
				//ptr_segment = (struct  segment *) calloc (L, sizeof(struct segment));
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
					printf("Chain: %s in structure: %s removed due to too long %d 'th segment\n", chainInStr.ptr_chainInStr[chainNr].chainId, I_measures.fileName,  segTooLong_b-1);
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
					//returnVal = aggrAndW_raw(ptr_segment, chainLen, order, full_b, &I_measures_raw);
					}
				if (incl_abs_b == 0 && split_b == 0 && closed_loops_b == 0){
					returnVal = aggrAndW_ExAbs(ptr_segment, chainLen, order, full_b, I_measures);
					}
				if (incl_abs_b == 1 && split_b == 1 && closed_loops_b == 0){
					returnVal = wAll(ptr_segment, chainLen, I_measures.wVal);
					returnVal = aggr(chainLen, order, full_b, I_measures.wVal, I_measures);
				}

				endComp = clock();

				compTime = ((double) ((endComp - startComp))/CLOCKS_PER_SEC);
				//compTime = ((double) (endComp - startComp));
				if(print_b == 1){
					printf("Comp time: %f\n",compTime);
					//printf("CPU time spend for measures on this file (only computation/aggregation):%lf \n", compTime);
				}
				if (compTime > max_compTime){
					max_compTime = compTime;
				}
				compTimeComplete += compTime;

				/*Code for searching through the DB starts here*/
				/**************************************************************************************/
				/*Scanning/matching based on window pairs*/
				/**************************************************************************************/

				if (write_chains_b != 0){
					returnVal = writeChain(fileNameQueryChains, I_measures.fileName, I_measures.structureName, I_measures.chainId, ptr_chain, *I_measures.chainLen);
				}

				returnVal = getInvariantsOnWindowPairs(&I_windowPairs, order, windowCoveringType, windowLength, stepSize, L, I_measures); 

				/*normalize the pairs' I-values if desired. The single window values (I_windows) and the window
				pair values (I_windowPairs) are normalized differently (for single windows see above)*/
				if(normalize_b == 1){

					//printf("Before normalizing I_windowPairs: I12: %lf; Ia12: %lf; I1234_full: %lf; I1324_full: %lf;; I1423: %lf\n",  I_windowPairs.I12[10][10], I_windowPairs.Ia12[10][10], I_windowPairs.I1234_full[10][10], I_windowPairs.I1324_full[10][10], I_windowPairs.I1423[10][10]);
					normalizeQueryWindowPairsPtr(I_windowPairs, I_windowPairs_order2_meanStddev_DB, nrOfEntriesForPairs);
					//printf("After normalizing I_windowPairs: I12: %lf; Ia12: %lf; I1234_full: %lf; I1324_full: %lf;; I1423: %lf\n",  I_windowPairs.I12[10][10], I_windowPairs.Ia12[10][10], I_windowPairs.I1234_full[10][10], I_windowPairs.I1324_full[10][10], I_windowPairs.I1423[10][10]);

				}
									
				/*reset*/
				cntPairs = -1;
				cntPairsConsidered = 0;
				absMutVal = -1e10;
				maxAbsMutValQuery = -1e10;
				maxPosMutValQuery = -1e10;
				maxNegMutValQuery = -1e10;

				optionalValue = -1e10;
				optionalValue_pos = -1e10;
				optionalValue_neg = 1e10;


				for(k = 0; k < I_windowPairs.nrOfWindows; k++){

					for (l = k+1; l < I_windowPairs.nrOfWindows; l++){ //k+1: if l=k the two windows are id

						cntPairs += 1; /*indexes the window pairs*/

						mutVal = I_windowPairs.I12[k][l];
						absMutVal = fabs(mutVal);

						/*if(absMutVal > 40){
							printf("(saaledes stadig ganske hemmelig) , str name: %s k:%d l: %d\n",  I_windowPairs.structureName, k,l);
							//getchar();
						}*/
								
						/*reset*/
						if(useOnlyMaxMutVal_b == 1 && signed_b == 1){
							ptr_Irarity_windowPairs_pos[cntPairs].rarityScore = 0.0;
							ptr_Irarity_windowPairs_neg[cntPairs].rarityScore = 0.0;
						}
						else{
							ptr_Irarity_windowPairs[cntPairs].rarityScore = 0.0;
						}

						/*Skip the pair if it's not disjoint:*/
						if(onlyDisjointPairs_b ==1 && I_windowPairs.windowPair[k][l].segIndices_1[1] > I_windowPairs.windowPair[k][l].segIndices_2[0]){
							continue;
						}

						//cntPairs += 1; /*counts the disjoint pairs*/

						/*If we only consider the max-abs mutual writhe, skip the pair if the abs-mutual value is less than the max so far; else skip the pair if the 
						threshold on mutual value is not exceded*/
						if(useOnlyMaxMutVal_b == 1){

							/*for counting the number of pairs with mutVal aove the threshold:*/
							if (absMutVal >= thresholdMutual){
								cntMutAboveThreshold +=1;
							}

							/*skip this window (nr: l) if mutual value lower than previous ones, else collect/replace the window info (k,l); obs: the p-value is just a dummy, pre-set to 1:*/ 
							if(signed_b ==0){
								if(absMutVal> maxAbsMutValQuery){
									maxAbsMutValQuery = absMutVal;
									optionalValue = mutVal; //for writing the actual mutual value to file
									/*overwrite the Irarity info, copying over the info from the Iwindows struct*/
									collect_Irarity_windowPairs(ptr_Irarity_windowPairs, &I_windowPairs, 0, k,l, 0 );
								}
								else{
									continue;
								}
							}
							else{ //if signed_b != 0

								if(mutVal > 0 && absMutVal> maxPosMutValQuery){
									maxPosMutValQuery = absMutVal;
									optionalValue_pos = mutVal; //for writing the actual mutual value to file
									/*overwrite the Irarity info, copying over the info from the Iwindows struct*/
									collect_Irarity_windowPairs(ptr_Irarity_windowPairs_pos, &I_windowPairs, 0, k,l, 0 );
								}
								else if(mutVal < 0 && absMutVal > maxNegMutValQuery){
									maxNegMutValQuery = absMutVal;
									optionalValue_neg = mutVal; //for writing the actual mutual value to file
									/*overwrite the Irarity info, copying over the info from the Iwindows struct*/
									collect_Irarity_windowPairs(ptr_Irarity_windowPairs_neg, &I_windowPairs, 0, k,l, 0 );
								}
								else{
									continue;
								}

							}

						}
						else if(useOnlyMaxMutVal_b == 0){

							/*skip this window (nr: l) if mutual value too low:*/ 
							if (absMutVal < thresholdMutual){
								continue;
							}
							else{
								cntMutAboveThreshold +=1;
							}
						
							/*find the index of the invariant's absolute value in the background distribution*/
							hitAbsMutValPair = bisectionalSearchValueL(absMutVal, ptr_absMutualWritheDistr, nrOfDBrecordsLoaded);

							//printf("S... hitAbsMutValPair:%d nrOfDBrecordsLoaded:%d \n", hitAbsMutValPair, nrOfDBrecordsLoaded);

							pValueAbsMutVal = (1.0 - (double)(hitAbsMutValPair)/nrOfDBrecordsLoaded); //hitDistValue is the left end of the interval in which distMin sits; we add 1 to avoid numm disturbance
							//printf("... absMutValPValue: %lf\n",pValueAbsMutVal);
							//getchar();

							/*re-init the Irarity info, copying over the info from the Iwindows struct*/
							collect_Irarity_windowPairs(ptr_Irarity_windowPairs, &I_windowPairs, cntPairs, k,l, -logbase(pValueAbsMutVal + 1e-20, 10) ); //cntPairsConsidered rather than cntPairs?? No cntPairs corr's to pair k,l

							cntPairsConsidered +=1;

						} /*end of if(useOnlyMaxMutVal_b == 0)*/

					/////////////////////////

							
					} //end query l-window loop

				} //end query k-window loop

				//getchar();

				//printf("... cntPairs: %d cntPairsConsidered:%d\n", cntPairs, cntPairsConsidered)


				/*if we use the scoring by the max abs mutual (for the query) we find its p-value and score here
				and else we sum the scores from the individual pairs that we found in the preceding loop:*/ 
				if(useOnlyMaxMutVal_b == 1){

					if(signed_b == 0){

						/*find the index of the invariant's absolute value in the background distribution*/
						hitMaxAbsMutValQuery = bisectionalSearchValueL(maxAbsMutValQuery, ptr_maxAbsMutualWritheDistr, nrOfDBstructuresLoaded);

						//printf("... maxAbsMutValQuery: %lf hitMaxAbsMutValQuery:%d nrOfDBstructuresLoaded:%d \n", maxAbsMutValQuery, hitMaxAbsMutValQuery, nrOfDBstructuresLoaded);
						//getchar();

						pValueQuery = (1 - (double)(hitMaxAbsMutValQuery)/nrOfDBstructuresLoaded); //hitDistValue is the left end of the interval in which distMin sits; 
						//printf("Langmodigt venter bysvalen  ... p val: %lf\n", pValueQuery);
						score = -logbase(pValueQuery + 1e-20, 10);

						/*get the info on the window pair at which the maxAbsMutVal was found:*/
						strcpy(qWindowsStringAid, "\0");
						strcpy(qWindowsString, "[(");

						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[0].windowPair.windowNr_1);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ",");
						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[0].windowPair.segIndices_1[0]);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ",");
						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[0].windowPair.segIndices_1[1]);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, "),(");

						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[0].windowPair.windowNr_2);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ",");
						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[0].windowPair.segIndices_2[0]);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ",");
						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[0].windowPair.segIndices_2[1]);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ")]");

						/*collect the results (if use max-abs-mut_writhe the optionalValue is maxAbsMutValQuery, else it's just a default -9999):*/
						collect_queryRawScore(ptr_queryRawScore, 1, I_windows_dummy, I_windowPairs, qWindowsString, score, pValueQuery, cntScoredQueries, optionalValue);

					}
					else{

						/*We handle the positive and negative parts similarly but separately; positive first:*/

						/*find the index of the invariant's absolute value in the background distribution*/
						hitMaxPosMutValQuery = bisectionalSearchValueL(maxPosMutValQuery, ptr_maxPosMutualWritheDistr, nrOfDBstructuresLoaded);

						pValueQuery = (1 - (double)(hitMaxPosMutValQuery)/nrOfDBstructuresLoaded); //hitDistValue is the left end of the interval in which distMin sits; 
						//printf("Sommeren mild oedsles kyssene nemmere pluds'ligt ... p val: %lf\n", pValueQuery);
						score = -logbase(pValueQuery + 1e-20, 10);

						/*get the info on the window pair at which the maxAbsMutVal was found:*/
						strcpy(qWindowsStringAid, "\0");
						strcpy(qWindowsString, "[(");

						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs_pos[0].windowPair.windowNr_1);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ",");
						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs_pos[0].windowPair.segIndices_1[0]);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ",");
						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs_pos[0].windowPair.segIndices_1[1]);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, "),(");

						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs_pos[0].windowPair.windowNr_2);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ",");
						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs_pos[0].windowPair.segIndices_2[0]);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ",");
						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs_pos[0].windowPair.segIndices_2[1]);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ")]");

						/*collect the results (if use max-abs-mut_writhe the optionalValue is maxAbsMutValQuery, else it's just a default -9999):*/
						collect_queryRawScore(ptr_queryRawScore_pos, 1, I_windows_dummy, I_windowPairs, qWindowsString, score, pValueQuery, cntScoredQueries, optionalValue_pos);

						/*Positive part now done, negative next: */

						/*find the index of the invariant's absolute value in the background distribution*/
						hitMaxNegMutValQuery = bisectionalSearchValueL(maxNegMutValQuery, ptr_maxNegMutualWritheDistr, nrOfDBstructuresLoaded);

						pValueQuery = (1 - (double)(hitMaxNegMutValQuery)/nrOfDBstructuresLoaded); //hitDistValue is the left end of the interval in which distMin sits; 
						//printf("(saaledes stadig ganske hemmelig) ... p val: %lf\n", pValueQuery);
						score = -logbase(pValueQuery + 1e-20, 10);

						/*get the info on the window pair at which the maxAbsMutVal was found:*/
						strcpy(qWindowsStringAid, "\0");
						strcpy(qWindowsString, "[(");

						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs_neg[0].windowPair.windowNr_1);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ",");
						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs_neg[0].windowPair.segIndices_1[0]);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ",");
						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs_neg[0].windowPair.segIndices_1[1]);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, "),(");

						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs_neg[0].windowPair.windowNr_2);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ",");
						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs_neg[0].windowPair.segIndices_2[0]);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ",");
						snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs_neg[0].windowPair.segIndices_2[1]);
						strcat(qWindowsString, qWindowsStringAid);
						strcat(qWindowsString, ")]");

						/*collect the results (if use max-abs-mut_writhe the optionalValue is maxAbsMutValQuery, else it's just a default -9999):*/
						collect_queryRawScore(ptr_queryRawScore_neg, 1, I_windows_dummy, I_windowPairs, qWindowsString, score, pValueQuery, cntScoredQueries, optionalValue_neg);

					}



				}
				else if(useOnlyMaxMutVal_b == 0){

					score = 0.0;
					strcpy(qWindowsStringAid, "\0");
					strcpy(qWindowsString, "\0");

					if(writeWindowInfo_b == 0){

						strcat(qWindowsString, "Getting scores per window pair was not called. Set writeWindowInfo_b (option: Y) to 1 to change.");

					}


					for(n = 0; n < cntPairs; n++){

						score += ptr_Irarity_windowPairs[n].rarityScore;

						/*get the window info where there was a contribution to the score:*/
						if(ptr_Irarity_windowPairs[n].rarityScore > 0){

							snprintf(rarityScoreString, 10, "%lf",ptr_Irarity_windowPairs[n].rarityScore); 

							if(writeWindowInfo_b != 0){
								strcat(qWindowsString, "[(");
								snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.windowNr_1);
								strcat(qWindowsString, qWindowsStringAid);
								strcat(qWindowsString, ",");
								snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.segIndices_1[0]);
								strcat(qWindowsString, qWindowsStringAid);
								strcat(qWindowsString, ",");
								snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.segIndices_1[1]);
								strcat(qWindowsString, qWindowsStringAid);
								strcat(qWindowsString, "),(");

								snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.windowNr_2);
								strcat(qWindowsString, qWindowsStringAid);
								strcat(qWindowsString, ",");
								snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.segIndices_2[0]);
								strcat(qWindowsString, qWindowsStringAid);
								strcat(qWindowsString, ",");
								snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.segIndices_2[1]);
								strcat(qWindowsString, qWindowsStringAid);
								strcat(qWindowsString, "),");
								strcat(qWindowsString, rarityScoreString);
								strcat(qWindowsString, "]");
							}


						}

					}

					//For averaging it could seem reasonable to use cntPairsConsidered rather than cntPairs: the back ground data base should 
					//be build with the same constraints: only disjoint pairs or not? and the same threshold on mutuals. Else, using the 
					//actual number of window pairs (of the the same disjoint-ness as in the data base), we get scores that reflect the hypothesis
					//that the number of rare configurations should be increasing (quadratically) with the length of the structure, and which are, to 
					//some extent, comparable across diff values of the threshold:
					if(getAvgScore_b == 1){

						if(cntPairs == -1){
							score = 0.0;
						}
						else{
							score /= (cntPairs+1); //+1 to avoid dividing by 0 
						}
					}

					/*if desired (and a file containing rarityScores is provided), compute the p-value of the rarity factor:*/
					if(absMutValScorePValues_b == 1){

						/*If the score is 0 we assign the most "conservative" p-value, and else we compute it
						by looking up the score value in the score-distribution*/
						if(score == 0.0){
							pValueQuery = 1;
						}
						else{

							hitAbsMutValScore = bisectionalSearchValueL(score, ptr_absMutValScoresDistr, sizeOfAbsMutValScoresDistr);

							//printf("Langmodig venter bysvalen .. ... score: %lf hitAbsMutValScore:%d sizeOfAbsMutValScoresDistr:%d \n", score, hitAbsMutValScore, sizeOfAbsMutValScoresDistr);

							pValueQuery = (1.0 - (double)(hitAbsMutValScore)/sizeOfAbsMutValScoresDistr) + 1e-20; //hitDistValue is the left end of the interval in which distMin sits; we add 1 to avoid numm disturbance
							//printf("... p val: %lf\n", pValueQuery);
							//getchar();
						}

					}

					/*collect the results (if use max-abs-mut_writhe the optionalValue is maxAbsMutValQuery, else it's just a default -9999):*/
					collect_queryRawScore(ptr_queryRawScore, 1, I_windows_dummy, I_windowPairs, qWindowsString, score, pValueQuery, cntScoredQueries, optionalValue);

				}
					

				cntScoredQueries += 1;

				if (write_windowPairs_b != 0){

					returnVal = writeInvariantsOnwindowPairs(ptr_fileNameOutwindowPairs, I_windowPairs, L, order, onlyDisjointPairs_b);
				}


				//Reset! The reason is that the pairs-loop may turn out to be empty -- no pairs fulfull the criteria in there; 
				//if omitting this reset the data from the most recent chain which had a pair fulfilling the criteria will be used in the ptr_queryRawScore/_pos/_neg!!
				if(useOnlyMaxMutVal_b == 1 && signed_b ==1 ){
					init_ptr_Irarity_windowPairs(&ptr_Irarity_windowPairs_pos, alreadyAllocatedTo_maxNrOfWindows);
					init_ptr_Irarity_windowPairs(&ptr_Irarity_windowPairs_neg, alreadyAllocatedTo_maxNrOfWindows);
				} else{
					init_ptr_Irarity_windowPairs(&ptr_Irarity_windowPairs, alreadyAllocatedTo_maxNrOfWindows);
				}


				endChain = clock();
				timeChain = ((double)(endChain - startChain))/CLOCKS_PER_SEC;
				//printf("CPU time spend for measures on this file incl load of chain and more:%lf \n", timeChain);
				if (timeChain > max_timeChain){
					max_timeChain = timeChain;
				}


				subStructureCnt += 1; /*will count and index the sub-structures for which we compute the invariants*/
			
			} /*end chainNr loop*/

			
			fclose(ptr_fileIn);

		} /*end fileNr loop*/


		//printf("... du er skoer, Topper!\n");

		/*sort the obtained scoring results:*/
		if(useOnlyMaxMutVal_b == 1 && signed_b == 1){
			heapSortRawScores(ptr_queryRawScore_pos, cntScoredQueries);
			heapSortRawScores(ptr_queryRawScore_neg, cntScoredQueries);
		}
		else{
			heapSortRawScores(ptr_queryRawScore, cntScoredQueries);
		}

		if(write_rarityScoresPairs_b == 1){

			if(useOnlyMaxMutVal_b == 1 && signed_b == 1){

				for(i=0;i<cntScoredQueries;i++){

					fprintf(ptr_fileOutScoresPairs_pos, "%s;%s;%s;%s;%s;%d;%d;%s;%lf;%lf;%lf;\n",
								DBName,
								ptr_queryRawScore_pos[i].fileName,
								ptr_queryRawScore_pos[i].structureName,
								ptr_queryRawScore_pos[i].chainId,
								ptr_queryRawScore_pos[i].classId,
								ptr_queryRawScore_pos[i].chainLen,
								ptr_queryRawScore_pos[i].nrOfWindows, //query nr of windows
								ptr_queryRawScore_pos[i].windowInfo,
								ptr_queryRawScore_pos[i].score,
								ptr_queryRawScore_pos[i].pValue,
								ptr_queryRawScore_pos[i].optionalValue
								);

					
					fprintf(ptr_fileOutScoresPairs_neg, "%s;%s;%s;%s;%s;%d;%d;%s;%lf;%lf;%lf;\n",
								DBName,
								ptr_queryRawScore_neg[i].fileName,
								ptr_queryRawScore_neg[i].structureName,
								ptr_queryRawScore_neg[i].chainId,
								ptr_queryRawScore_neg[i].classId,
								ptr_queryRawScore_neg[i].chainLen,
								ptr_queryRawScore_neg[i].nrOfWindows, //query nr of windows
								ptr_queryRawScore_neg[i].windowInfo,
								ptr_queryRawScore_neg[i].score,
								ptr_queryRawScore_neg[i].pValue,
								ptr_queryRawScore_neg[i].optionalValue
								);

				}

				/*close files for holding scores*/
				fclose(ptr_fileOutScoresPairs_pos);
				fclose(ptr_fileOutScoresPairs_neg);

			}
			else{

				for(i=0;i<cntScoredQueries;i++){

					fprintf(ptr_fileOutScoresPairs, "%s;%s;%s;%s;%s;%d;%d;%s;%lf;%lf;%lf;\n",
								DBName,
								ptr_queryRawScore[i].fileName,
								ptr_queryRawScore[i].structureName,
								ptr_queryRawScore[i].chainId,
								ptr_queryRawScore[i].classId,
								ptr_queryRawScore[i].chainLen,
								ptr_queryRawScore[i].nrOfWindows, //query nr of windows
								ptr_queryRawScore[i].windowInfo,
								ptr_queryRawScore[i].score,
								ptr_queryRawScore[i].pValue,
								ptr_queryRawScore[i].optionalValue
								);

				}

				/*close file for holding scores*/
				fclose(ptr_fileOutScoresPairs);

			}
					
		}

		//printf("Sommeren mild oesles kyssene nemmere pluds'ligt\n");

		endComplete = clock();
		//printf("start clocks:%d, end clocks: %d, clocks_per sec:%ld\n", endComplete, startComplete, CLOCKS_PER_SEC);
		timeComplete = ((double) (endComplete - startComplete))/CLOCKS_PER_SEC;
		printf("CPU time spend for all files incl loading data:%lf \n", timeComplete);
		printf("CPU time spend for computations on all files:%lf \n", compTimeComplete);
		printf("Max CPU time spend for measures across all files:%lf \n", max_timeChain);

		printf("Done number of files:%d \n", fileNr);
		printf("Done number of chains/sub-structures:%d \n", subStructureCnt);
		printf("%d too short and %d too long chains were skipped\n", chainsSkippedTooShort,chainsSkippedTooLong);
		printf("%d chains were skipped due to a too long segment (longer than sqrt of %d)\n",chainSkippedSegmTooLong, stdRealSegLength);


		printf("cntMutAboveThreshold: %d\n", cntMutAboveThreshold);

		/*If user wants to score another queries set without reloading the data base: 
		Prompt user a new queries set, else set keepGoing_b = 0 (will then leave the outer while-loop):*/
		keepGoing_b = 0;
		printf("Do you want to run another scan? (yes: 1, no: 0) ");
		fgets(keepGoingStr, 3, stdin);
		printf( "You entered: %s\n", keepGoingStr);
		keepGoing_b = atoi(keepGoingStr);
		//printf("keepGoing_b %d\n", keepGoing_b);
		if(keepGoing_b == 1){

			killroyWasHere_b = 1;

			/*Resets:*/
			strcpy(queriesDirPath, "\0");
			strcpy(queriesName, "\0");
			strcpy(outputPath, "\0");

			printf("Enter a full path to another queries set:");
			fgets(queriesDirPath, 1000, stdin);
			queriesDirPath[strcspn(queriesDirPath, "\n")] = 0; //removes trailing \n if any!
			printf("Queries path: %s", queriesDirPath );
			printf("Enter a name for this queries set:");
			fgets(queriesName, 100, stdin);
			queriesName[strcspn(queriesName, "\n")] = 0; //removes trailing \n if any!

			printf("Enter a full path to where the output should be placed:");
			fgets(outputPath, 1000, stdin);
			//strcat(outputPath,"\0" );
			outputPath[strcspn(outputPath, "\n")] = 0; //removes trainling \n if any!
			printf("Output path: %s", outputPath );

			printf("Is the new queries set from CATH, SCOP or none of these: (0,1 or 2)?");
			fgets(dataTypeIndicatorStr, 3, stdin);
			dataTypeIndicator = atoi(dataTypeIndicatorStr);
			use_cath_b = 0;
			use_scop_b = 0;
			dataTypeIndicator = 2;
			if(dataTypeIndicator == 0){
				use_cath_b = 1;
				printf("You have chosen to use CATH-data; if you want to only consider domains in a positive list, you can enter a full path to it here:\n");
				fgets(cathListPath, 1000, stdin );
				outputPath[strcspn(cathListPath, "\n")] = 0; //removes trainling \n if any!
			}
			else if(dataTypeIndicator == 1){
				use_scop_b = 1;
			}

			/*Reset file names:*/
		    strcpy(fileNameOutScoresPairs1, "");
			strcpy(fileNameOutScoresPairs1, "/RarityScan0_ScoresPairs_windowslgth_");
			strcpy(fileNameOutwindowPairs1,"/RarityScan0_MutualWrithe_Pairs_windowlgth_");	
			strcpy(fileNameOut2, "_order_");

			strcpy(fileNameOut3,"\0");

			strcpy(fileNameOutwindowPairs, "\0");
			strcpy(fileNameOutScoresPairs, "\0");
			strcpy(fileNameOutScoresPairs_pos, "\0");
			strcpy(fileNameOutScoresPairs_neg, "\0");

			/*Other resets*/
			if(useOnlyMaxMutVal_b == 1 && signed_b == 1){
				reinit_queryRawScore(ptr_queryRawScore_pos, cntScoredQueries);
				reinit_queryRawScore(ptr_queryRawScore_neg, cntScoredQueries);
			}
			else{
				reinit_queryRawScore(ptr_queryRawScore, cntScoredQueries);
			}
		

			cntMutAboveThreshold = 0;
			cntScoredQueries = 0; 
			fileNr = 0;
			chainNr = 0;
			maxChainLen = 0;
			subStructureCnt = 0;

			chainsSkippedTooShort = 0;
			chainsSkippedTooLong = 0;
			chainSkippedSegmTooLong = 0;


		}

	} /*end of outer while!*/


	/*free the memory allocated to this query set*/
	free_I_measures(I_measures, order, full_b, alreadyAllocatedTo_chainLen);
	free(ptr_segment);	
	free(ptr_chain);
	free(dirContent.ptr_dirList); //this is allocated in the ListDirectoryContents2 fct 
	free(dirContent.ptr_fileNameList); //this is allocated in the ListDirectoryContents2 fct 
	//Resets:
	dirContent.numberOfFiles = 0;
	subDirCnt = 0;
	//More freeing:
	for(i=0;i< alreadyAllocatedTo_subStructureCnt; i++){
		free(ptr_queryRawScore[i].fileName);
		free(ptr_queryRawScore[i].structureName);
		free(ptr_queryRawScore[i].classId);
		free(ptr_queryRawScore[i].chainId);
		free(ptr_queryRawScore[i].windowInfo);
	}
	free(ptr_queryRawScore);


	//getchar(); 


	returnVal = 1;

	return returnVal;
}


/*
Function for deciding if a structure (query) stands out as "rare" compared to a background (based on a set of PDB
files). The scoring is done by means of a "cross entropy" type. We often use the term "data base" for the back ground. The code 
allows running a set of queries against the data base; in fact it is set up so that one may run a series of query sets against the 
same back ground (withou reloading the back ground). 


The usage is:
1) Create the data base: this is done by running the main function here (SAonGISA_main_*) in the makeDB flavour. This wraps a call 
to the GI_windows function of the GISA_v*_unix-code. Must be run in "window pairs mode" (and) so that a file of mutual writhe values 
or, more generally, for the desired Gauss Integrals invariants, for all pairs of disjoint windows for each structure in the PDB-set 
is created. Set the order according to the set of invariants wanted (ie 1 if only the writhe and, possibly, the average crossing 
number are desired, else set the order to 2). It is though recommended only to use order 1, since the mutuals of the higher order 
invariants are not (all) true mutuals, which makes the results harder to interpret. Set windowing parameters as desired (e.g. windowLength 
20, stepSize 4); use windowType 0.
2) Now run the rarity scan with rawRarity1: the search method is aimed at detecting the frequency of data base window pairs 
with Gauss numbers in "epsilon-distance" to those of a given query pair (of windows); for a given query structure these 
frequencies are collected and a score is obtained by accumulation. The "epsilon" is set indirectly through binning of the invariant 
values. 

To obtain the p-values corresponding to the obtained scores output, it is necessary to first run rawRarity1 with query set = 
data base set and with rarityScorePValues_b = 0 (and else leave the settings). When done run rawRarity1 for the desired 
query set now with rarityScorePValues_b = 1.

It is possible to set a threshold filtering away all pairs with abs mutual value lower than this cutoff level. For speed it
is (very) desirable to set this threshold quite high; when the scan is only based on the writhe (and no other invariants) the 
threshold can be set to  e.g. 10 (and leave out normalization of the invariant values by setting the normalize_b parameter at 0). 
When running with more invariants (e.g. also the average crossing number) the values must be normalized (will be done to get mean 0 
and st dev 1) and the threshold must then be set differently (maybe to 0.5); this may need a little experimenting (e.g. have a look 
at the actual normalized invariant values -- these values can be written out).

An outline of the method:
a) The data base is converted to a "list of words": each data base element is a tuple/array of Gauss numbers (in number as many as the 
desired number of invariants for matching, and in any case limited the number of invariants of order 1 or order 2); each of these 
tuples is converted by binning the Gauss numbers into a tuple of integers (ie a "word"); after sorting these words lexicographically
the data base has the guise of a dictionary (though probably with many words repeated).
b) A given query is similarly translated, by the same binning, into a set of words (one for each window pair); the window pairs are 
now looped through and each is looked up in the data base dictionary (by setting a threshold as mentioned above so that only pairs 
having a mutual writhe (or mutual invariant) in absolute value above this threshold are considered). Either exact matching or matching 
allowing a set number of mismatches can be done. In both cases, the look-up gives a count of the number of matches, cntMatch, for 
each pair. The score is now the sum

Score = - Sum over query pairs log(cntMatch(pairs)/#data base) /#pairs in query

which equals

Score = - Sum over words (w)  ( #pairs in query of word w) /#pairs in query * log( #db-pairs of word w/#data base) 

If we write p_q (w ) = ( #pairs of word w) /#pairs in query and p_db (w) = #db-pairs of word w/#data base we have

Score = - Sum over words (w) p_q(w) * log(p_db(w)),

a cross entropy that is, and also showing the resemblance to the Kullback-Leibler relative entropy

KL = Sum over words (w) p_q(w) * log(p_q(w)) / p_db(w))

ie only the "q-idiosyncratic" term "Sum over words(w) p_q(w) * log(p_q(w))" is disregarded. 


The scoring method scans for structures having a "distribution of words" significantly different from that found in the background
distribution. One can also think of this as a way of determining whether a query has an unsual set of fragment pairs as
compared to the background. The threshold allows to focus on occurence of "rarer words".

Final scores are written to a file, in decreasing order ( ie with lowest p-value at the top). 
*/

int rawRarity1(  int halted_b, char rarityScorePairsDistrFilePath[1000], 
			    char rarityScorePairsDistrName[100], char DBName[100], char queriesName[100], 
			    char queriesDirPath[200], char DBResultsPairsFileName[1000], 
			    char fileNameQueryChains[2000], 
			    int onlyDisjointPairs_b,
			    double thresholdMutual,
				int rarityScorePValues_b, int use_scop_b, 
				int use_cath_b, char *cathListPath, int loadFromSubDirs_b, 
				int normalize_b, char *outputPath, int maxChainLength, 
				int nrOfEntriesForPairs, int binningType, 
				int nrOfBins, double compensationForRounding, int allowedNrOfMismatchesPairs,
				int nrOfWindowsForMatch, int chunkSize, int order, 
				int incl_abs_b, int full_b,  
				int write_matchWindowPairs_b, int windowCoveringType, int windowLength, 
				int stepSize, int write_windowPairs_b,  int write_rarityScoresPairs_b,
				int write_chains_b, int writeWindowInfo_b, int print_b, 
				int print_basic_b){

	int cntRuns = 0; //for testing code

	int matchWindowPairs_b = 1; /*SHOULD NOT BE CHANGED*/ 
	
	int use_lex = 0; /*keep this!*/

	int returnVal = 0;

	int split_b = 0;

	double rarityFactorPairs; //to measure how rare a GI-vector for pairs is (in the first filtering)
	int cntMutAboveThreshold = 0;
	//int windowLength = 16;
	//int write_windows_b = 1; /*for writing out I-values on windows for queries*/
	//int write_windowPairs_b = 0; /*at most one of write_windows_b and write_windowPairs_b should be set to 1; corr to get_windos_b resp. get_windowPairs_b*/
	/*for matching query to DB (set of structures in the set directory:*/
	int k,l, queryNrOfWindowPairs;
	int u;
	double maxMutVal, mutVal;
	//int nrOfWindowsForMatch = 3;
	
	double ***ptr_query_pairs, ***ptr_match_pairs; 

	//int strLengthCutOff = 10; /*chains shorter than this int will not be included in computation of invariants*/

	struct dirContent dirContent;
	int subDirCnt = 0;
	int numberOfFiles = 0;
	int subStructureCnt = 0; /*some structures are multi-mers; we compute the invariants for each sub-structure (each mono-mer)*/
	int maxNrOfWindows = 0; /*max nr of windows across loaded structures*/
	int maxNrOfWindowPairs = 0; /*max nr of window pairs across loaded structures*/
	const int numberOfPerturbations = 1; /*only there since needed when writing out the final values (preparation for perturbations code)*/
	const int pertNo = 0; /*ditto*/
	char **ptr_dirList;
	char **ptr_fileNameList;

	/*for "dynamically" generating the file names for output:*/
	char orderStr[10] = "-1"; //place holder value
	char inclAbsStr[10] = "-1"; //value is placeholder
	char nrInvsForMatchStr[10] = "-1"; //place holder value
	char nrWinsForMatchStr[10] = "-1"; //place holder value
	char allowedNrOfMismatchesPairsStr[10] = "-1"; //place holder value; 
	char binTypeStr[10] = "-1"; //place holder value
	char nrBinsStr[10] = "-1"; //place holder value
	char windowLengthStr[10] = "-1"; //place holder value
	char windowStepSizeStr[10] = "-1"; //value is placeholder
	char windowCovTypeStr[10] = "-1"; //value is placeholder
	char topNStr[10] = "-1"; //value is placeholder
	char thresholdMutualStr[10] = "-1"; //value is placeholder

	char fileNameOutScoresPairs1[1000] = "/RarityScan1_ScoresPairs_windowslgth_";
	char fileNameOutwindowPairs1[1000] = "/RarityScan1_Invariants_Pairs_windowlgth_";	
	char fileNameOut2[1000] = "_order_"; //; //; //"_minmax_";//"_incl_abs_"; //"_"; //

	char fileNameOut3[1000] = "\0";  //"_testQuery_top8000_5invs_norm.txt"; 
	char extensionName[50] = ".txt";

	char fileNameOutwindowPairs[2000] = "\0";
	char fileNameOutScoresPairs[2000] = "\0";


	char *ptr_fileNameOutwindowPairs; /*pointer to file for holding measures' values on pairs of sub-chains of pre-set windowlength*/
	FILE *ptr_fileOutwindowPairs;
	char *ptr_fileNameOutScoresPairs;
	FILE *ptr_fileOutScoresPairs;


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
	int structureNameLength = 7; /*string-length of structure name in DB, e.g SCOP, CATh*/
	char *currentStructureName;
	char *structureNameInit = "NN";
	char *currentClassId;
	char *classIdInit = "NN";
	char chainIdInit[2] = ">";
	int chainsSkippedTooShort = 0;
	int chainsSkippedTooLong = 0;
	int chainSkippedSegmTooLong = 0;
	int segTooLong_b = 0;
	//int resNr = 0; /*residue number*/
	//int resNrPrior = -999;
	int chainLen = 0; /*length of C-alpha chain*/
	int maxChainLen =0; /*max length of chain among all loaded*/
	int L = 0; /*length of segment chain = chainLen -1*/
	int simplexCnt = 0;/*size of simplex*/
	int maxSimplexCnt = 0; /*size of simplex corr to max chain length*/

	struct cAlpha * ptr_chain = NULL; /*for holding C-alpha chain */
	struct segment *ptr_segment = NULL; /* for holding chain of segments, ith C-alpha to jth C-alpha*/
	struct segment segCoords;	 

	int m = 0, n = 0, mMin, nMin;
	int i = 0;
	int j = 0;
	int cnt = 0;
	int fileNr = 0;
	int chainNr = 0;

	/* for holding meausure values (i.e. the core of the final output); 
	content to be received from aggr function:*/
	struct I_ptr I_measures;
	double ***I_measures_raw; /*triple pointer to replace the struct I_ptr use*/ 
	int alreadyAllocatedTo_chainLen  = 0;
	int alreadyAllocatedTo_maxNrOfWindows = 0;
	int alreadyAllocatedTo_maxNrOfWindowPairs = 0;
	int alreadyAllocatedTo_subStructureCnt = 0;

	/*pointer to hold the results (that this is a ptr of a ptr is due to that the write-out can also
	handle the perturbation case; we use pertNo = 0 here):*/
	struct I_values **ptr_I_values;

	/*struct for holding invariants on windows resp. window pairs for the query and for each DB entry:*/
	struct I_windowPairs_ptr I_windowPairs; //for query
	struct I_windows_ptr I_windows_dummy; //just placeholder for call to collect_queryRawScore
	double **bins;
	struct binned_I_windowPairs *ptr_binned_I_windowPairs_query; //for query, binned values
	struct I_windowPairs_ptr db_I_windowPairs;
	int nrOfDBstructuresLoaded = 0; 
	int nrOfDBrecordsLoaded = 0;
	int nrOfDBstructuresLoaded_pairsLoad = 0; 
	int nrOfDBrecordsLoaded_pairsLoad = 0;
	int maxNrOfWindowsDB = 0;
	int maxNrOfWindowPairsDB = 0;	
	/*pointer to hold DB results*/
	struct db_I_windowPairs_order2_ptr *ptr_dbResultWindowPairs;
	struct I_windowPairs_order2_meanStddev I_windowPairs_order2_meanStddev_DB;
	struct binned_I_counts *ptr_binned_I_counts;
	int wordCount = 0;
	int areId;
	int currentAllocSize_binned_I_counts_binned, newAllocSize_binned_I_counts_binned;
	int hitIndexQueryWord;


	/*for searching*/
	int hitLR[2];
	int *intArray2; //convenience array
	int *intArray3; //convenience array
	//struct binned_I_windows *ptr_binned_I_windows_DB;
	struct binned_I_windowPairs *ptr_binned_I_windowPairs_DB;
	int matchSetNrOfStructures = 1000;
	struct binned_I_windowPairs *ptr_binned_I_windowPairs_matchSet;
	struct binned_I_windowPairs **ptr_ptr_binned_I_windowPairs_matchSet;
	int db_fileNr;//for window pairs searching
	int db_wNr1;//for window pairs searching
	int db_wNr2;//for window pairs searching
	int cntMatch;
	int *ptr_cntMatch; /*to hold the cntMatch for the series of windows covering a query*/
	struct Irarity_windowPairs *ptr_Irarity_windowPairs;

	int nrOfMismatchesThisPair = 0;
	int cntMatch_12 = 0;
	int cntPairs = 0;
	int cntPairsConsidered = 0;


	int *ptr_matchIndicator;
	struct matchRange **ptr_newRange;
	struct matchRange range; 
	

	double score = 0.0; 
	int cntScoredQueries = 0; 
	double *ptr_rarityScorePairsDistr; //sorted list of dist-values (to be had from ...)
	FILE *ptr_fileRarityScorePairsDistr;
	int sizeOfRarityScorePairsDistr = -1;
	struct queryRawScore *ptr_queryRawScore;
	int hitRarityScorePairs;
	double pValuePairs = -1;
	char * qWindowsString;
	char * qWindowsStringAid;
	char * rarityScoreString;

	/*for closed loops finding, but only place holders here:*/
	int closed_loops_b = 0; //must be kept like this!
	struct twoSegmentIndex *ptr_closedLoopInd;
	int initValNrOfClosedLoops = 0;
	int nrOfClosedLoops = 0;
	int cntClosedLoop = 0;

	FILE *fp; //dummy file ptr used when checking file names

	/*params to control the execution*/
	int keepGoing_b = 1; 
	char keepGoingStr[10] = "-1"; //place holder value
	int dbLoaded_b = 0;
	char dataTypeIndicatorStr[3] = "-1"; //place holder value
	int dataTypeIndicator = -1;
	int killroyWasHere_b = 0;


	/*Declarations stop here*/


	/*The code for running one set of queries against the data base is wrapped in an outer while-loop, which allows the user to
	enter a new queris set, without having to reload the data base*/
	while(keepGoing_b == 1){

		/****************************************************************************************************************/
		/*Generation of output file names:*/
		/****************************************************************************************************************/

		/*we first build a string (fileNameOut3) shared by most file names, which contains various params values;
		fileNameOut3 will be concatenated into the file name strings below:*/

		sprintf(nrWinsForMatchStr, "%d", nrOfWindowsForMatch);
		strcat(fileNameOut3, nrWinsForMatchStr);
		strcat(fileNameOut3, "wins_");
		if(normalize_b == 1){
			strcat(fileNameOut3, "norm_");
		}
		else{
			strcat(fileNameOut3, "notNorm_");
		}
		sprintf(nrInvsForMatchStr, "%d", nrOfEntriesForPairs);
		strcat(fileNameOut3, nrInvsForMatchStr);
		strcat(fileNameOut3, "invs_");
		if(matchWindowPairs_b == 1){
			sprintf(allowedNrOfMismatchesPairsStr, "%d", allowedNrOfMismatchesPairs);
			strcat(fileNameOut3, allowedNrOfMismatchesPairsStr);
			strcat(fileNameOut3, "mmsPairs_");
		}
		sprintf(nrBinsStr, "%d", nrOfBins);
		strcat(fileNameOut3, nrBinsStr);
		strcat(fileNameOut3, "bins");
		sprintf(binTypeStr, "%d", binningType);
		strcat(fileNameOut3,binTypeStr);
		strcat(fileNameOut3, "_winCovType");
		sprintf(windowCovTypeStr, "%d", windowCoveringType);
		strcat(fileNameOut3, windowCovTypeStr);
		strcat(fileNameOut3, "_threshMut");
		sprintf(thresholdMutualStr, "%0.2f", thresholdMutual);
		strcat(fileNameOut3, thresholdMutualStr);


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
			strcat(fileNameOutwindowPairs, queriesName);
			strcat(fileNameOutwindowPairs, "_");
			strcat(fileNameOutwindowPairs, DBName);
			strcat(fileNameOutwindowPairs, "_");
			strcat(fileNameOutwindowPairs, fileNameOut3);
			strcat(fileNameOutwindowPairs, extensionName);
			printf("Results for window pairs will be written to: %s\n", fileNameOutwindowPairs);

			 
			ptr_fileNameOutwindowPairs = fileNameOutwindowPairs;
			
			/*clear the file by opening it in w-mode and then closing it again*/
			fp= fopen(ptr_fileNameOutwindowPairs, "w");
			if(fp == NULL){
				printf("Sorry, file name %s is invalid. Could the path be wrong? A blank space?", ptr_fileNameOutwindowPairs);
				getchar();
				exit(EXIT_FAILURE);
			}
			else{
				fclose(fp);
			}
		

		}



		if (write_rarityScoresPairs_b !=0){

			/*generate file name for writing out match values*/
			sprintf(orderStr, "%d", order);
			strcat(fileNameOutScoresPairs, outputPath);
			strcat(fileNameOutScoresPairs, fileNameOutScoresPairs1);
			sprintf(windowLengthStr, "%d", windowLength);
			strcat(fileNameOutScoresPairs, windowLengthStr);
			strcat(fileNameOutScoresPairs, "_");
			sprintf(windowStepSizeStr, "%d", stepSize);
			strcat(fileNameOutScoresPairs, windowStepSizeStr);
			strcat(fileNameOutScoresPairs, fileNameOut2);
			strcat(fileNameOutScoresPairs, orderStr);
			strcat(fileNameOutScoresPairs, "_");
			sprintf(inclAbsStr, "%d", incl_abs_b);
			strcat(fileNameOutScoresPairs, inclAbsStr);
			strcat(fileNameOutScoresPairs, "_");
			strcat(fileNameOutScoresPairs, queriesName);
			strcat(fileNameOutScoresPairs, "_");
			strcat(fileNameOutScoresPairs, DBName);
			strcat(fileNameOutScoresPairs, "_");
			strcat(fileNameOutScoresPairs, fileNameOut3);
	
			/*when getting the rarity score p-values, inlcude info on which distribution was used:*/
			if(rarityScorePValues_b == 1){
				strcat(fileNameOutScoresPairs, rarityScorePairsDistrName);
			}

			
			strcat(fileNameOutScoresPairs, extensionName);
			printf("Rarity scores will be written to: %s . Press enter to continue.\n", fileNameOutScoresPairs);
			
			ptr_fileNameOutScoresPairs = fileNameOutScoresPairs;
			
			/*clear the file by opening it in w-mode and then closing it again*/
			ptr_fileOutScoresPairs = fopen(ptr_fileNameOutScoresPairs, "w");

			fprintf(ptr_fileOutScoresPairs, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;\n",
				"DB name",
				"Query file",
				"structureName",
				"chainId",
				"classId",
				"chainLen",
				"nrOfWindows",
				"windowInfo",
				"rarityScore",
				"pValue"
				);  

			/*fprintf(ptr_fileOutScoresPairs, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;\n",
				"DB name",
				"Query file",
				"structureName",
				"chainId",
				"classId",
				"chainLen",
				"nrOfWindows",
				"qWinNr_1",
				"i1",
				"i2",
				"qWinNr_2",
				"j1",
				"j2",
				"rarityScore",
				"pValue"
				);  */
		

		}

		if(halted_b == 1){
			printf("If the output paths are fine, press enter to continue.\n");
			getchar();
		}

		/****************************************************************************************************************/
		/*Get content of the desired directory: number of files and a list of the file names:*/
		/****************************************************************************************************************/

		if(loadFromSubDirs_b ==0){
			dirContent = ListDirectoryContents(queriesDirPath);
		}
		else{
			returnVal = ListDirectoryContents2(queriesDirPath, &dirContent, &subDirCnt);
		}
		numberOfFiles = dirContent.numberOfFiles; 
		ptr_dirList = dirContent.ptr_dirList;
		ptr_fileNameList = dirContent.ptr_fileNameList;


		
		printf("Number of query files in directory:%d \n",numberOfFiles);
		printf("1st file: %s\n",ptr_fileNameList[0] );
		printf("Press enter to continue\n");
		if(halted_b == 1){
				getchar();
		}


		startComplete = clock();

		/****************************************************************************************************************/
		/*loop through the list of filenames in the directory to find chain info, eg the max chain length (for setting
		the size in the (global) mem alloc to I_measures's ptr's). First though we make a bulk allocation
		to the chain-info keeping pointer, using the pre-set maxNrOfChains:*/
		/****************************************************************************************************************/

		if(killroyWasHere_b == 0){
			returnVal = alloc_init_chainInStr(&chainInStr, maxNrOfChains);
		}


		/*if CATH data, read in CATH domain list. Note: we abuse the struct chainInStructure to accomodate
		domains rather than chains -- in CATH each chain may be split in several domains*/ 
		if(use_cath_b == 1){

			if(convToCLF20_b == 1){
				structureNameLength = lengthCLFformat20;
			}
			else{
				structureNameLength = 6;
			}

			/*first allocate a block of memory (if not large enough it will be extended when exectuing readCATHDomainList)*/
			ptr_cathDomain = (struct chainInfo *) malloc (maxNrOfDomains*sizeof(struct chainInfo));
			/*init*/
			for(i = 0; i< maxNrOfDomains; i++){
				ptr_cathDomain[i].chainId = (char *) malloc(sizeof(char)*2);
				ptr_cathDomain[i].structureName = (char *) malloc(sizeof(char)*(structureNameLength+2)); //to hold domain name
				ptr_cathDomain[i].classId = (char *) malloc(sizeof(char)*10);
				ptr_cathDomain[i].classIdConv = (int *) malloc(classifierSize*sizeof(int));

				strcpy(ptr_cathDomain[i].chainId, chainIdInit); /*init to unlikely value*/
				strcpy(ptr_cathDomain[i].structureName, structureNameInit);
				strcpy(ptr_cathDomain[i].classId, classIdInit);
				/*printf("Chain Id initial value:%c\n", *ptr_cathDomain[i].chainId);
				printf("Class Id initial value:%c\n", *ptr_cathDomain[i].classId);*/
				ptr_cathDomain[i].chainLength = -1;

				for(j=0; j< classifierSize; j++){ptr_cathDomain[i].classIdConv[j] = -1;}

			}

			/*open CATH file list:*/
			ptr_cath_domain_list_file = fopen(cathListPath, "r");

			if (!ptr_cath_domain_list_file){
				printf("Sorry, the file containing the positive-list of CATH domains was not found (%s was provided). Enter a valid file name and run again.\n", cathListPath);
				return 0;
			}

			/*read in contents of CATH domain list file*/
			nrOfCATHDomainsInList = readCATHDomainList(ptr_cath_domain_list_file, &ptr_cathDomain, convToCLF20_b );

			printf("nr of cath doms: %d\n", nrOfCATHDomainsInList);

			/*lex-sort the ptr_cath_domain_list_file on the domain name entry. Below we
			need to look up each file for which we have PDB-coords in this list. In other
			words, this list serves as a "positive list" of domains, for which to compute
			the invariants.*/
			heapSortCATHDomain(ptr_cathDomain, nrOfCATHDomainsInList, structureNameLength,0);

			/*for(i=0;i <nrOfCATHDomainsInList;i++){
				if(i%100 == 0){
					for(j=i;j <i+100;j++){
						printf("dom %d: %s\n", j, ptr_cathDomain[j].structureName);
					}
					getchar();
				}
			}*/
		}

		printf("Fetching chain info for data base structures ...\n");
		for(fileNr = 0; fileNr <numberOfFiles; fileNr ++){

			if (fileNr%100 == 0){printf(".. now fetched the chain info for %d structures\n", fileNr);}

			if (print_b == 1){
				printf("File name:%s path:%s fileNr:%d\n",ptr_fileNameList[fileNr], ptr_dirList[fileNr], fileNr);
			}

			/*if using CATH data: check if the file, ptr_fileNameList[fileNr], is in the
			domain list specified (and loaded and lex-sorted above)*/
			if(use_cath_b == 1){
				hitIndexDomain = bisectionalSearchCATHDomain(ptr_fileNameList[fileNr], ptr_cathDomain, nrOfCATHDomainsInList, lengthCLFformat20);
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
				continue;
			}

			/*get number of chains, their id's and lengths for this file:*/
			//chainInStr = readPDBChainStructure(ptr_fileIn);
			if(use_scop_b ==1){
				readPDBChainStructureSCOP(ptr_fileIn, &chainInStr );
			}
			else if(use_cath_b ==1){
				readPDBDomainStructureCATH(ptr_fileIn, &chainInStr );//will read chain length and chain id into chainInStr
			}
			else{
				readPDBChainStructure2(ptr_fileIn, &chainInStr );
			}


			fclose(ptr_fileIn);  


			//printf("chainInStr.numberOfChains: %d\n", chainInStr.numberOfChains);
			for (i=0;i <= chainInStr.numberOfChains -1;i++){
				//printf("chainInStr.ptr_chainInStr[i].chainLength: %d\n", chainInStr.ptr_chainInStr[i].chainLength);
				if(chainInStr.ptr_chainInStr[i].chainLength <= maxChainLength){ //we disregard very long chains for RAM availability reasons
					maxChainLen = max(chainInStr.ptr_chainInStr[i].chainLength,maxChainLen);
					if(chainInStr.ptr_chainInStr[i].chainLength >= strLengthCutOff){
										
						/*if using CATH data: we checked inside the fileNr-loop if the file, ptr_fileNameList[fileNr], is in the
						domain list specified (and loaded and lex-sorted above). We need another check to avoid loading in 
						several chains for the structureName under consideration:*/
						if(use_cath_b == 1){ //if not using CATH, chainInStr.ptr_chainInStr[i].classId was read in in the pdb-file reads above
							
							/*In case the pdb-file contains several chains: the domain is attached to a single chain (is a subset of 
							a single chain) and we want only to consider that chain. we therefore skip "wrong chains". Obs: hitIndexDomain
							was found above:*/
							if(chainInStr.numberOfChains > 1 && ptr_cathDomain[hitIndexDomain].structureName[4] != chainInStr.ptr_chainInStr[i].chainId[0]){

								if(print_basic_b > 0){
									printf("CATH domain: %s is regarded not to sit in this chain: %s\n", ptr_cathDomain[hitIndexDomain].structureName, chainInStr.ptr_chainInStr[i].chainId);
								}
								continue;
							}

						}

						subStructureCnt += 1; /*count and index the sub-structures for which the invariants will be computed*/
					
					}
				}
			}

		}

		printf("... done.\n");

		L = maxChainLen -1;


		printf("max chain len:%d\n", maxChainLen);
		printf("sub structure count:%d\n", subStructureCnt);
		if(halted_b == 1){
			getchar();
		}

		/****************************************************************************************************************/
		/*MEMORY ALLOCATION*/
		/****************************************************************************************************************/

		/*We allocate memory for ptr holding the I-values and a ptr for holding the 3d-coord's, both for the longest chain*/
		returnVal = alloc_init_I_measures(&I_measures, order, full_b, chainNr, maxChainLen, closed_loops_b, alreadyAllocatedTo_chainLen, &ptr_closedLoopInd);

		if(alreadyAllocatedTo_chainLen == 0){ //no allocation was done yet

			ptr_chain = (struct cAlpha *) malloc (maxChainLen*sizeof(struct cAlpha));
			/*init*/
			for ( n=0; n <= maxChainLen -1 ; n++ ){
				
				ptr_chain[n].coords.x = 0.0;
				ptr_chain[n].coords.y = 0.0;
				ptr_chain[n].coords.z = 0.0;

				ptr_chain[n].residueNr = -123;

			}
			ptr_segment = (struct  segment *) malloc (maxChainLen*sizeof(struct segment)); 
			/*init*/
			segCoords.s1.x = 0.0;
			segCoords.s1.y = 0.0;
			segCoords.s1.z = 0.0;
			segCoords.s2.x = 0.0;
			segCoords.s2.y = 0.0;
			segCoords.s2.z = 0.0;
			for ( n=0; n <= maxChainLen -2 ; n++ ){
					segCoords.s1 = ptr_chain[n].coords;
					segCoords.s2 = ptr_chain[n+1].coords;

					ptr_segment[n] = segCoords;

			}

		}
		else if(maxChainLen > alreadyAllocatedTo_chainLen){ //reallocation necc

			ptr_chain = realloc (ptr_chain, maxChainLen*sizeof(struct cAlpha));

			/*(re)init*/
			for ( n=0; n <= maxChainLen -1 ; n++ ){
				
				ptr_chain[n].coords.x = 0.0;
				ptr_chain[n].coords.y = 0.0;
				ptr_chain[n].coords.z = 0.0;

				ptr_chain[n].residueNr = -123;

			}

			ptr_segment = realloc (ptr_segment, maxChainLen*sizeof(struct segment));
			/*init*/
			segCoords.s1.x = 0.0;
			segCoords.s1.y = 0.0;
			segCoords.s1.z = 0.0;
			segCoords.s2.x = 0.0;
			segCoords.s2.y = 0.0;
			segCoords.s2.z = 0.0;
			for ( n=0; n <= maxChainLen -2 ; n++ ){

					segCoords.s1 = ptr_chain[n].coords;
					segCoords.s2 = ptr_chain[n+1].coords;

					ptr_segment[n] = segCoords;

			}


		}
		alreadyAllocatedTo_chainLen = max(alreadyAllocatedTo_chainLen, maxChainLen);


		/****************************************************************************************************************/
		/* Allocate mem to some ptr's for query: to hold I-values on window pairs, the binned version and for the scores	
		/****************************************************************************************************************/

		///*For query*/
		maxNrOfWindows =  (int) ceil((double) 1.1*(maxChainLen-windowLength)/stepSize); /*this should just suffice for the longest chain -- we're using overlapping windows*/
		
		/*We'll need I-values on the window pairs:*/
		maxNrOfWindowPairs = alloc_init_I_windowPairs_ptr(&I_windowPairs, order, full_b, maxChainLen, windowLength, stepSize, maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows);
		
		//printf("maxNrOfWindowPairs %d\n",maxNrOfWindowPairs );	
		
		///*allocate and initialize a ptr to hold the scores results for each query*/
		printf("subStructureCnt %d  alreadyAllocatedTo_subStructureCnt %d\n", subStructureCnt,alreadyAllocatedTo_subStructureCnt );
		printf("maxWindowInfoSize*maxNrOfWindowPairs %d maxWindowInfoSize*alreadyAllocatedTo_maxNrOfWindowPairs %d\n", maxWindowInfoSize*maxNrOfWindowPairs, maxWindowInfoSize*alreadyAllocatedTo_maxNrOfWindowPairs);
		returnVal = alloc_init_queryRawScore(&ptr_queryRawScore, subStructureCnt, alreadyAllocatedTo_subStructureCnt, maxWindowInfoSize*maxNrOfWindowPairs, maxWindowInfoSize*alreadyAllocatedTo_maxNrOfWindowPairs); /* maxWindowInfoSize = 50 ... should suffice*/
		
		/*allocate memory to struct binned_I_windowPairs and initialize*/
		returnVal = alloc_init_ptr_binned_I_windowPairs(&ptr_binned_I_windowPairs_query, maxNrOfWindowPairs, alreadyAllocatedTo_maxNrOfWindowPairs, nrOfEntriesForPairs);	

		/*for holding rarity scores to be written out:*/
		returnVal = alloc_init_ptr_Irarity_windowPairs(&ptr_Irarity_windowPairs, order, full_b, maxChainLen, windowLength, stepSize, maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows);

		/*convenience ptr for holding MUTUAL I values on window pairs, for query.
		On the diagonal these will though amount to the single window values:*/
		if(alreadyAllocatedTo_maxNrOfWindows == 0){ //indicates no alloc ws done yet
		
			ptr_query_pairs = (double ***) malloc(maxNrOfWindows*sizeof(double **));

			/*alloc and initialize*/
			for(k = 0; k < maxNrOfWindows; k++){
						
				ptr_query_pairs[k] = (double **) malloc(maxNrOfWindows*sizeof(double *));

				/*initialize*/
				for(l = 0; l < maxNrOfWindows; l++){
						
					ptr_query_pairs[k][l] = (double *) malloc(maxNrOfInvariants*sizeof(double));

					ptr_query_pairs[k][l][0] = 0.0;
					ptr_query_pairs[k][l][1] = 0.0;
						
					if(nrOfEntriesForPairs > 2){

						for(i= 2; i < maxNrOfInvariants;i++){
							ptr_query_pairs[k][l][i] = 0.0;

						}
					}
				}
			}

		}
		else if(maxNrOfWindows > alreadyAllocatedTo_maxNrOfWindows){ //need to alloc some more mem

			ptr_query_pairs = realloc(ptr_query_pairs, maxNrOfWindows*sizeof(double **));

			/*alloc and initialize*/
			for(k = 0; k < maxNrOfWindows; k++){
						
				ptr_query_pairs[k] = realloc(ptr_query_pairs[k], maxNrOfWindows*sizeof(double *));

				/*initialize*/
				for(l = 0; l < maxNrOfWindows; l++){
						
					ptr_query_pairs[k][l] = realloc(ptr_query_pairs[k][l], maxNrOfInvariants*sizeof(double));

					ptr_query_pairs[k][l][0] = 0.0;
					ptr_query_pairs[k][l][1] = 0.0;
						
					if(nrOfEntriesForPairs > 2){

						for(i= 2; i < maxNrOfInvariants;i++){
							ptr_query_pairs[k][l][i] = 0.0;

						}
					}
				}
			}

		}

		//Update
		alreadyAllocatedTo_maxNrOfWindows = max(maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows);
		alreadyAllocatedTo_maxNrOfWindowPairs = max(maxNrOfWindowPairs, alreadyAllocatedTo_maxNrOfWindowPairs);
		alreadyAllocatedTo_subStructureCnt = max(subStructureCnt, alreadyAllocatedTo_subStructureCnt);
		

		/****************************************************************************************************************/
		/*Load of DB for scoring/matching */
		/****************************************************************************************************************/

		if(dbLoaded_b == 0){
			
			/*We load in DB results for pairs; then bin the whole set (transforming the I-value tuples into "words"); sort lexicographically
			to obtain a "dictionary of the I-value tuples" */


			/*Load DB for the pair I-values:*/
			intArray3 = readDBresultsWindowPairsToPtr(&ptr_dbResultWindowPairs, DBResultsPairsFileName, chunkSize, &maxNrOfWindowsDB);

			nrOfDBstructuresLoaded_pairsLoad = intArray3[0];
			nrOfDBrecordsLoaded_pairsLoad =  intArray3[1]; //equals nr of window pairs loaded

			if(halted_b == 1){
				printf("Pairs load; no str's: %d pairs: %d\n", nrOfDBstructuresLoaded_pairsLoad, nrOfDBrecordsLoaded_pairsLoad);
				getchar();
			}

			/*for(m = 0; m < 5; m++){
				printf("db file at %d is %s, I12 at 0 0: %lf seg1 1: %d seg2 0: %d\n", m, ptr_dbResultWindowPairs[m].fileName, ptr_dbResultWindowPairs[m].I12[0][0], ptr_dbResultWindowPairs[m].segIndices_1[0][1], ptr_dbResultWindowPairs[m].segIndices_2[0][0]);
			}*/

			/*Convenience ptr for holding window MUTUAL values in ptr; we allocate some surplus:*/
			ptr_match_pairs = (double ***) malloc(maxNrOfWindowsDB*sizeof(double **));

			for (m=0; m<maxNrOfWindowsDB;m++){

				ptr_match_pairs[m] = (double **) malloc(maxNrOfWindowsDB*sizeof(double *));
				
				for (n=0; n<maxNrOfWindowsDB;n++){
					
					ptr_match_pairs[m][n] = (double *) malloc(maxNrOfInvariants*sizeof(double));

					ptr_match_pairs[m][n][0] = 0.0;
					ptr_match_pairs[m][n][1] = 0.0;	

					if(nrOfEntriesForPairs > 2){

						for(i = 2; i < maxNrOfInvariants; i++){
							ptr_match_pairs[m][n][i] = 0.0;

						}
					}

					
				}
			}

			nrOfDBstructuresLoaded = nrOfDBstructuresLoaded_pairsLoad; // intArray2[0];
			nrOfDBrecordsLoaded = nrOfDBrecordsLoaded_pairsLoad; // intArray2[1]; //equals nr of windows loaded

				
			/*normalize the DB results if desired*/
			if(normalize_b == 1){

				printf("normalizing DB result (pairs data) ..\n");

				if(matchWindowPairs_b == 1){

					I_windowPairs_order2_meanStddev_DB = normalizeDBresultsWindowPairsPtr(ptr_dbResultWindowPairs, nrOfDBstructuresLoaded_pairsLoad, nrOfDBrecordsLoaded_pairsLoad);

					if(print_basic_b == 1){
						printf("Pairs normalization: I12 move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_I12, I_windowPairs_order2_meanStddev_DB.stddev_I12);
						printf("Pairs normalization: Ia12 move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_Ia12, I_windowPairs_order2_meanStddev_DB.stddev_Ia12);
						printf("Pairs normalization: I1234_full move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_I1234_full, I_windowPairs_order2_meanStddev_DB.stddev_I1234_full);
						printf("Pairs normalization: I1324_full move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_I1324_full, I_windowPairs_order2_meanStddev_DB.stddev_I1324_full);
						printf("Pairs normalization: I1423 move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_I1423, I_windowPairs_order2_meanStddev_DB.stddev_I1423);
					}

				}

				printf(".. normalization done \n");
				if(halted_b == 1){
					getchar();
				}

			}


			/*Now transform the DB records (I-vectors) into binned versions and sort it all lexicographically*/
			/*First transform to array of binned I-results:*/
			printf("Binning DB results (pairs data) ... first generate bins\n");

			/*OBS: if binningType = 0 we use a custom "hand-held" binning defined in the generateBins fct; change
			of bin definitions is needed when nrOfBins change!!*/
			returnVal = generateBinsPairs(binningType, ptr_dbResultWindowPairs, nrOfDBstructuresLoaded, nrOfDBrecordsLoaded, nrOfEntriesForPairs, &bins, nrOfBins, compensationForRounding);

			printf("... done ... next the actual binning\n");


			/*OBS: the ptr ptr_binned_I_windows_DB is allocated to within the binDBresultsWindows fct! (unlike for ptr_binned_I_windows_query/binQueryresultsWindows)*/
			returnVal = binDBresultsWindowPairs(&ptr_binned_I_windowPairs_DB, ptr_dbResultWindowPairs , nrOfDBstructuresLoaded, nrOfDBrecordsLoaded, nrOfEntriesForPairs, bins, nrOfBins, onlyDisjointPairs_b);
			printf(".. done \n");


			/*for(i=0; i<nrOfDBrecordsLoaded;i++){
				if(i%50 == 0){getchar();}
				printf("rec nr %d: %d \n", i, ptr_binned_I_windowPairs_DB[i].binnedIvector[0]); //, ptr_binned_I_windowPairs_DB[i].binnedIvector[1], ptr_binned_I_windows_DB[i].binnedIvector[2], ptr_binned_I_windows_DB[i].binnedIvector[3], ptr_binned_I_windows_DB[i].binnedIvector[4]);
					//	if(i == 25){ break;}
			}*/

			/*Next sort lexicographically*/
			printf("sorting DB results (pairs data) ..\n");

			heapSortBinnedDBPairs(ptr_binned_I_windowPairs_DB, nrOfDBrecordsLoaded, nrOfEntriesForPairs);

			/*for(i=0; i<nrOfDBrecordsLoaded;i++){
				//if(i < 100){
					printf("rec nr %d: %d \n", i, ptr_binned_I_windowPairs_DB[i].binnedIvector[0]); //, ptr_binned_I_windowPairs_DB[i].binnedIvector[1], ptr_binned_I_windows_DB[i].binnedIvector[2], ptr_binned_I_windows_DB[i].binnedIvector[3], ptr_binned_I_windows_DB[i].binnedIvector[4]);
					//getchar();
					if(i%100 == 0){ getchar();}
				//}
			}
			getchar();*/

			printf(".. done \n");


			/****************************************************************************************************************/
			/* Create an array holding the number of occurences of/size of match set for each word in the data base; 
			/****************************************************************************************************************/
			
			printf("Create array of DB records holding the possible words and match set sizes ...\n");

			/*if we are not allowing mismatches, we populate the array right away (the created array is then automatically lex-sorted); this 
			will allow a very simple look-up of the size of the match set for a given query word (window pair). When we allow 
			mismatches we will record the size of the match set when handling a query (below); for each query (window pair) we'll 
			first look up in this array and only if it is not there we'll find the match set (and then add that word and its size 
			of match set to the array):*/
			
			currentAllocSize_binned_I_counts_binned = ((int) ceil(nrOfDBrecordsLoaded/10.0));
			//printf("I start by allocating: %d\n", currentAllocSize_binned_I_counts_binned);
			ptr_binned_I_counts = (struct binned_I_counts *) malloc(currentAllocSize_binned_I_counts_binned*sizeof(struct binned_I_counts)); /*just a guess; we reallocate more mem if necessary*/
			for(j=0;j<currentAllocSize_binned_I_counts_binned;j++){

				ptr_binned_I_counts[j].binnedIvector = (int *) malloc(nrOfEntriesForPairs*sizeof(int));
				/*init to unlikely values*/
				for(i=0;i<nrOfEntriesForPairs; i++){ ptr_binned_I_counts[j].binnedIvector[i] = 999999; /*high positive to retain lex-order*/}
				ptr_binned_I_counts[j].count = -1;
			}

			if(allowedNrOfMismatchesPairs == 0){

				wordCount = 0;
				//init:
				for(j=0; j< nrOfEntriesForPairs; j++){ptr_binned_I_counts[wordCount].binnedIvector[j] = ptr_binned_I_windowPairs_DB[0].binnedIvector[j];}
				//ptr_binned_I_counts[0]->binnedIvector = ptr_binned_I_windowPairs_DB[0]->binnedIvector;
				ptr_binned_I_counts[wordCount].count = 0;
				for(n=0;n<nrOfDBrecordsLoaded; n++){

					/*Is current word in ptr_binned_I_counts[wordCount] id to the previous one in the loop?:*/
					areId = 1;
					for(j=0; j< nrOfEntriesForPairs; j++){
						if(ptr_binned_I_counts[wordCount].binnedIvector[j] != ptr_binned_I_windowPairs_DB[n].binnedIvector[j]){
							areId = 0;
						}
					}
					

					if(areId == 1){

						ptr_binned_I_counts[wordCount].count +=1;

					}
					else{

						wordCount +=1;
						for(j=0; j< nrOfEntriesForPairs; j++){ptr_binned_I_counts[wordCount].binnedIvector[j] = ptr_binned_I_windowPairs_DB[n].binnedIvector[j];}
						//ptr_binned_I_counts[wordCount].binnedIvector = ptr_binned_I_windowPairs_DB[n].binnedIvector;
						ptr_binned_I_counts[wordCount].count = 1;

					}

					/*realloc if nec:*/
					if(wordCount > currentAllocSize_binned_I_counts_binned){			

						newAllocSize_binned_I_counts_binned = currentAllocSize_binned_I_counts_binned + ((int) ceil(nrOfDBrecordsLoaded/10.0)); 
						printf("I'm reallocating from %d to %d\n", currentAllocSize_binned_I_counts_binned, newAllocSize_binned_I_counts_binned);
						ptr_binned_I_counts = (struct binned_I_counts *) realloc(ptr_binned_I_counts, newAllocSize_binned_I_counts_binned*sizeof(struct binned_I_counts));
						for(j=currentAllocSize_binned_I_counts_binned;j<newAllocSize_binned_I_counts_binned;j++){

							ptr_binned_I_counts[j].binnedIvector = (int *) malloc(nrOfEntriesForPairs*sizeof(int));
							/*init to unlikely values*/
							for(i=0;i<nrOfEntriesForPairs; i++){ ptr_binned_I_counts[j].binnedIvector[i] = 999999; /*high positive to retain lex-order*/}
							ptr_binned_I_counts[j].count = 0;
						}

						currentAllocSize_binned_I_counts_binned = newAllocSize_binned_I_counts_binned;

					
					}



				}

				if(halted_b == 1){
					for(m = 0; m < wordCount; m++){
						printf("Word cnt: %d count, binnedIvector: %d,  reg'ed: %d", m, ptr_binned_I_counts[m].binnedIvector[0], ptr_binned_I_counts[m].count);
						getchar();
					}
				}

			}

			printf("... done\n");



			/****************************************************************************************************************/
			/* Allocate and initialize ptr to hold the match set returned from each "word search"; 
			/****************************************************************************************************************/

			printf("Create a match set array ... \n");
			
			/*we use this when looping through the match set so as to avoid recomputing multiple times for each matching structure.
			For assistance we declare a ptr holding a 0 or 1 for each DB record indicating whether the record
			is a match (1) or not (0).*/

			/*array of "GI words", one word for every possible match*/
			returnVal = alloc_init_ptr_binned_I_windowPairs(&ptr_binned_I_windowPairs_matchSet, nrOfDBrecordsLoaded, 0, nrOfEntriesForPairs);

			/*match indicator assistant --- gets reset for every new query (below)*/
			ptr_matchIndicator = (int *)malloc(nrOfDBrecordsLoaded*sizeof(int)); //(int *) malloc((nrOfDBrecordsLoaded+1)*sizeof(int));
			for(i = 0; i < nrOfDBrecordsLoaded; i++){ 
				ptr_matchIndicator[i] = 0;
			}

			printf("... done ... and DB load now completed.\n");

			dbLoaded_b = 1; 		

		}

		if(halted_b == 1){
			getchar();
		}

		/****************************************************************************************************************/
		/* Read in rarityScore-data and generate the distribution of these (just a sort):*/
		/****************************************************************************************************************/

		if(rarityScorePValues_b == 1){
				
			ptr_fileRarityScorePairsDistr = fopen(rarityScorePairsDistrFilePath, "r");

			if (!ptr_fileRarityScorePairsDistr){
				printf("Sorry, the file containing the back ground scores for the pairs, rarityScorePairsDistrFilePath, was not found (%s was provided). Enter a valid file name and run again.\n", rarityScorePairsDistrFilePath);
				return 0;
			}


			/*this will provide a sorted list (smallest to highest), ptr_distDistr, of dist-values.
			Will be used to score each window hit (as -log(p-value), where the p-value is obtained by a
			look up in this array: */
			sizeOfRarityScorePairsDistr = readSortRarityScoreData(ptr_fileRarityScorePairsDistr, &ptr_rarityScorePairsDistr);

			printf("sizeOfRarityScorePairsDistr %d\n", sizeOfRarityScorePairsDistr);
			/*for(i=0;i<200;i++){
				printf("rarity score at %d after sorting: %lf\n", i, ptr_rarityScorePairsDistr[i]);
				if(i%50 == 0){getchar();}
			}*/


		}
		
		/*Allocation to a pointer holding a match range for each entry (Gauss integral) to be matched (ie in number nrOfEntriesForPairs)*/  
		printf("nrOfEntries (pairs): %d\n", nrOfEntriesForPairs);
		printf("size of matchRange: %d\n", sizeof(struct matchRange));
		printf("size of matchRange ptr: %d\n", sizeof(struct matchRange *));

		ptr_newRange = (struct matchRange **) malloc(3*(nrOfEntriesForPairs+1)*sizeof(struct matchRange));
		
		if(ptr_newRange != NULL){
			for(i = 0; i < nrOfEntriesForPairs; i++){

				ptr_newRange[i] = (struct matchRange *) malloc(3*sizeof(struct matchRange));
				
				if(ptr_newRange[i] != NULL){
					for(j = 0; j < 3; j++){
						ptr_newRange[i][j].start = -1;
						ptr_newRange[i][j].end = -2;
						ptr_newRange[i][j].nrOfMismatches = 999;
					}
				}
				else{
					printf("newRange[%d] not allocated for\n", i);
				}
			}
		}
		else{
			printf("post ptr_newRange alloc, nrOfEntries (pairs) now: %d\n", nrOfEntriesForPairs);
			if(halted_b == 1){
				getchar();
			}
		}


		/*Finally a couple of allocations to strings for read-out of scores:*/
		if(write_rarityScoresPairs_b == 1){
			qWindowsStringAid = (char *) malloc(5*sizeof(char));
			qWindowsString = (char *) malloc(20*maxNrOfWindowPairs*sizeof(char));
			rarityScoreString = (char *) malloc(10*sizeof(char));

			/*init*/
			strcpy(qWindowsStringAid, "\0");
			strcpy(qWindowsString, "\0");
			strcpy(rarityScoreString, "\0");
		}



		/****************************************************************************************************************/
		/* MAIN LOOP */
		/****************************************************************************************************************/


		subStructureCnt = 0; /*reset*/

		//startComplete = clock();

		/*loop through the list of filenames in the queries' directory and compute the desired measures all along:*/
		for (fileNr = 0; fileNr <numberOfFiles; fileNr ++){

			if (fileNr%10 == 0){printf("Now at file nr:%d name: %s\n", fileNr, ptr_fileNameList[fileNr]);}

			//printf("Query: %s\n", ptr_dirList[fileNr]);
			//getchar();

			if (print_basic_b == 1){
				printf("Now handling file:%s\n fileNr:%d\n",ptr_dirList[fileNr], fileNr);
			}

			
			if(use_cath_b == 1){
				hitIndexDomain = bisectionalSearchCATHDomain(ptr_fileNameList[fileNr], ptr_cathDomain, nrOfCATHDomainsInList, lengthCLFformat20);
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
				printf("Sorry, the file %s could not be found. Execution stops here.\n", ptr_dirList[fileNr]);
				//getchar();
				return 0;
			}


			/*reinit chain info*/
			returnVal = reinit_chainInStr(chainInStr, maxNrOfChains);
					

			/*get chain info for the file (number of chains and their lengths for this file etc):*/
			if(use_scop_b ==1){
				readPDBChainStructureSCOP(ptr_fileIn, &chainInStr );
			}
			else if(use_cath_b ==1){
				readPDBDomainStructureCATH(ptr_fileIn, &chainInStr );//will read chain length and chain id into chainInStr; other values (structureName etc) will be read in below
				//printf("File: %s chainInStr.numberOfChains: %d \n",I_measures.fileName, chainInStr.numberOfChains);
				//getchar();
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
				}
			}

			/****************************************************************************************************************/
			/* INNER LOOP */
			/****************************************************************************************************************/
			/*loop through the chains in the current structure/file ... and compute ... :*/

			chainNr = 0; /*reset*/
			for (chainNr = 0; chainNr <= chainInStr.numberOfChains -1 ; chainNr ++){

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
							printf("CATH domain: %s is regarded not to sit in this chain: %s\n", ptr_cathDomain[hitIndexDomain].structureName, chainInStr.ptr_chainInStr[chainNr].chainId);
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

				I_measures.structureName = chainInStr.ptr_chainInStr[chainNr].structureName;
				I_measures.classId = chainInStr.ptr_chainInStr[chainNr].classId;
				I_measures.chainId = chainInStr.ptr_chainInStr[chainNr].chainId;
				I_measures.chainNr = &chainNr;

				chainLen = chainInStr.ptr_chainInStr[chainNr].chainLength;

				if (print_basic_b ==1){
					printf("sub str:%s\n",I_measures.structureName);
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

				
				if(use_scop_b ==1){
					if(strstr(chainInStr.ptr_chainInStr[chainNr].classId, "i") == chainInStr.ptr_chainInStr[chainNr].classId || 
						strstr(chainInStr.ptr_chainInStr[chainNr].classId, "j") == chainInStr.ptr_chainInStr[chainNr].classId ||
						strstr(chainInStr.ptr_chainInStr[chainNr].classId, "k") == chainInStr.ptr_chainInStr[chainNr].classId ||
						strstr(chainInStr.ptr_chainInStr[chainNr].classId, "h") == chainInStr.ptr_chainInStr[chainNr].classId){
						printf("Class %s of chain of %s was skipped since it is h, i, j or k\n", chainInStr.ptr_chainInStr[chainNr].chainId, ptr_dirList[fileNr], chainInStr.ptr_chainInStr[chainNr].classId);
						subStructureCnt += 1; //!
						continue;
					}
				}

				I_measures.chainLen = &chainLen;

				rewind(ptr_fileIn);
				returnVal = main_readPDB2(ptr_fileIn, ptr_chain, chainNr, chainLen);

				/*length of segments chain:*/
				L = chainLen - 1;
				
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
					printf("Chain: %s in structure: %s removed due to too long %d 'th segment\n", chainInStr.ptr_chainInStr[chainNr].chainId, I_measures.fileName,  segTooLong_b-1);
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
					//returnVal = aggrAndW_raw(ptr_segment, chainLen, order, full_b, &I_measures_raw);
					}
				if (incl_abs_b == 0 && split_b == 0 && closed_loops_b == 0){
					returnVal = aggrAndW_ExAbs(ptr_segment, chainLen, order, full_b, I_measures);
					}
				if (incl_abs_b == 1 && split_b == 1 && closed_loops_b == 0){
					returnVal = wAll(ptr_segment, chainLen, I_measures.wVal);
					returnVal = aggr(chainLen, order, full_b, I_measures.wVal, I_measures);
				}


				endComp = clock();

				compTime = ((double) ((endComp - startComp))/CLOCKS_PER_SEC);
				//compTime = ((double) (endComp - startComp));
				if(print_b == 1){
					printf("Comp time: %f\n",compTime);
					//printf("CPU time spend for measures on this file (only computation/aggregation):%lf \n", compTime);
				}
				if (compTime > max_compTime){
					max_compTime = compTime;
				}
				compTimeComplete += compTime;


				/*Code for searching through the DB starts here*/

				/**************************************************************************************/
				/*Scanning/matching (based on window pairs)*/
				/**************************************************************************************/

				if (write_chains_b != 0){
					returnVal = writeChain(fileNameQueryChains, I_measures.fileName, I_measures.structureName, I_measures.chainId, ptr_chain, *I_measures.chainLen);
				}

				returnVal = getInvariantsOnWindowPairs(&I_windowPairs, order, windowCoveringType, windowLength, stepSize, L, I_measures);

				/*normalize the pairs' I-values if desired. The single window values (I_windows) and the window
				pair values (I_windowPairs) are normalized differently*/
				if(normalize_b == 1){

					//printf("Before normalizing I_windowPairs: I12: %lf; Ia12: %lf; I1234_full: %lf; I1324_full: %lf;; I1423: %lf\n",  I_windowPairs.I12[10][10], I_windowPairs.Ia12[10][10], I_windowPairs.I1234_full[10][10], I_windowPairs.I1324_full[10][10], I_windowPairs.I1423[10][10]);
					normalizeQueryWindowPairsPtr(I_windowPairs, I_windowPairs_order2_meanStddev_DB, nrOfEntriesForPairs);
					//printf("After normalizing I_windowPairs: I12: %lf; Ia12: %lf; I1234_full: %lf; I1324_full: %lf;; I1423: %lf\n",  I_windowPairs.I12[10][10], I_windowPairs.Ia12[10][10], I_windowPairs.I1234_full[10][10], I_windowPairs.I1324_full[10][10], I_windowPairs.I1423[10][10]);

				}

				/*for later convenience load the query's I values for pairs (mutuals) into ptr's allocated for the purpose*/
				for (k = 0; k < I_windowPairs.nrOfWindows; k++){

					for (l = k; l < I_windowPairs.nrOfWindows; l++){

						ptr_query_pairs[k][l][0] = I_windowPairs.I12[k][l];
						//printf("q I12 at %d %d :%lf\n", k, l, I_windowPairs.I12[k][l]);
						if (nrOfEntriesForPairs >= 2){
							ptr_query_pairs[k][l][1] = I_windowPairs.Ia12[k][l];
						}
						if (nrOfEntriesForPairs > 2){
							ptr_query_pairs[k][l][2] = I_windowPairs.I1234_full[k][l];
							ptr_query_pairs[k][l][3] = I_windowPairs.I1324_full[k][l];
							ptr_query_pairs[k][l][4] = I_windowPairs.I1423[k][l];


							ptr_query_pairs[k][l][5] = I_windowPairs.Ia12a34_full[k][l];
							ptr_query_pairs[k][l][6] = I_windowPairs.Ia13a24_full[k][l];
							ptr_query_pairs[k][l][7] = I_windowPairs.Ia14a23[k][l];

							ptr_query_pairs[k][l][8] = I_windowPairs.Ia1234_full[k][l];
							ptr_query_pairs[k][l][9] = I_windowPairs.Ia1324_full[k][l];
							ptr_query_pairs[k][l][10] = I_windowPairs.Ia1423[k][l];

							ptr_query_pairs[k][l][11] = I_windowPairs.I12a34_full[k][l];
							ptr_query_pairs[k][l][12] = I_windowPairs.I13a24_full[k][l];
							ptr_query_pairs[k][l][13] = I_windowPairs.I14a23[k][l];

						}
					
					}
				
				}
						
				/*bin the query's I-window Pair values*/
				returnVal = binQueryResultsWindowPairs(ptr_binned_I_windowPairs_query, I_windowPairs, ptr_query_pairs, nrOfEntriesForPairs, bins, nrOfBins, onlyDisjointPairs_b);
		
				/*reset*/
				cntPairs = -1;
				cntPairsConsidered = 0;

				for(k = 0; k < I_windowPairs.nrOfWindows; k++){

					for (l = k + 1; l < I_windowPairs.nrOfWindows; l++){

						/*Skip the pair if it's not disjoint:*/
						if(I_windowPairs.windowPair[k][l].segIndices_1[1] > I_windowPairs.windowPair[k][l].segIndices_2[0]){
							continue;
						}

						//cntPairs += 1; /*counts disjoints pairs*/s

						cntPairs += 1; /*will index the window pairs as when binning the pairs' values (in fct binQueryResultsWindowPairs)*/
						/*reset*/
						rarityFactorPairs = 1.0;
						ptr_Irarity_windowPairs[cntPairs].rarityScore = 0.0;
						
						/*Skip the pair if the threshold on mutual value is not exceded:*/
						/*only if the mutual value ptr_query_pairs[k+u][l+u] is large enough do we care to find a match.
						But only when l != k, ie not for the pairs which actually are only a single window -- since we want to
						match also windows that have low I-values :*/
						maxMutVal = 0;
						if(l>k){
							for (u = 0; u <nrOfEntriesForPairs; u++){ 
									mutVal = ptr_query_pairs[k][l][u];
								if(fabs(mutVal) > maxMutVal){
									maxMutVal = fabs(mutVal);
								}
							}
						}
						/*skip this window (nr: l) if mutual value too low:*/ 
						if (l>k && maxMutVal < thresholdMutual){continue;}

						cntMutAboveThreshold +=1;
						
						//printf("Q: pairs nr %d: %d %d %d %d %d\n", cntPairs, ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[0], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[1], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[2], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[3], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[4]);
						//printf("Q: pairs nr %d: %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", cntPairs, ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[0], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[1], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[2], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[3], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[4], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[5], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[6], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[7], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[8], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[9], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[10], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[11], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[12], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[13]);
						//getchar();


						/*we look-up the word ptr_binned_I_windowPairs_query[cntPairs] in the array of word 
						counts (ptr_binned_I_counts). If we do not allow mismatches this will give us the
						size of the match set (the number of id words in the data base); if we do allow
						mismatches we find the match set (and thus its size) and record it in the array 
						of word counts:*/
						if(allowedNrOfMismatchesPairs == 0){

							hitIndexQueryWord = wordSearch_I_binned_I_counts(ptr_binned_I_windowPairs_query[cntPairs].binnedIvector, ptr_binned_I_counts, 0, wordCount, nrOfEntriesForPairs);
							cntMatch = ptr_binned_I_counts[hitIndexQueryWord].count;

						}
						else{

							hitIndexQueryWord = wordSearch_I_binned_I_counts(ptr_binned_I_windowPairs_query[cntPairs].binnedIvector, ptr_binned_I_counts, 0, wordCount+1, nrOfEntriesForPairs);

							//printf("langmodigt venter bysvalen  stadig  stædigt længes  ..., hit?: %d\n", hitIndexQueryWord);

							/*If we got a true hit we just look up the match set size; else we find the match set and
							record its size, so as to gather the info as we go through the queries:*/
							if(hitIndexQueryWord != -1){
								cntMatch = ptr_binned_I_counts[hitIndexQueryWord].count;

								//printf("... and cntMatch %d\n", cntMatch);
								//getchar();

							}
							else{

								/*realloc if nec:*/
								if(wordCount > currentAllocSize_binned_I_counts_binned){

									newAllocSize_binned_I_counts_binned = currentAllocSize_binned_I_counts_binned + ((int) ceil(nrOfDBrecordsLoaded/10.0)); 
									ptr_binned_I_counts = (struct binned_I_counts *) realloc(ptr_binned_I_counts, newAllocSize_binned_I_counts_binned);
									for(j=currentAllocSize_binned_I_counts_binned;j<newAllocSize_binned_I_counts_binned;j++){

										ptr_binned_I_counts[j].binnedIvector = (int *) malloc(nrOfEntriesForPairs*sizeof(int));
										/*init to unlikely values*/
										for(i=0;i<nrOfEntriesForPairs; i++){ ptr_binned_I_counts[j].binnedIvector[i] = 999999; /*high positive to retain lex-order*/}
										ptr_binned_I_counts[j].count = -1;
									}

									currentAllocSize_binned_I_counts_binned = newAllocSize_binned_I_counts_binned;

				
								}

								/*We find the match set:*/
								/*reset*/
								for(i=0; i <nrOfDBrecordsLoaded; i++){ ptr_matchIndicator[i] = 0;}
						
								range.start = 0;
								range.end = nrOfDBrecordsLoaded-1;
								range.nrOfMismatches = 0;
								returnVal = bisectionalSearchSingle_Pairs(ptr_matchIndicator, 0, ptr_newRange, ptr_binned_I_windowPairs_query[cntPairs], 0, range, ptr_binned_I_windowPairs_DB, allowedNrOfMismatchesPairs, nrOfEntriesForPairs);

								//getchar();

								cntMatch = 0;
								for(i=0; i <nrOfDBrecordsLoaded; i++){

									if(ptr_matchIndicator[i] == 1){
										//printf("cntMatch: %d at i: %d\n",cntMatch, i);
										ptr_binned_I_windowPairs_matchSet[cntMatch] =  ptr_binned_I_windowPairs_DB[i];
										//printf("rec nr %d: %d %d %d %d %d\n", i, ptr_binned_I_windowPairs_DB[i].binnedIvector[0], ptr_binned_I_windowPairs_DB[i].binnedIvector[1], ptr_binned_I_windowPairs_DB[i].binnedIvector[2], ptr_binned_I_windowPairs_DB[i].binnedIvector[3], ptr_binned_I_windowPairs_DB[i].binnedIvector[4]);
										//getchar();
										cntMatch +=1;
									}

								}

								/*Record the results and update the order:*/
								for(i=0;i<nrOfEntriesForPairs; i++){ptr_binned_I_counts[wordCount].binnedIvector[i] = ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[i];}
								ptr_binned_I_counts[wordCount].count = cntMatch;
								bubbleSort_binned_I_counts(ptr_binned_I_counts, wordCount+1, nrOfEntriesForPairs);

								/*for(m=0;m<20;m++){
									
									printf("(Sommeren mild ...), 1st letter in word %d in list of hits collected so far is: %d with cntMatch: %d\n", m, ptr_binned_I_counts[m].binnedIvector[0], ptr_binned_I_counts[m].count);
								}
								getchar();
								*/

								wordCount +=1;
								
								if(halted_b == 1){
									getchar();
								}
						
							}
						
						}


						//printf("Win pair %d %d Number of matches: %d\n", k, l, cntMatch);
						//getchar();

						/*We add a tiny constant to the "true" frequency so as to assign a huge score to cases having 
						no match in the data base*/
						rarityFactorPairs = ((float) cntMatch)/nrOfDBrecordsLoaded + 1e-20;
						//printf("Rarityfactor (pct) %lf, cntMatch %d, nrOfDBrecordsLoaded: %d\n", 100*rarityFactor, cntMatch, nrOfDBrecordsLoaded);
						//getchar();
					

						/*re-init the Irarity info, copying over the info from the Iwindows struct*/
						collect_Irarity_windowPairs(ptr_Irarity_windowPairs, &I_windowPairs, cntPairs, k,l, -logbase(rarityFactorPairs, 10) ); //cntPairsConsidered rather than cntPairs?? No cntPairs corr's to the pair k,l

						cntPairsConsidered += 1;

					/////////////////////////

							
					} //end query l-window loop

				} //end query k-window loop

				//printf("Langmodigt venter bysvalen ... cntPairs: %d cntPairsConsidered:%d\n", cntPairs, cntPairsConsidered);

				/*Alternative: do the heap sort and only sum up the first cntPairsConsidered */
				//heapSortRarityScoresPairs(ptr_Irarity_windowPairs, cntPairs);


				/*Sum up the scores obtained (for each window pair) and build a string to be written to file, containing all the window pair info of the scoring*/ 
				score = 0.0;
				strcpy(qWindowsStringAid, "\0");
				strcpy(qWindowsString, "\0");

				if(writeWindowInfo_b == 0){

						strcat(qWindowsString, "Getting scores per window pair was not called. Set writeWindowInfo_b (option: Y) to 1 to change.");

				}

				for(n = 0; n < cntPairs ; n++){ //cntPairsConsidered? no, cntPairs indexes the ptr_Irarity_windowPairs in collect_Irarity_windowPairs

					score += ptr_Irarity_windowPairs[n].rarityScore;


					/*get the window info where there was a contribution to the score:*/
					if(ptr_Irarity_windowPairs[n].rarityScore > 0){

						snprintf(rarityScoreString, 10, "%lf",ptr_Irarity_windowPairs[n].rarityScore); 

						if(writeWindowInfo_b != 0){
							strcat(qWindowsString, "[(");
							snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.windowNr_1);
							strcat(qWindowsString, qWindowsStringAid);
							strcat(qWindowsString, ",");
							snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.segIndices_1[0]);
							strcat(qWindowsString, qWindowsStringAid);
							strcat(qWindowsString, ",");
							snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.segIndices_1[1]);
							strcat(qWindowsString, qWindowsStringAid);
							strcat(qWindowsString, "),(");

							snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.windowNr_2);
							strcat(qWindowsString, qWindowsStringAid);
							strcat(qWindowsString, ",");
							snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.segIndices_2[0]);
							strcat(qWindowsString, qWindowsStringAid);
							strcat(qWindowsString, ",");
							snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.segIndices_2[1]);
							strcat(qWindowsString, qWindowsStringAid);
							strcat(qWindowsString, "),");
							strcat(qWindowsString, rarityScoreString);
							strcat(qWindowsString, "]");
						}

					}

				}

				//For averaging it could seem reasonable to use cntPairsConsidered rather than cntPairs: the back ground data base should 
				//be build with the same constraints: only disjoint pairs or not? and the same threshold on mutuals. However, using the 
				//actual number of windows (of the the same disjoint-ness as in the data base), we get scores that are comparable across
				//diff values of the threshold:
				if(cntPairs == -1){
					score = 0.0;
				}
				else{
					score /= (cntPairs+1); //+1 to avoid dividing by 0 
				}
				
				/*if desired (and a file containing rarityScores is provided), compute the p-value of the rarity factor:*/
				if(rarityScorePValues_b == 1){

					/*If the score is 0 we assign the most "conservative" p-value, and else we compute it
					by looking up the score value on the score-distribution*/
					if(score == 0.0){
						pValuePairs = 1;
					}
					else{

						hitRarityScorePairs = bisectionalSearchValueL(score, ptr_rarityScorePairsDistr, sizeOfRarityScorePairsDistr);

						//printf("(saaledes stadig ganske hemmelig)... score: %lf hitRarityScorePairs:%d sizeOfRarityScorePairsDistr:%d \n", score, hitRarityScorePairs, sizeOfRarityScorePairsDistr);

						/*we add a tiny constant to reveal cases having a score not in that of the range of scores in the data base:*/
						pValuePairs = (1.0 - (double)(hitRarityScorePairs +1)/sizeOfRarityScorePairsDistr) + 1e-20; //hitDistValue is the left end of the interval in which distMin sits; we add 1 to avoid numm disturbance
						//printf("Sommeren mild oedsles  ... p val: %lf\n", pValuePairs);
						//getchar();
					}

				}
						
				/*collect the results:*/
				collect_queryRawScore(ptr_queryRawScore, 1, I_windows_dummy, I_windowPairs, qWindowsString, score, pValuePairs, cntScoredQueries, -9999);

				cntScoredQueries += 1;

				if (write_windowPairs_b != 0){

					returnVal = writeInvariantsOnwindowPairs(ptr_fileNameOutwindowPairs, I_windowPairs, L, order, onlyDisjointPairs_b);
				}

				//Reset! The reason is that the pairs-loop may turn out to be empty -- no pairs fulfull the criteria in there; 
				//if omitting this reset the data from the most recent chain which had a pair fulfilling the criteria will be used in the ptr_queryRawScore!!
				init_ptr_Irarity_windowPairs(&ptr_Irarity_windowPairs, alreadyAllocatedTo_maxNrOfWindows);


				endChain = clock();
				timeChain = ((double)(endChain - startChain))/CLOCKS_PER_SEC;
				//printf("CPU time spend for measures on this file incl load of chain and more:%lf \n", timeChain);
				if (timeChain > max_timeChain){
					max_timeChain = timeChain;
				}


				subStructureCnt += 1; /*will count and index the sub-structures for which we compute the invariants*/
			
			} /*end chainNr loop*/

			
			fclose(ptr_fileIn);

		} /*end MAIN LOOP, fileNr loop*/

		/*sort the obtained scoring results:*/
		heapSortRawScores(ptr_queryRawScore, cntScoredQueries);

		if(write_rarityScoresPairs_b == 1){

			for(i=0;i<cntScoredQueries;i++){

				fprintf(ptr_fileOutScoresPairs, "%s;%s;%s;%s;%s;%d;%d;%s;%lf;%lf;\n",
							DBName,
							ptr_queryRawScore[i].fileName,
							ptr_queryRawScore[i].structureName,
							ptr_queryRawScore[i].chainId,
							ptr_queryRawScore[i].classId,
							ptr_queryRawScore[i].chainLen,
							ptr_queryRawScore[i].nrOfWindows, //query nr of windows
							ptr_queryRawScore[i].windowInfo,
							ptr_queryRawScore[i].score,
							ptr_queryRawScore[i].pValue
							);

			}

			fclose(ptr_fileOutScoresPairs);
			
		}


		/*close file for holding scores*/
		/*if(write_rarityScoresPairs_b ==1){
			fclose(ptr_fileOutScoresPairs);
		}*/

		endComplete = clock();
		//printf("start clocks:%d, end clocks: %d, clocks_per sec:%ld\n", endComplete, startComplete, CLOCKS_PER_SEC);
		timeComplete = ((double) (endComplete - startComplete))/CLOCKS_PER_SEC;
		printf("CPU time spend for all files incl loading data:%lf \n", timeComplete);
		printf("CPU time spend for computations on all files:%lf \n", compTimeComplete);
		printf("Max CPU time spend for measures across all files:%lf \n", max_timeChain);

		printf("Done number of files:%d \n", fileNr);
		printf("Done number of chains/sub-structures:%d \n", subStructureCnt);
		printf("%d too short and %d too long chains were skipped\n", chainsSkippedTooShort,chainsSkippedTooLong);
		printf("%d chains were skipped due to a too long segment (longer than sqrt of %d)\n",chainSkippedSegmTooLong, stdRealSegLength);


		printf("cntMutAboveThreshold: %d\n", cntMutAboveThreshold);

		/*If user wants to score another queries set without reloading the data base: 
		Prompt user a new queries set, else set keepGoing_b = 0 (will then leave the outer while-loop):*/
		keepGoing_b = 0;
		printf("Do you want to run another scan? (yes: 1, no: 0) ");
		fgets(keepGoingStr, 3, stdin);
		printf( "You entered: %s\n", keepGoingStr);
		keepGoing_b = atoi(keepGoingStr);
		//printf("keepGoing_b %d\n", keepGoing_b);
		if(keepGoing_b == 1){

			killroyWasHere_b = 1;

			/*Resets:*/
			strcpy(queriesDirPath, "\0");
			strcpy(queriesName, "\0");
			strcpy(outputPath, "\0");

			printf("Enter a full path to another queries set:");
			fgets(queriesDirPath, 1000, stdin);
			queriesDirPath[strcspn(queriesDirPath, "\n")] = 0; //removes trailing \n if any!
			printf("Queries path: %s", queriesDirPath );
			printf("Enter a name for this queries set:");
			fgets(queriesName, 100, stdin);
			queriesName[strcspn(queriesName, "\n")] = 0; //removes trailing \n if any!

			printf("Enter a full path to where the output should be placed:");
			fgets(outputPath, 1000, stdin);
			//strcat(outputPath,"\0" );
			outputPath[strcspn(outputPath, "\n")] = 0; //removes trainling \n if any!
			printf("Output path: %s", outputPath );

			printf("Is the new queries set from CATH, SCOP or none of these: (0,1 or 2)?");
			fgets(dataTypeIndicatorStr, 3, stdin);
			dataTypeIndicator = atoi(dataTypeIndicatorStr);
			use_cath_b = 0;
			use_scop_b = 0;
			if(dataTypeIndicator == 0){
				use_cath_b = 1;
			}
			else if(dataTypeIndicator == 1){
				use_scop_b = 1;
			}


			/*Reset file names:*/
		    strcpy(fileNameOutScoresPairs1, "");
			strcpy(fileNameOutScoresPairs1, "/RarityScan1_ScoresPairs_windowslgth_");
			strcpy(fileNameOutwindowPairs1,"/RarityScan1_Invariants_Pairs_windowlgth_");	
			strcpy(fileNameOut2, "_order_");

			strcpy(fileNameOut3,"\0");

			strcpy(fileNameOutwindowPairs, "\0");
			strcpy(fileNameOutScoresPairs, "\0");

			/*Other resets*/
			reinit_queryRawScore(ptr_queryRawScore, cntScoredQueries);

			subStructureCnt = 0;
			cntMutAboveThreshold = 0;
			cntScoredQueries = 0; 
			fileNr = 0;
			chainNr = 0;
			maxChainLen = 0;

			chainsSkippedTooShort = 0;
			chainsSkippedTooLong = 0;
			chainSkippedSegmTooLong = 0;			 

		}


	}


	/*free the memory allocated to this query set*/
	free_I_measures(I_measures, order, full_b, alreadyAllocatedTo_chainLen);
	free(ptr_segment);	
	free(ptr_chain);
	free(dirContent.ptr_dirList); //this is allocated in the ListDirectoryContents2 fct 
	free(dirContent.ptr_fileNameList); //this is allocated in the ListDirectoryContents2 fct 
	//Resets:
	dirContent.numberOfFiles = 0;
	subDirCnt = 0;
	//More freeing:
	for(i=0;i< subStructureCnt; i++){
		free(ptr_queryRawScore[i].fileName);
		free(ptr_queryRawScore[i].structureName);
		free(ptr_queryRawScore[i].classId);
		free(ptr_queryRawScore[i].chainId);
		free(ptr_queryRawScore[i].windowInfo);
	}
	free(ptr_queryRawScore);

	
	returnVal = 1;

	return returnVal;
}


/*
Function for deciding if a structure (query) stands out as "rare" compared to a background (based on a set of PDB
files). The scoring is done by means of a "cross entropy" type. We often use the term "data base" for the back ground. The code 
allows running a set of queries against the data base; in fact it is set up so that one may run a series of query sets against the 
same back ground (withou reloading the back ground). 

!Similar to rawRarity1; the data base is though here "dictionaried" based on the Gauss numbers of the windows and not the window pairs!
It is (therefore) possible to run a scan either in single window mode or have the pairs considered too.

The usage is:
1) Create the data base: this is done by running the main function here (SAonGISA_main_*) in the makeDB flavour. This wraps a call 
to the GI_windows function of the GISA_v*_unix-code. To run the scan in "window pairs mode" a data base for the pairs of 
windows must be created (run SAonGISA_main_* in flavour makeDB); this must not only contain results for the disjoint 
windows, but rather all combinations (since the single window results will be taken from this pairs data base). To run in "single 
window mode" a data base for the single windows results must be created (as for the pairs by running SAonGISA_main_* in flavour makeDB, but
here of course calling for having the single window values computed and stored). In sjort: in either case have the data base created for the 
desired Gauss Integrals invariants. For this set the order according to the set of invariants wanted (ie set order to 1 if only the writhe 
and average crossing number are desired, else 2); for the pairs-based scan it is recommended only to use order 1 (nrOfEntriesForPairs can 
be set to 1, while nrOfEntries can be left higher so as to allow to have higher order Gauss numbers included in the single windows matching). 
Set windowing parameters as desired (e.g. windowLength 20, stepSize 4).
2) Now run the rarity scan with rawRarity2: the search method is aimed at detecting the frequency of data base windows (and, possibly, 
pairs) with Gauss numbers in "epsilon-distance" to those of a given query window (pairs); for a given query structure these 
frequencies are collected and a score is obtained by accumulation. The "epsilon" is set indirectly through binning of the invariant 
values.  

To obtain the p-values corresponding to the obtained scores output, it is necessary to first run rawRarity1 with query set = 
data base set and with rarityScorePValues_b = 0 (and else leave the settings). When done run rawRarity2 for the desired 
query set now with rarityScorePValues_b = 1.

It is possible to set a threshold filtering out pairs with abs mutual value lower than this cutoff level. For speed it
is (very) desirable to set this threshold quite high; when the scan is only based on the writhe (and no other invariants) the 
threshold can be set to  e.g. 10 (and leave out normalization of the invariant values by setting the normalize_b parameter at 0). 
When running with more invariants (e.g. also the average crossing number) the values must be normalized (will be done to get mean 0 
and st dev 1) and the threshold must then be set differently (maybe to 0.5), but this may need a little experimenting (e.g. have a look 
at the actual normalized invariant values -- these values can be written out).

An outline of the method:
a) The data base is converted to a "list of words": each data base element is a tuple/array of Gauss number (in number as many as
many as there are invariants of the chosen order, ie 1 or 2); each of these tuples is converted by binning the Gauss numbers into a tuple of integers
(ie a "word"); after sorting these words lexicographically the data base has the guise of a dictionary (though with probably many words
repeated). As opposed to rawRarity1 this "dictionarying" is based on the Gauss numbers for the windows and not the pairs.
b) A given query is similarly translated into a set of words (one for each window); the windows are now looped through 
and each is looked up in the data base dictionary (for a pairs-based scan see more right below). This gives a count 
of the number of matches, cntMatch, for each window. The score is now the sum

Score = - Sum over query windows log(cntMatch(window)/#data base) /#windows in query

which equals

Score = - Sum over words (w)  ( #windows in query of word w) /#windows in query * log( #db-windows of word w/#data base) 

If we write p_q (w ) = ( #pairs of word w) /#windows in query and p_db (w) = #db-windows of word w/#data base we have

Score = - Sum over words (w) p_q(w) * log(p_db(w)),

a cross entropy that is, and also showing the resemblance to the Kullback-Leibler relative entropy

KL = Sum over words (w) p_q(w) * log(p_q(w)) / p_db(w))

ie only the "q-idiosyncratic" term "Sum over words(w) p_q(w) * log(p_q(w))" is disregarded. 


If a scan based on the pairs is also wanted, the pairs of windows in a query will subsequently be looped through (ie after 
all windows of the given query have been looped through and the matches are found and recorded); it is possible 
and desirable for speed to set a threshold so that only pairs having a mutual writhe (or mutual invariant) value above this 
threshold are considered). The scoring is verbatim the same as in the single window case, simply replace window by window pair (or see 
rawRarity1).

The scoring method scans for structures having a "distribution of words" significantly different from that found in the background
distribution. One can also think of this as a way of determining whether a query has an unusual set of fragments/fragment pairs as compared 
to the background. The threshold allows to focus on occurence of "more rare frament pairs" ("more rare words" in the pairs-based matching).

Final scores are written to a file, in decreasing order ( ie with lowest p-value at the top). 
*/
int rawRarity2( int halted_b, char rarityScoreDistrFilePath[1000], char rarityScoreDistrName[100], char rarityScorePairsDistrFilePath[1000], 
			   char rarityScorePairsDistrName[100], char DBName[100], char queriesName[100], 
			   char queriesDirPath[200], char DBResultsFilePath[1000],  char DBResultsPairsFileName[1000],
			   char fileNameQueryChains[2000], 
			   double thresholdMutual, int onlyDisjointPairs_b, 
			   int rarityScorePValues_b, int use_scop_b, 
			   int use_cath_b, char *cathListPath, int loadFromSubDirs_b, 
			   int normalize_b, char *outputPath, int maxChainLength, 
			   int nrOfEntries, int nrOfEntriesForPairs,  int binningType, 
			   int nrOfBins, double compensationForRounding, int allowedNrOfMismatches, 
			   int allowedNrOfMismatchesPairs, int nrOfWindowsForMatch,  int chunkSize, 
			   int order, int incl_abs_b, int full_b, 
			   int write_matchWindows_b, int matchWindowPairs_b, 
			   int write_matchWindowPairs_b, int windowCoveringType, int windowLength, 
			   int stepSize, int write_windows_b, int write_windowPairs_b, 
			   int write_rarityScores_b, int write_rarityScoresPairs_b, int write_chains_b, 
			   int writeWindowInfo_b, int print_b, int print_basic_b){

	int cntRuns = 0; //for testing code
	
	int use_lex = 0; /*keep this!*/

	int returnVal = 0;

	int split_b = 0;

	double rarityFactor; //to measure how rare a GI-vector is (in the first filtering)
	double rarityFactorPairs; //to measure how rare a GI-vector for pairs is (in the first filtering)

	//int windowLength = 16;
	//int write_windows_b = 1; /*for writing out I-values on windows for queries*/
	//int write_windowPairs_b = 0; /*at most one of write_windows_b and write_windowPairs_b should be set to 1; corr to get_windos_b resp. get_windowPairs_b*/
	/*for matching query to DB (set of structures in the set directory:*/
	int k,l, queryNrOfWindows, queryNrOfWindowPairs;
	int u;
	double maxMutVal, mutVal;
	//int nrOfWindowsForMatch = 3;
	int cntMutAboveThreshold = 0;
	
	double **ptr_query, **ptr_match; 
	double ***ptr_query_pairs, ***ptr_match_pairs; 

	//int strLengthCutOff = 10; /*chains shorter than this int will not be included in computation of invariants*/

	struct dirContent dirContent;
	int subDirCnt = 0;
	int numberOfFiles = 0;
	int subStructureCnt = 0; /*some structures are multi-mers; we compute the invariants for each sub-structure (each mono-mer)*/
	int maxNrOfWindows = 0; /*max nr of windows across loaded structures*/
	int maxNrOfWindowPairs = 0; /*max nr of window pairs across loaded structures*/
	const int numberOfPerturbations = 1; /*only there since needed when writing out the final values (preparation for perturbations code)*/
	const int pertNo = 0; /*ditto*/
	char **ptr_dirList;
	char **ptr_fileNameList;

	/*for "dynamically" generating the file names for output:*/
	char orderStr[10] = "-1"; //place holder value
	char inclAbsStr[10] = "-1"; //value is placeholder
	char nrInvsForMatchStr[10] = "-1"; //place holder value
	char nrWinsForMatchStr[10] = "-1"; //place holder value
	char allowedNrOfMismatchesStr[10] = "-1"; //place holder value; 
	char nrInvsForPairsMatchStr[10] = "-1"; //place holder value
	char allowedNrOfMismatchesPairsStr[10] = "-1"; //place holder value; 
	char binTypeStr[10] = "-1"; //place holder value
	char nrBinsStr[10] = "-1"; //place holder value
	char windowLengthStr[10] = "-1"; //place holder value
	char windowStepSizeStr[10] = "-1"; //value is placeholder
	char windowCovTypeStr[10] = "-1"; //value is placeholder
	char topNStr[10] = "-1"; //value is placeholder
	char thresholdMutualStr[10] = "-1"; //value is placeholder

	char fileNameOutScores1[1000] = "/RarityScan2_Scores_windowslgth_";
	char fileNameOutScoresPairs1[1000] = "/RarityScan2_ScoresPairs_windowslgth_";
	char fileNameOutWindows1[1000] = "/RarityScan2_Invariants_windowlgth_";
	char fileNameOutwindowPairs1[1000] = "/RarityScan2_Invariants_Pairs_windowlgth_";	
	char fileNameOut2[1000] = "_order_"; //; //; //"_minmax_";//"_incl_abs_"; //"_"; //

	char fileNameOut3[1000] = "\0";  //"_testQuery_top8000_5invs_norm.txt"; 
	char extensionName[50] = ".txt";

	char fileNameOutWindows[2000] = "\0";
	char fileNameOutwindowPairs[2000] = "\0";
	char fileNameOutScores[2000] = "\0";
	char fileNameOutScoresPairs[2000] = "\0";

	char *ptr_fileNameOutWindows; /*pointer to file for holding measures' values on sub-chains of pre-set windowlength*/
	FILE *ptr_fileOutWindows;
	char *ptr_fileNameOutwindowPairs; /*pointer to file for holding measures' values on pairs of sub-chains of pre-set windowlength*/
	FILE *ptr_fileOutwindowPairs;
	char *ptr_fileNameOutScores;
	FILE *ptr_fileOutScores;
	char *ptr_fileNameOutScoresPairs;
	FILE *ptr_fileOutScoresPairs;


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
	int structureNameLength = 7; /*string-length of structure name in DB, e.g SCOP, CATh*/
	char *currentStructureName;
	char *structureNameInit = "NN";
	char *currentClassId;
	char *classIdInit = "NN";
	char chainIdInit[2] = ">";
	int chainsSkippedTooShort = 0;
	int chainsSkippedTooLong = 0;
	int chainSkippedSegmTooLong = 0;
	int segTooLong_b = 0;
	//int resNr = 0; /*residue number*/
	//int resNrPrior = -999;
	int chainLen = 0; /*length of C-alpha chain*/
	int maxChainLen =0; /*max length of chain among all loaded*/
	int L = 0; /*length of segment chain = chainLen -1*/
	int simplexCnt = 0;/*size of simplex*/
	int maxSimplexCnt = 0; /*size of simplex corr to max chain length*/

	struct cAlpha * ptr_chain = NULL; /*for holding C-alpha chain */
	struct segment *ptr_segment = NULL; /* for holding chain of segments, ith C-alpha to jth C-alpha*/
	struct segment segCoords;	 

	int m = 0, n = 0, mMin, nMin;
	int i = 0;
	int j = 0;
	int cnt = 0;
	int fileNr = 0;
	int chainNr = 0;

	/* for holding meausure values (i.e. the core of the final output); 
	content to be received from aggr function:*/
	struct I_ptr I_measures;
	double ***I_measures_raw; /*triple pointer to replace the struct I_ptr use*/ 
	int alreadyAllocatedTo_chainLen  = 0;
	int alreadyAllocatedTo_maxNrOfWindows = 0;
	int alreadyAllocatedTo_maxNrOfWindowPairs = 0;
	int alreadyAllocatedTo_subStructureCnt = 0;

	/*pointer to hold the results (that this is a ptr of a ptr is due to that the write-out can also
	handle the perturbation case; we use pertNo = 0 here):*/
	struct I_values **ptr_I_values;

	/*struct for holding invariants on windows resp. window pairs for the query and for each DB entry:*/
	struct I_windows_ptr I_windows; //for query
	struct I_windowPairs_ptr I_windowPairs; //for query
	double **bins;
	double **binsPairs;
	struct binned_I_windows *ptr_binned_I_windows_query; //for query, binned values
	struct binned_I_windowPairs *ptr_binned_I_windowPairs_query; //for query, binned values
	struct I_windows_ptr db_I_windows;
	struct I_windowPairs_ptr db_I_windowPairs;
	int nrOfDBstructuresLoaded = 0; 
	int nrOfDBrecordsLoaded = 0;
	int nrOfDBstructuresLoaded_pairsLoad = 0; 
	int nrOfDBrecordsLoaded_pairsLoad = 0;
	int maxNrOfWindowsDB = 0;
	int maxNrOfWindowPairsDB = 0;	
	/*pointer to hold DB results*/
	struct db_I_windows_order2_ptr *ptr_dbResultWindows;
	struct db_I_windowPairs_order2_ptr *ptr_dbResultWindowPairs;
	struct I_windows_order2_meanStddev I_windows_order2_meanStddev_DB;
	struct I_windowPairs_order2_meanStddev I_windowPairs_order2_meanStddev_DB;

	/*for searching*/
	int hitLR[2];
	int *intArray2; //convenience array
	int *intArray3; //convenience array
	struct binned_I_windows *ptr_binned_I_windows_DB;
	int matchSetNrOfStructures = 1000;
	struct binned_I_windows *ptr_binned_I_windows_matchSet;
	struct binned_I_windows **ptr_ptr_binned_I_windows_matchSet;
	int db_fileNr;//for window pairs searching
	int db_wNr1;//for window pairs searching
	int db_wNr2;//for window pairs searching
	int cntMatch;
	int *ptr_cntMatch; /*to hold the cntMatch for the series of windows covering a query*/
	int *binnedIvector_matchPair; //for match, binned values
	struct Irarity_windows *ptr_Irarity_windows;
	struct Irarity_windows Irarity_windows_temp;
	struct Irarity_windowPairs *ptr_Irarity_windowPairs;

	int nrOfMismatchesThisPair = 0;
	int qualified_b = 1;
	int cntMatch_12 = 0;
	int cntPairs = 0;
	int cntPairsConsidered = 0;
	float diff;

	int *ptr_matchIndicator;
	struct matchRange **ptr_newRange;
	struct matchRange range; 
	//struct matchCandidate *ptr_matchCandidates; /*for holding the topN candidate single window matches to be written out (if desired)*/
	//struct matchScore *ptr_matchScores; /*for holding final scores for each query*/ 

	double score = 0.0;
	int cntScoredQueries = 0;
	double *ptr_rarityScoreDistr; //sorted list of dist-values (to be had from db-vs-db run)
	FILE *ptr_fileRarityScoreDistr;
	int sizeOfRarityScoreDistr = -1;
	int hitRarityScore;
	struct queryRawScore *ptr_queryRawScore;
	double pValue = -1;
	char * qWindowsString;
	char * qWindowsStringAid;

	double scorePairs = 0.0;
	int cntScoredQueriesPairs = 0; 
	double *ptr_rarityScorePairsDistr; //sorted list of dist-values (to be had from ...)
	FILE *ptr_fileRarityScorePairsDistr;
	int sizeOfRarityScorePairsDistr = -1;
	int hitRarityScorePairs;
	struct queryRawScore *ptr_queryRawScorePairs;
	double pValuePairs = -1;
	char * qWindowsStringPairs;
	char * qWindowsStringAidPairs;
	char * rarityScoreString;



	/*for closed loops finding, but only place holders here:*/
	int closed_loops_b = 0; //must be kept like this!
	struct twoSegmentIndex *ptr_closedLoopInd;
	int initValNrOfClosedLoops = 0;
	int nrOfClosedLoops = 0;
	int cntClosedLoop = 0;

	FILE *fp; //dummy file ptr used when checking file names

	/*params to control the execution*/
	int keepGoing_b = 1; 
	char keepGoingStr[10] = "-1"; //place holder value
	int dbLoaded_b = 0;
	char dataTypeIndicatorStr[3] = "-1"; //place holder value
	int dataTypeIndicator = -1;
	int killroyWasHere_b = 0;


	/*Declarations stop here*/

	/*The code for running one set of queries against the data base is wrapped in an outer while-loop, which allows the user to
	enter a new queries set, without having to reload the data base*/
	while(keepGoing_b == 1){

		/****************************************************************************************************************/
		/*Generation of output file names:*/
		/****************************************************************************************************************/

		/*we first build a string (fileNameOut3) shared by most file names, which contains various params values;
		fileNameOut3 will be concatenated into the file name strings below:*/

		//printf("(saaledes stadig ganske hemmelig) ?\n");

		sprintf(nrWinsForMatchStr, "%d", nrOfWindowsForMatch);
		strcat(fileNameOut3, nrWinsForMatchStr);
		strcat(fileNameOut3, "wins_");
		if(normalize_b == 1){
			strcat(fileNameOut3, "norm_");
		}
		else{
			strcat(fileNameOut3, "notnorm_");
		}
		sprintf(nrInvsForMatchStr, "%d", nrOfEntries);
		strcat(fileNameOut3, nrInvsForMatchStr);
		strcat(fileNameOut3, "invs_");
		sprintf(allowedNrOfMismatchesStr, "%d", allowedNrOfMismatches);
		strcat(fileNameOut3, allowedNrOfMismatchesStr);
		strcat(fileNameOut3, "mms_");
		if(matchWindowPairs_b == 1){
			sprintf(nrInvsForPairsMatchStr, "%d", nrOfEntriesForPairs);
			strcat(fileNameOut3, nrInvsForPairsMatchStr);
			strcat(fileNameOut3, "invsPairs_");
			sprintf(allowedNrOfMismatchesPairsStr, "%d", allowedNrOfMismatchesPairs);
			strcat(fileNameOut3, allowedNrOfMismatchesPairsStr);
			strcat(fileNameOut3, "mmsPairs_");
		}
		sprintf(nrBinsStr, "%d", nrOfBins);
		strcat(fileNameOut3, nrBinsStr);
		strcat(fileNameOut3, "bins");
		sprintf(binTypeStr, "%d", binningType);
		strcat(fileNameOut3, binTypeStr);
		sprintf(windowCovTypeStr, "%d", windowCoveringType);
		strcat(fileNameOut3,windowCovTypeStr);
		strcat(fileNameOut3, "_winCovType");
		strcat(fileNameOut3, windowCovTypeStr);
		if(matchWindowPairs_b == 1){
			strcat(fileNameOut3, "_threshMut");
			sprintf(thresholdMutualStr, "%0.2f", thresholdMutual);
			strcat(fileNameOut3, thresholdMutualStr);
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
			strcat(fileNameOutWindows, queriesName);
			strcat(fileNameOutWindows, "_");
			strcat(fileNameOutWindows, DBName);
			strcat(fileNameOutWindows, "_");
			strcat(fileNameOutWindows, fileNameOut3);
			strcat(fileNameOutWindows, extensionName);
			printf("Results for windows will be written to: %s\n", fileNameOutWindows);
			
			ptr_fileNameOutWindows = fileNameOutWindows;
			

			/*clear the file by opening it in w-mode and then closing it again*/
			fp= fopen(ptr_fileNameOutWindows, "w");
			if(fp == NULL){
				printf("Sorry, file name %s is invalid. Could the path be wrong? A blank space? --- execution stops here.", ptr_fileNameOutWindows);
				getchar();
				exit(EXIT_FAILURE);
			}
			else{
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
			strcat(fileNameOutwindowPairs, queriesName);
			strcat(fileNameOutwindowPairs, "_");
			strcat(fileNameOutwindowPairs, DBName);
			strcat(fileNameOutwindowPairs, "_");
			strcat(fileNameOutwindowPairs, fileNameOut3);
			strcat(fileNameOutwindowPairs, extensionName);
			printf("Results for window pairs will be written to: %s\n", fileNameOutwindowPairs);
			 
			ptr_fileNameOutwindowPairs = fileNameOutwindowPairs;
			
			/*clear the file by opening it in w-mode and then closing it again*/
			fp= fopen(ptr_fileNameOutwindowPairs, "w");
			if(fp == NULL){
				printf("Sorry, file name %s is invalid. Could the path be wrong? A blank space? --- execution stops here.", ptr_fileNameOutwindowPairs);
				getchar();
				exit(EXIT_FAILURE);
			}
			else{
				fclose(fp);
			}
		

		}



		
		if (write_rarityScores_b !=0){

			/*generate file name for writing out match values*/
			sprintf(orderStr, "%d", order);
			strcat(fileNameOutScores, outputPath);
			strcat(fileNameOutScores, fileNameOutScores1);
			sprintf(windowLengthStr, "%d", windowLength);
			strcat(fileNameOutScores, windowLengthStr);
			strcat(fileNameOutScores, "_");
			sprintf(windowStepSizeStr, "%d", stepSize);
			strcat(fileNameOutScores, windowStepSizeStr);
			strcat(fileNameOutScores, fileNameOut2);
			strcat(fileNameOutScores, orderStr);
			strcat(fileNameOutScores, "_");
			sprintf(inclAbsStr, "%d", incl_abs_b);
			strcat(fileNameOutScores, inclAbsStr);
			strcat(fileNameOutScores, "_");
			strcat(fileNameOutScores, queriesName);
			strcat(fileNameOutScores, "_");
			strcat(fileNameOutScores, DBName);
			strcat(fileNameOutScores, "_");
			strcat(fileNameOutScores, fileNameOut3);
	
			/*when getting the rarity score p-values, inlcude info on which distribution was used:*/
			if(rarityScorePValues_b == 1){
				strcat(fileNameOutScores, rarityScoreDistrName);
			}

			
			strcat(fileNameOutScores, extensionName);
			printf("Rarity scores will be written to: %s\n", fileNameOutScores);

			 
			ptr_fileNameOutScores = fileNameOutScores;
			
			/*clear the file by opening it in w-mode and then closing it again*/
			ptr_fileOutScores = fopen(ptr_fileNameOutScores, "w");

			fprintf(ptr_fileOutScores, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;\n",
				"DB name",
				"Query file",
				"structureName",
				"chainId",
				"classId",
				"chainLen",
				"nrOfWindows",
				"windowInfo",
				"rarityScore",
				"pValue"
				);  
		

		}

		if (write_rarityScoresPairs_b !=0){

			/*generate file name for writing out match values*/
			sprintf(orderStr, "%d", order);
			strcat(fileNameOutScoresPairs, outputPath);
			strcat(fileNameOutScoresPairs, fileNameOutScoresPairs1);
			sprintf(windowLengthStr, "%d", windowLength);
			strcat(fileNameOutScoresPairs, windowLengthStr);
			strcat(fileNameOutScoresPairs, "_");
			sprintf(windowStepSizeStr, "%d", stepSize);
			strcat(fileNameOutScoresPairs, windowStepSizeStr);
			strcat(fileNameOutScoresPairs, fileNameOut2);
			strcat(fileNameOutScoresPairs, orderStr);
			strcat(fileNameOutScoresPairs, "_");
			sprintf(inclAbsStr, "%d", incl_abs_b);
			strcat(fileNameOutScoresPairs, inclAbsStr);
			strcat(fileNameOutScoresPairs, "_");
			strcat(fileNameOutScoresPairs, queriesName);
			strcat(fileNameOutScoresPairs, "_");
			strcat(fileNameOutScoresPairs, DBName);
			strcat(fileNameOutScoresPairs, "_");
			strcat(fileNameOutScoresPairs, fileNameOut3);
	
			/*when getting the rarity score p-values, inlcude info on which distribution was used:*/
			if(rarityScorePValues_b == 1){
				strcat(fileNameOutScoresPairs, rarityScoreDistrName);
			}

			
			strcat(fileNameOutScoresPairs, extensionName);
			printf("Rarity scores will be written to: %s\n", fileNameOutScoresPairs);

			 
			ptr_fileNameOutScoresPairs = fileNameOutScoresPairs;
			
			/*clear the file by opening it in w-mode and then closing it again*/
			ptr_fileOutScoresPairs = fopen(ptr_fileNameOutScoresPairs, "w");

			fprintf(ptr_fileOutScoresPairs, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;\n",
				"DB name",
				"Query file",
				"structureName",
				"chainId",
				"classId",
				"chainLen",
				"nrOfWindows",
				"windowInfo",
				"rarityScore",
				"pValue"
				);  

			/*fprintf(ptr_fileOutScoresPairs, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;\n",
				"DB name",
				"Query file",
				"structureName",
				"chainId",
				"classId",
				"chainLen",
				"nrOfWindows",
				"qWinNr_1",
				"i1",
				"i2",
				"qWinNr_2",
				"j1",
				"j2",
				"rarityScore",
				"pValue"
				);  */
		

		}

		if(halted_b == 1){
			printf("If the output paths are fine, press enter to continue.\n");
			getchar();
		}

		/****************************************************************************************************************/
		/*Get content of the desired directory: number of files and a list of the file names:*/
		/****************************************************************************************************************/

		if(loadFromSubDirs_b ==0){
			dirContent = ListDirectoryContents(queriesDirPath);
		}
		else{
			returnVal = ListDirectoryContents2(queriesDirPath, &dirContent, &subDirCnt);
		}
		numberOfFiles = dirContent.numberOfFiles; 
		ptr_dirList = dirContent.ptr_dirList;
		ptr_fileNameList = dirContent.ptr_fileNameList;
		printf("1st file: %s\n",ptr_fileNameList[0] );

		printf("Number of query files in directory:%d \n",numberOfFiles);

		if(halted_b == 1){
			printf("Press enter to continue\n");
			getchar();
		}

		startComplete = clock();

		/****************************************************************************************************************/
		/*loop through the list of filenames in the directory to find chain info, eg the max chain length (for setting
		the size in the (global) mem alloc to I_measures's ptr's). First though we make a bulk allocation
		to the chain-info keeping pointer, using the pre-set maxNrOfChains:*/
		/****************************************************************************************************************/
		
		if(killroyWasHere_b == 0){
			returnVal = alloc_init_chainInStr(&chainInStr, maxNrOfChains);
		}


		//CALLING ALLOC_INIT_CHAININSTR IS THE SAME AS THIS, RIGHT?:

		/* chainInStr.numberOfChains = maxNrOfChains;
		chainInStr.ptr_chainInStr = (struct chainInfo *) malloc (maxNrOfChains*sizeof(struct chainInfo));
		//init
		for(i = 0; i< maxNrOfChains; i++){
			chainInStr.ptr_chainInStr[i].chainId = (char *) malloc(sizeof(char)*chainIdCharSize);
			chainInStr.ptr_chainInStr[i].structureName = (char *) malloc(sizeof(char)*structureNameCharSize);
			chainInStr.ptr_chainInStr[i].classId = (char *) malloc(sizeof(char)*classIdCharSize);
			chainInStr.ptr_chainInStr[i].classIdConv = (int *) malloc(classifierSize*sizeof(int));

			strcpy(chainInStr.ptr_chainInStr[i].chainId, chainIdInit); //init to unlikely value
			strcpy(chainInStr.ptr_chainInStr[i].structureName, structureNameInit);
			strcpy(chainInStr.ptr_chainInStr[i].classId, classIdInit);
			//printf("Chain Id initial value:%c\n", *chainInStr.ptr_chainInStr[i].chainId);
			//printf("Class Id initial value:%c\n", *chainInStr.ptr_chainInStr[i].classId);
			chainInStr.ptr_chainInStr[i].chainLength = -1;

			for(j=0; j< classifierSize; j++){chainInStr.ptr_chainInStr[i].classIdConv[j] = -1;}

		} */
		

		/*if CATH data, read in CATH domain list. Note: we abuse the struct chainInStructure to accomodate
		domains rather than chains -- in CATH each chain may be split in several domains*/ 
		if(use_cath_b == 1){

			if(convToCLF20_b == 1){
				structureNameLength = lengthCLFformat20;
			}
			else{
				structureNameLength = 6;
			}

			/*first allocate a block of memory (if not large enough it will be extended when exectuing readCATHDomainList)*/
			ptr_cathDomain = (struct chainInfo *) malloc (maxNrOfDomains*sizeof(struct chainInfo));
			/*init*/
			for(i = 0; i< maxNrOfDomains; i++){
				ptr_cathDomain[i].chainId = (char *) malloc(sizeof(char)*2);
				ptr_cathDomain[i].structureName = (char *) malloc(sizeof(char)*(structureNameLength+2)); //to hold domain name
				ptr_cathDomain[i].classId = (char *) malloc(sizeof(char)*10);
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

			if (!ptr_cath_domain_list_file){
				printf("Sorry, the file containing the positive-list of CATH domains was not found (%s was provided). Enter a valid file name and run again.\n", cathListPath);
				return 0;
			}

			/*read in contents of CATH domain list file*/
			nrOfCATHDomainsInList = readCATHDomainList(ptr_cath_domain_list_file, &ptr_cathDomain, convToCLF20_b );

			printf("nr of cath doms: %d\n", nrOfCATHDomainsInList);

			/*lex-sort the ptr_cath_domain_list_file on the domain name entry. Below we
			need to look up each file for which we have PDB-coords in this list. In other
			words, this list serves as a "positive list" of domains, for which to compute
			the invariants.*/
			heapSortCATHDomain(ptr_cathDomain, nrOfCATHDomainsInList, structureNameLength,0);

			/*for(i=0;i <nrOfCATHDomainsInList;i++){
				if(i%100 == 0){
					for(j=i;j <i+100;j++){
						printf("dom %d: %s\n", j, ptr_cathDomain[j].structureName);
					}
					getchar();
				}
			}*/
		}


		printf("Fetching chain info for data base structures ...\n");
		for(fileNr = 0; fileNr <numberOfFiles; fileNr ++){

			if (fileNr%100 == 0){printf(".. now fetched the chain info for %d structures\n", fileNr);}

			if (print_b == 1){
				printf("File name:%s path:%s fileNr:%d\n",ptr_fileNameList[fileNr], ptr_dirList[fileNr], fileNr);
			}

			/*if using CATH data: check if the file, ptr_fileNameList[fileNr], is in the
			domain list specified (and loaded and lex-sorted above)*/
			if(use_cath_b == 1){
				hitIndexDomain = bisectionalSearchCATHDomain(ptr_fileNameList[fileNr], ptr_cathDomain, nrOfCATHDomainsInList, lengthCLFformat20);
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
				continue;
			}

			/*get number of chains, their id's and lengths for this file:*/
			//chainInStr = readPDBChainStructure(ptr_fileIn);
			if(use_scop_b ==1){
				readPDBChainStructureSCOP(ptr_fileIn, &chainInStr );
			}
			else if(use_cath_b ==1){
				readPDBDomainStructureCATH(ptr_fileIn, &chainInStr );//will read chain length and chain id into chainInStr
			}
			else{
				readPDBChainStructure2(ptr_fileIn, &chainInStr );
			}


			fclose(ptr_fileIn);  

			for (i=0;i <= chainInStr.numberOfChains -1;i++){
				if(chainInStr.ptr_chainInStr[i].chainLength <= maxChainLength){ //we disregard very long chains for RAM availability reasons
					maxChainLen = max(chainInStr.ptr_chainInStr[i].chainLength,maxChainLen);
					if(chainInStr.ptr_chainInStr[i].chainLength >= strLengthCutOff){
										
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

						}

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

		/****************************************************************************************************************/
		/*MEMORY ALLOCATION*/
		/****************************************************************************************************************/

		/*We allocate memory for ptr holding the I-values and a ptr for holding the 3d-coord's, both for the longest chain*/
		returnVal = alloc_init_I_measures(&I_measures, order, full_b, chainNr, maxChainLen, closed_loops_b, alreadyAllocatedTo_chainLen, &ptr_closedLoopInd);

		if(alreadyAllocatedTo_chainLen == 0){ //no allocation was done yet

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
			for ( n=0; n <= maxChainLen -2 ; n++ ){
					segCoords.s1 = ptr_chain[n].coords;
					segCoords.s2 = ptr_chain[n+1].coords;

					ptr_segment[n] = segCoords;

			}

		}
		else if(maxChainLen > alreadyAllocatedTo_chainLen){ //reallocation necc

			ptr_chain = realloc (ptr_chain, maxChainLen*sizeof(struct cAlpha));

			/*(re)init*/
			for ( n=0; n <= maxChainLen -1 ; n++ ){
				
				ptr_chain[n].coords.x = 0.0;
				ptr_chain[n].coords.y = 0.0;
				ptr_chain[n].coords.z = 0.0;

				ptr_chain[n].residueNr = -123;


			}

			ptr_segment = realloc (ptr_segment, maxChainLen*sizeof(struct segment));
			/*init*/
			segCoords.s1.x = 0.0;
			segCoords.s1.y = 0.0;
			segCoords.s1.z = 0.0;
			segCoords.s2.x = 0.0;
			segCoords.s2.y = 0.0;
			segCoords.s2.z = 0.0;
			for ( n=0; n <= maxChainLen -2 ; n++ ){

					segCoords.s1 = ptr_chain[n].coords;
					segCoords.s2 = ptr_chain[n+1].coords;

					ptr_segment[n] = segCoords;

			}


		}
		alreadyAllocatedTo_chainLen = max(alreadyAllocatedTo_chainLen, maxChainLen);


		/****************************************************************************************************************/
		/* Allocate mem to some ptr's for query: to hold I-values on window pairs, the binned version and for the scores	
		/****************************************************************************************************************/


		/*Search on single windows, and possibly pairs: Allocate mem to some ptr's, for query, for read-out and for match; read in the DB results along the line*/

		/*For single windows*/
		maxNrOfWindows =  (int) ceil((double) 1.1*(maxChainLen-windowLength)/stepSize); /*this should just suffice for the longest chain -- we're using overlapping windows*/
		printf("maxNrOfWindows %d, alreadyAllocatedTo_maxNrOfWindows %d \n", maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows);
		/*allocate memory to struct I_windows and initialize*/
		//getchar();
      	maxNrOfWindows = alloc_init_I_windows_ptr(&I_windows, order, full_b, maxChainLen, windowLength, stepSize, maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows);

		/*for holding rarity scores to be written out:*/
		returnVal = alloc_init_ptr_Irarity_windows(&ptr_Irarity_windows, order, full_b, maxChainLen, windowLength, stepSize, maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows);


		/*allocate and initialize a ptr to hold the scores results for each query*/
		alloc_init_queryRawScore(&ptr_queryRawScore, subStructureCnt, alreadyAllocatedTo_subStructureCnt, maxWindowInfoSize*maxNrOfWindows, maxWindowInfoSize*alreadyAllocatedTo_maxNrOfWindows); /* maxWindowInfoSize = 50 ... should suffice*/

		/*allocate memory to struct binned_I_windows and initialize*/
		returnVal = alloc_init_ptr_binned_I_windows(&ptr_binned_I_windows_query, maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows, nrOfEntries);
		
		/*for pairs*/
		if(matchWindowPairs_b !=0){

			/*We'll need both the I-values on the single windows and on the pairs. The allocation for the single windows was just done, so only the pairs case needs attention:*/
			maxNrOfWindowPairs = alloc_init_I_windowPairs_ptr(&I_windowPairs, order, full_b, maxChainLen, windowLength, stepSize, maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows);
			
			/*allocate and initialize a ptr to hold the scores (on pairs) results for each query*/
			alloc_init_queryRawScore(&ptr_queryRawScorePairs, subStructureCnt, alreadyAllocatedTo_subStructureCnt, maxWindowInfoSize*maxNrOfWindowPairs, maxWindowInfoSize*alreadyAllocatedTo_maxNrOfWindowPairs); /* maxWindowInfoSize = 50 ... should suffice*/

			/*allocate memory to struct binned_I_windowPairs and initialize*/
			returnVal = alloc_init_ptr_binned_I_windowPairs(&ptr_binned_I_windowPairs_query, maxNrOfWindowPairs, alreadyAllocatedTo_maxNrOfWindowPairs, nrOfEntriesForPairs);

			/*for holding rarity scores to be written out:*/
			returnVal = alloc_init_ptr_Irarity_windowPairs(&ptr_Irarity_windowPairs, order, full_b, maxChainLen, windowLength, stepSize, maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows);
		
			/*allocate and init ptr_cntMatch -- to hold the cntMatch for a series of q-windows:*/
			if(alreadyAllocatedTo_maxNrOfWindows == 0){ //no allocation was done yet
				ptr_cntMatch = (int *) malloc(maxNrOfWindows*sizeof(int));
				/*init*/
				for(k=0; k< maxNrOfWindows; k++){
					ptr_cntMatch[k] = 0;
				}
			}
			else if(maxNrOfWindows > alreadyAllocatedTo_maxNrOfWindows){ //reallocation necc
				ptr_cntMatch = realloc(ptr_cntMatch, maxNrOfWindows*sizeof(int));
				/*init*/
				for(k=0; k< maxNrOfWindows; k++){
					ptr_cntMatch[k] = 0;
				}
			}

		}

		/*For match*/
		//printf("DBResultsFilePath:%s\n", DBResultsFilePath);

		/*convenience ptr for holding I values on windows, for query:*/
		if(alreadyAllocatedTo_maxNrOfWindows == 0){ //no allocation was done yet

			ptr_query = (double **) malloc(maxNrOfWindows*sizeof(double *)); 

			/*initialize*/
			for(k = 0; k < maxNrOfWindows; k++){
						
				ptr_query[k] = (double *) malloc(maxNrOfInvariants*sizeof(double));

				ptr_query[k][0] = 0;

				if(nrOfEntries > 1){
					ptr_query[k][1] = 0;
				}
				
				if(nrOfEntries > 2){

					for(i= 2; i < maxNrOfInvariants;i++){
						ptr_query[k][i] = 0;

					}
				}			
			
			}

		}
		else if(maxNrOfWindows > alreadyAllocatedTo_maxNrOfWindows){ //reallocation necc

			ptr_query = realloc(ptr_query, maxNrOfWindows*sizeof(double *)); 

			/*initialize*/
			for(k = 0; k < maxNrOfWindows; k++){
						
				if(k >= alreadyAllocatedTo_maxNrOfWindows){		
					ptr_query[k] = (double *) malloc(maxNrOfInvariants*sizeof(double));
				}

				ptr_query[k][0] = 0;

				if(nrOfEntries > 1){
					ptr_query[k][1] = 0;
				}
				
				if(nrOfEntries > 2){

					for(i= 2; i < maxNrOfInvariants;i++){
						ptr_query[k][i] = 0;

					}
				}			
			
			}
			
		}
		else { //just re-init

			for(k = 0; k < alreadyAllocatedTo_maxNrOfWindows; k++){

				ptr_query[k][0] = 0;

				if(nrOfEntries > 1){
					ptr_query[k][1] = 0;
				}
				
				if(nrOfEntries > 2){

					for(i= 2; i < maxNrOfInvariants;i++){
						ptr_query[k][i] = 0;

					}
				}			
			
			}
		
		
		}

		/*for pairs*/
		if(matchWindowPairs_b !=0){

			/*convenience ptr for holding MUTUAL I values on window pairs, for query.
			On the diagonal these will though amount to the single window values:*/
			if(alreadyAllocatedTo_maxNrOfWindows == 0){ //no allocation was done yet

				ptr_query_pairs = (double ***) malloc(maxNrOfWindows*sizeof(double **));

				/*alloc and initialize*/
				for(k = 0; k < maxNrOfWindows; k++){
							
					ptr_query_pairs[k] = (double **) malloc(maxNrOfWindows*sizeof(double *));

					/*initialize*/
					for(l = 0; l < maxNrOfWindows; l++){
							
						ptr_query_pairs[k][l] = (double *) malloc(maxNrOfInvariants*sizeof(double));

						ptr_query_pairs[k][l][0] = 0.0;
						
						if(nrOfEntriesForPairs > 1){
							
							ptr_query_pairs[k][l][1] = 0.0;
						}
						
						if(nrOfEntriesForPairs > 2){

							for(i= 2; i < maxNrOfInvariants;i++){
								ptr_query_pairs[k][l][i] = 0.0;

							}
						}

					}
				}
			}
			else if(maxNrOfWindows > alreadyAllocatedTo_maxNrOfWindows){ //reallocation necc

				ptr_query_pairs = realloc(ptr_query_pairs, maxNrOfWindows*sizeof(double **));

				/*alloc and initialize*/
				for(k = 0; k < maxNrOfWindows; k++){
							
					if(k < alreadyAllocatedTo_maxNrOfWindows){
						ptr_query_pairs[k] = realloc(ptr_query_pairs[k], maxNrOfWindows*sizeof(double *));
					}
					else if(k >= alreadyAllocatedTo_maxNrOfWindows){
						ptr_query_pairs[k] = (double **) malloc(maxNrOfWindows*sizeof(double *));
					}

					/*initialize*/
					for(l = 0; l < maxNrOfWindows; l++){
							
						if(l < alreadyAllocatedTo_maxNrOfWindows){
							ptr_query_pairs[k][l] = realloc(ptr_query_pairs[k][l], maxNrOfInvariants*sizeof(double));
						}
						else if(l >= alreadyAllocatedTo_maxNrOfWindows){
							ptr_query_pairs[k][l] = (double *) malloc(maxNrOfInvariants*sizeof(double));
						}

						ptr_query_pairs[k][l][0] = 0.0;
						
						if(nrOfEntriesForPairs > 1){
							
							ptr_query_pairs[k][l][1] = 0.0;
						}
						
						if(nrOfEntriesForPairs > 2){

							for(i= 2; i < maxNrOfInvariants;i++){
								ptr_query_pairs[k][l][i] = 0.0;

							}
						}

					}
				}
			}


		}


		//Update of alreadyAllocatedTo_maxNrOfWindows and alreadyAllocatedTo_maxNrOfWindowPairs is placed right before the main loop since
		//there one more place we need to handle (ie possibly reallocate based on the values so far) before we can do the update (see
		//call to aalloc_init_ptr_ptr_binned_I_windows below)
		

		/****************************************************/
		/*Load of DB for matching*/
		/****************************************************/

		if(dbLoaded_b == 0){

			/*We first load in DB results for pairs, if the "include pairs matching" is the desired search method. Subsequently we load the 
			/*DB for the single windows and bin these (the matching will go via the single  windows). In case the pairs-method is used this 
			/*load consists in populating the DB set by means of the pairs results; else the single windows DB is loaded from file.

			/**********************************************************************
			/*Search including window pairs: Allocate mem to ptr's for query, for read out and for match*/
			/**********************************************************************/
			if (matchWindowPairs_b !=0){
					

				/*Load DB for the pair I-values:*/
				intArray3 = readDBresultsWindowPairsToPtr(&ptr_dbResultWindowPairs, DBResultsPairsFileName, chunkSize, &maxNrOfWindowsDB);

				nrOfDBstructuresLoaded_pairsLoad = intArray3[0];
				nrOfDBrecordsLoaded_pairsLoad =  intArray3[1]; //equals nr of window pairs loaded

				if(halted_b == 1){
					printf("Pairs load; no str's: %d pairs: %d\n", nrOfDBstructuresLoaded_pairsLoad, nrOfDBrecordsLoaded_pairsLoad);
					getchar();
				}

				/*for(m = 0; m < 5; m++){
					printf("db file at %d is %s, I12 at 0 0: %lf\n", m, ptr_dbResultWindowPairs[m].fileName, ptr_dbResultWindowPairs[m].I12[0][0]);
				}*/

				/*NO more: Normalization of the pairs results is done below, after the single window values have been derived from the 
				pairs' values*/

				if(normalize_b == 1){

					I_windowPairs_order2_meanStddev_DB = normalizeDBresultsWindowPairsPtr(ptr_dbResultWindowPairs, nrOfDBstructuresLoaded_pairsLoad, nrOfDBrecordsLoaded_pairsLoad);

					printf("Pairs normalization: I12 move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_I12, I_windowPairs_order2_meanStddev_DB.stddev_I12);
					printf("Pairs normalization: Ia12 move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_Ia12, I_windowPairs_order2_meanStddev_DB.stddev_Ia12);
					printf("Pairs normalization: I1234_full move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_I1234_full, I_windowPairs_order2_meanStddev_DB.stddev_I1234_full);
					printf("Pairs normalization: I1324_full move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_I1324_full, I_windowPairs_order2_meanStddev_DB.stddev_I1324_full);
					printf("Pairs normalization: I1423 move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_I1423, I_windowPairs_order2_meanStddev_DB.stddev_I1423);

				}

				/*Compute the bins for the pairs data; these bins are used when binning the query's GI-pairs too. */
				printf("Generate bins based on the DB results for the pairs data) ... \n");

				/*OBS: if binningType = 0 we use a custom "hand-held" binning defined in the generateBins fct; change
				of bin definitions is needed when nrOfBins change!!*/
				returnVal = generateBinsPairs(binningType, ptr_dbResultWindowPairs, nrOfDBstructuresLoaded_pairsLoad, nrOfDBrecordsLoaded_pairsLoad, nrOfEntriesForPairs, &binsPairs, nrOfBins, compensationForRounding);

				printf("... done ... \n");

				/*Obs: we are not binning the entire pairs-db etc; we'll do the binning of the match case on the fly (each match case is a pair of single window matches from a single structure.
				We are though reserving a little convenience pointer for the match cases:*/

				/*For match*/

				/*Binned vector for holding the binned vector of each match to be considered:*/
				binnedIvector_matchPair = (int *) malloc(nrOfEntriesForPairs*sizeof(int));


				/*Convenience ptr for holding window MUTUAL values in ptr; we allocate some surplus:*/
				ptr_match_pairs = (double ***) malloc(maxNrOfWindowsDB*sizeof(double **));

				for (m=0; m<maxNrOfWindowsDB;m++){

					ptr_match_pairs[m] = (double **) malloc(maxNrOfWindowsDB*sizeof(double *));
					
					for (n=0; n<maxNrOfWindowsDB;n++){
						
						ptr_match_pairs[m][n] = (double *) malloc(maxNrOfInvariants*sizeof(double));

						ptr_match_pairs[m][n][0] = 0.0;
					
						if(nrOfEntriesForPairs > 1){
							
							ptr_match_pairs[m][n][1] = 0.0;
						}
						
						if(nrOfEntriesForPairs > 2){

							for(i= 2; i < maxNrOfInvariants;i++){
								ptr_match_pairs[m][n][i] = 0.0;
							}

						}

					}
				}

			} /*end matchWindowPairs_b != 0*/


			//printf("DBResultsFilePath:%s\n", DBResultsFilePath);

			/********************************************************************************************************************************************/
			/* Single windows DB matching/load of single window DB from pairs set:
			/********************************************************************************************************************************************/

			/*First read in DB results (single windows)*/
			printf("reading in the DB from a single windows source ...\n");
			intArray2 = readDBresultsWindowsToPtr(&ptr_dbResultWindows, DBResultsFilePath, chunkSize, &maxNrOfWindowsDB);
			printf("... done \n");

			/*  if(matchWindowPairs_b == 0){
				//First read in DB results (single windows)
				printf("reading in the DB from a single windows source ...\n");
				intArray2 = readDBresultsWindowsToPtr(&ptr_dbResultWindows, DBResultsFilePath, chunkSize, &maxNrOfWindowsDB);
				printf("... done \n");

			}
			else{ //populate the DB set, ptr_dbResultWindows, from the pairs results loaded from file above:
				if(matchWindowPairs_b != 0){ //we put an if-clause here just to be explicit about what's going on
					printf("reading in the DB from a window pairs source ...\n");
					intArray2 = loadDBWindowsFromDBPairs(&ptr_dbResultWindows, ptr_dbResultWindowPairs, nrOfDBstructuresLoaded_pairsLoad, chunkSize);
					printf("... done \n");
				}
			} */

			//getchar();

			nrOfDBstructuresLoaded = intArray2[0];
			nrOfDBrecordsLoaded=  intArray2[1]; //equals nr of windows loaded
			
			/*if(nrOfDBrecsLoaded < topN){
				
				printf("topN is set to nr of db records (%d) since this is lower than the set topN (%d)\n",  nrOfDBrecsLoaded, topN);
				topN = nrOfDBrecsLoaded;

			}*/

			printf("nr of DB records loaded: %d , nr of structures: %d,  maxNrOfWindowsDB: %d\n", nrOfDBrecordsLoaded, nrOfDBstructuresLoaded, maxNrOfWindowsDB);
			//getchar();


			/*convenience ptr for holding window values in ptr; we allocate some surplus:*/
			ptr_match = (double **) malloc(maxNrOfWindowsDB*sizeof(double *)); //malloc(nrOfEntries*maxNrOfWindowsDB*sizeof(double));

			for (n= 0; n< maxNrOfWindowsDB; n++){

				ptr_match[n] = (double *) malloc(maxNrOfInvariants*sizeof(double));

				/*init*/
				ptr_match[n][0] = 0;
				ptr_match[n][1] = 0;

				if(nrOfEntries > 2){

					for(i= 2; i < maxNrOfInvariants;i++){
						ptr_match[n][i] = 0;

					}
				}

			}


			/*normalize the DB results if desired*/
			if(normalize_b == 1){

				printf("normalizing DB results ..\n");

				//normalize the I-windows values in the DB; we get the means and std devs returned 
				//making normalization of the query's I-window values easy later: 

				printf("I12 before: %lf\n", ptr_dbResultWindows[0].I12[0]);
				if (nrOfEntries > 1){
					printf("Ia12 before: %lf\n", ptr_dbResultWindows[0].Ia12[0]);
				}
				if (nrOfEntries > 2){
					printf("I1234_full before: %lf\n", ptr_dbResultWindows[0].I1234_full[0]);
					printf("I1324_full before: %lf\n", ptr_dbResultWindows[0].I1324_full[0]);
					printf("I1423 before: %lf\n", ptr_dbResultWindows[0].I1423[0]);
				}

				I_windows_order2_meanStddev_DB = normalizeDBresultsWindowsPtr(ptr_dbResultWindows, nrOfDBstructuresLoaded);

				printf("I12 move and squeeze:%lf %lf\n", I_windows_order2_meanStddev_DB.mean_I12, I_windows_order2_meanStddev_DB.stddev_I12);
				printf("I12 after: %lf\n", ptr_dbResultWindows[0].I12[0]);
				if (nrOfEntries > 1){
					printf("Ia12 move and squeeze:%lf %lf\n", I_windows_order2_meanStddev_DB.mean_Ia12, I_windows_order2_meanStddev_DB.stddev_Ia12);
					printf("Ia12 after: %lf\n", ptr_dbResultWindows[0].Ia12[0]);
				}
				if (nrOfEntries > 2){
					printf("I1234_full move and squeeze:%lf %lf\n",  I_windows_order2_meanStddev_DB.mean_I1234_full, I_windows_order2_meanStddev_DB.stddev_I1234_full);
					printf("I1324_full move and squeeze:%lf %lf\n", I_windows_order2_meanStddev_DB.mean_I1324_full, I_windows_order2_meanStddev_DB.stddev_I1324_full);
					printf("I1423 move and squeeze:%lf %lf\n", I_windows_order2_meanStddev_DB.mean_I1423,  I_windows_order2_meanStddev_DB.stddev_I1423);
					printf("I1234_full after: %lf\n", ptr_dbResultWindows[0].I1234_full[0]);
					printf("I1324_full after: %lf\n", ptr_dbResultWindows[0].I1324_full[0]);
					printf("I1423 after: %lf\n", ptr_dbResultWindows[0].I1423[0]);
				}

				/*The normalization of the pairs results must be done AFTER the single window results have been derived from the 
				pairs' values:*/
				/* if(matchWindowPairs_b == 1){

					I_windowPairs_order2_meanStddev_DB = normalizeDBresultsWindowPairsPtr(ptr_dbResultWindowPairs, nrOfDBstructuresLoaded_pairsLoad, nrOfDBrecordsLoaded_pairsLoad);

					printf("Pairs normalization: I12 move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_I12, I_windowPairs_order2_meanStddev_DB.stddev_I12);
					printf("Pairs normalization: Ia12 move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_Ia12, I_windowPairs_order2_meanStddev_DB.stddev_Ia12);
					printf("Pairs normalization: I1234_full move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_I1234_full, I_windowPairs_order2_meanStddev_DB.stddev_I1234_full);
					printf("Pairs normalization: I1324_full move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_I1324_full, I_windowPairs_order2_meanStddev_DB.stddev_I1324_full);
					printf("Pairs normalization: I1423 move and squeeze:%lf %lf\n", I_windowPairs_order2_meanStddev_DB.mean_I1423, I_windowPairs_order2_meanStddev_DB.stddev_I1423);

				}*/

				printf(".. normalization done \n");

				if(halted_b == 1){
					getchar();
				}

			}


			/*for(m = 0; m < 5; m++){
				printf("db file at %d: %s, structure: %s, class: %s, chain: %s, I12 at 0: %lf\n", m, ptr_dbResultWindows[m].fileName,ptr_dbResultWindows[m].structureName, ptr_dbResultWindows[m].classId, ptr_dbResultWindows[m].chainId, ptr_dbResultWindows[m].I12[0]);
			}*/

			/*Now transform the DB records (I-vectors) into binned versions and sort it all lexicographically*/
			
			/*First transform to array of binned I-results:*/
			printf("binning DB results ..\n");

			/*OBS: if binningType = 0 we use a custom "hand-held" binning defined in the generateBins fct; change
			of bin definitions is needed when nrOfBins change!!*/
			returnVal = generateBins(binningType, ptr_dbResultWindows, nrOfDBstructuresLoaded, nrOfDBrecordsLoaded, nrOfEntries, &bins, nrOfBins, compensationForRounding);

			/*OBS: the ptr ptr_binned_I_windows_DB is allocated to within the binDBresultsWindows fct! (unlike for ptr_binned_I_windows_query/binQueryresultsWindows)*/
			returnVal = binDBresultsWindows(&ptr_binned_I_windows_DB, ptr_dbResultWindows, nrOfDBstructuresLoaded, nrOfDBrecordsLoaded, nrOfEntries, bins, nrOfBins);
			printf(".. done \n");

			/*Next sort lexicographically*/
			printf("sorting DB results ..");

			/*for(i=0; i<nrOfDBrecordsLoaded;i++){
				printf("rec nr %d: %d %d %d %d %d\n", i, ptr_binned_I_windows_DB[i].binnedIvector[0], ptr_binned_I_windows_DB[i].binnedIvector[1], ptr_binned_I_windows_DB[i].binnedIvector[2], ptr_binned_I_windows_DB[i].binnedIvector[3], ptr_binned_I_windows_DB[i].binnedIvector[4]);
				if(i == 25){ break;}
			}
			getchar();*/
			//bubbleSortLexicoBinnedDB_LR(ptr_binned_I_windows_DB, nrOfDBrecordsLoaded, nrOfEntries, allowedNrOfMismatches);
			heapSortBinnedDB(ptr_binned_I_windows_DB, nrOfDBrecordsLoaded, nrOfEntries);

			/*for(i=nrOfDBrecordsLoaded-1; i >= 0;i--){
				printf("fileNr: %d, rec nr %d: %d %d %d %d %d\n", ptr_binned_I_windows_DB[i].fileNr, i, ptr_binned_I_windows_DB[i].binnedIvector[0], ptr_binned_I_windows_DB[i].binnedIvector[1], ptr_binned_I_windows_DB[i].binnedIvector[2], ptr_binned_I_windows_DB[i].binnedIvector[3], ptr_binned_I_windows_DB[i].binnedIvector[4]);
				if(i == nrOfDBrecordsLoaded - 25){ break;}
			}
			*/
			printf(".. done \n");


			/****************************************************************************************************************/
			/* Allocate and initialize ptr to hold the match set returned from each "word search"; 
			/****************************************************************************************************************/

			printf("Create match set array ...\n");
			
			/*we use this when looping though the match set so as to avoid recomputing multiple times for each matching structure.
			For assistance we declare a ptr holding a 0 or 1 for each DB record indicating whether the record is a match (1) or not (0).*/

			/*array of "GI words", one word for every possible match*/
			//returnVal = alloc_init_ptr_binned_I_windows(&ptr_binned_I_windows_matchSet, nrOfDBrecordsLoaded, 0, nrOfEntries);

			/*match indicator assistant --- gets reset for every new query (below)*/
			ptr_matchIndicator = (int *)malloc(nrOfDBrecordsLoaded*sizeof(int)); 
			for(i = 0; i < nrOfDBrecordsLoaded; i++){ 
				ptr_matchIndicator[i] = 0;
			}
			
			printf(".. done \n");

			if(halted_b == 1){
				getchar();
			}

			dbLoaded_b = 1;

		}

		/*for(i=0; i<nrOfDBrecordsLoaded;i++){
			printf("rec nr %d: %d %d %d %d %d\n", i, ptr_binned_I_windows_DB[i].binnedIvector[0], ptr_binned_I_windows_DB[i].binnedIvector[1], ptr_binned_I_windows_DB[i].binnedIvector[2], ptr_binned_I_windows_DB[i].binnedIvector[3], ptr_binned_I_windows_DB[i].binnedIvector[4]);
			if(i == 25){ break;}
		}*/

		if (matchWindowPairs_b != 0){
			/*Pointer to array of "GI words"; one array (match set) for every window in the query (therefore: maxNrOfWindows and not maxNrOfWindowsDB!) 
			This serves to keep the matches found in the single windows matching for later use in the (possible) window pairs based matching (and is only
			used in that case)*/
			returnVal = alloc_init_ptr_ptr_binned_I_windows(&ptr_ptr_binned_I_windows_matchSet, nrOfDBrecordsLoaded, nrOfEntries, maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows);
		}

		/*for(i=0; i<nrOfDBrecordsLoaded;i++){
			printf("rec nr %d: %d %d %d %d %d\n", i, ptr_binned_I_windows_DB[i].binnedIvector[0], ptr_binned_I_windows_DB[i].binnedIvector[1], ptr_binned_I_windows_DB[i].binnedIvector[2], ptr_binned_I_windows_DB[i].binnedIvector[3], ptr_binned_I_windows_DB[i].binnedIvector[4]);
			if(i == 25){ break;}
		}*/


		/****************************************************************************************************************/
		/* Read in rarityScore-data and generate the distribution of these (just a sort):*/
		/****************************************************************************************************************/
		
		/*If p-values are desired (and made feasible by supplying the needed back ground distribution), read in 
		rarityScore-data and generate the distribution of these (just a sort):*/
		if(rarityScorePValues_b == 1){

			ptr_fileRarityScoreDistr = fopen(rarityScoreDistrFilePath, "r");

			if (!ptr_fileRarityScoreDistr){
				printf("Sorry, the file containing the back ground scores for the single windows, rarityScoreDistrFilePath, was not found (%s was provided). Enter a valid file name and run again.\n", rarityScoreDistrFilePath);
				return 0;
			}

			/*this will provide a sorted list (smallest to highest), ptr_distDistr, of dist-values.
			Will be used to score each window hit (as -log(p-value), where the p-value is obtained by a
			look up in this array: */
			sizeOfRarityScoreDistr = readSortRarityScoreData(ptr_fileRarityScoreDistr, &ptr_rarityScoreDistr);

			printf("sizeOfRarityScoreDistr %d\n", sizeOfRarityScoreDistr);
			if(halted_b == 1){
				getchar();
			}

			if(matchWindowPairs_b != 0){
				
				ptr_fileRarityScorePairsDistr = fopen(rarityScorePairsDistrFilePath, "r");

				if (!ptr_fileRarityScorePairsDistr){
					printf("Sorry, the file containing the back ground scores for the pairs, rarityScorePairsDistrFilePath, was not found (%s was provided). Enter a valid file name and run again.\n", rarityScorePairsDistrFilePath);
					return 0;
				}

				/*this will provide a sorted list (smallest to highest), ptr_distDistr, of dist-values.
				Will be used to score each window hit (as -log(p-value), where the p-value is obtained by a
				look up in this array: */
				sizeOfRarityScorePairsDistr = readSortRarityScoreData(ptr_fileRarityScorePairsDistr, &ptr_rarityScorePairsDistr);

				printf("sizeOfRarityScorePairsDistr %d\n", sizeOfRarityScorePairsDistr);
				if(halted_b == 1){
					getchar();
				}

			}

		}


		/*Allocation to a pointer holding a match range for each entry (Gauss integral) to be matched (ie in number nrOfEntries)*/  
		printf("nrOfEntries: %d\n", nrOfEntries);
		printf("size of matchRange: %d\n", sizeof(struct matchRange));
		printf("size of matchRange ptr: %d\n", sizeof(struct matchRange *));
		//returnVal= getchar();
		ptr_newRange = (struct matchRange **) malloc(3*(nrOfEntries+1)*sizeof(struct matchRange));
		
		if(ptr_newRange != NULL){
			for(i = 0; i < nrOfEntries; i++){

				ptr_newRange[i] = (struct matchRange *) malloc(3*sizeof(struct matchRange));
				
				if(ptr_newRange[i] != NULL){
					for(j = 0; j < 3; j++){
						ptr_newRange[i][j].start = -1;
						ptr_newRange[i][j].end = -2;
						ptr_newRange[i][j].nrOfMismatches = 999;
					}
				}
				else{
					printf("newRange[%d] not allocated for\n", i);
				}
			}
		}
		else{
			printf("post ptr_newRange alloc, nrOfEntries now: %d\n", nrOfEntries);
			if(halted_b == 1){
				getchar();
			}
		}


		/*Finally a couple of allocations to strings for read-out of scores:*/
		if(write_rarityScores_b == 1){
			qWindowsStringAid = (char *) malloc(5*sizeof(char));
			qWindowsString = (char *) malloc(20*maxNrOfWindows*sizeof(char));
			rarityScoreString = (char *) malloc(10*sizeof(char));

			/*init*/
			strcpy(qWindowsStringAid, "\0");
			strcpy(qWindowsString, "\0");
			strcpy(rarityScoreString, "\0");

			if(write_rarityScoresPairs_b == 1){
				qWindowsStringAidPairs = (char *) malloc(5*sizeof(char));
				qWindowsStringPairs = (char *) malloc(20*maxNrOfWindowPairs*sizeof(char));

				/*init*/
				strcpy(qWindowsStringAidPairs, "\0");
				strcpy(qWindowsStringPairs, "\0");

			}

		}

		//Update
		alreadyAllocatedTo_maxNrOfWindows = max(maxNrOfWindows, alreadyAllocatedTo_maxNrOfWindows);
		alreadyAllocatedTo_maxNrOfWindowPairs = max(maxNrOfWindowPairs, alreadyAllocatedTo_maxNrOfWindowPairs);
		alreadyAllocatedTo_subStructureCnt = max(alreadyAllocatedTo_subStructureCnt, subStructureCnt);
	
		//printf(" -- du er skoer, Topper\n");
		//getchar();


		/****************************************************************************************************************/
		/* MAIN LOOP */
		/****************************************************************************************************************/

		subStructureCnt = 0; /*reset*/
		//startComplete = clock();

		/*loop through the list of filenames in the queries' directory and compute the desired measures all along:*/
		for (fileNr = 0; fileNr <numberOfFiles; fileNr ++){

			if (fileNr%10 == 0){printf("Now at file nr:%d name: %s\n", fileNr, ptr_fileNameList[fileNr]);}
			
			//printf("Query: %s\n", ptr_dirList[fileNr]);
			//getchar();

			if (print_basic_b == 1){
				printf("Now handling file:%s\n fileNr:%d\n",ptr_dirList[fileNr], fileNr);
			}

			if(use_cath_b == 1){
				hitIndexDomain = bisectionalSearchCATHDomain(ptr_fileNameList[fileNr], ptr_cathDomain, nrOfCATHDomainsInList, lengthCLFformat20);
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

			//printf("p0: ..... %s\n", ptr_dbResultWindows[0].fileName);

			ptr_fileIn = fopen(ptr_dirList[fileNr], "r");

			//printf("p1: ..... %s\n", ptr_dbResultWindows[0].fileName);

			if (!ptr_fileIn){
				printf("Sorry, the file %s could not be found. Execution stops here.\n", ptr_dirList[fileNr]);
				//getchar();
				return 0;
			}
					
			I_measures.fileName = ptr_dirList[fileNr];

			/*reinit chain info*/
			returnVal = reinit_chainInStr(chainInStr, maxNrOfChains);

			/*get number of chains and their lengths for this file:*/
			//chainInStr = readPDBChainStructure(ptr_fileIn);
			if(use_scop_b ==1){
				readPDBChainStructureSCOP(ptr_fileIn, &chainInStr );
			}
			else if(use_cath_b ==1){
				readPDBDomainStructureCATH(ptr_fileIn, &chainInStr );//will read chain length and chain id into chainInStr; other values (structureName etc) where read in right above
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
				}
			}

			/****************************************************************************************************************/
			/* INNER LOOP */
			/****************************************************************************************************************/
			/*loop through the chains in the current structure/file ... and compute ... :*/
			chainNr = 0; /*reset*/
			for (chainNr = 0; chainNr <= chainInStr.numberOfChains -1 ; chainNr ++){

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

				I_measures.structureName = chainInStr.ptr_chainInStr[chainNr].structureName;
				I_measures.classId = chainInStr.ptr_chainInStr[chainNr].classId;
				I_measures.chainId = chainInStr.ptr_chainInStr[chainNr].chainId;
				I_measures.chainNr = &chainNr;


				//printf("Sommeren mild oedsles ...");
				chainLen = chainInStr.ptr_chainInStr[chainNr].chainLength;

				if (print_basic_b ==1){
					printf("sub str:%s\n",I_measures.structureName);
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
				
				if(use_scop_b ==1){
					if(strstr(chainInStr.ptr_chainInStr[chainNr].classId, "i") == chainInStr.ptr_chainInStr[chainNr].classId || 
						strstr(chainInStr.ptr_chainInStr[chainNr].classId, "j") == chainInStr.ptr_chainInStr[chainNr].classId ||
						strstr(chainInStr.ptr_chainInStr[chainNr].classId, "k") == chainInStr.ptr_chainInStr[chainNr].classId ||
						strstr(chainInStr.ptr_chainInStr[chainNr].classId, "h") == chainInStr.ptr_chainInStr[chainNr].classId){
						printf("Class %s of chain of %s was skipped since it is h, i, j or k\n", chainInStr.ptr_chainInStr[chainNr].chainId, ptr_dirList[fileNr], chainInStr.ptr_chainInStr[chainNr].classId);
						subStructureCnt += 1; //!
						continue;
					}
				}

				I_measures.chainLen = &chainLen;
				rewind(ptr_fileIn);
				returnVal = main_readPDB2(ptr_fileIn, ptr_chain, chainNr, chainLen);

				/*length of segments chain:*/
				L = chainLen - 1;
				/*allocate memory for array of segments*/
				//ptr_segment = (struct  segment *) calloc (L, sizeof(struct segment));
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
					printf("Chain: %s in structure: %s removed due to too long %d 'th segment\n", chainInStr.ptr_chainInStr[chainNr].chainId, I_measures.fileName,  segTooLong_b-1);
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
					//returnVal = aggrAndW_raw(ptr_segment, chainLen, order, full_b, &I_measures_raw);
					}
				if (incl_abs_b == 0 && split_b == 0 && closed_loops_b == 0){
					returnVal = aggrAndW_ExAbs(ptr_segment, chainLen, order, full_b, I_measures);
					}
				if (incl_abs_b == 1 && split_b == 1 && closed_loops_b == 0){
					returnVal = wAll(ptr_segment, chainLen, I_measures.wVal);
					returnVal = aggr(chainLen, order, full_b, I_measures.wVal, I_measures);
				}


				endComp = clock();

				compTime = ((double) (endComp - startComp))/CLOCKS_PER_SEC;
				//compTime = ((double) (endComp - startComp));
				if(print_b == 1){
					printf("Comp time: %f\n",compTime);
					//printf("CPU time spend for measures on this file (only computation/aggregation):%lf \n", compTime);
				}
				if (compTime > max_compTime){
					max_compTime = compTime;
				}
				compTimeComplete += compTime;

				/*Code for searching through the DB starts here*/
						
				/*********************************************************************************************/
				/*Scanning/matching based on single windows and, subsequently if desired, on window pairs*/
				/*********************************************************************************************/
				
				/*The method is: 
				1) get the query's results (per window), bin them, and look the "words" up in the likewise binned 
				database. If pairs-based matching is wanted, keep the match sets (in a pointer)
				2) Look-up each query-word (per window) in the (likewise binned) db and count the number of matches.
				3) If the additional pairs-matching is desired: the matches from the first aprt (on single windows) 
				are looped over twice (to get pairs) and matching is done (only for window pairs from one and the same 
				structure -- ie window pairs from two structures are disregarded)  
				*/
						
				if (write_chains_b != 0){
					returnVal = writeChain(fileNameQueryChains, I_measures.fileName, I_measures.structureName, I_measures.chainId, ptr_chain, *I_measures.chainLen);
				}

				returnVal = getInvariantsOnWindows(&I_windows, order, windowCoveringType, windowLength, stepSize, L, I_measures);

				/*normalize the I-values if desired*/
				if(normalize_b == 1){

					normalizeQueryWindowsPtr(I_windows, I_windows_order2_meanStddev_DB, nrOfEntries);

				}

				/*for later convenience load the query's I values into ptr's allocated for the purpose*/
				for (k = 0; k < I_windows.nrOfWindows; k++){

					ptr_query[k][0] = I_windows.I12[k];

					if (nrOfEntries > 1){
						ptr_query[k][1] = I_windows.Ia12[k];
					}
					
					if (nrOfEntries > 2){
						ptr_query[k][2] = I_windows.I1234_full[k];
						ptr_query[k][3] = I_windows.I1324_full[k];
						ptr_query[k][4] = I_windows.I1423[k];

						ptr_query[k][5] = I_windows.Ia12a34_full[k];
						ptr_query[k][6] = I_windows.Ia13a24_full[k];
						ptr_query[k][7] = I_windows.Ia14a23[k];

						ptr_query[k][8] = I_windows.Ia1234_full[k];
						ptr_query[k][9] = I_windows.Ia1324_full[k];
						ptr_query[k][10] = I_windows.Ia1423[k];

						ptr_query[k][11] = I_windows.I12a34_full[k];
						ptr_query[k][12] = I_windows.I13a24_full[k];
						ptr_query[k][13] = I_windows.I14a23[k];

					}

				}

				/*bin the query's I-window values*/
				returnVal = binQueryResultsWindows(ptr_binned_I_windows_query, I_windows, ptr_query, nrOfEntries, bins, nrOfBins);

				/*for(k=0;k<I_windows.nrOfWindows;k++){ 
					printf("win nr %d: %d %d %d %d %d\n", k, ptr_binned_I_windows_query[k].binnedIvector[0], ptr_binned_I_windows_query[k].binnedIvector[1], ptr_binned_I_windows_query[k].binnedIvector[2], ptr_binned_I_windows_query[k].binnedIvector[3], ptr_binned_I_windows_query[k].binnedIvector[4]);

				}*/

				/*reset*/
				if (matchWindowPairs_b != 0){
					for(k=0; k < I_windows.nrOfWindows; k++){
						ptr_cntMatch[k] = 0;
					}
				}

				for(k = 0; k < I_windows.nrOfWindows; k++){

					//printf("Q: win nr %d: %d %d %d %d %d\n", k, ptr_binned_I_windows_query[k].binnedIvector[0], ptr_binned_I_windows_query[k].binnedIvector[1], ptr_binned_I_windows_query[k].binnedIvector[2], ptr_binned_I_windows_query[k].binnedIvector[3], ptr_binned_I_windows_query[k].binnedIvector[4]);
					//printf("Q: win nr %d: %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", k, ptr_binned_I_windows_query[k].binnedIvector[0], ptr_binned_I_windows_query[k].binnedIvector[1], ptr_binned_I_windows_query[k].binnedIvector[2], ptr_binned_I_windows_query[k].binnedIvector[3], ptr_binned_I_windows_query[k].binnedIvector[4], ptr_binned_I_windows_query[k].binnedIvector[5], ptr_binned_I_windows_query[k].binnedIvector[6], ptr_binned_I_windows_query[k].binnedIvector[7], ptr_binned_I_windows_query[k].binnedIvector[8], ptr_binned_I_windows_query[k].binnedIvector[9], ptr_binned_I_windows_query[k].binnedIvector[10], ptr_binned_I_windows_query[k].binnedIvector[11], ptr_binned_I_windows_query[k].binnedIvector[12], ptr_binned_I_windows_query[k].binnedIvector[13]);
					//getchar();

					/*Obs: ptr_Irarity_windows is re-init'ed in the corr alloc_init fct;*/

					if(use_lex == 1){

						returnVal = bisectionalSearch2(hitLR, ptr_binned_I_windows_query[k], ptr_binned_I_windows_DB, nrOfDBrecordsLoaded, nrOfEntries, allowedNrOfMismatches);
						printf("left:%d right:%d\n", hitLR[0], hitLR[1]);
					
						/*loop over hits we just found; compute the distance between the invariants, sort
						and score the results. First we read in the matches to a ptr and sort acc to structure
						number; this allows us to avoid computing multiple times for each matching structure*/
						cntMatch = 0;
						for(n = hitLR[0]; n <= hitLR[1]; n++){
							printf("rec nr %d: %d %d %d %d %d\n", n, ptr_binned_I_windows_DB[n].binnedIvector[0], ptr_binned_I_windows_DB[n].binnedIvector[1], ptr_binned_I_windows_DB[n].binnedIvector[2], ptr_binned_I_windows_DB[n].binnedIvector[3], ptr_binned_I_windows_DB[n].binnedIvector[4]);
							//getchar();
							//ptr_binned_I_windows_matchSet[cntMatch] =  ptr_binned_I_windows_DB[n];
							cntMatch +=1;
						}

					}
					else{

						/*reset*/
						for(i=0; i <nrOfDBrecordsLoaded; i++){ ptr_matchIndicator[i] = 0;}
						
						range.start = 0;
						range.end = nrOfDBrecordsLoaded-1;
						range.nrOfMismatches = 0;

						returnVal = bisectionalSearchSingle(ptr_matchIndicator, 0, ptr_newRange, ptr_binned_I_windows_query[k], 0, range, ptr_binned_I_windows_DB, allowedNrOfMismatches, nrOfEntries);

						cntMatch = 0;
						for(i=0; i <nrOfDBrecordsLoaded; i++){

							if(ptr_matchIndicator[i] == 1){

								cntMatch +=1;

								//printf("cntMatch: %d at i: %d\n",cntMatch, i);
								//ptr_binned_I_windows_matchSet[cntMatch] = ptr_binned_I_windows_DB[i];
								//printf("rec nr %d: %d %d %d %d %d\n", i, ptr_binned_I_windows_DB[i].binnedIvector[0], ptr_binned_I_windows_DB[i].binnedIvector[1], ptr_binned_I_windows_DB[i].binnedIvector[2], ptr_binned_I_windows_DB[i].binnedIvector[3], ptr_binned_I_windows_DB[i].binnedIvector[4]);
								//getchar();
								/*if matching of pairs is to be done, collect the match sets for each window*/
								if (matchWindowPairs_b != 0){
									*(ptr_ptr_binned_I_windows_matchSet[k][cntMatch].binnedIvector) = *(ptr_binned_I_windows_DB[i].binnedIvector); //!!!!
									ptr_ptr_binned_I_windows_matchSet[k][cntMatch].fileNr = ptr_binned_I_windows_DB[i].fileNr;
									ptr_ptr_binned_I_windows_matchSet[k][cntMatch].windowNr = ptr_binned_I_windows_DB[i].windowNr;
									ptr_cntMatch[k] +=1;
								}
								
							}

						}
						
					}

					/*if(cntMatch < 100){

						for(i=0; i <nrOfDBrecordsLoaded; i++){
							if(ptr_matchIndicator[i] == 1){
								printf("rec nr %d: %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", i, ptr_binned_I_windows_DB[i].binnedIvector[0], ptr_binned_I_windows_DB[i].binnedIvector[1], ptr_binned_I_windows_DB[i].binnedIvector[2], ptr_binned_I_windows_DB[i].binnedIvector[3], ptr_binned_I_windows_DB[i].binnedIvector[4], ptr_binned_I_windows_DB[i].binnedIvector[5], ptr_binned_I_windows_DB[i].binnedIvector[6], ptr_binned_I_windows_DB[i].binnedIvector[7], ptr_binned_I_windows_DB[i].binnedIvector[8], ptr_binned_I_windows_DB[i].binnedIvector[9], ptr_binned_I_windows_DB[i].binnedIvector[10], ptr_binned_I_windows_DB[i].binnedIvector[11], ptr_binned_I_windows_DB[i].binnedIvector[12], ptr_binned_I_windows_DB[i].binnedIvector[13]);
								getchar();
							}
						}
					}*/

					//printf("Number of matches: %d\n", cntMatch);
					//getchar();

					/*We add a tiny constant to the "true" frequency so as to assign a huge score to cases having 
					no match in the data base*/
					rarityFactor = ((float) cntMatch)/nrOfDBrecordsLoaded + 1e-20; //+ 1e-20: "pseudo count" to avoid log(rarityFactor) "exploding"
					//printf("Rarityfactor (pct) %lf, cntMatch %d, nrOfDBrecordsLoaded: %d\n", 100*rarityFactor, cntMatch, nrOfDBrecordsLoaded);
					//getchar();
					/*re-init the Irarity info, copying over the info from the Iwindows struct*/ 
					collect_Irarity_windows(ptr_Irarity_windows, &I_windows, k, -logbase(rarityFactor, 10) );

				}
			
				if(write_rarityScores_b == 1){

					/*sort the scores:*/
					if(k != I_windows.nrOfWindows){
						printf("Warning: query windows seem not to be looped through properly!\n"); 
					}

					score = 0.0;
					strcpy(qWindowsStringAid, "\0");
					strcpy(qWindowsString, "\0");

					if(writeWindowInfo_b == 0){

						strcat(qWindowsString, "Getting scores per window pair was not called. Set writeWindowInfo_b (option: Y) to 1 to change.");

					}

					for(n = 0; n < I_windows.nrOfWindows; n++){

						score += ptr_Irarity_windows[n].rarityScore;

						//printf("(saaledes stadig ... at window %d score is: %lf\n", n, score);

						/*get the window info where there was a contribution to the score:*/
						if(ptr_Irarity_windows[n].rarityScore > 0){
							
							snprintf(rarityScoreString, 10, "%lf",ptr_Irarity_windows[n].rarityScore); 

							if(writeWindowInfo_b != 0){
								strcat(qWindowsString, "[(");
								snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windows[n].window.windowNr);
								strcat(qWindowsString, qWindowsStringAid);
								strcat(qWindowsString, ",");
								snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windows[n].window.segIndices[0]);
								strcat(qWindowsString, qWindowsStringAid);
								strcat(qWindowsString, ",");
								snprintf(qWindowsStringAid, sizeof(qWindowsStringAid)-1, "%d", ptr_Irarity_windows[n].window.segIndices[1]);
								strcat(qWindowsString, qWindowsStringAid);
								strcat(qWindowsString, "),");
								strcat(qWindowsString, rarityScoreString);
								strcat(qWindowsString, "]");
							}

						}

					}

					score /= (I_windows.nrOfWindows + 1); //+1 to avoid division by zero

					//printf("ganske hemmelig) ... scoren er: %lf\n", score);

					
					/*if desired (and a file containing rarityScores is provided), compute the p-value of the rarity factor:*/
					if(rarityScorePValues_b == 1){

						hitRarityScore = bisectionalSearchValueL(score, ptr_rarityScoreDistr, sizeOfRarityScoreDistr);

						/*we add a tiny constant to reveal cases having a score not in that of the range of scores
						in the data base:*/
						pValue = (1.0 -  (double)(hitRarityScore)/sizeOfRarityScoreDistr) + 1e-20; //hitDistValue is the left end of the interval in which distMin sits; we add 1e-20 to avoid numm disturbance

					}

					/*collect the results:*/
					collect_queryRawScore(ptr_queryRawScore, 0, I_windows, I_windowPairs, qWindowsString, score, pValue, cntScoredQueries, -9999);
					cntScoredQueries +=1;

				}


				/*write I-values on windows out if desired*/
				if (write_windows_b != 0){

					returnVal = writeInvariantsOnWindows(ptr_fileNameOutWindows, I_windows, L, order);
				}

				/**************************************************************************************/
				/*For scan based on pairs of windows*/
				/**************************************************************************************/

				if (matchWindowPairs_b != 0){


					returnVal = getInvariantsOnWindowPairs(&I_windowPairs, order, windowCoveringType, windowLength, stepSize, L, I_measures);

					if(I_windows.nrOfWindows != I_windowPairs.nrOfWindows){
							printf("For file: %s the nr of windows recorded by I_windows (%d) and by I_windowPairs (%d) are not equal! This is a no-go ...\n", I_measures.fileName, I_windows.nrOfWindows, I_windowPairs.nrOfWindows);
							return returnVal;
					}

					/*normalize the pairs' I-values if desired. The single window values (I_windows) and the window
					pair values (I_windowPairs) are normalized differently (for single windows see above)*/
					if(normalize_b == 1){

						//printf("Before normalizing I_windowPairs: I12: %lf; Ia12: %lf; I1234_full: %lf; I1324_full: %lf;; I1423: %lf\n",  I_windowPairs.I12[10][11], I_windowPairs.Ia12[10][11], I_windowPairs.I1234_full[10][11], I_windowPairs.I1324_full[10][11], I_windowPairs.I1423[10][11]);
						normalizeQueryWindowPairsPtr(I_windowPairs, I_windowPairs_order2_meanStddev_DB, nrOfEntriesForPairs);
						//printf("After normalizing I_windowPairs: I12: %lf; Ia12: %lf; I1234_full: %lf; I1324_full: %lf;; I1423: %lf\n",  I_windowPairs.I12[10][11], I_windowPairs.Ia12[10][11], I_windowPairs.I1234_full[10][11], I_windowPairs.I1324_full[10][11], I_windowPairs.I1423[10][11]);

					}


					/*for later convenience load the query's I values for pairs (mutuals) into ptr's allocated for the purpose*/
					for (k = 0; k < I_windows.nrOfWindows; k++){

						for (l = k+1; l < I_windows.nrOfWindows; l++){

							ptr_query_pairs[k][l][0] = I_windowPairs.I12[k][l];
							//printf("q I12 at %d %d :%lf\n", k, l, I_windowPairs.I12[k][l]);
							if (nrOfEntriesForPairs > 1){
								ptr_query_pairs[k][l][1] = I_windowPairs.Ia12[k][l];
							}
							if (nrOfEntriesForPairs > 2){

								ptr_query_pairs[k][l][2] = I_windowPairs.I1234_full[k][l];
								ptr_query_pairs[k][l][3] = I_windowPairs.I1324_full[k][l];
								ptr_query_pairs[k][l][4] = I_windowPairs.I1423[k][l];
								
								ptr_query_pairs[k][l][5] = I_windowPairs.Ia12a34_full[k][l];
								ptr_query_pairs[k][l][6] = I_windowPairs.Ia13a24_full[k][l];
								ptr_query_pairs[k][l][7] = I_windowPairs.Ia14a23[k][l];

								ptr_query_pairs[k][l][8] = I_windowPairs.Ia1234_full[k][l];
								ptr_query_pairs[k][l][9] = I_windowPairs.Ia1324_full[k][l];
								ptr_query_pairs[k][l][10] = I_windowPairs.Ia1423[k][l];

								ptr_query_pairs[k][l][11] = I_windowPairs.I12a34_full[k][l];
								ptr_query_pairs[k][l][12] = I_windowPairs.I13a24_full[k][l];
								ptr_query_pairs[k][l][13] = I_windowPairs.I14a23[k][l];
								
							}
					
						}
				
					}

									
					/*bin the query's I-window Pair values*/
					returnVal = binQueryResultsWindowPairs(ptr_binned_I_windowPairs_query, I_windowPairs, ptr_query_pairs, nrOfEntriesForPairs, binsPairs, nrOfBins, onlyDisjointPairs_b);
		

					cntPairs = -1;
					cntPairsConsidered = 0;
					/*Now the matching: double loop over the query's windows:*/
					for (k = 0; k < I_windows.nrOfWindows ; k++){

						for (l = k+1; l < I_windows.nrOfWindows; l++){

							/*Skip the pair if it's not disjoint:*/
							if(I_windows.window[k].segIndices[1] >= I_windows.window[l].segIndices[0]){
								continue;
							}
							/*if(I_windows.segIndices[k][1] >= I_windows.segIndices[l][0]){
								continue;
							}*/

							//cntPairs += 1; /*counts the disjoint pairs*/

							cntPairs += 1; /*will index the window pairs as when binning the pairs' values (in fct binQueryResultsWindowPairs)*/
							/*resets*/
							rarityFactorPairs = 1.0;
							ptr_Irarity_windowPairs[cntPairs].rarityScore = 0.0;

							//for(u=0; u< nrOfEntriesForPairs; u++){printf("ptr_binned_I_windowPairs_query at %d %lf \n", cntPairs, ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[u]);}


							/*skip the pair if the threshold on mutual value is not exceded:*/
							/*only if the mutual value ptr_query_pairs[k+u][l+u] is large enough do we care to find a match.
							But only when l != k, ie not for the pairs which actually are only a single window -- since we want to
							match also windows that have low I-values :*/
							maxMutVal = 0;
							if(l>k){
								for (u = 0; u < nrOfEntriesForPairs; u++){ 
										mutVal = ptr_query_pairs[k][l][u];
									if(fabs(mutVal) > maxMutVal){
										maxMutVal = fabs(mutVal);
									}
								}
							}
							/*skip this window (nr: l) if mutual value too low:*/ 
							if (l>k && maxMutVal < thresholdMutual){continue;}

							cntMutAboveThreshold += 1;

							/*Look up the candidates we found in the single window scan
							and loop over both these. We though only need to carry out the 
							matching if both matching windows belong to the same candidate:*/

							/*reset*/
							cntMatch_12 = 0;

							for(m=0; m < ptr_cntMatch[k]; m++){

								for(n=0; n < ptr_cntMatch[l]; n++){
						
									/*if the two matching windows do not belong to the same structure (the same candidate)
									we (must) skip it:*/
									if(ptr_ptr_binned_I_windows_matchSet[k][m].fileNr != ptr_ptr_binned_I_windows_matchSet[l][n].fileNr || ptr_ptr_binned_I_windows_matchSet[k][m].fileNr < 0){ 
										continue;
									}

									/*as we have the file nr of the match we can get the window nr's of the
									matching windows, and then just look up the invariant values needed for
									the scoring:*/
									db_fileNr = ptr_ptr_binned_I_windows_matchSet[k][m].fileNr;
									db_wNr1 = ptr_ptr_binned_I_windows_matchSet[k][m].windowNr;
									db_wNr2 =  ptr_ptr_binned_I_windows_matchSet[l][n].windowNr;


									/*for convenience load the I values into a ptr allocated for the purpose*/
									i = db_wNr1;
									j = db_wNr2;

									ptr_match_pairs[i][j][0] = ptr_dbResultWindowPairs[db_fileNr].I12[i][j];
									//printf("ptr_match_pairs at %d %d : %lf\n", i, j,ptr_match_pairs[i][j][0]);
									
									if (nrOfEntriesForPairs > 1){
										ptr_match_pairs[i][j][1] = ptr_dbResultWindowPairs[db_fileNr].Ia12[i][j];
									}
									if (nrOfEntriesForPairs > 2){

										ptr_match_pairs[i][j][2] = ptr_dbResultWindowPairs[db_fileNr].I1234_full[i][j];
										ptr_match_pairs[i][j][3] = ptr_dbResultWindowPairs[db_fileNr].I1324_full[i][j];
										ptr_match_pairs[i][j][4] = ptr_dbResultWindowPairs[db_fileNr].I1423[i][j];
									
										ptr_match_pairs[i][j][5] = ptr_dbResultWindowPairs[db_fileNr].Ia12a34_full[i][j];
										ptr_match_pairs[i][j][6] = ptr_dbResultWindowPairs[db_fileNr].Ia13a24_full[i][j];
										ptr_match_pairs[i][j][7] = ptr_dbResultWindowPairs[db_fileNr].Ia14a23[i][j];

										ptr_match_pairs[i][j][8] = ptr_dbResultWindowPairs[db_fileNr].Ia1234_full[i][j];
										ptr_match_pairs[i][j][9] = ptr_dbResultWindowPairs[db_fileNr].Ia1324_full[i][j];
										ptr_match_pairs[i][j][10] = ptr_dbResultWindowPairs[db_fileNr].Ia1423[i][j];

										ptr_match_pairs[i][j][11] = ptr_dbResultWindowPairs[db_fileNr].I12a34_full[i][j];
										ptr_match_pairs[i][j][12] = ptr_dbResultWindowPairs[db_fileNr].I13a24_full[i][j];
										ptr_match_pairs[i][j][13] = ptr_dbResultWindowPairs[db_fileNr].I14a23[i][j];
									
									}

									/*bin the match candidate's I-window pairs values*/
									//printf("ptr_match_pairs[i][j] %lf\n", ptr_match_pairs[i][j]);
									returnVal = binArray(binnedIvector_matchPair,  ptr_match_pairs[i][j], nrOfEntriesForPairs, binsPairs, nrOfBins);
									//for(u=0; u< nrOfEntriesForPairs; u++){printf("binnedIvector_matchPair for entry %d %d\n",u, binnedIvector_matchPair[u]); }
 
									nrOfMismatchesThisPair = 0;
									qualified_b = 1;
									for(u=0; u< nrOfEntriesForPairs; u++){

										//printf("Q: %s w1 %d w2 %d; db: %s w1 %d w2 %d\n", I_windowPairs.structureName, k, l, ptr_dbResultWindowPairs[db_fileNr].structureName,i, j);
										//printf("Sommeren mild oedsles k..... u:%d cntPairs: %d %d %d\n", u, cntPairs, binnedIvector_matchPair[u], ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[u]);
										diff = abs(binnedIvector_matchPair[u] - ptr_binned_I_windowPairs_query[cntPairs].binnedIvector[u]);
										if(diff < 2){

											nrOfMismatchesThisPair += diff;

										}
										else{
											
											qualified_b = 0;
										}
									}
									
									if(qualified_b == 1 && nrOfMismatchesThisPair <= allowedNrOfMismatchesPairs){

										cntMatch_12 += 1;

									}

									

								} //end n-loop (candidates to 2nd window)
				

							} //end m-loop (candidates to 1st window)
							
							//printf("cntMatch_12 %d\n", cntMatch_12 );
							rarityFactorPairs = ((float) cntMatch_12 )/nrOfDBrecordsLoaded_pairsLoad + 1e-20; //+ 1: pseudo count to avoid log(rarityFactorPairs) "exploding"
							/* if(rarityFactorPairs >0){
								printf("Rarityfactor pairs (pct) %lf, cntMatch_12 %d, nrOfDBrecordsLoaded pairs: %d\n",  100*rarityFactorPairs, cntMatch_12, nrOfDBrecordsLoaded_pairsLoad);
								//getchar();
							}*/

							/*re-init the Irarity info, copying over the info from the Iwindows struct*/ 
							collect_Irarity_windowPairs(ptr_Irarity_windowPairs, &I_windowPairs, cntPairs, k,l, -logbase(rarityFactorPairs, 10) ); //cntPairsConsidered rather than cntPairs?? No, cntPairs indexes the window pairs here

							/*increment here since the present (k,l)-pair has now been considered 
							(pairs not reaching this point have not)*/
							cntPairsConsidered += 1;
							
						} //end query l-window loop

					} //end query k-window loop

					//printf("Langmodigt venter bysvalen stadig stædigt længes");
					//getchar();


					if(write_rarityScoresPairs_b == 1){

						if(k != I_windows.nrOfWindows || l != I_windows.nrOfWindows){
							printf("Warning: query window pairs seem not to be looped through properly!\n"); 
						}

						/*Alternative???: do the heap sort and only sum up the first cntPairsConsidered */
						//heapSortRarityScoresPairs(ptr_Irarity_windowPairs, cntPairsConsidered);

						scorePairs = 0.0;
						strcpy(qWindowsStringAidPairs, "\0");
						strcpy(qWindowsStringPairs, "\0");

						if(writeWindowInfo_b == 0){

							strcat(qWindowsStringPairs, "Getting scores per window pair was not called. Set writeWindowInfo_b (option: Y) to 1 to change.");

						}

						/*looping up to cntPairsConsidered appears ok here, but we recorded ptr_Irarity_windowPairs in collect_Irarity_windowpairs using cntPairs*/
						for(n = 0; n < cntPairs  ; n++){ 

							scorePairs += ptr_Irarity_windowPairs[n].rarityScore;

							snprintf(rarityScoreString, 10, "%lf",ptr_Irarity_windowPairs[n].rarityScore); 

							/*get the window info where there was a contribution to the score:*/
							if(ptr_Irarity_windowPairs[n].rarityScore > 0){
							
								if(writeWindowInfo_b != 0){
									strcat(qWindowsStringPairs, "[(");
									snprintf(qWindowsStringAidPairs, sizeof(qWindowsStringAidPairs)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.windowNr_1);
									strcat(qWindowsStringPairs, qWindowsStringAidPairs);
									strcat(qWindowsStringPairs, ",");
									snprintf(qWindowsStringAidPairs, sizeof(qWindowsStringAidPairs)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.segIndices_1[0]);
									strcat(qWindowsStringPairs, qWindowsStringAidPairs);
									strcat(qWindowsStringPairs, ",");
									snprintf(qWindowsStringAidPairs, sizeof(qWindowsStringAidPairs)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.segIndices_1[1]);
									strcat(qWindowsStringPairs, qWindowsStringAidPairs);
									strcat(qWindowsStringPairs, "),(");

									snprintf(qWindowsStringAidPairs, sizeof(qWindowsStringAidPairs)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.windowNr_2);
									strcat(qWindowsStringPairs, qWindowsStringAidPairs);
									strcat(qWindowsStringPairs, ",");
									snprintf(qWindowsStringAidPairs, sizeof(qWindowsStringAidPairs)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.segIndices_2[0]);
									strcat(qWindowsStringPairs, qWindowsStringAidPairs);
									strcat(qWindowsStringPairs, ",");
									snprintf(qWindowsStringAidPairs, sizeof(qWindowsStringAidPairs)-1, "%d", ptr_Irarity_windowPairs[n].windowPair.segIndices_2[1]);
									strcat(qWindowsStringPairs, qWindowsStringAidPairs);
									strcat(qWindowsStringPairs, "),");
									strcat(qWindowsStringPairs, rarityScoreString);
									strcat(qWindowsStringPairs, "]");
								}

							}

						}


						//For averaging it could seem reasonable to use cntPairsConsidered rather than cntPairs: the back ground data base should 
						//be build with the same constraints: only disjoint pairs or not? and the same threshold on mutuals. However, using the 
						//actual number of windows (of the the same disjoint-ness as in the data base), we get scores that are comparable across
						//diff values of the threshold:
						if(cntPairs == -1){
							scorePairs = 0.0;
						}
						else{
							scorePairs /= (cntPairs+1); //+1 to avoid dividing by 0 
						}


						if(rarityScorePValues_b == 1){

							hitRarityScorePairs = bisectionalSearchValueL(scorePairs, ptr_rarityScorePairsDistr, sizeOfRarityScorePairsDistr);

							/*we add a tiny constant to reveal cases having a score not in that of the range of scores
							in the data base:*/
							pValuePairs = (1.0 -  (double)(hitRarityScorePairs)/sizeOfRarityScorePairsDistr) + 1e-20; //hitDistValue is the left end of the interval in which distMin sits; we add 1 to avoid numm disturbance

						}

						/*collect the results:*/
						collect_queryRawScore(ptr_queryRawScorePairs, 1, I_windows, I_windowPairs, qWindowsStringPairs, scorePairs, pValuePairs, cntScoredQueriesPairs, -9999);
						cntScoredQueriesPairs +=1;

					}

					if (write_windowPairs_b != 0){

						returnVal = writeInvariantsOnwindowPairs(ptr_fileNameOutwindowPairs, I_windowPairs, L, order, 0);
					}

				} //end if-matchWindowPairs_b

				//Reset! The reason is that the single-window-loop and/or pairs-loop may turn out to be empty -- no window*/pairs fulfull the criteria in there; 
				//if omitting these resets the data from the most recent chain which had a window/pair fulfilling the criteria will be used in the ptr_queryRawScore!!
				//*: for the single windows there is though currently no criteria within that loop, but the reset will not harm!
				init_ptr_Irarity_windows(&ptr_Irarity_windows, alreadyAllocatedTo_maxNrOfWindows);
				if(matchWindowPairs_b !=0){ init_ptr_Irarity_windowPairs(&ptr_Irarity_windowPairs, alreadyAllocatedTo_maxNrOfWindows);}



				endChain = clock();
				timeChain = ((double) (endChain - startChain))/CLOCKS_PER_SEC;
				//printf("CPU time spend for measures on this file incl load of chain and more:%lf \n", timeChain);
				if (timeChain > max_timeChain){
					max_timeChain = timeChain;
				}


				subStructureCnt += 1; /*will count and index the sub-structures for which we compute the invariants*/
			
			} /*end chainNr loop*/

			
			fclose(ptr_fileIn);

		} /*end fileNr loop*/


		/*sort the obtained scoring results:*/
		heapSortRawScores(ptr_queryRawScore, cntScoredQueries);
		heapSortRawScores(ptr_queryRawScorePairs, cntScoredQueriesPairs);

		if(write_rarityScores_b == 1){

			for(i=0;i<cntScoredQueries;i++){

				fprintf(ptr_fileOutScores, "%s;%s;%s;%s;%s;%d;%d;%s;%lf;%lf;\n",
							DBName,
							ptr_queryRawScore[i].fileName,
							ptr_queryRawScore[i].structureName,
							ptr_queryRawScore[i].chainId,
							ptr_queryRawScore[i].classId,
							ptr_queryRawScore[i].chainLen,
							ptr_queryRawScore[i].nrOfWindows, //query nr of windows
							ptr_queryRawScore[i].windowInfo,
							ptr_queryRawScore[i].score,
							ptr_queryRawScore[i].pValue
							);
			}

			fclose(ptr_fileOutScores);
		
		}

		if(write_rarityScoresPairs_b == 1){

			for(i=0;i<cntScoredQueriesPairs;i++){

				fprintf(ptr_fileOutScoresPairs, "%s;%s;%s;%s;%s;%d;%d;%s;%lf;%lf;\n",
							DBName,
							ptr_queryRawScorePairs[i].fileName,
							ptr_queryRawScorePairs[i].structureName,
							ptr_queryRawScorePairs[i].chainId,
							ptr_queryRawScorePairs[i].classId,
							ptr_queryRawScorePairs[i].chainLen,
							ptr_queryRawScorePairs[i].nrOfWindows, //query nr of windows
							ptr_queryRawScorePairs[i].windowInfo,
							ptr_queryRawScorePairs[i].score,
							ptr_queryRawScorePairs[i].pValue
							);
			}

			fclose(ptr_fileOutScoresPairs);
		
		}


		endComplete = clock();
		//printf("start clocks:%d, end clocks: %d, clocks_per sec:%d", endComplete, startComplete, CLOCKS_PER_SEC);
		timeComplete = ((double) (endComplete - startComplete))/CLOCKS_PER_SEC;
		//compTimeComplete = ((double) (compTimeComplete))/CLOCKS_PER_SEC;
		printf("CPU time spend for all files incl loading data:%lf \n", timeComplete);
		printf("CPU time spend for computations on all files:%lf \n", compTimeComplete);
		printf("Max CPU time spend for measures across all files:%lf \n", max_timeChain);

		printf("Done number of files:%d \n", fileNr);
		printf("Done number of chains/sub-structures:%d \n", subStructureCnt);
		printf("%d too short and %d too long chains were skipped\n", chainsSkippedTooShort,chainsSkippedTooLong);
		printf("%d chains were skipped due to a too long segment (longer than sqrt of %d)\n",chainSkippedSegmTooLong, stdRealSegLength);

		printf("cntMutAboveThreshold: %d\n", cntMutAboveThreshold);

		/*If user wants to score another queries set without reloading the data base: 
		Prompt user a new queries set, else set keepGoing_b = 0 (will then leave the outer while-loop):*/
		keepGoing_b = 0;
		printf("Do you want to run another scan? (yes: 1, no: 0) ");
		fgets(keepGoingStr, 3, stdin);
		//strcpy(keepGoingStr, "1");
		printf( "You entered: %s\n", keepGoingStr);
		keepGoing_b = atoi(keepGoingStr);
		//printf("keepGoing_b %d\n", keepGoing_b);
		//keepGoing_b = 1;
		if(keepGoing_b == 1){

			killroyWasHere_b = 1;

			/*Resets:*/
			strcpy(queriesDirPath, "\0");
			strcpy(queriesName, "\0");
			strcpy(outputPath, "\0");

			printf("Enter a full path to another queries set:");
			fgets(queriesDirPath, 1000, stdin);
			queriesDirPath[strcspn(queriesDirPath, "\n")] = 0; //removes trailing \n if any!
			printf("Queries path: %s", queriesDirPath );
			printf("Enter a name for this queries set:");
			fgets(queriesName, 100, stdin);
			queriesName[strcspn(queriesName, "\n")] = 0; //removes trailing \n if any!

			printf("Enter a full path to where the output should be placed:");
			fgets(outputPath, 1000, stdin);
			//strcat(outputPath,"\0" );
			outputPath[strcspn(outputPath, "\n")] = 0; //removes trainling \n if any!
			printf("Output path: %s", outputPath );

			printf("Is the new queries set from CATH, SCOP or none of these: (0,1 or 2)?");
			fgets(dataTypeIndicatorStr, 3, stdin);
			dataTypeIndicator = atoi(dataTypeIndicatorStr);
			//dataTypeIndicator = 2; 
			use_cath_b = 0;
			use_scop_b = 0;
			if(dataTypeIndicator == 0){
				use_cath_b = 1;
			}
			else if(dataTypeIndicator == 1){
				use_scop_b = 1;
			}


			/*Reset file names:*/
			strcpy(fileNameOutScores1, "");
			strcpy(fileNameOutScores1, "/RarityScan2_Scores_windowslgth_");
		    strcpy(fileNameOutScoresPairs1, "");
			strcpy(fileNameOutScoresPairs1, "/RarityScan2_ScoresPairs_windowslgth_");
			strcpy(fileNameOutwindowPairs1,"/RarityScan2_Invariants_Pairs_windowlgth_");	
			strcpy(fileNameOut2, "_order_");

			strcpy(fileNameOut3,"\0");

			strcpy(fileNameOutWindows, "\0");
			strcpy(fileNameOutScores, "\0");
			strcpy(fileNameOutwindowPairs, "\0");
			strcpy(fileNameOutScoresPairs, "\0");

			/*Other resets*/
			reinit_queryRawScore(ptr_queryRawScore, cntScoredQueries);
			reinit_queryRawScore(ptr_queryRawScorePairs, cntScoredQueriesPairs);

			subStructureCnt = 0;
			cntMutAboveThreshold = 0;
			cntScoredQueries = 0; 
			cntScoredQueriesPairs = 0;
			fileNr = 0;
			chainNr = 0;
			maxChainLen = 0;

			chainsSkippedTooShort = 0;
			chainsSkippedTooLong = 0;
			chainSkippedSegmTooLong = 0;


			/* 
			cntRuns += 1;

			if(cntRuns == 1){

				strcpy(queriesDirPath, "/isdata/kroghgrp/tkj375/data/structural/Kinemage/3strs");
				strcpy(queriesName, "3strs");
				strcpy(outputPath, "/isdata/kroghgrp/tkj375/GISA_unix/results/v9/top8000");

				//strcpy(queriesDirPath, "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH_test1");
				//strcpy(queriesName, "bysvalen");
				//strcpy(outputPath, "/isdata/kroghgrp/tkj375/GISA_unix/results/v9/CATH_test1");

			}
			else if(cntRuns == 2){

				strcpy(queriesDirPath, "/isdata/kroghgrp/tkj375/data/structural/Kinemage/3links");
				strcpy(queriesName, "3links");
				strcpy(outputPath, "/isdata/kroghgrp/tkj375/GISA_unix/results/v9/top100");


			}
			*/


			///isdata/kroghgrp/tkj375/data/structural/CATH/CATH_test1
			///isdata/kroghgrp/tkj375/GISA_unix/results/v9/CATH_test1

			///isdata/kroghgrp/tkj375/data/structural/Kinemage/3links/
			///isdata/kroghgrp/tkj375/data/structural/Kinemage/top100H/
			///isdata/kroghgrp/tkj375/GISA_unix/results/v9/top100


		}

	}

	/*free the memory allocated to this query set*/
	free_I_measures(I_measures, order, full_b, alreadyAllocatedTo_chainLen);
	free(ptr_segment);	
	free(ptr_chain);
	free(dirContent.ptr_dirList); //this is allocated in the ListDirectoryContents2 fct 
	free(dirContent.ptr_fileNameList); //this is allocated in the ListDirectoryContents2 fct 
	//Resets:
	dirContent.numberOfFiles = 0;
	subDirCnt = 0;
	//More freeing:
	//for single window scores
	for(i=0;i< subStructureCnt; i++){
		free(ptr_queryRawScore[i].fileName);
		free(ptr_queryRawScore[i].structureName);
		free(ptr_queryRawScore[i].classId);
		free(ptr_queryRawScore[i].chainId);
		free(ptr_queryRawScore[i].windowInfo);
	}
	free(ptr_queryRawScore);
	//for window pair scores
	for(i=0;i< subStructureCnt; i++){
		free(ptr_queryRawScorePairs[i].fileName);
		free(ptr_queryRawScorePairs[i].structureName);
		free(ptr_queryRawScorePairs[i].classId);
		free(ptr_queryRawScorePairs[i].chainId);
		free(ptr_queryRawScorePairs[i].windowInfo);
	}
	free(ptr_queryRawScorePairs);
	
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