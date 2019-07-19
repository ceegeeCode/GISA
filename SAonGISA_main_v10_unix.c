/* **********************************************************************************************************************************
*
*
* ***************************************   Wrapper of: Strucutral Analysis based on GISA  (SAonGISA) *******************************************************************************
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
#include "SAonGISA_v10_unix.h"
#include <error.h>
#include <argp.h> 


/*Arguments and options -- ie argp handling -- particular for a data base run*/ 
typedef struct arguments_db{
		
		/*full path to directory holding the pdb files (potentially in sub-directories)*/ 
		char dirPath[1000];  

		/*Optional name for data base*/
		char DBName[100];

		/*compute invariants up to and including this order:*/
		int order;
		/*include absolute value versions of order one and two invariants:*/
		int incl_abs_b;
		/* compute full order 2 measures across the simplex*/
		int full_b;

		int get_windows_b;
		int get_windowPairs_b; /*at most one of get_windows_b and get_windowPairs_b should be set to 1*/
		

} arguments_db;

static struct argp_option options_db[] = {

		{ "dirPath", 'd',"CHAR",0, "Full path to directory holding the pdb files (potentially in sub-directories). Mandatory in flavour makeDB."},
		{ "DBName", 'D', "CHAR",OPTION_ARG_OPTIONAL, "Name of the data base (your choice, default is '')."},
		{ "order", 'o', "NUM", OPTION_ARG_OPTIONAL, "Compute the Gauss integrals for orders up to and including this integer. Must be 0,1,2 or 3 (default: 1)."},
		{ "inclAbs_b", 'a', "NUM", OPTION_ARG_OPTIONAL, "Set to 1 to include computation of absolute value-version Gauss integrals (default: 0). Not relevant in flavour rar0."},
		{ "full_b", 'f', "NUM", OPTION_ARG_OPTIONAL, "Boolean. Set to 1 to include computation of Gauss integrals for all sub-chains (default: 0). Not relevant in flavour rar0."},
		{ "get_windows_b", 'g',  "NUM", OPTION_ARG_OPTIONAL, "Set this boolean to 1 to get the invariant values for all windows (default: 0). Do not set get_windowPairs_b to 1 too."},
		{ "get_windowPairs_b", 'G',  "NUM", OPTION_ARG_OPTIONAL, "Set this boolean to 1 to get the invariant values for all windows pairs (default: 0). Do not set get_windows_b to 1 too."},
		{ 0 }
    
};

static error_t parse_opt_db (int key, char *arg, struct argp_state *state){

	struct arguments_db *arguments = state->input;


	switch (key){

		case 'd':{

			if (arg == NULL) {
				printf ("Path to directory holding PDB files is mandatory, but not provided. Exection stops ... change to get output");
				//getchar();
				return;
				}
			else{
				sprintf(arguments->dirPath, "%s", arg);
				printf ("Parameter dirPath is set to: %s\n", arguments->dirPath);
				//getchar();
			} 

			break;	
		}

		case 'D':{

			if (arg == NULL) {
				printf ("Name of data base is left at default (blank).");
				}
			else{
				sprintf(arguments->DBName, "%s", arg);
				printf ("Parameter DBName is set to: %s\n", arguments->DBName);
			} 

			break;	
		}

		case 'o':{

			if (arg == NULL) {
				printf ("Parameter order is left at default: %d\n", arguments->order);
				//getchar();
				}
			else{
				arguments->order = atoi(arg); 
				printf ("Parameter order is set to: %d\n", arguments->order);
				//getchar();
			} 	

			break;
		
		}

		case 'f':{

			if (arg == NULL) {
				printf ("Boolean full_b is left at default: %d\n", arguments->full_b);
				}
			else{
				arguments->full_b = atoi(arg); 
				printf ("Boolean full_b is set to: %d\n", arguments->full_b);
			} 

			break;
		}

		case 'a':{

			if (arg == NULL) {
				printf ("Boolean incl_abs_b controlling whether to incl abs-value versions of invariants is left at default: %d\n", arguments->incl_abs_b);
				}
			else{
				arguments->incl_abs_b = atoi(arg); 
				printf ("Boolean incl_abs_b whether to incl abs-value versions of invariants is set to: %d\n", arguments->incl_abs_b);
			} 

			break;
		}


		case 'g':{

			if (arg == NULL) {
				printf ("Invariant values on windows not fetched since get_windows_b set to default: %d\n", arguments->get_windows_b);
				}
			else{
				arguments->get_windows_b = atoi(arg); 
				printf ("Invariant values on windows are fetched if this parameter (get_windows_b) is set 1. Is here set to: %d\n", arguments->get_windows_b);
			} 

			break;
		}

		case 'G':{

			if (arg == NULL) {
				printf ("Invariant values on windows pairs not fetched since get_windowPairs_b set to default: %d\n", arguments->get_windowPairs_b);
				}
			else{
				arguments->get_windowPairs_b = atoi(arg); 
				printf ("Invariant values on window pairs are fetched if this parameter (get_windowPairs_b) is set 1. Is here set to: %d\n", arguments->get_windowPairs_b);
			} 

			break;
		}
		
	}

	return 0;
}

static struct argp argp_db = { options_db, parse_opt_db };





/*Arguments and options -- ie argp handling -- particular for a rawRarity0 run (rar0)*/ 
typedef struct arguments_rar0{
		
		/*full path to directory holding the pdb files (potentially in sub-directories)*/ 
		char queriesDirPath[1000];
		/*Name of queries set*/  
		char queriesName[100];

		char DBResultsPairsFilePath[1000];

		int useOnlyMaxMutVal_b;
		int signed_b;

		int rarityScorePValues_b; /*flavour rarity 0: if set to 1 and useOnlyMaxMutVal_b = 0, provide a absMutValScoresDistr*/
		
		int getAvgScore_b; 


} arguments_rar0;

static struct argp_option options_rar0[] = {

		{ "queriesDirPath", 'q',"CHAR",0, "Full path to directory holding the pdb files of the queries (potentially in sub-directories). Mandatory."},
		{ "queriesName", 'Q',"CHAR",OPTION_ARG_OPTIONAL, "Name for queries set (recommended)."},
		{ "pairInvsDirPath", 'I',"CHAR",0, "Full path to file holding the Gauss invariants for the window pairs. Mandatory."},
		{ "useOnlyMaxMutVal_b", 'u', "NUM", OPTION_ARG_OPTIONAL, "Boolean (default: 0). If set to 1 the scan is based on the max mutual writhe values per structure."},
		{ "signed_b", 'U', "NUM", OPTION_ARG_OPTIONAL, "Boolean (default: 1). If the scan is based on the max mututal writhe values (-u1), setting this to 1 positive and negative writhe values are separated; else the scan will be based only on abs value of writhe."},
		{ "rarityScorePValues_b", 'P', "NUM", OPTION_ARG_OPTIONAL, "Boolean. If set to 1 (and, in flavour rar0, if useOnlyMaxMutVal_b = 0) p-values are computed based on the back ground scores distribution, which must therefore be provided (default: 0))."},
		{ "getAvgScore_b", 'A', "NUM", OPTION_ARG_OPTIONAL, "Boolean. If set to 1, the score is obtained as an average over all pairs in the structure (default: 1)."},
		{ 0 }
    
};

static error_t parse_opt_rar0 (int key, char *arg, struct argp_state *state){

	struct arguments_rar0 *arguments = state->input;


	switch (key){

		case 'q':{

			if (arg == NULL) {
				printf ("Path to directory holding the queries PDB files is mandatory, but not provided. Execution stops ... change to get output. \n");
				return;
				}
			else{
				sprintf(arguments->queriesDirPath, "%s", arg);
				printf ("Parameter queriesDirPath is set to: %s\n", arguments->queriesDirPath);
			} 

			break;	
		}

		case 'Q':{

			if (arg == NULL) {
				printf ("Name of queries set is left at default (blank).\n");
				}
			else{
				sprintf(arguments->queriesName, "%s", arg);
				printf ("Name of queries set is set to: %s\n", arguments->queriesName);
			} 

			break;	
		}


		case 'I':{

			if (arg == NULL) {
				printf ("Path to file holding the Gauss invariants for the window pairs is mandatory, but not provided. Execution stops ... change to get output.\n");
				return;
				}
			else{
				sprintf(arguments->DBResultsPairsFilePath, "%s", arg);
				printf ("Path to file holding the Gauss invariants for the window pairs is set to: %s\n", arguments->DBResultsPairsFilePath);
			} 

			break;	
		}

		case 'u':{

			if (arg == NULL) {
				printf ("Parameter useOnlyMaxMutVal_b is left at default: %d\n", arguments->useOnlyMaxMutVal_b);
				}
			else{
				arguments->useOnlyMaxMutVal_b = atoi(arg); 
				printf ("Parameter useOnlyMaxMutVal_b is set to: %d\n", arguments->useOnlyMaxMutVal_b);
			} 	

			break;
		
		}

		case 'U':{

			if (arg == NULL) {
				printf ("Parameter signed_b is left at default: %d\n", arguments->signed_b);
				}
			else{
				arguments->signed_b = atoi(arg); 
				printf ("Parameter signed_b is set to: %d\n", arguments->signed_b);
			} 	

			break;
		
		}

		case 'P':{

			if (arg == NULL) {
				printf ("Parameter rarityScorePValues_b is left at default: %d\n", arguments->rarityScorePValues_b);
				}
			else{
				arguments->rarityScorePValues_b = atoi(arg); 
				printf ("Parameter rarityScorePValues_b is set to: %d\n", arguments->rarityScorePValues_b);
			} 	

			break;
		
		}
		
		case 'A':{

			if (arg == NULL) {
				printf ("Parameter getAvgScore_b is left at default: %d\n", arguments->getAvgScore_b);
				}
			else{
				arguments->getAvgScore_b = atoi(arg); 
				printf ("Parameter getAvgScore_b is set to: %d\n", arguments->getAvgScore_b);
			} 	

			break;
		
		}
		
	}

	return 0;
}

static struct argp argp_rar0 = { options_rar0, parse_opt_rar0};







/*Arguments and options -- ie argp handling -- particular for a rarity2 scan run*/ 
typedef struct arguments_rar2{

		char backgroundDistrSingleFilePath[1000];

		char DBResultsFilePath[1000];

		int nrOfEntries;

		int matchWindowPairs_b;

		int allowedNrOfMismatches;
		

} arguments_rar2;

static struct argp_option options_rar2[] = {

		{ "backgroundDistrSingleFilePath", 'b',"CHAR", 0, "Path to file holding the back ground distribution of scores based on single windows. Mandatory." },
		{ "invsDirPath", 'i',"CHAR",0, "Full path to file holding the Gauss invariants single windows. Mandatory."},
		{ "nrOfEntries", 'e', "NUM", OPTION_ARG_OPTIONAL, "Number of invariants (Gauss integrals) to use for matching of single window values; 1: writhe, 2: abs-version of writhe too, 3 or more: all order 1 and order 2 invariants (default 1)"},
		{ "matchWindowPairs_b", 'M', "NUM", OPTION_ARG_OPTIONAL, "Boolean. Whether to include match on window pair values (1) or not (default: 0). Matching on single window values is always done."},
		{ "allowedNrOfMismatches", 'x', "NUM", OPTION_ARG_OPTIONAL, "Allowed number of mismatches in words corresponding to single window I-values (default 0)."},
		{ 0 }
    
};

static error_t parse_opt_rar2 (int key, char *arg, struct argp_state *state){

	struct arguments_rar2 *arguments = state->input;


	switch (key){

		
		case 'b':{

				if (arg == NULL) {
					printf ("Path to file holding the back ground scores (single windows) is mandatory, but not provided. Exection stops ... change to get output");
					return;
					}
				else{
					sprintf(arguments->backgroundDistrSingleFilePath, "%s", arg);
					printf ("Path to file holding the back ground scores (single windows) is set to: %s\n", arguments->backgroundDistrSingleFilePath);
				} 

				break;	
		}

		case 'i':{

			if (arg == NULL) {
				printf ("Path to file holding the Gauss invariants for single windows is mandatory, but not provided. Execution stops here ... change to get output");
				return;
				}
			else{
				sprintf(arguments->DBResultsFilePath, "%s", arg);
				printf ("Path to file holding the Gauss invariants for single windows is set to: %s\n", arguments->DBResultsFilePath);
			} 

			break;	
		}

		case 'e':{

			if (arg == NULL) {
				printf ("Number of invariants (Gauss integrals) to use for matching of single window values is left at default: %d\n", arguments-> nrOfEntries);
				}
			else{
				arguments->nrOfEntries = atoi(arg); 
				if(arguments->nrOfEntries != 0){
					printf ("Number of invariants (Gauss integrals) to use for matching of single window values is set to: %d\n", arguments->nrOfEntries);
				}
				else{
					arguments->nrOfEntries = 1; 
					printf ("Number of invariants (Gauss integrals) to use for matching of single window values was set to zero, so replaced by default: %d\n", arguments->nrOfEntries);
				}
				
			} 

			break;
		}

		case 'M':{

			if (arg == NULL) {
				printf ("Boolean matchWindowPairs_b is left at default: %d\n", arguments->matchWindowPairs_b);
				}
			else{
				arguments->matchWindowPairs_b = atoi(arg); 
				printf ("Boolean matchWindowPairs_b is set to: %d\n", arguments->matchWindowPairs_b);
			} 

			break;
		}

		case 'x':{

			if (arg == NULL) {
				printf ("Allowed number of mismatches in words corresponding to single window I-values is left at default: %d\n", arguments->allowedNrOfMismatches);
				}
			else{
				arguments->allowedNrOfMismatches = atoi(arg); 
				printf ("Allowed number of mismatches in words corresponding to single window I-values is set to: %d\n", arguments->allowedNrOfMismatches);
			} 

			break;
		}
		
	}

	return 0;
}

static struct argp argp_rar2 = { options_rar2, parse_opt_rar2 };



/*Arguments and options -- ie argp handling -- shared by all types of run ("flavours": data base run, rar0, rar 1 and rar2):*/ 
typedef struct arguments_common{

	char *argz;
	size_t argz_len;

	int halted_b; 

	int maxChainLength;

	int chunkSize;

	char flavour[10];

	char outputPath[1000];

	int use_scop_b; //if using data from SCOP (particular data structure when loading files)
	
	int use_cath_b;//if using data from CATH (particular data structure when loading files)
	/*"Positive list" of CAThj domains to consider; only relevant if use_cath_b is 1*/
	char CATHListPath[1000];

	int windowCoveringType; /*which method to use for deriving a "good" window length; used in setWindowCharacteristics*/
	int windowLength;
	int stepSize;

	//int write_windows_b;
	//int write_windowPairs_b; /*at most one of write_windows_b and write_windowPairs_b should be set to 1; corr to get_windos_b resp. get_windowPairs_b*/
	int onlyDisjointPairs_b;

	int write_chains_b;

	double thresholdMutual;
	char backgroundDistrPairsFilePath[1000];

	int nrOfEntriesForPairs;

	int allowedNrOfMismatchesPairs;

	int write_matchWindows_b; 
	int write_matchWindowPairs_b;

	int normalize_b;  /*only in its place when using more than one invariant for the scan/matching; */

	int nrOfBins;

	int writeWindowInfo_b;

	struct arguments_db arg_db;

	struct arguments_rar0 arg_rar0;

	struct arguments_rar2 arg_rar2;

} arguments_common;

static struct argp_option options_common[] = {

	{"halted_b", 'H', "NUM", 0, "Boolean. If set to 1 the execution will halt at a number of places, printing statements to screen thereby leaving a possible stop-go decision (default 0). "},
	{"maxChainLength", 'K',  "NUM", OPTION_ARG_OPTIONAL, "The maximal allowed chain-length; chains longer than this length will be disregarded (default: 1500)."},
	{"flavour", 'F', "CHAR", OPTION_ARG_OPTIONAL, "Flavour: set this according to the step in the pipeline: first run in flavour makeDB, then one of rar0, rar1 and rar2. If only the Gauss invariants are desired the flavour GI must be used. Mandatory."},
	{"outputPath", 'O',"CHAR",OPTION_ARG_OPTIONAL, "Full path to directory to which output is written."},
	{"CATHdata", 'C',"NUM",OPTION_ARG_OPTIONAL, "Whether the PDB-files are CATH domains or not (boolean, default 0)"},
	{"CATHListPath", 'p',"CHAR",OPTION_ARG_OPTIONAL, "Path to a positive list of CATH domains to be considered (only relevant when C/CATHdata is set to 1)."},
	{"SCOPdata", 'S',"NUM",OPTION_ARG_OPTIONAL, "Whether the PDB-files carry SCOP-labels or not (boolean, default 0)"},
	{"windowCoveringType", 't', "NUM", OPTION_ARG_OPTIONAL, "Which method to use for deriving a 'good' window length (0,1 or 2; default: 0)."},
	{"windowLength", 'l',  "NUM", OPTION_ARG_OPTIONAL, "Length of windows (default: 16)."},
	{"stepSize", 's',  "NUM", OPTION_ARG_OPTIONAL, "Step size (stride) for generating windows (default: 2)."},
	{"write_matchWindows_b", 'w',  "NUM", OPTION_ARG_OPTIONAL, "(Flavours: rar1/2). Boolean for whether to save obtained scores on single windows or not (default: 0)."},
	{"write_matchWindowPairs_b", 'W',  "NUM", OPTION_ARG_OPTIONAL, "(Flavours: rar1/2). Boolean for whether to save obtained scores on window pairs or not (default: 0)."},
	{"onlyDisjointPairs_b", 'j',  "NUM", OPTION_ARG_OPTIONAL, "Boolean (default: 1) for whether to (flavour: makeDB) save only invariant values on disjoint window pairs or not, or (flavour rar0/1/2) to only consider disjoint pairs in the scoring. In flavour makeDB only relevant to set to 1 when G/get_windowPairs_b is 1."},
	{"write_chains_b", 'c',  "NUM", OPTION_ARG_OPTIONAL, "Boolean (default 0) for whether to save the chains derived from the PDB-files (3d-coordianates of alpha Cs)."},
	{"thresholdMutual", 'T', "NUM", OPTION_ARG_OPTIONAL, "(Not relevant in flavour makeDB). Threshold on the absolute mutual writhe (default 5.00)."},
	{"backgroundDistrPairsFilePath", 'B',"CHAR", 0, "(Not relevant in flavour makeDB). Path to file holding the back ground distribution of scores based on window pairs. Mandatory ." },
	{"nrOfEntriesForPairs", 'E', "NUM", OPTION_ARG_OPTIONAL, "Number of invariants (Gauss integrals) to use for matching of window pair values (1: writhe, 2: abs-version of writhe too, 3 or more: all order 1 and order 2 invariants). Not relevant in flavour rar0. (default 1)"},
	{"allowedNrOfMismatchesPairs", 'X', "NUM", OPTION_ARG_OPTIONAL, "Allowed number of mismatches in words corresponding to mutual I-values in window pairs (default 0)."},
	{"normalize_b", 'z', "NUM", OPTION_ARG_OPTIONAL, "Boolean. If set to 1 the invariant values will be normalized to have mean 0 and std dev 1 (default: 0). Only in its place when using more than one invariant for the scan/matching. Not invoked (no effect) in flavour rar0."},
	{"nrOfBins", 'n', "NUM", OPTION_ARG_OPTIONAL, "Number of bins used for translation of invariant values into 'letters', ie discretization (default: 20)."},
	{"writeWindowInfo_b", 'Y', "NUM", OPTION_ARG_OPTIONAL, "Boolean. If set to 1, the score contribution from each window or window pair (depending on the scan-flavour) is written to file (default: 0). Set this to 1 only in combination with a high threshold on the writhe (T)."},
	{"chunkSize", 'k', "NUM", OPTION_ARG_OPTIONAL, "Technical parameter. Memory allocation to data base is done in chunks of this size (default: 100)."},
	{ 0 }
	
}; 


static error_t parse_opt_common (int key, char *arg, struct argp_state *state){

	struct arguments_common *arguments = (struct arguments_common*) state->input;


	switch (key){


		case 'H':{

			if (arg == NULL) {
				printf ("Boolean halted_b is left at default: %d\n", arguments->halted_b);
				}
			else{
				arguments->halted_b = atoi(arg); 
				printf ("Boolean halted_b is set to: %d\n", arguments->halted_b);
			} 

			break;
		}

		
		case 'K':{

			if (arg == NULL) {
				printf ("The maximal allowed chain-length is left at default: %d\n", arguments->maxChainLength);
				}
			else{
				arguments->maxChainLength = atoi(arg); 
				printf ("The maximal allowed chain-length is set to: %d\n", arguments->maxChainLength);
			} 

			break;
		}
	
		case 'F':{

			if (arg == NULL) {
				printf ("Flavour is left blank -- must be set. Possible values: makeDB, rar0, rar1 and rar2.\n");
				//getchar();
				//return;
				}
			else{
				sprintf(arguments->flavour, "%s", arg);
				//strcpy(arguments->flavour, arg);
				printf ("Parameter flavour is set to: %s\n", arguments->flavour);
				//getchar();
			} 

			break;	
		}



		case 'O':{

			if (arg == NULL) {
				printf ("Path to directory to which output is written is left blank ... change if you want to save the output\n");
				}
			else{
				sprintf(arguments->outputPath, "%s", arg);
				printf ("Parameter outputPath is set to: %s\n", arguments->outputPath);
			} 

			break;	
		}

		case 'C':{

			if (arg == NULL) {
				printf ("Boolean use_cath_b is left at default: %d\n", arguments->use_cath_b);
				}
			else{
				arguments->use_cath_b = atoi(arg); 
				printf ("Boolean use_cath_b is set to: %d\n", arguments->use_cath_b);
			} 

			break;
		}


		case 'p':{

			if (arg == NULL) {
				printf ("Path to a CATH positive list of domains was not provided (null).");
				}
			else{
				sprintf(arguments->CATHListPath, "%s", arg);
				printf ("Path to CATH positive list of domains is set to: %s\n", arguments->CATHListPath);
			} 

			break;
		}

		case 'S':{

			if (arg == NULL) {
				printf ("Boolean use_scop_b is left at default: %d\n", arguments->use_scop_b);
				//getchar();
				}
			else{
				arguments->use_scop_b = atoi(arg); 
				printf ("Boolean use_scop_b is set to: %d\n", arguments->use_scop_b);
				//getchar();
			} 

			break;
		}
		

		case 't':{

			if (arg == NULL) {
				printf ("Method to use for deriving a 'good' window length is left at default: %d\n", arguments->windowCoveringType);
				}
			else{
				arguments->windowCoveringType = atoi(arg); 
				printf ("Method to use for deriving a 'good' window length is set to: %d\n", arguments->windowCoveringType);
			} 

			break;
		}
		

		case 'l':{

			if (arg == NULL) {
				printf ("Window length (windowLength) is left at default: %d\n", arguments->windowLength);
				}
			else{
				arguments->windowLength = atoi(arg); 
				printf ("Window length (windowLength) is set to: %d\n", arguments->windowLength);
			} 

			break;
		}
		

		case 's':{

			if (arg == NULL) {
				printf ("Step size (stepSize) is left at default: %d\n", arguments->stepSize);
				}
			else{
				arguments->stepSize = atoi(arg); 
				printf ("Step size (stepSize) is set to: %d\n", arguments->stepSize);
			} 

			break;
		}

		case 'w':{

			if (arg == NULL) {
				printf ("Scores per window are NOT written to file (write_matchWindows_b is left at default: %d)\n", arguments->write_matchWindows_b);
				}
			else{
				arguments->write_matchWindows_b = atoi(arg); 
				if(arguments->write_matchWindows_b != 0){
					printf ("Scores per window are written to file (write_matchWindows_b is set to: %d)\n", arguments->write_matchWindows_b);
				}
				else{
					printf ("Scores per window are NOT written to file (write_matchWindows_b is set to: %d)\n", arguments->write_matchWindows_b);
				}
			} 

			break;
		}


		case 'W':{

			if (arg == NULL) {
				printf ("Scores per window pair are NOT written to file (write_matchWindowPairs_b is left at default: %d)\n", arguments->write_matchWindowPairs_b);
				}
			else{
				arguments->write_matchWindowPairs_b = atoi(arg); 
				if(arguments->write_matchWindowPairs_b != 0){
					printf ("Scores per window are written to file (write_matchWindowPairs_b is set to: %d)\n", arguments->write_matchWindowPairs_b);
				}
				else{
					printf ("Scores per window pair are NOT written to file (write_matchWindowPairs_b is set to: %d)\n", arguments->write_matchWindowPairs_b);
				}
			} 

			break;
		}


		/*case 'w':{

			if (arg == NULL) {
				printf ("Invariants on windows are NOT written to file (write_windows_b is left at default: %d)\n", arguments->write_windows_b);
				}
			else{
				arguments->write_windows_b = atoi(arg); 
				if(arguments->write_windows_b != 0){
					printf ("Invariants on windows are written to file (write_windows_b is set to: %d)\n", arguments->write_windows_b);
				}
				else{
					printf ("Invariants on windows are NOT written to file (write_windows_b is set to: %d)\n", arguments->write_windows_b);
				}
			} 

			break;
		}*/

		
		/*case 'W':{

			if (arg == NULL) {
				printf ("Invariants on window pairs are NOT written to file (write_windowPairs_b is left at default: %d)\n", arguments->write_windowPairs_b);
				}
			else{
				arguments->write_windowPairs_b = atoi(arg); 
				if(arguments->write_windowPairs_b != 0){
					printf ("Invariants on window pairs are written to file (write_windowPairs_b is set to: %d)\n", arguments->write_windowPairs_b);
				}
				else{
					printf ("Invariants on window pairs are NOT written to file (write_windowPairs_b is set to: %d)\n", arguments->write_windowPairs_b);
				}
			}

			break;
		}*/
		

		case 'j':{

			if (arg == NULL) {
				printf ("Only invariants on DISJOINT window pairs are written to file (writeOnlyDisjointPairs_b is left at default: %d)\n", arguments->onlyDisjointPairs_b);
				}
			else{
				arguments->onlyDisjointPairs_b = atoi(arg); 
				if(arguments->onlyDisjointPairs_b == 0){
					printf ("Invariants on ALL window pairs are consider for scoring/written to file (onlyDisjointPairs_b is set to: %d)\n", arguments->onlyDisjointPairs_b);
				}
				else
				{
					printf ("Only invariants on DISJOINT window pairs are considered for scoring/written to file (onlyDisjointPairs_b is set to: %d)\n", arguments->onlyDisjointPairs_b);
				}
				
			} 

			break;
		}

		case 'c':{

			if (arg == NULL) {
				printf ("Coordinates of the structures' alpha C-traces are NOT written to file (write_chains_b is left at default: %d)\n", arguments->write_chains_b);
				}
			else{
				arguments->write_chains_b = atoi(arg); 
				if(arguments->write_chains_b != 0){
					printf ("Coordinates of the structures' alpha C-traces are written to file (write_chains_b is set to: %d)\n", arguments->write_chains_b);
				}
				else
				{
					printf ("Coordinates of the structures' alpha C-traces are NOT written to file (write_chains_b is set to: %d)\n", arguments->write_chains_b);
				}
				
			} 

			break;
		}


		case 'T':{

			if (arg == NULL) {
				printf ("Threshold on absolute mutual writhe is left at default: %lf\n", arguments->thresholdMutual);
				}
			else{
				arguments->thresholdMutual = atof(arg); 
				printf ("Threshold on absolute mutual writhe is set to: %lf\n", arguments->thresholdMutual);
			} 

			break;
		}

		
		case 'B':{

			if (arg == NULL) {
				printf ("Path to file holding the back ground scores (window pairs) is mandatory, but not provided. Exection stops ... change to get output");
				return;
				}
			else{
				sprintf(arguments->backgroundDistrPairsFilePath, "%s", arg);
				printf ("Path to file holding the back ground scores (window pairs) is set to: %s\n", arguments->backgroundDistrPairsFilePath);
			} 

			break;	
		}

		case 'E':{

			if (arg == NULL) {
				printf ("Number of invariants (Gauss integrals) to use for matching of window pair values is left at default: %d)\n", arguments-> nrOfEntriesForPairs);
				}
			else{
				arguments->nrOfEntriesForPairs = atoi(arg); 
				if(arguments->nrOfEntriesForPairs != 0){
					printf ("Number of invariants (Gauss integrals) to use for matching of window pair values is set to: %d)\n", arguments->nrOfEntriesForPairs);
				}
				else
				{
					arguments->nrOfEntriesForPairs = 1; 
					printf ("Number of invariants (Gauss integrals) to use for matching of window pair values was set to zero, so replaced by default: %d)\n", arguments->nrOfEntriesForPairs);
				}
				
			} 

			break;
		}

		case 'X':{

			if (arg == NULL) {
				printf ("Allowed number of mismatches in words corresponding to mutual I-values in window pairs is left at default: %d\n", arguments->allowedNrOfMismatchesPairs);
				}
			else{
				arguments->allowedNrOfMismatchesPairs = atoi(arg); 
				printf ("Allowed number of mismatches in words corresponding to mutual I-values in window pairs is set to: %d\n", arguments->allowedNrOfMismatchesPairs);
			} 

			break;
		}

		case 'z':{

			if (arg == NULL) {
				printf ("Boolean normalize_b is left at default: %d\n", arguments->normalize_b);
				}
			else{
				arguments->normalize_b = atoi(arg); 
				printf ("Boolean normalize_b is set to: %d\n", arguments->normalize_b);
			} 

			break;
		}


		case 'n':{

			if (arg == NULL) {
				printf ("Parameter nrOfBins setting the number of bins used for translating of invariant values into 'letters', ie discretization, is left at default: %d\n", arguments->nrOfBins);
				}
			else{
				arguments->nrOfBins = atoi(arg); 
				printf ("Parameter nrOfBins setting the number of bins used for translating of invariant values into 'letters', ie discretization, is set to: %d\n", arguments->nrOfBins);
			} 

			break;
		}	

		case 'k':{

			if (arg == NULL) {
				printf ("Memory allocation to data base is done in chunks of (default): %d\n", arguments->chunkSize);
				}
			else{
				arguments->chunkSize = atoi(arg); 
				printf ("Memory allocation to data base is done in chunks of: %d\n", arguments->chunkSize);
			} 

			break;
		}	


		case 'Y':{

			if (arg == NULL) {
				printf ("Score contributions per window or window pair (depending on scan flavour) not written to file (writeWindowInfo_b option is left at default) %d\n", arguments->writeWindowInfo_b);
				}
			else{
				arguments->writeWindowInfo_b = atoi(arg); 
				printf ("Option for writing score contributions per window or window pair (depending on scan flavour) is set to: %d\n", arguments->writeWindowInfo_b);
			} 

			break;
		}	


		case ARGP_KEY_ARG:{
			argz_add (&arguments->argz, &arguments->argz_len, arg); 
			break;
		}

		case ARGP_KEY_INIT:{

			state->input = &arguments;

			arguments->argz = NULL;
			arguments->argz_len = 0;

			
			arguments->halted_b = 0;
			arguments->maxChainLength = 1500; //chains longer than this will be discarded; can be set acc to RAM availability
			arguments->chunkSize = 100;
			sprintf(arguments->flavour,"%s", "");
			strcpy(arguments->outputPath, "");	
			arguments->use_scop_b = 0;
			arguments->use_cath_b = 0;
			strcpy(arguments->CATHListPath, "");
			arguments->windowCoveringType = 0;
			arguments->windowLength = 16;
			arguments->stepSize = 2;
			arguments->onlyDisjointPairs_b = 1;
			arguments->write_chains_b = 0;
			arguments->thresholdMutual = 5;
			strcpy(arguments->backgroundDistrPairsFilePath, "");
			arguments->nrOfEntriesForPairs = 1;
			arguments->allowedNrOfMismatchesPairs = 0;
			arguments->write_matchWindows_b  = 0;
			arguments->write_matchWindowPairs_b  = 0;
			arguments->normalize_b = 0;
			arguments->nrOfBins = 20;
			arguments->writeWindowInfo_b = 0;

			strcpy(arguments->arg_db.dirPath , "");
			strcpy(arguments->arg_db.DBName, "DB");

			arguments->arg_db.order = 1;
			arguments->arg_db.incl_abs_b = 0;
			arguments->arg_db.full_b = 0;
			arguments->arg_db.get_windows_b = 0;
			arguments->arg_db.get_windowPairs_b = 0;
		
			//Default values for rar0 specific variables:
			strcpy(arguments->arg_rar0.queriesDirPath , "");
			strcpy(arguments->arg_rar0.queriesName , "");
			strcpy(arguments->arg_rar0.DBResultsPairsFilePath, "");

			arguments->arg_rar0.useOnlyMaxMutVal_b = 1;
			arguments->arg_rar0.signed_b = 1;
			arguments->arg_rar0.rarityScorePValues_b = 0;
			arguments->arg_rar0.getAvgScore_b = 1;


			//Default values for rar1/2 specific variables:
			strcpy(arguments->arg_rar2.backgroundDistrSingleFilePath, "");
			strcpy(arguments->arg_rar2.DBResultsFilePath, "");

			arguments->arg_rar2.nrOfEntries = 1;
			arguments->arg_rar2.matchWindowPairs_b = 0;
			arguments->arg_rar2.allowedNrOfMismatches = 0;
			

			//assign
			state->child_inputs[0] = &arguments->arg_db;
			state->child_inputs[1] = &arguments->arg_rar0;
			state->child_inputs[2] = &arguments->arg_rar2;
			break;
		}
	
		//break;

		// case ARGP_KEY_ARG:{
		//   printf("arg num: %d\n", state->arg_num);
		//   if (state->arg_num >= 2){
		//     /* Too many arguments. */
		// 	printf("Too many args\n");
		//     argp_usage (state);
		//   }
		//   arguments->args[state->arg_num] = arg;
		// }
		// break;

		// case ARGP_KEY_END:{
		// 	printf("arg num: %d\n", state->arg_num);
		//   	if (state->arg_num < 2){
		//     /* Not enough arguments. */
		// 		printf("Not enough args\n");
		// 	    argp_usage(state);
		// 	}
		// }
		// break;

		default:
			return ARGP_ERR_UNKNOWN;


	}


	
	return 0;
}


static struct argp_child children_parsers[] ={

		{&argp_db, 0, "Options specific for creation of a data base of Gauss invariant values on single windows or window pairs for a set of structures (flavour: makeDB):", 7 },
		{&argp_rar0, 0, "Options specific for a rarity scan of type 0 (flavour: rar0):", 7 },
		{&argp_rar2, 0, "Options specific for a rarity scan type 2 (flavour: rar2):", 7 },
		{0}

};

static struct argp argp_common = { options_common, parse_opt_common, 0, 0, children_parsers};


int main (int argc, char **argv){


	int returnVal = 0;
	
	//int maxChainLength = 1500; /*chains longer than this will be discarded; can be set acc to RAM availability*/
	
	
	/*Params for detecting and recording closed loops; serve as placeholders here:*/
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
	//int use_scop_b = 0; //if using data from SCOP (particular data structure when loading files)
	//int use_cath_b = 0;//if using data from CATH (particular data structure when loading files)

	int loadFromSubDirs_b = 1; //if PDB files sit in sub-directories at some level of the "root" directory, dirName, use loadFromSubDirs_b = 1; else 0.

	int omitSingletons_b = 0;

	//char DBName[100] = "top100";

	//char dirPath[1000] = "/home/tkj375/masters_study_pgm/projects/knots/code/top100H";

	//char CATHListPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH-v2_4_Kolodny/subset_list_web.txt";

	//char outputPath[1000]  = "/isdata/kroghgrp/tkj375/GISA_unix/results/v8/top100";

	char fileNameChains[1000] = "/chains_top100.txt";

	/*For getting and writing out the invariant values for a set windowlength:*/
	//int windowCoveringType = 0; /*which method to use for deriving a "good" window length; used in setWindowCharacteristics*/
	//int get_windows_b = 1;
	//int get_windowPairs_b = 0; /*at most one of get_windows_b and get_windowPairs_b should be set to 1*/
	//int windowLength = 16;
	//int stepSize = 2;
	int write_windows_b = 0;
	int write_windowPairs_b = 0; /*at most one of write_windows_b and write_windowPairs_b should be set to 1; corr to get_windos_b resp. get_windowPairs_b*/
	int maxNrOfWindows = 1;
	//int write_chains_b = 0;
	char pathChainsDB[2000] = "\0";
	char fileNameOutChainsDB[1000] =  "/chains_";
	char fileNameOutExtensionChains[20] = ".txt";

	/*Parameters specific to GI_windows:*/
	int writeOnlyDisjointPairs_b;


	/*Parameters specific to rarity0:*/
	//char queriesName[100] = "top8000"; //"CATH_test1";
	//char queriesDirPath[1000] = "/isdata/kroghgrp/tkj375/data/structural/CATH/CATH_test1";

	int chunkSize = 100; /* load in chunks of this size (nr of structures) */

	//char DBResultsPairsFilePath[1000] = "/isdata/kroghgrp/tkj375/GISA_unix/results/v9/top8000/Invariants_Pairs_windowlgth_20_4_order_1_computeGI_windows_top8000_winCovType0_onlyDisjointPairs1.txt";

	//double thresholdMutual = 5.0;

	//int useOnlyMaxMutVal_b = 1;

	int absMutValScorePValues_b = 0; /*default: if set to 1 and useOnlyMaxMutVal_b = 0, provide a absMutValScoresDistr*/

	char absMutValScoresDistrName[100] = "_pValuesIncl"; //_top100Matchtop100";
	char absMutValScoresDistrFilePath[1000] = ""; 

	char pathChainsQueries[2000] = "\0";
	
	//int onlyDisjointPairs_b = 1;

	int write_scoresPairs_b = 0; //set to 1 if scores should be output
	
	/*for rarity0: END*/	


	/*Parameters specific for rarity1/2:*/
	char rarityScorePairsDistrName[100] = "_pValuesIncl"; //_top100Matchtop100";
	char rarityScorePairsDistrFilePath[1000];

	int nrOfWindowsForMatch = 1; /*just leave it at 1, has no influence apart from appearing in file names*/

	//int allowedNrOfMismatchesPairs = 1; /*allowed nr of mismatches in word matching, window pairs case*/
	//int nrOfEntriesForPairs = 1; /*  setting this to 1 means that only the mutual writhe is considered; >1: mutuals of higher order inv's too*/

	/*For generating a set of bins*/
	int binningType = 1; /*don't use type 2 since this makes all bins equally populated (not what is wanted in a rarity search); don't use type 0 either, unless you want to go directly int the code and change the hardcoded values ...*/
	//int nrOfBins = 20;
	double compensationForRounding = 0.0000009; //When using equi-percentile binning (binningType 2) bin interval end points will be shifted to the right by this amount; to prevent separation of e.g. self-hits

	int write_matchWindowPairs_b = 0;/*write the match results for each window pair to file or not (0) */

	/*for rarity1: END*/	


	/*Parameters specific for rarity2:*/

	char rarityScoreDistrName[100] = "_pValuesIncl"; //_top100Matchtop100";
	char rarityScoreDistrFilePath[1000];

	char DBResultsFilePath[1000];

	//int allowedNrOfMismatches; /*allowed nr of mismatches in word matching, single window case*/

	int write_matchWindows_b = 0; /*write the match results for each window to file or not (0) */

	int write_scores_b = 1; //set to 1 if scores should be output

	/*for rarity2: END*/	


	struct arguments_common arguments;

	struct argp argp_common = { options_common, parse_opt_common, 0, 0, children_parsers};

	/*Init:
	strcpy(arguments.arg_db.dirPath , "");
	strcpy(arguments.arg_rar0.queriesDirPath , "");
	strcpy(arguments.arg_db.DBName, "DB");
	strcpy(arguments.arg_rar0.queriesName, "QS");
	strcpy(arguments.arg_rar0.DBResultsPairsFilePath, "");

	strcpy(arguments.backgroundDistrPairsFilePath , "");

	//Default values for shared variables:
	//strcpy(arguments.flavour, "(saaledes stadig ganske hemmelig)");
	arguments.halted_b = 0;
	arguments.maxChainLength = 1500; //chains longer than this will be discarded; can be set acc to RAM availability
	arguments.chunkSize = 100;
	strcpy(arguments.outputPath, "");
	strcpy(arguments.CATHListPath, "");
	arguments.use_cath_b = 0;
	arguments.use_scop_b = 0;
	arguments.arg_db.order = 1;
	arguments.arg_db.incl_abs_b = 0;
	arguments.arg_db.full_b = 0;
	arguments.windowCoveringType = 0;
	arguments.arg_db.get_windows_b = 0;
	arguments.arg_db.get_windowPairs_b = 0;
	arguments.windowLength = 16;
	arguments.stepSize = 2;
	//arguments.write_windows_b = 0;
	//arguments.write_windowPairs_b = 0;
	arguments.write_matchWindows_b  = 0;
	arguments.write_matchWindowPairs_b  = 0;
	arguments.onlyDisjointPairs_b = 1;
	arguments.write_chains_b = 0;
	arguments.thresholdMutual = 5;
	arguments.normalize_b = 0;
	arguments.nrOfBins = 20;
	

	//Default values for rar0 specific variables:
	arguments.arg_rar0.useOnlyMaxMutVal_b = 0;
	arguments.arg_rar0.rarityScorePValues_b = 0;


	//Default values for rar1/2 specific variables:
	arguments.arg_rar2.nrOfEntries = 1;
	arguments.arg_rar2.matchWindowPairs_b = 0;
	arguments.arg_rar2.allowedNrOfMismatches = 0;
	*/



	printf("----------------------------------------------------------------------------------------------------------------------------\n");
	printf("Parsing the command line arguments .... \n");
	printf("... the collected parameters are shown for you to check that the parser has picked up everything correctly ...\n");
	printf("... i.e that the settings are as you wanted (the actual parameter values may (will) be repeated below)\n");
	printf("... The actual parameter values will be repeated below so you have another chance to check it all.\n");
	printf("----------------------------------------------------------------------------------------------------------------------------\n");
	printf("\n");

	argp_parse(&argp_common, argc, argv, 0, 0, &arguments);

	printf("\n");



	printf("----------------------------------------------------------------------------------------------------------------------------\n");
	printf("General settings:\n");
	printf("----------------------------------------------------------------------------------------------------------------------------\n");
	printf("\n");

	printf("The flavour is: %s\n", arguments.flavour);
	
	printf("maxChainLength is set to %d; chains longer than this will be disregarded (can be set acc to RAM availability).\n", arguments.maxChainLength);
	printf("\n");
	printf("Memory allocation to data base is done in chunks of %d\n", arguments.chunkSize);
	printf("\n");


	printf("----------------------------------------------------------------------------------------------------------------------------\n");
	printf("Common parameters are set as follows:\n");
	printf("----------------------------------------------------------------------------------------------------------------------------\n");
	printf("\n");
	
	printf("Output path: %s\n", arguments.outputPath);

	printf("use_cath_b: %d\n", arguments.use_cath_b);
	if (arguments.CATHListPath != ""){
		printf("Path to CATH positive list: %s\n", arguments.CATHListPath);
		if(arguments.use_cath_b == 0){
			printf("As the use_cath_b  parameter is set to 0, the positive list is not used.\n");
		}
	}
	printf("use_scop_b: %d\n", arguments.use_scop_b);

	printf("windowCoveringType: %d\n", arguments.windowCoveringType);
	
	printf("windowLength: %d\n", arguments.windowLength);
	printf("stepSize: %d\n", arguments.stepSize);
	printf("write_windows_b: %d\n", write_windows_b);
	printf("write_windowPairs_b: %d\n", write_windowPairs_b);
	
	printf("write_chains_b is: %d\n", arguments.write_chains_b);

	
	if (arguments.write_chains_b !=0){
		/*generate file name for writing out the data base chains:*/
		strcat(pathChainsDB , arguments.outputPath);
		strcat(pathChainsDB , fileNameOutChainsDB);
		if (arguments.arg_db.DBName != ""){
			strcat(pathChainsDB , arguments.arg_db.DBName );
		}
		strcat(pathChainsDB , fileNameOutExtensionChains);

		printf("DB chains will be written to: %s\n", pathChainsDB);

		//clear the file by opening it in w-mode and then closing it again
		fclose(fopen(pathChainsDB , "w"));
	}

	//printf("If parameter values are fine, press any key to carry on.\n");
	//getchar();

	
	if (strcmp(arguments.flavour, "makeDB") == 0){

		printf("\n");
		printf("----------------------------------------------------------------------------------------------------------------------------\n");
		printf("I'll be running a DB-run with the following parameter values.\n"); 
		printf("----------------------------------------------------------------------------------------------------------------------------\n");
		printf("\n");
	
		if( strcmp(arguments.arg_db.dirPath,"") == 0){
		
				printf ("Path to directory holding the queries PDB files is mandatory, but not provided. Execution stops ... change to get output.\n");
				returnVal = 1;
				return;
			}
			else{
				printf("PDB-files directory (dirPath) is: %s\n", arguments.arg_db.dirPath);
		}
		

		if (strcmp(arguments.arg_db.DBName, "") != 0){

			printf("Name of data base is set to: %s\n", arguments.arg_db.DBName);
			printf("\n");
		}
		else{

			strcpy(arguments.arg_db.DBName, "DB");
			printf("Name of data base was set to '', but has been replaced by the default: %s\n", arguments.arg_db.DBName);
			printf("\n");
		
		}

		printf("order: %d\n", arguments.arg_db.order);
		printf("incl_abs_b: %d\n", arguments.arg_db.incl_abs_b);
		printf("full_b: %d\n", arguments.arg_db.full_b);


		printf("get_windows_b is: %d\n", arguments.arg_db.get_windows_b);
		printf("get_windowPairs_b is: %d\n", arguments.arg_db.get_windowPairs_b);

		if(arguments.arg_db.order > 2 && arguments.arg_db.get_windowPairs_b ==1){
			printf("With get_windowPairs_b set to 1, the order parameter can be at most 2. Execution stops here.\n");
			return;
		}
		
		if(arguments.arg_db.get_windows_b == arguments.arg_db.get_windowPairs_b ){
			printf("Set exactly one of get_windows_b and get_windowpairs_b to 1. Execution stops here.\n");
			return;
		}
		else if(arguments.arg_db.get_windows_b == 1){

			printf("Invariants' values on single windows will be computed and written to file.\n");
			write_windows_b = 1;
		
		}
		else if(arguments.arg_db.get_windowPairs_b == 1){

			printf("Invariants' values on window pairs will be computed and written to file.\n");
			write_windowPairs_b = 1;
		
		}
		
		

		/*if(arguments.write_windows_b == 1 && arguments.arg_db.get_windows_b == 0){
			printf("Set get_windows_b to 1 when write_windows_b is set to 1 (else nothing to write to file).");
		}*/

		/*if(arguments.write_windowPairs_b == 1 && arguments.arg_db.get_windowPairs_b == 0){
			printf("Set get_windowPairs_b to 1 when write_windowPairs_b is set to 1 (else nothing to write to file).");
		}*/

		writeOnlyDisjointPairs_b = arguments.onlyDisjointPairs_b;
		printf("writeOnlyDisjointPairs_b: %d\n", writeOnlyDisjointPairs_b);
		if(writeOnlyDisjointPairs_b == 1 && write_windowPairs_b == 0){
			printf("Value of writeOnlyDisjointPairs_b has no effect since write_windowPairs_b is set to 0.");
		}

		printf("If parameter values are fine, press any key to carry on.\n");
		getchar();

		returnVal = computeGI_windows(arguments.arg_db.DBName, arguments.arg_db.dirPath, arguments.use_scop_b, 
					arguments.use_cath_b, arguments.CATHListPath, loadFromSubDirs_b, 
					arguments.outputPath, arguments.maxChainLength, omitSingletons_b ,
					arguments.arg_db.order, arguments.arg_db.incl_abs_b, arguments.arg_db.full_b, 
					split_b, write_final_b, write_all_b, 
					closed_loops_b, closedLoopLength, closedLoopDist, 
					pokeLength, invValSubChainPairs_b, writeSubChainPairsAll_b, 
					subChainLength, arguments.windowCoveringType, arguments.arg_db.get_windows_b,
					arguments.arg_db.get_windowPairs_b, arguments.windowLength, arguments.stepSize, 
					write_windows_b, write_windowPairs_b, writeOnlyDisjointPairs_b, 
					maxNrOfWindows, arguments.write_chains_b, pathChainsDB,
					print_b, print_basic_b);

	}
	else if (strcmp(arguments.flavour, "rar0") == 0){
			
		printf("\n");
		printf("----------------------------------------------------------------------------------------------------------------------------\n");
		printf("I'll be running a rarity scan, type 0. This is based on mutual writhe values for pairs of windows.\n");
		printf("----------------------------------------------------------------------------------------------------------------------------\n");
		printf("\n");

		if(arguments.halted_b == 1){
			printf("The execution will be halted in a number of places (halted_b is set to 1).\n");
			printf("\n");
		}

		if(arguments.arg_rar0.useOnlyMaxMutVal_b == 1){

			printf("Only the segment/window pair having highest mutual writhe pair (per structure) will be used (useOnlyMaxMutVal_b is set to 1).\n");
			if(arguments.arg_rar0.signed_b == 1){
			
				printf("The positive and negative writhe cases will be handled separately (signed_b is set to 1)");

			}
			else{

				printf("Only the absolute value of the mutual writhe will be considered, i.e. the positive and negative writhe cases will NOT be handled separately (signed_b is set to diff from 1)");

			}

		}
		else if(arguments.arg_rar0.useOnlyMaxMutVal_b == 0){

			printf("All segment/window pairs having abs mutual writhe above the set threshold will be used (useOnlyMaxMutVal_b is set to 0).\n");

		}
			
		if(arguments.onlyDisjointPairs_b == 1){
			printf("Only disjoint pairs of windows are used.\n");
		}
		else{
			printf("All window pairs are considered, ie not only the disjoint pairs.\n");
		}
	

		if (strcmp(arguments.arg_db.DBName, "") != 0){

			printf("Name of data base is set to: %s\n", arguments.arg_db.DBName);
			printf("\n");
		}
		else{

			strcpy(arguments.arg_db.DBName, "DB");
			printf("Name of data base was set to '', but has been replaced by the default: %s\n", arguments.arg_db.DBName);
			printf("\n");
		
		}

		printf("DB results directory (DBResultsPairsFilePath) is: %s\n", arguments.arg_rar0.DBResultsPairsFilePath);

		printf("Queries' PDB-files directory (queriesDirPath) is: %s\n", arguments.arg_rar0.queriesDirPath);

		if (strcmp(arguments.arg_rar0.queriesName, "") != 0){
			printf("Name of queries is set to: %s\n", arguments.arg_rar0.queriesName);
			printf("\n");
		}
		else{
			strcpy(arguments.arg_rar0.queriesName, "QS");
			printf("Name of queries was set to '', but has been replaced by the default: %s\n", arguments.arg_rar0.queriesName);
			printf("\n");

		}

		printf("thresholdMutual: %lf\n", arguments.thresholdMutual);


		write_scoresPairs_b = 1; //scores should always be written to file
	
		absMutValScorePValues_b = arguments.arg_rar0.rarityScorePValues_b;
		printf("absMutValScorePValues_b: %d\n", arguments.arg_rar0.rarityScorePValues_b);

		if(absMutValScorePValues_b == 1 && strcmp(arguments.backgroundDistrPairsFilePath, "") == 0 ){

			printf("absMutValScorePValues_b is 1 while the back ground file for the scores is ''. Change to get output (execution stops here).\n");
			return 1;

		}
		else{

			strcpy(absMutValScoresDistrFilePath, arguments.backgroundDistrPairsFilePath);
			printf("The back ground scores file is: %s\n", absMutValScoresDistrFilePath);

		}

		if(arguments.arg_rar0.getAvgScore_b == 0){

			printf("Scores are NOT averaged over number of pairs in the query structure.\n");

		}
		else{

			printf("Scores are obtained by averaging over number of pairs in the query structure.\n");
		
		}


		if(arguments.writeWindowInfo_b == 0){

			printf("Score contributions per window pair are not written to file (writeWindowInfo_b is set to 0).\n"); 

		}
		else{

			printf("Score contributions per window pair are written to file (writeWindowInfo_b is not 0).\n"); 

		}

		printf("Only the mutual writhe is considered! (the parameter order is therefore set to 1 and incl_abs_b is set to 0).\n"); 

		if(arguments.normalize_b == 1)
			printf("The boolean normalize_b is set to 1, so the invariant (the writhe) will be normalized to mean 0 and std dev 1.\n");
		else{
			printf("The boolean normalize_b is set to 0, so the invariant (the writhe) will NOT be normalized.\n");
		}


		
		printf("If all settings are fine, press any key to carry on.\n"); 
		getchar();

		returnVal = rawRarity0(arguments.halted_b, arguments.arg_db.DBName, arguments.arg_rar0.queriesName, 
						arguments.arg_rar0.queriesDirPath, arguments.arg_rar0.DBResultsPairsFilePath, pathChainsQueries,
						absMutValScoresDistrName, absMutValScoresDistrFilePath, arguments.onlyDisjointPairs_b, 
						arguments.arg_rar0.useOnlyMaxMutVal_b, arguments.arg_rar0.signed_b, absMutValScorePValues_b, 
						arguments.thresholdMutual, arguments.use_scop_b, arguments.use_cath_b, 
						arguments.CATHListPath, loadFromSubDirs_b, arguments.normalize_b, 
						arguments.outputPath, arguments.maxChainLength, arguments.chunkSize, 
						arguments.windowCoveringType, arguments.windowLength, arguments.stepSize, 
						arguments.arg_rar0.getAvgScore_b, write_windowPairs_b, write_scoresPairs_b, 
						arguments.write_chains_b, arguments.writeWindowInfo_b, print_b, 
						print_basic_b);


	} 
	else if (strcmp(arguments.flavour, "rar1") == 0){

		printf("\n");
		printf("----------------------------------------------------------------------------------------------------------------------------\n");
		printf("I'll be running a rarity scan, type 1, for pairs of windows.\n"); 
		printf("----------------------------------------------------------------------------------------------------------------------------\n");

		if(arguments.halted_b == 1){
			printf("The execution will be halted in a number of places (halted_b is set to 1).\n");
			printf("\n");
		}

		if (strcmp(arguments.arg_db.DBName, "") != 0){

			printf("Name of data base is set to: %s\n", arguments.arg_db.DBName);
			printf("\n");
		}
		else{

			strcpy(arguments.arg_db.DBName, "DB");
			printf("Name of data base was set to '', but has been replaced by the default: %s\n", arguments.arg_db.DBName);
			printf("\n");
		
		}

		printf("DB results directory (DBResultsPairsFilePath) is: %s\n", arguments.arg_rar0.DBResultsPairsFilePath);


		if(arguments.arg_rar0.rarityScorePValues_b == 1){

			if(strcmp(arguments.backgroundDistrPairsFilePath, "") == 0 ){	
				printf("Obtaining p-values is called (rarityScorePValues_b is 1), but the path for the back ground file for the scores is left blank ''. Change to get output (execution stops here).\n");
				return 1;

			}
			else{

				strcpy(rarityScorePairsDistrFilePath, arguments.backgroundDistrPairsFilePath);
				printf("Obtaining p-values will be based on the back ground scores file: %s\n", rarityScorePairsDistrFilePath);

			}
		}

		printf("thresholdMutual: %lf\n", arguments.thresholdMutual);

		write_scoresPairs_b = 1;

		printf("Number of invariants used in  matching (nrOfEntriesForPairs): %d\n", arguments.nrOfEntriesForPairs);
		printf("The allowed number of mismatches is set to: %d\n", arguments.allowedNrOfMismatchesPairs);
		printf("The order, incl_abs_b and full_b parameters are therefore set as follows:\n");
		if (arguments.nrOfEntriesForPairs == 1){

			arguments.arg_db.order = 1;
			arguments.arg_db.full_b = 0;
			arguments.arg_db.incl_abs_b = 0;

		}
		else if(arguments.nrOfEntriesForPairs == 2){

			arguments.arg_db.order = 1;
			arguments.arg_db.full_b = 0;
			arguments.arg_db.incl_abs_b = 1;


		}
		else if(arguments.nrOfEntriesForPairs > 2){

			arguments.arg_db.order = 2;
			arguments.arg_db.full_b = 1;
			arguments.arg_db.incl_abs_b = 1;

		}
		printf("order: %d\n", arguments.arg_db.order);
		printf("incl_abs_b: %d\n", arguments.arg_db.incl_abs_b);
		printf("full_b: %d\n", arguments.arg_db.full_b);



		if(arguments.normalize_b == 1)
			printf("The boolean normalize_b is set to 1, so the invariants will be normalized to mean 0 and std dev 1.\n");
		else{
			printf("The boolean normalize_b is set to 0, so the invariants will NOT be normalized.\n");
		}

		printf("Parameter nrOfBins setting the number of bins used for translating of invariant values into 'letters', ie discretization, is: %d\n ", arguments.nrOfBins);

		if(arguments.writeWindowInfo_b == 0){

			printf("Score contributions per window pair are not written to file (writeWindowInfo_b is set to 0).\n"); 

		}
		else{

			printf("Score contributions per window pair are written to file (writeWindowInfo_b is not 0).\n"); 

		}

		if(arguments.write_matchWindowPairs_b == 1){

			printf("Score per window pair will be written to file.\n");

		}
		else{

			printf("Score per window pair will not be written to file.\n");

		}


		printf("If all settings are fine, press any key to carry on.\n"); 
		getchar();

		returnVal = rawRarity1(arguments.halted_b, rarityScorePairsDistrFilePath, rarityScorePairsDistrName, 
						arguments.arg_db.DBName, arguments.arg_rar0.queriesName, arguments.arg_rar0.queriesDirPath, 
						arguments.arg_rar0.DBResultsPairsFilePath, pathChainsQueries, arguments.onlyDisjointPairs_b,
						arguments.thresholdMutual, arguments.arg_rar0.rarityScorePValues_b, arguments.use_scop_b, 
						arguments.use_cath_b, arguments.CATHListPath, loadFromSubDirs_b, 
						arguments.normalize_b, arguments.outputPath, arguments.maxChainLength, 
						arguments.nrOfEntriesForPairs, binningType, arguments.nrOfBins, 
						compensationForRounding, arguments.allowedNrOfMismatchesPairs, nrOfWindowsForMatch, 
						arguments.chunkSize, arguments.arg_db.order, arguments.arg_db.incl_abs_b, 
						arguments.arg_db.full_b, arguments.write_matchWindowPairs_b, arguments.windowCoveringType, 
						arguments.windowLength, arguments.stepSize, write_windowPairs_b, 
						write_scoresPairs_b, arguments.write_chains_b, arguments.writeWindowInfo_b, 
						print_b, print_basic_b);




			
	}
	else if (strcmp(arguments.flavour, "rar2") == 0){

		printf("\n");
		printf("----------------------------------------------------------------------------------------------------------------------------\n");
		printf("I'll be running a rarity scan, type 2, for single/pairs of windows.\n"); 
		printf("----------------------------------------------------------------------------------------------------------------------------\n");
		printf("\n");

		if(arguments.halted_b == 1){
			printf("The execution will be halted in a number of places (halted_b is set to 1).\n");
			printf("\n");
		}

		
		if (strcmp(arguments.arg_db.DBName, "") != 0){

			printf("Name of data base is set to: %s\n", arguments.arg_db.DBName);
			printf("\n");
		}
		else{

			strcpy(arguments.arg_db.DBName, "DB");
			printf("Name of data base was set to '', but has been replaced by the default: %s\n", arguments.arg_db.DBName);
			printf("\n");
		
		}


		printf("DB results directory for the single windows (DBResultsFilePath) is: %s\n", arguments.arg_rar2.DBResultsFilePath);

		write_scores_b = 1;
		printf("Number of invariants used in single-window matching (nrOfEntries): %d\n", arguments.arg_rar2.nrOfEntries);
		printf("The allowed number of mismatches in the single-window matching is set to: %d\n", arguments.arg_rar2.allowedNrOfMismatches);
		
		printf("The order, incl_abs_b and full_b parameters are therefore set as follows:\n");
		if (arguments.arg_rar2.nrOfEntries == 1){

			arguments.arg_db.order = 1;
			arguments.arg_db.full_b = 0;
			arguments.arg_db.incl_abs_b = 0;

		}
		else if(arguments.arg_rar2.nrOfEntries == 2){

			arguments.arg_db.order = 1;
			arguments.arg_db.full_b = 0;
			arguments.arg_db.incl_abs_b = 1;


		}
		else if(arguments.arg_rar2.nrOfEntries > 2){

			arguments.arg_db.order = 2;
			arguments.arg_db.full_b = 1;
			arguments.arg_db.incl_abs_b = 1;

		}
		printf("order: %d\n", arguments.arg_db.order);
		printf("incl_abs_b: %d\n", arguments.arg_db.incl_abs_b);
		printf("full_b: %d\n", arguments.arg_db.full_b);

		if(arguments.arg_rar2.matchWindowPairs_b != 0){
			write_scoresPairs_b = 1;
			printf("Matching based on pairs of windows will be done too!");
			printf("thresholdMutual: %lf\n", arguments.thresholdMutual);
			printf("DB results directory for the pairs (DBResultsPairsFilePath) is: %s\n", arguments.arg_rar0.DBResultsPairsFilePath);
			printf("Number of invariants used in matching window pair values (nrOfEntriesForPairs): %d\n", arguments.nrOfEntriesForPairs);
			printf("The allowed number of mismatches in the pair-window matching is set to: %d\n", arguments.allowedNrOfMismatchesPairs);


			if(arguments.nrOfEntriesForPairs > arguments.arg_rar2.nrOfEntries){

				printf("You have set the nrOfEntriesForPairs (%d) larger than nrOfEntries (%d), which is not viable.\n", arguments.nrOfEntriesForPairs, arguments.arg_rar2.nrOfEntries);
				printf("The nrOfEntriesForPairs is therefore set equal to the set nrOfEntries!\n");
				arguments.nrOfEntriesForPairs = arguments.arg_rar2.nrOfEntries;

			}

			if(arguments.write_matchWindowPairs_b == 1){

				printf("Score per window pair will be written to file.\n");

			}
			else{

				printf("Score per window pair will not be written to file.\n");

			}

		}
		else{
			
			printf("Matching (and scoring) will only be based on single windows.\n");

			if(arguments.write_matchWindows_b == 1){

				printf("Score per window will be written to file.\n");

			}
			else{

				printf("Score per window will not be written to file.\n");

			}

		}
		


		if(arguments.arg_rar0.rarityScorePValues_b == 1){

			if(arguments.arg_rar2.matchWindowPairs_b != 0){
				if(strcmp(arguments.backgroundDistrPairsFilePath, "") == 0 || strcmp(arguments.arg_rar2.backgroundDistrSingleFilePath, "") == 0 ){	
					
					printf("Obtaining p-values is called (rarityScorePValues_b is 1), but at least one of the paths to the back ground files for the scores is left blank ''. Change to get output (execution stops here).\n");
					return 1;
				
				}
				else{
					
					strcpy(rarityScoreDistrFilePath, arguments.arg_rar2.backgroundDistrSingleFilePath);
					printf("Obtaining p-values for the single window matching will be based on the back ground scores file: %s\n", rarityScoreDistrFilePath);

					strcpy(rarityScorePairsDistrFilePath, arguments.backgroundDistrPairsFilePath);
					printf("Obtaining p-values for the pair window matching will be based on the back ground scores file: %s\n", rarityScorePairsDistrFilePath);

				}
			}
			else { //only single window matching

					if(strcmp(arguments.arg_rar2.backgroundDistrSingleFilePath, "") == 0 ){	
						
						printf("Obtaining p-values is called (rarityScorePValues_b is 1), but the path to the back ground files for the scores is left blank ''. Change to get output (execution stops here).\n");
						return 1;
					
					}
					else{
						
						strcpy(rarityScoreDistrFilePath, arguments.arg_rar2.backgroundDistrSingleFilePath);
						printf("Obtaining p-values for the single window matching will be based on the back ground scores file: %s\n", rarityScoreDistrFilePath);

					}

			}
		}
		
		if(arguments.normalize_b == 1){

			printf("The boolean normalize_b is set to 1, so the invariants will be normalized to mean 0 and std dev 1.\n");

		}
		else{

			printf("The boolean normalize_b is set to 0, so the invariants will NOT be normalized.\n");

		}

		printf("Parameter nrOfBins setting the number of bins used for translating of invariant values into 'letters', ie discretization, is: %d\n ", arguments.nrOfBins);

		if(arguments.writeWindowInfo_b == 0){

			printf("Score contributions per window pair are not written to file (writeWindowInfo_b is set to 0).\n"); 

		}
		else{

			printf("Score contributions per window pair are written to file (writeWindowInfo_b is not 0).\n"); 

		}


		printf("If all settings are fine, press any key to carry on.\n"); 
		getchar();

		returnVal = rawRarity2(arguments.halted_b, rarityScoreDistrFilePath, rarityScoreDistrName, 
								rarityScorePairsDistrFilePath, rarityScorePairsDistrName, arguments.arg_db.DBName, 
								arguments.arg_rar0.queriesName, arguments.arg_rar0.queriesDirPath, arguments.arg_rar2.DBResultsFilePath, 
								arguments.arg_rar0.DBResultsPairsFilePath, pathChainsQueries, arguments.thresholdMutual, arguments.onlyDisjointPairs_b, 
								arguments.arg_rar0.rarityScorePValues_b, arguments.use_scop_b, arguments.use_cath_b, 
								arguments.CATHListPath, loadFromSubDirs_b, arguments.normalize_b, 
								arguments.outputPath,arguments.maxChainLength, arguments.arg_rar2.nrOfEntries, 
								arguments.nrOfEntriesForPairs, binningType, arguments.nrOfBins, 
								compensationForRounding, arguments.arg_rar2.allowedNrOfMismatches, arguments.allowedNrOfMismatchesPairs, 
								nrOfWindowsForMatch, arguments.chunkSize, arguments.arg_db.order, 
								arguments.arg_db.incl_abs_b, arguments.arg_db.full_b, arguments.write_matchWindows_b, 
								arguments.arg_rar2.matchWindowPairs_b, arguments.write_matchWindowPairs_b, arguments.windowCoveringType, 
								arguments.windowLength, arguments.stepSize, write_windows_b, 
								write_windowPairs_b, write_scores_b, write_scoresPairs_b, 
								arguments.write_chains_b, arguments.writeWindowInfo_b, 
								print_b, print_basic_b);

		//getchar();

	}

	//printf("Sille, jaj ælsker daj ... Du er skoer, Topper\n");
	//getchar();

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