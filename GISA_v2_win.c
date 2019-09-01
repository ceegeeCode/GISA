/* **********************************************************************************************************************************
*
*
* ***************************************   GISA *******************************************************************************
*
*
*
*     Present version: v2, 8th May, 2017.
*
**************************************************************************************************************************************
* 
*     C-implementation of algorithm for computing Gauss integral based invariants for protein fold desciption. Code aimed for 
*     searching for "rare geometries" is included.
*
*     Author: Christian Grønbæk
*
*     The code is available under the GNU Public License v3 (the license text is appended at the end of the coding section below). 
*
*     If using this code cite this source:
*     C.Grønbæk., T.Hamelryck, P.Røgen: "GISA: Using Gauss Integrals to identify rare conformations in protein structures", TO APPEAR.
*     
*     Notes:
*     The function to use for computing the invariants is main_GISA, which in turn calls the core 
*     function, computeGI. Global memory allocation is used (memory is allocated up front); the code allows 
*     running the computations across a set (directory) of PDB-files; memory is then allocated to cope with 
*     the longest chain in the set. 
*     To use main_GISA set the input parameters as desired in the code (of the function); compile 
*     the code and, if compilation was successfull, run it (the function "main" consists in a call to 
*     main_GISA). Names of output files are generated inside computeGI, by concatention of some hard-coded
*     strings and values of external parameters; changing the hard-coded strings you need to modify the code
*     inside computeGI (see also below).

*     The code for main_GISA is preceded by a description/documentation. An outline of the code 
*     can be found in the Supplementary Data accompanying the Application Note mentioned above ("GIsquared: 
*     Using Gauss integrals to find rare local geometries in protein structures", Bioinformatics, 2017.) 
*     
*     This version of the code compiles in Windows (the author used Microsoft Visual Studio 12); for Unix 
*     compilation, the version GISA_v2_unix.c should be used.
*
*     Else, at its present state the C-code is admittedly not overly user friendly. Extending the code 
*     with a shell so as to call the compiled code from a command prompt is lacking. Also some checks 
*     and guiding response could be added (e.g. that a memory allocation was successful or that a parameter 
*     value was set improperly for the run). As for how the names of the files to which results are written 
*     out (if desired) are created, it is probably most easy to simply read the code: the names are put
*     togehter from set strings and values of parameters. As an example check e.g. the naming of the file
*     used if a write-all results are chosen; for this see the if-clause "if (write_all_b == 1){...}" in 
*     computeGI.
*
*	  OBS: pointers are allocated memory with explicit cast for the sake of avoiding troubles with missing 
*     casts if code is augmented with e.g. CUDA C for gpu-computing. 

**********************************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<windows.h>


#define pi acos(-1.0)
#define pathLength 2048
#define fileNameLength 100


/*struct for holding content of directory*/
typedef struct dirContent{

	int numberOfFiles;
	char ** ptr_dirList;
	char ** ptr_fileNameList;

} dirContent;

/*struct for keeping rough info about possible multi-mer PDB-file: the number of chains in the structure
and an array indexed by chain number and with entries (chainId, chainlength), ptr_chainInStr: */

typedef struct chainInfo{

	char *chainId; /*we shall assume that the chain id is always a single character (may be a blank)*/ 
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
	double ** I1234_full_aid;
	double ** I1324_full_aid;
    double ** I1234_full;
    double ** I1324_full;
	/* assisting the order 3 case:*/
	double ** I1324_full2_aid;
	double ** I1324_full2;
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


/**************************************************
************ Prototype functions ******************
***************************************************/

/*function for getting the contents of a named directory*/
struct dirContent ListDirectoryContents(char *);

/*functions for reading content of a PDB-file and outputting the array of C-alpha coordinates
The first, readPDBChainStructure, gets the number of chains in the file (structure) and the length of 
each if the chains and returns this info in an array {...,(chainNr, chainLength), ...}. The 
second, main_readPDB, reads in the array of C-alpha coordinates in a given file (structure) 
along with chain length and chain number.*/
struct chainInStructure readPDBChainStructure(FILE *);

struct cAlpha * main_readPDB(FILE *, int, int);


/*function computing the square of the distance between two C-alphas; used for finding closed loops
in a structure; we compute the square dist rather than the dist simply to avoid a sqrt (which is
unneccessary for our purposes)*/
double distCalphas(struct cAlphaCoords, struct cAlphaCoords);


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

int writheSubChainPairs(char *, struct segment *, int , struct I_ptr, int, char *, struct cAlpha *, int, double);

int genInvariantSubChainPairs(char *, struct segment *, int, struct I_ptr, double **, char *, int, int);

/*Functions for writing out the results to .txt files*/

int collectIvalues(struct I_values **, struct I_ptr, int, int, double, double);


int writeChain(char *, char *,  char *, struct cAlpha *, int);

int writeIvaluesToFile(char *, struct I_values **, int, int);

int writeAllIvaluesToFile(char *, int, struct I_ptr, int);


/*Code blocks for making code in main function (below) more transparent; take care of memory 
allocation, initialization and freeing memeory:*/

int alloc_init_I_values(struct I_values ***, int, int);

int alloc_init_I_measures(struct I_ptr *, int, int, int, int, int, struct twoSegmentIndex ** );


int alloc_init_segment_ptr(struct Segment_ptr *, int, int);

int init_segment_ptr(struct Segment_ptr, struct segment *, int, int *);


int init_I_measures(struct I_ptr, int, int, int);


/*
** functions for freeing memory allocated to pointers
*/
void freeDblArray(double **, size_t);

void free_I_measures(struct I_ptr, int, int, int);

void free_I_values(struct I_values **, int );

/*main program blocks*/

int computeGI(char , char *, char *, char *, int , int , int , int , int , int , 
			  int , int, int , int , double, double, int , int , int , double,  
			  int , int );



/****************************************************
*****************Coding section *********************
*****************************************************/


/*Directory walk for Windows OS
Parts inspired heavily by a post on http://stackoverflow.com by "NTDLS".
Function for getting the contents of a named directory*/
struct dirContent ListDirectoryContents(char *sDir){
    
	WIN32_FIND_DATA fdFile;
    HANDLE hFind = "*.*";

	int fileNumber = 0;
	char * sPath; /*for making the code more readable; will point to address with value of fileName*/
	char * fileName; /*to contain the file name, without the path in front, that is.*/
	char * fileNameWExt;
	char * ext;
	char ** ptr_dirList; /*pointer to contain path-filename list*/
	char ** ptr_fileNameList; /*pointer to contain filename list*/
	int i = 0;

	/*for holding output:*/
	struct dirContent output; 
	/*initialize:*/
	output.numberOfFiles = 0;
	output.ptr_dirList = NULL;
	output.ptr_fileNameList = NULL;

	sPath = (char *)calloc(pathLength, sizeof(char));
	fileName = (char *)calloc(fileNameLength, sizeof(char));
	fileNameWExt = (char *)calloc(fileNameLength, sizeof(char));
	ext = (char *)calloc(4, sizeof(char));

	//sprintf(sPath, "%s//*.*", sDir); //first directory "is" sDir//*.* 
	sprintf(sPath, "%s\\*.*", sDir); //first directory "is" sDir//*.* 


    if((hFind = FindFirstFile(sPath, &fdFile)) == INVALID_HANDLE_VALUE){
        printf("Sorry - directory not found: %s\n", sDir);
        return output;
    }

	/*loop through the directory and count the number of files we want to consider*/
    do{
        //Find first file will always return "."
        //    and ".." as the first two directories.
		// strcmp: ==0 when strings match
		// Also: avoid .txt files
        if(strcmp(fdFile.cFileName, ".") != 0
                && strcmp(fdFile.cFileName, "..") != 0){
            //Build up our file path using the passed in
            //  [sDir] and the file/foldername we just found:
            //sprintf(sPath, "%s//%s", sDir, fdFile.cFileName); //current sPath string "is" sDir//fdFile.cFileName            
			sprintf(sPath, "%s\\%s", sDir, fdFile.cFileName); //current sPath string "is" sDir//fdFile.cFileName 

            //Is the entity a file or a folder?
            if(fdFile.dwFileAttributes &FILE_ATTRIBUTE_DIRECTORY){
                printf("Pick files also in sub-directory: %s\n", sPath);
                ListDirectoryContents(sPath); //Recursion; note that, when the recursion halts, control returns to where the do-while where the recursion was initiated!
				//getchar();
            }
            else{
                //printf("File: %s\n", sPath);
				fileNumber += 1;
            }
        }
    }
    while(FindNextFile(hFind, &fdFile)); //Find the next file.

	FindClose(hFind);

	printf("Number of files found in this dir walk: %d\n", fileNumber);

	/*loop through the directory again and collect the file names*/
	/*first allocate memory to our output pointer:*/
	ptr_dirList = (char **)malloc(fileNumber*sizeof(sPath));
	ptr_fileNameList = (char **)malloc(fileNumber*sizeof(fileName));

	/*init*/
	for (i= 0; i < fileNumber; i++){
			 ptr_dirList[i] = (char *)calloc(pathLength, sizeof(char));
			 ptr_fileNameList[i] = (char *)calloc(fileNameLength, sizeof(char));
	}


	/*Loop through the directory again and record data in ptr*/

	//Reset.
	sPath = ptr_dirList[0]; /* let sPath use the first available address of ptr_dirList */
	fileName = ptr_fileNameList[0]; /* let fileName use the first available address of ptr_fileNameList */
	//sprintf(sPath, "%s//*.*", sDir); //sPath now becomes the directory sDir//*.*, where sDir is the directory lead to in the do-while above via the recursion!
	sprintf(sPath, "%s\\*.*", sDir); //sPath now becomes the directory sDir\*.*, where sDir is the directory lead to in the do-while above via the recursion!

	//printf("sPath: %s\n", sPath);

    hFind = FindFirstFile(sPath, &fdFile);

	fileNumber = 0;
	/*In this do-while the if(fdFile.dwFileAttributes &FILE_ATTRIBUTE_DIRECTORY)-clause is actually
	obsolete, since sPath is now a directory in which only files are found (if any)*/
	do{
        //Find first file will always return "."
        //    and ".." as the first two directories.
        if(strcmp(fdFile.cFileName, ".") != 0
                && strcmp(fdFile.cFileName, "..") != 0){
            //Build up our file path using the passed in
            //  [sDir] and the file/foldername we just found:
            //sprintf(sPath, "%s\\%s", sDir, fdFile.cFileName);
			//sprintf(sPath, "%s//%s", sDir, fdFile.cFileName);
			sprintf(sPath, "%s\\%s", sDir, fdFile.cFileName);
			sprintf(fileNameWExt, "%s",fdFile.cFileName);
			sscanf(fileNameWExt, "%s.%s",fileName, ext);
			//printf("fileName: %s ext: %s\n", fileName, ext);

            //Is the entity a File or Folder?
            if(fdFile.dwFileAttributes &FILE_ATTRIBUTE_DIRECTORY){
                printf("Pick files also in sub-directory: %s\n", sPath);
                ListDirectoryContents(sPath); //Recursion
            }
            else{
                /*printf("File: %s\n", sPath);*/
				fileNumber += 1;
				sPath = ptr_dirList[fileNumber]; /*get sPath a new address: the next in line*/
				fileName = ptr_fileNameList[fileNumber]; /*get fileName a new address: the next in line*/
            }
        }
    }
    while(FindNextFile(hFind, &fdFile)); //Find the next file.

	FindClose(hFind);

	/*
	for (i=0;i<10;i++){
		printf("check :%s ",ptr_dirList[i]);
	}
	*/

	output.numberOfFiles = fileNumber;
	output.ptr_dirList = ptr_dirList;
	output.ptr_fileNameList = ptr_fileNameList;

	return output; /*true;*/
}

struct dirContent ListDirectoryContents_old(char *sDir){
    
	WIN32_FIND_DATA fdFile;
    HANDLE hFind = "*.*";

	int fileNumber = 0;
	char * sPath; /*for making the code more readable; will point to address with value of fileName*/
	char ** ptr_dirList; /*pointer to contain filename list*/
	int i = 0;

	/*for holding output:*/
	struct dirContent output; 
	/*initialize:*/
	output.numberOfFiles = 0;
	output.ptr_dirList = NULL;

	sPath = (char *)calloc(2048, sizeof(char));

	sprintf(sPath, "%s//*.*", sDir);


    if((hFind = FindFirstFile(sPath, &fdFile)) == INVALID_HANDLE_VALUE){
        printf("Sorry - directory not found: %s\n", sDir);
        return output;
    }

	/*loop through the directory and count the number of files we want to consider*/
    do{
        //Find first file will always return "."
        //    and ".." as the first two directories.
		// strcmp: ==0 when strings match
		// Also: avoid .txt files
        if(strcmp(fdFile.cFileName, ".") != 0
                && strcmp(fdFile.cFileName, "..") != 0){
            //Build up our file path using the passed in
            //  [sDir] and the file/foldername we just found:
            sprintf(sPath, "%s//%s", sDir, fdFile.cFileName);

            //Is the entity a file or a folder?
            if(fdFile.dwFileAttributes &FILE_ATTRIBUTE_DIRECTORY){
                printf("Pick files also in sub-directory: %s\n", sPath);
                ListDirectoryContents(sPath); //Recursion
            }
            else{
                /*printf("File: %s\n", sPath);*/
				fileNumber += 1;
            }
        }
    }
    while(FindNextFile(hFind, &fdFile)); //Find the next file.

	FindClose(hFind);

	//printf("number of files found in dir walk: %d", fileNumber);

	/*loop through the directory again and collect the file names*/
	/*first allocate memory to our output pointer:*/
	ptr_dirList = (char **)malloc(fileNumber*sizeof(sPath));
	for (i=0; i < fileNumber; i++){
			 ptr_dirList[i] = (char *)calloc(2048, sizeof(char));
	}

	//Reset.
	sPath = ptr_dirList[0]; /* let sPath use the first available address of ptr_dirList */
	sprintf(sPath, "%s//*.*", sDir);


    hFind = FindFirstFile(sPath, &fdFile);

	/*loop through the directory and count the number of files we want to consider*/
	fileNumber = 0;
	do{
        //Find first file will always return "."
        //    and ".." as the first two directories.
        if(strcmp(fdFile.cFileName, ".") != 0
                && strcmp(fdFile.cFileName, "..") != 0){
            //Build up our file path using the passed in
            //  [sDir] and the file/foldername we just found:
            //sprintf(sPath, "%s\\%s", sDir, fdFile.cFileName);
			sprintf(sPath, "%s//%s", sDir, fdFile.cFileName);


            //Is the entity a File or Folder?
            if(fdFile.dwFileAttributes &FILE_ATTRIBUTE_DIRECTORY){
                printf("Pick files also in sub-directory: %s\n", sPath);
                ListDirectoryContents(sPath); //Recursion
            }
            else{
                /*printf("File: %s\n", sPath);*/
				fileNumber += 1;
				sPath = ptr_dirList[fileNumber]; /*get sPath a new address: the next in line*/
            }
        }
    }
    while(FindNextFile(hFind, &fdFile)); //Find the next file.

	FindClose(hFind);

	/*
	for (i=0;i<10;i++){
		printf("check :%s ",ptr_dirList[i]);
	}
	*/

	output.numberOfFiles = fileNumber;
	output.ptr_dirList = ptr_dirList;

	return output; /*true;*/
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

		/*allocate memory for chain structure pointer*/
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

struct cAlpha * main_readPDB(FILE * ptr_file, int chainNumber, int chainLength){
    
		char buf[1000];
		char *chainId;
		char chainIdInit[3] = "Ø"; /*some "unlikely" value for init of the chain identifier*/
		int chainNr = -1;
		int  n =0;
		int resNr = 0;
		int resNrPrior = -999;

		struct cAlpha *ptr_chain = NULL;

		chainId = (char *) malloc(sizeof(char)*3);

		/*allocate memory for array of coordinates*/
		ptr_chain = (struct cAlpha *) calloc (chainLength, sizeof(struct cAlpha));

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
				//printf("Res nr:%d\n", resNr);
				if (strstr(buf, "ATOM") == buf && buf[13] == 'C' && buf[14] == 'A' && (resNr != resNrPrior|| buf[26] == 'A' || buf[26] == 'B')){ /*Only read-in CA's; we also include "additional" residues marked with postfix A or B to the resNr */
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
			}
		}

		//printf("Populated chain Length: %d\n", n);


		/*printf("Chain: \n");
		for ( n=0; n<10; n++ ){
			printf("x: %lf y:%lf z:%lf \n", ptr_chain[n].x, ptr_chain[n].y, ptr_chain[n].z );
		}*/

    	return ptr_chain;

	}

/*function computing the square of the distance between two Calphas; used for finding closed loops
in a structure; we compute the square dist rather than the dist simply to avoid a sqrt (which is
unneccessary for our purposes)*/
double distCalphas(struct cAlphaCoords cA1, struct cAlphaCoords cA2){

	double d = 0;

	d = (cA1.x - cA2.x)*(cA1.x - cA2.x) + (cA1.y - cA2.y)*(cA1.y - cA2.y) + (cA1.z - cA2.z)*(cA1.z - cA2.z);

	return d;

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

	/* check if segments are adjacent, i.e. that vectors v has length "zero":*/
	if (l13 < 1e-16){
			return 0;
	}
	if (l23 < 1e-16){
			return 0;
	}
	if (l24 < 1e-16){
			return 0;
	}
	if (l14 < 1e-16){
			return 0;
	}

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
    
	/*call the computation*/
	return w(s101, s102, s103, s111, s112, s113, s201, s202, s203, s211, s212, s213);
	
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
	double wVal2;

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
						I1234_full_aid[i][j] += I12[b+1][j]*wVal1;
                        I1324_full_aid[i][j] += I12[b][j]*wVal1;
					}
                    /*I1234_full add up terms and recursion: I1234(i,j)*/
                    I1234_full[i][j] = I1234_full_aid[i][j] -I1234_full_aid[i][j-1];
                    I1234_full[i][j] += I1234_full[i][j-1] + I1234_full[i+1][j] - I1234_full[i+1][j-1];
                    /*I1324_full add up terms recursion: I1324(i,j)*/
                    I1324_full[i][j] = I1324_full_aid[i][j-1] - I1324_full_aid[i][j];
                    I1324_full[i][j] += (I12[i][j-1] - I12[i+1][j-1])*(I12[i+1][j] - I12[i+1][j-1]) + I1324_full[i][j-1] + I1324_full[i+1][j] - I1324_full[i+1][j-1];
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
	they are merely there for convenience -- readaility:*/

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


/* As aggrAndW but with a search for closed loops added: inside the traversal of the simplex (taking part of the 
recursion for computing the invariants) a 1d-loop is inserted traversing indices above the current "i" to search 
for the first place where a closed-loop distance criterium is fulfilled (if at any index). 
We neglect loops in which a segment has "too long" length ie which consists in an artificial peptide bond introduced 
by the way we load structures (as one long chain disregarding if it consists of several sub-chains (this can be 
avoided if changing the load procedure). A standard length of sqare root of 20 is used for this purpose. Also we 
neglect loops which are "too short": the minLoopLength allows varying this.*/
int aggrAndW_wClosedLoops(struct segment *ptr_segment, int chainLen, int order, int full_b, struct I_ptr I_measures, struct twoSegmentIndex *ptr_closedLoopInd, int closedLoopLength, double closedLoopDist, int * ptr_nrOfClosedLoops){

	int L = 0;
	int i = 0;
	int j = 0;
	int b = 0;
	int k = 0;
	int c = 0;
	int minLoopLength = 7;
	int loopTooLongInd = 0;
	double stdRealSegLength = 20; /*close to the square of 4.5: to avoid sub-chains with segments that are "strangely long" */
	int n = 0;
	int last_j = 0;
	
	int cnt = 0; /*for counting and indexing closed loops*/ 

	/* the parts of the output structure to hold the output:*/
	double wValue;
	double abs_wValue; 
	double wVal1;
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
			
			/*Find possible closed loops starting at first coord of segment[i].
			We discard loops of minLoopLength by only checking/admitting closedness for 
			segment index (n) at or above i+minLoopLength; from i to i+minLoopLength-1 
			we only check if the segments are "too long". 
			Also, we want to find the "longest irreducible loop". Since we reverse 
			in the i-loop, we try to extend a loop by adding to its "i head": if we 
			for i have found a closed loop terminating at n, we will not search further at 
			i (as we have then set last_j = n) and in the next step (at i-1) we will ultimately 
			check if (i-1,n-1) is also closed (but this will only happen if we have not hit 
			a loop (i-1,n') for some n' lower than n*/
			for (n = i; n <= min(i+ closedLoopLength,last_j-1); n++){
				/*if segment n is "too long" we skip the rest of the search at this i:*/
				if (distCalphas(ptr_segment[n].s1, ptr_segment[n].s2) > stdRealSegLength){
					break;
				} 
				//printf("i:%d n:%d last_j:%d\n", i,n, last_j);
				if (n >= i + minLoopLength && distCalphas(ptr_segment[i].s1, ptr_segment[n].s2) < closedLoopDist){
					ptr_closedLoopInd[cnt].idx1 = i;
					ptr_closedLoopInd[cnt].idx2 = n;
					//printf("Seg i: %d 1st pt is x y z : %lf %lf %lf\n", i, ptr_segment[i].s1.x,ptr_segment[i].s1.y, ptr_segment[i].s1.z );
					//printf("Seg n: %d 2nd pt is x y z : %lf %lf %lf\n", n, ptr_segment[n].s2.x,ptr_segment[n].s2.y, ptr_segment[n].s2.z );
					//printf("At i:%d n:%d dist is: %lf\n", i,n, distCalphas(ptr_segment[i].s1, ptr_segment[n].s2));
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
			
			/*Find possible closed loops starting at first coord of segment[i].
			We discard loops of minLoopLength by only checking/admitting closedness for 
			segment index (n) at or above i+minLoopLength; from i to i+minLoopLength-1 
			we only check if the segments are "too long". 
			Also, we want to find the "longest irreducible loop". Since we reverse 
			in the i-loop, we try to extend a loop by adding to its "i head": if we 
			for i have found a closed loop terminating at n, we will not search further at 
			i (as we have then set last_j = n) and in the next step (at i-1) we will ultimately 
			check if (i-1,n-1) is also closed (but this will only happen if we have not hit 
			a loop (i-1,n') for some n' lower than n*/
			for (n=i; n <= min(i+ closedLoopLength,last_j-1); n++){
				/*if segment n is "too long" we skip the rest of the search at this i:*/
				if (distCalphas(ptr_segment[n].s1, ptr_segment[n].s2) > stdRealSegLength){
					break;
				}
				//printf("i:%d n:%d last_j:%d\n", i,n, last_j);
				if (n >= i+minLoopLength && distCalphas(ptr_segment[i].s1, ptr_segment[n].s2) < closedLoopDist){
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
			
			/*Find possible closed loops starting at first coord of segment[i].
			We discard loops of minLoopLength by only checking/admitting closedness for 
			segment index (n) at or above i+minLoopLength; from i to i+minLoopLength-1 
			we only check if the segments are "too long". 
			Also, we want to find the "longest irreducible loop". Since we reverse 
			in the i-loop, we try to extend a loop by adding to its "i head": if we 
			for i have found a closed loop terminating at n, we will not search further at 
			i (as we have then set last_j = n) and in the next step (at i-1) we will ultimately 
			check if (i-1,n-1) is also closed (but this will only happen if we have not hit 
			a loop (i-1,n') for some n' lower than n*/
			for (n=i; n <= min(i+ closedLoopLength,last_j-1); n++){
				/*if segment n is "too long" we skip the rest of the search at this i:*/
				if (distCalphas(ptr_segment[n].s1, ptr_segment[n].s2) > stdRealSegLength){
					break;
				}
				//printf("i:%d n:%d last_j:%d\n", i,n, last_j);
				if (n >= i+minLoopLength && distCalphas(ptr_segment[i].s1, ptr_segment[n].s2) < closedLoopDist){
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
	double stdRealSegLength = 20; /*close to square of 4.5: to avoid sub-chains with segments that are "strangely long" */
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
		

		//printf("idx 1: %d,%d\n", i_1, j_1);

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
		poke is longer than the std length of a segment (~square of about 4-5 Ångstrøm). Obs:
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

		returnVal = writeChain(ptr_fileNameOutChains, I_measures.fileName, I_measures.chainId, ptr_chain, *I_measures.chainLen);

	}

	fclose(ptr_fileOut);

	return returnVal;
}

/*Simple fct which computes the mutual writhe of all (disjoint) pairs of sub-chains; the sub-chains 
considered all have the same fixed length, subChainLength (typically between 15 and 30). Writes for 
each pair the mutual writhe value to the file defined by ptr_fileNameOut, but only if writhe in absolute 
value is above the set threshold, thresholdWrithe.
Obs: unlike the restricted search (using closed loops, see examineClosedLoops) we here do not filter 
out "artificial" cases in which a segment of a sub-chain is "too long" (such cases do exist in 
pdb-files, but are rare)*/
int writheSubChainPairs(char *ptr_fileNameOut, struct segment *ptr_segment, int chainLength, struct I_ptr I_measures, int subChainLength, char *ptr_fileNameOutChains, struct cAlpha *ptr_chain, int writeAll_b, double thresholdWrithe){

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

	int writeChain_b = 0; 

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

			if(writeAll_b == 1 && fabs(I_subChainPair) > thresholdWrithe){
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

				/*when pair of sub-chains exists (as here) we want to write out the chain too:*/
				writeChain_b +=1;

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
		if(fabs(Imin) > thresholdWrithe){
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

			/*when pair of sub-chains exists (as here) we want to write out the chain too:*/
			writeChain_b +=1;

		}

		/*write Imax result to file:*/
		if(fabs(Imax) > thresholdWrithe){
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

			/*when pair of sub-chains exists (as here) we want to write out the chain too:*/
			writeChain_b +=1;
		}
	}

	/*write out the chain when called for:*/
	if(writeChain_b > 0){
		returnVal = writeChain(ptr_fileNameOutChains, I_measures.fileName, I_measures.chainId, ptr_chain, *I_measures.chainLen);
	}

	fclose(ptr_fileOut);

	return returnVal;
};

/*Version of writheSubChainPairs that does the same but for a generic invariant.
Additional input: a specific measure, e.g. I_measure.I1324, as "I_measure" and its name, I_measure_name, which in the 
example would be the string "I1324"; if the specified measure is I12 the fct works as writheSubChainPairs except that 
it does not write out the chains. And:
Obs: there is no set threshold for the value of the invariant limiting which cases to be written out (so all are written out)*/
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

	returnVal = 1;

	return returnVal;
};


/*Functions for writing out to files:*/

/*convenience fct's for writing out values:*/
int collectIvalues(struct I_values **ptr_I_values, struct I_ptr I_measures, int subStructureNr, int perturbationNumber, double chainTime, double compTime){

	int order = *I_measures.order;
	int L = *I_measures.chainLen - 1;

	ptr_I_values[subStructureNr][perturbationNumber].fileName = I_measures.fileName;
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

int writeChain(char *ptr_fileNameOut, char *fileName, char *chainId, struct cAlpha * ptr_chain, int chainLen){

	FILE *ptr_fileOut;

	int i;

	/*results will be written to this file. We use a-mode since we are calling this writing-fct
	repeatedly (in main):*/
	ptr_fileOut = fopen(ptr_fileNameOut, "a");

	fprintf(ptr_fileOut, "%s;%s;%s;%s\n",
		"PDB_file", 
		"chain_id",
		fileName,
		chainId
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
			"chain id",
			"chain nr",
		    "str length",
			"order",
			"perturbation number",
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
			"PDB_file", 
			"chain id",
			"chain nr",
		    "str length",
			"order",
			"perturbation nr");
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
int alloc_init_I_measures(struct I_ptr *ptr_I_measures, int order, int full_b, int chainNr, int chainLen, int closed_loops_b, struct twoSegmentIndex ** ptr_ptr_closedLoopInd){

	int i = 0;
	int j = 0;
	int L = chainLen - 1;
	int cnt = 0; /*only used if closed_loops_b == 1*/

	struct I_ptr I_measures;

	/*if closed_loops_b ==1 this ptr is in use:*/
	struct twoSegmentIndex *ptr_closedLoopInd;

	I_measures.fileName = (char *) malloc (sizeof(char));
	I_measures.fileName = "NN";

	I_measures.chainId = (char *) malloc (sizeof(char));
	I_measures.chainId = "NN";
	I_measures.chainNr = &chainNr;
	I_measures.chainLen = &chainLen;
	I_measures.order = &order;
	 
	if (order == 0){

		/*allocate memory for double array of w-values */
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


	if (order == 1){

		/*allocate memory for double arrays of w and I12 values (order <= 1 measures)*/
		I_measures.wVal = (double **) malloc ((L+1)*sizeof(double *));
		I_measures.I12 = (double **) malloc ((L+1)*sizeof(double *));
			
		I_measures.Ia12 = (double **) malloc ((L+1)*sizeof(double *));
			
		for (i=0; i <= L; i++){
			I_measures.wVal[i] = (double *)malloc((L+1) * sizeof(double));
			I_measures.I12[i] = (double *)malloc((L+1) * sizeof(double));
			I_measures.Ia12[i] = (double *) malloc ((L+1)*sizeof(double));
			/* Initialize to zero everywhere: */
			for (j=0; j <= L; j++){
				I_measures.wVal[i][j] = 0.0;
				I_measures.I12[i][j] = 0.0;
				I_measures.Ia12[i][j] = 0.0;
				cnt +=1;

			}
		}
	}


	if (order == 2){

		/*allocate memory for double arrays of for measures of order <= 2: */
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
				}

				cnt +=1;
			}
		}
	}


	if (order == 3){

		/*allocate memory for double arrays of for measures of order <= 3: */
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

			I_measures.I1324_full2_aid[i] = (double *)malloc((L+1) * sizeof(double));
			I_measures.I1324_full2[i] = (double *)malloc((L+1) * sizeof(double));
			I_measures.I1423_full0[i] = (double *)malloc((L+1) * sizeof(double));
			I_measures.I1423_full2_aid[i] = (double *)malloc((L+1) * sizeof(double));
			I_measures.I1423_full2[i] = (double *)malloc((L+1) * sizeof(double));


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

	*ptr_I_measures = I_measures;

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

	if (order == 0){
		for (i = 0; i < L+1; i++){
			free(I_measures.wVal[i]);
		}
		free(I_measures.wVal);
	}

	if (order == 1){
		for (i = 0; i < L+1; i++){
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
}

void free_I_values(struct I_values **ptr_I_values, int numberOfFiles){

	int i = 0;
	int j = 0;

	for (i=0; i< numberOfFiles; i++){
		free(ptr_I_values[i]);
	}
	free(ptr_I_values);
}


/*Main program blocks. Change the main by changing the fct called, chosing between
the main_* main-functions below (currently only one).*/

int main(){

	int returnVal = 0;

	returnVal = main_GISA();

	return returnVal;

}

/* This functions merely wraps a call to computGI; conveneint for setting and changing parameter values
For the meaning of each parameter see the comment block preceding computGI:*/
int main_GISA(){

	int returnVal = 0;
	
	int i = 0;

	/*compute invariants up to and including this order:*/
	int order = 1;
	/*inlcude absolute value versions of order one and two invariants:*/
	int incl_abs_b = 1;
	/* compute full order 2 measures across the simplex*/
	int full_b = 0;

	/*Detect and record closed loops; if chosen I12 values corr 
	to the loop will be stored for alter use (finding "pokes").
	Can only be used with order >= 1 (since I12 is needed):*/
	int closed_loops_b = 0;
	int closedLoopLength = 30;
	int pokeLength = 10;
	double closedLoopDist = 7*7; /*square of distance in nr of Ångstrøm -- we use square of dist*/
	double thresholdLinks = 0; //0.9*4*pi;
	double thresholdPokes = 0; //0.8*4*pi;
	/*For computing and recording the mutual invariant values (e.g writhe) of pairs of sub-chains (1) or not (0); aimed
	at more generally searching for particular geometric shapes than with the closed-loops
	search:*/
	int invValSubChainPairs_b = 0; /*if set to 1: computes the mutual writhe; if set to 2 computes the mutual invariant value for order 1 and order 2 "full" invariant (takes setting order = 3 or full_b = 1); if set to 3: all invariants (takes setting order = 3)*/
	int writeSubChainPairsAll_b = 0; /*if set to 1 the results for all sub-chain pairs will be attempted to be written out (may be lots of data!)*/
	int subChainLength = 15; /* the lenght of the sub-chains in the pairs considered*/
	double thresholdWrithe = 0; //4*pi;

	/*split versions: w computation split away from aggregation:*/
	int split_b = 0;
	/* print various intermediates on the screen (1):*/
	int print_b = 0;
	/*print put basic stuff to the screen:*/
	int print_basic_b = 0;
	/*write out to .txt file the values of the measures on the given chain (i.e. the "top corner" value for each chain)*/
	int write_final_b = 1;
	/*write out to .txt file the values of the measures on the given chain
	and all its sub-chains. Obs: this can become very large data sets*/
	int write_all_b = 0;


	/*name of folder holding the PDB files (the input):*/
	char DBName[100] = "top100";
	//char DBName[100] = "top8000";
	//char DBName[100] = "top8000_aug";
	//char dirName[1000] = "C:\\Users\\Christian\\BioInformatics\\projects\\knots\\code\\1ifcH";
	//char dirName[1000] = "C:\\Users\\Christian\\BioInformatics\\projects\\knots\\code\\2er7H";
	//char dirName[1000] = "C:\\Users\\Christian\\BioInformatics\\projects\\knots\\code\\3chains";
	//char dirName[1000] = "C:\\Users\\Christian\\BioInformatics\\masters_study_pgm\\projects\\knots\\code\\SMLX";
	//char dirName[1000] = "C:\\Users\\Christian\\BioInformatics\\projects\\knots\\code\\testSet22";
	char dirPath[1000] = "C:\\Users\\Christian\\BioInformatics\\masters_study_pgm\\projects\\knots\\code\\top100H";
	//char dirName[1000] = "E:\\Thesis\\data\\top8000_chains_70";
	//char dirPath[1000] = "C:\\Users\\Christian\\BioInformatics\\masters_study_pgm\\projects\\knots\\code\\top8000_chains_70";
	//char dirPath[1000] = "C:\\Users\\Christian\\BioInformatics\\masters_study_pgm\\projects\\knots\\code\\top8000_chains_70_aug";

	char outputPath[1000]  = "C:\\Users\\Christian\\Bioinformatics\\papers\\geometric_anomalies_in_folds\\c_code\\results_C_5th\\double_precision";

	//char fileNameChains[1000] = "\\chains_subChainPairs_top100.txt";
	char fileNameChains[1000] = "\\chains_closed_loops_top100.txt";

	//char fileNameChains[1000] = "\\chains_closed_loops_top8000.txt";
	//char fileNameChains[1000] = "\\chains_subChainPairs_top8000.txt";
	//char fileNameChains[1000] = "\\chains_subChainPairs_top8000_aug.txt";
	//char fileNameChains[1000] = "\\chains_subChainPairs_top8000_thrL09P08.txt";
	//char fileNameChains[1000] = "\\chains_closed_loops_top8000_thrL09P08.txt";

	for(i=2; i<3; i++){

		if(i==0){
			order = 1; 
			closed_loops_b = 0; 
			invValSubChainPairs_b = 0;
		}
		if(i==1){
			order = 2; 
			closed_loops_b = 1; 
			invValSubChainPairs_b = 0;
		}
		if(i==2){
			order = 3; 
			closed_loops_b = 0; 
			invValSubChainPairs_b = 1;
		}

		returnVal = computeGI(DBName, dirPath, outputPath, fileNameChains, order, incl_abs_b, full_b, split_b,write_final_b, write_all_b, 
					closed_loops_b, closedLoopLength, closedLoopDist, pokeLength, thresholdLinks, thresholdPokes, invValSubChainPairs_b, writeSubChainPairsAll_b, subChainLength, thresholdWrithe, 
					print_b, print_basic_b);
	
	}

	getchar();

	return returnVal;

}
 
/*Main function, computeGI: 

Computes all invariants up to and including the desired order (e.g. order = 2) for 
all PDB-files in a directory (defined by dirPath). Results are written out as desired: 
for each PDB-file, all results (write_all_b = 1) meaning the invariants' values for all 
vertices in the simplex or only final values (write_final_b =1) meaning the invariants' values 
of the given structure (PDB-file), i.e. the value at the top left corner of the simplex. As an 
additional feature it is possible to envoke a particular version of the computation of the 
invariants in which closed loops are detected and recorded, along with the writhe value (I12) 
of the loop. For running this set closed_loops_b to 1. A more general version of this is also 
callable in which the mutual writhe (I12) of all (disjoint) sub-chains of a specified length 
is computed and recorded (: written to a file).

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


Input presently only changeable by making a change inside the code block (ie parameters for which 
the value is hard-coded inside the block). Presently this only goes for a set of pre-set strings 
used for automatic generation of output file names:

fileNameOut*: define parts of file names of the files to be written out; see the code 
for how the file names are generated from these parts (or just try it out)

The remaining parameters of the function computGI are all internal, ie their value are set when the code
is run, and so keep values that may change dynamically.


Output:
the function returns an integer (=1 if succesful). The function writes results to .txt files
as set with the inputs: write_all_b, write_final_b, closed_loops_b and writheSubChainPairs_b. 
If desired the function can write excerpts of results to the screen (use print_b or print_basic_b). 
The function also writes certain time consumption values to the screen at the end of the run.

Notes:
This version of the function uses "global memory allocation": a block of memory sufficient for
holding the invariants' values of the longest chain among all in the given directory (of PDB-files) 
is allocated up front. The directory is looped through and, for each structure (chain or, possibly, 
sub-chain of the given PDB-file) the invariants are computed after having being reinitialized to zero 
everywhere relevant (i.e. set to zero in the simplex of size given by the current structure). For each 
structure (chain or, possibly, sub-chain of the given PDB-file) memory is allocated for the chain of 
segments and initialized with values as read in from the PDB-file. This memory is freed again at the 
completition of the computations/write-outs of the structure (chain or sub-chain of PDB-file).

The globally allocated memory (for the invariants' values) is freed after the completion of the last 
structure in the directory. 
*/
int computeGI(char DBName[100], char *dirPath, char *outputPath, char *fileNameChains, int order, int incl_abs_b, int full_b, int split_b, int write_final_b, int write_all_b, 
			  int closed_loops_b, int closedLoopLength, int closedLoopDist, int pokeLength, double thresholdLinks, double thresholdPokes, int invValSubChainPairs_b, int writeSubChainPairsAll_b, int subChainLength, double thresholdWrithe, 
			  int print_b, int print_basic_b){

	int returnval = 0;

	int strLengthCutOff = 10; /*chains shorter than this int will not be included in computation of invariants*/

	char* pokeLengthStr; /*for holding pokeLength converted to a string*/
	char* closedLoopLengthStr; /*for holding closedLoopLength converted to a string*/
	char* subChainLengthStr; /*for holding subChainLength converted to a string*/

	struct dirContent dirContent;
	int numberOfFiles = 0;
	int subStructureCnt = 0; /*some structures are multi-mers; we compute the invariants for each sub-structure (each mono-mer)*/
	const int numberOfPerturbations = 1; /*only there since needed when writing out the final values (preparation for perturbations code)*/
	const int pertNo = 0; /*ditto*/
	char ** ptr_dirList;
	
	/*for "dynamically" generating the file names for output:*/
	char orderStr[10] = "-1"; //value is placeholder

	char fileNameOut1[1000] = "\\Ivalues_order_";
	char fileNameOutAll1[1000] = "\\All_Ivalues_order_";
	char fileNameOutClosedLoops1[1000] = "\\ClosedLoopChars"; 
	char fileNameOutSubChainPairs1[1000] = "\\SubChainPairChars"; 
	char fileNameOut2[1000] = "_minmax_"; //; //; //"_incl_abs_"; //"_"; //

	//char fileNameOut3[1000] = "_mem_alloc_top100.txt"; 
	//char fileNameOutChains1[1000] = "\\chains_closed_loops_top8000.txt";
	char fileNameOutChains1[1000]; 
	char fileNameOut3[1000] = "computeGI_"; 
	char fileNameOut4[1000] = "_loops";
	char fileNameOut5[1000] = "_pokeLength";
	char fileNameOut6[1000] = "_subChainLength";
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

	//clock_t startChain, endChain, startComp, endComp, startComplete, endComplete;
	clock_t startChain, endChain, startComplete, endComplete;
    LARGE_INTEGER startComp,  endComp, ElapsedMicroseconds;
    LARGE_INTEGER Frequency; 
	long double timeChain = 0, compTime = 0, compTimeComplete = 0, timeComplete = 0, max_compTime = 0, max_timeChain = 0;

	FILE *ptr_fileIn;

	struct chainInStructure chainInStr;
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

	/*for closed loops finding:*/
	struct twoSegmentIndex *ptr_closedLoopInd;
	int initValNrOfClosedLoops = 0;
	int nrOfClosedLoops = 0;
	int cntClosedLoop = 0;

	FILE *fp; //dummy file ptr used when checking file names

	/*Declarations: done*/

	/*Generation of file names for output (as desired):*/

	if (write_all_b == 1){

		sprintf(orderStr, "%d", order);
		strcat(fileNameOutAll, outputPath);
		strcat(fileNameOutAll, fileNameOutAll1);
		strcat(fileNameOutAll, orderStr);
		strcat(fileNameOutAll, fileNameOut2);
		strcat(fileNameOutAll, fileNameOut3);
		strcat(fileNameOutAll, DBName);
		strcat(fileNameOutAll, extensionName);
		printf("Results for all vertices will be written to: %s\n", fileNameOutAll);

		//ptr_fileNameOutAll = "C:\\Users\\Christian\\Bioinformatics\\Thesis\\c_code\\results_C\\test\\All_Ivalues_order_X_excl_abs_global_mem_alloc.txt";
		ptr_fileNameOutAll = fileNameOutAll;
		
		/*clear the file by opening it in w-mode and then closing it again*/
		fp= fopen(ptr_fileNameOutAll, "w");
		if(fp == NULL){
			printf("Sorry, file name %s is invalid. Could the path be wrong? A blank space?", fileNameOutAll);
			getchar();
			exit(EXIT_FAILURE);
		}
		else{
		    fclose(fopen(ptr_fileNameOutAll, "w"));
		}

	}


	if (closed_loops_b == 1){

		pokeLengthStr = (char *) malloc(5*sizeof(char)); /*5 should suffice: can hold a poke length of 1000 ...*/
		closedLoopLengthStr = (char *) malloc(5*sizeof(char)); /*5 should suffice: can hold a closed loop length of 1000 ...*/

		sprintf(pokeLengthStr,"%d",pokeLength);
		sprintf(closedLoopLengthStr,"%d",closedLoopLength);

		/*generate file name for writing out closed loops/poke details:*/
		strcat(fileNameOutClosedLoops, outputPath);
		strcat(fileNameOutClosedLoops, fileNameOutClosedLoops1);
		strcat(fileNameOutClosedLoops, fileNameOut2);
		strcat(fileNameOutClosedLoops, fileNameOut3);
		strcat(fileNameOutClosedLoops, DBName);
		strcat(fileNameOutClosedLoops, fileNameOut4);
		strcat(fileNameOutClosedLoops, closedLoopLengthStr);
		strcat(fileNameOutClosedLoops, fileNameOut5);
		strcat(fileNameOutClosedLoops, pokeLengthStr);
		strcat(fileNameOutClosedLoops, extensionName);
		printf("Results for closed loops will be written to: %s\n", fileNameOutClosedLoops);

		//ptr_fileNameOutAll = "C:\\Users\\Christian\\Bioinformatics\\Thesis\\c_code\\results_C\\test\\All_Ivalues_order_X_excl_abs_global_mem_alloc.txt";
		ptr_fileNameOutClosedLoops = fileNameOutClosedLoops;

		/*clear the file by opening it in w-mode and then closing it again*/
		fp = fopen(ptr_fileNameOutClosedLoops, "w");
		if(fp == NULL){
			printf("Sorry, file name %s is invalid. Could the path be wrong? A blank space?", fileNameOutClosedLoops);
			getchar();
			exit(EXIT_FAILURE);
		}
		else{
		    fclose(fopen(ptr_fileNameOutClosedLoops, "w"));
		}


		/*generate file name for writing out chains:*/
		strcpy(fileNameOutChains1,fileNameChains);
		strcat(fileNameOutChains, outputPath);
		strcat(fileNameOutChains, fileNameOutChains1);
		printf("Chains will be written to: %s\n", fileNameOutChains);

		ptr_fileNameOutChains = fileNameOutChains;
	
		/*clear the file by opening it in w-mode and then closing it again*/
		fp = fopen(ptr_fileNameOutChains, "w");
		if(fp == NULL){
			printf("Sorry, file name %s is invalid. Could the path be wrong? A blank space?", fileNameOutChains);
			getchar();
			exit(EXIT_FAILURE);
		}
		else{
		    fclose(fopen(ptr_fileNameOutChains, "w"));
		}

	}

	

	if (invValSubChainPairs_b != 0 ){

		subChainLengthStr = (char *) malloc(5*sizeof(char)); /*5 should suffice: can hold a sub chain length of 1000 ...*/

		sprintf(subChainLengthStr,"%d",subChainLength);

		/*generate file name for writing out writhe results for sub-chain pairs:*/
		strcat(fileNameOutSubChainPairs, outputPath);
		strcat(fileNameOutSubChainPairs, fileNameOutSubChainPairs1);
		strcat(fileNameOutSubChainPairs, fileNameOut2);
		strcat(fileNameOutSubChainPairs, fileNameOut3);
		strcat(fileNameOutSubChainPairs, DBName);
		strcat(fileNameOutSubChainPairs, fileNameOut6);
		strcat(fileNameOutSubChainPairs, subChainLengthStr);
		strcat(fileNameOutSubChainPairs, extensionName);
		printf("Results for sub-chain pairs will be written to: %s\n", fileNameOutSubChainPairs);

		//ptr_fileNameOutAll = "C:\\Users\\Christian\\Bioinformatics\\Thesis\\c_code\\results_C\\test\\All_Ivalues_order_X_excl_abs_global_mem_alloc.txt";
		ptr_fileNameOutSubChainPairs = fileNameOutSubChainPairs;
		
		/*clear the file by opening it in w-mode and then closing it again*/
		fp = fopen(ptr_fileNameOutSubChainPairs, "w");
		if(fp == NULL){
			printf("Sorry, file name %s is invalid. Could the path be wrong? A blank space?", fileNameOutSubChainPairs);
			getchar();
			exit(EXIT_FAILURE);
		}
		else{
		    fclose(fopen(ptr_fileNameOutSubChainPairs, "w"));
		}

		/*generate file name for writing out chains:*/
		strcpy(fileNameOutChains1,fileNameChains);
		strcat(fileNameOutChains, outputPath);
		strcat(fileNameOutChains, fileNameOutChains1);
		printf("Chains will be written to: %s\n", fileNameOutChains);

		ptr_fileNameOutChains = fileNameOutChains;
	
		/*clear the file by opening it in w-mode and then closing it again*/
		fp = fopen(ptr_fileNameOutChains, "w");
		if(fp == NULL){
			printf("Sorry, file name %s is invalid. Could the path be wrong? A blank space?", fileNameOutChains);
			getchar();
			exit(EXIT_FAILURE);
		}
		else{
		    fclose(fopen(ptr_fileNameOutChains, "w"));
		}

	}
	

	//if (write_all_b == 1){

	//	sprintf(orderStr, "%d", order);
	//	strcat(fileNameOutAll, outputPath);
	//	strcat(fileNameOutAll, fileNameOutAll1);
	//	strcat(fileNameOutAll, orderStr);
	//	strcat(fileNameOutAll, fileNameOut2);
	//	strcat(fileNameOutAll, fileNameOut3);
	//	strcat(fileNameOutAll, DBName);
	//	strcat(fileNameOutAll, extensionName);
	//	printf("Results for all vertices will be written to: %s\n", fileNameOutAll);

	//	//ptr_fileNameOutAll = "C:\\Users\\Christian\\Bioinformatics\\Thesis\\c_code\\results_C\\test\\All_Ivalues_order_X_excl_abs_global_mem_alloc.txt";
	//	ptr_fileNameOutAll = fileNameOutAll;
	//	
	//	/*clear the file by opening it in w-mode and then closing it again*/
	//    fclose(fopen(ptr_fileNameOutAll, "w"));

	//}


	//if (closed_loops_b == 1){

	//	pokeLengthStr = (char *) malloc(5*sizeof(char)); /*5 should suffice: can hold a poke length of 1000 ...*/
	//	closedLoopLengthStr = (char *) malloc(5*sizeof(char)); /*5 should suffice: can hold a closed loop length of 1000 ...*/

	//	sprintf(pokeLengthStr,"%d",pokeLength);
	//	sprintf(closedLoopLengthStr,"%d",closedLoopLength);

	//	/*generate file name for writing out closed loops/poke details:*/
	//	strcat(fileNameOutClosedLoops, outputPath);
	//	strcat(fileNameOutClosedLoops, fileNameOutClosedLoops1);
	//	strcat(fileNameOutClosedLoops, fileNameOut2);
	//	strcat(fileNameOutClosedLoops, fileNameOut3);
	//	strcat(fileNameOutClosedLoops, DBName);
	//	strcat(fileNameOutClosedLoops, fileNameOut4);
	//	strcat(fileNameOutClosedLoops, closedLoopLengthStr);
	//	strcat(fileNameOutClosedLoops, fileNameOut5);
	//	strcat(fileNameOutClosedLoops, pokeLengthStr);
	//	strcat(fileNameOutClosedLoops, extensionName);
	//	printf("Results for closed loops will be written to: %s\n", fileNameOutClosedLoops);

	//	//ptr_fileNameOutAll = "C:\\Users\\Christian\\Bioinformatics\\Thesis\\c_code\\results_C\\test\\All_Ivalues_order_X_excl_abs_global_mem_alloc.txt";
	//	ptr_fileNameOutClosedLoops = fileNameOutClosedLoops;
	//	
	//	/*clear the file by opening it in w-mode and then closing it again*/
	//    fclose(fopen(ptr_fileNameOutClosedLoops, "w"));

	//	/*generate file name for writing out chains:*/
	//	strcpy(fileNameOutChains1,fileNameChains);
	//	strcat(fileNameOutChains, outputPath);
	//	strcat(fileNameOutChains, fileNameOutChains1);

	//	ptr_fileNameOutChains = fileNameOutChains;
	//
	//	/*clear the file by opening it in w-mode and then closing it again*/
	//    fclose(fopen(ptr_fileNameOutChains, "w"));

	//}
	//

	//if (invValSubChainPairs_b != 0 ){

	//	subChainLengthStr = (char *) malloc(5*sizeof(char)); /*5 should suffice: can hold a sub chain length of 1000 ...*/

	//	sprintf(subChainLengthStr,"%d",subChainLength);

	//	/*generate file name for writing out writhe results for sub-chain pairs:*/
	//	strcat(fileNameOutSubChainPairs, outputPath);
	//	strcat(fileNameOutSubChainPairs, fileNameOutSubChainPairs1);
	//	strcat(fileNameOutSubChainPairs, fileNameOut2);
	//	strcat(fileNameOutSubChainPairs, fileNameOut3);
	//	strcat(fileNameOutSubChainPairs, DBName);
	//	strcat(fileNameOutSubChainPairs, fileNameOut6);
	//	strcat(fileNameOutSubChainPairs, subChainLengthStr);
	//	strcat(fileNameOutSubChainPairs, extensionName);
	//	printf("Results for sub-chain pairs will be written to: %s\n", fileNameOutSubChainPairs);

	//	//ptr_fileNameOutAll = "C:\\Users\\Christian\\Bioinformatics\\Thesis\\c_code\\results_C\\test\\All_Ivalues_order_X_excl_abs_global_mem_alloc.txt";
	//	ptr_fileNameOutSubChainPairs = fileNameOutSubChainPairs;
	//	
	//	/*clear the file by opening it in w-mode and then closing it again*/
	//    fclose(fopen(ptr_fileNameOutSubChainPairs, "w"));

	//	/*generate file name for writing out chains:*/
	//	strcpy(fileNameOutChains1,fileNameChains);
	//	strcat(fileNameOutChains, outputPath);
	//	strcat(fileNameOutChains, fileNameOutChains1);

	//	ptr_fileNameOutChains = fileNameOutChains;
	//
	//	/*clear the file by opening it in w-mode and then closing it again*/
	//    fclose(fopen(ptr_fileNameOutChains, "w"));

	//}
	
	/*Get content of the desired directory: number of files and a list of the file names:*/
	dirContent = ListDirectoryContents(dirPath);

	numberOfFiles = dirContent.numberOfFiles; 
	ptr_dirList = dirContent.ptr_dirList;

	//printf("pi: %lf", pi);

	printf("Number of files in directory:%d \n",numberOfFiles);

	QueryPerformanceFrequency(&Frequency); //for measuring time consumption in microseconds

	startComplete = clock();

	/*loop through the list of filenames in the directory to find the max chain length (for setting
	the size in the (global) mem alloc to I_measures's ptr's):*/
	for (fileNr = 0; fileNr <numberOfFiles; fileNr ++){

		if (print_b == 1){
			printf("File name:%s fileNr:%d\n",ptr_dirList[fileNr], fileNr);
		}

		ptr_fileIn = fopen(ptr_dirList[fileNr], "r");


		if (!ptr_fileIn){
			printf("Sorry, the file: %s, fileNr: %d could not be read\n", ptr_dirList[fileNr], fileNr);
			continue;
		}

		/*get number of chains, their id's and lengths for this file:*/
		chainInStr = readPDBChainStructure(ptr_fileIn);

		fclose(ptr_fileIn); 

		for (i=0;i <= chainInStr.numberOfChains -1;i++){
			maxChainLen = max(chainInStr.ptr_chainInStr[i].chainLength,maxChainLen);
			if(chainInStr.ptr_chainInStr[i].chainLength >= strLengthCutOff){
				subStructureCnt += 1; /*count and index the sub-structures for which the invariants will be computed*/
			}
		}

	}

	L = maxChainLen -1;

	printf("max chain len:%d\n", maxChainLen);
	printf("sub structure count:%d\n", subStructureCnt);

	/*MEMORY ALLOCATION*/
	/*We allocate to longest chain*/
	returnVal = alloc_init_I_measures(&I_measures, order, full_b, chainNr, maxChainLen, closed_loops_b, &ptr_closedLoopInd);

	//ptr_segment = (struct  segment *) calloc (L, sizeof(struct segment)); /* for the case that re-allocating of mem to ptr_segment is used instead*/

	if (write_final_b == 1){

		/*allocate memory to pointer ptr_I_values and initialize (obs: numberOfPerturbations is preset to 1):*/
		returnVal = alloc_init_I_values(&ptr_I_values, subStructureCnt, numberOfPerturbations);

	}

	subStructureCnt =0; /*reset*/
	

	//startComplete = clock();

	/*Main loop: loop through the list of filenames in the directory and compute the desired measures all along:*/
	for (fileNr = 0; fileNr <numberOfFiles; fileNr ++){

		if (fileNr%100 == 0){printf("Now at file nr:%d\n", fileNr);}

		if (print_basic_b == 1){
			printf("Now handling file:%s\nfileNr:%d\n",ptr_dirList[fileNr], fileNr);
		}


		ptr_fileIn = fopen(ptr_dirList[fileNr], "r");

		I_measures.fileName = ptr_dirList[fileNr];


		if (!ptr_fileIn){
			printf("Sorry, the file could not be found\n");
			getchar();
			return 0;
		}


		/*get number of chains and their lengths for this file:*/
		chainInStr = readPDBChainStructure(ptr_fileIn);

		/*loop through the chains in the current structure/file ... and compute ... :*/
		chainNr = 0; /*reset*/
		for (chainNr = 0; chainNr <= chainInStr.numberOfChains -1 ; chainNr ++){

			startChain = clock();

			I_measures.chainId = chainInStr.ptr_chainInStr[chainNr].chainId;
			I_measures.chainNr = &chainNr;

			chainLen = chainInStr.ptr_chainInStr[chainNr].chainLength;

			if (print_basic_b ==1){
				printf("sub str cnt:%d\n",subStructureCnt);
				printf("Chain id:%s\n", chainInStr.ptr_chainInStr[chainNr].chainId);
				printf("Chain nr:%d\n", chainNr);
				printf("Chain Length: %d\n", chainLen);
			}

			/*if the chain length is less than 6 we skip the computation*/
			if(chainLen < strLengthCutOff){
				if (print_basic_b ==1){
					printf("Chain skipped since length is < %d\n",strLengthCutOff);
				}
				continue;
			}


			I_measures.chainLen = &chainLen;

			rewind(ptr_fileIn);
			ptr_chain = main_readPDB(ptr_fileIn, chainNr, chainLen);

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

			//startComp = clock();
	   	    QueryPerformanceCounter(&startComp); 

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


			/*If desired compute and record the mutual writhe of pairs of closed-loops and 
			pairs (closed loops, segment) at set lengths. Aim is to search for particular geometries 
			(links and "pokes")*/
			if (closed_loops_b == 1){
				returnVal = aggrAndW_wClosedLoops(ptr_segment, chainLen, order, full_b, I_measures, ptr_closedLoopInd, closedLoopLength, closedLoopDist, &nrOfClosedLoops);

				/*examine the closed loops for links and "pokes" and store/write out the results:*/
				returnVal = examineClosedLoops(ptr_fileNameOutClosedLoops, ptr_closedLoopInd, nrOfClosedLoops, I_measures, chainLen, ptr_segment, ptr_fileNameOutChains, ptr_chain, pokeLength, thresholdLinks, thresholdPokes);

			}

			//printf("Nr of closed loops found in this structure: %d\n", nrOfClosedLoops);


			/*If desired compute and record the mutual writhe/invariant value of pairs of sub-chains of 
			a set length. Aim is a more general search for particular geometries (e.g. links)*/
			if (invValSubChainPairs_b != 0){

				if(chainLen > 2*subChainLength + 1){

					if (invValSubChainPairs_b == 1){ /*compute only mutual writhe*/ 
						returnVal = writheSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, subChainLength, ptr_fileNameOutChains, ptr_chain, writeSubChainPairsAll_b, thresholdWrithe);
					} 

					if (invValSubChainPairs_b == 2){ /*compute only mutual invariant values for order 1 and order 2; this takes setting the order parameter to 3 or the full_b parameter to 1*/ 

						if (order > 2 || full_b ==1){

							/*genInvariantSubChainPairs(char *ptr_fileNameOut, struct segment *ptr_segment, int chainLength, struct I_ptr I_measures, double **I_measure, char **I_measure_name, int subChainLength, char *ptr_fileNameOutChains, struct cAlpha *ptr_chain, int writeAll_b)*/

							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures, I_measures.I12, "I12", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I1234_full, "I1234", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I1324_full, "I1324", subChainLength, writeSubChainPairsAll_b);
							returnVal = genInvariantSubChainPairs(ptr_fileNameOutSubChainPairs, ptr_segment, chainLen, I_measures,I_measures.I1423 , "I1423", subChainLength, writeSubChainPairsAll_b);

							/*write out the chain:*/
							returnVal = writeChain(ptr_fileNameOutChains, I_measures.fileName, I_measures.chainId, ptr_chain, *I_measures.chainLen);

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
							returnVal = writeChain(ptr_fileNameOutChains, I_measures.fileName, I_measures.chainId, ptr_chain, *I_measures.chainLen);

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


			QueryPerformanceCounter(&endComp); 
			//endComp = clock();

			//compTime = ((double) (endComp - startComp))/CLOCKS_PER_SEC;
			//compTime = ((long double) (endComp - startComp)/QueryPerformanceFrequency);
			ElapsedMicroseconds.QuadPart = endComp.QuadPart - startComp.QuadPart;
			if(print_b == 1){
				printf("Elapsed microsecs:%d\n", ElapsedMicroseconds.QuadPart);
			}
			ElapsedMicroseconds.QuadPart *= 1000000;
			ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
			compTime = ElapsedMicroseconds.QuadPart;

			//compTime = ((double) (endComp - startComp))/ CLOCKS_PER_SEC;
			//compTime = ((double) (endComp - startComp));
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


			free(ptr_chain);
			free(ptr_segment);	

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

	} /*end fileNr loop (main loop)*/


	/*free the globally allocated memory*/
	free_I_measures(I_measures, order, full_b, maxChainLen);

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
		//ptr_fileNameOut = "C:\\Users\\Christian\\Bioinformatics\\Thesis\\c_code\\results_C\\test\\Ivalues_order_X_excl_abs_global_mem_alloc.txt";

		/*clear the file by opening it in w-mode and then closing it again*/
	    fclose(fopen(ptr_fileNameOut, "w"));

		/*call writing the output (obs: numberOfPerturbations is preset to 1)*/
		returnVal = writeIvaluesToFile(ptr_fileNameOut, ptr_I_values, subStructureCnt, numberOfPerturbations);

		/*free the memory used for the purpose:*/
		free(ptr_I_values);
		 
	}

	endComplete = clock();
	//printf("start clocks:%d, end clocks: %d, clocks_per sec:%d", endComplete, startComplete, CLOCKS_PER_SEC);
	timeComplete = ((double) (endComplete - startComplete)) / CLOCKS_PER_SEC;
	//compTimeComplete = ((double) (compTimeComplete)) / CLOCKS_PER_SEC;
	printf("CPU time spend for all files incl loading data:%lf \n", timeComplete);
	printf("CPU time spend for computations on all files:%lf \n", compTimeComplete/=1000000);
	printf("Max CPU time spend for measures across all files:%lf \n", max_timeChain);

	printf("Done number of files:%d \n", fileNr);
	printf("Done number of chains/sub-structures:%d \n", subStructureCnt);

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