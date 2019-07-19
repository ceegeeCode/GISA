/* **********************************************************************************************************************************
*
*
* ***************************************  Sorting utils for GISA_v* and SAonGISA_v* *******************************************************************************
*
* See main code (GISA_v*) for do and how to cite.
*

Author: Christian Grønbæk
Version: v1, 3rd July 2017; v2 Nov 2017; _v3 Jan 2019

This version: _v4, > 1st Jan, 2019.

***********************************************************************************************/




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

void bubbleSort_binned_I_counts(struct binned_I_counts *, int , int );

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

int wordSearch_I_binned_I_counts(int *, struct binned_I_counts *, int , int, int );


/****************************************************
*****************Coding section *********************
*****************************************************/



/*bubble sort of match candidates, going from left to right;
in first iteration the largest element will be placed at the 
end of the array*/
void bubbleSort(struct matchCandidate *ptr_matchCandidate, int topN){

	int swapped_b = 0; //boolean to check if a swap was done or not
	struct matchCandidate temp;
	int i,j;

	// loop through all positions
   for(i = 0; i < topN-1; i++) { 

      swapped_b = 0;
		
      /*loop from left to right*/
      for(j = 0; j < topN-1-i; j++) {

         // check if next number is lesser than current no
         //   swap the numbers. 
         //  (Bubble up the lowest number, ie to prior index)
			
         if(ptr_matchCandidate[j].distMin > ptr_matchCandidate[j+1].distMin) {
            temp = ptr_matchCandidate[j];
            ptr_matchCandidate[j] = ptr_matchCandidate[j+1];
            ptr_matchCandidate[j+1] = temp;

            swapped_b = 1;
        
         }
      }

      // if no number was swapped that means 
      //   array is sorted now, break the loop. 
      if(swapped_b == 0) {
         break;
      }

   }


}

/*bubble sort of match candidates, going from right to left;
in first iteration the smallest element will be placed at the 
beginning of the array*/
void bubbleSort2(struct matchCandidate *ptr_matchCandidate, int topN){

	int swapped_b = 0; //boolean to check if a swap was done or not
	struct matchCandidate temp;
	int i,j;

	// loop through all positions
   for(i = 0; i < topN-1; i++) { 

      swapped_b = 0;

	  /*loop from right to left*/
      for(j = topN-1; j > i ; j--) {

         // check if next number is lesser than current no
         //   swap the numbers. 
         //  (Bubble up the lowest number, ie to prior index)
			
         if(ptr_matchCandidate[j].distMin < ptr_matchCandidate[j-1].distMin) {
            temp = ptr_matchCandidate[j-1];
            ptr_matchCandidate[j-1] = ptr_matchCandidate[j];
            ptr_matchCandidate[j] = temp;

            swapped_b = 1;
        
         }
		 
      }

	  //printf("distMin at j: %d is %f\n", j, ptr_matchCandidate[j].distMin );
      // if no number was swapped that means 
      //   array is sorted now, break the loop. 
      if(swapped_b == 0) {
		 //printf("At %d iterations no further swaps made\n", i);
         break;
      }

   }

	/*
	for(int x=0; x<n; x++)

	{

		for(int y=0; y<n-1; y++)

		{

			if(array[y]>array[y+1])

			{

				int temp = array[y+1];

				array[y+1] = array[y];

				array[y] = temp;

			}

		}

	}
	*/

}

/*As bubblesort2 but for match candidates using window pairs*/
void bubbleSort2_windowPairs(struct matchCandidate_windowPair *ptr_matchCandidate_windowPair, int topN){

	int swapped_b = 0; //boolean to check if a swap was done or not
	struct matchCandidate_windowPair temp;
	int i,j;

	// loop through all positions
   for(i = 0; i < topN-1; i++) { 

      swapped_b = 0;

	  /*loop from right to left*/
      for(j = topN-1; j > i ; j--) {

         // check if next number is lesser than current no
         //   swap the numbers. 
         //  (Bubble up the lowest number, ie to prior index)
			
         if(ptr_matchCandidate_windowPair[j].distMin < ptr_matchCandidate_windowPair[j-1].distMin) {
            temp = ptr_matchCandidate_windowPair[j-1];
            ptr_matchCandidate_windowPair[j-1] = ptr_matchCandidate_windowPair[j];
            ptr_matchCandidate_windowPair[j] = temp;

            swapped_b = 1;
        
         }
		 
      }

	  //printf("distMin at j: %d is %f\n", j, ptr_matchCandidate_windowPair[j].distMin );
      // if no number was swapped that means 
      //   array is sorted now, break the loop. 
      if(swapped_b == 0) {
		 //printf("At %d iterations no further swaps made\n", i);
         break;
      }

   }


}

/*scores bubble sort, left to rigth*/
void bubbleSortScoresLR(struct matchScore *ptr_matchScore, int nrOfDbRec){

	int swapped_b = 0; //boolean to check if a swap was done or not
	struct matchScore temp;
	int i,j;

	// loop through all positions
   for(i = 0; i < nrOfDbRec-1; i++) { 

      swapped_b = 0;
		
      /*loop from left to right*/
      for(j = 0; j < nrOfDbRec-1 -i  ; j++) {

         // check if next number is lesser than current no
         //   swap the numbers. 
         //  (Bubble up the highest number, ie to prior index)
			
         if(ptr_matchScore[j].score < ptr_matchScore[j+1].score) {
            temp = ptr_matchScore[j];
            ptr_matchScore[j] = ptr_matchScore[j+1];
            ptr_matchScore[j+1] = temp;

            swapped_b = 1;
        
         }
      }

      // if no number was swapped that means 
      //   array is sorted now, break the loop. 
      if(swapped_b == 0) {
         break;
      }

   }

}

void bubbleSortScoresRL(struct matchScore *ptr_matchScore, int nrOfDbRec){

	int swapped_b = 0; //boolean to check if a swap was done or not
	struct matchScore temp;
	int i,j;

	// loop through all positions
   for(i = 0; i < nrOfDbRec-1; i++) { 

      swapped_b = 0;
		
      /*loop from right to left*/
      for(j = nrOfDbRec-1; j > i ; j--) {

         // check if next number is lesser than current no
         //   swap the numbers. 
         //  (Bubble up the highest number, ie to prior index)
			
         if(ptr_matchScore[j].score < ptr_matchScore[j+1].score) {
            temp = ptr_matchScore[j];
            ptr_matchScore[j] = ptr_matchScore[j+1];
            ptr_matchScore[j+1] = temp;

            swapped_b = 1;
        
         }
      }

      // if no number was swapped that means 
      //   array is sorted now, break the loop. 
      if(swapped_b == 0) {
         break;
      }

   }

}

/*DOESN'T WORK! 
modified version of the lexicographic ordering allowing
a pre-set total nr of mismatches; THE INTENTION WAS ONLY TO ALLOW a mismatch 
of 1 per letter, but that results in an "order" which isn't transitive!
The same holds for the present version (which just allows some nr of mismatches) 
--- it's not a proper order.
*/
int lexicographicMod(int *v, int *w, int *nrOfMismatches, int allowedNrOfMismatches, int nrOfEntries){

	int returnVal = 0;
	int nrOfMms = 0;

	nrOfMms = *nrOfMismatches;
	//printf("nrOfMms: %d\n",nrOfMms);

	if(nrOfEntries >= 1){

		//printf("v0: %d w0: %d\n", v[0], w[0]);

		if(nrOfMms == allowedNrOfMismatches){ //proceed as for std lexicograhic

			if(v[0] < w[0]){returnVal = 1;}

			else if(v[0] > w[0]){returnVal = -1;}

			else{//recurse if v[0] == w[0]

				//printf("recurse");
	
				returnVal = lexicographicMod(v+1,w+1, &nrOfMms, allowedNrOfMismatches, nrOfEntries -1);

			}


		}
		else if(nrOfMms < allowedNrOfMismatches){ // allow a mismatch

			if(v[0] != w[0]){nrOfMms +=1;} 

			returnVal = lexicographicMod(v+1,w+1, &nrOfMms, allowedNrOfMismatches, nrOfEntries -1);

		}

	}

	return returnVal;
}


void bubbleSortLexicoBinnedDB_LR(struct binned_I_windows *ptr_binned_I_windows, int nrOfDbrecords, int nrOfEntries, int allowedNrOfMismatches){

	int nrOfMismatches = 0; 

	int swapped_b = 0; //boolean to check if a swap was done or not
	struct binned_I_windows temp;
	int i,j;

	int k, cnt; //for testing

	int lex;

	printf("nrOfDbrecords: %d\n", nrOfDbrecords);

	//getchar();

	// loop through all positions
   for(i = 0; i < nrOfDbrecords-1; i++) { 

	  if(i%100 == 0){ printf("First %d places done\n", i); }

      swapped_b = 0;

	  /*if(i < 10){
		for(k = 0; k < nrOfEntries; k++){ printf("Before. %d 'th rec, bin: %d\n", i, ptr_binned_I_windows[i].binnedIvector[k]);}
	  }*/

	  cnt =0;

      /*loop from left to right*/
      for(j = 0; j < nrOfDbrecords-1 -i  ; j++) {

		  //if(j%1000 ==0){printf("j:%d\n", j);}

         // check if next number is lesser than current no
         //   swap the numbers. 
         //  (Bubble up the lexicographically lowest binned I-vector, ie to prior index)

		 nrOfMismatches = 0; //reset;
		 //lex = lexicographicMod(ptr_binned_I_windows[j].binnedIvector, ptr_binned_I_windows[j+1].binnedIvector, &nrOfMismatches, allowedNrOfMismatches, nrOfEntries);
		 lex = lexicographic(ptr_binned_I_windows[j].binnedIvector, ptr_binned_I_windows[j+1].binnedIvector, nrOfEntries);

		 //printf("lex:%d\n", lex);

         if( lex == -1) {

			/*printf("Before: rec nr %d: %d %d %d %d %d\n", j, ptr_binned_I_windows[j].binnedIvector[0], ptr_binned_I_windows[j].binnedIvector[1], ptr_binned_I_windows[j].binnedIvector[2], ptr_binned_I_windows[j].binnedIvector[3], ptr_binned_I_windows[j].binnedIvector[4]);
			printf("Before: rec nr %d: %d %d %d %d %d\n", j+1, ptr_binned_I_windows[j+1].binnedIvector[0], ptr_binned_I_windows[j+1].binnedIvector[1], ptr_binned_I_windows[j+1].binnedIvector[2], ptr_binned_I_windows[j+1].binnedIvector[3], ptr_binned_I_windows[j+1].binnedIvector[4]);*/


            temp = ptr_binned_I_windows[j];
            ptr_binned_I_windows[j] = ptr_binned_I_windows[j+1];
            ptr_binned_I_windows[j+1] = temp;

            swapped_b = 1;

			cnt += 1;

			/*printf("After: rec nr %d: %d %d %d %d %d\n", j, ptr_binned_I_windows[j].binnedIvector[0], ptr_binned_I_windows[j].binnedIvector[1], ptr_binned_I_windows[j].binnedIvector[2], ptr_binned_I_windows[j].binnedIvector[3], ptr_binned_I_windows[j].binnedIvector[4]);
			printf("After: rec nr %d: %d %d %d %d %d\n", j+1, ptr_binned_I_windows[j+1].binnedIvector[0], ptr_binned_I_windows[j+1].binnedIvector[1], ptr_binned_I_windows[j+1].binnedIvector[2], ptr_binned_I_windows[j+1].binnedIvector[3], ptr_binned_I_windows[j+1].binnedIvector[4]);

			getchar();*/

         }
		 
		 //if(cnt%1000 ==0){printf("cnt:%d\n", cnt);}

      }

	  //printf("cnt:%d\n", cnt);

      // if no number was swapped that means 
      //   array is sorted now, break the loop. 
      if(swapped_b == 0) {
         break;
      }

   }


}

/*Sort binned DB results by lexicographic ordering;
nrOfDBrecords: the number of all windows in the DB (for which there is a set of I-values), summed up over
all structures (file, chain) that is*/
void bubbleSortLexicoBinnedDB_RL(struct binned_I_windows *ptr_binned_I_windows, int nrOfDbrecords, int nrOfEntries, int allowedNrOfMismatches){

	int nrOfMismatches = 0; 

	int swapped_b = 0; //boolean to check if a swap was done or not
	struct binned_I_windows temp;
	int i,j;

	int k, cnt; //for testing

	int lex;

	printf("nrOfDbrecords: %d\n", nrOfDbrecords);

	//getchar();

	// loop through all positions
   for(i = 0; i < nrOfDbrecords-1; i++) { 

	  if(i%100 == 0){ printf("First %d places done\n", i); }

      swapped_b = 0;

	  /*if(i < 10){
		for(k = 0; k < nrOfEntries; k++){ printf("Before. %d 'th rec, bin: %d\n", i, ptr_binned_I_windows[i].binnedIvector[k]);}
	  }*/

	  //getchar();
		

	  cnt =0;

      /*loop from right to left*/
      for(j = nrOfDbrecords-1; j > i ; j--) {

		  //if(j%1000 ==0){printf("j:%d\n", j);}

         // check if next number is lesser than current no
         //   swap the numbers. 
         //  (Bubble up the lexicographicaly lowest binned I-vector, ie to prior index)
		 nrOfMismatches = 0; //reset;
		 //lex = lexicographicMod(ptr_binned_I_windows[j].binnedIvector, ptr_binned_I_windows[j-1].binnedIvector, &nrOfMismatches, allowedNrOfMismatches, nrOfEntries);
		 lex = lexicographic(ptr_binned_I_windows[j].binnedIvector, ptr_binned_I_windows[j-1].binnedIvector, nrOfEntries);

		 //printf("lex:%d\n", lex);

         if( lex== 1) {
            temp = ptr_binned_I_windows[j-1];
            ptr_binned_I_windows[j-1] = ptr_binned_I_windows[j];
            ptr_binned_I_windows[j] = temp;

            swapped_b = 1;

			cnt += 1;
        
         }
		 
		 //if(cnt%1000 ==0){printf("cnt:%d\n", cnt);}

      }

	  /*if(i < 10){
		for(k = 0; k< nrOfEntries; k++){ printf("After. %d 'th rec, bin: %d\n", i, ptr_binned_I_windows[i].binnedIvector[k]);}
	  }*/

	  //getchar();

      // if no number was swapped that means 
      //   array is sorted now, break the loop. 
      if(swapped_b == 0) {
         break;
      }

   }

}

void bubbleSort_binned_I_counts(struct binned_I_counts *ptr_binned_I_counts, int nrOfRecords, int nrOfEntries){

	int swapped_b = 0; //boolean to check if a swap was done or not
	struct binned_I_counts temp;
	int i,j;

	int lex;

	// loop through all positions
   for(i = 0; i < nrOfRecords-1; i++) { 

	  //if(i%100 == 0){ printf("First %d places done\n", i); }

      swapped_b = 0;

	  /*if(i < 10){
		for(k = 0; k < nrOfEntries; k++){ printf("Before. %d 'th rec, bin: %d\n", i, ptr_binned_I_windows[i].binnedIvector[k]);}
	  }*/


      /*loop from left to right*/
      for(j = 0; j < nrOfRecords-1 -i  ; j++) {

         // check if next number is lesser than current
         //   swap the numbers. 
         //  (Bubble up the lexicographically lowest binned I-vector, ie to prior index)

				lex = lexicographic(ptr_binned_I_counts[j].binnedIvector, ptr_binned_I_counts[j+1].binnedIvector, nrOfEntries);

				//printf("lex:%d\n", lex);

         if( lex == -1) {

            temp = ptr_binned_I_counts[j];
            ptr_binned_I_counts[j] = ptr_binned_I_counts[j+1];
            ptr_binned_I_counts[j+1] = temp;

            swapped_b = 1;

         }
		 
      }

      // if no number was swapped that means 
      //   array is sorted now, break the loop. 
      if(swapped_b == 0) {
         break;
      }

   }


}


/*heap sort match candidates (on the dist-value "distMin")*/
void siftDownMatchCandidates(struct matchCandidate *ptr_matchCandidate, int start, int end){

	int root = start;
	int child;

	struct matchCandidate temp;

	while(2*root + 1 <= end){ //while root has a child

		child = 2*root + 1; //"left" child
		
		//printf("in sift down ... child:%d\n", child);

		/*if right child exists and its value is larger than left child's: move to right child instead*/
		if(child + 1 <= end && ptr_matchCandidate[child].distMin < ptr_matchCandidate[child+1].distMin){
			child = child +1;
		}
		/*swap root and child if value at child > value at root*/
		if(ptr_matchCandidate[root].distMin < ptr_matchCandidate[child].distMin){
			
			temp = ptr_matchCandidate[root];
			ptr_matchCandidate[root] =  ptr_matchCandidate[child];
			ptr_matchCandidate[child] = temp;

			root = child;
		}
		else{
			return;
		}

	}

};

void heapifyMatchCandidates(struct matchCandidate *ptr_matchCandidate, int topN){

	int start;

	/*start at last possible root:*/
	if(topN %2 == 0){
		
	start = (int) (topN - 2)/2;
	}
	else{
		start = (int) (topN - 1)/2;
	}

	while(start >= 0){

		//printf("in heapify .. start:%d\n", start);

		/*sift down the element at index start until all nodes
		below that index are in heap order*/
		siftDownMatchCandidates(ptr_matchCandidate, start, topN - 1);

		start = start - 1;
		/*when root has been sifted down all elts are in heap order*/ 
	}
};

void heapSortMatchCandidates(struct matchCandidate *ptr_matchCandidate, int topN){
	
	int end;

	struct matchCandidate temp;

	/*first pt whole array in max-heap order*/
	heapifyMatchCandidates(ptr_matchCandidate, topN);

	end = topN -1;
	while(end > 0){

		//printf("end:%d\n" , end);

		/*swap the present root and the last value of the heap
		(the root is then placed at the end of the heap/array)*/
		temp = ptr_matchCandidate[end];
		ptr_matchCandidate[end] = ptr_matchCandidate[0];
		ptr_matchCandidate[0] = temp;

		/*cut the array size by one so as to leave the final elt (which just came from the root)
		out of next step*/
		end = end -1;

		/*sift down the new root so as to heapify the array terminating at the (new) end*/
		siftDownMatchCandidates(ptr_matchCandidate, 0, end);

	}

};


/*heap sort for binned DB*/
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
void siftDownBinnedDB(struct binned_I_windows *ptr_binned_I_windows, int start, int end, int nrOfEntries){

	int nrOfMismatches = 0; //(int *)malloc(2*sizeof(int));

	int root = start;
	int child;

	struct binned_I_windows temp;

	while(2*root + 1 <= end){ //while root has a child

		child = 2*root + 1; //"left" child
		
		//printf("in sift down ... child:%d\n", child);

		/*if right child exists and its value is larger than left child's: move to right child instead*/
		nrOfMismatches = 0; //reset
		if(child + 1 <= end && lexicographic(ptr_binned_I_windows[child].binnedIvector, ptr_binned_I_windows[child +1].binnedIvector, nrOfEntries) == 1){
			child = child +1;
			nrOfMismatches = 0; //reset
		}
		nrOfMismatches = 0; //reset
		/*swap root and child if value at child > value at root*/
		if(lexicographic(ptr_binned_I_windows[root].binnedIvector, ptr_binned_I_windows[child].binnedIvector, nrOfEntries) == 1){
			
			temp = ptr_binned_I_windows[root];
			ptr_binned_I_windows[root] =  ptr_binned_I_windows[child];
			ptr_binned_I_windows[child] = temp;

			root = child;
			nrOfMismatches = 0; //reset
		}
		else{
			nrOfMismatches = 0; //reset
			return;
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
void heapifyBinnedDB(struct binned_I_windows *ptr_binned_I_windows, int nrOfDBrecords, int nrOfEntries){

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
		siftDownBinnedDB(ptr_binned_I_windows, start, nrOfDBrecords - 1, nrOfEntries);

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
void heapSortBinnedDB(struct binned_I_windows *ptr_binned_I_windows, int nrOfDBrecords, int nrOfEntries){
	
	int end;

	struct binned_I_windows temp;

	/*first pt whole array in max-heap order*/
	heapifyBinnedDB(ptr_binned_I_windows, nrOfDBrecords, nrOfEntries);

	end = nrOfDBrecords -1;
	while(end > 0){

		//printf("end:%d\n" , end);

		/*swap the present root and the last value of the heap
		(the root is then placed at the end of the heap/array)*/
		temp = ptr_binned_I_windows[end];
		ptr_binned_I_windows[end] =  ptr_binned_I_windows[0];
		ptr_binned_I_windows[0] = temp;

		/*cut the array size by one so as to leave the final elt (which just came from the root)
		out of next step*/
		end = end -1;

		/*sift down the new root so as to heapify the array terminating at the (new) end*/
		siftDownBinnedDB(ptr_binned_I_windows, 0, end, nrOfEntries);

	}

}


/*heap sort for binned DB*/
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
void siftDownBinnedDBPairs(struct binned_I_windowPairs *ptr_binned_I_windowPairs, int start, int end, int nrOfEntries){

	int nrOfMismatches = 0; //(int *)malloc(2*sizeof(int));

	int root = start;
	int child;

	struct binned_I_windowPairs temp;

	while(2*root + 1 <= end){ //while root has a child

		child = 2*root + 1; //"left" child
		
		//printf("in sift down ... child:%d\n", child);

		/*if right child exists and its value is larger than left child's: move to right child instead*/
		nrOfMismatches = 0; //reset
		if(child + 1 <= end && lexicographic(ptr_binned_I_windowPairs[child].binnedIvector, ptr_binned_I_windowPairs[child +1].binnedIvector, nrOfEntries) == 1){
			child = child +1;
			nrOfMismatches = 0; //reset
		}
		nrOfMismatches = 0; //reset
		/*swap root and child if value at child > value at root*/
		if(lexicographic(ptr_binned_I_windowPairs[root].binnedIvector, ptr_binned_I_windowPairs[child].binnedIvector, nrOfEntries) == 1){
			
			temp = ptr_binned_I_windowPairs[root];
			ptr_binned_I_windowPairs[root] =  ptr_binned_I_windowPairs[child];
			ptr_binned_I_windowPairs[child] = temp;

			root = child;
			nrOfMismatches = 0; //reset
		}
		else{
			nrOfMismatches = 0; //reset
			return;
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
void heapifyBinnedDBPairs(struct binned_I_windowPairs *ptr_binned_I_windowPairs, int nrOfDBrecords, int nrOfEntries){

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
		siftDownBinnedDBPairs(ptr_binned_I_windowPairs, start, nrOfDBrecords - 1, nrOfEntries);

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
void heapSortBinnedDBPairs(struct binned_I_windowPairs *ptr_binned_I_windowPairs, int nrOfDBrecords, int nrOfEntries){
	
	int end;

	struct binned_I_windowPairs temp;

	/*first pt whole array in max-heap order*/
	heapifyBinnedDBPairs(ptr_binned_I_windowPairs, nrOfDBrecords, nrOfEntries);

	end = nrOfDBrecords -1;
	while(end > 0){

		//printf("end:%d\n" , end);

		/*swap the present root and the last value of the heap
		(the root is then placed at the end of the heap/array)*/
		temp = ptr_binned_I_windowPairs[end];
		ptr_binned_I_windowPairs[end] =  ptr_binned_I_windowPairs[0];
		ptr_binned_I_windowPairs[0] = temp;

		/*cut the array size by one so as to leave the final elt (which just came from the root)
		out of next step*/
		end = end -1;

		/*sift down the new root so as to heapify the array terminating at the (new) end*/
		siftDownBinnedDBPairs(ptr_binned_I_windowPairs, 0, end, nrOfEntries);

	}

}


/*heap sort rarity scores*/
void siftDownRarityScores(struct Irarity_windows *ptr_Irarity_windows, int start, int end){

	int root = start;
	int child;

	struct Irarity_windows temp;

	while(2*root + 1 <= end){ //while root has a child

		child = 2*root + 1; //"left" child
		
		//printf("in sift down ... child:%d\n", child);

		/*if right child exists and its value is smaller than left child's: move to right child instead*/
		if(child + 1 <= end && ptr_Irarity_windows[child].rarityScore > ptr_Irarity_windows[child+1].rarityScore){
			child = child +1;
		}
		/*swap root and child if value at child > value at root*/
		if(ptr_Irarity_windows[root].rarityScore > ptr_Irarity_windows[child].rarityScore){
			
			temp = ptr_Irarity_windows[root];
			ptr_Irarity_windows[root] =  ptr_Irarity_windows[child];
			ptr_Irarity_windows[child] = temp;

			root = child;
		}
		else{
			return;
		}

	}

};

void heapifyRarityScores(struct Irarity_windows *ptr_Irarity_windows, int nrOfWindows){

	int start;

	/*start at last possible root:*/
	if(nrOfWindows%2 == 0){
		
	start = (int) (nrOfWindows - 2)/2;
	}
	else{
		start = (int) (nrOfWindows - 1)/2;
	}

	while(start >= 0){

		//printf("in heapify .. start:%d\n", start);

		/*sift down the element at index start until all nodes
		below that index are in heap order*/
		siftDownRarityScores(ptr_Irarity_windows, start, nrOfWindows - 1);

		start = start - 1;
		/*when root has been sifted down all elts are in heap order*/ 
	}
};

void heapSortRarityScores(struct Irarity_windows *ptr_Irarity_windows, int nrOfWindows){
	
	int end;

	struct Irarity_windows temp;

	/*first pt whole array in max-heap order*/
	heapifyRarityScores(ptr_Irarity_windows, nrOfWindows);

	end = nrOfWindows -1;
	while(end > 0){

		//printf("end:%d\n" , end);

		/*swap the present root and the last value of the heap
		(the root is then placed at the end of the heap/array)*/
		temp = ptr_Irarity_windows[end];
		ptr_Irarity_windows[end] = ptr_Irarity_windows[0];
		ptr_Irarity_windows[0] = temp;

		/*cut the array size by one so as to leave the final elt (which just came from the root)
		out of next step*/
		end = end -1;

		/*sift down the new root so as to heapify the array terminating at the (new) end*/
		siftDownRarityScores(ptr_Irarity_windows, 0, end);

	}

};


/*heap sort rarity scores, for pairs*/
void siftDownRarityScoresPairs(struct Irarity_windowPairs *ptr_Irarity_windowPairs, int start, int end){

	int root = start;
	int child;

	struct Irarity_windowPairs temp;

	while(2*root + 1 <= end){ //while root has a child

		child = 2*root + 1; //"left" child
		
		//printf("in sift down ... child:%d\n", child);

		/*if right child exists and its value is smaller than left child's: move to right child instead*/
		if(child + 1 <= end && ptr_Irarity_windowPairs[child].rarityScore > ptr_Irarity_windowPairs[child+1].rarityScore){
			child = child +1;
		}
		/*swap root and child if value at child > value at root*/
		if(ptr_Irarity_windowPairs[root].rarityScore > ptr_Irarity_windowPairs[child].rarityScore){
			
			temp = ptr_Irarity_windowPairs[root];
			ptr_Irarity_windowPairs[root] =  ptr_Irarity_windowPairs[child];
			ptr_Irarity_windowPairs[child] = temp;

			root = child;
		}
		else{
			return;
		}

	}

};

void heapifyRarityScoresPairs(struct Irarity_windowPairs *ptr_Irarity_windowPairs, int nrOfWindows){

	int start;

	/*start at last possible root:*/
	if(nrOfWindows%2 == 0){
		
	start = (int) (nrOfWindows - 2)/2;
	}
	else{
		start = (int) (nrOfWindows - 1)/2;
	}

	while(start >= 0){

		//printf("in heapify .. start:%d\n", start);

		/*sift down the element at index start until all nodes
		below that index are in heap order*/
		siftDownRarityScoresPairs(ptr_Irarity_windowPairs, start, nrOfWindows - 1);

		start = start - 1;
		/*when root has been sifted down all elts are in heap order*/ 
	}
};

void heapSortRarityScoresPairs(struct Irarity_windowPairs *ptr_Irarity_windowPairs, int nrOfWindows){
	
	int end;

	struct Irarity_windowPairs temp;

	/*first pt whole array in max-heap order*/
	heapifyRarityScoresPairs(ptr_Irarity_windowPairs, nrOfWindows);

	end = nrOfWindows -1;
	while(end > 0){

		//printf("end:%d\n" , end);

		/*swap the present root and the last value of the heap
		(the root is then placed at the end of the heap/array)*/
		temp = ptr_Irarity_windowPairs[end];
		ptr_Irarity_windowPairs[end] = ptr_Irarity_windowPairs[0];
		ptr_Irarity_windowPairs[0] = temp;

		/*cut the array size by one so as to leave the final elt (which just came from the root)
		out of next step*/
		end = end -1;

		/*sift down the new root so as to heapify the array terminating at the (new) end*/
		siftDownRarityScoresPairs(ptr_Irarity_windowPairs, 0, end);

	}

};


/*heap sort array of raw scores (over set of queries)*/
void siftDownRawScores(struct queryRawScore *ptr_queryRawScore, int start, int end){

	int root = start;
	int child;

	struct queryRawScore temp;

	while(2*root + 1 <= end){ //while root has a child

		child = 2*root + 1; //"left" child
		
		//printf("in sift down ... child:%d\n", child);

		/*if right child exists and its value is smaller than left child's: move to right child instead*/
		if(child + 1 <= end && ptr_queryRawScore[child].score > ptr_queryRawScore[child+1].score){
			child = child +1;
		}
		/*swap root and child if value at child > value at root*/
		if(ptr_queryRawScore[root].score > ptr_queryRawScore[child].score){
			
			temp = ptr_queryRawScore[root];
			ptr_queryRawScore[root] =  ptr_queryRawScore[child];
			ptr_queryRawScore[child] = temp;

			root = child;
		}
		else{
			return;
		}

	}

};

void heapifyRawScores(struct queryRawScore *ptr_queryRawScore, int nrOfStructures){

	int start;

	/*start at last possible root:*/
	if(nrOfStructures%2 == 0){
		
	start = (int) (nrOfStructures - 2)/2;
	}
	else{
		start = (int) (nrOfStructures - 1)/2;
	}

	while(start >= 0){

		//printf("in heapify .. start:%d\n", start);

		/*sift down the element at index start until all nodes
		below that index are in heap order*/
		siftDownRawScores(ptr_queryRawScore, start, nrOfStructures - 1);

		start = start - 1;
		/*when root has been sifted down all elts are in heap order*/ 
	}
};

void heapSortRawScores(struct queryRawScore *ptr_queryRawScore, int nrOfStructures){
	
	int end;

	struct queryRawScore temp;

	/*first pt whole array in max-heap order*/
	heapifyRawScores(ptr_queryRawScore, nrOfStructures);

	end = nrOfStructures -1;
	while(end > 0){

		//printf("end:%d\n" , end);

		/*swap the present root and the last value of the heap
		(the root is then placed at the end of the heap/array)*/
		temp = ptr_queryRawScore[end];
		ptr_queryRawScore[end] = ptr_queryRawScore[0];
		ptr_queryRawScore[0] = temp;

		/*cut the array size by one so as to leave the final elt (which just came from the root)
		out of next step*/
		end = end -1;

		/*sift down the new root so as to heapify the array terminating at the (new) end*/
		siftDownRawScores(ptr_queryRawScore, 0, end);

	}

};


/*heap sort for match set*/
void siftDownMatchSet(struct binned_I_windows *ptr_binned_I_windows, int start, int end, int nrOfEntries){

	int root = start;
	int child;

	struct binned_I_windows temp;

	while(2*root + 1 <= end){ //while root has a child

		child = 2*root + 1; //"left" child
		
		//printf("in sift down ... child:%d\n", child);

		/*if right child exists and its value is larger than left child's: move to right child instead*/
		if(child + 1 <= end && ptr_binned_I_windows[child].fileNr < ptr_binned_I_windows[child +1].fileNr){
			child = child +1;
		}
		/*swap root and child if value at child > value at root*/
		if(ptr_binned_I_windows[root].fileNr < ptr_binned_I_windows[child].fileNr){
			
			temp = ptr_binned_I_windows[root];
			ptr_binned_I_windows[root] =  ptr_binned_I_windows[child];
			ptr_binned_I_windows[child] = temp;

			root = child;
		}
		else{
			return;
		}

	}

}

void heapifyMatchSet(struct binned_I_windows *ptr_binned_I_windows, int nrOfDBrecords, int nrOfEntries){

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
		siftDownMatchSet(ptr_binned_I_windows, start, nrOfDBrecords - 1, nrOfEntries);

		start = start - 1;
		/*when root has been sifted down all elts are in heap order*/ 
	}
}

void heapSortMatchSet(struct binned_I_windows *ptr_binned_I_windows, int nrOfDBrecords, int nrOfEntries){
	
	int end;

	struct binned_I_windows temp;

	/*first pt whole array in max-heap order*/
	heapifyMatchSet(ptr_binned_I_windows, nrOfDBrecords, nrOfEntries);

	end = nrOfDBrecords -1;
	while(end > 0){

		//printf("end:%d\n" , end);

		/*swap the present root and the last value of the heap
		(the root is then placed at the end of the heap/array)*/
		temp = ptr_binned_I_windows[end];
		ptr_binned_I_windows[end] =  ptr_binned_I_windows[0];
		ptr_binned_I_windows[0] = temp;

		/*cut the array size by one so as to leave the final elt (which just came from the root)
		out of next step*/
		end = end -1;

		/*sift down the new root so as to heapify the array terminating at the (new) end*/
		siftDownMatchSet(ptr_binned_I_windows, 0, end, nrOfEntries);

	}

}


/*heap sort for match-scores array (set)*/
/*Input: if on_scores_b = 1 the sorting will be on the scores; else it will be on the structure name*/
void siftDownMatchScores(struct matchScore *ptr_matchScore, int start, int end, int nrOfEntries, int on_scores_b){

	int root = start;
	int child;

	struct matchScore temp;

	if(on_scores_b == 1){

		while(2*root + 1 <= end){ //while root has a child

			child = 2*root + 1; //"left" child
		
			//printf("in sift down ... child:%d\n", child);

			/*if right child exists and its value is larger than left child's: move to right child instead*/
			if(child + 1 <= end && ptr_matchScore[child].score > ptr_matchScore[child+1].score){
				child = child +1;
			}
			/*swap root and child if value at child > value at root*/
			if(ptr_matchScore[root].score > ptr_matchScore[child].score){
			
				temp = ptr_matchScore[root];
				ptr_matchScore[root] =  ptr_matchScore[child];
				ptr_matchScore[child] = temp;

				root = child;
			}
			else{
				return;
			}

		}

	} 
	
	else{ //on_scores_b != 1
		
		while(2*root + 1 <= end){ //while root has a child

			child = 2*root + 1; //"left" child
		
			//printf("in sift down ... child:%d\n", child);

			/*if right child exists and its value is larger than left child's: move to right child instead*/
			if(child + 1 <= end && lexicographicChar(ptr_matchScore[child].matchStructureName, ptr_matchScore[child+1].matchStructureName, nrOfEntries) == 1){
				child = child +1;
			}
			/*swap root and child if value at child > value at root*/
			if(lexicographicChar(ptr_matchScore[root].matchStructureName, ptr_matchScore[child].matchStructureName, nrOfEntries) == 1){
			
				temp = ptr_matchScore[root];
				ptr_matchScore[root] =  ptr_matchScore[child];
				ptr_matchScore[child] = temp;

				root = child;
			}
			else{
				return;
			}

		}
	} 


}

void heapifyMatchScores(struct matchScore *ptr_matchScore, int nrOfDBrecords, int nrOfEntries, int on_scores_b){

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
		siftDownMatchScores(ptr_matchScore, start, nrOfDBrecords - 1, nrOfEntries, on_scores_b);

		start = start - 1;
		/*when root has been sifted down all elts are in heap order*/ 
	}
}

void heapSortMatchScores(struct matchScore *ptr_matchScore, int nrOfDBrecords, int nrOfEntries, int on_scores_b){
	
	int end;

	struct matchScore temp;

	/*first pt whole array in max-heap order*/
	heapifyMatchScores(ptr_matchScore, nrOfDBrecords, nrOfEntries, on_scores_b);

	end = nrOfDBrecords -1;
	while(end > 0){

		//printf("end:%d\n" , end);

		/*swap the present root and the last value of the heap
		(the root is then placed at the end of the heap/array)*/
		temp = ptr_matchScore[end];
		ptr_matchScore[end] =  ptr_matchScore[0];
		ptr_matchScore[0] = temp;

		/*cut the array size by one so as to leave the final elt (which just came from the root)
		out of next step*/
		end = end -1;

		/*sift down the new root so as to heapify the array terminating at the (new) end*/
		siftDownMatchScores(ptr_matchScore, 0, end, nrOfEntries, on_scores_b);

	}

}

/*Search through the (lexicographically pre-sorted) match scores array (for a given structure name) 
returns: index in match scores array of hit*/
int bisectionalSearchMatchScores(char *structureName, struct matchScore *ptr_matchScore, int nrOfDbrecords, int nrOfEntries){

	int returnVal = 0;

	int left = 0;
	int right = nrOfDbrecords -1;
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
		lex = lexicographicChar(structureName, ptr_matchScore[mid].matchStructureName, nrOfEntries);

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
		if(lexicographicChar(structureName, ptr_matchScore[left].matchStructureName, nrOfEntries) == 0){hit = left;}
	}
	else if(lex == -1){
		if(lexicographicChar(structureName, ptr_matchScore[right].matchStructureName, nrOfEntries) == 0){hit = right;}
	}

	return hit;

}


/*Various bi-sectional searches:*/

/*Search through the (lexicographically pre-sorted) DB of binned I-vectors for a given
query (an I-window vector repr as binned I-vector)
returns: first and last index in sorted list (ptr_binned_I_windows_DB) matching the query*/
int bisectionalSearch(int *hitLR, struct binned_I_windows binned_I_windows_query, struct binned_I_windows *ptr_binned_I_windows_DB, int nrOfDbrecords, int nrOfEntries, int allowedNrOfMismatches){

	int returnVal = 0;

	int nrOfMismatches = 0; //(int *)malloc(2*sizeof(int));

	int iter = 0;

	int left = 0;
	int right = nrOfDbrecords -1;
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
		nrOfMismatches = 0; //reset
		lex = lexicographicMod(binned_I_windows_query.binnedIvector, ptr_binned_I_windows_DB[mid].binnedIvector, &nrOfMismatches, allowedNrOfMismatches, nrOfEntries);

		if( lex == 1){

			hit = left;


		}
		else if( lex == -1){

			hit = right;

		}
		else{ //lex == 0 so we have a hit "already" 
			 hit = mid;
			 //printf("break\n");
			 break;
		}

		//midNew = (int) ceil( (double) (right + left)/2 );

	}

	if(lex == 1){
		if(lexicographicMod(binned_I_windows_query.binnedIvector, ptr_binned_I_windows_DB[left].binnedIvector, &nrOfMismatches, allowedNrOfMismatches, nrOfEntries) == 0){hit = left;}
	}
	else if(lex == -1){
		if(lexicographicMod(binned_I_windows_query.binnedIvector, ptr_binned_I_windows_DB[right].binnedIvector, &nrOfMismatches, allowedNrOfMismatches, nrOfEntries) == 0){hit = right;}
	}


	/*extend the left and right from the hit if possible*/
	m = hit;
	//printf("hit:%d\n", hit);
	//printf("DB: %d %d %d %d %d\n", ptr_binned_I_windows_DB[hit].binnedIvector[0], ptr_binned_I_windows_DB[hit].binnedIvector[1], ptr_binned_I_windows_DB[hit].binnedIvector[2], ptr_binned_I_windows_DB[hit].binnedIvector[3], ptr_binned_I_windows_DB[hit].binnedIvector[4]);

	hitLR[0] = m;
	while(m >= 0 && lexicographicMod(binned_I_windows_query.binnedIvector, ptr_binned_I_windows_DB[m].binnedIvector, &nrOfMismatches, allowedNrOfMismatches, nrOfEntries) == 0) { 
		//printf("extend left, m:%d\n", m);
		hitLR[0] = m;
		nrOfMismatches = 0; //reset
		m -=1; 
	}

	n = hit;
	hitLR[1] = n;
	while(n < nrOfDbrecords && lexicographicMod(binned_I_windows_query.binnedIvector, ptr_binned_I_windows_DB[n].binnedIvector, &nrOfMismatches, allowedNrOfMismatches, nrOfEntries) == 0) { 
		//printf("extend right, n:%d\n", n);
		hitLR[1] = n;
		nrOfMismatches = 0; //reset
		n +=1; 
	}

	/*hitLR[0] = m+1;
	hitLR[1] = n-1;*/


	return returnVal;

}

/*older version*/
int bisectionalSearch_1(int *hitLR, struct binned_I_windows binned_I_windows_query, struct binned_I_windows *ptr_binned_I_windows_DB, int nrOfDbrecords, int nrOfEntries, int allowedNrOfMismatches){

	int returnVal = 0;

	int nrOfMismatches = 0; //(int *)malloc(2*sizeof(int));

	int iter = 0;

	int left = 0;
	int right = nrOfDbrecords -1;
	int mid = 0;
	int midNew = -1;

	int lex = -2;

	int hit; 

	int m, n;

	while(right -left > 1){

		mid = (int) ceil( (double) (right + left)/2 );

		//printf("mid:%d\n", mid);
		nrOfMismatches = 0; //reset
		lex = lexicographicMod(binned_I_windows_query.binnedIvector, ptr_binned_I_windows_DB[mid].binnedIvector, &nrOfMismatches, allowedNrOfMismatches, nrOfEntries);

		if( lex == 1){

			hit = left;

			right = mid;


		}
		else if( lex == -1){

			hit = right;

			left = mid;
		}
		else{ //lex == 0 so we have a hit "already" 
			 hit = mid;
			 //printf("break\n");
			 break;
		}

		//midNew = (int) ceil( (double) (right + left)/2 );

	}

	/*extend the left and right from the hit if possible*/
	m = hit;
	//printf("hit:%d\n", hit);
	//printf("DB: %d %d %d %d %d\n", ptr_binned_I_windows_DB[hit].binnedIvector[0], ptr_binned_I_windows_DB[hit].binnedIvector[1], ptr_binned_I_windows_DB[hit].binnedIvector[2], ptr_binned_I_windows_DB[hit].binnedIvector[3], ptr_binned_I_windows_DB[hit].binnedIvector[4]);

	hitLR[0] = m;
	while(m >= 0 && lexicographicMod(binned_I_windows_query.binnedIvector, ptr_binned_I_windows_DB[m].binnedIvector, &nrOfMismatches, allowedNrOfMismatches, nrOfEntries) == 0) { 
		//printf("extend left, m:%d\n", m);
		hitLR[0] = m;
		nrOfMismatches = 0; //reset
		m -=1; 
	}

	n = hit;
	hitLR[1] = n;
	while(n < nrOfDbrecords && lexicographicMod(binned_I_windows_query.binnedIvector, ptr_binned_I_windows_DB[n].binnedIvector, &nrOfMismatches, allowedNrOfMismatches, nrOfEntries) == 0) { 
		//printf("extend right, n:%d\n", n);
		hitLR[1] = n;
		nrOfMismatches = 0; //reset
		n +=1; 
	}

	/*hitLR[0] = m+1;
	hitLR[1] = n-1;*/


	return returnVal;

}

/*version of bisectional search which finds the range of hits directly, i.e. without out 
the final left- and right extensions.
29/3 2017: changed by adding "&& mid != lastRight" in first if-clause */
int bisectionalSearch2(int *hitLR, struct binned_I_windows binned_I_windows_query, struct binned_I_windows *ptr_binned_I_windows_DB, int nrOfDbrecords, int nrOfEntries, int allowedNrOfMismatches){

	int returnVal = 0;

	int nrOfMismatches =  0; //(int *)malloc(2*sizeof(int));

	int iter = 0;

	int left = 0;
	int right = nrOfDbrecords -1;
	int lastLeft = left;
	int lastRight = right;
	int mid = 0;
	int midNew = -1;

	int lex = -2;

	int hit; 

	while(right -left > 1 && mid != lastRight){

		lastLeft = left;
		lastRight = right;

		mid = (int) ceil( (double) (right + left)/2 );

		//printf("mid:%d\n", mid);
		nrOfMismatches = 0; //reset
		lex = lexicographicMod(binned_I_windows_query.binnedIvector, ptr_binned_I_windows_DB[mid].binnedIvector, &nrOfMismatches, allowedNrOfMismatches, nrOfEntries);

		if( lex == 1){

			hit = left;

			right = mid;


		}
		else if( lex == -1){

			hit = right;

			left = mid;
		}
		else{ //lex == 0 so we have a hit "already" 
			 hit = mid;
			 //printf("break\n");
			 break;
		}

		//midNew = (int) ceil( (double) (right + left)/2 );

	}

	/*extend the left and right from the hit if possible.
	This is here done in two further bisections:*/
	//printf("hit:%d\n", hit);
	//printf("DB: %d %d %d %d %d\n", ptr_binned_I_windows_DB[hit].binnedIvector[0], ptr_binned_I_windows_DB[hit].binnedIvector[1], ptr_binned_I_windows_DB[hit].binnedIvector[2], ptr_binned_I_windows_DB[hit].binnedIvector[3], ptr_binned_I_windows_DB[hit].binnedIvector[4]);

	/*left-ward extension:*/
	hitLR[0] = hit;
	left = lastLeft;
	right = hit;
	while(right -left > 1) {

		mid = (int) ceil( (double) (right + left)/2 );

		//printf("mid:%d\n", mid);
		nrOfMismatches = 0; //reset
		lex = lexicographicMod(binned_I_windows_query.binnedIvector, ptr_binned_I_windows_DB[mid].binnedIvector, &nrOfMismatches, allowedNrOfMismatches, nrOfEntries);

		
		if( lex == 0){

			right = mid;

		}
		else{

			left = mid;

		}

		hit = right;

		//printf("extend right, n:%d\n", n);
		hitLR[0] = hit;

	}

	/*right-ward extension:*/
	hitLR[1] = hit;
	left = hit;
	right = lastRight;
	while(right -left > 1) {

		mid = (int) ceil( (double) (right + left)/2 );

		//printf("mid:%d\n", mid);
		nrOfMismatches = 0; //reset
		lex = lexicographicMod(binned_I_windows_query.binnedIvector, ptr_binned_I_windows_DB[mid].binnedIvector, &nrOfMismatches, allowedNrOfMismatches, nrOfEntries);

		
		if( lex == 0){

			left = mid;

		}
		else{

			right = mid;

		}

		hit = left;

		//printf("extend right, n:%d\n", n);
		hitLR[1] = hit;

	}

	return returnVal;

}

/*bisectional search for "letter" at position position of query, binned_I_windows_query, among the letters at the same position
in the range of the array ptr_binned_I_windows_DB from  ptr_range.start to ptr_range.end.
Returns: index of the hit; updates ptr_range.start and .end with the last left resp. right interval bounds of the search*/
int bisectSingle(struct binned_I_windows binned_I_windows_query, int position, int *ptr_range, struct binned_I_windows *ptr_binned_I_windows_DB){

	int returnVal = 0;

	int iter = 0;

	int left = ptr_range[0];
	int right = ptr_range[1];

	//int lastLeft = left;
	int lastRight = right;
	int mid = -117;

	int hit = -1; 

	/*if the search range dow not contain the value, binned_I_windows_query.binnedIvector[position], there is no
	hit to be found, and hit = -1 is returned */
	if(binned_I_windows_query.binnedIvector[position] < ptr_binned_I_windows_DB[left].binnedIvector[position] ||
		binned_I_windows_query.binnedIvector[position] > ptr_binned_I_windows_DB[right].binnedIvector[position] ) {

			return hit; 

	}	/*If the range does not contain the value binned_I_windows_query.binnedIvector[position] there is no need to carry on*/


	/*The following if-clause seems superfluous by the preceding if-clause and since the while-loop to follow
	is now entered when left = right; however the hit value must be updated from -1 to the index (=left)*/
	//If the range is actually a single value and we have a hit we return that// 
	if(left == right){
		if(binned_I_windows_query.binnedIvector[position] == ptr_binned_I_windows_DB[left].binnedIvector[position]){
			hit = left;
		}
		//else: hit =-1 by the set initial value

		return hit;

	}
	
	//.. And else we search://
	while(mid != lastRight && left < right){//while(right -left > 1){ //while(mid  != lastRight && left <= right)

		//lastLeft = left;
		lastRight = right;

		mid = (int) ceil( (double) (right + left)/2 );

		//printf("mid:%d\n", mid);

		if( binned_I_windows_query.binnedIvector[position] < ptr_binned_I_windows_DB[mid].binnedIvector[position]){

			right = mid;


		}
		else if( binned_I_windows_query.binnedIvector[position] > ptr_binned_I_windows_DB[mid].binnedIvector[position] ){

			left = mid;
		}
		else{ //lex == 0 so we have a hit
				hit = mid;
				//printf("break\n");
				break; //break while-loop
		}

		/*there remains the possibility that the hit sits at the left end point (added 7th March '18)*/
		if (binned_I_windows_query.binnedIvector[position] == ptr_binned_I_windows_DB[left].binnedIvector[position]) {
			hit = left;
		}

	}

	ptr_range[0] = left;
	ptr_range[1] = right;

	return hit;
}


/*before 7th March '18 a hit != -1 would be returned even if there was no hit:*/
int bisectSingle_old_3(struct binned_I_windows binned_I_windows_query, int position, int *ptr_range, struct binned_I_windows *ptr_binned_I_windows_DB) {

	int returnVal = 0;

	int iter = 0;

	int left = ptr_range[0];
	int right = ptr_range[1];

	//int lastLeft = left;
	int lastRight = right;
	int mid = -117;

	int hit = -1;

	/*if the search range dow not contain the value, binned_I_windows_query.binnedIvector[position], there is no
	hit to be found, and hit = -1 is returned */
	if (binned_I_windows_query.binnedIvector[position] < ptr_binned_I_windows_DB[left].binnedIvector[position] ||
		binned_I_windows_query.binnedIvector[position] > ptr_binned_I_windows_DB[right].binnedIvector[position]) {

		return hit;

	}	/*If the range does not contain the value binned_I_windows_query.binnedIvector[position] there is no need to carry on*/


		/*The following if-clause seems superfluous by the preceding if-clause and since the while-loop to follow
		is now entered when left = right; however the hit value must be updated from -1 to the index (=left)*/
		//If the range is actually a single value and we have a hit we return that// 
	if (left == right) {
		if (binned_I_windows_query.binnedIvector[position] == ptr_binned_I_windows_DB[left].binnedIvector[position]) {
			hit = left;
		}
		//else: hit =-1 by the set initial value

		return hit;

	}

	//.. And else we search://
	while (mid != lastRight && left < right) {//while(right -left > 1){ //while(mid  != lastRight && left <= right)

											  //lastLeft = left;
		lastRight = right;

		mid = (int)ceil((double)(right + left) / 2);

		//printf("mid:%d\n", mid);

		if (binned_I_windows_query.binnedIvector[position] < ptr_binned_I_windows_DB[mid].binnedIvector[position]) {

			hit = left;

			right = mid;


		}
		else if (binned_I_windows_query.binnedIvector[position] > ptr_binned_I_windows_DB[mid].binnedIvector[position]) {

			hit = right;

			left = mid;
		}
		else { //lex == 0 so we have a hit "already" 
			hit = mid;
			//printf("break\n");
			break; //break while-loop
		}

		//midNew = (int) ceil( (double) (right + left)/2 );

	}

	ptr_range[0] = left;
	ptr_range[1] = right;

	return hit;
}

/*Before Jan 22, '18 (new version includes a trailing if-clause!)*/
int bisectSingle_old_2(struct binned_I_windows binned_I_windows_query, int position, int *ptr_range, struct binned_I_windows *ptr_binned_I_windows_DB){

	int returnVal = 0;

	int iter = 0;

	int left = ptr_range[0];
	int right = ptr_range[1];

	int lastLeft = left;
	int lastRight = right;
	int mid = -117;

	int hit = -1; 

	//If the range is actually a single value and we have a hit we return that// 
	if(left == right){
		if(binned_I_windows_query.binnedIvector[position] == ptr_binned_I_windows_DB[left].binnedIvector[position]){
			hit = left;
		}
		//else: hit =-1 by the set initial value

		return hit;

	}
	
	//.. And else we search://
	while(mid != lastRight && left < right){//while(right -left > 1){ //while(mid  != lastRight && left <= right)

		lastLeft = left;
		lastRight = right;

		mid = (int) ceil( (double) (right + left)/2 );

		//printf("mid:%d\n", mid);

		if( binned_I_windows_query.binnedIvector[position] < ptr_binned_I_windows_DB[mid].binnedIvector[position]){

			hit = left;

			right = mid;


		}
		else if( binned_I_windows_query.binnedIvector[position] > ptr_binned_I_windows_DB[mid].binnedIvector[position] ){

			hit = right;

			left = mid;
		}
		else{ //lex == 0 so we have a hit "already" 
				hit = mid;
				//printf("break\n");
				break; //break while-loop
		}

		//midNew = (int) ceil( (double) (right + left)/2 );

	}

	ptr_range[0] = left;
	ptr_range[1] = right;

	return hit;
}

/*Before april 4, '17:*/
int bisectSingle_old(struct binned_I_windows binned_I_windows_query, int position, int *ptr_range, struct binned_I_windows *ptr_binned_I_windows_DB){

	int returnVal = 0;

	int iter = 0;

	int left = ptr_range[0];
	int right = ptr_range[1];

	int lastLeft = left;
	int lastRight = right;
	int mid = -1;

	int hit = -1; 
	
	while(mid != lastRight && left <= right){//while(right -left > 1){ //while(mid  != lastRight && left <= right)

		lastLeft = left;
		lastRight = right;

		mid = (int) ceil( (double) (right + left)/2 );

		printf("mid:%d\n", mid);

		if( binned_I_windows_query.binnedIvector[position] < ptr_binned_I_windows_DB[mid].binnedIvector[position]){

			hit = left;

			right = mid;


		}
		else if( binned_I_windows_query.binnedIvector[position] > ptr_binned_I_windows_DB[mid].binnedIvector[position] ){

			hit = right;

			left = mid;
		}
		else{ //lex == 0 so we have a hit "already" 
				hit = mid;
				//printf("break\n");
				break;
		}

		//midNew = (int) ceil( (double) (right + left)/2 );

	}

	ptr_range[0] = lastLeft;
	ptr_range[1] = lastRight;

	return hit;
}


/*for extending left or right from a hit (with index hitIndex); searches for extension left from ptr_range.start and up to 
hitIndex, and right from hitIndex and to ptr_range.end. In particular hitIndex must sit in the interval [ptr_range.start,
ptr_range.end] for this to work.
Input: leftRight_b -- if 0 both extensions are done; if -1 only left, if 1 only right.
Transforms: updates ptr_range .start and .end to left resp right indexes found.*/
int bisectSingleExtend(int hitIndex, struct binned_I_windows binned_I_windows_query, int position, int *ptr_range, struct binned_I_windows *ptr_binned_I_windows_DB, int leftRight_b, int discrepancy){

	int returnVal = 0;

	int left;
	int lastLeft;
	int right;
	int lastRight;
	int mid = -1;
	//int hit = -1; ch'ed 7 March '18

	int hitLR[2];

	hitLR[0] = hitIndex;
	hitLR[1] = hitIndex;
	
	/*only extend if we have a hit*/
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
		
				if( binned_I_windows_query.binnedIvector[position] == ptr_binned_I_windows_DB[mid].binnedIvector[position] - discrepancy){

					right = mid;

				}
				else{

					left = mid;

				}

				//hit = right; ch'ed 7 March '18

				//printf("extend right, n:%d\n", n);
				hitLR[0] = right; // hit; ch'ed 7 March '18

			}

		}


		if(leftRight_b >= 0){
			/*right-ward extension:*/
			left = hitIndex;
			right = ptr_range[1];
			lastRight = right;
			mid = -1;
			while(lastRight - mid > 0) {

				lastRight =  right;
				mid = (int) ceil( (double) (right + left)/2 );

				//printf("mid:%d\n", mid);
		
				if( binned_I_windows_query.binnedIvector[position] == ptr_binned_I_windows_DB[mid].binnedIvector[position] + discrepancy){

					left = mid;

				}
				else{

					right = mid;

				}

				//hit = left; ch'ed 7 March '18

				//printf("extend right, n:%d\n", n);
				hitLR[1] = left; // hit; ch'ed 7 March '18

			}
		}
	}
	
	ptr_range[0] = hitLR[0];
	ptr_range[1] = hitLR[1];

	return returnVal;

}


/*Before april 4, '17:*/
int bisectSingleExtend_old(int hitIndex, struct binned_I_windows binned_I_windows_query, int position, int *ptr_range, struct binned_I_windows *ptr_binned_I_windows_DB, int leftRight_b){

	int returnVal = 0;

	int left;
	int lastLeft;
	int right;
	int lastRight;
	int mid = -1;
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
		
				if( binned_I_windows_query.binnedIvector[position] == ptr_binned_I_windows_DB[mid].binnedIvector[position]){

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
			while(lastRight - mid > 0) {

				lastRight =  right;
				mid = (int) ceil( (double) (right + left)/2 );

				//printf("mid:%d\n", mid);
		
				if( binned_I_windows_query.binnedIvector[position] == ptr_binned_I_windows_DB[mid].binnedIvector[position]){

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


/*7th March '18: returned to the "Before Jan 22, '18"-version; the problem with the misisng extension from the ends
seemed to be rather that a hit !=-1 was returned from bisectSingle so in the extension part the hit == -1 was virtually 
never entered */
int bisectionalSearchSingle_old_2(int *matchIndicator, int depth, struct matchRange **ptr_newRange, struct binned_I_windows binned_I_windows_query, int position, struct matchRange range, struct binned_I_windows *ptr_binned_I_windows_DB, int allowedNrOfMismatches, int nrOfEntries){

	int returnVal = 0;

	int left = range.start;
	int right = range.end;

	int lastLeft;
	int lastRight;
	int mid = 0;
	int midNew = -1;

	int hit = -1; 
	int hitLR[2];

	//struct matchRange *ptr_range_temp;

	int i, j;


	if(position < nrOfEntries){

		/*reset*/
		for(j = 0; j < 3; j++){

			//printf("Sommeren mild oedsles kyssene nemmere pluds'ligt .. position:%d depth: %d j : %d\n", position, depth, j);

			ptr_newRange[depth][j].start = -1;
			ptr_newRange[depth][j].end = -2;
			ptr_newRange[depth][j].nrOfMismatches = 999;
		}

		/*start the search for the "letter" at the given position*/
		//printf("Range left: %d right:%d\n", range.start, range.end);
		hitLR[0] = left;
		hitLR[1] = right;
		//printf("BisectSingle input; left: %d right:%d\n", left, right);
	
		hit = bisectSingle(binned_I_windows_query, position, hitLR, ptr_binned_I_windows_DB);

		/*extend left and right from the hit if possible.
		This is here done in two further bisections:*/
		/*printf("hit:%d\n", hit);
		if(hit !=-1){
			printf("hit in DB: %d %d %d %d %d\n", ptr_binned_I_windows_DB[hit].binnedIvector[0], ptr_binned_I_windows_DB[hit].binnedIvector[1], ptr_binned_I_windows_DB[hit].binnedIvector[2], ptr_binned_I_windows_DB[hit].binnedIvector[3], ptr_binned_I_windows_DB[hit].binnedIvector[4]);
		}*/

		returnVal = bisectSingleExtend(hit, binned_I_windows_query, position, hitLR, ptr_binned_I_windows_DB, 0, 0);

		/*We now have the start and end indices of the match of the query letter (ie without mismatches).
		So:*/

		ptr_newRange[depth][1].start = hitLR[0]; //range_temp.start;
		ptr_newRange[depth][1].end = hitLR[1]; //range_temp.end;
		ptr_newRange[depth][1].nrOfMismatches = range.nrOfMismatches;
		/*printf("New range 1 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);
		if(hit !=-1){
			if(hitLR[0]-1 >= range.start){ 
				printf("DB at %d (left endpt -1): %d %d %d %d %d\n", hitLR[0]-1, ptr_binned_I_windows_DB[hitLR[0]-1].binnedIvector[0], ptr_binned_I_windows_DB[hitLR[0]-1].binnedIvector[1], ptr_binned_I_windows_DB[hitLR[0]-1].binnedIvector[2], ptr_binned_I_windows_DB[hitLR[0]-1].binnedIvector[3], ptr_binned_I_windows_DB[hitLR[0]-1].binnedIvector[4]);
			}
			printf("DB at %d (left endpt): %d %d %d %d %d\n", hitLR[0], ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[0], ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[1], ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[2], ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[3], ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[4]);
			printf("DB at %d (right endpt): %d %d %d %d %d\n", hitLR[1], ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[0], ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[1], ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[2], ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[3], ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[4]);
			if(hitLR[1] +1 <= range.end){ 
				printf("DB at %d (right endpt +1) : %d %d %d %d %d\n", hitLR[1] +1, ptr_binned_I_windows_DB[hitLR[1]+1].binnedIvector[0], ptr_binned_I_windows_DB[hitLR[1]+1].binnedIvector[1], ptr_binned_I_windows_DB[hitLR[1]+1].binnedIvector[2], ptr_binned_I_windows_DB[hitLR[1]+1].binnedIvector[3], ptr_binned_I_windows_DB[hitLR[1]+1].binnedIvector[4]);
			}
			getchar();
		}*/

		/*Next, if there is still room for mismatches, search left/right of the matching range we just found
		allowing a difference between the query letter and the match of 1:*/
		if(range.nrOfMismatches < allowedNrOfMismatches){

			/*left-ward*/
			/*init*/
			/*take one step left, if possible; if difference to the query is one go on else skip*/
			/*OBS: If hit =-1 (returned from bisectSingle) then hitLR (returned from bisectSingleExtend) is 
			[-1, -2]; in this case, if there is a discrepancy of 1 between query and db at range.end, we try
			to extend this mismatch into the range*/
			if(ptr_newRange[depth][1].start > range.start || hit == -1 ){
				if(hit != -1){//take one step left, since there was a hit at ptr_newRange[depth][1].start
					hitLR[1] = ptr_newRange[depth][1].start - 1; //right interval bound
				}
				else{//hit =-1, so there was no hit at ptr_newRange[depth][1].start; we therefore start from
					//the end of the db-range, check if there is only a mismatch of 1 and then, in next if-clause, 
					//extend from there
					hitLR[1] = range.end;
				}
				if(ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[position] == binned_I_windows_query.binnedIvector[position] - 1 ){//+1: there can be no hit at hitLR[1]; we check if the next db entry, down in lex-order, differ by -1 
		
					/*extend the range left as far as possible*/
					hitLR[0] = range.start; //left interval bound

					/*call left-ward extension*/
					returnVal = bisectSingleExtend(hitLR[1], binned_I_windows_query, position, hitLR, ptr_binned_I_windows_DB, -1, -1);

					/*record the result*/
					ptr_newRange[depth][0].start = hitLR[0]; //range_temp.start;
					ptr_newRange[depth][0].end = hitLR[1]; //range_temp.end;
					ptr_newRange[depth][0].nrOfMismatches = range.nrOfMismatches + 1;
					//printf("new range 0 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);


				}

			}

			/*right-ward:*/
			/*take one step right, if possible; if difference to the query is one go on else skip*/
			if(0 <= ptr_newRange[depth][1].end  && ptr_newRange[depth][1].end < range.end || hit == -1 ){
				if(hit != -1){//take one step right, since there was a hit at ptr_newRange[depth][1].end
					hitLR[0] = ptr_newRange[depth][1].end + 1; //left interval bound
				}
				else{//hit =-1, so there was no hit at range.start; we therefore start from
					//the start of the db-range, check if there is only a mismatch of 1 there and then, in 
					//next if-clause, extend from there
					hitLR[0] = range.start;
				}
				if(ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[position] == binned_I_windows_query.binnedIvector[position] + 1 ){ //+1: there can be no hit at hitLR[0]; we check if the next db entry, up in lex-order, differ by 1  
		
					/*extend the range left as far as possible*/
					hitLR[1] = range.end; //right interval bound

					/*call left-ward extension*/
					returnVal = bisectSingleExtend(hitLR[0], binned_I_windows_query, position, hitLR, ptr_binned_I_windows_DB, 1, -1);

					/*record the result*/
					ptr_newRange[depth][2].start = hitLR[0]; //range_temp.start;
					ptr_newRange[depth][2].end = hitLR[1]; //range_temp.end;
					ptr_newRange[depth][2].nrOfMismatches = range.nrOfMismatches + 1;
					//printf("new range 2 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);

				}

			}
	
		}

		//getchar();

		/*Update the matchIndicator if at end of recursion*/
		if(depth == nrOfEntries-1){
			for(i=0; i < 3; i++){
				for(j= ptr_newRange[depth][i].start; j <= ptr_newRange[depth][i].end; j++){

					matchIndicator[j] = 1;

				}
				//printf("Recorded matches for a range %d from %d to %d\n", i, ptr_newRange[depth][i].start, ptr_newRange[depth][i].end); 
			}
			//getchar();
		}
		else{/*Recurse*/

			for(i=0; i < 3; i++){

				if(ptr_newRange[depth][i].start <= ptr_newRange[depth][i].end){
					//printf("Recurse from range %d to depth:%d position: %d range starts: %d ends: %d\n", i, depth+1, position+1, ptr_newRange[depth][i].start, ptr_newRange[depth][i].end);
					//getchar();
					returnVal = bisectionalSearchSingle(matchIndicator, depth+1, ptr_newRange, binned_I_windows_query, position + 1 , ptr_newRange[depth][i], ptr_binned_I_windows_DB, allowedNrOfMismatches, nrOfEntries);
				}
			}
		}

	}

	return returnVal;

}


/*Before Jan 22, '18; extension from the ends in case no hit was found was missing 
7th March '18: returned to this version; the problem seemed to be rather that a hit !=-1 was
returned from bisectSingle so in the extension part the hit == -1 was virtually never entered */
int bisectionalSearchSingle(int *matchIndicator, int depth, struct matchRange **ptr_newRange, struct binned_I_windows binned_I_windows_query, int position, struct matchRange range, struct binned_I_windows *ptr_binned_I_windows_DB, int allowedNrOfMismatches, int nrOfEntries){

	int returnVal = 0;

	int left = range.start;
	int right = range.end;

	int lastLeft;
	int lastRight;
	int mid = 0;
	int midNew = -1;

	int hit = -1; 
	int hitLR[2];

	struct matchRange *ptr_range_temp;

	int i, j;


	if(position < nrOfEntries){

		/*reset*/
		for(j = 0; j < 3; j++){
				ptr_newRange[depth][j].start = -1;
				ptr_newRange[depth][j].end = -2;
				ptr_newRange[depth][j].nrOfMismatches = 999;
			}

		/*start the search for the "letter" at the given position*/

		hitLR[0] = left;
		hitLR[1] = right;
		//printf("BisectSingle input; left: %d right:%d\n", left, right);
	
		hit = bisectSingle(binned_I_windows_query, position, hitLR, ptr_binned_I_windows_DB);

		/*extend left and right from the hit if possible.
		This is here done in two further bisections:*/
		/*printf("hit:%d\n", hit);
		if(hit !=-1){
			printf("hit in DB: %d %d %d %d %d\n", ptr_binned_I_windows_DB[hit].binnedIvector[0], ptr_binned_I_windows_DB[hit].binnedIvector[1], ptr_binned_I_windows_DB[hit].binnedIvector[2], ptr_binned_I_windows_DB[hit].binnedIvector[3], ptr_binned_I_windows_DB[hit].binnedIvector[4]);
		}*/

		returnVal = bisectSingleExtend(hit, binned_I_windows_query, position, hitLR, ptr_binned_I_windows_DB, 0, 0);

		/*We now have the start and end indices of the match of the query letter (ie without mismatches).
		So:*/

		ptr_newRange[depth][1].start = hitLR[0]; //range_temp.start;
		ptr_newRange[depth][1].end = hitLR[1]; //range_temp.end;
		ptr_newRange[depth][1].nrOfMismatches = range.nrOfMismatches;
		/*printf("New range 1 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);
		if(hit !=-1){
			printf("DB at %d (left endpt -1): %d %d %d %d %d\n", hitLR[0]-1, ptr_binned_I_windows_DB[hitLR[0]-1].binnedIvector[0], ptr_binned_I_windows_DB[hitLR[0]-1].binnedIvector[1], ptr_binned_I_windows_DB[hitLR[0]-1].binnedIvector[2], ptr_binned_I_windows_DB[hitLR[0]-1].binnedIvector[3], ptr_binned_I_windows_DB[hitLR[0]-1].binnedIvector[4]);
			printf("DB at %d (left endpt): %d %d %d %d %d\n", hitLR[0], ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[0], ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[1], ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[2], ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[3], ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[4]);
			printf("DB at %d (right endpt): %d %d %d %d %d\n", hitLR[1], ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[0], ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[1], ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[2], ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[3], ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[4]);
			printf("DB at %d (right endpt +1) : %d %d %d %d %d\n", hitLR[1] +1, ptr_binned_I_windows_DB[hitLR[1]+1].binnedIvector[0], ptr_binned_I_windows_DB[hitLR[1]+1].binnedIvector[1], ptr_binned_I_windows_DB[hitLR[1]+1].binnedIvector[2], ptr_binned_I_windows_DB[hitLR[1]+1].binnedIvector[3], ptr_binned_I_windows_DB[hitLR[1]+1].binnedIvector[4]);

		}*/

		/*Next, if there is still room for mismatches, search left/right of the matching range we just found
		allowing a difference between the query letter and the match of 1:*/
		if(range.nrOfMismatches < allowedNrOfMismatches){

			/*left-ward*/
			/*init*/
			/*take one step left, if possible; if difference to the query is one go on else skip*/
			if(ptr_newRange[depth][1].start > range.start){
				if(hit != -1){//take one step left, since there was a hit at ptr_newRange[depth][1].start
					hitLR[1] = ptr_newRange[depth][1].start - 1; //right interval bound
				}
				else{//hit =-1, so there was no hit at ptr_newRange[depth][1].start; we therefore start from there
					hitLR[1] = ptr_newRange[depth][1].start;
				}
				if(ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[position] == binned_I_windows_query.binnedIvector[position] - 1 ){//+1: there can be no hit at hitLR[1]; we check if the next db entry, down in lex-order, differ by -1 
		
					/*extend the range left as far as possible*/
					hitLR[0] = range.start; //left interval bound

					/*call left-ward extension*/
					returnVal = bisectSingleExtend(hitLR[1], binned_I_windows_query, position, hitLR, ptr_binned_I_windows_DB, -1, -1);

					/*record the result*/
					ptr_newRange[depth][0].start = hitLR[0]; //range_temp.start;
					ptr_newRange[depth][0].end = hitLR[1]; //range_temp.end;
					ptr_newRange[depth][0].nrOfMismatches = range.nrOfMismatches + 1;
					//printf("new range 0 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);


				}

			}

			/*right-ward:*/
			/*take one step right, if possible; if difference to the query is one go on else skip*/
			if(0 <= ptr_newRange[depth][1].end  && ptr_newRange[depth][1].end < range.end){
				if(hit != -1){//take one step right, since there was a hit at ptr_newRange[depth][1].end
					hitLR[0] = ptr_newRange[depth][1].end + 1; //left interval bound
				}
				else{//hit =-1, so there was no hit at ptr_newRange[depth][1].end; we therefore start from there
					hitLR[0] = ptr_newRange[depth][1].end;
				}
				if(ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[position] == binned_I_windows_query.binnedIvector[position] + 1 ){ //+1: there can be no hit at hitLR[0]; we check if the next db entry, up in lex-order, differ by 1  
		
					/*extend the range left as far as possible*/
					hitLR[1] = range.end; //right interval bound

					/*call left-ward extension*/
					returnVal = bisectSingleExtend(hitLR[0], binned_I_windows_query, position, hitLR, ptr_binned_I_windows_DB, 1, -1);

					/*record the result*/
					ptr_newRange[depth][2].start = hitLR[0]; //range_temp.start;
					ptr_newRange[depth][2].end = hitLR[1]; //range_temp.end;
					ptr_newRange[depth][2].nrOfMismatches = range.nrOfMismatches + 1;
					//printf("new range 2 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);

				}

			}
	
		}

		//getchar();

		/*Update the matchIndicator if at end of recursion*/
		if(depth == nrOfEntries-1){
			for(i=0; i < 3; i++){
				for(j= ptr_newRange[depth][i].start; j <= ptr_newRange[depth][i].end; j++){

					matchIndicator[j] = 1;

				}
				//printf("Recorded matches for a range %d from %d to %d\n", i, ptr_newRange[depth][i].start, ptr_newRange[depth][i].end); 
			}
			//getchar();
		}
		else{/*Recurse*/

			for(i=0; i < 3; i++){

				if(ptr_newRange[depth][i].start <= ptr_newRange[depth][i].end){
					//printf("Recurse from range %d to depth:%d position: %d range starts: %d ends: %d\n", i, depth+1, position+1, ptr_newRange[depth][i].start, ptr_newRange[depth][i].end);
					returnVal = bisectionalSearchSingle(matchIndicator, depth+1, ptr_newRange, binned_I_windows_query, position + 1 , ptr_newRange[depth][i], ptr_binned_I_windows_DB, allowedNrOfMismatches, nrOfEntries);
				}
			}
		}

	}

	return returnVal;

}


/*Before april 4, '17:*/
int bisectionalSearchSingle_old(int *matchIndicator, int depth, struct matchRange **ptr_newRange, struct binned_I_windows binned_I_windows_query, int position, struct matchRange range, struct binned_I_windows *ptr_binned_I_windows_DB, int allowedNrOfMismatches, int nrOfEntries){

	int returnVal = 0;

	int left = range.start;
	int right = range.end;

	int lastLeft;
	int lastRight;
	int mid = 0;
	int midNew = -1;

	int hit = -1; 
	int hitLR[2];

	struct matchRange *ptr_range_temp;

	int i, j;


	if(position < nrOfEntries){

		/*reset*/
		for(j = 0; j < 3; j++){
				ptr_newRange[depth][j].start = -1;
				ptr_newRange[depth][j].end = -2;
				ptr_newRange[depth][j].nrOfMismatches = 999;
			}

		/*start the search for the "letter" at the given position*/

		hitLR[0] = left;
		hitLR[1] = right;
		//printf("BisectSingle input; left: %d right:%d\n", left, right);
	
		hit = bisectSingle(binned_I_windows_query, position, hitLR, ptr_binned_I_windows_DB);

		/*extend left and right from the hit if possible.
		This is here done in two further bisections:*/
		/*printf("hit:%d\n", hit);
		if(hit !=-1){
			printf("DB: %d %d %d %d %d\n", ptr_binned_I_windows_DB[hit].binnedIvector[0], ptr_binned_I_windows_DB[hit].binnedIvector[1], ptr_binned_I_windows_DB[hit].binnedIvector[2], ptr_binned_I_windows_DB[hit].binnedIvector[3], ptr_binned_I_windows_DB[hit].binnedIvector[4]);
		}*/

		returnVal = bisectSingleExtend_old(hit, binned_I_windows_query, position, hitLR, ptr_binned_I_windows_DB, 0);

		/*We now have the start and end indices of the match of the query letter (ie without mismatches).
		So:*/

		ptr_newRange[depth][1].start = hitLR[0]; //range_temp.start;
		ptr_newRange[depth][1].end = hitLR[1]; //range_temp.end;
		ptr_newRange[depth][1].nrOfMismatches = range.nrOfMismatches;
		/*printf("new range 1 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);
		if(hit !=-1){
			printf("DB at left endpt: %d %d %d %d %d\n", ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[0], ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[1], ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[2], ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[3], ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[4]);
			printf("DB at right endpt: %d %d %d %d %d\n", ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[0], ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[1], ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[2], ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[3], ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[4]);
		}*/

		/*Next, if there is still room for mismatches, search left/right of the matching range we just found
		allowing a difference between the query letter and the match of 1:*/
		if(range.nrOfMismatches < allowedNrOfMismatches){

			/*left-ward*/
			/*init*/
			/*take one step left, if possible; if difference to the query is one go on else skip*/
			if(ptr_newRange[depth][1].start > range.start){
				hitLR[1] = ptr_newRange[depth][1].start - 1; //right interval bound
				if(ptr_binned_I_windows_DB[hitLR[1]].binnedIvector[position] == binned_I_windows_query.binnedIvector[position] - 1 ){
		
					/*extend the range left as far as possible*/
					hitLR[0] = range.start; //left interval bound

					/*call left-ward extension*/
					returnVal = bisectSingleExtend_old(hitLR[1], binned_I_windows_query, position, hitLR, ptr_binned_I_windows_DB, -1);

					/*record the result*/
					ptr_newRange[depth][0].start = hitLR[0]; //range_temp.start;
					ptr_newRange[depth][0].end = hitLR[1]; //range_temp.end;
					ptr_newRange[depth][0].nrOfMismatches = range.nrOfMismatches + 1;
					//printf("new range 0 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);


				}

			}

			/*right-ward:*/
			/*take one step right, if possible; if difference to the query is one go on else skip*/
			if(0 <= ptr_newRange[depth][1].end  && ptr_newRange[depth][1].end < range.end){
				hitLR[0] = ptr_newRange[depth][1].end + 1; //left interval bound
				if(ptr_binned_I_windows_DB[hitLR[0]].binnedIvector[position] == binned_I_windows_query.binnedIvector[position] + 1 ){
		
					/*extend the range left as far as possible*/
					hitLR[1] = range.end; //right interval bound

					/*call left-ward extension*/
					returnVal = bisectSingleExtend_old(hitLR[0], binned_I_windows_query, position, hitLR, ptr_binned_I_windows_DB, 1);

					/*record the result*/
					ptr_newRange[depth][2].start = hitLR[0]; //range_temp.start;
					ptr_newRange[depth][2].end = hitLR[1]; //range_temp.end;
					ptr_newRange[depth][2].nrOfMismatches = range.nrOfMismatches + 1;
					//printf("new range 2 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);

				}

			}
	
		}

		//getchar();

		/*Update the matchIndicator if at end of recursion*/
		if(depth == nrOfEntries-1){
			for(i=0; i < 3; i++){
				for(j= ptr_newRange[depth][i].start; j <= ptr_newRange[depth][i].end; j++){

					matchIndicator[j] = 1;

				}
			}
		}
		else{/*Recurse*/

			for(i=0; i < 3; i++){

				returnVal = bisectionalSearchSingle(matchIndicator, depth+1, ptr_newRange, binned_I_windows_query, position + 1 , ptr_newRange[depth][i], ptr_binned_I_windows_DB, allowedNrOfMismatches, nrOfEntries);

			}
		}

	}

	return returnVal;

}



/*bisectional search for "letter" at position position of query, binned_I_windows_query, among the letters at the same position
in the range of the array ptr_binned_I_windows_DB from  ptr_range.start to ptr_range.end.
Returns: index of the hit; updates ptr_range.start and .end with the last left resp. right interval bounds of the search*/
int bisectSingle_Pairs(struct binned_I_windowPairs binned_I_windowPairs_query, int position, int *ptr_range, struct binned_I_windowPairs *ptr_binned_I_windowPairs_DB){

	int returnVal = 0;

	int iter = 0;

	int left = ptr_range[0];
	int right = ptr_range[1];

	//int lastLeft = left;
	int lastRight = right;
	int mid = -117;

	int hit = -1; 

	/*if the search range dow not contain the value, binned_I_windowPairs_query.binnedIvector[position], there is no
	hit to be found, and hit = -1 is returned (and the input range is left as it is)*/
	if(binned_I_windowPairs_query.binnedIvector[position] < ptr_binned_I_windowPairs_DB[left].binnedIvector[position] ||
		binned_I_windowPairs_query.binnedIvector[position] > ptr_binned_I_windowPairs_DB[right].binnedIvector[position] ) {

			return hit; 

	}


	/*The folowing if-clause seems superfluous by the preceding if-clause and since the while-loop to follow
	is now entered when left = right; however the hit value must be updated from -1 to the index (=left)*/
	//If the range is actually a single value and we have a hit we return that// 
	if(left == right){
		if(binned_I_windowPairs_query.binnedIvector[position] == ptr_binned_I_windowPairs_DB[left].binnedIvector[position]){
			hit = left;
		}
		//else: hit =-1 by the set initial value

		return hit;

	}
	
	//.. And else we search://
	while(mid != lastRight && left < right){//while(right -left > 1){ //while(mid  != lastRight && left <= right)

		//lastLeft = left;
		lastRight = right;

		mid = (int) ceil( (double) (right + left)/2 );

		//printf("mid:%d\n", mid);

		if( binned_I_windowPairs_query.binnedIvector[position] < ptr_binned_I_windowPairs_DB[mid].binnedIvector[position]){

			right = mid;

		}
		else if( binned_I_windowPairs_query.binnedIvector[position] > ptr_binned_I_windowPairs_DB[mid].binnedIvector[position] ){

			left = mid;
		}
		else{ //lex == 0 so we have a hit
				hit = mid;
				//printf("break\n");
				break; //break while-loop
		}

		/*there remains the possibility that the hit sits at the left end point (added 7th March '18)*/
		if (binned_I_windowPairs_query.binnedIvector[position] == ptr_binned_I_windowPairs_DB[left].binnedIvector[position]) {
			hit = left;
		}

	}

	ptr_range[0] = left;
	ptr_range[1] = right;

	return hit;
}

/*Before 7th March '18; this version returns a hit != -1 even if there is none (the sought value sits either outside
the range sought in, or within a length-1 interval returned)*/
int bisectSingle_Pairs_old(struct binned_I_windowPairs binned_I_windowPairs_query, int position, int *ptr_range, struct binned_I_windowPairs *ptr_binned_I_windowPairs_DB) {

	int returnVal = 0;

	int iter = 0;

	int left = ptr_range[0];
	int right = ptr_range[1];

	//int lastLeft = left;
	int lastRight = right;
	int mid = -117;

	int hit = -1;

	/*if the search range dow not contain the value, binned_I_windowPairs_query.binnedIvector[position], there is no
	hit to be found, and hit = -1 is returned */
	if (binned_I_windowPairs_query.binnedIvector[position] < ptr_binned_I_windowPairs_DB[left].binnedIvector[position] ||
		binned_I_windowPairs_query.binnedIvector[position] > ptr_binned_I_windowPairs_DB[right].binnedIvector[position]) {

		return hit;

	}	/*If the range does not contain the value binned_I_windowPairs_query.binnedIvector[position] there is no


		/*The folowing if-clause seems superfluous by the preceding if-clause and since the while-loop to follow
		is now entered when left = right; however the hit value must be updated from -1 to the index (=left)*/
		//If the range is actually a single value and we have a hit we return that// 
	if (left == right) {
		if (binned_I_windowPairs_query.binnedIvector[position] == ptr_binned_I_windowPairs_DB[left].binnedIvector[position]) {
			hit = left;
		}
		//else: hit =-1 by the set initial value

		return hit;

	}

	//.. And else we search://
	while (mid != lastRight && left < right) {//while(right -left > 1){ //while(mid  != lastRight && left <= right)

											  //lastLeft = left;
		lastRight = right;

		mid = (int)ceil((double)(right + left) / 2);

		//printf("mid:%d\n", mid);

		if (binned_I_windowPairs_query.binnedIvector[position] < ptr_binned_I_windowPairs_DB[mid].binnedIvector[position]) {

			hit = left;

			right = mid;


		}
		else if (binned_I_windowPairs_query.binnedIvector[position] > ptr_binned_I_windowPairs_DB[mid].binnedIvector[position]) {

			hit = right;

			left = mid;
		}
		else { //lex == 0 so we have a hit "already" 
			hit = mid;
			//printf("break\n");
			break; //break while-loop
		}

		//midNew = (int) ceil( (double) (right + left)/2 );

	}

	ptr_range[0] = left;
	ptr_range[1] = right;

	return hit;
}

/*for extending left or right from a hit (with index hitIndex); searches for extension left from ptr_range.start and up to 
hitIndex, and right from hitIndex and to ptr_range.end. In particular hitIndex must sit in the interval [ptr_range.start,
ptr_range.end] for this to work.
Input: leftRigth_b -- if 0 both extensions are done; if -1 only left, if 1 only right.
Transforms: updates ptr_range .start and .end to left resp right indexes found.*/
int bisectSingleExtend_Pairs(int hitIndex, struct binned_I_windowPairs binned_I_windowPairs_query, int position, int *ptr_range, struct binned_I_windowPairs *ptr_binned_I_windowPairs_DB, int leftRight_b, int discrepancy){

	int returnVal = 0;

	int left;
	int lastLeft;
	int right;
	int lastRight;
	int mid = -1;
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
		
				if( binned_I_windowPairs_query.binnedIvector[position] == ptr_binned_I_windowPairs_DB[mid].binnedIvector[position] - discrepancy){

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
			while(lastRight - mid > 0) {

				lastRight =  right;
				mid = (int) ceil( (double) (right + left)/2 );

				//printf("mid:%d\n", mid);
		
				if( binned_I_windowPairs_query.binnedIvector[position] == ptr_binned_I_windowPairs_DB[mid].binnedIvector[position] + discrepancy){

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


/*7th March '18: changed similarly to single windows version (see above)*/
int bisectionalSearchSingle_Pairs(int *matchIndicator, int depth, struct matchRange **ptr_newRange, struct binned_I_windowPairs binned_I_windowPairs_query, int position, struct matchRange range, struct binned_I_windowPairs *ptr_binned_I_windowPairs_DB, int allowedNrOfMismatches, int nrOfEntries){

	int returnVal = 0;

	int left = range.start;
	int right = range.end;

	int lastLeft;
	int lastRight;
	int mid = 0;
	int midNew = -1;

	int hit = -1; 
	int hitLR[2];

	//struct matchRange *ptr_range_temp;

	int i, j;


	if(position < nrOfEntries){

		/*reset*/
		for(j = 0; j < 3; j++){

			//printf("Langmodigt venter bysvalen stædigt stadig længes .. position:%d depth: %d j : %d\n", position, depth, j);

			ptr_newRange[depth][j].start = -1;
			ptr_newRange[depth][j].end = -2;
			ptr_newRange[depth][j].nrOfMismatches = 999;
		}

		/*start the search for the "letter" at the given position*/
		//printf("Range left: %d right:%d\n", range.start, range.end);
		hitLR[0] = left;
		hitLR[1] = right;
		//printf("BisectSingle input; left: %d right:%d\n", left, right);
	
		hit = bisectSingle_Pairs(binned_I_windowPairs_query, position, hitLR, ptr_binned_I_windowPairs_DB);

		/*extend left and right from the hit if possible.
		This is here done in two further bisections:*/
		/*printf("hit:%d\n", hit);
		if(hit !=-1){
			printf("hit in DB: %d %d %d %d %d\n", ptr_binned_I_windowPairs_DB[hit].binnedIvector[0], ptr_binned_I_windowPairs_DB[hit].binnedIvector[1], ptr_binned_I_windowPairs_DB[hit].binnedIvector[2], ptr_binned_I_windowPairs_DB[hit].binnedIvector[3], ptr_binned_I_windowPairs_DB[hit].binnedIvector[4]);
		}*/

		returnVal = bisectSingleExtend_Pairs(hit, binned_I_windowPairs_query, position, hitLR, ptr_binned_I_windowPairs_DB, 0, 0);

		/*We now have the start and end indices of the match of the query letter (ie without mismatches).
		So:*/

		ptr_newRange[depth][1].start = hitLR[0]; //range_temp.start;
		ptr_newRange[depth][1].end = hitLR[1]; //range_temp.end;
		ptr_newRange[depth][1].nrOfMismatches = range.nrOfMismatches;
		/*printf("New range 1 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);
		if(hit !=-1){
			if(hitLR[0]-1 >= range.start){ 
				printf("DB at %d (left endpt -1): %d %d %d %d %d\n", hitLR[0]-1, ptr_binned_I_windowPairs_DB[hitLR[0]-1].binnedIvector[0], ptr_binned_I_windowPairs_DB[hitLR[0]-1].binnedIvector[1], ptr_binned_I_windowPairs_DB[hitLR[0]-1].binnedIvector[2], ptr_binned_I_windowPairs_DB[hitLR[0]-1].binnedIvector[3], ptr_binned_I_windowPairs_DB[hitLR[0]-1].binnedIvector[4]);
			}
			printf("DB at %d (left endpt): %d %d %d %d %d\n", hitLR[0], ptr_binned_I_windowPairs_DB[hitLR[0]].binnedIvector[0], ptr_binned_I_windowPairs_DB[hitLR[0]].binnedIvector[1], ptr_binned_I_windowPairs_DB[hitLR[0]].binnedIvector[2], ptr_binned_I_windowPairs_DB[hitLR[0]].binnedIvector[3], ptr_binned_I_windowPairs_DB[hitLR[0]].binnedIvector[4]);
			printf("DB at %d (right endpt): %d %d %d %d %d\n", hitLR[1], ptr_binned_I_windowPairs_DB[hitLR[1]].binnedIvector[0], ptr_binned_I_windowPairs_DB[hitLR[1]].binnedIvector[1], ptr_binned_I_windowPairs_DB[hitLR[1]].binnedIvector[2], ptr_binned_I_windowPairs_DB[hitLR[1]].binnedIvector[3], ptr_binned_I_windowPairs_DB[hitLR[1]].binnedIvector[4]);
			if(hitLR[1] +1 <= range.end){ 
				printf("DB at %d (right endpt +1) : %d %d %d %d %d\n", hitLR[1] +1, ptr_binned_I_windowPairs_DB[hitLR[1]+1].binnedIvector[0], ptr_binned_I_windowPairs_DB[hitLR[1]+1].binnedIvector[1], ptr_binned_I_windowPairs_DB[hitLR[1]+1].binnedIvector[2], ptr_binned_I_windowPairs_DB[hitLR[1]+1].binnedIvector[3], ptr_binned_I_windowPairs_DB[hitLR[1]+1].binnedIvector[4]);
			}
			getchar();
		}*/

		/*Next, if there is still room for mismatches, search left/right of the matching range we just found
		allowing a difference between the query letter and the match of 1:*/
		if(range.nrOfMismatches < allowedNrOfMismatches){

			/*left-ward*/
			/*init*/
			/*take one step left, if possible; if difference to the query is one go on else skip*/
			/*OBS: If hit =-1 (returned from bisectSingle) then hitLR (returned from bisectSingleExtend) is 
			[-1, -2]; in this case, if there is a discrepancy of 1 between query and db at range.end, we try
			to extend this mismatch into the range*/
			if(ptr_newRange[depth][1].start > range.start || hit == -1 ){
				if(hit != -1){//take one step left, since there was a hit at ptr_newRange[depth][1].start
					hitLR[1] = ptr_newRange[depth][1].start - 1; //right interval bound
				}
				else{//hit =-1, so there was no hit at ptr_newRange[depth][1].start; we therefore start from
					//ptr_newRange[depth][1].start, check in next if-clause if there is only a mismatch of -1 and if so 
					//extend from there. Obs: it is possible that the sought value does not sit in the DB-range at all
					//and then ptr_newRange[depth][1].start = range.start why no extension will be possible (but we still
					//should catch a case of a discrepancy of -1 at that range start)
					hitLR[1] = max(ptr_newRange[depth][1].start, range.start); //take the max since ptr_newRange[depth][1].start can be < 0 (so < range.start), in which case there's no db-element at ptr_newRange[depth][1].start
					//printf("(saaledes stadig ganske hemmelig) ... %d %d %d\n", ptr_newRange[depth][1].start, range.start, ptr_binned_I_windowPairs_DB[hitLR[1]].binnedIvector[position]);
					//getchar();
				}
				if(ptr_binned_I_windowPairs_DB[hitLR[1]].binnedIvector[position] == binned_I_windowPairs_query.binnedIvector[position] - 1 ){//+1: there can be no hit at hitLR[1]; we check if the next db entry, down in lex-order, differ by -1 
		
					/*extend the range left as far as possible*/
					hitLR[0] = range.start; //left interval bound

					/*call left-ward extension*/
					returnVal = bisectSingleExtend_Pairs(hitLR[1], binned_I_windowPairs_query, position, hitLR, ptr_binned_I_windowPairs_DB, -1, -1);

					/*record the result*/
					ptr_newRange[depth][0].start = hitLR[0]; //range_temp.start;
					ptr_newRange[depth][0].end = hitLR[1]; //range_temp.end;
					ptr_newRange[depth][0].nrOfMismatches = range.nrOfMismatches + 1;
					//printf("new range 0 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);


				}

			}

			/*right-ward:*/
			/*take one step right, if possible; if difference to the query is one go on else skip*/
			if(0 <= ptr_newRange[depth][1].end  && ptr_newRange[depth][1].end < range.end || hit == -1 ){
				if(hit != -1){//take one step right, since there was a hit at ptr_newRange[depth][1].end
					hitLR[0] = ptr_newRange[depth][1].end + 1; //left interval bound
				}
				else{//hit =-1, so there was no hit at ptr_newRange[depth][1].and; we therefore start from
					//ptr_newRange[depth][1].end, check in next if-clause if there is only a mismatch of +1 and if so 
					//extend from there. Obs: it is possible that the sought value does not sit in the DB-range at all
					//and then ptr_newRange[depth][1].end = range.end why no extension will be possible (but we still
					//should catch a case of a discrepancy of +1 at that range end)
					hitLR[0] = max(ptr_newRange[depth][1].end, range.end); //take the max since ptr_newRange[depth][1].end can be < 0, in which case there's no db-element at ptr_newRange[depth][1].end
					//printf("(saaledes stadig ganske hemmelig) ... %d %d %d\n", ptr_newRange[depth][1].end, range.end, ptr_binned_I_windowPairs_DB[hitLR[1]].binnedIvector[position]);
					//getchar();
				}
				if(ptr_binned_I_windowPairs_DB[hitLR[0]].binnedIvector[position] == binned_I_windowPairs_query.binnedIvector[position] + 1 ){ //+1: there can be no hit at hitLR[0]; we check if the next db entry, up in lex-order, differ by 1  
		
					/*extend the range left as far as possible*/
					hitLR[1] = range.end; //right interval bound

					/*call left-ward extension*/
					returnVal = bisectSingleExtend_Pairs(hitLR[0], binned_I_windowPairs_query, position, hitLR, ptr_binned_I_windowPairs_DB, 1, -1);

					/*record the result*/
					ptr_newRange[depth][2].start = hitLR[0]; //range_temp.start;
					ptr_newRange[depth][2].end = hitLR[1]; //range_temp.end;
					ptr_newRange[depth][2].nrOfMismatches = range.nrOfMismatches + 1;
					//printf("new range 2 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);

				}

			}
	
		}

		//getchar();

		/*Update the matchIndicator if at end of recursion*/
		if(depth == nrOfEntries-1){
			for(i=0; i < 3; i++){
				for(j= ptr_newRange[depth][i].start; j <= ptr_newRange[depth][i].end; j++){

					matchIndicator[j] = 1;

				}
				//printf("Recorded matches for a range %d from %d to %d\n", i, ptr_newRange[depth][i].start, ptr_newRange[depth][i].end); 
			}
			//getchar();
		}
		else{/*Recurse*/

			for(i=0; i < 3; i++){

				if(ptr_newRange[depth][i].start <= ptr_newRange[depth][i].end){
					//printf("Recurse from range %d to depth:%d position: %d range starts: %d ends: %d\n", i, depth+1, position+1, ptr_newRange[depth][i].start, ptr_newRange[depth][i].end);
					//getchar();
					returnVal = bisectionalSearchSingle_Pairs(matchIndicator, depth+1, ptr_newRange, binned_I_windowPairs_query, position + 1 , ptr_newRange[depth][i], ptr_binned_I_windowPairs_DB, allowedNrOfMismatches, nrOfEntries);
				}
			}
		}

	}

	return returnVal;

}

/*7th March '18: changed similarly to single windows version (see above)*/
int bisectionalSearchSingle_Pairs_old(int *matchIndicator, int depth, struct matchRange **ptr_newRange, struct binned_I_windowPairs binned_I_windowPairs_query, int position, struct matchRange range, struct binned_I_windowPairs *ptr_binned_I_windowPairs_DB, int allowedNrOfMismatches, int nrOfEntries) {

	int returnVal = 0;

	int left = range.start;
	int right = range.end;

	int lastLeft;
	int lastRight;
	int mid = 0;
	int midNew = -1;

	int hit = -1;
	int hitLR[2];

	//struct matchRange *ptr_range_temp;

	int i, j;


	if (position < nrOfEntries) {

		/*reset*/
		for (j = 0; j < 3; j++) {

			//printf("Sommeren mild oedsles .. position:%d depth: %d j : %d\n", position, depth, j);

			ptr_newRange[depth][j].start = -1;
			ptr_newRange[depth][j].end = -2;
			ptr_newRange[depth][j].nrOfMismatches = 999;
		}

		/*start the search for the "letter" at the given position*/
		//printf("Range left: %d right:%d\n", range.start, range.end);
		hitLR[0] = left;
		hitLR[1] = right;
		//printf("BisectSingle input; left: %d right:%d\n", left, right);

		hit = bisectSingle_Pairs(binned_I_windowPairs_query, position, hitLR, ptr_binned_I_windowPairs_DB);

		/*extend left and right from the hit if possible.
		This is here done in two further bisections:*/
		/*printf("hit:%d\n", hit);
		if(hit !=-1){
		printf("hit in DB: %d %d %d %d %d\n", ptr_binned_I_windowPairs_DB[hit].binnedIvector[0], ptr_binned_I_windowPairs_DB[hit].binnedIvector[1], ptr_binned_I_windowPairs_DB[hit].binnedIvector[2], ptr_binned_I_windowPairs_DB[hit].binnedIvector[3], ptr_binned_I_windowPairs_DB[hit].binnedIvector[4]);
		}*/

		returnVal = bisectSingleExtend_Pairs(hit, binned_I_windowPairs_query, position, hitLR, ptr_binned_I_windowPairs_DB, 0, 0);

		/*We now have the start and end indices of the match of the query letter (ie without mismatches).
		So:*/

		ptr_newRange[depth][1].start = hitLR[0]; //range_temp.start;
		ptr_newRange[depth][1].end = hitLR[1]; //range_temp.end;
		ptr_newRange[depth][1].nrOfMismatches = range.nrOfMismatches;
		/*printf("New range 1 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);
		if(hit !=-1){
		if(hitLR[0]-1 >= range.start){
		printf("DB at %d (left endpt -1): %d %d %d %d %d\n", hitLR[0]-1, ptr_binned_I_windowPairs_DB[hitLR[0]-1].binnedIvector[0], ptr_binned_I_windowPairs_DB[hitLR[0]-1].binnedIvector[1], ptr_binned_I_windowPairs_DB[hitLR[0]-1].binnedIvector[2], ptr_binned_I_windowPairs_DB[hitLR[0]-1].binnedIvector[3], ptr_binned_I_windowPairs_DB[hitLR[0]-1].binnedIvector[4]);
		}
		printf("DB at %d (left endpt): %d %d %d %d %d\n", hitLR[0], ptr_binned_I_windowPairs_DB[hitLR[0]].binnedIvector[0], ptr_binned_I_windowPairs_DB[hitLR[0]].binnedIvector[1], ptr_binned_I_windowPairs_DB[hitLR[0]].binnedIvector[2], ptr_binned_I_windowPairs_DB[hitLR[0]].binnedIvector[3], ptr_binned_I_windowPairs_DB[hitLR[0]].binnedIvector[4]);
		printf("DB at %d (right endpt): %d %d %d %d %d\n", hitLR[1], ptr_binned_I_windowPairs_DB[hitLR[1]].binnedIvector[0], ptr_binned_I_windowPairs_DB[hitLR[1]].binnedIvector[1], ptr_binned_I_windowPairs_DB[hitLR[1]].binnedIvector[2], ptr_binned_I_windowPairs_DB[hitLR[1]].binnedIvector[3], ptr_binned_I_windowPairs_DB[hitLR[1]].binnedIvector[4]);
		if(hitLR[1] +1 <= range.end){
		printf("DB at %d (right endpt +1) : %d %d %d %d %d\n", hitLR[1] +1, ptr_binned_I_windowPairs_DB[hitLR[1]+1].binnedIvector[0], ptr_binned_I_windowPairs_DB[hitLR[1]+1].binnedIvector[1], ptr_binned_I_windowPairs_DB[hitLR[1]+1].binnedIvector[2], ptr_binned_I_windowPairs_DB[hitLR[1]+1].binnedIvector[3], ptr_binned_I_windowPairs_DB[hitLR[1]+1].binnedIvector[4]);
		}
		getchar();
		}*/

		/*Next, if there is still room for mismatches, search left/right of the matching range we just found
		allowing a difference between the query letter and the match of 1:*/
		if (range.nrOfMismatches < allowedNrOfMismatches) {

			/*left-ward*/
			/*init*/
			/*take one step left, if possible; if difference to the query is one go on else skip*/
			/*OBS: If hit =-1 (returned from bisectSingle) then hitLR (returned from bisectSingleExtend) is
			[-1, -2]; in this case, if there is a discrepancy of 1 between query and db at range.end, we try
			to extend this mismatch into the range*/
			if (ptr_newRange[depth][1].start > range.start || hit == -1) {
				if (hit != -1) {//take one step left, since there was a hit at ptr_newRange[depth][1].start
					hitLR[1] = ptr_newRange[depth][1].start - 1; //right interval bound
				}
				else {//hit =-1, so there was no hit at ptr_newRange[depth][1].start; we therefore start from
					  //the end of the db-range, check if there is only a mismatch of 1 and then, in next if-clause, 
					  //extend from there
					hitLR[1] = range.end;
				}
				if (ptr_binned_I_windowPairs_DB[hitLR[1]].binnedIvector[position] == binned_I_windowPairs_query.binnedIvector[position] - 1) {//+1: there can be no hit at hitLR[1]; we check if the next db entry, down in lex-order, differ by -1 

																																			  /*extend the range left as far as possible*/
					hitLR[0] = range.start; //left interval bound

											/*call left-ward extension*/
					returnVal = bisectSingleExtend_Pairs(hitLR[1], binned_I_windowPairs_query, position, hitLR, ptr_binned_I_windowPairs_DB, -1, -1);

					/*record the result*/
					ptr_newRange[depth][0].start = hitLR[0]; //range_temp.start;
					ptr_newRange[depth][0].end = hitLR[1]; //range_temp.end;
					ptr_newRange[depth][0].nrOfMismatches = range.nrOfMismatches + 1;
					//printf("new range 0 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);


				}

			}

			/*right-ward:*/
			/*take one step right, if possible; if difference to the query is one go on else skip*/
			if (0 <= ptr_newRange[depth][1].end  && ptr_newRange[depth][1].end < range.end || hit == -1) {
				if (hit != -1) {//take one step right, since there was a hit at ptr_newRange[depth][1].end
					hitLR[0] = ptr_newRange[depth][1].end + 1; //left interval bound
				}
				else {//hit =-1, so there was no hit at range.start; we therefore start from
					  //the start of the db-range, check if there is only a mismatch of 1 there and then, in 
					  //next if-clause, extend from there
					hitLR[0] = range.start;
				}
				if (ptr_binned_I_windowPairs_DB[hitLR[0]].binnedIvector[position] == binned_I_windowPairs_query.binnedIvector[position] + 1) { //+1: there can be no hit at hitLR[0]; we check if the next db entry, up in lex-order, differ by 1  

																																			   /*extend the range left as far as possible*/
					hitLR[1] = range.end; //right interval bound

										  /*call left-ward extension*/
					returnVal = bisectSingleExtend_Pairs(hitLR[0], binned_I_windowPairs_query, position, hitLR, ptr_binned_I_windowPairs_DB, 1, -1);

					/*record the result*/
					ptr_newRange[depth][2].start = hitLR[0]; //range_temp.start;
					ptr_newRange[depth][2].end = hitLR[1]; //range_temp.end;
					ptr_newRange[depth][2].nrOfMismatches = range.nrOfMismatches + 1;
					//printf("new range 2 at depth %d: %d to %d\n", depth, hitLR[0], hitLR[1]);

				}

			}

		}

		//getchar();

		/*Update the matchIndicator if at end of recursion*/
		if (depth == nrOfEntries - 1) {
			for (i = 0; i < 3; i++) {
				for (j = ptr_newRange[depth][i].start; j <= ptr_newRange[depth][i].end; j++) {

					matchIndicator[j] = 1;

				}
				//printf("Recorded matches for a range %d from %d to %d\n", i, ptr_newRange[depth][i].start, ptr_newRange[depth][i].end); 
			}
			//getchar();
		}
		else {/*Recurse*/

			for (i = 0; i < 3; i++) {

				if (ptr_newRange[depth][i].start <= ptr_newRange[depth][i].end) {
					//printf("Recurse from range %d to depth:%d position: %d range starts: %d ends: %d\n", i, depth+1, position+1, ptr_newRange[depth][i].start, ptr_newRange[depth][i].end);
					//getchar();
					returnVal = bisectionalSearchSingle_Pairs(matchIndicator, depth + 1, ptr_newRange, binned_I_windowPairs_query, position + 1, ptr_newRange[depth][i], ptr_binned_I_windowPairs_DB, allowedNrOfMismatches, nrOfEntries);
				}
			}
		}

	}

	return returnVal;

}



/*For searching a word in an array of binned_I_counts*/
int wordSearch_I_binned_I_counts(int *binnedIvector, struct binned_I_counts *ptr_binned_I_counts, int start, int end, int nrOfEntries){

	int returnVal = 0;

	int left = start;
	int right = end;

	//int lastLeft = left;
	int lastRight = right;
	int mid = -117;

	int hit = -1; 

	/*if the search range does not contain the value, binned_I_windows_query.binnedIvector[position], there is no
	hit to be found, and hit = -1 is returned */
	if(lexicographic(binnedIvector, ptr_binned_I_counts[left].binnedIvector, nrOfEntries) == 1 ||
		lexicographic(binnedIvector, ptr_binned_I_counts[right].binnedIvector, nrOfEntries) == -1) {

			return hit; 

	}

	/*The following if-clause seems superfluous by the preceding if-clause and since the while-loop to follow
	is not entered when left = right; however the hit value must be updated from -1 to the index (=left)*/
	//If the range is actually a single value and we have a hit we return that// 
	if(left == right){
		if(lexicographic(binnedIvector, ptr_binned_I_counts[left].binnedIvector, nrOfEntries) == 0){
			hit = left;
		}
		//else: hit =-1 by the set initial value

		return hit;

	}
	
	//.. And else we search://
	while(mid != lastRight && left < right){//while(right -left > 1){ //while(mid  != lastRight && left <= right)

		//lastLeft = left;
		lastRight = right;

		mid = (int) ceil( (double) (right + left)/2 );

		//printf("mid:%d\n", mid);

		if( lexicographic(binnedIvector, ptr_binned_I_counts[mid].binnedIvector, nrOfEntries) == 1 ){

			right = mid;

		}
		else if( lexicographic(binnedIvector, ptr_binned_I_counts[mid].binnedIvector, nrOfEntries) == -1 ){

			left = mid;
		}
		else{ //lex == 0 so we have a hit "already" 
				hit = mid;
				//printf("break\n");
				break; //break while-loop
		}

		/*there remains the possibility that the hit sits at the left end point (added 7th March '18)*/
		if (lexicographic(binnedIvector, ptr_binned_I_counts[left].binnedIvector, nrOfEntries) == 0){
			hit = left;
		}

	}

	return hit;
}
