#include <stdio.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include <limits.h>
#include <math.h>
#include "sorting.h"

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1000000
#endif

#define MAXFILELEN   (50)   // length of string (path and filename)

void Save_seq1 (char *Filename, int N){
	FILE *fPtr; // File pointer.
	unsigned int element;
	int *seq1 = NULL;
	int p2 = 0, p3 = 0, u2=0, u3=0;
	int i=0, j=0;
  
	seq1 = (int *) malloc((N) * sizeof(int));
	seq1[0]= 1;
	for(i=1;i<=N;i++){
		if (seq1[p2] * 2 == seq1[i - 1])
			p2 += 1;
		if (seq1[p3] * 3 == seq1[i - 1])
			p3 += 1;
		u2 = seq1[p2] * 2;
		u3 = seq1[p3] * 3;
		if (u2 < u3){
			if (u2>=N) break;
			else{
				p2 += 1;
				seq1[i] = u2;
			}
		}
		else{
			if(u3>=N) break;
			else{
				p3 += 1;
				seq1[i] = u3;
			}
		}
	}
  
  fPtr = fopen(Filename, "w");  // Open file <file_name> with write mode.
  fprintf(fPtr, "%d\n", i);
  for (j=0; j<i; j++){
  	element=seq1[j];
  	fprintf(fPtr, "%d\n", element);
	}
  fclose(fPtr); // Close the output file.
  free(seq1); // Release the memory space of the data buffer
}
	
void Save_seq2 (char *Filename, int N){
	FILE *fPtr; // File pointer.
	unsigned int element;
	int *seq2 = NULL;
	int i=0, j=0, gap_n=0, gap=0;
	
  gap=N;
	while (gap!=1){
		gap_n++;
		gap=floor(gap/1.3);
		if (gap==9 || gap==10) gap=11;
	}
	seq2 = (int *) malloc((gap_n) * sizeof(int));
	gap=N;
	for (i=gap_n-1; i>=0; i--){
			gap=floor(gap/1.3);
			if (gap==9 || gap==10) gap=11;
			seq2[i]=gap;
	}
  
  fPtr = fopen(Filename, "w");  // Open file <file_name> with write mode.
  fprintf(fPtr, "%d\n", gap_n);
  for (j=0; j<gap_n; j++){
  	element=seq2[j];
  	fprintf(fPtr, "%d\n", element);
	}
  fclose(fPtr); // Close the output file.
  free(seq2); // Release the memory space of the data buffer
}

long *Load_File(char *Filename, int *Size){
	FILE *fPtr; // File pointer.
  long *buffer_ptr; // Pointer to the data buffer.
  int i=0, flag=0, fSize=0;
  long element;
  
  fPtr = fopen(Filename, "r"); 
  fscanf(fPtr,"%d", &fSize);//Read the first integer as the size of the data to be sorted.
  *Size = fSize;
  buffer_ptr = (long *) malloc((fSize)*sizeof(long)); // Allocate buffer memory.
  flag=0;
  i=0;
  while(flag==0){//Read other integers as data to be sorted.
  	fscanf(fPtr,"%d", &element);
  	if(i<fSize){
  		//element!='\n'&& element!='\r'&& element!=EOF && 
  		buffer_ptr[i] = element;
  		i++;
		}
		else if (element==EOF || i>=fSize)
			flag=1;
	}
  fclose(fPtr); // Close the input file.
  return buffer_ptr; // Return the memory address of the data buffer.
}

int Save_File(char *Filename, long *Array, int Size){
	FILE *fPtr; // File pointer.
	long element;
	int i=0;
  
  fPtr = fopen(Filename, "w");  // Open file <file_name> with write mode.
  fprintf(fPtr, "%d\n", Size);
  for (i=0; i<Size; i++){
  	element=Array[i];
  	fprintf(fPtr, "%ld\n", element);
	}
  fclose(fPtr); // Close the output file.
  return Size;
}

void Shell_Insertion_Sort(long *Array, int Size, double *N_Comp, double *N_Move){
	int *seq1 = NULL;
	double Ccount=0, Mcount=0;
	int p2 = 0, p3 = 0, u2=0, u3=0;
	int i=0, j=0, k=0, l=0, m=0, gap=0, fsize=0;
	long tmp=0;
	
	fsize=Size;
	seq1 = (int *) malloc((fsize) * sizeof(int));
	seq1[0]= 1;
	for(i=1;i<=fsize;i++){
		if (seq1[p2] * 2 == seq1[i - 1])
			p2 += 1;
		if (seq1[p3] * 3 == seq1[i - 1])
			p3 += 1;
		u2 = seq1[p2] * 2;
		u3 = seq1[p3] * 3;
		if (u2 < u3){
			if (u2>=fsize) break;
			else{
				p2 += 1;
				seq1[i] = u2;
			}
		}
		else{
			if(u3>=fsize) break;
			else{
				p3 += 1;
				seq1[i] = u3;
			}
		}
	}
	Ccount=0;
	Mcount=0;
	for(j=i-1; j>=0; j--){
		gap=seq1[j];
		for(k = 0; k < gap; k++){
      for(l = gap; l < fsize; l=l+gap){    
        for(m = l; m >= gap ; m -=gap){
        	Ccount++;
          if (Array[m] < Array[m -gap]){
          	tmp=Array[m -gap];
          	Mcount++;
          	Array[m -gap]= Array[m ];
          	Mcount++;
          	Array[m ] = tmp;
          	Mcount++;
					}
        }
    	}
    	if (l>=fsize) break;
		}
	}
	*N_Comp=Ccount;
	*N_Move=Mcount;
}

void Improved_Bubble_Sort(long *Array, int Size, double *N_Comp, double *N_Move){
	int *seq2 = NULL;
	double Ccount=0, Mcount=0;
	int i=0, j=0, k=0, l=0, m=0, gap_n=0, gap=0, fsize=0;
	long tmp=0;
	char Filename[MAXFILELEN] = "";
	
	fsize=Size;
	gap=fsize;
	while (gap!=1){
		gap_n++;
		gap=floor(gap/1.3);
		if (gap==9 || gap==10) gap=11;
	}
	seq2 = (int *) malloc((gap_n) * sizeof(int));
	gap=fsize;
	for (i=gap_n-1; i>=0; i--){
			gap=floor(gap/1.3);
			if (gap==9 || gap==10) gap=11;
			seq2[i]=gap;
	}
	Ccount=0;
	Mcount=0;
	for(j=gap_n-1; j>=0; j--){
		gap=seq2[j];
		for(k = 0; k < gap; k++){
      for(l = gap*(ceil(fsize/gap)-1); l>= gap; l=l-gap){  
        for(m = k; m < l ; m +=gap){
        	Ccount++;
          if (Array[m]>Array[m+gap]){
          	tmp=Array[m ];
          	Mcount++;
          	Array[m ] = Array[m +gap];
          	Mcount++;
          	Array[m +gap]=tmp;
          	Mcount++;
					}
        }
        if (l>=fsize) break;
    	}
    	if (k+gap>fsize) break;
		}
	}
	*N_Comp=Ccount;
	*N_Move=Mcount;
}
