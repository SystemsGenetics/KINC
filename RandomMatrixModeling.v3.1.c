// RandomMatrixModeling.v3.1.c
// uses RMT method to determine proper threshold for correlation matrix
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include "RandomMatrix.h"
#include <time.h>

/*============GLOBAL VARIABLES =================================*/
int* usedFlag; //will hold an array flagging use of indicies in the matrix
int* index1; //had to rename because index was already taken
int numGenes;
char* inputFileName;
char* inputDir;
int numLinesPerFile;
RMTParameters rmtParameter;
FILE* timing;
// FILE* runInfo;
int verbose; //a boolean for whether or not to print timing data

//set some global params for the RMM stepping
float thresholdStart; //the threshold to start generating cutMatrices with (actually, +thresholdStep)
float thresholdStep;  //the threshold step for each iteration of generation
float chiSoughtValue; //the chiValue the while loop will end on

/*==============FUNCTION DECLARATIONS BEGIN!====================*/

//writes a value to a string, then returns the pointer.  Only supports base 10 numbers.
char* itoa(int val, char* ptr){
	sprintf(ptr,"%d",val);
	return ptr;
} 

FILE* fileOpenHelper(char* extension){
	char* FNbuffer;//filename buffer
	FILE* fileObject; 
	int len;
	len = strlen(inputFileName);
	len+=strlen(extension);
	len++;
	FNbuffer = (char*)malloc(sizeof(char)*len);
	memset(FNbuffer,'\0', len);
	strcat(FNbuffer,inputFileName);
	strcat(FNbuffer,extension);
	fileObject = fopen(FNbuffer, "w");
	if (fileObject==NULL){
		printf("\nUnable to create file '%s' . Please check input filename for validity.\nExiting.\n", FNbuffer);
		exit(0);
	}
	free(FNbuffer);
	return fileObject;
}

int findNumUsed(){
	int i,sum=0;
	for(i=0;i<numGenes;i++){
		if(usedFlag[i]!=0) sum++;
	}
	return sum;
}

void setIndexArray(){
	int i,j=0;
	for(i=0;i<numGenes;i++){
		if(usedFlag[i]!=0){
			index1[i]=j;
			j++;
		}
	}
	return;
}

float* readPearsonCorrelationMatrix(float th, int *size){
	memset(index1, -1, sizeof(int)*(numGenes));
	float* rowj = (float*) malloc(sizeof(float)*numGenes);
	char* filename;
	char num[4];
	int len = strlen(inputFileName);
	len+=strlen(inputDir);
	len+=strlen(" Pearson Correlation ");
	len+=strlen("xxx.bin");
	len++;
	filename = (char*)malloc(sizeof(char)*len);
	memset(filename, '\0', len);
	int i,j, used, z, k, h, limit, junk;
	FILE* in;
	z = (numGenes-1)/numLinesPerFile;
	for(i=0;i<=z;i++){
		memset(filename, '\0', len);
		strcat(filename, inputDir);
		strcat(filename, inputFileName);
		strcat(filename, " Pearson Correlation ");
		strcat(filename, itoa(i, num));
		strcat(filename, ".bin");
		in = fopen(filename, "rb");
		fread(&junk, sizeof(int), 1, in);//numGenes
		fread(&junk, sizeof(int), 1, in);//numLinesPerFile
		if(i!=z){limit = (i+1)*numLinesPerFile;}
		else{limit = numGenes;}
		for(j = i*numLinesPerFile ; j< limit ;j++){
			fread(rowj, sizeof(float), j+1, in);
			for(k=0;k<j+1;k++){
				if(k!=j && fabs(rowj[k])>th){
					usedFlag[k]=1;
					usedFlag[j]=1;
				}
			}
		}	
		fclose(in);
	}
	
	used = findNumUsed();
	setIndexArray();
	//allocate memory for new cut matrix
	float* cutM = (float*)calloc(used*used, sizeof(float));	
	for(i=0;i<used;i++){	
		cutM[i+i*used]=1;  //initialize the diagonal to 1 
	}

	//tell determinePearson the size of the new matrix
	*size = used;

	//copy all eligible values to cut matrix
	//N.B. In this loop rowj[] actually corresponds to row i
	for(h=0;h<z;h++){
		memset(filename,'\0', len );
		strcat(filename, inputDir);
		strcat(filename, inputFileName);
		strcat(filename, " Pearson Correlation ");
		strcat(filename, itoa(h, num));
		strcat(filename, ".bin");
		in = fopen(filename, "rb");
		fread(&junk, sizeof(int), 1, in);//numGenes
		fread(&junk, sizeof(int), 1, in);//numLinesPerFile
		if(i!=z){limit = (h+1)*numLinesPerFile;}
		else{limit = numGenes;}
		for(i=h*numLinesPerFile ; i<limit;i++){
			fread(rowj, sizeof(float), i+1, in);//actually rowi in this loop
			for(j=0;j<i+1;j++){
				if(i!=j && fabs(rowj[j])>th){
					cutM[index1[j] + used*index1[i]] = rowj[j];//actually rowi in this loop
				}
			}
		}
		fclose(in);
	}

	free(rowj);
	free(filename);
	return cutM;
}
	
int determinePearsonCorrelationThreshold_LargeMatrix(){
	float* newM;
	int size;
	time_t start, end;
	// FILE* eigenF, *chiF;
	// eigenF = fileOpenHelper(".eigens.txt");
	// chiF = fileOpenHelper(".chiVals.txt");

	float th = thresholdStart;
	double finalTH = 0.0;
	double finalChi = 10000.0;
	double minTH = 1.0;
	double minChi = 10000.0;
	double maxTH = 0.0;
	double maxChi = 0.0;
	float* E;  //array for eigenvalues
	int i = 0;
	double chi;

	do{
		th = th - thresholdStep;
//		time(&start);
		newM = readPearsonCorrelationMatrix(th, &size);
//		time(&end);
//		if(verbose){
//			fprintf(timing,"Minutes for cut matrix: %f\n", (end-start)/60.0);
//		}
//		fprintf(runInfo,"%f\t%d\n",th,size);
//		fflush(runInfo);
		if(size >= 100){
//			time(&start);
			E = calculateEigen(newM, size);
//			time(&end);

			free(newM);
//			if(verbose){
//				fprintf(timing,"Minutes calculating Eigens: %f\n", (end-start)/60.0);
//			}
			/* print out eigenvalues to file */
			// fprintf(eigenF, "%f\t", th);
			// for(i=0 ; i<size ; i++){
				// fprintf(eigenF, "%f\t", E[i]);
			// }
			// fprintf(eigenF,"\n");
//			time(&start);
			chi = chiSquareTestUnfoldingNNSDWithPoisson(E, size, \
				rmtParameter.nnsdHistogramBin, \
				rmtParameter.minUnfoldingPace, \
				rmtParameter.maxUnfoldingPace);
//			time(&end);
//			if(verbose){
//				fprintf(timing,"Minutes Chi Square Testing: %f\n", (end-start)/60.0);
//			}
//			fflush(timing);
			// fprintf(chiF, "%f\t%f\n", th, chi);
			// fflush(chiF);
			free(E);
			if(chi < minChi){
				minChi = chi;
				minTH = th;
			}
			if (chi < rmtParameter.chiSquareTestThreshold){
				finalTH = th;
				finalChi = chi;
			}
			if (finalChi < rmtParameter.chiSquareTestThreshold && chi > finalChi && th < finalTH){
				maxChi = chi;
				maxTH = th;
			}
		}
		else{
			free(newM);
		}
	}while(maxChi < chiSoughtValue || size == numGenes);
	
	/*
	If finalChi is still greater than threshold, check the small scale
	*/

	if (finalChi > rmtParameter.chiSquareTestThreshold){
	// fprintf(runInfo,"checking small scale\n");
	// fprintf(chiF, "checking small scale\n");
		th = (float)minTH + 0.2;
		for(i = 0 ; i <= 40 ; i++){
			th = th - thresholdStep*i;
//			time(&start);
			newM = readPearsonCorrelationMatrix(th, &size);
//			time(&end);
//			if(verbose){
//				fprintf(timing,"Minutes for cut matrix: %f\n", (end-start)/60.0);
//			}
			// fprintf(runInfo,"%f\t%d\n",th,size);			
			if(size >= 100){
//				time(&start);
				E = calculateEigen(newM, size);
//				time(&end);
//				if(verbose){
//					fprintf(timing,"Minutes for cut matrix: %f\n", (end-start)/60.0);
//				}
				free(newM);
				/* print out eigenvalues to file */
				// fprintf(eigenF, "%f\t", th);
				// for(i=0 ; i<size ; i++){
					// fprintf(eigenF, "%f\t", E[i]);
				// }
				// fprintf(eigenF,"\n");
//				time(&start);
				chi = chiSquareTestUnfoldingNNSDWithPoisson(E, size, \
					rmtParameter.nnsdHistogramBin, \
					rmtParameter.minUnfoldingPace, \
					rmtParameter.maxUnfoldingPace);
//				time(&end);
//				if(verbose){
//					fprintf(timing,"Minutes Chi Square Testing: %f\n", (end-start)/60.0);
//				}
				// fprintf(chiF, "%f\t%f\n", th, chi);  //not done in this loop in original code
				free(E);
				if(chi < minChi){
					minChi = chi;
					minTH = th;
				}
				if (chi < rmtParameter.chiSquareTestThreshold){
					finalTH = th;
					finalChi = chi;
				}
			}//end if size >= 100
			else{
				free(newM);
			}
		}//end for 1 -> 40 loop
	}//end if finalChi > rmt... 

	//close the chi and eigen files now that results are written
	// fclose(chiF);
	// fclose(eigenF);
	/*
	Set the Properties file according to success or failure
	*/		
	// fprintf(runInfo, "==================================\n");
	if(finalChi < rmtParameter.chiSquareTestThreshold){
		finalTH = ceil(finalTH*10000)/10000.0;
		// fprintf(runInfo,"Step SUCCEEDED\n");
		// fprintf(runInfo, "Pearson Correlation Threshold = %f\n", finalTH);
		//fprintf(runInfo, "final chi val = %f\n", finalChi);
		FILE* th;
		th = fileOpenHelper("th.txt");
		fprintf(th, "%f", finalTH);
		fclose(th);
		return 0;
	}
	else{
		finalTH = ceil(finalTH*10000)/10000.0;
		// fprintf(runInfo, "Step FAILED.\n");
		// fprintf(runInfo, "Pearson Correlation Threshold = %f\n", finalTH);
		//fprintf(runInfo, "final chi val = %f\n", finalChi);
		return -2;
	}
}

int main(int argc, char** argv){
	int i, size;

	if(argc==1){
		printf("The arguments for this fucntion are:\n\n\
'-i':\tThe input file name. Same as used in previous step.\n\
\tMust be the same as the name used in the matrix binary\n\
\tfiles.  This name will also be used to create files \n\
\tassociated with this run of the program.\n\n\
Optional:\n\n\
'-b':\tThe initial threshold(+1*step) value that will be used.\n\
\t [b=Begin value] Default: 0.9200\n\
'-s':\tThe threshold step size used each iteration. Default: 0.001\n\
'-c':\tThe chi-square test value that the loop will stop on.\n\
\tDefault: 200\n\n\
'-v':\tSet the performance collection. Has two values possible values,\n\
ON/OFF . [v=Verbose] Default: ON\n\
Examples:\n\
<executable> -i <input.file.name> \n\
<exec> -i <input.file.name> -s <0.0001> -v ON\n\n\
Note: Order is not important, but spaces are required between the \n\
flag and its value.  Also note that erroneous inputs could cause the\n\
program to crash without a helpful error message...so get it right.\n\n\
			Exiting.\n");
		exit(0);
	}
	
	inputDir = "Pearson/";
	//initialize default values, which may be overwritten in the command
	verbose=1;
	thresholdStart = 0.92; 
	thresholdStep = 0.001;  
	chiSoughtValue = 200; 
	
	//parse the command line flags
	for(i=1;i<argc;i+=2){
		if(argv[i][1]=='i'){
			inputFileName = argv[i+1];
		}
		else if(argv[i][1]=='b'){
			thresholdStart = atof(argv[i+1]);
		}
		else if(argv[i][1]=='s'){
			thresholdStep = atof(argv[i+1]);
		}
		else if(argv[i][1]=='c'){
			chiSoughtValue = atof(argv[i+1]);
		}
		else if(argv[i][1]=='v'){
			if(strcmp(argv[i+1],"ON")==0 || strcmp(argv[i+1],"on")==0 || strcmp(argv[i+1],"On")==0 || strcmp(argv[i+1],"oN")==0){
				verbose = 1;
			}
			else{
				verbose = 0;
			}
		}
		else{
			printf("Flag '%s' not recognized, exiting.", argv[i]);
			exit(0);
		}
	}

	//read from first matrix file to determine some stuff
	char* filename;
	char num[4];
	int len = strlen(inputFileName);
	len+=strlen(inputDir);
	len+=strlen(" Pearson Correlation ");
	len+=strlen("xxx.bin");
	len++;
	filename = (char*)malloc(sizeof(char)*len);
	memset(filename, '\0', len);
	strcat(filename, inputDir);
	strcat(filename, inputFileName);
	strcat(filename, " Pearson Correlation ");
	strcat(filename, itoa(0, num));
	strcat(filename, ".bin");
	FILE* info;
	info = fopen(filename, "rb");
	fread(&numGenes, sizeof(int), 1, info);
	fread(&numLinesPerFile, sizeof(int), 1, info);
	fclose(info);
	free(filename);

	//initialize some global file pointers that will be used throughout	
	// runInfo = fileOpenHelper(".runInfo.txt");
	
	if(verbose == 1){
		timing = fileOpenHelper(".timing.txt");
	}

	//print some preliminary run info to the runInfo file
	// fprintf(runInfo, "number of genes: %d\n", numGenes);
	// fprintf(runInfo, "threshold start: %f\n", thresholdStart);
	// fprintf(runInfo, "threshold step:  %f\n", thresholdStep);
	// fprintf(runInfo, "sought chi val:  %f\n\n", chiSoughtValue);
	// fprintf(runInfo, "Below are the thresholds and corresponding cut matrix sizes\n");
	
	//initialize the global RMTParameters struct
	rmtParameter.nnsdHistogramBin = 0.05;
	rmtParameter.chiSquareTestThreshold = 99.607;
	rmtParameter.minUnfoldingPace=10;
	rmtParameter.maxUnfoldingPace=41;
	rmtParameter.mimiModuleSize = 4;
	rmtParameter.edHistogramBin = 0.1;

	//allocate memory for usedFlags
	usedFlag = (int*) malloc(numGenes*sizeof(int));

	//allocate memory for index 
	index1 = (int*) malloc(numGenes*sizeof(int));

	time_t start, end;
	fflush(timing);
	time(&start);
	int rc = determinePearsonCorrelationThreshold_LargeMatrix();
	time(&end);
	if(verbose){
		fprintf(timing,"Minutes in determining: %f\n", (end-start)/60.0);
	}
	free(index1);
	free(usedFlag);
	// fclose(runInfo);
	if(verbose){
		fclose(timing);
	}
	exit(rc);
}
