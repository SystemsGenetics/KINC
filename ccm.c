// Calculate Correlation Matrix ccm.c
// calculates a pearson correlation coefficient matrix based on input data
// Scott Gibson
// October 2011

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ROWS_PER_OUTPUT_FILE 10000
#define COUNT_POINTS 101

int main(int argc, char *argv[])
{
	if(argc != 4 && argc != 5)
	{
		printf("Usage: ./ccm <input file (without .txt)> <rows> <cols> <turn off perf monitor with '0'>\n");
		exit(-1);
	}	

	FILE *infile, *outfile;

	char string[50], outfilename[50], infilename[50];

	int i, j, k, m, rows = atoi(argv[2]), cols = atoi(argv[3]), perf = 1, corrdata = 0;

	float **data, r, one = 1.0;

	double *sum, *sumsquared, xy, x, y, x2, y2;

	time_t start_time, end_time;

	// if all arguments are present and the final one is '0', disable performance monitoring
	if(argc == 5 && atoi(argv[4]) == 0)
	{
		perf = 0;
	}

	if(perf) time(&start_time);

	// add ".txt" to the input file name
	sprintf(infilename, "%s.txt", argv[1]);

	// space for data to be read in
	data = (float**) malloc(sizeof(float *) * rows);
	for(i = 0; i < rows; i++)
	{
		data[i] = malloc(sizeof(float) * cols);
	}

	// sum of elements and sum of squares of elements for each row
	sum = (double*) malloc(sizeof(double) * rows);
	sumsquared = (double*) malloc(sizeof(double) * rows);

	// read in data from user-specified file
	infile = fopen(infilename, "r");
	for(i = 0; i < rows; i++)
	{	// the first entry on every line is a label string - read that in before the numerical data
		fscanf(infile, "%s", string);
		for(j = 0; j < cols; j++)
		{
			if(fscanf(infile, "%f", &data[i][j]) == EOF)
			{
				printf("Error: EOF reached early. Exiting.\n");
				exit(-1);
			}
		}
	}

	// calculate the sum and sumsquared of each row's data (to be used in the Pearson correlation formula)
	for(i = 0; i < rows; i++)
	{
		sum[i] = sumsquared[i] = 0;

		for(j = 0; j < cols; j++)
		{
			sum[i] = sum[i] + data[i][j];
			sumsquared[i] = sumsquared[i] + data[i][j] * data[i][j];
		}
	}

	// output a maximum of ROWS_PER_OUTPUT_FILE rows and then start a new file
	int z = (rows-1) / ROWS_PER_OUTPUT_FILE, j_limit;
	
	// keep track of the distribution of coefficients 
	// (if 'corrdata' is not FALSE, this option is used for debugging to check coefficient distribution)
	int count[COUNT_POINTS];
	
	if(corrdata) for(m = 0; m < COUNT_POINTS; m++) count[m] = 0;

	// each iteration of m is a new output file
	for(m = 0; m <= z; m++)
	{
		//printf("m = %d / %d\n", m, z);

		// the output file will be located in the Pearson directory and named based on the input file info
		sprintf(outfilename, "./Pearson/%s Pearson Correlation %d.bin", argv[1], m);

		outfile = fopen(outfilename, "wb");

		// calculate the limit on the rows to output based on where we are in the calculation
		if(m < z) j_limit = (m + 1) * ROWS_PER_OUTPUT_FILE;
		else j_limit = rows;

		// output the size of the overall matrix for the next step to use
		fwrite(&rows, sizeof(rows), 1, outfile);

		// determine and output the number of rows that will be stored in the current file
		i = j_limit - m * ROWS_PER_OUTPUT_FILE;
		fwrite(&i, sizeof(i), 1, outfile);

		for(j = m * ROWS_PER_OUTPUT_FILE; j < j_limit; j++)
		{	// lower triangular symmetric matrix (stop when col# > row#) 
			for(k = 0; k <= j; k++)
			{
				if(j == k)
				{	// correlation of an element with itself is 1
					fwrite(&one, sizeof(one), 1, outfile);
				}
				else
				{
					xy = 0;

					// calculate sum of XY for this pairwise correlation
					for(i = 0; i < cols; i++)
					{
						xy = xy + data[j][i] * data[k][i];
					}

					// Pearson's correlation formula
					r = (xy - sum[j] * sum[k] / cols) / sqrt( (sumsquared[j] - sum[j] * sum[j] / cols) * (sumsquared[k] - sum[k] * sum[k] / cols) );

					fwrite(&r, sizeof(r), 1, outfile);
					
					if(corrdata)
					{
						if(r < 0) r = -r;
						count[(int)(r * (COUNT_POINTS-1))]++;
					}
				}
			}
		}
		fclose(outfile);
	}
	
	if(corrdata)
	{
		outfile = fopen("corrdata.txt", "w");
		
		for(m = 0; m < COUNT_POINTS; m++)
		{
			fprintf(outfile, "%lf\t%d\n", 1.0 * m / (COUNT_POINTS-1), count[m]);
		}
		
		fclose(outfile);
	}

	if(perf)
	{
		time(&end_time);

		outfile = fopen("timingdata.txt", "a");
		
		fprintf(outfile, "CCM Runtime with %d x %d %s input dataset: %.2lf min\n", rows, cols, infilename, difftime(end_time, start_time)/60.0);
	}

	return 0;
}
