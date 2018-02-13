





// Fetch and build array of expressions for both genes, skipping any expressions that are missing
// for either gene. Also builds ordered rank list used for spearman algorithm.
//
// @param expressions Array of all expressions for all genes to generate gene lists from.
// @param size Size of lists for both genes.
// @param workSize Size of new work arrays for genes and rank list.
// @param vector Index into expression list for gene A.
// @param sampleMask Array of values denoting membership in a cluster.
// @param a New array of expressions for gene A that this function builds.
// @param b New array of expressions for gene B that this function builds.
// @param rankList New array that is initialized to start at 1 and increment by one for each
// successive element.
// @return Returns size of newly generated arrays which excludes any missing expression values.
int fetchData(
   __global const float *expressions,
   int size, int workSize,
   int2 vector,
   __global const char *sampleMask,
   __global float *a,
   __global float *b,
   __global int *rankList)
{
   // initialize counters and indexes
   __global const float *gene1 = &expressions[vector.x * size];
   __global const float *gene2 = &expressions[vector.y * size];

   // populate a and b with shared expressions of gene pair
   int j = 0;

   if ( sampleMask[0] != -1 )
   {
      // add samples that are in the cluster
      for ( int i = 0; i < size; ++i )
      {
         if ( sampleMask[i] == 1 )
         {
            a[j] = gene1[i];
            b[j] = gene2[i];
            rankList[j] = j + 1;
            ++j;
         }
      }
   }
   else
   {
      // add samples that are valid
      for ( int i = 0; i < size; ++i )
      {
         if ( !isnan(gene1[i]) && !isnan(gene2[i]) )
         {
            a[j] = gene1[i];
            b[j] = gene2[i];
            rankList[j] = j + 1;
            ++j;
         }
      }
   }

   // set any remaining values in work arrays to infinity or zero
   int numSamples = j;

   for ( int i = j; i < workSize; ++i )
   {
      a[i] = INFINITY;
      b[i] = INFINITY;
      rankList[i] = 0;
   }

   // return new size for generated lists
   return numSamples;
}






// Swap two floating point values.
//
// @param a First floating point value.
// @param b Second floating point value.
void swapF(__global float* a, __global float* b)
{
   float c = *a;
   *a = *b;
   *b = c;
}






// Sort a given list using the bitonic algorithm along with rearranging a second list with the same
// operations that are done to sort the sorted list.
//
// @param size Size of lists. It MUST be a power of 2.
// @param sortList The list to be sorted.
// @param extraList The extra list that will be rearranged identical to the sorted list.
void bitonicSortFF(int size, __global float* sortList, __global float* extraList)
{
   // initialize all variables
   int bsize = size/2;
   int ob,ib,i,dir,a,b,t;

   // bitonic algorithm, starting with an outer block of 2 and working up to total size of list
   for (ob = 2; ob <= size ;ob *= 2)
   {
      for (ib = ob; ib >= 2 ;ib /= 2)
      {
         t = ib/2;
         for (i = 0; i < bsize ;++i)
         {
            dir = -((i/(ob/2))&0x1);
            a = (i/t)*ib+(i%t);
            b = a+t;
            if ( ( ( sortList[a] > sortList[b] ) && !dir )
                 || ( ( sortList[a] < sortList[b] ) && dir ) )
            {
               swapF(&sortList[a],&sortList[b]);
               swapF(&extraList[a],&extraList[b]);
            }
         }
      }
   }
}






// Swap two integer values.
//
// @param a First integer value.
// @param b Second integer value.
void swapI(__global int* a, __global int* b)
{
   int c = *a;
   *a = *b;
   *b = c;
}






// Sort a given list using the bitonic algorithm along with rearranging a second list with the same
// operations that are done to sort the sorted list.
//
// @param size Size of lists. It MUST be a power of 2.
// @param sortList The list to be sorted.
// @param extraList The extra list that will be rearranged identical to the sorted list.
void bitonicSortFI(int size, __global float* sortList, __global int* extraList)
{
   // initialize all variables
   int bsize = size/2;
   int ob,ib,i,dir,a,b,t;

   // bitonic algorithm, starting with an outer block of 2 and working up to total size of list
   for (ob = 2; ob <= size ;ob *= 2)
   {
      for (ib = ob; ib >= 2 ;ib /= 2)
      {
         for (i = 0; i < bsize ;++i)
         {
            dir = -((i/(ob/2))&0x1);
            t = ib/2;
            a = (i/t)*ib+(i%t);
            b = a+t;
            if ( ( ( sortList[a] > sortList[b] ) && !dir )
                 || ( ( sortList[a] < sortList[b] ) && dir ) )
            {
               swapF(&sortList[a],&sortList[b]);
               swapI(&extraList[a],&extraList[b]);
            }
         }
      }
   }
}






// Make the final calculation of the spearman coefficient with the given presorted spearman ranking
// list.
//
// @param size Size of the ranking list.
// @param rankList Presorted spearman ranking list.
// @return Returns floating point spearman coefficient.
float calculateSpearman(int size, __global int* rankList)
{
   // declare and initialize all variables
   int i;
   long tmp;
   long difference = 0;

   // go through spearman sorted rank list and calculate difference from 1,2,3,... list
   for (i = 0; i < size ;++i)
   {
      tmp = (i+1)-rankList[i];
      difference += tmp*tmp;
   }

   // calculate and return spearman coefficient
   return 1.0-(6.0*(float)difference/((float)size*(((float)size*(float)size)-1.0)));
}






// Calculate a block of spearman coefficients given a block of gene pairs.
//
// @param expressions Row first 2 dimensional array of all gene expressions.
// @param size The size of the expressions/samples per gene.
// @param workSize The power of 2 work size, MUST be a power of 2.
// @param pairs Array of gene pairs.
// @param sampleMasks Array of sample masks for each gene pair.
// @param workLists Work space to be used for spearman calculations.
// @param rankLists Work space to be used for spearman calculations.
// @param results Array of output spearman coefficients for each gene pair.
__kernel void calculateSpearmanBlock(
   __global const float *expressions,
   int size, int workSize,
   __global const int2 *pairs,
   __global const char *sampleMasks,
   int minSamples,
   __global float *workLists,
   __global int *rankLists,
   __global float *results)
{
   int i = get_global_id(0);

   if ( pairs[i].x == 0 && pairs[i].y == 0 )
   {
      return;
   }

   // initialize workspace variables
   __global const char *sampleMask = &sampleMasks[i * size];
   __global float *a = &workLists[(2*i+0) * workSize];
   __global float *b = &workLists[(2*i+1) * workSize];
   __global int *rankList = &rankLists[2*i * workSize];

   // fetch a and b arrays from expression matrix
   int numSamples = fetchData(
      expressions, size, workSize,
      pairs[i],
      sampleMask,
      a, b, rankList
   );

   // compute correlation only if there are enough samples
   results[i] = NAN;

   if ( numSamples >= minSamples )
   {
      // get new power of 2 floor size
      int pow2Size = 2;
      while ( pow2Size < numSamples )
      {
         pow2Size *= 2;
      }

      // execute two bitonic sorts that is beginning of spearman algorithm
      bitonicSortFF(pow2Size, a, b);
      bitonicSortFI(pow2Size, b, rankList);

      // calculate spearman coefficient from rearranged rank list and save to result list
      results[i] = calculateSpearman(numSamples, rankList);
   }
}
