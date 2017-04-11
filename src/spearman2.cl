/* BASIC PLAN FOR NEW SPEARMAN KERNEL..

Prune lists... ignoring infinities
IF pruned list is big enough... (else just make result NaN)
use bitonic sort with for loop to do 2 required sorts for spearman
make summation of rank differences
compute the spearman coefficient


MEMORY REQUIREMENTS...
will need TWO float arrays of size n
one int array of size n
where n is the number of samples per gene
in turn, this will be multilpied by the total number of threads
so size = ( 2*sizeof(float) + sizeof(int) )*thread_size

*/






// Fetch and build array of expressions for both genes, skipping any expressions that are missing
// for either gene. Also builds ordered rank list used for spearman algorithm.
//
// @param indexA Index into expression list for gene A.
// @param indexB Index into expression list for gene B.
// @param size Size of lists for both genes.
// @param workSize Size of new work arrays for genes and rank list.
// @param listA New array of expressions for gene A that this function builds.
// @param listB New array of expressions for gene B that this function builds.
// @param rankList New array that is initialized to start at 1 and increment by one for each
// successive element.
// @param expressions Array of all expressions for all genes to generate gene lists from.
// @return Returns size of newly generated arrays which excludes any missing expression values.
int fetchLists(int indexA, int indexB, int size, int workSize, __global float* listA
               , __global float* listB, __global int* rankList, __global float* expressions)
{
   // initialize counters
   int i;
   int j = 0;
   int newSize = 2;

   // go through expression list with given indexes, generating new lists from it
   for (i = 0; i < size ;++i)
   {
      if ( !isnan(expressions[indexA+i]) && !isnan(expressions[indexB+i]) )
      {
         // if both expressions exist add expressions to new lists, next rank, and increment
         listA[j] = expressions[indexA+i];
         listB[j] = expressions[indexB+i];
         rankList[j] = j+1;
         j++;
      }
   }

   // set new size of generated lists and set unused end of lists to infinity or zero
   newSize = j;
   for (i = j; i < size;++i)
   {
      listA[i] = INFINITY;
      listB[i] = INFINITY;
      rankList[i] = 0;
   }

   // set any remaining values in work arrays to infinity or zero
   for ( i = size; i < workSize ;++i)
   {
      listA[i] = INFINITY;
      listB[i] = INFINITY;
      rankList[i] = 0;
   }

   // return new size for generated lists
   return newSize;
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






// Calculate a bulk selection of spearman coefficients with different genes.
//
// @param size The size of the expressions/samples per gene.
// @param workSize The power of 2 work size, MUST be a power of 2.
// @param targetList Array containing all gene targets to sort.
// @param expressions Row first 2 dimensional array of all gene expressions.
// @param workLists Work space to be used for spearman calculations.
// @param rankLists Work space to be used for spearman calculations.
// @param resultList Array that will contain all completed spearman coefficients of all gene
// correlations specified in target list.
__kernel void calculateSpearmanBlock(int size, int workSize, __global int* targetList
                                     , __global float* expressions, __global float* workLists
                                     , __global int* rankLists, __global float* resultList)
{
   // initialize all variables and get global id
   int newSize,pow2Size;
   int i = get_global_id(0);
   __global float* listA = &workLists[i*workSize];
   __global float* listB = &workLists[(i+1)*workSize];
   __global int* rankList = &rankLists[i*workSize];

   // fetch gene expressions lists for both genes from target list
   newSize = fetchLists(targetList[i],targetList[i+1],size,workSize,listA,listB,rankList
                        ,expressions);

   // get new power of 2 floor size
   pow2Size = 2;
   while ( pow2Size < newSize )
   {
      pow2Size *= 2;
   }

   // execute two bitonic sorts that is beginning of spearman algorithm
   bitonicSortFF(pow2Size,listA,listB);
   bitonicSortFI(pow2Size,listB,rankList);

   // calculate spearman coefficient from rearranged rank list and save to result list
   resultList[i] = calculateSpearman(newSize,rankList);
}
