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
// for either gene.
//
// @param indexA Index into expression list for gene A.
// @param indexB Index into expression list for gene B.
// @param size Size of lists for both genes.
// @param listA New array of expressions for gene A that this function builds.
// @param listB New array of expressions for gene B that this function builds.
// @param expressions Array of all expressions for all genes to generate gene lists from.
// @return Returns size of newly generated arrays which excludes any missing expression values.
int fetchLists(int indexA, int indexB, int size, __global float* listA, __global float* listB
               , __global float* expressions)
{
   // initialize counters
   int i;
   int j = 0;
   int newSize;

   // go through expression list with given indexes, generating new lists from it
   for (i = 0; i < size ;++i)
   {
      if ( !isnan(expressions[indexA+i]) && !isnan(expressions[indexB+i]) )
      {
         // if both expressions exist add expressions to new lists and increment
         listA[j] = expressions[indexA+i];
         listB[j] = expressions[indexB+i];
         ++j;
      }
   }

   // set new size of generated lists and set unused end of lists to infinity
   newSize = j;
   for (i = j; i < size;++i)
   {
      listA[i] = INFINITY;
      listB[i] = INFINITY;
   }

   // return new size for generated lists
   return newSize;
}
