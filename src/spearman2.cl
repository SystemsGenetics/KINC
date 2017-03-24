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
