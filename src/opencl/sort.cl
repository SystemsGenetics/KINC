





/*!
 * Swap two values
 *
 * @param a
 * @param b
 */
void swapF(__global float* a, __global float* b)
{
   float c = *a;
   *a = *b;
   *b = c;
}






/*!
 * Swap two values
 *
 * @param a
 * @param b
 */
void swapI(__global int* a, __global int* b)
{
   int c = *a;
   *a = *b;
   *b = c;
}






void siftDown(__global float *array, int start, int end)
{
   int root = start;

   while ( 2 * root + 1 <= end )
   {
      int child = 2 * root + 1;
      int swp = root;

      if ( array[swp] < array[child] )
      {
         swp = child;
      }

      if ( child + 1 <= end && array[swp] < array[child + 1] )
      {
         swp = child + 1;
      }

      if ( swp == root )
      {
         return;
      }
      else
      {
         swapF(&array[root], &array[swp]);
         root = swp;
      }
   }
}






void heapify(__global float *array, int n)
{
   int start = ((n-1) - 1) / 2;

   while ( start >= 0 )
   {
      siftDown(array, start, n - 1);
      start -= 1;
   }
}






/*!
 * Sort an array using heapsort.
 *
 * @param array
 * @param n
 */
void heapSort(__global float *array, int n)
{
   heapify(array, n);

   int end = n - 1;
   while ( end > 0 )
   {
      swapF(&array[end], &array[0]);
      end -= 1;

      siftDown(array, 0, end);
   }
}






/*!
 * Sort a list using bitonic sort, while also applying the same swap operations
 * to a second list of the same size. The lists should have a size which is a
 * power of two.
 *
 * @param size
 * @param sortList
 * @param extraList
 */
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






/*!
 * Sort a list using bitonic sort, while also applying the same swap operations
 * to a second list of the same size. The lists should have a size which is a
 * power of two.
 *
 * @param size
 * @param sortList
 * @param extraList
 */
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
               swapI(&extraList[a],&extraList[b]);
            }
         }
      }
   }
}
