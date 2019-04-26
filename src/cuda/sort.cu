





/*!
 * Compute the next power of 2 which occurs after a number.
 *
 * @param n
 */
__device__
int nextPower2(int n)
{
   int pow2 = 2;
   while ( pow2 < n )
   {
      pow2 *= 2;
   }

   return pow2;
}






/*!
 * Swap two values
 *
 * @param a
 * @param b
 */
__device__
void swapF(float *a, float *b)
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
__device__
void swapI(int *a, int *b)
{
   int c = *a;
   *a = *b;
   *b = c;
}






/*!
 * Sort an array using bitonic sort. The array should have a size which is a
 * power of two.
 *
 * @param array
 * @param size
 */
__device__
void bitonicSort(float *array, int size)
{
   int bsize = size / 2;
   int dir, a, b, t;

   for ( int ob = 2; ob <= size; ob *= 2 )
   {
      for ( int ib = ob; ib >= 2; ib /= 2 )
      {
         t = ib/2;
         for ( int i = 0; i < bsize; ++i )
         {
            dir = -((i/(ob/2)) & 0x1);
            a = (i/t) * ib + (i%t);
            b = a + t;
            if ( ((array[a] > array[b]) && !dir) || ((array[a] < array[b]) && dir) )
            {
               swapF(&array[a], &array[b]);
            }
         }
      }
   }
}






/*!
 * Sort an array using bitonic sort, while also applying the same swap operations
 * to a second array of the same size. The arrays should have a size which is a
 * power of two.
 *
 * @param size
 * @param array
 * @param extra
 */
__device__
void bitonicSortFF(int size, float *array, float *extra)
{
   int bsize = size / 2;
   int dir, a, b, t;

   for ( int ob = 2; ob <= size; ob *= 2 )
   {
      for ( int ib = ob; ib >= 2; ib /= 2 )
      {
         t = ib/2;
         for ( int i = 0; i < bsize; ++i )
         {
            dir = -((i/(ob/2)) & 0x1);
            a = (i/t) * ib + (i%t);
            b = a + t;
            if ( ((array[a] > array[b]) && !dir) || ((array[a] < array[b]) && dir) )
            {
               swapF(&array[a], &array[b]);
               swapF(&extra[a], &extra[b]);
            }
         }
      }
   }
}






/*!
 * Sort an array using bitonic sort, while also applying the same swap operations
 * to a second array of the same size. The arrays should have a size which is a
 * power of two.
 *
 * @param size
 * @param array
 * @param extra
 */
__device__
void bitonicSortFI(int size, float *array, int *extra)
{
   int bsize = size / 2;
   int dir, a, b, t;

   for ( int ob = 2; ob <= size; ob *= 2 )
   {
      for ( int ib = ob; ib >= 2; ib /= 2 )
      {
         t = ib/2;
         for ( int i = 0; i < bsize; ++i )
         {
            dir = -((i/(ob/2)) & 0x1);
            a = (i/t) * ib + (i%t);
            b = a + t;
            if ( ((array[a] > array[b]) && !dir) || ((array[a] < array[b]) && dir) )
            {
               swapF(&array[a], &array[b]);
               swapI(&extra[a], &extra[b]);
            }
         }
      }
   }
}
