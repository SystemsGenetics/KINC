





/*!
 * Compute the next power of 2 which occurs after a number.
 *
 * @param n
 */
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
void swap(__global float *a, __global float *b)
{
   float c = *a;
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
void bitonicSort(__global float *array, int size)
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
            if ( (!dir && (array[a] > array[b])) || (dir && (array[a] < array[b])) )
            {
               swap(&array[a], &array[b]);
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
void bitonicSortFF(int size, __global float *array, __global float *extra)
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
            if ( (!dir && (array[a] > array[b])) || (dir && (array[a] < array[b])) )
            {
               swap(&array[a], &array[b]);
               swap(&extra[a], &extra[b]);
            }
         }
      }
   }
}






/*!
 * Compute the rank of a sorted vector in place. In the event of ties,
 * the ranks are corrected using fractional ranking.
 *
 * @param array
 * @param n
 */
void computeRank(__global float *array, int n)
{
   int i = 0;

   while ( i < n - 1 )
   {
      float a_i = array[i];

      if ( a_i == array[i + 1] )
      {
         int j = i + 2;
         int k;
         float rank = 0;

         // we have detected a tie, find number of equal elements
         while ( j < n && a_i == array[j] )
         {
            ++j;
         }

         // compute rank
         for ( k = i; k < j; ++k )
         {
            rank += k;
         }

         // divide by number of ties
         rank /= (float) (j - i);

         for ( k = i; k < j; ++k )
         {
            array[k] = rank;
         }

         i = j;
      }
      else
      {
         // no tie - set rank to natural ordered position
         array[i] = i;
         ++i;
      }
   }

   if ( i == n - 1 )
   {
      array[n - 1] = (float) (n - 1);
   }
}
