#include "pairwise_spearman.h"



using namespace Pairwise;






/*!
 * Construct a Spearman correlation model.
 *
 * @param emx
 */
Spearman::Spearman(ExpressionMatrix* emx)
{
   // pre-allocate workspace
   _x_rank.resize(emx->sampleSize());
   _y_rank.resize(emx->sampleSize());
}






/*!
 * Compute the Spearman correlation of a cluster in a pairwise data array.
 *
 * @param x
 * @param y
 * @param labels
 * @param cluster
 * @param minSamples
 */
float Spearman::computeCluster(
   const float *x,
   const float *y,
   const QVector<qint8>& labels,
   qint8 cluster,
   int minSamples)
{
   // extract samples in pairwise cluster
   int n = 0;

   for ( int i = 0; i < labels.size(); ++i )
   {
      if ( labels[i] == cluster )
      {
         _x_rank[n] = x[i];
         _y_rank[n] = y[i];
         ++n;
      }
   }

   // compute correlation only if there are enough samples
   float result = NAN;

   if ( n >= minSamples )
   {
      // compute rank of x
      heapSort(_x_rank, _y_rank, n);
      computeRank(_x_rank, n);

      // compute rank of y
      heapSort(_y_rank, _x_rank, n);
      computeRank(_y_rank, n);

      // compute correlation of rank arrays
      float sumx = 0;
      float sumy = 0;
      float sumx2 = 0;
      float sumy2 = 0;
      float sumxy = 0;

      for ( int i = 0; i < n; ++i )
      {
         float x_i = _x_rank[i];
         float y_i = _y_rank[i];

         sumx += x_i;
         sumy += y_i;
         sumx2 += x_i * x_i;
         sumy2 += y_i * y_i;
         sumxy += x_i * y_i;
      }

      result = (n*sumxy - sumx*sumy) / sqrt((n*sumx2 - sumx*sumx) * (n*sumy2 - sumy*sumy));
   }

   return result;
}






void Spearman::siftDown(QVector<float>& array, QVector<float>& extra, int start, int end)
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
         std::swap(array[root], array[swp]);
         std::swap(extra[root], extra[swp]);
         root = swp;
      }
   }
}






/*!
 * Sort an array using heapsort, while also applying the same swap operations
 * to a second array of the same size.
 *
 * @param array
 * @param extra
 * @param n
 */
void Spearman::heapSort(QVector<float>& array, QVector<float>& extra, int n)
{
   // heapify the array
   int start = ((n-1) - 1) / 2;

   while ( start >= 0 )
   {
      siftDown(array, extra, start, n - 1);
      start -= 1;
   }

   // sort the array
   int end = n - 1;
   while ( end > 0 )
   {
      std::swap(array[end], array[0]);
      std::swap(extra[end], extra[0]);
      end -= 1;

      siftDown(array, extra, 0, end);
   }
}






/*!
 * Compute the rank of a sorted vector in place. In the event of ties,
 * the ranks are corrected using fractional ranking.
 *
 * @param array
 * @param n
 */
void Spearman::computeRank(QVector<float>& array, int n)
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
