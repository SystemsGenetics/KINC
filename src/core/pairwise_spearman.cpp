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
   _x.resize(emx->sampleSize());
   _y.resize(emx->sampleSize());
   _rank.resize(emx->sampleSize());
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
         _x[n] = x[i];
         _y[n] = y[i];
         _rank[n] = n;
         ++n;
      }
   }

   // compute correlation only if there are enough samples
   float result = NAN;

   if ( n >= minSamples )
   {
      // execute two sorts that are the beginning of the spearman algorithm
      heapSort(_x, _y, n);
      heapSort(_y, _rank, n);

      // go through spearman sorted rank list and calculate difference from 1,2,3,... list
      int diff = 0;

      for ( int i = 0; i < n; ++i )
      {
         int tmp = i - _rank[i];
         diff += tmp*tmp;
      }

      // compute spearman coefficient
      result = 1.0f - 6.0f * diff / (n * (n*n - 1));
   }

   return result;
}






template<typename T>
void Spearman::siftDown(QVector<float>& array, QVector<T>& extra, int start, int end)
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
template<typename T>
void Spearman::heapSort(QVector<float>& array, QVector<T>& extra, int n)
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
