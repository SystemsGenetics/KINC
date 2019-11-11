#include "pairwise_pearson.h"



using namespace Pairwise;



/*!
 * Compute the Pearson correlation of a cluster in a pairwise data array.
 *
 * @param x
 * @param y
 * @param labels
 * @param cluster
 * @param minSamples
 */
float Pearson::computeCluster(
   const float *x,
   const float *y,
   const QVector<qint8>& labels,
   qint8 cluster,
   int minSamples)
{
   // compute intermediate sums
   int n = 0;
   float sumx = 0;
   float sumy = 0;
   float sumx2 = 0;
   float sumy2 = 0;
   float sumxy = 0;

   for ( int i = 0; i < labels.size(); ++i )
   {
      if ( labels[i] == cluster )
      {
         float x_i = x[i];
         float y_i = y[i];

         sumx += x_i;
         sumy += y_i;
         sumx2 += x_i * x_i;
         sumy2 += y_i * y_i;
         sumxy += x_i * y_i;

         ++n;
      }
   }

   // compute correlation only if there are enough samples
   float result = NAN;

   if ( n >= minSamples )
   {
      result = (n*sumxy - sumx*sumy) / sqrtf((n*sumx2 - sumx*sumx) * (n*sumy2 - sumy*sumy));
   }

   return result;
}
