#include "genepair_pearson.h"



using namespace GenePair;






double Pearson::compute(
   ExpressionMatrix* input,
   Vector vector,
   const CCMatrix::Pair& pair, int cluster,
   int minSamples,
   int minExpression)
{
   // pre-allocate workspace
   _x.reserve(input->getSampleSize());
   _y.reserve(input->getSampleSize());

   // fetch a and b arrays from expression matrix
   fetchData(input, vector, pair, cluster, minExpression, _x, _y);

   // compute correlation only if there are enough samples
   const int n = _x.size();

   float result = NAN;

   if ( n >= minSamples )
   {
      // compute intermediate sums
      float sumx = 0;
      float sumy = 0;
      float sumx2 = 0;
      float sumy2 = 0;
      float sumxy = 0;

      for ( int i = 0; i < n; ++i )
      {
         sumx += _x[i];
         sumy += _y[i];
         sumx2 += _x[i] * _x[i];
         sumy2 += _y[i] * _y[i];
         sumxy += _x[i] * _y[i];
      }

      // compute Pearson correlation coefficient
      result = (n*sumxy - sumx*sumy) / sqrt((n*sumx2 - sumx*sumx) * (n*sumy2 - sumy*sumy));
   }

   return result;
}
