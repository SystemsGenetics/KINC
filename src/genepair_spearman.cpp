#include <gsl/gsl_statistics.h>

#include "genepair_spearman.h"



using namespace GenePair;






double Spearman::compute(
   ExpressionMatrix* input,
   Vector vector,
   const CCMatrix::Pair& pair, int cluster,
   int minSamples,
   int minExpression)
{
   // pre-allocate workspace
   _x.reserve(input->getSampleSize());
   _y.reserve(input->getSampleSize());
   _work.resize(input->getSampleSize() * 2);

   // fetch a and b arrays from expression matrix
   fetchData(input, vector, pair, cluster, minExpression, _x, _y);

   // compute correlation only if there are enough samples
   float result = NAN;

   if ( _x.size() >= minSamples )
   {
      result = gsl_stats_spearman(_x.data(), 1, _y.data(), 1, _x.size(), _work.data());
   }

   return result;
}
