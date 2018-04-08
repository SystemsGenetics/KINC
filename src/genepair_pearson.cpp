#include <ace/core/metadata.h>

#include "genepair_pearson.h"



using namespace GenePair;






void Pearson::initialize(ExpressionMatrix* input, CorrelationMatrix* output)
{
   // pre-allocate workspace
   _X.reserve(input->getSampleSize());

   // initialize correlation matrix
   EMetadata correlations(EMetadata::Array);
   EMetadata* name {new EMetadata(EMetadata::String)};
   *(name->toString()) = "pearson";
   correlations.toArray()->append(name);

   output->initialize(input->getGeneNames(), correlations);
}






float Pearson::compute(
   const QVector<Vector2>& data,
   const QVector<qint8>& labels, qint8 cluster,
   int minSamples)
{
   // extract samples in gene pair cluster
   _X.clear();

   for ( int i = 0; i < data.size(); ++i )
   {
      if ( labels[i] == cluster )
      {
         _X.append(data[i]);
      }
   }

   // compute correlation only if there are enough samples
   const int n = _X.size();

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
         sumx += _X[i].s[0];
         sumy += _X[i].s[1];
         sumx2 += _X[i].s[0] * _X[i].s[0];
         sumy2 += _X[i].s[1] * _X[i].s[1];
         sumxy += _X[i].s[0] * _X[i].s[1];
      }

      // compute Pearson correlation coefficient
      result = (n*sumxy - sumx*sumy) / sqrt((n*sumx2 - sumx*sumx) * (n*sumy2 - sumy*sumy));
   }

   return result;
}
