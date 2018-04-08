#include <ace/core/metadata.h>

#include "genepair_pearson.h"



using namespace GenePair;






void Pearson::initialize(ExpressionMatrix* input, CorrelationMatrix* output)
{
   // pre-allocate workspace
   _x.resize(input->getSampleSize());
   _y.resize(input->getSampleSize());

   // initialize correlation matrix
   EMetadata correlations(EMetadata::Array);
   EMetadata* name {new EMetadata(EMetadata::String)};
   *(name->toString()) = "pearson";
   correlations.toArray()->append(name);

   output->initialize(input->getGeneNames(), correlations);
}






float Pearson::computeCluster(
   const QVector<Vector2>& data,
   const QVector<qint8>& labels,
   qint8 cluster,
   int minSamples)
{
   // extract samples in gene pair cluster
   int n = 0;

   for ( int i = 0; i < data.size(); ++i )
   {
      if ( labels[i] == cluster )
      {
         _x[n] = data[i].s[0];
         _y[n] = data[i].s[1];
         ++n;
      }
   }

   // compute correlation only if there are enough samples
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
