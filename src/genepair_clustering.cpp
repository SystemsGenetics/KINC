#include "genepair_clustering.h"



using namespace GenePair;






void Clustering::initialize(ExpressionMatrix* input)
{
   _input = input;

   // pre-allocate workspace
   _X.resize(_input->getSampleSize());
   _labels.resize(_input->getSampleSize());
   _bestLabels.resize(_input->getSampleSize());
}






void Clustering::compute(
   Index index,
   int minSamples,
   int minExpression,
   qint8 minClusters,
   qint8 maxClusters,
   Criterion criterion,
   bool removePreOutliers,
   bool removePostOutliers)
{
   // fetch data matrix X from expression matrix
   int numSamples = fetchData(index, minExpression, _X, _bestLabels);

   // remove pre-clustering outliers
   if ( removePreOutliers )
   {
      markOutliers(_X, numSamples, 0, _bestLabels, 0, -7);
      markOutliers(_X, numSamples, 1, _bestLabels, 0, -7);
   }

   // perform clustering only if there are enough samples
   _bestK = 0;

   if ( numSamples >= minSamples )
   {
      float bestValue = INFINITY;

      for ( qint8 K = minClusters; K <= maxClusters; ++K )
      {
         // run each clustering model
         bool success = fit(_X, numSamples, K, _labels);

         if ( !success )
         {
            continue;
         }

         // evaluate model
         float value = INFINITY;

         switch (criterion)
         {
         case Criterion::BIC:
            value = computeBIC(K, logLikelihood(), numSamples, 2);
            break;
         case Criterion::ICL:
            value = computeICL(K, logLikelihood(), numSamples, 2, entropy());
            break;
         }

         // save the best model
         if ( value < bestValue )
         {
            _bestK = K;
            bestValue = value;

            for ( int i = 0, j = 0; i < numSamples; ++i )
            {
               if ( _bestLabels[i] >= 0 )
               {
                  _bestLabels[i] = _labels[j];
                  ++j;
               }
            }
         }
      }
   }

   if ( _bestK > 1 )
   {
      // remove post-clustering outliers
      if ( removePostOutliers )
      {
         for ( qint8 k = 0; k < _bestK; ++k )
         {
            markOutliers(_X, numSamples, 0, _bestLabels, k, -8);
            markOutliers(_X, numSamples, 1, _bestLabels, k, -8);
         }
      }
   }
}






int Clustering::fetchData(Index index, int minExpression, QVector<Vector2>& X, QVector<qint8>& labels)
{
   // read in gene expressions
   ExpressionMatrix::Gene gene1(_input);
   ExpressionMatrix::Gene gene2(_input);

   gene1.read(index.getX());
   gene2.read(index.getY());

   // populate X with shared expressions of gene pair
   int numSamples = 0;

   for ( int i = 0; i < _input->getSampleSize(); ++i )
   {
      if ( std::isnan(gene1.at(i)) || std::isnan(gene2.at(i)) )
      {
         labels[i] = -9;
      }
      else if ( gene1.at(i) < minExpression || gene2.at(i) < minExpression )
      {
         labels[i] = -6;
      }
      else
      {
         X[numSamples] = { gene1.at(i), gene2.at(i) };
         numSamples++;

         labels[i] = 0;
      }
   }

   // return size of X
   return numSamples;
}






void Clustering::markOutliers(const QVector<Vector2>& X, int N, int j, QVector<qint8>& labels, qint8 cluster, qint8 marker)
{
   // compute x_sorted = X[:, j], filtered and sorted
   QVector<float> x_sorted;
   x_sorted.reserve(N);

   for ( int i = 0; i < N; i++ )
   {
      if ( labels[i] == cluster || labels[i] == marker )
      {
         x_sorted.append(X[i].s[j]);
      }
   }

   if ( x_sorted.size() == 0 )
   {
      return;
   }

   std::sort(x_sorted.begin(), x_sorted.end());

   // compute quartiles, interquartile range, upper and lower bounds
   const int n = x_sorted.size();

   float Q1 = x_sorted[n * 1 / 4];
   float Q3 = x_sorted[n * 3 / 4];

   float T_min = Q1 - 1.5f * (Q3 - Q1);
   float T_max = Q3 + 1.5f * (Q3 - Q1);

   // mark outliers
   for ( int i = 0; i < N; ++i )
   {
      if ( labels[i] == cluster && (X[i].s[j] < T_min || T_max < X[i].s[j]) )
      {
         labels[i] = marker;
      }
   }
}






float Clustering::computeBIC(int K, float logL, int N, int D)
{
   int p = K * (1 + D + D * D);

   return log(N) * p - 2 * logL;
}






float Clustering::computeICL(int K, float logL, int N, int D, float E)
{
   int p = K * (1 + D + D * D);

   return log(N) * p - 2 * logL - 2 * E;
}
