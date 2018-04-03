#include "genepair_clustering.h"



using namespace GenePair;






void Clustering::compute(
   ExpressionMatrix* input,
   Vector vector,
   int minSamples,
   int minExpression,
   qint8 minClusters,
   qint8 maxClusters,
   Criterion criterion,
   bool removePreOutliers,
   bool removePostOutliers)
{
   // pre-allocate workspace
   _X.reserve(input->getSampleSize());
   _labels.resize(input->getSampleSize());
   _bestLabels.resize(input->getSampleSize());

   // fetch data matrix X from expression matrix
   fetchData(input, vector, minExpression);

   // remove pre-clustering outliers
   if ( removePreOutliers )
   {
      markOutliers(_X, 0, _bestLabels, 0, -7);
      markOutliers(_X, 1, _bestLabels, 0, -7);
   }

   // perform clustering only if there are enough samples
   _bestK = 0;

   if ( _X.size() >= minSamples )
   {
      float bestValue = INFINITY;

      for ( qint8 K = minClusters; K <= maxClusters; ++K )
      {
         // run each clustering model
         bool success = fit(_X, K, _labels);

         if ( !success )
         {
            continue;
         }

         // evaluate model
         float value = INFINITY;

         switch (criterion)
         {
         case Criterion::BIC:
            value = computeBIC(K, logLikelihood(), _X.size(), 2);
            break;
         case Criterion::ICL:
            value = computeICL(K, logLikelihood(), _X.size(), 2, entropy());
            break;
         }

         // save the best model
         if ( value < bestValue )
         {
            _bestK = K;
            bestValue = value;

            for ( int i = 0, j = 0; i < _X.size(); ++i )
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
            markOutliers(_X, 0, _bestLabels, k, -8);
            markOutliers(_X, 1, _bestLabels, k, -8);
         }
      }
   }
}






void Clustering::fetchData(ExpressionMatrix* input, Vector vector, int minExpression)
{
   // read in gene expressions
   ExpressionMatrix::Gene gene1(input);
   ExpressionMatrix::Gene gene2(input);

   gene1.read(vector.geneX());
   gene2.read(vector.geneY());

   // populate X with shared expressions of gene pair
   _X.clear();

   for ( int i = 0; i < input->getSampleSize(); ++i )
   {
      if ( std::isnan(gene1.at(i)) || std::isnan(gene2.at(i)) )
      {
         _labels[i] = -9;
      }
      else if ( gene1.at(i) < minExpression || gene2.at(i) < minExpression )
      {
         _labels[i] = -6;
      }
      else
      {
         _X.append({ gene1.at(i), gene2.at(i) });

         _labels[i] = 0;
      }
   }
}






void Clustering::markOutliers(const QVector<Vector2>& X, int j, QVector<qint8>& labels, qint8 cluster, qint8 marker)
{
   // compute x_sorted = X[:, j], filtered and sorted
   QVector<float> x_sorted;
   x_sorted.reserve(X.size());

   for ( int i = 0; i < X.size(); i++ )
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
   for ( int i = 0; i < X.size(); ++i )
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
