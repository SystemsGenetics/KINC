#include "genepair_correlation.h"



using namespace GenePair;






void Correlation::fetchData(ExpressionMatrix* input, Vector vector, const CCMatrix::Pair& pair, int cluster, int minExpression, QVector<double>& x, QVector<double>& y)
{
   // read in gene expressions
   ExpressionMatrix::Gene gene1(input);
   ExpressionMatrix::Gene gene2(input);

   gene1.read(vector.geneX());
   gene2.read(vector.geneY());

   // populate x and y with shared expressions of gene pair
   x.clear();
   y.clear();

   if ( pair.clusterSize() > 0 )
   {
      // add samples that are in the cluster
      for ( int i = 0; i < input->getSampleSize(); ++i )
      {
         if ( pair.at(cluster, i) == 1 )
         {
            x.append(gene1.at(i));
            y.append(gene2.at(i));
         }
      }
   }
   else
   {
      // add samples that are valid
      for ( int i = 0; i < input->getSampleSize(); ++i )
      {
         if ( !std::isnan(gene1.at(i)) && !std::isnan(gene2.at(i)) && minExpression <= gene1.at(i) && minExpression <= gene2.at(i) )
         {
            x.append(gene1.at(i));
            y.append(gene2.at(i));
         }
      }
   }
}
