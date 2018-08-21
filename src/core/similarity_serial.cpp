#include "similarity_serial.h"
#include "similarity_resultblock.h"
#include "similarity_workblock.h"
#include "expressionmatrix_gene.h"



using namespace std;






Similarity::Serial::Serial(Similarity* parent):
   EAbstractAnalytic::Serial(parent),
   _base(parent)
{
   // initialize clustering model
   if ( _base->_clusMethod != ClusteringMethod::None )
   {
      _base->_clusModel->initialize(_base->_input);
   }

   // initialize correlation model
   _base->_corrModel->initialize(_base->_input);
}






std::unique_ptr<EAbstractAnalytic::Block> Similarity::Serial::execute(const EAbstractAnalytic::Block* block)
{
   // cast block to work block
   const WorkBlock* workBlock {block->cast<WorkBlock>()};

   // initialize result block
   ResultBlock* resultBlock {new ResultBlock(workBlock->index(), workBlock->start())};

   // initialize workspace
   QVector<Pairwise::Vector2> data(_base->_input->sampleSize());
   QVector<qint8> labels(_base->_input->sampleSize());

   // iterate through all pairs
   Pairwise::Index index {workBlock->start()};

   for ( int i = 0; i < workBlock->size(); ++i )
   {
      // fetch pairwise input data
      int numSamples = fetchPair(index, data, labels);

      // remove pre-clustering outliers
      if ( _base->_removePreOutliers )
      {
         markOutliers(data, numSamples, labels, 0, -7);
      }

      // compute clusters
      qint8 K {1};

      if ( _base->_clusMethod != ClusteringMethod::None )
      {
         K = _base->_clusModel->compute(
            data,
            numSamples,
            labels,
            _base->_minSamples,
            _base->_minClusters,
            _base->_maxClusters,
            _base->_criterion
         );
      }

      // remove post-clustering outliers
      if ( K > 1 && _base->_removePostOutliers )
      {
         for ( qint8 k = 0; k < K; ++k )
         {
            markOutliers(data, numSamples, labels, k, -8);
         }
      }

      // compute correlations
      QVector<float> correlations = _base->_corrModel->compute(
         data,
         K,
         labels,
         _base->_minSamples
      );

      // save pairwise output data
      Pair pair;
      pair.K = K;

      if ( K > 1 )
      {
         pair.labels = labels;
      }

      if ( K > 0 )
      {
         pair.correlations = correlations;
      }

      resultBlock->append(pair);

      // increment to next pair
      ++index;
   }

   // return result block
   return unique_ptr<EAbstractAnalytic::Block>(resultBlock);
}






/*!
 * Extract pairwise data from an expression matrix given a pairwise index. Samples
 * with missing values and samples that fall below the expression threshold are
 * excluded. The number of extracted samples is returned.
 *
 * @param index
 * @param data
 * @param labels
 */
int Similarity::Serial::fetchPair(Pairwise::Index index, QVector<Pairwise::Vector2>& data, QVector<qint8>& labels)
{
   // read in gene expressions
   ExpressionMatrix::Gene gene1(_base->_input);
   ExpressionMatrix::Gene gene2(_base->_input);

   gene1.read(index.getX());
   gene2.read(index.getY());

   // extract pairwise samples
   int numSamples = 0;

   for ( int i = 0; i < _base->_input->sampleSize(); ++i )
   {
      // exclude samples with missing values
      if ( std::isnan(gene1.at(i)) || std::isnan(gene2.at(i)) )
      {
         labels[i] = -9;
      }

      // exclude samples which fall below the expression threshold
      else if ( gene1.at(i) < _base->_minExpression || gene2.at(i) < _base->_minExpression )
      {
         labels[i] = -6;
      }

      // include any remaining samples
      else
      {
         data[numSamples] = { gene1.at(i), gene2.at(i) };
         numSamples++;

         labels[i] = 0;
      }
   }

   // return number of extracted samples
   return numSamples;
}






/*!
 * Mark outliers in a vector of pairwise data. Outliers are detected independently
 * on each axis using the Tukey method, and marked with the given marker. Only those
 * samples in the given cluster are used in outlier detection. For unclustered data,
 * all samples should be labeled as 0, so a cluster value of 0 should be used.
 *
 * @param data
 * @param N
 * @param labels
 * @param cluster
 * @param marker
 */
void Similarity::Serial::markOutliers(const QVector<Pairwise::Vector2>& data, int N, QVector<qint8>& labels, qint8 cluster, qint8 marker)
{
   // extract univariate data from the given cluster
   QVector<float> x_sorted;
   QVector<float> y_sorted;

   x_sorted.reserve(N);
   y_sorted.reserve(N);

   for ( int i = 0; i < N; i++ )
   {
      if ( labels[i] == cluster )
      {
         x_sorted.append(data[i].s[0]);
         y_sorted.append(data[i].s[1]);
      }
   }

   // return if the given cluster is empty
   if ( x_sorted.size() == 0 || y_sorted.size() == 0 )
   {
      return;
   }

   std::sort(x_sorted.begin(), x_sorted.end());
   std::sort(y_sorted.begin(), y_sorted.end());

   // compute interquartile range and thresholds for each axis
   const int n = x_sorted.size();

   float Q1_x = x_sorted[n * 1 / 4];
   float Q3_x = x_sorted[n * 3 / 4];
   float T_x_min = Q1_x - 1.5f * (Q3_x - Q1_x);
   float T_x_max = Q3_x + 1.5f * (Q3_x - Q1_x);

   float Q1_y = y_sorted[n * 1 / 4];
   float Q3_y = y_sorted[n * 3 / 4];
   float T_y_min = Q1_y - 1.5f * (Q3_y - Q1_y);
   float T_y_max = Q3_y + 1.5f * (Q3_y - Q1_y);

   // mark outliers
   for ( int i = 0; i < N; ++i )
   {
      if ( labels[i] == cluster )
      {
         // mark samples that are outliers on either axis
         bool outlier_x = (data[i].s[0] < T_x_min || T_x_max < data[i].s[0]);
         bool outlier_y = (data[i].s[1] < T_y_min || T_y_max < data[i].s[1]);

         if ( outlier_x || outlier_y )
         {
            labels[i] = marker;
         }
      }
   }
}
