#include "similarity_serial.h"
#include "similarity_resultblock.h"
#include "similarity_workblock.h"
#include "expressionmatrix_gene.h"
#include "pairwise_gmm.h"
#include "pairwise_pearson.h"
#include "pairwise_spearman.h"
#include <ace/core/elog.h>



using namespace std;






/*!
 * Construct a new serial object with the given analytic as its parent.
 *
 * @param parent
 */
Similarity::Serial::Serial(Similarity* parent):
   EAbstractAnalytic::Serial(parent),
   _base(parent)
{
   EDEBUG_FUNC(this,parent);

   // initialize clustering model
   switch ( _base->_clusMethod )
   {
   case ClusteringMethod::None:
      _clusModel = nullptr;
      break;
   case ClusteringMethod::GMM:
      _clusModel = new Pairwise::GMM(_base->_input);
      break;
   }

   // initialize correlation model
   switch ( _base->_corrMethod )
   {
   case CorrelationMethod::Pearson:
      _corrModel = new Pairwise::Pearson();
      break;
   case CorrelationMethod::Spearman:
      _corrModel = new Pairwise::Spearman(_base->_input);
      break;
   }
}






/*!
 * Read in the given work block and save the results in a new result block. This
 * implementation takes the starting pairwise index and pair size from the work
 * block and processes those pairs.
 *
 * @param block
 */
std::unique_ptr<EAbstractAnalytic::Block> Similarity::Serial::execute(const EAbstractAnalytic::Block* block)
{
   EDEBUG_FUNC(this,block);

   if ( ELog::isActive() )
   {
      ELog() << tr("Executing(serial) work index %1.\n").arg(block->index());
   }

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
         numSamples = removeOutliers(data, numSamples, labels, 1, -7);
      }

      // compute clusters
      qint8 K {1};

      if ( _base->_clusMethod != ClusteringMethod::None )
      {
         K = _clusModel->compute(
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
      if ( _base->_removePostOutliers )
      {
         numSamples = removeOutliers(data, numSamples, labels, K, -8);
      }

      // compute correlations
      QVector<float> correlations = _corrModel->compute(
         data,
         K,
         labels,
         _base->_minSamples
      );

      // save pairwise output data
      Pair pair;
      pair.K = K;

      if ( K > 0 )
      {
         pair.labels = labels;
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
int Similarity::Serial::fetchPair(const Pairwise::Index& index, QVector<Pairwise::Vector2>& data, QVector<qint8>& labels)
{
   EDEBUG_FUNC(this,&index,&data,&labels);

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
 * Remove outliers from a vector of pairwise data. Outliers are detected independently
 * on each axis using the Tukey method, and marked with the given marker. Only the
 * samples in the given cluster are used in outlier detection. For unclustered data,
 * all samples are labeled as 0, so a cluster value of 0 should be used. The data
 * array should only contain samples that have a non-negative label.
 *
 * @param data
 * @param labels
 * @param cluster
 * @param marker
 */
int Similarity::Serial::removeOutliersCluster(QVector<Pairwise::Vector2>& data, QVector<qint8>& labels, qint8 cluster, qint8 marker)
{
   EDEBUG_FUNC(this,&data,&labels,cluster,marker);

   // extract univariate data from the given cluster
   QVector<float> x_sorted;
   QVector<float> y_sorted;

   x_sorted.reserve(labels.size());
   y_sorted.reserve(labels.size());

   for ( int i = 0, j = 0; i < labels.size(); i++ )
   {
      if ( labels[i] >= 0 )
      {
         if ( labels[i] == cluster )
         {
            x_sorted.append(data[j].s[0]);
            y_sorted.append(data[j].s[1]);
         }

         j++;
      }
   }

   // return if the given cluster is empty
   if ( x_sorted.size() == 0 || y_sorted.size() == 0 )
   {
      return 0;
   }

   // sort samples for each axis
   std::sort(x_sorted.begin(), x_sorted.end());
   std::sort(y_sorted.begin(), y_sorted.end());

   // compute quartiles and thresholds for each axis
   const int n = x_sorted.size();

   float Q1_x = x_sorted[n * 1 / 4];
   float Q3_x = x_sorted[n * 3 / 4];
   float T_x_min = Q1_x - 1.5f * (Q3_x - Q1_x);
   float T_x_max = Q3_x + 1.5f * (Q3_x - Q1_x);

   float Q1_y = y_sorted[n * 1 / 4];
   float Q3_y = y_sorted[n * 3 / 4];
   float T_y_min = Q1_y - 1.5f * (Q3_y - Q1_y);
   float T_y_max = Q3_y + 1.5f * (Q3_y - Q1_y);

   // remove outliers
   int numSamples = 0;

   for ( int i = 0, j = 0; i < labels.size(); i++ )
   {
      if ( labels[i] >= 0 )
      {
         // mark samples in the given cluster that are outliers on either axis
         if ( labels[i] == cluster && (data[j].s[0] < T_x_min || T_x_max < data[j].s[0] || data[j].s[1] < T_y_min || T_y_max < data[j].s[1]) )
         {
            labels[i] = marker;
         }

         // preserve all other non-outlier samples in the data array
         else
         {
            data[numSamples] = data[j];
            numSamples++;
         }

         j++;
      }
   }

   // return number of remaining samples
   return numSamples;
}






/*!
 * Perform outlier removal on each cluster in a parwise data array.
 *
 * @param data
 * @param numSamples
 * @param labels
 * @param clusterSize
 * @param marker
 */
int Similarity::Serial::removeOutliers(QVector<Pairwise::Vector2>& data, int numSamples, QVector<qint8>& labels, qint8 clusterSize, qint8 marker)
{
   EDEBUG_FUNC(this,&data,numSamples,&labels,clusterSize,marker);

   // do not perform post-clustering outlier removal if there is only one cluster
   if ( marker == -8 && clusterSize <= 1 )
   {
      return numSamples;
   }

   // perform outlier removal on each cluster
   for ( qint8 k = 0; k < clusterSize; ++k )
   {
      numSamples = removeOutliersCluster(data, labels, k, marker);
   }

   return numSamples;
}
