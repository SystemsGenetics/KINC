#include "cluster_filter_serial.h"
#include "cluster_filter_resultblock.h"
#include "cluster_filter_workblock.h"
#include "expressionmatrix_gene.h"
#include "pairwise_gmm.h"
#include "pairwise_pearson.h"
#include "pairwise_spearman.h"
#include <ace/core/elog.h>
#include "stats.hpp"


using namespace std;


/*!
 * Construct a new serial object with the given analytic as its parent.
 *
 * @param parent
 */
ClusterFilter::Serial::Serial(ClusterFilter* parent):
   EAbstractAnalyticSerial(parent),
   _base(parent)
{
   EDEBUG_FUNC(this,parent);

}


/*!
 * Read in the given work block and save the results in a new result block. This
 * implementation takes the starting pairwise index and pair size from the work
 * block and processes those pairs.
 *
 * @param block
 */
std::unique_ptr<EAbstractAnalyticBlock> ClusterFilter::Serial::execute(const EAbstractAnalyticBlock* block)
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

   // iterate through all pairs
   Pairwise::Index index {workBlock->start()};

   // create iterators for the CCM and CMX data objects.
   CCMatrix::Pair ccmPair = CCMatrix::Pair(_base->_ccm);
   CorrelationMatrix::Pair cmxPair = CorrelationMatrix::Pair(_base->_cmx);

   // move to the location in the CCM/CMX matrix where the work block starts.
   ccmPair.read(index);
   cmxPair.read(index);

   // iterate through the elements in the workblock.
   for ( int i = 0; i < workBlock->size(); ++i )
   {
       // Get the number of samples and clusters.
       int num_clusters = cmxPair.clusterSize();
       int num_samples = _base->_emx->sampleSize();

       // How many clusters will remain?
       qint8 num_final_K = 0;

       // Get the list of correlations and labels for each cluster.
       QVector<float> correlations = cmxPair.correlations();
       QVector<qint8> labels(num_samples, 0);

       // Initialize new correlation and labels lists
       QVector<float> new_correlations;
       QVector<qint8> new_labels;
       QVector<int> k_num_samples(num_clusters, 0);

       // Rebuild the pair labels.
       for ( qint8 k = 0; k < num_clusters; k++ ) {
           for ( int j = 0; j < num_samples; j++ ) {
             qint8 val = ccmPair.at(k, j);
             if (val == 1) {
               labels[j] = k;
               k_num_samples[k]++;
             }
             else {
               labels[j] = -val;
             }
           }
       }

       // iterate through the clusters.
       for ( qint8 k = 0; k < num_clusters; k++ ) {

           // Do the correlation power test. The following code is modeled
           // after the `pwr.r.test` function of the `pwr` package for R. Here
           // we calculate the power given the signficance level (alpha) provided
           // by the user, the number of samples in the cluster and we
           // compare the calculated power to that expected by the user.
           // This code uses functions from the Keith OHare StatsLib at
           // https://www.kthohr.com/statslib.html
           int n = k_num_samples[k];;
           double r = correlations[k];
           double sig_level = _base->_powerThresholdAlpha;
           double ttt = stats::qt(sig_level / 2, n - 2);
           double ttt_2 = pow(ttt,2);
           double rc = sqrt(ttt_2/(ttt_2 + n - 2));
           double zr = atanh(r) + r / (2 * (n - 1));
           double zrc = atanh(rc);
           double power = stats::pnorm((zr - zrc) * sqrt(n - 3), 0.0, 1) +
                          stats::pnorm((-zr - zrc) * sqrt(n - 3), 0.0, 1);

           // If the calculated power is >= the expected power then we
           // can keep this cluster.
           if (power >= _base->_powerThresholdPower) {
             new_correlations.append(correlations[k]);
             new_labels = labels;
             num_final_K++;
           }
       }

       // save pairwise output data
       Pair pair;
       pair.K = num_final_K;
       if (num_final_K > 0) {
           pair.correlations = new_correlations;
           pair.labels = new_labels;
       }
       resultBlock->append(pair);

       // move to the location in the CCM/CMX matrix where the work block starts.
       ccmPair.readNext();
       cmxPair.readNext();
   }
   // return result block
   return unique_ptr<EAbstractAnalyticBlock>(resultBlock);
}

