#include "corrpower_serial.h"
#include "corrpower_resultblock.h"
#include "corrpower_workblock.h"
#include <ace/core/elog.h>
#include "stats.hpp"



using namespace std;



/*!
 * Construct a new serial object with the given analytic as its parent.
 *
 * @param parent
 */
CorrPowerFilter::Serial::Serial(CorrPowerFilter* parent):
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
std::unique_ptr<EAbstractAnalyticBlock> CorrPowerFilter::Serial::execute(const EAbstractAnalyticBlock* block)
{
    EDEBUG_FUNC(this,block);

    if ( ELog::isActive() )
    {
        ELog() << tr("Executing(serial) work index %1.\n").arg(block->index());
    }

    // Cast block to work block.
    const WorkBlock* workBlock {block->cast<WorkBlock>()};

    // Initialize result block.
    ResultBlock* resultBlock {new ResultBlock(workBlock->index())};

    // Create iterators for the CCM and CMX data objects.
    CCMatrix::Pair ccmPair(_base->_ccm);
    CorrelationMatrix::Pair cmxPair(_base->_cmx);

    // Read the first pair of the work block.
    Pairwise::Index index(workBlock->start());

    ccmPair.read(index);
    cmxPair.read(index);

    // Iterate through each pair in the work block.
    for ( qint64 wbIndex = 0; wbIndex < workBlock->size(); ++wbIndex )
    {

        // Print warning if iterator indices do not match
        if ( ccmPair.index() != cmxPair.index() )
        {
            QString source = _base->_cmx->geneNames().at(cmxPair.index().getX()).toString();
            QString target = _base->_cmx->geneNames().at(cmxPair.index().getY()).toString();
            qInfo() << "warning: ccm and cmx files are out of sync at cmx coordinate:"
                    << source << "," << target << " ("
                    << cmxPair.index().getX() << "," << cmxPair.index().getY() <<").";
        }

        // Get the number of samples and clusters.
        int num_clusters = ccmPair.clusterSize();
        int num_samples = _base->_ccm->sampleSize();

        // Initialize new correlation and labels lists
        QVector<float> new_correlations;
        QVector<qint8> new_labels;
        QVector<int> k_num_samples(num_clusters, 0);
        QVector<int> k_keep;

        // Rebuild the pair labels from the pair sample mask. These
        // pair labels are the same as when the similarity analytic makes them.
        QVector<qint8> labels(num_samples, -128);

        for ( qint8 k = 0; k < num_clusters; k++ )
        {
            for ( int j = 0; j < num_samples; j++ )
            {
                qint8 val = ccmPair.at(k, j);

                // If the sample belongs to this cluster the value is 1
                if ( val == 1 )
                {
                    labels[j] = k;
                    k_num_samples[k]++;
                }

                // If the sample is not in another cluster then it was removed
                // or missing.
                else if ( val != 0 )
                {
                    labels[j] = -val;
                }
            }
        }

        // Iterate through the clusters and perform a correlation power analysis.
        for ( qint8 k = 0; k < num_clusters; k++ )
        {
            float r = cmxPair.correlations().at(k);

            // Perform the power analysis test.
            double power = pwr_r_test(
                abs(r),
                k_num_samples[k],
                _base->_powerThresholdAlpha);

            // If the calculated power is >= the expected power then we
            // can keep this cluster.  We keep it by adding the correlation
            // and the labels from the original cluster into new variables
            // for the correlation and labels.
            if ( power >= _base->_powerThresholdPower )
            {
                new_correlations.append(r);
                k_keep.append(k);
                new_labels = labels;
            }
        }

        // append pair to result block
        if ( new_correlations.size() > 0 )
        {
            resultBlock->append(Pair {
                cmxPair.index(),
                new_labels,
                new_correlations,
                k_keep
            });
        }

        // read the next pair
        cmxPair.readNext();
        ccmPair.read(cmxPair.index());

    }

    // We're done! Return the result block.
    return unique_ptr<EAbstractAnalyticBlock>(resultBlock);
}



/*!
 * Performs the correlation power test.
 *
 * @param r The correlation score
 * @param n The number of samples in the cluster
 * @param sig_level The desired signficance level (alpha) value.
 * @return The power from the test.
 */
double CorrPowerFilter::Serial::pwr_r_test(double r, int n, double sig_level)
{
    // The following code is modeled after the `pwr.r.test` function
    // of the `pwr` package for R. Here we calculate the power given the
    // signficance level (alpha) provided by the user, the number of
    // samples in the cluster and we compare the calculated power to that
    // expected by the user. This code uses functions from the Keith OHare StatsLib at
    // https://www.kthohr.com/statslib.html
    double ttt = stats::qt(sig_level / 2.0, n - 2.0);
    double ttt_2 = pow(ttt, 2);
    double rc = sqrt(ttt_2/(ttt_2 + (n - 2.0)));
    double zr = atanh(r) + r / (2.0 * (n - 1.0));
    double zrc = atanh(rc);
    double pnzrc = stats::pnorm((-zr - zrc) * sqrt(n - 3.0), 0.0, 1.0);
    double power = stats::pnorm((zr - zrc) * sqrt(n - 3.0) + pnzrc, 0.0, 1.0);

    return power;
}
