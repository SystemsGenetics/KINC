#include "corrpower.h"
#include "corrpower_input.h"
#include "corrpower_serial.h"
#include "corrpower_resultblock.h"
#include "corrpower_workblock.h"
#include "datafactory.h"
#include "expressionmatrix_gene.h"
#include "ccmatrix_pair.h"
#include "correlationmatrix_pair.h"
#include <ace/core/ace_qmpi.h>
#include <ace/core/elog.h>



using namespace std;



/*!
 * Return the total number of blocks this analytic must process as steps
 * or blocks of work.
 */
int CorrPowerFilter::size() const
{
    EDEBUG_FUNC(this);

    return (static_cast<qint64>(_ccm->size()) + _workBlockSize - 1) / _workBlockSize;
}



/*!
 * Create and return a work block for this analytic with the given index. This
 * implementation creates a work block with a start index and size denoting the
 * number of pairs to process.
 *
 * @param index
 */
std::unique_ptr<EAbstractAnalyticBlock> CorrPowerFilter::makeWork(int index) const
{
    EDEBUG_FUNC(this,index);

    if ( ELog::isActive() )
    {
        ELog() << tr("Making work index %1 of %2.\n").arg(index).arg(size());
    }

    qint64 start {index * static_cast<qint64>(_workBlockSize)};
    qint64 size {min(_cmx->size() - start, static_cast<qint64>(_workBlockSize))};

    return unique_ptr<EAbstractAnalyticBlock>(new WorkBlock(index, start, size));
}



/*!
 * Create an empty and uninitialized result block.
 */
std::unique_ptr<EAbstractAnalyticBlock> CorrPowerFilter::makeResult() const
{
    EDEBUG_FUNC(this);

    return unique_ptr<EAbstractAnalyticBlock>(new ResultBlock);
}



/*!
 * Create an empty and uninitialized work block.
 */
std::unique_ptr<EAbstractAnalyticBlock> CorrPowerFilter::makeWork() const
{
    EDEBUG_FUNC(this);

    return unique_ptr<EAbstractAnalyticBlock>(new WorkBlock);
}



/*!
 * Process the given index with a possible block of results if this analytic
 * produces work blocks. This implementation uses only the index of the result
 * block to determine which piece of work to do.
 *
 * @param result
 */
void CorrPowerFilter::process(const EAbstractAnalyticBlock* result)
{
    EDEBUG_FUNC(this);

    if ( ELog::isActive() )
    {
        ELog() << tr("Processing result %1 of %2.\n").arg(result->index()).arg(size());
    }

    // Iterate through the result block pairs.
    const ResultBlock* resultBlock {result->cast<ResultBlock>()};

    for ( auto& pair : resultBlock->pairs() )
    {
        // Create pair objects for both output data files.
        CCMatrix::Pair ccmPair(_ccmOut);
        CorrelationMatrix::Pair cmxPair(_cmxOut);

        // Iterate through the clusters in the pair.
        for ( qint8 k = 0; k < pair.correlations.size(); ++k )
        {
            // The pair.K indicates how many clusters remain in the pair
            // but the pair.labels are still the same as alwasy, and
            // we only want to create the sample string for kept clusters.
            // the pair.keep array lists the clusters that are kept.
            int ki = pair.keep[k];

            // determine whether correlation is within thresholds
            float corr = pair.correlations[k];

            // save sample string
            ccmPair.addCluster();

            // add each cluster sample string to the pair.
            for ( int i = 0; i < _ccm->sampleSize(); ++i )
            {
                qint8 val = pair.labels[i];

                if ( ki == val )
                {
                    val = 1;
                }

                // A value of -128 exists if the sample belong to
                // another cluster but that cluster is not present
                // in the input files.
                else if ( val == -128 )
                {
                    val = 0;
                }
                else if ( val > 0 )
                {
                    val = 0;
                }
                else
                {
                    val = -val;
                }

                ccmPair.at(k, i) = val;
            }

            // save correlation
            cmxPair.addCluster();
            cmxPair.at(k) = corr;
        }

        // write pairwise data to output files
        if ( ccmPair.clusterSize() > 0 )
        {
            ccmPair.write(pair.index);
            cmxPair.write(pair.index);
        }
    }
}



/*!
 * Make a new serial object and return its pointer.
 */
EAbstractAnalyticSerial* CorrPowerFilter::makeSerial()
{
    EDEBUG_FUNC(this);

    return new Serial(this);
}



/*!
 * Make a new input object and return its pointer.
 */
EAbstractAnalyticInput* CorrPowerFilter::makeInput()
{
    EDEBUG_FUNC(this);

    return new Input(this);
}



/*!
 * Initialize this analytic. This implementation checks to make sure the input
 * data objects and output file have been set.
 */
void CorrPowerFilter::initialize()
{
    EDEBUG_FUNC(this);

    // get MPI instance
    auto& mpi {Ace::QMPI::instance()};

    // only the master process needs to validate arguments
    if ( !mpi.isMaster() )
    {
        return;
    }

    // make sure input data is valid
    if ( !_ccm )
    {
        E_MAKE_EXCEPTION(e);
        e.setTitle(tr("Invalid Argument"));
        e.setDetails(tr("Did not get a valid input CCM data object."));
        throw e;
    }

    if ( !_cmx )
    {
        E_MAKE_EXCEPTION(e);
        e.setTitle(tr("Invalid Argument"));
        e.setDetails(tr("Did not get a valid input CMX data object."));
        throw e;
    }

    // initialize work block size
    if ( _workBlockSize == 0 )
    {
        int numWorkers = max(1, mpi.size() - 1);

        _workBlockSize = min(32768LL, _ccm->size() / numWorkers);
    }
}



/*!
 * Initialize the output data objects of this analytic.
 */
void CorrPowerFilter::initializeOutputs()
{
    EDEBUG_FUNC(this);

    // make sure output data is valid
    if ( !_ccmOut )
    {
        E_MAKE_EXCEPTION(e);
        e.setTitle(tr("Invalid Argument"));
        e.setDetails(tr("Did not get a valid CCM output data object."));
        throw e;
    }

    if ( !_cmxOut )
    {
        E_MAKE_EXCEPTION(e);
        e.setTitle(tr("Invalid Argument"));
        e.setDetails(tr("Did not get a valid CMX output data object."));
        throw e;
    }

    // initialize cluster matrix
    _ccmOut->initialize(_ccm->geneNames(), _ccm->maxClusterSize(), _ccm->sampleNames());

    // initialize correlation matrix
    _cmxOut->initialize(_cmx->geneNames(), _cmx->maxClusterSize(), _cmx->correlationName());
}
