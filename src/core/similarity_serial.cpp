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
    EAbstractAnalyticSerial(parent),
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
        _clusModel = new Pairwise::GMM(_base->_input, _base->_maxClusters);
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

    // initialize expression matrix
    _expressions = _base->_input->dumpRawData();
}



/*!
 * Read in the given work block and save the results in a new result block. This
 * implementation takes the starting pairwise index and pair size from the work
 * block and processes those pairs.
 *
 * @param block
 */
std::unique_ptr<EAbstractAnalyticBlock> Similarity::Serial::execute(const EAbstractAnalyticBlock* block)
{
    EDEBUG_FUNC(this,block);

    if ( ELog::isActive() )
    {
        ELog() << tr("Executing(serial) work index %1.\n").arg(block->index());
    }

    // cast block to work block
    const WorkBlock* workBlock {block->cast<WorkBlock>()};

    // initialize result block
    ResultBlock* resultBlock {new ResultBlock(workBlock->index())};

    // initialize workspace
    QVector<qint8> labels(_base->_input->sampleSize());

    // iterate through all pairs
    Pairwise::Index index {workBlock->start()};

    for ( int i = 0; i < workBlock->size(); ++i )
    {
        // fetch pairwise input data
        int numSamples = fetchPair(index, labels);

        // remove pre-clustering outliers
        if ( _base->_removePreOutliers )
        {
            numSamples = removeOutliers(index, numSamples, labels, 1, -7);
        }

        // compute clusters
        qint8 K {1};

        if ( _base->_clusMethod != ClusteringMethod::None )
        {
            K = _clusModel->compute(
                _expressions,
                index,
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
            numSamples = removeOutliers(index, numSamples, labels, K, -8);
        }

        // compute correlations
        QVector<float> correlations = _corrModel->compute(
            _expressions,
            index,
            K,
            labels,
            _base->_minSamples
        );

        // determine whether the pair contains any valid correlations
        bool valid = false;

        for ( qint8 k = 0; k < K; ++k )
        {
            // determine whether correlation is within thresholds
            float r = correlations[k];

            if ( !isnan(r) && _base->_minCorrelation <= abs(r) && abs(r) <= _base->_maxCorrelation )
            {
                valid = true;
                break;
            }
        }

        // save pair if it has any valid correlations
        if ( valid )
        {
            resultBlock->append(Pair {
                index,
                labels,
                correlations
            });
        }

        // increment to next pair
        ++index;
    }

    // return result block
    return unique_ptr<EAbstractAnalyticBlock>(resultBlock);
}



/*!
 * Compute the initial labels for a gene pair in an expression matrix. Samples
 * with missing values and samples that fall below the expression threshold are
 * labeled as such, all other samples are labeled as cluster 0. The number of
 * clean samples is returned.
 *
 * @param index
 * @param labels
 */
int Similarity::Serial::fetchPair(const Pairwise::Index& index, QVector<qint8>& labels)
{
    EDEBUG_FUNC(this,&index,&labels);

    // index into gene expressions
    const float *x = &_expressions[index.getX() * _base->_input->sampleSize()];
    const float *y = &_expressions[index.getY() * _base->_input->sampleSize()];

    // label the pairwise samples
    int numSamples = 0;

    for ( int i = 0; i < _base->_input->sampleSize(); ++i )
    {
        // label samples with missing values
        if ( std::isnan(x[i]) || std::isnan(y[i]) )
        {
            labels[i] = -9;
        }

        // label samples which are below the minimum expression threshold
        else if ( x[i] < _base->_minExpression || y[i] < _base->_minExpression )
        {
            labels[i] = -6;
        }

        // label samples which are above the maximum expression threshold
        else if ( x[i] > _base->_maxExpression || y[i] > _base->_maxExpression )
        {
            labels[i] = -6;
        }

        // label any remaining samples as cluster 0
        else
        {
            numSamples++;
            labels[i] = 0;
        }
    }

    // return number of clean samples
    return numSamples;
}



/*!
 * Remove outliers from a vector of pairwise data. Outliers are detected independently
 * on each axis using the Tukey method, and marked with the given marker. Only the
 * samples in the given cluster are used in outlier detection. For unclustered data,
 * all samples are labeled as 0, so a cluster value of 0 should be used.
 *
 * This function returns the number of clean samples remaining in the data array,
 * including samples in other clusters.
 *
 * @param x
 * @param y
 * @param labels
 * @param cluster
 * @param marker
 */
int Similarity::Serial::removeOutliersCluster(const float *x, const float *y, QVector<qint8>& labels, qint8 cluster, qint8 marker)
{
    EDEBUG_FUNC(this,x,y,&labels,cluster,marker);

    // extract samples from the given cluster into separate arrays
    QVector<float> x_sorted;
    QVector<float> y_sorted;

    x_sorted.reserve(labels.size());
    y_sorted.reserve(labels.size());

    for ( int i = 0; i < labels.size(); i++ )
    {
        if ( labels[i] == cluster )
        {
            x_sorted.append(x[i]);
            y_sorted.append(y[i]);
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

    // mark outliers
    int numSamples = 0;

    for ( int i = 0; i < labels.size(); i++ )
    {
        // mark samples in the given cluster that are outliers on either axis
        if ( labels[i] == cluster && (x[i] < T_x_min || T_x_max < x[i] || y[i] < T_y_min || T_y_max < y[i]) )
        {
            labels[i] = marker;
        }

        // count the number of remaining samples in the entire data array
        else if ( labels[i] >= 0 )
        {
            numSamples++;
        }
    }

    // return number of remaining samples
    return numSamples;
}



/*!
 * Perform outlier removal on each cluster in a parwise data array.
 *
 * @param index
 * @param numSamples
 * @param labels
 * @param clusterSize
 * @param marker
 */
int Similarity::Serial::removeOutliers(const Pairwise::Index& index, int numSamples, QVector<qint8>& labels, qint8 clusterSize, qint8 marker)
{
    EDEBUG_FUNC(this,&index,numSamples,&labels,clusterSize,marker);

    // index into gene expressions
    const float *x = &_expressions[index.getX() * _base->_input->sampleSize()];
    const float *y = &_expressions[index.getY() * _base->_input->sampleSize()];

    // do not perform post-clustering outlier removal if there is only one cluster
    if ( marker == -8 && clusterSize <= 1 )
    {
        return numSamples;
    }

    // perform outlier removal on each cluster
    for ( qint8 k = 0; k < clusterSize; ++k )
    {
        numSamples = removeOutliersCluster(x, y, labels, k, marker);
    }

    return numSamples;
}
