#include "similarity_cuda_worker.h"
#include "similarity_resultblock.h"
#include "similarity_workblock.h"
#include <ace/core/elog.h>



using namespace std;



/*!
 * Construct a new CUDA worker with the given parent analytic, CUDA object,
 * CUDA context, and CUDA program.
 *
 * @param base
 * @param baseCuda
 * @param program
 */
Similarity::CUDA::Worker::Worker(Similarity* base, Similarity::CUDA* baseCuda, ::CUDA::Program* program):
    _base(base),
    _baseCuda(baseCuda),
    _kernel(program)
{
    EDEBUG_FUNC(this,base,baseCuda,program);

    // initialize buffers
    int W {_base->_globalWorkSize};
    int N {_base->_input->sampleSize()};
    int N_pow2 {nextPower2(N)};
    int K {_base->_maxClusters};

    _buffers.in_index            = ::CUDA::Buffer<int2>   (W * 1);
    _buffers.work_x              = ::CUDA::Buffer<float>  (W * N_pow2, false);
    _buffers.work_y              = ::CUDA::Buffer<float>  (W * N_pow2, false);
    _buffers.work_gmm_data       = ::CUDA::Buffer<float2> (W * N, false);
    _buffers.work_gmm_labels     = ::CUDA::Buffer<qint8>  (W * N, false);
    _buffers.work_gmm_pi         = ::CUDA::Buffer<float>  (W * K, false);
    _buffers.work_gmm_mu         = ::CUDA::Buffer<float2> (W * K, false);
    _buffers.work_gmm_sigma      = ::CUDA::Buffer<float4> (W * K, false);
    _buffers.work_gmm_sigmaInv   = ::CUDA::Buffer<float4> (W * K, false);
    _buffers.work_gmm_normalizer = ::CUDA::Buffer<float>  (W * K, false);
    _buffers.work_gmm_MP         = ::CUDA::Buffer<float2> (W * K, false);
    _buffers.work_gmm_counts     = ::CUDA::Buffer<int>    (W * K, false);
    _buffers.work_gmm_logpi      = ::CUDA::Buffer<float>  (W * K, false);
    _buffers.work_gmm_gamma      = ::CUDA::Buffer<float>  (W * N * K, false);
    _buffers.out_K               = ::CUDA::Buffer<qint8>  (W * 1);
    _buffers.out_labels          = ::CUDA::Buffer<qint8>  (W * N);
    _buffers.out_correlations    = ::CUDA::Buffer<float>  (W * K);
}



/*!
 * Read in the given work block, execute the algorithms necessary to produce
 * results using CUDA acceleration, and save those results in a new result
 * block whose pointer is returned.
 *
 * @param block
 */
std::unique_ptr<EAbstractAnalyticBlock> Similarity::CUDA::Worker::execute(const EAbstractAnalyticBlock* block)
{
    EDEBUG_FUNC(this,block);

    if ( ELog::isActive() )
    {
        ELog() << tr("Executing(CUDA) work index %1.\n").arg(block->index());
    }

    // cast block to work block
    const WorkBlock* workBlock {block->cast<const WorkBlock>()};

    // initialize result block
    ResultBlock* resultBlock {new ResultBlock(workBlock->index())};

    // iterate through all pairs
    Pairwise::Index baseIndex {workBlock->start()};

    for ( int i = 0; i < workBlock->size(); i += _base->_globalWorkSize )
    {
        // initialize local index
        Pairwise::Index index {baseIndex};

        // write input buffers to device
        int numPairs {static_cast<int>(min(static_cast<qint64>(_base->_globalWorkSize), workBlock->size() - i))};

        for ( int j = 0; j < numPairs; ++j )
        {
            _buffers.in_index[j] = { index.getX(), index.getY() };
            ++index;
        }

        _buffers.in_index.write(_stream);

        // execute similiarity kernel
        _kernel.execute(
            _stream,
            _base->_globalWorkSize,
            _base->_localWorkSize,
            (int) _base->_clusMethod,
            (int) _base->_corrMethod,
            _base->_removePreOutliers,
            _base->_removePostOutliers,
            numPairs,
            &_baseCuda->_expressions,
            _base->_input->sampleSize(),
            &_buffers.in_index,
            _base->_minExpression,
            _base->_maxExpression,
            _base->_minSamples,
            _base->_minClusters,
            _base->_maxClusters,
            (int) _base->_criterion,
            &_buffers.work_x,
            &_buffers.work_y,
            &_buffers.work_gmm_data,
            &_buffers.work_gmm_labels,
            &_buffers.work_gmm_pi,
            &_buffers.work_gmm_mu,
            &_buffers.work_gmm_sigma,
            &_buffers.work_gmm_sigmaInv,
            &_buffers.work_gmm_normalizer,
            &_buffers.work_gmm_MP,
            &_buffers.work_gmm_counts,
            &_buffers.work_gmm_logpi,
            &_buffers.work_gmm_gamma,
            &_buffers.out_K,
            &_buffers.out_labels,
            &_buffers.out_correlations
        );

        // read results from device
        _buffers.out_K.read(_stream);
        _buffers.out_labels.read(_stream);
        _buffers.out_correlations.read(_stream);

        // wait for everything to finish
        _stream.wait();

        // save results
        index = baseIndex;

        for ( int j = 0; j < numPairs; ++j )
        {
            // extract output data for this pair
            qint8 K = _buffers.out_K.at(j);
            const qint8 *labels = &_buffers.out_labels.at(j * _base->_input->sampleSize());
            const float *correlations = &_buffers.out_correlations.at(j * _base->_maxClusters);

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
                    ResultBlock::makeVector(labels, _base->_input->sampleSize()),
                    ResultBlock::makeVector(correlations, K)
                });
            }

            // increment to next pair
            ++index;
        }

        // update base index
        baseIndex = index;
    }

    // return result block
    return unique_ptr<EAbstractAnalyticBlock>(resultBlock);
}
