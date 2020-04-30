#include "similarity_opencl_worker.h"
#include "similarity_resultblock.h"
#include "similarity_workblock.h"
#include <ace/core/elog.h>



using namespace std;



/*!
 * Construct a new OpenCL worker with the given parent analytic, OpenCL object,
 * OpenCL context, and OpenCL program.
 *
 * @param base
 * @param baseOpenCL
 * @param context
 * @param program
 */
Similarity::OpenCL::Worker::Worker(Similarity* base, Similarity::OpenCL* baseOpenCL, ::OpenCL::Context* context, ::OpenCL::Program* program):
    _base(base),
    _baseOpenCL(baseOpenCL),
    _queue(new ::OpenCL::CommandQueue(context, context->devices().first(), this)),
    _kernel(new OpenCL::Kernel(program, this))
{
    EDEBUG_FUNC(this,base,baseOpenCL,context,program);

    // initialize buffers
    int W {_base->_globalWorkSize};
    int N {_base->_input->sampleSize()};
    int N_pow2 {nextPower2(N)};
    int K {_base->_maxClusters};

    _buffers.in_index            = ::OpenCL::Buffer<cl_int2>   (context, W * 1);
    _buffers.work_x              = ::OpenCL::Buffer<cl_float>  (context, W * N_pow2);
    _buffers.work_y              = ::OpenCL::Buffer<cl_float>  (context, W * N_pow2);
    _buffers.work_gmm_data       = ::OpenCL::Buffer<cl_float2> (context, W * N);
    _buffers.work_gmm_labels     = ::OpenCL::Buffer<cl_char>   (context, W * N);
    _buffers.work_gmm_pi         = ::OpenCL::Buffer<cl_float>  (context, W * K);
    _buffers.work_gmm_mu         = ::OpenCL::Buffer<cl_float2> (context, W * K);
    _buffers.work_gmm_sigma      = ::OpenCL::Buffer<cl_float4> (context, W * K);
    _buffers.work_gmm_sigmaInv   = ::OpenCL::Buffer<cl_float4> (context, W * K);
    _buffers.work_gmm_normalizer = ::OpenCL::Buffer<cl_float>  (context, W * K);
    _buffers.work_gmm_MP         = ::OpenCL::Buffer<cl_float2> (context, W * K);
    _buffers.work_gmm_counts     = ::OpenCL::Buffer<cl_int>    (context, W * K);
    _buffers.work_gmm_logpi      = ::OpenCL::Buffer<cl_float>  (context, W * K);
    _buffers.work_gmm_gamma      = ::OpenCL::Buffer<cl_float>  (context, W * N * K);
    _buffers.out_K               = ::OpenCL::Buffer<cl_char>   (context, W * 1);
    _buffers.out_labels          = ::OpenCL::Buffer<cl_char>   (context, W * N);
    _buffers.out_correlations    = ::OpenCL::Buffer<cl_float>  (context, W * K);
}



/*!
 * Read in the given work block, execute the algorithms necessary to produce
 * results using OpenCL acceleration, and save those results in a new result
 * block whose pointer is returned.
 *
 * @param block
 */
std::unique_ptr<EAbstractAnalyticBlock> Similarity::OpenCL::Worker::execute(const EAbstractAnalyticBlock* block)
{
    EDEBUG_FUNC(this,block);

    if ( ELog::isActive() )
    {
        ELog() << tr("Executing(OpenCL) work index %1.\n").arg(block->index());
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

        _buffers.in_index.mapWrite(_queue).wait();

        for ( int j = 0; j < numPairs; ++j )
        {
            _buffers.in_index[j] = { index.getX(), index.getY() };
            ++index;
        }

        _buffers.in_index.unmap(_queue);

        // execute similiarity kernel
        _kernel->execute(
            _queue,
            _base->_globalWorkSize,
            _base->_localWorkSize,
            (cl_int) _base->_clusMethod,
            (cl_int) _base->_corrMethod,
            _base->_removePreOutliers,
            _base->_removePostOutliers,
            numPairs,
            &_baseOpenCL->_expressions,
            _base->_input->sampleSize(),
            &_buffers.in_index,
            _base->_minExpression,
            _base->_maxExpression,
            _base->_minSamples,
            _base->_minClusters,
            _base->_maxClusters,
            (cl_int) _base->_criterion,
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
        _buffers.out_K.mapRead(_queue);
        _buffers.out_labels.mapRead(_queue);
        _buffers.out_correlations.mapRead(_queue);

        // wait for everything to finish
        _queue->wait();

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

        _buffers.out_K.unmap(_queue);
        _buffers.out_labels.unmap(_queue);
        _buffers.out_correlations.unmap(_queue);

        // update base index
        baseIndex = index;
    }

    // return result block
    return unique_ptr<EAbstractAnalyticBlock>(resultBlock);
}
