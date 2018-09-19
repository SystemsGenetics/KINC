#include "similarity_opencl_worker.h"
#include "similarity_opencl_fetchpair.h"
#include "similarity_opencl_gmm.h"
#include "similarity_opencl_pearson.h"
#include "similarity_opencl_spearman.h"
#include "similarity_resultblock.h"
#include "similarity_workblock.h"
#include <ace/core/elog.h>
#include "pairwise_spearman.h"



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
   _queue(new ::OpenCL::CommandQueue(context, context->devices().first(), this))
{
   // initialize kernels
   _kernels.fetchPair = new OpenCL::FetchPair(program, this);
   _kernels.gmm = new OpenCL::GMM(program, this);
   _kernels.pearson = new OpenCL::Pearson(program, this);
   _kernels.spearman = new OpenCL::Spearman(program, this);

   // initialize buffers
   int kernelSize {_base->_kernelSize};
   int N {_base->_input->sampleSize()};
   int N_pow2 {Pairwise::Spearman::nextPower2(N)};
   int K {_base->_maxClusters};

   _buffers.in_index = ::OpenCL::Buffer<cl_int2>(context, 1 * kernelSize);
   _buffers.work_X = ::OpenCL::Buffer<Pairwise::Vector2>(context, N * kernelSize);
   _buffers.work_N = ::OpenCL::Buffer<cl_int>(context, 1 * kernelSize);
   _buffers.work_x = ::OpenCL::Buffer<cl_float>(context, N_pow2 * kernelSize);
   _buffers.work_y = ::OpenCL::Buffer<cl_float>(context, N_pow2 * kernelSize);
   _buffers.work_labels = ::OpenCL::Buffer<cl_char>(context, N * kernelSize);
   _buffers.work_components = ::OpenCL::Buffer<Pairwise::GMM::Component>(context, K * kernelSize);
   _buffers.work_MP = ::OpenCL::Buffer<Pairwise::Vector2>(context, K * kernelSize);
   _buffers.work_counts = ::OpenCL::Buffer<cl_int>(context, K * kernelSize);
   _buffers.work_logpi = ::OpenCL::Buffer<cl_float>(context, K * kernelSize);
   _buffers.work_loggamma = ::OpenCL::Buffer<cl_float>(context, N * K * kernelSize);
   _buffers.work_logGamma = ::OpenCL::Buffer<cl_float>(context, K * kernelSize);
   _buffers.work_rank = ::OpenCL::Buffer<cl_int>(context, N_pow2 * kernelSize);
   _buffers.out_K = ::OpenCL::Buffer<cl_char>(context, 1 * kernelSize);
   _buffers.out_labels = ::OpenCL::Buffer<cl_char>(context, N * kernelSize);
   _buffers.out_correlations = ::OpenCL::Buffer<cl_float>(context, K * kernelSize);
}






/*!
 * Read in the given work block, execute the algorithms necessary to produce
 * results using OpenCL acceleration, and save those results in a new result
 * block whose pointer is returned.
 *
 * @param block
 */
std::unique_ptr<EAbstractAnalytic::Block> Similarity::OpenCL::Worker::execute(const EAbstractAnalytic::Block* block)
{
   if ( ELog::isActive() )
   {
      ELog() << tr("Executing(OpenCL) work index %1.\n").arg(block->index());
   }

   // cast block to work block
   const WorkBlock* workBlock {block->cast<const WorkBlock>()};

   // initialize result block
   ResultBlock* resultBlock {new ResultBlock(workBlock->index(), workBlock->start())};

   // iterate through all pairs
   Pairwise::Index index {workBlock->start()};

   for ( int i = 0; i < workBlock->size(); i += _base->_kernelSize )
   {
      // write input buffers to device
      int steps {min(_base->_kernelSize, (int)workBlock->size() - i)};

      _buffers.in_index.mapWrite(_queue).wait();

      for ( int j = 0; j < steps; ++j )
      {
         _buffers.in_index[j] = { index.getX(), index.getY() };
         ++index;
      }

      for ( int j = steps; j < _base->_kernelSize; ++j )
      {
         _buffers.in_index[j] = { 0, 0 };
      }

      _buffers.in_index.unmap(_queue).wait();

      // execute fetch-pair kernel
      _kernels.fetchPair->execute(
         _queue,
         _base->_kernelSize,
         &_baseOpenCL->_expressions,
         _base->_input->sampleSize(),
         &_buffers.in_index,
         _base->_minExpression,
         &_buffers.work_X,
         &_buffers.work_N,
         &_buffers.out_labels
      ).wait();

      // execute clustering kernel
      if ( _base->_clusMethod == ClusteringMethod::GMM )
      {
         _kernels.gmm->execute(
            _queue,
            _base->_kernelSize,
            _base->_input->sampleSize(),
            _base->_minSamples,
            _base->_minClusters,
            _base->_maxClusters,
            _base->_criterion,
            _base->_removePreOutliers,
            _base->_removePostOutliers,
            &_buffers.work_X,
            &_buffers.work_N,
            &_buffers.work_x,
            &_buffers.work_y,
            &_buffers.work_labels,
            &_buffers.work_components,
            &_buffers.work_MP,
            &_buffers.work_counts,
            &_buffers.work_logpi,
            &_buffers.work_loggamma,
            &_buffers.work_logGamma,
            &_buffers.out_K,
            &_buffers.out_labels
         ).wait();
      }
      else
      {
         // set cluster size to 1 if clustering is disabled
         _buffers.out_K.mapWrite(_queue).wait();

         for ( int i = 0; i < _base->_kernelSize; ++i )
         {
            _buffers.out_K[i] = 1;
         }

         _buffers.out_K.unmap(_queue).wait();
      }

      // execute correlation kernel
      if ( _base->_corrMethod == CorrelationMethod::Pearson )
      {
         _kernels.pearson->execute(
            _queue,
            _base->_kernelSize,
            &_buffers.work_X,
            _base->_maxClusters,
            &_buffers.out_labels,
            _base->_input->sampleSize(),
            _base->_minSamples,
            &_buffers.out_correlations
         );
      }
      else if ( _base->_corrMethod == CorrelationMethod::Spearman )
      {
         _kernels.spearman->execute(
            _queue,
            _base->_kernelSize,
            &_buffers.work_X,
            _base->_maxClusters,
            &_buffers.out_labels,
            _base->_input->sampleSize(),
            _base->_minSamples,
            &_buffers.work_x,
            &_buffers.work_y,
            &_buffers.work_rank,
            &_buffers.out_correlations
         );
      }

      // read results from device
      auto e1 {_buffers.out_K.mapRead(_queue)};
      auto e2 {_buffers.out_labels.mapRead(_queue)};
      auto e3 {_buffers.out_correlations.mapRead(_queue)};

      e1.wait();
      e2.wait();
      e3.wait();

      // save results
      for ( int j = 0; j < steps; ++j )
      {
         // get pointers to the cluster labels and correlations for this pair
         const qint8 *labels = &_buffers.out_labels.at(j * _base->_input->sampleSize());
         const float *correlations = &_buffers.out_correlations.at(j * _base->_maxClusters);

         Pair pair;

         // save the number of clusters
         pair.K = _buffers.out_K.at(j);

         // save the cluster labels (if more than one cluster was found)
         if ( pair.K > 1 )
         {
            pair.labels = ResultBlock::makeVector(labels, _base->_input->sampleSize());
         }

         // save the correlations (if the pair was able to be processed)
         if ( pair.K > 0 )
         {
            pair.correlations = ResultBlock::makeVector(correlations, _base->_maxClusters);
         }

         resultBlock->append(pair);
      }

      auto e4 {_buffers.out_K.unmap(_queue)};
      auto e5 {_buffers.out_labels.unmap(_queue)};
      auto e6 {_buffers.out_correlations.unmap(_queue)};

      e4.wait();
      e5.wait();
      e6.wait();
   }

   // return result block
   return unique_ptr<EAbstractAnalytic::Block>(resultBlock);
}
