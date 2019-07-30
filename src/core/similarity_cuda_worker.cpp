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
   _kernels({
      .fetchPair = CUDA::FetchPair(program),
      .gmm = CUDA::GMM(program),
      .outlier = CUDA::Outlier(program),
      .pearson = CUDA::Pearson(program),
      .spearman = CUDA::Spearman(program)
   })
{
   EDEBUG_FUNC(this,base,baseCuda,program);

   // initialize buffers
   int W {_base->_globalWorkSize};
   int N {_base->_input->sampleSize()};
   int N_pow2 {nextPower2(N)};
   int K {_base->_maxClusters};

   _buffers.in_index = ::CUDA::Buffer<int2>(1 * W);
   _buffers.work_N = ::CUDA::Buffer<int>(1 * W, false);
   _buffers.work_xy = ::CUDA::Buffer<float>(2 * N_pow2 * W, false);
   _buffers.work_labels = ::CUDA::Buffer<qint8>(N * W, false);
   _buffers.work_components = ::CUDA::Buffer<cu_component>(K * W, false);
   _buffers.work_MP = ::CUDA::Buffer<float2>(K * W, false);
   _buffers.work_counts = ::CUDA::Buffer<int>(K * W, false);
   _buffers.work_logpi = ::CUDA::Buffer<float>(K * W, false);
   _buffers.work_gamma = ::CUDA::Buffer<float>(N * K * W, false);
   _buffers.out_K = ::CUDA::Buffer<qint8>(1 * W);
   _buffers.out_labels = ::CUDA::Buffer<qint8>(N * W);
   _buffers.out_correlations = ::CUDA::Buffer<float>(K * W);
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
   ResultBlock* resultBlock {new ResultBlock(workBlock->index(), workBlock->start())};

   // iterate through all pairs
   Pairwise::Index index {workBlock->start()};

   for ( int i = 0; i < workBlock->size(); i += _base->_globalWorkSize )
   {
      // write input buffers to device
      int globalWorkSize {static_cast<int>(min(static_cast<qint64>(_base->_globalWorkSize), workBlock->size() - i))};

      for ( int j = 0; j < globalWorkSize; ++j )
      {
         _buffers.in_index[j] = { index.getX(), index.getY() };
         ++index;
      }

      _buffers.in_index.write(_stream);

      // execute fetch-pair kernel
      _kernels.fetchPair.execute(
         _stream,
         globalWorkSize,
         _base->_localWorkSize,
         &_baseCuda->_expressions,
         _base->_input->sampleSize(),
         &_buffers.in_index,
         _base->_minExpression,
         &_buffers.work_N,
         &_buffers.out_labels
      );

      // execute outlier kernel (pre-clustering)
      if ( _base->_removePreOutliers )
      {
         _kernels.outlier.execute(
            _stream,
            globalWorkSize,
            _base->_localWorkSize,
            &_baseCuda->_expressions,
            _base->_input->sampleSize(),
            &_buffers.in_index,
            &_buffers.work_N,
            &_buffers.out_labels,
            &_buffers.out_K,
            -7,
            &_buffers.work_xy
         );
      }

      // execute clustering kernel
      if ( _base->_clusMethod == ClusteringMethod::GMM )
      {
         _kernels.gmm.execute(
            _stream,
            globalWorkSize,
            _base->_localWorkSize,
            &_baseCuda->_expressions,
            _base->_input->sampleSize(),
            &_buffers.in_index,
            _base->_minSamples,
            _base->_minClusters,
            _base->_maxClusters,
            (int) _base->_criterion,
            &_buffers.work_xy,
            &_buffers.work_N,
            &_buffers.work_labels,
            &_buffers.work_components,
            &_buffers.work_MP,
            &_buffers.work_counts,
            &_buffers.work_logpi,
            &_buffers.work_gamma,
            &_buffers.out_K,
            &_buffers.out_labels
         );
      }
      else
      {
         // set cluster size to 1 if clustering is disabled
         for ( int i = 0; i < globalWorkSize; ++i )
         {
            _buffers.out_K[i] = 1;
         }

         _buffers.out_K.write(_stream);
      }

      // execute outlier kernel (post-clustering)
      if ( _base->_removePostOutliers )
      {
         _kernels.outlier.execute(
            _stream,
            globalWorkSize,
            _base->_localWorkSize,
            &_baseCuda->_expressions,
            _base->_input->sampleSize(),
            &_buffers.in_index,
            &_buffers.work_N,
            &_buffers.out_labels,
            &_buffers.out_K,
            -8,
            &_buffers.work_xy
         );
      }

      // execute correlation kernel
      if ( _base->_corrMethod == CorrelationMethod::Pearson )
      {
         _kernels.pearson.execute(
            _stream,
            globalWorkSize,
            _base->_localWorkSize,
            &_baseCuda->_expressions,
            _base->_input->sampleSize(),
            &_buffers.in_index,
            _base->_maxClusters,
            &_buffers.out_labels,
            _base->_minSamples,
            &_buffers.out_correlations
         );
      }
      else if ( _base->_corrMethod == CorrelationMethod::Spearman )
      {
         _kernels.spearman.execute(
            _stream,
            globalWorkSize,
            _base->_localWorkSize,
            &_baseCuda->_expressions,
            _base->_input->sampleSize(),
            &_buffers.in_index,
            _base->_maxClusters,
            &_buffers.out_labels,
            _base->_minSamples,
            &_buffers.work_xy,
            &_buffers.out_correlations
         );
      }

      // read results from device
      _buffers.out_K.read(_stream);
      _buffers.out_labels.read(_stream);
      _buffers.out_correlations.read(_stream);

      // wait for everything to finish
      _stream.wait();

      // save results
      for ( int j = 0; j < globalWorkSize; ++j )
      {
         // get pointers to the cluster labels and correlations for this pair
         const qint8 *labels = &_buffers.out_labels.at(j * _base->_input->sampleSize());
         const float *correlations = &_buffers.out_correlations.at(j * _base->_maxClusters);

         Pair pair;

         // save the number of clusters
         pair.K = _buffers.out_K.at(j);

         // save the cluster labels and correlations (if the pair was able to be processed)
         if ( pair.K > 0 )
         {
            pair.labels = ResultBlock::makeVector(labels, _base->_input->sampleSize());
            pair.correlations = ResultBlock::makeVector(correlations, _base->_maxClusters);
         }

         resultBlock->append(pair);
      }
   }

   // return result block
   return unique_ptr<EAbstractAnalyticBlock>(resultBlock);
}
