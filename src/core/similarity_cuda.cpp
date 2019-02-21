#include "similarity_cuda.h"
#include "similarity_cuda_worker.h"



using namespace std;






/*!
 * Construct a new CUDA object with the given analytic as its parent.
 *
 * @param parent
 */
Similarity::CUDA::CUDA(Similarity* parent):
   EAbstractAnalytic::CUDA(parent),
   _base(parent)
{
   EDEBUG_FUNC(this,parent);
}






/*!
 * Create and return a new CUDA worker for the analytic.
 */
std::unique_ptr<EAbstractAnalytic::CUDA::Worker> Similarity::CUDA::makeWorker()
{
   EDEBUG_FUNC(this);

   return unique_ptr<EAbstractAnalytic::CUDA::Worker>(new Worker(_base, this, _program));
}






/*!
 * Initializes all CUDA resources used by this object's implementation.
 */
void Similarity::CUDA::initialize()
{
   EDEBUG_FUNC(this);

   // create list of cuda source files
   QStringList paths {
      ":/cuda/linalg.cu",
      ":/cuda/fetchpair.cu",
      ":/cuda/sort.cu",
      ":/cuda/outlier.cu",
      ":/cuda/gmm.cu",
      ":/cuda/pearson.cu",
      ":/cuda/spearman.cu"
   };

   // create program
   _program = new ::CUDA::Program(paths, this);

   // create buffer for expression data
   std::vector<float> rawData = _base->_input->dumpRawData();
   _expressions = ::CUDA::Buffer<float>(rawData.size());

   // copy expression data to device
   for ( size_t i = 0; i < rawData.size(); ++i )
   {
      _expressions[i] = rawData[i];
   }

   _expressions.write().wait();
}
