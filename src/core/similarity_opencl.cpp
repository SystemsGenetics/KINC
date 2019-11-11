#include "similarity_opencl.h"
#include "similarity_opencl_worker.h"



using namespace std;



/*!
 * Construct a new OpenCL object with the given analytic as its parent.
 *
 * @param parent
 */
Similarity::OpenCL::OpenCL(Similarity* parent):
   EAbstractAnalyticOpenCL(parent),
   _base(parent)
{
   EDEBUG_FUNC(this,parent);
}



/*!
 * Create and return a new OpenCL worker for the analytic.
 */
std::unique_ptr<EAbstractAnalyticOpenCLWorker> Similarity::OpenCL::makeWorker()
{
   EDEBUG_FUNC(this);

   return unique_ptr<EAbstractAnalyticOpenCLWorker>(new Worker(_base, this, _context, _program));
}



/*!
 * Initializes all OpenCL resources used by this object's implementation.
 *
 * @param context
 */
void Similarity::OpenCL::initialize(::OpenCL::Context* context)
{
   EDEBUG_FUNC(this,context);

   // create list of opencl source files
   QStringList paths {
      ":/opencl/linalg.cl",
      ":/opencl/fetchpair.cl",
      ":/opencl/sort.cl",
      ":/opencl/outlier.cl",
      ":/opencl/gmm.cl",
      ":/opencl/pearson.cl",
      ":/opencl/spearman.cl",
      ":/opencl/similarity.cl"
   };

   // create program
   _context = context;
   _program = new ::OpenCL::Program(context, paths, this);

   // create command queue
   _queue = new ::OpenCL::CommandQueue(context, context->devices().first(), this);

   // create buffer for expression data
   std::vector<float> rawData = _base->_input->dumpRawData();
   _expressions = ::OpenCL::Buffer<cl_float>(context,rawData.size());

   // copy expression data to device
   _expressions.mapWrite(_queue).wait();

   memcpy(_expressions.data(), rawData.data(), rawData.size() * sizeof(float));

   _expressions.unmap(_queue).wait();
}
