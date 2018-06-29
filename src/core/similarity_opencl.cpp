#include "similarity_opencl.h"
#include "similarity_opencl_worker.h"



using namespace std;






Similarity::OpenCL::OpenCL(Similarity* parent):
   EAbstractAnalytic::OpenCL(parent),
   _base(parent)
{
}






std::unique_ptr<EAbstractAnalytic::OpenCL::Worker> Similarity::OpenCL::makeWorker()
{
   return unique_ptr<EAbstractAnalytic::OpenCL::Worker>(new Worker(_base, this, _context, _program));
}






void Similarity::OpenCL::initialize(::OpenCL::Context* context)
{
   // create list of opencl source files
   QStringList paths {
      ":/opencl/linalg.cl",
      ":/opencl/fetchpair.cl",
      ":/opencl/sort.cl",
      ":/opencl/outlier.cl",
      ":/opencl/gmm.cl",
      ":/opencl/kmeans.cl",
      ":/opencl/pearson.cl",
      ":/opencl/spearman.cl"
   };

   // create program
   _context = context;
   _program = new ::OpenCL::Program(context, paths, this);

   // create command queue
   _queue = new ::OpenCL::CommandQueue(context, context->devices().first(), this);

   // create buffer for expression data
   _expressions = ::OpenCL::Buffer<cl_float>(context, _base->_input->getRawSize());

   unique_ptr<ExpressionMatrix::Expression> rawData(_base->_input->dumpRawData());
   ExpressionMatrix::Expression* rawDataRef {rawData.get()};

   // copy expression data to device
   _expressions.mapWrite(_queue).wait();

   for ( int i = 0; i < _base->_input->getRawSize(); ++i )
   {
      _expressions[i] = rawDataRef[i];
   }

   _expressions.unmap(_queue).wait();
}
