#ifndef SIMILARITY_OPENCL_H
#define SIMILARITY_OPENCL_H
#include <ace/core/openclxx.h>

#include "similarity.h"



class Similarity::OpenCL : public EAbstractAnalytic::OpenCL
{
   Q_OBJECT
public:
   class FetchPair;
   class GMM;
   class KMeans;
   class Pearson;
   class Spearman;
   class Worker;
   explicit OpenCL(Similarity* parent);
   virtual std::unique_ptr<EAbstractAnalytic::OpenCL::Worker> makeWorker() override final;
   virtual void initialize(::OpenCL::Context* context) override final;
private:
   Similarity* _base;
   ::OpenCL::Context* _context {nullptr};
   ::OpenCL::Program* _program {nullptr};
   ::OpenCL::CommandQueue* _queue {nullptr};

   ::OpenCL::Buffer<cl_float> _expressions;
};



#endif
