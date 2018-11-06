#ifndef SIMILARITY_OPENCL_H
#define SIMILARITY_OPENCL_H
#include <ace/core/openclxx.h>

#include "similarity.h"



/*!
 * This class implements the base OpenCL class of the similarity analytic.
 */
class Similarity::OpenCL : public EAbstractAnalytic::OpenCL
{
   Q_OBJECT
public:
   class FetchPair;
   class GMM;
   class Outlier;
   class Pearson;
   class Spearman;
   class Worker;
   explicit OpenCL(Similarity* parent);
   virtual std::unique_ptr<EAbstractAnalytic::OpenCL::Worker> makeWorker() override final;
   virtual void initialize(::OpenCL::Context* context) override final;
private:
   /*!
    * Pointer to the base analytic for this object.
    */
   Similarity* _base;
   /*!
    * Pointer to this object's base OpenCL context used to create all other resources.
    */
   ::OpenCL::Context* _context {nullptr};
   /*!
    * Pointer to this object's OpenCL program.
    */
   ::OpenCL::Program* _program {nullptr};
   /*!
    * Pointer to this object's OpenCL command queue.
    */
   ::OpenCL::CommandQueue* _queue {nullptr};
   /*!
    * Pointer to this object's OpenCL buffer for the expression matrix.
    */
   ::OpenCL::Buffer<cl_float> _expressions;
};



#endif
