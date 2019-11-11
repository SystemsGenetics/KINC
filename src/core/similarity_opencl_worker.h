#ifndef SIMILARITY_OPENCL_WORKER_H
#define SIMILARITY_OPENCL_WORKER_H
#include "similarity_opencl.h"
#include "similarity_opencl_kernel.h"



/*!
 * This class implements the OpenCL worker of the similarity analytic.
 */
class Similarity::OpenCL::Worker : public EAbstractAnalyticOpenCLWorker
{
   Q_OBJECT
public:
   explicit Worker(Similarity* base, Similarity::OpenCL* baseOpenCL, ::OpenCL::Context* context, ::OpenCL::Program* program);
   virtual std::unique_ptr<EAbstractAnalyticBlock> execute(const EAbstractAnalyticBlock* block) override final;
private:
   /*!
    * Pointer to the base analytic.
    */
   Similarity* _base;
   /*!
    * Pointer to the base OpenCL object.
    */
   Similarity::OpenCL* _baseOpenCL;
   /*!
    * Pointer to this worker's unique and private command queue.
    */
   ::OpenCL::CommandQueue* _queue;
   /*!
    * This worker's kernel.
    */
   OpenCL::Kernel* _kernel;
   /*!
    * Structure of this worker's buffers.
    */
   struct
   {
      ::OpenCL::Buffer<cl_int2> in_index;
      ::OpenCL::Buffer<cl_float> work_x;
      ::OpenCL::Buffer<cl_float> work_y;
      ::OpenCL::Buffer<cl_float2> work_gmm_data;
      ::OpenCL::Buffer<cl_char> work_gmm_labels;
      ::OpenCL::Buffer<cl_float> work_gmm_pi;
      ::OpenCL::Buffer<cl_float2> work_gmm_mu;
      ::OpenCL::Buffer<cl_float4> work_gmm_sigma;
      ::OpenCL::Buffer<cl_float4> work_gmm_sigmaInv;
      ::OpenCL::Buffer<cl_float> work_gmm_normalizer;
      ::OpenCL::Buffer<cl_float2> work_gmm_MP;
      ::OpenCL::Buffer<cl_int> work_gmm_counts;
      ::OpenCL::Buffer<cl_float> work_gmm_logpi;
      ::OpenCL::Buffer<cl_float> work_gmm_gamma;
      ::OpenCL::Buffer<cl_char> out_K;
      ::OpenCL::Buffer<cl_char> out_labels;
      ::OpenCL::Buffer<cl_float> out_correlations;
   } _buffers;
};



#endif
