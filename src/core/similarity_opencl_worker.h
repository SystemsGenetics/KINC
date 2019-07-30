#ifndef SIMILARITY_OPENCL_WORKER_H
#define SIMILARITY_OPENCL_WORKER_H
#include "similarity_opencl.h"
#include "similarity_opencl_fetchpair.h"
#include "similarity_opencl_gmm.h"
#include "similarity_opencl_outlier.h"
#include "similarity_opencl_pearson.h"
#include "similarity_opencl_spearman.h"



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
    * Structure of this worker's kernels.
    */
   struct
   {
      OpenCL::FetchPair* fetchPair;
      OpenCL::GMM* gmm;
      OpenCL::Outlier* outlier;
      OpenCL::Pearson* pearson;
      OpenCL::Spearman* spearman;
   } _kernels;
   /*!
    * Structure of this worker's buffers.
    */
   struct
   {
      ::OpenCL::Buffer<cl_int2> in_index;
      ::OpenCL::Buffer<cl_int> work_N;
      ::OpenCL::Buffer<cl_float2> work_X;
      ::OpenCL::Buffer<cl_char> work_labels;
      ::OpenCL::Buffer<cl_component> work_components;
      ::OpenCL::Buffer<cl_float2> work_MP;
      ::OpenCL::Buffer<cl_int> work_counts;
      ::OpenCL::Buffer<cl_float> work_logpi;
      ::OpenCL::Buffer<cl_float> work_gamma;
      ::OpenCL::Buffer<cl_float> work_x;
      ::OpenCL::Buffer<cl_float> work_y;
      ::OpenCL::Buffer<cl_char> out_K;
      ::OpenCL::Buffer<cl_char> out_labels;
      ::OpenCL::Buffer<cl_float> out_correlations;
   } _buffers;
};



#endif
