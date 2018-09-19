#ifndef SIMILARITY_OPENCL_WORKER_H
#define SIMILARITY_OPENCL_WORKER_H
#include "similarity_opencl.h"
#include "pairwise_gmm.h"



/*!
 * This class implements the OpenCL worker of the similarity analytic.
 */
class Similarity::OpenCL::Worker : public EAbstractAnalytic::OpenCL::Worker
{
   Q_OBJECT
public:
   explicit Worker(Similarity* base, Similarity::OpenCL* baseOpenCL, ::OpenCL::Context* context, ::OpenCL::Program* program);
   virtual std::unique_ptr<EAbstractAnalytic::Block> execute(const EAbstractAnalytic::Block* block) override final;
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
      OpenCL::Pearson* pearson;
      OpenCL::Spearman* spearman;
   } _kernels;
   /*!
    * Structure of this worker's buffers.
    */
   struct
   {
      ::OpenCL::Buffer<cl_int2> in_index;
      ::OpenCL::Buffer<Pairwise::Vector2> work_X;
      ::OpenCL::Buffer<cl_int> work_N;
      ::OpenCL::Buffer<cl_float> work_x;
      ::OpenCL::Buffer<cl_float> work_y;
      ::OpenCL::Buffer<cl_char> work_labels;
      ::OpenCL::Buffer<Pairwise::GMM::Component> work_components;
      ::OpenCL::Buffer<Pairwise::Vector2> work_MP;
      ::OpenCL::Buffer<cl_int> work_counts;
      ::OpenCL::Buffer<cl_float> work_logpi;
      ::OpenCL::Buffer<cl_float> work_loggamma;
      ::OpenCL::Buffer<cl_float> work_logGamma;
      ::OpenCL::Buffer<cl_int> work_rank;
      ::OpenCL::Buffer<cl_char> out_K;
      ::OpenCL::Buffer<cl_char> out_labels;
      ::OpenCL::Buffer<cl_float> out_correlations;
   } _buffers;
};



#endif
