#ifndef SIMILARITY_CUDA_OUTLIER_H
#define SIMILARITY_CUDA_OUTLIER_H
#include "similarity_cuda.h"



/*!
 * This class implements the outlier removal kernel for the similarity analytic.
 */
class Similarity::CUDA::Outlier : public ::CUDA::Kernel
{ 
public:
   /*!
    * Defines the arguments passed to the CUDA kernel.
    */
   enum Argument
   {
      GlobalWorkSize
      ,InData
      ,InN
      ,InLabels
      ,SampleSize
      ,InK
      ,Marker
      ,WorkX
      ,WorkY
   };
   explicit Outlier(::CUDA::Program* program);
   ::CUDA::Event execute(
      const ::CUDA::Stream& stream,
      int globalWorkSize,
      int localWorkSize,
      ::CUDA::Buffer<float2>* in_data,
      ::CUDA::Buffer<int>* in_N,
      ::CUDA::Buffer<qint8>* in_labels,
      int sampleSize,
      ::CUDA::Buffer<qint8>* in_K,
      qint8 marker,
      ::CUDA::Buffer<float>* work_x,
      ::CUDA::Buffer<float>* work_y
   );
};



#endif
