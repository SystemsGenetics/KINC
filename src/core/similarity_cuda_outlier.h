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
      NumPairs
      ,Expressions
      ,SampleSize
      ,InIndex
      ,InN
      ,InLabels
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
      int numPairs,
      ::CUDA::Buffer<float>* expressions,
      int sampleSize,
      ::CUDA::Buffer<int2>* in_index,
      ::CUDA::Buffer<int>* in_N,
      ::CUDA::Buffer<qint8>* in_labels,
      ::CUDA::Buffer<qint8>* in_K,
      qint8 marker,
      ::CUDA::Buffer<float>* work_x,
      ::CUDA::Buffer<float>* work_y
   );
};



#endif
