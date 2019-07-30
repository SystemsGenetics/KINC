#ifndef SIMILARITY_CUDA_FETCHPAIR_H
#define SIMILARITY_CUDA_FETCHPAIR_H
#include "similarity_cuda.h"



/*!
 * This class implements the fetch-pair kernel for the similarity analytic. This
 * kernel takes a list of pairwise indices and computes the pairwise data, the
 * number of clean samples, and the initial sample labels for each pair.
 */
class Similarity::CUDA::FetchPair : public ::CUDA::Kernel
{
public:
   /*!
    * Defines the arguments passed to the OpenCL kernel.
    */
   enum Argument
   {
      NumPairs
      ,Expressions
      ,SampleSize
      ,InIndex
      ,MinExpression
      ,OutN
      ,OutLabels
   };
   explicit FetchPair(::CUDA::Program* program);
   ::CUDA::Event execute(
      const ::CUDA::Stream& stream,
      int globalWorkSize,
      int localWorkSize,
      int numPairs,
      ::CUDA::Buffer<float>* expressions,
      int sampleSize,
      ::CUDA::Buffer<int2>* in_index,
      float minExpression,
      ::CUDA::Buffer<int>* out_N,
      ::CUDA::Buffer<qint8>* out_labels
   );
};



#endif
