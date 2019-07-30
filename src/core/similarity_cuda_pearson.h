#ifndef SIMILARITY_CUDA_PEARSON_H
#define SIMILARITY_CUDA_PEARSON_H
#include "similarity_cuda.h"



/*!
 * This class implements the Pearson kernel for the similarity analytic. This
 * kernel takes a list of pairwise data arrays (with cluster labels) and computes
 * the Pearson correlation for each cluster in each pair.
 */
class Similarity::CUDA::Pearson : public ::CUDA::Kernel
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
      ,ClusterSize
      ,InLabels
      ,MinSamples
      ,OutCorrelations
   };
   explicit Pearson(::CUDA::Program* program);
   ::CUDA::Event execute(
      const ::CUDA::Stream& stream,
      int globalWorkSize,
      int localWorkSize,
      int numPairs,
      ::CUDA::Buffer<float>* expressions,
      int sampleSize,
      ::CUDA::Buffer<int2>* in_index,
      char clusterSize,
      ::CUDA::Buffer<qint8>* in_labels,
      int minSamples,
      ::CUDA::Buffer<float>* out_correlations
   );
};



#endif
