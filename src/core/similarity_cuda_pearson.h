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
      GlobalWorkSize
      ,InData
      ,ClusterSize
      ,InLabels
      ,SampleSize
      ,MinSamples
      ,OutCorrelations
   };
   explicit Pearson(::CUDA::Program* program);
   ::CUDA::Event execute(
      const ::CUDA::Stream& stream,
      int globalWorkSize,
      int localWorkSize,
      ::CUDA::Buffer<float2>* in_data,
      char clusterSize,
      ::CUDA::Buffer<qint8>* in_labels,
      int sampleSize,
      int minSamples,
      ::CUDA::Buffer<float>* out_correlations
   );
};



#endif
