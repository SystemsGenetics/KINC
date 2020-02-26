#ifndef SIMILARITY_CUDA_KERNEL_H
#define SIMILARITY_CUDA_KERNEL_H
#include "similarity_cuda.h"



/*!
 * This class implements the CUDA kernel for the similarity analytic.
 */
class Similarity::CUDA::Kernel : public ::CUDA::Kernel
{
public:
    /*!
     * Defines the arguments passed to the CUDA kernel.
     */
    enum Argument
    {
        ClusMethod
        ,CorrMethod
        ,RemovePreOutliers
        ,RemovePostOutliers
        ,NumPairs
        ,Expressions
        ,SampleSize
        ,InIndex
        ,MinExpression
        ,MinSamples
        ,MinClusters
        ,MaxClusters
        ,Criterion
        ,WorkX
        ,WorkY
        ,WorkGmmData
        ,WorkGmmLabels
        ,WorkGmmPi
        ,WorkGmmMu
        ,WorkGmmSigma
        ,WorkGmmSigmaInv
        ,WorkGmmNormalizer
        ,WorkGmmMP
        ,WorkGmmCounts
        ,WorkGmmLogPi
        ,WorkGmmGamma
        ,OutK
        ,OutLabels
        ,OutCorrelations
    };
    explicit Kernel(::CUDA::Program* program);
    ::CUDA::Event execute(
        const ::CUDA::Stream& stream,
        int occupancy,
        int blockSize,
        int clusMethod,
        int corrMethod,
        bool removePreOutliers,
        bool removePostOutliers,
        int numPairs,
        ::CUDA::Buffer<float>* expressions,
        int sampleSize,
        ::CUDA::Buffer<int2>* in_index,
        int minExpression,
        int minSamples,
        char minClusters,
        char maxClusters,
        int criterion,
        ::CUDA::Buffer<float>* work_x,
        ::CUDA::Buffer<float>* work_y,
        ::CUDA::Buffer<float2>* work_gmm_data,
        ::CUDA::Buffer<qint8>* work_gmm_labels,
        ::CUDA::Buffer<float>* work_gmm_pi,
        ::CUDA::Buffer<float2>* work_gmm_mu,
        ::CUDA::Buffer<float4>* work_gmm_sigma,
        ::CUDA::Buffer<float4>* work_gmm_sigmaInv,
        ::CUDA::Buffer<float>* work_gmm_normalizer,
        ::CUDA::Buffer<float2>* work_gmm_MP,
        ::CUDA::Buffer<int>* work_gmm_counts,
        ::CUDA::Buffer<float>* work_gmm_logpi,
        ::CUDA::Buffer<float>* work_gmm_gamma,
        ::CUDA::Buffer<qint8>* out_K,
        ::CUDA::Buffer<qint8>* out_labels,
        ::CUDA::Buffer<float>* out_correlations
    );
};



#endif
