#include "similarity_cuda_kernel.h"

#define MAJOR_VERSION KINC_MAJOR_VERSION
#define MINOR_VERSION KINC_MINOR_VERSION
#define REVISION KINC_REVISION
#include <ace/core/ace_settings.h>



using namespace std;



/*!
 * Construct a new CUDA kernel object with the given CUDA program.
 *
 * @param program
 */
Similarity::CUDA::Kernel::Kernel(::CUDA::Program* program):
    ::CUDA::Kernel(program, "Similarity_compute")
{
    EDEBUG_FUNC(this,program);
}



/*!
 * Execute this kernel object's CUDA kernel using the given CUDA stream
 * and kernel arguments, returning the CUDA event associated with the
 * kernel execution.
 *
 * @param stream
 * @param occupancy
 * @param blockSize
 * @param clusMethod
 * @param corrMethod
 * @param removePreOutliers
 * @param removePostOutliers
 * @param numPairs
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param minExpression
 * @param minSamples
 * @param minClusters
 * @param maxClusters
 * @param criterion
 * @param work_x
 * @param work_y
 * @param work_gmm_data
 * @param work_gmm_labels
 * @param work_gmm_pi
 * @param work_gmm_mu
 * @param work_gmm_sigma
 * @param work_gmm_sigmaInv
 * @param work_gmm_normalizer
 * @param work_gmm_MP
 * @param work_gmm_counts
 * @param work_gmm_logpi
 * @param work_gmm_gamma
 * @param out_K
 * @param out_labels
 * @param out_correlations
 */
::CUDA::Event Similarity::CUDA::Kernel::execute(
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
)
{
    EDEBUG_FUNC(this,
        &stream,
        occupancy,
        blockSize,
        clusMethod,
        corrMethod,
        removePreOutliers,
        removePostOutliers,
        numPairs,
        expressions,
        sampleSize,
        in_index,
        minExpression,
        minSamples,
        minClusters,
        maxClusters,
        criterion,
        work_x,
        work_y,
        work_gmm_data,
        work_gmm_labels,
        work_gmm_pi,
        work_gmm_mu,
        work_gmm_sigma,
        work_gmm_sigmaInv,
        work_gmm_normalizer,
        work_gmm_MP,
        work_gmm_counts,
        work_gmm_logpi,
        work_gmm_gamma,
        out_K,
        out_labels,
        out_correlations);

    // set kernel arguments
    setArgument(ClusMethod, clusMethod);
    setArgument(CorrMethod, corrMethod);
    setArgument(RemovePreOutliers, removePreOutliers);
    setArgument(RemovePostOutliers, removePostOutliers);
    setArgument(NumPairs, numPairs);
    setBuffer(Expressions, expressions);
    setArgument(SampleSize, sampleSize);
    setBuffer(InIndex, in_index);
    setArgument(MinExpression, minExpression);
    setArgument(MinSamples, minSamples);
    setArgument(MinClusters, minClusters);
    setArgument(MaxClusters, maxClusters);
    setArgument(Criterion, criterion);
    setBuffer(WorkX, work_x);
    setBuffer(WorkY, work_y);
    setBuffer(WorkGmmData, work_gmm_data);
    setBuffer(WorkGmmLabels, work_gmm_labels);
    setBuffer(WorkGmmPi, work_gmm_pi);
    setBuffer(WorkGmmMu, work_gmm_mu);
    setBuffer(WorkGmmSigma, work_gmm_sigma);
    setBuffer(WorkGmmSigmaInv, work_gmm_sigmaInv);
    setBuffer(WorkGmmNormalizer, work_gmm_normalizer);
    setBuffer(WorkGmmMP, work_gmm_MP);
    setBuffer(WorkGmmCounts, work_gmm_counts);
    setBuffer(WorkGmmLogPi, work_gmm_logpi);
    setBuffer(WorkGmmGamma, work_gmm_gamma);
    setBuffer(OutK, out_K);
    setBuffer(OutLabels, out_labels);
    setBuffer(OutCorrelations, out_correlations);

    // select a grid size to achieve a given occupancy
    auto device {Ace::Settings::instance().cudaDevicePointer()};
    int numSMs {device->getAttribute(CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT)};
    int numBlocksPerSM {getMaxActiveBlocksPerMultiprocessor(blockSize)};

    int gridSize = numSMs * numBlocksPerSM * occupancy / 100;

    // set work sizes
    setSizes(gridSize, blockSize);

    // execute kernel
    return ::CUDA::Kernel::execute(stream);
}
