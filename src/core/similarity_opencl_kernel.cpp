#include "similarity_opencl_kernel.h"



using namespace std;



/*!
 * Construct a new OpenCL kernel object with the given OpenCL program and qt parent.
 *
 * @param program
 * @param parent
 */
Similarity::OpenCL::Kernel::Kernel(::OpenCL::Program* program, QObject* parent):
    ::OpenCL::Kernel(program, "Similarity_compute", parent)
{
    EDEBUG_FUNC(this,program,parent);
}



/*!
 * Execute this kernel object's OpenCL kernel using the given OpenCL command
 * queue and kernel arguments, returning the OpenCL event associated with the
 * kernel execution.
 *
 * @param queue
 * @param globalWorkSize
 * @param localWorkSize
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
::OpenCL::Event Similarity::OpenCL::Kernel::execute(
    ::OpenCL::CommandQueue* queue,
    int globalWorkSize,
    int localWorkSize,
    cl_int clusMethod,
    cl_int corrMethod,
    cl_int removePreOutliers,
    cl_int removePostOutliers,
    cl_int numPairs,
    ::OpenCL::Buffer<cl_float>* expressions,
    cl_int sampleSize,
    ::OpenCL::Buffer<cl_int2>* in_index,
    cl_int minExpression,
    cl_int minSamples,
    cl_char minClusters,
    cl_char maxClusters,
    cl_int criterion,
    ::OpenCL::Buffer<cl_float>* work_x,
    ::OpenCL::Buffer<cl_float>* work_y,
    ::OpenCL::Buffer<cl_float2>* work_gmm_data,
    ::OpenCL::Buffer<cl_char>* work_gmm_labels,
    ::OpenCL::Buffer<cl_float>* work_gmm_pi,
    ::OpenCL::Buffer<cl_float2>* work_gmm_mu,
    ::OpenCL::Buffer<cl_float4>* work_gmm_sigma,
    ::OpenCL::Buffer<cl_float4>* work_gmm_sigmaInv,
    ::OpenCL::Buffer<cl_float>* work_gmm_normalizer,
    ::OpenCL::Buffer<cl_float2>* work_gmm_MP,
    ::OpenCL::Buffer<cl_int>* work_gmm_counts,
    ::OpenCL::Buffer<cl_float>* work_gmm_logpi,
    ::OpenCL::Buffer<cl_float>* work_gmm_gamma,
    ::OpenCL::Buffer<cl_char>* out_K,
    ::OpenCL::Buffer<cl_char>* out_labels,
    ::OpenCL::Buffer<cl_float>* out_correlations
)
{
    EDEBUG_FUNC(this,
        queue,
        globalWorkSize,
        localWorkSize,
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

    // acquire lock for this kernel
    Locker locker {lock()};

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

    // set work sizes
    setSizes(0, globalWorkSize, localWorkSize);

    // execute kernel
    return ::OpenCL::Kernel::execute(queue);
}
