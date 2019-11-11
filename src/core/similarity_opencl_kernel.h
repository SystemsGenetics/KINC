#ifndef SIMILARITY_OPENCL_KERNEL_H
#define SIMILARITY_OPENCL_KERNEL_H
#include "similarity_opencl.h"



/*!
 * This class implements the OpenCL kernel for the similarity analytic.
 */
class Similarity::OpenCL::Kernel : public ::OpenCL::Kernel
{
public:
    /*!
     * Defines the arguments passed to the OpenCL kernel.
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
    explicit Kernel(::OpenCL::Program* program, QObject* parent = nullptr);
    ::OpenCL::Event execute(
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
    );
};



#endif
