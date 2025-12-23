#ifndef SIMILARITY_CUDA_WORKER_H
#define SIMILARITY_CUDA_WORKER_H
#include "similarity_cuda.h"
#include "similarity_cuda_kernel.h"



/*!
 * This class implements the CUDA worker of the similarity analytic.
 */
class Similarity::CUDA::Worker : public EAbstractAnalyticCUDAWorker
{
    Q_OBJECT
public:
    explicit Worker(Similarity* base, Similarity::CUDA* baseCuda, ::CUDA::Program* program);
    virtual std::unique_ptr<EAbstractAnalyticBlock> execute(const EAbstractAnalyticBlock* block) override final;
private:
    /*!
     * Pointer to the base analytic.
     */
    Similarity* _base;
    /*!
     * Pointer to the base CUDA object.
     */
    Similarity::CUDA* _baseCuda;
    /*!
     * Pointer to this worker's unique and private stream.
     */
    ::CUDA::Stream _stream;
    /*!
     * This worker's kernel.
     */
     CUDA::Kernel _kernel;
    /*!
     * Structure of this worker's buffers.
     */
    struct
    {
        ::CUDA::Buffer<int2> in_index;
        ::CUDA::Buffer<float> work_x;
        ::CUDA::Buffer<float> work_y;
        ::CUDA::Buffer<float2> work_gmm_data;
        ::CUDA::Buffer<qint8> work_gmm_labels;
        ::CUDA::Buffer<float> work_gmm_pi;
        ::CUDA::Buffer<float2> work_gmm_mu;
        ::CUDA::Buffer<float4> work_gmm_sigma;
        ::CUDA::Buffer<float4> work_gmm_sigmaInv;
        ::CUDA::Buffer<float> work_gmm_normalizer;
        ::CUDA::Buffer<float2> work_gmm_MP;
        ::CUDA::Buffer<int> work_gmm_counts;
        ::CUDA::Buffer<float> work_gmm_logpi;
        ::CUDA::Buffer<float> work_gmm_gamma;
        ::CUDA::Buffer<qint8> out_K;
        ::CUDA::Buffer<qint8> out_labels;
        ::CUDA::Buffer<float> out_correlations;
    } _buffers;
};



#endif
