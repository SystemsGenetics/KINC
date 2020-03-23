#ifndef SIMILARITY_CUDA_H
#define SIMILARITY_CUDA_H
#include <ace/core/cudaxx.h>

#include "similarity.h"



/*!
 * This class implements the base CUDA class of the similarity analytic.
 */
class Similarity::CUDA : public EAbstractAnalyticCUDA
{
    Q_OBJECT
public:
    class Kernel;
    class Worker;
    explicit CUDA(Similarity* parent);
    virtual std::unique_ptr<EAbstractAnalyticCUDAWorker> makeWorker() override final;
    virtual void initialize() override final;
private:
    /*!
     * Pointer to the base analytic for this object.
     */
    Similarity* _base;
    /*!
     * Pointer to this object's CUDA program.
     */
    ::CUDA::Program* _program {nullptr};
    /*!
     * CUDA buffer for this object's expression matrix.
     */
    ::CUDA::Buffer<float> _expressions;
};



#endif
