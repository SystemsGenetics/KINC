#ifndef SIMILARITY_OPENCL_OUTLIER_H
#define SIMILARITY_OPENCL_OUTLIER_H
#include "similarity_opencl.h"



/*!
 * This class implements the outlier removal kernel for the similarity analytic.
 */
class Similarity::OpenCL::Outlier : public ::OpenCL::Kernel
{
   Q_OBJECT
public:
   /*!
    * Defines the arguments passed to the OpenCL kernel.
    */
   enum Argument
   {
      GlobalWorkSize
      ,InData
      ,InN
      ,InLabels
      ,SampleSize
      ,InK
      ,Marker
      ,WorkXY
   };
   explicit Outlier(::OpenCL::Program* program, QObject* parent = nullptr);
   ::OpenCL::Event execute(
      ::OpenCL::CommandQueue* queue,
      int globalWorkSize,
      int localWorkSize,
      ::OpenCL::Buffer<cl_float2>* in_data,
      ::OpenCL::Buffer<cl_int>* in_N,
      ::OpenCL::Buffer<cl_char>* in_labels,
      cl_int sampleSize,
      ::OpenCL::Buffer<cl_char>* in_K,
      cl_char marker,
      ::OpenCL::Buffer<cl_float>* work_xy
   );
};



#endif
