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
      NumPairs
      ,Expressions
      ,SampleSize
      ,InIndex
      ,InN
      ,InLabels
      ,InK
      ,Marker
      ,WorkX
      ,WorkY
   };
   explicit Outlier(::OpenCL::Program* program, QObject* parent = nullptr);
   ::OpenCL::Event execute(
      ::OpenCL::CommandQueue* queue,
      int globalWorkSize,
      int localWorkSize,
      int numPairs,
      ::OpenCL::Buffer<cl_float>* expressions,
      cl_int sampleSize,
      ::OpenCL::Buffer<cl_int2>* in_index,
      ::OpenCL::Buffer<cl_int>* in_N,
      ::OpenCL::Buffer<cl_char>* in_labels,
      ::OpenCL::Buffer<cl_char>* in_K,
      cl_char marker,
      ::OpenCL::Buffer<cl_float>* work_x,
      ::OpenCL::Buffer<cl_float>* work_y
   );
};



#endif
