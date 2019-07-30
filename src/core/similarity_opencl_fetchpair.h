#ifndef SIMILARITY_OPENCL_FETCHPAIR_H
#define SIMILARITY_OPENCL_FETCHPAIR_H
#include "similarity_opencl.h"



/*!
 * This class implements the fetch-pair kernel for the similarity analytic. This
 * kernel takes a list of pairwise indices and computes the pairwise data, the
 * number of clean samples, and the initial sample labels for each pair.
 */
class Similarity::OpenCL::FetchPair : public ::OpenCL::Kernel
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
      ,MinExpression
      ,OutN
      ,OutLabels
   };
   explicit FetchPair(::OpenCL::Program* program, QObject* parent = nullptr);
   ::OpenCL::Event execute(
      ::OpenCL::CommandQueue* queue,
      int globalWorkSize,
      int localWorkSize,
      int numPairs,
      ::OpenCL::Buffer<cl_float>* expressions,
      cl_int sampleSize,
      ::OpenCL::Buffer<cl_int2>* in_index,
      cl_float minExpression,
      ::OpenCL::Buffer<cl_int>* out_N,
      ::OpenCL::Buffer<cl_char>* out_labels
   );
};



#endif
