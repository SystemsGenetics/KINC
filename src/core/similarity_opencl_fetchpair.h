#ifndef SIMILARITY_OPENCL_FETCHPAIR_H
#define SIMILARITY_OPENCL_FETCHPAIR_H
#include "similarity_opencl.h"



class Similarity::OpenCL::FetchPair : public ::OpenCL::Kernel
{
   Q_OBJECT
public:
   enum Argument
   {
      Expressions
      ,SampleSize
      ,InIndex
      ,MinExpression
      ,OutX
      ,OutN
      ,OutLabels
   };
   explicit FetchPair(::OpenCL::Program* program, QObject* parent = nullptr);
   ::OpenCL::Event execute(
      ::OpenCL::CommandQueue* queue,
      int kernelSize,
      ::OpenCL::Buffer<cl_float>* expressions,
      cl_int sampleSize,
      ::OpenCL::Buffer<cl_int2>* in_index,
      cl_int minExpression,
      ::OpenCL::Buffer<Pairwise::Vector2>* out_X,
      ::OpenCL::Buffer<cl_int>* out_N,
      ::OpenCL::Buffer<cl_char>* out_labels
   );
};



#endif
