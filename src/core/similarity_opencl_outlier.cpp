#include "similarity_opencl_outlier.h"



using namespace std;






/*!
 * Construct a new Outlier kernel object with the given OpenCL program and qt parent.
 *
 * @param program
 * @param parent
 */
Similarity::OpenCL::Outlier::Outlier(::OpenCL::Program* program, QObject* parent):
   ::OpenCL::Kernel(program, "removeOutliers", parent)
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
 * @param numPairs
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param in_N
 * @param in_labels
 * @param in_K
 * @param marker
 * @param work_x
 * @param work_y
 */
::OpenCL::Event Similarity::OpenCL::Outlier::execute(
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
)
{
   EDEBUG_FUNC(this,
      queue,
      globalWorkSize,
      localWorkSize,
      numPairs,
      expressions,
      sampleSize,
      in_index,
      in_N,
      in_labels,
      in_K,
      marker,
      work_x,
      work_y);

   // acquire lock for this kernel
   Locker locker {lock()};

   // set kernel arguments
   setArgument(NumPairs, numPairs);
   setBuffer(Expressions, expressions);
   setArgument(SampleSize, sampleSize);
   setBuffer(InIndex, in_index);
   setBuffer(InN, in_N);
   setBuffer(InLabels, in_labels);
   setBuffer(InK, in_K);
   setArgument(Marker, marker);
   setBuffer(WorkX, work_x);
   setBuffer(WorkY, work_y);

   // set work sizes
   setSizes(0, globalWorkSize / localWorkSize, localWorkSize);

   // execute kernel
   return ::OpenCL::Kernel::execute(queue);
}
