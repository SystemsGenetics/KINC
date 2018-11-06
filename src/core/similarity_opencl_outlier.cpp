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
 * @param in_data
 * @param in_N
 * @param in_labels
 * @param sampleSize
 * @param in_K
 * @param marker
 * @param work_x
 * @param work_y
 */
::OpenCL::Event Similarity::OpenCL::Outlier::execute(
   ::OpenCL::CommandQueue* queue,
   int globalWorkSize,
   int localWorkSize,
   ::OpenCL::Buffer<cl_float2>* in_data,
   ::OpenCL::Buffer<cl_int>* in_N,
   ::OpenCL::Buffer<cl_char>* in_labels,
   cl_int sampleSize,
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
      in_data,
      in_N,
      in_labels,
      sampleSize,
      in_K,
      marker,
      work_x,
      work_y);

   // acquire lock for this kernel
   Locker locker {lock()};

   // set kernel arguments
   setArgument(GlobalWorkSize, globalWorkSize);
   setBuffer(InData, in_data);
   setBuffer(InN, in_N);
   setBuffer(InLabels, in_labels);
   setArgument(SampleSize, sampleSize);
   setBuffer(InK, in_K);
   setArgument(Marker, marker);
   setBuffer(WorkX, work_x);
   setBuffer(WorkY, work_y);

   // set work sizes
   if ( localWorkSize == 0 )
   {
      localWorkSize = min(globalWorkSize, maxWorkGroupSize(queue->device()));
   }

   int numWorkgroups = (globalWorkSize + localWorkSize - 1) / localWorkSize;

   setSizes(0, numWorkgroups * localWorkSize, localWorkSize);

   // execute kernel
   return ::OpenCL::Kernel::execute(queue);
}
