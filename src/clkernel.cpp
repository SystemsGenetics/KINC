#include "clkernel.h"



CLKernel::~CLKernel()
{
   if (_isAlive)
   {
      clReleaseKernel(_id);
   }
}



CLKernel::CLKernel(const CLKernel& copy):
   _isAlive(true),
   _id(copy._id),
   _did(copy._did),
   _dims(copy._dims)
{
   for (int i=0;i<_dims;++i)
   {
      _gSizes[i] = copy._gSizes[i];
      _lSizes[i] = copy._lSizes[i];
   }
}



CLKernel& CLKernel::operator=(const CLKernel& copy)
{
   if (_isAlive)
   {
      clReleaseKernel(_id);
   }
   _isAlive = true;
   _id = copy._id;
   _did = copy._did;
   _dims = copy._dims;
   for (int i=0;i<_dims;++i)
   {
      _gSizes[i] = copy._gSizes[i];
      _lSizes[i] = copy._lSizes[i];
   }
}



CLKernel::CLKernel(CLKernel&& move):
   _isAlive(true),
   _id(move._id),
   _did(move._did),
   _dims(move._dims)
{
   for (int i=0;i<_dims;++i)
   {
      _gSizes[i] = move._gSizes[i];
      _lSizes[i] = move._lSizes[i];
   }
   move._isAlive = false;
}



CLKernel& CLKernel::operator=(CLKernel&& move)
{
   if (_isAlive)
   {
      clReleaseKernel(_id);
   }
   _isAlive = true;
   _id = move._id;
   _did = move._did;
   _dims = move._dims;
   for (int i=0;i<_dims;++i)
   {
      _gSizes[i] = move._gSizes[i];
      _lSizes[i] = move._lSizes[i];
   }
   move._isAlive = false;
}



void CLKernel::set_swarm_dims(cl_uint dims)
{
   assert<NotAlive>(_isAlive,__FILE__,__LINE__);
   assert<TooManyDims>(dims<=_maxDims,__FILE__,__LINE__);
   _dims = dims;
}



void CLKernel::set_swarm_size(int dim, cl_uint gSize, cl_uint lSize)
{
   assert<NotAlive>(_isAlive,__FILE__,__LINE__);
   assert<DimOutOfRange>(dim<_dims,__FILE__,__LINE__);
   _gSizes[dim] = gSize;
   _lSizes[dim] = lSize;
}



size_t CLKernel::get_wg_size()
{
   assert<NotAlive>(_isAlive,__FILE__,__LINE__);
   size_t ret;
   cl_int err = clGetKernelWorkGroupInfo(_id,_did,CL_KERNEL_WORK_GROUP_SIZE,
                                         sizeof(size_t),&ret,NULL);
   assert<CannotGetInfo>(err==CL_SUCCESS,__FILE__,__LINE__,err);
   return ret;
}



size_t CLKernel::get_wg_multiple()
{
   assert<NotAlive>(_isAlive,__FILE__,__LINE__);
   size_t ret;
   cl_int err = clGetKernelWorkGroupInfo(_id,_did,
                                   CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                         sizeof(size_t),&ret,NULL);
   assert<CannotGetInfo>(err==CL_SUCCESS,__FILE__,__LINE__,err);
   return ret;
}



CLKernel::CLKernel(cl_kernel id, cl_device_id did):
   _isAlive(true),
   _id(id),
   _did(did)
{}
