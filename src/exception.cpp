#include "exception.h"



Exception::Exception(const char* file, int line, const char* what):
   _file(file),
   _line(line),
   _what(what)
{}



const char* Exception::file()
{
   return _file;
}



const int Exception::line()
{
   return _line;
}



const char* Exception::what()
{
   return _what;
}



DataException::DataException(const char* file, int line, Data* who, const char* what):
   Exception(file,line,what),
   _who(who)
{}



Data* DataException::who()
{
   return _who;
}



AnalyticException::AnalyticException(const char* file, int line, Analytic* who, const char* what):
   Exception(file,line,what),
   _who(who)
{}



Analytic* AnalyticException::who()
{
   return _who;
}



SystemError::SystemError(const char* file, int line, const char* system):
   Exception(file,line,"SystemError"),
   _system(system)
{}



const char* SystemError::system()
{
   return _system;
}



// List of all possible OpenCL errors.
const char* OpenCLError::c_clDescErrors[] = {
   "CL_SUCCESS",
   "CL_DEVICE_NOT_FOUND",
   "CL_DEVICE_NOT_AVAILABLE",
   "CL_COMPILER_NOT_AVAILABLE",
   "CL_MEM_OBJECT_ALLOCATION_FAILURE",
   "CL_OUT_OF_RESOURCES",
   "CL_OUT_OF_HOST_MEMORY",
   "CL_PROFILING_INFO_NOT_AVAILABLE",
   "CL_MEM_COPY_OVERLAP",
   "CL_IMAGE_FORMAT_MISMATCH",
   "CL_IMAGE_FORMAT_NOT_SUPPORTED",
   "CL_BUILD_PROGRAM_FAILURE",
   "CL_MAP_FAILURE",
   "CL_MISALIGNED_SUB_BUFFER_OFFSET",
   "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST",
   "CL_UNKNOWN",
   "CL_UNKNOWN",
   "CL_UNKNOWN",
   "CL_UNKNOWN",
   "CL_UNKNOWN",
   "CL_UNKNOWN",
   "CL_UNKNOWN",
   "CL_UNKNOWN",
   "CL_UNKNOWN",
   "CL_UNKNOWN",
   "CL_UNKNOWN",
   "CL_UNKNOWN",
   "CL_UNKNOWN",
   "CL_UNKNOWN",
   "CL_UNKNOWN",
   "CL_INVALID_VALUE",
   "CL_INVALID_DEVICE_TYPE",
   "CL_INVALID_PLATFORM",
   "CL_INVALID_DEVICE",
   "CL_INVALID_CONTEXT",
   "CL_INVALID_QUEUE_PROPERTIES",
   "CL_INVALID_COMMAND_QUEUE",
   "CL_INVALID_HOST_PTR",
   "CL_INVALID_MEM_OBJECT",
   "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR",
   "CL_INVALID_IMAGE_SIZE",
   "CL_INVALID_SAMPLER",
   "CL_INVALID_BINARY",
   "CL_INVALID_BUILD_OPTIONS",
   "CL_INVALID_PROGRAM",
   "CL_INVALID_PROGRAM_EXECUTABLE",
   "CL_INVALID_KERNEL_NAME",
   "CL_INVALID_KERNEL_DEFINITION",
   "CL_INVALID_KERNEL",
   "CL_INVALID_ARG_INDEX",
   "CL_INVALID_ARG_VALUE",
   "CL_INVALID_ARG_SIZE",
   "CL_INVALID_KERNEL_ARGS",
   "CL_INVALID_WORK_DIMENSION",
   "CL_INVALID_WORK_GROUP_SIZE",
   "CL_INVALID_WORK_ITEM_SIZE",
   "CL_INVALID_GLOBAL_OFFSET",
   "CL_INVALID_EVENT_WAIT_LIST",
   "CL_INVALID_EVENT",
   "CL_INVALID_OPERATION",
   "CL_INVALID_GL_OBJECT",
   "CL_INVALID_BUFFER_SIZE",
   "CL_INVALID_MIP_LEVEL",
   "CL_INVALID_GLOBAL_WORK_SIZE",
   "CL_INVALID_PROPERTY"
};



OpenCLError::OpenCLError(const char* file, int line, cl::Error& e):
   Exception(file,line,"OpenCLError"),
   _clFunc(e.what()),
   _code(e.err())
{}



const char* OpenCLError::clFunc()
{
   return _clFunc;
}



cl_int OpenCLError::code()
{
   // Return OpenCL error code.
   return _code;
}



const char* OpenCLError::code_str()
{
   // Make sure OpenCL code is within range.
   if (_code<=0&&_code>=-64)
   {
      // Return pointer to string that describes error code.
      return c_clDescErrors[-_code];
   }
   else
   {
      // Unknown error code, return NULL.
      return NULL;
   }
}



InvalidInput::InvalidInput(const char* file, int line):
   Exception(file,line,"InvalidInput")
{}



InvalidUse::InvalidUse(const char* file, int line):
   Exception(file,line,"InvalidUse")
{}