#include "exception.h"



/// C-Style list of all possible OpenCL errors.
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



/// Initializes base exception object.
///
/// @param file Name of file where exception is thrown.
/// @param line Line number in file where exception is thrown.
/// @param what Description of exception thrown.
Exception::Exception(const char* file, int line, const char* what):
   _file(file),
   _line(line),
   _what(what)
{}



/// Gets file where exception was thrown.
///
/// @return Name of file.
const char* Exception::file()
{
   return _file;
}



/// Gets line number where exception was thrown.
///
/// @return File line number.
const int Exception::line()
{
   return _line;
}



/// Gets description from exception that was thrown.
///
/// @return Description of exception.
const char* Exception::what()
{
   return _what;
}



/// Initializes data plugin exception object.
///
/// @param file Name of file where exception is thrown.
/// @param line Line number in file where exception is thrown.
/// @param what Description of type of exception thrown.
/// @param lvl Severity level of exception.
DataException::DataException(const char* file, int line,
                             const char* what, Level lvl):
   Exception(file,line,what),
   _level(lvl)
{}



/// Get severity level of data exception that was thrown.
///
/// @return Severity of exception.
DataException::Level DataException::level()
{
   return _level;
}



/// Initializes analytic plugin exception object.
///
/// @param file Name of file where exception is thrown.
/// @param line Line number in file where exception is thrown.
/// @param what Description of type of exception thrown.
AnalyticException::AnalyticException(const char* file, int line,
                                     const char* what):
   Exception(file,line,what)
{}



/// Initializes system error exception object.
///
/// @param file Name of file where exception is thrown.
/// @param line Line number in file where exception is thrown.
/// @param system Name of system call that failed.
SystemError::SystemError(const char* file, int line, const char* system):
   Exception(file,line,"SystemError"),
   _system(system)
{}



/// Get name of system call that failed.
///
/// @return Name of system call.
const char* SystemError::system()
{
   return _system;
}



/// Initializes OpenCL error exception object.
///
/// @param file Name of file where exception is thrown.
/// @param line Line number in file where exception is thrown.
/// @param e OpenCL error object.
OpenCLError::OpenCLError(const char* file, int line, const char* func,
                         cl_int e):
   Exception(file,line,"OpenCLError"),
   _clFunc(func),
   _code(e)
{}



/// Get name of OpenCL function that failed.
///
/// @return Name of OpenCL function.
const char* OpenCLError::clFunc()
{
   return _clFunc;
}



/// Get OpenCL error code.
///
/// @return OpenCL error code.
cl_int OpenCLError::code()
{
   return _code;
}



/// Get textual description of OpenCL error code.
///
/// @return Description of OpenCL error code.
const char* OpenCLError::code_str()
{
   if (_code<=0&&_code>=-64)
   {
      return c_clDescErrors[-_code];
   }
   else
   {
      return c_clDescErrors[15];
   }
}
