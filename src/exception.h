/*
 * See GPL.txt for copyright information.
 *
 * Author: Joshua Burns
 *
 */
#ifndef EXCEPTION_H
#define EXCEPTION_H
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.h>
#define DATA_EXCEPTION(N,X,L) struct X : public DataException\
                              {\
                                 X(const char* file, int line):\
                                    DataException(file,line,#N"::"#X,\
                                                  Level::L)\
                                 {}\
                              };
#define ANALYTIC_EXCEPTION(N,X) struct X : public AnalyticException\
                                {\
                                   X(const char* file, int line):\
                                      AnalyticException(file,line,#N"::"#X)\
                                   {}\
                                };
#define OPENCL_EXCEPTION(X,F) struct X : public OpenCLError\
                              {\
                                 X(const char* file, int line, cl_int e):\
                                    OpenCLError(file,line,#F,e)\
                                 {}\
                              };
#define ACE_EXCEPTION(N,X) struct X : public Exception\
                           {\
                              X(const char* file, int line):\
                                 Exception(file,line,#N"::"#X)\
                              {}\
                           };



/// Tests assertion and throws basic Exception if it fails.
///
/// @param cond Condition to be tested.
/// @param file File where assertion is being tested.
/// @param line Line in file where assertion is being tested.
template<class T> inline void assert(bool cond, const char* file, int line)
{
   if (!cond)
   {
      throw T(file,line);
   }
}
/// Tests assertion and throws SystemError Exception if it fails.
///
/// @param cond Condition to be tested.
/// @param file File where assertion is being tested.
/// @param line Line in file where assertion is being tested.
/// @param system System function that failed.
template<class T> inline void assert(bool cond, const char* file, int line,
                                     const char* system)
{
   if (!cond)
   {
      throw T(file,line,system);
   }
}
/// Tests assertion and throws OpenCL Exception if it fails.
///
/// @param cond Condition to be tested.
/// @param file File where assertion is being tested.
/// @param line Line in file where assertion is being tested.
/// @param e Reference to OpenCL Error.
template<class T> inline void assert(bool cond, const char* file, int line,
                                     cl_int e)
{
   if (!cond)
   {
      throw T(file,line,e);
   }
}



/// @brief Base exception class.
///
/// Base class for any exception that this program will throw. It provides
/// information about what file and what line threw an exception along with
/// a textual field describing what type of exception was thrown.
///
/// @author Josh Burns
/// @date 21 March 2016
class Exception
{
public:
   // *
   // * BASIC METHODS
   // *
   Exception(const char*,int,const char*);
   // *
   // * FUNCTIONS
   // *
   const char* file();
   const int line();
   const char* what();
private:
   // *
   // * VARIABLES
   // *
   /// Source file where exception was thrown.
   const char* _file;
   /// Line in source file where exception was thrown.
   const int _line;
   /// Description of type of exception thrown.
   const char* _what;
};



/// @ingroup dataplugin
/// @brief Base exception class for Data objects.
///
/// Base exception class that any data implementation throws if an error has
/// occured within its object.
///
/// @author Josh Burns
/// @date 21 March 2016
class DataException : public Exception
{
public:
   /// Classifies severity of exception thrown.
   enum class Level
   {
      general, ///< General information, lowest severity.
      warning, ///< Warning with mild severity.
      caution, ///< Moderate severity issue.
      severe, ///< High severity issue.
      fatal ///< Fatal error occured, highest severity.
   };
   // *
   // * BASIC METHODS
   // *
   DataException(const char*,int,const char*,Level);
   // *
   // * FUNCTIONS
   // *
   Level level();
private:
   // *
   // * VARIABLES
   // *
   /// Severity level of exception.
   Level _level;
};



/// @brief Base exception class for analytic objects.
///
/// Base exception class that any analytic implementation throws if an error has
/// occured within its object.
///
/// @author Josh Burns
/// @date 21 March 2016
class AnalyticException : public Exception
{
public:
   // *
   // * BASIC METHODS
   // *
   AnalyticException(const char*,int,const char*);
};



/// An error occured after calling a system function.
///
/// @author Josh Burns
/// @date 21 March 2016
class SystemError : public Exception
{
public:
   // *
   // * BASIC METHODS
   // *
   SystemError(const char*,int,const char*);
   // *
   // * FUNCTIONS
   // *
   const char* system();
private:
   // *
   // * VARIABLES
   // *
   /// Name of system call that failed.
   const char* _system;
};



/// An error occured after calling an OpenCL function.
///
/// @author Josh Burns
/// @date 21 March 2016
class OpenCLError : public Exception
{
public:
   // *
   // * BASIC METHODS
   // *
   OpenCLError(const char*,int,const char*,cl_int);
   // *
   // * FUNCTIONS
   // *
   const char* clFunc();
   cl_int code();
   const char* code_str();
private:
   // *
   // * VARIABLES
   // *
   /// Static C-Style list of strings describing every possible OpenCL error.
   static const char* c_clDescErrors[];
   /// Name of OpenCL function that failed.
   const char* _clFunc;
   /// Internal OpenCL error code.
   cl_int _code;
};



#endif
