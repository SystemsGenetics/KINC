/*
 * See GPL.txt for copyright information.
 *
 * Author: Joshua Burns
 *
 */
#ifndef EXCEPTION_H
#define EXCEPTION_H
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>



class DataPlugin;
class Analytic;



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
                                     cl::Error& e)
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
   const char* _file;
   const int _line;
   const char* _what;
};



/// @brief Base exception class for Data objects.
///
/// Base exception class that any data implementation throws if an error has
/// occured within its object.
class DataException : public Exception
{
public:
   enum class Level
   {
      general,
      warning,
      caution,
      severe,
      fatal
   };
   // *
   // * BASIC METHODS
   // *
   DataException(const char*,int,DataPlugin*,const char*,Level);
   // *
   // * FUNCTIONS
   // *
   DataPlugin* who();
   Level level();
private:
   // *
   // * VARIABLES
   // *
   DataPlugin* _who;
   Level _level;
};



/// @brief Base exception class for analytic objects.
///
/// Base exception class that any analytic implementation throws if an error has
/// occured within its object.
class AnalyticException : public Exception
{
public:
   // *
   // * BASIC METHODS
   // *
   AnalyticException(const char*,int,Analytic*,const char*);
   // *
   // * FUNCTIONS
   // *
   Analytic* who();
private:
   // *
   // * VARIABLES
   // *
   Analytic* _who;
};



/// An error occured after calling a system function.
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
   const char* _system;
};



/// An error occured after calling an OpenCL function.
class OpenCLError : public Exception
{
public:
   // *
   // * BASIC METHODS
   // *
   OpenCLError(const char*,int,cl::Error&);
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
   static const char* c_clDescErrors[];
   const char* _clFunc;
   cl_int _code;
};



#endif
