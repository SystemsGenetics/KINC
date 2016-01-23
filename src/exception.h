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



class Data;
class Analytic;



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
   // *
   // * BASIC METHODS
   // *
   DataException(const char*,int,Data*,const char*);
   // *
   // * STATIC FUNCTIONS
   // *
   static void assert(bool,const char*,int,Data*,const char*);
   // *
   // * FUNCTIONS
   // *
   Data* who();
private:
   // *
   // * VARIABLES
   // *
   Data* _who;
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
   // * STATIC FUNCTIONS
   // *
   static void assert(bool,const char*,int,Analytic*,const char*);
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
   // * STATIC FUNCTIONS
   // *
   static void assert(bool,const char*,int,const char*);
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
   // * STATIC FUNCTIONS
   // *
   static void assert(bool,const char*,int,cl::Error&);
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



/// A function was given an argument that is invalid.
struct InvalidInput : public Exception
{
   // *
   // * BASIC METHODS
   // *
   InvalidInput(const char*,int);
   // *
   // * STATIC FUNCTIONS
   // *
   static void assert(bool,const char*,int);
};



/// An object was used in a way that is invalid.
struct InvalidUse : public Exception
{
   // *
   // * BASIC METHODS
   // *
   InvalidUse(const char*,int);
   // *
   // * STATIC FUNCTIONS
   // *
   static void assert(bool,const char*,int);
};



/// A function was given a numeric argument that is out of range.
struct OutOfRange : public Exception
{
   // *
   // * BASIC METHODS
   // *
   OutOfRange(const char*,int);
   // *
   // * STATIC FUNCTIONS
   // *
   static void assert(bool,const char*,int);
};



#endif
