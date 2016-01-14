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



class Exception
/*
 * Base class for any exception that this program will throw. It provides
 * information about what file and what line threw an exception along with
 * a textual field describing what type of exception was thrown.
 */
{
public:
   // ****************************** Functions ******************************
   Exception(const char*,int,const char*);
   const char* file();
   const int line();
   const char* what();
private:
   // ****************************** Variables ******************************
   const char* _file;
   const int _line;
   const char* _what;
};



class DataException : public Exception
/*
 * Base exception class that any data implementation throws if an error has
 * occured within its object.
 */
{
public:
   // ****************************** Functions ******************************
   static void assert(bool,const char*,int,Data*,const char*);
   DataException(const char*,int,Data*,const char*);
   Data* who();
private:
   // ****************************** Variables ******************************
   Data* _who;
};



class AnalyticException : public Exception
/*
 * Base exception class that any analytic implementation throws if an error has
 * occured within its object.
 */
{
public:
   // ****************************** Functions ******************************
   static void assert(bool,const char*,int,Analytic*,const char*);
   AnalyticException(const char*,int,Analytic*,const char*);
   Analytic* who();
private:
   // ****************************** Variables ******************************
   Analytic* _who;
};



class SystemError : public Exception
/*
 * An error occured after calling a system function.
 */
{
public:
   // ****************************** Functions ******************************
   static void assert(bool,const char*,int,const char*);
   SystemError(const char*,int,const char*);
   const char* system();
private:
   // ****************************** Variables ******************************
   const char* _system;
};



class OpenCLError : public Exception
/*
 * An error occured after calling an OpenCL function.
 */
{
public:
   // ****************************** Functions ******************************
   static void assert(bool,const char*,int,cl::Error&);
   OpenCLError(const char*,int,cl::Error&);
   const char* clFunc();
   cl_int code();
   const char* code_str();
private:
   // ****************************** Variables ******************************
   static const char* c_clDescErrors[];
   const char* _clFunc;
   cl_int _code;
};



struct InvalidInput : public Exception
/*
 * A function was given an argument that is invalid.
 */
{
   static void assert(bool,const char*,int);
   InvalidInput(const char*,int);
};



struct InvalidUse : public Exception
/*
 * An object was used in a way that is invalid.
 */
{
   static void assert(bool,const char*,int);
   InvalidUse(const char*,int);
};



struct OutOfRange : public Exception
/*
 * A function was given a numeric argument that is out of range.
 */
{
   static void assert(bool,const char*,int);
   OutOfRange(const char*,int);
};



#endif
