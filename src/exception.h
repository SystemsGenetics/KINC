#ifndef EXCEPTION_H
#define EXCEPTION_H
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>



class Data;
class Analytic;



class Exception
{
private:
   const char* _file;
   const int _line;
   const char* _what;
public:
   Exception(const char*,int,const char*);
   const char* file();
   const int line();
   const char* what();
};



class DataException : public Exception
{
private:
   Data* _who;
public:
   DataException(const char*,int,Data*,const char*);
   Data* who();
};



class AnalyticException : public Exception
{
private:
   Analytic* _who;
public:
   AnalyticException(const char*,int,Analytic*,const char*);
   Analytic* who();
};



class SystemError : public Exception
{
private:
   const char* _system;
public:
   SystemError(const char*,int,const char*);
   const char* system();
};



class OpenCLError : public Exception
{
private:
   static const char* c_clDescErrors[];
   const char* _clFunc;
   cl_int _code;
public:
   OpenCLError(const char*,int,cl::Error&);
   const char* clFunc();
   cl_int code();
   const char* code_str();
};



struct InvalidInput : public Exception
{
   InvalidInput(const char*,int);
};



struct InvalidUse : public Exception
{
   InvalidUse(const char*,int);
};



#endif
