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



DataException::DataException(const char* file, int line, Data* who,
                             const char* what):
   Exception(file,line,what),
   _who(who)
{}



Data* DataException::who()
{
   return _who;
}



AnalyticException::AnalyticException(const char* file, int line, Analytic* who,
                                     const char* what):
   Exception(file,line,what),
   _who(who)
{}



Analytic* AnalyticException::who()
{
   return _who;
}
