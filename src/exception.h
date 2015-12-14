#ifndef EXCEPTION_H
#define EXCEPTION_H



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



#endif
