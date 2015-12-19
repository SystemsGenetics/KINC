#ifndef CONSOLESTREAM_H
#define CONSOLESTREAM_H
#include <string>



class ConsoleStream
{
public:
   enum {out,warning,error};
private:
   int _which;
public:
   ConsoleStream(int);
   void print(short);
   void print(unsigned short);
   void print(int);
   void print(unsigned int);
   void print(long);
   void print(unsigned long);
   void print(float);
   void print(double);
   void print(const char*);
   void print(const std::string&);
   void flush();
   ConsoleStream& operator<<(short);
   ConsoleStream& operator<<(unsigned short);
   ConsoleStream& operator<<(int);
   ConsoleStream& operator<<(unsigned int);
   ConsoleStream& operator<<(long);
   ConsoleStream& operator<<(unsigned long);
   ConsoleStream& operator<<(float);
   ConsoleStream& operator<<(double);
   ConsoleStream& operator<<(const char*);
   ConsoleStream& operator<<(const std::string&);
};



#endif
