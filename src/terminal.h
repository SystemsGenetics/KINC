#ifndef TERMINAL_H
#define TERMINAL_H
#include <string>



class Terminal
{
protected:
   enum class Ops {newline,flush,general,warning,error};
public:
   virtual ~Terminal() {}
   virtual void header(const std::string&) = 0;
   virtual void set_ops(Ops) = 0;
   static Terminal& newline(Terminal&);
   static Terminal& endl(Terminal&);
   static Terminal& flush(Terminal&);
   static Terminal& general(Terminal&);
   static Terminal& warning(Terminal&);
   static Terminal& error(Terminal&);
   virtual void precision(int) = 0;
   Terminal& print(short);
   Terminal& print(unsigned short);
   Terminal& print(int);
   Terminal& print(unsigned int);
   Terminal& print(long);
   Terminal& print(unsigned long);
   Terminal& print(float);
   Terminal& print(double);
   Terminal& print(const char*);
   Terminal& print(const std::string&);
   virtual Terminal& operator<<(short) = 0;
   virtual Terminal& operator<<(unsigned short) = 0;
   virtual Terminal& operator<<(int) = 0;
   virtual Terminal& operator<<(unsigned int) = 0;
   virtual Terminal& operator<<(long) = 0;
   virtual Terminal& operator<<(unsigned long) = 0;
   virtual Terminal& operator<<(float) = 0;
   virtual Terminal& operator<<(double) = 0;
   virtual Terminal& operator<<(const char*) = 0;
   virtual Terminal& operator<<(const std::string&) = 0;
   virtual void operator>>(std::string&) = 0;
   Terminal& operator<<(Terminal& (*pf)(Terminal&));
};



#endif
