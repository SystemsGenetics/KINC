#include "consolestream.h"
#include <iostream>
#include "exception.h"



ConsoleStream::ConsoleStream(int which):
   _which(which)
{
   if (which!=out&&which!=warning&&which!=error)
   {
      throw InvalidInput(__FILE__,__LINE__);
   }
}



void ConsoleStream::print(short val)
{
   switch(_which)
   {
   case out:
      std::cout << val;
      break;
   case warning:
      std::cout << val;
      break;
   case error:
      std::cerr << val;
      break;
   }
}



void ConsoleStream::print(unsigned short val)
{
   switch(_which)
   {
   case out:
      std::cout << val;
      break;
   case warning:
      std::cout << val;
      break;
   case error:
      std::cerr << val;
      break;
   }
}



void ConsoleStream::print(int val)
{
   switch(_which)
   {
   case out:
      std::cout << val;
      break;
   case warning:
      std::cout << val;
      break;
   case error:
      std::cerr << val;
      break;
   }
}



void ConsoleStream::print(unsigned int val)
{
   switch(_which)
   {
   case out:
      std::cout << val;
      break;
   case warning:
      std::cout << val;
      break;
   case error:
      std::cerr << val;
      break;
   }
}



void ConsoleStream::print(long val)
{
   switch(_which)
   {
   case out:
      std::cout << val;
      break;
   case warning:
      std::cout << val;
      break;
   case error:
      std::cerr << val;
      break;
   }
}



void ConsoleStream::print(unsigned long val)
{
   switch(_which)
   {
   case out:
      std::cout << val;
      break;
   case warning:
      std::cout << val;
      break;
   case error:
      std::cerr << val;
      break;
   }
}



void ConsoleStream::print(float val)
{
   switch(_which)
   {
   case out:
      std::cout << val;
      break;
   case warning:
      std::cout << val;
      break;
   case error:
      std::cerr << val;
      break;
   }
}



void ConsoleStream::print(double val)
{
   switch(_which)
   {
   case out:
      std::cout << val;
      break;
   case warning:
      std::cout << val;
      break;
   case error:
      std::cerr << val;
      break;
   }
}



void ConsoleStream::print(const char* val)
{
   switch(_which)
   {
   case out:
      std::cout << val;
      break;
   case warning:
      std::cout << val;
      break;
   case error:
      std::cerr << val;
      break;
   }
}



void ConsoleStream::print(const std::string& val)
{
   switch(_which)
   {
   case out:
      std::cout << val;
      break;
   case warning:
      std::cout << val;
      break;
   case error:
      std::cerr << val;
      break;
   }
}



void ConsoleStream::flush()
{
   switch(_which)
   {
   case out:
      std::cout << std::flush;
      break;
   case warning:
      std::cout << std::flush;
      break;
   case error:
      std::cerr << std::flush;
      break;
   }
}



ConsoleStream& ConsoleStream::operator<<(short val)
{
   print(val);
   return *this;
}



ConsoleStream& ConsoleStream::operator<<(unsigned short val)
{
   print(val);
   return *this;
}



ConsoleStream& ConsoleStream::operator<<(int val)
{
   print(val);
   return *this;
}



ConsoleStream& ConsoleStream::operator<<(unsigned int val)
{
   print(val);
   return *this;
}



ConsoleStream& ConsoleStream::operator<<(long val)
{
   print(val);
   return *this;
}



ConsoleStream& ConsoleStream::operator<<(unsigned long val)
{
   print(val);
   return *this;
}



ConsoleStream& ConsoleStream::operator<<(float val)
{
   print(val);
   return *this;
}



ConsoleStream& ConsoleStream::operator<<(double val)
{
   print(val);
   return *this;
}



ConsoleStream& ConsoleStream::operator<<(const char* val)
{
   print(val);
   return *this;
}



ConsoleStream& ConsoleStream::operator<<(const std::string& val)
{
   print(val);
   return *this;
}
