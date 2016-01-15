#include "terminal.h"



Terminal& Terminal::newline(Terminal& term)
{
   term.set_ops(Ops::newline);
}



Terminal& Terminal::endl(Terminal& term)
{
   term.set_ops(Ops::newline);
   term.set_ops(Ops::flush);
}



Terminal& Terminal::flush(Terminal& term)
{
   term.set_ops(Ops::flush);
}



Terminal& Terminal::general(Terminal& term)
{
   term.set_ops(Ops::general);
}



Terminal& Terminal::warning(Terminal& term)
{
   term.set_ops(Ops::warning);
}



Terminal& Terminal::error(Terminal& term)
{
   term.set_ops(Ops::error);
}



Terminal& Terminal::print(short n)
{
   return *this << n;
}



Terminal& Terminal::print(unsigned short n)
{
   return *this << n;
}



Terminal& Terminal::print(int n)
{
   return *this << n;
}



Terminal& Terminal::print(unsigned int n)
{
   return *this << n;
}



Terminal& Terminal::print(long n)
{
   return *this << n;
}



Terminal& Terminal::print(unsigned long n)
{
   return *this << n;
}



Terminal& Terminal::print(float n)
{
   return *this << n;
}



Terminal& Terminal::print(double n)
{
   return *this << n;
}



Terminal& Terminal::print(const char* n)
{
   return *this << n;
}



Terminal& Terminal::print(const std::string& n)
{
   return *this << n;
}



Terminal& Terminal::operator<<(Terminal& (*pf)(Terminal&))
{
   return pf(*this);
}
