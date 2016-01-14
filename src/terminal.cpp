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



void Terminal::print(short n)
{
   *this << n;
}



void Terminal::print(unsigned short n)
{
   *this << n;
}



void Terminal::print(int n)
{
   *this << n;
}



void Terminal::print(unsigned int n)
{
   *this << n;
}



void Terminal::print(long n)
{
   *this << n;
}



void Terminal::print(unsigned long n)
{
   *this << n;
}



void Terminal::print(float n)
{
   *this << n;
}



void Terminal::print(double n)
{
   *this << n;
}



void Terminal::print(const char* n)
{
   *this << n;
}



void Terminal::print(const std::string& n)
{
   *this << n;
}



Terminal& Terminal::operator<<(Terminal& (*pf)(Terminal&))
{
   return pf(*this);
}
