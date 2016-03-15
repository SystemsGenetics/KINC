#include "terminal.h"



/// Sets a newline operator to Terminal output.
Terminal& Terminal::newline(Terminal& term)
{
   term.set_ops(Ops::newline);
   return term;
}



/// Sets a newline operator then a flush operator to Terminal output.
Terminal& Terminal::endl(Terminal& term)
{
   term.set_ops(Ops::newline);
   term.set_ops(Ops::flush);
   return term;
}



/// Sets a flush operator to Terminal output.
Terminal& Terminal::flush(Terminal& term)
{
   term.set_ops(Ops::flush);
   return term;
}



/// Sets a general message state operator to Terminal output.
Terminal& Terminal::general(Terminal& term)
{
   term.set_ops(Ops::general);
   return term;
}



/// Sets a warning message state operator to Terminal output.
Terminal& Terminal::warning(Terminal& term)
{
   term.set_ops(Ops::warning);
   return term;
}



/// Sets a error message state operator to Terminal output.
Terminal& Terminal::error(Terminal& term)
{
   term.set_ops(Ops::error);
   return term;
}



// Self-referencing print functions that redirect to output operators.
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



/// Takes pointers to functions that can manipulate Terminal operators.
Terminal& Terminal::operator<<(Terminal& (*pf)(Terminal&))
{
   return pf(*this);
}
