#include <unistd.h>
#include <termios.h>
#include <iostream>
#include <sys/ioctl.h>
#include <cstring>
#include <cstdio>
#include "linuxterm.h"
#include "exception.h"



bool LinuxTerm::_lock {false};
bool LinuxTerm::_cooked {true};



/// @brief Sets terminal attributes to raw.
///
/// Grabs the terminal attributes and sets them to not echo and not be
/// canonical. Changes terminal status to raw mode.
///
/// @pre Terminal status is in cooked mode.
void LinuxTerm::stty_raw()
{
   assert<InvalidUse>(_cooked,__FILE__,__LINE__);
   struct termios term = {0};
   bool cond;
   cond = tcgetattr(0,&term)>=0;
   assert<SystemError>(cond,__FILE__,__LINE__,"tcgetattr");
   term.c_lflag &= ~ICANON;
   term.c_lflag &= ~ECHO;
   term.c_cc[VMIN] = 1;
   term.c_cc[VTIME] = 0;
   cond = tcsetattr(0,TCSANOW,&term)>=0;
   assert<SystemError>(cond,__FILE__,__LINE__,"tcsetattr");
   _cooked = false;
}



/// @brief Sets terminal attributes to cooked.
///
/// Grabs the terminal attributes and sets them to echo and canoncial. Changes
/// terminal status to cooked mode.
///
/// @pre Terminal status is in raw mode.
/// @pre There can be no current instance of this function's class.
void LinuxTerm::stty_cooked()
{
   assert<InvalidUse>(!_cooked,__FILE__,__LINE__);
   assert<InvalidUse>(!_lock,__FILE__,__LINE__);
   struct termios term = {0};
   bool cond;
   cond = tcgetattr(0,&term)>=0;
   assert<SystemError>(cond,__FILE__,__LINE__,"tcgetattr");
   term.c_lflag |= ICANON;
   term.c_lflag |= ECHO;
   cond = tcsetattr(0,TCSANOW,&term)>=0;
   assert<SystemError>(cond,__FILE__,__LINE__,"tcsetattr");
   _cooked = true;
}



/// @brief Initializes terminal.
///
/// Grabs the width of the program's terminal in characters and locks terminal
/// for exclusive use.
///
/// @pre There can be no current instance of this constructor's class.
/// @pre The terminal must be in raw mode.
LinuxTerm::LinuxTerm():
   _i {_line.end()},
   _chCount {0}
{
   assert<InvalidUse>(!_lock,__FILE__,__LINE__);
   assert<InvalidUse>(!_cooked,__FILE__,__LINE__);
   _lock = true;
#ifdef TIOCGSIZE
   struct ttysize ts;
   ioctl(STDIN_FILENO,TIOCGSIZE,&ts);
   int _cols = ts.ts_cols;
#elif defined(TIOCGWINSZ)
   struct winsize ts;
   ioctl(STDIN_FILENO,TIOCGWINSZ,&ts);
   _cols = ts.ws_col;
#endif
}



/// Release lock on terminal.
LinuxTerm::~LinuxTerm()
{
   _lock = false;
}



/// Set new header for terminal input.
void LinuxTerm::header(const std::string& header)
{
   _header = header;
}



/// Set floating point precision.
///
/// @todo Need to finish this function.
void LinuxTerm::precision(int)
{}



// Output operators, self-explanatory.
LinuxTerm& LinuxTerm::operator<<(short n)
{
   std::cout << n;
   return *this;
}
LinuxTerm& LinuxTerm::operator<<(unsigned short n)
{
   std::cout << n;
   return *this;
}
LinuxTerm& LinuxTerm::operator<<(int n)
{
   std::cout << n;
   return *this;
}
LinuxTerm& LinuxTerm::operator<<(unsigned int n)
{
   std::cout << n;
   return *this;
}
LinuxTerm& LinuxTerm::operator<<(long n)
{
   std::cout << n;
   return *this;
}
LinuxTerm& LinuxTerm::operator<<(unsigned long n)
{
   std::cout << n;
   return *this;
}
LinuxTerm& LinuxTerm::operator<<(float n)
{
   std::cout << n;
   return *this;
}
LinuxTerm& LinuxTerm::operator<<(double n)
{
   std::cout << n;
   return *this;
}
LinuxTerm& LinuxTerm::operator<<(const char* n)
{
   std::cout << n;
   return *this;
}
LinuxTerm& LinuxTerm::operator<<(const std::string& n)
{
   std::cout << n;
   return *this;
}



/// @brief Read single line of input from user.
///
/// Reads a single line from the user on the terminal, blocking execution of the
/// program until the line has finished being read from the user.
///
/// @param buffer String that the line of new input will be written to.
///
/// @pre The terminal must be in raw mode and the cursor must be on a new and
/// blank line.
void LinuxTerm::operator>>(std::string& buffer)
{
   assert<InvalidUse>(!_cooked,__FILE__,__LINE__);
   buffer.clear();
   _line.clear();
   _i = _line.end();
   reprint(false);
   char input;
   while ((input = getchar())!='\n')
   {
      if (input>=' '&&input<='~')
      {
         _i = _line.insert(_i,input);
         _i++;
      }
      else
      {
         switch (input)
         {
         case _backspaceCh:
            if (_i!=_line.begin())
            {
               _i = _line.erase(--_i);
            }
            break;
         case _escapeCh:
            if (getchar()=='[')
            {
               switch(getchar())
               {
               case '3':
                  if (getchar()=='~'&&_i!=_line.end())
                  {
                     _i = _line.erase(_i);
                  }
                  break;
               case _arrowRightCh:
                  if (_i!=_line.end())
                  {
                     _i++;
                  }
                  break;
               case _arrowLeftCh:
                  if (_i!=_line.begin())
                  {
                     _i--;
                  }
                  break;
               }
            }
            break;
         }
      }
      reprint(true);
   }
   for (auto i:_line)
   {
      buffer += i;
   }
}



/// @brief Setting special control types.
///
/// Internal callback function for handling special control types enumerated by
/// Terminal interface class.
///
/// @todo Still need to finish this function.
void LinuxTerm::set_ops(Ops op)
{
   switch (op)
   {
   case Ops::newline:
      std::cout << "\n";
      break;
   case Ops::flush:
      std::cout << std::flush;
      break;
   }
}



/// @brief Move cursor up.
///
/// Move cursor up as many lines needed to get back to origional line from
/// number of characters that has been typed.
///
/// @param chCount Number of characters that has been typed.
void LinuxTerm::reset_cursor(int chCount)
{
   for (int i=0;i<(chCount/_cols);i++)
   {
      std::cout << _cursorUpStr;
   }
   std::cout << "\r";
}



/// @brief Prints current user input line.
///
/// If specified, moves cursor back to beginning of current input line and
/// header. Either way prints out header and input line along with one extra
/// white space. If the cursor resets itself, this effectively overwrites the
/// previous input line printed.
///
/// @param rewind Specifies if cursor will be reset to beginning of line.
void LinuxTerm::reprint(bool rewind)
{
   if (rewind)
   {
      reset_cursor(_chCount-1);
   }
   std::cout << _header;
   for (auto i:_line)
   {
      std::cout << i;
   }
   std::cout << " ";
   reset_cursor(_line.size()+_header.size());
   _chCount = _header.size();
   std::cout << _boldTextStr;
   std::cout << _header;
   std::cout << _normTextStr;
   for (auto i=_line.begin();i!=_i;i++)
   {
      std::cout << *i;
      _chCount++;
   }
   std::cout << std::flush;
}
