#include <unistd.h>
#include <termios.h>
#include <iostream>
#include <sys/ioctl.h>
#include <cstring>
#include <cstdio>
#include "linuxterm.h"
#include "exception.h"

//TODO: Change these defines into cosntexpr char variables.
#define LOWER_LIMIT 32
#define UPPER_LIMIT 126
#define BACKSPACE 127
#define SPEC_KEY 27
#define SPEC_KEY2 91
#define SPEC_KEY3 51
#define DELETE 126
#define ARROW_RIGHT 67
#define ARROW_LEFT 68



// Make sure only one instance of this class exists.
bool LinuxTerm::_lock = false;
// Tracks if terminal is currently in cooked or raw mode.
bool LinuxTerm::_cooked = true;



void LinuxTerm::stty_raw()
/*
 * Uses Linux system calls to turn off echoing and canonical for terminal this
 * program is using.
 *
 * PRECONDITION:
 * 1. System currently in cooked mode.
 */
{
   InvalidUse::assert(_cooked,__FILE__,__LINE__);
   struct termios term = {0};
   bool cond;
   cond = tcgetattr(0,&term)>=0;
   SystemError::assert(cond,__FILE__,__LINE__,"tcgetattr");
   term.c_lflag &= ~ICANON;
   term.c_lflag &= ~ECHO;
   term.c_cc[VMIN] = 1;
   term.c_cc[VTIME] = 0;
   cond = tcsetattr(0,TCSANOW,&term)>=0;
   SystemError::assert(cond,__FILE__,__LINE__,"tcsetattr");
   _cooked = false;
}



void LinuxTerm::stty_cooked()
/*
 * Uses Linux system calls to turn on echoing and canonical for terminal this
 * program is using.
 *
 * PRECONDITION:
 * 1. System currently in raw mode.
 */
{
   InvalidUse::assert(!_cooked,__FILE__,__LINE__);
   struct termios term = {0};
   bool cond;
   cond = tcgetattr(0,&term)>=0;
   SystemError::assert(cond,__FILE__,__LINE__,"tcgetattr");
   term.c_lflag |= ICANON;
   term.c_lflag |= ECHO;
   cond = tcsetattr(0,TCSANOW,&term)>=0;
   SystemError::assert(cond,__FILE__,__LINE__,"tcsetattr");
   _cooked = true;
}



LinuxTerm::LinuxTerm()
/*
 * Grabs the width of the program's terminal in characters.
 *
 * PRECONDITIONS:
 * 1. _lock variable is false.
 */
   :_i(_line.end()),
   _chCount(0)
{
   InvalidUse::assert(!_lock,__FILE__,__LINE__);
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



LinuxTerm::~LinuxTerm()
{
   _lock = false;
}



void LinuxTerm::header(const std::string& header)
{
   _header = header;
}



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



void LinuxTerm::precision(int)
{}



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



void LinuxTerm::operator>>(std::string& buffer)
/*
 * Uses object's list of characters to build input line from user.
 *
 * buffer: string that user input will be written to, overwrites anything
 *         currently in the string.
 *
 * PRECONDITIONS:
 * 1. _cooked variable is false.
 */
{
   InvalidUse::assert(!_cooked,__FILE__,__LINE__);
   buffer.clear();
   _line.clear();
   _i = _line.end();
   reprint(false);
   char input;
   while ((input = getchar())!='\n')
   {
      if (input>=LOWER_LIMIT&&input<=UPPER_LIMIT)
      {
         _i = _line.insert(_i,input);
         _i++;
      }
      else
      {
         switch (input)
         {
         case BACKSPACE:
            if (_i!=_line.begin())
            {
               _i = _line.erase(--_i);
            }
            break;
         case SPEC_KEY:
            if (getchar()==SPEC_KEY2)
            {
               switch(getchar())
               {
               case SPEC_KEY3:
                  if (getchar()==DELETE&&_i!=_line.end())
                  {
                     _i = _line.erase(_i);
                  }
                  break;
               case ARROW_RIGHT:
                  if (_i!=_line.end())
                  {
                     _i++;
                  }
                  break;
               case ARROW_LEFT:
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



void LinuxTerm::reset_cursor(int chCount)
/*
 * Rewinds cursor of program's terminal to where the user's input line began
 * along with the header string, if any, set by user of class. The number of
 * lines to go back is computed by width of terminal.
 */
{
   for (int i=0;i<(chCount/_cols);i++)
   {
      // TODO; make constexpr of this and explain! (move cursor up)
      std::cout << (char)27 << "[A";
   }
   std::cout << "\r";
}



void LinuxTerm::reprint(bool rewind)
/*
 * Erases the current line the user is inputing and move cursor to beginning if
 * desired. (Re)Print the line the user is currently inputing.
 *
 * rewind: specifies if the line has already been printed once and should be
 *         erased and the cursor rewound to the beginning to printing out the
 *         user line again.
 */
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
   // TODO; make constexpr of these and explain them!
   // https://en.wikipedia.org/wiki/ANSI_escape_code
   // Make text bold
   std::cout << (char)27 << "[1m";
   std::cout << _header;
   // Reset to normal text
   std::cout << (char)27 << "[0m";
   for (auto i=_line.begin();i!=_i;i++)
   {
      std::cout << *i;
      _chCount++;
   }
   std::cout << std::flush;
}
