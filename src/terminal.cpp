#include "terminal.h"
#include <unistd.h>
#include <termios.h>
#include "exception.h"
#include <iostream>
#include <sys/ioctl.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>

#define BACKSPACE 127
#define SPEC_KEY 27
#define SPEC_KEY2 91
#define SPEC_KEY3 51
#define DELETE 126
#define ARROW_RIGHT 67
#define ARROW_LEFT 68



void Terminal::reset_cursor(int chCount)
{
   for (int i=0;i<(chCount/_cols);i++)
   {
      std::cout << (char)27 << "[A";
   }
   std::cout << "\r";
}



void Terminal::reprint(bool rewind)
{
   if (rewind)
   {
      reset_cursor(_chCount-1);
   }
   std::cout << _header;
   for (std::list<char>::iterator i=_line.begin();i!=_line.end();i++)
   {
      std::cout << *i;
   }
   std::cout << " ";
   reset_cursor(_line.size()+_header.size());
   _chCount = _header.size();
   std::cout << _header;
   for (std::list<char>::iterator i=_line.begin();i!=_i;i++)
   {
      std::cout << *i;
      _chCount++;
   }
   std::cout << std::flush;
}



Terminal::Terminal():
   _i(_line.end()),
   _chCount(0)
{
#ifdef TIOCGSIZE
   struct ttysize ts;
   ioctl(STDIN_FILENO, TIOCGSIZE, &ts);
   int _cols = ts.ts_cols;
#elif defined(TIOCGWINSZ)
   struct winsize ts;
   ioctl(STDIN_FILENO, TIOCGWINSZ, &ts);
   _cols = ts.ws_col;
#endif
}



void Terminal::reset_line()
{
   reprint(false);
}



std::string Terminal::get()
{
   std::string ret;
   for (std::list<char>::iterator i = _line.begin();i!=_line.end();i++)
   {
      ret += *i;
   }
   return ret;
}



void Terminal::clear(const char* header)
{
   _header = header;
   _line.clear();
   reprint(false);
}



Terminal& Terminal::operator<<(char input)
{
   if (input>=32&&input<=126)
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



void Terminal::stty_raw()
{
   struct termios term = {0};
   if (tcgetattr(0,&term)<0)
   {
      throw SystemError(__FILE__,__LINE__,"tcgetattr");
   }
   term.c_lflag &= ~ICANON;
   term.c_lflag &= ~ECHO;
   term.c_cc[VMIN] = 1;
   term.c_cc[VTIME] = 0;
   if (tcsetattr(0,TCSANOW,&term)<0)
   {
      throw SystemError(__FILE__,__LINE__,"tcsetattr");
   }
}



void Terminal::stty_cooked()
{
   struct termios term = {0};
   if (tcgetattr(0,&term)<0)
   {
      throw SystemError(__FILE__,__LINE__,"tcgetattr");
   }
   term.c_lflag |= ICANON;
   term.c_lflag |= ECHO;
   if (tcsetattr(0,TCSANOW,&term)<0)
   {
      throw SystemError(__FILE__,__LINE__,"tcsetattr");
   }
}
