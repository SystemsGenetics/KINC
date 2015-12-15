#include "console.h"
#include <iostream>
#include <unistd.h>
#include <termios.h>
#include "exception.h"
#include <stdio.h>

Console g_console;



void Console::stty_raw(struct termios* term)
{
   if (tcgetattr(0,term)<0)
   {
      throw SystemError(__FILE__,__LINE__,"tcgetattr");
   }
   term->c_lflag &= ~ICANON;
   term->c_lflag &= ~ECHO;
   term->c_cc[VMIN] = 1;
   term->c_cc[VTIME] = 0;
   if (tcsetattr(0,TCSANOW,term)<0)
   {
      throw SystemError(__FILE__,__LINE__,"tcsetattr");
   }
}



void Console::stty_cooked(struct termios* term)
{
   term->c_lflag |= ICANON;
   term->c_lflag |= ECHO;
   if (tcsetattr(0,TCSANOW,term)<0)
   {
      throw SystemError(__FILE__,__LINE__,"tcsetattr");
   }
}



void Console::command()
{
   struct termios term = {0};
   stty_raw(&term);
   char input;
   int i = 0;
   std::string line;
   bool alive = true;
   while (alive)
   {
      input = getchar();
      if (input>=32&&input<=126) // regular input
      {
         if (i==line.size())
         {
            i++;
            line += input;
         }
         else
         {
            line[i++] = input;
         }
         std::cout << input;
      }
      else
      {
         switch(input)
         {
         case '\n': // newline
            alive = false;
            break;
         case 127: // backspace key
            if (i>0)
            {
               i--;
               line.erase(i,i+1);
               std::cout << "\b \b";
            }
            break;
         case 27: // special key hit
            if (getchar()==91) // Arrow key
            {
               switch(getchar())
               {
               case 67: //right
                  if (i<line.size())
                  {
                     std::cout << line[i++];
                  }
                  break;
               case 68: //left
                  if (i>0)
                  {
                     i--;
                     std::cout << "\b";
                  }
                  break;
               }
            }
            break;
         }
      }
   }
   std::cout << "\n" << line << std::endl;
   stty_cooked(&term);
}



Console::Console():
   out(ConsoleStream::out),
   warn(ConsoleStream::warning),
   err(ConsoleStream::error)
{}



void Console::run(int argc, char* argv[])
{
   command();
}
