#include "console.h"
#include "linuxterm.h"
#include "exception.h"
#include <iostream>



int main(int argc, char* argv[])
{
   try
   {
      LinuxTerm::stty_raw();
      LinuxTerm terminal;
      Console console(argc,argv,terminal);
      console.run();
   }
   catch (...)
   {
      LinuxTerm::stty_cooked();
      throw;
   }
   LinuxTerm::stty_cooked();
   return 0;
}
