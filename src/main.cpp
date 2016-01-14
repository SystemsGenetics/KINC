#include "console.h"
#include "linuxterm.h"
#include <iostream>



int main(int argc, char* argv[])
{
   try
   {
      LinuxTerm terminal;
      LinuxTerm::stty_raw();
      Console console(argc,argv,terminal);
      console.run();
      LinuxTerm::stty_cooked();
   }
   catch (...)
   {
      LinuxTerm::stty_cooked();
      throw;
   }
   return 0;
}
