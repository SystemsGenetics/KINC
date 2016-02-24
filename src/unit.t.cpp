#include <iostream>
#include "unit.h"
#include "console.h"
#include "linuxterm.h"



const char* unit::headerStr {nullptr};



int main(int argc, char* argv[])
{
   if (unit::getopts::main()&&
       unit::filemem::main()&&
       unit::histitem::main())
   {
      std::cout << "ALL PASSED." << std::endl;
   }
   else
   {
      std::cout << "FAILURE." << std::endl;
      exit(1);
   }
   try
   {
      DataMap dataMap;
      LinuxTerm::stty_raw();
      LinuxTerm terminal;
      Console console(argc,argv,terminal,dataMap);
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
