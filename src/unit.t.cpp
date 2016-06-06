#include <iostream>
#include "unit.h"
#include "console.h"
#include "linuxterm.h"
#include "plugins/plugins.h"



const char* unit::headerStr {nullptr};
int unit::numTestsDone {0};



AnalyticPlugin* KINCPlugins::new_analytic(const std::string& type)
{
   AnalyticPlugin* ret = nullptr;
   return ret;
}



int main(int argc, char* argv[])
{
   unit::initiate();
   if (unit::getopts::main()&&
       unit::filemem::main()&&
       unit::fstring::main()&&
       unit::histitem::main()&&
       unit::history::main()&&
       unit::kincfile::main()&&
       unit::datamap::main())
   {
      unit::complete();
   }
   else
   {
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
