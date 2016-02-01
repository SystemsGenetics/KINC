#include <iostream>
#include "unit.h"



int main(int argc, char* argv[])
{
   if (unit::linuxfile::main())
   {
      std::cout << "ALL PASSED." << std::endl;
   }
   else
   {
      std::cout << "FAILURE." << std::endl;
   }
   return 0;
}
