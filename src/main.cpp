#include "linker.cpp"
#include <iostream>



int main(int argc, char* argv[])
{
   KINCFactory factory;
   try
   {
      Ace::DataMap dataMap(factory);
      Ace::LinuxTerm::stty_raw();
      Ace::LinuxTerm terminal;
      Ace::Console console(argc,argv,terminal,factory,dataMap,"TEST");
      console.command("gpu set 0:0");
      console.command("spearman");
   }
   catch(Ace::Exception e)
   {
      std::cout << "Fatal Exception Caught!\n";
      std::cout << e.what() << "(" << e.function() << ":" << e.line() <<  ")";
      if (e.detail())
      {
         std::cout << ": " << e.detail();
      }
      std::cout << "\n";
      Ace::LinuxTerm::stty_cooked();
      return -1;
   }
   catch (std::exception& e)
   {
      std::cout << "Fatal STD Exception Caught!\n";
      std::cout << e.what() << "\n";
      Ace::LinuxTerm::stty_cooked();
      return -1;
   }
   catch (...)
   {
      Ace::LinuxTerm::stty_cooked();
      throw;
   }
   Ace::LinuxTerm::stty_cooked();
   return 0;
}
