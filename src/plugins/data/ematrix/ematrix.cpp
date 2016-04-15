#include "ematrix.h"
#include <fstream>
#include <sstream>



void ematrix::load(GetOpts &ops, Terminal &tm)
{
   assert<NotNewFile>(head()==FileMem::nullPtr,__FILE__,__LINE__);
   std::ifstream f(ops.com_front());
   assert<CannotOpen>(f.is_open(),__FILE__,__LINE__);
   std::string line;
   std::getline(f,line);
   std::istringstream s(line);
   int count = 0;
   std::string buf;
   while (s >> buf) ++count;
   tm << count << "\n";
}
