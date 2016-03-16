#ifndef UNIT_H
#define UNIT_H
#include <iostream>



namespace unit
{
   namespace getopts
   {
      bool main();
      bool init1();
      bool com_get1();
      bool com_get2();
      bool com_front1();
      bool com_pop1();
      bool com_empty1();
      bool size1();
      bool empty1();
      bool iterate1();
      bool erase1();
      bool key1();
      bool iter_empty1();
      bool operate1();
      bool operate2();
      bool operate3();
      bool operate4();
      bool operate5();
      bool operate6();
      bool operate7();
      bool operate8();
      bool operate9();
      bool operate10();
      bool operate11();
      bool operate12();
      bool operate13();
      bool operate14();
      bool operate15();
      bool operate16();
      bool operate17();
      bool operate18();
   }
   namespace filemem
   {
      bool main();
      bool static_();
      bool object();
      bool construct();
      bool reserve();
      bool expand();
      bool capacity();
      bool allocate();
      bool allot();
      bool clear();
      bool addr();
      bool sync();
   }
   namespace fstring
   {
      bool main();
      bool construct();
      bool move();
      bool addr();
      bool operat_deref();
      bool operat_fp();
      bool operat_set();
      bool final();
   }
   namespace histitem
   {
      bool main();
      bool construct();
      bool move();
      bool allocate();
      bool sync();
      bool timestamp();
      bool filename();
      bool object();
      bool command();
      bool next();
      bool childhead();
      bool operat_set();
      bool copy_from();
      bool mem();
      bool final();
   }
   namespace history
   {
      bool main();
      bool construct();
      bool operat_fp();
      bool add_child();
      bool iterate();
      bool operat_deref();
      bool iter_operat_deref();
      bool iter_operat_fp();
      bool iter_childhead();
   }
   namespace kincfile
   {
      bool main();
      bool construct();
      bool history();
      bool ident();
      bool head();
   }
   namespace datamap
   {
      bool main();
      bool construct();
      bool open();
      bool close();
      bool select();
      bool load();
      bool dump();
      bool query();
      bool find();
      bool iterate();
      bool iter_file();
      bool iter_type();
   }
   inline void initiate();
   inline void complete();
   inline void header(const char*);
   inline void end();
   inline void start();
   inline bool finish(bool,const char*);
   extern int numTestsDone;
   extern const char* headerStr;
}



inline void unit::initiate()
{
   numTestsDone = 0;
   std::cout << "Initiating Unit Tests >>>" << std::endl;
}



inline void unit::complete()
{
   std::cout << numTestsDone << " unit test(s) passed <<<" << std::endl;
}



inline void unit::header(const char* hdr)
{
   headerStr = hdr;
   std::cout << headerStr;
}



inline void unit::end()
{
   std::cout << std::endl;
}



inline void unit::start()
{
   numTestsDone++;
   std::cout << "." << std::flush;
}



inline bool unit::finish(bool cond, const char* name)
{
   if (!cond)
   {
      std::cout << std::endl << "unit::" << headerStr << "::" << name
                << " FAILED." << std::endl;
   }
   return cond;
}



#endif
