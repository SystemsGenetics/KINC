#ifndef UNIT_H
#define UNIT_H
#include <iostream>



namespace unit
{
   namespace getopts
   {
      bool main();
      bool construct();
      bool orig();
      bool com_size();
      bool com_empty();
      bool com_get();
      bool com_front();
      bool com_pop();
      bool com_empty();
      bool size();
      bool empty();
      bool has_opt();
      bool iterate();
      bool erase();
      bool iter_key();
      bool iter_value();
      bool iter_is_key();
      bool iter_val_empty();
      bool iter_operat_equal();
   }
   namespace filemem
   {
      bool main();
      bool static_();
      bool object();
      bool construct();
      bool reserve();
      bool expand();
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
      bool add_child();
      bool has_child();
      bool iterate();
      bool iter_load();
      bool iter_has_child();
      bool iter_child();
   }
   namespace kincfile
   {
      bool main();
      bool construct();
      bool history();
      bool clear();
      bool is_new();
      bool ident();
      bool head();
   }
   namespace datamap
   {
      bool main();
      bool construct();
      bool open();
      bool close();
      bool unselect();
      bool load();
      bool dump();
      bool query();
      bool iterate();
      bool selected();
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
   std::cout << "\nInitiating Unit Tests >>>\n\n";
}



inline void unit::complete()
{
   std::cout << numTestsDone << " unit test(s) passed <<<\n" << std::endl;
}



inline void unit::header(const char* hdr)
{
   headerStr = hdr;
   std::cout << headerStr << "\n";
}



inline void unit::end()
{
   std::cout << std::endl;
}



inline void unit::start()
{
   numTestsDone++;
   std::cout << ".";
}



inline bool unit::finish(bool cond, const char* name)
{
   std::cout << "  " << name << "\n";
   if (!cond)
   {
      std::cout << std::endl << "unit::" << headerStr << "::" << name
                << " FAILED." << std::endl;
   }
   return cond;
}



#endif
