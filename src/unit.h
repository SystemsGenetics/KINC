#ifndef UNIT_H
#define UNIT_H
#include <iostream>



namespace unit
{
   namespace linuxfile
   {
      bool main();
      bool init1();
      bool init2();
      bool init3();
      bool reserve1();
      bool reserve2();
      bool allocate1();
      bool write_read1();
      bool write_read2();
      bool allocate2();
      bool write1();
      bool read1();
      bool clear1();
      bool allocate3();
   }
   namespace filemem
   {
      bool main();
      bool object1();
      bool object2();
      bool dobject1();
      bool dobject2();
      bool dobject3();
      bool dobject4();
      bool dobject5();
      bool init1();
      bool init2();
      bool init3();
      bool init4();
      bool init5();
      bool init6();
      bool addr1();
      bool raw1();
      bool raw2();
      bool save1();
      bool operat1();
      bool operat2();
      bool operat3();
      bool operat4();
      bool operat5();
      bool operat6();
      bool operat7();
      bool operat8();
      bool operat9();
   }
   inline void header(const char*);
   inline void end();
   inline void start();
   inline bool finish(bool,const char*);
   static const char* headerStr {nullptr};
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
