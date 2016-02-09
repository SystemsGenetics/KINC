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
      bool clear1();
      bool reserve1();
      bool allocate1();
      bool available1();
      bool head1();
      bool read1();
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
      bool sync1();
      bool sync2();
      bool operat1();
   }
   inline void header(const char*);
   inline void end();
   inline void start();
   inline bool finish(bool,const char*);
   extern const char* headerStr;
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
