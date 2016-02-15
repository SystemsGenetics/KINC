#ifndef UNIT_H
#define UNIT_H
#include <iostream>



namespace unit
{
   namespace filemem
   {
      bool main();
      bool static1();
      bool object1();
      bool object2();
      bool object3();
      bool object4();
      bool object5();
      bool init1();
      bool init2();
      bool init3();
      bool init4();
      bool reserve1();
      bool capacity1();
      bool allocate1();
      bool allocate2();
      bool allocate3();
      bool clear1();
      bool addr1();
      bool sync1();
      bool sync2();
      bool sync3();
      bool sync4();
      bool sync5();
      bool sync6();
   }
   namespace histitem
   {
      bool main();
      bool init1();
      bool allocate1();
      bool init2();
      bool allocate2();
      bool init3();
      bool init4();
      bool timestamp1();
      bool timestamp2();
      bool filename1();
      bool filename2();
      bool filename3();
      bool filename4();
      bool object1();
      bool object2();
      bool object3();
      bool object4();
      bool command1();
      bool command2();
      bool command3();
      bool command4();
      bool next1();
      bool next2();
      bool next3();
      bool next4();
      bool childhead1();
      bool childhead2();
      bool childhead3();
      bool childhead4();
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
