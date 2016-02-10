#include <cstdlib>
#include <iostream>
#include "linuxfile.h"
#include "unit.h"



class PublicLinuxFile : public LinuxFile
{
public:
   using LinuxFile::LinuxFile;
   using LinuxFile::write;
   using LinuxFile::read;
};



namespace unit
{
   namespace linuxfile
   {
      constexpr const char* tmpFile = "testfiles/memfiletmp";
      constexpr const char* validFile = "testfiles/memfile";
      constexpr const char* invalidFile = "testfiles/notmemfile";
      const char* wData = "1234567890123456";
   }
}



bool unit::linuxfile::main()
{
   bool ret = false;
   header("LinuxFile");
   try
   {
      ret = init1()&&
            init2()&&
            init3()&&
            reserve1()&&
            clear1()&&
            available1()&&
            head1()&&
            read1();
   }
   catch (...)
   {
      system("rm -f testfiles/memfiletmp");
      end();
      throw;
   }
   system("rm -f testfiles/memfiletmp");
   end();
   return ret;
}



bool unit::linuxfile::init1()
{
   start();
   LinuxFile t(tmpFile);
   bool ret = t.size()==0;
   return finish(ret,"init1");
}



bool unit::linuxfile::init2()
{
   start();
   LinuxFile t(validFile);
   bool ret = t.size()==0;
   return finish(ret,"init2");
}



bool unit::linuxfile::init3()
{
   start();
   bool ret = false;
   try
   {
      LinuxFile t(invalidFile);
      t.size();
   }
   catch (LinuxFile::InvalidFile)
   {
      ret = true;
   }
   return finish(ret,"init3");
}



bool unit::linuxfile::reserve1()
{
   start();
   bool ret = false;
   {
      LinuxFile t(tmpFile);
      ret = t.reserve(128);
      ret = ret&&t.size()==128&&t.capacity()==128;
   }
   {
      LinuxFile t(tmpFile);
      ret = ret&&t.size()==128&&t.capacity()==128;
   }
   return finish(ret,"reserve1");
}



bool unit::linuxfile::clear1()
{
   start();
   bool ret = false;
   {
      LinuxFile t(tmpFile);
      t.allocate(64);
      t.clear();
      ret = t.size()==128&&t.capacity()==128;
   }
   {
      LinuxFile t(tmpFile);
      ret = ret&&t.size()==128&&t.capacity()==128;
   }
   return finish(ret,"reserve1");
}



bool unit::linuxfile::available1()
{
   start();
   LinuxFile t(tmpFile);
   t.allocate(64);
   bool ret = t.size()==128&&t.capacity()==64;
   return finish(ret,"available1");
}



bool unit::linuxfile::head1()
{
   start();
   LinuxFile t(tmpFile);
   t.clear();
   bool ret = t.head()==t.allocate(64);
   return finish(ret,"head1");
}



bool unit::linuxfile::read1()
{
   start();
   bool ret = false;
   {
      PublicLinuxFile t(tmpFile);
      t.write(wData,0,16);
   }
   {
      char buf[17] {'\0'};
      PublicLinuxFile t(tmpFile);
      t.read(buf,0,16);
      ret = strcmp(buf,wData)==0;
   }
   return finish(ret,"read1");
}
