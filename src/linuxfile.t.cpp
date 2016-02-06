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
   using LinuxFile::allocate;
};



const char* writeData = "1234567890123456";



bool unit::linuxfile::main()
{
   bool ret = false;
   std::cout << "LinuxFile";
   try
   {
      ret = init1()&&
            init2()&&
            init3()&&
            reserve1()&&
            reserve2()&&
            allocate1()&&
            write_read1()&&
            write_read2()&&
            allocate2()&&
            write1()&&
            read1()&&
            clear1()&&
            allocate3();
   }
   catch (...)
   {
      system("rm -f testfiles/kincdat2");
      std::cout << std::endl;
      throw;
   }
   system("rm -f testfiles/kincdat2");
   std::cout << std::endl;
   return ret;
}



bool unit::linuxfile::init1()
{
   std::cout << "." << std::flush;
   LinuxFile t("testfiles/kincdat1");
   bool ret = t.size()==0&&t.available()==0;
   if (!ret)
   {
      std::cout << std::endl << "unit::linuxfile::init1 FAILED." << std::endl;
   }
   return ret;
}



bool unit::linuxfile::init2()
{
   std::cout << "." << std::flush;
   LinuxFile t("testfiles/kincdat2");
   bool ret = t.size()==0&&t.available()==0;
   if (!ret)
   {
      std::cout << std::endl << "unit::linuxfile::init2 FAILED." << std::endl;
   }
   return ret;
}



bool unit::linuxfile::init3()
{
   std::cout << "." << std::flush;
   bool ret = false;
   try
   {
      LinuxFile t("testfiles/notkincdat");
      t.size();
   }
   catch (FileMem::InvalidFile)
   {
      ret = true;
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::linuxfile::init3 FAILED." << std::endl;
   }
   return ret;
}



bool unit::linuxfile::reserve1()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      LinuxFile t("testfiles/kincdat2");
      ret = t.reserve(128);
      ret = ret&&t.size()==128&&t.available()==128;
   }
   {
      LinuxFile t("testfiles/kincdat2");
      ret = ret&&t.size()==128&&t.available()==128;
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::linuxfile::reserve1 FAILED."
                << std::endl;
   }
   return ret;
}



bool unit::linuxfile::reserve2()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      LinuxFile t("testfiles/kincdat2");
      ret = t.reserve(128);
      ret = ret&&t.size()==256&&t.available()==256;
   }
   {
      LinuxFile t("testfiles/kincdat2");
      ret = ret&&t.size()==256&&t.available()==256;
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::linuxfile::reserve2 FAILED."
                << std::endl;
   }
   return ret;
}



bool unit::linuxfile::allocate1()
{
   std::cout << "." << std::flush;
   bool ret = false;
   FileMem::VPtr ptr;
   {
      PublicLinuxFile t("testfiles/kincdat2");
      ptr = t.allocate(16);
      ret = ptr==t.head()&&t.available()==(256 - 16);
   }
   {
      PublicLinuxFile t("testfiles/kincdat2");
      ret = ret&&t.allocate(16)==(t.head() + 16)&&t.available()==(256 - 32);
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::linuxfile::allocate1 FAILED."
                << std::endl;
   }
   return ret;
}



bool unit::linuxfile::write_read1()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      char buffer[17] {'\0'};
      PublicLinuxFile t("testfiles/kincdat2");
      t.write(writeData,0,16);
      t.read(buffer,0,16);
      ret = strcmp(writeData,buffer)==0;
   }
   {
      char buffer[17] {'\0'};
      PublicLinuxFile t("testfiles/kincdat2");
      t.read(buffer,0,16);
      ret = ret&&strcmp(writeData,buffer)==0;
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::linuxfile::write_read1 FAILED."
                << std::endl;
   }
   return ret;
}



bool unit::linuxfile::write_read2()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      char buffer[17] {'\0'};
      PublicLinuxFile t("testfiles/kincdat2");
      t.write(writeData,16,16);
      t.read(buffer,16,16);
      ret = strcmp(writeData,buffer)==0;
   }
   {
      char buffer[17] {'\0'};
      PublicLinuxFile t("testfiles/kincdat2");
      t.read(buffer,16,16);
      ret = ret&&strcmp(writeData,buffer)==0;
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::linuxfile::write_read2 FAILED."
                << std::endl;
   }
   return ret;
}



bool unit::linuxfile::allocate2()
{
   std::cout << "." << std::flush;
   bool ret = false;
   try
   {
      PublicLinuxFile t("testfiles/kincdat2");
      t.allocate(225);
   }
   catch (FileMem::OutOfMemory)
   {
      ret = true;
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::linuxfile::allocate2 FAILED."
                << std::endl;
   }
   return ret;
}



bool unit::linuxfile::write1()
{
   std::cout << "." << std::flush;
   bool ret = false;
   try
   {
      PublicLinuxFile t("testfiles/kincdat2");
      t.write(writeData,241,16);
   }
   catch (FileMem::FileSegFault)
   {
      ret = true;
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::linuxfile::write1 FAILED."
                << std::endl;
   }
   return ret;
}



bool unit::linuxfile::read1()
{
   std::cout << "." << std::flush;
   bool ret = false;
   try
   {
      char buffer[16];
      PublicLinuxFile t("testfiles/kincdat2");
      t.read(buffer,241,16);
   }
   catch (FileMem::FileSegFault)
   {
      ret = true;
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::linuxfile::read1 FAILED."
                << std::endl;
   }
   return ret;
}



bool unit::linuxfile::clear1()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      LinuxFile t("testfiles/kincdat2");
      t.clear();
      ret = t.size()==256&&t.available()==256;
   }
   {
      LinuxFile t("testfiles/kincdat2");
      ret = t.size()==256&&t.available()==256;
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::linuxfile::clear1 FAILED."
                << std::endl;
   }
   return ret;
}



bool unit::linuxfile::allocate3()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      PublicLinuxFile t("testfiles/kincdat2");
      ret = t.allocate(16)==t.head();
   }
   {
      PublicLinuxFile t("testfiles/kincdat2");
      ret = ret&&t.allocate(16)==(t.head() + 16);
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::linuxfile::allocate3 FAILED."
                << std::endl;
   }
   return ret;
}
