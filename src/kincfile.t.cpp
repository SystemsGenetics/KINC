#include "kincfile.h"
#include "unit.h"



namespace unit
{
   namespace kincfile
   {
      class PublicKincFile : public KincFile
      {
      public:
         using KincFile::KincFile;
         using KincFile::ident;
         using KincFile::head;
      };
      constexpr auto tmpFile = "kincfile.tmp";
      constexpr auto tmpFile2 = "kincfile2.tmp";
      constexpr auto tmpFile3 = "kincfile3.tmp";
      constexpr auto invalidFile = "notkincfile.tmp";
      constexpr auto invalidFile2 = "notkincfile2.tmp";
      constexpr auto invalidFile3 = "notkincfile3.tmp";
      constexpr const char* identStr = "1234567890123456";
      constexpr int dataPtr = 9999;
      constexpr int tStamp = 8888;
   }
}



void construct_kincfiles()
{
   {
      FileMem tf(unit::kincfile::tmpFile);
      KincFileData::Header header;
      tf.allot(header);
      History hist(tf);
      hist.timeStamp(unit::kincfile::tStamp);
      hist.sync();
      header.histHead() = hist.addr();
      header.dataHead() = unit::kincfile::dataPtr;
      FString ident(&tf);
      ident = unit::kincfile::identStr;
      header.ident() = ident.addr();
      memcpy(header.idString(),KincFileData::idString,KincFileData::idSz);
      tf.sync(header,FileSync::write);
   }
   {
      FileMem tf(unit::kincfile::invalidFile);
      KincFileData::Header header;
      tf.allot(header);
      History hist(tf);
      header.histHead() = hist.addr();
      char buffer[KincFileData::idSz] {'\0'};
      memcpy(header.idString(),buffer,KincFileData::idSz);
      tf.sync(header,FileSync::write);
   }
   {
      FileMem tf(unit::kincfile::invalidFile2);
      KincFileData::Header header;
      tf.allot(header);
      header.histHead() = FileMem::nullPtr;
      memcpy(header.idString(),KincFileData::idString,KincFileData::idSz);
      tf.sync(header,FileSync::write);
   }
   {
      FileMem tf(unit::kincfile::invalidFile3);
      FString data(&tf);
      data = "invalid";
   }
}



bool unit::kincfile::main()
{
   bool ret = false;
   header("KincFile");
   try
   {
      construct_kincfiles();
      ret = construct()&&
            history()&&
            clear()&&
            is_new()&&
            ident()&&
            head();
   }
   catch (...)
   {
      system("rm -f *.tmp");
      unit::end();
      throw;
   }
   system("rm -f *.tmp");
   unit::end();
   return ret;
}



bool unit::kincfile::construct()
{
   bool cont = true;
   {
      start();
      KincFile t(tmpFile2);
      bool test = t.is_new();
      cont = cont&&finish(test,"construct1");
   }
   if (cont)
   {
      start();
      KincFile t(tmpFile2);
      bool test = !t.is_new();
      cont = cont&&finish(test,"construct2");
   }
   if (cont)
   {
      start();
      KincFile t(tmpFile);
      bool test = !t.is_new();
      cont = cont&&finish(test,"construct3");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         KincFile t(invalidFile);
      }
      catch (KincFile::InvalidFile)
      {
         test = true;
      }
      cont = cont&&finish(test,"construct4");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         KincFile t(invalidFile2);
      }
      catch (KincFile::InvalidFile)
      {
         test = true;
      }
      cont = cont&&finish(test,"construct5");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         KincFile t(invalidFile3);
      }
      catch (KincFile::InvalidFile)
      {
         test = true;
      }
      cont = cont&&finish(test,"construct6");
   }
   return cont;
}



bool unit::kincfile::history()
{
   start();
   KincFile t(tmpFile);
   bool test = t.history().timeStamp()==tStamp;
   return finish(test,"history1");
}



bool unit::kincfile::clear()
{
   start();
   KincFile t(tmpFile2);
   t.history().timeStamp(6666);
   t.history().sync();
   t.clear();
   bool test = t.history().timeStamp()==0;
   return finish(test,"clear1");
}



bool unit::kincfile::is_new()
{
   bool cont {true};
   {
      start();
      KincFile t(tmpFile3);
      bool test = t.is_new();
      cont = cont&&finish(test,"is_new1");
   }
   if (cont)
   {
      start();
      KincFile t(tmpFile3);
      bool test = !t.is_new();
      cont = cont&&finish(test,"is_new2");
   }
   if (cont)
   {
      start();
      KincFile t(tmpFile3);
      t.clear();
      bool test = t.is_new();
      cont = cont&&finish(test,"is_new3");
   }
   return cont;
}



bool unit::kincfile::ident()
{
   std::string ident(identStr);
   bool cont = true;
   {
      start();
      PublicKincFile t(tmpFile);
      bool test = t.ident()==ident;
      cont = cont&&finish(test,"ident1");
   }
   if (cont)
   {
      start();
      {
         PublicKincFile t(tmpFile2);
         t.ident(ident);
      }
      PublicKincFile t(tmpFile2);
      bool test = t.ident()==ident;
      cont = cont&&finish(test,"ident2");
   }
   if (cont)
   {
      start();
      PublicKincFile t(tmpFile2);
      bool test = false;
      try
      {
         t.ident(ident);
      }
      catch (KincFile::AlreadySet)
      {
         test = true;
      }
      cont = cont&&finish(test,"ident3");
   }
   return cont;
}



bool unit::kincfile::head()
{
   bool cont = true;
   {
      start();
      PublicKincFile t(tmpFile);
      bool test = t.head()==dataPtr;
      cont = cont&&finish(test,"head1");
   }
   if (cont)
   {
      start();
      {
         PublicKincFile t(tmpFile2);
         t.head(dataPtr);
      }
      PublicKincFile t(tmpFile2);
      bool test = t.head()==dataPtr;
      cont = cont&&finish(test,"head2");
   }
   return cont;
}
