#ifndef KINCFILE_H
#define KINCFILE_H
#include <string>
#include "filemem.h"
#include "history.h"



namespace KincFileData
{
   constexpr auto idString = "\030\031\032\113\111\116\103\032\031\030";
   constexpr FileMem::SizeT idSz = 10;
   constexpr FileMem::SizeT hdrSz = idSz+16;
   struct Header;
}

struct KincFileData::Header : FileMem::Static<hdrSz>
{
   using Static<hdrSz>::Static;
   char* idString() { &get<char>(0); }
   FileMem::Ptr& histHead() { get<FileMem::Ptr>(idSz); }
   FileMem::Ptr& dataHead() { get<FileMem::Ptr>(idSz+8); }
};



class KincFile
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct InvalidFile;
   // *
   // * DECLERATIONS
   // *
   friend class Console;
   // *
   // * BASIC METHODS
   // *
   KincFile(const KincFile&) = delete;
   KincFile(KincFile&&) = delete;
   KincFile& operator=(const KincFile&) = delete;
   KincFile& operator=(KincFile&&) = delete;
   KincFile(const std::string&);
   ~KincFile();
protected:
   // *
   // * FUNCTIONS
   // *
   FileMem::Ptr head() const;
   void head(FileMem::Ptr);
private:
   // *
   // * DECLERATIONS
   // *
   using Header = KincFileData::Header;
   // *
   // * FUNCTIONS
   // *
   History& history();
   // *
   // * CONSTANTS
   // *
   constexpr static auto _idString = KincFileData::idString;
   constexpr static auto _idSz = KincFileData::idSz;
   constexpr static auto _hdrSz = KincFileData::hdrSz;
   // *
   // * VARIABLES
   // *
   FileMem* _mem;
   History* _hist;
   mutable Header _header;
};



struct KincFile::Exception : public ::Exception
{
   using ::Exception::Exception;
};

struct KincFile::InvalidFile : public KincFile::Exception
{
   InvalidFile(const char* file, int line):
      Exception(file,line,"KincFile::InvalidFile")
   {}
};



#endif
