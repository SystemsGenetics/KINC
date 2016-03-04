#ifndef KINCFILE_H
#define KINCFILE_H
#include <string>
#include "filemem.h"
#include "history.h"



namespace KincFileData
{
   struct Header;
   struct String;
   constexpr FileMem::SizeT idSz = 10;
   constexpr FileMem::SizeT hdrSz = idSz+28;
   constexpr auto idString = "\030\031\032\113\111\116\103\032\031\030";
}

struct KincFileData::Header : FileMem::Static<hdrSz>
{
   using FPtr = FileMem::Ptr;
   using FSizeT = FileMem::SizeT;
   using Static<hdrSz>::Static;
   char* idString() { &get<char>(0); }
   FPtr& ident() { get<FPtr>(idSz); }
   uint32_t identSize() { get<uint32_t>(idSz+8); }
   FPtr& histHead() { get<FPtr>(idSz+12); }
   FPtr& dataHead() { get<FPtr>(idSz+20); }
};

struct KincFileData::String : FileMem::Object
{
   using FSizeT = FileMem::SizeT;
   using Object::Object;
   String(FSizeT size): Object(size) {}
   char* c_str() { &get<char>(0); }
};



class KincFile
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct InvalidFile;
   struct AlreadySet;
   // *
   // * DECLERATIONS
   // *
   using string = std::string;
   using FPtr = FileMem::Ptr;
   // *
   // * BASIC METHODS
   // *
   KincFile(const string&);
   ~KincFile();
   // *
   // * COPY METHODS
   // *
   KincFile(const KincFile&) = delete;
   KincFile& operator=(const KincFile&) = delete;
   // *
   // * MOVE METHODS
   // *
   KincFile(KincFile&&) = delete;
   KincFile& operator=(KincFile&&) = delete;
   // *
   // * FUNCTIONS
   // *
   bool is_new();
   History& history();
protected:
   // *
   // * FUNCTIONS
   // *
   string ident() const;
   void ident(const string&);
   FPtr head() const;
   void head(FPtr);
private:
   // *
   // * DECLERATIONS
   // *
   using Header = KincFileData::Header;
   // *
   // * CONSTANTS
   // *
   constexpr static auto _idString {KincFileData::idString};
   constexpr static auto _idSz {KincFileData::idSz};
   constexpr static auto _hdrSz {KincFileData::hdrSz};
   // *
   // * VARIABLES
   // *
   bool _new {true};
   FileMem* _mem {nullptr};
   History* _hist {nullptr};
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

struct KincFile::AlreadySet : public KincFile::Exception
{
   AlreadySet(const char* file, int line):
      Exception(file,line,"KincFile::AlreadySet")
   {}
};



#endif
