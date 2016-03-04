#ifndef FSTRING_H
#define FSTRING_H
#include "filemem.h"



namespace FStringData
{
   struct Header;
   struct String;
   constexpr FileMem::SizeT hdrSz = 3;
   constexpr uint8_t strip = 170;
}

struct FStringData::Header : FileMem::Static<hdrSz>
{
   using Static<hdrSz>::Static;
   uint8_t& stripe() { get<uint8_t>(0); }
   uint16_t& sSize() { get<uint16_t>(1); }
};

struct FStringData::String : FileMem::Object
{
   using FPtr = FileMem::Ptr;
   using FSizeT = FileMem::SizeT;
   String(FSizeT size, FPtr ptr = FileMem::nullPtr): Object(size,ptr) {}
   char* c_str() { &get<char>(0); }
};



class FString
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct InvalidPtr;
   struct AlreadySet;
   // *
   // * DECLERATIONS
   // *
   using string = std::string;
   using FPtr = FileMem::Ptr;
   // *
   // * BASIC METHODS
   // *
   FString(FileMem*,FPtr = FileMem::nullPtr);
   // *
   // * COPY METHODS
   // *
   FString(const FString&) = delete;
   FString& operator=(const FString&) = delete;
   // *
   // * MOVE METHODS
   // *
   FString(FString&&);
   FString& operator=(FString&&);
   // *
   // * OPERATORS
   // *
   const string& operator*();
   const string* operator->();
   FString& operator=(const string&);
private:
   // *
   // * DECLERATIONS
   // *
   using Header = FStringData::Header;
   using String = FStringData::String;
   // *
   // * VARIABLES
   // *
   FileMem* _mem;
   Header _hdr;
   string _str;
};



struct FString::Exception : public ::Exception
{
   using ::Exception::Exception;
};

struct FString::InvalidPtr : public FString::Exception
{
   InvalidPtr(const char* file, int line):
      Exception(file,line,"FString::InvalidPtr")
   {}
};

struct FString::AlreadySet : public FString::Exception
{
   AlreadySet(const char* file, int line):
      Exception(file,line,"FString::AlreadySet")
   {}
};



#endif
