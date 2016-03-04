#ifndef FSTRING_H
#define FSTRING_H
#include "filemem.h"



namespace FStringData
{
   struct Header;
   struct String;
   constexpr FileMem::SizeT hdrSz = 3;
}

struct FStringData::Header : FileMem::Static<hdrSz>
{
   using Static<hdrSz>::Static;
   char& stripe() { get<char>(0); }
   uint16_t& size() { get<uint16_t>(1); }
};

struct FStringData::String : FileMem::Object
{
   using FSizeT = FileMem::SizeT;
   using Object::Object;
   String(FSizeT size): Object(size) {}
   char* c_str() { &get<char>(0); }
};



class FString
{
public:
   struct Exception;
   struct AlreadySet;
   using string = std::string;
   using FPtr = FileMem::Ptr;
   FString(const FString&) = delete;
   FString& operator=(const FString&) = delete;
   //
   FString(FString&&);
   FString& operator=(FString&&);
   FString(FileMem*);
   FString(FileMem*,FPtr);
   void set(const string&);
   const string& operator*();
   const string* operator->();
   FString& operator=(const string&);
private:
   using Header = FStringData::Header;
   using String = FStringData::String;
   FileMem* _mem;
   Header _hdr;
   string _str;
};



#endif
