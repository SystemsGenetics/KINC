#ifndef HISTITEM_H
#define HISTITEM_H
#include <string>
#include "filemem.h"
#include "exception.h"



namespace HistItemData
{
   struct Item;
   struct String;
   constexpr static int nodeLen = 60;
}

struct HistItemData::Item : FileMem::Static<nodeLen>
{
   using Static<nodeLen>::Static;
   int64_t& timeStamp() { get<int64_t>(0); }
   FileMem::Ptr& fileNamePtr() { get<FileMem::Ptr>(8); }
   int32_t& fileNameSize() { get<int32_t>(16); }
   FileMem::Ptr& objectPtr() { get<FileMem::Ptr>(20); }
   int32_t& objectSize() { get<int32_t>(28); }
   FileMem::Ptr& commandPtr() { get<FileMem::Ptr>(32); }
   int32_t& commandSize() { get<int32_t>(40); }
   FileMem::Ptr& childHead() { get<FileMem::Ptr>(44); }
   FileMem::Ptr& next() { get<FileMem::Ptr>(52); }
};

struct HistItemData::String : FileMem::Object
{
   using Object::Object;
   String(FileMem::SizeT size): Object(size) {}
   char* c_str() { &get<char>(0); }
};



class HistItem
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct AlreadySet;
   struct IsAllocated;
   struct IsNullPtr;
   // *
   // * BASIC METHODS
   // *
   HistItem(HistItem&&) = delete;
   HistItem& operator=(HistItem&&) = delete;
   HistItem(const HistItem&) = delete;
   HistItem& operator=(const HistItem&);
   HistItem(FileMem&,FileMem::Ptr = FileMem::nullPtr);
   void operator=(FileMem::Ptr);
   // *
   // * FUNCTIONS
   // *
   void allocate();
   void sync();
   void timeStamp(int64_t);
   int64_t timeStamp() const;
   void fileName(const std::string&);
   std::string fileName() const;
   void object(const std::string&);
   std::string object() const;
   void command(const std::string&);
   std::string command() const;
   void next(FileMem::Ptr);
   FileMem::Ptr next() const;
   void childHead(FileMem::Ptr);
   FileMem::Ptr childHead() const;
   FileMem::Ptr addr() const;
private:
   // *
   // * DECLERATIONS
   // *
   using Item = HistItemData::Item;
   using String = HistItemData::String;
   // *
   // * FUNCTIONS
   // *
   void load_item();
   std::string get_string(FileMem::Ptr,FileMem::SizeT);
   FileMem::Ptr set_string(const std::string&);
   FileMem::Ptr rec_add_item(FileMem&,FileMem::Ptr);
   // *
   // * CONSTANTS
   // *
   constexpr static int _nodeLen = HistItemData::nodeLen;
   // *
   // * VARIABLES
   // *
   FileMem& _mem;
   mutable Item _item;
   std::string _fileName;
   std::string _object;
   std::string _command;
};



struct HistItem::Exception : public ::Exception
{
   using ::Exception::Exception;
};

struct HistItem::AlreadySet : public HistItem::Exception
{
   AlreadySet(const char* file, int line):
      Exception(file,line,"History::AlreadySet")
   {}
};

struct HistItem::IsAllocated : public HistItem::Exception
{
   IsAllocated(const char* file, int line):
      Exception(file,line,"History::IsAllocated")
   {}
};

struct HistItem::IsNullPtr : public HistItem::Exception
{
   IsNullPtr(const char* file, int line):
      Exception(file,line,"History::IsNullPtr")
   {}
};



#endif
