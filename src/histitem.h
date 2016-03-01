#ifndef HISTITEM_H
#define HISTITEM_H
#include <string>
#include "filemem.h"
#include "exception.h"



namespace HistItemData
{
   struct Item;
   struct String;
   constexpr FileMem::SizeT nodeSz = 60;
}

struct HistItemData::Item : FileMem::Static<nodeSz>
{
   using FPtr = FileMem::Ptr;
   using Static<nodeSz>::Static;
   int64_t& timeStamp() { get<int64_t>(0); }
   FPtr& fileNamePtr() { get<FPtr>(8); }
   int32_t& fileNameSize() { get<int32_t>(16); }
   FPtr& objectPtr() { get<FPtr>(20); }
   int32_t& objectSize() { get<int32_t>(28); }
   FPtr& commandPtr() { get<FPtr>(32); }
   int32_t& commandSize() { get<int32_t>(40); }
   FPtr& childHead() { get<FPtr>(44); }
   FPtr& next() { get<FPtr>(52); }
};

struct HistItemData::String : FileMem::Object
{
   using FSizeT = FileMem::SizeT;
   using Object::Object;
   String(FSizeT size): Object(size) {}
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
   struct InvalidItem;
   // *
   // * DECLERATIONS
   // *
   using string = std::string;
   using FPtr = FileMem::Ptr;
   using TimeT = int64_t;
   // *
   // * BASIC METHODS
   // *
   HistItem(FileMem&,FPtr = FileMem::nullPtr);
   // *
   // * COPY METHODS
   // *
   HistItem(const HistItem&) = delete;
   HistItem& operator=(const HistItem&) = delete;
   // *
   // * MOVE METHODS
   // *
   HistItem(HistItem&&) = delete;
   HistItem& operator=(HistItem&&) = delete;
   // *
   // * FUNCTIONS
   // *
   void allocate();
   void copy_from(const HistItem&);
   void sync();
   void timeStamp(TimeT);
   TimeT timeStamp() const;
   void fileName(const string&);
   string fileName() const;
   void object(const string&);
   string object() const;
   void command(const string&);
   string command() const;
   void next(FPtr);
   FileMem::Ptr next() const;
   void childHead(FPtr);
   FPtr childHead() const;
   FPtr addr() const;
   // *
   // * FUNCTIONS
   // *
   void operator=(FPtr);
private:
   // *
   // * DECLERATIONS
   // *
   using Item = HistItemData::Item;
   using String = HistItemData::String;
   using FSizeT = FileMem::SizeT;
   // *
   // * FUNCTIONS
   // *
   inline void load_item();
   inline string get_string(FPtr,FSizeT);
   inline FPtr set_string(const string&);
   FPtr rec_add_item(FileMem&,FPtr);
   // *
   // * VARIABLES
   // *
   FileMem& _mem;
   mutable Item _item;
   string _fileName;
   string _object;
   string _command;
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

struct HistItem::InvalidItem : public HistItem::Exception
{
   InvalidItem(const char* file, int line):
      Exception(file,line,"History::InvalidItem")
   {}
};



#endif
