#ifndef HISTITEM_H
#define HISTITEM_H
#include <string>
#include "filemem.h"
#include "fstring.h"
#include "exception.h"



namespace HistItemData
{
   struct Skim;
   struct Item;
   constexpr FileMem::SizeT skimSz = 16;
   constexpr FileMem::SizeT nodeSz = 48;
}

struct HistItemData::Skim : FileMem::Static<skimSz>
{
   using FPtr = FileMem::Ptr;
   using Static<skimSz>::Static;
   FPtr& childHead() { get<FPtr>(0); }
   FPtr& next() { get<FPtr>(8); }
};

struct HistItemData::Item : FileMem::Static<nodeSz>
{
   using FPtr = FileMem::Ptr;
   using Static<nodeSz>::Static;
   FPtr& childHead() { get<FPtr>(0); }
   FPtr& next() { get<FPtr>(8); }
   int64_t& timeStamp() { get<int64_t>(16); }
   FPtr& fileNamePtr() { get<FPtr>(24); }
   FPtr& objectPtr() { get<FPtr>(32); }
   FPtr& commandPtr() { get<FPtr>(40); }
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
   HistItem(FileMem*,FPtr = FileMem::nullPtr);
   // *
   // * COPY METHODS
   // *
   HistItem(const HistItem&) = delete;
   HistItem& operator=(const HistItem&) = delete;
   // *
   // * MOVE METHODS
   // *
   HistItem(HistItem&&);
   HistItem& operator=(HistItem&&);
   // *
   // * FUNCTIONS
   // *
   void allocate();
   void copy_from(const HistItem&);
   void sync();
   void timeStamp(TimeT);
   TimeT timeStamp() const;
   void fileName(const string&);
   const string& fileName() const;
   void object(const string&);
   const string& object() const;
   void command(const string&);
   const string& command() const;
   void next(FPtr);
   FileMem::Ptr next() const;
   void childHead(FPtr);
   FPtr childHead() const;
   FPtr addr() const;
   FileMem* mem() const;
   // *
   // * OPERATORS
   // *
   void operator=(FPtr);
private:
   // *
   // * DECLERATIONS
   // *
   using Item = HistItemData::Item;
   using FSizeT = FileMem::SizeT;
   // *
   // * FUNCTIONS
   // *
   void load_item();
   FPtr rec_add_item(FileMem*,FPtr);
   // *
   // * VARIABLES
   // *
   FileMem* _mem;
   mutable Item _item;
   FString _fileName;
   FString _object;
   FString _command;
};



inline HistItem::HistItem(FileMem& mem, FileMem::Ptr ptr):
   HistItem(&mem,ptr)
{}



inline FileMem::Ptr HistItem::addr() const
{
   return _item.addr();
}



inline FileMem* HistItem::mem() const
{
   return _mem;
}



struct HistItem::Exception : public ::Exception
{
   using ::Exception::Exception;
};

struct HistItem::AlreadySet : public HistItem::Exception
{
   AlreadySet(const char* file, int line):
      Exception(file,line,"HistItem::AlreadySet")
   {}
};

struct HistItem::IsAllocated : public HistItem::Exception
{
   IsAllocated(const char* file, int line):
      Exception(file,line,"HistItem::IsAllocated")
   {}
};

struct HistItem::IsNullPtr : public HistItem::Exception
{
   IsNullPtr(const char* file, int line):
      Exception(file,line,"HistItem::IsNullPtr")
   {}
};

struct HistItem::InvalidItem : public HistItem::Exception
{
   InvalidItem(const char* file, int line):
      Exception(file,line,"HistItem::InvalidItem")
   {}
};



#endif
