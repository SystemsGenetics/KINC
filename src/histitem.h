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
   const FPtr& childHead() const { get<FPtr>(0); }
   FPtr& next() { get<FPtr>(8); }
   const FPtr& next() const { get<FPtr>(8); }
   int64_t& timeStamp() { get<int64_t>(16); }
   const int64_t& timeStamp() const { get<int64_t>(16); }
   FPtr& fileNamePtr() { get<FPtr>(24); }
   FPtr& objectPtr() { get<FPtr>(32); }
   FPtr& commandPtr() { get<FPtr>(40); }
};



/// @brief History item in file memory.
///
/// Represents a single history item in file memory. Has ability to create a new
/// history item or recursively make a copy of a given history item, adding
/// additional history items of all of the copied items children.
///
/// @author Josh Burns
/// @date 24 March 2016
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
   /// File memory object where history item is located.
   FileMem* _mem;
   /// File memory chunk of history item data.
   Item _item;
   /// File name of item.
   FString _fileName;
   /// Object that created item, if any.
   FString _object;
   /// Command that created item.
   FString _command;
};



/// Passes initialization to primary constructor.
///
/// @param mem File memory that will be used for item.
/// @param ptr Location where history item is located or nullptr if item to be
/// created.
inline HistItem::HistItem(FileMem& mem, FileMem::Ptr ptr):
   HistItem(&mem,ptr)
{}



/// Get file memory location of history item.
///
/// @return Location of item.
inline FileMem::Ptr HistItem::addr() const
{
   return _item.addr();
}



/// Get pointer to file memory object where item is located.
///
/// @return Pointer to file memory instance.
inline FileMem* HistItem::mem() const
{
   return _mem;
}



/// Generic base exception class for all exceptions thrown in HistItem class.
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
