#ifndef HISTORY_H
#define HISTORY_H
#include <string>
#include "histitem.h"
#include "exception.h"



class History : private HistItem
{
public:
   // *
   // * DECLERATIONS
   // *
   class Iterator;
   using FPtr = FileMem::Ptr;
   using HistItem::sync;
   using HistItem::timeStamp;
   using HistItem::fileName;
   using HistItem::object;
   using HistItem::command;
   using HistItem::addr;
   // *
   // * BASIC METHODS
   // *
   History(FileMem&,FPtr = FileMem::nullPtr);
   // *
   // * COPY METHODS
   // *
   History(const History&) = delete;
   History& operator=(const History&) = delete;
   // *
   // * MOVE METHODS
   // *
   History(History&&) = delete;
   History& operator=(History&&) = delete;
   // *
   // * FUNCTIONS
   // *
   void add_child(const History&);
   bool has_child() const;
   Iterator begin();
   Iterator end();
};



class History::Iterator
{
public:
   // *
   // * DECLERATIONS
   // *
   friend class History;
   using FPtr = History::FPtr;
   // *
   // * FUNCTIONS
   // *
   Iterator child();
   bool has_child() const;
   // *
   // * OPERATORS
   // *
   HistItem load();
   void operator++();
   bool operator!=(const Iterator&);
private:
   using Skim = HistItemData::Skim;
   // *
   // * BASIC METHODS
   // *
   Iterator(FileMem*,FPtr = FileMem::nullPtr);
   // *
   // * VARIABLES
   // *
   FileMem* _mem;
   mutable Skim _skim;
};



inline bool History::has_child() const
{
   return childHead()!=FileMem::nullPtr;
}



inline History::Iterator History::begin()
{
   return {mem(),childHead()};
}



inline History::Iterator History::end()
{
   return {mem()};
}



inline History::Iterator History::Iterator::child()
{
   return {_mem,_skim.childHead()};
}



inline bool History::Iterator::has_child() const
{
   return _skim.childHead()!=FileMem::nullPtr;
}



inline HistItem History::Iterator::load()
{
   return {_mem,_skim.addr()};
}



inline void History::Iterator::operator++()
{
   if (_skim.addr()!=FileMem::nullPtr)
   {
      _skim = _skim.next();
      if (_skim.addr()!=FileMem::nullPtr)
      {
         _mem->sync(_skim,FileSync::read);
      }
   }
}



inline bool History::Iterator::operator!=(const Iterator& cmp)
{
   return _skim.addr()!=cmp._skim.addr();
}



inline History::Iterator::Iterator(FileMem* mem, FPtr ptr):
   _mem(mem),
   _skim(ptr)
{
   if (ptr!=FileMem::nullPtr)
   {
      _mem->sync(_skim,FileSync::read);
   }
}



#endif
