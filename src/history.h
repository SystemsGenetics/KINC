#ifndef HISTORY_H
#define HISTORY_H
#include <string>
#include "histitem.h"
#include "exception.h"



class History
{
public:
   // *
   // * DECLERATIONS
   // *
   class Iterator;
   using FPtr = FileMem::Ptr;
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
   FileMem::Ptr addr();
   void add_child(const History&);
   Iterator begin();
   Iterator end();
   // *
   // * OPERATORS
   // *
   HistItem& operator*();
   HistItem* operator->();
private:
   // *
   // * VARIABLES
   // *
   FileMem& _mem;
   HistItem _head;
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
   Iterator childHead();
   // *
   // * OPERATORS
   // *
   HistItem& operator*();
   HistItem* operator->();
   void operator++();
   bool operator!=(const Iterator&);
private:
   // *
   // * BASIC METHODS
   // *
   Iterator(FileMem*,FPtr = FileMem::nullPtr);
   // *
   // * VARIABLES
   // *
   HistItem _item;
};



inline FileMem::Ptr History::addr()
{
   return _head.addr();
}



inline History::Iterator History::begin()
{
   return {&_mem,_head.childHead()};
}



inline History::Iterator History::end()
{
   return {&_mem};
}



inline HistItem& History::operator*()
{
   return _head;
}



inline HistItem* History::operator->()
{
   return &_head;
}



inline History::Iterator History::Iterator::childHead()
{
   return {_item.mem(),_item.childHead()};
}



inline HistItem& History::Iterator::operator*()
{
   return _item;
}



inline HistItem* History::Iterator::operator->()
{
   return &_item;
}



inline void History::Iterator::operator++()
{
   if (_item.addr()!=FileMem::nullPtr)
   {
      _item = _item.next();
   }
}



inline bool History::Iterator::operator!=(const Iterator& cmp)
{
   return _item.addr()!=cmp._item.addr();
}



inline History::Iterator::Iterator(FileMem* mem, FPtr ptr):
   _item(mem,ptr)
{}



#endif
