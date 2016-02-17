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
   // *
   // * BASIC METHODS
   // *
   History(const History&) = delete;
   History(History&&) = delete;
   History& operator=(const History&) = delete;
   History& operator=(History&&) = delete;
   History(FileMem&,FileMem::Ptr = FileMem::nullPtr);
   // *
   // * FUNCTIONS
   // *
   FileMem::Ptr addr();
   void timeStamp(int64_t);
   void fileName(const std::string&);
   void object(const std::string&);
   void command(const std::string&);
   void add_child(const History&);
   Iterator begin();
   Iterator end();
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
   // *
   // * FUNCTIONS
   // *
   Iterator childHead();
   // *
   // * OPERATORS
   // *
   HistItem& operator*();
   void operator++();
   bool operator!=(const Iterator&);
private:
   // *
   // * BASIC METHODS
   // *
   Iterator(FileMem&,FileMem::Ptr = FileMem::nullPtr);
   // *
   // * VARIABLES
   // *
   FileMem& _mem;
   HistItem _item;
};



#endif
