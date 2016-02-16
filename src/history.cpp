#include "history.h"



History::History(FileMem& mem, FileMem::Ptr ptr):
   _mem(mem),
   _head(mem,ptr)
{
   if (ptr==FileMem::nullPtr)
   {
      _head.allocate();
      _head.sync();
   }
}



void History::timeStamp(int64_t ts)
{
   _head.timeStamp(ts);
}



void History::fileName(const std::string& fileName)
{
   _head.fileName(fileName);
}



void History::object(const std::string& object)
{
   _head.object(object);
}



void History::command(const std::string& command)
{
   _head.command(command);
}



void History::add_child(const History& child)
{
   HistItem nChild(_mem);
   nChild = child._head;
   if (_head.childHead()==FileMem::nullPtr)
   {
      _head.childHead(nChild.addr());
      _head.sync();
   }
   else
   {
      HistItem tmp(_mem,_head.childHead());
      while (tmp.next()!=FileMem::nullPtr)
      {
         tmp = tmp.next();
      }
      tmp.next(nChild.addr());
      tmp.sync();
   }
}



History::Iterator History::begin()
{
   return {_mem,_head.childHead()};
}



History::Iterator History::end()
{
   return {_mem};
}



History::Iterator History::Iterator::childHead()
{
   return {_mem,_item.childHead()};
}



HistItem& History::Iterator::operator*()
{
   return _item;
}



void History::Iterator::operator++()
{
   if (_item.addr()!=FileMem::nullPtr)
   {
      _item = _item.next();
   }
}



bool History::Iterator::operator!=(const Iterator& cmp)
{
   return _item.addr()!=cmp._item.addr();
}



History::Iterator::Iterator(FileMem& mem, FileMem::Ptr ptr):
   _mem(mem),
   _item(mem,ptr)
{}
