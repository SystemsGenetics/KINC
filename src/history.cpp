#include "history.h"



History::History(FileMem& mem, FPtr ptr):
   HistItem(mem,ptr)
{
   if (ptr==FileMem::nullPtr)
   {
      allocate();
      sync();
   }
}



void History::add_child(const History& child)
{
   HistItem nChild(mem());
   nChild.copy_from(child);
   if (childHead()==FileMem::nullPtr)
   {
      childHead(nChild.addr());
      sync();
   }
   else
   {
      HistItem tmp(mem(),childHead());
      while (tmp.next()!=FileMem::nullPtr)
      {
         tmp = tmp.next();
      }
      tmp.next(nChild.addr());
      tmp.sync();
   }
}
