#include "history.h"



/// @brief Initialize history of file.
///
/// Initialize history data for file memory object from given file location. If
/// nullptr is given for location then a new history data structure is created
/// for the file.
///
/// @param mem File memory object.
/// @param ptr Location in file memory where history data is location, or
/// nullptr is history structure is to be created.
History::History(FileMem& mem, FPtr ptr):
   HistItem(mem,ptr)
{
   if (ptr==FileMem::nullPtr)
   {
      allocate();
      sync();
   }
}



/// @brief Add history child.
///
/// Add a new history child for this object by making a recursive copy of the
/// given history object.
///
/// @param child History object to make a child copy of.
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



/// Iterate to next history item in list of children.
void History::Iterator::operator++()
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



/// Initialize history iterator.
///
/// Initialize history item iterator from given file memory object at given
/// file location.
///
/// @param mem Pointer to file memory object.
/// @param ptr File location of history item iterator will point to.
History::Iterator::Iterator(FileMem* mem, FPtr ptr):
   _mem(mem),
   _skim(ptr)
{
   if (ptr!=FileMem::nullPtr)
   {
      _mem->sync(_skim,FileSync::read);
   }
}
