#include "histitem.h"



HistItem::HistItem(FileMem* mem, FileMem::Ptr ptr):
   _mem(mem),
   _item(ptr),
   _fileName(mem),
   _object(mem),
   _command(mem)
{
   if (_item.addr()!=FileMem::nullPtr)
   {
      load_item();
   }
}



HistItem::HistItem(HistItem&& tmp):
   _mem(tmp._mem),
   _item(tmp._item),
   _fileName(std::move(tmp._fileName)),
   _object(std::move(tmp._object)),
   _command(std::move(tmp._command))
{
   tmp._item = FileMem::nullPtr;
}



HistItem& HistItem::operator=(HistItem&& tmp)
{
   _mem = tmp._mem;
   _item = tmp._item;
   _fileName = std::move(tmp._fileName);
   _object = std::move(tmp._object);
   _command = std::move(tmp._command);
   tmp._item = FileMem::nullPtr;
}



void HistItem::allocate()
{
   bool cond = _item.addr()==FileMem::nullPtr;
   assert<IsAllocated>(cond,__FILE__,__LINE__);
   _mem->allot(_item);
   _item.timeStamp() = 0;
   _item.fileNamePtr() = FileMem::nullPtr;
   _item.objectPtr() = FileMem::nullPtr;
   _item.commandPtr() = FileMem::nullPtr;
   _item.childHead() = FileMem::nullPtr;
   _item.next() = FileMem::nullPtr;
}



void HistItem::copy_from(const HistItem& hist)
{
   bool cond = _item.addr()==FileMem::nullPtr;
   assert<IsAllocated>(cond,__FILE__,__LINE__);
   cond = hist.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   _mem->allot(_item);
   _fileName = hist.fileName();
   _object = hist.object();
   _command = hist.command();
   _item.timeStamp() = hist.timeStamp();
   _item.fileNamePtr() = _fileName.addr();
   _item.objectPtr() = _object.addr();
   _item.commandPtr() = _command.addr();
   _item.childHead() = rec_add_item(hist._mem,hist.childHead());
   _item.next() = rec_add_item(hist._mem,hist.next());
   _mem->sync(_item,FileSync::write);
}



void HistItem::sync()
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   _mem->sync(_item,FileSync::write);
}



void HistItem::timeStamp(int64_t n)
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   _item.timeStamp() = n;
}



int64_t HistItem::timeStamp() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return _item.timeStamp();
}



void HistItem::fileName(const std::string& str)
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   cond = _item.fileNamePtr()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _fileName = str;
   _item.fileNamePtr() = _fileName.addr();
}



const HistItem::string& HistItem::fileName() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return *_fileName;
}



void HistItem::object(const std::string& str)
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   cond = _item.objectPtr()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _object = str;
   _item.objectPtr() = _object.addr();
}



const HistItem::string& HistItem::object() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return *_object;
}



void HistItem::command(const std::string& str)
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   cond = _item.commandPtr()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _command = str;
   _item.commandPtr() = _command.addr();
}



const HistItem::string& HistItem::command() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return *_command;
}



void HistItem::next(FileMem::Ptr ptr)
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   cond = _item.next()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _item.next() = ptr;
}



FileMem::Ptr HistItem::next() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return _item.next();
}



void HistItem::childHead(FileMem::Ptr ptr)
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   cond = _item.childHead()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _item.childHead() = ptr;
}



FileMem::Ptr HistItem::childHead() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return _item.childHead();
}



void HistItem::operator=(FileMem::Ptr ptr)
{
   _item = ptr;
   _fileName.addr(FileMem::nullPtr);
   _object.addr(FileMem::nullPtr);
   _command.addr(FileMem::nullPtr);
   if (_item.addr()!=FileMem::nullPtr)
   {
      load_item();
   }
}



inline void HistItem::load_item()
{
   _mem->sync(_item,FileSync::read);
   if (_item.fileNamePtr()!=FileMem::nullPtr)
   {
      _fileName.addr(_item.fileNamePtr());
   }
   if (_item.objectPtr()!=FileMem::nullPtr)
   {
      _object.addr(_item.objectPtr());
   }
   if (_item.commandPtr()!=FileMem::nullPtr)
   {
      _command.addr(_item.commandPtr());
   }
}



FileMem::Ptr HistItem::rec_add_item(FileMem* mem, FileMem::Ptr ptr)
{
   FileMem::Ptr ret = ptr;
   if (ptr!=FileMem::nullPtr)
   {
      HistItem from(mem,ptr);
      HistItem to(_mem);
      to.copy_from(from);
      ret = to.addr();
   }
   return ret;
}
