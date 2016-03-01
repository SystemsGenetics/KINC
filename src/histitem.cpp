#include "histitem.h"



HistItem::HistItem(FileMem& mem, FileMem::Ptr ptr):
   _mem(mem),
   _item(ptr)
{
   if (_item.addr()!=FileMem::nullPtr)
   {
      load_item();
   }
}



void HistItem::allocate()
{
   bool cond = _item.addr()==FileMem::nullPtr;
   assert<IsAllocated>(cond,__FILE__,__LINE__);
   _mem.allot(_item);
   _item.timeStamp() = 0;
   _item.fileNamePtr() = FileMem::nullPtr;
   _item.fileNameSize() = 0;
   _item.objectPtr() = FileMem::nullPtr;
   _item.objectSize() = 0;
   _item.commandPtr() = FileMem::nullPtr;
   _item.commandSize() = 0;
   _item.childHead() = FileMem::nullPtr;
   _item.next() = FileMem::nullPtr;
}



void HistItem::copy_from(const HistItem& hist)
{
   bool cond = _item.addr()==FileMem::nullPtr;
   assert<IsAllocated>(cond,__FILE__,__LINE__);
   cond = hist.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   _mem.allot(_item);
   _fileName = hist.fileName();
   _object = hist.object();
   _command = hist.command();
   _item.timeStamp() = hist.timeStamp();
   _item.fileNamePtr() = set_string(_fileName);
   _item.fileNameSize() = _fileName.size()+1;
   _item.objectPtr() = set_string(_object);
   _item.objectSize() = _object.size()+1;
   _item.commandPtr() = set_string(_command);
   _item.commandSize() = _command.size()+1;
   _item.childHead() = rec_add_item(hist._mem,hist.childHead());
   _item.next() = rec_add_item(hist._mem,hist.next());
   _mem.sync(_item,FileSync::write);
}



void HistItem::sync()
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   _mem.sync(_item,FileSync::write);
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
   _item.fileNamePtr() = set_string(str);
   _item.fileNameSize() = str.size()+1;
}



std::string HistItem::fileName() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return _fileName;
}



void HistItem::object(const std::string& str)
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   cond = _item.objectPtr()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _object = str;
   _item.objectPtr() = set_string(str);
   _item.objectSize() = str.size()+1;
}



std::string HistItem::object() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return _object;
}



void HistItem::command(const std::string& str)
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   cond = _item.commandPtr()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _command = str;
   _item.commandPtr() = set_string(str);
   _item.commandSize() = str.size()+1;
}



std::string HistItem::command() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return _command;
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



FileMem::Ptr HistItem::addr() const
{
   return _item.addr();
}



void HistItem::operator=(FileMem::Ptr ptr)
{
   _item = ptr;
   _fileName.clear();
   _object.clear();
   _command.clear();
   if (_item.addr()!=FileMem::nullPtr)
   {
      load_item();
   }
}



inline void HistItem::load_item()
{
   _mem.sync(_item,FileSync::read);
   if (_item.fileNamePtr()!=FileMem::nullPtr)
   {
      _fileName = get_string(_item.fileNamePtr(),_item.fileNameSize());
   }
   if (_item.objectPtr()!=FileMem::nullPtr)
   {
      _object = get_string(_item.objectPtr(),_item.objectSize());
   }
   if (_item.commandPtr()!=FileMem::nullPtr)
   {
      _command = get_string(_item.commandPtr(),_item.commandSize());
   }
}



inline std::string HistItem::get_string(FileMem::Ptr ptr, FileMem::SizeT size)
{
   String str(size,ptr);
   _mem.sync(str,FileSync::read);
   bool cond = str.c_str()[size-1]=='\0';
   assert<InvalidItem>(cond,__FILE__,__LINE__);
   return {str.c_str()};
}



inline FileMem::Ptr HistItem::set_string(const std::string& newStr)
{
   String str(newStr.size()+1);
   _mem.allot(str);
   memcpy(str.c_str(),newStr.c_str(),newStr.size()+1);
   _mem.sync(str,FileSync::write);
   return str.addr();
}



FileMem::Ptr HistItem::rec_add_item(FileMem& mem, FileMem::Ptr ptr)
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
