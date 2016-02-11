#include "history.h"



History::History(const std::string& fileName):
   _mem(fileName)
{
   if (_mem.size()==0)
   {
      _mem.reserve(_emptyLen);
      _mem.allocate(_header);
      _mem.allocate(_info);
      _header.head() = _info.addr();
      _header.dataHead() = FileMem::nullPtr;
      _mem.sync(_header,FileSync::write);
      _info.timeStamp() = 0;
      _info.fileNamePtr() = FileMem::nullPtr;
      _info.fileNameSize() = 0;
      _info.objectPtr() = FileMem::nullPtr;
      _info.objectSize() = 0;
      _info.commandPtr() = FileMem::nullPtr;
      _info.commandSize() = 0;
      _info.childHead() = FileMem::nullPtr;
      _info.next() = FileMem::nullPtr;
      _mem.sync(_info,FileSync::write);
   }
   else
   {
      bool cond = _mem.size()>=_emptyLen;
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      _header = _mem.head();
      _mem.sync(_header,FileSync::read);
      cond = true;
      for (int i=0;i<_idLen;++i)
      {
         if (_header.ident()[i]!=_identString[i])
         {
            cond = false;
         }
      }
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      _info = _header.head();
      _mem.sync(_info,FileSync::read);
   }
}



void History::set_timestamp(int64_t tStamp)
{
   _info.timeStamp() = tStamp;
   _mem.sync(_info,FileSync::write);
}



void History::set_filename(const std::string& name)
{
   bool cond = _info.fileNamePtr()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _info.fileNamePtr() = write_string(name);
   _info.fileNameSize() = name.size();
   _mem.sync(_info,FileSync::write);
}



void History::set_object(const std::string& object)
{
   bool cond = _info.objectPtr()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _info.objectPtr() = write_string(object);
   _info.objectSize() = object.size();
   _mem.sync(_info,FileSync::write);
}



void History::set_command(const std::string& command)
{
   bool cond = _info.commandPtr()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _info.commandPtr() = write_string(command);
   _info.commandSize() = command.size();
   _mem.sync(_info,FileSync::write);
}



void History::add_child(History& child)
{
   if (_info.childHead()==FileMem::nullPtr)
   {
      _info.childHead() = rec_add_child(child.begin());
      _mem.sync(_info,FileSync::write);
   }
   else
   {
      Node i = _info.childHead();
      _mem.sync(i,FileSync::read);
      while (i.next()!=FileMem::nullPtr)
      {
         i = i.next();
         _mem.sync(i,FileSync::read);
      }
      i.next() = rec_add_child(child.begin());
      _mem.sync(i,FileSync::write);
   }
}



History::Iterator History::begin()
{
   return {_info.childHead(),_mem};
}



History::Iterator History::end()
{
   return {FileMem::nullPtr,_mem};
}



FileMem::Ptr History::write_string(const std::string& str)
{
   int rLen = str.size() + 1;
   if (_mem.capacity()<rLen)
   {
      _mem.reserve(rLen - _mem.capacity());
   }
   String string(rLen);
   memcpy(string.c_str(),str.c_str(),rLen);
   _mem.allocate(string);
   _mem.sync(string,FileSync::write);
}



FileMem::Ptr History::rec_add_child(Iterator i)
{
   if (i!=end())
   {
      Node child = *i;
      _mem.allocate(child);
      copy_child(child,i);
      child.childHead() = rec_add_child(i.childHead());
      child.next() = rec_add_child(++i);
      _mem.sync(child,FileSync::write);
      return child.addr();
   }
   else
   {
      return FileMem::nullPtr;
   }
}



void History::copy_child(Node& child, Iterator& i)
{
   bool cond = child.fileNamePtr()!=FileMem::nullPtr&&
               child.objectPtr()!=FileMem::nullPtr&&
               child.commandPtr()!=FileMem::nullPtr;
   assert<NotSet>(cond,__FILE__,__LINE__);
   child.fileNamePtr() = write_string(i.filename());
   child.objectPtr() = write_string(i.object());
   child.commandPtr() = write_string(i.command());
}



int64_t History::Iterator::timestamp()
{
   return _info.timeStamp();
}



std::string History::Iterator::filename()
{
   std::string ret;
   if (_ptr!=FileMem::nullPtr&&_info.fileNamePtr()!=FileMem::nullPtr)
   {
      String string(_info.fileNameSize(),_info.fileNamePtr());
      _mem.sync(string,FileSync::read);
      ret = string.c_str();
   }
   return ret;
}



std::string History::Iterator::object()
{
   std::string ret;
   if (_ptr!=FileMem::nullPtr&&_info.objectPtr()!=FileMem::nullPtr)
   {
      String string(_info.objectSize(),_info.objectPtr());
      _mem.sync(string,FileSync::read);
      ret = string.c_str();
   }
   return ret;
}



std::string History::Iterator::command()
{
   std::string ret;
   if (_ptr!=FileMem::nullPtr&&_info.commandPtr()!=FileMem::nullPtr)
   {
      String string(_info.commandSize(),_info.commandPtr());
      _mem.sync(string,FileSync::read);
      ret = string.c_str();
   }
   return ret;
}



History::Iterator History::Iterator::childHead()
{
   return {_info.childHead(),_mem};
}



History::Iterator& History::Iterator::operator++()
{
   _ptr = _info.next();
   if (_ptr!=FileMem::nullPtr)
   {
      _info = _ptr;
      _mem.sync(_info,FileSync::read);
   }
   return *this;
}



bool History::Iterator::operator!=(const Iterator& cmp)
{
   return _ptr!=cmp._ptr;
}



History::Iterator::Iterator(FileMem::Ptr ptr, FileMem& mem):
   _ptr(ptr),
   _mem(mem)
{
   _info.next() = FileMem::nullPtr;
   if (_ptr!=FileMem::nullPtr)
   {
      _mem.sync(_info,FileSync::read);
   }
}



History::Node& History::Iterator::operator*()
{
   return _info;
}
