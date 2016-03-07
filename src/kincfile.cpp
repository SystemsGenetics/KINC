#include <cstring>
#include "kincfile.h"



KincFile::KincFile(const std::string& fileName):
   _mem(fileName),
   _ident(&_mem)
{
   if (_mem.size()==0)
   {
      _mem.allot(_hdr);
      _hist = hptr(new History(_mem));
      _hdr.histHead() = _hist->addr();
      _hdr.dataHead() = FileMem::nullPtr;
      _hdr.ident() = FileMem::nullPtr;
      memcpy(_hdr.idString(),_idString,_idSz);
      _mem.sync(_hdr,FileSync::write);
   }
   else
   {
      bool cond = _mem.size()>=_hdrSz;
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      _hdr = _mem.head();
      _mem.sync(_hdr,FileSync::read);
      cond = _hdr.histHead()!=FileMem::nullPtr;
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      cond = strncmp(_hdr.idString(),_idString,_idSz)==0;
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      _hist = hptr(new History(_mem,_hdr.histHead()));
      _ident.addr(_hdr.ident());
      _new = false;
   }
}



void KincFile::ident(const string& id)
{
   try
   {
      _ident = id;
   }
   catch (FString::AlreadySet)
   {
      throw AlreadySet(__FILE__,__LINE__);
   }
   _hdr.ident() = _ident.addr();
   _mem.sync(_hdr,FileSync::write);
}



void KincFile::head(FileMem::Ptr ptr)
{
   _hdr.dataHead() = ptr;
   _mem.sync(_hdr,FileSync::write);
}
